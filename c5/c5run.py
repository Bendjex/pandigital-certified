"""C5 by top-level split, with the independent check of the hits in the same pass (03/11, "1'").
For every top-level split i of S = {1..n} with i % J == k and every constant: motor X's query of that split in provenance
mode (the computation of QX.split_hits, the same as C5's; Theorem X: every value of NE(S) in the window produced
by this split is covered by a returned hit) and then every value of every hit evaluated by lvl (hitcheck.check_hits).
The split is closed when its hits leave 0 undecided values.  One JSON line per (constant, split) is appended to the
output (flushed and synced): a relaunch skips the pairs already there, so the run is resumable split by split.
With C5RUN_CHUNKS=1 (10/11) a split is taken piece by piece (plan), the constants together, and every piece is
recorded when done: a restart repeats only the piece that was running.
(11/11) With E11_<const> (the optimum of the variant without unary minus, C11), the constant is queried with the wider
window, radius max(E*, E11), and each value is checked against both: "undecided" for E* (C7), "undecided_c11" for E11.
usage: c5run.py n ord consts k J out.jsonl        (E_<const> for constants without a certified record)"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hitcheck                                # its import sets the provenance tables (Ent.uniq keeps every derivation)
import paths, qx, cells_qpar
from cells_x import window
from flint import arb


CHUNKS = os.environ.get("C5RUN_CHUNKS") == "1"
C11 = {}          # (11/11) const -> radius E11 of the variant without unary minus (env E11_<const>); its window is
EQ = {}           # checked in the same run.  EQ: the radius each constant is queried with, max(E*, E11) (x (1 + 1e-9))
REDO = os.environ.get("C5RUN_REDO_UNDECIDED") == "1"    # repeat the pieces and splits recorded with undecided values


def plan(q, S1, S2, deep=None):
    """(10/11) the pieces of the top-level split S1 | S2, each a list of paths (part, only): part is one of the three
    parts of QX.split_hits (Lemma X1) and only = {set: pred(split number, part)} restricts the N queries of the sets
    below (QX.Q).  The pieces partition the parts of the split and, inside its part VN (V(S1), NE(S2)) when NE(S2) is
    a recursion, the parts of the splits of S2; for |S2| = 8 also the splits of its 7-digit subsets reached through a
    single digit.  So the hits of the pieces add up to those of QX.split_hits(S1, S2).
      |S2| <= 5 (tables): NU, VN, VV.
      |S2| = 6, 7: NU, VV, and inside VN one piece per split x | rest of S2 with a single digit x, and one per size
        of the table side for the others (r2, r3).
      |S2| = 8 (deep; any order): NU, VV, and inside VN
        T<digits>, one per subset T of S2 = T u {x, y}: everything that queries NE(T) as the query side:
          (V(xy), NE(T)), and (V(x), NE(S2 - x)) restricted below to (V(y), NE(T)), and the same with x, y exchanged;
        <x>, one per digit: the split x | S2 - x without the pairs (V(y), NE(T)) of the pieces T;
        r2: the parts NU and VV of the splits with a two-digit table side;  r3, r4: the splits of table side 3, 4.
    A piece queries each NE(T) with all its windows together, once per path."""
    D = lambda T: "".join(map(str, T))
    pieces = [("NU", [("NU", None)]), ("VV", [("VV", None)])]
    if q.tables_ok(S2, "N"):
        return pieces[:1] + [("VN", [("VN", None)])] + pieces[1:]
    sides = q.table_side(S2)
    size = [len(t[0]) for t in sides]
    num = {t[0]: i for i, t in enumerate(sides)}
    if deep is None:
        deep = len(S2) == 8 and not q.ord
    if not deep:
        for i, t in enumerate(sides):
            if size[i] == 1:
                pieces.append((f"VN.{D(t[0])}", [("VN", {S2: (lambda si, part, i=i: si == i)})]))
        for k in sorted(set(size) - {1}):
            pieces.append((f"VN.r{k}", [("VN", {S2: (lambda si, part, k=k: size[si] == k)})]))
        return pieces
    assert not q.ord and len(S2) >= 7
    import itertools
    sub = {}                                           # the 7-digit subsets: their splits with a single digit
    for x in S2:
        U = tuple(d for d in S2 if d != x)
        su = q.table_side(U)
        sub[x] = (U, {t[0][0]: i for i, t in enumerate(su) if len(t[0]) == 1}, [len(t[0]) for t in su])
    for T in itertools.combinations(S2, len(S2) - 2):
        x, y = [d for d in S2 if d not in T]
        paths = [("VN", {S2: (lambda si, part, i=num[(x, y)]: si == i and part == "VN")})]
        for a, b in ((x, y), (y, x)):
            U, one, _ = sub[a]
            paths.append(("VN", {S2: (lambda si, part, i=num[(a,)]: si == i and part == "VN"),
                                 U: (lambda si, part, i=one[b]: si == i and part == "VN")}))
        pieces.append((f"VN.T{D(T)}", paths))
    for x in S2:
        U, one, sz = sub[x]
        pieces.append((f"VN.{x}", [("VN", {S2: (lambda si, part, i=num[(x,)]: si == i),
                                           U: (lambda si, part, sz=sz: not (sz[si] == 1 and part == "VN"))})]))
    pieces.append(("VN.r2", [("VN", {S2: (lambda si, part: size[si] == 2 and part != "VN")})]))
    for k in sorted(set(size) - {1, 2}):
        pieces.append((f"VN.r{k}", [("VN", {S2: (lambda si, part, k=k: size[si] == k)})]))
    return pieces


def path_hits(q, S1, S2, forms, part, only, lo, hi):
    """the hits of one path of a piece: the part of QX.split_hits (the same calls, top level: c2 None) with the N
    queries below restricted by only"""
    q.only = only
    try:
        if part == "NU":
            if len(S1) >= 3 or q.tables_ok(S1, "N"):
                return q.combine(S1, "N", S2, "U", forms, lo, hi, want="any", c2=None)
            return qx.Hits.cat([])
        if part == "VN":
            return q.combine(S1, "V", S2, "N", forms, lo, hi, want="any", c2=None)
        assert qx.CLIP_OUT
        L, H, idx = qx.out_parts(lo, hi)
        if not len(L):
            return qx.Hits.cat([])
        h = q.combine(S1, "V", S2, "V", forms, L, H, want="out", c2=None)
        if len(h):
            h.w = idx[h.w]
        return h
    finally:
        q.only = None


def chunked_split(q, C, consts, i, side, wins, vws, out, done_chunks):
    """(10/11) one split for the constants consts (their windows together), piece by piece (plan): each piece is
    queried and checked on its own and recorded when done, one record per constant, so a restart repeats only the
    piece that was running; returns the totals per constant"""
    S1, S2, forms = side
    pieces = plan(q, S1, S2)
    tot = {c: [0, 0, [], 0, 0] for c in consts}
    for n_, (name, paths) in enumerate(pieces):
        todo = [c for c in consts if (c, i, name) not in done_chunks]
        if todo:
            lo = np.array([wins[c][0] for c in todo]); hi = np.array([wins[c][1] for c in todo])
            qx.PROV = qx._Prov()
            t1, c1 = time.time(), time.process_time()
            H = qx.Hits.cat([path_hits(q, S1, S2, forms, part, only, lo, hi) for (part, only) in paths])
            tq = time.process_time() - c1
            for k, c in enumerate(todo):
                Hc = H.sel(H.w == k)
                tot_k, det = hitcheck.check_hits_multi(C, Hc, c, vws[c])
                r = {"const": c, "split": i, "chunk": name, "hits": len(Hc), "undecided": tot_k["c7"][0],
                     "cpu": round(time.process_time() - c1, 1), "cpu_query": round(tq, 1), "seconds": round(time.time() - t1, 1),
                     "together": todo, "Eq": EQ[c], "undecided_detail": [d for d in det if d["undecided"]][:20]}
                if "c11" in tot_k:
                    r["undecided_c11"], r["inside_c11"] = tot_k["c11"]
                with open(out, "a", encoding="utf-8") as fo:
                    fo.write(json.dumps(r) + "\n"); fo.flush(); os.fsync(fo.fileno())
                done_chunks[(c, i, name)] = r
            del H
            if len(C.memo) > 200000:
                C.memo.clear()
        for c in consts:
            r = done_chunks[(c, i, name)]
            tot[c][0] += r["hits"]; tot[c][1] += r["undecided"]; tot[c][2] += r.get("undecided_detail", [])
            tot[c][3] += r.get("undecided_c11", 0); tot[c][4] += r.get("inside_c11", 0)
        r = done_chunks[(consts[0], i, name)]
        print(f"[chunks] split {i} ({','.join(consts)}): {n_ + 1}/{len(pieces)} chunks done ({name}: hits "
              f"{'/'.join(str(done_chunks[(c, i, name)]['hits']) for c in consts)}, undecided "
              f"{'/'.join(str(done_chunks[(c, i, name)]['undecided']) for c in consts)}"
              + (f", C11 {'/'.join(str(done_chunks[(c, i, name)].get('undecided_c11', '-')) for c in consts)}" if C11 else "")
              + f", {r['cpu']} s)", flush=True)
    return tot


def main():
    n, o = int(sys.argv[1]), int(sys.argv[2])
    consts = sys.argv[3].split(",")
    k, J = int(sys.argv[4]), int(sys.argv[5])
    out = sys.argv[6]
    done = set()
    done_chunks = {}
    split_eq = {}                                  # (const, split) -> query radius of its split record
    if os.path.exists(out) and os.path.getsize(out) > 0:
        with open(out, "rb") as f:                # a line cut by a power failure: close it, so that the next record
            f.seek(-1, 2)                          # starts on a line of its own (the cut line is ignored below)
            cut = f.read(1) != b"\n"
        if cut:
            with open(out, "a", encoding="utf-8") as f:
                f.write("\n"); f.flush(); os.fsync(f.fileno())
    if os.path.exists(out):
        for l in open(out, encoding="utf-8"):
            try:
                r = json.loads(l)
            except ValueError:
                continue                      # a line cut by a power failure: that split is simply done again
            if REDO and (r.get("undecided") or r.get("undecided_c11")):   # undecided values: made again (later lines win)
                done_chunks.pop((r["const"], r["split"], r.get("chunk")), None); done.discard((r["const"], r["split"]))
                continue
            if "chunk" in r:
                done_chunks[(r["const"], r["split"], r["chunk"])] = r
            elif "split" in r:
                done.add((r["const"], r["split"])); split_eq[(r["const"], r["split"])] = r.get("Eq")
    cert = {}
    for l in open(os.path.join(paths.RESULTS, "certified.jsonl")):
        r = json.loads(l); cert[(r["const"], r["n"], r["ord"])] = r
    qx.DEBUG = True
    tv, tn = (5, 5) if o else (4, 5)
    q = qx.QX(n, o, os.path.join(os.environ.get("RESDIR", paths.SRC), f"res_{n}_{o}.txt"), tv, tn)
    S = tuple(range(1, n + 1))
    sides = q.table_side(S)
    mine = [i for i in range(len(sides)) if i % J == k]
    C = hitcheck.Checker(q)
    wins, vws = {}, {}
    for c in consts:
        r = cert.get((c, n, o))
        E = (float(arb(r["optimal"][0]["error"]).upper()) if r else float(os.environ["E_" + c])) * (1 + 1e-9)
        vws[c] = {"c7": hitcheck.value_window(c, E)}
        Eq = E
        if os.environ.get("E11_" + c):
            # (11/11) the window of the variant without unary minus (C11), radius E11 >= E*: the query is made with the
            # wider of the two (its hits contain those of the narrower), and every value is tested against both
            E11 = float(os.environ["E11_" + c]) * (1 + 1e-9)
            # (11/11) E11 may be below the radius used for C7 when that radius is only an upper bound of E* (n = 9 in any
            # order: the windows of C7; gamma: C11 certified 1.98e-15 < 2.12e-15): the query keeps the wider radius
            C11[c] = E11
            vws[c]["c11"] = hitcheck.value_window(c, E11)
            Eq = max(E, E11)
        EQ[c] = Eq
        wins[c] = window(c, Eq)
    # records made with another query radius (another E11) are not used: their pieces are made again
    # (records of before 11/11 carry no radius: they count as made with E*, valid only for a constant without E11)
    plain = lambda c: None if c in C11 else EQ.get(c)
    for key in [k_ for k_, r_ in done_chunks.items() if r_.get("Eq", plain(k_[0])) != EQ.get(k_[0])]:
        del done_chunks[key]
    for key in [k_ for k_ in list(done) if (split_eq.get(k_) if split_eq.get(k_) is not None else plain(k_[0])) != EQ.get(k_[0])]:
        done.discard(key)

    def record(c, i, nh, und, det, t0, c0, tq, und11=None, ins11=None):
        rec = {"const": c, "n": n, "ord": o, "split": i, "S1": "".join(map(str, sides[i][0])), "hits": nh,
               "undecided": und, "cpu_query": round(tq, 1), "cpu": round(time.process_time() - c0, 1),
               "seconds": round(time.time() - t0, 1), "motorB_unresolved": q.motorb_unresolved(), "Eq": EQ[c],
               "undecided_detail": [d for d in det if d["undecided"]][:50]}
        if c in C11:
            rec["E11"] = C11[c]; rec["undecided_c11"] = und11; rec["inside_c11"] = ins11
        with open(out, "a", encoding="utf-8") as f:
            f.write(json.dumps(rec) + "\n"); f.flush(); os.fsync(f.fileno())
        print(f"RESULT {c} split {i} S1={rec['S1']}: hits {nh}, undecided {und}"
              + (f", C11 undecided {und11} (certainly inside {ins11})" if c in C11 else "") + f", cpu {rec['cpu']}s", flush=True)

    if CHUNKS:
        # (10/11) the constants together: their windows go down the recursion in the same batches
        for i in mine:
            todo = [c for c in consts if (c, i) not in done]
            if not todo:
                continue
            t0, c0 = time.time(), time.process_time()
            q.top = S
            tot = chunked_split(q, C, todo, i, sides[i], wins, vws, out, done_chunks)
            for c in todo:
                record(c, i, tot[c][0], tot[c][1], tot[c][2], t0, c0, 0.0, tot[c][3], tot[c][4])
    else:
        for c in consts:
            wlo, whi = wins[c]
            for i in mine:
                if (c, i) in done:
                    continue
                t0, c0 = time.time(), time.process_time()
                q.top = S
                qx.PROV = qx._Prov()               # provenance of this split only
                H = cells_qpar.part_hits(q, S, np.array([wlo]), np.array([whi]), 0, 1, select=lambda j, i=i: j == i, label=c)
                tq = time.process_time() - c0
                tot_k, det = hitcheck.check_hits_multi(C, H, c, vws[c])
                record(c, i, len(H), tot_k["c7"][0], det, t0, c0, tq, *(tot_k["c11"] if "c11" in tot_k else (None, None)))
                del H
                if len(C.memo) > 200000:
                    C.memo.clear()
    print("DONE", n, o, ",".join(consts), f"part {k}/{J}", flush=True)


if __name__ == "__main__":
    main()
