"""
C5 for the large cells, split into independent processes (no fork: also for Windows).
Part k of J handles the top-level splits (S1, S2) of {1..n} with index = k (mod J), exactly as QX.Q does
for kind "N" (qx.py), and writes its hits. The hits of the cell are the union over k = 0..J-1:
the cell is closed when every part reports 0 hits (or its hits are cleared by hitcheck.py).
usage: cells_qpar.py n ord tV tN const[,const...] k J [outfile]
"""
import sys, json, time, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qx, ivl
from cells_x import window
from flint import arb

import paths
RESDIR = os.environ.get("RESDIR", paths.SRC)
CERT = os.environ.get("CERT", os.path.join(paths.RESULTS, "certified.jsonl"))


LAST_CPU = [0.0]          # CPU seconds of the last split (process time; the wall time is in 'seconds')


def part_hits(q, S, lo, hi, k, J, select=None, on_split=None, label=""):
    """hits of the top-level splits idx with idx % J == k (or select(idx) True), as QX.Q for kind "N";
    on_split(idx, hits, seconds) is called after each split (checkpoint)"""
    res = []
    sides = q.table_side(S)
    mine = [i for i in range(len(sides)) if (select(i) if select else i % J == k)]
    t00 = time.time()
    for c_, idx in enumerate(mine):
        S1, S2, forms = sides[idx]
        print(f"[progress] {label} part {k}/{J}: split {c_}/{len(mine)} ({100 * c_ / max(1, len(mine)):.1f}%), "
              f"{time.time() - t00:.0f}s elapsed", flush=True)
        t1 = time.time(); c1 = time.process_time()
        h = qx.Hits.cat(q.split_hits(S1, S2, forms, lo, hi))
        LAST_CPU[0] = time.process_time() - c1
        res.append(h)
        if on_split:
            on_split(idx, h, time.time() - t1)
    print(f"[progress] {label} part {k}/{J}: split {len(mine)}/{len(mine)} (100.0%), {time.time() - t00:.0f}s elapsed", flush=True)
    return qx.Hits.cat(res)


if __name__ == "__main__":
    n, o, tv, tn = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    consts = sys.argv[5].split(",")
    k, J = int(sys.argv[6]), int(sys.argv[7])
    outf = sys.argv[8] if len(sys.argv) > 8 else f"cells_qpar_{n}_{o}_{k}of{J}.jsonl"
    os.makedirs(os.path.dirname(os.path.abspath(outf)), exist_ok=True)
    ck = outf + ".splits.jsonl"                      # checkpoint: one line per finished split
    done = {}
    if os.path.exists(ck):
        for l in open(ck):
            try:
                r = json.loads(l)
            except ValueError:
                break
            done[(r["const"], r["split"])] = r
    cert = {}
    for l in open(CERT):
        r = json.loads(l); cert[(r["const"], r["n"], r["ord"])] = r
    q = qx.QX(n, o, os.path.join(RESDIR, f"res_{n}_{o}.txt"), tv, tn)
    S = tuple(range(1, n + 1))
    nsplits = len(q.table_side(S))
    for c in consts:
        r = cert.get((c, n, o))
        E = float(arb(r["optimal"][0]["error"]).upper()) if r else float(os.environ["E_" + c])
        wlo, whi = window(c, E * (1 + 1e-9))
        T0 = time.time()

        def on_split(idx, h, sec, c=c):
            rec = {"const": c, "split": idx, "hits": len(h), "seconds": round(sec, 1), "cpu": round(LAST_CPU[0], 1),
                   "motorb_unresolved": q.motorb_unresolved()}
            if len(h):
                rec["sample"] = [[float(a), float(b)] for a, b in zip(h.A[:5], h.B[:5])]
            with open(ck, "a") as f:
                f.write(json.dumps(rec) + "\n")
            done[(c, idx)] = rec

        only = set(int(x) for x in os.environ.get("QPAR_ONLY", "").split(",") if x)   # tests: only these splits
        sel = lambda i, c=c: i % J == k and (c, i) not in done and (not only or i in only)
        part_hits(q, S, np.array([wlo]), np.array([whi]), k, J, select=sel, on_split=on_split, label=c)
        mine = [i for i in range(nsplits) if i % J == k]
        recs = [done[(c, i)] for i in mine if (c, i) in done]
        res = {"const": c, "n": n, "ord": o, "part": f"{k}/{J}", "E": E, "splits": len(mine),
               "splits_done": len(recs), "hits": sum(x["hits"] for x in recs),
               "splits_with_hits": [x["split"] for x in recs if x["hits"]],
               "seconds": round(sum(x["seconds"] for x in recs), 1),
               "motorb_unresolved": max([x.get("motorb_unresolved", 0) for x in recs] + [q.motorb_unresolved()])}
        print(json.dumps(res), flush=True)
        with open(outf, "a") as f:
            f.write(json.dumps(res) + "\n")
