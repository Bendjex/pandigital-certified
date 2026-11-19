"""combine() with the recursive sub-queries batched (QX_BATCHREC, 10/11: the windows of every form and branch of a
chunk in at most two queries) against one query per form and branch, and with each distinct window answered once
(QX_DEDUP): for random windows on a set above the table sizes, some of them repeated, the three modes must return
exactly the same hits (window index and every field of the enclosure), as multisets.
Also (QX.only): the hits of the query restricted to each split of the set, or to each part (NU, VN, VV) of its splits,
put together, are those of the whole query.
Windows: narrow, in range, far out of range on both sides, near 1 and at simple values.
usage: t_batchrec.py n digits kind nwin [seed]      e.g. t_batchrec.py 9 234689 N 40"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, qx, ivl
if os.environ.get("QX_ALLPROV") == "1":          # the tables of the provenance mode (as c5run.py)
    import hitcheck
    qx.DEBUG = True

n = int(sys.argv[1]); S = tuple(int(c) for c in sys.argv[2]); kind = sys.argv[3]; NW = int(sys.argv[4])
rng = np.random.default_rng(int(sys.argv[5]) if len(sys.argv) > 5 else 5)
q = qx.QX(n, 0, os.path.join(os.environ.get("RESDIR", paths.SRC), f"res_{n}_0.txt"), 4, 5, verbose=False)
assert not q.tables_ok(S, kind)
k = rng.integers(0, 5, NW)
c = np.where(k == 0, rng.normal(0, 3, NW),
    np.where(k == 1, rng.choice([-1, 1], NW) * rng.uniform(700, 4000, NW),
    np.where(k == 2, rng.choice([-1, 1], NW) * 10.0 ** rng.uniform(-14, -3, NW),
    np.where(k == 3, np.log(rng.choice([0.5, 2, 3, 0.25, 1.5, 10, 0.1, 8, 0.8], NW)), rng.uniform(-690, 690, NW)))))
w = np.where(k == 0, 1e-9, np.where(k == 1, 10.0 ** rng.uniform(-9, -6, NW), np.where(k == 2, np.abs(c) * 1e-6, np.where(k == 3, 3e-14, 1e-9))))
lo, hi = c - w, c + w


def rows(H):
    """the hits as sorted rows of bytes: window index and every field of the enclosure"""
    if not len(H):
        return []
    cols = [np.ascontiguousarray(H.w).view(np.uint8).reshape(len(H), -1)]
    cols += [np.ascontiguousarray(getattr(H, f)).view(np.uint8).reshape(len(H), -1) for f in qx.Hits.F]
    return sorted(bytes(r) for r in np.concatenate(cols, 1))


# some windows twice (bit for bit), as the preimage windows of the recursion are
dup = rng.integers(0, NW, max(2, NW // 4))
lo = np.concatenate([lo, lo[dup]]); hi = np.concatenate([hi, hi[dup]])
res = {}
for mode in (False, True, "dedup"):
    qx.BATCHREC = bool(mode); qx.DEDUP = mode == "dedup"
    qx.PROV = qx._Prov()
    t = time.time()
    H = q.Q(S, kind, lo, hi)
    res[mode] = rows(H)
    print(f"QX_BATCHREC={int(bool(mode))} QX_DEDUP={int(mode == 'dedup')}: {len(H)} hits, {time.time() - t:.1f} s", flush=True)
bad = int(res[False] != res[True]) + int(res[False] != res["dedup"])
# the query split by split (QX.only), batched
bad2 = 0
if kind == "N":
    allr = []
    ns = len(q.table_side(S))
    for j in range(ns):
        q.only = {S: (lambda si, part, j=j: si == j)}
        qx.PROV = qx._Prov()
        allr += rows(q.Q(S, kind, lo, hi))
    for part in ("NU", "VN", "VV"):                  # and part by part
        q.only = {S: (lambda si, p_, part=part: p_ == part)}
        qx.PROV = qx._Prov()
        allp = (allp if part != "NU" else []) + rows(q.Q(S, kind, lo, hi))
    q.only = None
    bad2 = int(sorted(allr) != res[True]) + int(sorted(allp) != res[True])
    print(f"split by split ({ns} splits): {len(allr)} hits; part by part: {len(allp)} hits")
print(f"t_batchrec {kind}{S}, {len(lo)} windows: {len(res[True])} hits; the three modes differ: {bad}; "
      f"split by split differs: {bad2}; motor B unresolved {q.motorb_unresolved()}")
sys.exit(1 if bad or bad2 or q.motorb_unresolved() else 0)
