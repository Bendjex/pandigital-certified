"""Reviewer's end-to-end check of c5run.plan (the piece partition used for the chunked runs): with tn = 4 the NE
query of a 5-digit S2 is a recursion, so plan() cuts the split S1 | S2 into pieces; their hits together must equal
those of QX.split_hits (as multisets of window index + every field)."""
import sys, os, time
PC = sys.argv[1]
sys.path.insert(0, os.path.join(PC, "c5")); sys.path.insert(1, os.path.join(PC, "src"))
import numpy as np
import qx, ivl, c5run, motorx
n, o, tv, tn = 6, 0, 4, 4
q = qx.QX(n, o, os.path.join(PC, "src", f"res_{n}_{o}.txt"), tv, tn)
S = tuple(range(1, n + 1))
sides = q.table_side(S)
rng = np.random.default_rng(5)
def key(H):
    cols = [H.w.astype(np.float64)] + [np.asarray(getattr(H, k)).astype(np.float64) for k in ivl.Q7.FIELDS]
    M = np.stack(cols, 1)
    M = np.where(np.isnan(M), 1e999, M)
    rows = [tuple(r) for r in M]
    rows.sort()
    return rows
bad = 0
tested = 0
for si in [i for i, t in enumerate(sides) if len(t[1]) == 5][:int(sys.argv[2])]:
    S1, S2, forms = sides[si]
    # windows: one around pi, one wide in-range, one near 1, one near a huge value
    lo = np.array([np.log(3.1415926), -1e-3, np.log(2.0), 50.0])
    hi = np.array([np.log(3.1415927), 1e-3, np.log(2.0001), 50.01])
    t0 = time.time()
    whole = qx.Hits.cat(q.split_hits(S1, S2, forms, lo, hi))
    t1 = time.time()
    pieces = c5run.plan(q, S1, S2)
    parts = []
    for name, paths in pieces:
        for (part, only) in paths:
            parts.append(c5run.path_hits(q, S1, S2, forms, part, only, lo, hi))
    tog = qx.Hits.cat(parts)
    same = key(whole) == key(tog)
    tested += 1
    print(f"split {si} S1={S1} S2={S2}: pieces {[p[0] for p in pieces]}; whole {len(whole)} hits ({t1 - t0:.0f}s), pieces {len(tog)} hits ({time.time() - t1:.0f}s); equal multisets: {same}", flush=True)
    bad += not same
print("splits tested", tested, "PROBLEMS", bad)
