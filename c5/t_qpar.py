import sys, numpy as np, time
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import paths, os
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qx, ivl, cells_qpar
cells_qpar.T0 = time.time()
for (n, o, tv, tn) in [(6, 1, 5, 5), (6, 0, 4, 5)]:
    q = qx.QX(n, o, os.path.join(paths.SRC, f"res_{n}_{o}.txt"), tv, tn)
    S = tuple(range(1, n + 1))
    lo, hi = 1 - 1e-9, 1 + 1e-9
    wl = np.array([float(ivl.ln_lo(np.array([lo]))[0])]); wh = np.array([float(ivl.ln_hi(np.array([hi]))[0])])
    full = q.Q(S, "N", wl, wh)
    key = lambda H: sorted(zip(H.A.tolist(), H.B.tolist()))
    for J in (2, 3):
        parts = [cells_qpar.part_hits(q, S, wl, wh, k, J) for k in range(J)]
        u = qx.Hits.cat(parts)
        print(n, o, "full", len(full), "J", J, "union", len(u), "same multiset:", key(full) == key(u))
