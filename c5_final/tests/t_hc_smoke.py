import sys, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import paths, os
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hitcheck, qx, ivl, lvl
qx.DEBUG = True
q = qx.QX(5, 1, os.path.join(paths.SRC, "res_5_1.txt"), 5, 5)
lo, hi = 1 - 1e-9, 1 + 1e-9
wl = float(ivl.ln_lo(np.array([lo]))[0]); wh = float(ivl.ln_hi(np.array([hi]))[0])
H = q.Q((1, 2, 3, 4, 5), "N", np.array([wl]), np.array([wh]))
print("hits", len(H))
C = hitcheck.Checker(q)
tot = ins = und = 0
for i in range(len(H)):
    for v, e in C.hit(int(H.pid[i])):
        tot += 1
        if v is None: und += 1; print("FAIL", e[:200]); continue
        if not lvl.outside(v, lo, hi): ins += 1
        if tot <= 4: print(repr(v), e[:160])
print("values", tot, "inside-or-undecided", ins, "fail", und)
