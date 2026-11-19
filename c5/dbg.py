import sys, json, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import paths, os
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qx, ivl
from cells_x import window
from flint import arb
qx.DEBUG = True
n, o, c, E = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], float(sys.argv[4])
q = qx.QX(n, o, os.path.join(paths.SRC, f"res_{n}_{o}.txt"), 5, 5)
wlo, whi = window(c, E * (1 + 1e-9))
H = q.Q(tuple(range(1, n + 1)), "N", np.array([wlo]), np.array([whi]))
print("hits", len(H))
for i in range(min(len(H), int(sys.argv[5]) if len(sys.argv) > 5 else 3)):
    print("HIT", H.A[i], H.B[i]); qx.explain(q, int(H.pid[i]))
