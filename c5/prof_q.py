import sys, numpy as np, time
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import paths, os
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qx
from cells_x import window
n, o, c, E, tv, tn = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], float(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
q = qx.QX(n, o, os.path.join(paths.SRC, f"res_{n}_{o}.txt"), tv, tn)
wlo, whi = window(c, E * (1 + 1e-9))
t0 = time.time()
H = q.Q(tuple(range(1, n + 1)), "N", np.array([wlo]), np.array([whi]))
print("hits", len(H), "time", round(time.time() - t0, 1))
for k, v in sorted(q.cstats.items(), key=lambda t: -t[1][1])[:25]:
    print(k, v)
