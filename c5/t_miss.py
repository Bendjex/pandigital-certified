import sys, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import paths, os
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import measure, motorx
S = tuple("12345"); measure.build(S)
mx = motorx.MX(5, 0, os.path.join(paths.SRC, "res_5_0.txt"), 4)
T = (1,2,3,4,5)
a = measure.A[S]; ne = a[~np.isin(measure.key(a), measure.key(measure.R[S]))]
e = mx.N(T)
order = np.argsort(e.A); A, B = e.A[order], e.B[order]; cm = np.maximum.accumulate(B)
miss = []
for v in ne[np.isfinite(ne)]:
    t = 1e-9 * (1 + abs(v)); i = np.searchsorted(A, v + t, side="right")
    if i == 0 or cm[i - 1] < v - t: miss.append(v)
miss = np.array(miss)
print(len(miss), "missing; sample", miss[:10])
print("abs<1e-6:", (np.abs(miss) < 1e-6).sum(), " near +-1e18:", (np.abs(miss) > 1e17).sum())
