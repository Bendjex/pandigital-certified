"""coverage test: every approximate NE value of measure.py (A \\ R, float) lies in some enclosure of
motor X's Nt (tolerance 1e-9 (1 + |L|) for the float approximations of measure.py)"""
import sys, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import paths, os
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import measure, motorx
n, o = int(sys.argv[1]), int(sys.argv[2])
assert o == 0
S = tuple(str(d) for d in range(1, n + 1))
measure.build(S)
mx = motorx.MX(n, o, os.path.join(paths.SRC, f"res_{n}_{o}.txt"), n)
for k in range(3, n + 1):
    for T in mx.sets(k):
        Ts = tuple(str(d) for d in T)
        a = measure.A[Ts]
        ne = a[~np.isin(measure.key(a), measure.key(measure.R[Ts]))]
        e = mx.N(T)
        # sort enclosures; for each approx value find a covering enclosure
        tol = 1e-9 * (1 + np.abs(ne))
        A, B = e.A, e.B
        order = np.argsort(A); A, B = A[order], B[order]
        cummaxB = np.maximum.accumulate(B)
        miss = 0
        fin = np.isfinite(ne)
        for v, t in zip(ne[fin], tol[fin]):
            i = np.searchsorted(A, v + t, side="right")      # entries with A <= v + t
            if i == 0 or cummaxB[i - 1] < v - t:
                miss += 1
        zeros = (~fin).sum()
        print("".join(map(str, T)), "approx NE", len(ne), "enclosures", len(e), "MISSING", miss, "(zeros/inf skipped", zeros, ")")
