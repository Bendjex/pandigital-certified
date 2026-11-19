import sys, numpy as np, itertools
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import measure
S = tuple(sys.argv[2]) if len(sys.argv) > 2 else tuple("12345"); measure.build(S)
target = float(sys.argv[1])
F = ["x*y", "x/y", "y/x", "x+y", "|x-y|", "x^y", "x^-y", "y^x", "y^-x"]
for r in (1, 2):
    for S1 in itertools.combinations(S, r):
        S2 = tuple(d for d in S if d not in S1)
        a, b = measure.A[S1], measure.A[S2]
        res = measure.combine(a, b)          # order: forms x 11? see combine
        n = len(a) * len(b)
        nf = len(res) // n if n else 0
        idx = np.where(np.abs(res - target) < 1e-9 * (1 + abs(target)))[0]
        for k in idx[:3]:
            f, rem = divmod(k, n); i, j = divmod(rem, len(b))
            print(S1, S2, "form#", f, "x L =", a[i], "y L =", b[j], "->", res[k])
