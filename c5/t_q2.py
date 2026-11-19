import sys, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import paths, os
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qx, ivl
n, o = 5, 0
q = qx.QX(n, o, os.path.join(paths.SRC, f"res_{n}_{o}.txt"), n, n, verbose=False)
S = tuple(range(1, n + 1))
rng = np.random.default_rng(7)
full = q.N(S); q.tn = n - 1
fin = np.isfinite(full.A)
c = full.A[fin][rng.integers(0, fin.sum(), 300)] + rng.normal(0, 1e-9, 300)
w = rng.choice([1e-12, 1e-6, 1e-3], 300)
lo, hi = c - w, c + w
H = q.Q(S, "N", lo, hi)
F = ivl.FORMS
for i in range(len(lo)):
    m = np.where((full.B >= lo[i]) & (full.A <= hi[i]) & ~full.Z)[0]
    hw = H.w == i
    for k in m:
        a, b = full.A[k], full.B[k]
        if not ((H.A[hw] <= b) & (H.B[hw] >= a)).any():
            p = full.prov[k]
            S1, S2, forms = q.splits(S)[p[0]]
            print("window", lo[i], hi[i], "uncovered", a, b, "prov split", S1, S2, F[p[1]], "I", p[2], "J", p[3], "kinds", p[4], p[5])
            e1 = q.N(S1) if p[4] == 1 else q.V(S1)
            print("  s:", e1.A[p[2]], e1.B[p[2]], e1.Z[p[2]], e1.O[p[2]])
            if p[5] == 2:
                V2 = q.V(S2); e2 = V2 if p[3] < len(V2) else q.N(S2); j = p[3] if p[3] < len(V2) else p[3] - len(V2)
            else:
                e2 = q.N(S2) if p[5] == 1 else q.V(S2); j = p[3]
            print("  l:", e2.A[j], e2.B[j], e2.Z[j], e2.O[j])
            wins = ivl.pre(p[1], np.array([e1.A[p[2]]]), np.array([e1.B[p[2]]]), np.array([e1.Z[p[2]]]), np.array([e1.O[p[2]]]), np.array([lo[i]]), np.array([hi[i]]))
            print("  pre windows:", [(x[0], y[0]) for x, y in wins])
