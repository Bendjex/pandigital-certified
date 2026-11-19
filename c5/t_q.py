"""consistency of the recursive queries: for random windows, every materialised NE (or V) enclosure
of a set meeting a window must be covered by an enclosure returned by the recursion for that window"""
import sys, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import paths, os
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qx, ivl
n, o, kind = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
q = qx.QX(n, o, os.path.join(paths.SRC, f"res_{n}_{o}.txt"), n, n, verbose=False)
S = tuple(range(1, n + 1))
rng = np.random.default_rng(7)
if kind == "N":
    full = q.N(S); q.tn = n - 1          # recursion for S
else:
    full = q.V(S); q.tv = n - 1
fin = np.isfinite(full.A)
c = full.A[fin][rng.integers(0, fin.sum(), 300)] + rng.normal(0, 1e-9, 300)
w = rng.choice([1e-12, 1e-6, 1e-3], 300)
lo, hi = c - w, c + w
q.idx.clear()
H = q.Q(S, kind, lo, hi)
miss = 0; tot = 0
for i in range(len(lo)):
    m = (full.B >= lo[i]) & (full.A <= hi[i]) & ~full.Z
    hw = H.w == i
    for a, b in zip(full.A[m], full.B[m]):
        tot += 1
        if not ((H.A[hw] <= a + 1e-12 * (1 + abs(a))) & (H.B[hw] >= b - 1e-12 * (1 + abs(b)))).any():
            # covered means: some returned enclosure contains the materialised one (up to rounding)
            if not ((H.A[hw] <= b) & (H.B[hw] >= a)).any():
                miss += 1
print(f"n={n} ord={o} kind={kind}: windows {len(lo)}, materialised matches {tot}, recursion hits {len(H)}, UNCOVERED {miss}")
print({k: (v[0], v[1], round(v[2], 2)) for k, v in q.stats.items()})
