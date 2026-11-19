import sys, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ivl
rng = np.random.default_rng(1)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
for f in range(9):
    Ls = rng.choice([-1e3, -50, -5, -0.5, -1e-6, 0, 1e-6, .5, 5, 50, 700, 1e4], N) * rng.random(N)
    Ll = rng.choice([-1e3, -50, -5, -0.5, -1e-6, 0, 1e-6, .5, 5, 50, 700, 1e4], N) * rng.random(N)
    z1 = rng.random(N) < .02; z2 = rng.random(N) < .02
    o1 = (Ls == 0) & ~z1; o2 = (Ll == 0) & ~z2
    Ls = np.where(z1, -np.inf, Ls); Ll = np.where(z2, -np.inf, Ll)
    r, ok = ivl.fwd(f, ivl.Q7(Ls, Ls, z1, o1), ivl.Q7(Ll, Ll, z2, o2))
    A, B, Z = r.A, r.B, r.Z
    c = rng.choice([-1e3, -10, -1, 0, 1, 10, 1e3], N) * rng.random(N)
    w = rng.choice([1e-12, 1e-6, 1e-2, 1, 10], N)
    al, be = c - w, c + w
    al = np.where(rng.random(N) < .05, -np.inf, al)
    meets = ok & ((Z & (al == -np.inf)) | (~Z & (B >= al) & (A <= be)))
    wins = ivl.pre(f, Ls, Ls, z1, o1, al, be)
    inside = np.zeros(N, bool)
    for lo, hi in wins:
        inside |= (Ll >= lo) & (Ll <= hi) | (z2 & (lo == -np.inf))
    bad = meets & ~inside
    print("form", f, ivl.FORMS[f], "meets", meets.sum(), "MISSED", bad.sum())
