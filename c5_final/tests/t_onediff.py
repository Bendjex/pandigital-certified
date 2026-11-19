"""(16/11) qx.one_diff_c2: the second level of l for |s - l| = u with s exactly 1 (split 20 of n = 9).
For random windows [al, be] of ln u (be <= -1, from -1 down to -1e7, narrow and wide) and random u inside them
(ln u uniform in the window), the true l = 1 + u (branch 0) and l = 1 - u (branch 1) must satisfy the constraint
of their branch: sign of ln l, and ln|ln l| in [lower, upper] (mpmath, 60 digits on the ln|ln l| scale).  Also: the
rule gives no constraint where it does not apply (s not exactly 1, s = 0, a window with u = 0 (al = -inf), be > -1).
usage: t_onediff.py [trials] [seed]"""
import sys, os
import numpy as np, mpmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qx, ivl
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 7)
mpmath.mp.dps = 60


class X7:
    def __init__(self, n, O, Z):
        self.A = np.zeros(n); self.B = np.zeros(n); self.O = O; self.Z = Z


bad = ok = 0
for t in range(N):
    k = rng.integers(0, 4)
    c = -10.0 ** rng.uniform(0, 7) if k else -rng.uniform(1, 5)        # centre of ln u
    w = abs(c) * 10.0 ** rng.uniform(-16, -1) if rng.random() < 0.7 else rng.uniform(0, 3)
    al, be = c - w, min(c + w, -1.0)
    if al > be:
        continue
    x = X7(1, np.array([True]), np.array([False]))
    r_ = rng.random()
    lnu = mpmath.mpf(al) if r_ < 0.15 else (mpmath.mpf(be) if r_ < 0.3 else
          mpmath.mpf(al) + (mpmath.mpf(be) - mpmath.mpf(al)) * mpmath.mpf(rng.random()))   # the ends too
    u = mpmath.exp(lnu)
    for bi in (0, 1):
        sg, m2, M2 = qx.one_diff_c2(bi, x, np.array([al]), np.array([be]))
        if bi == 0:
            lnl = mpmath.log1p(u); want = 1
        else:
            lnl = mpmath.log1p(-u); want = -1
        t2 = mpmath.log(abs(lnl))
        tol = mpmath.mpf(10) ** -40 * max(1, abs(t2))
        good = sg[0] == want and mpmath.mpf(m2[0]) - tol <= t2 <= mpmath.mpf(M2[0]) + tol
        if not good:
            bad += 1
            print("MISS", bi, al, be, mpmath.nstr(lnu, 20), int(sg[0]), m2[0], M2[0], mpmath.nstr(t2, 25))
        else:
            ok += 1
# where the rule must not apply
cases = [(np.array([False]), np.array([False]), -100.0, -99.0, "s not 1"),
         (np.array([True]), np.array([True]), -100.0, -99.0, "s = 0"),
         (np.array([True]), np.array([False]), -np.inf, -99.0, "u = 0 admitted"),
         (np.array([True]), np.array([False]), -3.0, -0.5, "u up to e^-0.5")]
for O, Z, al, be, name in cases:
    r = qx.one_diff_c2(0, X7(1, O, Z), np.array([al]), np.array([be]))
    if r is not None and r[0][0] != 0:
        bad += 1; print("APPLIED WHERE IT MUST NOT:", name)
print(f"t_onediff: {ok} branch constraints contain the true l; PROBLEMS {bad}")
sys.exit(1 if bad else 0)
