"""(17/11) lvl: integer powers of a huge exact power plus a constant (split 55 of n = 9: |3^(10^7) - (9^(10^6) +- c)^5|).
X = B^k with X of more than 600 000 bits (kept as its prime vector), q = X + s c, and T = X^m written in another way
(B0^(a k m) with B = B0^a, or B^(k m)).  |T - q^m| = |sum_{j >= 1} C(m, j) (s c)^j X^(m - j)|, whose ln is
ln(m c) + (m - 1) k ln B up to a relative c/X (below 10^-100000 here): the enclosure of lvl must contain it (mpmath,
300 digits) and the value must be decided outside a window at pi.  Also q^m against itself written twice: exact 0.
usage: t_lvlbin.py [trials] [seed]"""
import sys, os
import numpy as np, mpmath
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lvl
N = int(sys.argv[1]) if len(sys.argv) > 1 else 24
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
mpmath.mp.dps = 300
E = lvl.exact
BASES = [(3, 2), (2, 2), (2, 3), (5, 2), (3, 3), (7, 1), (3, 1)]    # B = B0^a
bad = ok = zeros = 0
for t in range(N):
    B0, a = BASES[int(rng.integers(len(BASES)))]
    B = B0 ** a
    k = int(rng.integers(int(620000 / (a * np.log2(B0))) + 1, int(1500000 / (a * np.log2(B0)))))
    m = int(rng.integers(2, 7))
    c = Fr(int(rng.integers(1, 40)), int(rng.choice([1, 5, 10])))
    s = 1 if rng.random() < 0.5 else -1
    try:
        X = lvl.pw(E(Fr(B)), E(k), 1)
        q = lvl.add(E(c), X) if s > 0 else lvl.sub(X, E(c))
        p = lvl.pw(q, E(m), 1)
        T = lvl.pw(E(Fr(B0)), E(a * k * m), 1) if rng.random() < 0.5 else lvl.pw(E(Fr(B)), E(k * m), 1)
        u = lvl.sub(T, p)
        L = lvl.lam_iv(u)
        tr = mpmath.log(m * c.numerator) - mpmath.log(c.denominator) + (m - 1) * k * mpmath.log(B)
        tol = mpmath.mpf(10) ** -250 * abs(tr)
        if not (mpmath.mpf(L.a) - tol <= tr <= mpmath.mpf(L.b) + tol):
            bad += 1; print("ENCLOSURE MISSES", B0, a, k, m, c, s, mpmath.nstr(tr, 30), L)
        elif not lvl.outside(u, 3.14159, 3.1416):
            bad += 1; print("NOT OUTSIDE A WINDOW AT PI", B0, a, k, m, c, s)
        else:
            ok += 1
        # the same value written twice: exact 0
        q2 = lvl.add(E(c), lvl.pw(E(Fr(B0)), E(a * k), 1)) if s > 0 else lvl.sub(lvl.pw(E(Fr(B0)), E(a * k), 1), E(c))
        z = lvl.sub(p, lvl.pw(q2, E(m), 1))
        if z.q == 0:
            zeros += 1
        else:
            bad += 1; print("NOT AN EXACT ZERO", B0, a, k, m, c, s)
    except lvl.Fail as ex:
        bad += 1; print("FAIL", B0, a, k, m, c, s, ex)
print(f"t_lvlbin: {ok} differences decided with the true ln inside, {zeros} exact zeros; PROBLEMS {bad}")
sys.exit(1 if bad else 0)
