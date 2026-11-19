"""lvl.py, integer powers of huge exact rationals close to a small one (14/11, lvl.anchor_of_q).
The split 100 of n = 9 had |3^128 - (9 +- 10^-65536)^64|: 9 +- 10^-65536 is an exact rational of 217 000 bits whose
64th power is too large to expand and cannot be factored, so lvl had only a plain enclosure and failed.  Now such a
q is written exactly as M (1 + d) (M of at most 31 bits in the denominator, d = q/M - 1 exactly) and the anchored
rules of pw apply.  Random q = M + s c 10^-K (K large enough that q^k is not expanded), k = 2..80:
  - q^k must contain the true value, and |M^k - q^k| must be decided with the true ln (mpmath, 60 digits on the
    logarithm scale: ln|q^k - M^k| = k ln M + ln|expm1(k log1p(d))|);
  - the result must be outside a window at pi.
usage: t_lvlanq.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lvl
N = int(sys.argv[1]) if len(sys.argv) > 1 else 150
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 100)
E = lvl.exact
mpmath.mp.dps = 60
bad = ok = 0
MS = [Fr(9), Fr(3, 2), Fr(7), Fr(9, 4), Fr(2), Fr(5, 3), Fr(10)]
for trial in range(N):
    M = MS[int(rng.integers(len(MS)))]
    k = int(rng.integers(2, 81))
    K = int(rng.integers(60000, 120000))
    c = Fr(int(rng.integers(1, 99)))
    sg = 1 if rng.random() < 0.5 else -1
    tiny = lvl.pw(E(Fr(1, 10)), E(K), 1)                     # 10^-K exact
    base = lvl.add(E(M), lvl.mul(E(c), tiny)) if sg > 0 else lvl.sub(E(M), lvl.mul(E(c), tiny))
    try:
        x = lvl.pw(base, E(k), 1)
        u = lvl.sub(E(M ** k), x)
    except lvl.Fail as ex:
        bad += 1; print("FAIL", M, k, K, c, sg, ex); continue
    # the truth: d = sg c 10^-K / M, ln|q^k - M^k| = k ln M + ln|expm1(k log1p(d))|
    lnd = mpmath.log(c) - K * mpmath.log(10) - (mpmath.log(M.numerator) - mpmath.log(M.denominator))
    d = sg * mpmath.exp(lnd)
    lnu = k * (mpmath.log(M.numerator) - mpmath.log(M.denominator)) + mpmath.log(abs(mpmath.expm1(k * mpmath.log1p(d))))
    try:
        L = lvl.lam_iv(u)
        if not (mpmath.mpf(L.a) - mpmath.mpf(10) ** -40 * abs(lnu) <= lnu <= mpmath.mpf(L.b) + mpmath.mpf(10) ** -40 * abs(lnu)):
            bad += 1; print("ENCLOSURE MISSES", M, k, K, c, sg, mpmath.nstr(lnu, 25), L)
        else:
            ok += 1
    except lvl.Fail as ex:
        bad += 1; print("NO ENCLOSURE", M, k, K, c, sg, ex)
    if not lvl.outside(u, 3.14159, 3.1416):
        bad += 1; print("NOT DECIDED OUTSIDE A WINDOW AT PI", M, k, K)
# (14/11) q = 2^K + c with K of millions (split 117: |4^(3^14) - (c + 2^(3^14))^2|): q kept as X (1 + d), X = 2^K
# by its vector; |4^K - q^2| = 2 c X + c^2, ln = ln(2c) + K ln 2 + ln(1 + c / (2 X))
ok2 = 0
for trial in range(max(3, N // 25)):
    K = int(rng.integers(3_000_000, 4_500_000))
    c = Fr(int(rng.integers(1, 50)), int(rng.choice([1, 10])))
    try:
        X = lvl.pw(E(Fr(2)), E(K), 1)
        q = lvl.add(E(c), X)
        u = lvl.sub(lvl.pw(E(Fr(4)), E(K), 1), lvl.pw(q, E(2), 1))
        L = lvl.lam_iv(u)
        lnu = mpmath.log(2 * c.numerator) - mpmath.log(c.denominator) + K * mpmath.log(2)
        if not (mpmath.mpf(L.a) - mpmath.mpf(10) ** -30 * lnu <= lnu <= mpmath.mpf(L.b) + mpmath.mpf(10) ** -30 * lnu):
            bad += 1; print("ENCLOSURE MISSES (2^K + c)^2", K, c, mpmath.nstr(lnu, 25), L)
        else:
            ok2 += 1
    except lvl.Fail as ex:
        bad += 1; print("FAIL (2^K + c)^2", K, c, ex)
print(f"t_lvlanq: {ok} powers of huge rationals near small ones, {ok2} squares of 2^K + c; PROBLEMS {bad}")
sys.exit(1 if bad else 0)
