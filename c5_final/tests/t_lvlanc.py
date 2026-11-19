"""lvl.py, first-order anchors at exact powers (03/11, mechanism B; lvl's own code, independent of ivl.py).  Random
values whose difference with an exact value is tiny but not zero, the residual hits of n = 9:
  F1  q^(a/r) with q = M^r + c (exact rational, many bits):     |M^a - q^(a/r)|          (near_root)
  F2  (B0 + s eps)^(a/r), then ^k (B0^(a k / r) rational):         |B0^(a k/r) - ...|      (anchored base, apf)
  F3  B^(c + s eps), then ^k (B^(c k) rational):                    |B^(c k) - ...|         (anchored exponent)
  F4  Q (m^(1/r) + s eps), then ^r:                                  |Q^r m - ...|           (exact power + tiny, mul, apf)
with eps = (p1/p2)^K tiny (ln eps down to -3e6).  lvl must decide every case at its base precision, the result must
contain the true value (mpmath; ln domain, F1 at the precision of q) and never be a false exact zero.
usage: t_lvlanc.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lvl
N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
E = lvl.exact


def contains_ln(v, L):
    """the lvl value v (a magnitude) contains e^L"""
    if v.q is not None:
        return v.q != 0 and abs(mpmath.log(mpmath.mpf(v.q.numerator)) - mpmath.log(mpmath.mpf(v.q.denominator)) - L) < 1e-25 * (1 + abs(L))
    la = lvl.lam_iv(v)
    tol = mpmath.mpf(10) ** -25 * (1 + abs(L))
    return la.a - tol <= L <= la.b + tol


def tiny():
    p1, p2 = int(rng.integers(2, 10)), int(rng.integers(2, 10))
    if p1 >= p2:
        p1, p2 = p2 - 1 if p2 > 2 else 1, p1 + 1
    lr = np.log2(p2)
    if rng.random() < 0.5:     # an exact rational of up to ~8000 bits (|ln eps| >= 800)
        K = int(rng.uniform(800 / np.log(p2 / p1), max(800 / np.log(p2 / p1) + 1, 8000 / lr)))
    else:                      # beyond lvl's exact size: kept as a vector (ln eps down to ~ -3e6)
        K = int(rng.uniform(8.1e6 / lr, 1.3e7 / lr))
    global EPS
    EPS = (p1, p2, K)
    return lvl.pw(E(Fr(p1, p2)), E(K), 1), K * (mpmath.log(p1) - mpmath.log(p2))     # value, ln eps (exact formula)


bad = dec = tot = 0
EPS = None
fam = {1: [0, 0], 2: [0, 0], 3: [0, 0], 4: [0, 0]}
for t in range(N):
    F = int(rng.integers(1, 5))
    DETAIL = None
    s = 1 if rng.random() < 0.5 else -1
    try:
        mpmath.mp.dps = 60
        if F == 1:
            M = Fr(int(rng.integers(2, 10)), int(rng.choice([1, 1, 1, 2, 3, 7])))
            if M == 1:
                M = Fr(3, 2)
            lb = max(1.0, float(abs(mpmath.log(M.numerator) - mpmath.log(M.denominator)) / mpmath.log(2)))
            r = max(2, int(rng.uniform(450, 6000) / lb))
            a = int(rng.choice([-3, -2, -1, 1, 2, 3]))
            c = int(rng.integers(1, 50)) * s
            q = M ** r + c
            if q <= 0:
                continue                      # (lvl takes magnitudes: only positive q)
            v = lvl.pw(E(q), E(Fr(abs(a), r)), 1 if a > 0 else -1)
            ref = M ** a
            z = lvl.sub(E(ref), v)
            mpmath.mp.dps = int(r * lb * 0.31) + 80
            tv = mpmath.mpf(q.numerator) / q.denominator
            L = mpmath.log(abs(mpmath.mpf(ref.numerator) / ref.denominator - tv ** (mpmath.mpf(a) / r)))
        elif F == 2:
            m = int(rng.integers(2, 8)); r = int(rng.integers(2, 6))
            B0 = Fr(m) ** r if rng.random() < 0.5 else Fr(m)          # B0^(1/r) rational, or irrational (apf path)
            a = int(rng.choice([1, 2, 3])); k = r                      # exponent a/r, then ^r: B0^a, rational
            ev, le = tiny()
            DETAIL = (int(B0), r, a, s, float(le))
            b = lvl.add(E(B0), ev) if s > 0 else lvl.sub(E(B0), ev)
            v = lvl.pw(lvl.pw(b, E(Fr(a, r)), 1), E(k), 1)
            ref = B0 ** a
            z = lvl.sub(E(ref), v)
            x = s * mpmath.exp(le - mpmath.log(int(B0)))      # eps / B0, sign s
            L = mpmath.log(int(ref)) + mpmath.log(abs(mpmath.expm1(a * mpmath.log1p(x))))
        elif F == 3:
            B = int(rng.integers(2, 10)); r = int(rng.integers(2, 6)); c = Fr(int(rng.integers(1, 4)), r)
            ev, le = tiny()
            DETAIL = (B, r, str(c), s, float(le))
            ex_ = lvl.add(E(c), ev) if s > 0 else lvl.sub(E(c), ev)
            v = lvl.pw(lvl.pw(E(B), ex_, 1), E(r), 1)                 # B^(r c + r s eps)
            ref = Fr(B) ** int(c * r)
            z = lvl.sub(E(ref), v)
            w = r * s * mpmath.exp(le) * mpmath.log(B)
            L = mpmath.log(int(ref)) + mpmath.log(abs(mpmath.expm1(w)))
        else:
            m = int(rng.choice([2, 3, 5, 6, 7, 10])); r = int(rng.integers(2, 6)); Q = Fr(int(rng.integers(1, 6)), int(rng.integers(1, 4)))
            P = lvl.pw(E(m), E(Fr(1, r)), 1)                          # m^(1/r), irrational: a vector
            ev, le = tiny()
            b = lvl.add(P, ev) if s > 0 else lvl.sub(P, ev)
            if Q != 1:
                b = lvl.mul(E(Q), b)
            v = lvl.pw(b, E(r), 1)
            ref = Q ** r * m
            z = lvl.sub(E(ref), v)
            x = s * mpmath.exp(le - mpmath.log(m) / r)
            L = mpmath.log(mpmath.mpf(ref.numerator) / ref.denominator) + mpmath.log(abs(mpmath.expm1(r * mpmath.log1p(x))))
    except lvl.Fail as f_:
        tot += 1; fam[F][0] += 1
        print("undecided:", F, f_, DETAIL, EPS)
        continue
    tot += 1; fam[F][0] += 1
    if z.q is not None and z.q == 0:
        bad += 1; print("FALSE ZERO", F); continue
    if not contains_ln(z, L):
        bad += 1; print("MISSES the true value", F, z, mpmath.nstr(L, 30)); continue
    dec += 1; fam[F][1] += 1
# (10/11) F1 with large numbers (near_root up to NEAR_BITS = 4e6 bits): q = M^r + c with r up to 2^21, as in the hits
# |5 - (5^262144 + 7)^(1/262144)| of n = 9.  The reference is the closed form: with x = c / M^r,
# |M^a - q^(a/r)| = M^a |(1 + x)^(a/r) - 1| = M^a |a/r| |x| (1 + O(x)), and |x| < 2^-500000
big = [0, 0]
mpmath.mp.dps = 60
for (M, r, c, a) in ((Fr(5), 262144, 7, 1), (Fr(5), 262144, -3, 1), (Fr(3), 1048576, 6, 1), (Fr(2), 2097152, -1, 2),
                     (Fr(7, 2), 300000, 5, -1), (Fr(9), 524288, Fr(7, 10), 3), (Fr(6), 78125 * 8, -Fr(1, 3), 1),
                     # a large M (the hits |8^15625 - (4^78125 + 1)^(3/10)| of n = 8: M = 2^15625, r = 10, a = 3)
                     (Fr(2) ** 15625, 10, 1, 3), (Fr(2) ** 15625, 10, -Fr(1, 10), 3), (Fr(3 ** 2000, 2 ** 1000), 7, -5, 2),
                     (Fr(5) ** 30000, 3, 2, -1), (Fr(7 ** 500, 10 ** 400), 12, Fr(3, 7), 5)):
    q = M ** r + c
    tot += 1; big[0] += 1
    try:
        v = lvl.pw(E(q), E(Fr(abs(a), r)), 1 if a > 0 else -1)
        z = lvl.sub(E(M ** a), v)
    except lvl.Fail as f_:
        print("undecided (large F1): M of", M.numerator.bit_length(), "bits, r, c, a =", r, c, a, f_)
        continue
    lM = mpmath.log(M.numerator) - mpmath.log(M.denominator)
    L = a * lM + mpmath.log(mpmath.mpf(abs(a)) / r) + mpmath.log(mpmath.mpf(abs(c).numerator) / abs(c).denominator) - r * lM
    if z.q is not None and z.q == 0:
        bad += 1; print("FALSE ZERO (large F1)", M.numerator.bit_length(), r, c, a); continue
    if not contains_ln(z, L):
        bad += 1; print("MISSES the true value (large F1)", M.numerator.bit_length(), r, c, a, mpmath.nstr(L, 30)); continue
    dec += 1; big[1] += 1
print("families (cases, decided):", fam, "large F1:", big)
print(f"t_lvlanc: {tot} cases, decided {dec}, problems {bad}")
sys.exit(1 if bad or dec != tot else 0)
