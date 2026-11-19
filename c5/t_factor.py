"""the exact factored-power channel of ivl.fwd (FP, FQ, FE: u = (FP/FQ)^FE): random exact bases q = p/s raised to
exact integers too large to expand, products, quotients and differences of such powers, including equal values
written differently (2^(2k) = 4^k, (8/27)^k = (2/3)^(3k)).  Checked exactly (integer arithmetic on the triples
and Fractions for the bases): the triple equals the value, |x - y| is the exact zero exactly when x = y, x/y is
exactly 1 exactly when x = y, and ln u stays in [A, B] (mpmath).  usage: t_factor.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ivl
mpmath.mp.dps = 50
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 23)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000


def exact_q(p, s):
    g = gcd(p, s); p //= g; s //= g
    L = mpmath.log(mpmath.mpf(p) / s)
    x = ivl.Q7(np.array([np.nextafter(float(L), -np.inf)]), np.array([np.nextafter(float(L), np.inf)]),
               np.array([False]), np.array([p == s]))
    if p == s:
        x.A[:] = 0; x.B[:] = 0
    x.AP = np.array([p], np.int64); x.AQ = np.array([s], np.int64)
    return Fr(p, s), x


def power(qv, q7, k, neg=False):
    kv, k7 = exact_q(k, 1)
    r, ok = ivl.fwd(6 if neg else 5, q7, k7)
    return r


def value_of(r):
    """(base Fraction, exponent) of the triple, or None"""
    if r.FE[0] == 0:
        return None
    return Fr(int(r.FP[0]), int(r.FQ[0])), int(r.FE[0])


def same_value(b1, e1, b2, e2):
    """exact test of b1^e1 == b2^e2 for positive Fractions (by prime exponents of numerator and denominator)"""
    from sympy import factorint
    def pf(b, e):
        d = {}
        for p_, a in factorint(b.numerator).items(): d[p_] = d.get(p_, 0) + a * e
        for p_, a in factorint(b.denominator).items(): d[p_] = d.get(p_, 0) - a * e
        return {p_: a for p_, a in d.items() if a}
    return pf(b1, e1) == pf(b2, e2)


bad = tested = zeros = ones = gaps = 0
bases = [(2, 1), (4, 1), (8, 1), (16, 1), (3, 1), (9, 1), (27, 1), (2, 3), (4, 9), (8, 27), (5, 1), (25, 1), (6, 1), (36, 1),
         (1, 2), (1, 4), (1, 8), (3, 2), (9, 4), (7, 1), (49, 1), (10, 1), (100, 1), (12, 1), (1, 3), (1, 9), (1, 16), (3, 2),
         (27, 8), (1, 10), (1, 100), (4, 1)]
for t in range(N):
    p1, s1 = bases[rng.integers(len(bases))]; p2, s2 = bases[rng.integers(len(bases))]
    k1 = int(rng.integers(30, 10 ** 6))
    # often choose k2 so that the two values are equal when the bases allow it
    q1, x1 = exact_q(p1, s1); q2, x2 = exact_q(p2, s2)
    b1 = ivl.primitive(p1, s1); b2 = ivl.primitive(p2, s2)
    inv = (b1[0], b1[1]) == (b2[1], b2[0]) and b1[0] != b1[1]
    if rng.random() < 0.5 and ((b1[0], b1[1]) == (b2[0], b2[1]) or inv) and (b1[2] * k1) % b2[2] == 0:
        k2 = b1[2] * k1 // b2[2]
    else:
        k2 = int(rng.integers(30, 10 ** 6))
    neg1, neg2 = bool(rng.random() < 0.3), bool(rng.random() < 0.3)
    if (b1[0], b1[1]) == (b2[0], b2[1]) and rng.random() < 0.5:
        neg2 = neg1
    if inv and rng.random() < 0.7:
        neg2 = not neg1
    u1 = power(q1, x1, k1, neg1); u2 = power(q2, x2, k2, neg2)
    for (q, k, neg, u) in ((q1, k1, neg1, u1), (q2, k2, neg2, u2)):
        tested += 1
        L = (-1 if neg else 1) * k * mpmath.log(mpmath.mpf(q.numerator) / q.denominator)
        if not (u.A[0] <= L <= u.B[0]):
            bad += 1; print("BAD enclosure", q, k, neg, u.A[0], u.B[0])
        v = value_of(u)
        if q != 1:
            if v is None:
                if u.AQ[0] == 0:
                    bad += 1; print("MISSING triple", q, k)
            elif not same_value(v[0], v[1], q, (-1 if neg else 1) * k):
                bad += 1; print("BAD triple", q, k, neg, v)
    eq = (q1 != 1 and q2 != 1) and same_value(q1, (-1 if neg1 else 1) * k1, q2, (-1 if neg2 else 1) * k2)
    d, okd = ivl.fwd(1, u1, u2)
    qq, okq = ivl.fwd(3, u1, u2)
    if eq:
        zeros += 1
        if not d.Z[0] or not qq.O[0]:
            bad += 1; print("EQUAL VALUES NOT DETECTED", q1, k1, neg1, q2, k2, neg2, d.Z[0], qq.O[0])
    else:
        if d.Z[0] or qq.O[0]:
            bad += 1; print("FALSE EQUALITY", q1, k1, neg1, q2, k2, neg2)
        elif q1 != 1 and q2 != 1:
            # the enclosure of |x - y| must contain ln|x - y| (mpmath, enough digits to separate the logs)
            mpmath.mp.dps = 80
            l1 = (-1 if neg1 else 1) * k1 * mpmath.log(mpmath.mpf(q1.numerator) / q1.denominator)
            l2 = (-1 if neg2 else 1) * k2 * mpmath.log(mpmath.mpf(q2.numerator) / q2.denominator)
            hi_, lo_ = max(l1, l2), min(l1, l2)
            Ld = hi_ + mpmath.log(-mpmath.expm1(lo_ - hi_))
            mpmath.mp.dps = 50
            if not (d.A[0] <= Ld <= d.B[0]):
                bad += 1; print("BAD difference enclosure", q1, k1, neg1, q2, k2, neg2, d.A[0], d.B[0], mpmath.nstr(Ld, 15))
            gaps += d.A[0] > -np.inf
print(f"tested {tested} powers, {zeros} equal pairs, {gaps} different pairs with a lower bound, bad {bad}")
