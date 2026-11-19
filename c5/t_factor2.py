"""the factored channel with RATIONAL exponents (03/11): u = (FP/FQ)^(FE/FD).  Random exact rational bases raised to
exact rational exponents (often chains: (q^t1)^t2, and equal values written differently such as (2^1000)^(2/5) =
(2^25)^16 or (3^.2)^10000 = 43046721^125), compared with |x - y| and x/y, multiplied with exact rationals of the same
base.  Checked exactly with prime-exponent vectors (Fractions): the triple equals the value; |x - y| is the exact zero
exactly when x = y; x/y is exactly 1 exactly when x = y; the enclosure of |x - y| contains the true value (mpmath,
120 digits) when x != y.  usage: t_factor2.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ivl
from sympy import factorint
mpmath.mp.dps = 120
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 29)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000


def exact_q(p, s):
    g = gcd(p, s); p //= g; s //= g
    L = mpmath.log(mpmath.mpf(p) / s)
    x = ivl.Q7(np.array([np.nextafter(float(L), -np.inf)]), np.array([np.nextafter(float(L), np.inf)]),
               np.array([False]), np.array([p == s]))
    if p == s:
        x.A[:] = 0; x.B[:] = 0
    x.AP = np.array([p], np.int64); x.AQ = np.array([s], np.int64)
    return x


def pf(q, t=Fr(1)):
    d = {}
    for p_, a in factorint(q.numerator).items(): d[p_] = d.get(p_, 0) + a * t
    for p_, a in factorint(q.denominator).items(): d[p_] = d.get(p_, 0) - a * t
    return {p_: a for p_, a in d.items() if a}


def value(r):
    """prime-exponent vector of the triple of r, or None"""
    if r.FE[0] == 0:
        return None
    return pf(Fr(int(r.FP[0]), int(r.FQ[0])), Fr(int(r.FE[0]), int(r.FD[0])))


def lnval(vec):
    return sum(mpmath.mpf(e.numerator) / e.denominator * mpmath.log(p_) for p_, e in vec.items())


bases = [(2, 1), (4, 1), (8, 1), (16, 1), (32, 1), (3, 1), (9, 1), (27, 1), (81, 1), (2, 3), (4, 9), (8, 27), (5, 1), (25, 1),
         (125, 1), (6, 1), (36, 1), (1, 2), (1, 4), (1, 8), (3, 2), (9, 4), (7, 1), (49, 1), (10, 1), (100, 1), (1, 10),
         (33554432, 1), (43046721, 1), (1, 3), (1, 9)]
expos = [(1, 1), (2, 1), (3, 1), (1, 2), (1, 3), (1, 5), (2, 5), (5, 2), (3, 4), (4, 3), (8, 3), (40, 3), (80, 3), (16, 1),
         (125, 1), (1000, 1), (10000, 1), (25, 1), (7, 6), (1, 100)]
bad = tested = zeros = ones = gaps = mixed = 0


def power(q7, t, neg):
    e7 = exact_q(*t)
    r, ok = ivl.fwd(6 if neg else 5, q7, e7)
    return r


for trial in range(N):
    # x = ((b1)^t1)^t2 and y = (b2)^t3, often built to be equal
    b1 = bases[rng.integers(len(bases))]; b2 = bases[rng.integers(len(bases))]
    t1 = expos[rng.integers(len(expos))]; t2 = expos[rng.integers(len(expos))]
    n1, n2 = bool(rng.random() < 0.3), bool(rng.random() < 0.3)
    vx = pf(Fr(*b1), Fr(*t1) * (-1 if n1 else 1)); vx = {p_: e * Fr(*t2) * (-1 if n2 else 1) for p_, e in vx.items()}
    x = power(power(exact_q(*b1), t1, n1), t2, n2)
    # y: the same value from another base when possible (equal pairs), else random
    if rng.random() < 0.5 and vx:
        b2 = bases[rng.integers(len(bases))]
        v2 = pf(Fr(*b2))
        # find t3 = vx / v2 (componentwise, the same ratio for every prime) if it exists
        ratios = {vx[p_] / v2[p_] for p_ in vx if p_ in v2} if set(vx) == set(v2) else set()
        if len(ratios) == 1:
            t3 = ratios.pop()
            n3 = t3 < 0; t3 = abs(t3)
            if t3.numerator < 2 ** 62 and t3.denominator < 2 ** 62:
                y = power(exact_q(*b2), (t3.numerator, t3.denominator), n3)
                vy = {p_: e * t3 * (-1 if n3 else 1) for p_, e in v2.items()}
            else:
                y = None
        else:
            y = None
    else:
        y = None
    if y is None:
        t3 = expos[rng.integers(len(expos))]; n3 = bool(rng.random() < 0.3)
        y = power(exact_q(*b2), t3, n3)
        vy = pf(Fr(*b2), Fr(*t3) * (-1 if n3 else 1))
    for (r, v) in ((x, vx), (y, vy)):
        tested += 1
        L = lnval(v) if v else mpmath.mpf(0)
        if not (r.A[0] <= L <= r.B[0]) and not (r.O[0] and not v):
            bad += 1; print("BAD enclosure", v, r.A[0], r.B[0])
        tv = value(r)
        if tv is not None and tv != v:
            bad += 1; print("BAD triple", v, tv)
        if r.AQ[0] > 0 and r.DS[0] == 0:                 # (R3) an exact rational must be the value itself
            if pf(Fr(int(r.AP[0]), int(r.AQ[0]))) != v:
                bad += 1; print("BAD exact rational", v, r.AP[0], r.AQ[0])
    eq = vx == vy
    d, okd = ivl.fwd(1, x, y)
    qq, okq = ivl.fwd(3, x, y)
    if eq:
        zeros += 1
        if (x.FE[0] != 0 or y.FE[0] != 0) and not (d.Z[0] and qq.O[0]):
            bad += 1; print("EQUAL VALUES NOT DETECTED", b1, t1, t2, b2, t3, d.A[0], d.B[0])
    else:
        if d.Z[0] or qq.O[0]:
            bad += 1; print("FALSE EQUALITY", b1, t1, t2, b2, t3)
        else:
            Lx, Ly = lnval(vx), lnval(vy)
            hi_, lo_ = max(Lx, Ly), min(Lx, Ly)
            Ld = hi_ + mpmath.log(-mpmath.expm1(lo_ - hi_))
            if not (d.A[0] <= Ld <= d.B[0]):
                bad += 1; print("BAD difference enclosure", b1, t1, t2, b2, t3, d.A[0], d.B[0], mpmath.nstr(Ld, 15))
            gaps += d.A[0] > -np.inf
    # product with an exact rational of the same primitive base, compared with the exact result
    if x.FE[0] != 0 and rng.random() < 0.5:
        a, b = int(x.FP[0]), int(x.FQ[0]); k = int(rng.integers(1, 4))
        z = exact_q(a ** k, b ** k)
        pr, okp = ivl.fwd(2, x, z)
        vp = dict(vx)
        for p_, e in pf(Fr(a ** k, b ** k)).items():
            vp[p_] = vp.get(p_, 0) + e
        vp = {p_: e for p_, e in vp.items() if e}
        mixed += 1
        tv = value(pr)
        if not vp:
            if not pr.O[0]:
                bad += 1; print("PRODUCT = 1 NOT DETECTED", vx, a, b, k)
        elif tv is not None and tv != vp:
            bad += 1; print("BAD product triple", vp, tv)
print(f"tested {tested} powers, {zeros} equal pairs, {gaps} different pairs with a lower bound, {mixed} mixed products, bad {bad}")
