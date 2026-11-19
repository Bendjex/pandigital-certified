"""soundness of the extension of rule C.4.6(8) to exact RATIONAL exponents (03/11) and of the exact 1/1 of the flag O:
random bases q = P/Q that are perfect powers (P^p and Q^p exact r-th powers) or not, exact or anchored q(1 + d) with
|d| from 1e-9 down to 1e-600 (both signs), exponents t = p/r, forms s^l, s^-l, l^s, l^-s; mpmath at 1400 digits.
Checked: when the rule applies, the anchor rational is q^(+-t) exactly, sign(delta') and ln|delta'| in [DM, DN]; ln u in
[A, B]; and a result with the flag O carries the exact rational 1/1.  usage: t_ancrat.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ivl
mpmath.mp.dps = 1400


def enc(v):
    a = float(v)
    return np.nextafter(a, -np.inf), np.nextafter(a, np.inf)


def q7_exact(p, r):
    L = mpmath.log(mpmath.mpf(p) / r)
    lo, hi = enc(L)
    x = ivl.Q7(np.array([lo]), np.array([hi]), np.array([False]), np.array([p == r]))
    if p == r:
        x.A[:] = 0; x.B[:] = 0
    x.AP = np.array([p], np.int64); x.AQ = np.array([r], np.int64)
    return mpmath.mpf(p) / r, x


def q7_anchored(p, r, d):
    v = mpmath.mpf(p) / r * (1 + d)
    L = mpmath.log(v)
    lo, hi = enc(L)
    x = ivl.Q7(np.array([lo]), np.array([hi]), np.array([False]), np.array([False]))
    x.AP = np.array([p], np.int64); x.AQ = np.array([r], np.int64)
    x.DS = np.array([1 if d > 0 else -1], np.int8)
    a, b = enc(mpmath.log(abs(d))); x.DM = np.array([a]); x.DN = np.array([b])
    if p == r:                                           # anchored at 1: also the second level, as in practice
        m = mpmath.log(abs(L)); a, b = enc(m)
        x.S = np.array([1 if L > 0 else -1], np.int8); x.M = np.array([a]); x.N = np.array([b])
    return v, x


rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 71)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000

roots = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12]
bad = tested = applied = exact_res = ones = 0
for trial in range(N):
    r_ = int(rng.choice([2, 3, 4, 5, 6]))
    a, b = int(rng.choice(roots)), int(rng.choice(roots))
    if rng.random() < 0.8:
        P, Q = a ** r_, b ** r_                    # a perfect r-th power: q^(1/r) rational
    else:
        P, Q = int(rng.integers(2, 60)), int(rng.integers(1, 30))
    g = gcd(P, Q); P //= g; Q //= g
    p_ = int(rng.integers(1, 6)); g = gcd(p_, r_)
    if g > 1:
        p_ //= g; rr = r_ // g
    else:
        rr = r_
    t = Fr(p_, rr)
    if t.denominator == 1:
        continue                                   # the integer case is rule 8 itself (t_ancpow.py)
    if rng.random() < 0.3:
        bv, base = q7_exact(P, Q); dd = 0
    else:
        dd = mpmath.mpf(int(rng.choice([-1, 1]))) * mpmath.mpf(10) ** float(-rng.uniform(9, 600))
        bv, base = q7_anchored(P, Q, dd)
    if bv == 1:
        continue
    _, e7 = q7_exact(t.numerator, t.denominator)
    f = int(rng.choice([5, 6, 7, 8]))
    x, y = (base, e7) if f in (5, 6) else (e7, base)
    sg = -1 if f in (6, 8) else 1
    r, ok = ivl.fwd(f, x, y)
    if not ok[0]:
        continue
    tested += 1
    u = bv ** (sg * mpmath.mpf(t.numerator) / t.denominator)
    L = mpmath.log(u)
    if not (r.A[0] <= L <= r.B[0]):
        bad += 1; print("BAD enclosure", P, Q, t, f, mpmath.nstr(dd, 5), r.A[0], r.B[0]); continue
    if r.O[0] and not (r.AP[0] == 1 and r.AQ[0] == 1 and r.DS[0] == 0):
        bad += 1; print("O WITHOUT 1/1", P, Q, t, f)
    if r.AQ[0] > 0:
        qa = mpmath.mpf(int(r.AP[0])) / int(r.AQ[0])
        if r.DS[0] == 0:
            exact_res += 1
            if abs(u - qa) > abs(u) * mpmath.mpf(10) ** -1300:
                bad += 1; print("BAD exact", P, Q, t, f, r.AP[0], r.AQ[0], mpmath.nstr(u, 20))
        else:
            applied += 1
            dl = u / qa - 1
            if dl == 0 or (dl > 0) != (r.DS[0] > 0):
                bad += 1; print("BAD sign", P, Q, t, f, mpmath.nstr(dd, 5), r.DS[0]); continue
            ld = mpmath.log(abs(dl))
            if not (r.DM[0] <= ld <= r.DN[0]):
                bad += 1; print("BAD delta", P, Q, t, f, mpmath.nstr(dd, 5), mpmath.nstr(ld, 15), r.DM[0], r.DN[0])
    # the flag O: 1^t and q^0-like results
    o, ok1 = ivl.fwd(5, q7_exact(1, 1)[1], e7)
    if o.O[0]:
        ones += 1
        if not (o.AP[0] == 1 and o.AQ[0] == 1 and o.DS[0] == 0):
            bad += 1; print("O WITHOUT 1/1 for 1^t", t)
print(f"tested {tested}, anchored results {applied}, exact results {exact_res}, exact ones {ones}, bad {bad}")

# (03/11) rule 7 with an anchored rational exponent t(1 + d) and an exact base q with q^t rational (R4)
bad4 = tested4 = applied4 = 0
for trial in range(N):
    r_ = int(rng.choice([2, 3, 4, 5]))
    a, b = int(rng.choice(roots[1:])), int(rng.choice(roots))
    P, Q = a ** r_, b ** r_
    g = gcd(P, Q); P //= g; Q //= g
    if P == Q:
        continue
    p_ = int(rng.integers(1, 7)); g = gcd(p_, r_); p_ //= g; rr = r_ // g
    if rr == 1:
        continue
    t = Fr(p_, rr)
    dd = mpmath.mpf(int(rng.choice([-1, 1]))) * mpmath.mpf(10) ** float(-rng.uniform(9, 600))
    bv, base = q7_exact(P, Q)
    tv, e7 = q7_anchored(t.numerator, t.denominator, dd)           # the exponent t (1 + d), anchored at t
    f = int(rng.choice([5, 6, 7, 8]))
    x, y = (base, e7) if f in (5, 6) else (e7, base)
    sg = -1 if f in (6, 8) else 1
    r, ok = ivl.fwd(f, x, y)
    if not ok[0]:
        continue
    tested4 += 1
    u = bv ** (sg * tv)
    if not (r.A[0] <= mpmath.log(u) <= r.B[0]):
        bad4 += 1; print("R4 BAD enclosure", P, Q, t, f); continue
    if r.AQ[0] > 0 and r.DS[0] != 0:
        applied4 += 1
        qa = mpmath.mpf(int(r.AP[0])) / int(r.AQ[0])
        dl = u / qa - 1
        if dl == 0 or (dl > 0) != (r.DS[0] > 0):
            bad4 += 1; print("R4 BAD sign", P, Q, t, f, mpmath.nstr(dd, 5)); continue
        ld = mpmath.log(abs(dl))
        if not (r.DM[0] <= ld <= r.DN[0]):
            bad4 += 1; print("R4 BAD delta", P, Q, t, f, mpmath.nstr(dd, 5), mpmath.nstr(ld, 15), r.DM[0], r.DN[0])
print(f"R4: tested {tested4}, anchored results {applied4}, bad {bad4}")
