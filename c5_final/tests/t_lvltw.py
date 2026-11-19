"""lvl.py, the canonical form of exact powers with exact irrational exponents (11/11, lvl.tower; lvl's own code).
The residual hits of n = 9 with a 4-digit table side on both sides are differences of two such powers that are the same
number written twice: 3^(2^(7/.4)) and 9^(8^(6-.5)), 2^(3^8.4) and (1/2)^-((9^7)^.6).  Random pairs
    u1 = A^(B^x)          u2 = (A^(B^t))^((B^i)^((x - t) / i))        (equal: A^(B^t B^(x - t)))
with A, B small rationals, x a non-integer rational, t >= 0 an integer and i >= 1, each built with lvl's operations:
  - equal pairs: |u1 - u2| must be the exact 0;
  - pairs made different (x, t, A or the outer sign changed): the result must never be the exact 0; when lvl decides
    it, it must contain the true value (mpmath: ln|u1 - u2| from the two logarithms).
usage: t_lvltw.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lvl
N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 6)
E = lvl.exact
mpmath.mp.dps = 400
SMALL = [Fr(2), Fr(3), Fr(5), Fr(6), Fr(7), Fr(8), Fr(9), Fr(4), Fr(1, 2), Fr(1, 5), Fr(3, 2), Fr(9, 4), Fr(1, 3)]


def lnq(q):
    return mpmath.log(q.numerator) - mpmath.log(q.denominator)


def build(A, B, x, t, i, neg):
    """(value of lvl, ln of the value) of (A^(B^t))^(+-(B^i)^((x - t)/i)), t = None: A^(+-B^x) directly"""
    if t is None:
        e = lvl.pw(E(B), E(abs(x)), 1 if x > 0 else -1)
        return lvl.pw(E(A), e, -1 if neg else 1), (-1 if neg else 1) * mpmath.exp(x * lnq(B)) * lnq(A)
    base = lvl.pw(E(A), E(B ** t), 1)
    y = (x - t) / i
    e = lvl.pw(E(B ** i), E(abs(y)), 1 if y > 0 else -1)
    return lvl.pw(base, e, -1 if neg else 1), (-1 if neg else 1) * mpmath.exp((x - t) * lnq(B)) * (B ** t) * lnq(A)


bad = eq = eq0 = df = dfdec = 0
for trial in range(N):
    A = SMALL[rng.integers(len(SMALL))]; B = SMALL[rng.integers(9)]          # B > 1
    x = Fr(int(rng.integers(5, 400)), int(rng.choice([2, 3, 4, 5, 10]))) * (1 if rng.random() < 0.85 else -1)
    if x.denominator == 1:
        x += Fr(1, 2)
    t = int(rng.integers(0, 3)); i = int(rng.integers(1, 4)); neg = rng.random() < 0.3
    if B ** t * abs(lnq(A)) > 10 ** 6 or abs(float(x) * float(lnq(B))) > 600:
        continue
    same = rng.random() < 0.5
    try:
        u1, L1 = build(A, B, x, None, 1, neg)
        if same:
            u2, L2 = build(A, B, x, t, i, neg)
        else:
            k = int(rng.integers(0, 4))
            x2 = x + Fr(int(rng.choice([-1, 1])), int(rng.choice([1, 7, 1000, 10 ** 9]))) if k == 0 else x
            t2 = t + 1 if k == 1 else t
            A2 = SMALL[(SMALL.index(A) + 1) % len(SMALL)] if k == 2 else A
            u2, L2 = build(A2, B, x2, t2, i, (not neg) if k == 3 else neg)
            if k == 1:
                L2 = L2                       # (A^(B^(t+1)))^(B^(x-t-1)) is again A^(B^x): an equal pair after all
                same = True
        z = lvl.sub(u1, u2)
    except lvl.Fail:
        if same:
            bad += 1; print("EQUAL pair not recognised:", A, B, x, t, i, neg)
        else:
            df += 1
        continue
    except lvl.Invalid:
        continue
    zero = z.q is not None and z.q == 0
    if same:
        eq += 1; eq0 += zero
        if not zero:
            bad += 1; print("EQUAL pair, not the exact 0:", A, B, x, t, i, neg, z)
        continue
    df += 1
    if zero:
        bad += 1; print("FALSE ZERO:", A, B, x, t, i, neg); continue
    if z.far is not None:
        dfdec += 1; continue
    # the true ln|u1 - u2| from L1, L2 (the larger one dominates unless they are close)
    hi, lo = max(L1, L2), min(L1, L2)
    if hi - lo >= 10 ** 4:
        Lt = hi                                          # (the relative correction is below e^-10000)
    else:
        Lt = hi + mpmath.log(-mpmath.expm1(lo - hi))     # ln(e^hi - e^lo), also for values within 1e-150 of 1
    try:
        la = lvl.lam_iv(z)
    except lvl.Fail:
        dfdec += 1; continue
    tol = mpmath.mpf(10) ** -25 * (1 + abs(Lt))
    if not (la.a - tol <= Lt <= la.b + tol):
        bad += 1; print("MISSES the true value:", A, B, x, t, i, neg, mpmath.nstr(Lt, 25))
    else:
        dfdec += 1
print(f"t_lvltw: {eq} equal pairs ({eq0} exact zeros), {df} different pairs ({dfdec} decided, none a false zero); problems {bad}")
sys.exit(1 if bad or eq != eq0 or not eq else 0)
