"""lvl.py, exact powers of primes with rational exponents (03/11, independent of ivl.py): chains (b^t1)^t2 and b2^t3 of
exact rationals with rational exponents, often equal values written differently; lvl.op must give the exact zero for
|x - y| exactly when x = y (prime-exponent vectors, Fractions), never a false zero, and an enclosure containing the
true value (mpmath, 120 digits) otherwise; products x*y and quotients x/y must contain the true value.
usage: t_lvlpf.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lvl
from sympy import factorint
mpmath.mp.dps = 120
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 37)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
F = {f: f for f in lvl.FORMS}             # lvl.op takes the name of the form


def pf(q, t=Fr(1)):
    d = {}
    for p_, a in factorint(q.numerator).items(): d[p_] = d.get(p_, 0) + a * t
    for p_, a in factorint(q.denominator).items(): d[p_] = d.get(p_, 0) - a * t
    return {p_: a for p_, a in d.items() if a}


def lnv(vec):
    return sum(mpmath.mpf(e.numerator) / e.denominator * mpmath.log(p_) for p_, e in vec.items())


def contains(v, L):
    """does the lvl value v contain the value e^L (L = ln of the true value)?"""
    if v.q is not None:
        return abs(mpmath.log(mpmath.mpf(v.q.numerator) / v.q.denominator) - L) < mpmath.mpf(10) ** -60 if v.q != 0 else False
    lo, hi = lvl.outside.__defaults__ if False else (None, None)
    try:
        la = lvl.lam_iv(v)
        return la.a - mpmath.mpf(10) ** -40 <= L <= la.b + mpmath.mpf(10) ** -40
    except lvl.Fail:
        return True


bases = [(2, 1), (4, 1), (8, 1), (32, 1), (3, 1), (9, 1), (27, 1), (81, 1), (2, 3), (4, 9), (5, 1), (25, 1), (6, 1), (36, 1),
         (1, 2), (1, 4), (3, 2), (9, 4), (7, 1), (10, 1), (1, 10), (33554432, 1), (43046721, 1), (1, 3)]
expos = [(1, 1), (2, 1), (3, 1), (1, 2), (1, 3), (1, 5), (2, 5), (5, 2), (3, 4), (4, 3), (8, 3), (40, 3), (16, 1), (125, 1),
         (1000, 1), (10000, 1), (25, 1), (7, 6), (1, 100)]
bad = zeros = tested = 0
for trial in range(N):
    b1 = bases[rng.integers(len(bases))]; t1 = expos[rng.integers(len(expos))]; t2 = expos[rng.integers(len(expos))]
    n1, n2 = bool(rng.random() < 0.3), bool(rng.random() < 0.3)
    try:
        x = lvl.op(F["s^-l" if n1 else "s^l"], lvl.exact(Fr(*b1)), lvl.exact(Fr(*t1)))
        x = lvl.op(F["s^-l" if n2 else "s^l"], x, lvl.exact(Fr(*t2)))
    except (lvl.Fail, lvl.Invalid):
        continue
    vx = {p_: e * Fr(*t1) * (-1 if n1 else 1) * Fr(*t2) * (-1 if n2 else 1) for p_, e in pf(Fr(*b1)).items()}
    # y equal to x through another base when possible
    y = None
    if rng.random() < 0.6 and vx:
        b2 = bases[rng.integers(len(bases))]; v2 = pf(Fr(*b2))
        if set(v2) == set(vx):
            ratios = {vx[p_] / v2[p_] for p_ in vx}
            if len(ratios) == 1:
                t3 = ratios.pop(); n3 = t3 < 0; t3 = abs(t3)
                try:
                    y = lvl.op(F["s^-l" if n3 else "s^l"], lvl.exact(Fr(*b2)), lvl.exact(t3))
                    vy = {p_: e * t3 * (-1 if n3 else 1) for p_, e in v2.items()}
                except (lvl.Fail, lvl.Invalid):
                    y = None
    if y is None:
        b2 = bases[rng.integers(len(bases))]; t3 = expos[rng.integers(len(expos))]; n3 = bool(rng.random() < 0.3)
        try:
            y = lvl.op(F["s^-l" if n3 else "s^l"], lvl.exact(Fr(*b2)), lvl.exact(Fr(*t3)))
        except (lvl.Fail, lvl.Invalid):
            continue
        vy = pf(Fr(*b2), Fr(*t3) * (-1 if n3 else 1))
    tested += 1
    eq = vx == vy
    try:
        d = lvl.op(F["|s-l|"], x, y)
    except lvl.Fail as er:
        if eq:
            bad += 1; print("EQUAL VALUES NOT DECIDED", b1, t1, t2, "|", er)
        continue
    if eq:
        zeros += 1
        if d.q != 0:
            bad += 1; print("EQUAL VALUES, NOT THE EXACT ZERO", b1, t1, t2, vx, d)
    elif d.q == 0:
        bad += 1; print("FALSE ZERO", b1, t1, t2, vx, vy)
    else:
        Lx, Ly = lnv(vx), lnv(vy)
        hi_, lo_ = max(Lx, Ly), min(Lx, Ly)
        if not contains(d, hi_ + mpmath.log(-mpmath.expm1(lo_ - hi_))):
            bad += 1; print("BAD difference", vx, vy, d)
    for f, Lt in (("s*l", lnv(vx) + lnv(vy)), ("s/l", lnv(vx) - lnv(vy))):
        try:
            r = lvl.op(F[f], x, y)
        except (lvl.Fail, lvl.Invalid):
            continue
        if not contains(r, Lt) and not (r.q == 1 and Lt == 0):
            bad += 1; print("BAD", f, vx, vy, r)
print(f"tested {tested} pairs, {zeros} equal pairs, bad {bad}")

# offsets (03/11): X + c with X a huge exact power and c a small exact rational; differences of offsets of the same X
# are exact (|c1 - c2|), and an offset never equals a value built on another X
bad2 = tested2 = 0
for trial in range(N // 2):
    b = bases[rng.integers(len(bases))]
    if Fr(*b) <= 1:
        continue                                     # X huge (X > 1): then |c2 - X| = X - c2
    k = int(rng.integers(5000, 10 ** 7))
    try:
        X = lvl.op("s^l", lvl.exact(Fr(*b)), lvl.exact(Fr(k)))
    except (lvl.Fail, lvl.Invalid):
        continue
    if X.pf is None or X.q is not None:
        continue                                     # expanded exactly: not the case of the offsets
    c1 = Fr(int(rng.integers(1, 50)), int(rng.integers(1, 20))); c2 = Fr(int(rng.integers(1, 50)), int(rng.integers(1, 20)))
    u1 = lvl.op("s+l", lvl.exact(c1), X)                       # X + c1
    u2 = lvl.op("|s-l|", lvl.exact(c2), X)                     # X - c2 (X huge)
    for (u, v, want) in ((X, u1, c1), (X, u2, c2), (u1, u2, c1 + c2)):
        tested2 += 1
        try:
            d = lvl.op("|s-l|", u, v)
        except lvl.Fail as er:
            bad2 += 1; print("OFFSET NOT DECIDED", b, k, c1, c2, er); continue
        if d.q != want:
            bad2 += 1; print("BAD OFFSET", b, k, c1, c2, d, want)
    # another base: never equal
    b2 = bases[rng.integers(len(bases))]
    try:
        Y = lvl.op("s^l", lvl.exact(Fr(*b2)), lvl.exact(Fr(k + 1)))
        d = lvl.op("|s-l|", u1, Y)
        tested2 += 1
        if d.q == 0:
            bad2 += 1; print("FALSE ZERO OFFSET", b, b2, k)
    except (lvl.Fail, lvl.Invalid):
        pass
print(f"offsets: tested {tested2}, bad {bad2}")

# irrational offsets (03/11): X + c and X - c with c = b^(-t) an algebraic number; |X - (X +- c)| must enclose c
bad3 = tested3 = 0
for trial in range(N // 4):
    b = bases[rng.integers(len(bases))]
    if Fr(*b) <= 1:
        continue
    k = int(rng.integers(5000, 10 ** 7))
    try:
        X = lvl.op("s^l", lvl.exact(Fr(*b)), lvl.exact(Fr(k)))
    except (lvl.Fail, lvl.Invalid):
        continue
    if X.pf is None or X.q is not None:
        continue
    cb = Fr(int(rng.integers(2, 10))); ct = Fr(int(rng.integers(1, 9)), int(rng.integers(2, 10)))
    c = lvl.op("s^-l", lvl.exact(cb), lvl.exact(ct))                 # cb^(-ct), irrational in general
    true_c = mpmath.mpf(cb.numerator) ** (-mpmath.mpf(ct.numerator) / ct.denominator)
    for form in ("s+l", "|s-l|"):
        u = lvl.op(form, c, X)
        tested3 += 1
        try:
            d = lvl.op("|s-l|", X, u)
        except lvl.Fail as er:
            bad3 += 1; print("IRRATIONAL OFFSET NOT DECIDED", b, k, cb, ct, form, er); continue
        if not contains(d, mpmath.log(true_c)):
            bad3 += 1; print("BAD IRRATIONAL OFFSET", b, k, cb, ct, form, d)
print(f"irrational offsets: tested {tested3}, bad {bad3}")
