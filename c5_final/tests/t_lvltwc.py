"""lvl.py, powers of canonical forms (13/11, lvl.pw with base.tw and an exact irrational exponent; lvl.tower_exact).
The 128 undecided values of n = 9, split 10, piece VN.7 were one expression, (|1/5 - (5^-(8^.6))^(4^-.9)|^-7)^(1/1000):
(5^-(8^.6))^(4^-.9) = 5^-(2^(9/5) 2^(-9/5)) = 1/5 exactly, the difference is 0 and 0^-7 has no value.
Random values
    u = (A^(s1 B^x))^(s2 C^y)        C = B^j (j = +-1, +-2), x and y non-integer rationals, B^x and C^y irrational
built with lvl's operations.  Their exact logarithm is ln u = s1 s2 B^(x + j y) ln A.
  - when x + j y makes B^(x + j y) rational (cancelling pairs, y = (k - x)/j): lvl must give an exact value (rational,
    or an exact power of primes) equal to A^(s1 s2 B^(x + j y)), and |u - that value| must be the exact 0 when rational;
  - otherwise (y perturbed): lvl must not give an exact value, and when it decides it, its enclosure of ln u must
    contain the true one (mpmath, 400 digits).
usage: t_lvltwc.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lvl
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 13)
E = lvl.exact
mpmath.mp.dps = 400
AS = [Fr(2), Fr(3), Fr(5), Fr(7), Fr(1, 2), Fr(1, 5), Fr(3, 2), Fr(9, 4), Fr(1, 3), Fr(6), Fr(10)]
BS = [Fr(2), Fr(3), Fr(4), Fr(5), Fr(6), Fr(8), Fr(9), Fr(1, 2), Fr(1, 4), Fr(3, 2), Fr(2, 3), Fr(10)]


def lnq(q):
    return mpmath.log(q.numerator) - mpmath.log(q.denominator)


def vec(q):
    from sympy import factorint
    out = {}
    for p_, a in factorint(q.numerator).items():
        out[int(p_)] = out.get(int(p_), 0) + a
    for p_, a in factorint(q.denominator).items():
        out[int(p_)] = out.get(int(p_), 0) - a
    return out


def rational_power(B, t):
    """B^t is rational"""
    return all((Fr(a) * t).denominator == 1 for a in vec(B).values())


def lvl_ln(v):
    """an enclosure (a, b) of ln v from lvl's result, or None when lvl gives none"""
    if v.q is not None:
        if v.q <= 0:
            return None
        l = lnq(v.q)
        return l, l
    try:
        L = lvl.lam_iv(v)
        return mpmath.mpf(L.a), mpmath.mpf(L.b)
    except lvl.Fail:
        return None


bad = canc = canc_exact = pert = pert_dec = skipped = 0
for trial in range(N):
    A = AS[rng.integers(len(AS))]; B = BS[rng.integers(len(BS))]
    j = int(rng.choice([1, 2, -1, -2])); C = B ** j
    x = Fr(int(rng.integers(1, 60)), int(rng.choice([2, 3, 5, 7, 10]))) * (1 if rng.random() < 0.7 else -1)
    k = int(rng.integers(-2, 4))
    y = Fr(k - x, j)
    cancel = rng.random() < 0.5
    if not cancel:
        y += Fr(int(rng.choice([-1, 1])), int(rng.choice([7, 11, 1000, 10 ** 9])))
    if x.denominator == 1 or y.denominator == 1 or y == 0 or rational_power(B, x) or rational_power(C, y):
        skipped += 1; continue
    t = x + j * y                                        # B^x C^y = B^t
    if abs(float(t) * float(lnq(B))) > 30 or abs(float(x) * float(lnq(B))) > 300:
        skipped += 1; continue
    s1 = 1 if rng.random() < 0.6 else -1; s2 = 1 if rng.random() < 0.6 else -1
    try:
        ex = lvl.pw(E(B), E(abs(x)), 1 if x > 0 else -1)
        b = lvl.pw(E(A), ex, s1)
        ey = lvl.pw(E(C), E(abs(y)), 1 if y > 0 else -1)
        u = lvl.pw(b, ey, s2)
    except lvl.Fail:
        skipped += 1; continue
    L = s1 * s2 * mpmath.exp(t * lnq(B)) * lnq(A)       # the true ln u
    exact_u = u.q is not None or u.pf is not None
    if rational_power(B, t):
        canc += 1
        if not exact_u:
            bad += 1; print("NOT EXACT:", A, B, j, x, y, s1, s2, u); continue
        canc_exact += 1
        Bt = Fr(1)
        for p_, a in vec(B).items():
            Bt *= Fr(p_) ** int(Fr(a) * t)
        # the exact value A^(s1 s2 B^t), as lvl's exact power of A's primes with the rational exponents
        w = {p_: Fr(a) * s1 * s2 * Bt for p_, a in vec(A).items()}
        want = lvl.pf_value(w)
        if want is not None:
            if u.q != want or lvl.sub(u, E(want)).q != 0:
                bad += 1; print("WRONG EXACT VALUE:", A, B, j, x, y, s1, s2, u, want)
        else:
            if u.pf is None or lvl.pf_key(u.pf) != lvl.pf_key(w):
                bad += 1; print("WRONG EXACT VECTOR:", A, B, j, x, y, s1, s2, u.pf, w)
        e_ = lvl_ln(u)
        if e_ is not None and not (e_[0] - mpmath.mpf(10) ** -300 <= L <= e_[1] + mpmath.mpf(10) ** -300):
            bad += 1; print("LN OUTSIDE (exact):", A, B, j, x, y, s1, s2, e_, mpmath.nstr(L, 20))
    else:
        pert += 1
        if u.q is not None or (u.pf is not None and lvl.pf_value(u.pf) is not None):
            bad += 1; print("EXACT RATIONAL FOR AN IRRATIONAL VALUE:", A, B, j, x, y, s1, s2, u); continue
        e_ = lvl_ln(u)
        if e_ is not None:
            pert_dec += 1
            if not (e_[0] <= L <= e_[1]):
                bad += 1; print("LN OUTSIDE:", A, B, j, x, y, s1, s2, e_, mpmath.nstr(L, 20))
# the expression of the lab: (|1/5 - (5^-(8^.6))^(4^-.9)|^-7)^(1/1000) has no value
z = lvl.pw(lvl.pw(E(5), lvl.pw(E(8), E(Fr(3, 5)), 1), -1), lvl.pw(E(4), E(Fr(9, 10)), -1), 1)
d = lvl.sub(E(Fr(1, 5)), z)
lab = "no value (0^-7)"
try:
    lvl.pw(d, E(7), -1); lab = "A VALUE"; bad += 1
except lvl.Invalid:
    pass
print(f"t_lvltwc: {canc} cancelling ({canc_exact} exact), {pert} perturbed ({pert_dec} decided), {skipped} skipped; "
      f"lab expression: {lab}; PROBLEMS {bad}")
sys.exit(1 if bad else 0)
