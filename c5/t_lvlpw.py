"""lvl.pw with an exponent (or base) whose exact rational value is known only through its prime-exponent vector (03/11):
e.g. 5^((10^(1/2))^6) = 5^1000, 8^((3 (9/10)^(-3/5))^5) = 8^(1000/3) = 2^1000 (the undecided values of split 45 of n = 9).
Random cases B^Y with Y = (m^(1/r))^(r j) * c^k ... built in lvl from exact rationals and rational powers (so Y carries a
vector, not a rational), against B2^K2 with B2^K2 = B^K exactly (equal pairs) or perturbed (different pairs):
  equal pairs: |x - y| must be the exact zero;  different pairs: never a false zero, and |x - y|, x, y must contain
  the true value (mpmath, 200 digits, ln domain).  Also the old behaviour where no exact value is known (Y irrational).
usage: t_lvlpw.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lvl
mpmath.mp.dps = 200
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 11)
TOL = mpmath.mpf(10) ** -40


def lnq(q):
    return mpmath.log(mpmath.mpf(q.numerator)) - mpmath.log(mpmath.mpf(q.denominator))


def contains(v, L):
    """the lvl value v contains e^L (L = -inf: the value 0)"""
    if v.q is not None:
        if v.q == 0:
            return L == -mpmath.inf
        return L != -mpmath.inf and abs(lnq(v.q) - L) < TOL * (1 + abs(L))
    if L == -mpmath.inf:
        return True                    # an enclosure may contain 0 (checked: it is not a false exact zero)
    try:
        la = lvl.lam_iv(v)
    except lvl.Fail:
        return True
    return la.a - TOL * (1 + abs(L)) <= L <= la.b + TOL * (1 + abs(L))


def exponent(K):
    """an lvl value equal to the rational K > 0, built so that it carries only a prime-exponent vector (no .q):
    ((m^(1/r))^(r) ...) written as a product of rational powers.  Returns (value, ln K)"""
    kind = rng.integers(0, 3)
    r = int(rng.integers(2, 7))
    if kind == 0:                       # (K^(1/r))^r
        y = lvl.pw(lvl.exact(K), lvl.exact(Fr(1, r)), 1)
        y = lvl.pw(y, lvl.exact(r), 1)
    elif kind == 1:                     # (c * (c^r / K)^(-1/r))^r  = c^r * K / c^r = K
        c = Fr(int(rng.integers(2, 12)), int(rng.integers(1, 5)))
        inner = lvl.pw(lvl.exact(c ** r / K), lvl.exact(Fr(1, r)), -1)
        y = lvl.pw(lvl.mul(lvl.exact(c), inner), lvl.exact(r), 1)
    else:                               # (K^(a/r))^(r/a) through two steps: (K^(1/r))^a, then ^(r/a)... kept integer: a | r
        a = int(rng.integers(1, 4))
        y = lvl.pw(lvl.exact(K), lvl.exact(Fr(a, r)), 1)
        y = lvl.pw(y, lvl.exact(r), 1)          # = K^a
        y = lvl.pw(y, lvl.exact(Fr(1, a)), 1) if a > 1 else y   # (rational exponent 1/a of a vector): K again
    return y


bad = eq_n = eqz = diff_n = old_n = 0
for t in range(N):
    p = int(rng.choice([2, 3, 5, 6, 7, 10]))
    d = int(rng.choice([1, 2, 3]))             # B = p^d^? : the two writings B^K = B2^K2 with B2 = p^d, K2 = K / d
    K = Fr(int(rng.integers(1, 1200)), int(rng.choice([1, 1, 2, 3, 5])))
    sg = 1 if rng.random() < 0.7 else -1        # s^l or s^-l
    try:
        y = exponent(K)
        x1 = lvl.pw(lvl.exact(Fr(p)), y, sg)                     # p^(sg K)
        equal = rng.random() < 0.6
        K2 = K / d if equal else K / d + Fr(int(rng.integers(1, 9)), int(rng.integers(2, 50)))
        x2 = lvl.pw(lvl.exact(Fr(p ** d)), lvl.exact(K2), sg)    # (p^d)^(sg K2)
        L1 = sg * K * mpmath.log(p); L2 = sg * K2 * d * mpmath.log(p)
        z = lvl.op("|s-l|", x1, x2)
    except lvl.Fail as f_:
        if equal:
            bad += 1; print("FAIL on an equal pair:", p, d, K, sg, f_)
        continue
    if not contains(x1, L1) or not contains(x2, L2):
        bad += 1; print("enclosure misses an operand", p, d, K, K2, sg)
    if equal:
        eq_n += 1
        if z.q == 0:
            eqz += 1
        else:
            bad += 1; print("NOT an exact zero:", p, d, K, sg, z)
    else:
        diff_n += 1
        if z.q is not None and z.q == 0:
            bad += 1; print("FALSE ZERO:", p, d, K, K2, sg)
        a, b = (L1, L2) if L1 >= L2 else (L2, L1)
        Lz = a + mpmath.log(-mpmath.expm1(b - a))
        if not contains(z, Lz):
            bad += 1; print("difference misses the true value", p, d, K, K2, sg, z)
# the old path: an irrational exponent (no exact value) still gives an enclosure
for t in range(200):
    p = int(rng.choice([2, 3, 5])); r = int(rng.integers(2, 7))
    m = int(rng.integers(2, 40))
    y = lvl.pw(lvl.exact(Fr(m)), lvl.exact(Fr(1, r)), 1)     # irrational unless m is a perfect r-th power
    x = lvl.pw(lvl.exact(Fr(p)), y, 1)
    if not contains(x, mpmath.mpf(m) ** (mpmath.mpf(1) / r) * mpmath.log(p)):
        bad += 1; print("irrational exponent: enclosure misses", p, m, r, x)
    old_n += 1
print(f"t_lvlpw: {eq_n} equal pairs ({eqz} exact zeros), {diff_n} different pairs, {old_n} irrational exponents; problems {bad}")
