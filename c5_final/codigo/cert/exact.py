"""
Exact / rigorous real arithmetic for motor B (C3/C4).

A value is an instance of X with one of two kinds:

  'S'  exact algebraic value: a finite sum  sum_r  c_r * r  where each r is a NORMALISED radical
       r = prod p^(e_p), p prime, e_p rational in (0, 1)  (r = () is the rational part), and c_r a
       nonzero Fraction.  The empty sum is 0.  Normal form: distinct radicals, nonzero coefficients.
       Zero test is exact: by the theorem of Besicovitch-Mordell-Siegel (positive real radicals
       whose pairwise quotients are irrational are linearly independent over Q), a normal-form sum
       with at least one term is nonzero.  Two distinct normalised radicals have irrational quotient
       (unique factorisation), so the theorem applies.
  'N'  numeric value: an expression tree over X values, evaluated with Arb balls at any precision.
       Used for transcendental values (x^y with y irrational) and for operations whose exact
       result is not kept (inverse of a sum, non-integer power of a sum, huge coefficients).

Only nonnegative bases are ever raised to powers (Lemma 2).  Every operation returns None when the
result is undefined (division by an exact zero, 0^y with y <= 0).
"""
from fractions import Fraction as Fr
import sys as _sys
_sys.set_int_max_str_digits(0)        # canonical forms print huge rationals
from functools import lru_cache
import math
from flint import arb, ctx
from sympy import factorint, isprime
import sympy.ntheory.factor_ as _sympy_factor
# sympy 1.14: with limit=, factorint's Fermat step can pass composite factors to its prime-factor cache,
# whose validation then raises ValueError. The cache only speeds up; factor_int checks every factor itself.
_sympy_factor.factor_cache.add = lambda n, factors: None

MAXCOEF_BITS = 60000          # coefficients larger than this are not kept exactly (-> numeric)
MAXTERMS = 24                 # sums with more terms are not kept exactly (-> numeric)


class Unresolved(Exception):
    pass


# ---------------------------------------------------------------- factorisation (cached, checked)
_SMALLP = None


@lru_cache(maxsize=None)
def factor_int(n):
    """prime factorisation of a positive integer as a tuple of (p, k); None if it cannot be obtained
    with certainty.  Every prime factor is either found by trial division (<= 10^6) or is a probable
    prime below 2^64, where the BPSW test used by sympy is proven exact.  The result is checked by
    multiplication."""
    global _SMALLP
    if n == 1:
        return ()
    if n.bit_length() > 400:                      # large: only numbers made of small primes (2^640, 6^200, ...)
        if _SMALLP is None:
            from sympy import primerange
            _SMALLP = list(primerange(2, 10 ** 4))
        m, out = n, []
        for p in _SMALLP:
            if m % p == 0:
                k = 0
                while m % p == 0:
                    m //= p
                    k += 1
                out.append((p, k))
                if m == 1:
                    break
        if m != 1:
            return None
        prod = 1
        for p, k in out:
            prod *= p ** k
        if prod != n:
            raise AssertionError("bad factorisation")
        return tuple(out)
    f = factorint(n, limit=10 ** 6, use_trial=True, use_rho=False, use_pm1=False, use_ecm=False)
    out = []
    prod = 1
    for p, k in f.items():
        if p > 10 ** 6 and (p >= 2 ** 64 or not isprime(p)):
            return None
        out.append((int(p), int(k)))
        prod *= int(p) ** int(k)
    if prod != n:                       # checked by multiplication (trust nothing from the factoriser)
        raise AssertionError("bad factorisation")
    return tuple(sorted(out))


# ---------------------------------------------------------------- radicals
def rad_mul(r1, r2):
    """product of two normalised radicals: (integer carry as Fraction, normalised radical)"""
    d = dict(r1)
    for p, e in r2:
        d[p] = d.get(p, Fr(0)) + e
    carry = 1
    out = []
    for p in sorted(d):
        e = d[p]
        k = e.numerator // e.denominator          # floor
        f = e - k
        if k:
            carry *= p ** k                        # k is 0 or 1 here (both e < 1)
        if f:
            out.append((p, f))
    return Fr(carry), tuple(out)


def rad_inv(r):
    """1/r = coefficient * normalised radical"""
    carry = Fr(1)
    out = []
    for p, e in r:
        carry /= p
        out.append((p, 1 - e))
    return carry, tuple(out)


def bits(q):
    return max(abs(q.numerator).bit_length(), q.denominator.bit_length())


# ---------------------------------------------------------------- the value class
class X:
    __slots__ = ("kind", "terms", "tree", "_key", "_ball")

    def __init__(self, kind, terms=None, tree=None):
        self.kind = kind
        self.terms = terms            # dict radical -> Fraction (kind 'S')
        self.tree = tree              # (op, args...) (kind 'N')
        self._key = None
        self._ball = {}

    # -- constructors
    @staticmethod
    def rat(q):
        q = Fr(q)
        return X("S", {(): q} if q else {})

    @staticmethod
    def from_terms(t):
        t = {r: c for r, c in t.items() if c}
        if len(t) > MAXTERMS or any(bits(c) > MAXCOEF_BITS for c in t.values()):
            return None
        return X("S", t)

    @staticmethod
    def num(op, *args):
        return X("N", tree=(op,) + args)

    # -- predicates
    def is_zero(self):
        return self.kind == "S" and not self.terms

    def is_rational(self):
        return self.kind == "S" and all(r == () for r in self.terms)

    def rational(self):
        return self.terms.get((), Fr(0)) if self.kind == "S" else None

    def is_monomial(self):
        return self.kind == "S" and len(self.terms) == 1

    def key(self):
        """exact values: their normal form; numeric values: the structure of their expression tree
        (equal trees give equal keys; different trees with equal values are not merged here)"""
        if self._key is None:
            if self.kind == "S":
                self._key = ("S", frozenset(self.terms.items()))
            else:
                self._key = ("N", self.tree[0]) + tuple(c.key() for c in self.tree[1:])
        return self._key

    # -- Arb enclosure
    def ball(self, prec):
        b = self._ball.get(prec)
        if b is not None:
            return b
        old = ctx.prec
        ctx.prec = prec
        try:
            if self.kind == "S":
                b = arb(0)
                for r, c in self.terms.items():
                    t = arb(c.numerator) / c.denominator
                    if r:
                        t *= radical_ball(r, prec)
                    b += t
            else:
                b = eval_tree(self.tree, prec)
        finally:
            ctx.prec = old
        self._ball[prec] = b
        return b

    def __repr__(self):
        if self.kind == "S":
            return "S" + repr(sorted(self.terms.items()))
        return "N" + repr(self.tree[0])


_RAD = {}


def radical_ball(r, prec):
    """Arb enclosure of a normalised radical prod p^e (cached; ctx.prec must already be prec)"""
    k = (r, prec)
    b = _RAD.get(k)
    if b is None:
        b = arb(1)
        for p, e in r:
            b *= (arb(p).log() * (arb(e.numerator) / e.denominator)).exp()
        if len(_RAD) > 200000:
            _RAD.clear()
        _RAD[k] = b
    return b


def eval_tree(t, prec):
    op = t[0]
    if op == "add":
        return t[1].ball(prec) + t[2].ball(prec)
    if op == "neg":
        return -t[1].ball(prec)
    if op == "mul":
        return t[1].ball(prec) * t[2].ball(prec)
    if op == "inv":
        b = t[1].ball(prec)
        if b.contains(0):
            raise Unresolved("inverse of a ball containing 0")
        return 1 / b
    if op == "pow":                      # base > 0 (checked when built)
        x = t[1].ball(prec)
        y = t[2].ball(prec)
        if not x > 0:
            raise Unresolved("power base not certified positive")
        return (x.log() * y).exp()
    raise ValueError(op)


ZERO = X.rat(0)
ONE = X.rat(1)
PRECS = (128, 256, 512, 1024, 2048, 4096, 8192)


# ---------------------------------------------------------------- sign / comparisons
def sign(x):
    """exact sign of x (Arb with increasing precision when needed); raises Unresolved"""
    if x.kind == "S":
        if not x.terms:
            return 0
        if len(x.terms) == 1:
            return 1 if next(iter(x.terms.values())) > 0 else -1
        for prec in PRECS:              # nonzero by the radical independence theorem: terminates
            b = x.ball(prec)
            if b > 0:
                return 1
            if b < 0:
                return -1
        raise Unresolved("sign of an algebraic sum not separated")
    for prec in PRECS:
        try:
            b = x.ball(prec)
        except Unresolved:
            continue
        if b > 0:
            return 1
        if b < 0:
            return -1
        if prec == 512:                 # possibly an exact zero: try the canonical form (canon.py)
            import canon
            if canon.is_zero(x):
                return 0
    import canon
    if canon.is_zero(x):
        return 0
    raise Unresolved("sign of a numeric value not decided (possible exact zero)")


# ---------------------------------------------------------------- arithmetic
def neg(x):
    if x.kind == "S":
        return X("S", {r: -c for r, c in x.terms.items()})
    return X.num("neg", x)


def add(x, y):
    if x.kind == "S" and y.kind == "S":
        t = dict(x.terms)
        for r, c in y.terms.items():
            t[r] = t.get(r, Fr(0)) + c
        z = X.from_terms(t)
        if z is not None:
            return z
    if x.is_zero():
        return y
    if y.is_zero():
        return x
    return X.num("add", x, y)


def sub(x, y):
    return add(x, neg(y))


def absval(x):
    sg = sign(x)
    if sg == 0:
        return ZERO                     # an exact zero (possibly found through canon.py)
    return neg(x) if sg < 0 else x


def mul(x, y):
    if x.is_zero() or y.is_zero():
        return ZERO
    if x.kind == "S" and y.kind == "S" and len(x.terms) * len(y.terms) <= MAXTERMS:
        t = {}
        for r1, c1 in x.terms.items():
            for r2, c2 in y.terms.items():
                carry, r = rad_mul(r1, r2)
                t[r] = t.get(r, Fr(0)) + c1 * c2 * carry
        z = X.from_terms(t)
        if z is not None:
            return z
    return X.num("mul", x, y)


def inv(x):
    if x.is_zero():
        return None
    if x.kind == "S" and len(x.terms) == 1:
        (r, c), = x.terms.items()
        carry, ri = rad_inv(r)
        z = X.from_terms({ri: carry / c})
        if z is not None:
            return z
    if sign(x) == 0:
        return None
    return X.num("inv", x)


def div(x, y):
    iy = inv(y)
    return None if iy is None else mul(x, iy)


def rat_pow_monomial(c, r, q):
    """(c * r)^q for c > 0 rational, r a normalised radical, q rational; exact or None"""
    if q.denominator == 1:
        k = q.numerator
        if abs(k) * max(bits(c), 1) > MAXCOEF_BITS:
            return None
        cc = c ** k
        carry = Fr(1)
        rr = ()
        # r^k: multiply exponents, renormalise
        d = {}
        for p, e in r:
            d[p] = e * k
    else:
        fn = factor_int(c.numerator)
        fd = factor_int(c.denominator)
        if fn is None or fd is None:
            return None
        d = {}
        for p, e in r:
            d[p] = d.get(p, Fr(0)) + e
        for p, k in fn:
            d[p] = d.get(p, Fr(0)) + k
        for p, k in fd:
            d[p] = d.get(p, Fr(0)) - k
        d = {p: e * q for p, e in d.items()}
        cc = Fr(1)
    carry = Fr(1)
    out = []
    for p in sorted(d):
        e = d[p]
        k = e.numerator // e.denominator
        f = e - k
        if k:
            if abs(k) * p.bit_length() > MAXCOEF_BITS:
                return None
            carry *= Fr(p) ** k
        if f:
            out.append((p, f))
    return X.from_terms({tuple(out): cc * carry})


def power(x, y):
    """x^y with x >= 0 (Lemma 2); None if undefined"""
    sy = None
    if y.is_zero():
        return None if x.is_zero() else ONE
    if x.is_zero():
        sy = sign(y)
        return ZERO if sy > 0 else None
    if x.kind == "S" and x.terms == {(): Fr(1)}:
        return ONE
    if y.is_rational():
        q = y.rational()
        if x.kind == "S" and len(x.terms) == 1:
            (r, c), = x.terms.items()
            if c < 0:
                raise ValueError("negative base")
            z = rat_pow_monomial(c, r, q)
            if z is not None:
                return z
        elif x.kind == "S" and q.denominator == 1 and 0 < q.numerator <= 6:
            z = x
            for _ in range(q.numerator - 1):
                z = mul(z, x)
            return z
    return X.num("pow", x, y)


# ---------------------------------------------------------------- range
LN_RANGE = None


def ln_range():
    global LN_RANGE
    if LN_RANGE is None:
        old = ctx.prec
        ctx.prec = 256
        LN_RANGE = arb(10).log() * 300
        ctx.prec = old
    return LN_RANGE


def in_range(x):
    """0 or 1e-300 <= |x| <= 1e300; values whose enclosure straddles a bound count as in range.
    Operands of every node are themselves in range (they were checked when produced), so a
    non-finite Arb ball can only come from a power x^y with |y ln x| huge: decided from y*ln x."""
    if x.is_zero():
        return True
    lr = ln_range()
    for prec in (128, 512):
        b = x.ball(prec)
        if not b.is_finite() or b.contains(0):
            # a power x^y (x > 0, so never 0) whose ball is useless: decide from y*ln x
            if x.kind == "N" and x.tree[0] == "pow":
                L = x.tree[2].ball(prec) * x.tree[1].ball(prec).log()
                if L > lr or L < -lr:
                    return False
            continue
        lb = abs(b).log()
        if lb > lr or lb < -lr:
            return False
        return True
    return True


# ---------------------------------------------------------------- adaptive enclosures
TPRECS = (256, 1024, 4096, 16384)


def tight_ball(x, rel_bits=100):
    """an enclosure of x at the lowest precision of TPRECS whose relative radius is <= 2^-rel_bits
    (the best one obtained otherwise).  Deep power towers such as (a^(b^-k))^(c^m) need far more
    working precision than their final accuracy, so a fixed precision can return useless balls."""
    best = None
    for prec in TPRECS:
        try:
            b = x.ball(prec)
        except Unresolved:
            continue
        best = b
        if x.is_zero() or (b.is_finite() and not b.contains(0)
                           and b.rad() <= abs(b.mid()) * (arb(2) ** -rel_bits)):
            return b
    if best is None:
        raise Unresolved("no enclosure at any precision")
    return best


def may_lie(x, lo, hi):
    """False only if x is certainly below lo or certainly above hi (lo, hi: Arb balls or None);
    the precision is raised while the enclosure straddles a bound"""
    for prec in TPRECS:
        try:
            b = x.ball(prec)
        except Unresolved:
            continue
        if (lo is not None and b < lo) or (hi is not None and b > hi):
            return False
        if (lo is None or b > lo) and (hi is None or b < hi):
            return True
    return True
