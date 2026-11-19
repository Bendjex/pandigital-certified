"""
Canonical forms for motor B values (exact equality of differently written values).

Every value X (kinds 'S' and 'N' of exact.py) is mapped to a LINEAR FORM: a finite sum of
coefficient * monomial, with nonzero Fraction coefficients and distinct monomials; a MONOMIAL is a
finite product of powers  base^exponent  with distinct bases and nonzero exponents, where

    base      ('p', prime)          a prime number
              ('r', q)              a positive rational that could not be factored (huge)
              ('s', linear form)    a positive sum with >= 2 terms, normalised to leading coefficient +-1
    exponent  a linear form (any real value)

Only identities that hold for positive bases and real exponents are used:
    (a b)^Y = a^Y b^Y,   (a^E)^Y = a^(E Y),   a^E a^F = a^(E+F),   distributivity,
and the integer part of the constant term of a prime's exponent is moved into the coefficient
(p^(k + E) = p^k p^E).  Hence two values with EQUAL canonical forms are EQUAL real numbers (the forms
are themselves expressions denoting the values).  Different forms do not imply different values.
"""
from fractions import Fraction as Fr
import exact as ex

MAXK = 4096            # integer parts of exponents moved to the coefficient only below this size
MAXEXPAND = 4          # integer powers of sums expanded up to this exponent


class Unsupported(Exception):
    pass


# ------------------------------------------------------------------ linear forms
# a linear form is a tuple of (monomial, Fraction) sorted by monomial; a monomial is a tuple of
# (base, exponent) sorted by base; this makes them hashable and comparable

def lin(d):
    return tuple(sorted(((m, c) for m, c in d.items() if c != 0), key=lambda t: repr(t[0])))


def lin_const(q):
    return lin({(): Fr(q)})


ZERO = ()
ONE = lin_const(1)


def lin_add(a, b):
    d = dict(a)
    for m, c in b:
        d[m] = d.get(m, Fr(0)) + c
    return lin(d)


def lin_scale(a, q):
    return lin({m: c * q for m, c in a})


def const_of(a):
    """the value if the form is a rational constant, else None"""
    if a == ZERO:
        return Fr(0)
    if len(a) == 1 and a[0][0] == ():
        return a[0][1]
    return None


def mono_mul(m1, m2):
    """product of two monomials: (coefficient, monomial)"""
    d = dict(m1)
    for b, e in m2:
        d[b] = lin_add(d[b], e) if b in d else e
    return normalise_mono(d)


def normalise_mono(d):
    coef = Fr(1)
    out = {}
    for b, e in d.items():
        if e == ZERO:
            continue
        if b[0] == 'p':
            c0 = const_of_part(e)
            k = c0.numerator // c0.denominator
            if k and abs(k) <= MAXK:
                coef *= Fr(b[1]) ** k
                e = lin_add(e, lin_const(-k))
                if e == ZERO:
                    continue
        out[b] = e
    return coef, tuple(sorted(out.items(), key=lambda t: repr(t[0])))


def const_of_part(e):
    for m, c in e:
        if m == ():
            return c
    return Fr(0)


def lin_mul(a, b):
    d = {}
    for m1, c1 in a:
        for m2, c2 in b:
            k, m = mono_mul(m1, m2)
            d[m] = d.get(m, Fr(0)) + c1 * c2 * k
    return lin(d)


# ------------------------------------------------------------------ powers
def rat_pow_form(q, Y):
    """q^Y for a positive rational q and an exponent form Y"""
    if q == 1:
        return ONE
    fn = ex.factor_int(q.numerator) if q.numerator > 1 else ()
    fd = ex.factor_int(q.denominator) if q.denominator > 1 else ()
    if fn is None or fd is None:
        k, m = normalise_mono({('r', q): Y})
        return lin({m: k})
    d = {}
    for p, k in fn:
        d[('p', p)] = lin_scale(Y, k)
    for p, k in fd:
        d[('p', p)] = lin_add(d.get(('p', p), ZERO), lin_scale(Y, -k))
    k, m = normalise_mono(d)
    return lin({m: k})


def mono_pow(m, Y):
    """m^Y for a monomial m (positive): exponents multiplied by Y"""
    d = {b: lin_mul(e, Y) for b, e in m}
    return normalise_mono(d)


def sum_base(a):
    """a positive sum with >= 2 terms as (coefficient > 0, normalised base)"""
    lead = a[0][1]
    c = abs(lead)
    return c, ('s', lin_scale(a, 1 / c))


def lin_pow(a, Y):
    """a^Y for a positive value with form a and an exponent form Y"""
    if a == ZERO:
        raise Unsupported("0^Y")
    yc = const_of(Y)
    if yc == 0:
        return ONE
    if len(a) == 1:
        m, c = a[0]
        if c < 0:
            raise Unsupported("negative single term as a base")
        out = rat_pow_form(c, Y)
        k, mm = mono_pow(m, Y)
        return lin_mul(out, lin({mm: k}))
    if yc is not None and yc.denominator == 1 and 1 <= yc <= MAXEXPAND:
        out = ONE
        for _ in range(int(yc)):
            out = lin_mul(out, a)
        return out
    c, B = sum_base(a)
    out = rat_pow_form(c, Y)
    k, m = normalise_mono({B: Y})
    return lin_mul(out, lin({m: k}))


def expand(a):
    """distribute sum-bases raised to small positive integer constants: (..) * B^k -> (..) * B * ... * B"""
    for i, (m, c) in enumerate(a):
        for b, e in m:
            k = const_of(e)
            if b[0] == 's' and k is not None and k.denominator == 1 and 1 <= k <= MAXEXPAND:
                rest = dict(m)
                del rest[b]
                kk, mm = normalise_mono(rest)
                t = lin({mm: c * kk})
                for _ in range(int(k)):
                    t = lin_mul(t, b[1])
                others = lin({m2: c2 for j, (m2, c2) in enumerate(a) if j != i})
                return expand(lin_add(others, t))
    return a


# ------------------------------------------------------------------ values
_cache = {}


def form(x):
    """canonical linear form of an exact.X value"""
    k = id(x)
    hit = _cache.get(k)
    if hit is not None and hit[0] is x:
        return hit[1]
    if x.kind == "S":
        d = {}
        for r, c in x.terms.items():
            mono = {('p', p): lin_const(e) for p, e in r}
            kk, m = normalise_mono(mono)
            d[m] = d.get(m, Fr(0)) + c * kk
        f = lin(d)
    else:
        t = x.tree
        op = t[0]
        if op == "add":
            f = lin_add(form(t[1]), form(t[2]))
        elif op == "neg":
            f = lin_scale(form(t[1]), -1)
        elif op == "mul":
            f = lin_mul(form(t[1]), form(t[2]))
        elif op == "inv":                    # sum bases are kept POSITIVE: use the true sign
            a = form(t[1])
            try:
                sg = ex.sign(t[1])
            except ex.Unresolved:
                raise Unsupported("sign of an inverted value")
            if sg > 0:
                f = lin_pow(a, lin_const(-1))
            else:
                f = lin_scale(lin_pow(lin_scale(a, -1), lin_const(-1)), -1)
        elif op == "pow":
            f = lin_pow(form(t[1]), form(t[2]))
        else:
            raise Unsupported(op)
    f = expand(f)
    if len(_cache) > 100000:            # bounded: the cache keeps its X objects alive
        _cache.clear()
    _cache[k] = (x, f)
    return f


def equal(x, y):
    """True if x and y are provably equal (equal canonical forms); False if not proved"""
    try:
        return form(x) == form(y)
    except (Unsupported, RecursionError, ZeroDivisionError):
        return False


def is_zero(x):
    try:
        return form(x) == ZERO
    except (Unsupported, RecursionError, ZeroDivisionError):
        return False
