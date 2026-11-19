"""
Independent multi-level interval evaluation of magnitudes built with + |-| * / ^ from exact or
enclosed leaves, for the residual hits of motor X (C5). mpmath intervals (unbounded exponents, outward
rounding) at PREC bits. A positive magnitude u is described by any of
  q    exact rational (Fraction);
  x    an interval for u (whenever |ln u| < ~1e4);
  s, Lam   ln u = s * exp(Lam), s in {-1, 0, +1} (s = 0: u = 1);
  anc  (aq, sig, D): u = aq * (1 + sig * exp(D)), aq exact > 0, exp(D) < 2^-(PREC-60)  (anchored);
  apf  (P, sig, D): u = P * (1 + sig * exp(D)), P = prod p^a an exact power of primes (rational exponents, P irrational),
       exp(D) < 2^-(PREC-60)  (anchored at an exact power; 03/11);
  far  'near1' | 'huge' | 'tiny': only a classification (too deep for any further operation).
Every rule is a rigorous enclosure; when no rule applies, Fail is raised (the hit stays undecided).
"""
from fractions import Fraction as Fr
import math
import mpmath
from mpmath import iv

PREC = 400
iv.prec = PREC
mpmath.mp.prec = PREC
LN2 = 0.6931471805599453
NEAR = -(PREC - 60) * LN2        # Lam below this: |ln u| < 2^-(PREC-60): anchored at 1
XMAX = 9.2                       # Lam below this (|ln u| < ~9900): direct interval x
LCAP = mpmath.mpf("1e600")       # Lam below this: ln u = s e^Lam computable as an interval (mpmath exponents are unbounded; checked up to e^(1e488))
TINYW = iv.mpf(2) ** -(PREC - 58)


def set_prec(p):
    """working precision of every interval (and of the constants that depend on it).  hitcheck.py repeats an
    undecided evaluation with more bits: the same arithmetic, so a hit decided at any precision is decided rigorously"""
    global PREC, NEAR, TINYW
    PREC = int(p)
    iv.prec = PREC
    mpmath.mp.prec = PREC
    NEAR = -(PREC - 60) * LN2
    TINYW = iv.mpf(2) ** -(PREC - 58)
BIGL = mpmath.mpf("1e100000")       # |Lam_base| below this is negligible against |ln y| > e^LCAP
EXACT_BITS = 8_000_000              # exact rationals q^k are expanded up to this size (Python integers; 03/11: was 200 000)


class Fail(Exception):
    pass


class Invalid(Exception):
    pass


def I(a, b=None):
    return iv.mpf([a, a if b is None else b])


def ivq(q):
    return iv.mpf(q.numerator) / iv.mpf(q.denominator)


def lnq(q):
    return iv.log(ivq(q))


class Val:
    __slots__ = ("q", "x", "s", "Lam", "anc", "far", "pf", "off", "apf", "tw")

    def __init__(self, q=None, x=None, s=None, Lam=None, anc=None, far=None, pf=None, off=None, apf=None):
        self.q, self.x, self.s, self.Lam, self.anc, self.far = q, x, s, Lam, anc, far
        self.tw = None        # (11/11) an exact power with an exact irrational exponent: ln u = sum_p sg_p C_p ln p, each
                              # C_p > 0 an exact power of primes; canonical tuple of (p, sg_p, key of C_p) (see tower)
        self.apf = apf        # (03/11) u = P (1 + sig e^D), P an exact power of primes given by its vector (see above)
        self.off = off        # (03/11) u = X + c exactly: (key of the exact power X (pf vector), rational c)
        self.pf = pf          # exact power of primes: canonical {prime: exponent}, exponents integers or Fractions
                              # (03/11: rational exponents, e.g. 3^(1/5); unique since the ln p are independent over Q)

    def __repr__(self):
        if self.far:
            return f"<far {self.far}>"
        if self.q is not None:
            if self.q.denominator < 10**20:
                return f"<q {self.q}>"
            try:
                return f"<q {float(self.q)}>"
            except OverflowError:            # (13/11) beyond the doubles: only the order of magnitude, from the bit lengths
                return f"<q ~2^{self.q.numerator.bit_length() - self.q.denominator.bit_length()}>"
        if self.anc is not None:
            aq, sg, D = self.anc
            return f"<{aq} (1 {'+' if sg > 0 else '-'} e^[{mpmath.nstr(D.a, 8)}, {mpmath.nstr(D.b, 8)}])>"
        if self.x is not None:
            return f"<x [{mpmath.nstr(self.x.a, 12)}, {mpmath.nstr(self.x.b, 12)}]>"
        return f"<ln u = {self.s} e^[{mpmath.nstr(self.Lam.a, 10)}, {mpmath.nstr(self.Lam.b, 10)}]>"


def set_lam(v):
    lam = iv.log(v.x)
    if lam.a > 0:
        v.s, v.Lam = 1, iv.log(lam)
    elif lam.b < 0:
        v.s, v.Lam = -1, iv.log(-lam)
    return v


def exact(q):
    q = Fr(q)
    if q < 0:
        q = -q
    if q == 0:
        return Val(q=q, x=I(0))
    v = Val(q=q, x=ivq(q))
    if q == 1:
        v.s = 0
        return v
    return set_lam(v)


def from_x(x):
    if not (x.a > 0):
        raise Fail("enclosure not positive")
    return set_lam(Val(x=x))


def from_lam(s, Lam):
    """u = exp(s * exp(Lam))"""
    if s == 0:
        return exact(1)
    if Lam.b < NEAR:
        # u - 1 = lam (1 + th), |th| <= |lam| < 2^-(PREC-60): ln|u-1| = Lam + [-2|lam|, 2|lam|]
        return Val(x=1 + I(-TINYW, TINYW), s=s, Lam=Lam, anc=(Fr(1), s, Lam + I(-TINYW, TINYW)))
    v = Val(s=s, Lam=Lam)
    if Lam.b < XMAX:
        v.x = iv.exp(s * iv.exp(Lam))
    return v


def from_lnu(lu):
    if lu.a > 0:
        return from_lam(1, iv.log(lu))
    if lu.b < 0:
        return from_lam(-1, iv.log(-lu))
    if lu.a == 0 and lu.b == 0:
        return exact(1)
    raise Fail("ln u straddles 0")


def make_anc(aq, sg, D):
    """u = aq (1 + sg exp(D))"""
    aq = Fr(aq)
    if D.b >= NEAR:
        return from_x(ivq(aq) * (1 + sg * iv.exp(D)))
    v = Val(anc=(aq, sg, D), x=ivq(aq) * (1 + I(-TINYW, TINYW)))
    if aq == 1:
        v.s, v.Lam = sg, D + I(-TINYW, TINYW)        # ln(1 + t) = t (1 + th), |th| <= |t|
        return v
    return set_lam(v)


# ---- (03/11) first-order anchors at exact powers: rigorous bounds, lvl's own arithmetic (unique factorisation and the
# expansions of ln(1 + x) and e^z - 1); used for the residual hits of n = 9 (tiny differences such as 3 - (3^32768 - 6)^(1/32768))
def lnpf(pf):
    """ln P for P = prod p^a (interval)"""
    r = I(0)
    for p_, a in pf.items():
        a = Fr(a)
        r = r + iv.mpf(a.numerator) / a.denominator * iv.log(iv.mpf(p_))
    return r


def make_apf(pf, sg, D):
    """u = P (1 + sg e^D), P = prod p^a: an exact rational P gives the ordinary anchor; otherwise ln u = ln P + ln(1 + x),
    |ln(1 + x)| <= 2|x| for |x| <= 1/2"""
    pf = {p_: Fr(a) for p_, a in pf.items() if a != 0}
    if not pf:
        return make_anc(Fr(1), sg, D)
    q = pf_value(pf)
    if q is not None:
        return make_anc(q, sg, D)
    lp = lnpf(pf)
    if D.b >= NEAR:
        return from_lnu(lp + iv.log(1 + sg * iv.exp(D)))
    xb = iv.exp(I(D.b)).b
    v = from_lnu(lp + I(-2 * xb, 2 * xb))
    v.apf = (pf, sg, D)
    return v


ANC_MAX = mpmath.mpf(2) ** -20          # first-order rules only while |x| and |z| are below this


def one_pow(sg, D, t):
    """(1 + x)^t = 1 + sg' e^D', x = sg e^D, t a nonzero rational: with |x| <= X <= 2^-20 and z = t ln(1 + x),
    ln(1 + x) = x (1 + th1), |th1| <= |x|;  e^z - 1 = z (1 + th2), |th2| <= |z| <= 2|t|X;  so (1 + x)^t - 1 =
    t x (1 + th1)(1 + th2): sg' = sg sign(t), D' = D + ln|t| + [-2X - 4|t|X, 2X + 4|t|X] (|ln(1 + th)| <= 2|th|).
    None when the bounds are not small enough."""
    X = iv.exp(I(D.b)).b
    tb = abs(Fr(t))
    zb = 2 * mpmath.mpf(tb.numerator) / tb.denominator * X
    if not (X <= ANC_MAX and zb <= ANC_MAX):
        return None
    w = 2 * X + 2 * zb
    return sg * (1 if t > 0 else -1), D + lnq(tb) + I(-w, w)


def iroot(n, r):
    """floor(n^(1/r)) for integers n >= 0, r >= 1 (Newton's method in integers)"""
    if n < 2 or r == 1:
        return n
    x = 1 << -(-n.bit_length() // r)
    while True:
        y = ((r - 1) * x + n // x ** (r - 1)) // r
        if y >= x:
            break
        x = y
    while x ** r > n:
        x -= 1
    while (x + 1) ** r <= n:
        x += 1
    return x


NEAR_BITS = 4_000_000        # near_root: exact rationals up to this size (10/11: was 400 000; 5^262144 + 7 has 608 698 bits)
ANC_HUGE_BITS = 400_000      # add, sub: above this size an exact q +- tiny becomes the anchor q (1 +- tiny/q) (_anc_not_huge)
EXP_GATE_BITS = 400_000      # pw: an exponent known through its vector or offset is used as a rational up to this size
NEAR_MBITS = 64              # the anchor M (or B0) must be a small rational: at most this many bits
NEAR_E_BITS = 2_000_000      # (10/11) exp_near_small: exact exponents up to this size ((4 + 2^-279936)^2 has 559 872 bits)


def near_root(q, t):
    """q^t for an exact rational q > 0 and a non-integer rational t, when q is close to a simple exact value:
      (a) q = M^r (1 + rho), r the denominator of t, M a rational (nearest r-th roots of numerator and
          denominator, or of the integer part of q or of 1/q): q^t = M^a (1 + rho)^t, M^a exact;
      (b) q = B0 (1 + rho), B0 the small rational nearest q (continued fractions): q^t = B0^t (1 + rho)^t, B0^t an
          exact power of primes (make_apf);
    rho tiny (exp(ln|rho|) < 2^-(PREC-60)), (1 + rho)^t by one_pow.  Both are identities; None if not applicable.
    The small anchors (M or B0 of at most NEAR_MBITS bits, the only ones before 10/11) are tried first: (9 + t)^(2/5)
    must be read as 9^(2/5) (1 + x), not as a 5000-bit rational to the 5th.  The others (10/11: M of any size with M^r
    and M^a within the exact sizes, e.g. (4^78125 + 1)^(3/10) with M = 2^15625) only if none of them applies.
    rho = N / Dn is kept as two integers, not reduced (a gcd of numbers of 10^6 bits costs far more than the rest);
    only its sign, whether it is 0, and an enclosure of ln|rho| are used."""
    a, r = t.numerator, t.denominator
    nb = max(q.numerator.bit_length(), q.denominator.bit_length())
    if r == 1 or nb < PREC - 64 or nb > NEAR_BITS:
        return None

    def root_cand(M):                                    # q / M^r - 1
        Mn, Md = M.numerator ** r, M.denominator ** r
        return ("root", M, q.numerator * Md - q.denominator * Mn, q.denominator * Mn)

    def attempt(cands):
        # by log2|rho|, roughly (bit lengths)
        for kind, M, N_, Dn in sorted(cands, key=lambda c: -10 ** 9 if c[2] == 0 else abs(c[2]).bit_length() - c[3].bit_length()):
            if N_ == 0:
                if kind == "root":
                    return exact(M ** a if a > 0 else 1 / M ** (-a))
                continue                                 # q = B0 exactly: the exact rules apply (pf of B0)
            D = iv.log(iv.mpf(abs(N_)) / iv.mpf(Dn))
            if D.b >= NEAR:
                continue
            op_ = one_pow(1 if N_ > 0 else -1, D, t)
            if op_ is None:
                continue
            if kind == "root":
                return make_anc(M ** a if a > 0 else 1 / M ** (-a), op_[0], op_[1])
            vec = pfactor(M)
            if vec is None:
                continue
            return make_apf({p_: Fr(x) * t for p_, x in vec.items()}, op_[0], op_[1])
        return None

    small, large = [], []
    if r <= nb:                                          # (an r-th root below 2 is no anchor; also bounds M^r below)
        rn, rd = iroot(q.numerator, r), iroot(q.denominator, r)
        pairs = [(num, den) for num in (rn, rn + 1) for den in (rd, rd + 1)]
        # (10/11) also the r-th roots of the integer part of q, or of 1/q: q = M^r + c with c not an integer
        # (4^78125 + 1/10 = (10 * 2^156250 + 1) / 10, whose numerator and denominator are no 10th powers)
        ri = iroot(q.numerator // q.denominator, r) if q > 1 else 0
        rj = iroot(q.denominator // q.numerator, r) if q < 1 else 0
        pairs += [(ri, 1), (ri + 1, 1)] if ri else []
        pairs += [(1, rj), (1, rj + 1)] if rj else []
        for num, den in dict.fromkeys(pairs):
            if num and den:
                M = Fr(num, den)
                mb = max(M.numerator.bit_length(), M.denominator.bit_length())
                if r * mb <= 2 * NEAR_BITS and abs(a) * mb <= EXACT_BITS:
                    if mb <= NEAR_MBITS:
                        small.append(root_cand(M))
                    else:
                        large.append(M)
    B0 = q.limit_denominator(2 ** NEAR_MBITS)
    if B0 > 0 and B0.numerator.bit_length() <= NEAR_MBITS:                    # q / B0 - 1
        small.append(("near", B0, q.numerator * B0.denominator - q.denominator * B0.numerator, q.denominator * B0.numerator))
    r_ = attempt(small)
    if r_ is None and large:
        r_ = attempt([root_cand(M) for M in dict.fromkeys(large)])
    return r_


def exp_near_small(base, e, sgn):
    """B^(sgn e) for an exact base B (rational or exact power, B != 1) and an exact rational exponent e of many bits that
    is close to a small rational c: e = c (1 + rho), rho tiny.  The value is the exact power prod p^(a sgn e) (kept as
    its vector) and also, by the anchored-exponent rule, B^(sgn c) (1 + x) with the bounds of pow_anchored."""
    q = Fr(e.q)
    nb = max(q.numerator.bit_length(), q.denominator.bit_length())
    if nb < PREC - 64 or nb > NEAR_E_BITS:
        return None
    c = q.limit_denominator(2 ** NEAR_MBITS)
    if c == 0 or c.numerator.bit_length() > NEAR_MBITS:
        return None
    rho = q / c - 1
    if rho == 0:
        return None
    D = lnq(abs(rho))
    if D.b >= NEAR:
        return None
    ea = Val(anc=(c, 1 if rho > 0 else -1, D))         # e as an anchored value (only its anchor is read)
    r_ = pow_anchored(base, ea, sgn)
    if r_ is None:
        return None
    vec = pf_of(base)
    if vec is not None and r_.q is None:
        r_.pf = {p_: Fr(x) * q * sgn for p_, x in vec.items()}   # the exact value, as a vector (equalities, same_pf)
    return r_


def pow_anchored(base, e, sgn):
    """the first-order rules of pw: an anchored base (anc or apf) to an exact rational, an exact base (rational or
    exact power) to an anchored exponent, an exact rational base close to an exact r-th power (near_root).  None if no
    rule applies (the caller then uses the plain enclosures)."""
    if e.q is not None and e.q != 0:
        t = Fr(e.q) * sgn
        if base.anc is not None or base.apf is not None:
            if base.anc is not None:
                aq, sg, D = base.anc
                vec = pfactor(aq) if t.denominator != 1 else None
                if t.denominator == 1:
                    k = int(t)
                    if abs(k) * max(aq.numerator.bit_length(), aq.denominator.bit_length()) >= EXACT_BITS:
                        return None
                    op_ = one_pow(sg, D, t)
                    return None if op_ is None else make_anc(aq ** k, op_[0], op_[1])
                if vec is None:
                    return None
            else:
                vec, sg, D = base.apf
            op_ = one_pow(sg, D, t)
            if op_ is None:
                return None
            return make_apf({p_: Fr(x) * t for p_, x in vec.items()}, op_[0], op_[1])
        if base.q is not None and t.denominator != 1:
            return near_root(base.q, t)
        return None
    if e.anc is not None and (base.q is not None or base.pf is not None) and not (base.q is not None and base.q == 1):
        # B^(c (1 + x)) = B^c e^w, w = c x ln B, |w| <= W <= 2^-20: e^w - 1 = w (1 + th), |th| <= |w|
        vec = pf_of(base)
        if vec is None:
            return None
        aq, sg, D = e.anc
        c = Fr(aq) * sgn
        lb = lnpf(vec)
        if lb.a <= 0 <= lb.b:
            return None
        X = iv.exp(I(D.b)).b
        W = mpmath.mpf(abs(c.numerator)) / c.denominator * X * abs(lb).b
        if not (W <= ANC_MAX):
            return None
        sgw = sg * (1 if c > 0 else -1) * (1 if lb.a > 0 else -1)
        D2 = D + lnq(abs(c)) + iv.log(abs(lb)) + I(-2 * W, 2 * W)
        return make_apf({p_: Fr(x) * c for p_, x in vec.items()}, sgw, D2)
    return None


def lam_iv(v):
    """ln u as an interval"""
    if v.far:
        raise Fail("far value used as an operand")
    if v.q is not None and v.q == 1:
        return I(0)
    if v.s is None:
        raise Fail("sign of ln u unknown")
    if v.Lam.b > LCAP:
        raise Fail("ln u too large for an interval")
    return v.s * iv.exp(v.Lam)


def is_tiny(v):
    return v.s == -1 and v.Lam is not None and v.Lam.a > XMAX and v.q is None


def chk(*vs):
    for v in vs:
        if v.far:
            raise Fail("far value used as an operand")


def ln1mexp(d):
    """ln(1 - exp(-d)) for an interval d > 0"""
    if not (d.a > 0):
        raise Fail("ln1mexp: d not positive")
    if d.a > 0.5:
        return iv.log(1 - iv.exp(-d))
    lo = iv.log(I(d.a)) + iv.log(1 - I(d.a) / 2)    # 1 - e^-d >= d - d^2/2
    hi = iv.log(1 - iv.exp(-I(d.b))) if d.b > 0.5 else iv.log(I(d.b))
    return I(lo.a, hi.b)


def by_logs(la, lb, plus):
    """ln(e^la + e^lb) (any la, lb), or ln|e^la - e^lb| (la, lb separated)"""
    if plus:
        hi, lo = (la, lb) if la.b >= lb.b else (lb, la)
        return hi + iv.log(1 + iv.exp(lo - hi))          # identity, evaluated in interval arithmetic
    if la.a > lb.b:
        hi, lo = la, lb
    elif lb.a > la.b:
        hi, lo = lb, la
    else:
        raise Fail("difference of values with overlapping logs")
    return hi + ln1mexp(hi - lo)


def add(a, b):
    chk(a, b)
    if a.q == 0:
        return b
    if b.q == 0:
        return a
    if a.q is not None and b.q is not None:
        r_ = _anc_not_huge(a, b, a.q, b.q, 1)          # (03/11) q + t, t tiny: the anchor, not a huge rational
        if r_ is not None:
            return r_
        return exact(a.q + b.q)
    for e, t in ((a, b), (b, a)):
        ft, fe = off_form(t), off_form(e)
        if ft is not None and (fe is None or fe[0] != ft[0]):
            # (X + sg c) + e = X + (sg c + e), exactly (03/11): for sg = +1 the new offset is c + e; for sg = -1 it is
            # known when c and e are exact rationals
            k, sg, c = ft
            r = _add_core(a, b)
            if sg > 0:
                r.off = (k, 1, add(c, e))
            elif c.q is not None and e.q is not None:
                d = e.q - c.q
                r.off = (k, 1 if d >= 0 else -1, exact(abs(d)))
            return r
    return _add_core(a, b)


def _apf_tiny(e, t):
    """D = ln(t / P) when e is an irrational exact power P (pf) and t is tiny relative to it (exp(D) < 2^-(PREC-60)); else None"""
    if e.pf is None or e.q is not None or t.q is not None and t.q == 0:
        return None
    try:
        D = lam_iv(t) - lnpf(e.pf)
    except Fail:
        return None
    return D if D.b < NEAR else None


def _rel_tiny(q, t):
    """D = ln(t / q) when t is tiny relative to the exact rational q > 0 (exp(D) < 2^-(PREC-60)); else None (03/11: the
    anchor q (1 +- t/q) applies whenever t/q is that small, not only when |ln t| > e^XMAX)"""
    if q == 0 or (t.q is not None) or t.anc is not None or t.apf is not None:
        return None
    try:
        D = lam_iv(t) - lnq(q)
    except Fail:
        return None
    return D if D.b < NEAR else None


def _add_core(a, b):
    for e, t in ((a, b), (b, a)):
        D = _apf_tiny(e, t)
        if D is not None:
            return make_apf(e.pf, 1, D)                  # P + t = P (1 + t/P)
        if e.q is None and e.apf is not None and t.q is not None and t.q == 0:
            return e
    for e, t in ((a, b), (b, a)):
        if e.q is not None:
            if t.anc is not None:
                aq, sg, D = t.anc
                q2 = e.q + aq
                return make_anc(q2, sg, D + lnq(aq) - lnq(q2))
            if is_tiny(t):
                return make_anc(e.q, 1, lam_iv(t) - lnq(e.q))
            D = _rel_tiny(e.q, t)
            if D is not None:
                return make_anc(e.q, 1, D)
    if a.x is not None and b.x is not None:
        return from_x(a.x + b.x)
    return from_lnu(by_logs(lam_iv(a), lam_iv(b), True))


_SMALLP = None


def _small_factor(n):
    """{prime: exponent} of an integer n >= 1 made of primes below 10^4 only (trial division, exact); None otherwise"""
    global _SMALLP
    if _SMALLP is None:
        from sympy import primerange
        _SMALLP = [int(p_) for p_ in primerange(2, 10 ** 4)]
    out = {}
    for p_ in _SMALLP:
        if n == 1:
            break
        if n % p_ == 0:
            k = 0
            while n % p_ == 0:                           # divide by p^(2^j) while possible: few big divisions
                pk, e = p_, 1
                while n % (pk * pk) == 0:
                    pk *= pk; e *= 2
                n //= pk; k += e
            out[p_] = k
    return out if n == 1 else None


_PFC = {}


def pfactor(q):
    """{prime: exponent} of a positive rational (cached)"""
    key = (q.numerator, q.denominator) if q > 0 else None
    if key is None:
        return None
    if key not in _PFC:
        if len(_PFC) > 20000:
            _PFC.clear()
        _PFC[key] = _pfactor(q)
    r = _PFC[key]
    return None if r is None else dict(r)


def _pfactor(q):
    """{prime: exponent} of a positive rational: sympy for numerator and denominator up to 64 bits; above, only
    numbers made of primes below 10^4 (trial division; 03/11, for large exact powers such as 3^125); None otherwise"""
    if q <= 0:
        return None
    out = {}
    for n, sg in ((q.numerator, 1), (q.denominator, -1)):
        if n.bit_length() <= 64:
            from sympy import factorint
            f = {int(p_): int(a) for p_, a in factorint(n).items()}
        else:
            f = _small_factor(n)
            if f is None:
                return None
        for p_, a in f.items():
            out[p_] = out.get(p_, 0) + sg * a
    return {p_: a for p_, a in out.items() if a}


def exact_of(v):
    """the exact rational value of v when it is known and expandable: v.q, an exact power (pf), or X +- c with X
    expandable and c exact; None otherwise"""
    if v.q is not None:
        return v.q
    if v.pf is not None:
        return pf_value(v.pf)
    if v.off is not None:
        k, sg, c = v.off
        if c.q is None:
            return None
        qx = pf_value(dict(k))
        return None if qx is None else qx + sg * c.q
    return None


_PFV = {}


def pf_value(pf):
    """the exact rational of a prime-exponent vector with integer exponents (cached; None if a fractional exponent,
    or larger than EXACT_BITS)"""
    key = pf_key(pf)
    if key not in _PFV:
        if len(_PFV) > 2000:
            _PFV.clear()
        _PFV[key] = _pf_value(pf)
    return _PFV[key]


def _pf_value(pf):
    """the exact rational of a prime-exponent vector with integer exponents (None if a fractional exponent, or larger
    than EXACT_BITS)"""
    if any(Fr(a).denominator != 1 for a in pf.values()):
        return None
    if any(abs(int(a)).bit_length() > 40 for a in pf.values()):
        return None                                      # (03/11) far beyond EXACT_BITS; also avoids float overflow
    if sum(abs(int(a)) * math.log2(p_) for p_, a in pf.items()) > EXACT_BITS - 64:
        return None
    q = Fr(1)
    for p_, a in pf.items():
        q *= Fr(p_) ** int(a)
    return q


def same_pf(a, b):
    """True/False when both values are exact powers of primes (a vector pf, or an exact rational); None otherwise.
    A vector with a fractional exponent is irrational, so it differs from every rational; a vector with integer
    exponents is compared with a large rational through its exact value."""
    fa = a.pf if a.pf is not None else (pfactor(a.q) if a.q is not None else None)
    fb = b.pf if b.pf is not None else (pfactor(b.q) if b.q is not None else None)
    if fa is not None and fb is not None:
        return {p_: Fr(x) for p_, x in fa.items()} == {p_: Fr(x) for p_, x in fb.items()}
    for v, f in ((a, fb), (b, fa)):                      # one side a large exact rational, the other a vector
        if v.q is not None and v.pf is None and f is not None:
            qv = pf_value(f)
            if qv is None:
                return False if any(Fr(x).denominator != 1 for x in f.values()) else None
            return qv == v.q
    return None


def pf_of(v):
    """prime-exponent vector of an exact value (its pf, or the factorisation of a small exact rational), else None"""
    if v.pf is not None:
        return {p_: Fr(x) for p_, x in v.pf.items()}
    if v.q is not None and v.q > 0:
        f = pfactor(v.q)
        return None if f is None else {p_: Fr(x) for p_, x in f.items()}
    return None


def pf_combine(fa, fb, sg):
    out = dict(fa)
    for p_, x in fb.items():
        out[p_] = out.get(p_, Fr(0)) + sg * x
    return {p_: x for p_, x in out.items() if x != 0}


def pf_key(pf):
    return tuple(sorted((p_, Fr(a)) for p_, a in pf.items()))


ZERO_VAL = None


def off_form(v):
    """(key of X, sg, c) with v = X + sg c exactly: X an exact power kept as a vector (not expanded), sg = +-1 and c a
    value of lvl (a magnitude, exact or an enclosure); None otherwise"""
    global ZERO_VAL
    if v.off is not None:
        return v.off
    if v.pf is not None and v.q is None:
        if ZERO_VAL is None:
            ZERO_VAL = exact(0)
        return pf_key(v.pf), 1, ZERO_VAL
    return None


def pf_result(pf, approx):
    """the value with exact vector pf: an exact rational if it has integer exponents and fits, else approx with pf"""
    if not pf:
        return exact(1)
    q = pf_value(pf)
    if q is not None:
        return exact(q)
    approx.pf = pf
    return approx


def _bits(q):
    return max(q.numerator.bit_length(), q.denominator.bit_length())


def _anc_not_huge(a, b, qa, qb, sg=-1):
    """|a - b| (sg = -1) or a + b (sg = +1) of two exact values whose result would be a rational of more than ANC_HUGE_BITS bits, one of them (q, of at
    most NEAR_MBITS bits) and the other (t) tiny relative to it: the anchor q (1 - t/q) instead of the huge rational,
    with the exact statement kept as the offset X - t, X = q as a vector (03/11).  None otherwise."""
    if max(qa.numerator.bit_length(), qa.denominator.bit_length(), qb.numerator.bit_length(), qb.denominator.bit_length()) <= ANC_HUGE_BITS:
        return None
    for (e, qe), (t, qt) in (((a, qa), (b, qb)), ((b, qb), (a, qa))):
        if qe == 0 or max(qe.numerator.bit_length(), qe.denominator.bit_length()) > NEAR_MBITS:
            continue
        lt = lnq(qt) if qt != 0 else None
        if lt is None:
            continue
        D = lt - lnq(qe)
        if D.b >= NEAR:
            continue
        r = make_anc(qe, sg, D)
        f = pfactor(qe)
        if f is not None and f:
            r.off = (pf_key(f), sg, t)
        return r
    return None


def sub(a, b):
    """|a - b|"""
    chk(a, b)
    if (a.pf is not None or b.pf is not None) and same_pf(a, b):
        return exact(0)                                  # equal exact rationals (unique factorisation)
    if a.tw is not None and a.tw == b.tw:
        return exact(0)                                  # (11/11) the same exact power written twice (tower)
    if not (a.q is not None and b.q is not None) and (a.pf is not None or b.pf is not None or a.off is not None or b.off is not None):
        # both values exactly known and small enough to expand (an exact rational, an exact power, or X +- c with c
        # exact): the exact difference (03/11)
        qa, qb = exact_of(a), exact_of(b)
        if qa is not None and qb is not None:
            r_ = _anc_not_huge(a, b, qa, qb)
            if r_ is not None:
                return r_
            return exact(abs(qa - qb))
    fa, fb = off_form(a), off_form(b)
    if fa is not None and fb is not None and fa[0] == fb[0]:
        # |(X + s1 c1) - (X + s2 c2)| = |s1 c1 - s2 c2|: |c1 - c2| for equal signs, c1 + c2 otherwise (03/11)
        return sub(fa[2], fb[2]) if fa[1] == fb[1] else add(fa[2], fb[2])
    for e, t, f in ((a, b, fb), (b, a, fa)):
        fe_ = off_form(e)
        if f is not None and (fe_ is None or fe_[0] != f[0]) and not (e.q is not None and e.q == 0):
            # |e - (X + sg c)| with X + sg c > e certainly (by the logs): X + sg c - e
            lt, le = lam_iv(t), lam_iv(e)
            if lt.a > le.b + 2:
                k, sg, c = f
                D_ = _apf_tiny(t, e) if (sg > 0 and c.q is not None and c.q == 0) else None
                r = make_apf(t.pf, -1, D_) if D_ is not None else from_lnu(by_logs(lt, le, False))   # (03/11) also the anchor
                if sg < 0:
                    r.off = (k, -1, add(c, e))
                elif c.q is not None and c.q == 0:
                    r.off = (k, -1, e)                         # X - e
                elif c.q is not None and e.q is not None:
                    d = c.q - e.q
                    r.off = (k, 1 if d >= 0 else -1, exact(abs(d)))
                return r
    if a.q == 0:
        return b
    if b.q == 0:
        return a
    for e, t in ((a, b), (b, a)):
        D = _apf_tiny(e, t)
        if D is not None:                                # P - t = P (1 - t/P), t tiny relative to P (03/11)
            r = make_apf(e.pf, -1, D)
            r.off = (pf_key(e.pf), -1, t)
            return r
    for e, t in ((a, b), (b, a)):
        if t.apf is not None and (e.pf is not None and e.q is None) and pf_key(e.pf) == pf_key(t.apf[0]):
            return from_lnu(lnpf(e.pf) + t.apf[2])       # |P - P (1 + x)| = P e^D (03/11)
    if a.apf is not None and b.apf is not None and pf_key(a.apf[0]) == pf_key(b.apf[0]):
        P, s1, D1 = a.apf
        _, s2, D2 = b.apf
        return from_lnu(lnpf(P) + by_logs(D1, D2, s1 != s2))
    if a.q is not None and b.q is not None:
        r_ = _anc_not_huge(a, b, a.q, b.q)
        if r_ is not None:
            return r_
        return exact(abs(a.q - b.q))
    for e, t in ((a, b), (b, a)):
        if e.q is not None:
            if t.anc is not None:
                aq, sg, D = t.anc
                if e.q == aq:
                    return from_lnu(lnq(aq) + D)                  # aq e^D
                d = e.q - aq
                q2 = abs(d)
                return make_anc(q2, -sg if d > 0 else sg, D + lnq(aq) - lnq(q2))
            if is_tiny(t):
                return make_anc(e.q, -1, lam_iv(t) - lnq(e.q))
            D = _rel_tiny(e.q, t)
            if D is not None:
                return make_anc(e.q, -1, D)
    if a.anc is not None and b.anc is not None and a.anc[0] == b.anc[0]:
        aq, s1, D1 = a.anc
        _, s2, D2 = b.anc
        return from_lnu(lnq(aq) + by_logs(D1, D2, s1 != s2))
    if a.x is not None and b.x is not None:
        d = a.x - b.x
        if d.a > 0:
            return from_x(d)
        if d.b < 0:
            return from_x(-d)
    return from_lnu(by_logs(lam_iv(a), lam_iv(b), False))


def dominant(a, b, sgn):
    """ln u = ln a + sgn ln b when one log is beyond intervals (Lam > LCAP): ln u = A (1 + th), |th| <= r"""
    for big, small, bigsign in ((a, b, 1), (b, a, sgn)):
        if big.s is not None and big.Lam is not None and big.Lam.a > LCAP:
            ls = abs(lam_iv(small))
            if ls.b == 0:
                return from_lam(bigsign * big.s, big.Lam)
            r = iv.exp(iv.log(ls) - big.Lam)
            if r.b < 0.1:
                return from_lam(bigsign * big.s, big.Lam + iv.log(1 + I(-r.b, r.b)))
    raise Fail("product of two very large logs")


def mul(a, b):
    chk(a, b)
    if a.q == 0 or b.q == 0:
        return exact(0)
    if a.q is not None and b.q is not None:
        return exact(a.q * b.q)
    for e, t in ((a, b), (b, a)):
        if e.q is not None and t.anc is not None:
            aq, sg, D = t.anc
            return make_anc(e.q * aq, sg, D)
        if t.apf is not None and (e.q is not None or e.pf is not None) and t.q is None:
            fe = pf_of(e)
            if fe is not None:                            # (03/11) Q P (1 + x) = (Q P)(1 + x)
                P, sg, D = t.apf
                return make_apf(pf_combine(P, fe, 1), sg, D)
    pfr = None
    if a.pf is not None or b.pf is not None:            # exact powers of primes: the vectors add
        fa, fb = pf_of(a), pf_of(b)
        if fa is not None and fb is not None:
            pfr = pf_combine(fa, fb, 1)
            if not pfr or pf_value(pfr) is not None:
                return pf_result(pfr, None)              # exactly 1, or an exact rational that fits
    try:
        r = from_lnu(lam_iv(a) + lam_iv(b))
    except Fail:
        r = dominant(a, b, 1)
    return r if pfr is None else pf_result(pfr, r)


def div(a, b):
    """a / b"""
    chk(a, b)
    if b.q == 0:
        raise Invalid("division by 0")
    if a.q == 0:
        return exact(0)
    if a.q is not None and b.q is not None:
        return exact(a.q / b.q)
    if a.anc is not None and b.q is not None:
        aq, sg, D = a.anc
        return make_anc(aq / b.q, sg, D)
    if a.apf is not None and (b.q is not None or b.pf is not None):
        fb_ = pf_of(b)
        if fb_ is not None:                               # (03/11) P (1 + x) / Q
            P, sg, D = a.apf
            return make_apf(pf_combine(P, fb_, -1), sg, D)
    if a.q is not None and b.anc is not None:
        aq, sg, D = b.anc
        return make_anc(a.q / aq, -sg, D + I(-TINYW, TINYW))      # 1/(1+t) = 1 - t(1+th)
    pfr = None
    if a.pf is not None or b.pf is not None:            # exact powers of primes: the vectors subtract
        fa, fb = pf_of(a), pf_of(b)
        if fa is not None and fb is not None:
            pfr = pf_combine(fa, fb, -1)
            if not pfr or pf_value(pfr) is not None:
                return pf_result(pfr, None)              # exactly 1, or an exact rational that fits
    try:
        r = from_lnu(lam_iv(a) - lam_iv(b))
    except Fail:
        r = dominant(a, b, -1)
    return r if pfr is None else pf_result(pfr, r)


def tower(vec, E, sgn):
    """(11/11) the canonical form of P^(sgn E) for P = prod p^(a_p) (vec) and E = prod r^(x_r) > 0 (a vector of prime
    exponents, the x_r rational): ln of it is sum_p (sgn a_p E) ln p, and each |a_p| E is again a product of prime
    powers (the factorisation of the rational |a_p| added to E).  Tuple of (p, sign, key of |a_p| E), sorted; None
    when some |a_p| cannot be factored.  Two values with the same tuple have the same ln, written term by term with
    the same exact numbers, so they are equal (the converse is not used)."""
    out = []
    for p_, a in vec.items():
        a = Fr(a)
        if a == 0:
            continue
        fa = pfactor(abs(a))
        if fa is None:
            return None
        out.append((p_, sgn * (1 if a > 0 else -1), pf_key(pf_combine(E, {r_: Fr(x) for r_, x in fa.items()}, 1))))
    return tuple(sorted(out)) if out else None


def pw(base, e, sgn):
    """base ^ (sgn * e), base and e magnitudes.  (11/11) The result of an exact power of primes raised to an exact
    power of primes that is not a usable rational (2^(35/2), 3^8.4) also carries its canonical form (tower); and a
    value with such a form raised to an exact rational k keeps one (its ln is multiplied by sgn k)."""
    r = _pw(base, e, sgn)
    if r.q is None and r.far is None and r.tw is None:
        try:
            if e.q is None and e.pf is not None:
                vec = pf_of(base)
                if vec:
                    r.tw = tower(vec, {p_: Fr(x) for p_, x in e.pf.items()}, sgn)
            elif e.q is not None and e.q > 0 and base.tw is not None:
                fk = pfactor(Fr(e.q))
                if fk is not None:
                    fk = {p_: Fr(x) for p_, x in fk.items()}
                    r.tw = tuple(sorted((p_, sg * sgn, pf_key(pf_combine(dict(k_), fk, 1))) for (p_, sg, k_) in base.tw))
        except Exception:
            r.tw = None                                  # the form is optional: without it nothing else changes
    return r


def _pw(base, e, sgn):
    """base ^ (sgn * e), base and e magnitudes"""
    chk(base, e)
    # (03/11) an exponent or base whose exact rational value is known only through its prime-exponent vector (pf with
    # integer exponents, e.g. (10^(1/2))^6 = 1000) or its exact offset is used as that rational (exact_of: the value
    # is exactly that number; None when unknown or too large, and then nothing changes)
    if e.q is None:
        q_ = exact_of(e)
        if q_ is not None and _bits(q_) <= EXP_GATE_BITS:
            e = exact(q_)
    if base.q is None:
        q_ = exact_of(base)
        if q_ is not None and _bits(q_) <= NEAR_BITS and base.anc is None and base.apf is None:
            base = exact(q_)                             # (an anchored base keeps its anchor)
    if e.q == 0:
        if base.q == 0:
            raise Invalid("0^0")
        return exact(1)
    if base.q == 0:
        if sgn > 0:
            return exact(0)
        raise Invalid("0^-x")
    if base.q == 1:
        return exact(1)
    if base.q is not None and e.q is not None and e.q.denominator == 1:
        k = int(e.q)
        if k * max(base.q.numerator.bit_length(), base.q.denominator.bit_length()) < EXACT_BITS:
            return exact(base.q ** (sgn * k))
        pf = pfactor(base.q)
        if pf is not None:                               # exact, kept factorised: q^(sgn k)
            v = from_lam(1 if sgn * k * lnq(base.q).a > 0 else -1, iv.log(abs(sgn * k * lnq(base.q))))
            v.pf = {p_: a * sgn * k for p_, a in pf.items()}
            return v
    if (base.anc is not None or base.apf is not None) and e.q is not None and e.q != 0:
        # (03/11) an anchored base: the first-order rule, keeping the exact vector too when the base has one
        r_ = pow_anchored(base, e, sgn)
        if r_ is not None:
            if base.pf is not None and r_.q is None and r_.pf is None:
                r_.pf = {p_: Fr(a) * Fr(e.q) * sgn for p_, a in base.pf.items()}
            return r_
    if base.pf is not None and e.q is not None and e.q.denominator == 1:
        k = int(e.q)
        v = from_lnu(sgn * k * lam_iv(base))
        v.pf = {p_: a * sgn * k for p_, a in base.pf.items()}
        return v
    if e.q is not None and e.q.denominator != 1 and (base.pf is not None or base.q is not None):
        r_ = exp_near_small(base, e, sgn)               # (03/11) e = c (1 + rho), c a small rational, rho tiny
        if r_ is not None:
            return r_
    if e.q is not None and e.q.denominator != 1 and (base.pf is not None or base.q is not None):
        # (03/11) an exact power of primes to an exact rational t: (prod p^a)^t = prod p^(a t), real and positive
        fb = pf_of(base)
        if fb is not None:
            t = Fr(e.q) * sgn
            v = from_lnu(t.numerator * lam_iv(base) / t.denominator) if base.pf is not None else None
            if v is None:
                if base.s is None:
                    raise Fail("base near 1 without an anchor")
                v = from_lam(sgn * base.s, lam_iv(e) + base.Lam) if base.Lam is not None else from_lnu(t.numerator * lam_iv(base) / t.denominator)
            return pf_result({p_: a * t for p_, a in fb.items()}, v)
    r_ = pow_anchored(base, e, sgn)                     # (03/11) first-order anchors
    if r_ is not None:
        return r_
    if base.s is None:
        raise Fail("base near 1 without an anchor")
    s_u = sgn * base.s
    try:
        ly = lam_iv(e)                                   # ln y
    except Fail:
        # |ln y| > e^LCAP while |ln|ln base|| < BIGL: ln|ln u| = ln y + Lam_base is dominated by ln y
        if e.s is not None and e.Lam is not None and e.Lam.a > LCAP and base.Lam.a > -BIGL and base.Lam.b < BIGL:
            if e.s == 1:
                return Val(far="huge" if s_u > 0 else "tiny")
            return Val(far="near1")
        raise
    return from_lam(s_u, ly + base.Lam)


FORMS = ["s+l", "|s-l|", "s*l", "s/l", "l/s", "s^l", "s^-l", "l^s", "l^-s"]


def op(f, s, l):
    if f == "s+l": return add(s, l)
    if f == "|s-l|": return sub(s, l)
    if f == "s*l": return mul(s, l)
    if f == "s/l": return div(s, l)
    if f == "l/s": return div(l, s)
    if f == "s^l": return pw(s, l, 1)
    if f == "s^-l": return pw(s, l, -1)
    if f == "l^s": return pw(l, s, 1)
    if f == "l^-s": return pw(l, s, -1)
    raise ValueError(f)


def outside(v, lo, hi):
    """True if the magnitude v is certainly outside [lo, hi] (0 < lo < 1 < ... or any positive window
    with |ln| < 10); False if undecided"""
    lo, hi = I(lo), I(hi)
    if v.far:
        if v.far == "near1":
            return not (lo.a <= 1.0000001 and hi.b >= 0.9999999)
        return True
    if v.q is not None:
        return ivq(v.q).b < lo.a or ivq(v.q).a > hi.b
    if v.x is not None and (v.x.b < lo.a or v.x.a > hi.b):
        return True
    if v.s is not None and v.Lam is not None and v.Lam.a > 2.31 and lo.a > 1e-4 and hi.b < 1e4:
        return True                                        # |ln u| > 10
    return False
