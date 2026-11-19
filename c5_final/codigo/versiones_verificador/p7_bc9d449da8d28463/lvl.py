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
    __slots__ = ("q", "x", "s", "Lam", "anc", "far", "pf", "off", "apf", "tw", "lin", "bnd")

    def __init__(self, q=None, x=None, s=None, Lam=None, anc=None, far=None, pf=None, off=None, apf=None):
        self.q, self.x, self.s, self.Lam, self.anc, self.far = q, x, s, Lam, anc, far
        self.bnd = None       # (17/11) far == "bound": (ln of a lower bound or None for 0, ln of an upper bound) of u (see bound_diff)
        self.lin = None       # (13/11) an exact linear form: u = c0 + sum_k a_k P_k (see lin_of)
        self.tw = None        # (11/11) an exact power with an exact irrational exponent: ln u = sum_p sg_p C_p ln p, each
                              # C_p > 0 an exact power of primes; canonical tuple of (p, sg_p, key of C_p) (see tower)
        self.apf = apf        # (03/11) u = P (1 + sig e^D), P an exact power of primes given by its vector (see above)
        self.off = off        # (03/11) u = X + c exactly: (key of the exact power X (pf vector), rational c)
        self.pf = pf          # exact power of primes: canonical {prime: exponent}, exponents integers or Fractions
                              # (03/11: rational exponents, e.g. 3^(1/5); unique since the ln p are independent over Q)

    def __repr__(self):
        if self.far == "bound":
            lo, hi = self.bnd
            return f"<bound: u in [{'0' if lo is None else 'e^' + mpmath.nstr(lo, 8)}, e^{mpmath.nstr(hi, 8)}]>"
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
        # (17/11) the r-th power part of a denominator (or numerator) made of small primes, and the r-th root of the
        # rest: (10/3)^78125 + 4/5 = (5 10^78125 + 4 3^78125) / (5 3^78125) is M^5 (1 + rho) with M = (10/3)^15625
        # (splits 146 and 177), which none of the roots above gives
        for big, small_, num_side in ((q.denominator, q.numerator, False), (q.numerator, q.denominator, True)):
            if big.bit_length() > 64:
                fb = _small_factor(big)
                if fb:
                    pr = 1
                    for p_, x in fb.items():
                        pr *= p_ ** (x // r)
                    if pr > 1:
                        rest = iroot(small_ * pr ** r // big, r)
                        if rest:
                            pairs += [(pr, rest), (pr, rest + 1)] if num_side else [(rest, pr), (rest + 1, pr)]
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


def anchor_of_q(q):
    """(14/11) an exact rational q of many bits (too large to raise to a power, and not factorable) written exactly as
    M (1 + sg e^D): M = the best rational approximation of q with denominator below 2^31, d = q / M - 1 computed exactly,
    D = ln|d| (an interval).  Only when d is tiny (D below NEAR), so that the anchored rules apply; None otherwise.
    Split 100 of n = 9: (9 + 10^-65536)^64 against 3^128 = 9^64."""
    if q <= 0 or max(q.numerator.bit_length(), q.denominator.bit_length()) <= NEAR_MBITS:
        return None
    e2 = q.numerator.bit_length() - q.denominator.bit_length() - 63
    if e2 > 64:
        # (14/11) q huge and dominated by X = m 2^e (m of 63 bits): q = X (1 + d) with X kept as its vector (never
        # expanded) and d = (q - X)/X exactly; (2^4782969 + c)^2 of split 117
        m = round(q / (Fr(2) ** e2))
        if m > 0:
            fm = pfactor(Fr(m))
            if fm is not None:
                vec = {p_: Fr(x) for p_, x in fm.items()}
                vec[2] = vec.get(2, Fr(0)) + e2
                diff = q - Fr(m) * Fr(2) ** e2
                if diff == 0:
                    return pf_result(vec, from_lnu(lnpf(vec)))
                try:
                    D = lnq(abs(diff)) - lnpf(vec)
                except Exception:
                    D = None
                if D is not None and D.b < NEAR:
                    return make_apf(vec, 1 if diff > 0 else -1, D)
    M = Fr(q).limit_denominator(2 ** 31)
    if M <= 0:
        return None
    d = q / M - 1
    if d == 0:
        return exact(M)
    try:
        D = lnq(abs(d))
    except Exception:
        return None
    if not (D.b < NEAR):
        return None
    return make_anc(M, 1 if d > 0 else -1, D)


def _root_of_offset(base, t):
    """(17/11) (a P + c)^t, t not an integer, for an exact rational base = a P + c given by its exact linear form: one
    exact power P (prime vector), a rational a > 0, a rational c != 0 tiny relative to a P; a P + c is checked to be
    the exact rational of the base (exact integers).  Then base = a P (1 + x),
    x = c / (a P) = sg e^D, and (1 + x)^t = 1 + sg' e^D' (one_pow): an anchor at the exact power (a P)^t, which stays
    exact under later integer powers and differences (((3^390625 + 3/5)^(1/10))^4 against 9^78125, splits 92 and 123),
    where near_root anchors at a rational close to the irrational root.  None when it does not apply."""
    f = lin_of(base)
    if f is None or len(f[1]) != 1 or f[0] == 0:
        return None
    c0 = Fr(f[0])
    (key, a), = f[1].items()
    a = Fr(a)
    if a <= 0:
        return None
    fa = pfactor(a)
    if fa is None:
        return None
    vec = pf_combine({p_: Fr(x) for p_, x in fa.items()}, {p_: Fr(x) for p_, x in key}, 1)
    if not vec:
        return None
    P = pf_value(vec)
    if P is None or base.q is None or base.q != c0 + P:
        return None                                      # the form must give the exact value of the base, checked exactly
    D = lnq(abs(c0)) - lnpf(vec)
    if not (D.b < NEAR):
        return None
    op_ = one_pow(1 if c0 > 0 else -1, D, t)
    if op_ is None:
        return None
    return make_apf({p_: Fr(x) * t for p_, x in vec.items()}, op_[0], op_[1])


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
                        # (14/11) aq^k too large to expand: aq by its vector (an anchor at an exact power, apf)
                        vec = pfactor(aq)
                        if vec is None:
                            return None
                    else:
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
            r_ = _root_of_offset(base, t)                # (17/11)
            if r_ is not None:
                return r_
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


def _add(a, b):
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


def _sub(a, b):
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
        eq_ = e.q
        if eq_ is None and t.anc is not None and (e.pf is not None or e.off is not None):
            eq_ = exact_of(e)                            # (17/11) an exact power known by its vector (2^4782969, split 92)
        if eq_ is not None and t.anc is not None:
            aq, sg, D = t.anc
            if eq_ == aq:
                return from_lnu(lnq(aq) + D)                  # aq e^D
            d = eq_ - aq
            q2 = abs(d)
            sg2 = -sg if d > 0 else sg
            try:
                return make_anc(q2, sg2, D + lnq(aq) - lnq(q2))
            except Fail:
                if sg2 > 0:
                    raise
                # (17/11) |e - aq (1 + sg e^D)| = |q2 - aq e^D|, two tiny exact pieces that almost cancel: bounds from
                # their logs (ln q2 exact up to rounding, ln aq + D), much tighter than those of e and t themselves
                return bound_diff(exact(q2), from_lnu(lnq(aq) + D))
        if e.q is not None:
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


def bound_diff(a, b):
    """(17/11) |a - b| when no rule gives it (its logs overlap: a and b equal or almost equal): only bounds, from the logs.
    For A >= B > 0, d = ln A - ln B >= 0:  A - B = A (1 - e^-d), and d e^-d <= 1 - e^-d <= d, so
        e^(min ln) |d|min e^(-|d|max) <= |a - b| <= e^(max ln) |d|max,
    with min ln / max ln the lower / upper end of the enclosures of ln a and ln b and |d| in the enclosure of |ln a - ln b|
    (the lower bound only when that enclosure excludes 0; otherwise 0).  The result can only be compared with a window
    (outside); as an operand it is refused (far)."""
    chk(a, b)
    la, lb = lam_iv(a), lam_iv(b)
    d = la - lb
    w_hi = max(abs(d.a), abs(d.b))
    lmax = max(la.b, lb.b)
    if w_hi == 0:
        hi = -mpmath.inf                                  # equal logs: a = b
    else:
        hi = (I(lmax) + iv.log(I(w_hi))).b
    lo = None
    if d.a > 0 or d.b < 0:
        w_lo = min(abs(d.a), abs(d.b))
        lo = (I(min(la.a, lb.a)) + iv.log(I(w_lo)) - I(w_hi)).a
    v = Val(far="bound")
    v.bnd = (lo, hi)
    return v


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


def _mul(a, b):
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


def _div(a, b):
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


# ---------------------------------------------------------------------------------------------------------------------
# (13/11) exact linear forms.  A value may carry  lin = (c0, terms):  u = c0 + sum_k a_k P_k  exactly, with c0 and the
# a_k nonzero rationals and each P_k = prod p^(f_p) an exact power of primes whose exponents f_p are in (0, 1) (the
# integer parts of the exponents are moved into the coefficient: 2^1.8 = 2 * 2^0.8, so equal powers written
# differently have the same key).  Sums, differences, products and quotients by a monomial combine the forms term by
# term; when every term cancels, the value is the exact rational c0.  Only that direction is used (a form that does
# not cancel says nothing more); the n = 9 residual hits of split 14 (pieces VN.2 and VN.5) are exactly such sums,
# e.g. (4^-.9) + (5 - (8^-.6)) = 2^-1.8 + 5 - 2^-1.8 = 5.
LIN_TERMS = 16            # at most this many terms
LIN_POW_MAX = 8           # (17/11) integer powers of a form expanded up to this exponent
LIN_BITS = 2048           # coefficients of the forms up to this size
LIN_C0_BITS = 600_000     # (14/11) the constant c0 up to this size (3^32768 of split 117 has 52 000 bits)


def _lin_ok(c, bits=None):
    return max(c.numerator.bit_length(), c.denominator.bit_length()) <= (LIN_BITS if bits is None else bits)


def lin_mono(pf):
    """(coefficient, key) of the exact power prod p^(e_p): exponents split into integer and fractional parts"""
    coef, key = Fr(1), []
    for p_, e in pf.items():
        e = Fr(e)
        fl = e.numerator // e.denominator
        fr = e - fl
        if fl:
            # (13/11) the size test on the bit length first: fl can be a huge integer (an exact power such as 3^(2^100000)),
            # and a huge int times a float overflows
            if abs(fl).bit_length() > 32 or abs(fl) * math.log2(p_) > LIN_BITS:
                # (17/11) a huge integer power of p (3^(2 10^6) of split 55): kept whole in the key, never expanded, so
                # that forms built from it cancel exactly (9^(10^6))^5 = 3^(10^7); the decision depends only on
                # this exponent, so a given power always gets the same key
                if abs(e.numerator).bit_length() > 64:
                    return None
                key.append((p_, e))
                continue
            coef *= Fr(p_) ** fl
        if fr:
            key.append((p_, fr))
    return coef, tuple(sorted(key))


def lin_of(v):
    """the exact linear form of v (c0, {key: a}), or None"""
    if v.q is not None and _lin_ok(v.q, LIN_C0_BITS):
        return (v.q, {})
    # (17/11) an exact rational too large for the constant of a form: its form or its exact power, when it has one
    if v.lin is not None:
        return v.lin[0], dict(v.lin[1])
    if v.pf is not None:
        m = lin_mono(v.pf)
        if m is None:
            return None
        coef, key = m
        if not key:
            return (coef, {}) if _lin_ok(coef) else None
        return Fr(0), {key: coef}
    return None


def _lin_norm(c0, terms):
    terms = {k: a for k, a in terms.items() if a != 0}
    if len(terms) > LIN_TERMS or not _lin_ok(c0, LIN_C0_BITS) or not all(_lin_ok(a) for a in terms.values()):
        return None
    return c0, terms


def _lin_value(f):
    """(14/11) the value c0 + sum a_k P_k of a form, in interval arithmetic (c0 is exact, so the cancellations
    between huge rationals have already happened exactly); (value, sign) with the magnitude as a value of lvl, or
    raises Fail when the enclosure contains 0"""
    c0, terms = f
    x = ivq(c0)
    for key, a in terms.items():
        lp = I(0)
        for (p_, fr) in key:
            lp = lp + iv.mpf(fr.numerator) / fr.denominator * iv.log(iv.mpf(p_))
        x = x + ivq(a) * iv.exp(lp)
    if x.a > 0:
        return from_x(x), 1
    if x.b < 0:
        return from_x(-x), -1
    raise Fail("linear form: enclosure contains 0")


def _lin_lin(la, lb, sg):
    """la + sg lb"""
    t = dict(la[1])
    for k, a in lb[1].items():
        t[k] = t.get(k, Fr(0)) + sg * a
    return _lin_norm(la[0] + sg * lb[0], t)


def _lin_mul(la, lb):
    c0 = la[0] * lb[0]
    t = {}
    A = [((), la[0])] + list(la[1].items())
    B = [((), lb[0])] + list(lb[1].items())
    if len(A) * len(B) > 4 * LIN_TERMS:
        return None
    for ka, a in A:
        for kb, b in B:
            if not ka and not kb:
                continue
            m = lin_mono(pf_combine(dict(ka), dict(kb), 1))
            if m is None:
                return None
            coef, key = m
            if key:
                t[key] = t.get(key, Fr(0)) + a * b * coef
            else:
                c0 += a * b * coef
    return _lin_norm(c0, t)


def _lin_inv(lb):
    """the form of 1/b when b is a rational or a monomial, else None"""
    c0, t = lb
    if not t:
        return (1 / c0, {}) if c0 != 0 else None
    if c0 != 0 or len(t) != 1:
        return None
    (k, a), = t.items()
    m = lin_mono({p_: -x for p_, x in k})
    if m is None:
        return None
    coef, key = m
    return Fr(0), {key: coef / a}


def _gt(a, b):
    """True when a > b certainly, False when b > a certainly, None when lvl's enclosures do not tell"""
    if a.q is not None and b.q is not None:
        return None if a.q == b.q else a.q > b.q
    if a.x is not None and b.x is not None:
        if a.x.a > b.x.b:
            return True
        if b.x.a > a.x.b:
            return False
    try:
        la, lb = lam_iv(a), lam_iv(b)
    except Fail:
        return None
    if la.a > lb.b:
        return True
    if lb.a > la.b:
        return False
    return None


def _lin_exact(c0, encl):
    """the exact rational c0 of a cancelled form, checked against an enclosure (an interval of the value, or None)"""
    if c0 < 0:
        raise Fail("linear form: negative magnitude")                 # cannot happen for a correct form
    if encl is not None and not (encl.a <= ivq(c0).b and ivq(c0).a <= encl.b):
        raise Fail("linear form: exact value outside the enclosure")  # cannot happen for a correct form
    return exact(c0)


def _form_agrees(r, f):
    """(17/11, parche 7) a safety check of every form before it is attached: its value (interval arithmetic) must be
    positive and its ln must meet the enclosure of ln r; a form that the enclosures cannot check (its value encloses 0,
    or no ln of r) is accepted only if r has no enclosure either.  A form that fails is a defect of a rule: Fail, so the
    value stays undecided (never a form that contradicts the value)"""
    try:
        lr = lam_iv(r)
    except Fail:
        lr = None
    try:
        v, sg = _lin_value(f)
    except Fail:
        return True                                      # nothing to compare: the form keeps the rule that made it
    if sg < 0:
        return False
    if lr is None:
        return True
    try:
        lf = lam_iv(v)
    except Fail:
        return True
    return lf.a <= lr.b and lr.a <= lf.b


def _attach(r, f):
    if f is not None and (r.q is None or not _lin_ok(r.q, LIN_C0_BITS)) and r.far is None and f[1]:
        if not _form_agrees(r, f):
            raise Fail("linear form disagrees with the value")
        r.lin = (f[0], tuple(sorted(f[1].items())))
        if r.pf is None and f[0] == 0 and len(f[1]) == 1:
            # (14/11) a form that reduced to one term a P (a > 0 rational, P an exact power) is the exact power a P:
            # its vector, so that integer powers of it are exact ((4^-.9 + 3/8^.6)^5 = (4 2^-1.8)^5 = 2)
            (k, a), = f[1].items()
            fa = pfactor(a) if a > 0 else None
            if fa is not None:
                r.pf = pf_combine({p_: Fr(x) for p_, x in fa.items()}, dict(k), 1)
    return r


def _lin_ratio(la, lb):
    """(14/11) the rational r with la = r lb (the same powers, every coefficient and the constant in the ratio r), or
    None: then a / b = r exactly ((6 - .5^.8) / (3 - 4^-.9) = 2)"""
    if not la[1] or set(la[1]) != set(lb[1]):
        return None
    r = None
    for k, x in la[1].items():
        y = lb[1][k]
        if y == 0:
            return None
        q = Fr(x) / Fr(y)
        if r is None:
            r = q
        elif q != r:
            return None
    if la[0] != r * lb[0] or r <= 0:
        return None
    return r


def add(a, b):
    """a + b (magnitudes); (13/11) with the exact linear forms"""
    la, lb = lin_of(a), lin_of(b)
    f = None
    if la is not None and lb is not None and (la[1] or lb[1]):
        f = _lin_lin(la, lb, 1)
        if f is not None and not f[1]:
            return _lin_exact(f[0], a.x + b.x if a.x is not None and b.x is not None else None)
    try:
        return _attach(_add(a, b), f)
    except Fail:
        if f is None:
            raise
        v, sg = _lin_value(f)                            # (14/11) the exact form when the enclosures fail
        return _attach(v, f if sg > 0 else None)


def sub(a, b):
    """|a - b|; (13/11) with the exact linear forms (the sign of a - b from the enclosures when the form does not
    cancel)"""
    la, lb = lin_of(a), lin_of(b)
    f = None
    if la is not None and lb is not None and (la[1] or lb[1]):
        f = _lin_lin(la, lb, -1)
        if f is not None and not f[1]:
            d = f[0]
            encl = None
            if a.x is not None and b.x is not None:
                dx = a.x - b.x
                encl = I(0, max(abs(dx.a), abs(dx.b))) if dx.a <= 0 <= dx.b else (dx if dx.a > 0 else -dx)
            return _lin_exact(abs(d), encl)
    try:
        r = _sub(a, b)
    except Fail:
        if f is None:
            return bound_diff(a, b)                      # (17/11) only bounds of |a - b| (raises Fail without logs)
        try:
            v, sg = _lin_value(f)                        # (14/11) |a - b| from the exact form of a - b
        except Fail:
            return bound_diff(a, b)                      # (17/11) (split 115)
        return _attach(v, f if sg > 0 else (-f[0], {k: -x for k, x in f[1].items()}))
    if r.far == "bound" and f is not None:
        # (17/11, parche 7) only bounds from the rules (bound_diff), but an exact form of a - b: its value is better
        # (|0.1 + 9^16384 - (6^.2 + 3^32768)| = 6^.2 - 0.1, split 173); the bounds only if the form cannot be evaluated
        try:
            v, sg = _lin_value(f)
        except Fail:
            return r
        return _attach(v, f if sg > 0 else (-f[0], {k: -x for k, x in f[1].items()}))
    if f is not None and (r.q is None or not _lin_ok(r.q, LIN_C0_BITS)):
        # the form of |a - b| is that of a - b or of b - a: by the order of a and b.  (17/11, parche 7) also when r is a
        # huge exact rational, whose form is attached too (since the parche 6): before, the form of a - b was attached
        # as is, so |c - P| (c < P, P a huge exact power) got the form c - P, the opposite of its value
        g = _gt(a, b)
        if g is None:
            f = None
        elif not g:
            f = (-f[0], {k: -x for k, x in f[1].items()})
    return _attach(r, f)


def mul(a, b):
    """a * b; (13/11) with the exact linear forms"""
    la, lb = lin_of(a), lin_of(b)
    f = None
    if la is not None and lb is not None and (la[1] or lb[1]):
        f = _lin_mul(la, lb)
        if f is not None and not f[1]:
            return _lin_exact(f[0], a.x * b.x if a.x is not None and b.x is not None else None)
    return _attach(_mul(a, b), f)


def _off_ratio(a, b):
    """(14/11) a / b when one of them is an exact power X (vector pf, not expanded) and the other is X + sg c exactly
    (offset form with the same X, c tiny relative to X): with t = c / X,
        X / (X + sg c) = 1 / (1 + sg t) = 1 - sg t th,  th = 1 / (1 + sg t) in [1 - 2t, 1 + 2t],
        (X + sg c) / X = 1 + sg t  exactly;
    as anchors at 1 (make_anc) with ln t = ln c - ln X.  The split 191 of n = 9 has 3^(2^49) / (c + 9^(2^48)), where the
    two logarithms (about 6e14) agree to 10^(-3e14) and no precision separates them.  None otherwise."""
    for x, y, xnum in ((a, b, True), (b, a, False)):
        if x.pf is None or x.q is not None or y.off is None or y.q is not None:
            continue
        k, sg, c = y.off
        if k != pf_key(x.pf):
            continue
        if c.q is not None and c.q == 0:
            return exact(1)
        try:
            D = lam_iv(c) - lnpf(x.pf)
        except Fail:
            return None
        if D.b >= NEAR:
            return None
        if xnum:
            return make_anc(Fr(1), -sg, D + I(-TINYW, TINYW))
        return make_anc(Fr(1), sg, D)
    return None


def div(a, b):
    """a / b; (13/11) with the exact linear forms when b is a rational or a monomial; (14/11) X / (X +- c)"""
    r_ = _off_ratio(a, b)
    if r_ is not None:
        return r_
    la, lb = lin_of(a), lin_of(b)
    f = None
    if la is not None and lb is not None and la[1] and lb[1]:
        rr = _lin_ratio(la, lb)
        if rr is not None:
            return _lin_exact(rr, a.x / b.x if a.x is not None and b.x is not None and b.x.a > 0 else None)
    if la is not None and lb is not None and (la[1] or lb[1]) and b.q != 0:
        ib = _lin_inv(lb)
        if ib is not None:
            f = _lin_mul(la, ib)
            if f is not None and not f[1]:
                return _lin_exact(f[0], a.x / b.x if a.x is not None and b.x is not None and b.x.a > 0 else None)
    return _attach(_div(a, b), f)


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


TW_EXACT_BITS = 256       # (13/11) a rational coefficient C_p of a canonical form is expanded up to this size


def tower_exact(tw):
    """(13/11) the exact value of a canonical form whose coefficients C_p are all rational (integer exponents in their
    keys): ln u = sum_p sg_p C_p ln p with rational coefficients, so u = prod p^(sg_p C_p), an exact power of primes
    (a rational when the exponents are integers).  None when some C_p is irrational or too large."""
    vec = {}
    for (p_, sg, k_) in tw:
        if any(Fr(x).denominator != 1 for (_, x) in k_):
            return None
        if any(abs(int(x)).bit_length() > 32 for (_, x) in k_) or \
                sum(abs(int(x)) * math.log2(r_) for (r_, x) in k_) > TW_EXACT_BITS:
            return None                                  # (13/11) bit length first: no float overflow
        c = Fr(1)
        for (r_, x) in k_:
            c *= Fr(r_) ** int(x)
        vec[p_] = vec.get(p_, Fr(0)) + sg * c
    vec = {p_: x for p_, x in vec.items() if x != 0}
    if not vec:
        return exact(1)
    q = pf_value(vec)
    if q is not None:
        return exact(q)
    try:
        lu = I(0)
        for p_, x in vec.items():
            lu = lu + ivq(Fr(x)) * lnq(Fr(p_))
        return pf_result(vec, from_lnu(lu))
    except Fail:
        return None


def pw(base, e, sgn):
    """base ^ (sgn * e), base and e magnitudes.  (11/11) The result of an exact power of primes raised to an exact
    power of primes that is not a usable rational (2^(35/2), 3^8.4) also carries its canonical form (tower); and a
    value with such a form raised to an exact rational k keeps one (its ln is multiplied by sgn k)."""
    if base.tw is not None and e.q is None and e.pf is not None:
        # (13/11) a value with a canonical form (ln base = sum_p sg_p C_p ln p) raised to an exact power of primes
        # F = prod r^(x_r) > 0 that is not rational: ln of the result = sum_p (sg_p sgn) (C_p F) ln p, with each C_p F
        # again an exact power of primes.  When every C_p F is rational the result is the exact power prod p^(sg_p
        # sgn C_p F): e.g. (5^-(8^.6))^(4^-.9) = 5^-(2^(9/5) 2^(-9/5)) = 1/5 (split 10 of n = 9, piece VN.7)
        try:
            F = {p_: Fr(x) for p_, x in e.pf.items()}
            tw = tuple(sorted((p_, sg * sgn, pf_key(pf_combine(dict(k_), F, 1))) for (p_, sg, k_) in base.tw))
        except Exception:
            tw = None
        if tw is not None:
            v = tower_exact(tw)
            if v is not None:
                return v
            r = _pw(base, e, sgn)
            if r.q is None and r.far is None and r.pf is None:
                r.tw = tw
            return r
    r = _pw(base, e, sgn)
    if sgn > 0 and e.q is not None and e.q == 1 and r.q is None and r.far is None:
        # (17/11) x^1 = x exactly: the exact forms of x (an offset X + c, a linear form, a prime vector, a tower) are
        # those of the result too; the evaluation above keeps its enclosure but may lose them
        # ((3^4194304 + 1/5)^1 - 9^2097152 = 1/5, split 128)
        for k_ in ("off", "lin", "pf", "tw"):
            if getattr(r, k_) is None and getattr(base, k_) is not None:
                setattr(r, k_, getattr(base, k_))
    if r.q is not None and r.pf is None and not _lin_ok(r.q, LIN_C0_BITS) and base.q is not None and base.q > 0 \
            and _lin_ok(base.q) and e.q is not None and e.q.denominator == 1:
        # (17/11) a huge exact rational b^k (b a small rational, k an integer): also its prime vector, so that it can be a
        # term of a form (9^(10^6) = 3^(2 10^6), split 55)
        fb = pfactor(base.q)
        if fb is not None:
            r.pf = {p_: Fr(x) * sgn * e.q for p_, x in fb.items()}
    if sgn > 0 and e.q is not None and e.q.denominator == 1 and 2 <= e.q <= LIN_POW_MAX and r.far is None \
            and r.pf is None and (r.q is None or not _lin_ok(r.q, LIN_C0_BITS)):
        # (17/11) an integer power k of a sum with an exact form: the form expanded (binomial products of the terms,
        # each an exact power): (9^(10^6) + c)^5 = sum C(5, j) c^(5 - j) 3^(2 10^6 j), whose difference with 3^(10^7)
        # then cancels exactly (split 55)
        la = lin_of(base)
        if la is not None and la[1]:
            f = la
            for _ in range(int(e.q) - 1):
                f = _lin_mul(f, la) if f is not None else None
            if f is not None:
                _attach(r, f)
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
        a_ = anchor_of_q(base.q)
        if a_ is not None:                               # (14/11) q = M (1 + d), d tiny: the anchored rules
            r_ = pow_anchored(a_, e, sgn)
            if r_ is not None:
                return r_
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
        # (14/11) an exact rational when the exponents become integers ((2^(1/5))^5 = 2)
        return pf_result({p_: a * sgn * k for p_, a in base.pf.items()}, v)
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
    if v.far == "bound":                                   # (17/11) bounds of a difference (bound_diff)
        llo, lhi = iv.log(lo), iv.log(hi)
        return v.bnd[1] < llo.a or (v.bnd[0] is not None and v.bnd[0] > lhi.b)
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
