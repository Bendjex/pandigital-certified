"""
A3b - rigorous error bounds for the e records.

For each record formula (1 + L)^R this script
  1. checks the digit usage (each digit once; for in-order records also the order),
  2. checks EXACTLY that the formula has the shape (1 + L)^R with L = b_L^e_L, R = b_R^e_R,
     b rational, e an exact (possibly huge) integer or rational, and L * R = 1
     (via prime-exponent vectors, no floating point),
  3. computes with Arb ball arithmetic (python-flint, rigorous) the bounds of
         Lemma.  For real N >= 1:   e/(2N) * (1 - 11/(12N))  <  e - (1+1/N)^N  <  e/(2N)
     and the number k of correct decimals, defined as the largest k with |x - e| < 0.5 * 10^-k.

Proof of the lemma (t = 1/N in (0,1]):
  f = N ln(1+t) = 1 - t/2 + t^2/3 - t^3/4 + ...  (alternating, terms strictly decreasing)
  so 1 - t/2 < f < 1 - t/2 + t^2/3, i.e. delta := 1 - f lies in (t/2 - t^2/3, t/2).
  e - (1+1/N)^N = e (1 - e^-delta) and delta - delta^2/2 < 1 - e^-delta < delta for delta > 0.
  Upper: < e*delta < e t/2.   Lower: > e (delta - delta^2/2) >= e (t/2 - t^2/3 - t^2/8) = e t/2 (1 - 11t/12)
  (g(d) = d - d^2/2 is increasing on (0,1)).
"""
import re, sys
from fractions import Fraction as F
from flint import arb, ctx

RECORDS = [
    # (label, digits that must be used, in_order, formula)
    ("1-7 any order (ours)",      "1234567",    False, "(1+4^(-36))^(.5^(-72))"),
    ("1-8 any order (ours)",      "12345678",   False, "(1+2^(-8^7))^(4^(.5^(-6/.3)))"),
    ("1-9 any order (Bamberger)", "123456789",  False, "(1+.2^(9^(6*7)))^(5^(3^84))"),
    ("1-8 in order (ours)",       "12345678",   True,  "(1+.2^34)^(5^(6*7-8))"),
    ("1-9 in order (ours)",       "123456789",  True,  "(1+.2^(3^4))^(5^(6*(7+8)*.9))"),
    ("0-9 any order (ours)",      "0123456789", False, "(1+.2^(9^70))^(5^(3^(84/.6)))"),
]

# ---------------- tiny exact parser ----------------
TOK = re.compile(r"\s*(\d*\.\d+|\d+|[-+*/^()])")
def tokenize(s):
    pos, out = 0, []
    while pos < len(s):
        m = TOK.match(s, pos)
        if not m: raise ValueError("bad token at " + s[pos:])
        out.append(m.group(1)); pos = m.end()
    return out

class P:
    """recursive-descent parser -> AST tuples: ('num', Fraction, text) | ('neg', a) | (op, a, b)"""
    def __init__(self, s): self.t = tokenize(s); self.i = 0
    def peek(self): return self.t[self.i] if self.i < len(self.t) else None
    def eat(self, x=None):
        tok = self.t[self.i]
        if x is not None and tok != x: raise ValueError(f"expected {x}, got {tok}")
        self.i += 1; return tok
    def expr(self):
        a = self.term()
        while self.peek() in ("+", "-"):
            op = self.eat(); a = (op, a, self.term())
        return a
    def term(self):
        a = self.unary()
        while self.peek() in ("*", "/"):
            op = self.eat(); a = (op, a, self.unary())
        return a
    def unary(self):
        if self.peek() == "-": self.eat(); return ("neg", self.unary())
        return self.power()
    def power(self):
        a = self.atom()
        if self.peek() == "^": self.eat(); return ("^", a, self.unary())   # right-assoc, unary minus allowed in exponent
        return a
    def atom(self):
        tok = self.peek()
        if tok == "(":
            self.eat("("); a = self.expr(); self.eat(")"); return a
        self.eat()
        if "." in tok:
            ip, fp = tok.split(".")
            return ("num", F(int(ip or "0") * 10 ** len(fp) + int(fp), 10 ** len(fp)), tok)
        return ("num", F(int(tok)), tok)

def ev(a):
    """exact value (Fraction) of a subtree whose value is rational and not astronomically large"""
    k = a[0]
    if k == "num": return a[1]
    if k == "neg": return -ev(a[1])
    x, y = ev(a[1]), ev(a[2])
    if k == "+": return x + y
    if k == "-": return x - y
    if k == "*": return x * y
    if k == "/": return x / y
    if k == "^":
        if y.denominator != 1: raise ValueError("non-integer exponent in exact subtree")
        if abs(y) > 10 ** 6: raise ValueError("exponent too large for explicit evaluation")
        return x ** int(y)
    raise ValueError(k)

def factor(n):
    n = abs(n); out = {}; p = 2
    while p * p <= n:
        while n % p == 0: out[p] = out.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: out[n] = out.get(n, 0) + 1
    return out

def vec(q):
    """prime-exponent vector of a positive rational"""
    v = dict(factor(q.numerator))
    for p, e in factor(q.denominator).items(): v[p] = v.get(p, 0) - e
    return {p: e for p, e in v.items() if e}

def power_parts(a):
    """a = ('^', base, exponent) possibly under unary minus -> (sign, base Fraction, exponent Fraction).
       The exponent may be a huge integer (computed exactly), the base must be a small rational."""
    sign = 1
    while a[0] == "neg": sign, a = -sign, a[1]
    if a[0] != "^": raise ValueError("not a power")
    base = ev(a[1])
    ex = a[2]; es = 1
    while ex[0] == "neg": es, ex = -es, ex[1]
    if ex[0] == "^":                       # exponent is itself a power b^c, evaluated exactly as a big integer
        b, c = ev(ex[1]), ev(ex[2])
        if c.denominator != 1: raise ValueError("bad exponent")
        e = F(b) ** int(c)
    else:
        e = ev(ex)
    return sign, base, es * e

def check(label, digits, in_order, s):
    used = re.findall(r"\d", s)
    assert sorted(used) == sorted(digits), f"{label}: digits used {''.join(used)}"
    if in_order: assert "".join(used) == digits, f"{label}: digits not in order"
    tree = P(s).expr()
    assert tree[0] == "^" and tree[1][0] == "+", "shape must be (U + L)^R"
    U, Lt, Rt = tree[1][1], tree[1][2], tree[2]
    assert ev(U) == 1, "U must be exactly 1"
    sL, bL, eL = power_parts(Lt)
    sR, bR, eR = power_parts(Rt)
    assert bL > 0 and bR > 0 and eL.denominator == 1 and eR.denominator == 1
    # L * R = 1 exactly  <=>  eL * vec(bL) + eR * vec(bR) = 0 (and signs)
    vL, vR = vec(bL), vec(bR)
    for p in set(vL) | set(vR):
        assert vL.get(p, 0) * eL + vR.get(p, 0) * eR == 0, f"{label}: L*R != 1"
    assert sL * sR == 1
    # N = R = bR^eR > 0
    ctx.prec = 4 * max(64, int(eR).bit_length() + 200) if eR.denominator == 1 else 2000
    lnN = arb(int(eR)) * arb(bR.numerator).log() - arb(int(eR)) * arb(bR.denominator).log()
    ln10 = arb(10).log()
    log10N = lnN / ln10
    E = arb(1).exp()
    # log10 of the bounds:  log10(e/2) - log10 N  (+ log10(1 - 11/(12N)) for the lower one, |.| < 1/N)
    up = (E / 2).log() / ln10 - log10N
    lo_corr = arb(11) / 12 * (-lnN).exp() * 2        # |log10(1 - z)| <= 2z/ln10 < 2z for z <= 1/2
    lo = up - lo_corr
    # k = largest k with err < 0.5*10^-k.  err < e/(2N) <= 0.5*10^-k  iff  k <= log10(N/e) =: u
    u = log10N - E.log() / ln10
    kf = u.floor().unique_fmpz()          # None if the ball for u straddles an integer
    assert kf is not None, "precision too low to decide floor(u)"
    # k = kf is valid (kf <= u) and kf+1 is not: kf+1 > u + delta, delta = -log10(1-11/(12N)) < 2*(11/12)/N
    assert bool(u >= arb(kf)) and bool(u + lo_corr < arb(kf) + 1)
    k_floor, ok = kf, True
    return dict(label=label, formula=s, N=f"({bR})^({eR if int(eR).bit_length() < 64 else str(int(eR).bit_length()) + '-bit integer'})",
                log10N=log10N, minus_log10_err=(-up, -lo), u=u, k=k_floor, kexact=ok, eR=eR, bR=bR)

if __name__ == "__main__":
    for lab, dg, io, s in RECORDS:
        r = check(lab, dg, io, s)
        print(f"== {r['label']}:  {r['formula']}")
        print(f"   identity L*R = 1 verified exactly; N = {r['N']}")
        print(f"   log10 N           = {r['log10N'].str(25, radius=False)}")
        a, b = r['minus_log10_err']
        print(f"   -log10(error)     in [{a.str(25, radius=False)}, {b.str(25, radius=False)}]")
        print(f"   correct decimals k = floor(log10(N/e)) = {r['k']}  (exactness of floor certified: {r['kexact']})")
