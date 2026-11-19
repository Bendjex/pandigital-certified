"""
B1b / B2b - exhaustive, rigorous search of the identity family F for e.

    F :  (U + s*A^X)^(t*C^Y)      s, t in {+1, -1}
         U, A, X, C, Y  expressions (rules below) on pairwise disjoint, non-empty digit blocks
         that together use every digit exactly once;  U = 1 exactly;  s*A^X * t*C^Y = 1 exactly.
    Then the value is (1 + 1/N)^N with N = t*C^Y, and the error is ~ e/(2|N|)  (A3b).
    Every two-level formula (U + A^(B^X))^(C^(D^Y)) is in F (take X' = B^X, Y' = D^Y), so F
    contains Sabey, Bamberger, Fastfaxr, McKernon's family and our records.

Expression rules (any-order mode): a block of k digits is either a literal (its digits in any
order, optionally with ONE decimal point before any digit: 12, 21, 1.2, .12, ...), or op(u, v)
with u, v expressions on a split of the block, op in {+, -, *, /, ^} (both orders), or -u.
x^y is defined iff x > 0, or x = 0 and y > 0, or x < 0 and y an integer.  x/0 is undefined.

Exactness.  Every value of a block of <= 2 digits is 0 or a signed radical s * prod p^(e_p),
e_p rational ("R").  Blocks of 3 digits may also give radical+rational sums or x^y with y an
irrational radical ("T", handled numerically with Arb balls).  R-values are exact.

Search.  For each candidate (U, side L = (A, X), side R = (C, Y)) the identity is
    mu_L + mu_R = 0,   mu_L = X ln|A|,  mu_R = Y ln|C|,     and N = t*C^Y with ln|N| = mu_R.
 * Screening: mu is enclosed in an Arb ball; two entries can form an identity only if
   mu_L + mu_R contains 0  (rigorous: if the balls are disjoint the identity is false).
 * Exact check: if A, C, X, Y are all R-values, the identity holds iff for every prime p
       X * alpha_p + Y * gamma_p = 0      (ln|A| = sum alpha_p ln p, ln|C| = sum gamma_p ln p);
   "if" is trivial; "only if" is Baker's theorem (logs of primes are linearly independent over
   the algebraic numbers).  Overlaps that fail the exact check, or involve T-values, are reported
   as UNRESOLVED (none may remain for the theorem to hold).
 * Size lemma: |mu| <= SB(k) for a side of k digits, computed below from the enumeration.
   Configurations whose smaller side has <= 3 digits are bounded by SB(3); they are skipped only
   if SB(3) < theta (the target), so blocks of more than 3 digits are never needed.
"""
import sys, itertools, math, time
from fractions import Fraction as F
from functools import lru_cache
from flint import arb, ctx

ctx.prec = 256
T0 = time.time()

# ------------------------------------------------------------ exact signed radicals
from sympy import factorint
@lru_cache(None)
def factor(n):
    return tuple(sorted(factorint(n).items())) if n > 1 else ()

ZERO = ("Z",)
FACTOR_LIMIT = 10 ** 18   # larger rationals are not factored: kept as numeric enclosures (kind "Q")

def R_from_frac(q):
    if q == 0: return ZERO
    if max(abs(q.numerator), q.denominator) > FACTOR_LIMIT:
        return ("T", 1 if q > 0 else -1, arb(abs(q.numerator)).log() - arb(q.denominator).log(), "Q")
    v = {}
    for p, e in factor(abs(q.numerator)): v[p] = v.get(p, 0) + e
    for p, e in factor(q.denominator): v[p] = v.get(p, 0) - e
    return ("R", 1 if q > 0 else -1, tuple(sorted((p, F(e)) for p, e in v.items() if e)))

def is_rat(x): return x[0] == "R" and all(e.denominator == 1 for _, e in x[2])
def to_frac(x):
    q = F(1)
    for p, e in x[2]: q *= F(p) ** int(e)
    return x[1] * q

ONE = R_from_frac(F(1))
LNP = {}
def lnp(p):
    if p not in LNP: LNP[p] = arb(p).log()
    return LNP[p]

def lnabs(x):
    """Arb enclosure of ln|x| (x nonzero)"""
    if x[0] == "R":
        s = arb(0)
        for p, e in x[2]: s += lnp(p) * arb(e.numerator) / e.denominator
        return s
    return x[2]

def as_arb(x):
    if x[0] == "Z": return arb(0)
    return x[1] * lnabs(x).exp()

def neg(x):
    if x[0] == "Z": return x
    return (x[0], -x[1]) + x[2:]

def T_from_arb(v, kind):
    if v.contains(0): raise ValueError("UNDECIDABLE ZERO")
    return ("T", 1 if v > 0 else -1, abs(v).log(), kind)

BIG = 20000.0     # |ln|v|| beyond which rationals are not expanded into explicit integers

def huge(x):
    return x[0] != "Z" and abs(float(lnabs(x).mid())) > BIG

def add(x, y):
    if x[0] == "Z": return y
    if y[0] == "Z": return x
    if huge(x) or huge(y):                                  # astronomically large/small: numeric enclosure only
        return add_numeric(x, y)
    if is_rat(x) and is_rat(y): return R_from_frac(to_frac(x) + to_frac(y))
    if x[0] == "R" and y[0] == "R" and x[2] == y[2]:        # same radical r: (+-r) + (+-r)
        return ZERO if x[1] == -y[1] else mul(x, R_from_frac(F(2)))
    # an irrational radical plus a rational (or two different radicals): numeric, provably != 0
    # (radical irrational + rational != 0; two distinct real radicals r1 + r2 = 0 iff r1 = -r2, excluded above)
    return add_numeric(x, y)

def add_numeric(x, y):
    """rigorous enclosure of x + y in log form: ln|big + small| = ln|big| + log1p(+-|small/big|)"""
    lx, ly = lnabs(x), lnabs(y)
    if ly < lx: big, sm, lb, ls = x, y, lx, ly
    elif lx < ly: big, sm, lb, ls = y, x, ly, lx
    else: return T_from_arb(as_arb(x) + as_arb(y), "A")    # comparable magnitudes: direct sum
    r = (ls - lb).exp()
    t = r if big[1] == sm[1] else -r
    if not (t > -0.5): return T_from_arb(as_arb(x) + as_arb(y), "A")
    return ("T", big[1], lb + t.log1p(), "A")

def mul(x, y):
    if x[0] == "Z" or y[0] == "Z": return ZERO
    if x[0] == "R" and y[0] == "R":
        v = dict(x[2])
        for p, e in y[2]: v[p] = v.get(p, F(0)) + e
        return ("R", x[1] * y[1], tuple(sorted((p, e) for p, e in v.items() if e)))
    return ("T", x[1] * y[1], lnabs(x) + lnabs(y), "A")

def inv(x):
    if x[0] == "Z": return None
    if x[0] == "R": return ("R", x[1], tuple((p, -e) for p, e in x[2]))
    return ("T", x[1], -x[2], x[3])

def div(x, y):
    iy = inv(y)
    return None if iy is None else mul(x, iy)

def pw(x, y):
    if y[0] == "Z": return None if x[0] == "Z" else ONE
    if x[0] == "Z": return ZERO if y[1] > 0 else None
    if is_rat(y) and huge(y):                               # exponent too large to expand exactly
        integer = all(e > 0 for _, e in y[2])
        if x[1] < 0 and not integer: return None
        even = integer and any(p == 2 for p, e in y[2])
        sg = 1 if x[1] > 0 or even else -1
        if x[0] == "R" and not x[2]: return ("R", sg, ())       # (+-1)^y
        return ("T", sg, as_arb(y) * lnabs(x), "H")
    if is_rat(y):
        q = to_frac(y)
        if x[1] < 0 and q.denominator != 1: return None
        sg = 1 if x[1] > 0 or q.numerator % 2 == 0 else -1
        if x[0] == "R": return ("R", sg, tuple((p, e * q) for p, e in x[2]))
        return ("T", sg, x[2] * (arb(q.numerator) / q.denominator), x[3])
    # irrational exponent: base must be positive
    if x[1] < 0: return None
    if x[0] == "R" and not x[2]: return ONE                        # 1^y
    return ("T", 1, as_arb(y) * lnabs(x), "T")

# ------------------------------------------------------------ block value sets
@lru_cache(None)
def V(S):
    """dict value -> one expression string, for the sorted digit tuple S"""
    out = {}
    def put(v, s):
        if v is None: return
        k = v if v[0] != "T" else ("T", s)
        if k not in out: out[k] = (v, s)
    for perm in set(itertools.permutations(S)):
        t = "".join(perm)
        put(R_from_frac(F(int(t))), t)
        for i in range(len(t)):
            put(R_from_frac(F(int(t), 10 ** (len(t) - i))), t[:i] + "." + t[i:])
    n = len(S)
    for r in range(1, n // 2 + 1):
        for S1 in itertools.combinations(S, r):
            S2 = tuple(d for d in S if d not in S1)
            if r * 2 == n and S1 > S2: continue
            for (u, su) in V(S1).values():
                for (v, sv) in V(S2).values():
                    for f, s in ((add(u, v), f"({su}+{sv})"), (add(u, neg(v)), f"({su}-{sv})"),
                                 (add(v, neg(u)), f"({sv}-{su})"), (mul(u, v), f"({su}*{sv})"),
                                 (div(u, v), f"({su}/{sv})"), (div(v, u), f"({sv}/{su})"),
                                 (pw(u, v), f"({su}^{sv})"), (pw(v, u), f"({sv}^{su})")):
                        put(f, s)
    for k, (v, s) in list(out.items()):
        put(neg(v), "-" + s)
    return out

# ------------------------------------------------------------ sides
@lru_cache(None)
def Vinfo(S):
    """per value: (value, expr, ln|v| as Arb, float ln|v|, is-integer flag); zero excluded"""
    out = []
    for (v, s) in V(S).values():
        if v[0] == "Z": continue
        L = lnabs(v)
        isint = v[0] == "R" and all(e.denominator == 1 and e > 0 for _, e in v[2])   # integer (no huge evaluation)
        out.append((v, s, L, float(L.mid()), isint))
    out.sort(key=lambda t: -t[3])
    return out

def defined_pow(a, xinfo):
    return a[1] > 0 or xinfo[4]           # a < 0 needs an integer exponent

def side_entries(T, theta):
    """all (ln|mu|, sign(mu), a, x, exprs) with A on SA, X on SX, SA u SX = T and ln|mu| >= theta possible"""
    res = []
    T = tuple(sorted(T))
    for r in range(1, len(T)):
        for SA in itertools.combinations(T, r):
            SX = tuple(d for d in T if d not in SA)
            VX = Vinfo(SX)
            for (a, sa, La, fa, _) in Vinfo(SA):
                if a[0] == "R" and not a[2]: continue                  # |A| = 1
                if La.contains(0): raise ValueError("cannot decide sign of ln|A|")
                lla = abs(La).log(); flla = float(lla.mid())
                for xi in VX:                                          # sorted by ln|x| decreasing
                    (x, sx, Lx, fx, _) = xi
                    if fx + flla < theta - 1e-3: break                 # float prefilter (margin >> rounding)
                    if not defined_pow(a, xi): continue
                    lmu = Lx + lla
                    if lmu.upper() < theta: continue
                    sgn = x[1] * (1 if La > 0 else -1)
                    res.append((lmu, sgn, a, x, sa, sx))
    return res

def side_bound(T):
    """rigorous upper bound of ln|mu| over all (A, X) on the digit set T"""
    best = -1e300
    T = tuple(sorted(T))
    for r in range(1, len(T)):
        for SA in itertools.combinations(T, r):
            SX = tuple(d for d in T if d not in SA)
            VX = Vinfo(SX)
            if not VX: continue                                        # X can only be 0: A^0 = 1, no side
            mx = VX[0][2]                                              # largest ln|x| (list sorted)
            for (a, _, La, _, _) in Vinfo(SA):
                if a[0] == "R" and not a[2]: continue
                if La.contains(0): raise ValueError("cannot decide sign of ln|A|")
                b = float((mx + abs(La).log()).upper()) * (1 + 1e-12) + 1e-12
                best = max(best, b)
    return best

def exact_identity(a, x, c, y):
    """True iff x ln|a| + y ln|c| = 0, for R-values (Baker); None if not decidable here"""
    if any(v[0] != "R" for v in (a, x, c, y)): return None
    al, ga = dict(a[2]), dict(c[2])
    if set(al) != set(ga): return False
    rho = div(x, y)                                    # need x/y = -gamma_p/alpha_p for all p
    targets = {-ga[p] / al[p] for p in al}
    if len(targets) != 1: return False
    return rho == R_from_frac(targets.pop())

def search(DIG, theta, units_max=3):
    DIG = tuple(DIG)
    hits, unresolved, skipped = [], [], []
    SB = {}
    for k in (2, 3):
        SB[k] = max(side_bound(T) for T in itertools.combinations(DIG, k))
    for ku in range(1, units_max + 1):
        for U in itertools.combinations(DIG, ku):
            if not any(v == ONE for (v, _) in V(U).values()): continue
            su = next(s for (v, s) in V(U).values() if v == ONE)
            rest = tuple(d for d in DIG if d not in U)
            for k1 in range(2, len(rest) - 1):
                for T1 in itertools.combinations(rest, k1):
                    T2 = tuple(d for d in rest if d not in T1)
                    m = min(len(T1), len(T2))
                    if m <= 3 and SB[m] < theta:            # lemma: cannot reach e^theta
                        continue
                    if max(len(T1), len(T2)) > 4:
                        skipped.append((U, T1, T2)); continue        # would need 4-digit blocks
                    EL = side_entries(T1, theta); ER = side_entries(T2, theta)
                    ER.sort(key=lambda e: float(e[0].mid()))
                    keys = [float(e[0].mid()) for e in ER]
                    import bisect
                    for (lm, sg, a, x, sa, sx) in EL:
                        f = float(lm.mid())
                        tol = 1e-9 * max(1.0, abs(f)) + 1e-6      # >> float error of the midpoints
                        i = bisect.bisect_left(keys, f - tol)
                        while i < len(ER) and keys[i] <= f + tol:
                            lm2, sg2, c, y, sc, sy = ER[i]; i += 1
                            if sg2 != -sg or sg2 < 0: continue            # need mu_R = -mu_L > 0
                            if not lm.overlaps(lm2): continue           # rigorous: not an identity
                            ok = exact_identity(a, x, c, y)
                            rec = (lm2, su, sa, sx, sc, sy)
                            if ok: hits.append(rec)
                            elif ok is None: unresolved.append(rec)
                            else: unresolved.append(rec + ("FAILS_EXACT",))
    return hits, unresolved, skipped, SB

def stories_of(expr):
    """number of stories of the exponent chain: a -> 0, a^b -> 1 + stories(b) (powers in the base do not count)"""
    def strip(e):
        while e.startswith("-"): e = e[1:]
        while e.startswith("(") and matching(e, 0) == len(e) - 1: e = e[1:-1]
        return e
    def matching(e, i):
        d = 0
        for j in range(i, len(e)):
            d += e[j] == "("; d -= e[j] == ")"
            if d == 0: return j
        return -1
    e = strip(expr)
    d = 0
    for j, ch in enumerate(e):                      # top-level '^' (fully parenthesised strings)
        d += ch == "("; d -= ch == ")"
        if ch == "^" and d == 0: return 1 + stories_of(e[j + 1:])
    return 0

if __name__ == "__main__":
    DIG = sys.argv[1]
    theta = float(sys.argv[2])            # search identities with ln(ln|N|) >= theta
    hits, unres, skipped, SB = search(DIG, theta)
    print(f"digits {DIG}  theta = {theta}  ({time.time()-T0:.0f}s)")
    print(f"side bounds (ln of max |mu|): SB(2) <= {SB[2]:.6f}  SB(3) <= {SB[3]:.6f}")
    print(f"configurations needing blocks > 3 digits (must be empty): {len(skipped)}")
    print(f"unresolved overlaps (must be empty): {len(unres)}")
    for u in unres[:10]: print("   ", u)
    best = sorted(hits, key=lambda h: -float(h[0].mid()))
    tops = {}
    for h in best:
        tops.setdefault(h[0].mid().str(15), []).append(h)
    print(f"verified identities: {len(hits)}, distinct ln ln N values: {len(tops)}")
    def fmt(h): return f"({h[1]} + {h[2]}^{h[3]})^({h[4]}^{h[5]})"
    for key, hs in list(tops.items())[:int(sys.argv[3]) if len(sys.argv) > 3 else 10]:
        print(f"  ln(ln N) = {hs[0][0].str(20)}   ({len(hs)} writings)")
        for h in hs[:int(sys.argv[4]) if len(sys.argv) > 4 else 2]:
            print(f"        {fmt(h)}   stories L={1+stories_of(h[3])} R={1+stories_of(h[5])}")
    by = {}
    for h in hits:
        st = max(1 + stories_of(h[3]), 1 + stories_of(h[5]))
        if st not in by or float(h[0].mid()) > float(by[st][0].mid()): by[st] = h
    for st in sorted(by):
        print(f"  best with {st} stories: ln(ln N) = {by[st][0].str(12)}   {fmt(by[st])}")
