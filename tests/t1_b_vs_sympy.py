"""Test 1 (sympy values filtered by the range of the final value; with <= 3 digits no
in-range value needs an out-of-range intermediate): motor B's complete value sets against an independent sympy enumerator (sets of <= 3 digits).
Values compared by sign and 25 significant digits of ln|v| (zero separately)."""
import sys, itertools, random
sys.path.insert(0, ".")
import sympy as sp, mpmath
import motorb
mpmath.mp.dps = 50

def lits(S, ordered):
    out = set()
    perms = [S] if ordered else set(itertools.permutations(S))
    for perm in perms:
        t = "".join(perm)
        out.add(sp.Integer(int(t)))
        for i in range(len(t)):
            out.add(sp.Rational(int(t), 10 ** (len(t) - i)))
    return out

def spow(x, y):
    if x == 0: return sp.Integer(0) if y > 0 else None
    if x > 0: return sp.Pow(x, y)
    if y.is_Integer: return sp.Pow(x, y)
    return None

memo = {}
def VS(S, ordered):
    key = (S, ordered)
    if key in memo: return memo[key]
    out = set(lits(S, ordered))
    n = len(S)
    if ordered:
        parts = [(S[:k], S[k:]) for k in range(1, n)]
    else:
        parts = []
        for r in range(1, n):
            for S1 in itertools.combinations(S, r):
                parts.append((S1, tuple(d for d in S if d not in S1)))
    for S1, S2 in parts:
        for u in VS(S1, ordered):
            for v in VS(S2, ordered):
                out |= {u + v, u - v, u * v}
                if v != 0: out.add(u / v)
                p = spow(u, v)
                if p is not None: out.add(p)
    out |= {-v for v in out}
    memo[key] = out
    return out

def k_sym(v):
    v = abs(v)
    if v == 0: return "0"
    return mpmath.nstr(mpmath.mpf(str(sp.log(v).evalf(60))), 25)

def k_b(x):
    if x.is_zero(): return "0"
    b = x.ball(256)
    return mpmath.nstr(mpmath.mpf(abs(b).log().mid().str(55, radius=False)), 25)

sets = []
for S in [("1",), ("7",), ("1","2"), ("3","8"), ("5","9"), ("1","2","3"), ("2","5","9"), ("3","4","8"), ("1","6","7"), ("4","5","6"), ("7","8","9")]:
    sets.append(S)
total = 0
for ordered in (False, True):
    for S in sets:
        b = motorb.B(ordered=ordered)
        A = set()
        for v in VS(S, ordered):
            k = k_sym(v)
            if k != "0" and abs(float(k)) > 690.7755278982137:   # out of range (final value)
                continue
            A.add(k)
        Bk = {k_b(v) for (v, e) in b.full(S).values()}
        bad = A ^ Bk
        total += len(bad)
        print(("ord " if ordered else "any ") + "".join(S), "sympy", len(A), "B", len(Bk), "diff", len(bad), "unresolved", len(b.unresolved), flush=True)
        for z in list(A - Bk)[:3]: print("   only sympy", z)
        for z in list(Bk - A)[:3]: print("   only B    ", z)
print("TOTAL DIFFERENCES", total)
