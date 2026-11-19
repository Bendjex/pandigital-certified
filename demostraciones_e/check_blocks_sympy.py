"""
Independent cross-check of the block value sets V(S) used by e_F.py.
A second, separately written enumerator builds every expression of the grammar with sympy
(exact symbolic numbers) and the two value sets are compared (values keyed by 30 significant
digits; zero separately).  Any difference is printed.
"""
import itertools, sys
import sympy as sp
import mpmath
mpmath.mp.dps = 50
sys.path.insert(0, ".")
import e_F

def lits(S):
    out = set()
    for perm in set(itertools.permutations(S)):
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
def VS(S):
    if S in memo: return memo[S]
    out = set(lits(S))
    for r in range(1, len(S)):
        for S1 in itertools.combinations(S, r):
            S2 = tuple(d for d in S if d not in S1)
            for u in VS(S1):
                for v in VS(S2):
                    out |= {u + v, u - v, u * v}
                    if v != 0: out.add(u / v)
                    p = spow(u, v)
                    if p is not None: out.add(p)
    out |= {-v for v in out}
    memo[S] = out
    return out

def key_sympy(v):
    if v == 0: return "0"
    sgn = "-" if v < 0 else "+"
    L = sp.log(abs(v)).evalf(60)                     # log avoids building huge integers twice
    return sgn + mpmath.nstr(mpmath.mpf(str(L)), 28)

def key_eF(v):
    if v[0] == "Z": return "0"
    sgn = "-" if v[1] < 0 else "+"
    return sgn + mpmath.nstr(mpmath.mpf(e_F.lnabs(v).mid().str(55, radius=False)), 28)

total_bad = 0
for S in [("2",), ("2", "3"), ("1", "4"), ("5", "9"), ("2", "5", "9"), ("3", "4", "8"), ("1", "6", "7")]:
    A = {key_sympy(v) for v in VS(S)}
    B = {key_eF(v) for (v, s) in e_F.V(S).values()}
    bad = len(A ^ B); total_bad += bad
    print(S, "sympy:", len(A), " e_F:", len(B), " differences:", bad, flush=True)
    for z in list(A - B)[:5]: print("    only sympy", z)
    for z in list(B - A)[:5]: print("    only e_F  ", z)
print("TOTAL DIFFERENCES:", total_bad)
