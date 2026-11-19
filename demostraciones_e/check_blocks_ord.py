"""Cross-checks of the in-order block sets Vo (e_F_ord.py):
 (1) runs 123, 234, 345, 456 against a separate in-order sympy enumerator (exact values, 28-digit keys);
 (2) for every 3-digit run, Vo(run) must be a subset of the any-order set V(sorted digits) of e_F.py."""
import itertools, sys
import sympy as sp, mpmath
mpmath.mp.dps = 50
sys.path.insert(0, ".")
import e_F, e_F_ord

memo = {}
def VSo(run):
    if run in memo: return memo[run]
    out = {sp.Integer(int(run))} | {sp.Rational(int(run), 10 ** (len(run) - i)) for i in range(len(run))}
    for m in range(1, len(run)):
        for u in VSo(run[:m]):
            for v in VSo(run[m:]):
                out |= {u + v, u - v, u * v}
                if v != 0: out.add(u / v)
                if u > 0 or (u == 0 and v > 0) or (u < 0 and v.is_Integer): out.add(sp.Pow(u, v))
    out |= {-v for v in out}
    memo[run] = out
    return out

def k_sym(v):
    if v == 0: return "0"
    return ("-" if v < 0 else "+") + mpmath.nstr(mpmath.mpf(str(sp.log(abs(v)).evalf(60))), 28)
def k_ef(v):
    if v[0] == "Z": return "0"
    return ("-" if v[1] < 0 else "+") + mpmath.nstr(mpmath.mpf(e_F.lnabs(v).mid().str(55, radius=False)), 28)

bad = 0
for run in ["123", "234", "345", "456"]:
    A = {k_sym(v) for v in VSo(run)}; B = {k_ef(v) for (v, s) in e_F_ord.Vo(run).values()}
    bad += len(A ^ B); print(run, "sympy", len(A), "e_F_ord", len(B), "differences", len(A ^ B), flush=True)
for i in range(7):
    run = "123456789"[i:i + 3]
    A = {k_ef(v) for (v, s) in e_F_ord.Vo(run).values()}
    B = {k_ef(v) for (v, s) in e_F.V(tuple(sorted(run))).values()}
    print(run, "in-order", len(A), "any-order", len(B), "in-order values missing from any-order:", len(A - B))
    bad += len(A - B)
print("TOTAL PROBLEMS:", bad)
