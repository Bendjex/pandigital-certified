"""
Independent re-derivation of the F result for digits 1-7 (the case where the brute force is small):
blocks from the separate sympy enumerator (check_blocks_sympy.VS), all splits U | (A,X) | (C,Y),
identity test X ln|A| + Y ln|C| = 0 at 80 digits with mpmath (tolerance 1e-60 relative), no
reuse of the matching code of e_F.py.  Prints every identity with ln N >= e^2.9 in the configurations the size lemma leaves (U = 1, sides 3 + 3).
"""
import itertools, sys
import sympy as sp, mpmath
mpmath.mp.dps = 80
sys.path.insert(0, ".")
from check_blocks_sympy import VS

DIG = tuple("1234567")
NUM = {}
def num(v):
    if v not in NUM: NUM[v] = mpmath.mpf(sp.N(v, 90).__str__())
    return NUM[v]

cache = {}
def side(T):
    T = tuple(sorted(T))
    if T in cache: return cache[T]
    out = []
    for r in range(1, len(T)):
        for SA in itertools.combinations(T, r):
            SX = tuple(d for d in T if d not in SA)
            for a in VS(SA):
                if a == 0 or abs(a) == 1: continue
                for x in VS(SX):
                    if x == 0: continue
                    if a < 0 and not x.is_Integer: continue
                    mu = num(x) * mpmath.log(abs(num(a)))
                    out.append((mu, a, x))
    cache[T] = out
    return out

found = set()
thr = mpmath.e ** 2.9
for ku in (1,):          # U = '1'; other U and side sizes are excluded by the size lemma (SB(2) < 2.9)
    for U in itertools.combinations(DIG, ku):
        if 1 not in VS(U): continue
        rest = tuple(d for d in DIG if d not in U)
        for k1 in (3,):
            for T1 in itertools.combinations(rest, k1):
                T2 = tuple(d for d in rest if d not in T1)
                R = sorted(side(T2), key=lambda t: t[0])
                keys = [t[0] for t in R]
                import bisect
                for (mL, a, x) in side(T1):
                    if mL >= 0: continue
                    target = -mL
                    if target < thr: continue
                    i = bisect.bisect_left(keys, target * (1 - mpmath.mpf(10) ** -60))
                    while i < len(R) and keys[i] <= target * (1 + mpmath.mpf(10) ** -60):
                        mR, c, y = R[i]; i += 1
                        found.add((U, str(a), str(x), str(c), str(y), mpmath.nstr(mpmath.log(mR), 15)))
for f in sorted(found, key=lambda z: z[-1]): print(f)
print("identities:", len(found), " distinct ln ln N:", sorted({f[-1] for f in found}))
