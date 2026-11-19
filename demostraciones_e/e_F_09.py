"""
B1 for the 0-9 variant (any order): is (1+.2^(9^70))^(5^(3^(84/.6))), N = 5^(3^140), optimal in F?

Grammar: as e_F.py, PERMISSIVE about zeros: any permutation of the block's digits is a literal,
leading zeros included ("05" = 5, "0.5", ".05").  A maximum under the permissive grammar is also the
maximum under any stricter convention in which the formula itself is valid (ours uses 70 and .6 only).

Size lemma: |mu| <= SB(3) ~ e^21.5 for a side of <= 3 digits, while our record has ln ln N ~ 154.28.
So U has 1 digit and the sides are 4+5 / 5+4, or U has 2 digits and the sides are 4+4; U with >= 3
digits leaves a side of <= 3 digits.  Blocks of up to 3 digits are enumerated exactly (e_F.V); for the
5-digit side the 4-digit blocks are enumerated only where they can matter:
   * X-role (A one digit): ln|X| >= theta - ln(ln 10) - 1
   * A-role (X one digit): ln|ln|A|| >= theta - ln 9 - 1
A 4-digit value is op(u, v) over a split of the block; the float pre-estimates below are upper bounds
(with margin 1) of the relevant magnitude, so nothing that can reach the threshold is discarded.
"""
import sys, itertools, bisect, math, time
from functools import lru_cache
import e_F
from e_F import V, Vinfo, add, neg, mul, div, pw, lnabs, ONE
from flint import arb
T0 = time.time()
LN_LN10 = math.log(math.log(10)); LN9 = math.log(9)

def fl(L):
    try: return float(L.mid())
    except OverflowError: return math.inf if L > 0 else -math.inf

@lru_cache(None)
def lninfo(S):
    """values of V(S) with float ln|v| (+-inf allowed) and float ln|ln|v|| (or -inf)"""
    out = []
    for (v, s, L, f, isint) in Vinfo(S):
        f = fl(L)
        lal = -math.inf if L.contains(0) else fl(abs(L).log())
        out.append((v, s, f, lal))
    return out

@lru_cache(None)
def V4big(S, mode, T):
    """4-digit block values that can matter. mode 'X': ln|v| >= T ; mode 'A': ln|ln|v|| >= T"""
    out = {}
    def put(v, s):
        if v is None or v[0] == "Z": return
        L = lnabs(v)
        if mode == "X":
            if L.upper() < T: return
        else:
            if L.contains(0) or abs(L).log().upper() < T: return
        k = v if v[0] != "T" else ("T", s)
        if k not in out: out[k] = (v, s)
    n = len(S)
    for r in range(1, n):
        for S1 in itertools.combinations(S, r):
            S2 = tuple(d for d in S if d not in S1)
            I1, I2 = lninfo(S1), lninfo(S2)
            for (u, su, fu, lu) in I1:
                for (v, sv, fv, lv) in I2:
                    # ---- power u^v : ln|u^v| = v ln|u| ;  ln|ln|u^v|| = ln|v| + ln|ln|u||
                    llp = fv + lu                                   # ln|ln|u^v||  (float, may be inf/-inf)
                    ok = (llp >= T - 1) if mode == "A" else (llp >= math.log(max(T - 1, 1e-300)) - 1 if T > 1 else True)
                    if ok: put(pw(u, v), f"({su}^{sv})")
                    # ---- sums, products, quotients: |ln| of the result <= |ln u| + |ln v| + ln 2
                    if mode == "X":
                        if abs(fu) + abs(fv) + 1 >= T - 1:              # |ln(u op v)| <= |ln u| + |ln v| + ln 2
                            for f, s in ((add(u, v), f"({su}+{sv})"), (add(u, neg(v)), f"({su}-{sv})"),
                                         (mul(u, v), f"({su}*{sv})"), (div(u, v), f"({su}/{sv})")):
                                put(f, s)
                    else:
                        bound = max(abs(fu), abs(fv)) + 1 if math.isfinite(fu) and math.isfinite(fv) else math.inf
                        if bound == math.inf or math.log(bound * 2 + 1) >= T - 1:
                            for f, s in ((add(u, v), f"({su}+{sv})"), (add(u, neg(v)), f"({su}-{sv})"),
                                         (mul(u, v), f"({su}*{sv})"), (div(u, v), f"({su}/{sv})")):
                                put(f, s)
    for k, (v, s) in list(out.items()):
        put(neg(v), "-" + s)
    return out

def info4(S, mode, T):
    res = []
    for (v, s) in V4big(S, mode, T).values():
        L = lnabs(v)
        isint = v[0] == "R" and all(e.denominator == 1 and e > 0 for _, e in v[2])
        res.append((v, s, L, fl(L), isint))
    res.sort(key=lambda t: -t[3])
    return res

def side_entries(T, theta):
    res = []
    T = tuple(sorted(T))
    for r in range(1, len(T)):
        for SA in itertools.combinations(T, r):
            SX = tuple(d for d in T if d not in SA)
            if len(SA) == 4: VA = info4(SA, "A", theta - LN9 - 1)
            elif len(SA) <= 3: VA = Vinfo(SA)
            else: raise ValueError("block too large")
            if len(SX) == 4: VX = info4(SX, "X", theta - LN_LN10 - 1)
            elif len(SX) <= 3: VX = Vinfo(SX)
            else: raise ValueError("block too large")
            if len(SA) == 4 and len(SX) != 1: continue        # A-role 4-digit blocks only pair with 1-digit X
            if len(SX) == 4 and len(SA) != 1: continue        # (other pairings have a <= 1-digit partner by size)
            for (a, sa, La, fa, _) in VA:
                if a[0] == "R" and not a[2]: continue
                if La.contains(0): raise ValueError("cannot decide sign of ln|A|")
                lla = abs(La).log(); flla = fl(lla)
                for (x, sx, Lx, fx, isint) in VX:
                    if fx + flla < theta - 1e-3: break
                    if not (a[1] > 0 or isint): continue
                    lmu = Lx + lla
                    if lmu.upper() < theta: continue
                    res.append((lmu, x[1] * (1 if La > 0 else -1), a, x, sa, sx))
    return res

def search(DIG, theta):
    DIG = tuple(DIG)
    SB3 = max(e_F.side_bound(T) for T in itertools.combinations(DIG, 3))
    assert SB3 < theta, "size lemma does not exclude 3-digit sides"
    hits, unres = [], []
    configs = []
    for ku in (1, 2):
        for U in itertools.combinations(DIG, ku):
            if not any(v == ONE for (v, _) in V(U).values()): continue
            su = next(s for (v, s) in V(U).values() if v == ONE)
            rest = tuple(d for d in DIG if d not in U)
            for k1 in range(4, len(rest) - 3):
                for T1 in itertools.combinations(rest, k1):
                    configs.append((su, T1, tuple(d for d in rest if d not in T1)))
    print(f"SB(3) <= e^{SB3:.4f};  configurations: {len(configs)}", flush=True)
    cache = {}
    def SE(T):
        T = tuple(sorted(T))
        if T not in cache:
            e = side_entries(T, theta); e.sort(key=lambda z: fl(z[0])); cache[T] = e
        return cache[T]
    for ci, (su, T1, T2) in enumerate(configs):
        EL, ER = SE(T1), SE(T2)
        keys = [fl(e[0]) for e in ER]
        for (lm, sg, a, x, sa, sx) in EL:
            f = fl(lm)
            if not math.isfinite(f): continue                 # ln|mu| beyond double: larger than any 4-digit side allows
            tol = 1e-9 * max(1.0, abs(f)) + 1e-6
            i = bisect.bisect_left(keys, f - tol)
            while i < len(ER) and keys[i] <= f + tol:
                lm2, sg2, c, y, sc, sy = ER[i]; i += 1
                if sg2 != -sg or sg2 < 0: continue
                if not lm.overlaps(lm2): continue
                ok = e_F.exact_identity(a, x, c, y)
                rec = (lm2, su, sa, sx, sc, sy)
                (hits if ok else unres).append(rec)
        if ci % 50 == 0: print(f"  {ci}/{len(configs)}  {time.time()-T0:.0f}s  hits {len(hits)} unresolved {len(unres)}", flush=True)
    return hits, unres

if __name__ == "__main__":
    theta = float(sys.argv[1]) if len(sys.argv) > 1 else 154.0
    hits, unres = search("0123456789", theta)
    print(f"digits 0-9, theta = {theta}  ({time.time()-T0:.0f}s)")
    print(f"unresolved overlaps (must be empty): {len(unres)}")
    for u in unres[:10]: print("   ", u)
    tops = {}
    for h in sorted(hits, key=lambda h: -fl(h[0])):
        tops.setdefault(h[0].mid().str(15, radius=False), []).append(h)
    print(f"verified identities: {len(hits)}, distinct ln ln N values: {len(tops)}")
    for key, hs in list(tops.items())[:8]:
        print(f"  ln(ln N) = {hs[0][0].str(20)}   ({len(hs)} writings)")
        for h in hs[:3]:
            print(f"        ({h[1]} + {h[2]}^{h[3]})^({h[4]}^{h[5]})")
