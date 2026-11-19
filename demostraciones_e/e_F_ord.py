"""
B1 for the IN-ORDER variant: exhaustive, rigorous search of F_ord for e with digits 1..n in increasing order.

    F_ord : (U + s*A^X)^(t*C^Y)   or   (s*A^X + U)^(t*C^Y),    s, t in {+1, -1},
            where the blocks are consecutive runs of the digit string, in the written order
            (U, A, X, C, Y) or (A, X, U, C, Y); U = 1 exactly; s*A^X * t*C^Y = 1 exactly.

In-order grammar: a block on the run d_i..d_j is a literal (the digits in this order, optionally with one
decimal point before any digit), or u op v with u on d_i..d_m and v on d_(m+1)..d_j (op in + - * / ^,
left operand first, no reordering), or -u.  Same domain rules as e_F.py.

Method as in e_F.py (exact radicals, Arb screening, exact verification, size lemma), with the side
bounds SB(k) computed for runs of k digits.  Blocks of up to 4 digits are enumerated (sides up to 5).
"""
import sys, itertools, bisect, time
from fractions import Fraction as F
from functools import lru_cache
import e_F
from e_F import R_from_frac, add, neg, mul, div, pw, lnabs, is_rat, ONE
from flint import arb
T0 = time.time()

@lru_cache(None)
def Vo(run):
    """value set of the in-order block on the digit string `run`"""
    out = {}
    def put(v, s):
        if v is None: return
        k = v if v[0] != "T" else ("T", s)
        if k not in out: out[k] = (v, s)
    t = run
    put(R_from_frac(F(int(t))), t)
    for i in range(len(t)):
        put(R_from_frac(F(int(t), 10 ** (len(t) - i))), t[:i] + "." + t[i:])
    for m in range(1, len(run)):
        for (u, su) in Vo(run[:m]).values():
            for (v, sv) in Vo(run[m:]).values():
                for f, s in ((add(u, v), f"({su}+{sv})"), (add(u, neg(v)), f"({su}-{sv})"),
                             (mul(u, v), f"({su}*{sv})"), (div(u, v), f"({su}/{sv})"),
                             (pw(u, v), f"({su}^{sv})")):
                    put(f, s)
    for k, (v, s) in list(out.items()):
        put(neg(v), "-" + s)
    return out

@lru_cache(None)
def Voinfo(run):
    out = []
    for (v, s) in Vo(run).values():
        if v[0] == "Z": continue
        L = lnabs(v)
        isint = v[0] == "R" and all(e.denominator == 1 and e > 0 for _, e in v[2])
        out.append((v, s, L, float(L.mid()), isint))
    out.sort(key=lambda t: -t[3])
    return out

@lru_cache(None)
def side_entries(Arun, Xrun, theta):
    res = []
    VX = Voinfo(Xrun)
    for (a, sa, La, fa, _) in Voinfo(Arun):
        if a[0] == "R" and not a[2]: continue
        if La.contains(0): raise ValueError("cannot decide sign of ln|A|")
        lla = abs(La).log(); flla = float(lla.mid())
        for xi in VX:
            (x, sx, Lx, fx, isint) = xi
            if fx + flla < theta - 1e-3: break
            if not (a[1] > 0 or isint): continue
            lmu = Lx + lla
            if lmu.upper() < theta: continue
            res.append((lmu, x[1] * (1 if La > 0 else -1), a, x, sa, sx))
    res.sort(key=lambda e: float(e[0].mid()))
    return res

@lru_cache(None)
def side_bound(Arun, Xrun):
    best = -1e300
    VX = Voinfo(Xrun); mx = VX[0][2]
    for (a, _, La, _, _) in Voinfo(Arun):
        if a[0] == "R" and not a[2]: continue
        if La.contains(0): raise ValueError("cannot decide sign of ln|A|")
        best = max(best, float((mx + abs(La).log()).upper()) * (1 + 1e-12) + 1e-12)
    return best

def SB_run(run):
    """max over splits of the run into A, X"""
    return max(side_bound(run[:m], run[m:]) for m in range(1, len(run)))

def search(n, theta, maxblock=4):
    D = "123456789"[:n]
    hits, unres, skipped = [], [], []
    layouts = []
    # cut points 0 < c1 < c2 < c3 < c4 < n : five consecutive non-empty runs
    for c in itertools.combinations(range(1, n), 4):
        r = [D[:c[0]], D[c[0]:c[1]], D[c[1]:c[2]], D[c[2]:c[3]], D[c[3]:]]
        layouts.append(("U+L", r[0], r[1], r[2], r[3], r[4]))      # U, A, X, C, Y
        layouts.append(("L+U", r[2], r[0], r[1], r[3], r[4]))      # A, X, U, C, Y  (stored as U, A, X, C, Y)
    for (kind, U, A, X, C, Y) in layouts:
        short, long_ = sorted([A + X, C + Y], key=len)
        if SB_run(short) < theta: continue                            # size lemma: the shorter side bounds |mu|
        if len(long_) > maxblock + 1 or max(len(A), len(X), len(C), len(Y)) > maxblock:
            skipped.append((kind, U, A, X, C, Y)); continue
        bL, bR = SB_run(A + X), SB_run(C + Y)
        if min(bL, bR) < theta: continue
        if len(U) > maxblock:
            skipped.append((kind, U, A, X, C, Y)); continue
        if not any(v == ONE for (v, _) in Vo(U).values()): continue
        su = next(s for (v, s) in Vo(U).values() if v == ONE)
        cap = min(bL, bR)                                             # |mu| of an identity is bounded by both sides
        EL = [e for e in side_entries(A, X, theta) if e[0].lower() <= cap]
        ER = [e for e in side_entries(C, Y, theta) if e[0].lower() <= cap]
        keys = [float(e[0].mid()) for e in ER]
        for (lm, sg, a, x, sa, sx) in EL:
            f = float(lm.mid()); tol = 1e-9 * max(1.0, abs(f)) + 1e-6
            i = bisect.bisect_left(keys, f - tol)
            while i < len(ER) and keys[i] <= f + tol:
                lm2, sg2, c, y, sc, sy = ER[i]; i += 1
                if sg2 != -sg or sg2 < 0: continue
                if not lm.overlaps(lm2): continue
                ok = e_F.exact_identity(a, x, c, y)
                rec = (lm2, kind, su, sa, sx, sc, sy)
                if ok: hits.append(rec)
                else: unres.append(rec + (ok,))
    return hits, unres, skipped

if __name__ == "__main__":
    n = int(sys.argv[1]); theta = float(sys.argv[2])
    hits, unres, skipped = search(n, theta)
    print(f"digits 1..{n} in order, theta = {theta}  ({time.time()-T0:.0f}s)")
    print(f"configurations needing blocks > 4 digits (must be empty): {len(skipped)}  {skipped[:3]}")
    print(f"unresolved overlaps (must be empty): {len(unres)}")
    for u in unres[:10]: print("   ", u)
    tops = {}
    for h in sorted(hits, key=lambda h: -float(h[0].mid())):
        tops.setdefault(h[0].mid().str(15, radius=False), []).append(h)
    print(f"verified identities: {len(hits)}, distinct ln ln N values: {len(tops)}")
    for key, hs in list(tops.items())[:int(sys.argv[3]) if len(sys.argv) > 3 else 8]:
        print(f"  ln(ln N) = {hs[0][0].str(20)}   ({len(hs)} writings)")
        for h in hs[:int(sys.argv[4]) if len(sys.argv) > 4 else 2]:
            k, su, sa, sx, sc, sy = h[1:7]
            L = f"{sa}^{sx}"
            print("        " + (f"({su} + {L})" if k == "U+L" else f"({L} + {su})") + f"^({sc}^{sy})")
