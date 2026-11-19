"""soundness of fwd() with anchors and the third level: random exact rationals q, anchored values
q(1+d) with tiny d, values near 1, deep values; exact result in log space with mpmath"""
import sys, numpy as np, mpmath
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ivl
mpmath.mp.dps = 520
rng = np.random.default_rng(5)
def enc(x):
    a = float(x); return np.nextafter(a, -np.inf), np.nextafter(a, np.inf)
def mk(kind):
    """returns (L = ln value as mpf, Q7 of one element)"""
    if kind == 0:     # exact rational
        p, q = int(rng.integers(1, 60)), int(rng.choice([1, 2, 3, 4, 5, 8, 10, 20, 100]))
        from math import gcd
        g = gcd(p, q); p //= g; q //= g
        L = mpmath.log(mpmath.mpf(p) / q)
        lo, hi = enc(L)
        if p == q:
            return L, ivl.Q7(np.array([0.0]), np.array([0.0]), np.array([False]), np.array([True]))
        x = ivl.Q7(np.array([lo]), np.array([hi]), np.array([False]), np.array([False]))
        x.AP = np.array([p], np.int64); x.AQ = np.array([q], np.int64)
        return L, x
    if kind == 1:     # anchored q(1+d), d = +-10^-e
        p, q = int(rng.integers(1, 30)), int(rng.choice([1, 2, 5, 10]))
        from math import gcd
        g = gcd(p, q); p //= g; q //= g
        sgn = int(rng.choice([-1, 1])); e = float(rng.uniform(25, 400))
        d = sgn * mpmath.mpf(10) ** (-e)
        L = mpmath.log(mpmath.mpf(p) / q) + mpmath.log1p(d)
        lo, hi = enc(L)
        x = ivl.Q7(np.array([lo]), np.array([hi]), np.array([False]), np.array([False]))
        x.AP = np.array([p], np.int64); x.AQ = np.array([q], np.int64); x.DS = np.array([sgn], np.int8)
        dl = mpmath.log(abs(d)); a, b = enc(dl); x.DM = np.array([a]); x.DN = np.array([b])
        if p == q:   # near 1: also the second level
            m = mpmath.log(abs(L)); a, b = enc(m)
            x.S = np.array([1 if L > 0 else -1], np.int8); x.M = np.array([a]); x.N = np.array([b])
        return L, x
    if kind == 2:     # tiny / huge moderate
        L = mpmath.mpf(float(rng.choice([-1, 1]) * 10 ** rng.uniform(-1, 3)))
        lo, hi = enc(L); x = ivl.Q7(np.array([lo]), np.array([hi]), np.array([False]), np.array([False]))
        return L, x
    # deep: ln value = +-10^e, e in [310, 400]
    sg = int(rng.choice([-1, 1])); e = float(rng.uniform(310, 400)); L = sg * mpmath.mpf(10) ** e
    lo, hi = enc(L); x = ivl.Q7(np.array([lo]), np.array([hi]), np.array([False]), np.array([False]))
    m = mpmath.log(abs(L)); a, b = enc(m); x.S = np.array([sg], np.int8); x.M = np.array([a]); x.N = np.array([b])
    return L, x
def lse(a, b):
    m = max(a, b); return m + mpmath.log1p(mpmath.exp(min(a, b) - m))
def ldf(a, b):
    m = max(a, b); return m + mpmath.log(-mpmath.expm1(min(a, b) - m))
bad = 0; tested = 0
for trial in range(int(sys.argv[1]) if len(sys.argv) > 1 else 3000):
    k1, k2 = int(rng.integers(0, 4)), int(rng.integers(0, 4))
    a, x = mk(k1); b, y = mk(k2)
    f = int(rng.integers(0, 9))
    R = (np.where(x.DS == 0, x.AP, 0), np.where(x.DS == 0, x.AQ, 0), np.where(y.DS == 0, y.AP, 0), np.where(y.DS == 0, y.AQ, 0)) if f == 1 else None
    r, ok = ivl.fwd(f, x, y, R)
    if not ok[0]:
        continue
    if f == 1 and a == b:
        if not (r.Z[0] or r.A[0] == -np.inf): bad += 1; print("zero missed")
        continue
    try:
        Lv = [lse(a, b), ldf(a, b), a + b, a - b, b - a, mpmath.exp(b) * a, -mpmath.exp(b) * a,
              mpmath.exp(a) * b, -mpmath.exp(a) * b][f]
    except Exception:
        continue
    tested += 1
    ok1 = r.A[0] <= Lv <= r.B[0]
    ok2 = True
    if ok1 and r.S[0] != 0 and Lv != 0:
        mu = mpmath.log(abs(Lv)); ok2 = ((Lv > 0) == (r.S[0] > 0)) and (r.M[0] <= mu <= r.N[0])
        if ok2 and r.T3[0]:
            ok2 = (mu < 0) and (r.M3[0] <= mpmath.log(-mu) <= r.N3[0])
    ok3 = True
    if r.AQ[0] > 0:
        qv = mpmath.mpf(int(r.AP[0])) / int(r.AQ[0])
        if r.DS[0] == 0:
            ok3 = abs(mpmath.exp(Lv) / qv - 1) < mpmath.mpf(10) ** -60 if abs(Lv) < 1e5 else True
        else:
            dd = mpmath.expm1(Lv - mpmath.log(qv))
            if dd == 0 and 3 in (k1, k2):
                ok3 = r.DN[0] <= -1e300     # delta below the test's precision: only its bound can be checked
            else:
                ok3 = (dd != 0) and ((dd > 0) == (r.DS[0] > 0)) and (r.DM[0] <= mpmath.log(abs(dd)) <= r.DN[0])
    if not (ok1 and ok2 and ok3):
        bad += 1
        print("BAD", ivl.FORMS[f], "kinds", k1, k2, "a", mpmath.nstr(a, 8), "b", mpmath.nstr(b, 8), "true L", mpmath.nstr(Lv, 12),
              "enc", r.A[0], r.B[0], "S", r.S[0], r.M[0], r.N[0], "anchor", r.AP[0], r.AQ[0], r.DS[0], r.DM[0], r.DN[0], ok1, ok2, ok3)
print("tested", tested, "bad", bad)
