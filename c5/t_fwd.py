"""soundness of fwd(): random point enclosures built from exact values (mpmath, 60 digits);
the exact result's ln must lie in [A, B], and ln|ln| in [M, N] when S != 0"""
import sys, numpy as np, mpmath
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ivl
mpmath.mp.dps = 60
rng = np.random.default_rng(3)
def enc(L):
    """enclosure of an exact ln value L (mpf): outward doubles"""
    a = float(L); lo = np.nextafter(a, -np.inf); hi = np.nextafter(a, np.inf)
    return lo, hi
bad = 0
for trial in range(int(sys.argv[1]) if len(sys.argv) > 1 else 4000):
    f = rng.integers(0, 9)
    vals = []
    for _ in range(2):
        kind = rng.integers(0, 5)
        if kind == 0:   L = mpmath.mpf(float(rng.normal(0, 3)))
        elif kind == 1: L = mpmath.mpf(int(rng.choice([-1, 1]))) * mpmath.mpf(10) ** float(rng.uniform(-300, -5))   # near 1
        elif kind == 2: L = mpmath.mpf(int(rng.choice([-1, 1]))) * mpmath.mpf(10) ** float(rng.uniform(3, 9))       # extreme
        elif kind == 4: L = mpmath.mpf(int(rng.choice([-1, 1]))) * mpmath.mpf(10) ** float(rng.uniform(310, 400))   # deep
        else:           L = mpmath.mpf(0)
        vals.append(L)
    x7 = []
    for L in vals:
        if L == 0:
            x7.append(ivl.Q7(np.array([0.0]), np.array([0.0]), np.array([False]), np.array([True])))
        else:
            lo, hi = enc(L)
            S = np.array([1 if L > 0 else -1], np.int8)
            m = mpmath.log(abs(L)); mlo, mhi = enc(m)
            x7.append(ivl.Q7(np.array([lo]), np.array([hi]), np.array([False]), np.array([False]), S, np.array([mlo]), np.array([mhi])))
    r, ok = ivl.fwd(int(f), x7[0], x7[1])
    if not ok[0]: continue
    a, b = vals
    def lse(x, y):
        m = max(x, y); return m + mpmath.log1p(mpmath.exp(min(x, y) - m))
    def ldf(x, y):
        m = max(x, y); return m + mpmath.log(-mpmath.expm1(min(x, y) - m))
    if f == 1 and a == b: continue
    Lv = [lse(a, b), ldf(a, b) if a != b else None, a + b, a - b, b - a,
          mpmath.exp(b) * a, -mpmath.exp(b) * a, mpmath.exp(a) * b, -mpmath.exp(a) * b][f]
    if not (r.A[0] <= Lv <= r.B[0]):
        bad += 1; print("BAD lambda form", ivl.FORMS[f], [mpmath.nstr(x, 8) for x in vals], "true", mpmath.nstr(Lv, 10), "enc", r.A[0], r.B[0], "S", r.S[0], r.M[0], r.N[0])
    elif r.S[0] != 0 and Lv != 0:
        mu = mpmath.log(abs(Lv))
        if (Lv > 0) != (r.S[0] > 0) or not (r.M[0] <= mu <= r.N[0]):
            bad += 1; print("BAD mu form", ivl.FORMS[f], [mpmath.nstr(x, 8) for x in vals], "true", mpmath.nstr(Lv, 10), mpmath.nstr(mu, 10), "S", r.S[0], r.M[0], r.N[0])
print("bad", bad)
