"""fwd() on values beyond the double range of ln (|ln v| > DBL_MAX, known by the second level only) whose second
levels are CLOSE (mu_a - mu_b in [0, 3]): the dominance rules of sums, differences, products and quotients (C.4.3,
C.4.4, C.4.2) with the strict separation mu_b < mu_a of 03/11 instead of mu_b < mu_a - 1.  The exact ln of the result
(mpmath, 80 digits) must lie in [A, B], and its ln|ln| in [M, N] when S != 0.  Also counts how often the result is
informative (S != 0).  usage: t_fwd2.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ivl
mpmath.mp.dps = 80
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 17)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000


def enc(L):
    a = float(L)
    return np.nextafter(a, -np.inf), np.nextafter(a, np.inf)


def q7(L):
    lo, hi = enc(L)
    m = mpmath.log(abs(L)); mlo, mhi = enc(m)
    return ivl.Q7(np.array([lo]), np.array([hi]), np.array([False]), np.array([False]),
                  np.array([1 if L > 0 else -1], np.int8), np.array([mlo]), np.array([mhi]))


def lse(x, y):
    m = max(x, y); return m + mpmath.log1p(mpmath.exp(min(x, y) - m))


def ldf(x, y):
    m = max(x, y); return m + mpmath.log(-mpmath.expm1(min(x, y) - m))


bad = informative = tested = 0
for t in range(N):
    f = int(rng.integers(0, 5))
    mu_a = mpmath.mpf(float(rng.uniform(709.0, 1500.0)))
    gap = mpmath.mpf(float(rng.choice([rng.uniform(0, 3), 10 ** rng.uniform(-14, 0)])))
    mu_b = mu_a - gap if rng.random() < 0.85 else mpmath.mpf(float(rng.uniform(0, 50)))
    a = (1 if rng.random() < 0.5 else -1) * mpmath.exp(mu_a)
    b = (1 if rng.random() < 0.5 else -1) * mpmath.exp(mu_b)
    if rng.random() < 0.5:
        a, b = b, a
    x, y = q7(a), q7(b)
    r, ok = ivl.fwd(f, x, y)
    if not ok[0] or (f == 1 and a == b):
        continue
    Lv = [lse(a, b), ldf(a, b), a + b, a - b, b - a][f]
    tested += 1
    if not (r.A[0] <= Lv <= r.B[0]):
        bad += 1; print("BAD lambda", ivl.FORMS[f], mpmath.nstr(mu_a, 12), mpmath.nstr(mu_b, 12), "true", mpmath.nstr(Lv, 10), "enc", r.A[0], r.B[0])
        continue
    if r.S[0] != 0:
        informative += 1
        if Lv == 0 or (Lv > 0) != (r.S[0] > 0):
            bad += 1; print("BAD sign", ivl.FORMS[f], mpmath.nstr(mu_a, 12), mpmath.nstr(mu_b, 12), "S", r.S[0]); continue
        mu = mpmath.log(abs(Lv))
        if not (r.M[0] <= mu <= r.N[0]):
            bad += 1; print("BAD mu", ivl.FORMS[f], mpmath.nstr(mu_a, 15), mpmath.nstr(mu_b, 15), "true", mpmath.nstr(mu, 15), "S", r.S[0], r.M[0], r.N[0])
print(f"tested {tested}, informative (S != 0) {informative}, bad {bad}")
