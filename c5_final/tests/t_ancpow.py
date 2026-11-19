"""soundness of the rule C.4.6(8) of ivl.anchors (exact or anchored base raised to an exact integer):
random bases q = p/r (exact, or anchored q(1 + d) with |d| from 1e-9 down to 1e-600, both signs), exact integer
exponents k (also near the 2^60 limit of ipow), forms s^l, s^-l, l^s, l^-s; the exact value is computed with
mpmath at 1400 digits.  Checked: the anchor rational is q^(+-k) exactly, sign(delta') and ln|delta'| in [DM, DN],
and ln u in [A, B].  Also checks the cases the rule must leave alone (q = 1, k = 0, overflow).
usage: t_ancpow.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ivl
mpmath.mp.dps = 1400
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 17)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 4000


def enc(v):
    a = float(v)
    return np.nextafter(a, -np.inf), np.nextafter(a, np.inf)


def q7_exact(p, r):
    L = mpmath.log(mpmath.mpf(p) / r)
    lo, hi = enc(L)
    x = ivl.Q7(np.array([lo]), np.array([hi]), np.array([False]), np.array([p == r]))
    if p == r:
        x.A[:] = 0; x.B[:] = 0
    x.AP = np.array([p], np.int64); x.AQ = np.array([r], np.int64)
    return mpmath.mpf(p) / r, x


def q7_anchored(p, r, d):
    v = mpmath.mpf(p) / r * (1 + d)
    L = mpmath.log(v)
    lo, hi = enc(L)
    x = ivl.Q7(np.array([lo]), np.array([hi]), np.array([False]), np.array([False]))
    x.AP = np.array([p], np.int64); x.AQ = np.array([r], np.int64)
    x.DS = np.array([1 if d > 0 else -1], np.int8)
    a, b = enc(mpmath.log(abs(d))); x.DM = np.array([a]); x.DN = np.array([b])
    if p == r:                                           # anchored at 1: also the second level, as in practice
        m = mpmath.log(abs(L)); a, b = enc(m)
        x.S = np.array([1 if L > 0 else -1], np.int8); x.M = np.array([a]); x.N = np.array([b])
    return v, x


bad = tested = anchored = exact = 0
for trial in range(N):
    p, r = int(rng.integers(1, 60)), int(rng.choice([1, 2, 3, 4, 5, 7, 8, 9, 10, 16, 20, 25, 100]))
    g = gcd(p, r); p //= g; r //= g
    kind = rng.random()
    if kind < 0.35:
        k = int(rng.integers(1, 8))
    elif kind < 0.7:
        k = int(rng.integers(8, 70))
    else:                                                # near the overflow limit 2^60 of ipow
        big = max(p, r, 2)
        k = max(1, int(60 / np.log2(big)) + int(rng.integers(-2, 2)))
    if rng.random() < 0.3:
        bv, base = q7_exact(p, r); dd = 0
    else:
        e = float(rng.uniform(9, 600)) if rng.random() < 0.8 else float(rng.uniform(8.8, 9.2))
        dd = int(rng.choice([-1, 1])) * mpmath.mpf(10) ** (-e)
        bv, base = q7_anchored(p, r, dd)
    ev, expo = q7_exact(k, 1)
    f = int(rng.choice([5, 6, 7, 8]))
    x, y = (base, expo) if f in (5, 6) else (expo, base)
    r7, ok = ivl.fwd(f, x, y)
    if not ok[0]:
        continue
    sg = -1 if f in (6, 8) else 1
    u = bv ** (sg * k)
    Lu = mpmath.log(u)
    tested += 1
    good = r7.A[0] <= Lu <= r7.B[0] if not (r7.O[0]) else u == 1
    if r7.AQ[0] > 0:
        qk = (Fr(p, r) ** (sg * k))
        if (int(r7.AP[0]), int(r7.AQ[0])) != (qk.numerator, qk.denominator):
            good = False; print("wrong anchor rational", p, r, k, f, r7.AP[0], r7.AQ[0])
        if r7.DS[0] == 0:
            exact += 1
            good &= (dd == 0)
        else:
            anchored += 1
            dl = u / (mpmath.mpf(qk.numerator) / qk.denominator) - 1
            if dl == 0:
                good = False
            else:
                good &= ((dl > 0) == (r7.DS[0] > 0)) and (r7.DM[0] <= mpmath.log(abs(dl)) <= r7.DN[0])
    if not good:
        bad += 1
        print("BAD form", ivl.FORMS[f], "base", p, "/", r, "d", mpmath.nstr(dd, 5), "k", k, "enc", r7.A[0], r7.B[0],
              "anchor", r7.AP[0], r7.AQ[0], r7.DS[0], r7.DM[0], r7.DN[0])
print(f"tested {tested} (anchored results {anchored}, exact results {exact}) bad {bad}")
