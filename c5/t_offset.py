"""the exact offset channel of ivl (03/11): u = X + c with X = (FP/FQ)^(FE/FD) a triple and c = OP/OQ.  Random exact
powers X (large and moderate, rational exponents), small exact rationals q, and chains X + q1, |q2 - (X + q1)|,
(X + q1) + q3; then |(X + c1) - (X + c2)| must be the exact rational |c1 - c2| (or the exact zero), and every
enclosure must contain the true value (mpmath, 400 digits); exactness is required for huge X (ln X > 40), where the channel is kept.  usage: t_offset.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
from fractions import Fraction as Fr
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ivl
mpmath.mp.dps = 400
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 53)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000


def exact_q(q):
    q = Fr(q); p, s = q.numerator, q.denominator
    L = mpmath.log(mpmath.mpf(p) / s) if p > 0 else None
    if p == 0:
        x = ivl.Q7(np.array([-np.inf]), np.array([-np.inf]), np.array([True]), np.array([False]))
    else:
        x = ivl.Q7(np.array([np.nextafter(float(L), -np.inf)]), np.array([np.nextafter(float(L), np.inf)]),
                   np.array([False]), np.array([p == s]))
        if p == s:
            x.A[:] = 0; x.B[:] = 0
    x.AP = np.array([p], np.int64); x.AQ = np.array([s], np.int64)
    return x


def lnv(v):
    """ln of an exact value given as ('pow', base, t, c): base^t + c (mpmath)"""
    _, b, t, c = v
    X = mpmath.exp(mpmath.mpf(t.numerator) / t.denominator * mpmath.log(mpmath.mpf(b.numerator) / b.denominator))
    return mpmath.log(abs(X + mpmath.mpf(c.numerator) / c.denominator))


def contains(r, L):
    return r.A[0] <= L <= r.B[0]


bases = [Fr(2), Fr(3), Fr(5), Fr(7), Fr(10), Fr(6), Fr(3, 2), Fr(9), Fr(4), Fr(1, 3), Fr(1, 10)]
expos = [Fr(1000), Fr(2000), Fr(32768), Fr(15625), Fr(100), Fr(40, 3), Fr(80, 3), Fr(1, 5), Fr(7, 2), Fr(3)]
bad = tested = exact_diffs = zeros = 0
for trial in range(N):
    b = bases[rng.integers(len(bases))]; t = expos[rng.integers(len(expos))]
    if rng.random() < 0.3:
        t = -t
    X = ivl.fwd(5, exact_q(b), exact_q(abs(t)))[0] if t > 0 else ivl.fwd(6, exact_q(b), exact_q(abs(t)))[0]
    if X.FE[0] == 0:
        continue                                     # small: an exact rational, not a triple
    qs = [Fr(int(rng.integers(0, 60)), int(rng.integers(1, 12))) for _ in range(3)]
    u1 = ivl.fwd(0, X, exact_q(qs[0]))[0]                    # X + q0
    u2 = ivl.fwd(0, exact_q(qs[1]), u1)[0]                   # (X + q0) + q1
    Xf = (b ** t if abs(t.numerator) * 64 < 10 ** 6 and t.denominator == 1 else None)
    big = X.A[0] > 50
    vals = [(X, ("pow", b, t, Fr(0))), (u1, ("pow", b, t, qs[0])), (u2, ("pow", b, t, qs[0] + qs[1]))]
    if big:
        u3 = ivl.fwd(1, exact_q(qs[2]), u1)[0]               # |q2 - (X + q0)| = X + q0 - q2 (X huge)
        vals.append((u3, ("pow", b, t, qs[0] - qs[2])))
    for (r, v) in vals:
        tested += 1
        if not contains(r, lnv(v)):
            bad += 1; print("BAD enclosure", b, t, v[3], r.A[0], r.B[0], mpmath.nstr(lnv(v), 15))
    # differences of two offsets of the same X
    for i in range(len(vals)):
        for j in range(i + 1, len(vals)):
            (ri, vi), (rj, vj) = vals[i], vals[j]
            d, ok = ivl.fwd(1, ri, rj)
            want = abs(vi[3] - vj[3])
            tested += 1
            if want == 0:
                zeros += 1
                # the channel is kept for huge X only (ln X > 40, where doubles cannot separate X + c1 from X + c2);
                # below, the enclosure must still contain the zero (A = -inf)
                if big and not d.Z[0]:
                    bad += 1; print("ZERO NOT DETECTED", b, t, vi[3], vj[3], d.A[0], d.B[0])
                elif not big and not (d.Z[0] or d.A[0] == -np.inf):
                    bad += 1; print("ZERO EXCLUDED", b, t, vi[3], vj[3], d.A[0], d.B[0])
            else:
                L = mpmath.log(mpmath.mpf(want.numerator) / want.denominator)
                if d.Z[0] or not contains(d, L):
                    bad += 1; print("BAD difference", b, t, vi[3], vj[3], d.Z[0], d.A[0], d.B[0], want)
                elif d.AQ[0] > 0 and d.DS[0] == 0:
                    exact_diffs += 1
                    if Fr(int(d.AP[0]), int(d.AQ[0])) != want:
                        bad += 1; print("BAD exact difference", d.AP[0], d.AQ[0], want)
print(f"tested {tested}, exact differences {exact_diffs}, exact zeros {zeros}, bad {bad}")
