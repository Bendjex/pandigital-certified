"""
Rigorous vectorised interval arithmetic on magnitudes, in the log domain (motor X, C5).

An enclosure is (A, B, Z, O):  Z: the value is exactly 0;  O: the value is exactly 1;  otherwise the
magnitude u satisfies  e^A <= u <= e^B  with e^-inf = 0 (A = -inf means "possibly 0 or tiny").
All bounds are outward: every + - of doubles is followed by pred/succ, exp and ln bounds come from
motor A's exp_lo/exp_hi/ln_lo/ln_hi (CORE-MATH pow; Lemma 4 of C4).  NaN results (inf - inf,
0 * inf) are replaced by -inf for lower bounds and +inf for upper bounds.

Forms (Lemma 3), with s from S1 and l from S2:
  0 s+l   1 |s-l|   2 s*l   3 s/l   4 l/s   5 s^l   6 s^-l   7 l^s   8 l^-s
fwd(f, s, l) gives an enclosure of every value of the form over the two enclosures (Lemma X2);
pre(f, s, w) gives windows (in ln l) containing every l with fwd(f, s, l) meeting w (Lemma X3).
"""
import numpy as np
from math import gcd
from fractions import Fraction as Fr
import xl

import os
INF = np.inf
pred, succ = xl.pred, xl.succ
ANCPOW = os.environ.get("IVL_ANCPOW", "1") == "1"      # rule C.4.6(8): exact/anchored base ^ exact integer
COUNT = {}                                             # instrumentation of the corrected rules (see bump)


def bump(key, k=1):
    COUNT[key] = COUNT.get(key, 0) + int(k)
FORMS = ["s+l", "|s-l|", "s*l", "s/l", "l/s", "s^l", "s^-l", "l^s", "l^-s"]
ORD_LEFT = (0, 1, 2, 3, 5, 6)      # in order, S1 is the left operand
ORD_RIGHT = (0, 1, 2, 4, 7, 8)     # in order, S1 is the right operand
# range: |v| in [1e-300, 1e300]; LN_HI <= ln(1e300) <= ..., safe thresholds (superset tests)
LN_HI_LO = float(xl.ln_lo(np.array([np.nextafter(1e300, 0)]))[0])     # <= ln 1e300
LN_HI_HI = float(xl.ln_hi(np.array([np.nextafter(1e300, INF)]))[0])   # >= ln 1e300
LN_LO_HI = float(xl.ln_hi(np.array([np.nextafter(1e-300, INF)]))[0])  # >= ln 1e-300
LN_LO_LO = float(xl.ln_lo(np.array([np.nextafter(1e-300, 0)]))[0])    # <= ln 1e-300


def fixlo(x):
    x = np.asarray(x, dtype=np.float64)
    return np.where(np.isnan(x), -INF, x)


def fixhi(x):
    x = np.asarray(x, dtype=np.float64)
    return np.where(np.isnan(x), INF, x)


def addlo(a, b):
    with np.errstate(invalid="ignore"):
        return fixlo(pred(fixlo(a + b)))


def addhi(a, b):
    with np.errstate(invalid="ignore"):
        return fixhi(succ(fixhi(a + b)))


def exp_lo(x):
    return xl.exp_lo(x)


def exp_hi(x):
    return xl.exp_hi(x)


def ln_lo(x):
    x = np.asarray(x, dtype=np.float64)
    out = np.full(x.shape, -INF)
    m = x > 0
    if m.any():
        out[m] = xl.ln_lo(x[m])
    # near 1: ln x >= u - u^2 with u = x - 1 exact (|u| < 2^-9)
    nr = (x > 1 - 2.0 ** -9) & (x < 1 + 2.0 ** -9)
    if nr.any():
        u = x[nr] - 1.0
        out[nr] = np.maximum(out[nr], pred(u - succ(u * u)))
    return out


def ln_hi(x):
    x = np.asarray(x, dtype=np.float64)
    out = np.full(x.shape, -INF)
    m = x > 0
    if m.any():
        out[m] = xl.ln_hi(x[m])
    out[np.isinf(x) & (x > 0)] = INF
    nr = (x > 1 - 2.0 ** -9) & (x < 1 + 2.0 ** -9)       # ln x <= x - 1 (exact)
    if nr.any():
        out[nr] = np.minimum(out[nr], x[nr] - 1.0)
    return out


def log1p_lo(t):
    """lower bound of ln(1+t), t >= 0 (t is a lower bound already)"""
    t = np.maximum(np.asarray(t, dtype=np.float64), 0.0)
    small = t < 2.0 ** -20
    out = np.empty(t.shape)
    # ln(1+t) >= t - t^2/2 >= t - t^2  (t >= 0)
    ts = t[small]
    out[small] = pred(pred(ts - succ(ts * ts)))
    tb = t[~small]
    out[~small] = ln_lo(pred(1.0 + tb))
    return out


def log1p_hi(t):
    """upper bound of ln(1+t), t >= 0 (t is an upper bound already)"""
    t = np.maximum(np.asarray(t, dtype=np.float64), 0.0)
    small = t < 2.0 ** -20
    out = np.empty(t.shape)
    out[small] = succ(t[small])                   # ln(1+t) <= t
    tb = t[~small]
    out[~small] = ln_hi(succ(1.0 + tb))
    return out


def log1m_lo(t):
    """lower bound of ln(1-t) for 0 <= t <= 1 (t is an upper bound of the true t); -inf if t >= 1"""
    t = np.asarray(t, dtype=np.float64)
    out = np.full(t.shape, -INF)
    small = (t >= 0) & (t < 2.0 ** -20)
    ts = t[small]
    # ln(1-t) >= -t - t^2 for 0 <= t <= 1/2
    out[small] = pred(-ts - succ(ts * ts))
    mid = (t >= 2.0 ** -20) & (t < 1)
    out[mid] = ln_lo(pred(1.0 - t[mid]))
    neg = t < 0
    out[neg] = 0.0
    return out


def log1m_hi(t):
    """upper bound of ln(1-t) for t >= 0 (t is a lower bound of the true t)"""
    t = np.maximum(np.asarray(t, dtype=np.float64), 0.0)
    out = np.empty(t.shape)
    small = t < 2.0 ** -20
    out[small] = pred(-t[small]) if False else -t[small]   # ln(1-t) <= -t (exact double negation)
    big = ~small
    out[big] = ln_hi(succ(1.0 - t[big]))
    return out


# ---------------------------------------------------------------- forward forms (Lemma X2)
class Q7:
    """a batch of enclosures: A <= ln u <= B, Z (u = 0), O (u = 1), and the second channel for values
    near 1:  S = sign of ln u (0: unknown), M <= ln|ln u| <= N (meaningful when S != 0)"""
    __slots__ = ("A", "B", "Z", "O", "S", "M", "N", "AP", "AQ", "DS", "DM", "DN", "T3", "M3", "N3", "FP", "FQ", "FE", "FD", "OP", "OQ")
    FIELDS = ("A", "B", "Z", "O", "S", "M", "N", "AP", "AQ", "DS", "DM", "DN", "T3", "M3", "N3", "FP", "FQ", "FE", "FD", "OP", "OQ")

    def __init__(self, A, B, Z, O, S=None, M=None, N=None, AP=None, AQ=None, DS=None, DM=None, DN=None,
                 T3=None, M3=None, N3=None, FP=None, FQ=None, FE=None, FD=None, OP=None, OQ=None):
        self.A, self.B, self.Z, self.O = A, B, Z, O
        if S is None:
            S, M, N = mu_default(A, B)
        self.S, self.M, self.N = S, M, N
        n = np.shape(A)
        # anchor: u = (AP/AQ) (1 + delta); AQ = 0: none; DS = 0 with AQ > 0: exact rational (delta = 0);
        # otherwise DS = sign(delta), DM <= ln|delta| <= DN
        self.AP = AP if AP is not None else np.zeros(n, np.int64)
        self.AQ = AQ if AQ is not None else np.zeros(n, np.int64)
        self.DS = DS if DS is not None else np.zeros(n, np.int8)
        self.DM = DM if DM is not None else np.full(n, -INF)
        self.DN = DN if DN is not None else np.full(n, -INF)
        # third level, for values 1 +- t with t below e^-(1e308): T3 = 1 when ln|ln u| = -e^nu with
        # M3 <= nu <= N3 (ln|ln u| itself is then far below the double range)
        self.T3 = T3 if T3 is not None else np.zeros(n, np.int8)
        self.M3 = M3 if M3 is not None else np.full(n, -INF)
        self.N3 = N3 if N3 is not None else np.full(n, INF)
        # (03/11) exact power of a rational: FE != 0 means u = (FP/FQ)^(FE/FD) exactly, with FP/FQ > 1 a reduced
        # rational that is not a perfect power (primitive) and FE/FD a reduced fraction, FD >= 1; two such values are
        # equal iff (FP, FQ, FE, FD) agree
        self.FP = FP if FP is not None else np.zeros(n, np.int64)
        self.FQ = FQ if FQ is not None else np.zeros(n, np.int64)
        self.FE = FE if FE is not None else np.zeros(n, np.int64)
        self.FD = FD if FD is not None else np.zeros(n, np.int64)
        # (03/11) exact offset of a triple: with FE != 0 and OQ > 0, u = (FP/FQ)^(FE/FD) + OP/OQ exactly (OP/OQ reduced,
        # either sign); OQ = 0: no offset, u is the triple itself
        self.OP = OP if OP is not None else np.zeros(n, np.int64)
        self.OQ = OQ if OQ is not None else np.zeros(n, np.int64)

    def take(self, idx):
        return Q7(*[getattr(self, k)[idx] for k in Q7.FIELDS])


def mu_default(A, B):
    A = np.asarray(A, dtype=np.float64); B = np.asarray(B, dtype=np.float64)
    S = np.zeros(A.shape, np.int8); M = np.full(A.shape, -INF); N = np.full(A.shape, INF)
    p = (A > 0)
    if p.any():
        S[p] = 1; M[p] = ln_lo(A[p]); N[p] = ln_hi(B[p])
    q = (B < 0)
    if q.any():
        S[q] = -1; M[q] = ln_lo(-B[q]); N[q] = ln_hi(-A[q])
    return S, M, N


def tighten(A, B, S, M, N):
    """intersect [A, B] with S*[e^M, e^N]"""
    p = S > 0
    if p.any():
        A[p] = np.maximum(A[p], exp_lo(M[p])); B[p] = np.minimum(B[p], exp_hi(N[p]))
    q = S < 0
    if q.any():
        A[q] = np.maximum(A[q], -exp_hi(N[q])); B[q] = np.minimum(B[q], -exp_lo(M[q]))
    return A, B


def fwd(f, x, y, R=None):
    """enclosure (Q7) of form f over the batches x (s) and y (l), elementwise; returns (Q7, ok)"""
    A1, B1, Z1, O1, S1, M1, N1 = x.A, x.B, x.Z, x.O, x.S, x.M, x.N
    A2, B2, Z2, O2, S2, M2, N2 = y.A, y.B, y.Z, y.O, y.S, y.M, y.N
    n = A1.shape
    Z = np.zeros(n, bool); O = np.zeros(n, bool); ok = np.ones(n, bool)
    S = None; M = None; N = None; pw_T3 = None
    with np.errstate(all="ignore"):
        if f == 2 or f in (3, 4):                         # s*l, s/l, l/s
            if f == 4:
                A1, B1, Z1, O1, S1, M1, N1, A2, B2, Z2, O2, S2, M2, N2 = A2, B2, Z2, O2, S2, M2, N2, A1, B1, Z1, O1, S1, M1, N1
            sg2 = 1 if f == 2 else -1                     # ln(result) = ln s + sg2 * ln l
            if f == 2:
                A = addlo(A1, A2); B = addhi(B1, B2)
                Z = Z1 | Z2
            else:
                A = addlo(A1, -B2); B = addhi(B1, -A2)
                ok = ~Z2; Z = Z1 & ~Z2
            O = O1 & O2
            S, M, N = mu_default(A, B)
            # exact 1 times a value near 1: same second channel
            c1 = O1 & ~O2 & (S2 != 0)
            S[c1] = sg2 * S2[c1]; M[c1] = M2[c1]; N[c1] = N2[c1]
            c2 = O2 & ~O1 & (S1 != 0)
            S[c2] = S1[c2]; M[c2] = M1[c2]; N[c2] = N1[c2]
            # two values near 1 whose logs have the same sign in the result: |ln| adds
            ss = ~O1 & ~O2 & (S1 != 0) & (S2 != 0) & (S1 == sg2 * S2)
            if ss.any():
                S[ss] = S1[ss]; M[ss] = lse_lo(M1[ss], M2[ss]); N[ss] = lse_hi(N1[ss], N2[ss])
            # opposite contributions with the second level of one above the other's: |ln u| = |a| - |b| with
            # mu_b <= N2 < M1 <= mu_a, so ln|ln u| = mu_a + ln(1 - e^(mu_b - mu_a)) >= M1 + ln(1 - e^(N2 - M1)) (03/11:
            # strict separation; it was N2 < M1 - 1, which the identity does not need)
            c1 = S1.astype(np.int16); c2 = (sg2 * S2).astype(np.int16)
            opp = ~O1 & ~O2 & (c1 != 0) & (c2 != 0) & (c1 == -c2)
            d1 = opp & (N2 < M1)
            if d1.any():
                S[d1] = c1[d1]; M[d1] = addlo(M1[d1], log1m_lo(exp_hi(addhi(N2[d1], -M1[d1])))); N[d1] = N1[d1]
            d2 = opp & (N1 < M2)
            if d2.any():
                S[d2] = c2[d2]; M[d2] = addlo(M2[d2], log1m_lo(exp_hi(addhi(N1[d2], -M2[d2])))); N[d2] = N2[d2]
        elif f == 0:                                      # s+l
            Mx = np.maximum(A1, A2); d1 = addlo(np.minimum(A1, A2), -Mx)
            A = addlo(Mx, log1p_lo(exp_lo(np.minimum(d1, 0.0))))
            My = np.maximum(B1, B2); d2 = addhi(np.minimum(B1, B2), -My)
            B = addhi(My, log1p_hi(exp_hi(d2)))
            A = np.where(Z1, A2, np.where(Z2, A1, A)); B = np.where(Z1, B2, np.where(Z2, B1, B))
            Z = Z1 & Z2
            O = (Z1 & O2) | (Z2 & O1)
            S, M, N = mu_default(A, B)
            # dominant operand with ln > 0 (value > 1): ln u = a + log1p(e^(b-a)) in [a, a + 1]... so
            # ln ln u in [ln a, ln a + 1/a]  (a = ln of the larger, b = ln of the smaller)
            for (Sa, Ma, Na, Ba_, Ab_, Sb, Nb, Aa, Bb_) in ((S1, M1, N1, B2, A1, S2, N2, A1, B1), (S2, M2, N2, B1, A2, S1, N1, A2, B2)):
                # the other is smaller: by ln (Ba_ < Ab_) or by ln ln (Nb < Ma, both > 1: b <= e^Nb < e^Ma <= a;
                # 03/11: it was Nb < Ma - 1, more than the bound ln u in [a, a + ln 2] needs)
                c = (Sa > 0) & (Ma > 2.0) & ((Ba_ < Ab_) | ((Sb > 0) & (Nb < Ma)))
                if c.any():
                    S[c] = 1; M[c] = np.maximum(M[c], Ma[c]); N[c] = np.minimum(N[c], addhi(Na[c], exp_hi(-Ma[c])))
                    A[c] = np.maximum(A[c], Aa[c]); B[c] = np.minimum(B[c], addhi(Bb_[c], 1.0))
            # 1 + t with t small: ln(1+t) in [t - t^2/2, t], so ln ln(1+t) in [ln t - t, ln t]
            T3 = np.zeros(n, np.int8); M3 = np.full(n, -INF); N3 = np.full(n, INF)
            for (Oa, Ab, Bb, Sb, Mb, Nb) in ((O1, A2, B2, S2, M2, N2), (O2, A1, B1, S1, M1, N1)):
                c = Oa & (Bb < -1.0)
                if c.any():
                    S[c] = 1; N[c] = Bb[c]; M[c] = addlo(Ab[c], -exp_hi(Bb[c]))
                    # t below the double range with ln t = -e^mu_t known: ln ln(1+t) = ln t + O(t) = -e^nu, nu = mu_t + O(...)
                    d = c & (Bb <= -1e308) & (Sb < 0) & np.isfinite(Mb) & np.isfinite(Nb)
                    if d.any():
                        T3[d] = 1; M3[d] = pred(Mb[d]); N3[d] = succ(Nb[d])
        elif f == 1:                                      # |s-l|
            Ba = addhi(B1, log1m_hi(exp_lo(np.minimum(addlo(A2, -B1), 0.0))))
            Bb = addhi(B2, log1m_hi(exp_lo(np.minimum(addlo(A1, -B2), 0.0))))
            A = np.full(n, -INF); B = np.maximum(Ba, Bb)
            g1 = B2 < A1
            if g1.any():
                A[g1] = addlo(A1[g1], log1m_lo(exp_hi(addhi(B2[g1], -A1[g1]))))
                B[g1] = addhi(B1[g1], log1m_hi(exp_lo(addlo(A2[g1], -B1[g1]))))
            g2 = B1 < A2
            if g2.any():
                A[g2] = addlo(A2[g2], log1m_lo(exp_hi(addhi(B1[g2], -A2[g2]))))
                B[g2] = addhi(B2[g2], log1m_hi(exp_lo(addlo(A1[g2], -B2[g2]))))
            Zr = np.zeros(n, bool)
            if R is not None:
                P1, Q1, P2, Q2 = R
                both = (Q1 > 0) & (Q2 > 0)
                eq = both & (P1 == P2) & (Q1 == Q2)
                ne_ = both & ~eq
                if ne_.any():
                    qq = succ(Q1[ne_].astype(np.float64) * Q2[ne_].astype(np.float64))
                    A[ne_] = np.maximum(A[ne_], -ln_hi(qq))
                Zr = eq
            # |e^a - e^b| = e^b |e^(a-b) - 1|; with |d| <= 1/2: ln|e^d - 1| in [ln|d| - 2|d|, ln|d| + 2|d|]
            def near(Ab, Bb, Sd, Md, Nd, base_lo, base_hi, m):
                # m: entries where d = a - b has known sign, ln|d| in [Md, Nd], |d| <= 1/2; result
                # ln = base + ln|e^d - 1| with base in [base_lo, base_hi]
                ed = 2 * exp_hi(Nd[m])
                A[m] = np.maximum(A[m], addlo(addlo(base_lo[m], Md[m]), -ed))
                B[m] = np.minimum(B[m], addhi(addhi(base_hi[m], Nd[m]), ed))
            # LH = fl(ln 1/2) exceeds ln 1/2 by 2.3e-17; the bound used by near(), |ln((e^d - 1)/d)| <= 2|d|,
            # holds for |d| <= 2, so N < LH is only a selection with ample slack (no rigour issue)
            LH = np.log(0.5)
            # s = 1 exactly, l near 1 (d = -ln l)
            m = O1 & ~O2 & (S2 != 0) & (N2 < LH)
            if m.any():
                near(None, None, S2, M2, N2, np.minimum(A2, 0.0), np.maximum(B2, 0.0), m)
            m = O2 & ~O1 & (S1 != 0) & (N1 < LH)
            if m.any():
                near(None, None, S1, M1, N1, np.minimum(A1, 0.0), np.maximum(B1, 0.0), m)
            # both near 1 with logs of opposite signs: |d| = |a| + |b|
            m = ~O1 & ~O2 & (S1 != 0) & (S2 != 0) & (S1 == -S2)
            if m.any():
                Md = np.full(n, -INF); Nd = np.full(n, INF)
                Md[m] = lse_lo(M1[m], M2[m]); Nd[m] = lse_hi(N1[m], N2[m])
                mm = m & (Nd < LH)
                if mm.any():
                    near(None, None, None, Md, Nd, np.minimum(A1, A2), np.maximum(B1, B2), mm)
            A = np.where(Z1, A2, np.where(Z2, A1, A)); B = np.where(Z1, B2, np.where(Z2, B1, B))
            Z = (Z1 & Z2) | (O1 & O2) | Zr
            O = (Z1 & O2) | (Z2 & O1)
            # 1 - t with t small (0 < t <= e^-1): -ln(1-t) in [t, t + t^2], so ln|ln| in [ln t, ln t + t]
            S, M, N = mu_default(A, B)
            T3 = np.zeros(n, np.int8); M3 = np.full(n, -INF); N3 = np.full(n, INF)
            for (Oa, Zb, Bb, Sb, Mb, Nb) in ((O1, Z2, B2, S2, M2, N2), (O2, Z1, B1, S1, M1, N1)):
                d = Oa & ~Zb & (Bb <= -1e308) & (Sb < 0) & np.isfinite(Mb) & np.isfinite(Nb)
                if d.any():
                    T3[d] = 1; M3[d] = pred(Mb[d]); N3[d] = succ(Nb[d])
            # dominant operand with ln a >= e^2 and the other at most e^(a-1): ln u in [a - 1, a], so
            # ln ln u in [ln(a - 1), ln a] subset [ln a - 2/a, ln a]
            for (Sa, Ma, Na, Bo, Aa, Sb, Nb, Ba2) in ((S1, M1, N1, B2, A1, S2, N2, B1), (S2, M2, N2, B1, A2, S1, N1, B2)):
                # the other is at most e^(a-1): by ln, or by ln ln: Nb < Ma and a lower bound >= 0 of
                # ln(e^Ma - e^Nb) give a - b >= e^Ma - e^Nb >= 1 (03/11; it was Nb < Ma - 1, i.e. b <= a/e <= a - 1)
                cl = (Sb > 0) & (Nb < Ma) & (Ma > 2.0)
                Lg = np.full(n, -INF)
                if cl.any():
                    Lg[cl] = addlo(Ma[cl], log1m_lo(exp_hi(addhi(Nb[cl], -Ma[cl]))))
                c = (Sa > 0) & (Ma > 2.0) & ((Bo < addlo(Aa, -1.0)) | (cl & (Lg >= 0.0)))
                if c.any():
                    S[c] = 1; M[c] = np.maximum(M[c], addlo(Ma[c], -2 * exp_hi(-Ma[c]))); N[c] = np.minimum(N[c], Na[c])
                    A[c] = np.maximum(A[c], addlo(Aa[c], -1.0)); B[c] = np.minimum(B[c], Ba2[c])
            for (Oa, Zb, Ab, Bb) in ((O1, Z2, A2, B2), (O2, Z1, A1, B1)):
                c = Oa & ~Zb & (Bb < -1.0)
                if c.any():
                    S[c] = -1; M[c] = Ab[c]; N[c] = addhi(Bb[c], exp_hi(Bb[c]))
                    A[c] = np.maximum(A[c], -exp_hi(N[c])); B[c] = np.minimum(B[c], -exp_lo(M[c]))
        else:                                             # powers
            if f in (7, 8):
                A1, B1, Z1, O1, S1, M1, N1, A2, B2, Z2, O2, S2, M2, N2 = A2, B2, Z2, O2, S2, M2, N2, A1, B1, Z1, O1, S1, M1, N1
            neg = f in (6, 8)
            y1 = exp_lo(A2); y2 = exp_hi(B2)
            c = [y1 * A1, y1 * B1, y2 * A1, y2 * B1]
            if neg:
                c = [-v for v in c]
            lo = fixlo(np.minimum.reduce([fixlo(v) for v in c]))
            hi = fixhi(np.maximum.reduce([fixhi(v) for v in c]))
            bad = np.isnan(c[0]) | np.isnan(c[1]) | np.isnan(c[2]) | np.isnan(c[3])
            A = np.where(bad, -INF, pred(lo)); B = np.where(bad, INF, succ(hi))
            # second channel: ln|ln result| = ln y + ln|ln x|, sign = +-sign(ln x)
            S = np.where(S1 != 0, (-S1 if neg else S1), 0).astype(np.int8)
            M = addlo(A2, M1); N = addhi(B2, N1)
            S = np.where(np.isnan(M) | np.isnan(N), 0, S).astype(np.int8)
            # base 1 +- t with ln|ln x| = -e^nu (third level) and exponent y = e^(e^mu_y) (y's ln above the
            # double range, mu_y = M2..N2 of the exponent): ln|ln result| = e^mu_y - e^nu
            T3b = x.T3 if f in (5, 6) else y.T3
            M3b = x.M3 if f in (5, 6) else y.M3
            N3b = x.N3 if f in (5, 6) else y.N3
            Sy = S2; My = M2; Ny = N2                        # exponent's second channel (after the swap)
            dd = (T3b > 0) & (S1 != 0) & (Sy > 0) & np.isfinite(My) & np.isfinite(Ny) & (A2 > 700)
            if dd.any():
                big = dd & (My > addhi(N3b, 1.0))            # e^mu_y - e^nu >= e^mu_y (1 - 1/e) > 0
                small = dd & (Ny < addlo(M3b, -1.0))         # e^mu_y - e^nu <= -e^nu (1 - 1/e) < 0
                if big.any():
                    # ln|ln result| = Lam in [e^My (1 - e^(N3-My)), e^Ny]
                    Llo = pred(exp_lo(My) * pred(1.0 - exp_hi(addhi(N3b, -My))))
                    Lhi = exp_hi(Ny)
                    alo = exp_lo(Llo); ahi = exp_hi(Lhi)
                    A = np.where(big, np.where(S > 0, alo, -ahi), A); B = np.where(big, np.where(S > 0, ahi, -alo), B)
                    M = np.where(big, Llo, M); N = np.where(big, Lhi, N)
                if small.any():
                    # ln|ln result| <= -e^M3 (1 - e^(Ny-M3)) =: -Lam
                    Lam = pred(exp_lo(M3b) * pred(1.0 - exp_hi(addhi(Ny, -M3b))))
                    ahi = exp_hi(-Lam)
                    A = np.where(small, np.where(S > 0, 0.0, -ahi), A); B = np.where(small, np.where(S > 0, ahi, 0.0), B)
                    M = np.where(small, -INF, M); N = np.where(small, -Lam, N)
                amb = dd & ~big & ~small
                S = np.where(amb, 0, S).astype(np.int8)
            # exponent below e^-(1e308) (ln y = -e^mu_y, mu_y known) and base with ln|ln x| = mu_x known:
            # ln|ln result| = ln y + mu_x = -(e^mu_y - mu_x) = -e^nu,  nu = ln(e^mu_y - mu_x), defined when
            # e^mu_y > mu_x (ensured by mu_y > ln max(|mu_x|, 1) + 1).  Bounds (mu_y in [M2, N2], mu_x in [M1, N1]):
            #   nu >= ln(e^M2 - N1^+) = M2 + ln(1 - N1^+ e^-M2),   nu <= ln(e^N2 + M1^-) = N2 + ln(1 + M1^- e^-N2)
            # and then ln|ln result| = -e^nu lies in [-e^N3, -e^M3] (finite or not).
            # (Corrected 02/11: the earlier bounds nu >= pred(M2), ln|ln u| <= -DBL_MAX fail for mu_y near 709.2
            # or |mu_x| > 1e295.)
            dt = (S2 < 0) & (B2 <= -1e308) & np.isfinite(M2) & np.isfinite(N2) & (S1 != 0) & np.isfinite(M1) & np.isfinite(N1)
            if dt.any():
                big_mu = np.maximum(np.abs(M1), np.abs(N1))
                ok3 = dt & (M2 > addhi(ln_hi(np.maximum(big_mu, 1.0)), 1.0))
                if ok3.any():
                    T3p = np.zeros(n, np.int8); M3p = np.full(n, -INF); N3p = np.full(n, INF)
                    tl = succ(np.maximum(N1[ok3], 0.0) * exp_hi(-M2[ok3]))      # >= N1^+ e^-M2  (< e^-1)
                    th = succ(np.maximum(-M1[ok3], 0.0) * exp_hi(-N2[ok3]))     # >= M1^- e^-N2  (< e^-1)
                    T3p[ok3] = 1
                    M3p[ok3] = addlo(M2[ok3], log1m_lo(tl)); N3p[ok3] = addhi(N2[ok3], log1p_hi(th))
                    M = np.where(ok3, -exp_hi(N3p), M); N = np.where(ok3, -exp_lo(M3p), N)
                    # instrumentation: would the bounds of the code before 02/11 have been invalid here?
                    bump("T3exp_fired", ok3.sum())
                    bump("T3exp_old_N_invalid", (N[ok3] > -np.finfo(float).max).sum())        # old: N = -DBL_MAX
                    bump("T3exp_old_nu_lo_invalid", (pred(M2[ok3] - 1e-300) > M3p[ok3]).sum())  # old: M3 = pred(M2)
                    pw_T3 = (T3p, M3p, N3p)
                else:
                    pw_T3 = None
            else:
                pw_T3 = None
            A, B = tighten(A, B, S, M, N)
            one = O1 | (Z2 & ~Z1)
            O = one.copy()
            zb = Z1 & ~Z2
            if neg:
                ok = ~zb & ~(Z1 & Z2)
            else:
                ok = ~(Z1 & Z2)
                Z = zb.copy(); O = O & ~zb
    A = np.where(Z, -INF, A); B = np.where(Z, -INF, B)
    A = np.where(O, 0.0, A); B = np.where(O, 0.0, B)
    if S is None:
        S, M, N = mu_default(A, B)
    r = Q7(A, B, Z, O, S, M, N)
    if f in (0, 1):
        r.T3, r.M3, r.N3 = T3, M3, N3
    elif f >= 5 and pw_T3 is not None:
        r.T3, r.M3, r.N3 = pw_T3
    anchors(f, x, y, r)
    if FACTOR:
        factored(f, x, y, r)
    r.S = np.where(r.Z | r.O, 0, r.S).astype(np.int8)
    # (03/11) a triple (FP/FQ)^FE with an integer exponent whose value fits is also that exact rational
    tr = (r.FE != 0) & (r.FD <= 1) & (r.OQ == 0) & ~r.Z & ~r.O & ((r.AQ == 0) | (r.DS != 0))
    if tr.any():
        rows = np.where(tr)[0]
        e_ = r.FE[rows].astype(np.int64); a_ = r.FP[rows].astype(np.int64); b_ = r.FQ[rows].astype(np.int64)
        k_ = np.abs(e_)
        pa, oka = ipow(a_, k_); pb, okb = ipow(b_, k_)
        okt = oka & okb
        if okt.any():
            rr = rows[okt]
            r.AP[rr] = np.where(e_[okt] > 0, pa[okt], pb[okt]); r.AQ[rr] = np.where(e_[okt] > 0, pb[okt], pa[okt])
            r.DS[rr] = 0; r.DM[rr] = -INF; r.DN[rr] = -INF
            bump("triple_rational", int(okt.sum()))
    # (03/11) an exact 1 is also the exact rational 1/1 (so that the exact differences of C.4.4(3) see it)
    one = r.O & ~r.Z & ((r.AQ == 0) | (r.DS != 0))
    if one.any():
        r.AP[one] = 1; r.AQ[one] = 1; r.DS[one] = 0; r.DM[one] = -INF; r.DN[one] = -INF
    r.A, r.B = tighten(r.A, r.B, r.S, r.M, r.N)
    return r, ok


FACTOR = os.environ.get("IVL_FACTOR", "1") == "1"
_PRIM = {}
_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61)


def _iroot(a, p):
    """the integer p-th root of a >= 1 if a is a perfect p-th power, else None (exact integer arithmetic)"""
    r = int(round(a ** (1.0 / p)))
    for c in (r - 1, r, r + 1):
        if c >= 1 and c ** p == a:
            return c
    return None


def primitive(P, Q):
    """(r_P, r_Q, m) with P/Q = (r_P/r_Q)^m and r_P/r_Q not a perfect power (P, Q > 0 reduced integers);
    exact: every root is checked by exponentiation"""
    key = (int(P), int(Q))
    if key not in _PRIM:
        a, b, m = key[0], key[1], 1
        changed = True
        while changed and not (a == 1 and b == 1):
            changed = False
            for p in _PRIMES:
                ra, rb = _iroot(a, p), _iroot(b, p)
                if ra is not None and rb is not None:
                    a, b, m = ra, rb, m * p; changed = True
                    break
        _PRIM[key] = (a, b, m)
    return _PRIM[key]


def _triple(z, i):
    """the triple (r_P, r_Q, e_num, e_den) of the value i of the batch z: its factored channel, or the primitive form of
    an exact rational (P/Q = r^m); None if neither (or the value is 1)"""
    if z.FE[i] != 0:
        if z.OQ[i] > 0:
            return None                                    # X + c: not a power
        return int(z.FP[i]), int(z.FQ[i]), int(z.FE[i]), max(int(z.FD[i]), 1)
    if z.AQ[i] > 0 and z.DS[i] == 0 and not z.Z[i] and not z.O[i] and z.AP[i] > 0:
        a, b, m = primitive(z.AP[i], z.AQ[i])
        if a == b:
            return None
        if a < b:
            a, b, m = b, a, -m
        return a, b, m, 1
    return None


def _set_triple(r, i, a, b, en, ed):
    """store u = (a/b)^(en/ed) in r[i] in canonical form (a/b > 1, en/ed reduced, ed >= 1); False if too large"""
    g = gcd(en, ed); en //= g; ed //= g
    if ed < 0:
        en, ed = -en, -ed
    if a < b:                                              # canonical orientation r > 1: r^e = (1/r)^(-e)
        a, b, en = b, a, -en
    if en == 0 or abs(en) >= 2 ** 62 or ed >= 2 ** 62:
        return False
    r.FP[i], r.FQ[i], r.FE[i], r.FD[i] = a, b, en, ed
    return True


OFFSET = os.environ.get("IVL_OFFSET", "1") == "1"


def _offset(z, i):
    """the exact offset c of z[i] = X + c (0 without one)"""
    return Fr(int(z.OP[i]), int(z.OQ[i])) if z.OQ[i] > 0 else Fr(0)


def _set_offset(r, t, i, c):
    """r[i] = X + c, X the triple of t[i] (the value itself when c = 0)"""
    if abs(c.numerator) >= 2 ** 62 or c.denominator >= 2 ** 62:
        return
    r.FP[i], r.FQ[i], r.FE[i], r.FD[i] = t.FP[i], t.FQ[i], t.FE[i], max(int(t.FD[i]), 1)
    if c == 0:
        r.OP[i], r.OQ[i] = 0, 0
    else:
        r.OP[i], r.OQ[i] = c.numerator, c.denominator
    bump("offset_set")


def _set_rational(r, i, d):
    """r[i] is exactly the positive rational d: exact anchor and its enclosure"""
    r.AP[i], r.AQ[i] = d.numerator, d.denominator
    r.DS[i] = 0; r.DM[i] = -INF; r.DN[i] = -INF
    lo = float(ln_lo(np.array([np.nextafter(d.numerator / d.denominator, 0)]))[0])
    hi = float(ln_hi(np.array([np.nextafter(d.numerator / d.denominator, INF)]))[0])
    r.A[i], r.B[i] = lo, hi
    r.Z[i] = False; r.O[i] = (d == 1)
    if d == 1:
        r.A[i] = r.B[i] = 0.0
    r.FE[i] = 0; r.OQ[i] = 0; r.T3[i] = 0


def factored(f, x, y, r):
    """the channel (FP, FQ, FE, FD) of r = f(x, y): u = (FP/FQ)^(FE/FD) exactly, with FP/FQ > 1 primitive and FE/FD a
    reduced fraction, so the triple of a value is unique: if r1^e1 = r2^e2 with r1, r2 > 1 primitive and e1, e2
    rational, multiplying the exponents by a common denominator D gives r1^(D e1) = r2^(D e2) with integer exponents,
    and unique factorisation with the primitivity of r1 and r2 forces r1 = r2 and e1 = e2.
      powers:    an exact rational q = r^m, or a triple r^e, raised to an exact positive rational +-p/q:
                 r^(+-m p/q), r^(+-e p/q), since (r^e)^t = r^(e t) for r > 0 and real t;
      products:  r^e1 r^e2 = r^(e1+e2);  quotients r^e1 / r^e2 = r^(e1-e2); exponent 0 means exactly 1;
      |s - l|:   the same triple gives exactly 0; different triples are different values, and the gap
                 ln|x - y| = max + ln(1 - e^-|L|), L = ln y - ln x != 0, is bounded in Arb (_factor_gap).
                 An exact rational on one side is compared through its own primitive triple.
    Everything is exact integer arithmetic; a triple is kept only while |FE|, FD < 2^62."""
    n = r.A.shape[0]
    fx = x.FE != 0; fy = y.FE != 0
    LIM = 2.0 ** 61                                        # int64 products are formed only below this (float check)
    if f >= 5:
        base, expo = (x, y) if f in (5, 6) else (y, x)
        sg = -1 if f in (6, 8) else 1
        # exponent: an exact positive rational p/q (q = 1: an integer), not itself a triple
        ek = (expo.AQ >= 1) & (expo.DS == 0) & (expo.AP >= 1) & ~expo.Z & (expo.FE == 0)
        cand = ek & (((base.FE != 0) & (base.OQ == 0)) | ((base.AQ > 0) & (base.DS == 0) & ~base.Z & ~base.O & (base.AP > 0))) & (r.AQ == 0)
        rows = np.where(cand)[0]
        if not len(rows):
            return
        ok, a, b, en, ed = _tri_vec(base, rows)
        p_ = expo.AP[rows].astype(np.int64); q_ = expo.AQ[rows].astype(np.int64)
        safe = ok & (np.abs(en).astype(np.float64) * p_.astype(np.float64) < LIM) & (ed.astype(np.float64) * q_.astype(np.float64) < LIM)
        num = np.where(safe, en * p_, 0) * sg; den = np.where(safe, ed * q_, 1)
        g = np.gcd(num, den); g = np.where(g == 0, 1, g)
        num //= g; den //= g
        st = safe & (num != 0)
        if st.any():
            rr = rows[st]
            r.FP[rr] = a[st]; r.FQ[rr] = b[st]; r.FE[rr] = num[st]; r.FD[rr] = den[st]
            bump("factor_pow", int(st.sum()))
        return
    if f == 0 and OFFSET:
        # X + c plus an exact rational q: X + (c + q), exactly; kept for huge X only (ln X > 40), where the double
        # enclosures cannot tell X + c1 from X + c2
        tA = np.where(fx, x.A, y.A)
        cand = (fx ^ fy) & ~r.Z & (r.FE == 0) & (tA > 40.0)
        rows = np.where(cand)[0]
        if not len(rows):
            return
        sx = fx[rows]
        g_ = lambda fld: np.where(sx, getattr(x, fld)[rows], getattr(y, fld)[rows])      # the X + c side
        h_ = lambda fld: np.where(sx, getattr(y, fld)[rows], getattr(x, fld)[rows])      # the other side
        OPt, OQt = g_("OP"), g_("OQ")
        oZ, oAP, oAQ, oDS = h_("Z"), h_("AP"), h_("AQ"), h_("DS")
        valid = oZ | ((oAQ > 0) & (oDS == 0))
        qn = np.where(oZ, 0, oAP).astype(np.int64); qd = np.where(oZ, 1, np.maximum(oAQ, 1)).astype(np.int64)
        cn = np.where(OQt > 0, OPt, 0).astype(np.int64); cd = np.where(OQt > 0, OQt, 1).astype(np.int64)
        safe = valid & (np.abs(cn).astype(np.float64) * qd + np.abs(qn).astype(np.float64) * cd < LIM) & (cd.astype(np.float64) * qd < LIM)
        num = np.where(safe, cn * qd + qn * cd, 0); den = np.where(safe, cd * qd, 1)
        gg = np.gcd(num, den); gg = np.where(gg == 0, 1, gg); num //= gg; den //= gg
        if safe.any():
            rr = rows[safe]; sf = safe
            r.FP[rr] = g_("FP")[sf]; r.FQ[rr] = g_("FQ")[sf]; r.FE[rr] = g_("FE")[sf]; r.FD[rr] = np.maximum(g_("FD")[sf], 1)
            r.OP[rr] = np.where(num[sf] != 0, num[sf], 0); r.OQ[rr] = np.where(num[sf] != 0, den[sf], 0)
            bump("offset_set", int(safe.sum()))
        return
    if f == 1 and OFFSET:
        # (a) X + c1 and X + c2 with the same X: |c1 - c2| exactly; (b) |q - (X + c)| = X + (c - q) when X + c > q;
        # for huge X only (ln X > 40), as for the sums
        big = ((x.A > 40.0) & fx) | ((y.A > 40.0) & fy)
        cand = ~r.Z & (fx | fy) & ((x.OQ > 0) | (y.OQ > 0) | (fx ^ fy)) & big
        for i in np.where(cand)[0]:
            if fx[i] and fy[i]:
                if (x.FP[i], x.FQ[i], x.FE[i], max(int(x.FD[i]), 1)) != (y.FP[i], y.FQ[i], y.FE[i], max(int(y.FD[i]), 1)):
                    continue
                d = abs(_offset(x, i) - _offset(y, i))
                if d == 0:
                    r.Z[i] = True; r.O[i] = False; r.A[i] = -INF; r.B[i] = -INF
                    r.AQ[i] = 0; r.DS[i] = 0; r.T3[i] = 0
                    bump("offset_zero")
                elif d.numerator < 2 ** 62 and d.denominator < 2 ** 62:
                    _set_rational(r, i, d)
                    bump("offset_diff")
            else:
                t, o = (x, y) if fx[i] else (y, x)
                if r.FE[i] != 0 or not t.A[i] > 40.0:
                    continue
                if o.Z[i]:
                    q = Fr(0)                              # |0 - (X + c)| = X + c
                elif o.AQ[i] > 0 and o.DS[i] == 0:
                    q = Fr(int(o.AP[i]), int(o.AQ[i]))
                else:
                    continue
                # X + c > q certainly: the enclosure of X + c is above ln q (A_t > ln_hi(q))
                if q == 0 or t.A[i] > float(ln_hi(np.array([float(q) * (1 + 2.0 ** -40)]))[0]):
                    _set_offset(r, t, i, _offset(t, i) - q)
    if f == 1:
        # rows with no lower bound where both sides are exact powers of rationals (a triple, or an exact rational)
        cand = ~r.Z & (r.A == -INF) & (fx | fy)
        rows = np.where(cand)[0]
        if not len(rows):
            return
        ok1, a1, b1, n1, d1 = _tri_vec(x, rows); ok2, a2, b2, n2, d2 = _tri_vec(y, rows)
        both = ok1 & ok2
        eq = both & (a1 == a2) & (b1 == b2) & (n1 == n2) & (d1 == d2)
        if eq.any():                                       # equal exact values: |s - l| = 0
            rr = rows[eq]
            r.Z[rr] = True; r.O[rr] = False; r.A[rr] = -INF; r.B[rr] = -INF
            r.AQ[rr] = 0; r.DS[rr] = 0; r.T3[rr] = 0
            bump("factor_zero", int(eq.sum()))
        for j in np.where(both & ~eq)[0]:
            _factor_gap((int(a1[j]), int(b1[j]), int(n1[j]), int(d1[j])), (int(a2[j]), int(b2[j]), int(n2[j]), int(d2[j])), r, rows[j])
        return
    if f in (2, 3, 4):
        # a triple times (or over) a triple or an exact rational of the same primitive base
        cand = (fx | fy) & ~r.Z & ~r.O & (r.FE == 0)
        rows = np.where(cand)[0]
        if not len(rows):
            return
        ok1, a1, b1, n1, d1 = _tri_vec(x, rows); ok2, a2, b2, n2, d2 = _tri_vec(y, rows)
        same = ok1 & ok2 & (a1 == a2) & (b1 == b2)
        if not same.any():
            return
        fa = np.abs(n1).astype(np.float64); fb = np.abs(n2).astype(np.float64)
        safe = same & (fa * d2 + fb * d1 < LIM) & (d1.astype(np.float64) * d2 < LIM)
        if f == 2:
            num = n1 * d2 + n2 * d1
        elif f == 3:
            num = n1 * d2 - n2 * d1
        else:
            num = n2 * d1 - n1 * d2
        num = np.where(safe, num, 0); den = np.where(safe, d1 * d2, 1)
        one = safe & (num == 0)
        if one.any():                                      # r^e / r^e = 1 exactly
            rr = rows[one]
            r.O[rr] = True; r.Z[rr] = False; r.A[rr] = 0.0; r.B[rr] = 0.0
            r.AQ[rr] = 0; r.DS[rr] = 0; r.T3[rr] = 0
            bump("factor_one", int(one.sum()))
        st = safe & (num != 0)
        if st.any():
            gg = np.gcd(num, den); gg = np.where(gg == 0, 1, gg); num //= gg; den //= gg
            rr = rows[st]
            r.FP[rr] = a1[st]; r.FQ[rr] = b1[st]; r.FE[rr] = num[st]; r.FD[rr] = den[st]


def _tri_vec(z, rows):
    """vectorised _triple for the given rows: (ok, a, b, en, ed) arrays with a/b > 1 primitive and en/ed reduced, ed >= 1:
    the factored channel of a pure power, or the primitive form of an exact rational (primitive() once per distinct
    rational)"""
    k = len(rows)
    ok = np.zeros(k, bool); a = np.zeros(k, np.int64); b = np.zeros(k, np.int64)
    en = np.zeros(k, np.int64); ed = np.ones(k, np.int64)
    FE = z.FE[rows]; OQ = z.OQ[rows]
    tm = (FE != 0) & (OQ == 0)
    if tm.any():
        a[tm] = z.FP[rows][tm]; b[tm] = z.FQ[rows][tm]; en[tm] = FE[tm]; ed[tm] = np.maximum(z.FD[rows][tm], 1); ok[tm] = True
    rm = (FE == 0) & (z.AQ[rows] > 0) & (z.DS[rows] == 0) & ~z.Z[rows] & ~z.O[rows] & (z.AP[rows] > 0)
    if rm.any():
        idx = np.where(rm)[0]
        pairs = np.stack([z.AP[rows][idx], z.AQ[rows][idx]], 1)
        uq, inv = np.unique(pairs, axis=0, return_inverse=True)
        inv = np.asarray(inv).ravel()
        pa = np.zeros(len(uq), np.int64); pb = np.zeros(len(uq), np.int64); pm = np.zeros(len(uq), np.int64); pok = np.zeros(len(uq), bool)
        for j in range(len(uq)):
            aa, bb, mm = primitive(int(uq[j, 0]), int(uq[j, 1]))
            if aa == bb:
                continue
            if aa < bb:
                aa, bb, mm = bb, aa, -mm
            pa[j], pb[j], pm[j], pok[j] = aa, bb, mm, True
        a[idx] = pa[inv]; b[idx] = pb[inv]; en[idx] = pm[inv]; ed[idx] = 1; ok[idx] = pok[inv]
    return ok, a, b, en, ed


def _factor_gap(t1, t2, r, i):
    """lower bound of ln|x - y| for different exact powers x = r1^e1, y = r2^e2 (triples t = (a, b, en, ed), e = en/ed;
    Arb balls at raised precision until L = e2 ln r2 - e1 ln r1 is separated from 0, which it is since x != y)"""
    from flint import arb, ctx
    old = ctx.prec
    try:
        for prec in (256, 1024, 4096, 16384):
            ctx.prec = prec
            l1 = arb(t1[2]) / arb(t1[3]) * (arb(t1[0]).log() - arb(t1[1]).log())
            l2 = arb(t2[2]) / arb(t2[3]) * (arb(t2[0]).log() - arb(t2[1]).log())
            L = l2 - l1
            if not L.contains(0):
                break
        else:
            return
        hi = l1 if l1 > l2 else l2                         # the larger logarithm (they are separated by L)
        lo_b = hi + (1 - (-abs(L)).exp()).log()
        a = float(lo_b.lower()); b = float(hi.upper())
        if not (np.isfinite(a) and np.isfinite(b)):
            return
        a = float(pred(np.array([a]))[0]); b = float(succ(np.array([b]))[0])
        r.A[i] = max(r.A[i], a); r.B[i] = min(r.B[i], b)
        bump("factor_gap")
    finally:
        ctx.prec = old


def possibly_out(A, B, Z):
    return ~Z & ((B > LN_HI_LO) | (A < LN_LO_HI))


def possibly_in(A, B, Z):
    return Z | ((A <= LN_HI_HI) & (B >= LN_LO_LO))


# ---------------------------------------------------------------- preimage windows (Lemma X3)
def divlo(a, b):
    with np.errstate(all="ignore"):
        return fixlo(pred(fixlo(np.asarray(a, dtype=np.float64) / b)))


def divhi(a, b):
    with np.errstate(all="ignore"):
        return fixhi(succ(fixhi(np.asarray(a, dtype=np.float64) / b)))


def ldiff_lo(x, y):
    """lower bound of ln(e^x - e^y) (x, y lower/upper bounds resp.); -inf if not x > y"""
    x = np.asarray(x, dtype=np.float64); y = np.broadcast_to(np.asarray(y, dtype=np.float64), x.shape)
    out = np.full(x.shape, -INF)
    m = (x > y) & np.isfinite(x)
    if m.any():
        out[m] = addlo(x[m], log1m_lo(exp_hi(addhi(y[m], -x[m]))))
    return out


def ldiff_hi(x, y):
    """upper bound of ln(e^x - e^y) (x upper, y lower bound); -inf if x <= y"""
    x = np.asarray(x, dtype=np.float64); y = np.broadcast_to(np.asarray(y, dtype=np.float64), x.shape)
    out = np.full(x.shape, -INF)
    m = x > y
    if m.any():
        xm, ym = x[m], y[m]
        r = np.where(np.isinf(xm), INF, 0.0)
        f = np.isfinite(xm)
        r[f] = addhi(xm[f], log1m_hi(exp_lo(np.minimum(addlo(ym[f], -xm[f]), 0.0))))
        out[m] = r
    return out


def lse_lo(x, y):
    M = np.maximum(x, y); d = addlo(np.minimum(x, y), -M)
    return addlo(M, log1p_lo(exp_lo(np.minimum(d, 0.0))))


def lse_hi(x, y):
    M = np.maximum(x, y); d = addhi(np.minimum(x, y), -M)
    return addhi(M, log1p_hi(exp_hi(d)))


def _ratio_window(al, be, A1, B1):
    """ln of {y > 0 : y*mu in [al, be] for some mu in [A1, B1]} as [lo, hi] (outward)"""
    n = A1.shape
    lo = np.full(n, INF); hi = np.full(n, -INF)           # empty by default
    pos = A1 > 0
    if pos.any():
        a, b, A, B = al[pos], be[pos], A1[pos], B1[pos]
        ok = be[pos] > 0
        ylo = np.where(a > 0, divlo(np.maximum(a, 0), B), 0.0)
        yhi = divhi(b, A)
        l = np.where(ylo > 0, ln_lo(ylo), -INF); h = ln_hi(yhi)
        lo[pos] = np.where(ok, l, INF); hi[pos] = np.where(ok, h, -INF)
    neg = B1 < 0
    if neg.any():
        a, b, A, B = al[neg], be[neg], A1[neg], B1[neg]
        ok = a < 0
        ylo = np.where(b < 0, divlo(-b, -A), 0.0)
        yhi = divhi(-a, -B)
        l = np.where(ylo > 0, ln_lo(ylo), -INF); h = ln_hi(yhi)
        lo[neg] = np.where(ok, l, INF); hi[neg] = np.where(ok, h, -INF)
    mid = ~pos & ~neg
    lo[mid] = -INF; hi[mid] = INF
    return lo, hi


def pre(f, A1, B1, Z1, O1, al, be):
    out = _pre(f, A1, B1, Z1, O1, al, be)
    A1, B1, Z1, O1, al, be = [np.asarray(x) for x in np.broadcast_arrays(A1, B1, Z1, O1, al, be)]
    # a window [-inf, -inf] admits only l = 0: keep it only where f(s, 0) is defined and may lie in the window
    s_in = (B1 >= al) & (A1 <= be)
    has0 = al == -INF; has1 = (al <= 0) & (be >= 0)
    z_ok = {0: s_in, 1: s_in, 2: has0, 3: np.zeros(al.shape, bool), 4: has0 & ~Z1, 5: has1, 6: has1,
            7: has0, 8: np.zeros(al.shape, bool)}[f]
    res = []
    for lo, hi in out:
        only0 = (hi == -INF)
        lo = np.where(only0 & ~z_ok, INF, lo)
        res.append((lo, hi))
    return res


def _pre(f, A1, B1, Z1, O1, al, be):
    """windows [lo, hi] (ln l) containing every l with f(s, l) in [al, be] (ln), for s in (A1, B1, Z1, O1).
    Returns a list of (lo, hi) pairs (several branches); lo > hi means empty.  A window with lo = -inf
    also admits l = 0 (zero entries have A = B = -inf)."""
    A1, B1, Z1, O1, al, be = [np.asarray(x) for x in np.broadcast_arrays(A1, B1, Z1, O1, al, be)]
    n = A1.shape
    full_lo, full_hi = np.full(n, -INF), np.full(n, INF)
    has0 = al == -INF                                        # the window contains the value 0
    has1 = (al <= 0) & (be >= 0)                             # the window contains the value 1
    if f == 2:
        lo, hi = addlo(al, -B1), addhi(be, -A1)
        lo = np.where(Z1, np.where(has0, -INF, INF), lo); hi = np.where(Z1, np.where(has0, INF, -INF), hi)
        return [(lo, hi)]
    if f == 3:
        lo, hi = addlo(A1, -be), addhi(B1, -al)
        lo = np.where(Z1, np.where(has0, -INF, INF), lo); hi = np.where(Z1, np.where(has0, INF, -INF), hi)
        return [(lo, hi)]
    if f == 4:
        lo, hi = addlo(al, A1), addhi(be, B1)
        return [(np.where(Z1, INF, lo), np.where(Z1, -INF, hi))]
    if f == 0:
        hi = ldiff_hi(be, A1)
        lo = ldiff_lo(al, B1)
        hi = np.where(hi == -INF, -INF, hi)
        lo = np.where(Z1, al, lo); hi = np.where(Z1, be, hi)
        # l = 0 admissible when s itself may be in the window: lo = -inf covers it
        s_in = (B1 >= al) & (A1 <= be)
        lo = np.where(s_in, -INF, lo)
        hi = np.where(s_in & (hi < lo), -INF, hi) if False else hi
        return [(lo, np.where(s_in, np.maximum(hi, -INF), hi))]
    if f == 1:
        lo1, hi1 = lse_lo(A1, al), lse_hi(B1, be)
        hi2 = np.maximum(ldiff_hi(B1, al), ldiff_hi(be, A1))
        lo2 = np.where(B1 < al, ldiff_lo(al, B1), np.where(A1 > be, ldiff_lo(A1, be), -INF))
        lo1 = np.where(Z1, al, lo1); hi1 = np.where(Z1, be, hi1)
        lo2 = np.where(Z1, INF, lo2); hi2 = np.where(Z1, -INF, hi2)
        # l = 0: |s - 0| = s in the window
        s_in = (B1 >= al) & (A1 <= be)
        lo2 = np.where(s_in, -INF, lo2); hi2 = np.where(s_in, np.maximum(hi2, -INF), hi2)
        return [(lo1, hi1), (lo2, hi2)]
    if f in (5, 6):
        a, b = (al, be) if f == 5 else (-be, -al)
        lo, hi = _ratio_window(a, b, A1, B1)
        lo = np.where(has1, -INF, lo)                        # l = 0 gives s^0 = 1
        hi = np.where(has1 & (hi < -1e308), -INF, hi)
        lo = np.where(O1, np.where(has1, -INF, INF), lo); hi = np.where(O1, np.where(has1, INF, -INF), hi)
        lo = np.where(Z1, np.where(has0 & (f == 5), -INF, INF), lo); hi = np.where(Z1, np.where(has0 & (f == 5), INF, -INF), hi)
        return [(lo, hi)]
    if f in (7, 8):
        a, b = (al, be) if f == 7 else (-be, -al)
        tlo, thi = exp_lo(A1), exp_hi(B1)
        lo = np.where(a >= 0, divlo(a, thi), divlo(a, tlo))
        hi = np.where(b <= 0, divhi(b, thi), divhi(b, tlo))
        lo = np.where(Z1, np.where(has1, -INF, INF), lo); hi = np.where(Z1, np.where(has1, INF, -INF), hi)
        lo = np.where(has0 & (f == 7), -INF, lo)              # l = 0: 0^s = 0
        return [(lo, hi)]
    raise ValueError(f)



# ---------------------------------------------------------------- anchors (shadows of exact rationals)
DTINY = -20.0          # anchors are kept only while ln|delta| <= DTINY (|delta| < 2e-9)
_I62 = float(2 ** 62)


_LNQ = {}


def lnq_bounds(P, Q):
    """bounds of ln(P/Q) for positive int64 arrays (rigorous bounds computed once per distinct P/Q)"""
    if len(P) == 0:
        return np.empty(0), np.empty(0)
    key = np.stack([P, Q], 1)
    u, inv = np.unique(key, axis=0, return_inverse=True)
    inv = inv.ravel()
    lo = np.empty(len(u)); hi = np.empty(len(u))
    todo = []
    for i, (p, q) in enumerate(map(tuple, u)):
        c = _LNQ.get((p, q))
        if c is None:
            todo.append(i)
        else:
            lo[i], hi[i] = c
    if todo:
        t = np.array(todo)
        l, h = _lnq_raw(u[t, 0], u[t, 1])
        lo[t] = l; hi[t] = h
        for i, a, b in zip(t, l, h):
            _LNQ[(int(u[i, 0]), int(u[i, 1]))] = (a, b)
    return lo[inv], hi[inv]


def _lnq_raw(P, Q):
    p = P.astype(np.float64); q = Q.astype(np.float64)
    lo = addlo(ln_lo(pred(p)), -ln_hi(succ(q))); hi = addhi(ln_hi(succ(p)), -ln_lo(pred(q)))
    exactf = (np.abs(P) < 2 ** 53) & (Q < 2 ** 53)
    lo = np.where(exactf, addlo(ln_lo(p), -ln_hi(q)), lo); hi = np.where(exactf, addhi(ln_hi(p), -ln_lo(q)), hi)
    return lo, hi


def _roots_vec(P, r):
    """exact integer r-th roots of the int64 arrays P >= 1 (0 where P is not a perfect r-th power); the roots are
    computed once per distinct (P, r) with exact integer arithmetic (_iroot_big)"""
    out = np.zeros(len(P), np.int64)
    if not len(P):
        return out
    pairs = np.stack([P, r], 1)
    uq, inv = np.unique(pairs, axis=0, return_inverse=True)
    inv = np.asarray(inv).ravel()
    vals = np.zeros(len(uq), np.int64)
    for j in range(len(uq)):
        v = _iroot_big(int(uq[j, 0]), int(uq[j, 1])) if uq[j, 0] >= 1 else None
        vals[j] = v if v is not None else 0
    return vals[inv]


def _iroot_big(a, r):
    """the integer r-th root of a >= 1 if a is a perfect r-th power, else None (exact, Python integers)"""
    if a == 1:
        return 1
    lo, hi = 1, 1 << ((a.bit_length() + r - 1) // r + 1)
    while lo < hi:                                         # largest c with c^r <= a
        mid = (lo + hi + 1) // 2
        if mid ** r <= a:
            lo = mid
        else:
            hi = mid - 1
    return lo if lo ** r == a else None


def ipow(P, k):
    """P^k for int64 arrays P >= 1, k >= 1, exactly; ok False where P^k may reach 2^61 (then P^k is not computed)"""
    P = np.asarray(P, np.int64); k = np.asarray(k, np.int64)
    with np.errstate(all="ignore"):
        ok = (P >= 1) & (k >= 1) & (k.astype(np.float64) * np.log2(np.maximum(P, 1).astype(np.float64)) < 60.0)
    e = np.where(ok, k, 0); b = np.where(ok, P, 1); r = np.ones(P.shape, np.int64)
    while (e > 0).any():
        r = np.where((e & 1) == 1, r * b, r)
        e = e >> 1
        b = np.where(e > 0, b * b, b)          # b = P^(2^j) <= P^k < 2^60 while bits of k remain
    return r, ok


def rat_op(op, P1, Q1, P2, Q2):
    """exact rational op (add, sub, mul, div) on int64 arrays; ok False where it would overflow"""
    f = lambda a: a.astype(np.float64)
    if op == "add" or op == "sub":
        big = (np.abs(f(P1) * f(Q2)) > _I62 / 4) | (np.abs(f(P2) * f(Q1)) > _I62 / 4) | (f(Q1) * f(Q2) > _I62)
        P1s, Q1s, P2s, Q2s = [np.where(big, 1, a) for a in (P1, Q1, P2, Q2)]
        num = P1s * Q2s + (P2s * Q1s if op == "add" else -P2s * Q1s); den = Q1s * Q2s
    elif op == "mul":
        big = (np.abs(f(P1) * f(P2)) > _I62) | (f(Q1) * f(Q2) > _I62)
        P1s, Q1s, P2s, Q2s = [np.where(big, 1, a) for a in (P1, Q1, P2, Q2)]
        num = P1s * P2s; den = Q1s * Q2s
    else:
        big = (np.abs(f(P1) * f(Q2)) > _I62) | (np.abs(f(Q1) * f(P2)) > _I62) | (P2 == 0)
        P1s, Q1s, P2s, Q2s = [np.where(big, 1, a) for a in (P1, Q1, P2, Q2)]
        num = P1s * Q2s; den = Q1s * P2s
    sg = np.sign(num) * np.sign(den)
    num = np.abs(num); den = np.abs(den)
    g = np.gcd(num, den); g = np.where(g == 0, 1, g)
    num //= g; den //= g
    return num, den, sg, ~big


def _lnq(P, Q, m):
    """lnq_bounds evaluated only where m (elsewhere: -inf, inf)"""
    lo = np.full(P.shape, -INF); hi = np.full(P.shape, INF)
    if m.any():
        lo[m], hi[m] = lnq_bounds(np.maximum(P[m], 1), np.maximum(Q[m], 1))
    return lo, hi


def _virtual_anchor(x):
    """values near 1 (second channel, |ln u| tiny) seen as anchors 1(1 + delta), delta = e^ln u - 1:
    ln|delta| in ln|ln u| +- 2|ln u|"""
    m = (x.AQ == 0) & (x.S != 0) & ~x.O & ~x.Z & (x.N < DTINY)
    if not m.any():
        return x
    y = x.take(np.arange(len(x.A)))
    e2 = 2 * exp_hi(x.N[m])
    y.AP = x.AP.copy(); y.AQ = x.AQ.copy(); y.DS = x.DS.copy(); y.DM = x.DM.copy(); y.DN = x.DN.copy()
    y.AP[m] = 1; y.AQ[m] = 1; y.DS[m] = x.S[m]; y.DM[m] = addlo(x.M[m], -e2); y.DN[m] = addhi(x.N[m], e2)
    return y


def anchors(f, x, y, r):
    """set the anchor fields of r = f(x, y) (and tighten the value / the near-1 channel)"""
    n = r.A.shape
    x = _virtual_anchor(x); y = _virtual_anchor(y)
    if not ((x.AQ > 0).any() or (y.AQ > 0).any()):
        return
    with np.errstate(all="ignore"):
        ex1 = (x.AQ > 0) & (x.DS == 0) & ~x.Z; ex2 = (y.AQ > 0) & (y.DS == 0) & ~y.Z
        an1 = (x.AQ > 0) & (x.DS != 0); an2 = (y.AQ > 0) & (y.DS != 0)
        AP = np.zeros(n, np.int64); AQ = np.zeros(n, np.int64); DS = np.zeros(n, np.int8)
        DM = np.full(n, -INF); DN = np.full(n, -INF)

        def put(m, P, Q, ds, dm, dn):
            if m.any():
                AP[m] = P[m]; AQ[m] = Q[m]; DS[m] = ds[m] if np.ndim(ds) else ds
                DM[m] = dm[m]; DN[m] = dn[m]

        def setval(m, lo, hi):                      # the value itself: ln u in [lo, hi]
            if m.any():
                r.A[m] = np.maximum(r.A[m], lo[m]); r.B[m] = np.minimum(r.B[m], hi[m])

        if f in (0, 1):
            for (e1, a1, xx, e2, a2, yy) in ((ex1, an1, x, ex2, an2, y), (ex2, an2, y, ex1, an1, x)):
                cand = (e1 | a1) & ~(yy.AQ > 0) & ~yy.Z & ~yy.O & (yy.B < -40 + np.log(np.maximum(xx.AP, 1) / np.maximum(xx.AQ, 1)) + 1)
                ql, qh = _lnq(xx.AP, xx.AQ, cand)
                # exact q + tiny t (t not anchored, not exact): delta = +-t/q
                t = e1 & ~(yy.AQ > 0) & ~yy.Z & ~yy.O & (yy.B < ql - 40)
                put(t, xx.AP, xx.AQ, 1 if f == 0 else -1, addlo(yy.A, -qh), addhi(yy.B, -ql))
                # anchored q(1+d) + tiny t: same sign of d and +-t: |delta| = |d| + t/q
                sgn_t = 1 if f == 0 else -1
                t2 = a1 & ~(yy.AQ > 0) & ~yy.Z & ~yy.O & (yy.B < ql - 40) & (xx.DS == sgn_t)
                put(t2, xx.AP, xx.AQ, xx.DS, lse_lo(xx.DM, addlo(yy.A, -qh)), lse_hi(xx.DN, addhi(yy.B, -ql)))
            if f == 0:
                # exact + exact: exact rational
                m = ex1 & ex2
                if m.any():
                    P, Q, sg, okr = rat_op("add", x.AP, x.AQ, y.AP, y.AQ)
                    mm = m & okr
                    put(mm, P, Q, 0, DM, DN)
                # anchored + exact: (q1 + q2)(1 + q1 d1/(q1+q2))
                for (a1, xx, e2, yy) in ((an1, x, ex2, y), (an2, y, ex1, x)):
                    m = a1 & e2
                    if m.any():
                        P, Q, sg, okr = rat_op("add", xx.AP, xx.AQ, yy.AP, yy.AQ)
                        mm = m & okr & (P > 0)
                        l1l, l1h = _lnq(xx.AP, xx.AQ, mm)
                        l3l, l3h = _lnq(P, Q, mm)
                        put(mm, P, Q, xx.DS, addlo(addlo(l1l, -l3h), xx.DM), addhi(addhi(l1h, -l3l), xx.DN))
            else:
                # |exact - exact|: exact (zero handled by fwd)
                m = ex1 & ex2
                if m.any():
                    P, Q, sg, okr = rat_op("sub", x.AP, x.AQ, y.AP, y.AQ)
                    mm = m & okr & (P > 0)
                    put(mm, P, Q, 0, DM, DN)
                for (e1, xx, a2, yy) in ((ex1, x, an2, y), (ex2, y, an1, x)):
                    m = e1 & a2
                    if not m.any():
                        continue
                    same = m & (xx.AP == yy.AP) & (xx.AQ == yy.AQ)
                    if same.any():                     # |q - q(1+d)| = q |d|
                        ql, qh = _lnq(yy.AP, yy.AQ, same)
                        setval(same, addlo(ql, yy.DM), addhi(qh, yy.DN))
                        r.Z[same] = False
                    diff = m & ~same
                    if diff.any():                     # |D - q2 d2|, D = q1 - q2
                        P, Q, sg, okr = rat_op("sub", xx.AP, xx.AQ, yy.AP, yy.AQ)
                        mm = diff & okr & (P > 0)
                        l2l, l2h = _lnq(yy.AP, yy.AQ, mm)
                        lDl, lDh = _lnq(P, Q, mm)
                        ds = np.where(sg > 0, -yy.DS, yy.DS).astype(np.int8)
                        put(mm, P, Q, ds, addlo(addlo(l2l, -lDh), yy.DM), addhi(addhi(l2h, -lDl), yy.DN))
        elif f == 2:
            m = ex1 & ex2
            if m.any():
                P, Q, sg, okr = rat_op("mul", x.AP, x.AQ, y.AP, y.AQ)
                put(m & okr, P, Q, 0, DM, DN)
            for (e1, xx, a2, yy) in ((ex1, x, an2, y), (ex2, y, an1, x)):
                m = e1 & a2
                if m.any():
                    P, Q, sg, okr = rat_op("mul", xx.AP, xx.AQ, yy.AP, yy.AQ)
                    put(m & okr, P, Q, yy.DS, yy.DM, yy.DN)
        elif f >= 5:
            # q^(+-k(1+d)) = q^(+-k) e^(+-k d ln q): anchor at q^(+-k) (exact), delta' = e^(+-k d ln q) - 1,
            # ln|delta'| in ln k + ln|d| + ln|ln q| +- 2|k d ln q|, sign = sign(d) sign(ln q) (+-)
            base, expo = (x, y) if f in (5, 6) else (y, x)
            negx = f in (6, 8)
            eb = (base.AQ > 0) & (base.DS == 0) & ~base.Z & (base.AP != base.AQ)
            ae = (expo.AQ == 1) & (expo.DS != 0) & (expo.AP >= 1) & (expo.AP <= 6)
            m = eb & ae
            if m.any():
                for i in np.where(m)[0]:
                    p, q, k = int(base.AP[i]), int(base.AQ[i]), int(expo.AP[i])
                    pp, qq = (p ** k, q ** k) if not negx else (q ** k, p ** k)
                    if pp >= 2 ** 62 or qq >= 2 ** 62:
                        continue
                    lql, lqh = lnq_bounds(np.array([p], np.int64), np.array([q], np.int64))
                    lq_abs_lo = min(abs(lql[0]), abs(lqh[0])) if lql[0] * lqh[0] > 0 else 0.0
                    lq_abs_hi = max(abs(lql[0]), abs(lqh[0]))
                    if lq_abs_lo <= 0:
                        continue
                    sg = int(expo.DS[i]) * (1 if lql[0] > 0 else -1) * (-1 if negx else 1)
                    # outward sums (corrected 02/11: the plain double sums could round inwards by up to 1 ulp)
                    dm = float(addlo(addlo(np.array([expo.DM[i]]), ln_lo(np.array([lq_abs_lo]))), ln_lo(np.array([float(k)])))[0])
                    dn = float(addhi(addhi(np.array([expo.DN[i]]), ln_hi(np.array([lq_abs_hi]))), ln_hi(np.array([float(k)])))[0])
                    e2 = 2 * float(exp_hi(np.array([dn]))[0])
                    bump("ancexp_fired")
                    _odm = expo.DM[i] + float(ln_lo(np.array([lq_abs_lo]))[0]) + float(ln_lo(np.array([float(k)]))[0])
                    _odn = expo.DN[i] + float(ln_hi(np.array([lq_abs_hi]))[0]) + float(ln_hi(np.array([float(k)]))[0])
                    # was the bound of the code before 02/11 (plain double sums) above the exact combination of its
                    # own input bounds (lower side) or below it (upper side)?  Exact rational arithmetic on the doubles.
                    from fractions import Fraction as _Fr
                    _a = [float(expo.DM[i]), float(ln_lo(np.array([lq_abs_lo]))[0]), float(ln_lo(np.array([float(k)]))[0])]
                    _b = [float(expo.DN[i]), float(ln_hi(np.array([lq_abs_hi]))[0]), float(ln_hi(np.array([float(k)]))[0])]
                    _olo = float(pred(np.array([_odm - e2]))[0]); _ohi = float(succ(np.array([_odn + e2]))[0])
                    if all(np.isfinite(_a + _b + [e2, _olo, _ohi])):
                        bump("ancexp_old_invalid", (_Fr(_olo) > sum(map(_Fr, _a)) - _Fr(e2)) or (_Fr(_ohi) < sum(map(_Fr, _b)) + _Fr(e2)))
                    put(np.arange(n[0]) == i, np.full(n, pp, np.int64), np.full(n, qq, np.int64), sg,
                        np.full(n, float(addlo(np.array([dm]), np.array([-e2]))[0])),
                        np.full(n, float(addhi(np.array([dn]), np.array([e2]))[0])))
            # (03/11) the same with an anchored RATIONAL exponent t (1 + d), t = p/r (r > 1, p <= 6), when q^t is rational
            # (P^p and Q^p exact r-th powers):  q^(t(1+d)) = q^t q^(t d) = q^t (1 + delta'), delta' = e^z - 1 with
            # z = t d ln q, so sign(delta') = sign(d) sign(ln q) (+-), and ln|delta'| = ln|d| + ln|ln q| + ln t + O(|z|),
            # bounded by +-2 e^(DN') as in the integer case (DN' the upper bound of ln|z|)
            aer = (expo.AQ > 1) & (expo.AQ <= 64) & (expo.DS != 0) & (expo.AP >= 1) & (expo.AP <= 6)
            mr7 = eb & aer & ~(AQ > 0)
            rows7 = np.where(mr7)[0]
            if len(rows7):
                p7 = expo.AP[rows7].astype(np.int64); r7 = expo.AQ[rows7].astype(np.int64)
                Pb = base.AP[rows7].astype(np.int64); Qb = base.AQ[rows7].astype(np.int64)
                rp7 = _roots_vec(Pb, r7); rq7 = _roots_vec(Qb, r7)
                ok7 = (rp7 > 0) & (rq7 > 0)
                Pt7, okp7 = ipow(np.where(ok7, rp7, 1), np.where(ok7, p7, 1))
                Qt7, okq7 = ipow(np.where(ok7, rq7, 1), np.where(ok7, p7, 1))
                ok7 &= okp7 & okq7
                lql7, lqh7 = lnq_bounds(Pb, Qb)
                same_sign = (lql7 * lqh7) > 0
                ok7 &= same_sign
                if ok7.any():
                    rr = rows7[ok7]
                    labs_lo = np.minimum(np.abs(lql7), np.abs(lqh7))[ok7]; labs_hi = np.maximum(np.abs(lql7), np.abs(lqh7))[ok7]
                    tf = p7[ok7].astype(np.float64) / r7[ok7].astype(np.float64)
                    t_lo = pred(tf); t_hi = succ(tf)
                    dm = addlo(addlo(expo.DM[rr], ln_lo(labs_lo)), ln_lo(t_lo))
                    dn = addhi(addhi(expo.DN[rr], ln_hi(labs_hi)), ln_hi(t_hi))
                    e2 = 2 * exp_hi(dn)
                    PP = np.where(not negx, Pt7[ok7], Qt7[ok7]); QQ = np.where(not negx, Qt7[ok7], Pt7[ok7])
                    gg = np.gcd(PP, QQ); gg = np.where(gg == 0, 1, gg)
                    AP[rr] = PP // gg; AQ[rr] = QQ // gg
                    DS[rr] = (expo.DS[rr] * np.where(lql7[ok7] > 0, 1, -1) * (-1 if negx else 1)).astype(np.int8)
                    DM[rr] = addlo(dm, -e2); DN[rr] = addhi(dn, e2)
                    bump("ancexp_rat", int(ok7.sum()))
            # (added 02/11) base exact or anchored, q (1 + d) with q = P/Q, and exponent an exact integer k >= 1:
            #   (q (1 + d))^(+-k) = q^(+-k) (1 + d)^(+-k):  exact q^(+-k) when d = 0; otherwise an anchor at q^(+-k)
            #   with delta' = e^z - 1, z = +-k ln(1 + d), so sign(delta') = +-sign(d) and
            #   ln|delta'| = ln k + ln|d| + ln(ln(1+d)/d) + ln((e^z - 1)/z),  |ln(ln(1+d)/d)| <= 2|d|, |ln((e^z-1)/z)| <= |z| <= 2k|d|
            bq = (base.AQ > 0) & ~base.Z & ~base.O & (base.AP > 0)
            ek = (expo.AQ == 1) & (expo.DS == 0) & (expo.AP >= 1) & ~expo.Z
            m = bq & ek & ~(AQ > 0) & ANCPOW
            if m.any():
                bump("ancpow_rows", m.sum())
                kk = np.where(m, expo.AP, 1)
                P_, Q_ = (base.AP, base.AQ) if not negx else (base.AQ, base.AP)
                pp, okp = ipow(np.where(m, P_, 1), kk); qq, okq = ipow(np.where(m, Q_, 1), kk)
                mm = m & okp & okq
                put(mm & (base.DS == 0), pp, qq, 0, DM, DN)
                an = mm & (base.DS != 0)
                if an.any():
                    kf = kk.astype(np.float64)
                    e2 = np.full(n, INF); e2[an] = succ(2.0 * (kf[an] + 1.0) * exp_hi(base.DN[an]))
                    dm = np.full(n, -INF); dn = np.full(n, INF)
                    dm[an] = addlo(addlo(base.DM[an], ln_lo(kf[an])), -e2[an])
                    dn[an] = addhi(addhi(base.DN[an], ln_hi(kf[an])), e2[an])
                    put(an, pp, qq, (base.DS * (-1 if negx else 1)).astype(np.int8), dm, dn)
            # (03/11) the same with an exact rational exponent t = p/r (r > 1) when q^t is rational: q = P/Q with
            # P^p and Q^p exact r-th powers (integer roots, checked by exponentiation).  (q(1 + d))^t = q^t (1 + d)^t,
            # and the bound of rule 8 holds with k replaced by t: |ln(ln(1+d)/d)| <= 2|d|, |z| = |t ln(1 + d)| <= 2t|d|
            er = (expo.AQ > 1) & (expo.DS == 0) & (expo.AP >= 1) & ~expo.Z
            mr = bq & er & ~(AQ > 0) & ANCPOW & (expo.AQ <= 64)
            rows = np.where(mr)[0]
            if len(rows):
                p_ = expo.AP[rows].astype(np.int64); r_ = expo.AQ[rows].astype(np.int64)
                P0 = (base.AP if not negx else base.AQ)[rows].astype(np.int64)
                Q0 = (base.AQ if not negx else base.AP)[rows].astype(np.int64)
                rp = _roots_vec(P0, r_); rq = _roots_vec(Q0, r_)      # exact integer r-th roots (0: none)
                okr = (rp > 0) & (rq > 0)
                Pt, okp = ipow(np.where(okr, rp, 1), np.where(okr, p_, 1))
                Qt, okq = ipow(np.where(okr, rq, 1), np.where(okr, p_, 1))
                ok_ = okr & okp & okq
                if ok_.any():
                    gg = np.gcd(Pt, Qt); gg = np.where(gg == 0, 1, gg); Pt = Pt // gg; Qt = Qt // gg
                    bump("ancpow_rat", int(ok_.sum()))
                    ex_ = ok_ & (base.DS[rows] == 0)
                    if ex_.any():
                        rr = rows[ex_]
                        AP[rr] = Pt[ex_]; AQ[rr] = Qt[ex_]; DS[rr] = 0; DM[rr] = -INF; DN[rr] = -INF
                    an_ = ok_ & (base.DS[rows] != 0)
                    if an_.any():
                        rr = rows[an_]
                        tf = p_[an_].astype(np.float64) / r_[an_].astype(np.float64)
                        t_lo = pred(tf); t_hi = succ(tf)
                        e2 = succ(2.0 * (t_hi + 1.0) * exp_hi(base.DN[rr]))
                        AP[rr] = Pt[an_]; AQ[rr] = Qt[an_]
                        DS[rr] = (base.DS[rr] * (-1 if negx else 1)).astype(np.int8)
                        DM[rr] = addlo(addlo(base.DM[rr], ln_lo(t_lo)), -e2)
                        DN[rr] = addhi(addhi(base.DN[rr], ln_hi(t_hi)), e2)
        elif f in (3, 4):
            num, den = (x, y) if f == 3 else (y, x)
            e_n = (num.AQ > 0) & (num.DS == 0) & ~num.Z; e_d = (den.AQ > 0) & (den.DS == 0) & ~den.Z
            a_n = (num.AQ > 0) & (num.DS != 0); a_d = (den.AQ > 0) & (den.DS != 0)
            P, Q, sg, okr = rat_op("div", num.AP, num.AQ, den.AP, den.AQ)
            put(e_n & e_d & okr, P, Q, 0, DM, DN)
            put(a_n & e_d & okr, P, Q, num.DS, num.DM, num.DN)
            m = e_n & a_d & okr                         # 1/(1+d) - 1: sign flips, |.| = |d| (1 + O(d))
            if m.any():
                e2 = 2 * exp_hi(den.DN)
                put(m, P, Q, (-den.DS).astype(np.int8), addlo(den.DM, -e2), addhi(den.DN, e2))
        # drop anchors whose delta is not tiny
        bad = (AQ > 0) & (DS != 0) & ~(DN <= DTINY)
        AQ[bad] = 0; DS[bad] = 0
        r.AP, r.AQ, r.DS, r.DM, r.DN = AP, AQ, DS, DM, DN
        # anchor 1 with delta: the value is near 1, ln u = ln(1+d) has the sign of d, ln|ln u| in ln|d| +- 2|d|
        one = (AQ > 0) & (AP == AQ) & (DS != 0)
        if one.any():
            e2 = 2 * exp_hi(DN[one])
            r.S[one] = DS[one]; r.M[one] = addlo(DM[one], -e2); r.N[one] = addhi(DN[one], e2)
        # value from the anchor: ln u = ln q + ln(1+d), |ln(1+d)| <= 2|d|
        an = (AQ > 0) & (DS != 0)
        if an.any():
            ql, qh = lnq_bounds(np.maximum(AP[an], 1), np.maximum(AQ[an], 1))
            e2 = 2 * exp_hi(DN[an])
            idx = np.where(an)[0]
            r.A[idx] = np.maximum(r.A[idx], addlo(ql, -e2)); r.B[idx] = np.minimum(r.B[idx], addhi(qh, e2))
        exr = (AQ > 0) & (DS == 0)
        if exr.any():
            r.Z[exr] = False
            r.O[exr] = (AP[exr] == AQ[exr])
