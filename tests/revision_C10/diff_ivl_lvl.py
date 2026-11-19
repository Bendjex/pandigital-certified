"""Differential test (reviewer's own): random expression trees over digit literals, evaluated
 (a) with ivl.fwd chains (motor X enclosures, all channels), and
 (b) with lvl.op (the package's independent evaluator) and, when the value is an exact rational, exactly.
Any statement of the ivl enclosure contradicted by (b) is reported."""
import sys, os, random, math
PC = sys.argv[1]
sys.path.insert(0, os.path.join(PC, "c5")); sys.path.insert(1, os.path.join(PC, "src"))
import numpy as np, mpmath
from fractions import Fraction as Fr
import ivl, motorx, lvl
sys.set_int_max_str_digits(0)
from mpmath import iv, mp
mp.prec = 600
seed = int(sys.argv[2]); NT = int(sys.argv[3]); MAXLEAF = int(sys.argv[4]) if len(sys.argv) > 4 else 6
rnd = random.Random(seed)

def leaf_ivl(q):
    x = np.array([q.numerator / q.denominator])
    lo = ivl.ln_lo(np.nextafter(x, 0)); hi = ivl.ln_hi(np.nextafter(x, np.inf))
    P = np.array([q.numerator], np.int64); Q = np.array([q.denominator], np.int64)
    e = motorx.Ent(lo, hi, np.zeros(1, bool), P == Q, None, P, Q)
    if q == 1:
        e.A[:] = 0.0; e.B[:] = 0.0
        S, M, N = ivl.mu_default(e.A, e.B); e.S, e.M, e.N = S, M, N
    return e.q7()

LITS = []
for d in range(1, 10):
    LITS += [Fr(d), Fr(d, 10)]
for a in range(1, 10):
    for b in range(1, 10):
        if a != b:
            LITS += [Fr(10 * a + b), Fr(10 * a + b, 10), Fr(10 * a + b, 100)]

def rand_tree(nleaf):
    if nleaf == 1:
        q = rnd.choice(LITS if rnd.random() < 0.25 else LITS[:18])
        return ("L", q)
    k = rnd.randint(1, nleaf - 1)
    f = rnd.choices(range(9), weights=[2, 3, 2, 2, 2, 4, 4, 4, 4])[0]
    return ("O", f, rand_tree(k), rand_tree(nleaf - k))

class Undef(Exception):
    pass

def ev_ivl(t):
    if t[0] == "L":
        return leaf_ivl(t[1])
    _, f, a, b = t
    x = ev_ivl(a); y = ev_ivl(b)
    R = None
    if f == 1:
        R = (np.where(x.DS == 0, x.AP, 0), np.where(x.DS == 0, x.AQ, 0), np.where(y.DS == 0, y.AP, 0), np.where(y.DS == 0, y.AQ, 0))
    r, ok = ivl.fwd(f, x, y, R)
    if not ok[0]:
        raise Undef()
    return r

def ev_lvl(t):
    if t[0] == "L":
        return lvl.exact(t[1])
    _, f, a, b = t
    return lvl.op(lvl.FORMS[f], ev_lvl(a), ev_lvl(b))

def show(t):
    if t[0] == "L":
        return str(t[1])
    return f"[{ivl.FORMS[t[1]]} | {show(t[2])} | {show(t[3])}]"

def lnq_mp(q):
    """ln q for a positive Fraction, without cancellation near 1"""
    d = q - 1
    if d == 0:
        return mp.mpf(0)
    if abs(d) < Fr(1, 4):
        dm = mp.mpf(d.numerator) / mp.mpf(d.denominator)      # relative precision mp.prec
        return mp.log1p(dm)
    return mp.log(q.numerator) - mp.log(q.denominator)

stats = {"trees": 0, "lvl_fail": 0, "lvl_invalid": 0, "ivl_undef": 0, "checked": 0, "bad": 0, "nonstrictS": 0, "nonstrictDS": 0,
         "exactq": 0, "far": 0, "S_set": 0, "anc": 0, "T3": 0, "FE": 0, "off": 0, "Z": 0, "O": 0}
bad_list = []

def report(kind, t, r, v, extra=""):
    stats["bad"] += 1
    print("BAD-NOW", kind, "|", extra, "|", show(t), flush=True)
    if len(bad_list) < 40:
        bad_list.append((kind, show(t), {k: getattr(r, k)[0] for k in ivl.Q7.FIELDS}, repr(v), extra))

def check(t, r, v):
    """r: ivl Q7 (1 entry); v: lvl Val"""
    A, B = float(r.A[0]), float(r.B[0]); Z, O = bool(r.Z[0]), bool(r.O[0]); S = int(r.S[0]); M, N = float(r.M[0]), float(r.N[0])
    stats["checked"] += 1
    if v.far:
        stats["far"] += 1
        # only the classification: huge/tiny must not contradict a bounded enclosure
        if v.far == "huge" and B < 1e300: report("far-huge vs B", t, r, v)
        if v.far == "tiny" and A > -1e300: report("far-tiny vs A", t, r, v)
        if v.far == "near1" and (A > 1e-7 or B < -1e-7): report("far-near1", t, r, v)
        return
    if v.q is not None:
        stats["exactq"] += 1
        q = v.q
        if Z:
            stats["Z"] += 1
            if q != 0: report("Z but q != 0", t, r, v)
            return
        if O:
            stats["O"] += 1
            if q != 1: report("O but q != 1", t, r, v)
            return
        if q == 0:
            if A != -np.inf: report("q == 0 but A > -inf", t, r, v)
            if S > 0: report("q == 0 but S > 0", t, r, v)
            if r.AQ[0] > 0 and r.DS[0] == 0 and r.AP[0] != 0: report("q == 0 but exact anchor", t, r, v)
            return
        L = lnq_mp(q)
        if not (mp.mpf(A) <= L <= mp.mpf(B)): report("ln q outside [A,B]", t, r, v, f"ln={mp.nstr(L, 25)}")
        if S != 0:
            stats["S_set"] += 1
            if L == 0:
                stats["nonstrictS"] += 1
                if stats["nonstrictS"] <= 3: print("NONSTRICT-S example:", show(t), {k: getattr(r, k)[0] for k in ("A","B","S","M","N","AP","AQ","DS","DM","DN")})
                if M != -np.inf: report("S != 0, ln u = 0, M finite", t, r, v)
            else:
                if (L > 0) != (S > 0): report("sign S wrong", t, r, v)
                LL = mp.log(abs(L))
                if not (mp.mpf(M) <= LL <= mp.mpf(N)): report("ln|ln q| outside [M,N]", t, r, v, f"lnln={mp.nstr(LL, 25)}")
        if r.AQ[0] > 0:
            stats["anc"] += 1
            a = Fr(int(r.AP[0]), int(r.AQ[0]))
            if a == 0:
                report("anchor 0 for q != 0", t, r, v); return
            d = q / a - 1
            if r.DS[0] == 0:
                if d != 0: report("exact anchor wrong", t, r, v, f"anchor {a}")
            elif d == 0:
                stats["nonstrictDS"] += 1
                if stats["nonstrictDS"] <= 3: print("NONSTRICT-DS example:", show(t), {k: getattr(r, k)[0] for k in ("A","B","S","M","N","AP","AQ","DS","DM","DN")})
                if r.DM[0] != -np.inf: report("DS != 0, delta = 0, DM finite", t, r, v)
            else:
                if (d > 0) != (r.DS[0] > 0): report("sign DS wrong", t, r, v, f"anchor {a}")
                LD = lnq_mp(abs(d))
                if not (mp.mpf(float(r.DM[0])) <= LD <= mp.mpf(float(r.DN[0]))): report("ln|delta| outside [DM,DN]", t, r, v, f"lnd={mp.nstr(LD, 25)} anchor {a}")
        if r.FE[0] != 0:
            stats["FE"] += 1
            base = Fr(int(r.FP[0]), int(r.FQ[0])); e = Fr(int(r.FE[0]), max(int(r.FD[0]), 1))
            c = Fr(int(r.OP[0]), int(r.OQ[0])) if r.OQ[0] > 0 else Fr(0)
            if r.OQ[0] > 0: stats["off"] += 1
            X = q - c
            if X <= 0: report("factored: q - c <= 0", t, r, v)
            else:
                d = lnq_mp(X) - mp.mpf(e.numerator) / e.denominator * lnq_mp(base)
                if abs(d) > mp.mpf(2) ** -500 * (1 + abs(lnq_mp(X))): report("factored triple wrong", t, r, v, f"{base}^{e}+{c}")
        if r.T3[0]:
            stats["T3"] += 1
            report("T3 on exact rational (check by hand)", t, r, v)
        return
    # interval value
    try:
        lam = lvl.lam_iv(v)
    except (lvl.Fail, OverflowError, MemoryError):
        lam = None
    if Z: report("Z but lvl value not exact 0", t, r, v); return
    if O:
        if not (v.s == 0): report("O but lvl value not exact 1", t, r, v)
        return
    if v.s == 0:
        if not (A <= 0 <= B): report("lvl exact 1 outside [A,B]", t, r, v)
        return
    if lam is not None:
        if lam.b < A or lam.a > B: report("lvl ln u disjoint from [A,B]", t, r, v, f"lam=[{mpmath.nstr(lam.a, 20)}, {mpmath.nstr(lam.b, 20)}]")
    elif v.s is not None and v.Lam is not None:
        # ln u = s e^Lam beyond intervals: compare at the second level with [A, B]
        if v.s > 0 and np.isfinite(B) and B > 0 and v.Lam.a > mp.log(mp.mpf(B)): report("lvl huge above B", t, r, v)
        if v.s > 0 and B <= 0: report("lvl ln>0 but B<=0", t, r, v)
        if v.s < 0 and np.isfinite(A) and A < 0 and v.Lam.a > mp.log(mp.mpf(-A)): report("lvl tiny below A", t, r, v)
        if v.s < 0 and A >= 0 and not (A == 0 and False): report("lvl ln<0 but A>=0", t, r, v)
    if S != 0 and v.s is not None and v.Lam is not None:
        stats["S_set"] += 1
        if v.s != S: report("sign S vs lvl s", t, r, v)
        elif v.Lam.b < M or v.Lam.a > N: report("lvl Lam disjoint from [M,N]", t, r, v, f"Lam=[{mpmath.nstr(v.Lam.a, 20)}, {mpmath.nstr(v.Lam.b, 20)}]")
    if r.AQ[0] > 0:
        stats["anc"] += 1
        a = Fr(int(r.AP[0]), int(r.AQ[0]))
        if r.DS[0] == 0:
            # ivl states the value is exactly a; lvl only has an enclosure: it must contain a
            xa = lvl.ivq(a)
            if v.x is not None and (v.x.b < xa.a or v.x.a > xa.b): report("exact anchor outside lvl interval", t, r, v, f"anchor {a}")
            elif v.x is None and lam is not None:
                la = lvl.lnq(a)
                if lam.b < la.a or lam.a > la.b: report("exact anchor outside lvl ln interval", t, r, v, f"anchor {a}")
            stats["exact_vs_interval"] = stats.get("exact_vs_interval", 0) + 1
        else:
            # u = a (1 + delta): ln u - ln a = ln(1 + delta)
            if v.anc is not None and v.anc[0] == a:
                aq, sg, D = v.anc
                if sg != r.DS[0]: report("sign DS vs lvl anc", t, r, v)
                elif D.b < r.DM[0] or D.a > r.DN[0]: report("lvl D disjoint from [DM,DN]", t, r, v, f"D=[{mpmath.nstr(D.a, 20)}, {mpmath.nstr(D.b, 20)}]")
            elif lam is not None:
                la = lvl.lnq(a)
                dl = lam - la                      # ln(1 + delta)
                if dl.a > 0 and r.DS[0] < 0: report("sign DS vs ln u - ln a", t, r, v)
                if dl.b < 0 and r.DS[0] > 0: report("sign DS vs ln u - ln a", t, r, v)
                if dl.a > 0 or dl.b < 0:
                    ad = abs(iv.exp(dl) - 1)
                    if ad.a > 0:
                        LD = iv.log(ad)
                        if LD.b < r.DM[0] or LD.a > r.DN[0]: report("ln|delta| (from lvl ln) disjoint from [DM,DN]", t, r, v, f"LD=[{mpmath.nstr(LD.a, 20)}, {mpmath.nstr(LD.b, 20)}]")
    if r.T3[0]:
        stats["T3"] += 1
        if v.s is not None and v.Lam is not None:
            # ln|ln u| = -e^nu, nu in [M3, N3]
            if v.Lam.a >= 0: report("T3 but Lam >= 0", t, r, v)
            else:
                nu = iv.log(-v.Lam)
                if nu.b < r.M3[0] or nu.a > r.N3[0]: report("lvl nu disjoint from [M3,N3]", t, r, v, f"nu=[{mpmath.nstr(nu.a, 20)}, {mpmath.nstr(nu.b, 20)}]")
    if r.FE[0] != 0:
        stats["FE"] += 1
        if r.OQ[0] > 0: stats["off"] += 1
        base = Fr(int(r.FP[0]), int(r.FQ[0])); e = Fr(int(r.FE[0]), max(int(r.FD[0]), 1))
        lX = iv.mpf(e.numerator) / e.denominator * lvl.lnq(base)
        if r.OQ[0] == 0 and lam is not None:
            if lam.b < lX.a or lam.a > lX.b: report("factored triple vs lvl ln", t, r, v, f"{base}^{e}")
        if v.pf is not None and r.OQ[0] == 0:
            lp = lvl.lnpf(v.pf)
            if lp.b < lX.a or lp.a > lX.b: report("factored triple vs lvl pf", t, r, v, f"{base}^{e}")

import time as _t
_t0 = _t.time()
for it in range(NT):
    if it % 2000 == 0 and it: print("progress", it, round(_t.time() - _t0), stats["bad"], flush=True)
    t = rand_tree(rnd.randint(2, MAXLEAF))
    stats["trees"] += 1
    try:
        v = ev_lvl(t)
    except lvl.Invalid:
        stats["lvl_invalid"] += 1
        continue
    except lvl.Fail:
        stats["lvl_fail"] += 1
        continue
    except (OverflowError, MemoryError, ValueError, ZeroDivisionError) as er:
        stats["lvl_fail"] += 1
        continue
    try:
        r = ev_ivl(t)
    except Undef:
        stats["ivl_undef"] += 1
        report("ivl says undefined, lvl has a value", t, ivl.Q7(np.zeros(1), np.zeros(1), np.zeros(1, bool), np.zeros(1, bool)), v)
        continue
    try:
        check(t, r, v)
    except Exception as er:
        report("exception in check: " + repr(er), t, r, v)
print(stats)
for b in bad_list:
    print("BAD", b[0], "|", b[4]); print("   ", b[1]); print("   ", b[2]); print("   lvl:", b[3])
