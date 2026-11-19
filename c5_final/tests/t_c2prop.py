"""soundness of the propagation of second-level window constraints (03/11): qx.level2_lin (forms s+l, |s-l|, s*l, s/l,
l/s), qx.pre_level2 and qx.pre2 with a constraint c2 = (sign of w, bounds of ln|w|) on the windows w = ln u.  Random s
of moderate size (both signs of ln s), random u beyond the double range (|ln u| = e^mu, mu in [700, 3000], both
signs) or inside it; the exact l with f(s, l) = u is computed with mpmath; whenever the derived constraint is active,
l must satisfy it (sign of ln l, ln|ln l| in [m, M]), and for s^l, s^-l, ln l must lie in a window of pre2.
usage: t_c2prop.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ivl, qx
mpmath.mp.dps = 60
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 83)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
DM = 1.7976931348623157e308


def nx(v, d):
    return np.nextafter(float(v), d)


def q7_of(lam):
    """enclosure of s = e^lam (lam an mpf): first level and second level"""
    A = np.array([nx(lam, -np.inf)]); B = np.array([nx(lam, np.inf)])
    if lam == 0:
        return ivl.Q7(A, B, np.array([False]), np.array([False]))
    m = mpmath.log(abs(lam))
    return ivl.Q7(A, B, np.array([False]), np.array([False]), np.array([1 if lam > 0 else -1], np.int8),
                  np.array([nx(m, -np.inf)]), np.array([nx(m, np.inf)]))


bad = tested = active = 0
for t in range(N):
    lam = mpmath.mpf(float(rng.normal(0, 3))) if rng.random() < 0.8 else mpmath.mpf(float(rng.choice([-1, 1]) * 10 ** rng.uniform(0, 2.5)))
    x7 = q7_of(lam)
    sgw = int(rng.choice([-1, 1]))
    rr = rng.random()
    if rr < 0.65:
        mu = mpmath.mpf(float(rng.uniform(700, 3000)))
    elif rr < 0.8:
        mu = mpmath.mpf(float(rng.uniform(3, 700)))
    else:                                                   # |ln u| below the subnormals: u within 1e-324 of 1
        mu = mpmath.mpf(float(rng.uniform(-3000, -746)))
    w = sgw * mpmath.exp(mu)                                  # ln u
    eps = float(10 ** rng.uniform(-14, -2))
    c2 = (np.array([sgw], np.int8), np.array([float(mu) - eps]), np.array([float(mu) + eps]))
    if mu < -745:
        al, be = np.array([-5e-324]), np.array([5e-324])   # the first level straddles ln 1 = 0
    elif mu > 709.7:
        al, be = (np.array([DM]), np.array([np.inf])) if sgw > 0 else (np.array([-np.inf]), np.array([-DM]))
    else:
        al, be = np.array([nx(w * (1 - 1e-12), -np.inf)]), np.array([nx(w * (1 + 1e-12), np.inf)])
        if sgw < 0:
            al, be = np.array([nx(w * (1 + 1e-12), -np.inf)]), np.array([nx(w * (1 - 1e-12), np.inf)])
    f = int(rng.integers(0, 9))
    s = mpmath.exp(lam)
    # the exact l values (ln l) with f(s, l) = u
    lls = []
    if f == 2: lls = [w - lam]
    elif f == 3: lls = [lam - w]
    elif f == 4: lls = [w + lam]
    elif f == 0:
        if w > lam:                                         # l = u - s > 0
            lls = [w + mpmath.log(-mpmath.expm1(lam - w))]
    elif f == 1:
        lls = [w + mpmath.log1p(mpmath.exp(lam - w))]       # l = u + s
        if w > lam: lls.append(w + mpmath.log(-mpmath.expm1(lam - w)))     # l = s + u... l = u - s
        elif w < lam: lls.append(lam + mpmath.log(-mpmath.expm1(w - lam)))  # l = s - u
    elif f in (5, 6):
        if lam != 0:
            y = w / lam if f == 5 else -w / lam              # l = y must be > 0
            if y > 0: lls = [mpmath.log(y)]
    else:                                                     # l^s, l^-s: s ln l = +-w
        lls = [(w if f == 7 else -w) / s]                     # ln l
        lls = [("lnl", lls[0])]
    if not lls:
        continue
    tested += 1
    if f in (7, 8):
        sg, m2, M2 = qx.pre_level2(f, x7, al, be, c2)
        lnl = lls[0][1]
        if sg[0] != 0:
            active += 1
            ok = (lnl > 0) == (sg[0] > 0) and m2[0] <= mpmath.log(abs(lnl)) <= M2[0]
            if not ok:
                bad += 1; print("MISSED l^s", f, mpmath.nstr(lam, 8), sgw, mpmath.nstr(mu, 10), sg[0], m2[0], M2[0], mpmath.nstr(mpmath.log(abs(lnl)), 12))
        continue
    if f in (5, 6):
        wins = qx.pre2(f, x7, al, be, c2)
        for ll in lls:
            if not any(lo[0] <= ll <= hi[0] for lo, hi in wins):
                bad += 1; print("MISSED pre2", f, mpmath.nstr(lam, 8), sgw, mpmath.nstr(mu, 10), mpmath.nstr(ll, 12), [(lo[0], hi[0]) for lo, hi in wins])
        active += 1
        continue
    sg, m2, M2 = qx.level2_lin(f, x7, al, be, c2)
    if sg[0] != 0:
        active += 1
        for ll in lls:
            ok = ll != 0 and (ll > 0) == (sg[0] > 0) and m2[0] <= mpmath.log(abs(ll)) <= M2[0]
            if not ok:
                bad += 1; print("MISSED lin", f, mpmath.nstr(lam, 8), sgw, mpmath.nstr(mu, 10), sg[0], m2[0], M2[0], mpmath.nstr(ll, 12))
print(f"tested {tested}, constraint active {active}, missed {bad}")
