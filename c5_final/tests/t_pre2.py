"""soundness of qx.pre2 (preimage windows of s^l and s^-l intersected with the window from the second level of s):
random s given by ln|ln s| (moderate, near 1 and double-exponential bases, both signs of ln s), random l and
windows; whenever f(s, l) lies in the window (exactly, in mpmath), l must lie in one of the windows of pre2.
usage: t_pre2.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ivl, qx
mpmath.mp.dps = 40
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 31)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
bad = tested = inside = 0
for t in range(N):
    sg = int(rng.choice([-1, 1]))
    kind = rng.integers(0, 3)
    mu = mpmath.mpf(float(rng.uniform(-30, 5))) if kind == 0 else (mpmath.mpf(float(rng.uniform(5, 700))) if kind == 1 else mpmath.mpf(float(rng.uniform(700, 900))))
    lam = sg * mpmath.exp(mu)                        # ln s
    m_lo, m_hi = np.nextafter(float(mu), -np.inf), np.nextafter(float(mu), np.inf)
    a = float(lam)
    A = np.nextafter(a, -np.inf) if np.isfinite(a) else (-np.inf if a < 0 else np.nextafter(1.7976931348623157e308, -np.inf))
    B = np.nextafter(a, np.inf) if np.isfinite(a) else (np.inf if a > 0 else -np.nextafter(1.7976931348623157e308, -np.inf))
    x = ivl.Q7(np.array([A]), np.array([B]), np.array([False]), np.array([False]),
               np.array([sg], np.int8), np.array([m_lo]), np.array([m_hi]))
    f = int(rng.choice([5, 6]))
    # an l whose value lands in a window around a random target, or a random l
    target = mpmath.mpf(float(rng.normal(0, 5))) if rng.random() < 0.7 else mpmath.mpf(float(rng.choice([-1, 1]) * 10 ** rng.uniform(-3, 4)))
    if rng.random() < 0.6:
        ll = target / ((1 if f == 5 else -1) * lam) * (1 + mpmath.mpf(float(rng.normal(0, 1e-9))))
        if ll <= 0:
            continue
        lnl = mpmath.log(ll)
    else:
        lnl = mpmath.mpf(float(rng.uniform(-50, 50)))
    val = (1 if f == 5 else -1) * mpmath.exp(lnl) * lam        # ln f(s, l)
    w = 10.0 ** rng.uniform(-12, 0)
    al, be = float(target) - w, float(target) + w
    tested += 1
    if not (al <= val <= be):
        continue
    inside += 1
    wins = qx.pre2(f, x, np.array([al]), np.array([be]))
    ok = any(lo[0] <= lnl <= hi[0] for lo, hi in wins)
    if not ok:
        bad += 1
        print("MISSED", f, sg, mpmath.nstr(mu, 8), "ln l", mpmath.nstr(lnl, 12), "window", al, be, [(lo[0], hi[0]) for lo, hi in wins])
print(f"tested {tested}, with f(s, l) in the window {inside}, missed {bad}")
