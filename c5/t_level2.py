"""soundness of qx.pre_level2 (second-level constraint on the bases l of l^s and l^-s): random s (ln s in an
enclosure [A, B]: moderate, tiny, huge and double-exponential), random l (moderate, near 1 on both sides, huge and
tiny) and windows [al, be] for ln f(s, l) = +-s ln l.  Whenever f(s, l) lies in the window (exactly, in mpmath)
and the constraint is active (sign != 0), ln l must have that sign and ln|ln l| must lie in [m, M].
usage: t_level2.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ivl, qx
mpmath.mp.dps = 60
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 41)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
bad = tested = inside = active = 0
for t in range(N):
    # s = e^lam, lam in [A, B]
    kind = rng.integers(0, 4)
    if kind == 0:
        lam = mpmath.mpf(float(rng.normal(0, 3)))
    elif kind == 1:
        lam = mpmath.mpf(float(rng.choice([-1, 1]) * 10 ** rng.uniform(0, 2.8)))
    elif kind == 2:
        lam = mpmath.mpf(float(rng.choice([-1, 1]) * 10 ** rng.uniform(2.8, 15)))
    else:
        lam = mpmath.mpf(float(rng.choice([-1, 1]) * 10 ** rng.uniform(-12, 0)))
    w = abs(float(lam)) * float(10 ** rng.uniform(-16, -2)) if rng.random() < 0.8 else float(10 ** rng.uniform(-3, 1))
    A = np.nextafter(float(lam) - w * rng.random(), -np.inf); B = np.nextafter(float(lam) + w * rng.random(), np.inf)
    lam_true = lam if A <= lam <= B else mpmath.mpf(A + (B - A) * rng.random())
    x = ivl.Q7(np.array([A]), np.array([B]), np.array([False]), np.array([False]))
    f = int(rng.choice([7, 8]))
    s = mpmath.exp(lam_true)
    # l with ln l = mu
    k2 = rng.integers(0, 4)
    if k2 == 0:
        mu = mpmath.mpf(float(rng.normal(0, 3)))
    elif k2 == 1:
        mu = mpmath.mpf(float(rng.choice([-1, 1]) * 10 ** rng.uniform(-300, -1)))
    elif k2 == 2:
        mu = mpmath.mpf(float(rng.choice([-1, 1]) * 10 ** rng.uniform(1, 300)))
    else:                                   # aimed at a window around a target value of ln f
        target = float(rng.normal(0, 5))
        mu = mpmath.mpf(target) / ((1 if f == 7 else -1) * s) * (1 + mpmath.mpf(float(rng.normal(0, 1e-9))))
    val = (1 if f == 7 else -1) * s * mu      # ln f(s, l)
    c = float(val) if rng.random() < 0.8 else float(rng.normal(0, 5))
    ww = float(10 ** rng.uniform(-14, 1))
    al, be = c - ww * rng.random(), c + ww * rng.random()
    if rng.random() < 0.05:
        al = -np.inf
    tested += 1
    if not (al <= val <= be):
        continue
    inside += 1
    sg, m2, M2 = qx.pre_level2(f, x, np.array([al]), np.array([be]))
    if sg[0] == 0:
        continue
    active += 1
    ok = (mu > 0 if sg[0] == 1 else mu < 0)
    if ok:
        ll = mpmath.log(abs(mu))
        ok = m2[0] <= ll <= M2[0]
    if not ok:
        bad += 1
        print("MISSED", f, "lnS", A, B, "ln l", mpmath.nstr(mu, 12), "window", al, be, "constraint", sg[0], m2[0], M2[0])
print(f"tested {tested}, with f(s, l) in the window {inside}, constraint active {active}, missed {bad}")
