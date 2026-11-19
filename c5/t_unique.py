"""qx.pre2 and qx.pre_level2 evaluate the logs of the window bounds once per distinct value (qx._on_unique):
check that the windows are bit-identical to the direct evaluation on random batches with repeated bounds
(as in combine: every window repeated once per table entry).  usage: t_unique.py [trials] [seed]"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ivl, qx
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
bad = 0
for t in range(N):
    nw, n1 = int(rng.integers(1, 40)), int(rng.integers(1, 300))
    c = rng.normal(0, 3, nw) * 10.0 ** rng.integers(-300, 300, nw) * rng.choice([0, 1], nw, p=[0.05, 0.95])
    w = np.abs(rng.normal(0, 1, nw)) * 10.0 ** rng.integers(-15, 3, nw)
    lo = c - w; hi = c + w
    lo[rng.random(nw) < 0.05] = -np.inf; hi[rng.random(nw) < 0.05] = np.inf
    W = np.repeat(np.arange(nw), n1)
    n = len(W)
    A = rng.normal(0, 50, n) * 10.0 ** rng.integers(-5, 300, n)
    B = A + np.abs(rng.normal(0, 1, n)) * 10.0 ** rng.integers(-16, 2, n)
    S = np.sign(A).astype(np.int8); S[(A < 0) & (B > 0)] = 0
    with np.errstate(all="ignore"):
        M = np.log(np.abs(A)); Nn = np.log(np.abs(B))
    M, Nn = np.minimum(M, Nn), np.maximum(M, Nn)
    x = ivl.Q7(A, B, rng.random(n) < 0.02, rng.random(n) < 0.02, S, M, Nn)
    for f in (5, 6):
        got = qx.pre2(f, x, lo[W], hi[W])
        qx.ln_lo_u, qx.ln_hi_u = ivl.ln_lo, ivl.ln_hi
        ref = qx.pre2(f, x, lo[W], hi[W])
        qx.ln_lo_u, qx.ln_hi_u = (lambda v: qx._on_unique(ivl.ln_lo, v)), (lambda v: qx._on_unique(ivl.ln_hi, v))
        for (g1, g2), (r1, r2) in zip(got, ref):
            if not (np.array_equal(g1, r1, equal_nan=True) and np.array_equal(g2, r2, equal_nan=True)):
                bad += 1; print("pre2 differs", f, t)
    for f in (7, 8):
        got = qx.pre_level2(f, x, lo[W], hi[W])
        qx.ln_lo_u, qx.ln_hi_u = ivl.ln_lo, ivl.ln_hi
        ref = qx.pre_level2(f, x, lo[W], hi[W])
        qx.ln_lo_u, qx.ln_hi_u = (lambda v: qx._on_unique(ivl.ln_lo, v)), (lambda v: qx._on_unique(ivl.ln_hi, v))
        for g, r in zip(got, ref):
            if not np.array_equal(g, r, equal_nan=True):
                bad += 1; print("pre_level2 differs", f, t)
print(f"trials {N}: bit-identical windows, differences {bad}")
