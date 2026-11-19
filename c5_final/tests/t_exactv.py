"""exact_V (motor X recursion with exact leaves) against motor B's query: for random ln-windows on sets above
the table size, every value motor B finds certainly inside the window must also be found by exact_V
(numerically equal: overlapping 256-bit balls).  usage: t_exactv.py n ord nwin [seed]"""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, qx, ivl
from flint import arb, ctx
n, o, nwin = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
rng = np.random.default_rng(seed)
tv = 5 if o else 4
q = qx.QX(n, o, os.path.join(paths.SRC, f"res_{n}_{o}.txt"), tv, 5)
sets = [S for S in q.sets(tv + 1)]
miss = tot = 0; tB = tX = 0.0
for k in range(nwin):
    S = sets[rng.integers(len(sets))]
    kind = rng.integers(0, 4)
    if kind == 0:   c = rng.normal(0, 3); w = 10.0 ** rng.uniform(-12, -3)
    elif kind == 1: c = float(rng.choice([-1, 1]) * 10.0 ** rng.uniform(-12, -2)); w = abs(c) * 10.0 ** rng.uniform(-6, -1)   # near 1
    elif kind == 2: c = float(rng.uniform(-690, 690)); w = 10.0 ** rng.uniform(-9, -2)                         # extreme
    else:           c = float(np.log(rng.choice([0.5, 2, 3, 0.25, 1.5, 10, 0.1]))); w = 1e-10                   # simple values
    A, B = c - w, c + w
    t0 = time.time()
    lo = float(ivl.exp_lo(np.array([A]))[0]); hi = float(ivl.exp_hi(np.array([B]))[0])
    mb = list(q.B().query(tuple(str(d) for d in S), arb(np.nextafter(lo, 0)), arb(np.nextafter(hi, np.inf))).values())
    t1 = time.time()
    xv = q.exact_V(S, A, B)
    t2 = time.time(); tB += t1 - t0; tX += t2 - t1
    old = ctx.prec; ctx.prec = 256
    xb = [v.ball(256) for v, _ in xv]
    inner_lo, inner_hi = arb(hi) * 0 + arb(lo), arb(hi)
    for v, e in mb:
        b = v.ball(256)
        if not (b > arb(lo) and b < arb(hi)):
            continue                                   # not certainly inside: motor B keeps it as a superset
        tot += 1
        if not any(b.overlaps(x) for x in xb):
            miss += 1; print("MISSED by exact_V:", S, e, b, "window", A, B)
    ctx.prec = old
    if k % 10 == 0:
        print(f"win {k}: set {''.join(map(str, S))} kind {kind} motorB {len(mb)} exactV {len(xv)}  (B {t1 - t0:.2f}s, X {t2 - t1:.2f}s)", flush=True)
print(f"n={n} ord={o}: windows {nwin}, motor-B values certainly inside {tot}, MISSED {miss}; time motor B {tB:.1f}s, exact_V {tX:.1f}s")
