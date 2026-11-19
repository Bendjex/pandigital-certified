import sys, os
PC = sys.argv[1]
sys.path.insert(0, os.path.join(PC, "c5")); sys.path.insert(1, os.path.join(PC, "src"))
import numpy as np, mpmath
from mpmath import mp, mpf
import xl
mp.prec = 300
rng = np.random.default_rng(7)
N = 200000
z = np.concatenate([rng.normal(0, 1, N // 4), rng.uniform(-745, 710, N // 4), rng.normal(0, 1e-8, N // 4) , np.exp(rng.uniform(-700, 6.5, N // 4)) * rng.choice([-1, 1], N // 4)])
elo = xl.exp_lo(z); ehi = xl.exp_hi(z)
bad = 0
for i in range(0, len(z), 7):
    t = mp.exp(mpf(float(z[i])))
    if not (mpf(float(elo[i])) <= t <= (mpf(float(ehi[i])) if np.isfinite(ehi[i]) else mp.inf)): bad += 1; print("EXP BAD", z[i], elo[i], ehi[i])
x = np.concatenate([np.exp(rng.uniform(-700, 700, N // 2)), 1 + rng.normal(0, 1e-9, N // 4), rng.uniform(0.5, 2, N // 4)])
x = x[x > 0]
llo = xl.ln_lo(x); lhi = xl.ln_hi(x)
worst = 0
for i in range(0, len(x), 7):
    t = mp.log(mpf(float(x[i])))
    if not (mpf(float(llo[i])) <= t <= mpf(float(lhi[i]))): bad += 1; print("LN BAD", x[i], llo[i], lhi[i])
    w = float((mpf(float(lhi[i])) - mpf(float(llo[i]))) / max(abs(t), mpf(10) ** -300))
    worst = max(worst, w)
print("exp/ln bounds: checked", len(range(0, len(z), 7)) + len(range(0, len(x), 7)), "bad", bad, " worst relative width of [ln_lo, ln_hi]:", worst)
print("lambda1 < log2(e) < lambda2:", mpf(float.fromhex("0x1.71547652b82fep+0")) < 1 / mp.log(2) < mpf(float.fromhex("0x1.71547652b82ffp+0")))
# special values
for v in (0.0, 5e-324, 2.2250738585072014e-308, 1.0, np.nextafter(1.0, 0), np.nextafter(1.0, 2), 1.7976931348623157e308, np.inf):
    a = np.array([v])
    print("  x =", v, " ln_lo, ln_hi =", xl.ln_lo(a)[0] if v > 0 else None, xl.ln_hi(a)[0] if v > 0 else None)
for v in (-np.inf, -1e308, -745.2, -744.0, 0.0, 709.78, 709.79, 1e308, np.inf):
    a = np.array([v]); print("  z =", v, " exp_lo, exp_hi =", xl.exp_lo(a)[0], xl.exp_hi(a)[0])
