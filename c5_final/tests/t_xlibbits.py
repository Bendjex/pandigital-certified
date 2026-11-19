"""bit identity of xlib (xlib.so / xlib.dll): rigorous exp/ln bounds on 800 000 deterministic inputs and motor A's
tables of n = 5 (any order) dumped as enclosures; prints hashes to compare across machines.  usage: t_xlibbits.py res_5_0.txt"""
import sys, os, hashlib, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xl
rng = np.random.default_rng(20261103)
m = rng.random(200000); e = rng.integers(-1000, 1000, 200000)
v = m * np.sign(e + 0.5) * np.ldexp(1.0, (np.abs(e) % 60) * np.sign(e))
for name in ("exp_lo", "exp_hi", "ln_lo", "ln_hi"):
    x = np.clip(v, -740, 700) if name.startswith("exp") else np.abs(v) + 1e-300
    y = getattr(xl, name)(x)
    print(name, hashlib.sha256(y.tobytes()).hexdigest()[:16])
xl.init(5, 0, sys.argv[1], 4)
h = hashlib.sha256(); tot = 0
for mask in range(1, 32):
    if bin(mask).count("1") > 4:
        continue
    lo, hi, fl = xl.table(mask)
    h.update(lo.tobytes()); h.update(hi.tobytes()); h.update(fl.tobytes()); tot += len(lo)
print("tables_n5", tot, h.hexdigest()[:16])
