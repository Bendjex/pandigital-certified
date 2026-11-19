"""C5 bit-for-bit test: hits of an artificial window around 1 (many near-1 NE values), n = 5 in order (table
lookup) and n = 6 in order (top-level split loop); prints a hash of every enclosure field, to compare across machines"""
import sys, os, hashlib, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, qx, ivl
rows = []
for n in (5, 6):
    q = qx.QX(n, 1, os.path.join(paths.SRC, f"res_{n}_1.txt"), 5, 5)
    lo, hi = 1 - 1e-9, 1 + 1e-9
    wl = np.array([float(ivl.ln_lo(np.array([lo]))[0])]); wh = np.array([float(ivl.ln_hi(np.array([hi]))[0])])
    H = q.Q(tuple(range(1, n + 1)), "N", wl, wh)
    for i in range(len(H)):
        rows.append(f"{n} " + " ".join(float(getattr(H, f)[i]).hex() for f in qx.Hits.F))
rows.sort()
print("BITS", len(rows), hashlib.sha256("\n".join(rows).encode()).hexdigest())
