import sys, numpy as np, re
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import paths, os
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qx, ivl
from cells_x import window
qx.DEBUG = True
q = qx.QX(8, 1, os.path.join(paths.SRC, "res_8_1.txt"), 5, 5)
wlo, whi = window("pi", 5.279748721052567e-09)
# s + l in W with s = 1 (first digit, in order): l in W - 1
import math
e1 = q.V((8,))
for i in range(len(e1)):
    for (plo, phi) in ivl.pre(0, e1.A[i:i+1], e1.B[i:i+1], e1.Z[i:i+1], e1.O[i:i+1], np.array([wlo]), np.array([whi])):
        print("s", e1.A[i], "pre window", plo, phi)
        lo, hi = plo, phi
H = q.Q((1, 2, 3, 4, 5, 6, 7), "N", lo, hi)
print("hits", len(H))
w = H.B - H.A
if len(H): print("width quantiles", np.quantile(np.where(np.isfinite(w), w, 1e300), [0.1, 0.5, 0.9]))
print("A=-inf:", np.isinf(H.A).sum(), "B=inf:", np.isinf(H.B).sum())
import io, contextlib
cats = {}
rng = np.random.default_rng(0)
for i in rng.choice(len(H), 40, replace=False):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        qx.explain(q, int(H.pid[i]))
    t = buf.getvalue().strip().split("\n")
    key = " / ".join(re.sub(r"\[[^\]]*\]", "[]", l.strip())[:60] for l in t[-3:])
    cats[key] = cats.get(key, 0) + 1
for k, v in sorted(cats.items(), key=lambda t: -t[1])[:10]:
    print(v, k)
