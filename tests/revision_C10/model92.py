import numpy as np
from math import comb, factorial, log, exp
C = {"pi": np.pi, "gamma": 0.5772156649015329, "delta": 4.669201609102991}
DREL = 0.05
DW = {0: {2: (1, 0, 0), 3: (21, 14, 6), 4: (900, 669, 505), 5: (46260, 38739, 34511), 6: (4234516, 4175942, 3333557)},
      1: {2: (1, 0, 0), 3: (15, 8, 5), 4: (86, 70, 74), 5: (1492, 1615, 2017), 6: (37805, 60643, 32699), 7: (909808, 946739, 819132)}}
E = {0: {"pi": [1.4159e-1, 8.2593e-3, 3.2618e-5, 5.9270e-6, 1.3801e-8, 3.0184e-10, 1.5987e-12, 2.81216e-14],
         "gamma": [5.3742e-2, 1.3460e-4, 1.0408e-5, 2.7235e-7, 8.3812e-9, 2.5780e-11, 6.0839e-13, 1.13861e-15],
         "delta": [3.3080e-1, 2.7613e-2, 2.5349e-3, 2.9461e-5, 1.2702e-7, 3.0581e-11, 2.2255e-12, 1.37758e-14]},
     1: {"pi": [1.4159e-1, 2.0685e-2, 7.1145e-4, 1.3961e-4, 3.8200e-6, 1.8892e-7, 5.2797e-9, 2.3672e-10],
         "gamma": [5.3742e-2, 1.4880e-3, 6.4537e-4, 8.6432e-5, 2.7235e-7, 7.1303e-8, 1.4779e-9, 5.8208e-11],
         "delta": [3.3080e-1, 2.7613e-2, 2.5349e-3, 1.1726e-4, 8.1961e-6, 6.9319e-8, 1.2724e-8, 1.1591e-9]}}
def nexpr(n, ordered):
    A = {}
    for k in range(1, n + 1):
        lit = (k + 1) if ordered else factorial(k) * (k + 1)
        s = sum((1 if ordered else comb(k, i)) * 5 * (2 * A[i]) * (2 * A[k - i]) for i in range(1, k))
        A[k] = lit + s
    return {k: 2 * A[k] for k in A}
def model(o, ci):
    N = nexpr(9, o); r = {k: N[k] / N[k - 1] for k in range(2, 10)}
    d = {k: v[ci] for k, v in DW[o].items()}
    n0 = max(d); g = d[n0] / d[n0 - 1]; kappa = g / r[n0]
    dd = dict(d)
    for m in range(n0 + 1, 10):
        dd[m] = dd[m - 1] * kappa * r[m]
    return n0, dd, N
names = ["pi", "gamma", "delta"]
rows = []
for o in (0, 1):
    for ci, c in enumerate(names):
        n0, dd, N = model(o, ci)
        Delta = DREL * C[c]
        for n in range(3, 10):
            em = Delta / (dd[n] + 1); ratio = E[o][c][n - 2] / em
            rows.append((o, c, n, n <= n0, E[o][c][n - 2], em, ratio))
print("N_n any order:", nexpr(9, 0)); print("N_n in order :", nexpr(9, 1))
inb = [r for r in rows if 0.051 <= r[6] <= 3.0]
print("cells:", len(rows), " in band:", len(inb), " expected:", 0.9 * len(rows))
print("outside:", [(("any" if r[0] == 0 else "ord"), r[1], r[2], round(r[6], 2), "measured" if r[3] else "extrapolated") for r in rows if not (0.051 <= r[6] <= 3.0)])
gm = exp(np.mean([log(r[6]) for r in rows])); print("geometric mean of ratio (all 42):", round(gm, 3), "; exponential law predicts exp(-gamma) =", round(exp(-0.5772156649), 3))
meas = [r for r in rows if r[3]]; ext = [r for r in rows if not r[3]]
print("measured cells:", len(meas), "in band", sum(0.051 <= r[6] <= 3.0 for r in meas), "GM", round(exp(np.mean([log(r[6]) for r in meas])), 3),
      "| extrapolated:", len(ext), "in band", sum(0.051 <= r[6] <= 3.0 for r in ext), "GM", round(exp(np.mean([log(r[6]) for r in ext])), 3))
ml = np.mean([log(r[6]) for r in rows]); se = (np.pi / np.sqrt(6)) / np.sqrt(len(rows))
print("mean ln ratio:", round(ml, 3), " expected -0.577, s.e. (if independent)", round(se, 3), " z =", round((ml + 0.5772) / se, 2))
print("above 1:", sum(r[6] > 1 for r in rows), "of", len(rows), "(exponential: P(X>1) = 0.368 -> expected", round(0.368 * len(rows), 1), ")")
for o in (0, 1):
    for c in names:
        print(("any " if o == 0 else "ord "), c, " ".join(f"n{r[2]}:{r[6]:.2f}{'' if r[3] else '*'}" for r in rows if r[0] == o and r[1] == c))
print("model n=9 any order:", {r[1]: f"{r[5]:.2e}" for r in rows if r[0] == 0 and r[2] == 9})
