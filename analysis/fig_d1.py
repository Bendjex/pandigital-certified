"""D1: certified optimal error vs n against the random (uniform) model. Data: dp.py window counts."""
import numpy as np, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import comb, factorial

C = {"pi": np.pi, "gamma": 0.5772156649015329, "delta": 4.669201609102991}
DREL = 0.05
# distinct values within ±5% of c (dp.py; exact double-precision DP, merged at ~1e-12 relative)
DW = {0: {2: (1, 0, 0), 3: (21, 14, 6), 4: (900, 669, 505), 5: (46260, 38739, 34511), 6: (4234516, 4175942, 3333557)},
      1: {2: (1, 0, 0), 3: (15, 8, 5), 4: (86, 70, 74), 5: (1492, 1615, 2017), 6: (37805, 60643, 32699), 7: (909808, 946739, 819132)}}
# certified optima (Tables 2-3); any order n = 9 certified by C7 on 05/11/2026 (laboratory)
E = {0: {"pi": [1.4159e-1, 8.2593e-3, 3.2618e-5, 5.9270e-6, 1.3801e-8, 3.0184e-10, 1.5987e-12, 2.8122e-14],
         "gamma": [5.3742e-2, 1.3460e-4, 1.0408e-5, 2.7235e-7, 8.3812e-9, 2.5780e-11, 6.0839e-13, 1.1386e-15],
         "delta": [3.3080e-1, 2.7613e-2, 2.5349e-3, 2.9461e-5, 1.2702e-7, 3.0581e-11, 2.2255e-12, 1.3776e-14]},
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
    g1 = d[n0 - 1] / d[n0 - 2]; kappa1 = g1 / r[n0 - 1]
    dd, lo, hi = dict(d), dict(d), dict(d)
    for m in range(n0 + 1, 10):
        dd[m] = dd[m - 1] * kappa * r[m]
        alt = (lo[m - 1] if kappa1 < kappa else hi[m - 1]) * kappa1 * r[m]
        lo[m] = min(dd[m], alt) if kappa1 < kappa else dd[m]
        hi[m] = max(dd[m], alt) if kappa1 > kappa else dd[m]
        if kappa1 < kappa: hi[m] = dd[m]; lo[m] = lo[m - 1] * kappa1 * r[m]
        else: lo[m] = dd[m]; hi[m] = hi[m - 1] * kappa1 * r[m]
    return n0, dd, lo, hi

names = ["pi", "gamma", "delta"]; labels = {"pi": "π", "gamma": "γ", "delta": "δ"}
col = {"pi": "#2a78d6", "gamma": "#eb6834", "delta": "#1baf7a"}
mk = {"pi": "o", "gamma": "s", "delta": "^"}
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#888", "xtick.color": "#444", "ytick.color": "#444"})
fig, axs = plt.subplots(2, 2, figsize=(7.0, 5.2), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
out = {}
for j, (o, title) in enumerate([(0, "Digits 1–n in any order"), (1, "Digits 1–n in increasing order")]):
    ax, axr = axs[0, j], axs[1, j]
    for ci, c in enumerate(names):
        n0, dd, lo, hi = model(o, ci)
        Delta = DREL * C[c]
        ns = [n for n in range(3, 10)]
        em = {n: Delta / (dd[n] + 1) for n in ns}
        meas = [n for n in ns if n <= n0]; ext = [n for n in ns if n >= n0]
        ax.plot(meas, [em[n] for n in meas], "-", color=col[c], lw=1.2, alpha=0.9)
        ax.plot(ext, [em[n] for n in ext], "--", color=col[c], lw=1.2, alpha=0.9)
        ax.fill_between(ext, [Delta / (hi[n] + 1) for n in ext], [Delta / (lo[n] + 1) for n in ext], color=col[c], alpha=0.12, lw=0)
        errs = E[o][c]; nn = list(range(2, 10))
        pend = False
        ax.plot(nn[:-1] if pend else nn, errs[:-1] if pend else errs, mk[c], color=col[c], ms=5, label=labels[c], mec="white", mew=0.6)
        if pend: ax.plot([9], [errs[-1]], mk[c], color=col[c], ms=5, mfc="white", mew=1.2)
        rat = [errs[n - 2] / em[n] for n in ns]
        axr.plot(ns, rat, mk[c] + "-", color=col[c], ms=4, lw=0.8, mec="white", mew=0.5)
        out[f"{o}_{c}"] = {n: (errs[n - 2], em[n]) for n in ns}
    ax.set_yscale("log"); ax.set_title(title, fontsize=9.5, loc="left")
    ax.grid(axis="y", color="#e6e6e6", lw=0.6); ax.set_ylim(1e-16, 1)
    axr.set_yscale("log"); axr.axhspan(0.051, 3.0, color="#999", alpha=0.15, lw=0)
    axr.axhline(1, color="#888", lw=0.6); axr.set_ylim(0.01, 30); axr.set_xlabel("n")
    axr.grid(axis="y", color="#e6e6e6", lw=0.6)
axs[0, 0].set_ylabel("|error|"); axs[1, 0].set_ylabel("optimum / model")
axs[0, 0].legend(frameon=False, loc="lower left", title=None, handletextpad=0.3)
axs[1, 1].text(2.0, 0.07, "90% band of the model", fontsize=7, color="#666", va="bottom")
fig.tight_layout()
fig.savefig("./fig_error_vs_n.pdf"); fig.savefig("./fig_error_vs_n.png", dpi=180)
json.dump(out, open("./model_vs_opt.json", "w"), indent=0)
for k, v in out.items(): print(k, {n: f"{a:.2e}/{b:.2e}={a/b:.2f}" for n, (a, b) in v.items()})
