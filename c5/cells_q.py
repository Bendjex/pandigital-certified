"""C5 check of cells with recursive window queries (motor X): hits of NE({1..n}) in W(E*)"""
import sys, json, time, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import paths, os
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qx, ivl
from cells_x import window, const
cert = {}
for l in open(os.path.join(paths.RESULTS, "certified.jsonl")):
    r = json.loads(l); cert[(r["const"], r["n"], r["ord"])] = r
cert[("e", 7, 1)] = {"optimal": [{"error": "[2.165605521030721023076466e-8 +/- 3.08e-33]"}]}
n, o = int(sys.argv[1]), int(sys.argv[2])
tv, tn = int(sys.argv[3]), int(sys.argv[4])
consts = sys.argv[5].split(",") if len(sys.argv) > 5 else ["pi", "gamma", "delta", "e"]
from flint import arb
q = qx.QX(n, o, os.path.join(paths.SRC, f"res_{n}_{o}.txt"), tv, tn)
S = tuple(range(1, n + 1))
for c in consts:
    r = cert.get((c, n, o))
    if r is None:
        continue
    E = float(arb(r["optimal"][0]["error"]).upper())
    wlo, whi = window(c, E * (1 + 1e-9))
    t0 = time.time(); q.stats.clear(); q.top = S
    H = q.Q(S, "N", np.array([wlo]), np.array([whi]))
    res = {"const": c, "n": n, "ord": o, "E": E, "hits": len(H), "seconds": round(time.time() - t0, 1),
           "stats": {f"{k[0]}{k[1]}": [v[0], v[1], round(v[2], 1)] for k, v in q.stats.items()}}
    if len(H):
        res["sample"] = [[float(a), float(b)] for a, b in zip(H.A[:10], H.B[:10])]
    print(json.dumps(res), flush=True)
    with open(os.path.join(paths.C5, "cells_q.jsonl"), "a") as f:
        f.write(json.dumps(res) + "\n")
