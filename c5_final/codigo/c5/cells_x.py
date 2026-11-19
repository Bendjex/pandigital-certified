"""C5 check of cells by forward enumeration (motor X): hits of NE(S) in the window of radius E*"""
import sys, os, json, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import paths; sys.path.insert(0, paths.SRC)
import motorx, ivl, xl
from flint import arb, ctx
ctx.prec = 200
DELTA = "4.669201609102990671853203820466"

def const(c):
    return {"pi": arb.pi(), "e": arb(1).exp(), "gamma": arb.const_euler(), "delta": arb(DELTA) + arb(0, 1e-30)}[c]

def window(c, E):
    x = const(c)
    lo = float((x - arb(E)).lower()); hi = float((x + arb(E)).upper())
    lo = np.nextafter(lo, 0); hi = np.nextafter(hi, np.inf)
    return float(ivl.ln_lo(np.array([lo]))[0]), float(ivl.ln_hi(np.array([hi]))[0])

def _main():
    pass
cert = {}
for l in open(os.path.join(paths.RESULTS, "certified.jsonl")):
    r = json.loads(l); cert[(r["const"], r["n"], r["ord"])] = r
if __name__ == "__main__":
    n, o = int(sys.argv[1]), int(sys.argv[2])
    consts = sys.argv[3].split(",") if len(sys.argv) > 3 else ["pi", "gamma", "delta", "e"]
    mx = motorx.MX(n, o, os.path.join(paths.SRC, f"res_{n}_{o}.txt"), min(n - 1, 5))
    S = tuple(range(1, n + 1))
    out = []
    for c in consts:
        r = cert.get((c, n, o))
        if r is None and (c, n, o) == ("e", 7, 1):
            E = 2.165605521030721e-08
        elif r is None:
            continue
        else:
            E = float(arb(r["optimal"][0]["error"]).upper())
        wlo, whi = window(c, E * (1 + 1e-9))
        t0 = time.time()
        h = mx.cell_hits(S, wlo, whi)
        res = {"const": c, "n": n, "ord": o, "E": E, "hits": len(h), "seconds": round(time.time() - t0, 1)}
        if len(h):
            res["hit_enclosures"] = [[float(a), float(b), p.tolist()] for a, b, p in zip(h.A[:20], h.B[:20], h.prov[:20])]
        print(json.dumps(res), flush=True)
        out.append(res)
    with open(os.path.join(paths.C5, "cells_x.jsonl"), "a") as f:
        for r in out:
            f.write(json.dumps(r) + "\n")
    