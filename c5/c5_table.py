"""Table 6 of the article: C5 results per cell (last run of each cell in cells_q.jsonl)"""
import json, collections
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import paths, os
rows = collections.OrderedDict()
for l in open(os.path.join(paths.C5, "cells_q.jsonl")):
    r = json.loads(l); rows[(r["ord"], r["n"], r["const"])] = r
def fmt(x):
    m, e = f"{x:.2e}".split("e"); return f"{m}·10^{int(e)}"
print("| mode | n | constant | window radius E* | hits | time (s) |")
print("|---|---|---|---|---|---|")
for (o, n, c), r in sorted(rows.items(), key=lambda t: (-t[0][0], t[0][1], ["pi", "gamma", "delta", "e"].index(t[0][2]))):
    print(f"| {'in order' if o else 'any order'} | {n} | {c} | {fmt(r['E'])} | {r['hits']} | {r['seconds']} |")
