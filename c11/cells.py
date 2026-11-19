"""list of C11 cells: n ord const E_start (the certified with-minus optimum), cheapest first"""
import os
import json
cells = []
for l in open(os.environ.get("PANDI_RESULTS", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")) + "/certified.jsonl"):
    r = json.loads(l)
    E = float(r["optimal"][0]["error"].strip("[").split()[0])
    cells.append((r["n"], r["ord"], r["const"], E))
cells.append((7, 1, "e", 2.165605521030721e-08))
cost = lambda c: (c[0] if c[1] == 1 else c[0] + 2.5, c[2])
for c in sorted(cells, key=cost):
    print(*c)
