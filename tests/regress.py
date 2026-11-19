"""compare a results file with a reference results file: usage regress.py new.jsonl [reference.jsonl]"""
import json, sys
ref = {}
for l in open(sys.argv[2] if len(sys.argv) > 2 else "results/certified.jsonl"):
    r = json.loads(l); ref[(r["const"], r["n"], r["ord"])] = r
bad = 0; n = 0
for l in open(sys.argv[1]):
    r = json.loads(l); k = (r["const"], r["n"], r["ord"]); n += 1
    e1 = r["optimal"][0]["error"].split()[0]
    ok = all(c["valid_in_range"] and c["digits_ok"] for g in r["optimal"] for c in g["checked"])
    if k in ref:
        e0 = ref[k]["optimal"][0]["error"].split()[0]
        if e0 != e1: bad += 1; print("DIFFERENT", k, e0, e1)
    if not ok or r["unresolved"] or r["ties"]:
        bad += 1; print("PROBLEM", k, ok, r["unresolved"], r["ties"])
print(n, "cells,", bad, "problems")
