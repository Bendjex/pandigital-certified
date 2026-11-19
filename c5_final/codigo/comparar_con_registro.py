"""Compare a split computed again (the .jsonl written by c5run.py) with the record used in the closure.

usage:  python comparar_con_registro.py <new.jsonl> <resultados_n9/registros_usados_n9.jsonl>

For every (constant, split) and every piece of the new file: the hits must be the same as in the record used, and the
new run must leave 0 values undecided (C7 and C11), 0 values certainly inside the C11 window and 0 undecided
operations of motor B. The fingerprint of the checker of each side is printed (field "code"; records made before
13 November 2026 have none: checker of the first version, c90dbdc006cb5be9, or f70deff8e1e04e0e, which differs only
in printing).
"""
import json, sys


def load(path, wrapped):
    sp, ch = {}, {}
    recs = []
    for li, l in enumerate(open(path, encoding="utf-8")):
        try:
            r = json.loads(l)
        except ValueError:
            continue
        if wrapped and "registro" in r:
            r = r["registro"]
        if isinstance(r, dict) and "const" in r and "split" in r:
            recs.append((float(r.get("t", 0)), li, r))
    recs.sort(key=lambda x: x[:2])
    for _, _, r in recs:
        if "chunk" in r:
            ch[(r["const"], int(r["split"]), r["chunk"])] = r
        elif "S1" in r:
            sp[(r["const"], int(r["split"]))] = r
    return sp, ch


def main():
    nsp, nch = load(sys.argv[1], False)
    osp, och = load(sys.argv[2], True)
    bad = 0
    for k, r in sorted(nsp.items()):
        o = osp.get(k)
        und = (r.get("undecided") or 0) + (r.get("undecided_c11") or 0) + (r.get("inside_c11") or 0) + \
            (r.get("motorB_unresolved") or 0)
        same = o is not None and o["hits"] == r["hits"]
        bad += (not same) + (und > 0)
        print(f"{k[0]:6s} split {k[1]:3d}: hits {r['hits']} (record {o['hits'] if o else '-'}) "
              f"{'SAME' if same else 'DIFFERENT'}; undecided {und}; checker {r.get('code')} (record {o.get('code') if o else '-'})")
    pieces = [(k, r) for k, r in nch.items() if k in och]
    dp = [k for k, r in pieces if r["hits"] != och[k]["hits"]]
    bad += len(dp)
    print(f"pieces compared: {len(pieces)}; different hits: {len(dp)} {dp[:5]}")
    print("COMPARISON: " + ("OK" if not bad else "DIFFERENT"))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
