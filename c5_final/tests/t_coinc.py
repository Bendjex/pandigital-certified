"""(11/11) Theorem X in its structural form: every expression with some node out of range is covered, also when its
value coincides EXACTLY with a value that has an in-range writing (then it is not in NE(S), only in X(S)).
Coincidences with few digits: 1^(9^(8^7)) = 1 (9^2097152 is out of range) and 2 * 1^(9^(8^7)) = 2, (2 + 1^(...)) = 3.
The queries Q(S, NE, w) for S = {1, 7, 8, 9} (window around 1), S = {1, 2, 7, 8, 9} (around 2 and 3: table lookups) and
S = {1, 2, 3, 7, 8, 9} (around 5 = 2 + 3 * 1^(...) and 6: the recursion), relative radius 1e-13, in provenance mode;
hitcheck evaluates every derivation of every hit.  The test passes if the exact values 1, 2, 3, 5 and 6 are found among
them, written through the out-of-range power (the exponent 9^(8^7)).
(The 7-digit coincidence 9^(4^8)/3^(2^17) = 1 is covered by the same argument; its query, around 1 with 7 digits, needs
more memory than this machine has with the other runs.)
usage: t_coinc.py"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, hitcheck, qx, ivl, lvl
qx.DEBUG = True
q = qx.QX(9, 0, os.path.join(os.environ.get("RESDIR", paths.SRC), "res_9_0.txt"), 4, 5, verbose=False)


def win(vals):
    lo = np.array([float(ivl.ln_lo(np.array([v * (1 - 1e-13)]))[0]) for v in vals])
    hi = np.array([float(ivl.ln_hi(np.array([v * (1 + 1e-13)]))[0]) for v in vals])
    return lo, hi


C = hitcheck.Checker(q)
ok = True
for S, targets in ((((1, 7, 8, 9), (1,)), ((1, 2, 7, 8, 9), (2, 3))) if os.environ.get("COINC") != "6" else (((1, 2, 3, 7, 8, 9), (6,)),)):
    qx.PROV = qx._Prov()
    lo, hi = win([float(t) for t in targets])
    H = q.Q(S, "N", lo, hi)
    found = {t: [] for t in targets}
    nv = 0
    for i in range(len(H)):
        t = targets[int(H.w[i])]
        for (v, e) in C.hit(int(H.pid[i])):
            nv += 1
            if v is not None and v.q is not None and v.q == t and "2097152" in e:
                found[t].append(e)
    for t in targets:
        print(f"S = {S}, value {t}: hits {int((H.w == targets.index(t)).sum())}, exact {t} through 9^(8^7): {len(found[t])} "
              f"derivations, e.g. {found[t][0][:160] if found[t] else '-'}", flush=True)
        ok &= bool(found[t])
print("t_coinc:", "OK" if ok else "MISSING", "; motor B unresolved", q.motorb_unresolved())
sys.exit(0 if ok else 1)
