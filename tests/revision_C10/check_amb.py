"""¿Se dispara refine_near1 (clases con 0 en su enclosure de ln que no son el 1 exacto) en alguna tabla V?
uso: check_amb.py PC n ord tv RES"""
import sys, os, itertools, time
PC = sys.argv[1]
sys.path.insert(0, os.path.join(PC, "c5")); sys.path.insert(1, os.path.join(PC, "src"))
import numpy as np
import xl, motorx
n, o, tv = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
t0 = time.time()
xl.init(n, o, sys.argv[5], tv)
D = tuple(range(1, n + 1))
tot = amb_tot = sets = 0; worst = []
for k in range(1, tv + 1):
    Ss = [D[i:i + k] for i in range(n - k + 1)] if o else list(itertools.combinations(D, k))
    for S in Ss:
        A, B, fl = xl.table(motorx.bits(S))
        Z = (fl & 1) > 0; O = (fl & 2) > 0
        amb = (A <= 0) & (B >= 0) & ~O & ~Z
        tot += len(A); amb_tot += int(amb.sum()); sets += 1
        if amb.any(): worst.append((S, int(amb.sum()), [(float(a), float(b)) for a, b in zip(A[amb][:3], B[amb][:3])]))
print(f"n={n} ord={o} tv={tv}: {sets} conjuntos, {tot} clases; clases ambiguas (0 en [A,B], ni 1 ni 0 exactos): {amb_tot}  ({time.time() - t0:.0f}s)")
for w in worst[:10]: print("   ", w)
