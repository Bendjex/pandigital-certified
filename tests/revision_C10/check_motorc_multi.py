"""Lema 15, comprobación empírica del revisor: todo camino candidato que motor B encuentra no vacío debe estar en la
salida de motor C, en cada configuración.  uso (desde src/): check_motorc_multi.py n ord window cands kept1 kept2 ..."""
import sys, time
sys.path.insert(0, ".")
from flint import arb
import verify
n, o, wf, cf = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
kfs = sys.argv[5:]
a, b = [float.fromhex(x) for x in open(wf).read().split()]
kept = [set(l.strip() for l in open(k) if l.strip()) for k in kfs]
V = verify.Verifier(n, o, fullmax=4)
cands = list(dict.fromkeys(l.strip() for l in open(cf) if l.strip()))
bad = [0] * len(kfs); nonempty = 0; t0 = time.time()
for l in cands:
    vals = V.verify_line(l.split(" ", 1)[1], arb(a), arb(b))
    if vals:
        nonempty += 1
        for i, ks in enumerate(kept):
            if l not in ks:
                bad[i] += 1; print("DROPPED BUT NON-EMPTY (%s):" % kfs[i], l)
print("distinct candidates", len(cands), "non-empty for motor B", nonempty, "motor B unresolved", len(V.B.unresolved), "time %.0fs" % (time.time() - t0))
for k, ks, bd in zip(kfs, kept, bad): print("   ", k, "kept", len(ks), "dropped wrongly", bd)
