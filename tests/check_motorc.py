"""every candidate path that motor B finds non-empty must be kept by motor C. usage: n ord window cands kept"""
import sys
sys.path.insert(0, ".")
from flint import arb
import verify
n, o, wf, cf, kf = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
a, b = [float.fromhex(x) for x in open(wf).read().split()]
kept = set(l.strip() for l in open(kf) if l.strip())
V = verify.Verifier(n, o, fullmax=4)
bad = nonempty = 0
for l in open(cf):
    l = l.strip()
    if not l: continue
    vals = V.verify_line(l.split(" ", 1)[1], arb(a), arb(b))
    if vals:
        nonempty += 1
        if l not in kept:
            bad += 1; print("DROPPED BUT NON-EMPTY:", l)
print("non-empty paths", nonempty, "dropped wrongly", bad)
