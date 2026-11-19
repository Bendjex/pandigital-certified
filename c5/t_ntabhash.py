"""bit identity of motor X's N tables: N(S) for every set of 4 and 5 digits of {1..n} (any order), hash of every
field.  usage: t_ntabhash.py n ord"""
import sys, os, hashlib, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, motorx
n, o = int(sys.argv[1]), int(sys.argv[2])
q = motorx.MX(n, o, os.path.join(os.environ.get("RESDIR", paths.SRC), f"res_{n}_{o}.txt"), 5 if o else 4)
h = hashlib.sha256(); tot = 0
for k in (3, 4, 5):
    for S in q.sets(k):
        e = q.N(S)
        for f in motorx.Ent.F + ("prov",):
            h.update(np.ascontiguousarray(getattr(e, f)).tobytes())
        tot += len(e)
print("NTAB", tot, h.hexdigest())
