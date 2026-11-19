import sys, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import paths, os
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import motorx, ivl
mx = motorx.MX(5, 0, os.path.join(paths.SRC, "res_5_0.txt"), 4)
L = 18.82654565077563
for name, e in (("V", mx.V((1,3,4,5))), ("N", mx.N((1,3,4,5)))):
    m = np.where((e.A <= L + 1e-9) & (e.B >= L - 1e-9))[0]
    print(name, "entries containing y:", [(e.A[i], e.B[i], e.S[i], e.M[i], e.N[i]) for i in m])
    for i in m:
        x = mx.V((2,))
        r, ok = ivl.fwd(6, x.q7(np.array([1])), e.q7(np.array([i])))
        print("  fwd 2^-y:", r.A, r.B, ok, "possibly_out", ivl.possibly_out(r.A, r.B, r.Z))
N = mx.N((1,2,3,4,5)); t = -104013653.29732081
print("in N12345:", ((N.A <= t + 1) & (N.B >= t - 1)).sum())
