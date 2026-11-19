import sys, os
PC = sys.argv[1]
sys.path.insert(0, os.path.join(PC, "c5")); sys.path.insert(1, os.path.join(PC, "src")); sys.path.insert(2, os.path.join(PC, "c11"))
import numpy as np
import ivl
# (a) non-strict sign: s = exact 1, l = enclosure [-inf, -5] (possibly 0, not flagged Z)
one = ivl.Q7(np.array([0.0]), np.array([0.0]), np.array([False]), np.array([True]))
t = ivl.Q7(np.array([-np.inf]), np.array([-5.0]), np.array([False]), np.array([False]))
for f in (0, 1):
    r, ok = ivl.fwd(f, one, t)
    print(f"(a) fwd({ivl.FORMS[f]}; 1, t in [0, e^-5]):  S={int(r.S[0])} M={r.M[0]} N={r.N[0]:.4g}  A={r.A[0]:.4g} B={r.B[0]:.4g}  O={bool(r.O[0])}   <- u = 1 (t = 0) is in [A,B] but 'sign(ln u) = S' is false for it")
# (b) conversion of Q > 2^53
Q = 2 ** 53 + 1
qq = float(ivl.succ(np.array([float(Q) * float(Q)]))[0])
from fractions import Fraction as Fr
print("(b) Q1 = Q2 = 2^53+1:  succ(fl(Q1)*fl(Q2)) - Q1*Q2 =", Fr(qq) - Q * Q, " (negative: the 'upper bound' of q1 q2 is below it)")
# (c) an integer exponent that motor B does not see as an integer
import exact as ex, exprparse, canon
E = exprparse.to_x("((7^(2^.5))^(8^(3/6)))")
print("(c) exponent (7^(2^.5))^(8^(3/6)): kind", E.kind, " is_rational:", E.is_rational(), " ball:", E.ball(256).str(20), " canon == 2401:", canon.equal(E, ex.X.rat(2401)))
import nm
print("    nm.is_int(exponent) =", nm.is_int(E), " -> a negative base to this exponent is treated as invalid")
base = exprparse.to_x("(.14-.9)")
from mpmath import mp, mpf
mp.dps = 30
print("    value of (.14-.9)^2401 =", mp.nstr(-(mpf("0.76") ** 2401), 8), "(in range, minus-free, uses 1..9 once)")
