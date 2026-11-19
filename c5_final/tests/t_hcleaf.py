"""hitcheck.leaf with the writing (13/11): the members of the V tables whose writing has only literals, products,
quotients and powers get lvl's canonical form (tower) from their writing.
  (1) every member of V(S), |S| <= 3 (motor B's full enumeration), and of random 4-digit sets: when the leaf gets a form,
      lvl's evaluation of the writing contains motor B's value (ln, 12 000 bits);
  (2) soundness against motor B's exact arithmetic: two members with the same form are equal for motor B
      (exact.sub = 0), and no member with a form is a rational;
  (3) the case of n = 9, split 10, piece VN.7: (5^-(4^-.9)) of V(459) raised to (8^.6) of V(68) is exactly 1/5, so
      |1/5 - that| = 0 and the hit has no value.
usage: t_hcleaf.py [sets of 4 digits] [seed]"""
import sys, os, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, qx, hitcheck, lvl
import exact as ex
from fractions import Fraction as Fr
N4 = int(sys.argv[1]) if len(sys.argv) > 1 else 4
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
q = qx.QX(9, 0, os.path.join(os.environ.get("RESDIR", paths.SRC), "res_9_0.txt"), 4, 5, verbose=False)
import mpmath
mpmath.mp.dps = 1000
bad = 0; nform = 0; nval = 0; nlin = 0
D = tuple(range(1, 10))
sets = list(itertools.combinations(D, 2)) + list(itertools.combinations(D, 3))
s4 = list(itertools.combinations(D, 4))
sets += [s4[i] for i in rng.choice(len(s4), size=N4, replace=False)]
for S in sets:
    full = q.B().full(tuple(str(d) for d in S)) if len(S) <= 3 else None
    if full is not None:
        vals = list(full.values())
    else:
        e = q.V(S); vals = []
        for j in rng.choice(len(e), size=min(400, len(e)), replace=False):
            vals += q.members(S, e, int(j))
    byform = {}
    for (x, w) in vals:
        nval += 1
        v = hitcheck.leaf(x, w)
        if v.tw is None and v.lin is None:
            continue
        nform += 1
        if v.lin is not None:                                   # (13/11) the linear form against motor B's value
            nlin += 1
            c0, terms = v.lin
            val = mpmath.mpf(c0.numerator) / c0.denominator
            for key, a in terms:
                m = mpmath.mpf(a.numerator) / a.denominator
                for (p_, f) in key:
                    m *= mpmath.power(p_, mpmath.mpf(f.numerator) / f.denominator)
                val += m
            bx = hitcheck.arb_iv(ex.absval(x).ball(lvl.PREC + 64))      # motor B's value, an interval of lvl.PREC bits
            tol = mpmath.mpf(2) ** -(lvl.PREC - 20) * max(1, abs(val))
            if not (mpmath.mpf(bx.a) - tol <= val <= mpmath.mpf(bx.b) + tol):
                bad += 1; print("LINEAR FORM NOT MOTOR B'S VALUE:", S, w, v.lin)
        if ex.absval(x).is_rational():
            bad += 1; print("FORM FOR A RATIONAL:", S, w)
        vw = hitcheck.lvl_of_writing(w)
        if vw is None or vw.tw != v.tw or vw.lin != v.lin:
            bad += 1; print("FORM NOT FROM THE WRITING:", S, w)
        try:
            lx = lvl.lam_iv(lvl.from_x(hitcheck.arb_iv(ex.absval(x).ball(12000))))
            lw = lvl.lam_iv(vw)
            if not (lw.a <= lx.b and lx.a <= lw.b):
                bad += 1; print("WRITING DOES NOT CONTAIN MOTOR B'S VALUE:", S, w)
        except lvl.Fail:
            pass
        byform.setdefault((v.tw, v.lin), []).append((x, w))
    for f, xs in byform.items():
        for (x, w) in xs[1:]:
            if ex.sign(ex.sub(ex.absval(x), ex.absval(xs[0][0]))) != 0:
                bad += 1; print("SAME FORM, DIFFERENT VALUES:", S, w, xs[0][1])
print(f"(1)-(2) {len(sets)} sets, {nval} members, {nform} with a form ({nlin} linear), problems {bad}", flush=True)
# (4) the case of n = 9, split 14, pieces VN.2 and VN.5: (4^-.9) of V(49) + (5-(8^-.6)) of V(568) = 5 exactly
A2 = [(x, w) for (x, w) in q.B().full(("4", "9")).values() if w.replace(" ", "") == "(4^-.9)"]
B2 = [(x, w) for (x, w) in q.B().full(("5", "6", "8")).values() if w.replace(" ", "") == "(5-(8^-.6))"]
case2 = "not found"
if A2 and B2:
    z = lvl.add(hitcheck.leaf(*A2[0]), hitcheck.leaf(*B2[0]))
    case2 = f"(4^-.9) + (5-(8^-.6)) = {z}"
    if z.q != 5:
        bad += 1
else:
    print("writings found:", len(A2), len(B2)); bad += 1
print("split 14 case:", case2, flush=True)
# (3)
A = [(x, w) for (x, w) in q.B().full(("4", "5", "9")).values() if w.replace(" ", "") == "(5^-(4^-.9))"]
B = [(x, w) for (x, w) in q.B().full(("6", "8")).values() if w.replace(" ", "") == "(8^.6)"]
case = "not found"
if A and B:
    a = hitcheck.leaf(*A[0]); b = hitcheck.leaf(*B[0])
    z = lvl.pw(a, b, 1)
    case = f"(5^-(4^-.9))^(8^.6) = {z}"
    if z.q != Fr(1, 5):
        bad += 1
    try:
        lvl.pw(lvl.sub(lvl.exact(Fr(1, 5)), z), lvl.exact(7), -1)
        bad += 1; case += ", |1/5 - z|^-7 HAS A VALUE"
    except lvl.Invalid:
        case += ", |1/5 - z|^-7 has no value"
else:
    print("writings found:", len(A), len(B)); bad += 1
print(f"t_hcleaf: {nform} members with a form among {nval}; lab case: {case}; PROBLEMS {bad}")
sys.exit(1 if bad else 0)
