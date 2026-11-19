"""lvl.py, quotients of an exact power and the same power plus or minus a small value (14/11, lvl._off_ratio).
The split 191 of n = 9 had 3^(2^49) / (c + 9^(2^48)): 9^(2^48) = 3^(2^49), so the value is 1 / (1 + c/X), 1 minus
about 10^(-2.7e14); the two logarithms agree to 10^(-3e14) and lvl failed ("product of two very large logs").
Random X = B^k (k huge, so X stays an unexpanded exact power), written a second time as (B^j)^(k/j); c a small
rational; the four quotients X / (c + X'), X / |c - X'|, (c + X') / X, |c - X'| / X must be anchors at 1 whose
enclosure of ln|u - 1| contains the true value (mpmath, 300 digits) with the right sign; c = 0 gives the exact 1.
usage: t_lvloff.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lvl
N = int(sys.argv[1]) if len(sys.argv) > 1 else 400
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 191)
E = lvl.exact
mpmath.mp.dps = 300
bad = ok = 0
for trial in range(N):
    B = [Fr(2), Fr(3), Fr(5), Fr(7), Fr(6), Fr(10)][int(rng.integers(6))]
    j = int(rng.choice([2, 3, 4]))
    k = j * int(rng.integers(10 ** 5, 10 ** 13))
    c = Fr(int(rng.integers(0, 50)), int(rng.choice([1, 2, 5, 10])))
    X = lvl.pw(E(B), E(k), 1)
    Xp = lvl.pw(E(B ** j), E(k // j), 1)                  # the same power written differently
    lnt = (mpmath.log(c.numerator) - mpmath.log(c.denominator) - k * mpmath.log(B)) if c else None
    cases = [("X/(c+X)", lambda: lvl.div(X, lvl.add(E(c), Xp)), -1, "plus"),
             ("X/|c-X|", lambda: lvl.div(X, lvl.sub(E(c), Xp)), 1, "minus"),
             ("(c+X)/X", lambda: lvl.div(lvl.add(E(c), Xp), X), 1, "exact"),
             ("|c-X|/X", lambda: lvl.div(lvl.sub(E(c), Xp), X), -1, "exact")]
    for name, f, sg_want, kind in cases:
        try:
            u = f()
        except lvl.Fail as ex:
            bad += 1; print("FAIL", name, B, k, c, ex); continue
        if not c:
            if u.q != 1:
                bad += 1; print("NOT THE EXACT 1", name, B, k)
            else:
                ok += 1
            continue
        if u.anc is None or u.anc[0] != 1 or u.anc[1] != sg_want:
            bad += 1; print("NOT AN ANCHOR AT 1 WITH SIGN", sg_want, name, B, k, c, u); continue
        # the true ln|u - 1|: ln t - ln(1 + t) (plus), ln t - ln(1 - t) (minus), ln t (exact); |ln(1 +- t)| < 2t
        D = u.anc[2]
        slack = 2 * mpmath.exp(lnt) if kind != "exact" else mpmath.mpf(0)
        if not (mpmath.mpf(D.a) - slack <= lnt <= mpmath.mpf(D.b) + slack):
            bad += 1; print("ENCLOSURE MISSES", name, B, k, c, mpmath.nstr(lnt, 30), D)
        else:
            ok += 1
        if not lvl.outside(u, 3.14159, 3.1416):
            bad += 1; print("NOT DECIDED OUTSIDE A WINDOW AT PI", name)
print(f"t_lvloff: {ok} quotients checked; PROBLEMS {bad}")
sys.exit(1 if bad else 0)
