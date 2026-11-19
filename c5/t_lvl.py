"""self-test of lvl.py: random expressions (biased to deep and near-1 values) evaluated by lvl and by
direct mpmath at 3000 digits; every field lvl gives must contain the direct value"""
import sys, random
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
import mpmath
import lvl

mp = mpmath.mp.clone() if hasattr(mpmath.mp, "clone") else None
LEAVES = ["1", "2", "3", "4", "5", ".5", ".1", "6", "7", ".7", "9", "1.5", ".2", "8"]
FORMS = lvl.FORMS
W = {"s+l": 1, "|s-l|": 3, "s*l": 1, "s/l": 1, "l/s": 1, "s^l": 3, "s^-l": 3, "l^s": 2, "l^-s": 2}


def tree(d):
    if d == 0 or random.random() < 0.25:
        return random.choice(LEAVES)
    f = random.choices(FORMS, [W[f] for f in FORMS])[0]
    return (f, tree(d - 1), tree(d - 1))


def ev_lvl(t):
    if isinstance(t, str):
        return lvl.exact(Fr(t))
    return lvl.op(t[0], ev_lvl(t[1]), ev_lvl(t[2]))


def ev_mp(t):
    """direct value at 3000 digits; None if too large"""
    if isinstance(t, str):
        return mpmath.mpf(t)
    a, b = ev_mp(t[1]), ev_mp(t[2])
    if a is None or b is None:
        return None
    f = t[0]
    if f == "s+l": r = a + b
    elif f == "|s-l|": r = abs(a - b)
    elif f == "s*l": r = a * b
    elif f == "s/l": r = a / b if b != 0 else "inv"
    elif f == "l/s": r = b / a if a != 0 else "inv"
    else:
        base, e = (a, b) if f in ("s^l", "s^-l") else (b, a)
        if f.endswith("-l") or f.endswith("-s"):
            e = -e
        if base == 0:
            r = mpmath.mpf(0) if e > 0 else "inv"
            if e == 0: r = "inv"
        else:
            lu = e * mpmath.log(base)
            if abs(lu) > 20000:
                return None
            r = mpmath.exp(lu)
    if isinstance(r, str):
        raise lvl.Invalid
    if r != 0 and abs(mpmath.log(r)) > 20000:
        return None
    return r


def check(v, u):
    """direct value u inside every field of v"""
    if v.far:
        return v.far == "near1" and abs(u - 1) < mpmath.mpf(10) ** -100 or v.far == "huge" and u > 1e100 or v.far == "tiny" and u < 1e-100
    ok = True
    if v.q is not None:
        # the direct value is computed through exp/log at 3000 digits, so an exact zero of lvl can show up there
        # as a residue of order 10^-3000 (e.g. |.5 - 3^1/6|): the tolerance is relative, with a floor of 10^-2990
        qv = mpmath.mpf(v.q.numerator) / v.q.denominator
        return abs(u - qv) <= max(abs(u), abs(qv)) * mpmath.mpf(10) ** -2900 + mpmath.mpf(10) ** -2990
    if v.x is not None:
        ok &= v.x.a <= u <= v.x.b
    if v.s is not None and v.s != 0 and u != 0:
        lu = mpmath.log(u)
        if lu != 0 and abs(lu) > mpmath.mpf(10) ** -2900:
            ok &= (lu > 0) == (v.s > 0)
            L = mpmath.log(abs(lu))
            ok &= v.Lam.a - 1e-50 <= L <= v.Lam.b + 1e-50
    if v.anc is not None:
        aq, sg, D = v.anc
        t = u / (mpmath.mpf(aq.numerator) / aq.denominator) - 1
        if abs(t) > mpmath.mpf(10) ** -2900:
            ok &= (t > 0) == (sg > 0)
            ok &= D.a - 1e-50 <= mpmath.log(abs(t)) <= D.b + 1e-50
    return ok


if __name__ == "__main__":
    random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    if len(sys.argv) > 3:                       # working precision of lvl (hitcheck raises it for undecided hits)
        lvl.set_prec(int(sys.argv[3]))
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    mpmath.mp.dps = 3000
    bad = fails = cmp = skipped = 0
    kinds = {}
    for k in range(N):
        t = tree(4)
        try:
            v = ev_lvl(t)
        except lvl.Invalid:
            continue
        except lvl.Fail:
            fails += 1
            continue
        mpmath.mp.dps = 3000
        try:
            u = ev_mp(t)
        except (lvl.Invalid, ZeroDivisionError):
            continue
        finally:
            lvl.set_prec(lvl.PREC)
        mpmath.mp.dps = 3000
        if u is None:
            skipped += 1
            continue
        kinds[("q" if v.q is not None else "anc" if v.anc else "x" if v.x is not None else "far" if v.far else "lam")] = kinds.get(("q" if v.q is not None else "anc" if v.anc else "x" if v.x is not None else "far" if v.far else "lam"), 0) + 1
        cmp += 1
        if not check(v, u):
            bad += 1
            print("BAD", t, v, mpmath.nstr(u, 30))
    print("compared", cmp, "bad", bad, "lvl Fail", fails, "too large for direct", skipped, kinds)
