"""(17/11, parche 7) lvl: |a - b| of two huge values that are equal or almost equal (the undecided values of n = 9 left by
the parts and the cloud sessions on 17/11: splits 83, 88, 90, 92, 112, 115, 119, 121, 123, 128, 146, 168, 176, 177, 179,
180, 192).  Rules: bound_diff (bounds of |a - b| from the logs when no rule gives the difference), the anchored
difference |q - aq (1 + sg e^D)| = |(q - aq) - sg aq e^D| bounded from its two exact pieces, x^1 keeping the exact forms
of x, and an exact power known by its prime vector compared with an anchor.
  (1) bound_diff against the truth: random ln a, ln b (overlapping enclosures), random points x, y inside them:
      e^lo <= |e^x - e^y| <= e^hi (mpmath, 600 digits);
  (2) anchored near-cancellations with a known value: A = B^k, A' = (B^(k m) + c)^(1/m) = A (1 + c / (m A^m) + ...):
      ln|A - A'| = ln(c/m) - (m - 1) k ln B up to a relative c/A^m; the result (a value or bounds) must contain it and be
      decided outside a window at pi;
  (3) the recorded undecided values (one of each family, rational leaves and motor B writings): decided outside the
      windows of pi, gamma and delta (as hitcheck: the precisions in turn);
  (6) (18/11) towers written in two ways (3^(2^a) = 9^(8^b), a = 1 + 3b): (T + c) - T = c exactly, x 1 = x / 1 = x,
      1 / (9^-(8^b)) = 9^(8^b); two different towers are never made equal;
  (5) the exact forms of |a - b| and a + b with a huge exact power (the sign error of the parche 6, |c - P| with the form
      c - P): every form equals the exact value (Python fractions);
  (4) no false decision: values that are exactly a constant near the window are never placed outside it
      (|T - (T' + 355/113)| = 355/113 with T = T' written differently, |3^(8^4.8) - (1 + 9^(4^6.7))| = 1, a window around
      2 for |T - T'| + 2 written as a sum) at any precision.
usage: t_lvlbound.py [trials] [seed]"""
import sys, os
here = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault("PANDI_SRC", os.path.join(os.path.dirname(here), "cert"))
os.environ.setdefault("PANDI_RESULTS", os.path.join(os.path.dirname(os.path.dirname(here)), "data"))
sys.path.insert(0, here)
sys.path.insert(0, os.path.join(os.path.dirname(here), "cert"))
import numpy as np, mpmath
from fractions import Fraction as Fr
from mpmath import iv
import lvl
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
E = lvl.exact
PRECS = (400, 800, 1600, 3200, 6400, 12800, 25600, 51200)
bad = 0
counts = {"ineq": 0, "anc": 0, "rec": 0, "nofalse": 0}


def problem(*a):
    global bad
    bad += 1
    print("PROBLEM", *a, flush=True)


# ---- (1) the inequality of bound_diff
mpmath.mp.dps = 600
for t in range(N):
    lvl.set_prec(400)
    m = mpmath.mpf(float(rng.uniform(0.5, 6000)))
    wa = mpmath.mpf(10) ** float(-rng.uniform(2, 110)) * m
    wb = mpmath.mpf(10) ** float(-rng.uniform(2, 110)) * m
    off = (wa + wb) * float(rng.uniform(-0.9, 0.9))
    la = iv.mpf([m - wa, m + wa])
    lb = iv.mpf([m + off - wb, m + off + wb])
    try:
        a, b = lvl.from_lnu(la), lvl.from_lnu(lb)
        v = lvl.bound_diff(a, b)
    except lvl.Fail as ex:
        problem("bound_diff failed", ex); continue
    lo, hi = v.bnd
    A, B = lvl.lam_iv(a), lvl.lam_iv(b)
    for _ in range(6):
        x = mpmath.mpf(A.a) + (mpmath.mpf(A.b) - mpmath.mpf(A.a)) * float(rng.random())
        y = mpmath.mpf(B.a) + (mpmath.mpf(B.b) - mpmath.mpf(B.a)) * float(rng.random())
        if x == y:
            continue
        lt = max(x, y) + mpmath.log(-mpmath.expm1(-abs(x - y)))      # ln |e^x - e^y|
        tol = mpmath.mpf(10) ** -500 * (1 + abs(lt))
        if not (lt <= mpmath.mpf(hi) + tol and (lo is None or mpmath.mpf(lo) - tol <= lt)):
            problem("bound_diff misses", mpmath.nstr(lt, 20), lo, hi)
            break
    else:
        counts["ineq"] += 1
lvl.set_prec(400)

# ---- (2) anchored near-cancellations
mpmath.mp.dps = 300
BASES = [Fr(5, 3), Fr(10, 3), Fr(3), Fr(2), Fr(100, 9), Fr(9, 4), Fr(7, 2)]
for t in range(max(12, N // 10)):
    B = BASES[int(rng.integers(len(BASES)))]
    m = int(rng.choice([2, 4, 5, 10]))
    k = int(rng.integers(2000, 80000))
    c = Fr(int(rng.integers(1, 10)), int(rng.choice([1, 5, 10])))
    tr = mpmath.log(c.numerator) - mpmath.log(c.denominator) - mpmath.log(m) - (m - 1) * k * (mpmath.log(B.numerator) - mpmath.log(B.denominator))
    done = False
    for p in PRECS:
        lvl.set_prec(p)
        try:
            A = lvl.pw(E(B), E(k), 1)
            A2 = lvl.pw(lvl.add(E(c), lvl.pw(E(B), E(k * m), 1)), E(Fr(1, m)), 1)
            u = lvl.sub(A, A2)
        except lvl.Fail as ex:
            continue
        tol = mpmath.mpf(10) ** -250 * abs(tr)
        if u.far == "bound":
            lo, hi = u.bnd
            if not (tr <= mpmath.mpf(hi) + tol and (lo is None or mpmath.mpf(lo) - tol <= tr)):
                problem("anchored: bounds miss the truth", B, k, m, c, mpmath.nstr(tr, 20), lo, hi); done = True; break
        elif u.far is None:
            try:
                L = lvl.lam_iv(u)
            except lvl.Fail:
                L = None
            if L is not None and not (mpmath.mpf(L.a) - tol <= tr <= mpmath.mpf(L.b) + tol):
                problem("anchored: enclosure misses the truth", B, k, m, c, mpmath.nstr(tr, 20), L); done = True; break
        if lvl.outside(u, 3.14159, 3.1416):
            counts["anc"] += 1; done = True
            break
    if not done:
        problem("anchored: not decided outside a window at pi", B, k, m, c)
lvl.set_prec(400)

# ---- (5) the exact forms of differences and sums with a huge exact power (the sign error of the parche 6: |c - P| got
# the form c - P): random trees of +, |-| over small rationals and huge exact powers P = B^K (> 600 000 bits); every
# node that is an exact rational and carries a form: the form gives exactly that rational, and the rational is the
# true value (Python fractions)
counts["forms"] = 0
lvl.set_prec(400)
HUGE = [(Fr(3), 390625), (Fr(5), 262144), (Fr(9), 196608), (Fr(10, 3), 390625), (Fr(7), 215000)]
for t in range(max(30, N // 4)):
    B, K = HUGE[int(rng.integers(len(HUGE)))]
    P = lvl.pw(E(B), E(K), 1)
    PT = B ** K
    leaves = [(P, PT)] + [(E(c), c) for c in (Fr(int(rng.integers(1, 30)), int(rng.choice([1, 5, 10]))) for _ in range(3))]
    nodes = list(leaves)
    for _ in range(6):
        (a, ta), (b, tb) = nodes[int(rng.integers(len(nodes)))], nodes[int(rng.integers(len(nodes)))]
        try:
            if rng.random() < 0.5:
                v, tv = lvl.sub(a, b), abs(ta - tb)
            else:
                v, tv = lvl.add(a, b), ta + tb
        except lvl.Fail:
            continue
        if v.far:
            continue
        nodes.append((v, tv))
        if v.q is not None and v.q != tv:
            problem("forms: wrong exact value", B, K)
        if v.lin is not None:
            c0, terms = v.lin
            tot = Fr(c0)
            ok_ = True
            for key, a_ in terms:
                pv = lvl.pf_value(dict(key))
                if pv is None:
                    ok_ = False; break
                tot += Fr(a_) * pv
            if ok_ and tot != tv:
                problem("forms: the form of a value is not its exact value (sign or constant)", B, K, c0)
                break
            counts["forms"] += 1
lvl.set_prec(400)

# ---- (6) (18/11) powers with an irrational exponent (towers) written in two ways, T = 3^(2^a) = 9^(8^b) with
# a = 1 + 3b (the splits 187 and 252: 3^(2^(7/.4)) and 9^(8^(6-.5))): |T - (T' + c)| = c, |T - (T' - c)| = c,
# |T - 1 T'| = |T - T'/1| = |T - 1/(9^-(8^b))| = 0 exactly; and two towers that are NOT equal (a != 1 + 3b) are never
# made equal: their difference must enclose the truth (mpmath, small towers) and never be an exact 0 or c
counts["towers"] = 0
lvl.set_prec(400)
mpmath.mp.dps = 120
for t in range(max(30, N // 5)):
    small = rng.random() < 0.5
    bnum = int(rng.integers(1, 9)); bden = int(rng.choice([2, 4, 5, 10]))
    b_ = Fr(bnum, bden) if small else Fr(int(rng.integers(40, 120)), 10)
    a_ = 1 + 3 * b_
    c = Fr(int(rng.integers(1, 30)), int(rng.choice([1, 2, 10])))
    def T3(a):   # 3^(2^a)
        return lvl.pw(E(Fr(3)), lvl.pw(E(Fr(2)), E(a), 1), 1)
    def T9(b, sg=1):   # 9^(sg 8^b)
        return lvl.pw(E(Fr(9)), lvl.pw(E(Fr(8)), E(b), 1), sg)
    if a_.denominator == 1 or b_.denominator == 1:
        continue                                          # rational exponents are another rule
    try:
        A, B = T3(a_), T9(b_)
        cases = [("T - (T' + c)", lvl.sub(A, lvl.add(E(c), B)), c), ("T - (T' - c)", lvl.sub(A, lvl.sub(E(c), B)), c),
                 ("(T' + c) - T", lvl.sub(lvl.add(B, E(c)), A), c), ("T - 1 T'", lvl.sub(A, lvl.mul(E(Fr(1)), B)), Fr(0)),
                 ("T - T'/1", lvl.sub(A, lvl.div(B, E(Fr(1)))), Fr(0)),
                 ("T - 1/(9^-8^b)", lvl.sub(A, lvl.div(E(Fr(1)), T9(b_, -1))), Fr(0))]
    except lvl.Fail as ex:
        problem("towers: Fail", a_, b_, ex); continue
    if small:
        # small towers: the truth by mpmath (|T - (T' - c)| is c only when c < T'); the enclosure must contain it.
        # (18/11) the working precision must exceed the digits of T' (up to 10^1400 here): at a fixed 120 digits,
        # T' - (T' - c) loses c and the test's own "truth" was 0 (seed 5, b = 5/2: T' ~ 10^172.7)
        mpmath.mp.dps = 120 + int(float(8 ** (b_.numerator / b_.denominator)) * 0.9542425094) + 10
        xb = mpmath.mpf(9) ** (mpmath.mpf(8) ** (mpmath.mpf(b_.numerator) / b_.denominator))
        cc = mpmath.mpf(c.numerator) / c.denominator
        truth = {"T - (T' + c)": cc, "T - (T' - c)": abs(xb - abs(cc - xb)), "(T' + c) - T": cc}
    for nm, v, want in cases:
        if small:
            tv = truth.get(nm, mpmath.mpf(0))
            if v.q is not None:
                ok_ = abs(mpmath.mpf(v.q.numerator) / v.q.denominator - tv) < mpmath.mpf(10) ** -80
            elif v.x is not None:
                ok_ = mpmath.mpf(v.x.a) - mpmath.mpf(10) ** -80 <= tv <= mpmath.mpf(v.x.b) + mpmath.mpf(10) ** -80
            else:
                ok_ = True                                # no enclosure to compare (it stays undecided)
            if not ok_:
                problem("towers: the value misses the truth", nm, a_, b_, c, repr(v)[:120])
            else:
                counts["towers"] += 1
        elif v.q is None or v.q != want:
            problem("towers: not the exact value", nm, a_, b_, c, repr(v)[:120])
        else:
            counts["towers"] += 1
    # a different tower: a2 = a + 1/7 (not equal): never an exact 0 or c, and the truth inside when computable
    a2 = a_ + Fr(1, 7)
    try:
        A2 = T3(a2)
        v0 = lvl.sub(A2, B); vc = lvl.sub(A2, lvl.add(E(c), B))
    except lvl.Fail:
        continue
    for nm, v, bad_q in (("T2 - T'", v0, Fr(0)), ("T2 - (T' + c)", vc, c)):
        if v.q is not None and v.q == bad_q:
            problem("towers: two different towers made equal", nm, a2, b_)
        elif small and not v.far:
            x2 = mpmath.mpf(3) ** (mpmath.mpf(2) ** (mpmath.mpf(a2.numerator) / a2.denominator))
            xb = mpmath.mpf(9) ** (mpmath.mpf(8) ** (mpmath.mpf(b_.numerator) / b_.denominator))
            tr = abs(x2 - xb - (mpmath.mpf(c.numerator) / c.denominator if nm.endswith("c)") else 0))
            try:
                L = lvl.lam_iv(v)
                if not (mpmath.mpf(L.a) - 1e-60 <= mpmath.log(tr) <= mpmath.mpf(L.b) + 1e-60):
                    problem("towers: enclosure misses the truth", nm, a2, b_, c)
            except lvl.Fail:
                pass
mpmath.mp.dps = 120
lvl.set_prec(400)

# ---- (3) the recorded undecided values, one of each family
import hitcheck, exprparse
REC = [
    ("83", "[|s-l| | [s^l | 2/1 | 6561/1] | [s^-l | 1/2 | [l^s | 8/1 | [s+l | 3/1 | [s^l | 1/10 | 117649/1]]]]]"),
    ("83", "[|s-l| | [s^-l | 1/5 | 6561/1] | [s^l | 5/1 | [l^s | 8/1 | [|s-l| | 3/1 | [s^l | 1/10 | 117649/1]]]]]"),
    ("88", "[|s-l| | [s^-l | 3/5 | 78125/1] | [l^s | 1/10 | [s+l | 3/10 | [s^-l | 9/25 | 390625/1]]]]"),
    ("90/121", "[|s-l| | [s^l | 9/1 | 15625/1] | [l^-s | 4/1 | [l^-s | 1/10 | [|s-l| | 8/1 | [s^l | 3/1 | 78125/1]]]]]"),
    ("92", "[|s-l| | [s^l | 2/1 | 4782969/1] | [l^s | 1/2 | [s+l | 1/10 | [s^l | 4/1 | 4782969/1]]]]"),
    ("92/123", "[|s-l| | [s^l | 9/1 | 78125/1] | [l^s | 4/1 | [l^s | 1/10 | [s+l | 3/5 | [s^l | 3/1 | 390625/1]]]]]"),
    ("115", "[|s-l| | (9^64) | [s^l | 3/1 | [l^s | 7/1 | [s+l | 2/1 | [s^-l | 8/1 | 100000/1]]]]]"),
    ("119", "[|s-l| | [s^-l | 3/5 | 78125/1] | [l^s | 1/10 | [|s-l| | 3/10 | [l^-s | 390625/1 | 9/25]]]]"),
    ("123", "[|s-l| | [s^-l | 1/2 | 4782969/1] | [l^s | 1/1 | [s+l | 4/1 | [s^l | 2/1 | 4782969/1]]]]"),
    ("128", "[|s-l| | [s^l | 9/1 | 2097152/1] | [l^s | 1/1 | [s+l | 1/5 | [s^l | 3/1 | 4194304/1]]]]"),
    ("146/177", "[|s-l| | [s^-l | 9/100 | 15625/1] | [l^s | 2/5 | [s+l | 4/5 | [s^-l | 3/10 | 78125/1]]]]"),
    ("168", "[|s-l| | [s^l | 5/1 | (.1^-(4^.9))] | [s^-l | 1/5 | ((3+7)^(8^.6))]]"),
    ("176", "[|s-l| | [s^-l | 1/2 | (.1^-(8^.6))] | [s^l | 2/1 | ((3+7)^(4^.9))]]"),
    ("179", "[|s-l| | [s^-l | 9/100 | 78125/1] | [l^s | 2/5 | [|s-l| | 3/5 | [s^-l | 3/10 | 390625/1]]]]"),
    ("180", "[|s-l| | ((.9^.1)^-(8^5)) | [l^s | 1/5 | [s^-l | 9/10 | 16384/1]]]"),
    ("192", "[|s-l| | [s^l | 3/1 | (8^(5-.2))] | [s+l | 1/1 | [s^l | 9/1 | (4^6.7)]]]"),
    ("192", "[|s-l| | [s^l | 3/1 | (8^(5-.2))] | [s/l | 1/1 | [s^-l | 9/1 | (4^6.7)]]]"),
    ("187", "[|s-l| | [s^l | 3/1 | (2^(7/.4))] | [s+l | 1/10 | [s^l | 9/1 | (8^(6-.5))]]]"),
    ("187", "[|s-l| | [s^l | 3/1 | (2^(7/.4))] | [|s-l| | 1/1 | [s^l | 9/1 | (8^(6-.5))]]]"),
    ("187", "[|s-l| | [s^l | 3/1 | (2^(7/.4))] | [s/l | 1/1 | [s^-l | 9/1 | (8^(6-.5))]]]"),
    ("252", "[|s-l| | [s^l | 9/1 | (8^(6-.5))] | [s*l | 1/1 | [s^l | 3/1 | (2^(7/.4))]]]"),
    ("252", "[|s-l| | [s^l | 9/1 | (8^(6-.5))] | [l/s | 1/1 | [s^l | 3/1 | (2^(7/.4))]]]"),
]


def parse(s, i=0):
    """the value of a recorded writing of a hit (hitcheck: [form | A | B], rationals p/q, motor B writings (...)), or
    a string for a failure"""
    if s[i] == "[":
        f = next(f for f in sorted(lvl.FORMS, key=len, reverse=True) if s.startswith(f, i + 1))
        i = s.index("|", i + 1 + len(f)) + 1
        while s[i] == " ":
            i += 1
        a, i = parse(s, i)
        while s[i] in " |":
            i += 1
        b, i = parse(s, i)
        while s[i] == " ":
            i += 1
        if isinstance(a, str) or isinstance(b, str):
            return (a if isinstance(a, str) else b), i + 1
        try:
            return lvl.op(f, a, b), i + 1
        except lvl.Fail as ex:
            return "FAIL: " + str(ex), i + 1
    if s[i] == "(":
        d, j = 0, i
        while True:
            d += {"(": 1, ")": -1}.get(s[j], 0)
            if d == 0:
                break
            j += 1
        w = s[i:j + 1]
        return hitcheck.leaf(exprparse.to_x(w), w), j + 1
    j = i
    while j < len(s) and s[j] not in " |]":
        j += 1
    return E(Fr(s[i:j])), j


def first_outside(expr, lo, hi):
    for p in PRECS:
        lvl.set_prec(p)
        v, _ = parse(expr)
        if not isinstance(v, str) and lvl.outside(v, lo, hi):
            lvl.set_prec(400)
            return p
    lvl.set_prec(400)
    return None


mpmath.mp.dps = 50
CONST = {"pi": mpmath.pi, "gamma": mpmath.euler, "delta": mpmath.mpf("4.669201609102990671853203820466201617258185577475768632745651343004134330211314737138689744023948013817165984855189815134408627142027932522312442988890890859944935463236713411532481714219947455644365823793202009561058330575458617652222070385410646749494284981453391726200568755665952339875603825637225648")}
for tag, expr in REC:
    for c, x in CONST.items():
        lo, hi = float(x - mpmath.mpf("1e-11")), float(x + mpmath.mpf("1e-11"))
        p = first_outside(expr, lo, hi)
        if p is None:
            problem("recorded value still undecided", tag, c, expr)
        else:
            counts["rec"] += 1
    print(f"  split {tag}: decided outside the windows of pi, gamma and delta", flush=True)

# ---- (4) no false decision
NOFALSE = [
    ("|T - (T' + 355/113)| = 355/113", "[|s-l| | ((.9^.1)^-(8^5)) | [s+l | 355/113 | [l^s | 1/5 | [s^-l | 9/10 | 16384/1]]]]",
     Fr(355, 113)),
    ("|3^(8^4.8) - (1 + 9^(4^6.7))| = 1", "[|s-l| | [s^l | 3/1 | (8^(5-.2))] | [s+l | 1/1 | [s^l | 9/1 | (4^6.7)]]]", Fr(1)),
    ("|(5/3)^78125 - ((25/9)^390625 + 3/10)^(1/10) - 2| ~ 2", "[s+l | 2/1 | [|s-l| | [s^-l | 3/5 | 78125/1] | [l^s | 1/10 | [s+l | 3/10 | [s^-l | 9/25 | 390625/1]]]]]",
     Fr(2)),
    ("|3^(2^(7/.4)) - (9^(8^(6-.5)) + 355/113)| = 355/113", "[|s-l| | [s^l | 3/1 | (2^(7/.4))] | [s+l | 355/113 | [s^l | 9/1 | (8^(6-.5))]]]",
     Fr(355, 113)),
    ("|2^4782969 - (4^4782969 + 1/10)^(1/2) + 1/5^...| (|x - (x + 5/2)|)", "[|s-l| | [s^l | 2/1 | 4782969/1] | [s+l | 5/2 | [l^s | 1/2 | [s^l | 4/1 | 4782969/1]]]]",
     Fr(5, 2)),
]
for tag, expr, val in NOFALSE:
    x = float(val)
    for (lo, hi) in ((x - 1e-9, x + 1e-9), (x - 1e-3, x + 1e-3)):
        for p in PRECS:
            lvl.set_prec(p)
            v, _ = parse(expr)
            if not isinstance(v, str) and lvl.outside(v, lo, hi):
                problem("FALSE DECISION: a value equal to", val, "placed outside", (lo, hi), "at", p, "bits:", tag, repr(v)[:200])
                break
        else:
            counts["nofalse"] += 1
        lvl.set_prec(400)
    print(f"  no false decision: {tag}", flush=True)

print(f"t_lvlbound: {counts['ineq']} bounds of |e^x - e^y| hold, {counts['anc']} anchored near-cancellations decided "
      f"with the truth inside, {counts['forms']} exact forms equal to their values, {counts['towers']} tower identities exact, {counts['rec']} recorded values (and "
      f"constant) decided, {counts['nofalse']} windows around a true value never excluded; PROBLEMS {bad}")
sys.exit(1 if bad else 0)
