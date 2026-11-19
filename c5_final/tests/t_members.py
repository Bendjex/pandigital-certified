"""The exact members of the entries of the V tables (motorx.MX.class_members, dedup, refine_near1; 10/11), against
motor B's unfiltered values and an independent evaluation of their writings (exact rationals while possible, then mpmath at 12 000 bits).

What the exact decisions of motor X need (invariant I): every a in V(S) has an entry of the table of S whose
statement holds for a (the exact rational, the exact 1, or A <= ln a <= B) and whose list class_members contains a.

  (1) all the values of V(S) (motor B's full enumeration), |S| <= 3: invariant I, value by value;
  (2) |S| = 4, values of motor B in windows around random entries, around exact rationals and around 1: invariant I;
  (3) every value that the test of class_members removes from motor B's answer has its ln outside [A, B] (mpmath);
  (4) ln_bounds(v) contains ln v (mpmath);
  (5) dedup keeps exactly the distinct rows of the raw table.
usage: t_members.py [sets of 3 digits] [sets of 4 digits] [entries per set] [seed]"""
import sys, os, itertools
from fractions import Fraction
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, qx, ivl, motorx, xl
import exact as ex
from flint import arb
import mpmath
from mpmath import mp, mpf

N3 = int(sys.argv[1]) if len(sys.argv) > 1 else 10
N4 = int(sys.argv[2]) if len(sys.argv) > 2 else 6
NE = int(sys.argv[3]) if len(sys.argv) > 3 else 60
rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 10)
mp.prec = 12000
q = qx.QX(9, 0, os.path.join(paths.SRC, "res_9_0.txt"), 4, 5, verbose=False)
bad = 0
cnt = {"values": 0, "removed": 0, "lnb": 0, "rows": 0}


# ---------------------------------------------------------------- independent evaluation of a writing
def parse(s):
    """value (mpf) of a writing: literals (digits and a point), reduced rationals p/q, (a op b) with op in
    + - * / ^ ^-, and |a-b|"""
    pos = [0]

    def number():
        i = pos[0]
        while pos[0] < len(s) and (s[pos[0]].isdigit() or s[pos[0]] == "."):
            pos[0] += 1
        t = s[i:pos[0]]
        assert t, (s, i)
        v = Fraction(t if not t.startswith(".") else "0" + t)          # literals are exact decimals
        # a rational p/q: two integers joined by a slash
        if "." not in t and pos[0] < len(s) and s[pos[0]] == "/" and pos[0] + 1 < len(s) and s[pos[0] + 1].isdigit():
            j = pos[0] + 1
            while j < len(s) and s[j].isdigit():
                j += 1
            if j == len(s) or s[j] != ".":
                v = v / Fraction(s[pos[0] + 1:j]); pos[0] = j
        return v

    def m(x):
        return mpf(x.numerator) / x.denominator if isinstance(x, Fraction) else x

    def exact_ok(a, b):
        return isinstance(a, Fraction) and isinstance(b, Fraction)

    def expr():
        c = s[pos[0]]
        if c == "(":
            pos[0] += 1
            a = expr()
            op = s[pos[0]]; pos[0] += 1
            neg = False
            if op == "^" and s[pos[0]] == "-":
                neg = True; pos[0] += 1
            b = expr()
            assert s[pos[0]] == ")", (s, pos[0]); pos[0] += 1
            # (11/11) exact rational arithmetic while both operands are rationals (an exact 1 such as .4*(.5/(.9-.7))
            # must have ln exactly 0, which decimal floats at 12 000 bits do not give)
            if op in "+-*/" and exact_ok(a, b):
                return a + b if op == "+" else a - b if op == "-" else a * b if op == "*" else a / b
            if op == "^" and exact_ok(a, b) and b.denominator == 1 and abs(b.numerator) <= 4096 and a != 0:
                k = -b.numerator if neg else b.numerator
                if abs(k) * max(a.numerator.bit_length(), a.denominator.bit_length()) <= 200000:
                    return a ** k
            a, b = m(a), m(b)
            if op == "+":
                return a + b
            if op == "-":
                return a - b
            if op == "*":
                return a * b
            if op == "/":
                return a / b
            assert op == "^", (s, op)
            return mpmath.power(a, -b if neg else b)
        if c == "|":
            pos[0] += 1
            a = expr()
            assert s[pos[0]] == "-", (s, pos[0]); pos[0] += 1
            b = expr()
            assert s[pos[0]] == "|", (s, pos[0]); pos[0] += 1
            return abs(a - b) if exact_ok(a, b) else abs(m(a) - m(b))
        return number()

    v = expr()
    assert pos[0] == len(s), (s, pos[0])
    return v


_LN = {}


def ln_of(w):
    if w not in _LN:
        v = parse(w)
        if isinstance(v, Fraction):
            _LN[w] = mpf("-inf") if v == 0 else (mpf(0) if v == 1 else mpmath.log(v.numerator) - mpmath.log(v.denominator))
        else:
            _LN[w] = mpmath.log(v) if v > 0 else mpf("-inf")
    return _LN[w]


def holds(e, j, a, w):
    """the statement of the entry j of the table e, for the exact value a with writing w"""
    if e.Z[j]:
        return ex.sign(a) == 0
    if e.Q[j] > 0:
        return ex.sign(ex.sub(a, ex.X.rat(Fraction(int(e.P[j]), int(e.Q[j]))))) == 0
    if e.O[j]:
        return ex.sign(ex.sub(a, ex.ONE)) == 0
    l = ln_of(w)
    return mpf(float(e.A[j])) <= l <= mpf(float(e.B[j]))


def check_value(S, e, a, w, order):
    """invariant I for the value a of V(S)"""
    global bad
    cnt["values"] += 1
    if ex.sign(a) == 0:
        cand = np.where(e.Z)[0]
    else:
        l = float(ln_of(w))
        i0 = np.searchsorted(e.A[order], l + 1e-9, side="right")
        cand = order[max(0, i0 - 4000):i0]
        cand = cand[(e.B[cand] >= l - 1e-9)]
        cand = cand[~e.Z[cand]]
    k = a.key()
    for j in cand:
        j = int(j)
        if holds(e, j, a, w) and any(v.key() == k for (v, _) in q.members(S, e, j)):
            return
    bad += 1
    print(f"NO ENTRY for {w} of V{S}: {len(cand)} candidates,",
          [(int(j), bool(holds(e, int(j), a, w)), len(q.members(S, e, int(j)))) for j in cand[:6]], flush=True)


def window_of(e, i, widen=0):
    lo = float(ivl.exp_lo(np.array([e.A[i]]))[0]); hi = float(ivl.exp_hi(np.array([e.B[i]]))[0])
    lo, hi = float(np.nextafter(lo, 0)), float(np.nextafter(hi, np.inf))
    for _ in range(widen):
        lo, hi = float(np.nextafter(lo, 0)), float(np.nextafter(hi, np.inf))
    return lo, hi


def bq(S, lo, hi):
    return list(q.B().query(tuple(str(d) for d in S), arb(lo), arb(hi)).values())


def raw_table(S):
    A, B, fl = xl.table(motorx.bits(S)); P, Q = xl.rats(motorx.bits(S))
    return motorx.Ent(A, B, (fl & 1) > 0, (fl & 2) > 0, None, P, Q)


D = tuple(range(1, 10))
# ---------------------------------------------------------------- (1) full enumeration, |S| <= 3
sets = [(d,) for d in D] + list(itertools.combinations(D, 2))
s3 = list(itertools.combinations(D, 3))
sets += [s3[i] for i in rng.choice(len(s3), size=min(N3, len(s3)), replace=False)]
for S in sets:
    e = q.V(S)
    order = np.argsort(e.A, kind="stable")
    full = q.B().full(tuple(str(d) for d in S))
    for (a, w) in full.values():
        check_value(S, e, a, w, order)
    # (5) dedup
    r = raw_table(S)
    k1 = {bytes(x) for x in r.uniq_key()}; k2 = {bytes(x) for x in motorx.MX.dedup(r).uniq_key()}
    cnt["rows"] += len(r)
    if k1 != k2 or len(motorx.MX.dedup(r)) != len(k1):
        bad += 1; print("DEDUP changes the rows of", S)
print(f"(1) {len(sets)} sets of 1-3 digits: {cnt['values']} values, problems so far {bad}", flush=True)

# ---------------------------------------------------------------- (2)-(4) |S| = 4
s4 = list(itertools.combinations(D, 4))
for S in [s4[i] for i in rng.choice(len(s4), size=min(N4, len(s4)), replace=False)]:
    e = q.V(S)
    order = np.argsort(e.A, kind="stable")
    rat = np.where((e.Q > 0) & ~e.Z)[0]
    near1 = np.where((np.abs(e.A) < 1e-15) & ~e.Z)[0]
    A_, cts = np.unique(e.A, return_counts=True)
    picks = list(rng.choice(len(e), size=NE, replace=False)) + list(rng.choice(rat, size=min(NE, len(rat)), replace=False)) \
        + list(rng.choice(near1, size=min(NE, len(near1)), replace=False))
    seen_w = set()
    for i in picks:
        i = int(i)
        if e.Z[i]:
            continue
        lo, hi = window_of(e, i, widen=2)
        if (lo, hi) in seen_w:
            continue
        seen_w.add((lo, hi))
        vals = bq(S, lo, hi)
        for (a, w) in vals:
            check_value(S, e, a, w, order)                         # (2)
            l, h = motorx.ln_bounds(a)                              # (4)
            cnt["lnb"] += 1
            la = ln_of(w)
            if not (mpf(l) <= la <= mpf(h)):
                bad += 1; print("LN_BOUNDS misses ln of", w, l, h, mpmath.nstr(la, 25))
        # (3) what the entry i keeps and removes
        if e.Q[i] > 0:
            continue
        lo, hi = window_of(e, i)
        old = bq(S, lo, hi); new = q.members(S, e, i)
        newk = {v.key() for (v, _) in new}
        if not newk <= {v.key() for (v, _) in old}:
            bad += 1; print("MEMBERS not among motor B's values:", S, i)
        for (v, w) in old:
            if v.key() in newk:
                continue
            cnt["removed"] += 1
            la = ln_of(w)
            if mpf(float(e.A[i])) <= la <= mpf(float(e.B[i])):
                bad += 1; print("REMOVED a value whose ln is in [A, B]:", S, i, w, float(e.A[i]), float(e.B[i]))
    print(f"V{''.join(map(str, S))}: {len(e)} entries, {len(seen_w)} windows; totals: "
          f"{cnt['values']} values checked, {cnt['removed']} removed values, problems {bad}", flush=True)

print(f"t_members: {cnt['values']} values with invariant I, {cnt['removed']} removed values checked, "
      f"{cnt['lnb']} ln_bounds, {cnt['rows']} raw rows (dedup); motor B unresolved {q.motorb_unresolved()}; "
      f"PROBLEMS {bad}")
sys.exit(1 if bad or q.motorb_unresolved() else 0)
