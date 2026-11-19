"""lvl.py, exact linear forms (13/11: lvl.lin_of, add/sub/mul/div with the forms; 14/11: quotients of proportional
forms are the exact ratio, and a form reduced to one term is the exact power).
The residual hits of n = 9, split 14 (pieces VN.2 and VN.5) are sums in which exact irrational powers cancel:
(4^-.9) + (5 - (8^-.6)) = 2^-1.8 + 5 - 2^-1.8 = 5, then |3 - 5| = 2, |2 - 2| = 0 and 0^-(10^-7) has no value.
Random linear combinations  L = c0 + sum_i a_i B_i^(e_i)  (B_i small rationals, e_i non-integer rationals, a_i nonzero
rationals of both signs, c0 large enough that every partial sum is positive), each built twice with lvl's operations:
  - T1 and T2: the terms in different orders and the powers written differently ((B^k)^(e/k), B^(e+1)/B, B*B^(e-1),
    a*m as m/(1/a)), c0 split in two: |T1 - T2| must be the exact 0, and T1 must be the exact rational when L is one;
  - T3: T2 with one term perturbed (its coefficient or its exponent): |T1 - T3| must never be the exact 0;
  - products (T1 * T2') and quotients by a monomial: compared with the truth;
  - every value lvl gives (exact or an enclosure) must contain the true value (mpmath, 300 digits, the same magnitude
    semantics: |a - b| at every difference).
usage: t_lvllin.py [trials] [seed]"""
import sys, os, numpy as np, mpmath
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lvl
N = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 17)
E = lvl.exact
mpmath.mp.dps = 300
BASES = [Fr(2), Fr(3), Fr(5), Fr(6), Fr(1, 2), Fr(2, 3), Fr(3, 2), Fr(10), Fr(1, 3)]
COEF = [Fr(1), Fr(2), Fr(3), Fr(1, 2), Fr(1, 4), Fr(3, 10), Fr(5, 2), Fr(7)]
TOL = mpmath.mpf(10) ** -250


def mq(q):
    return mpmath.mpf(q.numerator) / q.denominator


class V:
    """a value built twice: with lvl (v) and with mpmath (t, the same magnitudes)"""
    def __init__(self, v, t):
        self.v, self.t = v, t


def lit(q):
    return V(E(q), mq(q))


def power(B, e, way):
    """B^e written in one of several ways"""
    if way == 1 and B.denominator == 1 and B.numerator in (2, 3) :
        k = int(rng.integers(2, 4))                         # (B^k)^(e/k)
        x = e / k
        return V(lvl.pw(E(B ** k), E(abs(x)), 1 if x > 0 else -1), mpmath.power(mq(B), mq(e)))
    if way == 2:                                            # B^(e+1) / B
        x = e + 1
        return V(lvl.div(lvl.pw(E(B), E(abs(x)), 1 if x > 0 else -1), E(B)), mpmath.power(mq(B), mq(e)))
    if way == 3:                                            # B * B^(e-1)
        x = e - 1
        return V(lvl.mul(E(B), lvl.pw(E(B), E(abs(x)), 1 if x > 0 else -1)), mpmath.power(mq(B), mq(e)))
    return V(lvl.pw(E(B), E(abs(e)), 1 if e > 0 else -1), mpmath.power(mq(B), mq(e)))


def scaled(a, m, way):
    """|a| * m"""
    a = abs(a)
    if way == 1 and a != 1:
        return V(lvl.div(m.v, E(1 / a)), m.t * mq(a))
    return V(lvl.mul(E(a), m.v), m.t * mq(a))


def plus(x, y):
    return V(lvl.add(x.v, y.v), x.t + y.t)


def minus(x, y):
    return V(lvl.sub(x.v, y.v), abs(x.t - y.t))


def build(c0, terms, rand_ways):
    """c0 + sum a_i B_i^(e_i), terms in the given order"""
    if rand_ways and rng.random() < 0.5:
        c1 = Fr(int(rng.integers(1, 4)), 10) * c0             # at most 0.3 c0: every partial sum stays positive
        acc = lit(c0 - c1)
        rest = [("c", c1)]
    else:
        acc = lit(c0)
        rest = []
    items = [("t", t) for t in terms] + rest
    if rand_ways:
        items = [items[i] for i in rng.permutation(len(items))]
    for kind, it in items:
        if kind == "c":
            acc = plus(acc, lit(it))
            continue
        a, B, e = it
        m = power(B, e, int(rng.integers(0, 4)) if rand_ways else 0)
        am = scaled(a, m, int(rng.integers(0, 2)) if rand_ways else 0)
        acc = plus(acc, am) if a > 0 else minus(acc, am)
    return acc


def check(x, what):
    """the value lvl gives must contain the truth"""
    global bad, exacts
    v, t = x.v, x.t
    if v.q is not None:
        exacts += 1
        if abs(mq(v.q) - t) > TOL * max(1, abs(t)):
            bad += 1; print("WRONG EXACT", what, v.q, mpmath.nstr(t, 30)); return False
        return True
    try:
        if v.x is not None:
            ok = v.x.a - TOL <= t <= v.x.b + TOL
        else:
            L = lvl.lam_iv(v); ok = L.a - TOL <= mpmath.log(t) <= L.b + TOL
    except lvl.Fail:
        return True
    if not ok:
        bad += 1; print("ENCLOSURE MISSES", what, v, mpmath.nstr(t, 30))
    return ok


bad = exacts = zeros = zeros_want = fakes = pert = 0
for trial in range(N):
    k = int(rng.integers(1, 5))
    terms = []
    for _ in range(k):
        B = BASES[int(rng.integers(len(BASES)))]
        e = Fr(int(rng.integers(1, 40)), int(rng.choice([2, 3, 5, 10]))) * (1 if rng.random() < 0.6 else -1)
        if e.denominator == 1:
            e += Fr(1, 2)
        a = COEF[int(rng.integers(len(COEF)))] * (1 if rng.random() < 0.5 else -1)
        terms.append((a, B, e))
    tot = sum(abs(mq(a)) * mpmath.power(mq(B), mq(e)) for (a, B, e) in terms)
    if tot > 1e6:
        continue
    c0 = Fr(2 * int(mpmath.ceil(tot)) + int(rng.integers(1, 5)), 1) * Fr(int(rng.choice([1, 3, 7])), 1) + Fr(int(rng.integers(0, 10)), 10)
    try:
        T1 = build(c0, terms, False)
        T2 = build(c0, terms, True)
    except lvl.Fail:
        continue
    check(T1, "T1"); check(T2, "T2")
    try:
        d = minus(T1, T2)
        zeros_want += 1
        if d.v.q is not None and d.v.q == 0:
            zeros += 1
        check(d, "T1-T2")
    except lvl.Fail:
        pass
    # a perturbed copy: never an exact 0
    i = int(rng.integers(len(terms)))
    a, B, e = terms[i]
    t3 = list(terms)
    t3[i] = (a + Fr(1, 7), B, e) if rng.random() < 0.5 else (a, B, e + Fr(1, 1000))
    if t3[i][0] == 0:
        continue
    try:
        T3 = build(c0, t3, True)
        d3 = minus(T1, T3)
        pert += 1
        if d3.v.q is not None and d3.v.q == 0:
            fakes += 1; bad += 1; print("FALSE EXACT ZERO", terms, t3)
        check(d3, "T1-T3")
        check(V(lvl.mul(T1.v, T3.v), T1.t * T3.t), "T1*T3")
        m = power(B, e, 0)
        check(V(lvl.div(T1.v, m.v), T1.t / m.t), "T1/m")
        check(plus(T1, T3), "T1+T3")
    except lvl.Fail:
        pass
# (14/11) proportional forms: (r L) / L = r exactly; a perturbed copy is never an exact rational
prop = prop_bad_free = 0
for trial in range(N // 3):
    k = int(rng.integers(1, 4)); terms = []
    for _ in range(k):
        B = BASES[int(rng.integers(len(BASES)))]
        e = Fr(int(rng.integers(1, 30)), int(rng.choice([2, 3, 5, 10]))) * (1 if rng.random() < 0.6 else -1)
        if e.denominator == 1:
            e += Fr(1, 2)
        terms.append((COEF[int(rng.integers(len(COEF)))] * (1 if rng.random() < 0.5 else -1), B, e))
    tot = sum(abs(mq(a)) * mpmath.power(mq(B), mq(e)) for (a, B, e) in terms)
    if tot > 1e6:
        continue
    c0 = Fr(2 * int(mpmath.ceil(tot)) + 1)
    r = Fr(int(rng.integers(1, 9)), int(rng.choice([1, 2, 3, 7])))
    try:
        T1 = build(c0, terms, False)
        T2 = build(c0 * r, [(a * r, B, e) for (a, B, e) in terms], True)
        q = lvl.div(T2.v, T1.v)
        prop += 1
        if q.q != r:
            bad += 1; print("PROPORTIONAL NOT EXACT", r, terms, q)
        i = int(rng.integers(len(terms)))
        t3 = [(a * r + (Fr(1, 7) if j == i else 0), B, e) for j, (a, B, e) in enumerate(terms)]
        if t3[i][0] != 0:
            T3 = build(c0 * r, t3, True)
            q3 = lvl.div(T3.v, T1.v)
            if q3.q is not None:
                bad += 1; print("FALSE EXACT RATIO", terms, t3)
            check(V(q3, T3.t / T1.t), "T3/T1")
    except lvl.Fail:
        pass
# (14/11) a sum reduced to one term a B^e is the exact power: (a B^e)^q exact when q e is an integer
mono = 0
for trial in range(N // 3):
    B = BASES[int(rng.integers(len(BASES)))]
    e = Fr(int(rng.integers(1, 12)), int(rng.choice([2, 3, 5])))
    if e.denominator == 1:
        continue
    a1, a2 = Fr(int(rng.integers(1, 5))), Fr(int(rng.integers(1, 5)), int(rng.choice([1, 2, 3])))
    try:
        m1 = power(B, e, 0); m2 = power(B, e, int(rng.integers(1, 4)))
        x = plus(scaled(a1, m1, 0), scaled(a2, m2, 1))
        k = e.denominator
        y = lvl.pw(x.v, E(k), 1)
        want = (a1 + a2) ** k * B ** int(e * k)
        mono += 1
        if y.q != want:
            bad += 1; print("SINGLE TERM POWER NOT EXACT", B, e, a1, a2, y, want)
    except lvl.Fail as ex:
        bad += 1; print("FAIL single term", B, e, ex)
# (14/11) huge exact rationals in the forms (split 117): |X - (c + |m - X|)| = |m - c| with X = B^k of tens of
# thousands of bits and m a small exact power; the enclosure must contain |m - c| (it fails without the exact form)
big = 0
for trial in range(N // 6):
    B = [Fr(3), Fr(9), Fr(5), Fr(7)][int(rng.integers(4))]
    k = int(rng.integers(5000, 40000))
    Bm = BASES[int(rng.integers(len(BASES)))]
    e = Fr(int(rng.integers(1, 10)), int(rng.choice([2, 5, 10]))) * (1 if rng.random() < 0.5 else -1)
    if e.denominator == 1:
        continue
    c = Fr(int(rng.integers(1, 20)), 10)
    try:
        X = V(lvl.pw(E(B), E(k), 1), None)
        m = power(Bm, e, int(rng.integers(0, 4)))
        y = lvl.sub(X.v, lvl.add(E(c), lvl.sub(m.v, X.v)))
        truth = abs(m.t - mq(c))
        if truth < mpmath.mpf(10) ** -30:
            continue
        big += 1
        check(V(y, truth), "|X - (c + |m - X|)|")
    except lvl.Fail as ex:
        bad += 1; print("FAIL big constant form", B, k, Bm, e, c, ex)
print(f"t_lvllin (14/11): {prop} proportional quotients exact, {mono} single-term powers exact, {big} forms with a huge constant")
# the expressions of the lab (split 14): no value
s = lvl.pw(E(4), E(Fr(9, 10)), -1); X = lvl.pw(E(8), E(Fr(3, 5)), -1)
z = lvl.add(s, lvl.sub(E(5), X))                                # (4^-.9) + (5 - (8^-.6)) = 5
w = lvl.sub(E(2), lvl.sub(E(3), z))
try:
    lvl.pw(w, E(Fr(1, 10 ** 7)), -1); bad += 1; lab = "A VALUE"
except lvl.Invalid:
    lab = "no value"
z2 = lvl.add(lvl.pw(E(4), E(Fr(9, 10)), -1), lvl.sub(E(Fr(3, 10)), X))   # (4^-.9) + (.3 - (8^-.6)) = .3
w2 = lvl.sub(E(Fr(1, 2)), lvl.add(E(Fr(1, 5)), z2))
try:
    lvl.pw(w2, E(Fr(1, 10 ** 7)), -1); bad += 1; lab += ", A VALUE"
except lvl.Invalid:
    lab += ", no value"
print(f"t_lvllin: {zeros}/{zeros_want} equal pairs give the exact 0, {pert} perturbed pairs ({fakes} false zeros), "
      f"{exacts} exact results checked; lab expressions: {lab}; PROBLEMS {bad}")
sys.exit(1 if bad or zeros < zeros_want else 0)
