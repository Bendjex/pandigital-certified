"""Reviewer's structural coverage test of motor X (the role of the missing t_coinc.py, done at random):
random EXPRESSION TREES over a digit set S are evaluated node by node with lvl (the package's independent
evaluator, exact rationals where possible).  If some node is certainly out of range ([1e-300, 1e300] u {0}), the
magnitude of the tree is a value of X(S) (structural definition, C.2), whatever its value is (it may coincide with
an in-range value).  Then
  (1) some entry of the materialised table N^(S) (motorx.N) must contain it, and
  (2) the window query Q(S, "N", w) with a narrow window w around it must return an enclosure containing it.
usage: struct_cover.py PC n ord tv tn seed ntrees digits [query]"""
import sys, os, random, time
PC = sys.argv[1]
sys.path.insert(0, os.path.join(PC, "c5")); sys.path.insert(1, os.path.join(PC, "src"))
import numpy as np, mpmath
from mpmath import mp, iv
from fractions import Fraction as Fr
import ivl, motorx, lvl, qx
mp.prec = 600
n, o, tv, tn, seed, NT = [int(x) for x in sys.argv[2:8]]
S = tuple(int(c) for c in sys.argv[8])
DOQ = len(sys.argv) > 9 and sys.argv[9] == "query"
rnd = random.Random(seed)
q = qx.QX(n, o, os.path.join(PC, "src", f"res_{n}_{o}.txt"), tv, tn)

LO, HI = Fr(1, 10 ** 300), Fr(10 ** 300)
LLO, LHI = iv.log(iv.mpf(10)) * (-300), iv.log(iv.mpf(10)) * 300

def rng_class(v):
    """'in', 'out' or '?' for a lvl magnitude (0 counts as in range)"""
    if v.far:
        return "out" if v.far in ("huge", "tiny") else "in"
    if v.q is not None:
        if v.q == 0: return "in"
        return "in" if LO <= v.q <= HI else "out"
    try:
        lam = lvl.lam_iv(v)
    except lvl.Fail:
        if v.s is not None and v.Lam is not None and v.Lam.a > 7:
            return "out"
        return "?"
    if lam.a > LHI.b or lam.b < LLO.a: return "out"
    if lam.b < LHI.a and lam.a > LLO.b: return "in"
    return "?"

def literals(T):
    out = []
    import itertools
    perms = [T] if o else itertools.permutations(T)
    for p in perms:
        t = "".join(map(str, p))
        for i in range(len(t) + 1):
            out.append(Fr(int(t), 10 ** (len(t) - i)))
    return out

def rand_tree(T):
    """random derivation over the digit tuple T (magnitude forms of Lemma 3)"""
    if len(T) == 1 or rnd.random() < 0.12:
        return ("L", rnd.choice(literals(T)))
    if o:
        i = rnd.randint(1, len(T) - 1)
        S1, S2 = T[:i], T[i:]
        forms = ivl.ORD_LEFT
    else:
        while True:
            S1 = tuple(d for d in T if rnd.random() < 0.5)
            if 0 < len(S1) < len(T): break
        S2 = tuple(d for d in T if d not in S1)
        forms = tuple(range(9))
    w = {0: 2, 1: 3, 2: 2, 3: 2, 4: 2, 5: 5, 6: 5, 7: 5, 8: 5}
    f = rnd.choices(forms, weights=[w[k] for k in forms])[0]
    return ("O", f, rand_tree(S1), rand_tree(S2))

class Skip(Exception):
    pass

def ev(t):
    """(lvl value, has an out-of-range node)"""
    if t[0] == "L":
        return lvl.exact(t[1]), False
    _, f, a, b = t
    va, oa = ev(a); vb, ob = ev(b)
    v = lvl.op(lvl.FORMS[f], va, vb)
    c = rng_class(v)
    if c == "?":
        raise Skip()
    return v, (oa or ob or c == "out")

def show(t):
    if t[0] == "L": return str(t[1])
    return f"[{ivl.FORMS[t[1]]} | {show(t[2])} | {show(t[3])}]"

def contains(e, i, v):
    """may the enclosure i of e contain the lvl value v?  (main statement A <= ln u <= B, flags)"""
    A, B = float(e.A[i]), float(e.B[i])
    if v.q is not None and v.q == 0:
        return bool(e.Z[i]) or A == -np.inf
    if e.Z[i]: return False
    if v.q is not None and v.q == 1:
        return bool(e.O[i]) or (A <= 0 <= B)
    if e.O[i]: return False
    try:
        lam = lvl.lam_iv(v)
        return not (lam.b < A or lam.a > B)
    except lvl.Fail:
        if v.far == "huge" or (v.s == 1 and v.Lam is not None): return B == np.inf
        if v.far == "tiny" or (v.s == -1 and v.Lam is not None): return A == -np.inf
        return True

def window(v):
    """narrow ln-window certainly containing the value (None if not available)"""
    if v.far or (v.q is not None and v.q == 0): return None
    try:
        lam = lvl.lam_iv(v)
    except lvl.Fail:
        return None
    a, b = float(lam.a), float(lam.b)
    if not (np.isfinite(a) and np.isfinite(b)): return None
    return np.nextafter(np.nextafter(a, -np.inf), -np.inf), np.nextafter(np.nextafter(b, np.inf), np.inf)

t0 = time.time()
Nt = q.N(S) if len(S) <= tn else None
print("N table of", S, ":", None if Nt is None else len(Nt), "entries", round(time.time() - t0, 1), "s", flush=True)
stats = {"trees": 0, "skip": 0, "invalid": 0, "fail": 0, "X": 0, "X_in_range_value": 0, "X_exact_rational": 0, "X_zero": 0, "X_one": 0,
         "tab_missed": 0, "q_done": 0, "q_missed": 0}
missed = []
for it in range(NT):
    t = rand_tree(S)
    stats["trees"] += 1
    try:
        v, out = ev(t)
    except Skip:
        stats["skip"] += 1; continue
    except lvl.Invalid:
        stats["invalid"] += 1; continue
    except (lvl.Fail, OverflowError, MemoryError, ValueError, ZeroDivisionError):
        stats["fail"] += 1; continue
    if not out:
        continue
    stats["X"] += 1
    c = rng_class(v)
    if c == "in": stats["X_in_range_value"] += 1
    if v.q is not None:
        stats["X_exact_rational"] += 1
        if v.q == 0: stats["X_zero"] += 1
        if v.q == 1: stats["X_one"] += 1
    if Nt is not None:
        # candidates by the main statement (vectorised prefilter), then the exact test
        ok = False
        if v.q is not None and v.q == 0:
            ok = bool((Nt.Z | (Nt.A == -np.inf)).any())
        else:
            try:
                lam = lvl.lam_iv(v); a, b = float(lam.a), float(lam.b)
                cand = np.where((Nt.B >= np.nextafter(a, -np.inf)) & (Nt.A <= np.nextafter(b, np.inf)) & ~Nt.Z)[0]
            except lvl.Fail:
                cand = np.arange(len(Nt))
            for i in cand:
                if contains(Nt, int(i), v): ok = True; break
        if not ok:
            stats["tab_missed"] += 1
            if len(missed) < 20: missed.append(("table", show(t), repr(v)))
    if DOQ and c == "in" and (it % 1 == 0):
        w = window(v)
        if w is not None:
            stats["q_done"] += 1
            H = q.Q(S, "N", np.array([w[0]]), np.array([w[1]]))
            ok = any(contains(H, i, v) for i in range(len(H)))
            if not ok:
                stats["q_missed"] += 1
                if len(missed) < 20: missed.append(("query", show(t), repr(v), w))
    if it % 2000 == 0 and it: print("progress", it, round(time.time() - t0), stats, flush=True)
print("FINAL", stats, round(time.time() - t0, 1), "s")
for m in missed: print("MISSED", m)
