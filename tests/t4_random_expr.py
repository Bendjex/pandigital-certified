"""Test 4 (Theorem B end to end on windows that certainly contain values): random valid expressions
over the digits 1..n are generated, each window is centred on the absolute value of one of them, and
the values found by motor A's candidate paths (verified by motor B) must equal the values motor B finds
directly; the centre expression's own value must be among them.
usage: t4_random_expr.py n ord nwin seed [relmin relmax [full|centre]]
  full:   also compare with motor B's direct query (expensive for n >= 8)
  centre: only check that the centre expression's value is found through motor A's paths"""
import sys, os, random, subprocess, time, math
sys.path.insert(0, ".")
import mpmath
from flint import arb
import motorb, verify
import exact as ex

mpmath.mp.dps = 60
n, ordered, nwin, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
relmin = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-12
relmax = float(sys.argv[6]) if len(sys.argv) > 6 else 1e-8
MODE = sys.argv[7] if len(sys.argv) > 7 else "full"
random.seed(seed)
LO, HI = mpmath.mpf(10) ** -300, mpmath.mpf(10) ** 300


class Bad(Exception):
    pass


def ok(x):
    if x == 0:
        return x
    if not (LO <= abs(x) <= HI):
        raise Bad()
    return x


def gen(D):
    """random expression over the digit tuple D: (string, mpmath value)"""
    if len(D) == 1 or (len(D) <= 3 and random.random() < 0.35):
        perm = list(D) if ordered else random.sample(list(D), len(D))
        t = "".join(perm)
        i = random.randint(0, len(t))
        s = t if i == len(t) else t[:i] + "." + t[i:]
        return s, ok(mpmath.mpf(s))
    if ordered:
        k = random.randint(1, len(D) - 1)
        L, R = D[:k], D[k:]
    else:
        m = list(D)
        random.shuffle(m)
        k = random.randint(1, len(D) - 1)
        L, R = tuple(m[:k]), tuple(m[k:])
    (a, x), (b, y) = gen(L), gen(R)
    if random.random() < 0.3:
        b, y = "-" + b, -y
    op = random.choice("+-*/^^")
    if op == "+": v = x + y
    elif op == "-": v = x - y
    elif op == "*": v = x * y
    elif op == "/":
        if y == 0: raise Bad()
        v = x / y
    else:
        if x == 0:
            if y <= 0: raise Bad()
            v = mpmath.mpf(0)
        elif x < 0:
            raise Bad()                       # negative bases left out (Lemma 2 covers them)
        else:
            if abs(y * mpmath.log(x)) > 700: raise Bad()
            v = mpmath.power(x, y)
    return f"({a}{op}{b})", ok(v)


digs = tuple(str(d) for d in range(1, n + 1))
exprs = []
while len(exprs) < nwin:
    try:
        s, v = gen(digs)
    except (Bad, ZeroDivisionError):
        continue
    v = abs(v)
    if v == 0 or not (1e-290 < v < 1e290):
        continue
    if abs(mpmath.log(v)) < 1e-6:          # windows around 1 hold millions of values (x^tiny): too costly here
        continue
    exprs.append((s, v))
wins = []
for s, v in exprs:
    r = 10 ** random.uniform(math.log10(relmin), math.log10(relmax))
    a = float(v * (1 - r))
    b = float(v * (1 + r))
    wins.append((a, b))
res = f"res_{n}_{ordered}.txt"
if not os.path.exists(res) or os.path.getsize(res) == 0:
    subprocess.run(["./build_res.sh", str(n), str(ordered), res], check=True)
wf, cf = f"wins4_{n}_{ordered}.txt", f"cand4_{n}_{ordered}.txt"
with open(wf, "w") as f:
    for a, b in wins:
        f.write(f"{float(a).hex()} {float(b).hex()}\n")
t0 = time.time()
subprocess.run(["./motora", "multi", str(n), str(ordered), res, wf, cf], check=True,
               stderr=open(cf + ".log", "w"))
print(f"motor A: {time.time() - t0:.1f}s", flush=True)
cands = {}
for line in open(cf):
    w, rest = line.split(" ", 1)
    cands.setdefault(int(w[1:]), []).append(rest.strip())
V = verify.Verifier(n, ordered)


def k(x):
    if x.is_zero(): return "0"
    b = ex.tight_ball(x)
    return mpmath.nstr(mpmath.mpf(abs(b).log().mid().str(50, radius=False)), 25)


bad = 0
for i, (a, b) in enumerate(wins):
    lo, hi = arb(a), arb(b)
    t1 = time.time()
    VA = {}
    for c in cands.get(i, []):
        VA.update(V.verify_line(c, lo, hi))
    t2 = time.time()
    VB = V.B.query(digs, lo, hi) if MODE == "full" else VA
    t3 = time.time()
    ka = {k(v) for (v, e) in VA.values()}
    kb = {k(v) for (v, e) in VB.values()}
    centre = mpmath.nstr(mpmath.log(exprs[i][1]), 25)
    miss = kb - ka
    extra = ka - kb
    has_c = any(abs(float(z) - float(centre)) <= 1e-20 * max(1, abs(float(centre))) for z in kb) if kb else False
    if miss or extra or not has_c: bad += 1
    print(f"window {i} centre {exprs[i][0]} = {mpmath.nstr(exprs[i][1], 12)}: A cand {len(cands.get(i, []))}, "
          f"values A {len(ka)} B {len(kb)}, missing {len(miss)}, extra {len(extra)}, centre found {has_c} "
          f"({t2-t1:.1f}s + {t3-t2:.1f}s)", flush=True)
    for z in list(miss)[:3]:
        e = [e for (v, e) in VB.values() if k(v) == z]
        print("     MISSING", z, e[:2])
print(f"n={n} ord={ordered}: {len(wins)} windows, {bad} with differences; unresolved B {len(V.B.unresolved)}")
