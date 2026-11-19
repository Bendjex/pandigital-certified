"""Test 2: motor B's window query (preimage windows + recursion) against complete enumeration.
For 4-digit sets (any order) and 4-5 digit runs (in order): for many windows, the values returned by
query() must be exactly the values of full() that may lie in the window."""
import sys, random
sys.path.insert(0, ".")
import mpmath
from flint import arb
import motorb
mpmath.mp.dps = 40
random.seed(12345)

def k(x):
    if x.is_zero(): return "0"
    return mpmath.nstr(mpmath.mpf(abs(x.ball(256)).log().mid().str(50, radius=False)), 25)

def windows(vals):
    out = []
    fl = [float(v.ball(64).mid()) for v in vals]
    fl = [f for f in fl if f == f]
    for _ in range(25):                                   # around attainable values, several widths
        f = random.choice(fl)
        for w in (1e-15, 1e-9, 1e-3):
            out.append((f * (1 - w), f * (1 + w)))
    for _ in range(10):                                   # log-uniform random windows
        f = 10 ** random.uniform(-40, 40)
        out.append((f, f * (1 + 10 ** random.uniform(-12, 0))))
    for c in (0.0, 1.0, 2.0, 0.5, 3.0, 10.0, 0.1):        # special points
        out.append((max(0.0, c - 1e-12), c + 1e-12))
        out.append((max(0.0, c - 1e-3), c + 1e-3))
    out.append((0.0, 1e-30))
    out.append((1e30, 1e300))
    return out

total_bad = 0
cases = [(False, ("1","2","3","4")), (False, ("2","5","7","9")), (False, ("1","3","6","8")),
         (True, ("1","2","3","4")), (True, ("5","6","7","8")), (True, ("1","2","3","4","5"))]
for ordered, S in cases:
    bf = motorb.B(ordered=ordered, fullmax=len(S))
    bq = motorb.B(ordered=ordered, fullmax=3)
    F = bf.full(S)
    vals = [v for (v, e) in F.values()]
    nb = 0
    for (a, b) in windows(vals):
        lo, hi = arb(a), arb(b)
        want = {k(v) for (v, e) in F.values() if motorb.B.may_lie(v.ball(256), lo, hi)}
        got = {k(v) for (v, e) in bq.query(S, lo, hi).values()}
        if want != got:
            nb += 1
            if nb <= 3:
                print("   MISMATCH window", a, b, "missing", list(want - got)[:3], "extra", list(got - want)[:3])
    total_bad += nb
    print(("ord " if ordered else "any ") + "".join(S), "values", len(F), "windows", len(windows(vals)),
          "mismatches", nb, "unresolved", len(bf.unresolved) + len(bq.unresolved), flush=True)
print("TOTAL MISMATCHES", total_bad)
