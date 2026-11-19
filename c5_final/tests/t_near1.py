"""near1(S), the exact members of V(S) with |ln v| <= DELTA1 (motor B, qx.near1), against motor X's own recursion:
  (1) for random windows [A, B] inside [-DELTA1, DELTA1] (around 0, one-sided, tiny and wide), every exact value
      that exact_V finds with the near1 shortcut OFF (Theorem X recursion with exact leaves) and whose ln may lie in
      the window is covered by a row of near1(S) meeting the window (its Arb ball of ln v overlaps the row), and
      the exact 1 is a row with the flag O;
  (2) IndexN1.look returns exactly the pairs of the plain test (B >= lo and A <= hi) on the rows;
  (3) every row: the sign and the second level agree with a fresh Arb evaluation of its writing.
usage: t_near1.py n digits [windows] [seed]      e.g. t_near1.py 5 12345 40"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, qx, ivl
import exact as ex
from flint import arb, ctx
n = int(sys.argv[1]); S = tuple(int(c) for c in sys.argv[2])
NW = int(sys.argv[3]) if len(sys.argv) > 3 else 40
rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 9)
q = qx.QX(n, 0, os.path.join(paths.SRC, f"res_{n}_0.txt"), 4, 5)
e = q.near1(S)
assert e is not None, "near1 not available (motor B left operations undecided)"
rows = q._near1_rows(S)
D = qx.DELTA1
print(f"near1{S}: {len(e)} rows (exact 1: {int(e.O.sum())}, undecided sign: {int((e.S == 0).sum())}), radius {D:g}")
bad = 0

# (3) rows against a fresh evaluation of their writing
ctx.prec = 512
for k in rng.choice(len(rows), size=min(len(rows), 400), replace=False):
    v, w = q.near1_value(S, int(k))
    if e.O[k]:
        if ex.sign(ex.sub(v, ex.ONE)) != 0:
            bad += 1; print("ROW flagged 1 is not 1:", w)
        continue
    b = ex.tight_ball(v).log()
    if not (e.A[k] <= b.mid() <= e.B[k]) and not b.contains(0):
        if not (b.upper() >= e.A[k] and b.lower() <= e.B[k]):
            bad += 1; print("ROW enclosure misses ln v:", w, e.A[k], e.B[k], b)
    if e.S[k] != 0 and not ((b > 0 and e.S[k] == 1) or (b < 0 and e.S[k] == -1) or b.contains(0)):
        bad += 1; print("ROW sign wrong:", w, e.S[k], b)
    if e.S[k] != 0 and not b.contains(0):
        lb = abs(b).log()
        if not (lb.upper() >= e.M[k] and lb.lower() <= e.N[k]):
            bad += 1; print("ROW second level misses ln|ln v|:", w, e.M[k], e.N[k], lb)
ctx.prec = 53

# (2) the index against the plain test
for t in range(200):
    c = float(rng.uniform(-D, D)) * float(10.0 ** -rng.integers(0, 300))
    w = abs(float(rng.normal(0, 1))) * float(10.0 ** -rng.integers(0, 320)) * D
    lo = np.array([c - w, -np.inf if rng.random() < 0.1 else c - 2 * w]); hi = np.array([c + w, c + 3 * w])
    ix = qx.IndexN1(e)
    W, K = ix.look(lo, hi)
    got = set(zip(W.tolist(), K.tolist()))
    ref = set()
    for j in range(len(lo)):
        m = np.where((e.B >= lo[j]) & (e.A <= hi[j]))[0]
        ref |= {(j, int(i)) for i in m}
    if got != ref:
        bad += 1; print("INDEX differs", lo, hi, len(got), len(ref))

# (1) the recursion (near1 shortcut off) against the rows
qx.NEAR1 = False
checked = covered = 0
for t in range(NW):
    kind = rng.integers(0, 4)
    if kind == 0:                    # around 0
        r = D * float(10.0 ** -rng.uniform(0, 12)); A, B = -r, r * float(rng.uniform(0.1, 1))
    elif kind == 1:                  # positive side
        a = D * float(10.0 ** -rng.uniform(0, 30)); A, B = a, min(D, a * float(rng.uniform(1, 50)))
    elif kind == 2:                  # negative side
        a = D * float(10.0 ** -rng.uniform(0, 30)); A, B = -min(D, a * float(rng.uniform(1, 50))), -a
    else:                            # the whole radius
        A, B = -D, D
    vals = q.exact_V(S, A, B)
    for (v, w) in vals:
        try:
            if ex.sign(ex.sub(v, ex.ONE)) == 0:
                if not (A <= 0 <= B):
                    continue             # exact_V returns a superset: this 1 is not in the window
                checked += 1
                if not (e.O & (e.A <= B) & (e.B >= A)).any():
                    bad += 1; print("MISSED exact 1", w)
                else:
                    covered += 1
                continue
        except ex.Unresolved:
            pass
        old = ctx.prec; ctx.prec = 1024
        try:
            lb = ex.tight_ball(v).log()
        finally:
            ctx.prec = old
        if lb.upper() < A or lb.lower() > B:
            continue                 # not in the window
        checked += 1
        lo_v, hi_v = float(lb.lower()), float(lb.upper())
        tol = max(abs(lo_v), abs(hi_v)) * 1e-15 + 1e-320          # symmetric slack (both signs of ln v)
        m = (e.A <= hi_v + tol) & (e.B >= lo_v - tol) & (e.B >= A) & (e.A <= B)
        if not m.any():
            bad += 1; print("MISSED", w, "ln v ~", lb.mid(), "window", A, B)
        else:
            covered += 1
print(f"windows {NW}: exact values of the recursion in the windows {checked}, covered by rows {covered}; problems {bad}")
