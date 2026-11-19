"""(03/11) exact_V(S, A, B, c2): the exact members of V(S) in a window with a second-level constraint c2 = (sign of ln v,
lower, upper of ln|ln v|), the query that refine_base1 makes for the base l of l^s, l^-s (QX_NARROW_REFINE): the window
is the base's enclosure intersected with its preimage window (Lemma X3) and c2 is the sub-query's constraint
(pre_level2).  The claim tested: every member v of V(S) with ln v in [A, B], sign(ln v) = sign and ln|ln v| in
[lower, upper] is among the values returned (the recursion with c2 only drops values that fail the constraint).
 (1) |S| = 5, against motor B: the rows of near1(S) (every exact member of V(S) with |ln v| <= 1e-11, motor B, Lemma 14).
     Windows and constraints are built around random rows (tight, with slack, touching the row's ends, and random
     ones inside the radius).  Every row whose Arb ball of ln v lies inside the window and whose sign and second
     level satisfy the constraint must be matched by a returned value (their 1024-bit balls of ln v overlap).  The
     near-1 shortcut of Q is off here, so the recursion of Theorem X runs (with the constraint passed down).
 (2) |S| = 6, against the same recursion without the constraint: for random windows [A, B] around values of V(S), the
     values of exact_V(S, A', B', c2) (A <= A' <= B' <= B, c2 from a returned value's second level, with slack)
     contain every value of exact_V(S, A, B) whose ball lies inside [A', B'] and satisfies c2.
usage: t_narrow.py [windows per set] [seed] [sets5 comma-separated] [set6]
       e.g. t_narrow.py 30 1 45678,56789 456789      (n = 9 tables: QX_NCACHE)"""
import sys, os, bisect, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, qx, ivl
import exact as ex
from flint import arb, ctx
NW = int(sys.argv[1]) if len(sys.argv) > 1 else 30
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
SETS5 = [tuple(int(c) for c in s) for s in (sys.argv[3] if len(sys.argv) > 3 else "45678,56789").split(",") if s]
SET6 = tuple(int(c) for c in (sys.argv[4] if len(sys.argv) > 4 else "456789")) if (len(sys.argv) <= 4 or sys.argv[4]) else None
q = qx.QX(9, 0, os.path.join(paths.SRC, "res_9_0.txt"), 4, 5)
bad = 0
INF = np.inf


def lnball(v, prec=1024):
    old = ctx.prec; ctx.prec = prec
    try:
        return ex.tight_ball(v).log()
    finally:
        ctx.prec = old


def balls(vals):
    """(sorted midpoints, list of (lo, hi, is_one)) of the ln balls of the returned values"""
    out = []
    for (v, w) in vals:
        try:
            if ex.sign(ex.sub(v, ex.ONE)) == 0:
                out.append((0.0, 0.0, True)); continue
        except ex.Unresolved:
            pass
        b = lnball(v)
        out.append((float(b.lower()), float(b.upper()), False))
    out.sort()
    return out


def covered(bl, lo, hi, one=False):
    """some returned ball overlaps [lo, hi] (exact 1: some returned exact 1)"""
    if one:
        return any(o for (_, _, o) in bl)
    i = bisect.bisect_left(bl, (lo - max(abs(lo), abs(hi)) * 1e-12 - 1e-300, -INF, False))
    for (a, b, o) in bl[max(0, i - 2):]:
        if a > hi:
            break
        if b >= lo and a <= hi:
            return True
    return False


def c2_ok(sg, m, M, c2):
    """the row's sign and second level [m, M] certainly satisfy the constraint (None: no constraint)"""
    if c2 is None:
        return True
    return sg == c2[0] and m >= c2[1] and M <= c2[2]


# (1) 5-digit sets against the near1 rows (motor B)
D = qx.DELTA1
for S in SETS5:
    e = q.near1(S)
    assert e is not None, "near1 not available"
    rows = q._near1_rows(S)
    nz = np.where(~e.O & (e.S != 0))[0]
    print(f"near1{S}: {len(e)} rows ({len(nz)} with a known sign and second level)", flush=True)
    qx.NEAR1 = False
    req_tot = cov_tot = 0
    t0 = time.time()
    for t in range(NW):
        kind = t % 5
        k = int(rng.choice(nz))
        a, b = float(e.A[k]), float(e.B[k])
        w = max(abs(a), abs(b))
        if kind == 0:                                     # the row's ball exactly
            A, B = a, b
        elif kind == 1:                                   # relative slack
            A, B = a - w * float(10.0 ** -rng.uniform(0, 14)), b + w * float(10.0 ** -rng.uniform(0, 14))
        elif kind == 2:                                   # wide around the row, still inside the radius
            A, B = max(-D, a - w * float(rng.uniform(0, 3))), min(D, b + w * float(rng.uniform(0, 3)))
        elif kind == 3:                                   # one end at the row's end
            A, B = (a, b + w * float(rng.uniform(0, 1))) if rng.random() < 0.5 else (a - w * float(rng.uniform(0, 1)), b)
        else:                                             # random inside the radius, one-signed
            x = D * float(10.0 ** -rng.uniform(0, 25)) * (1 if rng.random() < 0.5 else -1)
            A, B = sorted((x, x * float(rng.uniform(1, 30))))
        # the constraint: the row's sign and second level with random slack (sometimes none, sometimes a
        # constraint that excludes the row: then nothing is required of it)
        r = rng.random()
        if r < 0.15:
            c2 = None
        elif r < 0.85:
            c2 = (int(e.S[k]), float(e.M[k]) - float(10.0 ** -rng.uniform(0, 12)) * (rng.random() < 0.7),
                  float(e.N[k]) + float(10.0 ** -rng.uniform(0, 12)) * (rng.random() < 0.7))
        else:
            c2 = (int(e.S[k]), float(e.N[k]) + 0.01, float(e.N[k]) + 1.0)       # excludes the row k
        vals = q.exact_V(S, A, B, c2)
        bl = balls(vals)
        for j in range(len(e)):
            if e.O[j]:
                if A <= 0 <= B and c2 is None:
                    req_tot += 1
                    if covered(bl, 0.0, 0.0, one=True):
                        cov_tot += 1
                    else:
                        bad += 1; print("MISSED exact 1", rows[j][6], (A, B), flush=True)
                continue
            if not (e.A[j] >= A and e.B[j] <= B):        # the row's ball is not inside the window
                continue
            if e.S[j] == 0 or not c2_ok(int(e.S[j]), float(e.M[j]), float(e.N[j]), c2):
                continue
            req_tot += 1
            v, wr = q.near1_value(S, j)
            bj = lnball(v)
            if covered(bl, float(bj.lower()), float(bj.upper())):
                cov_tot += 1
            else:
                bad += 1; print("MISSED", wr, "ln v ~", bj.mid(), "window", A, B, "c2", c2, flush=True)
    qx.NEAR1 = True
    print(f"  {NW} windows: required rows {req_tot}, matched {cov_tot}, {time.time() - t0:.0f}s", flush=True)

# (2) a 6-digit set: the constrained narrowed query against the unconstrained one
if SET6:
    t0 = time.time()
    req_tot = cov_tot = 0
    nwin = max(2, NW // 5)
    for t in range(nwin):
        x = float(rng.uniform(-6, 6))
        if abs(x) < 1e-3:
            x = 0.5
        w = abs(x) * float(10.0 ** -rng.uniform(12, 14.5))
        A, B = x - w, x + w
        full = q.exact_V(SET6, A, B)
        bl_full = [(float(b.lower()), float(b.upper())) for b in (lnball(v) for (v, _) in full)]
        inside = [i for i, (a, b) in enumerate(bl_full) if a >= A and b <= B]
        if not inside:
            print(f"  window {t}: no value of V{SET6} inside [{A}, {B}]", flush=True)
            continue
        i = int(rng.choice(inside))
        a, b = bl_full[i]
        A2 = max(A, a - (a - A) * float(rng.uniform(0, 1))); B2 = min(B, b + (B - b) * float(rng.uniform(0, 1)))
        sg = 1 if a > 0 else -1
        lb = abs(lnball(full[i][0])).log()
        c2 = (sg, float(lb.lower()) - float(10.0 ** -rng.uniform(0, 12)), float(lb.upper()) + float(10.0 ** -rng.uniform(0, 12)))
        nar = q.exact_V(SET6, A2, B2, c2)
        bl = balls(nar)
        for j, (a_, b_) in enumerate(bl_full):
            if not (a_ >= A2 and b_ <= B2):
                continue
            lbj = abs(lnball(full[j][0])).log()
            if not ((a_ > 0) == (sg > 0) and float(lbj.lower()) >= c2[1] and float(lbj.upper()) <= c2[2]):
                continue
            req_tot += 1
            if covered(bl, a_, b_):
                cov_tot += 1
            else:
                bad += 1; print("MISSED (6 digits)", full[j][1], "window", A2, B2, "c2", c2, flush=True)
        print(f"  window {t}: [{A:.6g} +- {w:.2g}]: {len(full)} values, narrowed {len(nar)}", flush=True)
    print(f"{SET6}: {nwin} windows: required values {req_tot}, matched {cov_tot}, {time.time() - t0:.0f}s", flush=True)
print(f"t_narrow: problems {bad}")
