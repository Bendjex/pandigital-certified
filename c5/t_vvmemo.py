"""resolve_vv with the differences kept per pair of lists (qx.vv_hull, 10/11) against the code it replaces (one
computation per pair of hits, copied below as resolve_vv_ref): on every call made by a real query, the same pairs
are kept (meet) and the same enclosures are left (r.A, r.B), bit for bit.
The query: one split S1 | S2 of the recursion for NE(S), S a set of digits, for 17 generic in-range windows (those
of the forms of V(1) around pi).
usage: t_vvmemo.py n digits_of_S digits_of_S1      e.g. t_vvmemo.py 9 23456789 57"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, qx, ivl
if os.environ.get("QX_ALLPROV") == "1":          # the tables of the provenance mode (as c5run.py)
    import hitcheck
    qx.DEBUG = True
from cells_x import window

n = int(sys.argv[1]); S = tuple(int(c) for c in sys.argv[2]); S1 = tuple(int(c) for c in sys.argv[3])
q = qx.QX(n, 0, os.path.join(os.environ.get("RESDIR", paths.SRC), f"res_{n}_0.txt"), 4, 5, verbose=False)


def resolve_vv_ref(self, S1, S2, I, sub, amb, meet, wlo, whi, r, want="out"):
    """|s - l| of two overlapping in-range classes (enclosure [-inf, B]), decided with the exact members (motor B,
    Lemma 14: the members of each class are a superset of its values).
      want "out" (pairs (V, V) of the top level): only an out-of-range difference is a new value (Lemma X1(c)); an
        exact zero or an in-range difference is dropped, an out-of-range one replaces the enclosure by its ball;
      want "in" (pairs of a recursive V query, windows clipped to the range): only an in-range difference is a
        member of V(S); an exact zero cannot meet a clipped window (0 is returned by zero_hits before the
        clipping), an out-of-range difference is not in V(S), and an in-range one replaces the enclosure by its ball.
    A difference whose sign stays undecided keeps [-inf, B]; a pair whose l is not a table class, or whose exact
    members cannot be obtained (Unresolved), keeps its enclosure."""
    import exact as ex
    from exact import Unresolved
    from flint import ctx, arb
    meet = meet.copy()
    e1 = self.V(S1)
    TH = arb(10) ** -299
    keep_fn = ivl.possibly_out if want == "out" else ivl.possibly_in
    # (03/11) recursive l values that overlap the same class of s have nearly the same narrow enclosure (e.g. many
    # writings of 4 against s = 4): one exact query on the union of their enclosures (a superset for each of them,
    # Lemma 14 / exact_V), whose members are then filtered by their own Arb ball of ln, replaces one query per pair
    pre_ys = {}
    if want == "out" and qx.GROUP_VV:
        rec = [t for t in np.where(amb)[0] if int(sub.ti[t]) < 0 and np.isfinite(sub.A[t]) and np.isfinite(sub.B[t])
               and float(sub.B[t]) - float(sub.A[t]) < 1e-6]
        by_cls = {}
        for t in rec:
            by_cls.setdefault(int(I[t]), []).append(t)
        for i_cls, ts in by_cls.items():
            if len(ts) < 2:
                continue
            A0 = float(min(sub.A[t] for t in ts)); B0 = float(max(sub.B[t] for t in ts))
            if not (B0 - A0 < 1e-6):
                continue
            try:
                vals = self.exact_vals(S2, A0, B0)
            except Unresolved:
                continue
            balls = []
            for (v, e) in vals:
                try:
                    old_p = ctx.prec; ctx.prec = 256
                    try:
                        b = self._ball1200(v)
                        lb = abs(b).log() if not b.contains(0) else None
                    finally:
                        ctx.prec = old_p
                except Exception:
                    lb = None
                if lb is None:
                    balls.append((v, e, -np.inf, np.inf))  # undecided ln: kept for every pair (superset)
                else:
                    balls.append((v, e, float(np.nextafter(float(lb.lower()), -np.inf)), float(np.nextafter(float(lb.upper()), np.inf))))
            for t in ts:
                pre_ys[t] = [(v, e) for (v, e, lo_, hi_) in balls if hi_ >= sub.A[t] and lo_ <= sub.B[t]]
            self.grouped_vv = getattr(self, "grouped_vv", 0) + len(ts)
    for t in np.where(amb)[0]:
        j = int(sub.ti[t])
        if j < 0 and want == "in":
            continue                                   # a recursive value (not a table class): kept as it is
        try:
            if j >= 0:
                ys = self.members(S2, self.V(S2), j)
            elif t in pre_ys:
                ys = pre_ys[t]
            else:
                ys = self.exact_vals(S2, float(sub.A[t]), float(sub.B[t]))     # a recursive V value: its exact members
            xs = self.members(S1, e1, int(I[t]))
        except Unresolved:
            self.unresolved_vv = getattr(self, "unresolved_vv", 0) + 1
            continue                                   # undecided: the pair keeps its enclosure (a hit if it meets)
        lo_t, hi_t = np.inf, -np.inf
        for (x, xe) in xs:
            bx = self._ball1200(x)
            for (y, ye) in ys:
                # a difference certainly >= 1e-299 is in range: no exact work needed ("out": not a new value;
                # "in": its enclosure comes from the two 1200-bit balls)
                try:
                    db = (bx - self._ball1200(y))
                    big = db.abs_lower() > TH
                except Exception:
                    big = False
                if big:
                    if want == "out":
                        continue
                    old = ctx.prec; ctx.prec = 256
                    lb = abs(db).log(); ctx.prec = old
                    lo_t = min(lo_t, np.nextafter(float(lb.lower()), -np.inf)); hi_t = max(hi_t, np.nextafter(float(lb.upper()), np.inf))
                    continue
                d = ex.sub(x, y)
                try:
                    sg = ex.sign(d)
                except Unresolved:
                    sg = None
                if sg == 0:
                    continue
                if sg is None:
                    lo_t = -np.inf; hi_t = max(hi_t, float(r.B[t])); continue
                b = ex.tight_ball(ex.absval(d))
                old = ctx.prec; ctx.prec = 256
                lb = b.log(); ctx.prec = old
                a_, b_ = np.nextafter(float(lb.lower()), -np.inf), np.nextafter(float(lb.upper()), np.inf)
                if keep_fn(np.array([a_]), np.array([b_]), np.array([False]))[0]:
                    lo_t = min(lo_t, a_); hi_t = max(hi_t, b_)
        self.resolved_vv = getattr(self, "resolved_vv", 0) + 1
        if lo_t > hi_t or hi_t < wlo[t] or lo_t > whi[t]:
            meet[t] = False
        else:
            r.A[t], r.B[t] = max(lo_t, float(r.A[t])), min(hi_t, float(r.B[t]))
    return meet


new_fn = qx.QX.resolve_vv
stat = {"calls": 0, "pairs": 0, "dropped": 0, "narrowed": 0, "bad": 0}


def both(self, S1, S2, I, sub, amb, meet, wlo, whi, r, want="out"):
    A0, B0 = r.A.copy(), r.B.copy()
    m_ref = resolve_vv_ref(self, S1, S2, I, sub, amb, meet, wlo, whi, r, want)
    A_ref, B_ref = r.A.copy(), r.B.copy()
    r.A[:] = A0; r.B[:] = B0
    m_new = new_fn(self, S1, S2, I, sub, amb, meet, wlo, whi, r, want)
    stat["calls"] += 1; stat["pairs"] += int(amb.sum()); stat["dropped"] += int((meet & ~m_new).sum())
    stat["narrowed"] += int(((r.A != A0) | (r.B != B0)).sum())
    if not (np.array_equal(m_ref, m_new) and np.array_equal(A_ref, r.A) and np.array_equal(B_ref, r.B)):
        stat["bad"] += 1
        print("DIFFERENT:", S1, S2, want, int((m_ref != m_new).sum()), int((A_ref != r.A).sum()), int((B_ref != r.B).sum()))
    return m_new


qx.QX.resolve_vv = both
qx.QX.unresolved_vv = 0
e1 = q.V((1,)); k1 = len(e1)
wlo, whi = window("pi", 3.23e-14 * (1 + 1e-9))
x7 = e1.q7(np.arange(k1))
lo = np.full(k1, wlo); hi = np.full(k1, whi)
wins = [(a[a <= b], b[a <= b]) for f in range(9) for (a, b) in qx.pre2(f, x7, lo, hi, None)]
L = np.concatenate([a for a, b in wins]); H = np.concatenate([b for a, b in wins])
sides = {s[0]: s for s in q.table_side(S)}
S1_, S2_, forms = sides[S1]
res = qx.Hits.cat(q.split_hits(S1_, S2_, forms, L, H, None))
print(f"t_vvmemo {S} = {S1_} | {S2_}, {len(L)} windows: {len(res)} hits; resolve_vv calls {stat['calls']}, ambiguous pairs "
      f"{stat['pairs']} (dropped {stat['dropped']}, narrowed {stat['narrowed']}), exact queries {getattr(q, 'nxV', 0)}, "
      f"undecided {q.unresolved_vv}, motor B unresolved {q.motorb_unresolved()}; calls that differ: {stat['bad']}")
sys.exit(1 if stat["bad"] or q.motorb_unresolved() else 0)
