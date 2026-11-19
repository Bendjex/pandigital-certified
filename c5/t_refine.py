"""refine_base1 with the forward enclosures of all its hits computed in blocks (10/11) against the code it replaces (one
call of ivl.fwd per hit, copied below as refine_base1_ref): on every call made by a real query, the same hits are kept
(meet) and the same enclosures are left (r.A, r.B), bit for bit.
The query: one split S1 | S2 of the recursion for NE(S), S a set of digits, for 17 generic in-range windows (those
of the forms of V(1) around pi) and 12 windows far out of range.
usage: t_refine.py n digits_of_S digits_of_S1 [k: every k-th window]      e.g. t_refine.py 9 123456789 1234 4"""
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


def refine_base1_ref(self, f, S1, S2, x, sub, wide, meet, wlo, whi, r, lwin=None, c2l=None):
    """the hits `wide` of form f whose base is an in-range value possibly equal to 1 (its ln enclosure contains 0)
    and whose result is wide: the base is replaced by its exact members (exact_near1) and the hit is kept only
    if some member gives a value meeting the window.  For l^s, l^-s (f = 7, 8; base l from the query side),
    lwin = (lower, upper) per hit: the preimage window of l of the sub-query that returned it (Lemma X3: every l
    with f(s, l) in the window has ln l in it), and c2l = (sign, lower, upper) per hit: the second-level
    constraint of that sub-query (pre_level2).  The members are then only required inside the window and the
    constraint (03/11): a member outside them gives no value in the window for any s of the class."""
    meet = meet.copy()
    if not hasattr(self, "_n1arr_ref"):
        self._n1arr_ref = {}
    for t in np.where(wide)[0]:
        c2t = None
        if f in (5, 6):
            S, A, B = S1, x.A[t], x.B[t]
        else:
            S, A, B = S2, sub.A[t], sub.B[t]
            if lwin is not None:
                A2, B2 = max(float(A), float(lwin[0][t])), min(float(B), float(lwin[1][t]))
                if A2 <= B2:                     # (always: the hit meets its window)
                    A, B = A2, B2
                    pass
            if c2l is not None and int(c2l[0][t]) != 0:
                c2t = (int(c2l[0][t]), float(c2l[1][t]), float(c2l[2][t]))
        key = (tuple(S), float(A), float(B), qx._c2key(c2t))
        if key not in self._n1arr_ref:
            mem = self.exact_near1(S, A, B, c2t)
            if mem:
                one = np.array([m[0] for m in mem]); la = np.array([m[1] for m in mem]); lb = np.array([m[2] for m in mem])
                self._n1arr_ref[key] = ivl.Q7(la, lb, np.zeros(len(mem), bool), one)
            else:
                self._n1arr_ref[key] = None
        q = self._n1arr_ref[key]
        pass
        if q is None:
            meet[t] = False
            continue
        k = len(q.A)
        idx = np.full(k, t)
        other = sub.q7().take(idx) if f in (5, 6) else x.take(idx)
        rr, okk = ivl.fwd(f, q, other) if f in (5, 6) else ivl.fwd(f, other, q)
        m = okk & np.where(rr.Z, wlo[t] == -np.inf, (rr.B >= wlo[t]) & (rr.A <= whi[t]))
        if not m.any():
            meet[t] = False
        else:
            r.A[t], r.B[t] = rr.A[m].min(), rr.B[m].max()
    return meet


new_fn = qx.QX.refine_base1
stat = {"calls": 0, "hits": 0, "dropped": 0, "narrowed": 0, "bad": 0}


def both(self, f, S1, S2, x, sub, wide, meet, wlo, whi, r, lwin=None, c2l=None):
    A0, B0 = r.A.copy(), r.B.copy()
    m_ref = refine_base1_ref(self, f, S1, S2, x, sub, wide, meet, wlo, whi, r, lwin, c2l)
    A_ref, B_ref = r.A.copy(), r.B.copy()
    r.A[:] = A0; r.B[:] = B0
    m_new = new_fn(self, f, S1, S2, x, sub, wide, meet, wlo, whi, r, lwin, c2l)
    stat["calls"] += 1; stat["hits"] += int(wide.sum()); stat["dropped"] += int((meet & ~m_new).sum())
    stat["narrowed"] += int(((r.A != A0) | (r.B != B0)).sum())
    if not (np.array_equal(m_ref, m_new) and np.array_equal(A_ref, r.A, equal_nan=True) and np.array_equal(B_ref, r.B, equal_nan=True)):
        stat["bad"] += 1
        print("DIFFERENT:", f, S1, S2, int((m_ref != m_new).sum()), int((A_ref != r.A).sum()), int((B_ref != r.B).sum()))
    return m_new


qx.QX.refine_base1 = both
e1 = q.V((1,)); k1 = len(e1)
wlo, whi = window("pi", 3.23e-14 * (1 + 1e-9))
x7 = e1.q7(np.arange(k1))
lo = np.full(k1, wlo); hi = np.full(k1, whi)
wins = [(a[a <= b], b[a <= b]) for f in range(9) for (a, b) in qx.pre2(f, x7, lo, hi, None)]
rng = np.random.default_rng(4)
c = rng.choice([-1, 1], 12) * rng.uniform(700, 60000, 12)
L = np.concatenate([a for a, b in wins] + [c - 1e-8]); H = np.concatenate([b for a, b in wins] + [c + 1e-8])
if len(sys.argv) > 4:                                # fewer windows: every k-th
    L = L[::int(sys.argv[4])]; H = H[::int(sys.argv[4])]
sides = {s[0]: s for s in q.table_side(S)}
S1_, S2_, forms = sides[S1]
res = qx.Hits.cat(q.split_hits(S1_, S2_, forms, L, H, None))
print(f"t_refine {S} = {S1_} | {S2_}, {len(L)} windows: {len(res)} hits; refine_base1 calls {stat['calls']}, hits refined "
      f"{stat['hits']} (dropped {stat['dropped']}, narrowed {stat['narrowed']}), motor B unresolved {q.motorb_unresolved()}; "
      f"calls that differ: {stat['bad']}")
sys.exit(1 if stat["bad"] or not stat["hits"] or q.motorb_unresolved() else 0)
