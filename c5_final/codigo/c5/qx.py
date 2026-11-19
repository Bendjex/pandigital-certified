"""
Motor X, window queries (C5, Theorem X).

Q(S, kind, lo, hi): enclosures (with the index of the window they meet) covering
   kind 'V': V(S) (in-range values; motor A classes for |S| <= tv, recursion above);
   kind 'N': NE(S) (materialised for |S| <= tn, recursion above, by Lemma X1);
   kind 'U': both;
restricted to values that may lie in the window [lo_i, hi_i] (ln).  Recursion over the splits
S = S1 u S2 with S1 the table side (|S1| <= 4): for every s of the S1 table and every form, pre()
gives the windows of l (Lemma X3), the S2 query answers them, and fwd() keeps the pairs whose value
meets the original window (Lemma X2).  Theorem X: by induction on |S|, every value of the requested
kind lying in a window is covered by a returned enclosure for that window.
"""
import sys, os, itertools, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import paths
import ivl, motorx
from ivl import fwd, pre, INF
from motorx import Ent

CHUNK = int(os.environ.get("QX_CHUNK", "3000000"))     # pairs per chunk (memory: ~1 GB per million pairs)


class _Prov:
    """debug provenance stored by batches (arrays), PROV[pid] gives the description tuple"""
    def __init__(self):
        self.bases, self.batches, self.n = [], [], 0
    def __len__(self):
        return self.n
    def add(self, meta, *arrays):
        base = self.n
        self.bases.append(base); self.batches.append((meta, [np.asarray(a).copy() for a in arrays]))
        self.n += len(arrays[0])
        return base
    def __getitem__(self, pid):
        import bisect
        b = bisect.bisect_right(self.bases, pid) - 1
        meta, arrays = self.batches[b]
        o = pid - self.bases[b]
        if meta[0] == "tab":
            return meta + (int(arrays[0][o]),)
        if meta[0] == "zero":
            return ("zero",)
        if meta[0] == "near1":
            return ("near1", meta[1], int(arrays[0][o]))
        _, S1, k1, f, S2, k2 = meta
        return ("pair", S1, k1, int(arrays[0][o]), f, S2, k2, int(arrays[1][o]))


PROV = _Prov()     # debug provenance: PROV[pid] = description tuple
DEBUG = False
CLIP_OUT = os.environ.get("QX_CLIP_OUT", "1") == "1"     # (V, V) pairs: only the out-of-range parts of the windows
ZERO_V = os.environ.get("QX_ZERO_V", "1") == "1"         # recursive V queries also return the exact zero
SWAP = os.environ.get("QX_SWAP", "1") == "1"             # combine_swapped when the query side is a much smaller table
SWAP_RATIO = float(os.environ.get("QX_SWAP_RATIO", "4"))
SW = (0, 1, 2, 4, 3, 7, 8, 5, 6)                          # f(s, l) = f'(l, s) with f' = SW[f] (forms of ivl.py)
INDEX_NEW = os.environ.get("QX_INDEX", "new") == "new"   # IndexT (tiers by absolute width) instead of Index
BATCH = os.environ.get("QX_BATCH", "0") == "1"
BATCHREC = os.environ.get("QX_BATCHREC", "1") == "1"   # (10/11) recursive sub-queries: every form of a chunk at once
DEDUP = os.environ.get("QX_DEDUP", "1") == "1"         # (10/11) recursive queries: each distinct window answered once
GROUP_VV = os.environ.get("QX_GROUP_VV", "1") == "1"      # resolve_vv: one exact query per class of s
C2PROP = os.environ.get("QX_C2PROP", "1") == "1"         # second-level window constraints passed down the recursion
EXACT1_C2 = os.environ.get("QX_EXACT1_C2", "1") == "1"   # (14/11) pre2: base exactly 1 against a window with a sign
RESOLVE_IN = os.environ.get("QX_RESOLVE_IN", "1") == "1"   # recursive V queries: exact |s - l| of overlapping classes
NEAR1 = os.environ.get("QX_NEAR1", "1") == "1"           # windows within 1e-300 of 1: answered by near1(S)
NEAR1_MAXSET = int(os.environ.get("QX_NEAR1_MAXSET", "5"))
DELTA1 = float(os.environ.get("QX_NEAR1_R", "1e-11"))  # radius (in ln) of the near-1 rows near1(S): |ln v| <= DELTA1
MEMBERS_MAX = int(os.environ.get("QX_MEMBERS_MAX", "500000"))  # (04/11) exact values kept by the cache of members()
XV_MAX = int(os.environ.get("QX_XV_MAX", "100000"))   # (13/11) exact values kept by the cache of exact_V
NARROW_REFINE = os.environ.get("QX_NARROW_REFINE", "1") == "1"   # refine_base1 of l^s, l^-s: exact members of l only
                                                                    # inside l's preimage window and second level


def _on_unique(fn, x):
    """fn(x) evaluated once per distinct value of x: the window bounds of a batch repeat once per pair (W, I),
    and ivl.ln_lo/ln_hi are pure functions of the double, so the result is identical to fn(x)"""
    x = np.asarray(x, dtype=np.float64)
    if x.size < 64:
        return fn(x)
    u, inv = np.unique(x, return_inverse=True)
    if len(u) * 2 > x.size:
        return fn(x)
    return fn(u)[inv.reshape(x.shape)]


def ln_lo_u(x):
    return _on_unique(ivl.ln_lo, x)


def ln_hi_u(x):
    return _on_unique(ivl.ln_hi, x)


PRE2 = os.environ.get("QX_PRE2", "1") == "1"            # pre() also with the second level of s (forms s^l, s^-l)


def pre2(f, x7, al, be, c2p=None):
    """preimage windows of Lemma X3 for the batch x7 (Q7 of s), intersected, for the forms s^l and s^-l, with the
    window given by the second level of s: if sign(ln s) = sg and ln|ln s| in [M, N], then +-l ln s in [al, be]
    means l |ln s| in [a, b] with (a, b) = (al, be), (-be, -al), (-be, -al), (al, be) for (s^l, sg = +1), (s^l, -1),
    (s^-l, +1), (s^-l, -1); so b > 0 and ln l in [ln a - N, ln b - M] (lower end -inf if a <= 0).  Both windows
    contain every l with f(s, l) in the window, so their intersection does too.  Not applied where the window
    contains the value 1 (l = 0 gives s^0 = 1) or s is 0 or exactly 1."""
    wins = pre(f, x7.A, x7.B, x7.Z, x7.O, al, be)
    if c2p is not None and f in (5, 6) and EXACT1_C2:
        # (14/11) a base exactly 1: s^l = s^-l = 1 for every l.  When the second level of the window gives the sign of
        # ln w (so w != 1), no l gives a value in the window: the preimage is empty.  Without this, a window that
        # contains ln 1 = 0 only through outward rounding (ln l in [-5e-324, 6e-309] for l^s = pi with s ~ e^800,
        # split 100 of n = 9) made the window of l the whole line, and the recursion enumerated every value of a
        # 5-digit V (12 million) only for c2_drop to remove the results at the end (an exact 1 against a constraint
        # with a sign).  The same pairs are removed, before they are enumerated.
        one = x7.O & ~x7.Z
        if one.any():
            sgw, _, _, _ = _wlog(al, be, c2p)
            kill = one & (sgw != 0)
            if kill.any():
                wins = [(np.where(kill, INF, lo_), np.where(kill, -INF, hi_)) for (lo_, hi_) in wins]
    if not PRE2 or f not in (5, 6):
        return wins
    S, M, N = x7.S, x7.M, x7.N
    base_ok = (S != 0) & ~x7.Z & ~x7.O & np.isfinite(M) & np.isfinite(N)
    use = base_ok & ~((al <= 0) & (be >= 0))             # first-level windows that do not contain ln 1 = 0
    if c2p is not None:
        # (03/11) with the second level of the window, the sign of w is known, so w != 0 and l = 0 (s^0 = 1) is not
        # admissible: the window of the second level applies even when [al, be] contains 0
        sgw, lw, hw, empty = _wlog(al, be, c2p)
        use_c2 = base_ok & (sgw != 0)
    else:
        use_c2 = np.zeros(len(al), bool)
    if not (use.any() or use_c2.any()):
        return wins
    flip = (S < 0) if f == 5 else (S > 0)
    a = np.where(flip, -be, al); b = np.where(flip, -al, be)
    with np.errstate(all="ignore"):
        lo2 = np.where(a > 0, ivl.addlo(ln_lo_u(np.where(a > 0, a, 1.0)), -N), -INF)
        hi2 = np.where(b > 0, ivl.addhi(ln_hi_u(np.where(b > 0, b, 1.0)), -M), -INF)
    lo2 = np.where(b > 0, lo2, INF)                       # b <= 0: no l at all
    lo2 = np.where(use, lo2, -INF); hi2 = np.where(use, hi2, INF)    # no first-level bound where the window holds 0
    if c2p is not None:
        # l |ln s| = +-w with ln|w| in [lw, hw]: ln l in [lw - N, hw - M]; intersected with the windows above (each
        # contains every admissible l)
        need = np.where(flip, -1, 1)
        with np.errstate(all="ignore"):
            l3 = ivl.addlo(lw, -N); h3 = ivl.addhi(hw, -M)
        good = use_c2 & (sgw == need) & ~empty
        lo2 = np.where(good, np.maximum(lo2, l3), lo2); hi2 = np.where(good, np.minimum(hi2, h3), hi2)
        bad_ = use_c2 & ((sgw != need) | empty)
        lo2 = np.where(bad_, INF, lo2); hi2 = np.where(bad_, -INF, hi2)
    app = use | use_c2
    out = []
    for (lo, hi) in wins:
        out.append((np.where(app, np.maximum(lo, lo2), lo), np.where(app, np.minimum(hi, hi2), hi)))
    return out


LEVEL2 = os.environ.get("QX_LEVEL2", "1") == "1"        # second-level constraint for the bases of l^s, l^-s
C2STAT = [0]


def _wlog(al, be, c2p):
    """sign and bounds of ln|w| for the values w = ln u of the windows [al, be], intersected with the second-level
    constraint c2p = (sign of w, lower, upper of ln|w|) of the windows when given (03/11).  Returns (sg, lo, hi, empty):
    sg = 0 where the sign of w is unknown; empty marks windows whose two descriptions contradict each other."""
    n = len(al)
    sg = np.zeros(n, np.int8); lo = np.full(n, -INF); hi = np.full(n, INF)
    pos = al > 0; neg = be < 0
    with np.errstate(all="ignore"):
        if pos.any():
            sg[pos] = 1; lo[pos] = ln_lo_u(al[pos]); hi[pos] = ln_hi_u(be[pos])
        if neg.any():
            sg[neg] = -1; lo[neg] = ln_lo_u(-be[neg]); hi[neg] = ln_hi_u(-al[neg])
    lo = np.where(np.isnan(lo), -INF, lo); hi = np.where(np.isnan(hi), INF, hi)
    empty = np.zeros(n, bool)
    if c2p is not None:
        cs, cm, cM = c2p
        use = cs != 0
        empty = use & (sg != 0) & (sg != cs)
        ok = use & ~empty
        sg = np.where(ok, cs, sg).astype(np.int8)
        lo = np.where(ok, np.maximum(lo, cm), lo); hi = np.where(ok, np.minimum(hi, cM), hi)
        empty |= ok & (lo > hi)
    return sg, lo, hi, empty


def pre_level2(f, x7, al, be, c2p=None):
    """for l^s (f = 7) and l^-s (f = 8): +-s ln l = w with w in the window [al, be] (and its second-level constraint
    c2p, when given) and s = e^lam, lam in [A, B] (the batch x7 of s).  If the sign of w is known, ln l has the sign
    of +-w and ln|ln l| = ln|w| - lam lies in [lo(ln|w|) - B, hi(ln|w|) - A].  Returns (sign, lower, upper) arrays
    (sign 0: no constraint)."""
    sgw, lw, hw, empty = _wlog(al, be, c2p)
    A, B = x7.A, x7.B
    n = len(al)
    sg = np.zeros(n, np.int8); m2 = np.full(n, -INF); M2 = np.full(n, INF)
    ok = ~x7.Z & ~x7.O & np.isfinite(A) & (sgw != 0)
    with np.errstate(all="ignore"):
        if ok.any():
            sg[ok] = sgw[ok] * (1 if f == 7 else -1)
            m2[ok] = ivl.addlo(lw[ok], -B[ok])
            M2[ok] = ivl.addhi(hw[ok], -A[ok])
    m2 = np.where(np.isnan(m2), -INF, m2); M2 = np.where(np.isnan(M2), INF, M2)
    # a contradictory window admits no l: an empty second-level interval
    m2 = np.where(empty & ok, INF, m2); M2 = np.where(empty & ok, -INF, M2)
    return sg, m2, M2


def level2_lin(f, x7, al, be, c2p):
    """(03/11) second-level constraint on l for the forms s + l, |s - l|, s l, s/l, l/s when the windows w = ln u are
    beyond the double range and carry a second level (c2p), with s of moderate size (x7: ln s in [A, B], finite):
      s l, s/l, l/s: ln l = w - lam, lam - w, w + lam with |lam| <= L = max(|A|, |B|) <= |w|/2, so ln l has the sign of
        w (-w for s/l) and |ln l| is in [|w| - L, |w| + L]: ln|ln l| in [lo + ln(1 - t), hi + t], t = L e^-lo;
      s + l = u with w > 0 and ln s <= w - 1: l in [u - s, u] gives ln l in [w - 1, w];  with w < 0: l <= u, so ln l <= w
        < 0 and ln|ln l| >= lo (no upper bound);
      |s - l| = u with w > 0 and ln s <= w - 1: l in [u - s, u + s] gives ln l in [w - 1, w + 1];
    and for w in [w - 1, w + 1] with w >= e^lo > 1: ln(w - 1) >= lo + ln(1 - e^-lo), ln(w + 1) <= hi + e^-lo.
    Returns (sign, lower, upper) arrays (sign 0: no constraint)."""
    n = len(al)
    sg = np.zeros(n, np.int8); m2 = np.full(n, -INF); M2 = np.full(n, INF)
    if c2p is None:
        return sg, m2, M2
    sgw, lw, hw, empty = _wlog(al, be, c2p)
    A, B = x7.A, x7.B
    with np.errstate(all="ignore"):
        fin = np.isfinite(A) & np.isfinite(B) & ~x7.Z & (sgw != 0) & ~empty & np.isfinite(lw) & (lw > 2.0)
        if f in (2, 3, 4):
            L = np.maximum(np.abs(A), np.abs(B))
            t = ivl.succ(L * ivl.exp_hi(-lw))
            ok = fin & (t <= 0.5)
            if ok.any():
                sg[ok] = sgw[ok] * (-1 if f == 3 else 1)
                m2[ok] = ivl.addlo(lw[ok], ivl.log1m_lo(t[ok]))
                M2[ok] = ivl.addhi(hw[ok], t[ok])
        elif f in (0, 1):
            e = ivl.exp_hi(-lw)                            # >= 1/w
            small_s = fin & (ivl.addhi(B, 1.0) <= ivl.exp_lo(lw))
            ok = small_s & (sgw > 0)
            if ok.any():
                sg[ok] = 1
                m2[ok] = ivl.addlo(lw[ok], ivl.log1m_lo(e[ok]))
                M2[ok] = hw[ok] if f == 0 else ivl.addhi(hw[ok], e[ok])
            if f == 0:
                okn = ~x7.Z & (sgw < 0) & ~empty & np.isfinite(lw)
                if okn.any():
                    sg[okn] = -1; m2[okn] = lw[okn]; M2[okn] = INF
    m2 = np.where(np.isnan(m2), -INF, m2); M2 = np.where(np.isnan(M2), INF, M2)
    return sg, m2, M2


ONEDIFF_C2 = os.environ.get("QX_ONEDIFF_C2", "1") == "1"   # (16/11) |1 - l| = u: the second level of l per branch


def one_diff_c2(bi, x7, al, be):
    """(16/11) the second-level constraint on l for the form |s - l| = u with s exactly 1 (flag O, not 0), for the
    branch bi of ivl.pre (0: l = s + u, 1: l = |s - u|), where the window [al, be] of ln u excludes u = 0 (al finite)
    and u <= e^-1 (be <= -1):
      branch 0: l = 1 + u, ln l = log1p(u) > 0, and u (1 - u/2) <= log1p(u) <= u gives
                ln|ln l| in [ln u + ln(1 - u/2), ln u], contained in [al - e^be, be]  (ln(1 - x) >= -2x for x <= 1/2);
      branch 1: l = 1 - u (u < 1), ln l = log1p(-u) < 0, and u <= -log1p(-u) <= u / (1 - u) gives
                ln|ln l| in [ln u, ln u - ln(1 - u)], contained in [al, be + 2 e^be]  (-ln(1 - u) <= 2u for u <= 1/2).
    Split 20 of n = 9: a window u ~ e^-240779 gave for l = 1 +- u the window ln l in [-2e-323, 1.5e-323] (outward
    rounding: it contains ln 1 = 0) without a second level, and one level down the forms s^l admitted every tiny l
    (about 3e5 values per window) for windows that the true l (|ln l| ~ e^-240779) cannot meet.  Returns
    (sign, lower, upper) arrays (sign 0: no constraint), or None where the rule applies nowhere."""
    one = x7.O & ~x7.Z
    ok = one & np.isfinite(al) & (be <= -1.0)
    if not ok.any():
        return None
    n = len(al)
    sg = np.zeros(n, np.int8); m2 = np.full(n, -INF); M2 = np.full(n, INF)
    with np.errstate(all="ignore"):
        eps = ivl.exp_hi(be[ok])                           # >= u for every u of the window
        if bi == 0:
            sg[ok] = 1
            m2[ok] = ivl.addlo(al[ok], -eps)
            M2[ok] = be[ok]
        else:
            sg[ok] = -1
            m2[ok] = al[ok]
            M2[ok] = ivl.addhi(be[ok], ivl.succ(2.0 * eps))
    return sg, m2, M2


def c2_drop(r, widx, c2):
    """entries (Q7-like r, window index widx) that cannot satisfy the second-level constraint c2 of their window:
    a known sign of ln u different from the window's, a second level disjoint from it, an exact 1 (ln u = 0 has no
    sign) or an exact 0 (ln u = -inf: sign -1 and ln|ln u| = +inf) against a constraint that excludes it"""
    sg, m2, M2 = c2
    s_w = sg[widx]; lo_w = m2[widx]; hi_w = M2[widx]
    use = s_w != 0
    drop = use & (r.S != 0) & ((r.S != s_w) | (r.N < lo_w) | (r.M > hi_w))
    drop |= use & r.O & ~r.Z
    drop |= use & r.Z & ((s_w > 0) | (hi_w < INF))
    return drop


def c2_and(c2a, c2b):
    """the intersection of two second-level constraints on the same values (either may be None)"""
    if c2a is None:
        return c2b
    if c2b is None:
        return c2a
    sa, ma, Ma = c2a; sb, mb, Mb = c2b
    sg = np.where(sa != 0, sa, sb).astype(np.int8)
    both = (sa != 0) & (sb != 0)
    contra = both & (sa != sb)
    m = np.where(both, np.maximum(ma, mb), np.where(sa != 0, ma, mb))
    M = np.where(both, np.minimum(Ma, Mb), np.where(sa != 0, Ma, Mb))
    m = np.where(contra, INF, m); M = np.where(contra, -INF, M)
    return sg, m, M


def out_parts(lo, hi):
    """the out-of-range parts of the windows [lo, hi] (in ln): the tiny part [lo, min(hi, LN_LO_HI)] and the
    huge part [max(lo, LN_HI_LO), hi], as arrays (lo', hi', idx), idx the index of the original window.
    Every out-of-range magnitude u (u < 1e-300 or u > 1e300) with ln u in [lo, hi] lies in one of the parts,
    since LN_LO_HI >= ln 1e-300 and LN_HI_LO <= ln 1e300 (ivl.py)."""
    it = np.where(lo < ivl.LN_LO_HI)[0]
    ig = np.where(hi > ivl.LN_HI_LO)[0]
    L = np.concatenate([lo[it], np.maximum(lo[ig], ivl.LN_HI_LO)])
    H = np.concatenate([np.minimum(hi[it], ivl.LN_LO_HI), hi[ig]])
    return L, H, np.concatenate([it, ig]).astype(np.int64)


def zero_hits(widx):
    """the exact zero, for the windows widx (V queries)"""
    k = len(widx)
    q = ivl.Q7(np.full(k, -INF), np.full(k, -INF), np.ones(k, bool), np.zeros(k, bool))
    h = Hits(np.asarray(widx, np.int64), q)
    if DEBUG:
        base = PROV.add(("zero",), np.zeros(k, np.int64))
        h.pid = np.arange(base, base + k)
    return h


class Hits:
    """matches: window index w, enclosure (ivl.Q7 fields), provenance id (DEBUG)"""
    F = ivl.Q7.FIELDS

    def __init__(self, w, q, pid=None, ti=None):
        self.w = w
        for k in self.F:
            setattr(self, k, getattr(q, k))
        self.pid = pid if pid is not None else np.full(len(w), -1, np.int64)
        self.ti = ti if ti is not None else np.full(len(w), -1, np.int64)     # table entry index (V lookups)
        self.isv = np.zeros(len(w), bool)                                    # value from a V query

    def __len__(self):
        return len(self.w)

    def q7(self):
        return ivl.Q7(*[getattr(self, k) for k in self.F])

    @staticmethod
    def cat(lst):
        lst = [h for h in lst if len(h)]
        h = Hits.__new__(Hits)
        if not lst:
            e = ivl.Q7(np.empty(0), np.empty(0), np.empty(0, bool), np.empty(0, bool))
            return Hits(np.empty(0, np.int64), e)
        for k in Hits.F + ("w", "pid", "ti", "isv"):
            setattr(h, k, np.concatenate([getattr(x, k) for x in lst]))
        return h

    def sel(self, m):
        h = Hits.__new__(Hits)
        for k in Hits.F + ("w", "pid", "ti", "isv"):
            setattr(h, k, getattr(self, k)[m])
        return h


class TooMany(Exception):
    """a table lookup would test more candidates than its limit: the caller takes the windows in smaller groups"""


def _spend(bud, k):
    """count k more candidates against the budget bud = [remaining] of a lookup (None: no limit)"""
    if bud is not None:
        bud[0] -= k
        if bud[0] < 0:
            raise TooMany()


def _ranges(i0, i1, bud=None):
    """(window index, position) for the half-open ranges [i0[w], i1[w]) (empty when i1 <= i0)"""
    cnt = np.maximum(i1 - i0, 0)
    tot = int(cnt.sum())
    _spend(bud, tot)
    if not tot:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    w = np.repeat(np.arange(len(i0)), cnt)
    off = np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    return w, np.repeat(i0, cnt) + off


def _mk_tiers(idx, key, width, T):
    """the entries idx split in tiers of width (width in (T[t-1], T[t]]), each sorted by key: [(idx, key, wmax)]"""
    out = []
    lo_t = -INF
    for t in T:
        m = (width > lo_t) & (width <= t) if np.isfinite(t) else (width > lo_t)
        if m.any():
            o = np.argsort(key[m], kind="stable")
            out.append((idx[m][o], key[m][o], max(0.0, float(width[m].max()))))
        lo_t = t
    return out


def _slack(x):
    """2^-40 |x| + 2^-1000 (inf for infinite x): covers the rounding errors of np.log near x with a wide margin"""
    with np.errstate(invalid="ignore", over="ignore"):
        return np.abs(x) * 2.0 ** -40 + 2.0 ** -1000


class IndexT:
    """enclosure table for window lookups, by classes of entries.  look() returns exactly the pairs of the plain test
        zero entry: lo = -inf;   other entries: B >= lo and A <= hi
    (the same as the original Index); the classes only decide which candidates are tested:
      - finite one-signed enclosures, A > 0 or B < 0, keyed by ln|.|: (ln A, ln B) or (ln -B, ln -A), in tiers of
        the width of that key (relative width of the enclosure), each sorted by its first key;
      - the other finite enclosures (A <= 0 <= B), in tiers of absolute width, each sorted by A;
      - half-lines [-inf, B] (sorted by B), half-lines [A, +inf] (sorted by A), the whole line, the exact zeros.
    Version ixT1 had only the absolute tiers for every finite enclosure: entries with |A| ~ 1e72 and width ~ 1e60
    (relative 1e-12) became candidates of every window within 1e60 of them, e.g. near ln 1 (n = 9: 4e8 candidates for
    one batch).  Candidates, for a tier of width wmax and a window [lo, hi]:
      key a in [L - wmax (1 + 2^-40) - slack(L), H + slack(H)], rounded outward twice, where for A > 0: L = ln lo
      (-inf if lo <= 0), H = ln hi (windows with hi > 0 only); for B < 0: L = ln -hi (-inf if hi >= 0), H = ln -lo
      (windows with lo < 0 only).
    A match (B >= lo, A <= hi, A > 0) has ln A <= ln hi and ln B >= ln lo, so fl(ln A) <= H + errors and
    fl(ln A) >= fl(ln B) - wmax (1 + 2^-52) >= L - wmax (1 + 2^-52) - errors; the errors of np.log (a few ulps) are
    far below slack(.) and wmax 2^-40 (same for B < 0 with the mirrored keys).  So every match is a candidate, and
    the final test on (A, B) itself gives exactly the plain pairs."""
    TIERS = (2.0 ** -40, 2.0 ** -20, 2.0 ** -6, 1.0, 2.0 ** 6, 2.0 ** 20, 2.0 ** 60, 2.0 ** 200, 2.0 ** 600, INF)
    LTIERS = (2.0 ** -46, 2.0 ** -40, 2.0 ** -34, 2.0 ** -28, 2.0 ** -20, 2.0 ** -12, 2.0 ** -6, 1.0, 2.0 ** 4, INF)

    def __init__(self, e):
        self.e = e
        A, B, Z = e.A, e.B, e.Z
        nz = ~Z
        fa, fb = np.isfinite(A), np.isfinite(B)
        with np.errstate(invalid="ignore", over="ignore"):
            width = np.where(fa & fb, B - A, INF)
        fin = nz & fa & fb & ~np.isnan(width)
        pos = fin & (A > 0)
        neg = fin & (B < 0)
        mid = fin & ~pos & ~neg
        im = np.where(mid)[0]
        self.tiers = _mk_tiers(im, A[im], width[im], self.TIERS)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            ip = np.where(pos)[0]; ka = np.log(A[ip]); kb = np.log(B[ip])
            self.ptiers = _mk_tiers(ip, ka, kb - ka, self.LTIERS)
            ineg = np.where(neg)[0]; ka = np.log(-B[ineg]); kb = np.log(-A[ineg])
            self.ntiers = _mk_tiers(ineg, ka, kb - ka, self.LTIERS)
        m = nz & (A == -INF) & ~(B == INF) & ~np.isnan(B)      # [-inf, B] (B finite or -inf)
        idx = np.where(m)[0]; o = np.argsort(B[idx], kind="stable"); self.lh = (idx[o], B[idx[o]])
        m = nz & fa & (B == INF)                                  # [A, +inf]
        idx = np.where(m)[0]; o = np.argsort(A[idx], kind="stable"); self.rh = (idx[o], A[idx[o]])
        self.full = np.where(nz & (A == -INF) & (B == INF))[0]  # the whole line
        self.zero = np.where(Z)[0]
        known = np.zeros(len(A), bool)
        for idx, _, _ in self.tiers + self.ptiers + self.ntiers:
            known[idx] = True
        for idx in (self.lh[0], self.rh[0], self.full, self.zero):
            known[idx] = True
        self.other = np.where(~known)[0]                          # anything else (NaN, A = +inf): checked directly

    VERSION = "ixT2"
    PARTS = ("t", "p", "n")

    def _parts(self):
        return {"t": self.tiers, "p": self.ptiers, "n": self.ntiers}

    def save(self, d):
        """store the index arrays next to the cached table (memory-mapped by the other processes)"""
        import ncache
        if os.path.exists(os.path.join(d, "IX_" + self.VERSION)):
            return
        arr = {"ix_lh_i": self.lh[0], "ix_lh_B": self.lh[1], "ix_rh_i": self.rh[0], "ix_rh_A": self.rh[1],
               "ix_full": self.full, "ix_zero": self.zero, "ix_other": self.other}
        counts = []
        for pn, tiers in self._parts().items():
            arr[f"ix2_{pn}_wmax"] = np.array([w for _, _, w in tiers], np.float64)
            counts.append(len(tiers))
            for t, (idx, keys, _) in enumerate(tiers):
                arr[f"ix2_{pn}{t}_i"] = idx; arr[f"ix2_{pn}{t}_k"] = keys
        try:
            for k, a in arr.items():
                tmp = os.path.join(d, f"{k}.tmp{os.getpid()}.npy")
                np.save(tmp, np.ascontiguousarray(a)); os.replace(tmp, os.path.join(d, k + ".npy"))
            with open(os.path.join(d, "IX_" + self.VERSION), "w") as f:
                f.write(" ".join(map(str, counts)) + "\n")
        except OSError:
            pass

    @classmethod
    def load(cls, e, d):
        m = os.path.join(d, "IX_" + cls.VERSION)
        if not os.path.exists(m):
            return None
        try:
            counts = [int(x) for x in open(m).read().split()]
            L = lambda k: np.load(os.path.join(d, k + ".npy"), mmap_mode="r")
            ix = cls.__new__(cls)
            ix.e = e
            parts = []
            for pn, T in zip(cls.PARTS, counts):
                wm = L(f"ix2_{pn}_wmax")
                parts.append([(L(f"ix2_{pn}{t}_i"), L(f"ix2_{pn}{t}_k"), float(wm[t])) for t in range(T)])
            ix.tiers, ix.ptiers, ix.ntiers = parts
            ix.lh = (L("ix_lh_i"), L("ix_lh_B")); ix.rh = (L("ix_rh_i"), L("ix_rh_A"))
            ix.full, ix.zero, ix.other = L("ix_full"), L("ix_zero"), L("ix_other")
            return ix
        except Exception:
            return None

    @staticmethod
    def _search(tiers, L, H, wsel, lo, hi, e, W, E, bud=None):
        """candidates of the keyed tiers for the windows wsel with key windows [L, H]; appends the exact matches"""
        for idx, keys, wmax in tiers:
            lo2 = np.nextafter(np.nextafter(L - (wmax * (1 + 2.0 ** -40)) - _slack(L), -INF), -INF)
            hi2 = np.nextafter(np.nextafter(H + _slack(H), INF), INF)
            lo2 = np.where(np.isnan(lo2), -INF, lo2); hi2 = np.where(np.isnan(hi2), INF, hi2)
            w, pos = _ranges(np.searchsorted(keys, lo2, side="left"), np.searchsorted(keys, hi2, side="right"), bud)
            if len(w):
                k = idx[pos]; ww = wsel[w]
                m = (e.B[k] >= lo[ww]) & (e.A[k] <= hi[ww])
                W.append(ww[m]); E.append(k[m])

    def look(self, lo, hi, limit=None):
        """pairs (window index, entry index) with entry meeting [lo, hi].  limit (10/11): raises TooMany, before
        building them, if more than limit candidates would be tested (the pairs are otherwise the same)"""
        e = self.e
        W, E = [], []
        nw = len(lo)
        if not nw:
            return np.empty(0, np.int64), np.empty(0, np.int64)
        bud = None if limit is None else [int(limit)]
        with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
            for idx, As, wmax in self.tiers:
                # B >= lo and B - A <= width imply A >= lo - width: the candidates have A in [lo', hi] with
                # lo' <= lo - wmax (1 + 2^-40) - 2^-1000, rounded down twice (fl(B - A) may be below B - A)
                lo2 = np.nextafter(np.nextafter(lo - (wmax * (1 + 2.0 ** -40) + 2.0 ** -1000), -INF), -INF)
                lo2 = np.where(np.isnan(lo2), -INF, lo2)
                w, pos = _ranges(np.searchsorted(As, lo2, side="left"), np.searchsorted(As, hi, side="right"), bud)
                if len(w):
                    k = idx[pos]
                    m = (e.B[k] >= lo[w]) & (e.A[k] <= hi[w])
                    W.append(w[m]); E.append(k[m])
            if self.ptiers:                                       # A > 0: windows with hi > 0
                ws = np.where((hi > 0) & ~np.isnan(lo))[0]
                if len(ws):
                    L = np.where(lo[ws] > 0, np.log(np.where(lo[ws] > 0, lo[ws], 1.0)), -INF)
                    H = np.log(hi[ws])
                    self._search(self.ptiers, L, H, ws, lo, hi, e, W, E, bud)
            if self.ntiers:                                       # B < 0: windows with lo < 0
                ws = np.where((lo < 0) & ~np.isnan(hi))[0]
                if len(ws):
                    L = np.where(hi[ws] < 0, np.log(np.where(hi[ws] < 0, -hi[ws], 1.0)), -INF)
                    H = np.log(-lo[ws])
                    self._search(self.ntiers, L, H, ws, lo, hi, e, W, E, bud)
            # (a NaN bound fails the plain test: such windows get no pair from the half-lines and the whole line)
            idx, Bs = self.lh                                     # [-inf, B]: B >= lo (A = -inf <= hi unless hi is NaN)
            if len(idx):
                i0 = np.where(np.isnan(hi), len(idx), np.searchsorted(Bs, lo, side="left"))
                w, pos = _ranges(i0, np.full(nw, len(idx)), bud)
                W.append(w); E.append(idx[pos])
            idx, As = self.rh                                     # [A, +inf]: A <= hi (B = inf >= lo unless lo is NaN)
            if len(idx):
                i1 = np.where(np.isnan(lo) | np.isnan(hi), 0, np.searchsorted(As, hi, side="right"))
                w, pos = _ranges(np.zeros(nw, np.int64), i1, bud)
                W.append(w); E.append(idx[pos])
            if len(self.full):
                wf = np.where(~np.isnan(lo) & ~np.isnan(hi))[0]
                _spend(bud, len(wf) * len(self.full))
                w = np.repeat(wf, len(self.full)); W.append(w); E.append(np.tile(self.full, len(wf)))
            if len(self.zero):
                wz = np.where(lo == -INF)[0]
                if len(wz):
                    _spend(bud, len(wz) * len(self.zero))
                    W.append(np.repeat(wz, len(self.zero))); E.append(np.tile(self.zero, len(wz)))
            for k0 in self.other:
                m = (e.B[k0] >= lo) & (e.A[k0] <= hi)
                W.append(np.where(m)[0]); E.append(np.full(int(m.sum()), k0))
        if not W:
            return np.empty(0, np.int64), np.empty(0, np.int64)
        return np.concatenate(W).astype(np.int64), np.concatenate(E).astype(np.int64)


def _cand_counts(self, lo, hi):
    """(10/11) for each window, the number of candidate entries look() tests (the lengths of its key ranges, plus the
    half-lines, the whole line, the zeros and the unclassified entries): an upper bound of its number of pairs.  It
    is used only to cut a batch of windows into groups of bounded size (look_batches, QX.window_groups): it has no
    effect on the pairs returned."""
    nw = len(lo)
    c = np.zeros(nw, np.int64)
    if not nw:
        return c
    with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
        for idx, As, wmax in self.tiers:
            lo2 = np.nextafter(np.nextafter(lo - (wmax * (1 + 2.0 ** -40) + 2.0 ** -1000), -INF), -INF)
            lo2 = np.where(np.isnan(lo2), -INF, lo2)
            c += np.maximum(np.searchsorted(As, hi, side="right") - np.searchsorted(As, lo2, side="left"), 0)
        for tiers, ws, Lf, Hf in ((self.ptiers, np.where((hi > 0) & ~np.isnan(lo))[0], lambda w: np.where(lo[w] > 0, np.log(np.where(lo[w] > 0, lo[w], 1.0)), -INF), lambda w: np.log(hi[w])),
                                  (self.ntiers, np.where((lo < 0) & ~np.isnan(hi))[0], lambda w: np.where(hi[w] < 0, np.log(np.where(hi[w] < 0, -hi[w], 1.0)), -INF), lambda w: np.log(-lo[w]))):
            if not len(tiers) or not len(ws):
                continue
            L, H = Lf(ws), Hf(ws)
            for idx, keys, wmax in tiers:
                lo2 = np.nextafter(np.nextafter(L - (wmax * (1 + 2.0 ** -40)) - _slack(L), -INF), -INF)
                hi2 = np.nextafter(np.nextafter(H + _slack(H), INF), INF)
                lo2 = np.where(np.isnan(lo2), -INF, lo2); hi2 = np.where(np.isnan(hi2), INF, hi2)
                c[ws] += np.maximum(np.searchsorted(keys, hi2, side="right") - np.searchsorted(keys, lo2, side="left"), 0)
        idx, Bs = self.lh
        if len(idx):
            c += len(idx) - np.where(np.isnan(hi), len(idx), np.searchsorted(Bs, lo, side="left"))
        idx, As = self.rh
        if len(idx):
            c += np.where(np.isnan(lo) | np.isnan(hi), 0, np.searchsorted(As, hi, side="right"))
        c += len(self.full) + len(self.other)
        if len(self.zero):
            c += np.where(lo == -INF, len(self.zero), 0)
    return c


def _groups(counts, limit):
    """consecutive groups [i, j) of the windows 0..n-1, each with at least one window and a sum of counts <= limit
    (a single window above the limit is a group of its own)"""
    n = len(counts)
    if not n:
        return []
    cs = np.cumsum(counts)
    if int(cs[-1]) <= limit:
        return [(0, n)]
    out = []
    i = 0
    while i < n:
        base = int(cs[i - 1]) if i else 0
        j = max(int(np.searchsorted(cs, base + limit, side="right")), i + 1)
        out.append((i, j)); i = j
    return out


def _look_batches(self, lo, hi, limit):
    """look(lo, hi) by consecutive groups of windows of at most limit candidates (cand_counts): yields the pairs
    (window index, entry index) of each group; together they are the pairs of look(lo, hi)"""
    for (i, j) in _groups(self.cand_counts(lo, hi), limit):
        w, k = self.look(lo[i:j], hi[i:j])
        yield w + i, k


IndexT.cand_counts = _cand_counts
IndexT.look_batches = _look_batches
LOOK_MAX = int(os.environ.get("QX_LOOK_MAX", str(4 * CHUNK)))   # candidates of a table lookup handled at once
REC_PAIRS = int(os.environ.get("QX_REC_PAIRS", str(2 * CHUNK)))  # (14/11) (window, sub-hit) pairs of a recursion at once
REC_WIN0 = int(os.environ.get("QX_REC_WIN0", "4"))     # (14/11) distinct windows of the first batch of a recursion (some windows have ~3e5 sub-hits)


class IndexN1(IndexT):
    """IndexT for the near-1 rows (|ln v| <= DELTA1): the same look(), with width tiers fine enough for enclosures
    of ln v of size ~1e-11 and below (the tiers only decide which candidates are tested, never the answer)"""
    TIERS = tuple(2.0 ** -k for k in (1000, 600, 400, 300, 250, 200, 170, 150, 130, 110, 100, 90, 80, 70, 60, 50, 45, 40, 36)) + (INF,)


class Index:
    """sorted enclosure table for window lookups (the original version; QX_INDEX=old)"""

    def __init__(self, e):
        self.e = e
        fin = np.isfinite(e.A) & np.isfinite(e.B)
        width = np.where(fin, e.B - e.A, INF)
        narrow = fin & (width <= 1e-3 * (1 + np.abs(e.A)))
        self.ni = np.where(narrow)[0]
        order = np.argsort(e.A[self.ni], kind="stable")
        self.ni = self.ni[order]
        self.nA = e.A[self.ni]
        self.wmax = float(width[narrow].max()) if narrow.any() else 0.0
        self.wi = np.where(~narrow)[0]               # wide, infinite or zero entries: checked directly

    def look(self, lo, hi):
        """pairs (window index, entry index) with entry meeting [lo, hi]"""
        e = self.e
        W, E = [], []
        if len(self.ni) and len(lo):
            # entries with A in [lo - wmax*(1+|lo|) , hi]
            slack = self.wmax * (1 + np.abs(np.where(np.isfinite(lo), lo, 0)))
            i0 = np.searchsorted(self.nA, np.where(np.isfinite(lo), lo - slack, -INF), side="left")
            i1 = np.searchsorted(self.nA, hi, side="right")
            cnt = np.maximum(i1 - i0, 0)
            tot = int(cnt.sum())
            if tot:
                w = np.repeat(np.arange(len(lo)), cnt)
                off = np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt)
                k = self.ni[np.repeat(i0, cnt) + off]
                m = (e.B[k] >= lo[w]) & (e.A[k] <= hi[w])
                W.append(w[m]); E.append(k[m])
        if len(self.wi) and len(lo):
            for j0 in range(0, len(self.wi), 2000):
                wj = self.wi[j0:j0 + 2000]
                ww = np.repeat(np.arange(len(lo)), len(wj)); kk = np.tile(wj, len(lo))
                ezero = e.Z[kk]
                m = np.where(ezero, lo[ww] == -INF, (e.B[kk] >= lo[ww]) & (e.A[kk] <= hi[ww]))
                W.append(ww[m]); E.append(kk[m])
        if not W:
            return np.empty(0, np.int64), np.empty(0, np.int64)
        return np.concatenate(W), np.concatenate(E)


def _look_all(self, lo, hi, limit):
    """the original Index (QX_INDEX=old): one batch"""
    yield self.look(lo, hi)


Index.look_batches = _look_all


class QX(motorx.MX):
    def __init__(self, n, ordered, res, tv, tn, verbose=True):
        super().__init__(n, ordered, res, tv)
        self.tn = tn
        self.cstats = {}
        self.idx = {}
        self.stats = {}
        self.verbose = verbose

    def index(self, S, kind):
        key = (S, kind)
        if key not in self.idx:
            if kind == "V":
                e = self.V(S) if len(S) <= self.tv else self.literals(S)
            else:
                e = self.N(S)
            d = getattr(e, "cache_dir", None)
            ix = IndexT.load(e, d) if (INDEX_NEW and d) else None
            if ix is None:
                ix = (IndexT if INDEX_NEW else Index)(e)
                if INDEX_NEW and d:
                    ix.save(d)
            self.idx[key] = ix
        return self.idx[key]

    def literals(self, S):
        """enclosures of the literals of S (for V-queries above the table size)"""
        perms = [S] if self.ord else itertools.permutations(S)
        vals = set()
        for p in perms:
            t = "".join(map(str, p))
            for i in range(len(t) + 1):
                vals.add((int(t), 10 ** (len(t) - i)))
        from math import gcd
        vals = sorted(vals)
        num = np.array([float(a) for a, b in vals]); den = np.array([float(b) for a, b in vals])
        x = num / den
        lo = ivl.ln_lo(np.nextafter(x, 0)); hi = ivl.ln_hi(np.nextafter(x, INF))
        P = np.array([a // gcd(a, b) for a, b in vals], np.int64); Q = np.array([b // gcd(a, b) for a, b in vals], np.int64)
        return Ent(lo, hi, np.zeros(len(x), bool), P == Q, None, P, Q)

    def tables_ok(self, S, kind):
        return (kind == "V" and len(S) <= self.tv) or (kind == "N" and len(S) <= self.tn)

    def lookup(self, S, kind, lo, hi, c2=None):
        ix = self.index(S, kind)
        lim = getattr(self, "_look_limit", None)        # (10/11) set by combine() around a table sub-query
        w, k = ix.look(lo, hi, lim) if (lim is not None and isinstance(ix, IndexT)) else ix.look(lo, hi)
        e = ix.e
        if c2 is not None and len(k):
            # second-level constraint of the windows (pre_level2): sign(ln l) = sg and ln|ln l| in [m, M];
            # an entry with a known second level outside it holds no admissible l (entries with S = 0 are kept)
            class _V:                                   # the entries' fields used by c2_drop
                pass
            v = _V(); v.S = e.S[k]; v.N = e.N[k]; v.M = e.M[k]; v.O = e.O[k]; v.Z = e.Z[k]
            drop = c2_drop(v, w, c2)
            if drop.any():
                keep = ~drop
                w, k = w[keep], k[keep]
                C2STAT[0] += int(drop.sum())
        pid = None
        if DEBUG and len(k):
            base = PROV.add(("tab", S, kind), k)
            pid = np.arange(base, base + len(k))
        return Hits(w, e.q7(k), pid, k if (kind == "V" and len(S) <= self.tv) else None)

    def table_side(self, S):
        """splits (S1 table side, S2 query side, forms)"""
        out = []
        if self.ord:
            for i in range(1, len(S)):
                L, R = S[:i], S[i:]
                if len(L) <= len(R):
                    out.append((L, R, ivl.ORD_LEFT))
                else:
                    out.append((R, L, ivl.ORD_RIGHT))
        else:
            for r in range(1, len(S) // 2 + 1):
                for S1 in itertools.combinations(S, r):
                    S2 = tuple(d for d in S if d not in S1)
                    if 2 * r == len(S) and S1 > S2:
                        continue
                    out.append((S1, S2, tuple(range(9))))
        return out

    def Q(self, S, kind, lo, hi, c2=None):
        """c2: optional second-level constraint of the windows (used by table lookups; a recursive query may
        ignore it, which only returns more)"""
        S = tuple(S)
        lo = np.asarray(lo, dtype=np.float64); hi = np.asarray(hi, dtype=np.float64)
        if not len(lo):
            return Hits.cat([])
        if kind == "U":
            return Hits.cat([self.Q(S, "V", lo, hi, c2), self.Q(S, "N", lo, hi, c2)])
        if self.tables_ok(S, kind):
            h = self.lookup(S, kind, lo, hi, c2)
            h.isv[:] = (kind == "V")
            return h
        if DEDUP and len(lo) > 1:
            # (10/11) windows equal bit for bit (bounds and second-level constraint) are the same query: it is
            # answered once and its hits are given to each of them.  The preimage windows of different entries and
            # forms coincide often (s l = w and l / (1/s) = w, the forms with s = 1, ...)
            rows = [lo, hi] if c2 is None else [lo, hi, c2[0].astype(np.float64), np.asarray(c2[1], np.float64),
                                                np.asarray(c2[2], np.float64)]
            rows = np.ascontiguousarray(np.stack(rows, 1))
            _, first, inv = np.unique(rows.view(np.dtype((np.void, rows.dtype.itemsize * rows.shape[1]))).ravel(),
                                      return_index=True, return_inverse=True)
            if len(first) < len(lo):
                h = self.Q(S, kind, lo[first], hi[first], tuple(a[first] for a in c2) if c2 is not None else None)
                inv = inv.ravel()
                order = np.argsort(h.w, kind="stable")
                cnt = np.bincount(h.w, minlength=len(first))
                start = np.cumsum(cnt) - cnt
                w, pos = _ranges(start[inv], start[inv] + cnt[inv])
                h = h.sel(order[pos])
                h.w = w
                st = self.stats.setdefault((len(S), kind), [0, 0, 0.0])
                st[0] += len(lo) - len(first); st[1] += len(h) - int(cnt.sum())    # the windows and hits not recomputed
                return h
        t0 = time.time()
        res = []
        if kind == "V" and NEAR1 and len(S) <= NEAR1_MAXSET:
            # windows inside |ln u| <= DELTA1 (values within DELTA1 of 1 in ln) are answered by near1(S), the
            # exact members of V(S) there (computed once per set, motor B): no enumeration of the recursion
            inside = (lo >= -DELTA1) & (hi <= DELTA1)
            if inside.any() and self.near1(S) is None:      # motor B could not certify the rows: recursion
                inside[:] = False
            if inside.any():
                hn = self.near1_hits(S, lo, hi, np.where(inside)[0])
                if inside.all():
                    hn.isv[:] = True
                    return hn
                res.append(hn)
                lo = np.where(inside, INF, lo); hi = np.where(inside, -INF, hi)
        if kind == "V":
            # the exact zero (a value of V(S) for many S): admitted by the windows with lo = -inf.  It is
            # returned as a possible member (superset), since the recursion below only meets nonzero values
            if ZERO_V:
                zw = np.where(lo == -INF)[0]
                if len(zw):
                    res.append(zero_hits(zw))
            # in range only: clip the windows to the range
            lo = np.maximum(lo, ivl.LN_LO_LO); hi = np.minimum(hi, ivl.LN_HI_HI)
            ok = lo <= hi
            res.append(self.lookup(S, "V", np.where(ok, lo, INF), np.where(ok, hi, -INF)))   # literals
            # V(S): s in V(S1), l in V(S2)
            for S1, S2, forms in self.table_side(S):
                res.append(self.combine(S1, "V", S2, "V", forms, lo, hi, want="in", c2=c2 if C2PROP else None))
        else:
            sides = self.table_side(S)
            # (10/11) c5run takes a query in pieces: self.only = {set: pred}, and the N query of that set keeps the
            # part p ("NU", "VN", "VV") of its split number si only if pred(si, p).  Q(S, "N") is the concatenation
            # of the hits of the parts of its splits, so pieces whose predicates partition them add up to the query
            pred = (getattr(self, "only", None) or {}).get(S)
            for si_, (S1, S2, forms) in enumerate(sides):
                parts = None
                if pred is not None:
                    parts = tuple(p_ for p_ in ("NU", "VN", "VV") if pred(si_, p_))
                    if not parts:
                        continue
                if S == getattr(self, "top", None):             # progress of the top-level query
                    el = time.time() - t0
                    print(f"[progress] split {si_}/{len(sides)} ({100 * si_ / len(sides):.1f}%) "
                          f"{''.join(map(str, S1))}|{''.join(map(str, S2))}, {el:.0f}s elapsed", flush=True)
                res += self.split_hits(S1, S2, forms, lo, hi, c2 if C2PROP else None, parts)
        h = Hits.cat(res)
        h.isv[:] = (kind == "V")
        if c2 is not None and len(h):              # the second-level constraint on the recursive results too
            drop = c2_drop(h, h.w, c2)
            if drop.any():
                C2STAT[0] += int(drop.sum())
                h = h.sel(~drop)
        k = (len(S), kind)
        st = self.stats.setdefault(k, [0, 0, 0.0])
        st[0] += len(lo); st[1] += len(h); st[2] += time.time() - t0
        return h

    def near1(self, S):
        """the exact members of V(S) with |ln v| <= DELTA1, as an Ent (exact 1: flag O; others: enclosures from
        Arb balls of ln v at the precision that separates it from 0).  Motor B's query with the window
        [1 - 2 DELTA1, 1 + 2 DELTA1] contains every such value (Lemma 14); a value whose |ln v| is shown to exceed
        DELTA1 is dropped.  Cached in memory and on disk (QX_NCACHE/near1, keyed by the code of motor B)."""
        S = tuple(S)
        if not hasattr(self, "_near1"):
            self._near1 = {}
        if S in self._near1:
            return self._near1[S]
        import json, hashlib
        from exact import Unresolved
        fn = None
        if motorx.NCACHE:
            h = hashlib.sha256()
            for nm in ("motorb.py", "exact.py", "canon.py", "exprparse.py"):
                h.update(open(os.path.join(paths.SRC, nm), "rb").read().replace(b"\r\n", b"\n"))
            h.update(f"near1 v2 ord={int(self.ord)} S={S} delta={DELTA1!r}".encode())
            d = os.path.join(motorx.NCACHE, "near1"); os.makedirs(d, exist_ok=True)
            fn = os.path.join(d, f"V{int(self.ord)}_{''.join(map(str, S))}_{h.hexdigest()[:20]}.json")
            if os.path.exists(fn):
                try:
                    rows = json.load(open(fn))
                    self.__dict__.setdefault("_near1_rowcache", {})[S] = rows
                    self._near1[S] = self._near1_ent(rows)
                    return self._near1[S]
                except ValueError:
                    pass
        import exact as ex
        from flint import arb, ctx
        t0 = time.time()
        # |ln v| <= DELTA1 < 1 implies |v - 1| <= e^DELTA1 - 1 <= 2 DELTA1: the window holds every such value
        dl = arb(DELTA1)
        try:
            vals = self.bquery(S, 1 - 2 * dl, 1 + 2 * dl)
        except Unresolved as u:
            self.say(f"near1 V{''.join(map(str, S))}: not available ({u}); the recursion is used instead")
            self._near1[S] = None
            return None
        rows = []
        for (v, e) in vals.values():
            try:
                sg = ex.sign(ex.sub(v, ex.ONE))
            except Unresolved:
                sg = None
            if sg == 0:
                rows.append([0.0, 0.0, 1, 0, 0.0, 0.0, e]); continue
            lb = None
            for prec in (256, 1024, 4096, 16384, 65536):
                old = ctx.prec; ctx.prec = prec
                try:
                    lb = v.ball(prec).log()
                finally:
                    ctx.prec = old
                if not lb.contains(0):
                    break
            if lb is None or lb.contains(0):              # undecided: kept as possibly 1, without a second level
                rows.append([-DELTA1, DELTA1, 0, 0, -INF, INF, e]); continue
            lo_, hi_ = float(lb.lower()), float(lb.upper())
            lo_, hi_ = float(np.nextafter(lo_, -INF)), float(np.nextafter(hi_, INF))
            if lo_ > DELTA1 or hi_ < -DELTA1:              # certainly not within DELTA1 of ln 1 (ln 1 = 0)
                continue
            # the sign of ln v from the ball (it excludes 0); the rounded bound can cross 0 when |ln v| < 2^-1075
            # (float() of a lower bound below half the least subnormal is 0, and nextafter gives -2^-1074):
            # then the bound on that side is 0, still valid since ln v has the sign of the ball
            if lb > 0:
                sgn = 1; lo_ = max(lo_, 0.0)
            else:
                sgn = -1; hi_ = min(hi_, 0.0)
            ab = abs(lb).log()
            rows.append([lo_, hi_, 0, sgn, float(np.nextafter(float(ab.lower()), -INF)), float(np.nextafter(float(ab.upper()), INF)), e])
        if fn:
            tmp = fn + f".tmp{os.getpid()}"
            json.dump(rows, open(tmp, "w")); os.replace(tmp, fn)
        self.say(f"near1 V{''.join(map(str, S))}: {len(vals)} candidates of motor B, {len(rows)} with |ln v| <= {DELTA1:g}, {time.time() - t0:.1f}s")
        self.__dict__.setdefault("_near1_rowcache", {})[S] = rows
        self._near1[S] = self._near1_ent(rows)
        return self._near1[S]

    def _near1_rows(self, S):
        """the cached rows of near1(S) (hitcheck, _xvals); each row ends with the writing of its value.  (13/11) When
        the rows are not in memory (c5run.release empties the caches between pieces) near1(S) is made again: its rows
        are a deterministic function of S (read back from QX_NCACHE, or recomputed), so the row numbers stay valid"""
        S = tuple(S)
        if S not in self.__dict__.get("_near1_rowcache", {}):
            self.__dict__.get("_near1", {}).pop(S, None)
        self.near1(S)
        return self._near1_rowcache[S]

    def near1_value(self, S, k):
        """(exact value, writing) of the row k of near1(S): the writing is motor B's expression of the value,
        re-evaluated exactly (exprparse.to_x, motor B's arithmetic), as hitcheck does"""
        import exprparse
        e = self._near1_rows(S)[k][6]
        return exprparse.to_x(e), e

    def bquery(self, S, lo, hi):
        """motor B query (Lemma 14: a superset of V(S) in [lo, hi]), refused when motor B has skipped an operation it
        could not decide (motorb.B.unresolved not empty): then its answers need not be supersets"""
        from exact import Unresolved
        B = self.B()
        r = B.query(tuple(str(d) for d in S), lo, hi)
        if B.unresolved:
            raise Unresolved(f"motor B left {len(B.unresolved)} operations undecided (first: {B.unresolved[0]!r:.200})")
        return r

    def motorb_unresolved(self):
        """number of operations motor B has left undecided in this process (must be 0)"""
        return len(self._B.unresolved) if hasattr(self, "_B") else 0

    def members(self, S, ent, i):
        """the exact values the entry i of the V table of S may stand for, with their expressions
        (motorx.MX.class_members), with a BOUNDED cache (04/11): least recently used entries are dropped beyond
        MEMBERS_MAX and recomputed if needed again.  The list is a deterministic function of (S, i), so a
        recomputation gives the same values; only memory and time change.  Refused if motor B has left an operation
        undecided: its query is then not certainly a superset."""
        from exact import Unresolved
        cache = self.__dict__.get("_memlru")
        if cache is None:
            import collections
            cache = self._memlru = collections.OrderedDict()
            self._memlru_n = 0
        key = (S, i)
        if key in cache:
            cache.move_to_end(key)
            return cache[key]
        r = self.class_members(S, ent, i)
        if hasattr(self, "_B") and self._B.unresolved:
            raise Unresolved(f"motor B left {len(self._B.unresolved)} operations undecided")
        cache[key] = r
        self._memlru_n += len(r) + 1
        while self._memlru_n > MEMBERS_MAX and len(cache) > 1:
            _, old = cache.popitem(last=False)
            self._memlru_n -= len(old) + 1
            self.members_evicted = getattr(self, "members_evicted", 0) + 1
        return r

    @staticmethod
    def _near1_ent(rows):
        n = len(rows)
        A = np.array([r[0] for r in rows], np.float64); B = np.array([r[1] for r in rows], np.float64)
        O = np.array([bool(r[2]) for r in rows], bool); Sg = np.array([r[3] for r in rows], np.int8)
        M = np.array([r[4] for r in rows], np.float64); N = np.array([r[5] for r in rows], np.float64)
        e = Ent(A, B, np.zeros(n, bool), O, None, None, None, S=Sg, M=M, N=N) if n else Ent.cat([])
        return e

    def near1_hits(self, S, lo, hi, widx):
        """Hits of near1(S) for the windows widx: the rows meeting each window (the plain test B >= lo, A <= hi,
        through IndexN1, which returns exactly those pairs)"""
        e = self.near1(S)
        if not len(e) or not len(widx):
            return Hits.cat([])
        key = ("near1ix", tuple(S))
        if key not in self.idx:
            self.idx[key] = IndexN1(e)
        w, K = self.idx[key].look(lo[widx], hi[widx])
        W = widx[w]
        h = Hits(W, e.q7(K), None)
        if DEBUG and len(K):
            base = PROV.add(("near1", S), K)
            h.pid = np.arange(base, base + len(K))
        return h

    def split_hits(self, S1, S2, forms, lo, hi, c2=None, parts=None):
        """the NE(S) values of the windows produced by one split S = S1 u S2 (S1 the table side), by
        Lemma X1: (NE(S1), U(S2)) [part "NU"], (V(S1), NE(S2)) ["VN"] and (V(S1), V(S2)) with an out-of-range
        result ["VV"]; parts (10/11): only those parts (None: the three)"""
        res = []
        if (parts is None or "NU" in parts) and (len(S1) >= 3 or self.tables_ok(S1, "N")):
            res.append(self.combine(S1, "N", S2, "U", forms, lo, hi, want="any", c2=c2))
        if parts is None or "VN" in parts:
            res.append(self.combine(S1, "V", S2, "N", forms, lo, hi, want="any", c2=c2))
        if parts is not None and "VV" not in parts:
            return res
        if CLIP_OUT:
            # (V, V) pairs give a new value only if it is out of range: only the out-of-range parts of
            # the windows are queried (out_parts), each hit keeps the index of its original window
            L, H, idx = out_parts(lo, hi)
            if len(L):
                c2o = tuple(a[idx] for a in c2) if c2 is not None else None
                h = self.combine(S1, "V", S2, "V", forms, L, H, want="out", c2=c2o)
                if len(h):
                    h.w = idx[h.w]
                res.append(h)
        else:
            out_possible = (lo < ivl.LN_LO_HI) | (hi > ivl.LN_HI_LO)
            if out_possible.any():
                li = np.where(out_possible, lo, INF); hi_ = np.where(out_possible, hi, -INF)
                res.append(self.combine(S1, "V", S2, "V", forms, li, hi_, want="out"))
        return res

    def window_groups(self, S, kind, lo, hi, mi):
        """consecutive groups (g0, g1) of the windows mi (indices into lo, hi) for the query (S, kind): when the
        query is answered by table lookups, groups of at most LOOK_MAX candidate entries (IndexT.cand_counts, which
        bounds the number of pairs of each lookup); a single group otherwise.  The answer to a window does not
        depend on the other windows of its batch, so the groups only bound the memory in use."""
        kinds = ("V", "N") if kind == "U" else (kind,)
        if len(mi) < 2 or not all(self.tables_ok(S, k) for k in kinds):
            return [(0, len(mi))]
        ixs = [self.index(S, k) for k in kinds]
        if not all(hasattr(ix, "cand_counts") for ix in ixs):
            return [(0, len(mi))]
        l, h = lo[mi], hi[mi]
        c = ixs[0].cand_counts(l, h)
        for ix in ixs[1:]:
            c = c + ix.cand_counts(l, h)
        return _groups(c, LOOK_MAX)

    def qtable(self, S, kind):
        """the table answering the query (S, kind) when S is within the table sizes: (Ent, isv, ti, n_V),
        exactly the entries a lookup of (S, kind) can return (kind U: the V table, then the N table)"""
        key = ("qtable", S, kind)
        if key not in self.idx:
            parts = []
            if kind in ("V", "U"):
                parts.append(self.V(S))
            if kind in ("N", "U"):
                parts.append(self.N(S))
            e = Ent.cat(parts) if len(parts) > 1 else parts[0]
            nv = len(parts[0]) if kind in ("V", "U") else 0
            isv = np.zeros(len(e), bool); isv[:nv] = True
            ti = np.full(len(e), -1, np.int64); ti[:nv] = np.arange(nv)
            self.idx[key] = (e, isv, ti, nv)
        return self.idx[key]

    def combine(self, S1, k1, S2, k2, forms, lo, hi, want, c2=None):
        """values f(s, l), s in table k1 of S1, l in query k2 of S2, meeting the windows"""
        if k1 == "N":
            e1 = self.N(S1)
        else:
            e1 = self.V(S1)
        n1, nw = len(e1), len(lo)
        if not n1 or not nw:
            return Hits.cat([])
        if k2 == "U" and self.tables_ok(S2, "V") != self.tables_ok(S2, "N"):
            # (10/11) U(S2) with one kind in a table and the other a recursion (|S2| = 5): the two kinds one after the
            # other, so that the table lookup is bounded (LOOK_MAX) and the recursion gets its windows together.  The
            # query U is the query V followed by the query N, and each sub-hit is treated on its own: the same hits
            return Hits.cat([self.combine(S1, k1, S2, "V", forms, lo, hi, want, c2),
                             self.combine(S1, k1, S2, "N", forms, lo, hi, want, c2)])
        if SWAP and all(self.tables_ok(S2, k) for k in (("V", "N") if k2 == "U" else (k2,))):
            e2, isv2, ti2, nv2 = self.qtable(S2, k2)
            if len(e2) * SWAP_RATIO < n1:
                return self.combine_swapped(S1, k1, e1, S2, k2, e2, isv2, ti2, nv2, forms, lo, hi, want)
        valid = lo <= hi
        widx = np.where(valid)[0]
        if not len(widx):
            return Hits.cat([])
        out = []
        if BATCH:
            return self.combine_batched(S1, k1, e1, S2, k2, forms, lo, hi, want, widx)
        # (10/11) a recursive sub-query (S2 above the table sizes) receives the windows of every form and branch of
        # the chunk at once, in at most two queries (the forms with a second-level constraint on l, and the others),
        # instead of one query per form and branch: the answer to a window does not depend on the other windows of
        # its batch, so the sub-hits are the same.  Table lookups keep one query per form and branch.
        tab2 = all(self.tables_ok(S2, k) for k in (("V", "N") if k2 == "U" else (k2,)))
        rec = BATCHREC and not tab2
        step = max(1, CHUNK // (n1 * (len(forms) if rec else 1)))
        for c0 in range(0, len(widx), step):
            wv = widx[c0:c0 + step]
            W = np.repeat(wv, n1); I = np.tile(np.arange(n1), len(wv))
            a1, b1, z1, o1 = e1.A[I], e1.B[I], e1.Z[I], e1.O[I]
            x7 = e1.q7(I)
            c2W = tuple(a[W] for a in c2) if c2 is not None else None   # the windows' second level, per pair
            pair = (S1, k1, S2, k2, e1, W, I, x7, a1, b1, o1, lo, hi, want, c2)
            grp = ([], [])                                 # rec: the windows without / with a constraint on l
            for f in forms:
                if LEVEL2 and f in (7, 8):
                    c2all = pre_level2(f, x7, lo[W], hi[W], c2W)
                elif c2W is not None and f <= 4:
                    c2all = level2_lin(f, x7, lo[W], hi[W], c2W)
                else:
                    c2all = None
                c2f = c2all
                for bi, (plo, phi) in enumerate(pre2(f, x7, lo[W], hi[W], c2W)):
                    c2all = c2f
                    if ONEDIFF_C2 and f == 1:
                        # (16/11) |1 - l| = u: each branch of l with its second level (one_diff_c2)
                        c2b = one_diff_c2(bi, x7, lo[W], hi[W])
                        if c2b is not None:
                            c2all = c2_and(c2f, c2b)
                    mi_all = np.where(plo <= phi)[0]
                    if not len(mi_all):
                        continue
                    ck = (len(S1), len(S2), k1, k2, f)
                    cs = self.cstats.setdefault(ck, [0, 0, 0])
                    cs[0] += len(mi_all)
                    if rec:
                        c2m = tuple(a[mi_all] for a in c2all) if c2all is not None else None
                        grp[c2all is not None].append((f, mi_all, plo[mi_all], phi[mi_all], c2m))
                        continue
                    # (10/11) a sub-query answered by tables: the lookup counts its candidates and stops (TooMany) above
                    # LOOK_MAX; the windows then go in groups of at most LOOK_MAX candidates (window_groups).  The same
                    # sub-hits, a bounded number of them in memory at once
                    groups = [(0, len(mi_all))]
                    first = None
                    try:
                        self._look_limit = LOOK_MAX if tab2 else None
                        c2sub = tuple(a[mi_all] for a in c2all) if c2all is not None else None
                        first = self.Q(S2, k2, plo[mi_all], phi[mi_all], c2sub)
                    except TooMany:
                        groups = self.window_groups(S2, k2, plo, phi, mi_all)
                    finally:
                        self._look_limit = None
                    for (g0, g1) in groups:
                        mi = mi_all[g0:g1]
                        if first is not None:
                            sub = first
                        else:
                            c2sub = tuple(a[mi] for a in c2all) if c2all is not None else None   # the constraint on l
                            sub = self.Q(S2, k2, plo[mi], phi[mi], c2sub)
                        cs[1] += len(sub)
                        if not len(sub):
                            continue
                        pi = mi[sub.w]                              # index into the (W, I) pairs
                        h = self.finish_form(f, pair, pi, sub, plo[pi], phi[pi],
                                             tuple(a[pi] for a in c2all) if c2all is not None else None, cs)
                        if h is not None:
                            out.append(h)
            for g in grp:
                if not g:
                    continue
                PF = np.concatenate([np.full(len(t[1]), t[0], np.int8) for t in g])
                PP = np.concatenate([t[1] for t in g])
                PL = np.concatenate([t[2] for t in g]); PH = np.concatenate([t[3] for t in g])
                C2 = tuple(np.concatenate([t[4][k] for t in g]) for k in range(3)) if g[0][4] is not None else None
                for sub_all in self.rec_chunks(S2, k2, PL, PH, C2):
                    fw = PF[sub_all.w]
                    for f in forms:
                        sel = np.where(fw == f)[0]
                        if not len(sel):
                            continue
                        sub = sub_all.sel(sel)
                        cs = self.cstats[(len(S1), len(S2), k1, k2, f)]
                        cs[1] += len(sub)
                        j = sub.w                               # index into the windows of the group
                        h = self.finish_form(f, pair, PP[j], sub, PL[j], PH[j],
                                             tuple(a[j] for a in C2) if C2 is not None else None, cs)
                        if h is not None:
                            out.append(h)
                    del sub_all
        return Hits.cat(out)

    def rec_chunks(self, S2, k2, PL, PH, C2):
        """(14/11) the sub-hits of the recursive query Q(S2, k2) for the windows (PL, PH, C2), as successive Hits whose
        w index those windows, each with at most REC_PAIRS (window, sub-hit) pairs (more only when a single window has
        more).  Windows equal bit for bit (bounds and second-level constraint) are queried once, and their sub-hits are
        given to each of them chunk by chunk.  Before, Q (DEDUP) built every pair in one array: for split 100 of n = 9
        (piece NU, 18 entries of N(357) sharing their windows) 1.3e9 pairs, 9.7 GB at once.
        (14/11, b) The distinct windows are queried in batches whose size adapts to the sub-hits they return (about
        REC_PAIRS per batch, from REC_WIN0 windows, growing at most 8 times per batch): split 20 of n = 9 (piece VN.1)
        asked Q(S2) for all its distinct windows at once, tens of millions of sub-hits in one array.
        The pairs are the same (each window gets the sub-hits of its own query, which does not depend on the other
        windows of the batch); only their order changes, and finish_form treats each pair on its own, so combine
        returns the same hits (possibly in another order)."""
        rows = [PL, PH] if C2 is None else [PL, PH, C2[0].astype(np.float64), np.asarray(C2[1], np.float64),
                                            np.asarray(C2[2], np.float64)]
        rows = np.ascontiguousarray(np.stack(rows, 1))
        _, first, inv = np.unique(rows.view(np.dtype((np.void, rows.dtype.itemsize * rows.shape[1]))).ravel(),
                                  return_index=True, return_inverse=True)
        del rows
        inv = inv.ravel()
        nu = len(first)
        byu = np.argsort(inv, kind="stable")             # the windows grouped by distinct window, in order in each
        ucnt = np.bincount(inv, minlength=nu)
        uend = np.cumsum(ucnt)
        u0, nb = 0, min(nu, REC_WIN0)
        while u0 < nu:
            u1 = min(nu, u0 + nb)
            fu = first[u0:u1]
            hu = self.Q(S2, k2, PL[fu], PH[fu], tuple(a[fu] for a in C2) if C2 is not None else None)
            nh = len(hu)
            nb = (u1 - u0) * 8 if nh * 8 <= REC_PAIRS else max(1, int((u1 - u0) * REC_PAIRS / nh))
            if nh:
                wins = np.sort(byu[(uend[u0 - 1] if u0 else 0):uend[u1 - 1]])   # the windows of this batch
                il = inv[wins] - u0
                order = np.argsort(hu.w, kind="stable")
                cnt = np.bincount(hu.w, minlength=u1 - u0)
                start = np.cumsum(cnt) - cnt
                per = cnt[il]                             # sub-hits of each window
                cum = np.cumsum(per)
                n, a = len(wins), 0
                while a < n:
                    base = int(cum[a - 1]) if a else 0
                    b = max(int(np.searchsorted(cum, base + REC_PAIRS, side="right")), a + 1)
                    idx = np.arange(a, min(b, n))
                    w, pos = _ranges(start[il[idx]], start[il[idx]] + per[idx])
                    if len(w):
                        h = hu.sel(order[pos])
                        h.w = wins[idx[w]]
                        yield h
                        del h
                    a = int(idx[-1]) + 1
            del hu
            u0 = u1

    def finish_form(self, f, pair, pi, sub, plo, phi, c2l_all, cs):
        """the sub-hits of form f for the pairs pi of a chunk of combine(): forward enclosure (Lemma X2), the exact
        decisions and the filters; plo, phi (the preimage window of l) and c2l_all (its second-level constraint, or
        None), one per sub-hit.  Returns the Hits meeting their windows, or None."""
        S1, k1, S2, k2, e1, W, I, x7, a1, b1, o1, lo, hi, want, c2 = pair
        r, ok = fwd(f, x7.take(pi), sub.q7(),
                    (e1.P[I[pi]], e1.Q[I[pi]], np.where(sub.DS == 0, sub.AP, 0),
                     np.where(sub.DS == 0, sub.AQ, 0)) if f == 1 else None)
        A, B, Z, O = r.A, r.B, r.Z, r.O
        wo = W[pi]
        meet = ok & np.where(Z, lo[wo] == -INF, (B >= lo[wo]) & (A <= hi[wo]))
        if f >= 5:
            # a base that is an in-range value possibly equal to 1 (class or query enclosure with
            # 0 in its ln) gives no information with a huge exponent: decide it exactly (motor B)
            if f in (5, 6):
                bvz = (k1 == "V") & (a1[pi] <= 0) & (b1[pi] >= 0) & ~o1[pi]
            else:
                bvz = sub.isv & (sub.A <= 0) & (sub.B >= 0) & ~sub.O
            with np.errstate(invalid="ignore"):
                wide = meet & bvz & ~((B - A) < 1.0)
            if wide.any():
                # (03/11) for l^s, l^-s the exact members of the base l are only needed inside the
                # preimage window of l (Lemma X3) and its second-level constraint (pre_level2)
                lw_ = (plo, phi) if (NARROW_REFINE and f in (7, 8)) else None
                c2l = c2l_all if (NARROW_REFINE and f in (7, 8)) else None
                meet = self.refine_base1(f, S1, S2, x7.take(pi), sub, wide, meet, lo[wo], hi[wo], r, lw_, c2l)
        if f == 1 and k1 == "V" and k2 == "V" and (want == "out" or (want == "in" and RESOLVE_IN)):
            amb = meet & (A == -INF) & ~Z
            if amb.any():
                meet = self.resolve_vv(S1, S2, I[pi], sub, amb, meet, lo[wo], hi[wo], r, want)
        if want == "in":
            meet &= ivl.possibly_in(A, B, Z)
        elif want == "out":
            meet &= ivl.possibly_out(A, B, Z)
        if c2 is not None and meet.any():
            # (03/11) the windows' second-level constraint, applied as the hits are produced (the same
            # filter as at the end of Q; applying it early only keeps fewer candidates in memory)
            meet &= ~c2_drop(r, wo, c2)
        cs[2] += int(meet.sum())
        if not meet.any():
            return None
        pid = None
        if DEBUG:
            base = PROV.add(("pair", S1, k1, f, S2, k2), I[pi][meet], sub.pid[meet])
            pid = np.arange(base, base + int(meet.sum()))
        return Hits(wo[meet], r.take(meet), pid)

    def finish_pairs(self, f, S1, k1, S2, k2, e1, I, sub, wo, lo, hi, want, cs):
        """the pairs (e1[I], sub) of form f found for the windows wo: forward enclosure (Lemma X2), the exact
        decisions of C.5 and the filters, exactly as in combine(); returns the Hits meeting their windows"""
        x7 = e1.q7(I)
        r, ok = fwd(f, x7, sub.q7(),
                    (e1.P[I], e1.Q[I], np.where(sub.DS == 0, sub.AP, 0), np.where(sub.DS == 0, sub.AQ, 0)) if f == 1 else None)
        A, B, Z, O = r.A, r.B, r.Z, r.O
        meet = ok & np.where(Z, lo[wo] == -INF, (B >= lo[wo]) & (A <= hi[wo]))
        if f >= 5:
            # a base that is an in-range value possibly equal to 1 (class or query enclosure with 0 in its ln)
            # gives no information with a huge exponent: decide it exactly
            if f in (5, 6):
                bvz = (k1 == "V") & (e1.A[I] <= 0) & (e1.B[I] >= 0) & ~e1.O[I]
            else:
                bvz = sub.isv & (sub.A <= 0) & (sub.B >= 0) & ~sub.O
            with np.errstate(invalid="ignore"):
                wide = meet & bvz & ~((B - A) < 1.0)
            if wide.any():
                meet = self.refine_base1(f, S1, S2, x7, sub, wide, meet, lo[wo], hi[wo], r)
        if f == 1 and k1 == "V" and k2 == "V" and (want == "out" or (want == "in" and RESOLVE_IN)):
            amb = meet & (A == -INF) & ~Z
            if amb.any():
                meet = self.resolve_vv(S1, S2, I, sub, amb, meet, lo[wo], hi[wo], r, want)
        if want == "in":
            meet &= ivl.possibly_in(A, B, Z)
        elif want == "out":
            meet &= ivl.possibly_out(A, B, Z)
        cs[2] += int(meet.sum())
        if not meet.any():
            return None
        pid = None
        if DEBUG:
            base = PROV.add(("pair", S1, k1, f, S2, k2), I[meet], sub.pid[meet])
            pid = np.arange(base, base + int(meet.sum()))
        return Hits(wo[meet], r.take(meet), pid)

    def combine_batched(self, S1, k1, e1, S2, k2, forms, lo, hi, want, widx):
        """combine() with one recursive query per chunk: the preimage windows of every form and branch
        (Lemma X3) are concatenated and answered by a single Q(S2, k2, ...); each answer is then sent back to
        its form and pair.  The pairs and the treatment of each pair are those of combine()."""
        n1 = len(e1)
        out = []
        step = max(1, CHUNK // (n1 * max(1, len(forms))))
        for c0 in range(0, len(widx), step):
            wv = widx[c0:c0 + step]
            W = np.repeat(wv, n1); I = np.tile(np.arange(n1), len(wv))
            a1, b1, z1, o1 = e1.A[I], e1.B[I], e1.Z[I], e1.O[I]
            PL, PH, PF, PP = [], [], [], []
            x7b = e1.q7(I)
            for f in forms:
                for (plo, phi) in pre2(f, x7b, lo[W], hi[W]):
                    m = np.where(plo <= phi)[0]
                    if len(m):
                        PL.append(plo[m]); PH.append(phi[m]); PF.append(np.full(len(m), f, np.int8)); PP.append(m)
            if not PL:
                continue
            PL = np.concatenate(PL); PH = np.concatenate(PH); PF = np.concatenate(PF); PP = np.concatenate(PP)
            sub_all = self.Q(S2, k2, PL, PH)
            fw = PF[sub_all.w] if len(sub_all) else np.empty(0, np.int8)
            for f in forms:
                ck = (len(S1), len(S2), k1, k2, f)
                cs = self.cstats.setdefault(ck, [0, 0, 0])
                cs[0] += int((PF == f).sum())
                sel = np.where(fw == f)[0]
                cs[1] += len(sel)
                if not len(sel):
                    continue
                sub = sub_all.sel(sel)
                pi = PP[sub.w]                                  # index into the (W, I) pairs
                h = self.finish_pairs(f, S1, k1, S2, k2, e1, I[pi], sub, W[pi], lo, hi, want, cs)
                if h is not None:
                    out.append(h)
        return Hits.cat(out)

    def combine_swapped(self, S1, k1, e1, S2, k2, e2, isv2, ti2, nv2, forms, lo, hi, want):
        """combine() with the roles of the two tables exchanged, when the query side S2 is itself a table much
        smaller than the table of S1: for every window, every entry l of the S2 table and every form f, the
        preimage window of s is pre(SW[f], l, w) (Lemma X3 for the form f'(l, s) = f(s, l)), and the S1 table
        answers it by lookup (all forms and branches of a chunk in one lookup).  The pairs found contain every
        pair (s, l) with f(s, l) in w, and each pair is then treated exactly as in combine() (finish_pairs)."""
        ix1 = self.index(S1, k1)
        assert ix1.e is e1
        n2 = len(e2)
        widx = np.where(lo <= hi)[0]
        if not len(widx) or not n2:
            return Hits.cat([])
        out = []
        step = max(1, CHUNK // (n2 * max(1, len(forms))))
        for c0 in range(0, len(widx), step):
            wv = widx[c0:c0 + step]
            W = np.repeat(wv, n2); J = np.tile(np.arange(n2), len(wv))
            x7j = e2.q7(J)
            PL, PH, PF, PP = [], [], [], []
            for f in forms:
                # (10/11) pre2: for f' = SW[f] in (s'^l', s'^-l'), i.e. l^s and l^-s, the window of s is also the one
                # given by the second level of the base l (the entry of S2), as in combine() for the bases s
                for (plo, phi) in pre2(SW[f], x7j, lo[W], hi[W]):
                    m = np.where(plo <= phi)[0]
                    if len(m):
                        PL.append(plo[m]); PH.append(phi[m]); PF.append(np.full(len(m), f, np.int8)); PP.append(m)
            if not PL:
                continue
            PL = np.concatenate(PL); PH = np.concatenate(PH); PF = np.concatenate(PF); PP = np.concatenate(PP)
            for f in forms:
                cs = self.cstats.setdefault((len(S1), len(S2), k1, k2, f, "sw"), [0, 0, 0])
                cs[0] += int((PF == f).sum())
            # (10/11) the lookup by groups of windows of at most LOOK_MAX candidates, and the pairs of each form by
            # blocks of CHUNK: the pairs are the same, never more than a bounded number of them in memory at once
            try:
                batches = [ix1.look(PL, PH, LOOK_MAX)] if isinstance(ix1, IndexT) else [ix1.look(PL, PH)]
            except TooMany:
                batches = ix1.look_batches(PL, PH, LOOK_MAX)
            for ww, Iall in batches:                                # s = e1[I] meets the preimage window ww
                fw = PF[ww]
                for f in forms:
                    cs = self.cstats[(len(S1), len(S2), k1, k2, f, "sw")]
                    sel_all = np.where(fw == f)[0]
                    cs[1] += len(sel_all)
                    for b0 in range(0, len(sel_all), CHUNK):
                        sel = sel_all[b0:b0 + CHUNK]
                        pj = PP[ww[sel]]                          # index into the (W, J) pairs
                        I = Iall[sel]; Jl = J[pj]; wo = W[pj]
                        sub = Hits(np.zeros(len(Jl), np.int64), e2.q7(Jl), None, ti2[Jl])
                        sub.isv = isv2[Jl]
                        if DEBUG:
                            sub.pid = np.full(len(Jl), -1, np.int64)
                            for kind, sl, off in (("V", Jl < nv2, 0), ("N", Jl >= nv2, nv2)):
                                if sl.any():
                                    base = PROV.add(("tab", S2, kind), Jl[sl] - off)
                                    sub.pid[sl] = np.arange(base, base + int(sl.sum()))
                        h = self.finish_pairs(f, S1, k1, S2, k2, e1, I, sub, wo, lo, hi, want, cs)
                        if h is not None:
                            out.append(h)
        return Hits.cat(out)


def _ball1200(self, v):
    """Arb ball of an exact value at 1200 bits (cached per value object)"""
    from flint import ctx
    try:
        return v.ball(1200)
    except Exception:
        old = ctx.prec; ctx.prec = 1200
        try:
            return ex_tight(v)
        finally:
            ctx.prec = old


def ex_tight(v):
    import exact as ex
    return ex.tight_ball(v)


QX._ball1200 = _ball1200


def _resolve_vv(self, S1, S2, I, sub, amb, meet, wlo, whi, r, want="out"):
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
    if want == "out" and NEAR_X:
        # (16/11) a recursive l against a tiny window (ln u <= NEAR_X_LN): only the members of l within u of a member
        # x of the class of s can give a difference in the window (near_x_ys); one list per class of s, filtered below
        # by the ball of each l (as the group query of 03/11)
        rec = [t for t in np.where(amb)[0] if int(sub.ti[t]) < 0 and np.isfinite(whi[t]) and whi[t] <= NEAR_X_LN]
        by_cls = {}
        for t in rec:
            by_cls.setdefault(int(I[t]), []).append(t)
        for i_cls, ts in by_cls.items():
            uhi = float(max(whi[t] for t in ts))
            vals = self.near_x_ys(S1, e1, i_cls, S2, uhi)
            if vals is None:
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
                    balls.append((v, e, -np.inf, np.inf))
                else:
                    balls.append((v, e, float(np.nextafter(float(lb.lower()), -np.inf)), float(np.nextafter(float(lb.upper()), np.inf))))
            for t in ts:
                pre_ys[t] = (balls, "nearx", uhi)
            self.nearx_vv = getattr(self, "nearx_vv", 0) + len(ts)
    if want == "out" and GROUP_VV:
        rec = [t for t in np.where(amb)[0] if int(sub.ti[t]) < 0 and np.isfinite(sub.A[t]) and np.isfinite(sub.B[t])
               and float(sub.B[t]) - float(sub.A[t]) < 1e-6 and t not in pre_ys]
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
                pre_ys[t] = (balls, A0, B0)        # the values of the pair t: those whose ball meets its enclosure
            self.grouped_vv = getattr(self, "grouped_vv", 0) + len(ts)
    # (10/11) the enclosures of the differences depend only on the two lists of exact values, not on the window:
    # they are computed once per pair of lists (vv_hull) and kept; a list is identified by what it is computed from
    # (the class of s; the class of l, or the exact query of the group and the enclosure of l that filters it)
    memo = self.__dict__.setdefault("_vvmemo", {})
    for t in np.where(amb)[0]:
        j = int(sub.ti[t])
        if j < 0 and want == "in":
            continue                                   # a recursive value (not a table class): kept as it is
        i = int(I[t])
        if j >= 0:
            key = (want, S1, i, S2, j)
        elif t in pre_ys:
            key = (want, S1, i, S2, pre_ys[t][1], pre_ys[t][2], float(sub.A[t]), float(sub.B[t]))
        else:
            key = (want, S1, i, S2, None, None, float(sub.A[t]), float(sub.B[t]))
        h = memo.get(key)
        if h is None:
            try:
                if j >= 0:
                    ys = self.members(S2, self.V(S2), j)
                elif t in pre_ys:
                    ys = [(v, e) for (v, e, lo_, hi_) in pre_ys[t][0] if hi_ >= sub.A[t] and lo_ <= sub.B[t]]
                else:
                    ys = self.exact_vals(S2, float(sub.A[t]), float(sub.B[t]))     # a recursive V value: its exact members
                xs = self.members(S1, e1, i)
            except Unresolved:
                self.unresolved_vv = getattr(self, "unresolved_vv", 0) + 1
                continue                               # undecided: the pair keeps its enclosure (a hit if it meets)
            if len(memo) >= VV_MEMO_MAX:
                memo.clear()
            h = memo[key] = self.vv_hull(xs, ys, want)
        lo_t, hi_t, und = h
        if und:                                        # a difference of undecided sign: nothing known below
            lo_t = -np.inf; hi_t = max(hi_t, float(r.B[t]))
        self.resolved_vv = getattr(self, "resolved_vv", 0) + 1
        if lo_t > hi_t or hi_t < wlo[t] or lo_t > whi[t]:
            meet[t] = False
        else:
            r.A[t], r.B[t] = max(lo_t, float(r.A[t])), min(hi_t, float(r.B[t]))
    return meet


NEAR_X = os.environ.get("QX_NEAR_X", "1") == "1"      # (16/11) resolve_vv "out": the members of l near each x only
NEAR_X_LN = -690.0                                      # applied to windows with ln u <= this (u below ~1e-300)


def _near_x_ys(self, S1, e1, i, S2, uhi):
    """(16/11) the exact members y of V(S2) such that |x - y| = u with ln u <= uhi for some exact member x of the class
    i of V(S1) (a superset), or None when this cannot be given cheaply (the caller then keeps its general query).
    For each x (0 < u <= x/2): |ln y - ln x| = |ln(1 +- u/x)| <= 2u/x = dl, so ln y lies in [ln x - dl, ln x + dl]
    (Arb, 256 bits, rounded outward to doubles).  When that window does not contain 0 it is the exact query;
    when it does (x within ~1e-300 of 1): if x is exactly 1, y = 1 +- u gives ln y of either sign with
    |ln y| <= 2u = dl, i.e. two queries with the second level (+-1, -inf, ln dl) (y = 1 itself gives |x - y| = 0,
    never a value of the window); if x is not exactly 1, None.
    Split 20 of n = 9: x and l near 1 and u ~ e^-240779; the general query asked for every member of V(345789)
    within 2e-12 of 1 (hours, > 9 GB) while only those within ~e^-240779 of x can matter."""
    import exact as ex
    from exact import Unresolved
    from flint import ctx, arb
    key = (tuple(S1), int(i), tuple(S2), float(uhi).hex())
    cache = self.__dict__.setdefault("_nearx", {})
    if key in cache:
        return cache[key]
    xs = self.members(S1, e1, i)
    out = {}
    res = None
    old = ctx.prec; ctx.prec = 256
    try:
        u = arb(float(uhi)).exp()                         # ln u <= uhi: u <= this (arb, no underflow)
        for (x, xe) in xs:
            xb = self._ball1200(x)
            if not (xb > 2 * u):                          # u <= x/2 needed (and x > 0)
                break
            lnx = xb.log()
            dl = 2 * u / xb
            A_ = float(np.nextafter(float((lnx - dl).lower()), -np.inf))
            B_ = float(np.nextafter(float((lnx + dl).upper()), np.inf))
            if A_ > 0 or B_ < 0:
                vals = self.exact_vals(S2, A_, B_) if len(S2) <= self.tv else self.exact_V(S2, A_, B_)
            else:
                try:
                    is_one = ex.sign(ex.sub(x, ex.ONE)) == 0
                except Unresolved:
                    is_one = False
                if not is_one or len(S2) <= self.tv:
                    break
                M = float(np.nextafter(float(dl.log().upper()), np.inf))
                vals = list(self.exact_V(S2, A_, B_, (1, -np.inf, M))) + list(self.exact_V(S2, A_, B_, (-1, -np.inf, M)))
            for (v, e) in vals:
                out.setdefault(v.key(), (v, e))
        else:
            res = list(out.values())
    except Unresolved:
        res = None
    finally:
        ctx.prec = old
    cache[key] = res
    return res


QX.near_x_ys = _near_x_ys
VV_MEMO_MAX = int(os.environ.get("QX_VV_MEMO_MAX", "400000"))   # pairs of lists kept by resolve_vv


def _vv_hull(self, xs, ys, want):
    """(lo, hi, und) for the exact values xs of s and ys of l: the hull [lo, hi] of the enclosures of ln|x - y| over
    the pairs whose difference is nonzero and of the wanted kind (out of range for "out", in range for "in"), empty
    (lo > hi) if there is none, and und = True if the sign of some difference stayed undecided.  The computation
    resolve_vv made for every pair of hits, as a function of the two lists and of want only."""
    import exact as ex
    from exact import Unresolved
    from flint import ctx, arb
    TH = arb(10) ** -299
    keep_fn = ivl.possibly_out if want == "out" else ivl.possibly_in
    lo_t, hi_t, und = np.inf, -np.inf, False
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
                und = True; continue
            b = ex.tight_ball(ex.absval(d))
            old = ctx.prec; ctx.prec = 256
            lb = b.log(); ctx.prec = old
            a_, b_ = np.nextafter(float(lb.lower()), -np.inf), np.nextafter(float(lb.upper()), np.inf)
            if keep_fn(np.array([a_]), np.array([b_]), np.array([False]))[0]:
                lo_t = min(lo_t, a_); hi_t = max(hi_t, b_)
    return (float(lo_t), float(hi_t), und)


QX.vv_hull = _vv_hull
QX.resolve_vv = _resolve_vv


def _exact_vals(self, S, A, B):
    """exact values (with writings) of V(S) with ln in [A, B]: a superset.  |S| <= tv: motor B (Lemma 14);
    |S| > tv: motor X's own recursion with exact leaves (exact_V), which needs no motor B query on S"""
    if EXACT_V and len(S) > self.tv:
        return self.exact_V(S, A, B)
    key = (tuple(S), A, B)
    if not hasattr(self, "_xv"):
        self._xv = {}
    if key not in self._xv:
        from flint import arb
        lo = float(ivl.exp_lo(np.array([A]))[0]); hi = float(ivl.exp_hi(np.array([B]))[0])
        vals = self.bquery(S, arb(np.nextafter(lo, 0)), arb(np.nextafter(hi, np.inf)))   # Unresolved: caller keeps
        self._xv[key] = list(vals.values())
    return self._xv[key]


def exact_form(f, a, ea, b, eb):
    """[(value, writing)] of form f on the exact magnitudes a (s) and b (l): empty if f is undefined there or
    the value is not in range (V(S) only holds values with every intermediate 0 or in range; ex.in_range
    keeps the values whose enclosure straddles a range bound, so this is a superset).  The same exact
    operations as motor B's forms (motorb.B.forms)."""
    import exact as ex
    if f == 0:
        v, e = ex.add(a, b), f"({ea}+{eb})"
    elif f == 1:
        v, e = ex.absval(ex.sub(a, b)), f"|{ea}-{eb}|"
    elif f == 2:
        v, e = ex.mul(a, b), f"({ea}*{eb})"
    elif f == 3:
        v, e = ex.div(a, b), f"({ea}/{eb})"
    elif f == 4:
        v, e = ex.div(b, a), f"({eb}/{ea})"
    elif f == 5:
        v, e = ex.power(a, b), f"({ea}^{eb})"
    elif f == 6:
        v, e = ex.power(a, ex.neg(b)), f"({ea}^-{eb})"
    elif f == 7:
        v, e = ex.power(b, a), f"({eb}^{ea})"
    else:
        v, e = ex.power(b, ex.neg(a)), f"({eb}^-{ea})"
    if v is None or not ex.in_range(v):
        return []
    return [(v, e)]


def _c2key(c2):
    """hashable form of a single-window second-level constraint (sign, lower, upper of ln|ln v|) or None"""
    if c2 is None or int(c2[0]) == 0:
        return None
    return (int(c2[0]), float(c2[1]).hex(), float(c2[2]).hex())


def _exact_V(self, S, A, B, c2=None):
    """exact values, with writings, of every member of V(S) whose ln may lie in [A, B], for |S| > tv (a superset).
    With c2 = (sign, lower, upper) (03/11): only the members whose ln has that sign and ln|ln v| in [lower, upper]
    are required (the query passes the constraint down the recursion, which only removes members outside it).
    Motor X's recursion for kind V runs in provenance mode (Theorem X: every value of V(S) in the window is
    f(a, b) for a returned derivation, with a in a table class of S1 and b covered by the derivation of the
    sub-query), and each returned derivation is evaluated exactly: the exact members of its table classes
    (motorx.members, Lemma 14), its literals, and motor B's exact arithmetic (exact_form).  Cached by the exact
    query: the digit set and the two window bounds (hex of the doubles).  Raises Unresolved when an exact
    operation cannot be decided (the caller then keeps the enclosure)."""
    global PROV, DEBUG
    import exact as ex
    from exact import Unresolved
    S = tuple(S)
    ck = _c2key(c2)
    key = (S, float(A).hex(), float(B).hex(), ck)
    cache = self.__dict__.get("_xV")
    if cache is None:
        import collections
        cache = self._xV = collections.OrderedDict(); self._xV_n = 0
    if key in cache:
        cache.move_to_end(key)
        r = cache[key]
        if isinstance(r, Unresolved):
            raise r
        return r
    saved = (PROV, DEBUG)
    PROV, DEBUG = _Prov(), True
    try:
        c2a = None if ck is None else (np.array([ck[0]], np.int8), np.array([float(c2[1])]), np.array([float(c2[2])]))
        H = self.Q(S, "V", np.array([float(A)]), np.array([float(B)]), c2a)
        memo, out = {}, {}
        for pid in H.pid:
            for (v, e) in self._xvals(int(pid), memo):
                out.setdefault(v.key(), (v, e))
        res = list(out.values())
    except Unresolved as u:
        res = u
    finally:
        PROV, DEBUG = saved
    # (13/11) a bounded cache: the least recently used results are dropped beyond XV_MAX exact values (some queries
    # return 10^5 values and several GB).  A result is a deterministic function of its key: recomputing it gives the
    # same values, so only time changes
    nres = len(res) if isinstance(res, list) else 1
    if nres <= XV_MAX:
        cache[key] = res
        self._xV_n += nres
        while self._xV_n > XV_MAX and len(cache) > 1:
            _, old = cache.popitem(last=False)
            self._xV_n -= len(old) if isinstance(old, list) else 1
    self.nxV = getattr(self, "nxV", 0) + 1
    if isinstance(res, Unresolved):
        raise res
    return res


def _xvals(self, pid, memo):
    """exact values of the derivation PROV[pid] of a V query (provenance mode)"""
    if pid in memo:
        return memo[pid]
    import exact as ex
    from fractions import Fraction as Fr
    p = PROV[pid]
    if p[0] == "zero":
        out = [(ex.ZERO, "0")]
    elif p[0] == "near1":                                  # a row of near1(S): motor B's value, re-evaluated
        out = [self.near1_value(p[1], p[2])]
    elif p[0] == "tab":
        _, S, kind, j = p
        assert kind == "V"
        if len(S) <= self.tv:
            out = list(self.members(S, self.V(S), j))
        else:                                              # a literal of S (QX.literals)
            e = self.index(S, "V").e
            out = [(ex.X.rat(Fr(int(e.P[j]), int(e.Q[j]))), f"{int(e.P[j])}/{int(e.Q[j])}")]
    else:
        _, S1, k1, ii, f, S2, k2, cp = p
        assert k1 == "V" and k2 == "V"
        A_ = self.members(S1, self.V(S1), ii)
        B_ = self._xvals(cp, memo)
        d = {}
        for (a, ea) in A_:
            for (b, eb) in B_:
                for (v, e) in exact_form(f, a, ea, b, eb):
                    d.setdefault(v.key(), (v, e))
        out = list(d.values())
    memo[pid] = out
    return out


QX.exact_V = _exact_V
QX._xvals = _xvals
EXACT_V = os.environ.get("QX_EXACT_V", "1") == "1"      # exact values of V(S), |S| > tv, without motor B on S


QX.exact_vals = _exact_vals


def _exact_near1(self, S, A, B, c2=None):
    """exact members (motor B) of V(S) with ln in [A, B]: list of (is_one, ln-ball lo, hi).  c2 (sign, lower, upper
    of ln|ln v|, optional): only the members satisfying it are required (a superset of them is returned)"""
    key = (tuple(S), float(A), float(B), _c2key(c2))
    if not hasattr(self, "_n1"):
        self._n1 = {}
    if key in self._n1:
        return self._n1[key]
    import exact as ex
    from exact import Unresolved
    from flint import arb, ctx
    if NEAR1 and self.tv < len(S) <= NEAR1_MAXSET and -DELTA1 <= A and B <= DELTA1:
        # inside the radius of near1(S): its rows hold every member of V(S) with |ln v| <= DELTA1, with exact
        # decisions (exact 1) and Arb enclosures of ln v; the rows meeting [A, B] contain the members there
        e = self.near1(S)
        if e is not None:
            m = np.where((e.B >= A) & (e.A <= B))[0]
            self._n1[key] = [(bool(e.O[i]), float(e.A[i]), float(e.B[i])) for i in m]
            self.n1_rows_used = getattr(self, "n1_rows_used", 0) + 1
            return self._n1[key]
    if EXACT_V and len(S) > self.tv:
        try:
            vlist = self.exact_V(S, A, B, c2)
        except Unresolved:                     # undecided: the base keeps its enclosure (no refinement)
            self._n1[key] = [(False, float(A), float(B))]
            return self._n1[key]
    else:
        lo = float(ivl.exp_lo(np.array([A]))[0]); hi = float(ivl.exp_hi(np.array([B]))[0])
        try:
            vlist = self.bquery(S, arb(np.nextafter(lo, 0)), arb(np.nextafter(hi, np.inf))).values()
        except Unresolved:
            self._n1[key] = [(False, float(A), float(B))]
            return self._n1[key]
    out = []
    for (v, e) in vlist:
        d = ex.sub(v, ex.ONE)
        try:
            sg = ex.sign(d)
        except Unresolved:
            sg = None
        if sg == 0:
            out.append((True, 0.0, 0.0)); continue
        lb = None
        for prec in (256, 1024, 4096, 16384, 65536):
            old = ctx.prec; ctx.prec = prec
            try:
                lb = v.ball(prec).log()
            finally:
                ctx.prec = old
            if not lb.contains(0):
                break
        out.append((False, np.nextafter(float(lb.lower()), -np.inf), np.nextafter(float(lb.upper()), np.inf)))
    self._n1[key] = out
    return out


def _refine_base1(self, f, S1, S2, x, sub, wide, meet, wlo, whi, r, lwin=None, c2l=None):
    """the hits `wide` of form f whose base is an in-range value possibly equal to 1 (its ln enclosure contains 0)
    and whose result is wide: the base is replaced by its exact members (exact_near1) and the hit is kept only
    if some member gives a value meeting the window.  For l^s, l^-s (f = 7, 8; base l from the query side),
    lwin = (lower, upper) per hit: the preimage window of l of the sub-query that returned it (Lemma X3: every l
    with f(s, l) in the window has ln l in it), and c2l = (sign, lower, upper) per hit: the second-level
    constraint of that sub-query (pre_level2).  The members are then only required inside the window and the
    constraint (03/11): a member outside them gives no value in the window for any s of the class."""
    meet = meet.copy()
    if not hasattr(self, "_n1arr"):
        self._n1arr = {}
    ts, kid, mems = [], [], {}
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
                    self.narrowed_base = getattr(self, "narrowed_base", 0) + 1
            if c2l is not None and int(c2l[0][t]) != 0:
                c2t = (int(c2l[0][t]), float(c2l[1][t]), float(c2l[2][t]))
        key = (tuple(S), float(A), float(B), _c2key(c2t))
        if key not in self._n1arr:
            mem = self.exact_near1(S, A, B, c2t)
            if mem:
                self._n1arr[key] = (np.array([m[0] for m in mem]), np.array([m[1] for m in mem]), np.array([m[2] for m in mem]))
            else:
                self._n1arr[key] = None
        q = self._n1arr[key]
        self.refined_base = getattr(self, "refined_base", 0) + 1
        if q is None:
            meet[t] = False
            continue
        ts.append(t); kid.append(mems.setdefault(key, (len(mems), q))[0])
    if not ts:
        return meet
    # (10/11) the outcome of a hit depends only on its list of members, its other operand and its window.  Many hits
    # share the three (the derivations of the base with the same narrowed enclosure, against the same entry and
    # window): each distinct combination is computed once, in blocks (ivl.fwd works entry by entry), and given to all
    # its hits.  Per combination: the hull of the forward enclosures of the members that meet the window
    ts = np.array(ts); kid = np.array(kid, np.float64)
    lists = [None] * len(mems)
    for (i_, q) in mems.values():
        lists[i_] = q
    oth = sub.q7() if f in (5, 6) else x
    rows = np.concatenate([kid[:, None], wlo[ts, None].astype(np.float64), whi[ts, None].astype(np.float64)] +
                          [np.ascontiguousarray(getattr(oth, k)[ts]).astype(np.float64).reshape(len(ts), -1) for k in ivl.Q7.FIELDS], 1)
    rows = np.ascontiguousarray(rows)
    _, first, inv = np.unique(rows.view(np.dtype((np.void, rows.dtype.itemsize * rows.shape[1]))).ravel(),
                              return_index=True, return_inverse=True)
    inv = inv.ravel()
    tg = ts[first]; kg = kid[first].astype(np.int64)             # one representative hit per combination
    some = np.zeros(len(tg), bool); lo_ = np.full(len(tg), np.inf); hi_ = np.full(len(tg), -np.inf)
    b0 = 0
    while b0 < len(tg):
        b1, tot = b0, 0
        while b1 < len(tg) and (b1 == b0 or tot + len(lists[kg[b1]][0]) <= CHUNK):
            tot += len(lists[kg[b1]][0]); b1 += 1
        cnt = np.array([len(lists[k_][0]) for k_ in kg[b0:b1]])
        idx = np.repeat(tg[b0:b1], cnt)
        q = ivl.Q7(np.concatenate([lists[k_][1] for k_ in kg[b0:b1]]), np.concatenate([lists[k_][2] for k_ in kg[b0:b1]]),
                   np.zeros(tot, bool), np.concatenate([lists[k_][0] for k_ in kg[b0:b1]]))
        other = oth.take(idx)
        rr, okk = ivl.fwd(f, q, other) if f in (5, 6) else ivl.fwd(f, other, q)
        with np.errstate(invalid="ignore"):
            m = okk & np.where(rr.Z, wlo[idx] == -np.inf, (rr.B >= wlo[idx]) & (rr.A <= whi[idx]))
        st = np.cumsum(cnt) - cnt
        some[b0:b1] = np.logical_or.reduceat(m, st)
        lo_[b0:b1] = np.minimum.reduceat(np.where(m, rr.A, np.inf), st)
        hi_[b0:b1] = np.maximum.reduceat(np.where(m, rr.B, -np.inf), st)
        b0 = b1
    sm = some[inv]
    meet[ts[~sm]] = False
    r.A[ts[sm]] = lo_[inv][sm]; r.B[ts[sm]] = hi_[inv][sm]
    self.refine_groups = getattr(self, "refine_groups", 0) + len(tg)
    return meet


QX.exact_near1 = _exact_near1
QX.refine_base1 = _refine_base1


def explain(q, pid, depth=0):
    p = PROV[pid]
    ind = "  " * depth
    if p[0] == "zero":
        print(f"{ind}exact zero (recursive V query)")
        return
    if p[0] == "near1":
        print(f"{ind}near1 V{''.join(map(str, p[1]))}[{p[2]}]")
        return
    if p[0] == "tab":
        _, S, kind, j = p
        e = q.index(S, kind).e
        line = f"{ind}{kind}{''.join(map(str, S))}[{j}] = [{e.A[j]:.6g}, {e.B[j]:.6g}] Z={e.Z[j]} O={e.O[j]}"
        print(line)
        if kind == "N":
            pr = e.prov[j]
            S1, S2, _ = q.splits(S)[pr[0]]
            print(f"{ind}   (materialised: split {S1}|{S2} form {ivl.FORMS[pr[1]]} I={pr[2]} k1={pr[4]} J={pr[3]} k2={pr[5]})")
        return
    _, S1, k1, ii, f, S2, k2, cp = p
    e1 = q.N(S1) if k1 == "N" else q.V(S1)
    print(f"{ind}pair {ivl.FORMS[f]}: s = {k1}{''.join(map(str, S1))}[{ii}] = [{e1.A[ii]:.6g}, {e1.B[ii]:.6g}] Z={e1.Z[ii]} O={e1.O[ii]}; l from {k2}{''.join(map(str, S2))}:")
    explain(q, cp, depth + 1)


# ---------------------------------------------------------------- optional timing of the expensive steps (QX_PROF=1)
PROF = {}


def _timed(name, fn):
    def w(*a, **k):
        t = time.time()
        try:
            return fn(*a, **k)
        finally:
            r = PROF.setdefault(name, [0, 0.0]); r[0] += 1; r[1] += time.time() - t
    return w


if os.environ.get("QX_PROF"):
    for _nm in ("resolve_vv", "exact_vals", "members", "refine_base1", "exact_near1", "N", "refine_near1"):
        setattr(QX, _nm, _timed(_nm, getattr(QX, _nm)))
