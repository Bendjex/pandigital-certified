"""
Motor X (C5): values of expressions with intermediates OUT of range.

Notation (magnitudes, Lemmas 1-3 of C4): U(S) all valid expressions over S, V(S) those with every
intermediate 0 or in [1e-300, 1e300] (the set certified by C3), NE(S) = U(S) \\ V(S).
Lemma X1 (generation): u in NE(S) is f(u1, u2) over a split S = S1 u S2 with u1 in NE(S1), or
u2 in NE(S2), or u out of range with u1 in V(S1), u2 in V(S2).
Tables (enclosures, ivl.py):
  Vt[S]  motor A's classes of T(S) (Theorem A: they cover V(S)), |S| <= TV;
  Nt[S]  a list of enclosures covering NE(S), built with Lemma X1 from Vt and Nt of the parts.
"""
import sys, os, itertools, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import paths; sys.path.insert(1, paths.SRC)
import xl, ivl
from ivl import fwd, INF

CH = 4_000_000          # pairs per chunk
NCACHE = os.environ.get("QX_NCACHE", "")     # directory of the persistent N-table cache ("": no cache)
ALLPROV = [False]                            # hitcheck: N tables keep all the derivations of each entry
COUNT_OLD_UNIQ = os.environ.get("QX_COUNT_OLD_UNIQ", "1") == "1"


class AllProv:
    """all derivations of each N-table entry (hitcheck), stored as rows + offsets: ap[j] = list of 6-tuples"""
    def __init__(self, rows, off):
        self.rows, self.off = rows, off

    def __getitem__(self, j):
        return [tuple(int(t) for t in r) for r in self.rows[self.off[j]:self.off[j + 1]]]

    def __len__(self):
        return len(self.off) - 1

    def __iter__(self):
        return (self[j] for j in range(len(self)))
_FP = {}


def code_fingerprint():
    """sha256 of the source files the N tables depend on (line endings normalised)"""
    if "code" not in _FP:
        import hashlib
        h = hashlib.sha256()
        here = os.path.dirname(os.path.abspath(__file__))
        for d, names in ((here, ("ivl.py", "motorx.py", "xl.py", "xlib.c")),
                         (paths.SRC, ("motora.c", os.path.join("coremath", "pow.c"), "exact.py", "motorb.py", "canon.py",
                                      "exprparse.py"))):
            for nm in names:
                b = open(os.path.join(d, nm), "rb").read().replace(b"\r\n", b"\n")
                if nm == "motorx.py":            # the cache input/output itself is not part of the computation
                    import re as _re
                    b = _re.sub(rb"# <cache-io>.*?# </cache-io>", b"", b, flags=_re.S)
                h.update(nm.encode()); h.update(b)
        _FP["code"] = h.hexdigest()
    return _FP["code"]


class Ent:
    """a list of enclosures (ivl.Q7 fields) with provenance; P, Q: exact rational (Q = 0: none)"""
    F = ivl.Q7.FIELDS

    def __init__(self, A, B, Z, O, prov=None, P=None, Q=None, S=None, M=None, N=None, **kw):
        n = len(A)
        P = P if P is not None else np.zeros(n, np.int64)
        Q = Q if Q is not None else np.zeros(n, np.int64)
        q = ivl.Q7(A, B, Z, O, S, M, N, kw.get("AP"), kw.get("AQ"), kw.get("DS"), kw.get("DM"), kw.get("DN"))
        if kw.get("AQ") is None:                 # exact rationals become exact anchors
            q.AP = P.copy(); q.AQ = Q.copy()
        for k in self.F:
            setattr(self, k, getattr(q, k))
        self.prov = prov if prov is not None else np.zeros((n, 6), np.int32) - 1

    @property
    def P(self):
        return np.where(self.DS == 0, self.AP, 0)

    @property
    def Q(self):
        return np.where(self.DS == 0, self.AQ, 0)

    def __len__(self):
        return len(self.A)

    @staticmethod
    def from_q7(q, prov=None):
        e = Ent.__new__(Ent)
        for k in Ent.F:
            setattr(e, k, getattr(q, k))
        e.prov = prov if prov is not None else np.zeros((len(q.A), 6), np.int32) - 1
        return e

    @staticmethod
    def cat(lst):
        lst = [e for e in lst if len(e)]
        if not lst:
            return Ent(np.empty(0), np.empty(0), np.empty(0, bool), np.empty(0, bool), np.empty((0, 6), np.int32))
        e = Ent.__new__(Ent)
        for k in Ent.F + ("prov",):
            setattr(e, k, np.concatenate([getattr(x, k) for x in lst]))
        return e

    def q7(self, idx=None):
        if idx is None:
            return ivl.Q7(*[getattr(self, k) for k in self.F])
        return ivl.Q7(*[getattr(self, k)[idx] for k in self.F])

    def sel(self, m):
        e = Ent.__new__(Ent)
        for k in Ent.F + ("prov",):
            setattr(e, k, getattr(self, k)[m])
        return e

    def uniq_key(self):
        """exact key of an entry: the raw bytes of every field (corrected 03/11: the key of 01/11 left out the
        third level T3, M3, N3 and compared AP, AQ as doubles, so two different entries could be merged)"""
        n = len(self)
        return np.concatenate([np.ascontiguousarray(getattr(self, k)).view(np.uint8).reshape(n, -1) for k in self.F], 1)

    def uniq(self):
        if not len(self):
            return self
        _, idx = np.unique(self.uniq_key(), axis=0, return_index=True)
        if COUNT_OLD_UNIQ:                    # instrumentation: would the key of 01/11 have merged different entries?
            old = np.stack([self.A, self.B, self.Z.astype(float), self.O.astype(float), self.S.astype(float),
                            self.M, self.N, self.AP.astype(float), self.AQ.astype(float), self.DS.astype(float),
                            self.DM, self.DN], 1)
            _, iold = np.unique(old, axis=0, return_index=True)
            ivl.bump("uniq_old_lost_entries", len(idx) - len(iold))
        return self.sel(np.sort(idx))


def bits(S):
    m = 0
    for d in S:
        m |= 1 << (int(d) - 1)
    return m


def ln_bounds(v):
    """outward double bounds (lo, hi) of ln v for an exact value v >= 0 of motor B: (0, 0) for the exact 1,
    (-inf, -inf) for the exact 0; otherwise from an Arb ball of ln v, at the first precision whose relative radius
    is <= 2^-60 (the best one obtained otherwise), each bound rounded to the nearest double and moved one double
    outwards.  (-inf, inf) when nothing can be certified: such a value passes every test on its ln."""
    import exact as ex
    from flint import arb, ctx
    try:
        if ex.sign(v) == 0:
            return (-INF, -INF)
        if ex.sign(ex.sub(v, ex.ONE)) == 0:
            return (0.0, 0.0)
    except ex.Unresolved:
        return (-INF, INF)
    best = None
    for prec in (256, 1024, 4096, 16384):
        old = ctx.prec; ctx.prec = prec
        try:
            try:
                b = v.ball(prec)
            except ex.Unresolved:
                continue
            if not (b.is_finite() and b > 0):
                continue
            lb = b.log()
            if not lb.is_finite():
                continue
            best = lb
            if not lb.contains(0) and lb.rad() <= abs(lb.mid()) * (arb(2) ** -60):
                break
        finally:
            ctx.prec = old
    if best is None:
        return (-INF, INF)
    return (float(np.nextafter(float(best.lower()), -INF)), float(np.nextafter(float(best.upper()), INF)))


class MX:
    def __init__(self, n, ordered, res, tv):
        self.n, self.ord, self.tv = n, bool(ordered), tv
        t0 = time.time()
        xl.init(n, int(ordered), res, tv)
        self.Vt, self.Nt = {}, {}
        self.log = []
        self.say(f"motor A tables up to {tv}: {time.time() - t0:.1f}s")

    def say(self, s):
        self.log.append(s)
        print(s, flush=True)

    # ------------------------------------------------ sets and splits
    def sets(self, k):
        D = tuple(range(1, self.n + 1))
        if self.ord:
            return [D[i:i + k] for i in range(self.n - k + 1)]
        return list(itertools.combinations(D, k))

    def splits(self, S):
        """(S1, S2, forms): every split once, with the forms (s from S1, l from S2) of Lemma 3"""
        out = []
        if self.ord:
            for i in range(1, len(S)):
                out.append((S[:i], S[i:], ivl.ORD_LEFT))
        else:
            for r in range(1, len(S) // 2 + 1):
                for S1 in itertools.combinations(S, r):
                    S2 = tuple(d for d in S if d not in S1)
                    if 2 * r == len(S) and S1 > S2:
                        continue
                    out.append((S1, S2, tuple(range(9))))
        return out

    # ------------------------------------------------ tables
    def V(self, S):
        S = tuple(S)
        if S not in self.Vt:
            assert len(S) <= self.tv
            A, B, fl = xl.table(bits(S))
            P, Q = xl.rats(bits(S))
            e = self.dedup(Ent(A, B, (fl & 1) > 0, (fl & 2) > 0, None, P, Q))
            self.Vt[S] = self.refine_near1(S, e)
        return self.Vt[S]

    @staticmethod
    def dedup(e):
        """(10/11) one entry per distinct statement.  Motor A's classes are finer than the double enclosures of the
        table, so several classes can give entries equal in every field (62 of them at 4 + tiny in V(3469)).  Equal
        entries state the same thing, and the exact decisions depend only on the fields of the entry
        (class_members), so one of them stands for all the values the others stood for."""
        if len(e) < 2:
            return e
        _, idx = np.unique(e.uniq_key(), axis=0, return_index=True)
        if len(idx) == len(e):
            return e
        return e.sel(np.sort(idx))

    def refine_near1(self, S, e):
        """classes whose ln-enclosure contains 0 (value possibly 1) but are not the exact 1: replaced by
        their exact members (motor B): the exact 1, or precise enclosures with the sign of ln known"""
        amb = np.where((e.A <= 0) & (e.B >= 0) & ~e.O & ~e.Z)[0]
        if not len(amb):
            return e
        import exact as ex
        from exact import Unresolved
        from flint import arb, ctx
        keep = np.ones(len(e), bool); keep[amb] = False
        A, B, Z, O = [], [], [], []
        for i in amb:
            for (v, expr) in self.members(S, e, int(i)):
                d = ex.sub(v, ex.ONE)
                try:
                    sg = ex.sign(d)
                except Unresolved:
                    sg = None
                if sg == 0:
                    A.append(0.0); B.append(0.0); Z.append(False); O.append(True); continue
                lb = None
                for prec in (256, 1024, 4096, 16384):
                    old = ctx.prec; ctx.prec = prec
                    try:
                        lb = ex.tight_ball(v).log() if False else v.ball(prec).log()
                    finally:
                        ctx.prec = old
                    if not lb.contains(0):
                        break
                lo, hi = float(lb.lower()), float(lb.upper())
                A.append(np.nextafter(lo, -np.inf)); B.append(np.nextafter(hi, np.inf)); Z.append(False); O.append(False)
                if lb.contains(0):
                    self.say(f"near-1 class of {S} not separated from 1: {expr}")
        add = Ent(np.array(A), np.array(B), np.array(Z, bool), np.array(O, bool))
        self.say(f"V{''.join(map(str, S))}: {len(amb)} classes near 1 refined into {len(A)} exact members")
        return Ent.cat([e.sel(keep), add])

    def pairs(self, e1, e2, forms, keep, tag):
        """all f(u1, u2), u1 in e1, u2 in e2, f in forms; keep(A, B, Z, O) selects results"""
        out = []
        n1, n2 = len(e1), len(e2)
        if not n1 or not n2:
            return Ent.cat([])
        step = max(1, CH // max(1, n2))
        for i0 in range(0, n1, step):
            i1 = min(n1, i0 + step)
            I = np.repeat(np.arange(i0, i1), n2); J = np.tile(np.arange(n2), i1 - i0)
            for f in forms:
                r, ok = fwd(f, e1.q7(I), e2.q7(J), (e1.P[I], e1.Q[I], e2.P[J], e2.Q[J]) if f == 1 else None)
                m = ok & keep(r.A, r.B, r.Z, r.O)
                if m.any():
                    prov = np.stack([np.full(m.sum(), tag[0]), np.full(m.sum(), f), I[m], J[m],
                                     np.full(m.sum(), tag[1]), np.full(m.sum(), tag[2])], 1).astype(np.int32)
                    out.append(Ent.from_q7(r.take(m), prov))
        return Ent.cat(out)

    def ncache_path(self, S):
        """file of the persistent cache of N(S): the name holds a fingerprint of everything N(S) is computed
        from (the source code, the mode, tV and the V tables of the subsets of S), so a stored table is
        only ever used where the same deterministic computation would give it again"""
        import hashlib
        h = hashlib.sha256()
        h.update(f"{code_fingerprint()} ord={int(self.ord)} tv={self.tv} S={S} allprov={ALLPROV[0]} ancpow={ivl.ANCPOW}".encode())
        for k in range(1, len(S)):
            for T in itertools.combinations(S, k):
                if k <= self.tv:
                    e = self.V(T)
                    for f in ("A", "B", "Z", "O", "AP", "AQ", "DS"):
                        h.update(np.ascontiguousarray(getattr(e, f)).tobytes())
        tag = "H" if ALLPROV[0] else ""
        return os.path.join(NCACHE, f"N{tag}{int(self.ord)}_{''.join(map(str, S))}_{h.hexdigest()[:20]}.npz")

    def N(self, S):
        """materialised enclosures covering NE(S)"""
        S = tuple(S)
        if S in self.Nt:
            return self.Nt[S]
        t0 = time.time()
        # <cache-io>
        use_cache = bool(NCACHE) and len(S) >= 4
        if use_cache:
            import ncache
            fn = self.ncache_path(S)[:-4]                # a directory (ncache.py)
            try:
                z = ncache.load(fn, Ent.F + ("prov",), ("ap_rows", "ap_off"))
            except Exception as er:                     # a damaged table (e.g. power cut): recompute
                z = None
                self.say(f"N cache: {fn} unreadable ({er!r}), recomputed")
            if z is not None and (not ALLPROV[0] or "ap_rows" in z):
                e = Ent.__new__(Ent)
                for k in Ent.F + ("prov",):
                    setattr(e, k, z[k])
                if ALLPROV[0]:
                    e.allprov = AllProv(z["ap_rows"], z["ap_off"])
                e.cache_dir = fn
                self.Nt[S] = e
                self.say(f"N{''.join(map(str, S))}: {len(e)} enclosures (from the cache), {time.time() - t0:.1f}s")
                return e
        # </cache-io>
        parts = []
        if len(S) >= 2:
            for si, (S1, S2, forms) in enumerate(self.splits(S)):
                V1, V2 = self.V(S1), self.V(S2)
                N1, N2 = self.N(S1), self.N(S2)
                every = lambda A, B, Z, O: np.ones(len(A), bool)
                # NE on the left (with anything on the right) and NE on the right (with V on the left)
                parts.append(self.pairs(N1, Ent.cat([V2, N2]), forms, every, (si, 1, 2)))
                parts.append(self.pairs(V1, N2, forms, every, (si, 0, 1)))
                # V x V out of range; |s-l| of overlapping classes is decided exactly (motor B)
                fo = tuple(f for f in forms if f != 1)
                parts.append(self.pairs(V1, V2, fo, lambda A, B, Z, O: ivl.possibly_out(A, B, Z), (si, 0, 0)))
                if 1 in forms:
                    d = self.pairs(V1, V2, (1,), lambda A, B, Z, O: ivl.possibly_out(A, B, Z), (si, 0, 0))
                    amb = np.isinf(d.A)
                    parts.append(d.sel(~amb))
                    if amb.any():
                        parts.append(self.resolve_sub(S1, S2, d.prov[amb], si))
        e = Ent.cat(parts).uniq()
        self.Nt[S] = e
        # <cache-io>
        if use_cache:
            os.makedirs(NCACHE, exist_ok=True)
            arrays = {k: getattr(e, k) for k in Ent.F + ("prov",)}
            if ALLPROV[0]:
                rows = [r for g in e.allprov for r in g]
                arrays["ap_rows"] = np.array(rows, np.int32).reshape(-1, 6)
                arrays["ap_off"] = np.cumsum([0] + [len(g) for g in e.allprov]).astype(np.int64)
            ncache.save(fn, arrays)
            e.cache_dir = fn
        # </cache-io>
        if len(S) >= 3:
            self.say(f"N{''.join(map(str, S))}: {len(e)} enclosures "
                     f"({int(ivl.possibly_in(e.A, e.B, e.Z).sum())} possibly in range), {time.time() - t0:.1f}s")
        return e

    # ------------------------------------------------ exact resolution (motor B)
    def B(self):
        if not hasattr(self, "_B"):
            sys.path.insert(0, paths.SRC)
            import motorb
            self._B = motorb.B(ordered=self.ord, fullmax=int(os.environ.get("QX_BFULL", "3")))
            self.unresolved = []
        return self._B

    def members(self, S, ent, i):
        """exact values (motor B) of the entry i of Vt[S], with their expressions (class_members, cached)"""
        key = (S, i)
        if not hasattr(self, "_mem"):
            self._mem = {}
        if key not in self._mem:
            self._mem[key] = self.class_members(S, ent, i)
        return self._mem[key]

    def class_members(self, S, ent, i):
        """(10/11) the exact values, with expressions, that the entry i of the table ent of V(S) may stand for: a
        superset of M(i) = {a in V(S): the statement of the entry holds for a}.
          - an exact rational (Q > 0): that rational;
          - any other entry: motor B's values in the value window [e^A, e^B] rounded outwards (Lemma 14: every
            value of V(S) there), without those whose ln is certainly outside [A, B] (ln_bounds): the statement of
            the entry, A <= ln u <= B, is false for them.  Near 1 the value window is [pred(1), succ(1)] whatever
            A and B are, so without this test every class within 1e-16 of 1 returned all its neighbours.
        Every a in V(S) has an entry whose statement holds for a (Theorem A for the classes of motor A; an entry
        made by refine_near1 from a states the exact 1 or the ball of ln a), so a is in the list of that entry:
        what the exact decisions need."""
        if ent.Q[i] > 0:
            import exact as ex
            from fractions import Fraction as Fr
            return [(ex.X.rat(Fr(int(ent.P[i]), int(ent.Q[i]))), f"{ent.P[i]}/{ent.Q[i]}")]
        lo = float(ivl.exp_lo(np.array([ent.A[i]]))[0]); hi = float(ivl.exp_hi(np.array([ent.B[i]]))[0])
        lo, hi = float(np.nextafter(lo, 0)), float(np.nextafter(hi, np.inf))
        vals = self.bwindow(S, lo, hi)
        A, B = float(ent.A[i]), float(ent.B[i])
        return [(v, e) for (v, e, l, h) in vals if h >= A and l <= B]

    BWIN_MAX = int(os.environ.get("QX_BWIN_MAX", "400000"))     # exact values kept by the cache of bwindow()

    def bwindow(self, S, lo, hi):
        """motor B's values of V(S) in the value window [lo, hi] (doubles), each with outward double bounds of its
        ln: [(value, expression, ln_lo, ln_hi)].  Cached by the window (least recently used windows are dropped
        beyond BWIN_MAX values and recomputed if needed: a deterministic function of (S, lo, hi))."""
        import collections
        from flint import arb
        c = self.__dict__.get("_bwin")
        if c is None:
            c = self._bwin = collections.OrderedDict(); self._bwin_n = 0
        key = (S, lo, hi)
        if key in c:
            c.move_to_end(key)
            return c[key]
        r = list(self.B().query(tuple(str(d) for d in S), arb(lo), arb(hi)).values())
        out = [(v, e) + ln_bounds(v) for (v, e) in r]
        c[key] = out
        self._bwin_n += len(out) + 1
        while self._bwin_n > self.BWIN_MAX and len(c) > 1:
            _, old = c.popitem(last=False)
            self._bwin_n -= len(old) + 1
        return out

    def resolve_sub(self, S1, S2, prov, si):
        """|s - l| for overlapping classes of V(S1), V(S2): exact differences; only out-of-range nonzero
        differences are new values (an exact zero, or a difference in range, is in V)"""
        import exact as ex
        from exact import Unresolved
        from flint import arb, ctx
        V1, V2 = self.V(S1), self.V(S2)
        A, B, Z, O, P = [], [], [], [], []
        nz = 0
        for row in prov:
            i, j = int(row[2]), int(row[3])
            for (x, xe) in self.members(S1, V1, i):
                for (y, ye) in self.members(S2, V2, j):
                    d = ex.sub(x, y)
                    try:
                        sg = ex.sign(d)
                    except Unresolved:
                        sg = None
                    if sg == 0:
                        continue
                    if sg is None:
                        self.unresolved.append((S1, xe, S2, ye))
                        A.append(-INF); B.append(float(max(V1.B[i], V2.B[j]))); Z.append(False); O.append(False); P.append(row)
                        continue
                    b = ex.tight_ball(ex.absval(d))
                    old = ctx.prec; ctx.prec = 256
                    lb = b.log(); ctx.prec = old
                    lo, hi = float(lb.lower()), float(lb.upper())
                    lo, hi = np.nextafter(lo, -np.inf), np.nextafter(hi, np.inf)
                    if ivl.possibly_out(np.array([lo]), np.array([hi]), np.array([False]))[0]:
                        nz += 1
                        A.append(lo); B.append(hi); Z.append(False); O.append(False); P.append(row)
        if A:
            return Ent(np.array(A), np.array(B), np.array(Z), np.array(O), np.array(P, np.int32))
        return Ent.cat([])

    def cell_hits(self, S, wlo, whi):
        """enclosures covering NE(S) that meet the window [wlo, whi] (ln) - all splits of S, forward"""
        S = tuple(S)
        hits = []
        keepw = lambda A, B, Z, O: ~Z & (B >= wlo) & (A <= whi)
        for si, (S1, S2, forms) in enumerate(self.splits(S)):
            V1, V2 = self.V(S1), self.V(S2)
            N1, N2 = self.N(S1), self.N(S2)
            hits.append(self.pairs(N1, Ent.cat([V2, N2]), forms, keepw, (si, 1, 2)))
            hits.append(self.pairs(V1, N2, forms, keepw, (si, 0, 1)))
        return Ent.cat(hits)
