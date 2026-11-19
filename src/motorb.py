"""
Motor B: exact value sets and window queries (C3), any order and in order.

Values are nonnegative representatives |val E| (Lemma 1), powers have nonnegative bases (Lemma 2)
and the forms of Lemma 3 are used.  Only in-range values are kept (all subexpressions 0 or with
1e-300 <= |x| <= 1e300); values whose enclosure straddles a range bound are kept (superset).

  full(S)            every value of V+(S), exactly, for small S (deduplicated by exact value)
  query(S, lo, hi)   every value of V+(S) that may lie in [lo, hi] (a superset of V+(S) n [lo, hi])

Digit sets are tuples of characters; in the in-order mode (ORD=True) they are consecutive digits in
increasing order.  Windows are pairs of Arb balls (lo, hi); "may lie" means the value's ball is not
certainly below lo nor certainly above hi.
"""
import itertools, bisect
from fractions import Fraction as Fr
from flint import arb, ctx
import exact as ex
from exact import X, ZERO, ONE, Unresolved

PREC = 256
ctx.prec = PREC
BIG = arb(10) ** 300          # in-range bound
FULLMAX = 3                   # sets up to this size are enumerated completely
COMPACT = 4                   # fully enumerated sets of this size or more are stored compactly


def pb(e):
    """a base written with a leading unary minus needs parentheses: -x^y means -(x^y)"""
    return "(" + e + ")" if e.startswith("-") else e


def mn(a, b):
    return a if a < b else b


def mx(a, b):
    return a if a > b else b


class B:
    def __init__(self, ordered=False, fullmax=FULLMAX):
        self.ORD = ordered
        self.fullmax = fullmax
        self._full = {}
        self._sorted = {}
        self.unresolved = []

    # ------------------------------------------------------------ literals
    def literals(self, S):
        perms = [S] if self.ORD else sorted(set(itertools.permutations(S)))
        out = {}
        for perm in perms:
            t = "".join(perm)
            for i in range(len(t) + 1):
                if i == len(t):
                    v, s = X.rat(Fr(int(t))), t
                else:
                    v, s = X.rat(Fr(int(t), 10 ** (len(t) - i))), t[:i] + "." + t[i:]
                if ex.in_range(v):
                    out.setdefault(v.key(), (v, s))
        return out

    # ------------------------------------------------------------ splits
    def splits(self, S, smaller_first=True):
        """(S1, S2, s_is_left): S1 the part that plays 's'; for queries S1 is the smaller part"""
        n = len(S)
        out = []
        if self.ORD:
            for k in range(1, n):
                L, R = S[:k], S[k:]
                if smaller_first:
                    if len(L) <= len(R):
                        out.append((L, R, True))
                    else:
                        out.append((R, L, False))
                else:
                    out.append((L, R, True))
        else:
            for r in range(1, n // 2 + 1):
                for S1 in itertools.combinations(S, r):
                    S2 = tuple(d for d in S if d not in S1)
                    if 2 * r == n and S1 > S2:
                        continue
                    out.append((S1, S2, None))
        return out

    # ------------------------------------------------------------ forms of Lemma 3
    def forms(self, s, se, l, le, side):
        """values |op(+-s, +-l)| with s from the set S1; side True: S1 left, False: right, None: any"""
        out = []
        def put(v, e):
            if v is not None and ex.in_range(v):
                out.append((v, e))
        L = side in (True, None)
        R = side in (False, None)
        if side is False:
            put(ex.add(s, l), f"({le}+{se})")
        else:
            put(ex.add(s, l), f"({se}+{le})")
        d = ex.sub(s, l)
        sg = ex.sign(d)
        if side is False:
            put(ex.absval(d), f"({le}-{se})" if sg <= 0 else f"-({le}-{se})")
        elif side is True:
            put(ex.absval(d), f"({se}-{le})" if sg >= 0 else f"-({se}-{le})")
        else:
            put(ex.absval(d), f"({se}-{le})" if sg >= 0 else f"({le}-{se})")
        put(ex.mul(s, l), f"({le}*{se})" if side is False else f"({se}*{le})")
        if L:
            put(ex.div(s, l), f"({se}/{le})")
            put(ex.power(s, l), f"({pb(se)}^{le})")
            if not l.is_zero():
                put(ex.power(s, ex.neg(l)), f"({pb(se)}^-{le})")
        if R:
            put(ex.div(l, s), f"({le}/{se})")
            put(ex.power(l, s), f"({pb(le)}^{se})")
            if not s.is_zero():
                put(ex.power(l, ex.neg(s)), f"({pb(le)}^-{se})")
        return out

    # ------------------------------------------------------------ complete enumeration
    def full(self, S):
        S = tuple(S)
        if S in self._full:
            return self._full[S]
        out = dict(self.literals(S))
        if len(S) > 1:
            for S1, S2, side in self.splits(S, smaller_first=False):
                A, Bv = self.full(S1), self.full(S2)
                for (s, se) in A.values():
                    for (l, le) in Bv.values():
                        try:
                            for v, e in self.forms(s, se, l, le, True if self.ORD else None):
                                out.setdefault(v.key(), (v, e))
                        except Unresolved as u:
                            self.unresolved.append((S, se, le, str(u)))
        self._full[S] = out
        return out

    def sorted_full(self, S):
        S = tuple(S)
        if S not in self._sorted:
            items = []
            for (v, e) in self.full(S).values():
                b = ex.tight_ball(v)
                items.append((float(b.mid()), b, v, e))
            items.sort(key=lambda t: t[0])
            if len(S) >= COMPACT:
                # compact form for large sets: keep only key, ball and expression; the exact value is rebuilt
                # from the expression when the value is actually used (filter_full checks the two agree)
                wide = [t for t in items if t[1].rad() > abs(t[1].mid()) * 1e-12]
                setattr(self, "_wide_" + "".join(S), wide)
                items = [(k, b, None, e) for (k, b, v, e) in items]
                self._full.pop(S, None)
            self._sorted[S] = ([t[0] for t in items], items)
        return self._sorted[S]

    @staticmethod
    def rebuild(b, e):
        import exprparse
        v = exprparse.to_x(e)
        if ex.sign(v) < 0:
            v = ex.neg(v)
        w = ex.tight_ball(v)
        if not w.overlaps(b):
            raise AssertionError("compact value does not match its expression: " + e)
        return v

    # ------------------------------------------------------------ windows
    @staticmethod
    def may_lie(b, lo, hi):
        return not (b < lo) and not (b > hi)

    def filter_full(self, S, lo, hi):
        keys, items = self.sorted_full(S)
        flo = float(lo.lower()) if lo is not None else -1.0
        fhi = float(hi.upper()) if hi is not None else float("inf")
        i = bisect.bisect_left(keys, flo - abs(flo) * 1e-9 - 1e-300)
        out = {}
        while i < len(keys) and keys[i] <= fhi + abs(fhi) * 1e-9 + 1e-300:
            _, b, v, e = items[i]
            if v is None:
                if (lo is not None and b < lo) or (hi is not None and b > hi):
                    i += 1
                    continue
                v = self.rebuild(b, e)
            if ex.may_lie(v, lo, hi):
                out[v.key()] = (v, e)
            i += 1
        # values whose float midpoint is far from their ball (wide balls) are rechecked directly
        for (_, b, v, e) in self.wide(S):
            if ex.may_lie(v, lo, hi):
                out[v.key()] = (v, e)
        return out

    def wide(self, S):
        keys, items = self.sorted_full(S)
        attr = "_wide_" + "".join(S)
        if not hasattr(self, attr):
            setattr(self, attr, [t for t in items if t[1].rad() > abs(t[1].mid()) * 1e-12])
        return getattr(self, attr)

    def windows(self, fam, s, lo, hi):
        """list of (l_lo, l_hi) Arb windows for l in V+(S2), or ('const', value) entries"""
        zero = arb(0)
        out = []
        sb = ex.tight_ball(s)
        def W(a, b):
            a = a if a > zero else zero
            b = b if b < BIG else BIG
            if not (a > b):
                out.append((a, b))
        if fam == "add":
            W(lo - sb, hi - sb)
        elif fam == "sub":
            W(sb + lo, sb + hi)
            W(sb - hi, sb - lo)
        elif fam == "mul":
            if s.is_zero():
                out.append(("const", ZERO))
            else:
                W(lo / sb, hi / sb)
        elif fam == "sl":                                  # s / l
            if s.is_zero():
                out.append(("const", ZERO))
            else:
                W(sb / hi, sb / lo if lo > zero else BIG)
        elif fam == "ls":                                  # l / s
            if not s.is_zero():
                W(lo * sb, hi * sb)
        elif fam == "ps":                                  # s^(+-l)
            if s.is_zero():
                out.append(("const", ZERO))
            elif s.key() == ONE.key():
                out.append(("const", ONE))
            else:
                lam = sb.log()
                A = lo.log() if lo > zero else -BIG
                Bh = hi.log()
                t1, t2 = A / lam, Bh / lam
                tlo = mn(t1.lower(), t2.lower())
                thi = mx(t1.upper(), t2.upper())
                if tlo > zero:                             # l = |t|, t in [tlo, thi]
                    W(tlo, thi)
                elif thi < zero:
                    W(-thi, -tlo)
                else:
                    W(zero, mx(-tlo, thi))
        elif fam == "pl":                                  # l^(+-s)
            if s.is_zero():
                out.append(("const", ONE))
            else:
                inv = 1 / sb
                a1 = (lo.log() * inv).exp() if lo > zero else zero
                b1 = (hi.log() * inv).exp()
                W(a1, b1)
                a2 = (hi.log() * (-inv)).exp()
                b2 = (lo.log() * (-inv)).exp() if lo > zero else BIG
                W(a2, b2)
        return out

    def families(self, side):
        if side is True:
            return ("add", "sub", "mul", "sl", "ps")
        if side is False:
            return ("add", "sub", "mul", "ls", "pl")
        return ("add", "sub", "mul", "sl", "ls", "ps", "pl")

    def family_values(self, fam, s, se, l, le, side):
        """the values of family fam for the pair (s, l), with expressions"""
        sel = []
        if fam == "add":
            v = ex.add(s, l)
            sel.append((v, f"({le}+{se})" if side is False else f"({se}+{le})"))
        elif fam == "sub":
            d = ex.sub(s, l)
            sg = ex.sign(d)
            if side is False:
                sel.append((ex.absval(d), f"({le}-{se})" if sg <= 0 else f"-({le}-{se})"))
            elif side is True:
                sel.append((ex.absval(d), f"({se}-{le})" if sg >= 0 else f"-({se}-{le})"))
            else:
                sel.append((ex.absval(d), f"({se}-{le})" if sg >= 0 else f"({le}-{se})"))
        elif fam == "mul":
            sel.append((ex.mul(s, l), f"({le}*{se})" if side is False else f"({se}*{le})"))
        elif fam == "sl":
            sel.append((ex.div(s, l), f"({se}/{le})"))
        elif fam == "ls":
            sel.append((ex.div(l, s), f"({le}/{se})"))
        elif fam == "ps":
            sel.append((ex.power(s, l), f"({pb(se)}^{le})"))
            if not l.is_zero():
                sel.append((ex.power(s, ex.neg(l)), f"({pb(se)}^-{le})"))
        elif fam == "pl":
            sel.append((ex.power(l, s), f"({pb(le)}^{se})"))
            if not s.is_zero():
                sel.append((ex.power(l, ex.neg(s)), f"({pb(le)}^-{se})"))
        return [(v, e) for v, e in sel if v is not None and ex.in_range(v)]

    # ------------------------------------------------------------ the window query
    def query(self, S, lo, hi):
        S = tuple(S)
        if len(S) <= self.fullmax:
            return self.filter_full(S, lo, hi)
        out = {}
        for k, (v, e) in self.literals(S).items():
            if ex.may_lie(v, lo, hi):
                out[k] = (v, e)
        for S1, S2, side in self.splits(S):
            witness = None
            for (s, se) in self.full(S1).values():
                for fam in self.families(side):
                    try:
                        wins = self.windows(fam, s, lo, hi)
                    except Unresolved as u:
                        self.unresolved.append((S, se, fam, str(u)))
                        continue
                    for w in wins:
                        if w[0] == "const":
                            val = w[1]
                            if ex.may_lie(val, lo, hi):
                                if witness is None:
                                    witness = self.positive_literal(S2)
                                for v, e in self.family_values(fam, s, se, witness[0], witness[1], side):
                                    if v.key() == val.key():
                                        out.setdefault(v.key(), (v, e))
                            continue
                        sub = self.query(S2, w[0], w[1])
                        for (l, le) in sub.values():
                            try:
                                for v, e in self.family_values(fam, s, se, l, le, side):
                                    if ex.may_lie(v, lo, hi):
                                        out.setdefault(v.key(), (v, e))
                            except Unresolved as u:
                                self.unresolved.append((S, se, le, str(u)))
        return out

    def positive_literal(self, S):
        t = "".join(S)
        return (X.rat(Fr(int(t))), t)
