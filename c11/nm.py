"""
C11: the variant WITHOUT unary minus.

Grammar: literal | (E1 op E2), op in {+, -, *, /, ^}; same values and validity rules as C1
(x^y needs x > 0, or x = 0 and y > 0, or x < 0 and y an integer), same range rule as C3.
Without unary minus the value set is no longer symmetric, so every magnitude u in V+(S) (motor B)
carries two flags:
    POS: +u is the value of some minus-free expression,  NEG: -u is.
(0 gets both flags.)  The flags are propagated through exactly the same enumeration as motor B
(complete tables, window queries, path verification), OR-ed over every derivation of a value.
Because every minus-free expression is also an expression of C1, its magnitude is found by the
C3 machinery (motor A paths + motor B verification); the flags only decide which of the values
found are reachable with the right sign without unary minus.

Entry of a value dict: key -> [v, expr, flags, expr_pos, expr_neg]   (expr: a C1 writing of |v|)
"""
import sys, os
SRC = os.environ.get("PANDI_SRC", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
sys.path.insert(0, SRC)
import motorb, verify, exact as ex
from exact import Unresolved, ONE, ZERO
from motorb import pb

POS, NEG = 1, 2


def signs(fl):
    out = []
    if fl & POS:
        out.append(1)
    if fl & NEG:
        out.append(-1)
    return out


def is_int(x):
    if not x.is_rational():
        return None
    q = x.rational()
    return q.numerator if q.denominator == 1 else None


def par(e):
    return "(" + e + ")"


def signed_family(fam, s, sfl, sp, sn, l, lfl, lp, ln, side):
    """minus-free signed results of family fam for the pair (s, l).
    s, l: magnitudes; sfl/lfl: flags; sp/sn, lp/ln: minus-free expressions of +s/-s, +l/-l.
    side: True S1 is the left operand, False S1 is the right operand, None any order.
    Returns list of (magnitude, flagbit, minus-free expression)."""
    out = []
    S = [(sg, sp if sg > 0 else sn) for sg in signs(sfl)]
    L = [(tg, lp if tg > 0 else ln) for tg in signs(lfl)]
    orient = [True, False] if side is None else [side]
    def emit(mag, sgn, e):
        if mag is None or not ex.in_range(mag):
            return
        if mag.is_zero():
            out.append((mag, POS | NEG, e))
        else:
            out.append((mag, POS if sgn > 0 else NEG, e))
    for (sg, se) in S:
        for (tg, le) in L:
            for s_left in orient:
                A, ea, B, eb = (s, se, l, le) if s_left else (l, le, s, se)
                ag, bg = (sg, tg) if s_left else (tg, sg)
                if fam in ("add", "sub"):
                    for op in "+-":
                        eps = bg if op == "+" else -bg           # value = ag*A + eps*B
                        e = f"({ea}{op}{eb})"
                        if fam == "add" and ag == eps:
                            emit(ex.add(s, l), ag, e)
                        elif fam == "sub" and ag != eps:
                            d = ex.sub(s, l)                     # value = ag*(A - B)
                            sd = ex.sign(d)
                            sgn = ag * (sd if s_left else -sd)
                            emit(ex.absval(d), sgn, e)
                elif fam == "mul":
                    emit(ex.mul(s, l), ag * bg, f"({ea}*{eb})")
                elif fam in ("sl", "ls"):
                    # sl: magnitude s/l, needs s as numerator (left);  ls: l/s, needs l left
                    if (fam == "sl") == s_left:
                        mag = ex.div(s, l) if fam == "sl" else ex.div(l, s)
                        emit(mag, ag * bg, f"({ea}/{eb})")
                elif fam in ("ps", "pl"):
                    # ps: base s (left), exponent +-l;  pl: base l (left), exponent +-s
                    if (fam == "ps") == s_left:
                        base, bsg, expo, esg = (s, sg, l, tg) if fam == "ps" else (l, tg, s, sg)
                        mag = ex.power(base, expo if esg > 0 else ex.neg(expo)) if not (esg < 0 and expo.is_zero()) else None
                        if esg < 0 and expo.is_zero():
                            continue                             # x^-0 is the same as x^0 (other sign)
                        if mag is None:
                            continue
                        if bsg > 0 or base.is_zero():
                            emit(mag, 1, f"({pb(ea)}^{eb})")
                        else:
                            k = is_int(expo)
                            if k is None:
                                continue                         # negative base, non-integer exponent
                            emit(mag, 1 if k % 2 == 0 else -1, f"({ea}^{eb})")
    return out


class NB(motorb.B):
    """motor B with sign flags (minus-free reachability)"""
    WARN = []

    def const_partners(self, S2, sfl):
        """partners l for a constant family (s = 0 or 1, value independent of l): complete table
        for small S2; otherwise the literals, which suffice unless s = 1 is reachable only as -1
        (then the parity of an integer l matters): such cases are recorded in WARN"""
        if len(S2) <= self.fullmax:
            return list(self.full(S2).values())
        if not (sfl & POS):
            self.WARN.append(S2)
        return list(self.literals(S2).values())

    def literals(self, S):
        return {k: [v, e, POS, e, None] for k, (v, e) in super().literals(S).items()}

    def merge(self, out, k, v, e, fl, epos, eneg):
        cur = out.get(k)
        if cur is None:
            out[k] = [v, e, fl, epos if fl & POS else None, eneg if fl & NEG else None]
        else:
            if fl & POS and cur[3] is None:
                cur[3] = epos
            if fl & NEG and cur[4] is None:
                cur[4] = eneg
            cur[2] |= fl

    def pair_into(self, out, fams, s, se, sfl, sp, sn, l, le, lfl, lp, ln, side, lo=None, hi=None):
        for fam in fams:
            # with-minus writing (for reference) and the magnitudes come from motor B itself
            base = {}
            for v, e in self.family_values(fam, s, se, l, le, side):
                if lo is None or ex.may_lie(v, lo, hi):
                    base.setdefault(v.key(), (v, e))
            for k, (v, e) in base.items():
                self.merge(out, k, v, e, 0, None, None)
            if lo is not None and not base:
                # no magnitude of this family in the window: the signed results are among motor B's
                # magnitudes (asserted below whenever they are computed), so none of them is either
                continue
            for mag, bit, e in signed_family(fam, s, sfl, sp, sn, l, lfl, lp, ln, side):
                k = mag.key()
                if k in base:
                    self.merge(out, k, mag, base[k][1], bit, e, e)
                elif lo is None or ex.may_lie(mag, lo, hi):
                    # the magnitude must be one of motor B's values (same arithmetic calls)
                    raise AssertionError("signed result not among motor B values: " + e)

    def full(self, S):
        S = tuple(S)
        if S in self._full:
            return self._full[S]
        out = dict(self.literals(S))
        if len(S) > 1:
            for S1, S2, side in self.splits(S, smaller_first=False):
                A, Bv = self.full(S1), self.full(S2)
                sd = True if self.ORD else None
                for (s, se, sfl, sp, sn) in A.values():
                    for (l, le, lfl, lp, ln) in Bv.values():
                        try:
                            self.pair_into(out, self.families(sd),
                                           s, se, sfl, sp, sn, l, le, lfl, lp, ln, sd)
                        except Unresolved as u:
                            self.unresolved.append((S, se, le, str(u)))
        self._full[S] = out
        return out

    def sorted_full(self, S):
        S = tuple(S)
        if S not in self._sorted:
            items = []
            for ent in self.full(S).values():
                b = ex.tight_ball(ent[0])
                items.append((float(b.mid()), b, ent))
            items.sort(key=lambda t: t[0])
            wide = [t for t in items if t[1].rad() > abs(t[1].mid()) * 1e-12]
            self._sorted[S] = ([t[0] for t in items], items, wide)
        return self._sorted[S]

    def filter_full(self, S, lo, hi):
        import bisect
        keys, items, wide = self.sorted_full(S)
        flo = float(lo.lower()) if lo is not None else -1.0
        fhi = float(hi.upper()) if hi is not None else float("inf")
        i = bisect.bisect_left(keys, flo - abs(flo) * 1e-9 - 1e-300)
        out = {}
        while i < len(keys) and keys[i] <= fhi + abs(fhi) * 1e-9 + 1e-300:
            ent = items[i][2]
            if ex.may_lie(ent[0], lo, hi):
                out[ent[0].key()] = list(ent)
            i += 1
        for (_, b, ent) in wide:                      # wide balls rechecked directly
            if ex.may_lie(ent[0], lo, hi):
                out[ent[0].key()] = list(ent)
        return out

    def query(self, S, lo, hi):
        S = tuple(S)
        if len(S) <= self.fullmax:
            return self.filter_full(S, lo, hi)
        out = {}
        for k, ent in self.literals(S).items():
            if ex.may_lie(ent[0], lo, hi):
                out[k] = ent
        for S1, S2, side in self.splits(S):
            for (s, se, sfl, sp, sn) in self.full(S1).values():
                for fam in self.families(side):
                    try:
                        wins = self.windows(fam, s, lo, hi)
                    except Unresolved as u:
                        self.unresolved.append((S, se, fam, str(u)))
                        continue
                    for w in wins:
                        if w[0] == "const":
                            # s = 0 or s = 1: every l of S2 gives the same magnitude
                            for (l, le, lfl, lp, ln) in self.const_partners(S2, sfl):
                                self.pair_into(out, (fam,), s, se, sfl, sp, sn, l, le, lfl, lp, ln, side, lo, hi)
                            continue
                        sub = self.query(S2, w[0], w[1])
                        for (l, le, lfl, lp, ln) in sub.values():
                            try:
                                self.pair_into(out, (fam,), s, se, sfl, sp, sn, l, le, lfl, lp, ln, side, lo, hi)
                            except Unresolved as u:
                                self.unresolved.append((S, se, le, str(u)))
        return out


class NVerifier(verify.Verifier):
    def __init__(self, n, ordered, fullmax=motorb.FULLMAX):
        self.n = n
        self.B = NB(ordered=bool(ordered), fullmax=fullmax)
        self.cache = {}

    def leaf(self, digs, parts, lo, hi):
        what = parts[1]
        cw = verify.class_window(tuple(parts[2:]))
        out = {}
        if what == "lit":
            for k, ent in self.B.literals(digs).items():
                if ex.may_lie(ent[0], *cw) and ex.may_lie(ent[0], lo, hi):
                    out[k] = ent
            return out
        w = verify.isect(cw, (lo, hi))
        if w[0] > w[1]:
            return out
        return dict(self.B.query(digs, *w))

    def verify(self, digs, steps, leaf, lo, hi):
        if not steps:
            return self.leaf(digs, leaf, lo, hi)
        mask, fam, side = steps[0][:3]
        spec = tuple(steps[0][3:])
        side = verify.SIDE[side]
        S1 = verify.mask_digits(int(mask), self.n)
        S2 = tuple(d for d in digs if d not in S1)
        out = {}
        for (s, se, sfl, sp, sn) in self.class_values(S1, spec).values():
            if len(steps) == 1 and leaf[1] == "const":
                for (l, le, lfl, lp, ln) in self.B.const_partners(S2, sfl):
                    self.B.pair_into(out, (fam,), s, se, sfl, sp, sn, l, le, lfl, lp, ln, side, lo, hi)
                continue
            for win in self.B.windows(fam, s, lo, hi):
                if win[0] == "const":
                    continue
                sub = self.verify(S2, steps[1:], leaf, win[0], win[1])
                for (l, le, lfl, lp, ln) in sub.values():
                    self.B.pair_into(out, (fam,), s, se, sfl, sp, sn, l, le, lfl, lp, ln, side, lo, hi)
        return out
