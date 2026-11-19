"""
Motor B: verification of motor A's candidate paths (Theorem C).

A path is a sequence of steps  "mask:fam:side:o:v:code"  or  "mask:fam:side:l:llo:lhi"  (the set S1,
the family, the side of S1, and the class of T(S1)) followed by a leaf: "mask:leaf:o:v:code",
"mask:leaf:l:llo:lhi", "mask:lit:o:v:code" or "mask:const:value".

verify_path(B, digits, path, lo, hi) returns every value of V+(digits) covered by the path and lying
(possibly) in [lo, hi], computed exactly by motor B, as a dict key -> (value, expression).
"""
import math
from flint import arb
import motorb, exact as ex
from exact import Unresolved
from resolve import hull_of

SIDE = {"1": True, "0": False, "-1": None}


def mask_digits(mask, n):
    return tuple(str(b + 1) for b in range(n) if mask >> b & 1)


def class_window(spec):
    """spec: ('o', v, code) | ('l', m1, m2) | ('a', 'p/q', m1, m2)"""
    from resolve import window_of
    return window_of(":".join(spec))


def isect(w1, w2):
    lo = w1[0] if w1[0] > w2[0] else w2[0]
    hi = w1[1] if w1[1] < w2[1] else w2[1]
    return lo, hi


class Verifier:
    def __init__(self, n, ordered, fullmax=motorb.FULLMAX):
        self.n = n
        self.B = motorb.B(ordered=bool(ordered), fullmax=fullmax)
        self.cache = {}

    def class_values(self, digs, spec):
        key = (digs, spec)
        if key not in self.cache:
            if len(self.cache) > 20000:      # bounded memory on long runs
                self.cache.clear()
            self.cache[key] = self.B.query(digs, *class_window(spec))
        return self.cache[key]

    def leaf(self, digs, parts, lo, hi):
        mask, what = parts[0], parts[1]
        if what == "const":
            raise ValueError("const leaf without a step")
        cw = class_window(tuple(parts[2:]))
        out = {}
        if what == "lit":
            for k, (v, e) in self.B.literals(digs).items():
                if ex.may_lie(v, *cw) and ex.may_lie(v, lo, hi):
                    out[k] = (v, e)
            return out
        w = isect(cw, (lo, hi))
        if w[0] > w[1]:
            return out
        for k, (v, e) in self.B.query(digs, *w).items():
            out[k] = (v, e)
        return out

    def verify(self, digs, steps, leaf, lo, hi):
        if not steps:
            return self.leaf(digs, leaf, lo, hi)
        mask, fam, side = steps[0][:3]
        spec = tuple(steps[0][3:])
        side = SIDE[side]
        S1 = mask_digits(int(mask), self.n)
        S2 = tuple(d for d in digs if d not in S1)
        out = {}
        for (s, se) in self.class_values(S1, spec).values():
            if len(steps) == 1 and leaf[1] == "const":
                w = self.B.positive_literal(S2)
                for v, e in self.B.family_values(fam, s, se, w[0], w[1], side):
                    if ex.may_lie(v, lo, hi):
                        out.setdefault(v.key(), (v, e))
                continue
            for win in self.B.windows(fam, s, lo, hi):
                if win[0] == "const":
                    continue
                sub = self.verify(S2, steps[1:], leaf, win[0], win[1])
                for (l, le) in sub.values():
                    for v, e in self.B.family_values(fam, s, se, l, le, side):
                        if ex.may_lie(v, lo, hi):
                            out.setdefault(v.key(), (v, e))
        return out

    def verify_line(self, line, lo, hi):
        toks = line.split()
        if toks and toks[0] == "C":
            toks = toks[1:]
        steps = [tuple(t.split(":")) for t in toks[:-1]]
        leaf = toks[-1].split(":")
        digs = tuple(str(d) for d in range(1, self.n + 1))
        return self.verify(digs, steps, leaf, lo, hi)
