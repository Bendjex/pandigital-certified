"""
Motor B: resolution of motor A's singular derivations (Lemma 7).

Input: lines "KEY | kind form" written by `motora stage`, with
    KEY = "S form m1 op1 m2 op2",  op = "o:<v>:<code>" (ordinary class) or "l:<llo>:<lhi>" (L-class)
Output (appended to the resolution file): "KEY | item item ..." with items
    q:p/q        an exact rational value (|p|, q <= 2^53)
    x:lo:hi      an enclosure (outward doubles) not containing 0 or 1
    l:llo:lhi    an enclosure of ln x, for values with |ln x| < 2^-10 (x != 1)
covering f(s, l) for every s in V+(S1) within op1's enclosure and l in V+(S2) within op2's enclosure,
with f the form, whenever f(s, l) is defined and in range.  Values are computed exactly (motor B).
Anything that cannot be decided is reported as UNRESOLVED and makes the run fail.
"""
import sys, math
from fractions import Fraction as Fr
from flint import arb, ctx
import motorb, exact as ex
from exact import X, Unresolved

ctx.prec = 256
NEAR1 = 2.0 ** -10
RLIM = 2 ** 53
FORMS = ["x+y", "|x-y|", "x*y", "x/y", "y/x", "x^y", "x^-y", "y^x", "y^-x"]


def digits(mask, n):
    return tuple(str(b + 1) for b in range(n) if mask >> b & 1)


def hull_of(v, code):
    if code == 0:
        return v, v
    e = math.frexp(v)[1]
    rho = math.ldexp(1.0, e - 53 + code - 1)
    a = math.nextafter(v - rho, -math.inf)
    b = math.nextafter(v + rho, math.inf)
    return max(a, 0.0), b


def window_of(op):
    """Arb window [lo, hi] for an operand description: o:v:code, l:m1:m2 (e^m), a:p/q:m1:m2 (q e^m)"""
    parts = op.split(":")
    kind = parts[0]
    if kind == "o":
        lo, hi = hull_of(float.fromhex(parts[1]), int(parts[2]))
        return arb(lo), arb(hi)
    if kind == "l":
        q, m1, m2 = Fr(1), parts[1], parts[2]
    else:
        q, m1, m2 = Fr(parts[1]), parts[2], parts[3]
    Q = arb(q.numerator) / q.denominator
    lo = Q * arb(float.fromhex(m1)).exp()
    hi = Q * arb(float.fromhex(m2)).exp()
    return arb(lo.lower()), arb(hi.upper())


def apply_form(f, x, y):
    if f == 0: return ex.add(x, y)
    if f == 1: return ex.absval(ex.sub(x, y))
    if f == 2: return ex.mul(x, y)
    if f == 3: return ex.div(x, y)
    if f == 4: return ex.div(y, x)
    if f == 5: return ex.power(x, y)
    if f == 6: return None if y.is_zero() else ex.power(x, ex.neg(y))
    if f == 7: return ex.power(y, x)
    if f == 8: return None if x.is_zero() else ex.power(y, ex.neg(x))


def out_down(b):
    return math.nextafter(float(b), -math.inf)


def out_up(b):
    return math.nextafter(float(b), math.inf)


def item_of(v):
    """exact rational, anchored (rational q times e^mu, mu tiny), near 1, or an enclosure"""
    if v.is_zero():
        return "q:0/1"
    if v.is_rational():
        q = abs(v.rational())
        if abs(q.numerator) <= RLIM and q.denominator <= RLIM:
            return f"q:{q.numerator}/{q.denominator}"
    anchors = []
    if v.kind == "S" and () in v.terms:
        Q = abs(v.terms[()])
        if Q.numerator <= RLIM and Q.denominator <= RLIM:
            anchors.append(Q)
    for prec in (256, 1024, 4096, 16384):
        old = ctx.prec
        ctx.prec = prec
        try:
            b = v.ball(prec)
            if b.contains(0):
                continue
            for Q in anchors:                      # v = Q e^mu with mu tiny
                mu = (abs(b) / (arb(Q.numerator) / Q.denominator)).log()
                if mu.contains(0):
                    continue
                if abs(mu) < NEAR1:
                    if Q == 1:
                        return f"l:{out_down(mu.lower()).hex()}:{out_up(mu.upper()).hex()}"
                    return f"a:{Q.numerator}/{Q.denominator}:{out_down(mu.lower()).hex()}:{out_up(mu.upper()).hex()}"
            L = abs(b).log()
            if L.contains(0):
                continue
            if abs(L) < NEAR1:
                return f"l:{out_down(L.lower()).hex()}:{out_up(L.upper()).hex()}"
            return f"x:{out_down(b.lower()).hex()}:{out_up(b.upper()).hex()}"
        finally:
            ctx.prec = old
    raise Unresolved("item not separated from 0 or 1")


_qcache = {}


def cached_query(B, digs, op):
    key = (digs, op)
    if key not in _qcache:
        _qcache[key] = B.query(digs, *window_of(op))
    return _qcache[key]


def main():
    n, ordered, pend, resfile = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4]
    B = motorb.B(ordered=bool(ordered))
    lines = [l.rstrip("\n") for l in open(pend) if " | " in l]
    out = open(resfile, "a")
    nres = nunres = 0
    for line in lines:
        key = line.split(" | ")[0]
        S, form, m1, op1, m2, op2 = key.split(" ")
        form, m1, m2 = int(form), int(m1), int(m2)
        try:
            X1 = cached_query(B, digits(m1, n), op1)
            X2 = cached_query(B, digits(m2, n), op2)
            items = set()
            for (x, xe) in X1.values():
                for (y, ye) in X2.values():
                    v = apply_form(form, x, y)
                    if v is None or not ex.in_range(v):
                        continue
                    items.add(item_of(v))
            out.write(key + " | " + " ".join(sorted(items)) + "\n")
            nres += 1
        except Unresolved as u:
            print("UNRESOLVED", key, u)
            nunres += 1
    out.close()
    if B.unresolved:
        print("UNRESOLVED inside motor B:", B.unresolved[:5])
        nunres += len(B.unresolved)
    print(f"resolved {nres}, unresolved {nunres}")
    sys.exit(1 if nunres else 0)


if __name__ == "__main__":
    main()
