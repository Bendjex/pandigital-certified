"""brute-force check of the flags: all minus-free expressions for small digit sets, exact values"""
import os
import sys, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import nm, exact as ex, motorb
from fractions import Fraction as Fr
from exact import X

def lits(S, ordered):
    perms = [S] if ordered else set(itertools.permutations(S))
    out = []
    for p in perms:
        t = "".join(p)
        for i in range(len(t) + 1):
            v = X.rat(Fr(int(t))) if i == len(t) else X.rat(Fr(int(t), 10 ** (len(t) - i)))
            out.append(v)
    return out

def signed_vals(S, ordered, memo={}):
    """exact signed values (as (sign, magnitude key)) of minus-free expressions; in range only"""
    k = (S, ordered)
    if k in memo: return memo[k]
    out = {}
    def put(v):
        if v is None or not ex.in_range(ex.absval(v)): return
        out[(ex.sign(v), ex.absval(v).key())] = v
    for v in lits(S, ordered): put(v)
    n = len(S)
    splits = []
    if ordered:
        splits = [(S[:i], S[i:]) for i in range(1, n)]
    else:
        for r in range(1, n):
            for A in itertools.combinations(S, r):
                splits.append((A, tuple(d for d in S if d not in A)))
    for A, B in splits:
        for a in signed_vals(A, ordered).values():
            for b in signed_vals(B, ordered).values():
                put(ex.add(a, b)); put(ex.sub(a, b)); put(ex.mul(a, b))
                if not b.is_zero(): put(ex.div(a, b))
                sa = ex.sign(a)
                if sa > 0: put(ex.power(a, b))
                elif sa == 0:
                    if ex.sign(b) > 0: put(ex.ZERO)
                else:
                    k2 = nm.is_int(b)
                    if k2 is not None:
                        m = ex.power(ex.neg(a), b)
                        if m is not None: put(m if k2 % 2 == 0 else ex.neg(m))
    memo[k] = out
    return out

for ordered in (True, False):
    for S in [("1","2","3"), ("2","3","4"), ("1","2","3","4")] + ([("3","4","5","6")] if ordered else []):
        B = nm.NB(ordered=ordered, fullmax=4)
        full = B.full(S)
        ref = signed_vals(S, ordered)
        def nk(v):
            b = ex.tight_ball(v)
            return b.mid().str(40, radius=False)
        got = set()
        for key, ent in full.items():
            if ent[2] & nm.POS: got.add((1 if not ent[0].is_zero() else 0, nk(ent[0])))
            if ent[2] & nm.NEG and not ent[0].is_zero(): got.add((-1, nk(ent[0])))
        refs = set((sg, nk(ex.absval(v))) for (sg, _), v in ref.items())
        print("ord" if ordered else "any", "".join(S), "flags", len(got), "brute", len(refs), "missing", len(refs - got), "extra", len(got - refs))
