"""B2b: among the identities of F (1-9) with ln ln N >= 50, which ones admit a writing in which
X or Y is a three-story tower d1^(d2^d3) (single-digit blocks)?  Exact comparison of values."""
import itertools, e_F
from e_F import V, pw, search
hits, unres, skipped, SB = search("123456789", 50.0)
# values of d1^(d2^d3) on 3-digit sets
tower = {}
for S in itertools.combinations("23456789", 3):
    for perm in itertools.permutations(S):
        for (a, sa) in V((perm[0],)).values():
            for (b, sb) in V((perm[1],)).values():
                for (c, sc) in V((perm[2],)).values():
                    inner = pw(b, c)
                    if inner is None: continue
                    v = pw(a, inner)
                    if v is None or v[0] == "T": continue
                    for w in (v, e_F.neg(v)):
                        tower.setdefault((tuple(sorted(S)), w), f"{sa}^({sb}^{sc})")
def digits(s): return tuple(sorted(ch for ch in s if ch.isdigit()))
import re
res = {}
for h in hits:
    lm, su, sa, sx, sc, sy = h[:6]
    # recover the block values from e_F by expression string
    def val(expr, D):
        for (v, s) in V(D).values():
            if s == expr: return v
        return None
    X = val(sx, digits(sx)); Y = val(sy, digits(sy))
    tx = tower.get((digits(sx), X)); ty = tower.get((digits(sy), Y))
    key = lm.mid().str(12, radius=False)
    res.setdefault(key, set())
    if tx: res[key].add("X=" + tx)
    if ty: res[key].add("Y=" + ty)
for k in sorted(res, key=float, reverse=True): print("ln ln N =", k, "three-story writings:", sorted(res[k])[:4] or "none")
