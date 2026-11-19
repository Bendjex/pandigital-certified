"""
Upper bound for the minus-free optimum of a cell in any order (C11), by a meet-in-the-middle search in
double precision over expressions whose top operation splits the digits 4+4 or 3+5 (or k+(n-k) for the
sizes given). Rules as C1 without unary minus: literals are the digits of a set in any order with at
most one decimal point; operations + - * / ^; negative base only with an integer exponent; every
intermediate value 0 or of magnitude in [1e-300, 1e300].

The search is heuristic (doubles); only its output is used, and each output formula is re-evaluated
independently at 80 digits (certify.independent_eval: digits, validity, range), so the bound is rigorous.
usage: mitm.py n const[,const...] [splits e.g. 4,3] [keep]
"""
import os
import sys, itertools, json, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import certify, mpmath
import c11_certify as cc

LO, HI = 1e-300, 1e300
OPS = ["+", "-", "r-", "*", "/", "r/", "^", "r^"]          # r: operands swapped (b op a)
np.seterr(all="ignore")


def literals(S):
    out = []
    for p in itertools.permutations(S):
        t = "".join(p)
        out.append(t)
        for i in range(len(t)):
            out.append(t[:i] + "." + t[i:])
    return sorted(set(out))


def apply(op, a, b):
    """op on arrays a, b (broadcast); returns values with NaN where invalid or out of range"""
    if op == "+": r = a + b
    elif op == "-": r = a - b
    elif op == "r-": r = b - a
    elif op == "*": r = a * b
    elif op == "/": r = np.where(b != 0, a / np.where(b != 0, b, 1), np.nan)
    elif op == "r/": r = np.where(a != 0, b / np.where(a != 0, a, 1), np.nan)
    else:
        x, y = (a, b) if op == "^" else (b, a)
        x = x + 0 * y; y = y + 0 * x
        yi = np.round(y)
        isint = (y == yi) & (np.abs(y) < 2 ** 52)
        with np.errstate(all="ignore"):
            pos = np.where(x > 0, np.power(np.abs(x), y), np.nan)
            neg = np.where((x < 0) & isint, np.power(np.abs(x), y) * np.where(np.mod(yi, 2) == 0, 1.0, -1.0), np.nan)
            zer = np.where((x == 0) & (y > 0), 0.0, np.nan)
        r = np.where(x > 0, pos, np.where(x < 0, neg, zer))
    ar = np.abs(r)
    ok = np.isfinite(r) & ((r == 0) | ((ar >= LO) & (ar <= HI)))
    return np.where(ok, r, np.nan)


class Table:
    """values of one digit set with provenance: (kind, ...) per value"""
    def __init__(self, vals, prov):
        self.v = vals; self.p = prov          # prov: list-like of tuples or arrays


TAB = {}


def table(S):
    S = tuple(S)
    if S in TAB:
        return TAB[S]
    lits = literals(S)
    vals = [np.array([float(l) for l in lits])]
    provs = [("lit", lits)]
    if len(S) > 1:
        for r in range(1, len(S) // 2 + 1):
            for S1 in itertools.combinations(S, r):
                S2 = tuple(d for d in S if d not in S1)
                if 2 * r == len(S) and S1 > S2:
                    continue
                A, B = table(S1), table(S2)
                a = A.v[:, None]; b = B.v[None, :]
                for k, op in enumerate(OPS):
                    rr = apply(op, a, b)
                    I, J = np.nonzero(np.isfinite(rr))
                    vals.append(rr[I, J])
                    provs.append(("op", op, S1, S2, I.astype(np.int32), J.astype(np.int32)))
    allv = np.concatenate(vals)
    # provenance index: segment and offset
    seg = np.concatenate([np.full(len(v), i, dtype=np.int32) for i, v in enumerate(vals)])
    off = np.concatenate([np.arange(len(v), dtype=np.int32) for v in vals])
    u, first = np.unique(allv, return_index=True)
    t = Table(u, (provs, seg[first], off[first]))
    TAB[S] = t
    return t


def expr(S, i):
    t = table(S)
    provs, seg, off = t.p
    pv = provs[seg[i]]
    if pv[0] == "lit":
        return pv[1][off[i]]
    _, op, S1, S2, I, J = pv
    ea, eb = expr(S1, I[off[i]]), expr(S2, J[off[i]])
    if op.startswith("r"):
        ea, eb, op = eb, ea, op[1:]
    return f"({ea}{op}{eb})"


def best_pairs(a, b, c, keep):
    """indices (i, j, form) of the pairs whose combination is closest to c (heuristic)"""
    order = np.argsort(b); bs = b[order]
    res = []
    for op in OPS:
        with np.errstate(all="ignore"):
            if op == "+": t = c - a
            elif op == "-": t = a - c
            elif op == "r-": t = c + a
            elif op == "*": t = np.where(a != 0, c / a, np.nan)
            elif op == "/": t = a / c
            elif op == "r/": t = c * a
            elif op == "^": t = np.where((a > 0) & (a != 1), math.log(c) / np.log(np.where(a > 0, a, 1)), np.nan)
            else: t = np.where(a != 0, np.power(c, 1.0 / np.where(a != 0, a, 1)), np.nan)
        ok = np.isfinite(t)
        ii = np.nonzero(ok)[0]
        k = np.searchsorted(bs, t[ii])
        for d in (-1, 0):
            kk = np.clip(k + d, 0, len(bs) - 1)
            jj = order[kk]
            r = apply(op, a[ii], b[jj])
            e = np.abs(r - c)
            e = np.where(np.isfinite(e), e, np.inf)
            if len(e) > keep:
                sel = np.argpartition(e, keep)[:keep]
            else:
                sel = np.arange(len(e))
            for s in sel:
                if np.isfinite(e[s]):
                    res.append((float(e[s]), int(ii[s]), int(jj[s]), op))
    res.sort()
    return res[:keep]


def search(n, consts, sizes, keep=8):
    D = tuple(str(d) for d in range(1, n + 1))
    cand = {c: [] for c in consts}
    for k in sizes:
        for S1 in itertools.combinations(D, k):
            S2 = tuple(d for d in D if d not in S1)
            if 2 * k == n and S1 > S2:
                continue
            A, B = table(S1), table(S2)
            for c in consts:
                cv = float(certify.mp_constant(c))
                for e, i, j, op in best_pairs(A.v, B.v, cv, keep):
                    ea, eb = expr(S1, i), expr(S2, j)
                    if op.startswith("r"):
                        ea, eb, o2 = eb, ea, op[1:]
                    else:
                        o2 = op
                    cand[c].append((e, f"({ea}{o2}{eb})"))
            # free the large tables of this split
            for S in (S1, S2):
                if len(S) >= 5:
                    TAB.pop(S, None)
    out = {}
    for c in consts:
        best = None
        for e, f in sorted(cand[c])[:60]:
            dok, rok, val, margin = certify.independent_eval(f, n, 0)
            if not (dok and rok) or cc.has_unary_minus(f):
                continue
            with mpmath.workdps(80):
                er = abs(val - certify.mp_constant(c))
            if best is None or er < best[0]:
                best = (er, f)
        out[c] = {"error": mpmath.nstr(best[0], 20), "formula": best[1]}
        print(c, n, out[c], flush=True)
    return out


if __name__ == "__main__":
    n = int(sys.argv[1]); consts = sys.argv[2].split(",")
    sizes = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [n // 2]
    keep = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    res = search(n, consts, sizes, keep)
    print(json.dumps({"n": n, "sizes": sizes, "bounds": res}))
