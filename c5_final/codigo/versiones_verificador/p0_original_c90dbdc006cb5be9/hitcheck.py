"""
C5, residual hits of motor X: each hit is reconstructed explicitly and evaluated independently.

1. The query of cells_q.py is rerun in DEBUG mode (provenance of every hit), with Ent.uniq patched so
   that each N-table entry keeps ALL the derivations merged into it (uniq() keeps one entry per
   identical enclosure; the computation itself is unchanged).
2. A hit expands into the set of all values it stands for: V classes into their exact members
   (motor B, motorx.members), N entries into all their derivations, query pairs into their two parts.
3. Every value is evaluated with lvl.py (multi-level intervals, mpmath, independent of motor X's
   double-precision enclosures) and must lie outside the window [c - E', c + E'], E' = E (1 + 1e-9).
usage: hitcheck.py n ord const[,const...] [out.json]
"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, qx, ivl, motorx, lvl
from cells_x import const
from flint import arb
from mpmath import iv
import exact as ex


def uniq_all(self):
    if not len(self):
        self.allprov = []
        return self
    _, idx, inv = np.unique(self.uniq_key(), axis=0, return_index=True, return_inverse=True)
    inv = np.asarray(inv).ravel()
    order = np.argsort(idx)
    pos = np.empty(len(idx), np.int64); pos[order] = np.arange(len(idx))
    out = self.sel(np.sort(idx))
    groups = [[] for _ in range(len(idx))]
    for k in range(len(self)):
        groups[pos[inv[k]]].append(tuple(int(t) for t in self.prov[k]))
    out.allprov = groups
    return out


motorx.Ent.uniq = uniq_all
motorx.ALLPROV[0] = True      # N tables with all their derivations (cached apart, files N H...)


def arb_iv(b):
    m, e = b.mid().man_exp()
    rm, re_ = b.rad().man_exp()
    mid = iv.mpf(int(m)) * iv.mpf(2) ** int(e)
    rad = (iv.mpf(int(rm)) * iv.mpf(2) ** int(re_)).b if int(rm) else 0
    return mid + iv.mpf([-rad, rad])


def leaf(x):
    x = ex.absval(x)
    if x.is_rational():
        return lvl.exact(x.rational())
    v = lvl.from_x(arb_iv(x.ball(lvl.PREC + 64)))
    # (03/11) an exact monomial of motor B, c * prod p^e (kind 'S', one term, e rational in (0, 1)): lvl also gets its
    # vector of prime exponents, so that equal powers written differently are recognised (lvl.same_pf)
    if x.kind == "S" and len(x.terms) == 1:
        from fractions import Fraction as Fr
        (r, c), = x.terms.items()
        c = Fr(c)
        f = lvl.pfactor(c) if c > 0 else None
        if f is not None:
            pfv = {p_: Fr(a) for p_, a in f.items()}
            for p_, e in r:
                pfv[p_] = pfv.get(p_, Fr(0)) + Fr(e)
            pfv = {p_: a for p_, a in pfv.items() if a != 0}
            if pfv:
                v.pf = pfv
    return v


class Checker:
    def __init__(self, q):
        self.q = q
        self.memo = {}

    def tab(self, S, kind, j):
        key = (S, kind, j)
        if key in self.memo:
            return self.memo[key]
        q = self.q
        if kind in ("V", 0) and len(S) > q.tv:
            # a V lookup above the table size is a lookup of the literals of S (QX.index, QX.literals)
            from fractions import Fraction as Fr
            e = q.index(S, "V").e
            out = [(lvl.exact(Fr(int(e.P[j]), int(e.Q[j]))), f"{int(e.P[j])}/{int(e.Q[j])}")]
        elif kind in ("V", 0):
            out = [(leaf(x), e) for (x, e) in q.members(S, q.V(S), j)]
        else:
            ent = q.N(S)
            out = []
            for (si, f, I, J, k1, k2) in ent.allprov[j]:
                S1, S2, _ = q.splits(S)[si]
                A = self.tab(S1, "V" if k1 == 0 else "N", I)
                if k2 == 0:
                    B = self.tab(S2, "V", J)
                elif k2 == 1:
                    B = self.tab(S2, "N", J)
                else:
                    nv = len(q.V(S2))
                    B = self.tab(S2, "V", J) if J < nv else self.tab(S2, "N", J - nv)
                out += combine(ivl.FORMS[f], A, B)
        self.memo[key] = out
        return out

    def hit(self, pid):
        p = qx.PROV[pid]
        if p[0] == "zero":
            return [(lvl.exact(0), "0")]
        if p[0] == "near1":                    # an exact member of V(S) within 1e-300 of 1 (qx.near1): its writing
            import exprparse
            S = p[1]; e = self.q._near1_rows(S)[p[2]][6]
            v = exprparse.to_x(e)
            return [(leaf(v), e)]
        if p[0] == "tab":
            _, S, kind, j = p
            return self.tab(S, kind, j)
        _, S1, k1, ii, f, S2, k2, cp = p
        return combine(ivl.FORMS[f], self.tab(S1, k1, ii), self.hit(cp))


def combine(f, A, B):
    out = []
    for (a, ea) in A:
        for (b, eb) in B:
            e = f"[{f} | {ea} | {eb}]"
            if a is None or b is None:
                out.append((None, e))
                continue
            try:
                out.append((lvl.op(f, a, b), e))
            except lvl.Invalid:
                pass
            except lvl.Fail as er:
                out.append((None, e + " FAIL: " + str(er)))
    return out


PRECS = (400, 800, 1600, 3200, 6400, 12800, 25600, 51200)      # lvl working precisions tried for an undecided hit

# HC_SPLITS="pi:3,17;gamma:5": for these constants, rerun only those top-level splits (others: full query)
SPLITS = {}
for _part in filter(None, os.environ.get("HC_SPLITS", "").split(";")):
    _c, _l = _part.split(":")
    SPLITS[_c] = set(int(x) for x in _l.split(",") if x != "")


def check_hits(C, H, c, vw):
    """every value of every hit of H (constant c, value window vw) evaluated by lvl, with more bits for an undecided one:
    (number of undecided values, detail records); a value is undecided when lvl fails or does not place it outside vw"""
    res, undecided = [], 0
    for i in range(len(H)):
        vals = C.hit(int(H.pid[i]))
        bad = [(e, repr(v)) for (v, e) in vals if v is None or not lvl.outside(v, *vw)]
        prec = lvl.PREC
        for p in PRECS[1:]:                  # undecided: the same evaluation with more bits (lvl.set_prec)
            if not bad:
                break
            lvl.set_prec(p); C.memo.clear(); prec = p
            vals = C.hit(int(H.pid[i]))
            bad = [(e, repr(v)) for (v, e) in vals if v is None or not lvl.outside(v, *vw)]
        if lvl.PREC != PRECS[0]:
            lvl.set_prec(PRECS[0]); C.memo.clear()
        undecided += len(bad)
        res.append({"const": c, "enclosure_ln": [float(H.A[i]), float(H.B[i])], "values": len(vals), "prec_bits": prec,
                    "undecided": bad[:5], "sample": [[e[:400], repr(v)] for (v, e) in vals[:3]]})
        print(f"[progress] {c} hit {i + 1}/{len(H)} ({100 * (i + 1) / len(H):.1f}%) values {len(vals)} "
              f"undecided {len(bad)}", (bad[0][0][:300] if bad else ""), flush=True)
    return undecided, res


def certainly_inside(v, lo, hi):
    """True if the magnitude v is certainly inside [lo, hi] (an exact rational, or an enclosure, within it): such a value
    stays inside at any precision; used only to stop raising the precision and to report it"""
    lo, hi = I_(lo), I_(hi)
    if v is None or v.far:
        return False
    if v.q is not None:
        x = lvl.ivq(v.q)
        return x.a >= lo.b and x.b <= hi.a
    return v.x is not None and v.x.a >= lo.b and v.x.b <= hi.a


def I_(a):
    return lvl.I(a)


def check_hits_multi(C, H, c, vws):
    """(11/11) check_hits against several value windows at once, vws = {label: window}: every value of every hit is
    evaluated once by lvl and tested against each window, with more bits while some value is undecided for some
    window and not certainly inside it (a value certainly inside a window stays so).  The windows are nested in use
    (C7: radius E*; C11: radius E11 >= E*), and the hits are those of the widest one, so the hits of every narrower
    window are among them.  Returns ({label: (undecided, certainly inside)}, detail records)."""
    res = []
    tot = {k: [0, 0] for k in vws}
    for i in range(len(H)):
        prec = lvl.PREC
        for p in PRECS:
            if p != PRECS[0]:
                lvl.set_prec(p); C.memo.clear(); prec = p
            vals = C.hit(int(H.pid[i]))
            bad = {k: [(e, v) for (v, e) in vals if v is None or not lvl.outside(v, *w)] for k, w in vws.items()}
            if all(all(certainly_inside(v, *vws[k]) for (e, v) in b) for k, b in bad.items()):
                break
        if lvl.PREC != PRECS[0]:
            lvl.set_prec(PRECS[0]); C.memo.clear()
        rec = {"const": c, "enclosure_ln": [float(H.A[i]), float(H.B[i])], "values": len(vals), "prec_bits": prec,
               "sample": [[e[:400], repr(v)] for (v, e) in vals[:3]], "undecided": []}
        for k, b in bad.items():
            ins = sum(1 for (e, v) in b if certainly_inside(v, *vws[k]))
            tot[k][0] += len(b); tot[k][1] += ins
            rec["undecided_" + k] = [(e, repr(v)) for (e, v) in b[:5]]
            rec["inside_" + k] = ins
            if b:
                rec["undecided"] = rec["undecided"] or [(e, repr(v)) for (e, v) in b[:5]]
        res.append(rec)
        print(f"[progress] {c} hit {i + 1}/{len(H)} ({100 * (i + 1) / len(H):.1f}%) values {len(vals)} " +
              " ".join(f"undecided_{k} {len(b)}" for k, b in bad.items()), flush=True)
    return {k: tuple(v) for k, v in tot.items()}, res


def value_window(c, E):
    """the value window [c - E, c + E] (doubles, outward)"""
    x = const(c)
    return (np.nextafter(float((x - arb(E)).lower()), 0), np.nextafter(float((x + arb(E)).upper()), np.inf))


def main():
    n, o = int(sys.argv[1]), int(sys.argv[2])
    consts = sys.argv[3].split(",")
    outf = sys.argv[4] if len(sys.argv) > 4 else os.path.join(paths.C5, f"hitcheck_{n}_{o}.json")
    cert = {}
    for l in open(os.path.join(paths.RESULTS, "certified.jsonl")):
        r = json.loads(l); cert[(r["const"], r["n"], r["ord"])] = r
    from cells_x import window
    qx.DEBUG = True
    tv, tn = (5, 5) if o else (4, 5)
    q = qx.QX(n, o, os.path.join(os.environ.get("RESDIR", paths.SRC), f"res_{n}_{o}.txt"), tv, tn)
    S = tuple(range(1, n + 1))
    C = Checker(q)
    res = []
    undecided = 0
    nhits = 0
    for c in consts:                     # one window at a time: the debug provenance is cleared in between
        qx.PROV = qx._Prov()
        r = cert.get((c, n, o))
        E = (float(arb(r["optimal"][0]["error"]).upper()) if r else float(os.environ["E_" + c])) * (1 + 1e-9)
        wlo, whi = window(c, E)
        vw = value_window(c, E)
        t0 = time.time(); q.top = S
        sp = SPLITS.get(c)
        if sp is not None:                   # only the top-level splits that produced hits (cells_qpar)
            import cells_qpar
            H = cells_qpar.part_hits(q, S, np.array([wlo]), np.array([whi]), 0, 1, select=lambda i: i in sp, label=c)
        else:
            H = q.Q(S, "N", np.array([wlo]), np.array([whi]))
        print(f"{c}: hits {len(H)} ({time.time() - t0:.0f}s)", flush=True)
        nhits += len(H)
        und_c, det_c = check_hits(C, H, c, vw)
        undecided += und_c
        res += det_c
        with open(outf + ".partial.jsonl", "a") as f:     # results of this constant, written at once
            f.write(json.dumps({"const": c, "hits": len(H), "undecided": sum(len(r["undecided"]) for r in res if r["const"] == c),
                                "detail": [r for r in res if r["const"] == c]}) + "\n")
        del H
    json.dump({"n": n, "ord": o, "hits": nhits, "undecided": undecided, "detail": res}, open(outf, "w"), indent=1)
    print("TOTAL hits", nhits, "undecided values", undecided)


if __name__ == "__main__":
    main()
