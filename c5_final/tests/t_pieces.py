"""The pieces of c5run.plan (10/11) add up to the split.
 (1) Logic: the units of a top-level split S1 | S2 are its three parts and, inside VN when NE(S2) is a recursion, the
     parts of every split of S2 and, for a split x | U with a single digit x, the parts of every split of U inside its
     VN.  Every unit must be taken by exactly one path of one piece (the predicates of the paths evaluated on it):
     checked for every top-level split of {1..9} and {1..8} (any order) and of {1..9}, {1..8} in order.
 (2) End to end (n = 8, any order, the split 1 | 2345678 and one more): the hits of the pieces, with and without the
     deep pieces T, put together, are exactly those of QX.split_hits (window index and every field of the enclosure,
     as multisets), for the window of pi of the cell and for a window 1e6 times wider.
usage: t_pieces.py [0|1: also the end-to-end part]"""
import sys, os, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, qx, c5run
from cells_x import window

bad = 0
PARTS = ("NU", "VN", "VV")


class FakeQ:
    """the table sizes and splits of QX, without tables"""
    def __init__(self, n, o):
        self.ord = bool(o); self.tv, self.tn = (5, 5) if o else (4, 5)
    tables_ok = qx.QX.tables_ok
    table_side = qx.QX.table_side


def units(q, S1, S2, deep):
    """the units of the split: tuples (part,) or ("VN", si, part) or ("VN", si, "VN", sj, part)"""
    out = [("NU",), ("VV",)]
    if q.tables_ok(S2, "N"):
        return out + [("VN",)]
    for si, t in enumerate(q.table_side(S2)):
        for p in PARTS:
            if deep and len(t[0]) == 1 and p == "VN" and not q.tables_ok(t[1], "N"):
                out += [("VN", si, "VN", sj, p2) for sj in range(len(q.table_side(t[1]))) for p2 in PARTS]
            else:
                out.append(("VN", si, p))
    return out


def takes(q, S2, part, only, u):
    """does the path (part, only) take the unit u?"""
    if u[0] != part:
        return False
    if len(u) == 1:
        return True
    pred = (only or {}).get(S2)
    if pred is not None and not pred(u[1], u[2]):
        return False
    if len(u) == 3:
        return True
    U = q.table_side(S2)[u[1]][1]
    pred = (only or {}).get(U)
    return pred is None or bool(pred(u[3], u[4]))


nsplit = npieces = nunits = 0
for n, o in ((9, 0), (8, 0), (9, 1), (8, 1), (7, 0)):
    q = FakeQ(n, o)
    for (S1, S2, forms) in q.table_side(tuple(range(1, n + 1))):
        for deep in ((None,) if (o or len(S2) < 7) else (None, True, False)):
            pieces = c5run.plan(q, S1, S2, deep)
            dp = (len(S2) == 8 and not o) if deep is None else deep
            us = units(q, S1, S2, dp)
            names = [nm for nm, _ in pieces]
            if len(set(names)) != len(names):
                bad += 1; print("DUPLICATE piece names", n, o, S1)
            for u in us:
                k = sum(takes(q, S2, part, only, u) for _, paths in pieces for (part, only) in paths)
                if k != 1:
                    bad += 1
                    if bad < 10:
                        print(f"UNIT {u} of {S1}|{S2} (n={n}, ord={o}, deep={dp}) is taken {k} times")
            nsplit += 1; npieces += len(pieces); nunits += len(us)
print(f"(1) {nsplit} plans, {npieces} pieces, {nunits} units: each unit in exactly one path; problems {bad}", flush=True)

if len(sys.argv) > 1 and sys.argv[1] == "1":
    if os.environ.get("QX_ALLPROV") == "1":
        import hitcheck
        qx.DEBUG = True
    n = 8
    q = qx.QX(n, 0, os.path.join(os.environ.get("RESDIR", paths.SRC), f"res_{n}_0.txt"), 4, 5, verbose=False)
    S = tuple(range(1, n + 1))
    sides = q.table_side(S)

    def rows(H):
        if not len(H):
            return []
        cols = [np.ascontiguousarray(H.w).view(np.uint8).reshape(len(H), -1)]
        cols += [np.ascontiguousarray(getattr(H, f)).view(np.uint8).reshape(len(H), -1) for f in qx.Hits.F]
        return sorted(bytes(r) for r in np.concatenate(cols, 1))

    E = float(os.environ.get("E_pi", "1e-12"))
    w1 = window("pi", E); w2 = window("pi", E * 1e6)
    lo = np.array([w1[0], w2[0]]); hi = np.array([w1[1], w2[1]])
    for i in (0, 8 + 5):
        S1, S2, forms = sides[i]
        qx.PROV = qx._Prov()
        ref = rows(qx.Hits.cat(q.split_hits(S1, S2, forms, lo, hi)))
        for deep in (False, True) if len(S2) >= 7 else (False,):
            got = []
            for name, paths in c5run.plan(q, S1, S2, deep):
                qx.PROV = qx._Prov()
                for (part, only) in paths:
                    got += rows(c5run.path_hits(q, S1, S2, forms, part, only, lo, hi))
            ok = sorted(got) == ref
            bad += not ok
            print(f"(2) split {i} {S1}|{S2}, deep={deep}: {len(ref)} hits, pieces {len(got)} hits, equal: {ok}", flush=True)
    bad += q.motorb_unresolved()
print("t_pieces: PROBLEMS", bad)
sys.exit(1 if bad else 0)
