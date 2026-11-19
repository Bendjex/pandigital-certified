"""Check that the filtered 4-digit enumeration V4big of e_F_09.py loses nothing: compare it with the
FULL 4-digit value set e_F.V(S) on test blocks, for thresholds low enough that many values qualify."""
import sys, time
import e_F, e_F_09
from e_F import lnabs
bad = 0
for S, TX, TA in [(("0", "2", "3", "4"), 30.0, 3.0), (("2", "5", "7", "9"), 60.0, 4.0), (("0", "3", "6", "8"), 100.0, 4.5)]:
    t = time.time()
    full = e_F.V(S)
    def key(v, s): return v if v[0] != "T" else ("T", s)
    for mode, T in (("X", TX), ("A", TA)):
        need = set()
        for (v, s) in full.values():
            if v[0] == "Z": continue
            L = lnabs(v)
            if mode == "X":
                if L.lower() >= T: need.add(key(v, s))
            else:
                if not L.contains(0) and abs(L).log().lower() >= T: need.add(key(v, s))
        got = set(e_F_09.V4big(S, mode, T).keys())
        # T-values are keyed by expression; compare those by numeric value instead
        def num(k, table):
            v = table[k][0]
            return (v[1], e_F.lnabs(v).mid().str(25, radius=False))
        needR = {k for k in need if k[0] != "T"}; gotR = {k for k in got if k[0] != "T"}
        needT = {num(k, full) for k in need if k[0] == "T"}
        V4 = e_F_09.V4big(S, mode, T)
        gotT = {num(k, V4) for k in got if k[0] == "T"}
        miss = len(needR - gotR) + len(needT - gotT)
        bad += miss
        print(S, mode, T, "full-set values above threshold:", len(need), " in V4big:", len(got), " missing:", miss, f"({time.time()-t:.0f}s)", flush=True)
print("TOTAL MISSING:", bad)
