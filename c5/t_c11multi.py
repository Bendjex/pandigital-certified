"""hitcheck.check_hits_multi (11/11), the check of the hits of one query against several value windows (C7: E*, C11:
E11 >= E*), against hitcheck.check_hits window by window.  The hits: NE((1..5)) in order around 1 (the values of
t_hc_smoke.py, many of them inside).  For the window of the query, a nested narrower one and a window around 2 (no
value there): the undecided count of check_hits_multi must equal that of check_hits for each window, the values
certainly inside must be among the undecided ones, and around 2 there must be none.
usage: t_c11multi.py"""
import sys, numpy as np
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, hitcheck, qx, ivl, lvl
qx.DEBUG = True
q = qx.QX(5, 1, os.path.join(paths.SRC, "res_5_1.txt"), 5, 5, verbose=False)
lo, hi = 1 - 1e-9, 1 + 1e-9
wl = float(ivl.ln_lo(np.array([lo]))[0]); wh = float(ivl.ln_hi(np.array([hi]))[0])
H = q.Q((1, 2, 3, 4, 5), "N", np.array([wl]), np.array([wh]))
H = H.sel(np.arange(0, len(H), max(1, len(H) // 40)))     # 40 of them (check_hits raises every inside value to 51 200 bits)
C = hitcheck.Checker(q)
vws = {"wide": (lo, hi), "narrow": (1 - 1e-12, 1 + 1e-12), "two": (2 - 1e-9, 2 + 1e-9)}
tot, det = hitcheck.check_hits_multi(C, H, "one", vws)
bad = 0
for k, vw in vws.items():
    C.memo.clear()
    und, _ = hitcheck.check_hits(C, H, "one", vw)
    print(f"{k}: check_hits undecided {und}; check_hits_multi undecided {tot[k][0]}, certainly inside {tot[k][1]}")
    bad += (und != tot[k][0]) + (tot[k][1] > tot[k][0])
bad += tot["two"][0] != 0
bad += not (tot["wide"][0] >= tot["narrow"][0] > 0)
print(f"t_c11multi: {len(H)} hits; problems {bad}")
sys.exit(1 if bad else 0)
