"""c5run.release (13/11): emptying the caches between pieces changes no answer and leaves no cache half-emptied.
The lab's part 0 stopped (13/11 18:42, split 12, piece VN.3) with KeyError in QX._near1_rows: release emptied the rows
of near1 but not near1 itself, so near1(S) returned at once and its rows were missing.
  (1) the exact values near 1 of sets whose rows are on disk: near1_value(S, k) before and after release, equal;
  (2) the members of random entries of the V tables, exact_V of a window, and motor B's answer by window: equal
      before and after release;
  (3) release twice in a row, then everything again.
usage: t_release.py [sets] [seed]"""
import sys, os, glob, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, qx, motorx, hitcheck, c5run
N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
q = qx.QX(9, 0, os.path.join(os.environ.get("RESDIR", paths.SRC), "res_9_0.txt"), 4, 5, verbose=False)
C = hitcheck.Checker(q)
bad = 0

on_disk = sorted({tuple(int(c) for c in os.path.basename(f).split("_")[1])
                  for f in glob.glob(os.path.join(motorx.NCACHE or "", "near1", "V0_*.json"))})
n1sets = [on_disk[i] for i in rng.choice(len(on_disk), size=min(N, len(on_disk)), replace=False)] if on_disk else []
if (2, 4, 6, 8, 9) in on_disk and (2, 4, 6, 8, 9) not in n1sets:
    n1sets.append((2, 4, 6, 8, 9))                           # the set of the lab's KeyError
vsets = [tuple(sorted(s)) for s in (rng.choice(range(1, 10), size=4, replace=False) for _ in range(N))]


def snapshot():
    out = {}
    for S in n1sets:
        rows = q._near1_rows(S)
        ks = list(range(min(len(rows), 5)))
        out[("n1", S)] = [(q.near1_value(S, k)[1], q._near1_rows(S)[k][:6]) for k in ks]
    for S in vsets:
        e = q.V(S)
        js = [int(j) for j in np.random.default_rng(sum(S)).choice(len(e), size=20, replace=False)]
        out[("mem", S)] = [sorted(w for (_, w) in q.members(S, e, j)) for j in js]
        A, B = float(e.A[js[0]]), float(e.B[js[0]])
        out[("xV", S)] = sorted(w for (_, w) in q.exact_V(S, A - 1e-3, B + 1e-3))
    return out


print("near1 sets", n1sets, "V sets", vsets, flush=True)
s0 = snapshot()
for rnd in range(2):
    c5run.release(q, C)
    if rnd == 1:
        c5run.release(q, C)                                   # twice in a row
    s1 = snapshot()
    for k in s0:
        if s0[k] != s1[k]:
            bad += 1; print("DIFFERENT AFTER RELEASE:", k)
print(f"t_release: {len(n1sets)} sets near 1, {len(vsets)} V sets, {len(s0)} answers compared over 2 releases; PROBLEMS {bad}")
sys.exit(1 if bad else 0)
