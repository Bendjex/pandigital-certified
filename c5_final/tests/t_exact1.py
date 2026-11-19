"""qx.pre2, item 41 (14/11): a base exactly 1 against a window whose second-level constraint gives the sign of ln w.
For the entries of real V tables (with the exact 1 of V(1), V(19), ...) and random windows (ln bounds) with random
second-level constraints (sign -1, 0 or +1 and bounds of ln|ln w|), pre2 with the rule (QX_EXACT1_C2) and without it:
  - for an exact 1 (flag O), forms s^l and s^-l, constraint with a sign: the windows of l are empty with the rule
    (s^(+-l) = 1 for every l, and the constraint excludes 1);
  - every other (entry, form, window): the same windows bit for bit.
usage: t_exact1.py [windows] [seed]"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, qx
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 41)
q = qx.QX(9, 0, os.path.join(os.environ.get("RESDIR", paths.SRC), "res_9_0.txt"), 4, 5, verbose=False)
bad = changed = kept = ones = 0
for S in [(1,), (1, 9), (1, 2), (1, 5, 7), (2, 4), (3, 8, 9)]:
    e = q.V(S)
    I = rng.integers(0, len(e), size=N)
    one = np.where(e.O & ~e.Z)[0]
    if len(one):
        I[: N // 3] = rng.choice(one, size=N // 3)          # many exact 1
    x7 = e.q7(I)
    a = rng.uniform(-50, 50, N); w = rng.exponential(1.0, N) * rng.choice([1e-12, 1e-3, 1, 30], N)
    al, be = a, a + w
    near = rng.random(N) < 0.3                             # windows around ln 1 = 0, as the rounded ones
    al = np.where(near, -5e-324, al); be = np.where(near, rng.choice([6e-309, 1e-10], N), be)
    sg = rng.choice([-1, 0, 1], N).astype(np.int8)
    m2 = rng.uniform(-800, 5, N); M2 = m2 + rng.exponential(2, N)
    c2 = (sg, m2, M2)
    for f in range(9):
        qx.EXACT1_C2 = False
        off = qx.pre2(f, x7, al, be, c2)
        qx.EXACT1_C2 = True
        on = qx.pre2(f, x7, al, be, c2)
        sgw, _, _, _ = qx._wlog(al, be, c2)
        target = (f in (5, 6)) & x7.O & ~x7.Z & (sgw != 0)
        for (lo0, hi0), (lo1, hi1) in zip(off, on):
            same = (lo0 == lo1) & (hi0 == hi1)
            empty = lo1 > hi1
            if not np.all(same | target):
                bad += 1; print("CHANGED A WINDOW OUTSIDE THE RULE:", S, f)
            if not np.all(empty[target]):
                bad += 1; print("A WINDOW OF THE RULE NOT EMPTY:", S, f)
            changed += int((target & ~same).sum()); kept += int((~target).sum()); ones += int(target.sum())
print(f"t_exact1: {ones} (exact 1, s^l or s^-l, signed window) windows emptied ({changed} of them were not empty "
      f"before), {kept} other windows unchanged; PROBLEMS {bad}")
sys.exit(1 if bad else 0)
