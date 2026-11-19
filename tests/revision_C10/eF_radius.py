import sys, runpy, itertools
sys.argv = ["e_F.py", sys.argv[1], sys.argv[2]]
ns = runpy.run_path("e_F.py", run_name="not_main")
DIG = tuple(sys.argv[1]); theta = float(sys.argv[2])
side_entries = ns["side_entries"]
worst = 0.0; n = 0; kinds = {}
for k in (3, 4):
    for T in itertools.combinations(DIG, k):
        for (lm, sg, a, x, sa, sx) in side_entries(T, theta):
            n += 1
            rr = float(lm.rad()) / max(1.0, abs(float(lm.mid())))
            kinds[(a[0], x[0])] = kinds.get((a[0], x[0]), 0) + 1
            if rr > worst: worst = rr; wex = (sa, sx, lm.str(10))
print("digits", sys.argv[1], "theta", theta, ": side entries", n, "kinds", kinds, " worst relative radius of ln|mu|:", worst, wex if n else "")
