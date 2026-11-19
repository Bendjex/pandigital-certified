"""build (and store in the cache QX_NCACHE) the N tables of the sets of 4 and 5 digits of {1..n}, part k of J
(the sets with index = k mod J), so that the C5 processes find them ready.  usage: ntables.py n ord tV tN k J"""
import sys, os, time, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, motorx, qx
if os.environ.get("QX_ALLPROV") == "1":
    # (03/11) the tables of the provenance mode used by hitcheck.py / c5run.py (every derivation of each entry kept,
    # cached apart as files N H...): the same patch, applied by importing hitcheck
    import hitcheck
    assert motorx.ALLPROV[0]
n, o, tv, tn, k, J = map(int, sys.argv[1:7])
assert motorx.NCACHE, "QX_NCACHE must be set"
RESDIR = os.environ.get("RESDIR", paths.SRC)
# qx.QX: the same tables (motorx.MX.N), with motor B's answers refused if motor B left an operation undecided
q = qx.QX(n, o, os.path.join(RESDIR, f"res_{n}_{o}.txt"), tv, tn)
sets = [S for m in (4, 5) if m <= tn for S in q.sets(m)]
t0 = time.time()
mine = [S for i, S in enumerate(sets) if i % J == k]
for c, S in enumerate(mine):
    q.N(S)
    q.Nt = {T: e for T, e in q.Nt.items() if len(T) < 4}          # keep the memory low: big tables are on disk
    print(f"[progress] ntables part {k}/{J}: {c + 1}/{len(mine)} ({100 * (c + 1) / len(mine):.1f}%), {time.time() - t0:.0f}s", flush=True)
print("done", len(mine), "tables", round(time.time() - t0, 1), "s; motor B undecided operations:", q.motorb_unresolved())
