"""(16/11) qx.near_x_ys (resolve_vv "out" with a tiny window): soundness against the general query on real cases.
Runs motor X's query on top-level splits of n (all pieces, the windows of pi, gamma, delta as c5run) with resolve_vv
wrapped: for every recursive l (against a tiny window, the case of near_x_ys) the general list of exact members of
l (exact_vals on l's own enclosure, the code before) is computed too, and every pair (x, y) of a member x of the
class of s and a member y of that list whose exact difference is nonzero and has ln|x - y| meeting the window of the
pair must have y in the new list.  The hits with and without the rule are counted (the rule may narrow the enclosure of a hit, so
they are reported, not compared).
usage: t_nearx.py n splits(comma) [maxpairs] [pieces(comma)]   (T_NEARX_ONLY_NEW=1: no run without the rule)"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hitcheck, paths, qx, c5run
from cells_x import window
import exact as ex
from exact import Unresolved
n = int(sys.argv[1]); SPL = [int(x) for x in sys.argv[2].split(",")]
MAXP = int(sys.argv[3]) if len(sys.argv) > 3 else 200000
PIECES = sys.argv[4].split(",") if len(sys.argv) > 4 else None       # only these pieces of the plan
ONLY_NEW = os.environ.get("T_NEARX_ONLY_NEW") == "1"                   # the run with the rule only
MAXCASES = int(os.environ.get("T_NEARX_MAXCASES", "200"))              # cases compared with the general query
OLD_S = int(os.environ.get("T_NEARX_OLD_S", "120"))                    # seconds allowed to one general query
LOG = open(os.environ.get("T_NEARX_LOG", "t_nearx_cases.log"), "a", buffering=1)
q = qx.QX(n, 0, os.path.join(os.environ.get("RESDIR", paths.SRC), f"res_{n}_0.txt"), 4, 5, verbose=False)
S = tuple(range(1, n + 1)); sides = q.table_side(S)
E = {"pi": 3.23e-14, "gamma": 2.12e-15, "delta": 2.2e-14} if n == 9 else {"pi": 1.7754842022891951e-12, "gamma": 6.083860837112889e-13, "delta": 3.0580735760575984019e-11}
lo = np.array([window(c, E[c] * (1 + 1e-9))[0] for c in E]); hi = np.array([window(c, E[c] * (1 + 1e-9))[1] for c in E])
stats = {"groups": 0, "ts": 0, "pairs": 0, "in_window": 0, "missing": 0, "old_unres": 0}
_rvv = qx.QX.resolve_vv


def check(self, S1, S2, I, sub, amb, meet, wlo, whi, r, want="out"):
    if want == "out":
        e1 = self.V(S1)
        for t in np.where(amb)[0]:
            if int(sub.ti[t]) >= 0 or not (np.isfinite(whi[t]) and whi[t] <= qx.NEAR_X_LN):
                continue
            new = self.near_x_ys(S1, e1, int(I[t]), S2, float(whi[t]))
            if new is None:
                continue
            stats["ts"] += 1
            if float(sub.A[t]) <= 0 <= float(sub.B[t]):
                # l pegged to 1: the general query is the one that cannot be afforded (split 20); this case (x = 1,
                # two signed queries with |ln y| <= 2u) is covered by its proof and t_onediff.py
                stats["near1_skipped"] = stats.get("near1_skipped", 0) + 1
                continue
            if stats.get("checked", 0) >= MAXCASES:
                continue
            newk = {v.key() for (v, e) in new}
            import signal
            def _alarm(*_):
                raise TimeoutError
            old_h = signal.signal(signal.SIGALRM, _alarm); signal.alarm(OLD_S)
            try:
                old = self.exact_vals(S2, float(sub.A[t]), float(sub.B[t]))
                xs = self.members(S1, e1, int(I[t]))
            except Unresolved:
                stats["old_unres"] += 1; continue
            except TimeoutError:
                stats["old_timeout"] = stats.get("old_timeout", 0) + 1
                LOG.write(f"timeout old query S2={S2} [{float(sub.A[t])!r}, {float(sub.B[t])!r}]\n"); continue
            finally:
                signal.alarm(0); signal.signal(signal.SIGALRM, old_h)
            stats["checked"] = stats.get("checked", 0) + 1
            p0, w0, m0 = stats["pairs"], stats["in_window"], stats["missing"]
            for (x, xe) in xs:
                for (y, ye) in old:
                    if stats["pairs"] >= MAXP:
                        break
                    stats["pairs"] += 1
                    d = ex.sub(x, y)
                    try:
                        sg = ex.sign(d)
                    except Unresolved:
                        sg = None
                    if sg == 0:
                        continue
                    try:
                        b = ex.tight_ball(ex.absval(d))
                        from flint import ctx
                        o = ctx.prec; ctx.prec = 256; lb = b.log(); ctx.prec = o
                        a_, b_ = float(lb.lower()), float(lb.upper())
                    except Exception:
                        a_, b_ = -np.inf, np.inf
                    if b_ < wlo[t] or a_ > whi[t]:
                        continue
                    stats["in_window"] += 1
                    if y.key() not in newk:
                        stats["missing"] += 1
                        print("MISSING", S1, S2, xe, ye, a_, b_, wlo[t], whi[t], flush=True)
            LOG.write(f"case {stats['checked']}: S1={S1} S2={S2} l in [{float(sub.A[t])!r}, {float(sub.B[t])!r}] window ln u <= "
                      f"{float(whi[t])!r}: {len(xs)} x, {len(old)} old y, {len(new)} new y, pairs {stats['pairs'] - p0}, in window "
                      f"{stats['in_window'] - w0}, missing {stats['missing'] - m0}\n")
    return _rvv(self, S1, S2, I, sub, amb, meet, wlo, whi, r, want)


def run(rule):
    qx.NEAR_X = rule
    qx.QX.resolve_vv = check if rule else _rvv
    hits = []
    for sp in SPL:
        S1, S2, forms = sides[sp]
        for name, paths_ in c5run.plan(q, S1, S2):
            if PIECES and name not in PIECES:
                continue
            q.top = S
            H = qx.Hits.cat([c5run.path_hits(q, S1, S2, forms, part, only, lo, hi) for (part, only) in paths_])
            cols = [np.ascontiguousarray(H.w).view(np.uint8).reshape(len(H), -1)] if len(H) else []
            cols += [np.ascontiguousarray(getattr(H, f)).view(np.uint8).reshape(len(H), -1) for f in qx.Hits.F] if len(H) else []
            hits += sorted(bytes(x) for x in np.concatenate(cols, 1)) if len(H) else []
        print(f"rule={rule} split {sp}: {len(hits)} hits so far", flush=True)
    return sorted(hits)


qx.DEBUG = True
t0 = time.time()
h_new = run(True)
for k in ("_nearx", "_xV", "_xv", "_vvmemo"):
    self_ = q.__dict__.pop(k, None)
h_old = run(False) if not ONLY_NEW else h_new
import collections
cn, co = collections.Counter(h_new), collections.Counter(h_old)
extra = sum((cn - co).values()); dropped = sum((co - cn).values())
print(f"t_nearx: {stats}; hits with the rule {len(h_new)}, without {len(h_old)}, new hits not in the old {extra}, "
      f"old hits removed {dropped}; {time.time() - t0:.0f}s")
bad = stats["missing"]          # (a narrower hull changes the enclosure bytes of a hit: extra/dropped are reported only)
print(f"t_nearx: PROBLEMS {bad}")
sys.exit(1 if bad else 0)
