"""(14/11) qx.rec_chunks: the distinct windows of a recursive sub-query in adaptive batches (QX_REC_WIN0, QX_REC_PAIRS).
For random windows on a set above the table sizes (some repeated bit for bit), the query must return exactly the
same hits (window index and every field of the enclosure, as multisets) with one batch for everything (the code
before) and with tiny batches and tiny chunks (many batches at every level of the recursion).
Windows: narrow, in range, far out of range on both sides, near 1 and at simple values (as t_batchrec.py).
usage: t_recbatch.py n digits kind nwin [seed]      e.g. t_recbatch.py 8 234678 N 40"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import paths, qx

n = int(sys.argv[1]); S = tuple(int(c) for c in sys.argv[2]); kind = sys.argv[3]; NW = int(sys.argv[4])
rng = np.random.default_rng(int(sys.argv[5]) if len(sys.argv) > 5 else 5)
q = qx.QX(n, 0, os.path.join(os.environ.get("RESDIR", paths.SRC), f"res_{n}_0.txt"), 4, 5, verbose=False)
assert not q.tables_ok(S, kind)
k = rng.integers(0, 5, NW)
c = np.where(k == 0, rng.normal(0, 3, NW),
    np.where(k == 1, rng.choice([-1, 1], NW) * rng.uniform(700, 4000, NW),
    np.where(k == 2, rng.choice([-1, 1], NW) * 10.0 ** rng.uniform(-14, -3, NW),
    np.where(k == 3, np.log(rng.choice([0.5, 2, 3, 0.25, 1.5, 10, 0.1, 8, 0.8], NW)), rng.uniform(-690, 690, NW)))))
w = np.where(k == 0, 1e-9, np.where(k == 1, 10.0 ** rng.uniform(-9, -6, NW), np.where(k == 2, np.abs(c) * 1e-6, np.where(k == 3, 3e-14, 1e-9))))
lo, hi = c - w, c + w
dup = rng.integers(0, NW, max(2, NW // 4))
lo = np.concatenate([lo, lo[dup]]); hi = np.concatenate([hi, hi[dup]])


def rows(H):
    if not len(H):
        return []
    cols = [np.ascontiguousarray(H.w).view(np.uint8).reshape(len(H), -1)]
    cols += [np.ascontiguousarray(getattr(H, f)).view(np.uint8).reshape(len(H), -1) for f in qx.Hits.F]
    return sorted(bytes(r) for r in np.concatenate(cols, 1))


calls = [0]
_rc = qx.QX.rec_chunks
def rec_chunks(self, *a):
    for h in _rc(self, *a):
        calls[0] += 1
        yield h
qx.QX.rec_chunks = rec_chunks
res = {}
for (w0, rp) in ((10 ** 12, 10 ** 12), (1, 50), (7, 3000)):
    qx.REC_WIN0, qx.REC_PAIRS = w0, rp
    qx.PROV = qx._Prov()
    calls[0] = 0
    t = time.time()
    H = q.Q(S, kind, lo, hi)
    res[(w0, rp)] = rows(H)
    print(f"REC_WIN0={w0} REC_PAIRS={rp}: {len(H)} hits, {calls[0]} chunks, {time.time() - t:.1f} s", flush=True)
ref = res[(10 ** 12, 10 ** 12)]
bad = sum(int(v != ref) for v in res.values())
print(f"t_recbatch: {len(ref)} hits; modes different from one batch: {bad}; PROBLEMS {bad}")
sys.exit(1 if bad or not ref else 0)
