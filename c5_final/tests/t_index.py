"""IndexT against the original Index: identical sets of (window, entry) pairs, on random enclosure tables with
every kind of entry (narrow, wide, huge, half-lines, whole line, zeros) and random windows (also -inf/inf)"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qx, ivl
from motorx import Ent
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 11)
bad = 0
for trial in range(int(sys.argv[1]) if len(sys.argv) > 1 else 200):
    n = int(rng.integers(1, 3000))
    c = rng.choice([-1, 1], n) * 10.0 ** rng.uniform(-20, 305, n) * rng.random(n)
    c = np.where(rng.random(n) < 0.3, rng.normal(0, 3, n), c)
    w = 10.0 ** rng.uniform(-18, 300, n) * rng.random(n)
    w = np.where(rng.random(n) < 0.5, np.abs(c) * 10.0 ** rng.uniform(-17, -3, n), w)
    A, B = c - w, c + w
    k = rng.random(n)
    A = np.where(k < 0.05, -np.inf, A); B = np.where((k > 0.05) & (k < 0.1), np.inf, B)
    A = np.where((k > 0.1) & (k < 0.12), -np.inf, A); B = np.where((k > 0.1) & (k < 0.12), np.inf, B)
    Z = (k > 0.12) & (k < 0.15)
    A = np.where(Z, -np.inf, A); B = np.where(Z, -np.inf, B)
    B = np.where((k > 0.15) & (k < 0.16), -np.inf, B); A = np.where((k > 0.15) & (k < 0.16), -np.inf, A)
    e = Ent(A, B, Z, np.zeros(n, bool))
    m = int(rng.integers(1, 400))
    cw = rng.choice([-1, 1], m) * 10.0 ** rng.uniform(-20, 305, m)
    cw = np.where(rng.random(m) < 0.4, c[rng.integers(0, n, m)], cw)
    ww = 10.0 ** rng.uniform(-18, 300, m) * rng.random(m)
    lo, hi = cw - ww, cw + ww
    kk = rng.random(m)
    lo = np.where(kk < 0.15, -np.inf, lo); hi = np.where((kk > 0.15) & (kk < 0.25), np.inf, hi)
    lo = np.where((kk > 0.25) & (kk < 0.3), np.inf, lo); hi = np.where((kk > 0.25) & (kk < 0.3), -np.inf, hi)   # empty windows
    a = qx.Index(e).look(lo, hi); b = qx.IndexT(e).look(lo, hi)
    sa = set(zip(a[0].tolist(), a[1].tolist())); sb = set(zip(b[0].tolist(), b[1].tolist()))
    if sa != sb or len(a[0]) != len(sa) or len(b[0]) != len(sb):
        bad += 1
        print("DIFF trial", trial, "old", len(sa), "new", len(sb), "only old", list(sa - sb)[:3], "only new", list(sb - sa)[:3])
print("trials", trial + 1, "differences", bad)

# second family (ixT2, 03/11): one-signed entries of every magnitude with tiny relative width, entries touching 0,
# signed zeros and subnormals, windows near 0 and at huge magnitudes; IndexT against the plain test itself
bad2 = 0
T2 = int(sys.argv[3]) if len(sys.argv) > 3 else 300
for trial in range(T2):
    n = int(rng.integers(1, 4000))
    mag = 10.0 ** rng.uniform(-320, 308, n)
    c = rng.choice([-1, 1], n) * mag
    rel = 10.0 ** rng.uniform(-16, 0, n)
    w = np.abs(c) * rel * rng.random(n)
    A, B = c - w, c + w
    k = rng.random(n)
    A = np.where(k < 0.03, 0.0, A); B = np.where((k > 0.03) & (k < 0.06), -0.0, B)
    A = np.where((k > 0.06) & (k < 0.08), -0.0, A)
    A = np.where((k > 0.08) & (k < 0.1), np.nextafter(0, 1), A)
    B = np.maximum(A, B)
    A = np.where((k > 0.1) & (k < 0.12), -np.inf, A); B = np.where((k > 0.12) & (k < 0.14), np.inf, B)
    Z = (k > 0.14) & (k < 0.15)
    A = np.where(Z, -np.inf, A); B = np.where(Z, -np.inf, B)
    A = np.where((k > 0.15) & (k < 0.155), np.nan, A)
    e = Ent(A, B, Z, np.zeros(n, bool))
    m = int(rng.integers(1, 500))
    cw = rng.choice([-1, 1], m) * 10.0 ** rng.uniform(-320, 308, m)
    cw = np.where(rng.random(m) < 0.4, np.where(np.isfinite(A[rng.integers(0, n, m)]), A[rng.integers(0, n, m)], 0.0), cw)
    ww = np.abs(cw) * 10.0 ** rng.uniform(-17, 1, m) * rng.random(m)
    ww = np.where(rng.random(m) < 0.1, 10.0 ** rng.uniform(-300, 300, m), ww)
    lo, hi = cw - ww, cw + ww
    kk = rng.random(m)
    lo = np.where(kk < 0.1, -np.inf, lo); hi = np.where((kk > 0.1) & (kk < 0.2), np.inf, hi)
    lo = np.where((kk > 0.2) & (kk < 0.25), np.inf, lo); hi = np.where((kk > 0.2) & (kk < 0.25), -np.inf, hi)
    lo = np.where((kk > 0.25) & (kk < 0.28), 0.0, lo); hi = np.where((kk > 0.28) & (kk < 0.31), -0.0, hi)
    lo = np.where((kk > 0.31) & (kk < 0.32), np.nan, lo)
    b = qx.IndexT(e).look(lo, hi)
    sb = set(zip(b[0].tolist(), b[1].tolist()))
    ref = set()
    with np.errstate(invalid="ignore"):
        for j in range(m):
            mm = np.where(Z, lo[j] == -np.inf, (B >= lo[j]) & (A <= hi[j]))
            ref |= {(j, int(i)) for i in np.where(mm)[0]}
    if sb != ref or len(b[0]) != len(sb):
        bad2 += 1
        print("DIFF2 trial", trial, "index", len(sb), "plain", len(ref), "missing", list(ref - sb)[:3], "extra", list(sb - ref)[:3])
print("trials", T2, "(plain test, one-signed entries of every magnitude) differences", bad2)

# third family (10/11): the bounded batches.  cand_counts bounds the pairs of each window from above; the groups
# of qx._groups are consecutive, cover every window once and hold at most `limit` candidates (or a single window);
# look_batches returns, group by group, exactly the pairs of look
bad3 = 0
T3 = int(sys.argv[4]) if len(sys.argv) > 4 else 300
for trial in range(T3):
    n = int(rng.integers(1, 4000))
    mag = 10.0 ** rng.uniform(-320, 308, n)
    c = rng.choice([-1, 1], n) * mag
    c = np.where(rng.random(n) < 0.3, rng.normal(0, 3, n), c)
    w = np.abs(c) * 10.0 ** rng.uniform(-16, 0, n) * rng.random(n)
    w = np.where(rng.random(n) < 0.2, 10.0 ** rng.uniform(-18, 300, n), w)
    A, B = c - w, c + w
    k = rng.random(n)
    A = np.where(k < 0.05, -np.inf, A); B = np.where((k > 0.05) & (k < 0.1), np.inf, B)
    A = np.where((k > 0.1) & (k < 0.12), -np.inf, A); B = np.where((k > 0.1) & (k < 0.12), np.inf, B)
    Z = (k > 0.12) & (k < 0.15)
    A = np.where(Z, -np.inf, A); B = np.where(Z, -np.inf, B)
    A = np.where((k > 0.15) & (k < 0.155), np.nan, A)
    e = Ent(A, B, Z, np.zeros(n, bool))
    m = int(rng.integers(1, 500))
    cw = rng.choice([-1, 1], m) * 10.0 ** rng.uniform(-320, 308, m)
    cw = np.where(rng.random(m) < 0.5, np.where(np.isfinite(A[rng.integers(0, n, m)]), A[rng.integers(0, n, m)], 0.0), cw)
    ww = np.abs(cw) * 10.0 ** rng.uniform(-17, 1, m) * rng.random(m)
    ww = np.where(rng.random(m) < 0.2, 10.0 ** rng.uniform(-300, 300, m), ww)
    lo, hi = cw - ww, cw + ww
    kk = rng.random(m)
    lo = np.where(kk < 0.15, -np.inf, lo); hi = np.where((kk > 0.15) & (kk < 0.3), np.inf, hi)
    lo = np.where((kk > 0.3) & (kk < 0.35), np.inf, lo); hi = np.where((kk > 0.3) & (kk < 0.35), -np.inf, hi)
    lo = np.where((kk > 0.35) & (kk < 0.36), np.nan, lo); hi = np.where((kk > 0.36) & (kk < 0.37), np.nan, hi)
    ix = qx.IndexT(e)
    W, K = ix.look(lo, hi)
    cnt = ix.cand_counts(lo, hi)
    real = np.bincount(W, minlength=m)
    limit = int(rng.choice([1, 7, 100, 3000, 10 ** 6]))
    gr = qx._groups(cnt, limit)
    ok = (cnt >= real).all() and [g[0] for g in gr] == [0] + [g[1] for g in gr[:-1]] and gr[-1][1] == m \
        and all(j > i and (j == i + 1 or int(cnt[i:j].sum()) <= limit) for (i, j) in gr)
    pb = [p for p in ix.look_batches(lo, hi, limit)]
    Wb = np.concatenate([p[0] for p in pb]); Kb = np.concatenate([p[1] for p in pb])
    ok = ok and len(pb) == len(gr) and sorted(zip(W.tolist(), K.tolist())) == sorted(zip(Wb.tolist(), Kb.tolist())) \
        and all(int(real[i:j].sum()) == len(p[0]) for (i, j), p in zip(gr, pb))
    # look with a limit: the same pairs when the candidates fit, TooMany otherwise (never when limit >= candidates)
    tot = int(cnt.sum())
    for lim in (limit, tot, tot + 1, max(0, len(W) - 1)):
        try:
            Wl, Kl = ix.look(lo, hi, lim)
            ok = ok and np.array_equal(Wl, W) and np.array_equal(Kl, K) and len(W) <= lim
        except qx.TooMany:
            ok = ok and lim < tot
    if not ok:
        bad3 += 1
        print("DIFF3 trial", trial, "limit", limit, "pairs", len(W), "in batches", len(Wb), "count below real:", int((cnt < real).sum()))
print("trials", T3, "(bounded batches: cand_counts, groups, look_batches) differences", bad3)
sys.exit(1 if bad or bad2 or bad3 else 0)
