"""D1: distinct-value DP (double precision, heuristic only) for the random model.
Magnitudes only (value sets are symmetric under unary minus). In-range rule as Theorem 1."""
import numpy as np, sys, itertools, time, json
LO, HI = 1e-300, 1e300
MASK = np.uint64(~((1 << 12) - 1) & 0xFFFFFFFFFFFFFFFF)   # merge values equal to ~1e-12 relative

def uniq(v):
    v = v[np.isfinite(v)]
    v = v[(v == 0) | ((v >= LO) & (v <= HI))]
    b = (v.view(np.uint64) & MASK)
    _, idx = np.unique(b, return_index=True)
    return v[idx]

def literals(digs, ordered):
    out = []
    perms = [digs] if ordered else itertools.permutations(digs)
    for p in perms:
        s = ''.join(map(str, p))
        out.append(float(s))
        for i in range(len(s)):
            out.append(float(s[:i] + '.' + s[i:]))
    return np.array(out)

def combine(A, B, ordered=False):
    a = A[:, None]; b = B[None, :]
    with np.errstate(all='ignore'):
        r = [a + b, np.abs(a - b), a * b]
        r.append(np.where(b > 0, a / np.where(b > 0, b, 1), np.nan))
        if not ordered: r.append(np.where(a > 0, b / np.where(a > 0, a, 1), np.nan))
        for x, y in (((a, b),) if ordered else ((a, b), (b, a))):
            pos = x > 0
            xx = np.where(pos, x, 1.0)
            r.append(np.where(pos, xx ** y, np.where(y > 0, 0.0, np.nan)))
            r.append(np.where(pos, xx ** (-y), np.nan))
    return np.concatenate([q.ravel() for q in r])

def run(n, ordered, consts, delta_rel=0.05, maxsize=None):
    digs = list(range(1, n + 1))
    V = {}
    res = {}
    for k in range(1, n + 1):
        if ordered:
            sets = [tuple(range(i, i + k)) for i in range(1, n - k + 2)]
        else:
            sets = list(itertools.combinations(digs, k))
        for S in sets:
            parts = [literals(list(S), ordered)]
            if ordered:
                splits = [(S[:i], S[i:]) for i in range(1, k)]
            else:
                splits = []
                for r in range(1, k):
                    for S1 in itertools.combinations(S, r):
                        if S1 < tuple(x for x in S if x not in S1):
                            splits.append((S1, tuple(x for x in S if x not in S1)))
            for S1, S2 in splits:
                parts.append(combine(V[S1], V[S2], ordered))
                if ordered:   # in order, A op B and B op A are different only via sign symmetry; a/b and b/a: keep only S1 op S2
                    pass
            V[S] = uniq(np.concatenate(parts))
        top = tuple(digs[:k]) if ordered else tuple(digs[:k])
        if top in V:
            v = V[top]
            row = {"n": k, "distinct": int(v.size)}
            for name, c in consts.items():
                d = delta_rel * c
                row[name] = int(np.count_nonzero(np.abs(v - c) <= d))
                row[name + "_best"] = float(np.min(np.abs(v - c)))
            res[k] = row
            print(json.dumps(row), flush=True)
        if maxsize and k >= maxsize: break
    return res

if __name__ == "__main__":
    n = int(sys.argv[1]); ordered = sys.argv[2] == "1"
    consts = {"pi": np.pi, "gamma": 0.5772156649015329, "delta": 4.669201609102991}
    run(n, ordered, consts)

def window_top(n, ordered, consts, delta_rel=0.05):
    """V for all proper subsets (sizes <= n-1), then only the values of the full set inside each window."""
    digs = tuple(range(1, n + 1))
    V = {}
    for k in range(1, n):
        sets = [tuple(range(i, i + k)) for i in range(1, n - k + 2)] if ordered else list(itertools.combinations(digs, k))
        for S in sets:
            cur = uniq(literals(list(S), ordered))
            for S1, S2 in splits_of(S, ordered):
                A, B = V[S1], V[S2]
                step = max(1, int(1e7 // max(1, B.size)))
                for i in range(0, A.size, step):
                    cur = uniq(np.concatenate([cur, uniq(combine(A[i:i + step], B, ordered))]))
            V[S] = cur
    S = digs
    acc = {name: [] for name in consts}
    def add(vals):
        for name, c in consts.items():
            d = delta_rel * c
            acc[name].append(vals[np.abs(vals - c) <= d])
    add(literals(list(S), ordered))
    for S1, S2 in splits_of(S, ordered):
        A, B = V[S1], V[S2]
        if A.size < B.size and not ordered:
            A, B = B, A
        step = max(1, int(2e7 // max(1, B.size)))
        for i in range(0, A.size, step):
            add(combine(A[i:i + step], B, ordered) if not (ordered and A is not V[S1]) else combine(A[i:i+step], B, ordered))
    row = {"n": n}
    for name, c in consts.items():
        v = uniq(np.concatenate(acc[name]))
        row[name] = int(v.size); row[name + "_best"] = float(np.min(np.abs(v - c)))
    print(json.dumps(row), flush=True)
    return row

def splits_of(S, ordered):
    k = len(S)
    if ordered:
        return [(S[:i], S[i:]) for i in range(1, k)]
    out = []
    for r in range(1, k):
        for S1 in itertools.combinations(S, r):
            S2 = tuple(x for x in S if x not in S1)
            if S1 < S2: out.append((S1, S2))
    return out
