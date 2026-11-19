"""
C5 measurement (approximate, float log-domain): for digit sets S (any order), compare
  R(S): magnitudes |val E| of expressions with every intermediate 0 or in [1e-300, 1e300]
  A(S): magnitudes of all valid expressions (no range limit)
and classify the in-range values of A(S) that are not (approximately) in R(S).
Values are stored as L = ln|v| (float64); 0 is L = -inf.  |L| > DEEP is clipped to +-BIGL (deep).
Not rigorous: only to learn which mechanisms exist and how many values they produce.
"""
import sys, itertools, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import paths, os
LR = np.log(1e300)          # in range: |L| <= LR
DEEP = 1e14
BIGL = 1e18
KEYS = 2.0 ** 34


def lits(S):
    out = set()
    for perm in itertools.permutations(S):
        t = "".join(perm)
        for i in range(len(t) + 1):
            out.add(float(int(t)) / 10 ** (len(t) - i))
    return np.log(np.array(sorted(out)))


def clipL(L):
    L = np.where(L > DEEP, BIGL, L)
    L = np.where(L < -DEEP, -BIGL, L)
    return L


def combine(La, Lb):
    """all Lemma-3 forms of magnitudes; returns L array (NaN = undefined)"""
    A = La[:, None]
    B = Lb[None, :]
    A, B = np.broadcast_arrays(A, B)
    A = A.ravel(); B = B.ravel()
    out = []
    with np.errstate(all="ignore"):
        out.append(A + B)                          # x*y
        out.append(A - B)                          # x/y
        out.append(B - A)                          # y/x
        # x+y, |x-y| via log-sum
        m = np.maximum(A, B); d = -np.abs(A - B)
        s_add = m + np.log1p(np.exp(d))
        s_sub = m + np.log(-np.expm1(d))           # -inf if equal
        both_deep = (np.abs(A) >= BIGL) & (np.abs(B) >= BIGL)
        s_add = np.where(both_deep & (A != B), np.nan, s_add)
        s_sub = np.where(both_deep, np.nan, s_sub)
        out.append(s_add); out.append(s_sub)
        # powers: x^(+-y): L = +- y * Lx ; y = exp(B)
        for (X, Y) in ((A, B), (B, A)):
            # magnitude of exponent * ln x
            lx = X
            ly = Y
            mag = ly + np.log(np.abs(lx))          # ln(|y|*|Lx|)
            val = np.where(mag > np.log(DEEP), np.sign(lx) * BIGL, np.sign(lx) * np.exp(np.minimum(mag, 50)))
            val = np.where(lx == 0, 0.0, val)      # 1^y = 1
            val = np.where(ly == -np.inf, 0.0, val)  # x^0 = 1
            zero_base = lx == -np.inf
            p1 = np.where(zero_base, np.where(ly > -np.inf, -np.inf, np.nan), val)   # 0^y (y>0) = 0
            p2 = np.where(zero_base, np.nan, -val)                                     # 0^-y undefined
            p2 = np.where((~zero_base) & (ly == -np.inf), np.nan, p2)                 # x^-0 same as x^0: skip dup
            out.append(p1); out.append(p2)
    L = np.concatenate(out)
    L = L[~np.isnan(L)]
    return clipL(L)


def key(L):
    k = np.where(np.isinf(L), np.where(L > 0, 2 ** 62, -2 ** 62), np.round(np.clip(L, -2e8, 2e8) * KEYS))
    return k.astype(np.int64)


def uniq(L):
    k, idx = np.unique(key(L), return_index=True)
    return L[idx]


def inrange(L):
    return (np.abs(L) <= LR) | (L == -np.inf)


R, A = {}, {}


def build(S):
    S = tuple(S)
    if S in A:
        return
    lr = lits(S)
    ra, aa = [lr], [lr]
    for r in range(1, len(S) // 2 + 1):
        for S1 in itertools.combinations(S, r):
            S2 = tuple(d for d in S if d not in S1)
            if 2 * r == len(S) and S1 > S2:
                continue
            build(S1); build(S2)
            x = combine(R[S1], R[S2]); ra.append(x[inrange(x)])
            aa.append(combine(A[S1], A[S2]))
    R[S] = uniq(np.concatenate(ra))
    A[S] = uniq(np.concatenate(aa))


if __name__ == "__main__":
    S = tuple(sys.argv[1])
    build(S)
    for T in sorted(A, key=lambda t: (len(t), t)):
        if len(T) < 3:
            continue
        a = A[T]; a_in = a[inrange(a)]
        new = a_in[~np.isin(key(a_in), key(R[T]))]
        # shadows: approx equal to a value of a proper subset (or 1, or 0)
        sub = [np.array([0.0, -np.inf])]
        for r in range(1, len(T)):
            for U in itertools.combinations(T, r):
                sub.append(R[U])
        sub = np.concatenate(sub)
        sh = np.isin(key(new), key(sub))
        print("".join(T), "R", len(R[T]), "A_in", len(a_in), "extreme", int((~inrange(a)).sum()),
              "new", len(new), "shadow", int(sh.sum()), "exotic", int((~sh).sum()))
        if len(T) == len(S):
            ex = np.sort(new[~sh])
            np.save(os.path.join(paths.C5, "exotic_") + "".join(T) + ".npy", ex)
            print("  exotic sample (value):", [f"{np.exp(v):.6g}" for v in ex[:: max(1, len(ex) // 20)]])
