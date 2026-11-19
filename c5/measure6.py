r"""
C5 measurement for one 6-digit set (any order), approximate float log-domain.
New(S) = in-range values attainable only with extreme intermediates.  Only combinations where at
least one operand is new/extreme can produce new values, so we combine:
   NE(S1) x A(S2)  and  A(S1) x NE(S2),   NE = A \ R  (new or extreme)
and remove values (approximately) in R(S).
"""
import sys, itertools, numpy as np
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__))); import paths, os
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from measure import lits, combine, key, uniq, inrange, build, R, A

S = tuple(sys.argv[1])
consts = {"pi": np.log(np.pi), "gamma": np.log(0.5772156649015329), "delta": np.log(4.66920160910299),
          "e": 1.0}
for r in range(1, len(S)):
    for T in itertools.combinations(S, r):
        build(T)
NE = {}
for T in A:
    a = A[T]
    NE[T] = a[~np.isin(key(a), key(R[T]))]
cand = []
splits = []
for r in range(1, len(S) // 2 + 1):
    for S1 in itertools.combinations(S, r):
        S2 = tuple(d for d in S if d not in S1)
        if 2 * r == len(S) and S1 > S2:
            continue
        splits.append((S1, S2))
        for (P, Q) in ((NE[S1], A[S2]), (A[S1], NE[S2])):
            if len(P) and len(Q):
                for i in range(0, len(P), 20000):
                    y = combine(P[i:i + 20000], Q)
                    cand.append(uniq(y[inrange(y)]))
cand = uniq(np.concatenate(cand))
ck = np.sort(key(cand))
hit = np.zeros(len(ck), bool)
L0 = lits(S)
hit |= np.isin(ck, key(L0))
for S1, S2 in splits:
    P, Q = R[S1], R[S2]
    if len(P) > len(Q):
        P, Q = Q, P
    for i in range(0, len(P), 2000):
        y = combine(P[i:i + 2000], Q)
        hit |= np.isin(ck, key(y[inrange(y)]))
print("candidates", len(ck), "in R(S)", int(hit.sum()), flush=True)
new_keys = ck[~hit]
new = cand[np.isin(key(cand), new_keys)]
Rk = None
# shadows: approx equal to an in-range value of a subset with >= 3 digits fewer, or to 1 / 0
sub = [np.array([0.0, -np.inf])]
for r in range(1, len(S) - 2):
    for U in itertools.combinations(S, r):
        sub.append(R[U])
sh = np.isin(key(new), key(np.concatenate(sub)))
exo = np.sort(new[~sh])
print("".join(S), "new", len(new), "shadow", int(sh.sum()), "exotic", len(exo))
v = exo / np.log(10)
print("histogram log10:", np.histogram(v, bins=[-310, -100, -30, -10, -3, -1, 0, 1, 3, 10, 30, 100, 310])[0])
for c, lc in consts.items():
    i = np.searchsorted(exo, lc)
    near = exo[max(0, i - 3): i + 3]
    print(c, "closest exotic values:", [f"{np.exp(t):.12g} (err {abs(np.exp(t) - np.exp(lc)):.2e})" for t in near])
np.save(os.path.join(paths.C5, "exotic6_") + "".join(S) + ".npy", exo)
