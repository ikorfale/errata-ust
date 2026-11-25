#!/usr/bin/env python3
"""Fast screen (mod p): edge sets S with a quadratic identity supported on S not in span(family v3).
new identity on S  <=>  |idx| - rank(M[idx]) > rank(K) - rank(K[:, not idx])."""
import sys, itertools, numpy as np, sympy as sp
sys.path.insert(0, '.')
import check_rank_cutside as cc
from check_rank_family import rank_mod
cc.BUILD_ONLY = True
spec, kmin, kmax = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); p = 1000003
E, T, M, K, pairs, pi, m = cc.run(spec); nE = len(E)
rK = rank_mod(K, p); G = M @ M.T; hits = []
for k in range(kmin, kmax + 1):
    for S in itertools.combinations(range(nE), k):
        if any(set(h) <= set(S) for h in hits): continue
        idx = [0] + [1 + e for e in S] + [pi(a, b) for a, b in itertools.combinations(S, 2)]
        rest = [i for i in range(m) if i not in set(idx)]
        if len(idx) - rank_mod(G[np.ix_(idx, idx)], p) > rK - rank_mod(K[:, rest], p):
            hits.append(S); print(k, [E[e] for e in S], flush=True)
            for v in sp.Matrix(M[idx]).T.nullspace():
                v = v * sp.ilcm(*[x.q for x in v])
                lab = ['1'] + [f'x{E[e]}' for e in S] + [f'x{E[a]}x{E[b]}' for a, b in itertools.combinations(S, 2)]
                print('   ', ' '.join(f'{int(x):+d}*{l}' for x, l in zip(v, lab) if x != 0), flush=True)
    if hits: break
print('done', len(hits))
