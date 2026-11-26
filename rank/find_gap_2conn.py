#!/usr/bin/env python3
"""The one seed-1 gap left by family v5: a 2-vertex-connected graph (V=8, E=13). Smallest edge sets S carrying a
quadratic identity not in span(v5), with the identity printed over Q. Screen as in find_new_identity3.py."""
import sys, itertools, numpy as np, sympy as sp, networkx as nx
sys.argv += ['--v5']
sys.path.insert(0, '.')
import check_rank_cutside as cc
from check_rank_family import rank_mod
E0 = [(0, 5), (0, 7), (1, 2), (1, 3), (1, 5), (1, 6), (1, 7), (2, 3), (2, 4), (3, 7), (4, 7), (5, 6), (6, 7)]
G0 = nx.Graph(E0); print('2-vertex-connected:', nx.is_biconnected(G0), 'cut vertices:', list(nx.articulation_points(G0)))
print('degrees:', dict(G0.degree()))
cc.grid = lambda r, c: (8, list(E0)); cc.BUILD_ONLY = True
E, T, M, K, pairs, pi, m = cc.run('1x1'); nE = len(E); p = 1000003
rK = rank_mod(K, p); G = M @ M.T; hits = []
print('m', m, 'rank M', rank_mod(G, p), 'family upper', m - rK)
for k in range(2, nE + 1):
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
