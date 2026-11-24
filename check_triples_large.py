"""Claw triples past enumeration: at every interior vertex (degree 4) of an L x L grid, the four 'claws' (three of its
four edges). Exact P(claw in T) = det Y[S,S] (determinantal law, Burton-Pemantle). One observation per tree: the number
of claws present; t uses the empirical sd across trees. Usage: python3 check_triples_large.py L k seed sampler"""
import sys, json, itertools, numpy as np
from ust import sampler
from pairs import transfer, pair_classes
L, k, seed, name = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
n, E, nb, C = pair_classes(L); Y = transfer(n, E)
K = np.array([c for v in range(n) if len(nb[v]) == 4 for c in itertools.combinations([e for _, e in nb[v]], 3)])
ex = np.array([np.linalg.det(Y[np.ix_(s, s)]) for s in K]); R = np.diag(Y)
rng = np.random.default_rng(seed); s = sampler(name, n, E, nb)
X = np.zeros((k, len(E)), np.int8)
for i in range(k): X[i, list(s(rng))] = 1
per = X[:, K].prod(2).sum(1)
out = dict(grid=f'{L}x{L}', sampler=name, samples=k, seed=seed, claws=len(K), exact_per_claw=round(ex.mean(), 4),
           if_independent=round(float(R[K].prod(1).mean()), 4), observed_per_claw=round(per.mean() / len(K), 4),
           count_mean=float(per.mean()), count_sd=float(per.std(ddof=1)), exact_count=float(ex.sum()),
           t=round((per.mean() - ex.sum()) / (per.std(ddof=1) / np.sqrt(k)), 2))
print(json.dumps(out))
