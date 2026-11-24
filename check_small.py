"""Full chi-square on the 3x3 grid: all spanning trees enumerated (matrix-tree count must agree), k samples per sampler.
Usage: python3 check_small.py k seed"""
import sys, json, numpy as np
from scipy import stats
from ust import grid, all_trees, tree_count, sampler
k, seed = int(sys.argv[1]), int(sys.argv[2])
n, E, nb = grid(3)
T = all_trees(n, E); idx = {t: i for i, t in enumerate(T)}
assert len(T) == tree_count(n, E), (len(T), tree_count(n, E))
for name in ('wilson', 'aldous_broder', 'random_mst', 'random_dfs'):
    rng = np.random.default_rng(seed); s = sampler(name, n, E, nb)
    c = np.bincount([idx[s(rng)] for _ in range(k)], minlength=len(T))
    chi = stats.chisquare(c)
    print(json.dumps(dict(grid='3x3', trees=len(T), sampler=name, samples=k, seed=seed, chi2=round(chi.statistic, 1),
                          df=len(T) - 1, p=float(f'{chi.pvalue:.3g}'), max_over_min=round(c.max() / max(c.min(), 1), 2))))
