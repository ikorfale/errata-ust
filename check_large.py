"""Large grids, past enumeration: Kirchhoff says P(e in T) = R_eff(e). One observation per tree: for an edge class C
(corner edges, other boundary edges, interior), count |T ∩ C| per tree; its exact mean is sum_{e in C} R_eff(e);
a t-statistic uses the empirical sd across trees (edges inside one tree are dependent). Also max |z| over single
edges (binomial sd, Bonferroni). The total over all edges is n-1 for every tree (Foster): forced, so not reported as a test.
Usage: python3 check_large.py L k seed sampler"""
import sys, json, numpy as np
from scipy import stats
from ust import grid, resistances, sampler
L, k, seed, name = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
n, E, nb = grid(L); R = resistances(n, E)
assert abs(R.sum() - (n - 1)) < 1e-6
deg = np.array([len(x) for x in nb])
cls = np.array(['corner' if min(deg[a], deg[b]) == 2 else 'boundary' if (deg[a] < 4 and deg[b] < 4) else 'interior' for a, b in E])
rng = np.random.default_rng(seed); s = sampler(name, n, E, nb)
X = np.zeros((k, len(E)), np.int8)
for i in range(k): X[i, list(s(rng))] = 1
out = dict(grid=f'{L}x{L}', sampler=name, samples=k, seed=seed)
for c in ('corner', 'boundary', 'interior'):
    m = cls == c; per = X[:, m].sum(1)
    t = (per.mean() - R[m].sum()) / (per.std(ddof=1) / np.sqrt(k))
    out[c] = dict(edges=int(m.sum()), exact=round(R[m].sum() / m.sum(), 4), observed=round(per.mean() / m.sum(), 4), t=round(t, 2))
f = X.mean(0); z = (f - R) / np.sqrt(R * (1 - R) / k)
out['max_abs_edge_z'] = round(float(np.abs(z).max()), 2); out['bonferroni_p'] = float(f'{min(1, len(E) * 2 * stats.norm.sf(np.abs(z).max())):.3g}')
print(json.dumps(out))
