"""U triples past 3x3: for every unit square, s = number of its four sides in the tree (s <= 3, a tree has no cycle).
1{s=3} = 1 - s + C(s,2) - 1{s=0}, so the U triples are fixed by singles and pairs only where P(s=0) = 0 (squares touching
a degree-2 corner). Exact: P(U in T) = det Y[S,S]; P(s=0) = det(I - Y[S,S]) (determinantal law, complement form).
One observation per tree: number of squares with s=3 (interior squares only). Usage: python3 check_u_large.py L k seed sampler"""
import sys, json, numpy as np
from ust import sampler, grid
from pairs import transfer
L, k, seed, name = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
n, E, nb = grid(L); Y = transfer(n, E); pos = {e: i for i, e in enumerate(E)}
sq, cornered = [], []
for i in range(L - 1):
    for j in range(L - 1):
        a = i * L + j
        sq.append([pos[(a, a + 1)], pos[(a + L, a + L + 1)], pos[(a, a + L)], pos[(a + 1, a + 1 + L)]])
        cornered.append((i in (0, L - 2)) and (j in (0, L - 2)))
sq = np.array(sq); cornered = np.array(cornered)
p0 = np.array([np.linalg.det(np.eye(4) - Y[np.ix_(s, s)]) for s in sq])
p3 = np.array([sum(np.linalg.det(Y[np.ix_(t, t)]) for t in ([s[0], s[1], s[2]], [s[0], s[1], s[3]], [s[0], s[2], s[3]], [s[1], s[2], s[3]])) for s in sq])
inner = ~cornered
rng = np.random.default_rng(seed); smp = sampler(name, n, E, nb)
X = np.zeros((k, len(E)), np.int8)
for r in range(k): X[r, list(smp(rng))] = 1
S = X[:, sq].sum(2)                      # trees x squares
per3 = (S[:, inner] == 3).sum(1); per0 = (S[:, inner] == 0).sum(1)
t = lambda v, ex: round((v.mean() - ex) / (v.std(ddof=1) / np.sqrt(k)), 2)
print(json.dumps(dict(grid=f'{L}x{L}', sampler=name, samples=k, seed=seed, squares=int(inner.sum()),
    p0_corner_max=float(p0[cornered].max()), p0_inner_min=round(float(p0[inner].min()), 5), p0_inner_max=round(float(p0[inner].max()), 5),
    exact_U_count=round(float(p3[inner].sum()), 3), obs_U_count=round(float(per3.mean()), 3), t_U=t(per3, p3[inner].sum()),
    exact_s0_count=round(float(p0[inner].sum()), 3), obs_s0_count=round(float(per0.mean()), 3), t_s0=t(per0, p0[inner].sum()))))
