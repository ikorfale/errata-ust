"""3x3: pair law checked against all 192 trees (exact), then an adversarial sampler that matches every single-edge
marginal exactly but gets pairs wrong (an LP over the 192 trees), to show what single-edge tests miss.
Usage: python3 check_pairs_small.py k seed"""
import sys, json, numpy as np
from fractions import Fraction
from scipy import stats, optimize
from ust import all_trees, resistances, sampler
from pairs import transfer, pair_classes, exact_pairs
k, seed = int(sys.argv[1]), int(sys.argv[2])
n, E, nb, C = pair_classes(3)
T = all_trees(n, E); M = np.zeros((len(T), len(E)))
for i, t in enumerate(T): M[i, list(t)] = 1
Y = transfer(n, E); R = np.diag(Y)
# 1. exact pair law against enumeration, every pair of edges
emp = M.T @ M / len(T); ex = np.outer(R, R) - Y ** 2; np.fill_diagonal(ex, R)
print(json.dumps(dict(check='transfer-current vs enumeration, all pairs', trees=len(T), max_abs_err=float(f'{np.abs(emp-ex).max():.2e}'),
                      corner_pair=str(Fraction(round(emp[C['corner'][0][0], C['corner'][0][1]] * 192), 192)))))
# 2. adversarial law: same marginals as uniform, maximise the number of 'perp' pairs both present
w = np.zeros(len(T))
for e, f in C['perp']: w += M[:, e] * M[:, f]
res = optimize.linprog(-w, A_eq=np.vstack([M.T, np.ones(len(T))]), b_eq=np.r_[R, 1], bounds=(0, None), method='highs')
q = np.clip(res.x, 0, None); q /= q.sum()
print(json.dumps(dict(check='adversarial LP law', support=int((q > 1e-12).sum()), max_marginal_err=float(f'{np.abs(M.T@q-R).max():.1e}'),
                      perp_mean_uniform=round(w.mean(), 4), perp_mean_adversarial=round(float(w @ q), 4))))
# 3. sample: uniform (wilson), adversarial, random MST; single-edge max |z| vs pair-class t
rng0 = np.random.default_rng(seed)
samplers = {'wilson': sampler('wilson', n, E, nb), 'random_mst': sampler('random_mst', n, E, nb),
            'adversarial': lambda r: T[r.choice(len(T), p=q)]}
for name, s in samplers.items():
    rng = np.random.default_rng(seed); X = np.zeros((k, len(E)), np.int8)
    for i in range(k): X[i, list(s(rng))] = 1
    z = (X.mean(0) - R) / np.sqrt(R * (1 - R) / k)
    out = dict(sampler=name, samples=k, seed=seed, single_edge_max_abs_z=round(float(np.abs(z).max()), 2),
               single_edge_bonferroni_p=float(f'{min(1, len(E)*2*stats.norm.sf(np.abs(z).max())):.3g}'))
    for c, P in C.items():
        per = (X[:, P[:, 0]] & X[:, P[:, 1]]).sum(1); mu = exact_pairs(Y, P).sum()
        out[c] = dict(pairs=len(P), exact=round(mu, 4), observed=round(per.mean(), 4), t=round((per.mean() - mu) / (per.std(ddof=1) / np.sqrt(k)), 2))
    print(json.dumps(out))
