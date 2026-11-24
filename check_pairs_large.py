"""Pair-class test past enumeration (see pairs.py). The 'corner' class is forced by the single-edge marginals
(one of the two corner edges is always present, so P(both) = 2R - 1): it is reported but tests nothing new.
Covariance check (marginals removed): mean over the class of observed P(ef)-P(e)P(f) vs exact -Y_ef^2, sd by 200 bootstrap resamples of trees.
Usage: python3 check_pairs_large.py L k seed sampler"""
import sys, json, numpy as np
from ust import sampler
from pairs import transfer, pair_classes, exact_pairs
L, k, seed, name = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
n, E, nb, C = pair_classes(L); Y = transfer(n, E)
assert abs(np.trace(Y) - (n - 1)) < 1e-6
rng = np.random.default_rng(seed); s = sampler(name, n, E, nb)
X = np.zeros((k, len(E)), np.int8)
for i in range(k): X[i, list(s(rng))] = 1
out = dict(grid=f'{L}x{L}', sampler=name, samples=k, seed=seed)
for c, P in C.items():
    p = exact_pairs(Y, P); per = (X[:, P[:, 0]] & X[:, P[:, 1]]).sum(1)
    indep = (np.diag(Y)[P[:, 0]] * np.diag(Y)[P[:, 1]]).sum()
    Xa, Xb = X[:, P[:, 0]].astype(float), X[:, P[:, 1]].astype(float); Z = Xa * Xb
    cov = lambda w: ((w @ Z) / k - (w @ Xa) * (w @ Xb) / k ** 2).mean()   # w = multinomial resampling counts
    bs = np.random.default_rng(0); boot = [cov(bs.multinomial(k, np.full(k, 1 / k)).astype(float)) for _ in range(200)]
    exc = -(Y[P[:, 0], P[:, 1]] ** 2).mean(); obc = cov(np.ones(k))
    out[c] = dict(cov_exact=round(exc, 4), cov_observed=round(obc, 4), cov_z=round((obc - exc) / np.std(boot), 2),
                  pairs=len(P), exact=round(p.sum() / len(P), 4), if_independent=round(indep / len(P), 4),
                  observed=round(per.mean() / len(P), 4), t=round((per.mean() - p.sum()) / (per.std(ddof=1) / np.sqrt(k)), 2))
print(json.dumps(out))
