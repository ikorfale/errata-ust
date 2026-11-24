"""3x3: the uniform spanning tree is determinantal, P(S subset of T) = det Y[S,S] with Y the transfer-current matrix
(Burton-Pemantle). Step 1 checks this on all 220 edge triples against the 192 trees. Step 2 builds a law over the 192
trees that matches EVERY single-edge and EVERY pair probability of the uniform law exactly (an LP) and pushes the
'claw' triples (three of the four edges at the centre vertex) as far as it can, both ways; step 3 samples it and
shows that single-edge and pair tests pass while the claw count fails. The 'U' triples (three sides of one unit
square) are reported as a forced class: no law that matches the pairs can move their total (the tree structure
fixes it), just as the 'corner' pairs are forced by the singles. Found by an LP range scan over all 220 triples:
164 can move, the claws furthest (uniform 1/6 each, range 0.104 to 0.208).
Usage: python3 check_triples_small.py k seed"""
import sys, json, itertools, numpy as np
from scipy import stats, optimize
from ust import all_trees, sampler
from pairs import transfer, pair_classes, exact_pairs
k, seed = int(sys.argv[1]), int(sys.argv[2])
L = 3; n, E, nb, C = pair_classes(L)
T = all_trees(n, E); M = np.zeros((len(T), len(E)))
for i, t in enumerate(T): M[i, list(t)] = 1
Y = transfer(n, E); R = np.diag(Y)
tri = list(itertools.combinations(range(len(E)), 3))
emp = np.array([M[:, list(s)].prod(1).mean() for s in tri]); ex = np.array([np.linalg.det(Y[np.ix_(s, s)]) for s in tri])
print(json.dumps(dict(check='det law vs enumeration, all triples', triples=len(tri), trees=len(T), max_abs_err=float(f'{np.abs(emp-ex).max():.2e}'))))
pos = {e: i for i, e in enumerate(E)}
U = []
for i in range(L - 1):
    for j in range(L - 1):
        a = i * L + j; sq = [pos[(a, a + 1)], pos[(a + L, a + L + 1)], pos[(a, a + L)], pos[(a + 1, a + 1 + L)]]
        U += list(itertools.combinations(sq, 3))
U = np.array(U); wU = M[:, U].prod(2).sum(1)                      # U triples per tree (at most 1 per square)
ctr = [k for _, k in nb[n // 2]]; K = np.array(list(itertools.combinations(ctr, 3))); wK = M[:, K].prod(2).sum(1)  # claws at the centre
pairs = list(itertools.combinations(range(len(E)), 2))
A = np.vstack([M.T, np.array([M[:, a] * M[:, b] for a, b in pairs]), np.ones(len(T))])
b = np.r_[R, [R[a] * R[b] - Y[a, b] ** 2 for a, b in pairs], 1]
uni = np.full(len(T), 1 / len(T)); assert np.abs(A @ uni - b).max() < 1e-9
laws = {}
for sense, sgn in (('max', -1), ('min', 1)):
    rU = [optimize.linprog(g * wU, A_eq=A, b_eq=b, bounds=(0, None), method='highs').fun * g for g in (-1, 1)]
    r = optimize.linprog(sgn * wK, A_eq=A, b_eq=b, bounds=(0, None), method='highs')
    q = np.clip(r.x, 0, None); q /= q.sum(); laws['lp_' + sense] = q
    print(json.dumps(dict(check=f'LP law matching all singles and pairs, {sense} claws', U_range_any_such_law=[round(x, 4) for x in rU], status=r.status, support=int((q > 1e-12).sum()),
                          rank_constraints=int(np.linalg.matrix_rank(A)), max_constraint_err=float(f'{np.abs(A@q-b).max():.1e}'),
                          claw_mean_uniform=round(wK.mean(), 4), claw_mean_lp=round(float(wK @ q), 4))))
samplers = {'wilson': sampler('wilson', n, E, nb)}
for name, q in laws.items(): samplers[name] = (lambda q: lambda r: T[r.choice(len(T), p=q)])(q)
for name, s in samplers.items():
    rng = np.random.default_rng(seed); X = np.zeros((k, len(E)), np.int8)
    for i in range(k): X[i, list(s(rng))] = 1
    z = (X.mean(0) - R) / np.sqrt(R * (1 - R) / k)
    out = dict(sampler=name, samples=k, seed=seed, single_edge_max_abs_z=round(float(np.abs(z).max()), 2))
    for c, P in C.items():
        per = (X[:, P[:, 0]] & X[:, P[:, 1]]).sum(1); mu = exact_pairs(Y, P).sum()
        out[c + '_t'] = round((per.mean() - mu) / (per.std(ddof=1) / np.sqrt(k)), 2)
    allp = np.array([(X[:, a] & X[:, b]).mean() for a, b in pairs]); exp2 = b[len(E):len(E) + len(pairs)]
    out['all_66_pairs_max_abs_z'] = round(float(np.abs((allp - exp2) / np.sqrt(exp2 * (1 - exp2) / k)).max()), 2)
    perU = X[:, U].prod(2).sum(1); muU = sum(np.linalg.det(Y[np.ix_(s, s)]) for s in U)
    out['U'] = dict(triples=len(U), exact=round(muU, 4), observed=round(perU.mean(), 4), t=round((perU.mean() - muU) / (perU.std(ddof=1) / np.sqrt(k)), 2))
    perK = X[:, K].prod(2).sum(1); muK = sum(np.linalg.det(Y[np.ix_(s, s)]) for s in K)
    out['claw'] = dict(triples=len(K), exact=round(muK, 4), observed=round(perK.mean(), 4), t=round((perK.mean() - muK) / (perK.std(ddof=1) / np.sqrt(k)), 2))
    print(json.dumps(out))
