"""Does networkx.random_spanning_tree (weight=None) sample uniformly? 3x3 grid, 192 trees, chi-square over all trees,
plus single-edge and pair laws. Weighted modes: edge weights 1,2,3 assigned by a fixed rule; mult -> P(T) ~ product, add -> P(T) ~ sum.
Usage: python3 nx/check_nx_small.py k seed [uniform|mult|add]"""
import sys, json, itertools, time, numpy as np, networkx as nx
from scipy.stats import chi2
k, seed = int(sys.argv[1]), int(sys.argv[2]); mode = sys.argv[3] if len(sys.argv) > 3 else 'uniform'
G = nx.grid_2d_graph(3, 3); E = sorted(tuple(sorted(e)) for e in G.edges())
trees = [frozenset(c) for c in itertools.combinations(E, 8) if nx.is_tree(nx.Graph(list(c)))]
W = {e: 1 + (i % 3) for i, e in enumerate(E)}
nx.set_edge_attributes(G, {e: W[e] for e in E}, 'w')
w = np.array([1.0 if mode == 'uniform' else (np.prod([W[e] for e in t]) if mode == 'mult' else sum(W[e] for e in t)) for t in trees], float); w /= w.sum()
idx = {t: i for i, t in enumerate(trees)}; cnt = np.zeros(len(trees)); t0 = time.time()
import random; rng = random.Random(seed)
for _ in range(k):
    T = nx.random_spanning_tree(G, None if mode == 'uniform' else 'w', multiplicative=(mode != 'add'), seed=rng); cnt[idx[frozenset(tuple(sorted(e)) for e in T.edges())]] += 1
ex = k * w; x2 = float(((cnt - ex) ** 2 / ex).sum())
print(json.dumps(dict(mode=mode, networkx=nx.__version__, trees=len(trees), samples=k, seed=seed, chi2=round(x2, 1), df=len(trees) - 1,
    p=round(float(chi2.sf(x2, len(trees) - 1)), 4), worst=[(round(float(cnt[i]),0), round(float(ex[i]),1)) for i in np.argsort(-(cnt - ex) ** 2 / ex)[:3]], secs=round(time.time() - t0, 1))))
