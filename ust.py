"""Uniform spanning trees of a grid graph: samplers and exact checks.

Samplers
- wilson: Wilson's algorithm (loop-erased random walks), exactly uniform.
- aldous_broder: first-entrance edges of a random walk, exactly uniform (slow; a second positive control).
- random_mst: Kruskal on i.i.d. uniform edge weights. A spanning tree "at random", but NOT uniform (negative control).
- random_dfs: depth-first search with shuffled neighbours, also not uniform (negative control).

Exact laws
- Kirchhoff: P(e in T) = effective resistance across e (unit weights), from the Laplacian pseudoinverse, at any size.
- Foster: the sum of those resistances is n-1, so the mean number of edges is forced and tests nothing; the per-edge
  profile is not forced.
- Small grids: every spanning tree enumerated, a full chi-square against the uniform law.
"""
import numpy as np, itertools


def grid(L, M=None):
    M = M or L
    idx = lambda i, j: i * M + j
    E = [(idx(i, j), idx(i, j + 1)) for i in range(L) for j in range(M - 1)] + \
        [(idx(i, j), idx(i + 1, j)) for i in range(L - 1) for j in range(M)]
    n = L * M
    nb = [[] for _ in range(n)]
    for k, (a, b) in enumerate(E):
        nb[a].append((b, k)); nb[b].append((a, k))
    return n, E, nb


def wilson(n, nb, rng):
    intree = np.zeros(n, bool); nxt = -np.ones(n, int); nxte = -np.ones(n, int)
    intree[rng.integers(n)] = True
    for s in rng.permutation(n):
        u = s
        while not intree[u]:          # random walk, remembering the last exit (this erases loops)
            v, k = nb[u][rng.integers(len(nb[u]))]
            nxt[u], nxte[u] = v, k; u = v
        u = s
        while not intree[u]:
            intree[u] = True; u = nxt[u]
    return _wilson_edges(n, nxt, nxte, intree)


def _wilson_edges(n, nxt, nxte, intree):
    # the root keeps nxte = -1 only if it was never walked from; every other vertex's final exit edge is a tree edge
    return frozenset(int(k) for k in nxte if k >= 0)


def aldous_broder(n, nb, rng):
    seen = np.zeros(n, bool); u = rng.integers(n); seen[u] = True; left = n - 1; T = []
    while left:
        v, k = nb[u][rng.integers(len(nb[u]))]
        if not seen[v]: seen[v] = True; left -= 1; T.append(k)
        u = v
    return frozenset(T)


def random_mst(n, E, rng):
    par = list(range(n))
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    T = []
    for k in np.argsort(rng.random(len(E))):
        a, b = f(E[k][0]), f(E[k][1])
        if a != b: par[a] = b; T.append(int(k))
    return frozenset(T)


def random_dfs(n, nb, rng):
    seen = np.zeros(n, bool); s = int(rng.integers(n)); seen[s] = True; T = []; stack = [s]
    while stack:
        u = stack[-1]; opts = [(v, k) for v, k in nb[u] if not seen[v]]
        if not opts: stack.pop(); continue
        v, k = opts[rng.integers(len(opts))]; seen[v] = True; T.append(k); stack.append(v)
    return frozenset(T)


def resistances(n, E):
    Lap = np.zeros((n, n))
    for a, b in E:
        Lap[a, a] += 1; Lap[b, b] += 1; Lap[a, b] -= 1; Lap[b, a] -= 1
    P = np.linalg.pinv(Lap)
    return np.array([P[a, a] + P[b, b] - 2 * P[a, b] for a, b in E])


def all_trees(n, E):
    out = []
    for S in itertools.combinations(range(len(E)), n - 1):
        par = list(range(n)); ok = True
        def f(x):
            while par[x] != x: x = par[x]
            return x
        for k in S:
            a, b = f(E[k][0]), f(E[k][1])
            if a == b: ok = False; break
            par[a] = b
        if ok: out.append(frozenset(S))
    return out


def tree_count(n, E):
    Lap = np.zeros((n, n))
    for a, b in E:
        Lap[a, a] += 1; Lap[b, b] += 1; Lap[a, b] -= 1; Lap[b, a] -= 1
    return round(np.linalg.det(Lap[1:, 1:]))


def sampler(name, n, E, nb):
    return {'wilson': lambda r: wilson(n, nb, r), 'aldous_broder': lambda r: aldous_broder(n, nb, r),
            'random_mst': lambda r: random_mst(n, E, r), 'random_dfs': lambda r: random_dfs(n, nb, r)}[name]
