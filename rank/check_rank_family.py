#!/usr/bin/env python3
"""Rank of the single+pair statistics of the uniform spanning tree, across grids.

Question (board thread on errata-ust, 25.11): the linear identities among the m = 1 + E + C(E,2)
statistics {1, x_e, x_e x_f} over all spanning trees are claimed to be exactly
  Foster (1) + edge-pair (E: sum_{e!=f} x_ef = (V-2) x_f) + one per 2-edge cut {a,b}: 1 - x_a - x_b + x_ab = 0,
so rank(M) = min(T, m - rank(K_fam)).  Checked so far on 2x3 and 4x4 only.

Exact certificate per grid, all arithmetic modular (two primes):
  lower bound  rank_Q(M) >= rank_p(G), G = M M^T (a nonzero minor mod p is nonzero over Z);
  upper bound  rank_Q(M) <= m - rank_Q(K) <= m - rank_p(K), every row of K verified K@M == 0 over Z.
If they meet, rank_Q(M) is exact and the family K is complete.
"""
import itertools, sys, json
import numpy as np

PRIMES = [2147483647, 1000003]

def grid(r, c):
    V = r * c; E = []
    for i in range(r):
        for j in range(c):
            v = i * c + j
            if j + 1 < c: E.append((v, v + 1))
            if i + 1 < r: E.append((v, v + c))
    return V, E

def trees(V, E):
    """All spanning trees by edge-include/exclude backtracking with a small union-find copy."""
    out = []; n = len(E)
    def find(p, x):
        while p[x] != x: x = p[x]
        return x
    def connected_possible(k, chosen):
        # remaining edges k.. plus chosen must connect the graph
        p = list(range(V))
        for e in chosen + list(range(k, n)):
            a, b = find(p, E[e][0]), find(p, E[e][1])
            if a != b: p[a] = b
        return len({find(p, v) for v in range(V)}) == 1
    def rec(k, chosen, p):
        if len(chosen) == V - 1:
            out.append(tuple(chosen)); return
        if k == n or n - k < V - 1 - len(chosen): return
        a, b = find(p, E[k][0]), find(p, E[k][1])
        if a != b:
            q = p[:]; q[a] = b
            rec(k + 1, chosen + [k], q)
        if connected_possible(k + 1, chosen):
            rec(k + 1, chosen, p)
    rec(0, [], list(range(V)))
    return out

def rank_mod(A, p):
    A = (np.array(A, dtype=object) % p).astype(np.int64)
    A = A.copy(); rows, cols = A.shape; r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c] % p: piv = i; break
        if piv is None: continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        nz = np.nonzero(A[:, c])[0]
        for i in nz:
            if i != r:
                A[i] = (A[i] - A[i, c] * A[r]) % p
        r += 1
        if r == rows: break
    return r

def two_cuts(V, E):
    cuts = []
    for a, b in itertools.combinations(range(len(E)), 2):
        p = list(range(V))
        def find(x):
            while p[x] != x: x = p[x]
            return x
        for k, (u, v) in enumerate(E):
            if k in (a, b): continue
            x, y = find(u), find(v)
            if x != y: p[x] = y
        if len({find(v) for v in range(V)}) > 1: cuts.append((a, b))
    return cuts

def run(r, c):
    V, E = grid(r, c); nE = len(E)
    T = trees(V, E)
    pairs = list(itertools.combinations(range(nE), 2)); pidx = {pq: k for k, pq in enumerate(pairs)}
    m = 1 + nE + len(pairs)
    X = np.zeros((len(T), nE), dtype=np.int64)
    for t, tr in enumerate(T): X[t, list(tr)] = 1
    M = np.empty((m, len(T)), dtype=np.int64)
    M[0] = 1; M[1:1 + nE] = X.T
    for k, (a, b) in enumerate(pairs): M[1 + nE + k] = X[:, a] * X[:, b]
    G = M @ M.T  # entries <= T, exact in int64
    def P(a, b): return 1 + nE + pidx[(min(a, b), max(a, b))]
    K = []
    row = np.zeros(m, dtype=np.int64); row[1:1 + nE] = 1; row[0] = -(V - 1); K.append(row)
    for f in range(nE):
        row = np.zeros(m, dtype=np.int64); row[1 + f] = -(V - 2)
        for e in range(nE):
            if e != f: row[P(e, f)] += 1
        K.append(row)
    cuts = two_cuts(V, E)
    for a, b in cuts:
        row = np.zeros(m, dtype=np.int64); row[0] = 1; row[1 + a] = -1; row[1 + b] = -1; row[P(a, b)] = 1
        K.append(row)
    K = np.array(K)
    viol = int(np.abs(K @ M).max())
    lo = [rank_mod(G, p) for p in PRIMES]
    rk = [rank_mod(K, p) for p in PRIMES]
    up = m - max(rk)  # rank_p(K) <= rank_Q(K); use the largest as the tightest valid bound
    res = dict(grid=f"{r}x{c}", V=V, E=nE, T=len(T), m=m, two_cuts=len(cuts),
               family_rows=len(K), family_rank=rk, K_violations=viol,
               rank_lower=lo, rank_upper=up, exact=(min(lo) == max(lo) == up),
               predicted=min(len(T), m - (1 + nE + len(cuts))))
    print(json.dumps(res), flush=True)
    return res

if __name__ == "__main__":
    for spec in sys.argv[1:]:
        r, c = map(int, spec.split("x")); run(r, c)
