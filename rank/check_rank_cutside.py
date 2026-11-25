#!/usr/bin/env python3
"""Family v2 = Foster + edge-pair + 2-cut + cut-side rows, where for each 2-edge cut {a,b}, each side A (|A|>1)
and each cut edge a:  (1 - x_a) * (|A| - 1 - sum_{e in E(A)} x_e) = 0  (if a is absent, b is A's only exit, so the
tree restricted to A is connected). Certificate as in check_rank_family.py: rank_p(MM^T) <= rank_Q(M) <= m - rank_p(K).
Family v3 adds cut-rung rows: if edge r joins the endpoints u,w of the two cut edges inside A, then x_r * (|A|-1 - e_T(A)) = 0
(a disconnected T[A] has exactly two parts, one holding u and one holding w; r would join them).
Family v4 (--v4) adds rows for every vertex set A whose exit edges touch at most two vertices of A.
Graphs: grids RxC and ladders with extra chords given as 'RxC+u-v,...'."""
import sys, itertools, json, numpy as np
sys.path.insert(0, '.')
from check_rank_family import grid, trees, two_cuts, rank_mod, PRIMES

def side(V, E, a, b):
    p = list(range(V))
    def f(x):
        while p[x] != x: x = p[x]
        return x
    for k, (u, v) in enumerate(E):
        if k in (a, b): continue
        x, y = f(u), f(v)
        if x != y: p[x] = y
    r0 = f(0); A = {v for v in range(V) if f(v) == r0}
    return [A, set(range(V)) - A]

def run(spec):
    base, _, extra = spec.partition('+')
    r, c = map(int, base.split('x')); V, E = grid(r, c)
    for ch in filter(None, extra.split(',')):
        u, v = map(int, ch.split('-')); E.append((u, v))
    nE = len(E); T = trees(V, E)
    pairs = list(itertools.combinations(range(nE), 2)); P = {p: 1 + nE + k for k, p in enumerate(pairs)}
    pi = lambda a, b: P[(min(a, b), max(a, b))]; m = 1 + nE + len(pairs)
    X = np.zeros((len(T), nE), dtype=np.int64)
    for t, tr in enumerate(T): X[t, list(tr)] = 1
    M = np.empty((m, len(T)), dtype=np.int64); M[0] = 1; M[1:1 + nE] = X.T
    for k, (a, b) in enumerate(pairs): M[1 + nE + k] = X[:, a] * X[:, b]
    K = []
    row = np.zeros(m, np.int64); row[0] = -(V - 1); row[1:1 + nE] = 1; K.append(row)
    for f in range(nE):
        row = np.zeros(m, np.int64); row[1 + f] = -(V - 2)
        for e in range(nE):
            if e != f: row[pi(e, f)] += 1
        K.append(row)
    cuts = two_cuts(V, E)
    for a, b in cuts:
        row = np.zeros(m, np.int64); row[0] = 1; row[1 + a] = -1; row[1 + b] = -1; row[pi(a, b)] = 1; K.append(row)
    n_old = len(K); n_side = 0
    for a, b in cuts:
        for A in side(V, E, a, b):
            if len(A) < 2: continue
            EA = [k for k, (u, v) in enumerate(E) if u in A and v in A]
            for x in (a, b):   # (1 - x)(|A|-1 - sum EA) = (|A|-1) - (|A|-1) x - sum EA + sum x*EA
                row = np.zeros(m, np.int64); row[0] = len(A) - 1; row[1 + x] -= len(A) - 1
                for e in EA: row[1 + e] -= 1; row[pi(x, e)] += 1
                K.append(row); n_side += 1
    n_chain = 0
    if not NO_RUNG:
        for a, b in cuts:   # cut-rung rows: r joins the two cut endpoints inside side A  =>  x_r (|A|-1 - e_T(A)) = 0
            for A in side(V, E, a, b):
                if len(A) < 2: continue
                u = [w for w in E[a] if w in A][0]; w = [z for z in E[b] if z in A][0]
                if u == w: continue
                EA = [k for k, (p_, q_) in enumerate(E) if p_ in A and q_ in A]
                for rr in EA:
                    if set(E[rr]) == {u, w}:
                        row = np.zeros(m, np.int64); row[1 + rr] = len(A) - 1
                        for e in EA: row[pi(rr, e) if e != rr else 1 + rr] -= 1
                        K.append(row); n_chain += 1
    n_term = 0
    if TWO_TERMINAL:   # v4: every vertex set A whose exits touch at most two vertices of A (2 <= |A| <= V-1)
        adj_exit = {}
        for A_bits in range(1, 2 ** V - 1):
            A = {v for v in range(V) if A_bits >> v & 1}
            if len(A) < 2: continue
            ex = [k for k, (p_, q_) in enumerate(E) if (p_ in A) != (q_ in A)]
            term = {}
            for k in ex:
                t_ = E[k][0] if E[k][0] in A else E[k][1]; term.setdefault(t_, []).append(k)
            if len(term) > 2: continue
            EA = [k for k, (p_, q_) in enumerate(E) if p_ in A and q_ in A]
            def conn_row(coef_idx=None, minus=None):
                row = np.zeros(m, np.int64)
                if minus is None:   # x_r * (|A|-1 - e_T(A))
                    row[1 + coef_idx] = len(A) - 1
                    for e in EA: row[pi(coef_idx, e) if e != coef_idx else 1 + coef_idx] -= 1
                else:               # (1 - x_a) * (|A|-1 - e_T(A))
                    row[0] = len(A) - 1; row[1 + minus] -= len(A) - 1
                    for e in EA: row[1 + e] -= 1; row[pi(minus, e)] += 1
                return row
            if len(term) == 1:      # one terminal: T[A] always connected (linear rows)
                row = np.zeros(m, np.int64); row[0] = len(A) - 1
                for e in EA: row[1 + e] -= 1
                K.append(row); n_term += 1; continue
            (u, eu), (w, ew) = term.items()
            for exits in (eu, ew):
                if len(exits) == 1: K.append(conn_row(minus=exits[0])); n_term += 1
            for rr in EA:
                if set(E[rr]) == {u, w}: K.append(conn_row(coef_idx=rr)); n_term += 1
    if BUILD_ONLY: return E, T, M, np.array(K), pairs, pi, m
    K = np.array(K); viol = int(np.abs(K @ M).max())
    G = M @ M.T
    lo = [rank_mod(G, p) for p in PRIMES]
    old_up = m - max(rank_mod(K[:n_old], p) for p in PRIMES)
    up = m - max(rank_mod(K, p) for p in PRIMES)
    res = dict(graph=spec, V=V, E=nE, T=len(T), m=m, cuts=len(cuts), side_rows=n_side, cut_rung_rows=n_chain, two_terminal_rows=n_term, K_violations=viol,
               rank_lower=lo, upper_old_family=old_up, upper_new_family=up, exact=min(lo) == max(lo) == up)
    print(json.dumps(res), flush=True)

BUILD_ONLY = False
NO_RUNG = False
TWO_TERMINAL = '--v4' in sys.argv
if TWO_TERMINAL: sys.argv.remove('--v4')
if __name__ == '__main__':
    for s in sys.argv[1:]: run(s)
