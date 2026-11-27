"""Matroid version of check_sep_span.py (board #130777 claims the theorem for every matroid):
degree-<=2 identities on basis indicators == span of Foster q_E*{1,x_f}, 1-separation rows q_X*{1,x_f},
2-separation rows q_X(q_X+1). Matroids from a rank oracle: uniform, GF(p) matrices (Fano, non-Fano via GF(3),
random), Vamos, and random sparse paving matroids (circuit-hyperplanes pairwise meeting in <= r-2), most of
which are not representable. Exact mod P. 'no2' drops the 2-separation rows (control: gaps must appear)."""
import itertools, random, sys
from functools import lru_cache
P = 1000003
NO2 = 'no2' in sys.argv
def rank_mod(rows):
    basis = {}
    for r in rows:
        r = [v % P for v in r]
        for c, br in basis.items():
            if r[c]:
                t = r[c]; r = [(x - t * y) % P for x, y in zip(r, br)]
        piv = next((c for c, v in enumerate(r) if v), None)
        if piv is None: continue
        inv = pow(r[piv], P - 2, P); r = [v * inv % P for v in r]
        for c in basis:
            br = basis[c]
            if br[piv]:
                t = br[piv]; basis[c] = [(x - t * y) % P for x, y in zip(br, r)]
        basis[piv] = r
    return len(basis)
def gf_rank(cols, p):
    rows = [list(c) for c in cols]  # rank of vectors over GF(p)
    rows = [r[:] for r in rows]; rk = 0; n = len(rows[0]) if rows else 0
    for c in range(n):
        piv = next((i for i in range(rk, len(rows)) if rows[i][c] % p), None)
        if piv is None: continue
        rows[rk], rows[piv] = rows[piv], rows[rk]; inv = pow(rows[rk][c], p - 2, p)
        rows[rk] = [v * inv % p for v in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][c] % p:
                t = rows[i][c]; rows[i] = [(a - t * b) % p for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk
def vector_matroid(cols, p):
    @lru_cache(None)
    def rk(S): return gf_rank([cols[i] for i in S], p) if S else 0
    return len(cols), lambda S: rk(tuple(sorted(S)))
def uniform(r, n): return n, lambda S: min(len(set(S)), r)
def sparse_paving(n, r, H):  # H: list of frozensets of size r, pairwise |A&B| <= r-2
    H = [frozenset(h) for h in H]
    def rk(S):
        S = frozenset(S)
        if len(S) < r: return len(S)
        if len(S) == r and S in H: return r - 1
        return r
    return n, rk
def rand_sparse_paving(rng, n, r, tries=200):
    H = []
    for _ in range(tries):
        h = frozenset(rng.sample(range(n), r))
        if all(len(h & g) <= r - 2 for g in H): H.append(h)
    return H
def check(n, rk):
    k = n; pairs = list(itertools.combinations(range(k), 2)); pi = {pr: i for i, pr in enumerate(pairs)}; m = 1 + k + len(pairs)
    def vec(const, lin, quad):
        v = [0] * m; v[0] = const % P
        for e, c in lin.items(): v[1 + e] = (v[1 + e] + c) % P
        for (e, f), c in quad.items():
            j = 1 + k + pi[(min(e, f), max(e, f))]; v[j] = (v[j] + c) % P
        return v
    full = set(range(k)); r = rk(full)
    T = [S for S in itertools.combinations(range(k), r) if rk(S) == r]
    id_dim = m - rank_mod([vec(1, {e: 1 for e in t}, {pr: 1 for pr in itertools.combinations(t, 2)}) for t in T])
    def q_times(X, rX, f):
        lin = {f: -rX}; quad = {}
        for e in X:
            if e == f: lin[f] += 1
            else: quad[(e, f)] = quad.get((e, f), 0) + 1
        return vec(0, lin, quad)
    sep = []; sets = [(full, r, 'F')]
    for size in range(1, k):
        for X in itertools.combinations(range(k), size):
            rX = rk(X); lam = rX + rk(full - set(X)) - r
            if lam <= 1: sets.append((set(X), rX, lam))
    for X, rX, lam in sets:
        if lam in ('F', 0):
            sep.append(vec(-rX, {e: 1 for e in X}, {})); sep += [q_times(X, rX, f) for f in range(k)]
        elif not NO2:
            sep.append(vec(rX * rX - rX, {e: 2 - 2 * rX for e in X}, {pr: 2 for pr in itertools.combinations(sorted(X), 2)}))
    for s in sep:
        assert all((s[0] + sum(s[1 + e] for e in t) + sum(s[1 + k + pi[pr]] for pr in itertools.combinations(t, 2))) % P == 0 for t in T[:60])
    n2 = sum(1 for x in sets if x[2] == 1)
    return len(T), id_dim, rank_mod(sep), n2
def direct_sum(A, B):
    (na, ra), (nb, rb) = A, B
    return na + nb, lambda S: ra([e for e in S if e < na]) + rb([e - na for e in S if e >= na])
def two_sum(A, B):  # basepoint = last element of A and first of B; result drops both
    (na, ra), (nb, rb) = A, B; pa, pb = na - 1, 0
    # rank of 2-sum via the parallel connection: r(X) = r_P(X) on P(A,B) with identified point
    def rP(S):  # parallel connection rank: elements 0..na-1 of A (pa shared), then B\pb as na..na+nb-2
        SA = [e for e in S if e < na]; SB = [e - na + 1 for e in S if e >= na]
        a = ra(SA); b = rb(SB)
        return min(a + b, ra(SA + [pa]) + rb(SB + [pb]) - 1)
    N = na - 1 + nb - 1
    def rk(S): return rP([e if e < pa else e + 1 for e in S])
    return N, rk
if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 1)
    N = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 20
    fano = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)]
    vamos_H = [{0,1,2,3},{0,1,4,5},{2,3,4,5},{0,1,6,7},{2,3,6,7}]  # 5 dependent 4-sets of V8 (sparse paving)
    named = {'U(2,4)': uniform(2, 4), 'U(2,6)': uniform(2, 6), 'U(3,7)': uniform(3, 7), 'U(4,8)': uniform(4, 8),
             'F7': vector_matroid(fano, 2), 'nonFano (GF3)': vector_matroid(fano, 3),
             'Vamos': sparse_paving(8, 4, vamos_H), 'U(2,4)+U(1,3) direct': direct_sum(uniform(2, 4), uniform(1, 3)),
             'F7 (+)2 U(2,4)': two_sum(vector_matroid(fano, 2), uniform(2, 4)),
             'U(2,5) (+)2 U(1,4) (circuit leaf? no: U(3,4)=C4)': two_sum(uniform(2, 5), uniform(3, 4)),
             'U(2,5) (+)2 U(1,4) (bond leaf)': two_sum(uniform(2, 5), uniform(1, 4))}
    bad = 0
    for name, (n, rk) in named.items():
        t, d, s, n2 = check(n, rk); bad += d != s
        print(f'{name}: E={n} r={rk(range(n))} B={t} 2sep={n2} id_dim={d} sep_rank={s} {"OK" if d == s else "GAP"}', flush=True)
    for i in range(N):
        kind = rng.choice(['sp', 'sp', 'gf', 'sum'])
        if kind == 'sp':
            n = rng.randint(6, 10); r = rng.randint(2, min(5, n - 2)); M = sparse_paving(n, r, rand_sparse_paving(rng, n, r))
        elif kind == 'gf':
            p = rng.choice([2, 3, 5]); r = rng.randint(2, 4); n = rng.randint(r + 2, 10)
            while True:
                cols = [tuple(rng.randrange(p) for _ in range(r)) for _ in range(n)]
                if gf_rank(cols, p) == r: break
            M = vector_matroid(cols, p)
        else:
            a = sparse_paving(6, 3, rand_sparse_paving(rng, 6, 3)); b = rng.choice([uniform(2, 4), uniform(3, 4), uniform(1, 4), vector_matroid(fano, 2)])
            M = two_sum(a, b) if rng.random() < 0.7 else direct_sum(a, b)
        n, rk = M; t, d, s, n2 = check(n, rk); bad += d != s
        print(f'rand{i} {kind}: E={n} r={rk(range(n))} B={t} 2sep={n2} id_dim={d} sep_rank={s} {"OK" if d == s else "GAP"}', flush=True)
    print('gaps:', bad)
