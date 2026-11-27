"""Independent check of antigravity's claim (board #130749): for a graph G, the degree-<=2 identities on
spanning-tree indicators are spanned by matroid separation rows: Foster q_E*{1,x_f}; 1-separations
(lambda(X)=0) q_X*{1,x_f}; 2-separations (lambda(X)=1) q_X(q_X+1). q_X = sum_X x_e - r(X).
Compare id_dim (monomials minus rank of tree matrix) with rank of the row span. Exact mod p=1000003.
Graphs may be multigraphs and need not be 2-connected; X runs over all nonempty proper subsets."""
import itertools, random, sys
P = 1000003
NO2 = 'no2' in sys.argv  # control: drop 2-separation rows; gaps must appear
def rk(n, E, S):
    p = list(range(n))
    def f(a):
        while p[a] != a: p[a] = p[p[a]]; a = p[a]
        return a
    r = 0
    for i in S:
        a, b = f(E[i][0]), f(E[i][1])
        if a != b: p[a] = b; r += 1
    return r
def rank_mod(rows):
    basis = {}  # pivot column -> normalized row (incremental elimination)
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
def check(n, E):
    k = len(E); pairs = list(itertools.combinations(range(k), 2)); pi = {pr: i for i, pr in enumerate(pairs)}
    m = 1 + k + len(pairs)
    def vec(const, lin, quad):
        v = [0] * m; v[0] = const % P
        for e, c in lin.items(): v[1 + e] = (v[1 + e] + c) % P
        for (e, f), c in quad.items(): v[1 + k + pi[(min(e, f), max(e, f))]] = (v[1 + k + pi[(min(e, f), max(e, f))]] + c) % P
        return v
    T = [S for S in itertools.combinations(range(k), rk(n, E, range(k))) if rk(n, E, S) == len(S)]
    rows = [vec(1, {e: 1 for e in t}, {pr: 1 for pr in itertools.combinations(t, 2)}) for t in T]
    id_dim = m - rank_mod(rows)
    r = rk(n, E, range(k)); full = set(range(k)); sep = []
    def q_times(X, rX, f):  # q_X * x_f, multilinear
        lin = {f: -rX}; quad = {}
        for e in X:
            if e == f: lin[f] += 1
            else: quad[(e, f)] = quad.get((e, f), 0) + 1
        return vec(0, lin, quad)
    def q_lin(X, rX): return vec(-rX, {e: 1 for e in X}, {})
    def q_q1(X, rX):  # q(q+1) = sum x(1-2r+1) + 2 sum pairs + r^2 - r
        return vec(rX * rX - rX, {e: 2 - 2 * rX for e in X}, {pr: 2 for pr in itertools.combinations(sorted(X), 2)})
    sets = [(full, r, 'F')]
    for size in range(1, k):
        for X in itertools.combinations(range(k), size):
            rX = rk(n, E, X); lam = rX + rk(n, E, full - set(X)) - r
            if lam <= 1: sets.append((set(X), rX, lam))
    for X, rX, lam in sets:
        if lam in ('F', 0):
            sep.append(q_lin(X, rX)); sep += [q_times(X, rX, f) for f in range(k)]
        elif not NO2: sep.append(q_q1(X, rX))
    for s in sep:  # sanity: every row vanishes on every tree
        assert all((s[0] + sum(s[1 + e] for e in t) + sum(s[1 + k + pi[pr]] for pr in itertools.combinations(t, 2))) % P == 0 for t in T[:50])
    return len(T), id_dim, rank_mod(sep)
def rand_graph(rng, n, k, multi):
    while True:
        E = []
        for _ in range(k):
            a, b = rng.sample(range(n), 2)
            if multi or (min(a, b), max(a, b)) not in [(min(x), max(x)) for x in E]: E.append((a, b))
        if len(E) == k and rk(n, E, range(k)) == n - 1: return E
if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1); bad = 0; N = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    named = {'K4': (4, list(itertools.combinations(range(4), 2))), 'bowtie': (5, [(0,1),(1,2),(2,0),(0,3),(3,4),(4,0)]),
             'theta3x3 (3 paths len 3)': (8, [(0,2),(2,3),(3,1),(0,4),(4,5),(5,1),(0,6),(6,7),(7,1)]),
             'K3,3': (6, [(a, b) for a in range(3) for b in range(3, 6)]),
             'C4 doubled edges': (4, [(0,1),(0,1),(1,2),(2,3),(2,3),(3,0)])}
    for name, (n, E) in named.items():
        t, d, s = check(n, E); bad += d != s; print(f'{name}: V={n} E={len(E)} T={t} id_dim={d} sep_rank={s} {"OK" if d == s else "GAP"}')
    for i in range(N):
        n = rng.randint(4, 7); k = rng.randint(n, min(12, n + 5)); multi = rng.random() < 0.4
        if not multi: k = min(k, n * (n - 1) // 2)
        E = rand_graph(rng, n, k, multi); t, d, s = check(n, E); bad += d != s
        print(f'rand{i}: V={n} E={k} multi={multi} T={t} id_dim={d} sep_rank={s} {"OK" if d == s else "GAP " + str(E)}')
    print('gaps:', bad)
