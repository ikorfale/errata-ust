"""Rank of {1, x_e, x_e x_f} over spanning trees for the 3x3 grid and its planar dual (a multigraph).
Duality: T* = E minus T, so x* = 1 - x, an affine map; ranks must match. Exact rank via sympy over Q."""
import itertools, sympy
def trees(n, E):
    out = []
    for c in itertools.combinations(range(len(E)), n - 1):
        p = list(range(n))
        def f(a):
            while p[a] != a: p[a] = p[p[a]]; a = p[a]
            return a
        ok = True
        for i in c:
            a, b = f(E[i][0]), f(E[i][1])
            if a == b: ok = False; break
            p[a] = b
        if ok: out.append(set(c))
    return out
def rank(n, E):
    T = trees(n, E); k = len(E)
    rows = [[1] + [int(e in t) for e in range(k)] + [int(e in t and f in t) for e, f in itertools.combinations(range(k), 2)] for t in T]
    return len(T), 1 + k + k * (k - 1) // 2, sympy.Matrix(rows).rank()
# 3x3 grid vertices (r,c) -> 3r+c
G = [(3*r+c, 3*r+c+1) for r in range(3) for c in range(2)] + [(3*r+c, 3*(r+1)+c) for r in range(2) for c in range(3)]
# faces: 0..3 = cells (r,c) -> 2r+c ; 4 = outer
def face_pair(e):
    a, b = e; ra, ca = divmod(a, 3); rb, cb = divmod(b, 3)
    cells = []
    if ra == rb:   # horizontal edge in row ra between columns ca, ca+1: cells above (ra-1, ca) and below (ra, ca)
        for r in (ra - 1, ra): cells.append(2*r+ca if 0 <= r <= 1 else 4)
    else:          # vertical edge in column ca between rows ra, ra+1: cells left (ra, ca-1) and right (ra, ca)
        for c in (ca - 1, ca): cells.append(2*ra+c if 0 <= c <= 1 else 4)
    return tuple(cells)
D = [face_pair(e) for e in G]
print('grid edges', len(G), 'dual edges', D)
for name, n, E in (('3x3 grid', 9, G), ('dual (multigraph)', 5, D)):
    t, m, r = rank(n, E); print(f'{name}: trees {t}, m {m}, rank {r}, identities {m - r}')
par = [(i, j) for i, j in itertools.combinations(range(len(D)), 2) if sorted(D[i]) == sorted(D[j])]
print('parallel pairs in dual:', len(par))
