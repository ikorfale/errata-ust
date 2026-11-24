"""Pair law (Burton-Pemantle transfer-current theorem): with Y = B L^+ B^T (B the oriented edge-vertex incidence,
unit weights), P(e and f in T) = Y_ee Y_ff - Y_ef^2 for a uniform spanning tree. Exact at any size, no enumeration.
Pair classes on the grid, by geometry:
- corner: the two edges at a degree-2 vertex (one of them is forced, so P = 2R - 1)
- perp: two perpendicular edges sharing a vertex (not at a corner)
- straight: two collinear edges sharing a vertex
- face: opposite sides of one unit square (parallel, no shared vertex)
One observation per tree: the number of pairs of a class with both edges in the tree; exact mean = sum of P over the class;
t uses the empirical sd across trees."""
import numpy as np
from ust import grid


def transfer(n, E):
    B = np.zeros((len(E), n))
    for k, (a, b) in enumerate(E): B[k, a] = 1; B[k, b] = -1
    Lap = B.T @ B
    return B @ np.linalg.pinv(Lap) @ B.T


def pair_classes(L):
    n, E, nb = grid(L); deg = [len(x) for x in nb]
    horiz = [b - a == 1 for a, b in E]
    cls = {'corner': [], 'perp': [], 'straight': [], 'face': []}
    for v in range(n):
        ks = [k for _, k in nb[v]]
        for i in range(len(ks)):
            for j in range(i + 1, len(ks)):
                e, f = ks[i], ks[j]
                if deg[v] == 2: cls['corner'].append((e, f))
                elif horiz[e] != horiz[f]: cls['perp'].append((e, f))
                else: cls['straight'].append((e, f))
    pos = {e: k for k, e in enumerate(E)}
    for i in range(L - 1):
        for j in range(L - 1):
            a = i * L + j
            cls['face'].append((pos[(a, a + 1)], pos[(a + L, a + L + 1)]))   # top and bottom
            cls['face'].append((pos[(a, a + L)], pos[(a + 1, a + 1 + L)]))   # left and right
    return n, E, nb, {c: np.array(v) for c, v in cls.items()}


def exact_pairs(Y, P):
    e, f = P[:, 0], P[:, 1]
    return Y[e, e] * Y[f, f] - Y[e, f] ** 2
