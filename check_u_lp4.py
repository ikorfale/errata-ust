"""4x4: is the centre square's U total (P(s=3)) fixed by the single-edge and pair laws, as every square is on 3x3?
Enumerate all spanning trees, then an LP over laws matching every single-edge and pair probability of the uniform law
exactly; min and max of P(s=3) for the centre square (and the sum over all 9 squares)."""
import itertools, json, numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix, vstack
L = 4; n = L * L
E = [(a, a + 1) for a in range(n) if a % L < L - 1] + [(a, a + L) for a in range(n - L)]
def is_tree(c):
    p = list(range(n))
    def f(x):
        while p[x] != x: p[x] = p[p[x]]; x = p[x]
        return x
    for k in c:
        a, b = f(E[k][0]), f(E[k][1])
        if a == b: return False
        p[a] = b
    return True
T = [c for c in itertools.combinations(range(len(E)), n - 1) if is_tree(c)]
M = np.zeros((len(T), len(E)), np.int8)
for i, c in enumerate(T): M[i, list(c)] = 1
pos = {e: k for k, e in enumerate(E)}
sq = [[pos[(a, a + 1)], pos[(a + L, a + L + 1)], pos[(a, a + L)], pos[(a + 1, a + 1 + L)]] for a in [i * L + j for i in range(L - 1) for j in range(L - 1)]]
s3 = np.array([(M[:, s].sum(1) == 3) for s in sq]).astype(float)   # squares x trees
rows = [M.T.astype(float)] + [(M[:, e] * M[:, f]).astype(float)[None, :] for e, f in itertools.combinations(range(len(E)), 2)]
A = csr_matrix(np.vstack([np.ones((1, len(T)))] + rows)); u = np.full(len(T), 1 / len(T)); b = A @ u
out = dict(trees=len(T), constraints=A.shape[0], uniform_centre=round(float(s3[4] @ u), 6), uniform_all=round(float(s3.sum(0) @ u), 6))
for name, c in [('centre', s3[4]), ('all9', s3.sum(0))]:
    lo = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method='highs'); hi = linprog(-c, A_eq=A, b_eq=b, bounds=(0, None), method='highs')
    out[name] = [round(lo.fun, 6) if lo.status == 0 else lo.message, round(-hi.fun, 6) if hi.status == 0 else hi.message]
print(json.dumps(out))
