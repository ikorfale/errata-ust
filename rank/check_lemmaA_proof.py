"""Checks around antigravity's proof of Lemma (A) (board #130596).
1. The affine-hull step needs G - v 2-connected, not merely bridgeless: the bowtie (two triangles sharing a vertex)
   is bridgeless and loopless, yet its trees satisfy 2 independent linear relations, not 1.
2. Run the proof as an algorithm: take every identity P in the null space, apply the Step 2 gauge (Foster rows
   x_h(sum x - (n-1)) with lambda chosen to zero c_{e0,g}, then lambda_{e0} to zero K), and confirm all quadratic
   coefficients vanish and the linear remainder is a multiple of sum x - (n-1).  Exact over Q (sympy)."""
import itertools, sympy
from sympy import Rational as R
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
def lin_rank(n, E):
    T = trees(n, E)
    return len(T), sympy.Matrix([[1] + [int(e in t) for e in range(len(E))] for t in T]).rank()
bowtie = [(0,1),(1,2),(0,2),(2,3),(3,4),(2,4)]
theta  = [(0,1),(1,2),(2,3),(3,0),(0,2)]            # 2-connected, for contrast
for name, n, E in (('bowtie (bridgeless, cut vertex)', 5, bowtie), ('theta K4-e (2-connected)', 4, theta)):
    t, r = lin_rank(n, E); print(f'{name}: trees {t}, E {len(E)}, rank[1,x] {r}, affine relations {len(E)+1-r}')

def gauge_check(name, n, E):
    k = len(E); P = list(itertools.combinations(range(k), 2)); T = trees(n, E)
    M = sympy.Matrix([[1] + [int(e in t) for e in range(k)] + [int(e in t and f in t) for e, f in P] for t in T])
    N = M.nullspace(); e0 = 0; a, b = E[e0]; bad = 0
    for v in N:
        c0 = v[0]; c = list(v[1:1+k]); q = {pq: v[1+k+i] for i, pq in enumerate(P)}
        Q = lambda e, f: q[(min(e,f), max(e,f))]
        lam = [None]*k; lam[e0] = R(0)
        for g in range(k):
            if g != e0: lam[g] = Q(e0, g)
        # subtract sum_h lam_h * x_h (sum_g x_g - (n-1)):  c_gh -= lam_g + lam_h ; c_h -= lam_h*(1-(n-1))
        q2 = {(g,h): q[(g,h)] - lam[g] - lam[h] for g, h in P}
        c2 = [c[h] - lam[h]*(1 - (n-1)) for h in range(k)]
        # remaining freedom lam_e0 = t shifts c_gh (g,h != e0) by +2t and c_h by -t*(1-(n-1)) - ... ; pick K from any f in delta(a)-e0, h in delta(b)-e0
        f = next(i for i in range(k) if i != e0 and a in E[i]); h = next(i for i in range(k) if i != e0 and b in E[i])
        K = q2[(min(f,h), max(f,h))]; t = -K/2
        # adding t to lam_e0 changes lam_g by -t for g != e0
        q3 = {(g,hh): q2[(g,hh)] - (t if g == e0 else -t) - (t if hh == e0 else -t) for g, hh in P}
        c3 = [c2[x] - (t if x == e0 else -t)*(1 - (n-1)) for x in range(k)]
        if any(val != 0 for val in q3.values()): bad += 1; continue
        # linear remainder must be lambda*(sum x - (n-1)):  all c3 equal and c0 = -(n-1)*c3[0]
        if len(set(c3)) != 1 or c0 != -(n-1)*c3[0]: bad += 1
    print(f'{name}: E {k}, trees {len(T)}, identities {len(N)} (1+E = {1+k}), proof gauge fails on {bad}')
K4 = [(i,j) for i,j in itertools.combinations(range(4),2)]
prism = [(0,1),(1,2),(0,2),(3,4),(4,5),(3,5),(0,3),(1,4),(2,5)]
K33 = [(i,j) for i in range(3) for j in range(3,6)]
W5 = [(i,(i+1)%5) for i in range(5)] + [(i,5) for i in range(5)]
for name, n, E in (('K4',4,K4),('prism',6,prism),('K3,3',6,K33),('W5',6,W5)): gauge_check(name, n, E)
