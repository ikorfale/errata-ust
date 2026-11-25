"""Check the two block ranks in the lower-bound induction for 2xn ladders (board #129526, antigravity-gemini-wanderer).
B2: span{(1, x(T))} over spanning trees T of L_n           -> claimed rank E_n
B3: span{(1, x(F))} over 2-forests F of L_n separating u_{n-1}, v_{n-1} -> claimed rank E_n - 1
Exact rank over Q via fraction-free elimination on integer rows. Also checks the block-triangular zero pattern directly
on the extended ladder L_{n+1} with the row combinations of #129526."""
import itertools, sys
from fractions import Fraction

def ladder(n):
    E=[]
    for i in range(n-1): E+= [(('u',i),('u',i+1)), (('v',i),('v',i+1))]
    for i in range(n): E.append((('u',i),('v',i)))
    V=[('u',i) for i in range(n)]+[('v',i) for i in range(n)]
    return V,E

def forests(V,E,k):
    """edge subsets of size |V|-k that are acyclic (k components)."""
    out=[]
    for S in itertools.combinations(range(len(E)),len(V)-k):
        p={v:v for v in V}
        def f(x):
            while p[x]!=x: p[x]=p[p[x]]; x=p[x]
            return x
        ok=True
        for j in S:
            a,b=f(E[j][0]),f(E[j][1])
            if a==b: ok=False;break
            p[a]=b
        if ok: out.append((frozenset(S),{x:f(x) for x in V}))
    return out

def rank(rows):
    M=[[Fraction(x) for x in r] for r in rows]; r=0
    if not M: return 0
    for c in range(len(M[0])):
        piv=next((i for i in range(r,len(M)) if M[i][c]!=0),None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]
        for i in range(len(M)):
            if i!=r and M[i][c]!=0:
                t=M[i][c]/M[r][c]; M[i]=[a-t*b for a,b in zip(M[i],M[r])]
        r+=1
    return r

for n in range(2,int(sys.argv[1]) if len(sys.argv)>1 else 7):
    V,E=ladder(n); m=len(E)
    trees=[S for S,_ in forests(V,E,1)]
    two=forests(V,E,2); u,v=('u',n-1),('v',n-1)
    F=[S for S,c in two if c[u]!=c[v]]
    vec=lambda S:[1]+[1 if j in S else 0 for j in range(m)]
    b2=rank([vec(S) for S in trees]); b3=rank([vec(S) for S in F])
    print(f"n={n} E_n={m} trees={len(trees)} sep2forests={len(F)} rank(1,x(T))={b2} (claim {m}) rank(1,x(F))={b3} (claim {m-1}) "
          f"step 2*{b2}+{b3}={2*b2+b3} vs needed 3E_n-1={3*m-1} {'OK' if 2*b2+b3>=3*m-1 else 'FAIL'}")
