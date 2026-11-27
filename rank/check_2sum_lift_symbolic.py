# Symbolic check (sympy, ranks as symbols) of the 2-sum lift identities in board #130896, and theone #130908 correction.
import sympy as sp, itertools
r1,r2=sp.symbols('r1 r2')
X=sp.symbols('x1:4'); Y=sp.symbols('y1:4')
def e2(V): return sum(a*b for a,b in itertools.combinations(V,2))
def h(V,r): return e2(V)+(1-r)*sum(V)+r*(r-1)/2
sX,sY=sum(X),sum(Y); qX=sX-r1; qY=sY-r2
rM=r1+r2-1; qM=sX+sY-rM
xe=1+qX
g2=e2(Y)+xe*sY - r2*(r2-1)/2          # e2(Y u e) - C(r2,2)
hE=h(X+Y,rM); hX=h(X,r1)
print('LHS-(hE-hX) =', sp.factor(sp.expand(g2-(hE-hX))))
print('LHS-(hE-hX+(r2-1)qM) =', sp.expand(g2-(hE-hX+(r2-1)*qM)))
d=sp.expand(g2-(hE-hX))
print('r1=r2=2, sX=sY=2:', d.subs({r1:2,r2:2}).subs({X[0]:1,X[1]:1,X[2]:0,Y[0]:1,Y[1]:1,Y[2]:0}))
# identity 3: (1+qX)(qY+1) vs hE-hX-hY
hY=h(Y,r2)
print('id3 diff =', sp.factor(sp.expand((1+qX)*(qY+1)-(hE-hX-hY))))
ra=sp.symbols('ra'); A=Y[:2]
lhs=e2(A)+xe*sum(A)+(1-ra)*(sum(A)+xe)+ra*(ra-1)/2
rhs=h(X+A,r1+ra-1)-h(X,r1)
print('id1 general A diff =', sp.factor(sp.expand(lhs-rhs)))
# Foster x_a q rows: red(x_a q_M2)|sub vs red(x_a q_M)
a=Y[0]
def red(p):
    p=sp.expand(p)
    for v in X+Y: p=p.subs(v**2,v)
    return sp.expand(p)
qM2=sY+xe-r2
print('id2 diff =', sp.expand(red(a*qM2)-red(a*qM)))
