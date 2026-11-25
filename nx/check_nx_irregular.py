"""Probe: is networkx.random_spanning_tree uniform (and weight-proportional) on small graphs?
Enumerate every spanning tree, sample N, chi-square against the exact law."""
import sys, time, itertools, random, math
import networkx as nx
from scipy.stats import chisquare
def trees(G):
    E=list(G.edges()); n=G.number_of_nodes(); out=[]
    for c in itertools.combinations(E,n-1):
        H=nx.Graph(); H.add_nodes_from(G); H.add_edges_from(c)
        if nx.is_tree(H): out.append(frozenset(frozenset(e) for e in c))
    return out
def run(name,G,N,weight=None,mult=True,seed=1):
    T=trees(G); idx={t:i for i,t in enumerate(T)}
    if weight:
        w=[math.prod(G.edges[tuple(e)][weight] for e in t) if mult else sum(G.edges[tuple(e)][weight] for e in t) for t in T]
    else: w=[1]*len(T)
    s=sum(w); exp=[N*x/s for x in w]; obs=[0]*len(T)
    t0=time.time()
    for k in range(N):
        S=nx.random_spanning_tree(G,weight=weight,multiplicative=mult,seed=seed*10**6+k)
        obs[idx[frozenset(frozenset(e) for e in S.edges())]]+=1
    st,p=chisquare(obs,exp)
    print(f"{name}: trees={len(T)} N={N} weight={weight} mult={mult} chi2={st:.1f} df={len(T)-1} p={p:.3g} ({time.time()-t0:.0f}s)",flush=True)
case=sys.argv[1]
if case=="grid": run("grid3x3",nx.grid_2d_graph(3,3),20000)
if case=="ear":
    K=nx.complete_graph(4); K.add_edges_from([(3,4),(4,0)]); run("K4+ear",K,20000)
if case in("mult","add"):
    W=nx.wheel_graph(5); rnd=random.Random(7)
    for u,v in W.edges(): W.edges[u,v]['w']=rnd.choice([1,2,3,5])
    run("wheel5 weighted "+case,W,8000,weight='w',mult=(case=="mult"))
