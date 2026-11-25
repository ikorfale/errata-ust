"""networkx.random_spanning_tree returns non-trees when the graph has a self-loop.
A 4-cycle with a loop at node 0 has exactly 4 spanning trees (a loop is never in a tree).
Run: python -I selfloop_repro.py   (prints networkx version, invalid share, and the workaround's law)."""
import collections, networkx as nx
from scipy.stats import chisquare
N = 4000
def sample(G, **kw):
    bad = 0; law = collections.Counter()
    for k in range(N):
        T = nx.random_spanning_tree(G, seed=k, **kw)
        if not nx.is_tree(T) or T.number_of_nodes() != G.number_of_nodes(): bad += 1
        else: law[frozenset(frozenset(e) for e in T.edges())] += 1
    return bad, law
print("networkx", nx.__version__)
G = nx.cycle_graph(4); G.add_edge(0, 0)
for u, v in G.edges(): G.edges[u, v]["w"] = 2
for label, kw in [("uniform", {}), ("multiplicative w", {"weight": "w"}), ("additive w", {"weight": "w", "multiplicative": False})]:
    bad, law = sample(G, **kw)
    print(f"C4 + loop, {label}: {bad}/{N} outputs are not spanning trees ({100*bad/N:.0f}%); valid outputs: {len(law)} distinct")
H = G.copy(); H.remove_edges_from(list(nx.selfloop_edges(H)))
bad, law = sample(H)
print(f"workaround (loops removed first): {bad}/{N} invalid, {len(law)} distinct trees, counts {sorted(law.values())}, chi2 p={chisquare(list(law.values())).pvalue:.3g}")
