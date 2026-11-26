#!/usr/bin/env python3
"""v5 on larger graphs: Petersen, and seeded random bridgeless graphs with V in [VMIN,VMAX], E <= V+EX.
Usage: check_rank_big.py SEED N VMIN VMAX EX  (prints one line per graph, GAP lines when bounds don't meet)."""
import sys, random, json, io, contextlib
sys.path.insert(0, '.')
import check_rank_cutside as cc
import networkx as nx
cc.TWO_TERMINAL = True; cc.V5 = True
def one(V, E):
    cc.grid = lambda r, c, V=V, E=E: (V, list(E))
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf): cc.run('1x1')
    return json.loads(buf.getvalue())
seed, N, vmin, vmax, ex = map(int, sys.argv[1:6])
if seed == 0:
    G = nx.petersen_graph(); E = sorted(tuple(sorted(e)) for e in G.edges())
    r = one(10, E); print('petersen', json.dumps(r)); sys.exit()
rng = random.Random(seed); done = bad = 0; cutv = 0
while done < N:
    V = rng.randint(vmin, vmax); Emax = rng.randint(V + 1, V + ex)
    E = set(); perm = list(range(V)); rng.shuffle(perm)
    for i in range(1, V): E.add(tuple(sorted((perm[i], perm[rng.randrange(i)]))))
    while len(E) < Emax:
        a, b = rng.sample(range(V), 2); E.add(tuple(sorted((a, b))))
    E = sorted(E); G = nx.Graph(E)
    if any(True for _ in nx.bridges(G)): continue
    hasc = any(True for _ in nx.articulation_points(G)); cutv += hasc
    r = one(V, E); done += 1
    gap = not r['exact'] or r['rank_lower'][0] != r['rank_lower'][1]
    bad += gap
    print('GAP' if gap else 'ok', V, len(E), 'cutvertex' if hasc else '2conn', json.dumps(r) if gap else '', flush=True)
print('summary', done, 'graphs,', cutv, 'with a cut vertex,', bad, 'with a gap')
