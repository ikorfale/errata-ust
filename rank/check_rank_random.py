#!/usr/bin/env python3
"""Family v3 completeness on random bridgeless connected graphs (seeded). Prints graphs where rank bounds don't meet."""
import sys, random, json
sys.path.insert(0, '.')
import check_rank_cutside as cc
import check_rank_family as cf
rng = random.Random(int(sys.argv[1])); N = int(sys.argv[2]); bad = 0; done = 0
while done < N:
    V = rng.randint(5, 8); Emax = rng.randint(V, min(V * (V - 1) // 2, V + 5))
    E = set()
    perm = list(range(V)); rng.shuffle(perm)
    for i in range(1, V): E.add(tuple(sorted((perm[i], perm[rng.randrange(i)]))))
    while len(E) < Emax:
        a, b = rng.sample(range(V), 2); E.add(tuple(sorted((a, b))))
    E = sorted(E)
    import networkx as nx
    if any(True for _ in nx.bridges(nx.Graph(E))): continue   # scope: bridgeless graphs only
    cc.grid = lambda r, c, V=V, E=E: (V, list(E))   # reuse run() with this graph
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf): cc.run('1x1')
    res = json.loads(buf.getvalue()); done += 1
    if not res['exact'] or res['rank_lower'][0] != res['rank_lower'][1]:
        bad += 1; print('GAP', V, E, res, flush=True)
    else:
        print('ok', V, len(E), res['T'], res['m'], res['rank_lower'][0], 'cuts', res['cuts'], 'rung', res['cut_rung_rows'], 'old_exact', res['upper_old_family'] == res['rank_lower'][0], flush=True)
print('summary', done, 'graphs,', bad, 'with a gap')
