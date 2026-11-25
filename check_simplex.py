"""antigravity 129029: with P(s=4)=0, fixing E[s] and E[C(s,2)] of one square leaves (p0+d, p1-3d, p2+3d, p3-d);
so P(s=3) in [p3 - p1/3, p3 + p0] (when p0 <= p2/3). Check: 4x4 per-square laws by enumeration vs the LP in
results/u_lp4.json, and the 20x20 centre square's exact law from transfer-current determinants."""
import itertools, json, numpy as np
from ust import grid, all_trees
from pairs import transfer
from math import comb

def law_from_Y(Y, s):
    """P(s=k), k=0..4, by inclusion-exclusion over subsets of the square's 4 edges: P(all of S in tree) = det Y[S,S]."""
    q = {}
    for r in range(5):
        for S in itertools.combinations(s, r):
            q[S] = 1.0 if r == 0 else float(np.linalg.det(Y[np.ix_(S, S)]))
    return [sum((-1) ** (len(S) - k) * comb(len(S), k) * q[S] for S in q if len(S) >= k) for k in range(5)]

def bounds(p):
    return p[3] - p[1] / 3, p[3] + p[0], p[0] <= p[2] / 3

out = {}
for L in (4, 20):
    n, E, _ = grid(L)
    pos = {tuple(sorted(e)): k for k, e in enumerate(E)}
    sqs = [(i, j) for i in range(L - 1) for j in range(L - 1)]
    def edges(i, j):
        a = i * L + j
        return [pos[(a, a + 1)], pos[(a + L, a + L + 1)], pos[(a, a + L)], pos[(a + 1, a + 1 + L)]]
    Y = transfer(n, E)
    if L == 4:
        T = all_trees(n, E)
        res = {}
        for (i, j) in sqs:
            s = edges(i, j)
            cnt = [0] * 5
            for t in T:
                cnt[len(set(s) & set(t))] += 1
            p = [c / len(T) for c in cnt]
            pd = law_from_Y(Y, s)
            res[f"{i},{j}"] = dict(counts=cnt, det_maxerr=max(abs(a - b) for a, b in zip(p, pd)), lo=bounds(p)[0], hi=bounds(p)[1], cond=bounds(p)[2])
        out["4x4"] = dict(trees=len(T), squares=res, sum_lo=sum(r["lo"] for r in res.values()), sum_hi=sum(r["hi"] for r in res.values()),
                          lp=json.load(open("results/u_lp4.json")))
    else:
        c = (L // 2 - 1, L // 2 - 1)
        p = law_from_Y(Y, edges(*c))
        out["20x20_centre"] = dict(square=c, p=[round(x, 6) for x in p], lo=round(bounds(p)[0], 6), hi=round(bounds(p)[1], 6), cond=bounds(p)[2])
print(json.dumps(out, default=float)[:3000])
json.dump(out, open("results/simplex_bounds.json", "w"), default=float, indent=1)
