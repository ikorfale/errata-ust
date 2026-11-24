"""Per-edge z = (observed inclusion - effective resistance) / binomial sd, 20x20 grid, 5,000 trees, seed 1."""
import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from ust import grid, resistances, sampler
L, k = 20, 5000; n, E, nb = grid(L); R = resistances(n, E)
fig, ax = plt.subplots(1, 2, figsize=(12, 6.4))
for a, name, title in zip(ax, ('wilson', 'random_mst'), ("Wilson's algorithm (uniform)", 'Kruskal on random weights (not uniform)')):
    rng = np.random.default_rng(1); s = sampler(name, n, E, nb); f = np.zeros(len(E))
    for _ in range(k): f[list(s(rng))] += 1
    z = (f / k - R) / np.sqrt(R * (1 - R) / k)
    seg = [[(b % L, -(b // L)) for b in e] for e in E]
    lc = LineCollection(seg, cmap='RdBu_r', norm=matplotlib.colors.Normalize(-6, 6), linewidths=3.2)
    lc.set_array(np.clip(z, -6, 6)); a.add_collection(lc); a.set_xlim(-1, L); a.set_ylim(-L, 1); a.set_aspect('equal'); a.axis('off')
    a.set_title(f'{title}\nmax |z| = {np.abs(z).max():.1f}', fontsize=12)
cb = fig.colorbar(lc, ax=ax, shrink=0.7, label='z vs Kirchhoff (red: kept more often, blue: less; clipped ±6)')
fig.suptitle('Does each edge appear as often as its effective resistance? 20×20 grid, 5,000 trees each', fontsize=13)
fig.savefig('edge_z.png', dpi=110, bbox_inches='tight')
