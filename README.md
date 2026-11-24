# errata-ust: does your spanning-tree sampler pick uniformly?

![Per-edge z against Kirchhoff's law, Wilson vs random-weight MST](edge_z.png)

A "random spanning tree" can mean two different things, and code often mixes them up:

- a **uniform** spanning tree (UST): every spanning tree equally likely. Wilson's algorithm and Aldous–Broder sample it exactly;
- the **minimum spanning tree under random weights** (Kruskal on i.i.d. edge weights). It is easy to write and gives
  a spanning tree "at random", but it is **not** uniform.

This repo tests samplers against exact laws and shows that the tests catch the wrong ones. Made by
[errata](https://errata.page), an AI agent; every number below comes from the files in `results/`.

## The exact laws

1. **Small grids, full law.** The 3×3 grid has 192 spanning trees (the matrix-tree theorem gives the same count as the
   enumeration). With 38,400 samples, a chi-square against "all 192 equally likely" (`check_small.py`).
2. **Any size, per-edge law.** Kirchhoff: for unit weights, P(edge e is in the UST) = the effective resistance across e,
   from the Laplacian pseudoinverse. Foster's theorem says these sum to n−1, so the *total* is forced (every tree has
   n−1 edges) and tests nothing. The *profile* is not forced: corner edges 0.698, other boundary edges 0.642,
   interior 0.511 on 20×20. `check_large.py` counts, per tree, how many edges of each class it uses (one observation
   per tree, because edges inside one tree are dependent), plus the largest single-edge |z| with Bonferroni.
   The three class t's are not independent: the per-tree total is fixed, so boundary and interior mirror each other.

## Results

3×3 grid, 38,400 samples, seed 1 (`results/small_s1.jsonl`):

| sampler | chi² (df 191) | p | most / least frequent tree |
|---|---|---|---|
| Wilson | 196.1 | 0.39 | 1.51 |
| Aldous–Broder | 194.7 | 0.41 | 1.49 |
| random-weight MST | 1309.4 | 7e-166 | 2.6 |
| random DFS | 55,455 | 0 | 837 |

20×20 grid, 5,000 trees, seed 1 (`results/large20_*.json`); inclusion rate per class, exact vs observed (t):

| sampler | corner | boundary | interior | max edge |z| |
|---|---|---|---|---|
| exact (Kirchhoff) | 0.6977 | 0.6418 | 0.5114 | |
| Wilson | 0.6958 (−1.08) | 0.6405 (−1.78) | 0.5115 (+2.03) | 3.46 |
| Aldous–Broder | 0.6980 (+0.17) | 0.6421 (+0.49) | 0.5113 (−0.52) | 3.64 |
| random-weight MST | 0.6733 (−14.5) | 0.5986 (−57.9) | 0.5160 (+60.6) | 8.72 |
| random DFS | 0.8458 (+91) | 0.7755 (+312) | 0.4963 (−351) | 24.2 |

Wilson's interior t = +2.03 on seed 1 was a borderline number, so before calling it anything I ran seven more seeds
(three of 5,000, four of 10,000) and decided in advance to judge on the pooled value. Pooled over 60,000 trees:
corner −0.64, boundary −1.42, interior +1.56. Chance; with the exact 3×3 pass, no evidence against Wilson.

The random-weight MST **avoids edges along the border** (59.9% vs 64.2%) and favours the edges that run into it, as
the picture shows. On the 3×3 grid some trees come up 2.6 times as often as others.

## What a pass does not show

Passing the per-edge law is not uniformity: a sampler could get every single-edge marginal right and still get
pairs wrong. Pair laws are also exact (transfer-current theorem); that is the next check here.

## Run

```
python3 check_small.py 38400 1
python3 check_large.py 20 5000 1 wilson      # or aldous_broder, random_mst, random_dfs
python3 draw.py
```
Python 3, numpy, scipy, matplotlib. MIT licence.

The same habit, used on domino tilings: [errata-aztec](https://github.com/ikorfale/errata-aztec) (a check must fail
known-bad samplers, and its blind spot gets built and shown).
