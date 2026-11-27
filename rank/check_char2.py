"""Characteristic-2 follow-up to check_sep_span_matroid.py (board #130822, #130830, #130838).
Over GF(2) the 2-separation rows q_X(q_X+1) are 0 and Foster loses rank (U(1,3): det -2). theone (#130838) proposed
the integral halves h_X = q_X(q_X+1)/2 = e2(X) + (1-r(X)) s(X) + C(r(X),2), s = sum of x over X, e2 = pairs in X.
Modes: base   = the char-0 generator set reduced mod 2 (expect gaps);
       half2  = base + h_X for 2-separations;
       halfall= base + h_X for E, every 1-separation and every 2-separation;
       g      = base + only g = h_E (theone #130853: h_{E-f} = g mod Foster for non-coloop f);
       gsep   = base + g + h_X only for 2-separations with |X|, |E-X| >= 2 (no complement-of-singleton rows).
Each h_X is checked to vanish on every basis before use. Usage: check_char2.py MODE [SEED [N_RANDOM]] (positional; P=2 fixed; "seed 2" would silently mean seed 1, N=2)."""
import itertools, random, sys
import check_sep_span_matroid as C
C.P = 2
MODE = sys.argv[1] if len(sys.argv) > 1 else 'base'
def comb2(r): return r * (r - 1) // 2
def check(n, rk):
    k = n; pairs = list(itertools.combinations(range(k), 2)); pi = {pr: i for i, pr in enumerate(pairs)}; m = 1 + k + len(pairs)
    def vec(const, lin, quad):
        v = [0] * m; v[0] = const % 2
        for e, c in lin.items(): v[1 + e] = (v[1 + e] + c) % 2
        for (e, f), c in quad.items():
            j = 1 + k + pi[(min(e, f), max(e, f))]; v[j] = (v[j] + c) % 2
        return v
    full = set(range(k)); r = rk(full)
    T = [S for S in itertools.combinations(range(k), r) if rk(S) == r]
    id_dim = m - C.rank_mod([vec(1, {e: 1 for e in t}, {pr: 1 for pr in itertools.combinations(t, 2)}) for t in T])
    def q_times(X, rX, f):
        lin = {f: -rX}; quad = {}
        for e in X:
            if e == f: lin[f] += 1
            else: quad[(e, f)] = quad.get((e, f), 0) + 1
        return vec(0, lin, quad)
    def half(X, rX): return vec(comb2(rX), {e: 1 - rX for e in X}, {pr: 1 for pr in itertools.combinations(sorted(X), 2)})
    sep = []; sets = [(full, r, 'F')]
    for size in range(1, k):
        for X in itertools.combinations(range(k), size):
            rX = rk(X); lam = rX + rk(full - set(X)) - r
            if lam <= 1: sets.append((set(X), rX, lam))
    for X, rX, lam in sets:
        if lam in ('F', 0):
            sep.append(vec(-rX, {e: 1 for e in X}, {})); sep += [q_times(X, rX, f) for f in range(k)]
            if MODE == 'halfall' or (MODE in ('g', 'gsep') and lam == 'F'): sep.append(half(X, rX))
        else:
            sep.append(vec(rX * rX - rX, {e: 2 - 2 * rX for e in X}, {pr: 2 for pr in itertools.combinations(sorted(X), 2)}))
            if MODE in ('half2', 'halfall') or (MODE == 'gsep' and 2 <= len(X) <= k - 2): sep.append(half(X, rX))
    for s in sep:
        assert all((s[0] + sum(s[1 + e] for e in t) + sum(s[1 + k + pi[pr]] for pr in itertools.combinations(t, 2))) % 2 == 0 for t in T), 'row does not vanish'
    return len(T), id_dim, C.rank_mod(sep), sum(1 for x in sets if x[2] == 1)
C.check = check
if __name__ == '__main__':
    sys.argv = [sys.argv[0]] + sys.argv[2:]
    exec(open(C.__file__).read().split("if __name__ == '__main__':")[1].replace('\n    ', '\n'), {**C.__dict__, 'check': check, 'sys': sys})
