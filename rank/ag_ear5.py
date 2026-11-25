#!/usr/bin/env python3
"""antigravity #129563 example: K4 on {0,2,3,4} plus ear 0-1-2. Family v3 with and without the cut-rung row; plus F_n = T(L_{n+1}) - 3 T(L_n) on ladders."""
import sys, json, io, contextlib
sys.path.insert(0, '.')
import check_rank_cutside as cc
V = 5
E = [(0,1),(1,2),(0,3),(0,4),(2,3),(0,2),(2,4),(3,4)]   # K4 on {0,2,3,4} + ear 0-1-2
cc.grid = lambda r, c: (V, list(E))
for nr in (False, True):
    cc.NO_RUNG = nr
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf): cc.run('1x1')
    res = json.loads(buf.getvalue())
    print('NO_RUNG' if nr else 'v3', res)
