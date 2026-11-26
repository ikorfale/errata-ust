#!/usr/bin/env python3
"""antigravity #130010 example: V=7, E=9, cut vertex 6 holding triangle {2,4,6}. Family v4 vs v5."""
import sys, json, io, contextlib
sys.path.insert(0, '.')
import check_rank_cutside as cc
E = [(0,1),(0,5),(1,3),(1,6),(2,4),(2,6),(3,5),(4,6),(5,6)]
cc.grid = lambda r, c: (7, list(E))
for v5 in (False, True):
    cc.TWO_TERMINAL = True; cc.V5 = v5
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf): cc.run('1x1')
    print('v5' if v5 else 'v4', buf.getvalue().strip())
