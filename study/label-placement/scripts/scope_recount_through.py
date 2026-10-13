#!/usr/bin/env python3
"""Recount 'leader through another label' (information only) from NAME.scope-geom.json with the fixed test of
scope_run.py (_seg_poly) and update NAME.scope.json. usage: scope_recount_through.py DIR"""
import glob, json, math, pathlib, sys
from collections import Counter
src = (pathlib.Path(__file__).resolve().parent / 'scope_run.py').read_text()
ns = {'math': math}
exec(src[src.index('def _seg_poly'):src.index('# ---------- SVG offsets')], ns)
seg_poly = ns['_seg_poly']
for f in sorted(glob.glob(str(pathlib.Path(sys.argv[1]) / '*.scope-geom.json'))):
    g = json.load(open(f)); L = g['labels']
    boxes = [(min(p[0] for p in l['poly']), min(p[1] for p in l['poly']), max(p[0] for p in l['poly']), max(p[1] for p in l['poly'])) for l in L]
    order = sorted(range(len(L)), key=lambda j: boxes[j][0])
    xs = [boxes[j][0] for j in order]
    import bisect
    through = Counter(); tl = []
    for i, l in enumerate(L):
        for p, q in zip(l['lead'], l['lead'][1:]):
            sx0, sx1 = min(p[0], q[0]), max(p[0], q[0]); sy0, sy1 = min(p[1], q[1]), max(p[1], q[1])
            hi = bisect.bisect_right(xs, sx1)
            for k in range(hi):
                j = order[k]
                b = boxes[j]
                if j == i or b[2] < sx0 or sy1 < b[1] or sy0 > b[3]:
                    continue
                if seg_poly(p, q, L[j]['poly']):
                    key = f"{l['kind']}>{L[j]['kind']}"; through[key] += 1
                    if len(tl) < 2000:
                        tl.append([key, i, j])
                    break
    g['through'] = tl
    json.dump(g, open(f, 'w'))
    sf = f.replace('.scope-geom.json', '.scope.json')
    m = json.load(open(sf)); m['through_text'] = dict(through); json.dump(m, open(sf, 'w'), indent=1)
    if through:
        print(pathlib.Path(f).name, dict(through))
