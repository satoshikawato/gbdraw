#!/usr/bin/env python3
"""Where the leaders of a render are bent (round 3 of prototype P): usage breakdown_p.py NAME=DIR:SUFFIX ...
Per leader (labels with a drawn leader): straight (turn under 0.5 deg); slight bend (0.5 to 8 deg); then by cluster (labels of one
column and half, boxes at most 6.1 px apart, as work/J/lop.py): lone label bent 8 deg or more; bent at the outside of a cluster's straight
run (R3); bent in a cluster without any straight leader; bent between two straight ones (mixed)."""
import glob, json, os, sys
T = '/home/kawato/gbdraw-baselines/label-placement-20261011'
sys.path.insert(0, f'{T}/work/J'); import lop

def tally(path):
    rows = []
    for l in json.load(open(path))['labels']:
        if l.get('is_embedded') or not l.get('leader_points') or len(l['leader_points']) < 2: continue
        b = l['aabb_local']; rows.append(dict(t=l.get('proto_tag') or '?', half=l['feature_anchor_y'] < 0, y0=b[1], y1=b[3], s=abs(lop.turns(l['leader_points']))))
    c = dict(straight=0, slight=0, lone=0, outer=0, nostraight=0, between=0)
    for key in {(r['t'], r['half']) for r in rows}:
        col = sorted((r for r in rows if (r['t'], r['half']) == key), key=lambda r: r['y0']); cur, groups = [], []
        for r in col:
            if cur and r['y0'] - cur[-1]['y1'] <= 6.1: cur.append(r)
            else:
                if cur: groups.append(cur)
                cur = [r]
        groups.append(cur)
        for g in groups:
            flat = [r['s'] < 0.5 for r in g]
            first = flat.index(True) if True in flat else None; last = len(g) - 1 - flat[::-1].index(True) if first is not None else None
            for i, r in enumerate(g):
                if flat[i]: c['straight'] += 1
                elif r['s'] < 8.0: c['slight'] += 1
                elif len(g) == 1: c['lone'] += 1
                elif first is None: c['nostraight'] += 1
                elif first < i < last: c['between'] += 1
                else: c['outer'] += 1
    return c

for spec in sys.argv[1:]:
    name, rest = spec.split('=', 1); folder, suffix = rest.rsplit(':', 1)
    tot = {}
    for f in sorted(glob.glob(f'{T}/{folder}/*__{suffix}.labels.json')):
        for k, v in tally(f).items(): tot[k] = tot.get(k, 0) + v
    print(f"{name}: " + ' | '.join(f'{k} {v}' for k, v in tot.items()) + f" | all {sum(tot.values())}")
