#!/usr/bin/env python3
"""LP-D26: do the labels of a cluster follow the circumference?

The Owner: labels placed along the circumference look better; a jagged stack is not beautiful. For every cluster
(labels of one column and one half whose boxes are at most 6.1 px apart, as in work/J/lop.py, but also when all
leaders are straight) the radial distance of each leader end from the ring centre is taken in stack order.
`step` is the largest change of that distance between two neighbours (0 on an arc concentric with the ring),
`spread` is largest minus smallest. A cluster is "off the arc" when step > `--step` px (default 2.0).
Side clusters (every leader end within 60 degrees of the horizontal) and pole clusters are counted apart: rows
stacked above or below the ring cannot follow it.

LP-D30: with `--rowem EM` the two distances that mean "packed" (cluster reach 6.1 px, loose pair above 5.0 px) grow by
what the font-relative row rule adds to 3 px (max(0, EM * font size - 0.1 - 3.0)), so that packed rows of large
letters are not counted as loose.

LP-D31: `--rowmin PX` is the floor of that rule.

usage: arc_fit.py [-v] [--step PX] [--rowem EM] [--rowmin PX] NAME=DIR:SUFFIX ... | FILE.labels.json ...
"""
import glob, json, math, os, sys
sys.path.insert(0, '/home/kawato/gbdraw-baselines/label-placement-20261011/work/J')
import lop

ROWEM = 0.0
ROWMIN = 0.0


def extra(labels):
    fs = max((float(l.get('font_size', 0.0)) for l in labels if not l.get('is_embedded')), default=0.0)
    return max(0.0, max(ROWMIN, ROWEM * fs - 0.1 if ROWEM else 0.0) - 3.0)


def loose(path, lo=5.0, hi=20.0):
    """Round 4: pairs of neighbours (consecutive in a column, both leader ends within 60 degrees of the horizontal, box gap in (lo, hi] px) with at
    least one bent leader.  Returns (pairs with a straight neighbour, pairs of two bent leaders)."""
    rows = []; labels = json.load(open(path))['labels']; lo += extra(labels)
    for l in labels:
        pts = l.get('leader_points')
        if l.get('is_embedded') or not pts or len(pts) < 2 or not l.get('aabb_local') or abs(float(l.get('rotation_deg', 0.0))) > 1e-6:
            continue
        b = l['aabb_local']; e = pts[-1]
        rows.append(dict(t=l.get('proto_tag') or '?', half=l['feature_anchor_y'] < 0, y0=b[1], y1=b[3], bent=abs(lop.turns(pts)) >= 0.5,
                         steep=abs(math.degrees(math.atan2(abs(e[1]), abs(e[0]))))))
    mixed = both = 0
    for key in {(r['t'], r['half']) for r in rows}:
        col = sorted((r for r in rows if (r['t'], r['half']) == key), key=lambda r: r['y0'])
        for a, b in zip(col, col[1:]):
            if max(a['steep'], b['steep']) <= 60.0 and lo < b['y0'] - a['y1'] <= hi and (a['bent'] or b['bent']):
                mixed += not (a['bent'] and b['bent']); both += a['bent'] and b['bent']
    return mixed, both


def clusters(path):
    rows = []; labels = json.load(open(path))['labels']; reach = 6.1 + extra(labels)
    for l in labels:
        pts = l.get('leader_points')
        if l.get('is_embedded') or not pts or len(pts) < 2 or not l.get('aabb_local') or abs(float(l.get('rotation_deg', 0.0))) > 1e-6:
            continue
        b = l['aabb_local']; e = pts[-1]
        rows.append(dict(t=l.get('proto_tag') or '?', half=l['feature_anchor_y'] < 0, y0=b[1], y1=b[3], r=math.hypot(*e),
                         steep=abs(math.degrees(math.atan2(abs(e[1]), abs(e[0])))), name=l['label_text'], bent=abs(lop.turns(pts)) >= 0.5))
    out = []
    for key in sorted({(r['t'], r['half']) for r in rows}):
        col = sorted((r for r in rows if (r['t'], r['half']) == key), key=lambda r: r['y0']); cur = []
        for r in col:
            if cur and r['y0'] - cur[-1]['y1'] <= reach:
                cur.append(r)
            else:
                out.append(cur); cur = [r]
        out.append(cur)
    res = []
    for c in out:
        if len(c) < 2:
            continue
        rs = [r['r'] for r in c]
        res.append(dict(n=len(c), side=max(r['steep'] for r in c) <= 60.0, step=max(abs(b - a) for a, b in zip(rs, rs[1:])),
                        spread=max(rs) - min(rs), zz=sum(1 for a, b, c_ in zip(rs, rs[1:], rs[2:]) if (b - a > 2.0 and b - c_ > 2.0) or (a - b > 2.0 and c_ - b > 2.0)), names=[r['name'] for r in c], r=[round(v, 1) for v in rs]))
    return res


if __name__ == '__main__':
    args = sys.argv[1:]; verbose = '-v' in args; lim = 2.0
    if '--rowmin' in args:
        ROWMIN = float(args[args.index('--rowmin') + 1]); del args[args.index('--rowmin'):args.index('--rowmin') + 2]
    if '--rowem' in args:
        ROWEM = float(args[args.index('--rowem') + 1]); del args[args.index('--rowem'):args.index('--rowem') + 2]
    if '--step' in args:
        lim = float(args[args.index('--step') + 1]); del args[args.index('--step'):args.index('--step') + 2]
    for spec in (a for a in args if a != '-v'):
        if '=' in spec:
            name, rest = spec.split('=', 1); folder, suffix = rest.rsplit(':', 1)
            files = sorted(glob.glob(os.path.join(folder, f'*__{suffix}.labels.json')))
        else:
            name, files = os.path.basename(spec), [spec]
        side = [0, 0]; pole = [0, 0]; steps = []; lz = [0, 0]; zz = 0
        for f in files:
            m_, b_ = loose(f); lz[0] += m_; lz[1] += b_
            for c in clusters(f):
                k = side if c['side'] else pole
                zz += c['zz']
                k[0] += 1; k[1] += c['step'] > lim
                if c['side']:
                    steps.append(c['step'])
                if verbose and c['side'] and c['step'] > lim:
                    print(f"  {os.path.basename(f).split('__')[0]}: step {c['step']:.1f} spread {c['spread']:.1f} r {c['r'][:8]} {c['names'][:5]}")
        steps.sort()
        print(f"{name}: side clusters off the arc {side[1]} of {side[0]} (median step {steps[len(steps) // 2]:.1f} px, max {steps[-1]:.1f}); pole clusters {pole[1]} of {pole[0]}; zigzag labels {zz} (all clusters); loose pairs {sum(lz)} (next to a straight leader {lz[0]}, two bent {lz[1]})" if steps else f'{name}: no side clusters')
