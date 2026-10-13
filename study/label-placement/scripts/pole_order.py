#!/usr/bin/env python3
"""LP-D33: how orderly are the leaders, by distance of the feature from a pole?

The Owner: toward the top and bottom ends the leaders of the dense maps get disordered. For every outside label the
leader is reduced to its chord (feature anchor -> leader end): length L and lean = signed angle between the chord
and the outward radius at the feature (0 = radial). Labels are taken in angular order; a pair = two neighbours whose
features lie in the same bin of pole distance (angle between the feature's radius and the vertical: 0 = pole,
90 = equator) and in the same half. Per bin:
  n          labels
  len        mean leader length (px)
  lean       mean |lean| (deg)
  d_lean     mean |lean difference| between neighbours (deg): 0 in a fan whose leaders turn evenly
  jumps      pairs whose lean differs by more than --jump deg (default 15)
  d_len      mean |length difference| between neighbours (px)
  mixed      pairs with one straight and one bent leader (turn < 1 deg = straight)
  by_text    px of leader that run within --near px (default 5) of the box of another label, per label
  n_by       leaders with 15 px or more of that
  rough      outline roughness: the near end of each label box (the end facing its feature, at mid height) is taken in
             angular order; mean length of the second difference over three neighbours of one side (px): 0 when the
             ends lie evenly spaced on a line, small on a smooth outline with even rows
  breaks     neighbours of one side whose near ends are more than 2.5 row heights apart
  bundle     px over which a leader runs within 4 px of its neighbour's leader, more than 15 px away from both features, per pair

usage: pole_order.py [--jump DEG] [--bins 20,40,60,90] NAME=FILE.labels.json ...
"""
import json, math, sys
sys.path.insert(0, '/home/kawato/gbdraw-baselines/label-placement-20261011/scripts')
import lp_metrics
from gbdraw.labels import circular as lc

a = sys.argv[1:]
NEAR = 5.0; JUMP = 15.0; BINS = [20.0, 40.0, 60.0, 90.0]
while a and a[0].startswith('--'):
    if a[0] == '--jump': JUMP = float(a[1])
    if a[0] == '--near': NEAR = float(a[1])
    if a[0] == '--bins': BINS = [float(v) for v in a[1].split(',')]
    a = a[2:]


def samples(pts, step=1.0):
    out = []
    for p, q in zip(pts, pts[1:]):
        n = max(1, int(math.dist(p, q) / step))
        out += [(p[0] + (q[0] - p[0]) * k / n, p[1] + (q[1] - p[1]) * k / n) for k in range(n)]
    return out + [pts[-1]]


def box_dist(p, b):
    return math.hypot(max(b[0] - p[0], 0.0, p[0] - b[2]), max(b[1] - p[1], 0.0, p[1] - b[3]))


def rows(path):
    data = json.load(open(path)); total = data['total_length']
    labels = [l for l in data['labels'] if not l.get('is_embedded')]
    boxes = []
    for l in labels:
        x0, y0, x1, y1 = lp_metrics.label_box(l, total, lc)
        boxes.append((min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)))
    out = []
    for k, l in enumerate(labels):
        pts = lp_metrics.label_leader(l)
        ax, ay = pts[0]; ex, ey = pts[-1]
        if math.hypot(ex - ax, ey - ay) < 0.5:
            continue
        theta = math.degrees(math.atan2(ax, -ay)) % 360.0            # clockwise from 12 o'clock
        lean = (math.degrees(math.atan2(ex - ax, -(ey - ay))) - theta + 180.0) % 360.0 - 180.0
        length = sum(math.dist(p, q) for p, q in zip(pts, pts[1:]))
        pole = math.degrees(math.atan2(abs(ax), abs(ay)))
        sm = samples(pts)
        by = sum(1 for p in sm if any(box_dist(p, b) < NEAR for j, b in enumerate(boxes) if j != k))
        out.append(dict(theta=theta, lean=lean, len=length, pole=pole, top=ay < 0, turn=lp_metrics.leader_turns([pts])[0], text=l['label_text'], by=by, sm=sm, a=(ax, ay),
                        end=((boxes[k][0] if boxes[k][0] + boxes[k][2] > 2 * ax else boxes[k][2]), (boxes[k][1] + boxes[k][3]) / 2), h=boxes[k][3] - boxes[k][1], right=boxes[k][0] + boxes[k][2] > 2 * ax))
    return sorted(out, key=lambda r: r['theta'])


def which(v):
    return next(k for k, hi in enumerate(BINS) if v <= hi + 1e-9)


print(f'{"render":<14}{"pole dist":>10}{"n":>5}{"len":>7}{"lean":>7}{"d_lean":>8}{"jumps":>10}{"d_len":>7}{"mixed":>7}{"by_text":>9}{"n_by":>6}{"bundle":>8}{"rough":>7}{"breaks":>8}')
for spec in a:
    name, path = spec.split('=', 1)
    rs = rows(path); n = len(rs)
    agg = [dict(n=0, len=0.0, lean=0.0, pairs=0, dl=0.0, j=0, dn=0.0, mx=0, by=0, nby=0, bu=0, ro=0.0, tr=0, br=0) for _ in BINS]
    for i, r in enumerate(rs):
        b = agg[which(r['pole'])]; b['n'] += 1; b['len'] += r['len']; b['lean'] += abs(r['lean']); b['by'] += r['by']; b['nby'] += r['by'] >= 15
        s = rs[(i + 1) % n]; t = rs[(i + 2) % n]
        if which(s['pole']) == which(r['pole']) and s['top'] == r['top'] and s['right'] == r['right']:
            b['br'] += math.dist(r['end'], s['end']) > 2.5 * r['h']
            if which(t['pole']) == which(r['pole']) and t['top'] == r['top'] and t['right'] == r['right']:
                b['ro'] += math.hypot(r['end'][0] - 2 * s['end'][0] + t['end'][0], r['end'][1] - 2 * s['end'][1] + t['end'][1]); b['tr'] += 1
        if which(s['pole']) == which(r['pole']) and s['top'] == r['top'] and abs((s['theta'] - r['theta'] + 180) % 360 - 180) < 30:
            d = abs((s['lean'] - r['lean'] + 180) % 360 - 180)
            b['pairs'] += 1; b['dl'] += d; b['j'] += d > JUMP; b['dn'] += abs(s['len'] - r['len']); b['mx'] += (r['turn'] < 1.0) != (s['turn'] < 1.0)
            far = [q for q in s['sm'] if math.dist(q, s['a']) > 15 and math.dist(q, r['a']) > 15]
            b['bu'] += sum(1 for p in r['sm'] if math.dist(p, r['a']) > 15 and math.dist(p, s['a']) > 15 and any(math.dist(p, q) < 4.0 for q in far))
    lo = 0.0
    for b, hi in zip(agg, BINS):
        if b['n']:
            p = max(1, b['pairs'])
            print(f'{name:<14}{f"{lo:.0f}-{hi:.0f}":>10}{b["n"]:>5}{b["len"] / b["n"]:>7.0f}{b["lean"] / b["n"]:>7.0f}{b["dl"] / p:>8.1f}{str(b["j"]) + "/" + str(b["pairs"]):>10}{b["dn"] / p:>7.1f}{b["mx"]:>7}{b["by"] / b["n"]:>9.1f}{b["nby"]:>6}{b["bu"] / p:>8.1f}{b["ro"] / max(1, b["tr"]):>7.1f}{b["br"]:>8}')
        lo = hi
