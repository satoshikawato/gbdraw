#!/usr/bin/env python3
"""LP-D37 / LP-D38: do the outside labels stand on one arc, and how near 12 and 6 o'clock do the end labels stand?

Per render, from the near end of every outside label (the middle of the box side that faces the vertical axis):
  A, B     half width and half height of the ellipse around the ring centre that fits the near ends best
  off max  largest distance of a near end from that ellipse, along its ray from the centre, px; off rms likewise
  top, bottom  px between the near ends of the highest label of the left and of the right half (and the lowest):
           the gap the labels leave open at 12 and at 6 o'clock

usage: arc_shape.py NAME=FILE.labels.json ...
"""
import json, math, sys
import numpy as np
sys.path.insert(0, '/home/kawato/gbdraw-baselines/label-placement-20261011/scripts')
import lp_metrics
from gbdraw.labels import circular as lc

print(f'{"render":<30}{"A":>7}{"B":>7}{"off max":>9}{"off rms":>9}{"top":>7}{"bottom":>8}')
for spec in sys.argv[1:]:
    name, path = spec.split('=', 1)
    d = json.load(open(path)); ls = [l for l in d['labels'] if not l.get('is_embedded')]
    R = max(float(l['feature_outer_radius_px']) for l in ls)
    ends = []
    for l in ls:
        x0, y0, x1, y1 = lp_metrics.label_box(l, d['total_length'], lc)
        x0, x1, y0, y1 = min(x0, x1), max(x0, x1), min(y0, y1), max(y0, y1)
        if max(abs(x0), abs(x1)) ** 2 + max(abs(y0), abs(y1)) ** 2 > R * R:      # outside labels only
            ends.append((x0 if x0 + x1 > 0 else x1, 0.5 * (y0 + y1)))
    p = np.array(ends)
    (u, v), *_ = np.linalg.lstsq(p ** 2, np.ones(len(p)), rcond=None)
    off = np.hypot(p[:, 0], p[:, 1]) * (1.0 - 1.0 / np.sqrt(u * p[:, 0] ** 2 + v * p[:, 1] ** 2)) if u > 0 and v > 0 else np.full(len(p), np.nan)
    gap = []
    for pick in (np.argmin, np.argmax):
        right, left = p[p[:, 0] > 0], p[p[:, 0] <= 0]
        gap.append(right[pick(right[:, 1]), 0] - left[pick(left[:, 1]), 0])
    A, B = (1 / math.sqrt(u), 1 / math.sqrt(v)) if u > 0 and v > 0 else (float('nan'),) * 2
    print(f'{name:<30}{A:>7.1f}{B:>7.1f}{np.max(np.abs(off)):>9.1f}{math.sqrt(float(np.mean(off ** 2))):>9.1f}{gap[0]:>7.0f}{gap[1]:>8.0f}')
