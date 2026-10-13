#!/usr/bin/env python3
"""LP-D34: how far do the labels reach above the top and below the bottom of the map?

Per render: R = outer radius of the outermost labelled feature track; top / bottom = px by which the outside label
boxes pass R upward / downward (the "limit" today's placement keeps without naming it); side = the same sideways;
rows = outside labels on the right / left half; need = px of column height the fuller half needs at the row rule
(box height + max(--rowmin, --rowem * font size) ), have = 2 * (R + larger of top and bottom).

usage: pole_extent.py [--rowem 0.22] [--rowmin 4] NAME=FILE.labels.json ...
"""
import json, sys
sys.path.insert(0, '/home/kawato/gbdraw-baselines/label-placement-20261011/scripts')
import lp_metrics
from gbdraw.labels import circular as lc
a = sys.argv[1:]; ROWEM, ROWMIN = 0.22, 4.0
while a and a[0].startswith('--'):
    if a[0] == '--rowem': ROWEM = float(a[1])
    if a[0] == '--rowmin': ROWMIN = float(a[1])
    a = a[2:]
print(f'{"render":<34}{"R":>6}{"top":>6}{"bottom":>8}{"side":>6}{"rows R/L":>10}{"box h":>7}{"need":>7}{"have":>7}')
for spec in a:
    name, path = spec.split('=', 1)
    d = json.load(open(path)); ls = [l for l in d['labels'] if not l.get('is_embedded')]
    if not ls:
        continue
    R = max(float(l['feature_outer_radius_px']) for l in ls)
    bs = []
    for l in ls:
        x0, y0, x1, y1 = lp_metrics.label_box(l, d['total_length'], lc); bs.append((min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)))
    out = [b for b in bs if max(abs(b[0]), abs(b[2])) ** 2 + max(abs(b[1]), abs(b[3])) ** 2 > R * R]      # outside labels only
    top = -min(b[1] for b in out) - R; bot = max(b[3] for b in out) - R; side = max(max(-b[0], b[2]) for b in out) - R
    nr = sum(1 for b in out if b[0] + b[2] > 0); nl = len(out) - nr
    h = sorted(b[3] - b[1] for b in out)[len(out) // 2]; fs = max(float(l['font_size']) for l in ls)
    need = max(nr, nl) * (h + max(ROWMIN, ROWEM * fs)); have = 2 * (R + max(top, bot))
    print(f'{name:<34}{R:>6.0f}{top:>6.0f}{bot:>8.0f}{side:>6.0f}{f"{nr}/{nl}":>10}{h:>7.1f}{need:>7.0f}{have:>7.0f}')
