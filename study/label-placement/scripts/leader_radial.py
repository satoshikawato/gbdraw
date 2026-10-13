#!/usr/bin/env python3
"""LP-D40: does every leader cross the axis at a right angle?

The Owner (2026-10-13): 「リーダー線は一応axisとは直交してほしい。…直交したうえで曲がる分には構わんから。」 For every outside label the
first piece of the leader (feature anchor -> first vertex, or the whole leader when it is one line) is compared with the
outward radius at the anchor:
  n        outside labels with a leader
  dev_max  largest angle between the first piece and the radius (deg; 0 = perpendicular to the axis)
  n>2      leaders whose first piece leans more than 2 deg      (a 2 deg lean is not visible)
  n>5      ... more than 5 deg
  first    shortest first piece measured from the outer edge of the feature's track (px; <0: the bend lies inside the track)
  worst    the three labels with the largest lean

usage: leader_radial.py NAME=FILE.labels.json ...
"""
import json, math, sys


def report(name, path):
    d = json.load(open(path))
    rows = []
    for l in d['labels']:
        if l.get('is_embedded'):
            continue
        pts = l.get('leader_points')
        if not pts:                                         # today's engine: anchor -> elbow (middle) -> text start
            pts = [(l['feature_anchor_x'], l['feature_anchor_y']), (l['middle_x'], l['middle_y']), (l['leader_start_x'], l['leader_start_y'])]
        if len(pts) < 2:
            continue
        (ax, ay), (bx, by) = pts[0], pts[1]
        ra = math.hypot(ax, ay)
        dx, dy, L = bx - ax, by - ay, math.hypot(bx - ax, by - ay)
        if ra < 1e-9 or L < 1e-9:
            continue
        cosv = max(-1.0, min(1.0, (dx * ax + dy * ay) / (L * ra)))
        dev = math.degrees(math.acos(cosv))
        outer = l.get('feature_outer_radius_px') or ra
        rows.append((dev, math.hypot(bx, by) - outer, l['label_text']))
    rows.sort(reverse=True)
    n = len(rows)
    worst = ', '.join(f'{t} {v:.1f}' for v, _, t in rows[:3])
    print(f"{name:34} n={n:4d} dev_max={rows[0][0] if rows else 0:5.1f} n>2={sum(v > 2 for v, _, _ in rows):3d} "
          f"n>5={sum(v > 5 for v, _, _ in rows):3d} first={min(f for _, f, _ in rows) if rows else 0:6.1f}  worst: {worst}")


if __name__ == '__main__':
    for arg in sys.argv[1:]:
        name, path = arg.split('=', 1)
        report(name, path)
