#!/usr/bin/env python3
"""LP-D23: does every leader point at one of the 12 contact points of its label's bounding box?

The box is the ink box of the text (scripts/lp_ports.py: what Inkscape shows as the bounding box). A leader points
at a contact point when the line of its last piece passes within `--tol` px (default 0.75) of it and the point lies
0 to `--reach` px (default 5) ahead of the leader's end. Rotated and embedded labels are skipped.
Also reported: the arrival rule (which points a direction may use; round 2, after the Owner's picture): within 45
degrees of a short side's normal only the middle of that side; a vertex between 45 and 90 degrees, from outside
that corner (round Q3: up to 90, so that a steep leader near a pole may end at a vertex as in J2; the box then
stands beside its leader and a staircase of wide labels is possible); a long-side point within 45 degrees of
its normal. `--rule 2` is the rule of P and Q rounds 1 and 2 (vertex 45 to 70), `--rule 1` the first one (20 to 70).

usage: port_check.py [-v] [--tol PX] [--reach PX] NAME=DIR:SUFFIX ... | FILE.labels.json ...
"""
import glob, json, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lp_ports as P

VERTEX_FROM = 45.0
VERTEX_TO = 90.0
NORMAL = {'L': (-1, 0), 'R': (1, 0), 'T1': (0, -1), 'T2': (0, -1), 'T3': (0, -1), 'B1': (0, 1), 'B2': (0, 1), 'B3': (0, 1),
          'TL': (-1, -1), 'TR': (1, -1), 'BL': (-1, 1), 'BR': (1, 1)}


def arrival_ok(name, d):
    """d = unit direction of travel of the last piece; the leader must come from outside, against the normal."""
    nx, ny = NORMAL[name]; ax, ay = -d[0], -d[1]                # where the leader comes from, seen from the point
    if len(name) == 2 and name[0] in 'TB' and name[1] in 'LR':  # vertex: from outside the corner, 20 to 70 degrees
        if ax * nx <= 0 or ay * ny <= 0:
            return False
        deg = math.degrees(math.atan2(abs(ay), abs(ax)))      # 0 = along the normal of the short side (horizontal text)
        return VERTEX_FROM - 0.5 <= deg <= VERTEX_TO + 0.5
    n = math.hypot(nx, ny)
    return (ax * nx + ay * ny) / n >= math.cos(math.radians(45.5))


def check(path, tol=0.75, reach=5.0):
    rows = []
    for lab in json.load(open(path))['labels']:
        pts = lab.get('leader_points')
        if lab.get('is_embedded') or not pts or len(pts) < 2 or abs(float(lab.get('rotation_deg', 0.0))) > 1e-6:
            continue
        kept = [p for k, p in enumerate(pts) if k == 0 or math.dist(p, pts[k - 1]) >= 3.0]
        a, e = (kept[-2], pts[-1]) if len(kept) >= 2 else (pts[0], pts[-1])
        n = math.dist(a, e)
        if n < 1e-9:
            continue
        d = ((e[0] - a[0]) / n, (e[1] - a[1]) / n)
        box = P.ink_box(lab['label_text'], float(lab['text_x']), float(lab['text_y']), float(lab.get('font_size', 14.0)), lab.get('text_anchor', 'start'))
        best = None
        for name, (px, py) in P.ports(box).items():
            ahead = (px - e[0]) * d[0] + (py - e[1]) * d[1]; off = abs((px - e[0]) * d[1] - (py - e[1]) * d[0])
            if -0.5 <= ahead <= reach and off <= tol and (best is None or off < best[1]):
                best = (name, off, ahead)
        rows.append((lab['label_text'], best[0] if best else None, arrival_ok(best[0], d) if best else None, round(best[2], 1) if best else None))
    return rows


if __name__ == '__main__':
    args = sys.argv[1:]; verbose = '-v' in args; opt = {'--tol': 0.75, '--reach': 5.0}
    if '--rule' in args:
        rule = args[args.index('--rule') + 1]; del args[args.index('--rule'):args.index('--rule') + 2]
        VERTEX_FROM = 20.0 if rule == '1' else 45.0; VERTEX_TO = 70.0 if rule in ('1', '2') else 90.0
    for k in list(opt):
        if k in args:
            opt[k] = float(args[args.index(k) + 1]); del args[args.index(k):args.index(k) + 2]
    for spec in (a for a in args if a != '-v'):
        if '=' in spec:
            name, rest = spec.split('=', 1); folder, suffix = rest.rsplit(':', 1)
            files = sorted(glob.glob(os.path.join(folder, f'*__{suffix}.labels.json')))
        else:
            name, files = os.path.basename(spec), [spec]
        total = off = wrong = 0; use = {}
        for f in files:
            rows = check(f, opt['--tol'], opt['--reach'])
            total += len(rows); off += sum(r[1] is None for r in rows); wrong += sum(r[2] is False for r in rows)
            for r in rows:
                use[r[1] or 'none'] = use.get(r[1] or 'none', 0) + 1
            if verbose:
                bad = [f'{r[0]} ({r[1] or "no point"})' for r in rows if r[1] is None or r[2] is False]
                if bad:
                    print(f"  {os.path.basename(f).split('__')[0]}: " + '; '.join(bad[:12]))
        print(f'{name}: {total} leaders, off a contact point {off}, at a point its direction may not use {wrong}; ' + ' '.join(f'{k} {v}' for k, v in sorted(use.items(), key=lambda kv: -kv[1])))
