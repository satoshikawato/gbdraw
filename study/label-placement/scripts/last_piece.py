#!/usr/bin/env python3
"""LP-D24: how long is the last piece of a bent leader (from its last bend to its end)?

The Owner on tRNA-Trp: the horizontal piece after the bend may be shorter (22 px was too long). A leader is bent
from `--bent` degrees of turn (default 3.0); pieces under 3 px are dropped as in lp_metrics.leader_turns.

usage: last_piece.py [-v] [--bent DEG] NAME=DIR:SUFFIX ... | FILE.labels.json ...
"""
import glob, json, math, os, sys


def rows(path, bent=3.0):
    out = []
    for lab in json.load(open(path))['labels']:
        pts = lab.get('leader_points')
        if lab.get('is_embedded') or not pts or len(pts) < 3:
            continue
        kept = [p for k, p in enumerate(pts) if k == 0 or math.dist(p, pts[k - 1]) >= 3.0]
        if kept[-1] != pts[-1]:
            kept[-1:] = [pts[-1]] if len(kept) > 1 else [kept[0], pts[-1]]
        turn = 0.0
        for a, b, c in zip(kept, kept[1:], kept[2:]):
            u = (b[0] - a[0], b[1] - a[1]); v = (c[0] - b[0], c[1] - b[1]); nu, nv = math.hypot(*u), math.hypot(*v)
            if nu >= 3.0 and nv >= 3.0:
                turn += math.degrees(math.acos(max(-1.0, min(1.0, (u[0] * v[0] + u[1] * v[1]) / (nu * nv)))))
        if turn >= bent and len(kept) >= 3:
            out.append((lab['label_text'], math.dist(kept[-2], kept[-1]), turn))
    return out


if __name__ == '__main__':
    args = sys.argv[1:]; verbose = '-v' in args; bent = 3.0
    if '--bent' in args:
        bent = float(args[args.index('--bent') + 1]); del args[args.index('--bent'):args.index('--bent') + 2]
    for spec in (a for a in args if a != '-v'):
        if '=' in spec:
            name, rest = spec.split('=', 1); folder, suffix = rest.rsplit(':', 1)
            files = sorted(glob.glob(os.path.join(folder, f'*__{suffix}.labels.json')))
        else:
            name, files = os.path.basename(spec), [spec]
        allr = []
        for f in files:
            r = rows(f, bent); allr += r
            if verbose:
                for t, n, turn in sorted(r, key=lambda x: -x[1])[:6]:
                    print(f"  {os.path.basename(f).split('__')[0]}: {t} last piece {n:.1f} px, turn {turn:.1f}")
        ls = sorted(n for _, n, _ in allr)
        if ls:
            print(f'{name}: {len(ls)} bent leaders; last piece median {ls[len(ls) // 2]:.1f} px, max {ls[-1]:.1f}, over 15 px {sum(v > 15 for v in ls)}, under 6 px {sum(v < 6 for v in ls)}')
        else:
            print(f'{name}: no bent leaders')
