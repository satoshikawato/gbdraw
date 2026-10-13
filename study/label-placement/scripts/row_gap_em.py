#!/usr/bin/env python3
"""LP-D29: how far apart are stacked label rows, in px and relative to the font size?

Read from the SVG itself (works for any render, also a published Gallery figure): every horizontal <text> of the
`label_text` group gets its ink box (scripts/lp_ports.py: what Inkscape shows as the bounding box). A pair is
"stacked" when the two ink boxes overlap in x (or are less than `--side` px apart in x) and the vertical distance
between them is under `--reach` em. Reported per file: number of stacked pairs, smallest and median vertical ink
gap in px and in em (px / font size), and the pairs under `--tight` em.

usage: row_gap_em.py [-v] [--tight EM] [--reach EM] NAME=FILE.svg ...
"""
import os, re, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lp_ports as P


def texts(svg):
    t = open(svg, encoding='utf-8').read()
    i = t.index('id="label_text"'); s = t.rfind('<g', 0, i)
    depth, k = 0, s
    for m in re.finditer(r'<(/?)g\b', t[s:]):
        depth += -1 if m.group(1) else 1
        if depth == 0:
            k = s + m.start(); break
    out = []
    for m in re.finditer(r'<text\b([^>]*)>(.*?)</text>', t[s:k], re.S):
        a = dict(re.findall(r'([\w:-]+)="([^"]*)"', m.group(1)))
        rot = re.search(r'rotate\(\s*([-0-9.eE]+)', a.get('transform', ''))
        if rot and abs(float(rot.group(1))) > 1e-6:
            continue
        txt = re.sub(r'<[^>]+>', '', m.group(2))
        size = float(a.get('font-size', '14').replace('px', ''))
        try:
            box = P.ink_box(txt, float(a['x']), float(a['y']), size, a.get('text-anchor', 'start'))
        except KeyError:
            continue
        out.append((txt, size, box))
    return out


def pairs(labs, reach, side):
    res = []
    for i, (ta, sa, a) in enumerate(labs):
        for tb, sb, b in labs[i + 1:]:
            em = max(sa, sb)
            dx = max(a[0] - b[2], b[0] - a[2])
            dy = max(a[1] - b[3], b[1] - a[3])
            if dx < side and 0 <= dy < reach * em:
                res.append((dy, dy / em, ta, tb))
    return sorted(res)


if __name__ == '__main__':
    a = sys.argv[1:]; verbose = '-v' in a; a = [v for v in a if v != '-v']
    tight, reach, side = 0.2, 0.6, 0.0
    while a and a[0].startswith('--'):
        key, val = a[0], float(a[1]); a = a[2:]
        tight, reach, side = (val, reach, side) if key == '--tight' else (tight, val, side) if key == '--reach' else (tight, reach, val)
    for spec in a:
        name, _, svg = spec.rpartition('=')
        labs = texts(svg); ps = pairs(labs, reach, side)
        if not ps:
            print(f'{name or svg}: {len(labs)} horizontal labels, no stacked pair within {reach} em'); continue
        size = statistics.median(s for _, s, _ in labs)
        n_tight = sum(1 for p in ps if p[1] < tight)
        print(f'{name or svg}: {len(labs)} horizontal labels at {size:g} pt, stacked pairs {len(ps)}, ink gap min {ps[0][0]:.1f} px = {ps[0][1]:.2f} em, '
              f'median {statistics.median(p[0] for p in ps):.1f} px = {statistics.median(p[1] for p in ps):.2f} em, under {tight} em: {n_tight}')
        if verbose:
            for dy, em, ta, tb in ps:
                print(f'    {dy:5.1f} px {em:.2f} em  {ta} / {tb}')
