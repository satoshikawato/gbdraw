#!/usr/bin/env python3
"""Fit study table (LP-D41): smallest single change that makes today's placement fit, per definition.

Usage: fit_table.py [DIR=work/FIT]
D1: no overlapping label boxes (text_overlap_pairs 0). D2: D1 and no leader crossing. D3: D2 and no leader through text.
Prints, per case: today's numbers, then per axis (font, y, xy, spacing, font+y) the first value meeting D1 / D2 / D3,
with the canvas size there. A value later in the axis that fails again is flagged with '!'.
"""
import glob, json, os, re, sys

d = sys.argv[1] if len(sys.argv) > 1 else 'work/FIT'
rows = {}
for f in glob.glob(f'{d}/*.metrics.json'):
    case, var = os.path.basename(f)[:-len('.metrics.json')].split('__', 1)
    rows.setdefault(case, {})[var] = json.load(open(f))
DEF = {'D1': lambda m: m['text_overlap_pairs'] == 0,
       'D2': lambda m: m['text_overlap_pairs'] == 0 and m['leader_crossing_pairs'] == 0,
       'D3': lambda m: m['text_overlap_pairs'] == 0 and m['leader_crossing_pairs'] == 0 and m['leader_through_text'] == 0}


def axis(vars_, prefix):
    out = []
    for v in vars_:
        m = re.fullmatch(prefix + r'([0-9.]+)', v)
        if m:
            out.append((float(m.group(1)), v))
    return out


for case in sorted(rows):
    r = rows[case]
    cur = r.get('cur')
    if not cur:
        continue
    fsz = cur['font_size'][0] if cur.get('font_size') else None
    print(f"\n{case}: {cur['labels_external']} outer labels at {fsz} pt; overlap {cur['text_overlap_pairs']} "
          f"(overprint {cur['text_overprint_pairs']}), crossings {cur['leader_crossing_pairs']}, through text "
          f"{cur['leader_through_text']}, canvas {cur.get('canvas_px')}")
    for name, prefix, rev in (('font', 'f', True), ('x', 'x', False), ('y', 'y', False), ('xy', 'xy', False), ('spacing', 's', True)):
        ax = sorted(axis(r, prefix), reverse=rev)
        if not ax:
            continue
        parts = []
        for dn, ok in DEF.items():
            first = next((i for i, (_v, var) in enumerate(ax) if ok(r[var])), None)
            if first is None:
                parts.append(f'{dn}: none up to {ax[-1][1]}')
                continue
            again = any(not ok(r[var]) for _v, var in ax[first + 1:])
            m = r[ax[first][1]]
            parts.append(f"{dn}: {ax[first][1]}{'!' if again else ''} {m.get('canvas_px')}")
        print(f'  {name:8s}', ' | '.join(parts))
    combos = sorted((v for v in r if '+' in v), key=lambda v: (-float(v.split('+')[0][1:]), float(v.split('+')[1][1:])))
    if combos:
        fits = [v for v in combos if DEF['D2'](r[v])]
        print('  font+y  D2 fits:', ', '.join(f"{v} {r[v].get('canvas_px')}" for v in fits) or 'none')
