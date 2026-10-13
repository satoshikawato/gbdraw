#!/usr/bin/env python3
"""Lever table of the scope study: per failing case, one change at a time.
Columns: inner overlap pairs / inner crossing pairs / outer overlaps / outer crossings / inner labels over the
definition text / placement s / canvas. A row 'fits' (F) when it has no overlap and no crossing of any kind."""
import json, pathlib, re, sys
T = pathlib.Path(__file__).resolve().parent.parent
D = T / 'work' / 'SCOPE'
cases = ['c_t5_both_def', 'c_t5_both_f14', 'c_tobacco_both_f14', 'c_tobacco_both_def', 'c_p1_both_f14', 'c_plasmid37_both_def',
         'c_epifagus_both_def', 'c_mjenmv_both_def', 'c_lambda_both_def', 'c_wssv_both_def', 'c_tobacco_mid_both_def']
order = ['', 'f13', 'f12', 'f11', 'f10', 'f9', 'f8', 'f7.5', 'f7', 'f6.5', 'f6', 'ixy0.9', 'ixy0.95', 'ixy1.05', 'ixy1.1', 'ixy1.2',
         'iy1.1', 'ix1.1', 's0.5', 's1', 's2', 's5', 'out', 'radial', 'mid']


def get(m, kind, field):
    tot = 0
    for k, v in m.get(field, {}).items():
        a, b, *_ = k.split('|')
        if kind in (a, b):
            tot += v
    return tot


for c in cases:
    print(f'== {c}')
    print(f'  {"variant":8} {"in-ov":>6} {"in-x":>6} {"out-ov":>6} {"out-x":>6} {"r-x":>5} {"def":>4} {"place s":>8} canvas  fit')
    for v in order:
        f = D / (c + (f'__{v}' if v else '') + '.scope.json')
        if not f.exists():
            continue
        m = json.loads(f.read_text())
        inov, inx = get(m, 'c-inner', 'overlap_pairs'), get(m, 'c-inner', 'crossing_pairs')
        oov, ox = get(m, 'c-outer', 'overlap_pairs'), get(m, 'c-outer', 'crossing_pairs')
        rx = get(m, 'c-radial-inner', 'crossing_pairs') + get(m, 'c-radial-outer', 'crossing_pairs') + get(m, 'c-radial-inner', 'overlap_pairs') + get(m, 'c-radial-outer', 'overlap_pairs')
        de = sum(m.get('labels_over_definition', {}).values())
        t = m['timers'].get('circular_label_placement_s', 0)
        fit = 'F' if m['overlap_total'] == 0 and m['crossing_total'] == 0 else ('x-only' if m['overlap_total'] == 0 else '')
        print(f'  {v or "cur":8} {inov:6} {inx:6} {oov:6} {ox:6} {rx:5} {de:4} {t:8.1f} {m["canvas_px"]} {fit}')
