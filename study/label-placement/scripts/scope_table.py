#!/usr/bin/env python3
"""Table of the scope study renders: label counts per kind, overlapping pairs, crossing leader pairs, leaders through
another label's box (information only), label placement seconds, canvas size.
usage: scope_table.py [DIR=work/SCOPE] [PATTERN=*]"""
import glob, json, pathlib, sys

T = pathlib.Path(__file__).resolve().parent.parent
d = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else T / 'work' / 'SCOPE'
pat = sys.argv[2] if len(sys.argv) > 2 else '*'
AB = {'c-outer': 'out', 'c-inner': 'in', 'c-embedded': 'emb', 'c-radial-outer': 'r-out', 'c-radial-inner': 'r-in',
      'l-external': 'rows', 'l-embedded': 'emb', 'l-above': 'above'}


def fmt(counter, sep='|'):
    out = []
    for k, v in sorted(counter.items()):
        parts = k.split(sep) if sep in k else k.split('>')
        kinds = '/'.join(AB.get(p, p) for p in parts[:2])
        if sep == '|' and parts[2] == 'cross':
            kinds += '(x-rec)'
        out.append(f'{kinds} {v}')
    return ', '.join(out) or '0'


print(f'{"case":34} {"labels by kind":34} {"overlapping pairs":30} {"crossing leader pairs":34} {"leader through label":26} {"place s":>7} canvas')
for f in sorted(glob.glob(str(d / f'{pat}.scope.json'))):
    m = json.load(open(f))
    name = pathlib.Path(f).name[:-len('.scope.json')]
    kinds = ' '.join(f'{AB[k]} {v}' for k, v in sorted(m['labels_by_kind'].items()))
    t = m['timers']
    ps = t.get('circular_label_placement_s', t.get('linear_label_placement_s', 0.0))
    print(f'{name:34} {kinds:34} {fmt(m["overlap_pairs"]):30} {fmt(m["crossing_pairs"]):34} {fmt(m["through_text"], ">"):26} {ps:7.2f} {m["canvas_px"]}'
          + ('' if not m.get('unmatched_labels') else f' unmatched={m["unmatched_labels"]}'))
for f in sorted(glob.glob(str(d / f'{pat}.failed.json'))):
    m = json.load(open(f))
    print(f'{pathlib.Path(f).name[:-len(".failed.json")]:34} FAILED rc={m["rc"]}: {m["error"][:160]}')
