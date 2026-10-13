#!/usr/bin/env python3
"""Information only (not in DEF-CROSS): Circular labels whose box overlaps a line of the centre definition text by more
than 0.5 px in x and y. Definition lines are measured with the engine's text measure (bold not widened).
usage: PYTHONPATH=<worktree> python3.13 -P scope_definition.py DIR [PATTERN]"""
import glob, json, pathlib, re, sys
import xml.etree.ElementTree as ET
from collections import Counter
from gbdraw.core.text import calculate_bbox_dimensions
TR = re.compile(r"translate\(\s*([-0-9.eE]+)[ ,]+([-0-9.eE]+)\s*\)")
d = pathlib.Path(sys.argv[1]); pat = sys.argv[2] if len(sys.argv) > 2 else '[cr]_*'
for f in sorted(glob.glob(str(d / f'{pat}.scope-geom.json'))):
    svg = f.replace('.scope-geom.json', '.svg')
    root = ET.parse(svg).getroot()
    boxes = []
    for g in root.iter():
        if not g.get('id', '').endswith('_definition'):
            continue
        cx = sum(float(a) for a, _ in TR.findall(g.get('transform', ''))); cy = sum(float(b) for _, b in TR.findall(g.get('transform', '')))
        for t in g.iter():
            if t.tag.split('}')[1] != 'text':
                continue
            text = ''.join(t.itertext()).strip()
            if not text:
                continue
            fam = t.get('font-family', 'Liberation Sans').split(',')[0].strip("' ")
            w, h = calculate_bbox_dimensions(text, fam, float(t.get('font-size', '14')), 96)
            x = cx + float(t.get('x', 0)); y = cy + float(t.get('y', 0))
            boxes.append((x - w / 2, y - h / 2, x + w / 2, y + h / 2, text))
    g = json.load(open(f))
    hit = Counter(); texts = []
    for lab in g['labels']:
        xs = [p[0] for p in lab['poly']]; ys = [p[1] for p in lab['poly']]
        for b in boxes:
            if min(max(xs), b[2]) - max(min(xs), b[0]) > 0.5 and min(max(ys), b[3]) - max(min(ys), b[1]) > 0.5:
                hit[lab['kind']] += 1; texts.append(f"{lab['text']} / {b[4]}"); break
    m = json.load(open(f.replace('.scope-geom.json', '.scope.json')))
    m['labels_over_definition'] = dict(hit)
    json.dump(m, open(f.replace('.scope-geom.json', '.scope.json'), 'w'), indent=1)
    if hit:
        print(pathlib.Path(f).name.replace('.scope-geom.json', ''), dict(hit), texts[:4])
