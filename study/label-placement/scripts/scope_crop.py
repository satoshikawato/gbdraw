#!/usr/bin/env python3
"""Zoomed crop of a scope_run.py render with its defects marked.

usage: scope_crop.py OUT.png ZOOM DIR/NAME KIND [FILTER] [--region X0 Y0 X1 Y1] [--max N] [--all-boxes]
KIND: overlap | crossing | through. FILTER: a substring of the kind-pair key (for example 'c-inner|c-inner').
Red boxes: overlapping labels; red circles: leader crossings; orange box + thick line: a leader through a label.
Thin blue outlines (--all-boxes): every label box in the crop. Coordinates are SVG canvas px. Without --region the crop
is the bounding box of the first N marked defects plus a margin.
"""
import argparse, io, json, re, xml.etree.ElementTree as ET
import cairosvg
from PIL import Image, ImageDraw

Image.MAX_IMAGE_PIXELS = None
ap = argparse.ArgumentParser()
ap.add_argument('out'); ap.add_argument('zoom', type=float); ap.add_argument('base'); ap.add_argument('kind')
ap.add_argument('filter', nargs='?', default='')
ap.add_argument('--region', nargs=4, type=float); ap.add_argument('--max', type=int, default=6)
ap.add_argument('--all-boxes', action='store_true'); ap.add_argument('--pad', type=float, default=40)
ns = ap.parse_args()
g = json.load(open(ns.base + '.scope-geom.json'))
L = g['labels']
marks = []
if ns.kind == 'overlap':
    for key, i, j, d in g['overlaps']:
        if ns.filter in key:
            marks.append(('poly', L[i]['poly'])); marks.append(('poly', L[j]['poly']))
elif ns.kind == 'crossing':
    for key, i, j, pt in g['crossings']:
        if ns.filter in key:
            marks.append(('pt', pt)); marks.append(('lead', L[i]['lead'])); marks.append(('lead', L[j]['lead']))
elif ns.kind == 'through':
    for key, i, j in g['through']:
        if ns.filter in key:
            marks.append(('opoly', L[j]['poly'])); marks.append(('olead', L[i]['lead']))
per = {'overlap': 2, 'crossing': 3, 'through': 2}[ns.kind]
marks = marks[:ns.max * per]
if ns.region:
    box = ns.region
else:
    pts = [p for t, geo in marks for p in (geo if t != 'pt' else [geo])]
    box = (min(p[0] for p in pts) - ns.pad, min(p[1] for p in pts) - ns.pad, max(p[0] for p in pts) + ns.pad, max(p[1] for p in pts) + ns.pad)
x0, y0, x1, y1 = box
ET.register_namespace('', 'http://www.w3.org/2000/svg'); ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')
root = ET.parse(ns.base + '.svg').getroot()
root.set('viewBox', f'{x0} {y0} {x1 - x0} {y1 - y0}')
root.set('width', f'{(x1 - x0)}px'); root.set('height', f'{(y1 - y0)}px')
png = cairosvg.svg2png(bytestring=ET.tostring(root), background_color='white', scale=ns.zoom)
im = Image.open(io.BytesIO(png)).convert('RGB')
dr = ImageDraw.Draw(im, 'RGBA')
X = lambda x: (x - x0) * ns.zoom
Y = lambda y: (y - y0) * ns.zoom
if ns.all_boxes:
    for lab in L:
        p = lab['poly']
        if all(x0 - 50 < q[0] < x1 + 50 and y0 - 50 < q[1] < y1 + 50 for q in p):
            dr.polygon([(X(a), Y(b)) for a, b in p], outline=(40, 90, 220, 160))
for t, geo in marks:
    if t in ('poly', 'opoly'):
        dr.polygon([(X(a), Y(b)) for a, b in geo], outline=(220, 0, 0, 230) if t == 'poly' else (240, 130, 0, 230), width=2)
    elif t in ('lead', 'olead'):
        dr.line([(X(a), Y(b)) for a, b in geo], fill=(220, 0, 0, 120) if t == 'lead' else (240, 130, 0, 200), width=3)
    else:
        r = 8
        dr.ellipse([X(geo[0]) - r, Y(geo[1]) - r, X(geo[0]) + r, Y(geo[1]) + r], outline=(220, 0, 0, 240), width=3)
im.save(ns.out)
print(ns.out, im.size, len(marks), 'marks', [round(v) for v in box])
