#!/usr/bin/env python3
"""Labels and leaders against everything else drawn (LP-D06, LP-D10): ticks, tick labels, tracks, the ring.

The other objects come from a second render of the same command with `--labels none` (`<case>__bare.svg`), aligned on
the map centre. Only the map itself is used (out to the first empty ring beyond the axis), so the legend and a plot
title are not checked. Adds to each metrics.json:
  labels_on_objects    text boxes that touch another object (list in labels_on_objects_list)
  labels_near_objects  text boxes closer than `clear` px (default 3) to another object
  leaders_on_objects   leaders that touch another object more than 8 px away from the feature band

usage: lp_objects.py DIR [DIR ...]   (every <case>__<design>.metrics.json whose <case>__bare.svg is found in a DIR)
"""
import json, math, sys
from pathlib import Path
import numpy as np
from PIL import Image, ImageFilter

T = Path('/home/kawato/gbdraw-baselines/label-placement-20261011')
sys.path.insert(0, str(T / 'page')); sys.path.insert(0, str(T / 'scripts'))
import figs
import lp_metrics
from gbdraw.labels import circular as lc

CLEAR = 3


def bare_mask(stem, r_axis):
    panel = figs.load(stem)
    rgb = np.asarray(panel['img'].convert('RGB')).astype(int)
    mask = (rgb.min(axis=2) < 235)
    cx, cy = panel['center']
    yy, xx = np.indices(mask.shape)
    rr = np.hypot(xx - cx, yy - cy)
    profile = np.bincount(rr[mask].astype(int), minlength=int(rr.max()) + 2)
    r_end = int(r_axis)
    while r_end + 12 < len(profile) and profile[r_end:r_end + 12].any():
        r_end += 1
    mask &= rr <= r_end + 1
    near = np.asarray(Image.fromarray((mask * 255).astype('uint8')).filter(ImageFilter.MaxFilter(2 * CLEAR + 1))) > 0
    return mask, near, (cx, cy), r_end


def check(mj, bare):
    lj = Path(str(mj)[:-len('.metrics.json')] + '.labels.json')
    data = json.loads(lj.read_text()); met = json.loads(mj.read_text())
    ext = [l for l in data['labels'] if not l.get('is_embedded')]
    if not ext:
        return None
    r_out = max(float(l.get('feature_outer_radius_px', 0.0)) for l in ext)
    r_in = min(float(l.get('feature_inner_radius_px', r_out)) for l in ext)
    mask, near, (cx, cy), r_end = bare
    h, w = mask.shape
    on, close, lead = [], [], []
    for lab in ext:
        if abs(float(lab.get('rotation_deg') or 0.0)) > 1e-6:
            continue
        x0, y0, x1, y1 = lp_metrics.label_box(lab, data['total_length'], lc)
        a, b, c, d = int(math.floor(x0 + cx + 0.5)), int(math.floor(y0 + cy + 0.5)), int(math.ceil(x1 + cx - 0.5)), int(math.ceil(y1 + cy - 0.5))
        a, b, c, d = max(a, 0), max(b, 0), min(c, w), min(d, h)     # a box off the bare map's image touches nothing
        inside = c > a and d > b
        if inside and mask[b:d, a:c].any():
            on.append(str(lab.get('label_text', '')))
        elif inside and near[b:d, a:c].any():
            close.append(str(lab.get('label_text', '')))
        pts = lp_metrics.label_leader(lab)
        hit = False
        for p, q in zip(pts, pts[1:]):
            n = max(2, int(math.dist(p, q) / 0.5))
            for k in range(n + 1):
                x = p[0] + (q[0] - p[0]) * k / n; y = p[1] + (q[1] - p[1]) * k / n
                r = math.hypot(x, y)
                if r_in - 8 <= r <= r_out + 8:
                    continue
                xi, yi = int(round(x + cx)), int(round(y + cy))
                if 0 <= xi < w and 0 <= yi < h and mask[yi, xi]:
                    hit = True; break
            if hit:
                break
        if hit:
            lead.append(str(lab.get('label_text', '')))
    met.update(labels_on_objects=len(on), labels_on_objects_list=on[:40], labels_near_objects=len(close),
               labels_near_objects_list=close[:40], leaders_on_objects=len(lead), leaders_on_objects_list=lead[:40])
    mj.write_text(json.dumps(met, indent=1))
    return len(on), len(close), len(lead)


if __name__ == '__main__':
    dirs = [Path(d) for d in sys.argv[1:]]
    bares = {}
    for d in dirs:
        for svg in d.glob('*__bare.svg'):
            bares.setdefault(svg.name[:-len('__bare.svg')], svg.with_suffix(''))
    cache, done = {}, 0
    for d in dirs:
        for mj in sorted(d.glob('*.metrics.json')):
            case = mj.name.split('__')[0]
            if case not in bares:
                continue
            labels = json.loads(Path(str(mj)[:-len('.metrics.json')] + '.labels.json').read_text())['labels']
            ext = [l for l in labels if not l.get('is_embedded')]
            if not ext:
                continue
            if case not in cache:
                cache[case] = bare_mask(bares[case], max(float(l.get('feature_outer_radius_px', 0.0)) for l in ext))
            if check(mj, cache[case]) is not None:
                done += 1
    print('checked', done, 'renders against', len(cache), 'bare maps')
