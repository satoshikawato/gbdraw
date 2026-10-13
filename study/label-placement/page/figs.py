"""Shared figure helpers for the label-placement study: load a render, align panels on the map centre, crop, compose.

A render is the set <dir>/<case>__<design>.{svg,png,metrics.json,labels.json} written by scripts/run_case.py.
Coordinates called "local" are circle-local pixels (map centre = 0,0; y down), the same in every design, so one
local box shows the same part of the map whatever the canvas size.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

Image.MAX_IMAGE_PIXELS = None
FONT_DIR = Path('/home/kawato/gbdraw-work/gbdraw/data')
PAD = 24


def font(size, bold=False):
    return ImageFont.truetype(str(FONT_DIR / ('LiberationSans-Bold.ttf' if bold else 'LiberationSans-Regular.ttf')), size)


def find_render(dirs, case, design):
    for d in dirs:
        stem = Path(d) / f'{case}__{design}'
        if stem.with_suffix('.png').exists() and stem.with_suffix('.svg').exists():
            return stem
    return None


def center_of(svg_path):
    """Map centre in canvas pixels: the summed translate() of the Axis group."""
    head = Path(svg_path).read_text(encoding='utf-8')[:2_000_000]
    tag = re.search(r'<g\b[^>]*\bid="Axis"[^>]*>', head)
    if not tag:
        raise ValueError(f'no Axis group in {svg_path}')
    cx = cy = 0.0
    for x, y in re.findall(r'translate\(\s*([-0-9.eE]+)[ ,]+([-0-9.eE]+)\s*\)', tag.group(0)):
        cx += float(x); cy += float(y)
    return cx, cy


def rasterize(svg_path, pad=PAD):
    """PNG of the SVG on white with `pad` px added on every side (cached next to the SVG).

    gbdraw's own PNG has no margin, and cairosvg draws text a little wider than gbdraw measured it, so a label that
    ends at the canvas edge loses its first letters (OV-471). The margin keeps that export fault out of the comparison.
    """
    svg_path = Path(svg_path)
    cache = svg_path.with_suffix(f'.pad{pad}.png')
    if not cache.exists() or cache.stat().st_mtime < svg_path.stat().st_mtime:
        import cairosvg
        text = svg_path.read_text(encoding='utf-8')
        tag = re.search(r'<svg\b[^>]*>', text).group(0)
        w = float(re.search(r'\bwidth="([0-9.]+)', tag).group(1)); h = float(re.search(r'\bheight="([0-9.]+)', tag).group(1))
        new = re.sub(r'\bwidth="[^"]*"', f'width="{w + 2 * pad}px"', tag, count=1)
        new = re.sub(r'\bheight="[^"]*"', f'height="{h + 2 * pad}px"', new, count=1)
        new = re.sub(r'\bviewBox="[^"]*"', f'viewBox="{-pad} {-pad} {w + 2 * pad} {h + 2 * pad}"', new, count=1)
        cairosvg.svg2png(bytestring=text.replace(tag, new, 1).encode('utf-8'), write_to=str(cache), background_color='white')
    return Image.open(cache).convert('RGB')


def load(stem, pad=PAD):
    """Image (white ground, `pad` px margin), map centre in image pixels, metrics and drawn labels of one render."""
    stem = Path(stem)
    svg = stem.with_suffix('.svg')
    cx, cy = center_of(svg)
    read = lambda suffix: json.loads(Path(str(stem) + suffix).read_text()) if Path(str(stem) + suffix).exists() else {}
    metrics = read('.metrics.json')
    return {'stem': stem, 'img': rasterize(svg, pad), 'center': (cx + pad, cy + pad), 'metrics': metrics,
            'canvas': tuple(metrics.get('canvas_px') or ()), 'labels': read('.labels.json').get('labels', [])}


def load_case(dirs, case, designs, base='cur'):
    """Panels of one case by design key, and the designs left out because their SVG equals the base design's."""
    panels, same, base_bytes = {}, [], None
    for k in designs:
        stem = find_render(dirs, case, k)
        if not stem:
            continue
        data = stem.with_suffix('.svg').read_bytes()
        if k == base:
            base_bytes = data
        elif base_bytes is not None and data == base_bytes:
            same.append(k)
            continue
        panels[k] = load(stem)
    return panels, same


def union(panels):
    """Canvas that holds every panel with the map centres on one point: (width, height, [(x, y) offset per panel])."""
    left = max(p['center'][0] for p in panels); top = max(p['center'][1] for p in panels)
    right = max(p['img'].width - p['center'][0] for p in panels)
    bottom = max(p['img'].height - p['center'][1] for p in panels)
    return left + right, top + bottom, [(left - p['center'][0], top - p['center'][1]) for p in panels]


def crop_local(panel, box):
    x0, y0, x1, y1 = box
    cx, cy = panel['center']
    out = Image.new('RGB', (round(x1 - x0), round(y1 - y0)), 'white')
    out.paste(panel['img'], (round(-(cx + x0)), round(-(cy + y0))))
    return out


def crop_zoom(panel, box, zoom, pad=PAD):
    """The region of crop_local drawn `zoom` times larger from the SVG (to show differences of a few px)."""
    import io
    import cairosvg
    text = Path(panel['stem']).with_suffix('.svg').read_text(encoding='utf-8')
    tag = re.search(r'<svg\b[^>]*>', text).group(0)
    w = float(re.search(r'\bwidth="([0-9.]+)', tag).group(1)); h = float(re.search(r'\bheight="([0-9.]+)', tag).group(1))
    new = re.sub(r'\bwidth="[^"]*"', f'width="{w + 2 * pad}px"', tag, count=1)
    new = re.sub(r'\bheight="[^"]*"', f'height="{h + 2 * pad}px"', new, count=1)
    new = re.sub(r'\bviewBox="[^"]*"', f'viewBox="{-pad} {-pad} {w + 2 * pad} {h + 2 * pad}"', new, count=1)
    png = cairosvg.svg2png(bytestring=text.replace(tag, new, 1).encode('utf-8'), background_color='white', scale=zoom)
    big = Image.open(io.BytesIO(png)).convert('RGB')
    x0, y0, x1, y1 = box
    cx, cy = panel['center']
    out = Image.new('RGB', (round((x1 - x0) * zoom), round((y1 - y0) * zoom)), 'white')
    out.paste(big, (round(-(cx + x0) * zoom), round(-(cy + y0) * zoom)))
    return out


def label_points(labels):
    points = []
    for lab in labels:
        if lab.get('is_embedded'):
            continue
        if lab.get('aabb_local'):
            x0, y0, x1, y1 = lab['aabb_local']
            points.append(((x0 + x1) / 2, (y0 + y1) / 2))
        elif 'start_x' in lab:
            points.append((float(lab['start_x']), float(lab['start_y'])))
    return points


def densest_box(labels, w=520, h=380):
    """Local box of size w x h that holds the most external labels (centred on their mean)."""
    points = label_points(labels)
    if not points:
        return (-w / 2, -h / 2, w / 2, h / 2)
    best = max(points, key=lambda c: sum(abs(x - c[0]) <= w / 2 and abs(y - c[1]) <= h / 2 for x, y in points))
    inside = [(x, y) for x, y in points if abs(x - best[0]) <= w / 2 and abs(y - best[1]) <= h / 2]
    mx = sum(x for x, _ in inside) / len(inside); my = sum(y for _, y in inside) / len(inside)
    return (round(mx - w / 2), round(my - h / 2), round(mx + w / 2), round(my + h / 2))


def summary(panel):
    """Caption facts of one panel from its metrics: font, labels placed and dropped, canvas."""
    m = panel['metrics']
    proto = m.get('proto') or {}
    dropped = proto.get('dropped')
    dropped = len(dropped) if isinstance(dropped, list) else (dropped or 0)
    size = '/'.join(f'{v:g}' for v in m.get('font_size', [])) or '?'
    canvas = ' x '.join(str(v) for v in panel['canvas']) or '?'
    return f"{size} pt, {m.get('labels_external', '?')} labels, {dropped} dropped, {canvas}"


def summary_ja(panel):
    m = panel['metrics']
    dropped = (m.get('proto') or {}).get('dropped')
    dropped = len(dropped) if isinstance(dropped, list) else (dropped or 0)
    size = '/'.join(f'{v:g}' for v in m.get('font_size', [])) or '?'
    canvas = ' × '.join(str(v) for v in panel['canvas']) or '?'
    return f"{size} pt · ラベル {m.get('labels_external', '?')} 個 · 省略 {dropped} 個 · {canvas} px"


def grid(cells, cols, pad=14, caption_px=22):
    """Compose captioned images (list of (caption, PIL image)) into one grid on white; all cells share one size."""
    w = max(im.width for _, im in cells); h = max(im.height for _, im in cells)
    bar = caption_px + 12
    rows = -(-len(cells) // cols)
    out = Image.new('RGB', (cols * w + (cols + 1) * pad, rows * (h + bar) + (rows + 1) * pad), 'white')
    draw = ImageDraw.Draw(out)
    for i, (caption, im) in enumerate(cells):
        x = pad + (i % cols) * (w + pad); y = pad + (i // cols) * (h + bar + pad)
        head, _, rest = caption.partition('|')
        draw.text((x + 2, y), head.strip(), fill=(20, 60, 150), font=font(caption_px, bold=True))
        if rest:
            offset = draw.textlength(head.strip() + '  ', font=font(caption_px, bold=True))
            draw.text((x + 2 + offset, y + 3), rest.strip(), fill=(70, 78, 92), font=font(caption_px - 4))
        out.paste(im, (x, y + bar))
        draw.rectangle((x - 1, y + bar - 1, x + w, y + bar + h), outline=(205, 211, 221))
    return out


def aligned(panels, scale=1.0):
    """Every panel drawn on the shared canvas of union(panels), optionally scaled."""
    w, h, offsets = union(panels)
    out = []
    for p, (ox, oy) in zip(panels, offsets):
        canvas = Image.new('RGB', (round(w), round(h)), 'white')
        canvas.paste(p['img'], (round(ox), round(oy)))
        if scale != 1.0:
            canvas = canvas.resize((max(1, round(w * scale)), max(1, round(h * scale))), Image.LANCZOS)
        out.append(canvas)
    return out


MARK_COLORS = {'back': (226, 112, 0), 'far': (214, 30, 30), 'offmid': (20, 110, 220), 'unclear': (176, 38, 196), 'behind': (176, 38, 196)}


def draw_marks(im, panel, box, zoom, kinds):
    """Study marks drawn on a crop made by crop_local (zoom 1) or crop_zoom. `kinds`: any of
    'back' (labels whose box comes more than 1 px nearer to the ring than the end of their leader: box outline, a dot
    on the end, and the circle around the map centre through the end, so the part of the box inside the circle is
    the amount), 'far' (LP-D15: leaders that end on the far half of their box: ring on the end, box outline),
    'unclear' (ends nearer to another label's text than their own + 1 px: ring on the end), 'behind' (leaders that
    pass behind another label's box: the leader redrawn in colour, the box outlined). Returns the number of marks."""
    import math
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
    import lp_metrics
    labs = [l for l in panel['labels'] if not l.get('is_embedded') and l.get('aabb_local') and l.get('leader_points')
            and abs(float(l.get('rotation_deg') or 0.0)) < 1e-6]
    boxes = [tuple(float(v) for v in l['aabb_local']) for l in labs]
    leads = [[tuple(p) for p in l['leader_points']] for l in labs]
    px = lambda x, y: ((x - box[0]) * zoom, (y - box[1]) * zoom)
    inside = lambda b: b[2] > box[0] and b[0] < box[2] and b[3] > box[1] and b[1] < box[3]
    draw = ImageDraw.Draw(im)
    wide = max(2, round(0.7 * zoom + 1))
    count = 0

    def outline(b, color):
        (x0, y0), (x1, y1) = px(b[0] - 1.2, b[1] - 1.2), px(b[2] + 1.2, b[3] + 1.2)
        draw.rectangle((x0, y0, x1, y1), outline=color, width=max(1, round(zoom * 0.5)))

    def ring(point, color):
        x, y = px(*point); r = 4.2 * zoom if zoom > 1 else 5
        draw.ellipse((x - r, y - r, x + r, y + r), outline=color, width=wide)

    if 'back' in kinds:
        color = MARK_COLORS['back']
        for lab, b, lead in zip(labs, boxes, leads):
            if not inside(b):
                continue
            back = lp_metrics.back_px([lab], [b], [lead], [False])
            if not back or back[0][0] <= 1.0:
                continue
            count += 1
            outline(b, color)
            radius = math.hypot(*lead[-1])
            angles = [math.atan2(y, x) for x in (b[0], b[2]) for y in (b[1], b[3])] + [math.atan2(lead[-1][1], lead[-1][0])]
            base = angles[-1]
            rel = [(a - base + math.pi) % (2 * math.pi) - math.pi for a in angles]
            a0, a1 = base + min(rel) - 0.05, base + max(rel) + 0.05
            steps = max(8, int((a1 - a0) * radius * zoom / 6))
            for k in range(0, steps, 2):                     # dashed arc
                p = px(radius * math.cos(a0 + (a1 - a0) * k / steps), radius * math.sin(a0 + (a1 - a0) * k / steps))
                q = px(radius * math.cos(a0 + (a1 - a0) * (k + 1) / steps), radius * math.sin(a0 + (a1 - a0) * (k + 1) / steps))
                draw.line((p, q), fill=color, width=wide)
            x, y = px(*lead[-1]); r = 1.6 * zoom + 1
            draw.ellipse((x - r, y - r, x + r, y + r), fill=color)
    if 'far' in kinds:
        for lab, b, lead in zip(labs, boxes, leads):
            if inside(b) and lp_metrics.end_far_half([lab], [b], [lead], [False])[0]:
                count += 1
                outline(b, MARK_COLORS['far']); ring(lead[-1], MARK_COLORS['far'])
    if 'offmid' in kinds:                                # LP-D18: head-on leaders that end away from the middle of the box end
        for lab, b, lead in zip(labs, boxes, leads):
            if inside(b) and lp_metrics.end_off_middle([lab], [b], [lead], [False]):
                count += 1
                outline(b, MARK_COLORS['offmid']); ring(lead[-1], MARK_COLORS['offmid'])
    if 'unclear' in kinds:
        raster = lp_metrics.ink_raster(str(Path(panel['stem']).with_suffix('.svg')))
        gaps = lp_metrics.end_gaps(labs, boxes, leads, [False] * len(labs), raster)
        for (text, _b, own, other), b, lead in zip(gaps, boxes, leads):
            if inside(b) and own is not None and other is not None and other < own + 1.0:
                count += 1
                ring(lead[-1], MARK_COLORS['unclear'])
    if 'behind' in kinds:
        color = MARK_COLORS['behind']
        for i, lead in enumerate(leads):
            segs = list(zip(lead, lead[1:]))
            hit = [j for j, b in enumerate(boxes) if j != i and any(lp_metrics._seg_box(p, q, b, -0.5) for p, q in segs)]
            if not hit or not any(inside(boxes[j]) for j in hit):
                continue
            count += 1
            draw.line([px(*p) for p in lead], fill=color, width=wide)
            for j in hit:
                outline(boxes[j], color)
    return count
