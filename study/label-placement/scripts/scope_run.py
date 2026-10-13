#!/usr/bin/env python3
"""Scope study (LP-D48 Q10, SCOPE-ALL-LABELS): render one case with today's engine and count, on the drawn layout,
label overlaps and leader crossings for every label kind, per record and across records.

Usage:
  PYTHONPATH=<worktree> /home/kawato/micromamba/bin/python3.13 -P scope_run.py --out DIR --name NAME \
      [--formats svg] -- <circular|linear> <gbdraw args without -o/-f>

Kinds: c-outer, c-inner (horizontal Circular labels outside / inside the ring), c-embedded (Circular text on the
feature arc), c-radial-outer, c-radial-inner (--label_placement radial); l-external (Linear rows above / below),
l-embedded (Linear text inside the feature), l-above (Linear --label_placement above_feature).
Definitions (same as the outer-label fit check, LP-D43 / OD-F1): two label boxes overlap when the separating-axis
penetration depth is more than 0.5 px (axis-aligned boxes: both the x and the y overlap > 0.5 px); two leaders cross
when two of their segments properly cross. A Circular embedded label is the annular sector of its text (centre
radius +- half the text height, text width along the arc), as a convex polygon. Every label is the one handed to the
label drawer (the drawn layout); its canvas offset is the sum of the translate() of its SVG ancestors, found through
its data-label-feature-id. Writes DIR/NAME.svg and DIR/NAME.scope.json.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import pathlib
import re
import sys
import time
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

os.environ.setdefault("TMPDIR", "/home/kawato/gbdraw-baselines/tmp")
ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--name", required=True)
ap.add_argument("--formats", default="svg")
ap.add_argument("rest", nargs=argparse.REMAINDER)
ns = ap.parse_args()
rest = ns.rest[1:] if ns.rest and ns.rest[0] == "--" else ns.rest
out = pathlib.Path(ns.out).resolve()
out.mkdir(parents=True, exist_ok=True)

import gbdraw.labels.circular as lc
import gbdraw.render.drawers.circular.labels as cdraw
import gbdraw.render.drawers.linear.labels as ldraw
from gbdraw.labels.linear import _label_corners

drawn: list[tuple[str, dict, int]] = []      # (mode, label, record length)
timers: Counter = Counter()

_c_draw = cdraw.LabelDrawer.draw


def _c_capture(self, label, group, record_length, *a, **k):
    drawn.append(("c", label, int(record_length)))
    return _c_draw(self, label, group, record_length, *a, **k)


cdraw.LabelDrawer.draw = _c_capture
_l_draw = ldraw.LabelDrawer.draw


def _l_capture(self, label_entry, group):
    drawn.append(("l", label_entry, 0))
    return _l_draw(self, label_entry, group)


ldraw.LabelDrawer.draw = _l_capture


def _timed(module, name, key):
    orig = getattr(module, name)

    def wrapper(*a, **k):
        t0 = time.perf_counter()
        try:
            return orig(*a, **k)
        finally:
            timers[key] += time.perf_counter() - t0
            timers[key + "_calls"] += 1
    setattr(module, name, wrapper)


import gbdraw.diagrams.circular.assemble as casm
import gbdraw.diagrams.linear.precalc as lpre
_timed(casm, "prepare_label_list", "circular_label_placement_s")
_timed(lpre, "prepare_label_list_linear", "linear_label_placement_s")
try:
    import gbdraw.render.groups.linear.seq_record as lseq
    _timed(lseq, "prepare_label_list_linear", "linear_label_placement_s")
except Exception:
    pass

from gbdraw.cli import main

os.chdir(out)
sys.argv = ["gbdraw", *rest, "-o", ns.name, "-f", ns.formats, "--overwrite"]
t0 = time.time()
rc = main()
elapsed = time.time() - t0


# ---------- geometry ----------
def kind_of(mode, lab):
    if mode == "c":
        if lab.get("is_embedded"):
            return "c-embedded"
        if lab.get("placement") == "radial":
            return "c-radial-inner" if lab.get("is_inner") else "c-radial-outer"
        return "c-inner" if lab.get("is_inner") else "c-outer"
    if not lab.get("is_embedded"):
        return "l-external"
    return "l-above" if "above_feature_place_above" in lab else "l-embedded"


def sector_polygon(lab, total_length, steps=6):
    r = float(lab.get("feature_center_radius_px") or 0.0)
    w = float(lab["width_px"]); h = float(lab["height_px"])
    if r <= 0:
        return None
    mid = math.radians(360.0 * float(lab["middle"]) / total_length - 90.0)
    half = min(math.pi * 0.95, 0.5 * w / r)
    ro, ri = r + 0.5 * h, max(0.0, r - 0.5 * h)
    outer = [(ro * math.cos(mid - half + 2 * half * k / steps), ro * math.sin(mid - half + 2 * half * k / steps)) for k in range(steps + 1)]
    inner = [(ri * math.cos(mid + half - 2 * half * k / steps), ri * math.sin(mid + half - 2 * half * k / steps)) for k in range(steps + 1)]
    return outer + inner


def geometry(mode, lab, total_length):
    """(polygon, leader polyline) in the label's own SVG coordinates."""
    if mode == "c":
        if lab.get("is_embedded"):
            return sector_polygon(lab, total_length), []
        if lab.get("placement") == "radial" and lab.get("oriented_corners"):
            poly = [tuple(map(float, p)) for p in lab["oriented_corners"]]
        else:
            x0, y0, x1, y1 = lc._label_bbox_rect(lab, total_length, 0.0)
            poly = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        lead = [
            (float(lab.get("feature_anchor_x", lab["feature_middle_x"])), float(lab.get("feature_anchor_y", lab["feature_middle_y"]))),
            (float(lab["middle_x"]), float(lab["middle_y"])),
            (float(lab.get("leader_start_x", lab["start_x"])), float(lab.get("leader_start_y", lab["start_y"]))),
        ]
        return poly, lead
    poly = [tuple(map(float, p)) for p in _label_corners(lab)]
    if lab.get("leader_line"):
        lead = [(float(lab["leader_start_x"]), float(lab["leader_start_y"])), (float(lab["leader_end_x"]), float(lab["leader_end_y"]))]
    elif not lab.get("is_embedded"):
        ly, fy, hh = float(lab["middle_y"]), float(lab["feature_middle_y"]), float(lab["height_px"])
        edge = ly - 0.45 * hh if ly >= fy else ly + 0.45 * hh
        lead = [(float(lab["middle"]), fy), (float(lab["middle"]), edge)]
    else:
        lead = []
    return poly, lead


def sat_depth(a, b):
    """Smallest projection overlap of two convex polygons (> 0: they intersect that deep; <= 0: separated)."""
    best = float("inf")
    for poly in (a, b):
        n = len(poly)
        for i in range(n):
            (x0, y0), (x1, y1) = poly[i], poly[(i + 1) % n]
            ex, ey = x1 - x0, y1 - y0
            ln = math.hypot(ex, ey)
            if ln < 1e-9:
                continue
            ax, ay = -ey / ln, ex / ln
            pa = [x * ax + y * ay for x, y in a]
            pb = [x * ax + y * ay for x, y in b]
            ov = min(max(pa), max(pb)) - max(min(pa), min(pb))
            if ov < best:
                best = ov
                if best <= 0:
                    return best
    return best


def _cross(p, q, r, s):
    def o(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    d1, d2, d3, d4 = o(r, s, p), o(r, s, q), o(p, q, r), o(p, q, s)
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)) and all(abs(d) > 1e-6 for d in (d1, d2, d3, d4))


def _cross_point(p, q, r, s):
    den = (q[0] - p[0]) * (s[1] - r[1]) - (q[1] - p[1]) * (s[0] - r[0])
    t = ((r[0] - p[0]) * (s[1] - r[1]) - (r[1] - p[1]) * (s[0] - r[0])) / den
    return (p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1]))


def _seg_poly(p, q, poly, shrink=0.5):
    """Segment p-q passes through convex polygon `poly` more than `shrink` px deep: along every polygon edge normal
    the segment's projection overlaps the polygon's by more than `shrink`, and along the segment's own normal the
    polygon reaches more than `shrink` beyond the segment line on both sides."""
    n = len(poly)
    for i in range(n):
        (x0, y0), (x1, y1) = poly[i], poly[(i + 1) % n]
        ax, ay = -(y1 - y0), x1 - x0
        ln = math.hypot(ax, ay)
        if ln < 1e-9:
            continue
        ax, ay = ax / ln, ay / ln
        pa = [x * ax + y * ay for x, y in poly]
        ps = [p[0] * ax + p[1] * ay, q[0] * ax + q[1] * ay]
        if min(max(pa), max(ps)) - max(min(pa), min(ps)) <= shrink:
            return False
    ax, ay = -(q[1] - p[1]), q[0] - p[0]
    ln = math.hypot(ax, ay)
    if ln < 1e-9:
        return False
    ax, ay = ax / ln, ay / ln
    s0 = p[0] * ax + p[1] * ay
    pa = [x * ax + y * ay for x, y in poly]
    return min(pa) < s0 - shrink and max(pa) > s0 + shrink


# ---------- SVG offsets ----------
TR = re.compile(r"(\w+)\s*\(([^)]*)\)")


def parse_translate(transform):
    dx = dy = 0.0
    other = []
    for fn, argstr in TR.findall(transform or ""):
        nums = [float(v) for v in re.split(r"[\s,]+", argstr.strip()) if v]
        if fn == "translate":
            dx += nums[0]; dy += nums[1] if len(nums) > 1 else 0.0
        else:
            other.append(fn)
    return dx, dy, other


svg_path = out / f"{ns.name}.svg"
offsets: dict[str, list[tuple[float, float]]] = defaultdict(list)
non_translate = Counter()
canvas = None
if svg_path.exists():
    root = ET.parse(svg_path).getroot()
    canvas = [round(float(re.match(r"[0-9.]+", root.get("width")).group(0))), round(float(re.match(r"[0-9.]+", root.get("height")).group(0)))]

    def walk(el, dx, dy):
        for child in el:
            if child.get("data-gbdraw-label-binding-schema") == "1" and child.get("data-label-feature-id"):
                offsets[child.get("data-label-feature-id")].append((dx, dy))
                continue
            tx, ty, other = parse_translate(child.get("transform"))
            if other and len(child):
                non_translate.update(other)
            walk(child, dx + tx, dy + ty)
    walk(root, 0.0, 0.0)

# ---------- collect labels ----------
labels = []
seen = set()
used = Counter()
record_of_length: dict = {}
record_counter = 0
for mode, lab, tl in drawn:
    if id(lab) in seen:
        continue
    seen.add(id(lab))
    fid = str(lab.get("feature_id"))
    k = used[fid]; used[fid] += 1
    off = offsets.get(fid, [])
    dx, dy = off[k] if k < len(off) else (None, None)
    poly, lead = geometry(mode, lab, tl)
    if poly is None or dx is None:
        labels.append(None)
        continue
    poly = [(x + dx, y + dy) for x, y in poly]
    lead = [(x + dx, y + dy) for x, y in lead]
    rec = lab.get("record_index")
    labels.append({"kind": kind_of(mode, lab), "text": str(lab.get("label_text", "")), "poly": poly, "lead": lead,
                   "rec": (round(dx, 3), round(dy, 3)) if mode == "c" else None, "fid": fid})
missing = sum(1 for v in labels if v is None)
labels = [v for v in labels if v is not None]

# record identity: Circular = the record's translate (one per record); Linear = the record group of the label
if labels and labels[0]["rec"] is None and svg_path.exists():
    # Linear: map each label to its record by the record group that contains its text element.
    rec_of = defaultdict(list)

    def walk2(el, rec):
        for child in el:
            r = child.get("data-gbdraw-record-index", rec)
            if child.get("data-gbdraw-label-binding-schema") == "1" and child.get("data-label-feature-id"):
                rec_of[child.get("data-label-feature-id")].append(r)
                continue
            walk2(child, r)
    walk2(root, None)
    used2 = Counter()
    for v in labels:
        k = used2[v["fid"]]; used2[v["fid"]] += 1
        rr = rec_of.get(v["fid"], [])
        v["rec"] = rr[k] if k < len(rr) else None
recs = sorted({str(v["rec"]) for v in labels})
rec_index = {r: i for i, r in enumerate(recs)}
for v in labels:
    v["rec"] = rec_index[str(v["rec"])]

# ---------- pairs ----------
def aabb(pts):
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    return min(xs), min(ys), max(xs), max(ys)


n = len(labels)
boxes = [aabb(v["poly"]) for v in labels]
order = sorted(range(n), key=lambda i: boxes[i][0])
overlap = Counter(); overlap_list = []; depth_max = Counter()
active: list[int] = []
for i in order:
    x0 = boxes[i][0]
    active = [j for j in active if boxes[j][2] > x0 + 0.5]
    for j in active:
        if boxes[j][1] >= boxes[i][3] - 0.5 or boxes[i][1] >= boxes[j][3] - 0.5:
            continue
        d = sat_depth(labels[i]["poly"], labels[j]["poly"])
        if d > 0.5:
            a, b = sorted((labels[i]["kind"], labels[j]["kind"]))
            key = f"{a}|{b}|{'same' if labels[i]['rec'] == labels[j]['rec'] else 'cross'}"
            overlap[key] += 1
            depth_max[key] = max(depth_max[key], round(d, 1))
            if len(overlap_list) < 4000:
                overlap_list.append([key, i, j, round(d, 1)])
    active.append(i)

segs = [list(zip(v["lead"], v["lead"][1:])) for v in labels]
lbox = [aabb(v["lead"]) if v["lead"] else None for v in labels]
lorder = sorted([i for i in range(n) if lbox[i]], key=lambda i: lbox[i][0])
crossing = Counter(); crossing_list = []
active = []
for i in lorder:
    x0 = lbox[i][0]
    active = [j for j in active if lbox[j][2] >= x0]
    for j in active:
        if lbox[j][1] > lbox[i][3] or lbox[i][1] > lbox[j][3]:
            continue
        pts = [_cross_point(p, q, r, s) for p, q in segs[i] for r, s in segs[j] if _cross(p, q, r, s)]
        if pts:
            a, b = sorted((labels[i]["kind"], labels[j]["kind"]))
            key = f"{a}|{b}|{'same' if labels[i]['rec'] == labels[j]['rec'] else 'cross'}"
            crossing[key] += 1
            if len(crossing_list) < 4000:
                crossing_list.append([key, i, j, [round(pts[0][0], 2), round(pts[0][1], 2)]])
    active.append(i)

through = Counter(); through_list = []
for i in range(n):
    for p, q in segs[i]:
        sb = aabb([p, q])
        for j in range(n):
            if j == i:
                continue
            b = boxes[j]
            if sb[2] < b[0] or sb[0] > b[2] or sb[3] < b[1] or sb[1] > b[3]:
                continue
            if _seg_poly(p, q, labels[j]["poly"]):
                key = f"{labels[i]['kind']}>{labels[j]['kind']}"
                through[key] += 1
                if len(through_list) < 2000:
                    through_list.append([key, i, j])
                break

kinds = Counter(v["kind"] for v in labels)
per_rec = Counter((v["rec"], v["kind"]) for v in labels)
result = {
    "argv": rest, "rc": rc, "seconds": round(elapsed, 2), "canvas_px": canvas,
    "timers": {k: round(v, 3) if isinstance(v, float) else v for k, v in timers.items()},
    "labels_by_kind": dict(kinds), "records": len(recs),
    "labels_by_record_kind": {f"{r}:{k}": c for (r, k), c in sorted(per_rec.items())},
    "unmatched_labels": missing, "non_translate_ancestors": dict(non_translate),
    "overlap_pairs": dict(overlap), "overlap_depth_max": dict(depth_max),
    "crossing_pairs": dict(crossing), "through_text": dict(through),
    "overlap_total": sum(overlap.values()), "crossing_total": sum(crossing.values()),
}
(out / f"{ns.name}.scope.json").write_text(json.dumps(result, indent=1))
(out / f"{ns.name}.scope-geom.json").write_text(json.dumps({
    "labels": [{"kind": v["kind"], "text": v["text"], "rec": v["rec"], "poly": [[round(x, 2), round(y, 2)] for x, y in v["poly"]],
                "lead": [[round(x, 2), round(y, 2)] for x, y in v["lead"]]} for v in labels],
    "overlaps": overlap_list, "crossings": crossing_list, "through": through_list}))
print(json.dumps({k: result[k] for k in ("labels_by_kind", "records", "overlap_pairs", "crossing_pairs", "through_text",
                                         "unmatched_labels", "seconds", "canvas_px", "timers")}))
sys.exit(rc or 0)
