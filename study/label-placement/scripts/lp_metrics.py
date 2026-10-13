"""Metrics for one rendered circular label layout (study script; works for today's labels and prototypes)."""
from __future__ import annotations

import math
import re


def label_box(label, total_length, lc):
    if label.get("aabb_local") is not None:
        return tuple(float(v) for v in label["aabb_local"])
    return tuple(float(v) for v in lc._label_bbox_rect(label, total_length, 0.0))


def label_leader(label):
    if label.get("leader_points") is not None:
        return [tuple(p) for p in label["leader_points"]]
    return [
        (float(label.get("feature_anchor_x", label["feature_middle_x"])), float(label.get("feature_anchor_y", label["feature_middle_y"]))),
        (float(label["middle_x"]), float(label["middle_y"])),
        (float(label.get("leader_start_x", label["start_x"])), float(label.get("leader_start_y", label["start_y"]))),
    ]


def _cross(p, q, r, s):
    def o(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    d1, d2, d3, d4 = o(r, s, p), o(r, s, q), o(p, q, r), o(p, q, s)
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)) and all(abs(d) > 1e-6 for d in (d1, d2, d3, d4))


def _seg_box(p, q, box, pad=0.0):
    x0, y0, x1, y1 = box[0] - pad, box[1] - pad, box[2] + pad, box[3] + pad
    dx, dy = q[0] - p[0], q[1] - p[1]
    t0, t1 = 0.0, 1.0
    for d, lo, hi, st in ((dx, x0, x1, p[0]), (dy, y0, y1, p[1])):
        if abs(d) < 1e-12:
            if st <= lo or st >= hi:
                return False
            continue
        ta, tb = (lo - st) / d, (hi - st) / d
        if ta > tb:
            ta, tb = tb, ta
        t0, t1 = max(t0, ta), min(t1, tb)
        if t0 >= t1:
            return False
    return True


def _seg_min_radius(p, q):
    dx, dy = q[0] - p[0], q[1] - p[1]
    l2 = dx * dx + dy * dy
    t = 0.0 if l2 == 0 else max(0.0, min(1.0, -(p[0] * dx + p[1] * dy) / l2))
    return math.hypot(p[0] + t * dx, p[1] + t * dy)


def min_box_gap(boxes, skip=()):
    """Smallest distance between two label boxes (px); negative = they overlap by that depth. Rotated boxes skipped."""
    best = None
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            if i in skip or j in skip:
                continue
            a, b = boxes[i], boxes[j]
            dx = max(a[0] - b[2], b[0] - a[2]); dy = max(a[1] - b[3], b[1] - a[3])
            gap = math.hypot(dx, dy) if dx > 0 and dy > 0 else max(dx, dy)
            best = gap if best is None else min(best, gap)
    return None if best is None else round(best, 1)


def toward_axis(ext, boxes, leaders, rotated, tol=1.0):
    """Labels that run back toward the ring (LP-D09): the end of the text box away from the leader is closer to
    the ring than the end the leader meets. Returns the label texts. Rotated text is skipped."""
    out = []
    for lab, box, lead, rot in zip(ext, boxes, leaders, rotated):
        if rot or not lead:
            continue
        px, py = lead[-1]
        ring = float(lab.get("feature_center_radius_px", 0.0)) or math.hypot(*lead[0])
        dist = lambda x, y: abs(math.hypot(x, y) - ring)
        corners = [(x, y) for x in (box[0], box[2]) for y in (box[1], box[3])]
        near = min(corners, key=lambda c: (round(math.dist(c, (px, py)), 1), dist(*c)))
        far_x = box[2] if near[0] == box[0] else box[0]
        if min(dist(far_x, box[1]), dist(far_x, box[3])) < dist(*near) - tol:
            out.append(str(lab.get("label_text", "")))
    return out


def back_px(ext, boxes, leaders, rotated):
    """LP-D09 after LP-D13 (a leader may end anywhere on the edge of its label box): per label, how many px the
    point of its box nearest to the ring is closer to the ring than the end of its leader, and whether that is more
    than the box height (the corner rule of H2 allows up to one box height: the leader may meet either corner of
    the text end). Returns [(px, over the box height, text)], largest first. Rotated text is skipped."""
    out = []
    for lab, box, lead, rot in zip(ext, boxes, leaders, rotated):
        if rot or not lead:
            continue
        ring = float(lab.get("feature_center_radius_px", 0.0)) or math.hypot(*lead[0])
        end_r = math.hypot(*lead[-1])
        if end_r >= ring:                                  # outside label: the box point nearest to the centre
            back = end_r - math.hypot(min(max(0.0, box[0]), box[2]), min(max(0.0, box[1]), box[3]))
        else:                                              # inside label: the box corner farthest from the centre
            back = max(math.hypot(x, y) for x in (box[0], box[2]) for y in (box[1], box[3])) - end_r
        back = max(0.0, back)
        out.append((round(back, 1), back > (box[3] - box[1]) + 0.5, str(lab.get("label_text", ""))))
    return sorted(out, reverse=True)


def _arrival(pts, min_seg=3.0):
    """Unit direction in which a leader reaches its label: its last piece of at least `min_seg` px (shorter pieces
    are dropped as in leader_turns), or the whole leader when it has no such piece."""
    kept = [p for k, p in enumerate(pts) if k == 0 or math.dist(p, pts[k - 1]) >= min_seg]
    a, b = (kept[-2], kept[-1]) if len(kept) >= 2 else (pts[0], pts[-1])
    n = math.dist(a, b)
    return ((b[0] - a[0]) / n, (b[1] - a[1]) / n) if n > 1e-9 else (0.0, 0.0)


def end_far_half(ext, boxes, leaders, rotated, min_deg=10.0, tol=1.0, eps=0.05):
    """LP-D15: leaders that end on the far half of their label box. A leader that comes from above must end between
    the top vertex and the middle of the box side it meets, one from below between the bottom vertex and the middle,
    and the same along a top or bottom side for a leader that comes from the left or the right. An end off a corner is
    tested along both sides of that corner. No side is "far" when the leader reaches the side within `min_deg` degrees
    of the perpendicular, or when the end is within `tol` px of the middle. Rotated text is skipped.
    Returns (far, places): far = [(text, 'x'|'y', position along the side from -1 to 1 with 0 = middle and positive =
    far side, degrees off the perpendicular)], worst first; places = how many ends are off a corner, on a side, inside."""
    far, places = [], {"corner": 0, "side": 0, "inside": 0}
    for lab, box, lead, rot in zip(ext, boxes, leaders, rotated):
        if rot or len(lead) < 2:
            continue
        (ex, ey), d = lead[-1], _arrival(lead)
        out_x = ex <= box[0] + eps or ex >= box[2] - eps
        out_y = ey <= box[1] + eps or ey >= box[3] - eps
        places["corner" if out_x and out_y else "side" if out_x or out_y else "inside"] += 1
        worst = None
        for axis, on_side in ((1, out_x), (0, out_y)):      # along a left or right side: y; along a top or bottom side: x
            if not on_side:
                continue
            half = (box[axis + 2] - box[axis]) / 2.0
            off = (ex, ey)[axis] - (box[axis] + box[axis + 2]) / 2.0
            deg = math.degrees(math.asin(max(-1.0, min(1.0, abs(d[axis])))))
            if deg >= min_deg and off * d[axis] > 0 and abs(off) > tol:
                item = (str(lab.get("label_text", "")), "xy"[axis], round(min(1.5, abs(off) / half), 2), round(deg, 1))
                worst = item if worst is None or item[2] > worst[2] else worst
        if worst:
            far.append(worst)
    return sorted(far, key=lambda it: (-it[2], -it[3])), places


def end_off_middle(ext, boxes, leaders, rotated, max_deg=10.0, frac=0.35, eps=0.05):
    """LP-D18: leaders that reach the left or right end of their label box head-on (within `max_deg` degrees of the
    horizontal) but end away from the middle of that end: more than `frac` of the half height above or below it
    (an end off a corner counts as 1). Such a leader should come out of the middle of the box end. Top and bottom
    sides are not tested (the Owner spoke of the box end; open for the poles). Rotated text is skipped.
    Returns [(text, position from -1 = top to 1 = bottom, degrees off the horizontal)], worst first."""
    out = []
    for lab, box, lead, rot in zip(ext, boxes, leaders, rotated):
        if rot or len(lead) < 2:
            continue
        (ex, ey), d = lead[-1], _arrival(lead)
        if not (ex <= box[0] + eps or ex >= box[2] - eps):
            continue
        deg = math.degrees(math.asin(max(-1.0, min(1.0, abs(d[1])))))
        pos = (ey - (box[1] + box[3]) / 2.0) / ((box[3] - box[1]) / 2.0)
        if deg < max_deg and abs(d[0]) > 0.5 and abs(pos) > frac:
            out.append((str(lab.get("label_text", "")), round(max(-1.5, min(1.5, pos)), 2), round(deg, 1)))
    return sorted(out, key=lambda it: -abs(it[1]))


def ink_raster(svg_path, scale=3, pad=24):
    """The label text of a render alone, drawn `scale` times larger by cairosvg (the renderer of the Owner pages):
    (boolean array of ink pixels, function from circle-local px to array coordinates). None when the SVG has no
    label_text group."""
    import io
    import xml.etree.ElementTree as ET
    import cairosvg
    import numpy as np
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    root = ET.parse(svg_path).getroot()
    group = next((g for g in root if g.get("id") == "label_text"), None)
    if group is None:
        return None
    for child in list(root):
        if child is not group and not child.tag.endswith("defs"):
            root.remove(child)
    w = float(re.match(r"[0-9.]+", root.get("width")).group(0)); h = float(re.match(r"[0-9.]+", root.get("height")).group(0))
    root.set("width", f"{w + 2 * pad}px"); root.set("height", f"{h + 2 * pad}px")
    root.set("viewBox", f"{-pad} {-pad} {w + 2 * pad} {h + 2 * pad}")
    cx = cy = 0.0
    for x, y in re.findall(r"translate\(\s*([-0-9.eE]+)[ ,]+([-0-9.eE]+)\s*\)", group.get("transform", "")):
        cx += float(x); cy += float(y)
    png = cairosvg.svg2png(bytestring=ET.tostring(root), background_color="white", scale=scale)
    ink = np.asarray(Image.open(io.BytesIO(png)).convert("L")) < 128
    return ink, lambda x, y: ((x + cx + pad) * scale, (y + cy + pad) * scale), scale


def end_gaps(ext, boxes, leaders, rotated, raster=None, grow=1.5, reach=14.0):
    """LP-D15: how far the end of each leader stands from its label. Per label: distance in px to the text box, and,
    with `raster` (ink_raster), to the nearest ink of its own text and of any other label's text within `reach` px
    (None when there is none). A label's ink is the ink inside its box grown by `grow` px. Rotated text is skipped.
    Returns [(text, to box, to own ink or None, to other ink or None)]."""
    import numpy as np
    out = []

    def ink_distance(end, box):
        ink, to_px, scale = raster
        (x0, y0), (x1, y1) = to_px(box[0] - grow, box[1] - grow), to_px(box[2] + grow, box[3] + grow)
        c0, r0 = max(0, int(math.floor(x0))), max(0, int(math.floor(y0)))
        ys, xs = np.nonzero(ink[r0:int(math.ceil(y1)) + 1, c0:int(math.ceil(x1)) + 1])
        if not len(xs):
            return None
        px, py = to_px(*end)
        return float(np.sqrt(((xs + c0 + 0.5 - px) ** 2 + (ys + r0 + 0.5 - py) ** 2).min())) / scale

    def box_distance(end, box):
        return math.hypot(max(box[0] - end[0], 0.0, end[0] - box[2]), max(box[1] - end[1], 0.0, end[1] - box[3]))

    for i, (lab, box, lead, rot) in enumerate(zip(ext, boxes, leaders, rotated)):
        if rot or len(lead) < 2:
            continue
        end = lead[-1]
        own = other = None
        if raster is not None:
            own = ink_distance(end, box)
            near = [ink_distance(end, b) for j, b in enumerate(boxes) if j != i and not rotated[j] and box_distance(end, b) <= reach]
            near = [v for v in near if v is not None]
            other = min(near) if near else None
        out.append((str(lab.get("label_text", "")), box_distance(end, box), own, other))
    return out


def end_fields(ext, boxes, leaders, rotated, svg_path=None):
    """Metrics fields of LP-D15 (where a leader ends, and how far from its label)."""
    far, places = end_far_half(ext, boxes, leaders, rotated)
    off_mid = end_off_middle(ext, boxes, leaders, rotated)
    raster = None
    if svg_path:
        try:
            raster = ink_raster(svg_path)
        except Exception as error:                      # no cairosvg, or an SVG it cannot read: the box fields remain
            raster = None
            print("lp_metrics: no ink raster:", type(error).__name__, error)
    gaps = end_gaps(ext, boxes, leaders, rotated, raster)
    to_box = [g[1] for g in gaps]
    inked = [g for g in gaps if g[2] is not None]
    to_ink = sorted(g[2] for g in inked)
    unclear = [g for g in inked if g[3] is not None and g[3] < g[2] + 1.0]
    r1 = lambda v: round(v, 1)
    return {
        "leaders_end_far_half": len(far),
        "leaders_end_far_half_list": [f"{t}: {axis} {pos} {deg}" for t, axis, pos, deg in far[:40]],
        "leader_end_places": places,
        "leaders_end_off_middle": len(off_mid),
        "leaders_end_off_middle_list": [f"{t}: {pos} {deg}" for t, pos, deg in off_mid[:40]],
        "leader_end_gap_box_min": r1(min(to_box)) if to_box else None,
        "leader_end_gap_box_max": r1(max(to_box)) if to_box else None,
        "leader_end_gap_ink_min": r1(to_ink[0]) if to_ink else None,
        "leader_end_gap_ink_median": r1(to_ink[len(to_ink) // 2]) if to_ink else None,
        "leader_end_gap_ink_max": r1(to_ink[-1]) if to_ink else None,
        "leader_end_gap_ink_over5": sum(v > 5.0 for v in to_ink) if to_ink else None,
        "leader_end_gap_ink_list": [f"{t}: {r1(own)} (box {r1(b)})" for t, b, own, _o in sorted(inked, key=lambda g: -g[2])[:15]],
        "leader_ends_unclear": len(unclear) if raster is not None else None,
        "leader_ends_unclear_list": [f"{t}: own {r1(own)} other {r1(o)}" for t, _b, own, o in unclear[:40]],
    }


def acute_bends(leaders, min_seg=3.0, limit=89.0):
    """Leaders with a bend sharper than a right angle (LP-D10). The first piece (anchor to the ring edge) counts,
    so a leader that leaves the ring at a grazing angle is found too. Returns (indexes, smallest angle in degrees)."""
    found, smallest = [], 180.0
    for i, pts in enumerate(leaders):
        pts = [p for k, p in enumerate(pts) if k == 0 or math.dist(p, pts[k - 1]) >= min_seg]
        worst = 180.0
        for a, b, c in zip(pts, pts[1:], pts[2:]):
            u = (a[0] - b[0], a[1] - b[1]); v = (c[0] - b[0], c[1] - b[1])
            nu, nv = math.hypot(*u), math.hypot(*v)
            if nu < min_seg or nv < min_seg:
                continue
            worst = min(worst, math.degrees(math.acos(max(-1.0, min(1.0, (u[0] * v[0] + u[1] * v[1]) / (nu * nv))))))
        smallest = min(smallest, worst)
        if worst < limit:
            found.append(i)
    return found, round(smallest, 1)


def leader_turns(leaders, min_seg=3.0):
    """Per leader: how far the line turns in total, in degrees (0 = one straight line; LP-D11). Pieces shorter than
    `min_seg` are dropped as in acute_bends, so a label seated on its ring-edge point counts as straight."""
    out = []
    for pts in leaders:
        pts = [p for k, p in enumerate(pts) if k == 0 or math.dist(p, pts[k - 1]) >= min_seg]
        turn = 0.0
        for a, b, c in zip(pts, pts[1:], pts[2:]):
            u = (b[0] - a[0], b[1] - a[1]); v = (c[0] - b[0], c[1] - b[1])
            nu, nv = math.hypot(*u), math.hypot(*v)
            if nu < min_seg or nv < min_seg:
                continue
            turn += math.degrees(math.acos(max(-1.0, min(1.0, (u[0] * v[0] + u[1] * v[1]) / (nu * nv)))))
        out.append(turn)
    return out


def turn_fields(ext, leaders):
    drawn = [(i, t) for i, t in enumerate(leader_turns(leaders)) if len(leaders[i]) >= 2]
    turns = [t for _, t in drawn]
    worst = sorted(drawn, key=lambda it: -it[1])[:15]
    return {
        "leader_turn_deg_sum": round(sum(turns), 1),
        "leader_turn_deg_mean": round(sum(turns) / len(turns), 2) if turns else 0.0,
        "leader_turn_deg_max": round(max(turns), 1) if turns else 0.0,
        "leaders_straight": sum(t < 3.0 for t in turns),            # reads as one straight line
        "leaders_turn_ge30": sum(t >= 30.0 for t in turns),
        "leaders_drawn": len(turns),
        "leader_turn_worst": [[str(ext[i].get("label_text", "")), round(t, 1)] for i, t in worst],
    }


def compute(labels, total_length, lc, svg_path=None):
    ext = [lab for lab in labels if not lab.get("is_embedded")]
    boxes = [label_box(lab, total_length, lc) for lab in ext]
    leaders = [label_leader(lab) for lab in ext]
    n = len(ext)
    from gbdraw.layout.text_geometry import convex_polygons_intersect
    rotated = [abs(float(lab.get("rotation_deg", 0.0))) > 1e-6 and lab.get("oriented_corners") is not None for lab in ext]
    overlap_pairs = overprint_pairs = 0
    overprinted = set()
    for i in range(n):
        a = boxes[i]
        for j in range(i + 1, n):
            b = boxes[j]
            ix = min(a[2], b[2]) - max(a[0], b[0])
            iy = min(a[3], b[3]) - max(a[1], b[1])
            if ix > 0.5 and iy > 0.5 and (rotated[i] or rotated[j]):
                # rotated text (today's radial placement): exact oriented-box test, no area estimate
                if convex_polygons_intersect(ext[i]["oriented_corners"], ext[j]["oriented_corners"], touching=False):
                    overlap_pairs += 1
                continue
            if ix > 0.5 and iy > 0.5:
                overlap_pairs += 1
                small = min((a[2] - a[0]) * (a[3] - a[1]), (b[2] - b[0]) * (b[3] - b[1]))
                if small > 0 and ix * iy >= 0.2 * small:
                    overprint_pairs += 1
                    overprinted.update((i, j))
    crossings = 0
    through_text = 0
    segs = [list(zip(l, l[1:])) for l in leaders]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if j > i and any(_cross(p, q, r, s) for p, q in segs[i] for r, s in segs[j]):
                crossings += 1
            if not rotated[j] and any(_seg_box(p, q, boxes[j], -0.5) for p, q in segs[i]):
                through_text += 1
    lengths = [sum(math.dist(p, q) for p, q in s) for s in segs if s]
    feat_r = max((float(lab.get("feature_outer_radius_px", 0.0)) for lab in ext), default=0.0)
    over_features = sum(
        1 for s in segs if any(_seg_min_radius(p, q) < feat_r - 1.0 for p, q in s[1:])
    )
    far = []
    for lab, box in zip(ext, boxes):
        ax = float(lab.get("feature_middle_x", 0.0)); ay = float(lab.get("feature_middle_y", 0.0))
        nx = min(max(ax, box[0]), box[2]); ny = min(max(ay, box[1]), box[3])
        far.append(math.hypot(nx - ax, ny - ay))
    turned = toward_axis(ext, boxes, leaders, rotated)
    back = back_px(ext, boxes, leaders, rotated)
    sharp, min_angle = acute_bends(leaders)
    tight, _ = acute_bends(leaders, limit=110.0)
    result = {
        "labels_external": n,
        "labels_embedded": len(labels) - n,
        "font_size": sorted({float(lab.get("font_size", 0)) for lab in ext}),
        "text_overlap_pairs": overlap_pairs,
        "text_overprint_pairs": overprint_pairs,
        "labels_overprinted": len(overprinted),
        "leader_crossing_pairs": crossings,
        "leader_through_text": through_text,
        "leaders_over_features": over_features,
        "leader_len_mean": round(sum(lengths) / len(lengths), 1) if lengths else 0.0,
        "leader_len_max": round(max(lengths), 1) if lengths else 0.0,
        "leader_len_sum": round(sum(lengths), 1),
        "label_gap_min_px": min_box_gap(boxes, {i for i, r in enumerate(rotated) if r}),
        "labels_toward_axis": len(turned),
        "label_back_px_max": back[0][0] if back else 0.0,
        "labels_back_over_height": sum(1 for _b, over, _t in back if over),
        "label_back_list": [f"{t}: {b}" for b, _o, t in back[:10] if b > 0],
        **turn_fields(ext, leaders),
        **end_fields(ext, boxes, leaders, rotated, svg_path),
        "leader_acute_bends": len(sharp),
        "leader_acute_list": [str(ext[i].get("label_text", "")) for i in sharp][:40],
        "leader_min_angle_deg": min_angle,
        "leader_bends_lt110": len(tight),
        "leader_bends_lt110_list": [str(ext[i].get("label_text", "")) for i in tight][:40],
        "labels_toward_axis_list": turned[:40],
        "label_to_feature_px_mean": round(sum(far) / len(far), 1) if far else 0.0,
        "label_to_feature_px_max": round(max(far), 1) if far else 0.0,
        "labels_bbox": [round(min(b[0] for b in boxes), 1), round(min(b[1] for b in boxes), 1),
                        round(max(b[2] for b in boxes), 1), round(max(b[3] for b in boxes), 1)] if boxes else None,
    }
    if svg_path:
        head = open(svg_path, encoding="utf-8").read(200000)
        m = re.search(r'<svg[^>]*?\bheight="([0-9.]+)[a-z]*"', head)
        w = re.search(r'<svg[^>]*?\bwidth="([0-9.]+)[a-z]*"', head)
        if m and w:
            result["canvas_px"] = [round(float(w.group(1))), round(float(m.group(1)))]
    return result
