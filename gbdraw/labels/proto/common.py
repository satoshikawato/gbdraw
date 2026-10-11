"""Shared geometry for the label-placement study prototypes.

Coordinates are circle-local pixels (centre 0,0; y grows downwards).
Angles are radians clockwise from twelve o'clock, so unit(theta) = (sin, -cos).
A label box is an axis-aligned rectangle (x0, y0, x1, y1); text is horizontal.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Sequence

Point = tuple[float, float]
Box = tuple[float, float, float, float]


def unit(theta: float) -> Point:
    return math.sin(theta), -math.cos(theta)


def polar(radius: float, theta: float) -> Point:
    return radius * math.sin(theta), -radius * math.cos(theta)


def angle_of(x: float, y: float) -> float:
    """Clockwise-from-top angle (radians, 0..2pi) of a point."""
    return math.atan2(x, -y) % (2.0 * math.pi)


def angle_diff(a: float, b: float) -> float:
    """Signed smallest difference a-b in (-pi, pi]."""
    return (a - b + math.pi) % (2.0 * math.pi) - math.pi


@dataclass
class Item:
    index: int            # position in clockwise order of the feature midpoints
    src: dict             # label dict built by prepare_label_list (outer, not embedded)
    text: str
    w: float              # measured text width (px)
    h: float              # measured text height (px)
    theta: float          # feature midpoint angle
    anchor: Point         # where the leader meets the feature
    feat_outer_r: float   # outer radius of the feature's lane
    span_bp: float        # length of the labelled segment (a simple priority: longer = more important)


@dataclass
class Ctx:
    total_length: int
    radius: float         # feature anchor radius R (390 px by default)
    r_elbow: float        # radius of today's leader elbow (just outside the outermost ring)
    r_keep: float         # no label box may come inside this radius
    gap: float            # minimum gap between two label boxes (labels.spacing.circular, 3 px)
    font_size: float
    arc_x_radius: float   # today's label ellipse (the default label band)
    arc_y_radius: float
    options: dict[str, str] = field(default_factory=dict)
    dropped: list[str] = field(default_factory=list)   # texts of labels the design gave up on
    notes: dict[str, Any] = field(default_factory=dict)
    radial_layout: Any = None   # CircularRadialLayout (inner tracks, ticks, title radius); used by design G
    feature_colors: dict[str, Any] = field(default_factory=dict)   # stable_id -> feature colour; used by design G

    @classmethod
    def build(cls, items, *, total_length, radius, gap, font_size, arc_x_radius, arc_y_radius, options):
        elbows = sorted(math.hypot(float(i.src["middle_x"]), float(i.src["middle_y"])) for i in items)
        outer = max((i.feat_outer_r for i in items), default=radius)
        r_elbow = max(elbows[len(elbows) // 2] if elbows else radius, outer + 4.0)
        return cls(
            total_length=total_length, radius=radius, r_elbow=r_elbow, r_keep=r_elbow + 4.0,
            gap=float(gap), font_size=font_size, arc_x_radius=arc_x_radius, arc_y_radius=arc_y_radius,
            options=dict(options),
        )

    def opt(self, key: str, default: float) -> float:
        try:
            return float(self.options.get(key, default))
        except ValueError:
            return default



def make_items(labels: Sequence[dict], total_length: int) -> list[Item]:
    ordered = sorted(labels, key=lambda label: float(label["middle"]) % total_length)
    items = []
    for index, label in enumerate(ordered):
        theta = 2.0 * math.pi * ((float(label["middle"]) % total_length) / float(total_length))
        items.append(
            Item(
                index=index, src=label, text=str(label["label_text"]),
                w=float(label["width_px"]), h=float(label["height_px"]), theta=theta,
                anchor=(float(label["feature_anchor_x"]), float(label["feature_anchor_y"])),
                feat_outer_r=float(label.get("feature_outer_radius_px", 0.0)),
                span_bp=abs(float(label["end"]) - float(label["start"])),
            )
        )
    return items


def line_pitch(items: Sequence[Item], gap: float) -> float:
    return max((i.h for i in items), default=0.0) + gap


# ---------------------------------------------------------------- boxes and segments

def box_from_left(x0: float, yc: float, w: float, h: float) -> Box:
    return (x0, yc - 0.5 * h, x0 + w, yc + 0.5 * h)


def box_from_right(x1: float, yc: float, w: float, h: float) -> Box:
    return (x1 - w, yc - 0.5 * h, x1, yc + 0.5 * h)


def box_from_center(cx: float, cy: float, w: float, h: float) -> Box:
    return (cx - 0.5 * w, cy - 0.5 * h, cx + 0.5 * w, cy + 0.5 * h)


def box_on_ray(theta: float, r: float, w: float, h: float) -> Box:
    """Box centred on the ray theta and pushed out until it clears the circle of radius r."""
    sx, sy = unit(theta)
    ext = 0.5 * w * abs(sx) + 0.5 * h * abs(sy)
    return box_from_center((r + ext) * sx, (r + ext) * sy, w, h)


def boxes_overlap(a: Box, b: Box, gap: float = 0.0) -> bool:
    return a[0] < b[2] + gap and b[0] < a[2] + gap and a[1] < b[3] + gap and b[1] < a[3] + gap


def box_min_radius(box: Box) -> float:
    """Distance from the centre to the nearest point of the box."""
    dx = max(box[0], 0.0, -box[2])
    dy = max(box[1], 0.0, -box[3])
    return math.hypot(dx, dy)


def nearest_point_on_box(box: Box, p: Point, inset: float = 0.0) -> Point:
    """Nearest point of the box outline to p (p is assumed outside the box)."""
    x = min(max(p[0], box[0] + inset), box[2] - inset) if box[2] - box[0] > 2 * inset else 0.5 * (box[0] + box[2])
    y = min(max(p[1], box[1] + inset), box[3] - inset) if box[3] - box[1] > 2 * inset else 0.5 * (box[1] + box[3])
    return x, y


def segments_cross(p: Point, q: Point, r: Point, s: Point) -> bool:
    """Proper intersection of segments pq and rs (shared end points do not count)."""
    def orient(a: Point, b: Point, c: Point) -> float:
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    d1, d2, d3, d4 = orient(r, s, p), orient(r, s, q), orient(p, q, r), orient(p, q, s)
    return ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)) and all(abs(d) > 1e-9 for d in (d1, d2, d3, d4))


def segment_hits_box(p: Point, q: Point, box: Box, pad: float = 0.0) -> bool:
    """True when segment pq passes through the box grown by pad (Liang-Barsky)."""
    x0, y0, x1, y1 = box[0] - pad, box[1] - pad, box[2] + pad, box[3] + pad
    dx, dy = q[0] - p[0], q[1] - p[1]
    t0, t1 = 0.0, 1.0
    for d, lo, hi, start in ((dx, x0, x1, p[0]), (dy, y0, y1, p[1])):
        if abs(d) < 1e-12:
            if start <= lo or start >= hi:
                return False
            continue
        ta, tb = (lo - start) / d, (hi - start) / d
        if ta > tb:
            ta, tb = tb, ta
        t0, t1 = max(t0, ta), min(t1, tb)
        if t0 >= t1:
            return False
    return True


def segment_min_radius(p: Point, q: Point) -> float:
    """Smallest distance from the centre to segment pq (to test whether a leader cuts into the circle)."""
    dx, dy = q[0] - p[0], q[1] - p[1]
    length2 = dx * dx + dy * dy
    t = 0.0 if length2 == 0 else max(0.0, min(1.0, -(p[0] * dx + p[1] * dy) / length2))
    return math.hypot(p[0] + t * dx, p[1] + t * dy)


def polyline_segments(points: Sequence[Point]) -> list[tuple[Point, Point]]:
    return list(zip(points, points[1:]))


class Scene:
    """Boxes and leader polylines placed so far, with the collision queries the designs need."""

    def __init__(self, gap: float, leader_pad: float = 1.0) -> None:
        self.gap = gap
        self.leader_pad = leader_pad
        self.boxes: list[Box] = []
        self.leaders: list[list[Point]] = []

    def box_free(self, box: Box, gap: float | None = None) -> bool:
        gap = self.gap if gap is None else gap
        return not any(boxes_overlap(box, other, gap) for other in self.boxes)

    def box_clear_of_leaders(self, box: Box) -> bool:
        return not any(
            segment_hits_box(p, q, box, self.leader_pad)
            for leader in self.leaders for p, q in polyline_segments(leader)
        )

    def leader_free(self, points: Sequence[Point], *, check_crossings: bool = True, skip_last_box: Box | None = None) -> bool:
        """True when the polyline passes through no placed box and crosses no placed leader."""
        for p, q in polyline_segments(points):
            for other in self.boxes:
                if other is skip_last_box:
                    continue
                if segment_hits_box(p, q, other, self.leader_pad):
                    return False
            if check_crossings:
                for leader in self.leaders:
                    for r, s in polyline_segments(leader):
                        if segments_cross(p, q, r, s):
                            return False
        return True

    def add(self, box: Box, leader: Sequence[Point] | None) -> None:
        self.boxes.append(box)
        if leader:
            self.leaders.append(list(leader))


# ---------------------------------------------------------------- one-dimensional ordered packing

def pack_ordered(targets: Sequence[float], sizes: Sequence[float], lo: float | None = None, hi: float | None = None) -> list[float]:
    """Centres c_i in the given order, non-overlapping (c_{i+1}-c_i >= (s_i+s_{i+1})/2), as close to the
    targets as possible in the least-squares sense (pool-adjacent-violators), then clamped as one rigid
    stack into [lo, hi] when bounds are given. Returns the centres. If the stack is longer than the
    interval it is centred on it and overflows equally at both ends."""
    n = len(targets)
    if n == 0:
        return []
    offsets = [0.0] * n   # offset of centre i from centre 0 when packed tight
    for i in range(1, n):
        offsets[i] = offsets[i - 1] + 0.5 * (sizes[i - 1] + sizes[i])
    blocks: list[list[float]] = []   # [sum of (target-offset), count] -> start value = mean
    for i in range(n):
        blocks.append([targets[i] - offsets[i], 1.0])
        while len(blocks) > 1 and blocks[-2][0] / blocks[-2][1] > blocks[-1][0] / blocks[-1][1]:
            total, count = blocks.pop()
            blocks[-1][0] += total
            blocks[-1][1] += count
    centres: list[float] = []
    for total, count in blocks:
        start = total / count
        for _ in range(int(count)):
            centres.append(start + offsets[len(centres)])
    if lo is not None and hi is not None:
        top = centres[0] - 0.5 * sizes[0]
        bottom = centres[-1] + 0.5 * sizes[-1]
        if bottom - top > hi - lo:
            shift = 0.5 * (lo + hi) - 0.5 * (top + bottom)
            centres = [c + shift for c in centres]
        else:
            # clamp from the top, then from the bottom, keeping the order constraints
            centres[0] = max(centres[0], lo + 0.5 * sizes[0])
            for i in range(1, n):
                centres[i] = max(centres[i], centres[i - 1] + 0.5 * (sizes[i - 1] + sizes[i]))
            centres[-1] = min(centres[-1], hi - 0.5 * sizes[-1])
            for i in range(n - 2, -1, -1):
                centres[i] = min(centres[i], centres[i + 1] - 0.5 * (sizes[i] + sizes[i + 1]))
    return centres


# ---------------------------------------------------------------- output

def emit(item: Item, box: Box, leader: Sequence[Point] | None, *, align: str = "start",
         extra: Sequence[Sequence[Point]] = (), tag: str = "") -> dict:
    """Build the label dict the renderer draws: horizontal text in `box`, aligned start/middle/end,
    with the leader polyline `leader` (feature end first, label end last; None or [] for no leader).
    `extra` adds further strokes drawn like leaders (for example a bracket)."""
    label = dict(item.src)
    x0, y0, x1, y1 = box
    text_x = {"start": x0, "middle": 0.5 * (x0 + x1), "end": x1}[align]
    text_y = 0.5 * (y0 + y1)
    points = [(float(x), float(y)) for x, y in (leader or [])]
    centre_angle = math.degrees(angle_of(0.5 * (x0 + x1), text_y))
    label.update(
        placement="radial", rotation_deg=0.0, text_anchor=align, dominant_baseline="middle",
        text_x=text_x, text_y=text_y, start_x=text_x, start_y=text_y,
        aabb_local=(x0, y0, x1, y1), oriented_corners=((x0, y0), (x1, y0), (x1, y1), (x0, y1)),
        placed_angle_deg=centre_angle, leader_points=points,
        leader_extra=[[(float(x), float(y)) for x, y in line] for line in extra],
        proto_tag=tag,
    )
    if points:
        label.update(
            middle_x=points[min(1, len(points) - 1)][0], middle_y=points[min(1, len(points) - 1)][1],
            leader_start_x=points[-1][0], leader_start_y=points[-1][1],
            feature_anchor_x=points[0][0], feature_anchor_y=points[0][1],
        )
    else:
        label.update(middle_x=text_x, middle_y=text_y, leader_start_x=text_x, leader_start_y=text_y,
                     feature_anchor_x=text_x, feature_anchor_y=text_y)
    return label
