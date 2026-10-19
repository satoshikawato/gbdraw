"""What is drawn behind a label's letters: a halo, or a white rounded box.

The box geometry is one owner for region and feature labels: a white box with
rounded corners, ``BACKING_PADDING_EM`` font sizes beyond the letters on every
side and as the corner radius, which may stop at an edge on one side (a region
label's mark, a feature label's own leader end or underline) so it never
covers what it belongs to. A region label is drawn on it at
``BACKING_OPACITY`` (LP-D82 ``RGN-FRONT``); a feature label with
``protection = box`` at ``label_style.BOX_OPACITY``, bound to its label and
taking no pointer events.

A halo (``protection = halo``) is a ``<use>`` of the label's text, stroked
white ``label_style.HALO_WIDTH_EM`` font sizes wide with round joins and no
fill, drawn just before the letters. It is a reference, not a second
``<text>``, so it follows any edit of its text, and ``fill="none"`` keeps a
CairoSVG PDF from carrying its letters twice. CairoSVG has no
``paint-order``, so the halo cannot be a stroke on the text itself.

A Linear underline has a casing (LP-D103): a white line under it, as wide as
the line plus ``label_style.HALO_WIDTH_EM`` font sizes on each side, with butt
caps and each end pushed out by that margin and one line width. Drawn after
every leader and before every underline and letter, it cuts the leaders of
other labels where they cross the underline. The label's own leader ends where
it enters its casing, and the rest of it (its tail) is drawn after the
casings, so the casing never cuts its joint and no ink is drawn twice.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from gbdraw.labels.label_ink import FeatureLabelInk
from gbdraw.labels.label_style import BOX_OPACITY, HALO_COLOR, HALO_WIDTH_EM, UNDERLINE_GAP_EM, leader_width_px
from gbdraw.layout.spatial import Aabb
from gbdraw.layout.text_geometry import Point
from gbdraw.render.label_binding import bind_label_part
from gbdraw.svg.elements import Line, Rect, Text, Use

BACKING_PADDING_EM = 0.18
BACKING_FILL = "#ffffff"
BACKING_OPACITY = 0.7


def label_halo(text: Text, text_id: str) -> Use:
    """The halo of ``text``, which gets ``text_id`` as its ``id``; the caller draws it before the text."""

    text.attribs["id"] = text_id
    halo = Use(
        href=f"#{text_id}",
        fill="none",
        stroke=HALO_COLOR,
        stroke_width=HALO_WIDTH_EM * float(text.attribs["font-size"]),
        stroke_linejoin="round",
        stroke_linecap="round",
    )
    halo.attribs["aria-hidden"] = "true"
    halo.attribs["pointer-events"] = "none"
    return halo


@dataclass(frozen=True)
class UnderlineCasing:
    """The casing of one underline: a butt-capped white line from ``start`` to ``end``, ``width`` wide."""

    start: Point
    end: Point
    width: float

    def entry(self, outside: Point, inside: Point) -> Point:
        """Where the segment from ``outside`` to ``inside`` (a point within the casing) enters it; ``outside``
        when the whole segment lies within it."""

        length = math.dist(self.start, self.end)
        reach = math.dist(inside, outside)
        if length == 0.0 or reach == 0.0:
            return outside
        ux, uy = (self.end[0] - self.start[0]) / length, (self.end[1] - self.start[1]) / length
        dx, dy = (outside[0] - inside[0]) / reach, (outside[1] - inside[1]) / reach
        along = (inside[0] - self.start[0]) * ux + (inside[1] - self.start[1]) * uy
        across = (inside[1] - self.start[1]) * ux - (inside[0] - self.start[0]) * uy
        exit_at = math.inf
        for origin, step, low, high in (
            (along, dx * ux + dy * uy, 0.0, length),
            (across, dy * ux - dx * uy, -0.5 * self.width, 0.5 * self.width),
        ):
            if abs(step) > 1e-12:
                exit_at = min(exit_at, max((low - origin) / step, (high - origin) / step))
        if exit_at >= reach:
            return outside
        return (inside[0] + dx * exit_at, inside[1] + dy * exit_at)

    def element(self, feature_id: str) -> Line:
        """The drawn casing, bound to its label so that it is hidden with it."""

        line = Line(start=self.start, end=self.end, stroke=HALO_COLOR, stroke_width=self.width)
        line.attribs["stroke-linecap"] = "butt"
        line.attribs["aria-hidden"] = "true"
        line.attribs["pointer-events"] = "none"
        return bind_label_part(line, feature_id)


def underline_casing(start: Point, end: Point, font_size: float) -> UnderlineCasing:
    """The casing of the underline from ``start`` to ``end`` of a label of ``font_size`` px."""

    line_width = leader_width_px(font_size)
    margin = HALO_WIDTH_EM * float(font_size)
    reach = margin + line_width
    length = math.dist(start, end)
    ux, uy = ((end[0] - start[0]) / length, (end[1] - start[1]) / length) if length else (0.0, 0.0)
    return UnderlineCasing(
        start=(start[0] - ux * reach, start[1] - uy * reach),
        end=(end[0] + ux * reach, end[1] + uy * reach),
        width=line_width + 2.0 * margin,
    )


def backing_box(ink: Aabb, font_size: float, *, edge_x: float | None = None, edge_y: float | None = None) -> Aabb:
    """The box a backed label is drawn on: its ``ink`` padded on every side, stopped at ``edge_x`` /
    ``edge_y`` on the side of the letters where that edge lies."""

    box = ink.expanded(BACKING_PADDING_EM * font_size)
    min_x, min_y, max_x, max_y = box.min_x, box.min_y, box.max_x, box.max_y
    if edge_y is not None:
        if ink.max_y <= edge_y:
            max_y = min(max_y, edge_y)
        else:
            min_y = max(min_y, edge_y)
    if edge_x is not None:
        if ink.max_x <= edge_x:
            max_x = min(max_x, edge_x)
        else:
            min_x = max(min_x, edge_x)
    return Aabb(min_x, min_y, max_x, max_y)


def backing_rect(box: Aabb, font_size: float, *, opacity: float = BACKING_OPACITY) -> Rect:
    """The drawn backing of ``box``."""

    radius = BACKING_PADDING_EM * font_size
    return Rect(
        insert=(box.min_x, box.min_y),
        size=(box.width, box.height),
        rx=radius,
        ry=radius,
        fill=BACKING_FILL,
        fill_opacity=opacity,
        stroke="none",
    )


def feature_label_backing(label: FeatureLabelInk, leader: Sequence[tuple[float, float]]) -> Aabb:
    """A feature label's box in the frame of its ``ink``: stopped short of its underline by the line's gap
    when it has one, else at the end of its drawn ``leader`` on the side the leader arrives from."""

    line = label.underline
    if line is not None and line[1] != line[2]:
        line_y = label.to_label_frame(line[0])[1]
        gap = UNDERLINE_GAP_EM * label.font_size
        edge = line_y - gap if line_y > label.ink.max_y else line_y + gap
        return backing_box(label.ink, label.font_size, edge_y=edge)
    if not leader:
        return backing_box(label.ink, label.font_size)
    end_x, end_y = label.to_label_frame(leader[-1])
    ink = label.ink
    out_x = max(ink.min_x - end_x, end_x - ink.max_x)
    out_y = max(ink.min_y - end_y, end_y - ink.max_y)
    if out_x >= out_y:
        edge = max(end_x, ink.max_x) if end_x > 0.5 * (ink.min_x + ink.max_x) else min(end_x, ink.min_x)
        return backing_box(ink, label.font_size, edge_x=edge)
    edge = max(end_y, ink.max_y) if end_y > 0.5 * (ink.min_y + ink.max_y) else min(end_y, ink.min_y)
    return backing_box(ink, label.font_size, edge_y=edge)


def feature_label_box(label: FeatureLabelInk, leader: Sequence[tuple[float, float]], feature_id: str | None) -> Rect:
    """The drawn box of a feature label (``protection = box``), turned with its letters; ``leader`` is the
    leader as drawn. Like a halo or a casing it takes no pointer events, and it is bound to its label
    (``feature_id``, required: ``bind_label_part`` refuses none) so that it is hidden with it."""

    rect = backing_rect(feature_label_backing(label, leader), label.font_size, opacity=BOX_OPACITY)
    if label.rotation_deg != 0.0:
        rect.rotate(label.rotation_deg, center=label.center)
    rect.attribs["aria-hidden"] = "true"
    rect.attribs["pointer-events"] = "none"
    return bind_label_part(rect, feature_id)


__all__ = [
    "BACKING_FILL",
    "BACKING_OPACITY",
    "BACKING_PADDING_EM",
    "backing_box",
    "backing_rect",
    "feature_label_backing",
    "feature_label_box",
    "label_halo",
    "underline_casing",
    "UnderlineCasing",
]
