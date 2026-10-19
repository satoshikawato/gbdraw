"""The translucent backing of a label: one geometry for region and feature labels.

A backed label is drawn on a white box with rounded corners at 70 % opacity
(LP-D82 ``RGN-FRONT``), ``BACKING_PADDING_EM`` font sizes beyond its letters
on every side, with the same corner radius. A box may stop at an edge on one
side (a region label's mark, a feature label's own leader end), so it never
covers what it belongs to.

Study prototype (LP-D85): feature labels take the same backing where
``GBDRAW_LABEL_BACKING`` asks for it: ``through`` backs a label that another
label's leader passes through, ``all`` backs every label with a leader.
"""

from __future__ import annotations

import math
import os
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal, TypeAlias

from gbdraw.layout.spatial import Aabb
from gbdraw.layout.text_geometry import Point, TextAnchor
from gbdraw.render.drawers.annotation_marks import label_ink_offsets
from gbdraw.svg.elements import Rect

BACKING_PADDING_EM = 0.18
BACKING_FILL = "#ffffff"
BACKING_OPACITY = 0.7

FeatureLabelBackingMode: TypeAlias = Literal["off", "through", "all"]
BACKING_MODE_ENV = "GBDRAW_LABEL_BACKING"


def backing_box(ink: Aabb, font_size: float, *, edge_x: float | None = None, edge_y: float | None = None) -> Aabb:
    """The box a backed label is drawn on: its ``ink`` padded on every side, stopped at ``edge_x`` /
    ``edge_y`` on the side of the letters where that edge lies (past the letters' far side when the
    edge is not before them)."""

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


def backing_rect(box: Aabb, font_size: float) -> Rect:
    """The drawn backing of ``box``."""

    radius = BACKING_PADDING_EM * font_size
    return Rect(
        insert=(box.min_x, box.min_y),
        size=(box.width, box.height),
        rx=radius,
        ry=radius,
        fill=BACKING_FILL,
        fill_opacity=BACKING_OPACITY,
        stroke="none",
    )


def feature_label_backing_mode() -> FeatureLabelBackingMode:
    """The prototype switch (``GBDRAW_LABEL_BACKING``); ``off`` when unset or unknown."""

    value = os.environ.get(BACKING_MODE_ENV, "").strip().lower()
    if value == "through":
        return "through"
    if value == "all":
        return "all"
    return "off"


@dataclass(frozen=True)
class FeatureLabelInk:
    """One feature label with a leader, in its record's drawing frame.

    ``ink`` is the box of its letters before the ``rotation_deg`` turn about
    ``center`` (the point its text is drawn at); ``leader`` is its drawn leader
    polyline, which ends at its letters.
    """

    ink: Aabb
    font_size: float
    rotation_deg: float
    center: Point
    leader: tuple[Point, ...]

    def to_label_frame(self, point: Point) -> Point:
        """``point`` turned back by ``rotation_deg`` about ``center`` (the frame of ``ink``)."""

        if self.rotation_deg == 0.0:
            return point
        angle = math.radians(-self.rotation_deg)
        cx, cy = self.center
        dx, dy = point[0] - cx, point[1] - cy
        return (cx + dx * math.cos(angle) - dy * math.sin(angle), cy + dx * math.sin(angle) + dy * math.cos(angle))

    @property
    def backing(self) -> Aabb:
        """Its backing in the frame of ``ink``, stopped at its own leader's end on the side it arrives from."""

        if not self.leader:
            return backing_box(self.ink, self.font_size)
        end_x, end_y = self.to_label_frame(self.leader[-1])
        ink = self.ink
        out_x = max(ink.min_x - end_x, end_x - ink.max_x)
        out_y = max(ink.min_y - end_y, end_y - ink.max_y)
        if out_x >= out_y:
            edge = max(end_x, ink.max_x) if end_x > 0.5 * (ink.min_x + ink.max_x) else min(end_x, ink.min_x)
            return backing_box(ink, self.font_size, edge_x=edge)
        edge = max(end_y, ink.max_y) if end_y > 0.5 * (ink.min_y + ink.max_y) else min(end_y, ink.min_y)
        return backing_box(ink, self.font_size, edge_y=edge)

    def backing_rect(self) -> Rect:
        rect = backing_rect(self.backing, self.font_size)
        if self.rotation_deg != 0.0:
            rect.rotate(self.rotation_deg, center=self.center)
        return rect


def feature_label_ink(
    text: str,
    font_family: str,
    font_size: float,
    *,
    at: Point,
    anchor: TextAnchor,
    dominant_baseline: str,
    rotation_deg: float,
    leader: Sequence[Point],
) -> FeatureLabelInk:
    """A drawn feature label (``text`` drawn at ``at``, turned ``rotation_deg`` about it) as its backing sees it."""

    left, right, top, bottom = label_ink_offsets(
        text, font_family, font_size, anchor, dominant_baseline=dominant_baseline
    )
    x, y = at
    return FeatureLabelInk(
        ink=Aabb(x + left, y + top, x + right, y + bottom),
        font_size=float(font_size),
        rotation_deg=float(rotation_deg),
        center=(float(x), float(y)),
        leader=tuple((float(px), float(py)) for px, py in leader),
    )


def backing_rects(labels: Sequence[FeatureLabelInk], mode: FeatureLabelBackingMode) -> list[Rect]:
    """The backings of one record's feature labels, to draw before (under) all of their letters."""

    backed = backed_feature_labels(labels, mode)
    return [label.backing_rect() for index, label in enumerate(labels) if index in backed]


def _segment_through_box(start: Point, end: Point, box: Aabb) -> bool:
    """Whether the segment passes through the interior of ``box`` (a touch of its edge does not count)."""

    t0, t1 = 0.0, 1.0
    for origin, delta, low, high in (
        (start[0], end[0] - start[0], box.min_x, box.max_x),
        (start[1], end[1] - start[1], box.min_y, box.max_y),
    ):
        if abs(delta) < 1e-12:
            if origin <= low or origin >= high:
                return False
            continue
        ta, tb = (low - origin) / delta, (high - origin) / delta
        t0, t1 = max(t0, min(ta, tb)), min(t1, max(ta, tb))
        if t1 - t0 <= 1e-9:
            return False
    return True


def leaders_through_text(labels: Sequence[FeatureLabelInk]) -> list[tuple[int, int]]:
    """(leader of, text of) pairs: a label's leader drawn through another label's letters."""

    found: list[tuple[int, int]] = []
    for text_of, label in enumerate(labels):
        for leader_of, other in enumerate(labels):
            if leader_of == text_of or not other.leader:
                continue
            points = [label.to_label_frame(point) for point in other.leader]
            if any(_segment_through_box(a, b, label.ink) for a, b in zip(points, points[1:])):
                found.append((leader_of, text_of))
    return found


def backed_feature_labels(labels: Sequence[FeatureLabelInk], mode: FeatureLabelBackingMode) -> frozenset[int]:
    """The labels drawn on a backing: none, those another leader passes through, or every one with a leader."""

    if mode == "all":
        return frozenset(index for index, label in enumerate(labels) if label.leader)
    if mode == "through":
        return frozenset(text_of for _leader_of, text_of in leaders_through_text(labels))
    return frozenset()


__all__ = [
    "BACKING_FILL",
    "BACKING_MODE_ENV",
    "BACKING_OPACITY",
    "BACKING_PADDING_EM",
    "FeatureLabelBackingMode",
    "FeatureLabelInk",
    "backed_feature_labels",
    "backing_box",
    "backing_rect",
    "backing_rects",
    "feature_label_ink",
    "feature_label_backing_mode",
    "leaders_through_text",
]
