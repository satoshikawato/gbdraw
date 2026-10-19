"""A drawn feature label with a leader, as its backing, its underline and the fit count see it.

Study prototype (LP-D85, LP-D88). ``GBDRAW_LABEL_UNDERLINE=on`` draws a line
under the label's letters that continues its leader (the LoVis4u convention):
the leader ends where it meets the line, at the point of the line nearest the
leader's previous point (the near end for a leader that arrives from the side,
a T for a leader that arrives from straight below or above; a Circular leader
always meets it at the end where it ends today, LP-D90). The line runs
along the letters' ink, one ``UNDERLINE_GAP_EM`` beyond the font's descent
line, or beyond its ascent line when the leader would otherwise run through
the letters (it arrives from above, as for a Linear label below its feature). The drawers and the fit count
read the same polyline (``drawn``): the line is part of the label's leader for
crossings. ``draw`` draws the line without counting it (a study view).
"""

from __future__ import annotations

import math
import os
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from ..core.text import dominant_baseline_shift_em, get_font_vertical_metrics, text_ink_offsets
from ..layout.spatial import Aabb
from ..layout.text_geometry import Point, TextAnchor, normalize_text_anchor
from .circular_leaders import leader_polyline, rim_text_placement
from .circular_types import CircularLabel

UNDERLINE_ENV = "GBDRAW_LABEL_UNDERLINE"
UNDERLINE_GAP_EM = 0.06


def label_underline_on() -> bool:
    """The prototype switch (``GBDRAW_LABEL_UNDERLINE``): ``on`` draws and counts the underline,
    ``draw`` only draws it (the layout and the fit count stay as without one; a study view)."""

    return os.environ.get(UNDERLINE_ENV, "").strip().lower() in {"on", "draw"}


def label_underline_counted() -> bool:
    """Whether the fit count (and so the engine that reads it) sees the underline as part of the leader."""

    return os.environ.get(UNDERLINE_ENV, "").strip().lower() == "on"


def segment_through_box(start: Point, end: Point, box: Aabb) -> bool:
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


@dataclass(frozen=True)
class FeatureLabelInk:
    """One feature label with a leader, in its record's drawing frame.

    ``ink`` is the box of its letters before the ``rotation_deg`` turn about
    ``center`` (the point its text is drawn at); ``leader`` is the leader the
    engine places, which ends at the letters. ``under_y`` / ``over_y`` are the
    y (label frame) of an underline below or above the letters; None without one.
    ``join_at_end``: the leader meets the line at the end on the side where it
    ends today, never inside it (Circular, LP-D90: its last segment is re-aimed
    there, no bend is added).
    """

    ink: Aabb
    font_size: float
    rotation_deg: float
    center: Point
    leader: tuple[Point, ...]
    under_y: float | None = None
    over_y: float | None = None
    join_at_end: bool = False

    def _turned(self, point: Point, degrees: float) -> Point:
        if degrees == 0.0:
            return point
        angle = math.radians(degrees)
        cx, cy = self.center
        dx, dy = point[0] - cx, point[1] - cy
        return (cx + dx * math.cos(angle) - dy * math.sin(angle), cy + dx * math.sin(angle) + dy * math.cos(angle))

    def to_label_frame(self, point: Point) -> Point:
        """``point`` turned back by ``rotation_deg`` about ``center`` (the frame of ``ink``)."""

        return self._turned(point, -self.rotation_deg)

    def from_label_frame(self, point: Point) -> Point:
        return self._turned(point, self.rotation_deg)

    @property
    def underline(self) -> tuple[Point, Point, Point] | None:
        """(where the leader meets the line, one end, the other end) in the drawing frame; None without one."""

        if self.under_y is None or self.over_y is None or len(self.leader) < 2:
            return None
        previous = self.to_label_frame(self.leader[-2])
        meet_x = min(max(previous[0], self.ink.min_x), self.ink.max_x)
        near_x, far_x = (
            (self.ink.min_x, self.ink.max_x)
            if abs(meet_x - self.ink.min_x) <= abs(meet_x - self.ink.max_x)
            else (self.ink.max_x, self.ink.min_x)
        )
        if self.join_at_end:
            # The end on the side where the engine's leader ends today (the text anchor's side).
            end_x = self.to_label_frame(self.leader[-1])[0]
            near_x, far_x = (
                (self.ink.min_x, self.ink.max_x)
                if abs(end_x - self.ink.min_x) <= abs(end_x - self.ink.max_x)
                else (self.ink.max_x, self.ink.min_x)
            )
            meet_x = near_x
        # Under the letters, unless the leader would then run through them (it arrives from above): over them.
        y = self.over_y if segment_through_box(previous, (meet_x, self.under_y), self.ink) else self.under_y
        return (
            self.from_label_frame((meet_x, y)),
            self.from_label_frame((near_x, y)),
            self.from_label_frame((far_x, y)),
        )

    @property
    def drawn(self) -> tuple[Point, ...]:
        """The drawn leader as one polyline: the engine's leader, or with an underline the leader to the
        line, then the line (through its near end to its far end)."""

        line = self.underline
        if line is None:
            return self.leader
        meet, near, far = line
        points = [*self.leader[:-1], meet]
        if math.dist(meet, near) > 1e-9:
            points.append(near)
        points.append(far)
        return tuple(points)


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
    underline: bool = False,
    join_at_end: bool = False,
) -> FeatureLabelInk:
    """A drawn feature label (``text`` drawn at ``at``, turned ``rotation_deg`` about it)."""

    size = float(font_size)
    offsets = text_ink_offsets(text, font_family, size, dominant_baseline=dominant_baseline, text_anchor=anchor)
    if offsets is None:
        width = 0.6 * len(text) * size
        origin = -{"start": 0.0, "middle": 0.5, "end": 1.0}[anchor] * width
        offsets = (origin, origin + width, -0.5 * size, 0.5 * size)
    left, right, top, bottom = offsets
    x, y = float(at[0]), float(at[1])
    under_y = over_y = None
    if underline:
        metrics = get_font_vertical_metrics(font_family)
        shift = None if metrics is None else dominant_baseline_shift_em(dominant_baseline, metrics)
        if metrics is not None and shift is not None:
            under_y = y + (shift + metrics.descent + UNDERLINE_GAP_EM) * size
            over_y = y + (shift - metrics.ascent - UNDERLINE_GAP_EM) * size
        else:
            under_y, over_y = y + bottom + UNDERLINE_GAP_EM * size, y + top - UNDERLINE_GAP_EM * size
    return FeatureLabelInk(
        ink=Aabb(x + left, y + top, x + right, y + bottom),
        font_size=size,
        rotation_deg=float(rotation_deg),
        center=(x, y),
        leader=tuple((float(px), float(py)) for px, py in leader),
        under_y=under_y,
        over_y=over_y,
        join_at_end=join_at_end,
    )


def linear_label_ink(label: Mapping[str, object], leader: Sequence[Point], *, underline: bool) -> FeatureLabelInk:
    """A drawn Linear feature label (``dominant-baseline=central``)."""

    return feature_label_ink(
        str(label["label_text"]),
        str(label["font_family"]),
        float(str(label["font_size"])),
        at=(float(str(label["middle_x"])), float(str(label["middle_y"]))),
        anchor=normalize_text_anchor(label.get("text_anchor", "middle")),
        dominant_baseline="central",
        rotation_deg=float(str(label.get("rotation_deg", 0.0))),
        leader=leader,
        underline=underline,
    )


def circular_label_ink(
    label: CircularLabel, record_length: int, font_size: float, font_family: str, *, underline: bool
) -> FeatureLabelInk:
    """A drawn outer or inner Circular label."""

    (x, y), anchor, baseline, rotation_deg = rim_text_placement(label, record_length)
    return feature_label_ink(
        str(label["label_text"]),
        font_family,
        font_size,
        at=(x, y),
        anchor=normalize_text_anchor(anchor),
        dominant_baseline=baseline,
        rotation_deg=rotation_deg,
        leader=leader_polyline(label),
        underline=underline,
        join_at_end=True,
    )


__all__ = [
    "UNDERLINE_ENV",
    "UNDERLINE_GAP_EM",
    "segment_through_box",
    "FeatureLabelInk",
    "circular_label_ink",
    "feature_label_ink",
    "label_underline_counted",
    "label_underline_on",
    "linear_label_ink",
]
