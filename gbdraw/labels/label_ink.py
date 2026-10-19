"""A drawn feature label with a leader, as its backing, its underline and the fit count see it.

Study prototype (LP-D85, LP-D88). ``GBDRAW_LABEL_UNDERLINE=on`` draws a line
under the label's letters that continues its leader (the LoVis4u convention):
the leader ends where it meets the line, at the point of the line nearest the
leader's previous point (the near end for a leader that arrives from the side,
a T for a leader that arrives from straight below or above). The line runs
along the letters' ink, one ``UNDERLINE_GAP_EM`` beyond the font's descent
line, or beyond its ascent line when the leader arrives from straight above
the letters (a Linear label below its feature). The drawers and the fit count
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


@dataclass(frozen=True)
class FeatureLabelInk:
    """One feature label with a leader, in its record's drawing frame.

    ``ink`` is the box of its letters before the ``rotation_deg`` turn about
    ``center`` (the point its text is drawn at); ``leader`` is the leader the
    engine places, which ends at the letters. ``under_y`` / ``over_y`` are the
    y (label frame) of an underline below or above the letters; None without one.
    """

    ink: Aabb
    font_size: float
    rotation_deg: float
    center: Point
    leader: tuple[Point, ...]
    under_y: float | None = None
    over_y: float | None = None

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
        # Over the letters only for a leader that arrives from straight above them (a T); a leader from
        # below or from the side meets the line under them.
        from_above = previous[1] < self.ink.min_y and self.ink.min_x <= previous[0] <= self.ink.max_x
        y = self.over_y if from_above else self.under_y
        meet_x = min(max(previous[0], self.ink.min_x), self.ink.max_x)
        near_x, far_x = (
            (self.ink.min_x, self.ink.max_x)
            if abs(meet_x - self.ink.min_x) <= abs(meet_x - self.ink.max_x)
            else (self.ink.max_x, self.ink.min_x)
        )
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
    )


__all__ = [
    "UNDERLINE_ENV",
    "UNDERLINE_GAP_EM",
    "FeatureLabelInk",
    "circular_label_ink",
    "feature_label_ink",
    "label_underline_counted",
    "label_underline_on",
    "linear_label_ink",
]
