"""A drawn feature label with a leader: its letters' box, its underline, and where its leader ends.

The geometry owner that both drawers (underline, leader end, box backing)
and the Linear region-label placer (the underline as part of the letters) read
(LP-D88, LP-D90, LP-D98). The underline takes no room in a Linear label's row
box.
With an underline (``label_style.FeatureLabelStyle.draws_underline``: the
style's ``line`` is ``under`` and the label is not turned) the leader
ends where it meets the line, at the point of the line nearest the
leader's previous point (the near end for a leader that arrives from the side,
a T for one from straight below); a Circular leader meets it at the end where
it ends today, with its last segment re-aimed there and no added bend (LP-D90).
The line runs along the letters' ink, ``UNDERLINE_GAP_EM`` below the font's
descent line, measured to the line's edge. A leader that would reach the line
through the letters (it arrives from above, as for a Linear label below its
feature) draws no line and stops ``GAP_STOP_EM`` above the letters.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from ..core.text import dominant_baseline_shift_em, get_font_vertical_metrics, text_ink_offsets
from ..layout.spatial import Aabb
from ..layout.text_geometry import Point, TextAnchor, normalize_text_anchor
from .circular_leaders import leader_polyline, rim_text_placement
from .circular_types import CircularLabel
from .label_style import GAP_STOP_EM, UNDERLINE_GAP_EM, UNDERLINE_TRIM_EM, FeatureLabelStyle, leader_width_px


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
    engine places, which ends at the letters. ``under_y`` / ``stop_y`` are the
    y (label frame) of its underline and of its leader's stop above the
    letters; None without an underline. ``join_at_end``: the leader meets the
    line at the end on the side where it ends today, never inside it.
    """

    ink: Aabb
    font_size: float
    rotation_deg: float
    center: Point
    leader: tuple[Point, ...]
    under_y: float | None = None
    stop_y: float | None = None
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
        """(where the leader ends, one end of the line, its other end) in the drawing frame; the three points
        are equal when the leader stops above the letters (no line). None without an underline. The line is
        the letters' ink less ``UNDERLINE_TRIM_EM`` at each end (at most a quarter of it), and the leader
        ends on it (LP-D103)."""

        if self.under_y is None or self.stop_y is None or len(self.leader) < 2:
            return None
        previous = self.to_label_frame(self.leader[-2])
        meet_x = min(max(previous[0], self.ink.min_x), self.ink.max_x)
        side_x = self.to_label_frame(self.leader[-1])[0] if self.join_at_end else meet_x
        near_x, far_x = (
            (self.ink.min_x, self.ink.max_x)
            if abs(side_x - self.ink.min_x) <= abs(side_x - self.ink.max_x)
            else (self.ink.max_x, self.ink.min_x)
        )
        if self.join_at_end:
            # The end on the side where the engine's leader ends today (the text anchor's side).
            meet_x = near_x
        if segment_through_box(previous, (meet_x, self.under_y), self.ink):
            # The leader arrives from above: no line; it stops in the gap above the letters.
            stop = self.from_label_frame((meet_x, self.stop_y))
            return (stop, stop, stop)
        trim = min(UNDERLINE_TRIM_EM * self.font_size, 0.25 * abs(far_x - near_x))
        toward_far = math.copysign(trim, far_x - near_x)
        near_x, far_x = near_x + toward_far, far_x - toward_far
        meet_x = min(max(meet_x, min(near_x, far_x)), max(near_x, far_x))
        return (
            self.from_label_frame((meet_x, self.under_y)),
            self.from_label_frame((near_x, self.under_y)),
            self.from_label_frame((far_x, self.under_y)),
        )

    @property
    def underline_box(self) -> Aabb | None:
        """The drawn underline's box with its stroke, in the drawing frame; None when no line is drawn."""

        line = self.underline
        if line is None or line[1] == line[2]:
            return None
        (x1, y1), (x2, y2) = line[1], line[2]
        return Aabb(x1, y1, x2, y2).expanded(0.5 * leader_width_px(self.font_size))


def feature_label_ink(
    text: str,
    font_family: str,
    font_size: float,
    *,
    at: Point,
    anchor: TextAnchor,
    dominant_baseline: str,
    rotation_deg: float,
    turned: bool,
    leader: Sequence[Point],
    style: FeatureLabelStyle,
    join_at_end: bool = False,
) -> FeatureLabelInk:
    """A drawn feature label (``text`` drawn at ``at``, turned ``rotation_deg`` about it), with an underline
    when ``style`` draws one under it and the label is not ``turned``."""

    size = float(font_size)
    offsets = text_ink_offsets(text, font_family, size, dominant_baseline=dominant_baseline, text_anchor=anchor)
    if offsets is None:
        # The fallback width of calculate_bbox_dimensions.
        width = 0.6 * len(text) * size
        origin = -{"start": 0.0, "middle": 0.5, "end": 1.0}[anchor] * width
        offsets = (origin, origin + width, -0.5 * size, 0.5 * size)
    left, right, top, bottom = offsets
    x, y = float(at[0]), float(at[1])
    under_y = stop_y = None
    if style.draws_underline(turned):
        # The gap is measured to the line's edge, so its centre lies half its width further out.
        reach = UNDERLINE_GAP_EM * size + 0.5 * leader_width_px(size)
        metrics = get_font_vertical_metrics(font_family)
        shift = None if metrics is None else dominant_baseline_shift_em(dominant_baseline, metrics)
        under_y = (
            y + (shift + metrics.descent) * size + reach
            if metrics is not None and shift is not None
            else y + bottom + reach
        )
        stop_y = y + top - GAP_STOP_EM * size
    return FeatureLabelInk(
        ink=Aabb(x + left, y + top, x + right, y + bottom),
        font_size=size,
        rotation_deg=float(rotation_deg),
        center=(x, y),
        leader=tuple((float(px), float(py)) for px, py in leader),
        under_y=under_y,
        stop_y=stop_y,
        join_at_end=join_at_end,
    )


def linear_label_ink(label: Mapping[str, object], leader: Sequence[Point], style: FeatureLabelStyle) -> FeatureLabelInk:
    """A drawn Linear feature label (``dominant-baseline=central``); turned when given a rotation."""

    rotation_deg = float(str(label.get("rotation_deg", 0.0)))
    return feature_label_ink(
        str(label["label_text"]),
        str(label["font_family"]),
        float(str(label["font_size"])),
        at=(float(str(label["middle_x"])), float(str(label["middle_y"]))),
        anchor=normalize_text_anchor(label.get("text_anchor", "middle")),
        dominant_baseline="central",
        rotation_deg=rotation_deg,
        turned=rotation_deg % 360.0 != 0.0,
        leader=leader,
        style=style,
    )


def circular_label_ink(label: CircularLabel, record_length: int, style: FeatureLabelStyle) -> FeatureLabelInk:
    """A drawn outer or inner Circular label; a radial one is turned, even where it reads horizontally."""

    (x, y), anchor, baseline, rotation_deg = rim_text_placement(label, record_length)
    return feature_label_ink(
        str(label["label_text"]),
        str(label["font_family"]),
        float(label["font_size"]),
        at=(x, y),
        anchor=normalize_text_anchor(anchor),
        dominant_baseline=baseline,
        rotation_deg=rotation_deg,
        turned=label.get("placement") == "radial",
        leader=leader_polyline(label),
        style=style,
        join_at_end=True,
    )


__all__ = [
    "FeatureLabelInk",
    "circular_label_ink",
    "feature_label_ink",
    "linear_label_ink",
    "segment_through_box",
]
