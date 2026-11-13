"""Where a Linear label drawn inside its feature fits (Owner D-161, OV-3380, OV-3397).

A label inside its feature fits where the ink of its letters, with the caller's margin more at each end (half
the room Owner D-20 gives a label inside its feature), lies inside the drawn outline of the block it sits on
(``LinearBlockOutline``, the outline the block's path is drawn from) and clear of every part a later feature
draws over it. It stays at the middle of
its block when it fits there, else it takes the fitting x nearest to the middle; it has no x when none fits.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import NamedTuple, cast

from ..layout.spatial import Interval, IntervalIndex
from ..svg.arrows import LinearBlockOutline

_EPSILON = 1e-9


class InkBox(NamedTuple):
    """The box a label's letters paint, in px from the point it is drawn at (``text_ink_offsets``), y down."""

    left: float
    right: float
    top: float
    bottom: float


# A lane of a record (the y of its top, middle and bottom) and a wrap row in it.
LaneRow = tuple[tuple[float, float, float], int]


class DrawnPartSpans:
    """The x spans of the parts the features of one record draw, by lane and wrap row, with the order the
    features are drawn in."""

    def __init__(self, bucket_size: float) -> None:
        self._bucket_size = bucket_size
        self._indexes: dict[LaneRow, IntervalIndex] = {}
        self._spans: list[tuple[int, float, float]] = []

    def add(self, order: int, lane_row: LaneRow, start_x: float, end_x: float) -> None:
        """Add a part drawn by the ``order``-th feature in ``lane_row`` from ``start_x`` to ``end_x``."""

        low, high = min(start_x, end_x), max(start_x, end_x)
        index = self._indexes.get(lane_row)
        if index is None:
            index = self._indexes[lane_row] = IntervalIndex(self._bucket_size)
        index.insert(len(self._spans), Interval(low, high))
        self._spans.append((order, low, high))

    def drawn_after(self, order: int, lane_row: LaneRow, start_x: float, end_x: float) -> list[tuple[float, float]]:
        """The spans that features drawn after the ``order``-th one draw in ``lane_row`` over
        ``start_x``..``end_x``: their outlines are drawn over that feature."""

        index = self._indexes.get(lane_row)
        if index is None:
            return []
        spans = []
        for span_id in index.query(Interval(start_x, end_x)):
            later, low, high = self._spans[cast(int, span_id)]
            if later > order and low < end_x and high > start_x:
                spans.append((low, high))
        return spans


def _near_edge_range(outline: LinearBlockOutline, ink: InkBox, margin_px: float) -> tuple[float, float] | None:
    """The distances from the block's tail that the near end of the widened ink box may take, the label
    centred on the block's middle y: anywhere along a shaft tall enough for the ink, and in the head only
    as far as its half-height, which falls to 0 at the tip, still holds the ink."""

    length = abs(outline.tip_x - outline.tail_x)
    head = abs(outline.tip_x - outline.shoulder_x)
    above, below = -ink.top, ink.bottom
    shaft_holds = (
        above <= outline.middle_y - outline.shaft_top_y + _EPSILON
        and below <= outline.shaft_bottom_y - outline.middle_y + _EPSILON
    )
    if head > 0.0:
        head_above, head_below = outline.middle_y - outline.top_y, outline.bottom_y - outline.middle_y
        if head_above <= 0.0 or head_below <= 0.0:
            return None
        share = max(above / head_above, below / head_below)
        if share > 1.0 + _EPSILON:
            return None
        far = length - head * share
    elif shaft_holds:
        far = length
    else:
        return None
    near = 0.0 if shaft_holds else length - head
    last_near = far - (ink.right - ink.left + 2.0 * margin_px)
    if last_near < near - _EPSILON:
        return None
    return near, max(near, last_near)


def outline_anchor_range(outline: LinearBlockOutline, ink: InkBox, margin_px: float) -> tuple[float, float] | None:
    """The anchor x range (text-anchor middle) where the widened ink box lies inside ``outline``, or None."""

    near_range = _near_edge_range(outline, ink, margin_px)
    if near_range is None:
        return None
    first, last = near_range
    tail_x = outline.tail_x
    if outline.tip_x >= tail_x:
        return tail_x + first - ink.left + margin_px, tail_x + last - ink.left + margin_px
    return tail_x - last - ink.right - margin_px, tail_x - first - ink.right - margin_px


def embedded_label_anchor_x(
    outline: LinearBlockOutline,
    ink: InkBox,
    margin_px: float,
    centre_x: float,
    later_spans: Sequence[tuple[float, float]],
) -> float | None:
    """The anchor x of a label inside ``outline``: ``centre_x`` when it fits there, else the fitting x
    nearest to it, or None when no x fits. An x fits when the ink box widened by ``margin_px`` at each end
    lies inside the outline and meets none of ``later_spans``, the parts drawn over the block after it."""

    anchor_range = outline_anchor_range(outline, ink, margin_px)
    if anchor_range is None:
        return None
    low, high = anchor_range
    blocked = [(start - ink.right - margin_px, end - ink.left + margin_px) for start, end in later_spans]

    def fits(anchor_x: float) -> bool:
        return low - _EPSILON <= anchor_x <= high + _EPSILON and not any(
            first < anchor_x < last for first, last in blocked
        )

    if fits(centre_x):
        return centre_x
    candidates = [anchor_x for anchor_x in (low, high, *(x for span in blocked for x in span)) if fits(anchor_x)]
    return min(candidates, key=lambda anchor_x: abs(anchor_x - centre_x)) if candidates else None


__all__ = [
    "DrawnPartSpans",
    "InkBox",
    "LaneRow",
    "embedded_label_anchor_x",
    "outline_anchor_range",
]
