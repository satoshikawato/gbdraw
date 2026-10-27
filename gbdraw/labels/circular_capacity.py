"""Room for one record's horizontal Circular labels, read before they are placed (PD-OI-097).

Horizontal outer and inner labels are anchored on an ellipse (:class:`LabelArc`)
and their boxes run away from (outer) or toward (inner) the vertical axis.
Walking the labels around the ellipse, two neighbours that do not overlap are
apart either vertically by their heights or horizontally by a width. The
anchors travel ``4 * y_radius`` vertically and ``4 * x_radius`` horizontally
around the ellipse, so the labels need

    needed = sum(heights) - (heights the horizontal travel saves) - 2 * max(heights)

of the vertical room ``free = 4 * y_radius``. The saving is a fractional
knapsack of the horizontal travel over the labels of the smallest width per
height; ``2 * max(heights)`` allows the two crossings at the poles. With every
anchor on the ellipse this is a lower bound: ``ratio = needed / free`` over 1
cannot be placed without overlaps. The repair and lane clearance move anchors
off the ellipse (up to about 1.2 times in the calibration), which the stop
ratio allows for.

The lever is the same ratio after the strongest single change a fix row can
make: the font row at its floor, the outward rows at their limit (outer side),
and putting the inner labels outside (inner side). A side stops at once
(CAPACITY-STOP) when its ratio is over ``CAPACITY_STOP_RATIO`` and its lever over
``CAPACITY_LEVER_RATIO``; it warns before placement (CAPACITY-WARN) from
``CAPACITY_WARN_RATIO``. Radial and embedded labels are not read here.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from contextvars import ContextVar
import warnings
from dataclasses import dataclass
from typing import Literal, TypeAlias

from .circular_fit import FONT_LIMIT_PT, OFFSET_LIMIT
from .circular_types import LABEL_FIT_X_MARGIN_PX, LABEL_FIT_Y_MARGIN_PX, CircularLabel
from ..exceptions import LabelsMayNotFitWarning
from .fit import LabelCapacityRecord, LabelCapacityRemedy, label_capacity_sentence

# A side whose labels need more than this many times the room stops at once when its lever needs more than
# CAPACITY_LEVER_RATIO times too (CAPACITY-STOP); from CAPACITY_WARN_RATIO it warns before placement (CAPACITY-WARN).
CAPACITY_WARN_RATIO = 0.8
CAPACITY_STOP_RATIO = 1.5
CAPACITY_LEVER_RATIO = 1.0

CapacityLevel: TypeAlias = Literal["fits", "near", "over"]
HorizontalLabelKind: TypeAlias = Literal["circular_outer", "circular_inner"]

# Off while a fit row's settings are verified with the full assembly, so the rows are computed as before.
_CHECKED: ContextVar[bool] = ContextVar("gbdraw_circular_label_capacity_checked", default=True)


@dataclass(frozen=True)
class LabelArc:
    """The ellipse one side's horizontal labels are placed on, in px around the record's centre."""

    x_radius: float
    y_radius: float
    center_x: float


@dataclass(frozen=True)
class CircularLabelCapacity:
    """One side's room: the ratio of the current settings, the lever, and the level they give."""

    label_kind: HorizontalLabelKind
    record: LabelCapacityRecord
    lever: float
    level: CapacityLevel


class CircularLabelCapacityExceeded(Exception):
    """A record's horizontal labels need much more room than the map has, so they were not placed (CAPACITY-STOP)."""

    def __init__(self, capacities: Sequence[CircularLabelCapacity]) -> None:
        super().__init__("Circular labels need much more room than the map has.")
        self.capacities = tuple(capacities)


@contextmanager
def label_capacity_unchecked() -> Iterator[None]:
    """Place labels without the capacity check: a fit row's verification keeps the rows it found before the check."""

    token = _CHECKED.set(False)
    try:
        yield
    finally:
        _CHECKED.reset(token)


def label_capacity_checked() -> bool:
    """Whether a drawn placement reads its capacity first (not during a fit row's verification)."""

    return _CHECKED.get()


def _needed_px(labels: Sequence[CircularLabel], arc: LabelArc, scale: float) -> float:
    sizes = [(scale * float(label["width_px"]), scale * float(label["height_px"])) for label in labels]
    if not sizes:
        return 0.0
    heights = [height + LABEL_FIT_Y_MARGIN_PX for _, height in sizes]
    # Horizontal travel saves the heights of the labels that cost the least width per height.
    budget = max(0.0, 4.0 * abs(arc.x_radius))
    saved = 0.0
    for width, height in sorted(sizes, key=lambda size: (size[0] + LABEL_FIT_X_MARGIN_PX) / max(size[1], 1e-9)):
        if budget <= 0.0:
            break
        box_width = width + LABEL_FIT_X_MARGIN_PX
        take = min(1.0, budget / box_width)
        saved += take * (height + LABEL_FIT_Y_MARGIN_PX)
        budget -= take * box_width
    # The two neighbours that meet across the poles run away from each other.
    if len(heights) > 2:
        saved += 2.0 * max(heights)
    return max(0.0, sum(heights) - saved)


def _ratio(labels: Sequence[CircularLabel], arc: LabelArc, *, scale: float = 1.0) -> float:
    free = 4.0 * abs(arc.y_radius)
    needed = _needed_px(labels, arc, scale)
    if free <= 0.0:
        return float("inf") if needed > 0.0 else 0.0
    return needed / free


def _widened(arc: LabelArc, x_offset: float, y_offset: float) -> LabelArc:
    """The outer ellipse with both outer offsets at the outward rows' limit (the radii scale with the offsets)."""

    def factor(offset: float) -> float:
        return max(1.0, OFFSET_LIMIT / offset) if offset > 0.0 else 1.0

    return LabelArc(arc.x_radius * factor(x_offset), arc.y_radius * factor(y_offset), arc.center_x)


def _remedies(font_size: float) -> tuple[LabelCapacityRemedy, ...]:
    if font_size > FONT_LIMIT_PT:
        return ("filter", "label_text", "font", "linear")
    # The font row stops at its floor, so a smaller font is not offered from there.
    return ("filter", "label_text", "linear")


def _level(ratio: float, lever: float) -> CapacityLevel:
    if ratio > CAPACITY_STOP_RATIO and lever > CAPACITY_LEVER_RATIO:
        return "over"
    if ratio >= CAPACITY_WARN_RATIO:
        return "near"
    return "fits"


def _capacity(
    label_kind: HorizontalLabelKind, labels: Sequence[CircularLabel], ratio: float, lever: float, font_size: float
) -> CircularLabelCapacity:
    record = LabelCapacityRecord(len(labels), float(font_size), ratio, _remedies(font_size))
    return CircularLabelCapacity(label_kind, record, lever, _level(ratio, lever))


def circular_label_capacity(
    outer_labels: Sequence[CircularLabel],
    inner_labels: Sequence[CircularLabel],
    *,
    outer_arc: LabelArc,
    inner_arc: LabelArc,
    outer_offsets: tuple[float, float],
    font_size: float,
) -> tuple[CircularLabelCapacity, ...]:
    """The capacity of each side with labels: outer first, then inner.

    ``outer_offsets`` are the outer labels' X and Y radius offsets, which the
    outward rows raise to their limit.
    """

    font_scale = min(font_size, FONT_LIMIT_PT) / font_size if font_size > 0.0 else 1.0
    capacities: list[CircularLabelCapacity] = []
    if outer_labels:
        ratio = _ratio(outer_labels, outer_arc)
        lever = min(
            _ratio(outer_labels, outer_arc, scale=font_scale),
            _ratio(outer_labels, _widened(outer_arc, *outer_offsets)),
        )
        capacities.append(_capacity("circular_outer", outer_labels, ratio, lever, font_size))
    if inner_labels:
        ratio = _ratio(inner_labels, inner_arc)
        lever = min(
            _ratio(inner_labels, inner_arc, scale=font_scale),
            # All labels outside: the inner labels join the outer ones on the outer ellipse.
            _ratio([*outer_labels, *inner_labels], outer_arc),
        )
        capacities.append(_capacity("circular_inner", inner_labels, ratio, lever, font_size))
    return tuple(capacities)


def check_label_capacity(record_id: str, capacities: Sequence[CircularLabelCapacity]) -> None:
    """Stop a record whose labels need much more room than the map has, or warn for the sides near the limit.

    Raises :class:`CircularLabelCapacityExceeded` with every side over the limit
    (CAPACITY-STOP); otherwise warns ``LabelsMayNotFitWarning`` once per side near it
    (CAPACITY-WARN), before the labels are placed.
    """

    over = [capacity for capacity in capacities if capacity.level == "over"]
    if over:
        raise CircularLabelCapacityExceeded(over)
    for capacity in capacities:
        if capacity.level == "near":
            message = label_capacity_sentence(record_id, capacity.label_kind, capacity.record, placed=True)
            warnings.warn(
                LabelsMayNotFitWarning(
                    message, record_id=record_id, label_kind=capacity.label_kind, capacity=capacity.record
                ),
                stacklevel=2,
            )


__all__ = [
    "CAPACITY_LEVER_RATIO",
    "CAPACITY_STOP_RATIO",
    "CAPACITY_WARN_RATIO",
    "CapacityLevel",
    "CircularLabelCapacity",
    "CircularLabelCapacityExceeded",
    "HorizontalLabelKind",
    "LabelArc",
    "check_label_capacity",
    "circular_label_capacity",
    "label_capacity_checked",
    "label_capacity_unchecked",
]
