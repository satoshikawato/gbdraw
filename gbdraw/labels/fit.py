"""Whether drawn labels fit, the one decision point, and the report type.

Labels fit when no two label boxes overlap and no two leader polylines cross
(LP-D43 ``DEF-CROSS``). The counts take plain boxes and polylines, so every
label kind (each with its own adapter) shares them and the report.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal, TypeAlias

from ..layout.spatial import Aabb, candidate_aabb_pairs
from ..layout.text_geometry import Point, segments_intersect

LabelKind: TypeAlias = Literal["circular_outer"]
LabelFitOutcome: TypeAlias = Literal["fits", "warn", "stop"]
LabelFitRowKind: TypeAlias = Literal[
    "outward_both",
    "outward_vertical",
    "outward_horizontal",
    "smaller_font",
    "tighter_rows",
    "combination",
]
# ``found``: ``values`` fit; ``not_found``: no value within the lever's limit
# (``reason`` says which); ``not_computed``: the bounded search stopped first.
LabelFitRowStatus: TypeAlias = Literal["found", "not_found", "not_computed"]
Polyline: TypeAlias = Sequence[Point]


@dataclass(frozen=True)
class LabelFit:
    """Overlapping label-box pairs and crossing leader pairs; the labels fit when both are 0."""

    overlap_pairs: int
    crossing_pairs: int

    @property
    def fits(self) -> bool:
        return self.overlap_pairs == 0 and self.crossing_pairs == 0


def _closed_overlap(first: Aabb, second: Aabb) -> bool:
    return (
        first.min_x <= second.max_x
        and second.min_x <= first.max_x
        and first.min_y <= second.max_y
        and second.min_y <= first.max_y
    )


def label_box_overlap_count(boxes: Sequence[Aabb]) -> int:
    """Count box pairs that overlap or touch; callers include each kind's margins in the boxes."""

    extent = max((max(box.width, box.max_y - box.min_y) for box in boxes), default=16.0)
    bucket_size = max(16.0, min(256.0, extent))
    return sum(
        1
        for left, right in candidate_aabb_pairs(boxes, bucket_size=bucket_size)
        if _closed_overlap(boxes[left], boxes[right])
    )


def _polyline_aabb(polyline: Polyline) -> Aabb:
    return Aabb(
        min(point[0] for point in polyline),
        min(point[1] for point in polyline),
        max(point[0] for point in polyline),
        max(point[1] for point in polyline),
    )


def polyline_crossing_count(polylines: Sequence[Polyline]) -> int:
    """Count polyline pairs that cross (touching does not count)."""

    segments = [list(zip(polyline, polyline[1:])) for polyline in polylines]
    count = 0
    for left, right in candidate_aabb_pairs(
        [_polyline_aabb(polyline) for polyline in polylines], bucket_size=64.0
    ):
        if any(
            segments_intersect(first, second, touching=False)
            for first in segments[left]
            for second in segments[right]
        ):
            count += 1
    return count


def label_fit(boxes: Sequence[Aabb], leaders: Sequence[Polyline]) -> LabelFit:
    """The fit of one record's labels of one kind."""

    return LabelFit(
        overlap_pairs=label_box_overlap_count(boxes),
        crossing_pairs=polyline_crossing_count(leaders),
    )


def label_fit_outcome(fits: Sequence[LabelFit]) -> LabelFitOutcome:
    """The one decision point: overlaps stop, crossings alone warn (LP-D42, LP-D46)."""

    if any(fit.overlap_pairs > 0 for fit in fits):
        return "stop"
    if any(fit.crossing_pairs > 0 for fit in fits):
        return "warn"
    return "fits"


@dataclass(frozen=True)
class LabelFitRecord:
    """A record whose labels of one kind do not fit."""

    record_id: str
    label_kind: LabelKind
    label_count: int
    font_size: float
    fit: LabelFit


@dataclass(frozen=True)
class LabelFitRow:
    """One lever: its mildest fitting value, or why there is none.

    ``values`` maps config paths to the values to set and ``canvas_size`` is
    the (width, height) of the verified figure; both are set only when
    ``status`` is ``found``.
    """

    kind: LabelFitRowKind
    status: LabelFitRowStatus
    values: Mapping[str, float] | None = None
    canvas_size: tuple[float, float] | None = None
    reason: str | None = None


@dataclass(frozen=True)
class LabelFitReport:
    """Failing records and the rows, in the dialog order.

    ``drawn`` is True when the figure was drawn anyway (crossing leaders only).
    """

    records: tuple[LabelFitRecord, ...]
    rows: tuple[LabelFitRow, ...]
    drawn: bool

    @property
    def overlap_pairs(self) -> int:
        return sum(record.fit.overlap_pairs for record in self.records)

    @property
    def crossing_pairs(self) -> int:
        return sum(record.fit.crossing_pairs for record in self.records)


__all__ = [
    "LabelFit",
    "LabelFitOutcome",
    "LabelFitRecord",
    "LabelFitReport",
    "LabelFitRow",
    "LabelFitRowKind",
    "LabelFitRowStatus",
    "LabelKind",
    "Polyline",
    "label_box_overlap_count",
    "label_fit",
    "label_fit_outcome",
    "polyline_crossing_count",
]
