"""Surface-neutral serialization of resolved track-slot geometry.

Each record's entry also carries its ``recordKey`` and its drawn ``axis``
(``CircularAxisGeometry`` or ``LinearAxisGeometry``): where the axis lies and
which source base each drawn base shows, so the Web maps a pointer on the
drawing to a record position without measuring drawn shapes. The coordinates
are in the local frame of the record's drawn group (the innermost ``g`` that
carries ``data-gbdraw-record-index``, outside any record definition and other
than a drawn annotation's group), so a
Layout edit move of the record keeps them right.
"""

from __future__ import annotations

from collections.abc import MutableMapping
from dataclasses import dataclass
from typing import Literal, Mapping, Sequence, TypedDict

from Bio.SeqRecord import SeqRecord

from gbdraw.features.placement import FeaturePlacementSlot
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.record_transform import spans_may_cross_origin
from gbdraw.svg.circular_arcs import position_angle


class AxisSourceMap(TypedDict):
    """Which source base a drawn base shows.

    The drawn bases are numbered 0 to ``drawnLength - 1`` along the axis (the
    display index); display index ``i`` shows source base
    ``(anchorBase - 1 + step * i) mod sourceLength + 1`` (1-based, the input
    file's coordinates; ``step`` -1 on a reverse-complemented record).
    ``recordBase`` is the source base of the drawn record's first base (its
    feature coordinates count from it): record base ``k`` (0-based) is source
    base ``(recordBase - 1 + step * k) mod sourceLength + 1``.
    ``isCircular`` says whether a span may cross the record's origin (``spans_may_cross_origin``,
    the annotation resolver's rule): the record is drawn as circular (``RecordDisplayTransform.is_circular``,
    which the record planner sets from ``ResolvedRecordProvenance.is_circular``) and no region was
    applied to it.
    """

    anchorBase: int
    recordBase: int
    step: Literal[-1, 1]
    sourceLength: int
    drawnLength: int
    isCircular: bool


class CircularAxisGeometry(TypedDict):
    """A Circular record's axis: centred on the frame's origin.

    Display offset ``d`` (bp from the drawn start) lies at the angle
    ``startAngleDeg + degreesPerBp * d`` (degrees, SVG orientation). The band
    ``bandInnerPx``..``bandOuterPx`` is the radial extent of the drawn feature
    tracks.
    """

    kind: Literal["circular"]
    radiusPx: float
    bandInnerPx: float
    bandOuterPx: float
    startAngleDeg: float
    degreesPerBp: float
    source: AxisSourceMap


class LinearAxisRow(TypedDict):
    """One drawn row of a Linear record: display offsets ``displayStart``..``displayEnd``
    run from ``xStartPx`` to ``xEndPx`` on the axis at ``yPx``; the feature tracks
    lie between ``bandTopPx`` and ``bandBottomPx``."""

    xStartPx: float
    xEndPx: float
    yPx: float
    bandTopPx: float
    bandBottomPx: float
    displayStart: int
    displayEnd: int


class LinearAxisGeometry(TypedDict):
    """A Linear record's axis: one row, or one per Wrap Rows row, top to bottom."""

    kind: Literal["linear"]
    rows: list[LinearAxisRow]
    source: AxisSourceMap


class RadialExtent(TypedDict):
    """A drawn track's radial extent (inner and outer radius, px)."""

    innerPx: float
    outerPx: float


class CircularTrackSlotEntry(TypedDict):
    """One resolved track slot of a Circular record: its identity, renderer, side and radial size."""

    slotIndex: int
    slotId: str
    renderer: str
    side: str
    widthPx: float
    widthFactor: float | None
    radiusFactor: float
    innerGapPx: float
    outerGapPx: float
    compressed: bool
    source: Literal["resolved"]


class _CircularTrackSlotRecordFields(TypedDict):
    recordIndex: int
    recordId: str
    recordLabel: str
    recordKey: str
    axisRadiusPx: float
    axis: CircularAxisGeometry
    slots: list[CircularTrackSlotEntry]


class CircularTrackSlotRecord(_CircularTrackSlotRecordFields, total=False):
    """One drawn Circular record's track-slot geometry. ``featurePlacementTargets`` is present only when the
    drawing has a feature slot (``attach_feature_placement_targets`` writes it)."""

    featurePlacementTargets: list[dict[str, object]]


@dataclass(frozen=True)
class RecordAxisInput:
    """One drawn record as its axis is built from: the record, its identity in the
    Result (``recordKey``), and the transform it is drawn by (the plan's, or
    ``bind_record_transforms``' for a record no plan bound one to)."""

    record: SeqRecord
    record_key: str
    transform: RecordDisplayTransform

    @classmethod
    def of_records(
        cls,
        records: Sequence[SeqRecord],
        record_keys: Sequence[str],
        transforms: Sequence[RecordDisplayTransform],
    ) -> list[RecordAxisInput]:
        """One input per drawn record, in drawing order."""

        return [
            cls(record, str(record_key), transform)
            for record, record_key, transform in zip(records, record_keys, transforms, strict=True)
        ]

    def source(self) -> AxisSourceMap:
        """Which source base each drawn base of the record shows."""

        transform = self.transform
        return {
            "anchorBase": transform.display_index_to_source_base(0),
            "recordBase": transform.source_base,
            "step": transform.source_step,
            "sourceLength": int(transform.length),
            "drawnLength": len(self.record.seq),
            "isCircular": spans_may_cross_origin(self.record, transform),
        }


def circular_axis_geometry(
    drawn: RecordAxisInput,
    *,
    radius_px: float,
    feature_bands: Sequence[RadialExtent],
) -> CircularAxisGeometry:
    """The axis of one Circular record; the band spans ``feature_bands`` (the axis alone without one)."""

    source = drawn.source()
    drawn_length = source["drawnLength"]
    start_angle = position_angle(0, drawn_length)
    return {
        "kind": "circular",
        "radiusPx": float(radius_px),
        "bandInnerPx": min((band["innerPx"] for band in feature_bands), default=float(radius_px)),
        "bandOuterPx": max((band["outerPx"] for band in feature_bands), default=float(radius_px)),
        "startAngleDeg": start_angle,
        "degreesPerBp": position_angle(1, drawn_length) - start_angle,
        "source": source,
    }


def linear_axis_geometry(
    drawn: RecordAxisInput,
    *,
    rows: Sequence[tuple[float, float, float]],
    row_length_bp: int,
    px_per_bp: float,
) -> LinearAxisGeometry:
    """The axis of one Linear record from each row's ``(axis y, band top, band bottom)``.

    Row ``r`` shows the display offsets from ``r * row_length_bp``; every row
    starts at x 0 of the record's frame.
    """

    source = drawn.source()
    drawn_length = source["drawnLength"]
    payload: list[LinearAxisRow] = []
    for row_index, (axis_y, band_top, band_bottom) in enumerate(rows):
        start = row_index * row_length_bp
        end = min(drawn_length, start + row_length_bp)
        payload.append(
            {
                "xStartPx": 0.0,
                "xEndPx": float((end - start) * px_per_bp),
                "yPx": float(axis_y),
                "bandTopPx": float(band_top),
                "bandBottomPx": float(band_bottom),
                "displayStart": int(start),
                "displayEnd": int(end),
            }
        )
    return {"kind": "linear", "rows": payload, "source": source}


# The attribute a drawing carries its resolved track-slot geometry in.
TRACK_SLOT_GEOMETRY_ATTR = "_gbdraw_track_slot_geometry"


def attach_track_slot_geometry(canvas: object, geometry: Mapping[str, object]) -> None:
    """Record a drawing's resolved track-slot geometry on it, where ``collect_track_slot_geometry_records`` reads it."""

    setattr(canvas, TRACK_SLOT_GEOMETRY_ATTR, geometry)


def attach_feature_placement_targets(
    records: Sequence[CircularTrackSlotRecord] | Sequence[MutableMapping[str, object]],
    slot: FeaturePlacementSlot | None,
) -> None:
    """Give each record's geometry the placement targets its feature slot supports (no slot: none written)."""

    if slot is None:
        return
    targets = slot.supported_targets()
    for record_geometry in records:
        record_geometry["featurePlacementTargets"] = targets


def collect_track_slot_geometry_records(
    canvas: object,
    *,
    result_index: int,
    result_name: str,
) -> list[dict[str, object]]:
    """Copy track-slot geometry from a canvas and add result identity."""

    geometry = getattr(canvas, TRACK_SLOT_GEOMETRY_ATTR, None)
    if not isinstance(geometry, Mapping):
        return []
    records = geometry.get("records")
    if not isinstance(records, Sequence) or isinstance(records, (str, bytes)):
        return []
    copied_records: list[dict[str, object]] = []
    for record in records:
        if not isinstance(record, Mapping):
            continue
        copied_record = dict(record)
        copied_record["resultIndex"] = int(result_index)
        copied_record["resultName"] = str(result_name)
        slots = copied_record.get("slots")
        if isinstance(slots, Sequence) and not isinstance(slots, (str, bytes)):
            copied_record["slots"] = [
                dict(slot) for slot in slots if isinstance(slot, Mapping)
            ]
        else:
            copied_record["slots"] = []
        copied_records.append(copied_record)
    return copied_records


class TrackSlotGeometry(TypedDict):
    """A render's resolved track-slot geometry: one entry per drawn record of each Result."""

    schema: Literal[1]
    mode: Literal["circular", "linear"]
    source: Literal["resolved"]
    records: list[dict[str, object]]


def resolved_track_slot_geometry(
    *,
    mode: Literal["circular", "linear"],
    records: Sequence[Mapping[str, object]],
) -> TrackSlotGeometry | None:
    """The geometry of the collected records; ``None`` when no Result has any."""

    if not records:
        return None
    return {
        "schema": 1,
        "mode": mode,
        "source": "resolved",
        "records": [dict(record) for record in records],
    }


__all__ = [
    "AxisSourceMap",
    "CircularAxisGeometry",
    "LinearAxisGeometry",
    "LinearAxisRow",
    "RadialExtent",
    "RecordAxisInput",
    "TrackSlotGeometry",
    "circular_axis_geometry",
    "linear_axis_geometry",
    "collect_track_slot_geometry_records",
    "resolved_track_slot_geometry",
]
