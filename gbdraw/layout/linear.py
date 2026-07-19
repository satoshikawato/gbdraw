#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from dataclasses import dataclass, field
import math
from types import MappingProxyType
from typing import Iterable, Mapping, Sequence


@dataclass(frozen=True)
class VerticalBand:
    """One closed vertical interval in record-axis-local coordinates."""

    top_y: float
    bottom_y: float

    def __post_init__(self) -> None:
        top_y = float(self.top_y)
        bottom_y = float(self.bottom_y)
        if not math.isfinite(top_y) or not math.isfinite(bottom_y):
            raise ValueError("vertical band coordinates must be finite")
        if top_y > bottom_y:
            raise ValueError("vertical band top_y must not exceed bottom_y")
        object.__setattr__(self, "top_y", top_y)
        object.__setattr__(self, "bottom_y", bottom_y)

    @property
    def height(self) -> float:
        return self.bottom_y - self.top_y

    def translate(self, offset_y: float) -> VerticalBand:
        offset = float(offset_y)
        if not math.isfinite(offset):
            raise ValueError("vertical band translation must be finite")
        return VerticalBand(self.top_y + offset, self.bottom_y + offset)

    def expand(self, top: float = 0.0, bottom: float | None = None) -> VerticalBand:
        top_value = max(0.0, float(top))
        bottom_value = top_value if bottom is None else max(0.0, float(bottom))
        return VerticalBand(self.top_y - top_value, self.bottom_y + bottom_value)

    def union(self, *others: VerticalBand) -> VerticalBand:
        bands = (self, *others)
        return VerticalBand(
            min(band.top_y for band in bands),
            max(band.bottom_y for band in bands),
        )

    def intersects(self, other: VerticalBand, *, epsilon: float = 1e-9) -> bool:
        tolerance = max(0.0, float(epsilon))
        return (
            self.top_y < other.bottom_y - tolerance
            and self.bottom_y > other.top_y + tolerance
        )


def union_vertical_bands(
    bands: Iterable[VerticalBand],
    *,
    default: VerticalBand | None = None,
) -> VerticalBand:
    """Return the smallest band containing every input band."""

    collected = tuple(bands)
    if not collected:
        if default is None:
            raise ValueError("at least one vertical band is required")
        return default
    return VerticalBand(
        min(band.top_y for band in collected),
        max(band.bottom_y for band in collected),
    )


@dataclass(frozen=True)
class LinearFeatureLane:
    """Measured geometry for one rendered Linear feature lane."""

    strand_pool: str
    track_id: int
    top_y: float
    middle_y: float
    bottom_y: float
    band: VerticalBand

    @property
    def positions(self) -> tuple[float, float, float]:
        """Return the resolved feature path coordinates for this lane."""

        return self.top_y, self.middle_y, self.bottom_y


@dataclass(frozen=True)
class LinearFeatureLaneGeometry:
    """All feature lanes and their combined record-axis-local footprint."""

    lanes: tuple[LinearFeatureLane, ...]
    occupied_band: VerticalBand
    _lanes_by_identity: Mapping[tuple[str, int], LinearFeatureLane] = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        lanes_by_identity = {
            (lane.strand_pool, lane.track_id): lane
            for lane in self.lanes
        }
        if len(lanes_by_identity) != len(self.lanes):
            raise ValueError("Linear feature lane identities must be unique")
        object.__setattr__(
            self,
            "_lanes_by_identity",
            MappingProxyType(lanes_by_identity),
        )

    def lane_for(
        self,
        *,
        strand: str,
        track_id: int,
        separate_strands: bool,
    ) -> LinearFeatureLane:
        """Return the lane shared by measurement, labels, and SVG rendering."""

        strand_pool = str(strand) if separate_strands else "shared"
        resolved_track_id = int(track_id)
        identity = (strand_pool, resolved_track_id)
        try:
            return self._lanes_by_identity[identity]
        except KeyError:
            raise KeyError(
                "no measured Linear feature lane for "
                f"strand_pool={strand_pool!r}, track_id={resolved_track_id}"
            ) from None


def _legacy_middle_factors(strand: str, track_id: int, separate_strands: bool) -> list[float]:
    """Legacy factor calculation used for backward-compatible 'middle' layout."""
    initial_offset = 0.1
    feature_height = 0.5

    def calculate_track_offset(track_num: int) -> float:
        if separate_strands:
            return (
                track_num * -1.5 * initial_offset
                if strand == "positive"
                else (track_num - 1) * -1.5 * initial_offset
            )
        return track_num * -3 * initial_offset

    if not separate_strands:
        track_offset = calculate_track_offset(track_id)
        track = track_id
        return [
            -track - feature_height + track_offset,
            -track + track_offset,
            -track + feature_height + track_offset,
        ]

    # Keep middle layout vertically symmetric around axis in separated-strand mode.
    # Positive track 0 and negative track -1 should be mirrored with equal axis gaps.
    half_height = feature_height / 2.0
    track_step = 0.65
    axis_to_center = 0.35
    if strand == "negative":
        negative_index = max(0, abs(int(track_id)) - 1)
        middle = axis_to_center + (track_step * negative_index)
    else:
        positive_index = max(0, int(track_id))
        middle = -(axis_to_center + (track_step * positive_index))
    return [middle - half_height, middle, middle + half_height]


def _resolve_axis_gap_factor(
    *,
    separate_strands: bool,
    axis_gap_factor: float | None,
) -> float:
    """Resolve axis-to-feature edge gap in factor units (scaled by cds_height)."""
    if axis_gap_factor is not None:
        return max(0.0, float(axis_gap_factor))
    # Auto defaults tuned for above/below layouts to keep labels and features legible.
    if separate_strands:
        return 0.25
    return 0.30


def resolve_feature_axis_gap_linear(
    *,
    cds_height: float,
    separate_strands: bool,
    axis_gap: float | None,
) -> float:
    """Return the resolved feature-edge to axis gap in pixels."""
    if axis_gap is not None:
        return max(0.0, float(axis_gap))
    return max(0.0, float(cds_height)) * _resolve_axis_gap_factor(
        separate_strands=separate_strands,
        axis_gap_factor=None,
    )


def calculate_feature_position_factors_linear(
    strand: str,
    track_id: int,
    separate_strands: bool,
    track_layout: str = "middle",
    axis_gap_factor: float | None = None,
) -> list[float]:
    """
    Calculates feature position factors for linear track layouts.

    Returns:
        list[float]: Three position factors [top, middle, bottom] for feature placement.
    """
    layout = str(track_layout).strip().lower()
    if layout in {"spreadout", "above"}:
        layout = "above"
    elif layout in {"tuckin", "below"}:
        layout = "below"
    else:
        layout = "middle"

    if layout == "middle":
        return _legacy_middle_factors(strand, track_id, separate_strands)

    resolved_axis_gap = _resolve_axis_gap_factor(
        separate_strands=separate_strands,
        axis_gap_factor=axis_gap_factor,
    )

    if not separate_strands:
        track_index = max(0, abs(int(track_id)))
        track_spacing = 1.3
        half_height = 0.5
        middle = resolved_axis_gap + half_height + (track_spacing * track_index)
        if layout == "above":
            middle = -middle
        return [middle - half_height, middle, middle + half_height]

    # Separate-strands above/below modes keep positive strand above negative strand
    # while shifting both bands to the selected side of the axis.
    track_spacing = 0.65
    half_height = 0.25
    strand_band_gap = 0.6
    positive_index = max(0, int(track_id))
    negative_index = max(0, abs(int(track_id)) - 1)

    if layout == "above":
        if strand == "negative":
            middle = -(resolved_axis_gap + half_height + (track_spacing * negative_index))
        else:
            # Keep above-layout positive tracks displaced outward while aligning
            # band spacing with track spacing to avoid subtle track mismatch.
            middle = -(resolved_axis_gap + half_height + track_spacing + (track_spacing * positive_index))
    else:  # below
        if strand == "negative":
            middle = resolved_axis_gap + half_height + strand_band_gap + (track_spacing * negative_index)
        else:
            # Place displaced positive tracks from the default bottom band so
            # positive track i aligns with negative track -i for i >= 1.
            if positive_index == 0:
                middle = resolved_axis_gap + half_height + (track_spacing * positive_index)
            else:
                middle = (
                    resolved_axis_gap
                    + half_height
                    + strand_band_gap
                    + (track_spacing * (positive_index - 1))
                )

    return [middle - half_height, middle, middle + half_height]


def measure_linear_feature_lanes(
    feature_dict: Mapping[str, object],
    *,
    cds_height: float,
    separate_strands: bool,
    track_layout: str = "middle",
    axis_gap: float | None = None,
    stroke_width: float = 0.0,
) -> LinearFeatureLaneGeometry:
    """Measure feature lanes using the same factors consumed by the renderer."""

    height = max(0.0, float(cds_height))
    half_stroke = 0.5 * max(0.0, float(stroke_width))
    axis_gap_factor = (
        float(axis_gap) / height
        if axis_gap is not None and height > 0.0
        else None
    )
    lanes_by_identity: dict[tuple[str, int], LinearFeatureLane] = {}
    for feature in feature_dict.values():
        track_id = int(getattr(feature, "feature_track_id", 0))
        strand = str(getattr(feature, "strand", "undefined"))
        strand_pool = strand if separate_strands else "shared"
        identity = (strand_pool, track_id)
        if identity in lanes_by_identity:
            continue
        factors = calculate_feature_position_factors_linear(
            strand=strand,
            track_id=track_id,
            separate_strands=bool(separate_strands),
            track_layout=track_layout,
            axis_gap_factor=axis_gap_factor,
        )
        top_y = height * float(factors[0])
        middle_y = height * float(factors[1])
        bottom_y = height * float(factors[2])
        band = VerticalBand(top_y - half_stroke, bottom_y + half_stroke)
        lanes_by_identity[identity] = LinearFeatureLane(
            strand_pool=strand_pool,
            track_id=track_id,
            top_y=top_y,
            middle_y=middle_y,
            bottom_y=bottom_y,
            band=band,
        )

    lanes = tuple(
        sorted(
            lanes_by_identity.values(),
            key=lambda lane: (lane.band.top_y, lane.band.bottom_y, lane.strand_pool, lane.track_id),
        )
    )
    occupied_band = union_vertical_bands(
        (lane.band for lane in lanes),
        default=VerticalBand(0.0, 0.0),
    )
    return LinearFeatureLaneGeometry(lanes=lanes, occupied_band=occupied_band)


def measure_linear_label_band(
    labels: Sequence[Mapping[str, object]],
    *,
    leader_stroke_width: float = 0.0,
) -> VerticalBand | None:
    """Measure the vertical paint extent of prepared Linear labels and leaders."""

    # Imported lazily to keep the lane geometry module independent of label setup.
    from ..labels.linear import calculate_label_y_bounds

    if not labels:
        return None
    half_stroke = 0.5 * max(0.0, float(leader_stroke_width))
    bands: list[VerticalBand] = []
    for label in labels:
        top_y, bottom_y = calculate_label_y_bounds(label)
        leader_values = (
            label.get("leader_start_y"),
            label.get("leader_end_y"),
        )
        leader_y = [
            float(value)
            for value in leader_values
            if value is not None and math.isfinite(float(value))
        ]
        if leader_y:
            top_y = min(float(top_y), min(leader_y) - half_stroke)
            bottom_y = max(float(bottom_y), max(leader_y) + half_stroke)
        bands.append(VerticalBand(float(top_y), float(bottom_y)))
    return union_vertical_bands(bands)


__all__ = [
    "LinearFeatureLane",
    "LinearFeatureLaneGeometry",
    "VerticalBand",
    "calculate_feature_position_factors_linear",
    "measure_linear_feature_lanes",
    "measure_linear_label_band",
    "resolve_feature_axis_gap_linear",
    "union_vertical_bands",
]


