#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal, Mapping, Sequence, TypeAlias

from ._immutable import freeze_layout_mapping

if TYPE_CHECKING:
    from ..config.models import CircularRenderProfile


LAYOUT_EPSILON = 1e-6


@dataclass(frozen=True)
class RadialBand:
    inner_px: float
    outer_px: float

    def __post_init__(self) -> None:
        inner, outer = sorted((max(0.0, float(self.inner_px)), max(0.0, float(self.outer_px))))
        object.__setattr__(self, "inner_px", inner)
        object.__setattr__(self, "outer_px", outer)

    @property
    def width_px(self) -> float:
        return max(0.0, float(self.outer_px) - float(self.inner_px))

    @property
    def center_px(self) -> float:
        return (float(self.inner_px) + float(self.outer_px)) / 2.0

    def expanded(self, spacing_px: float) -> "RadialBand":
        spacing = max(0.0, float(spacing_px))
        return RadialBand(max(0.0, self.inner_px - spacing), self.outer_px + spacing)

    def expanded_sides(self, inner_gap_px: float, outer_gap_px: float) -> "RadialBand":
        inner_gap = max(0.0, float(inner_gap_px))
        outer_gap = max(0.0, float(outer_gap_px))
        return RadialBand(max(0.0, self.inner_px - inner_gap), self.outer_px + outer_gap)


@dataclass(frozen=True)
class CircularFeatureLane:
    track_id: int
    strand_group: str
    inner_px: float
    center_px: float
    outer_px: float

    @property
    def band_px(self) -> RadialBand:
        return RadialBand(self.inner_px, self.outer_px)


@dataclass(frozen=True)
class CircularFeatureStackMetrics:
    lane_count: int
    lane_width_px: float
    lane_spacing_px: float
    band_width_px: float
    center_radius_px: float
    inner_radius_px: float
    outer_radius_px: float
    lane_centers_by_track_id: Mapping[int, float]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "lane_centers_by_track_id",
            freeze_layout_mapping(self.lane_centers_by_track_id),
        )


@dataclass(frozen=True)
class CircularFeatureLayout:
    anchor_radius_px: float
    width_px: float
    lanes_by_track_id: Mapping[int, CircularFeatureLane]
    primary_band_px: RadialBand
    all_band_px: RadialBand
    stack_metrics: CircularFeatureStackMetrics | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "lanes_by_track_id",
            freeze_layout_mapping(self.lanes_by_track_id),
        )

    def lane_for_track_id(self, track_id: int) -> CircularFeatureLane:
        if int(track_id) in self.lanes_by_track_id:
            return self.lanes_by_track_id[int(track_id)]
        if 0 in self.lanes_by_track_id:
            return self.lanes_by_track_id[0]
        return next(iter(self.lanes_by_track_id.values()))


@dataclass(frozen=True)
class CircularAxisLayout:
    radius_px: float
    stroke_width_px: float


@dataclass(frozen=True)
class CircularTickLayout:
    anchor_radius_px: float
    tick_band_px: RadialBand
    label_band_px: RadialBand | None
    reserved_band_px: RadialBand
    track_preset: str
    label_side: str
    tick_side: str
    tick_length_px: float | None


@dataclass(frozen=True)
class CircularFeatureLabelLayout:
    side: str
    reserved_band_px: RadialBand
    labels: tuple[Mapping[str, Any], ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "labels",
            tuple(freeze_layout_mapping(label) for label in self.labels),
        )


CircularSlotPayload: TypeAlias = CircularFeatureLayout | CircularTickLayout | CircularFeatureLabelLayout
CircularFeatureLaneDirection: TypeAlias = Literal["inside", "outside", "split"]


@dataclass(frozen=True)
class CircularResolvedSlot:
    slot_index: int
    id: str
    renderer: str
    side: str
    z: int
    anchor_radius_px: float | None
    anchor_offset_px: float | None
    requested_width_px: float | None
    resolved_width_px: float | None
    packing_band_px: RadialBand | None
    draw_band_px: RadialBand | None
    reserved_band_px: RadialBand | None
    inner_gap_px: float
    outer_gap_px: float
    params: Mapping[str, Any]
    payload: CircularSlotPayload | None = None
    explicit_anchor: bool = False
    explicit_width: bool = False
    compressed: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "params", freeze_layout_mapping(self.params))

    @property
    def center_radius_px(self) -> float:
        if self.draw_band_px is not None:
            return float(self.draw_band_px.center_px)
        return float(self.anchor_radius_px or 0.0)

    @property
    def draw_inner_radius_px(self) -> float:
        return float(self.draw_band_px.inner_px) if self.draw_band_px is not None else 0.0

    @property
    def draw_outer_radius_px(self) -> float:
        return float(self.draw_band_px.outer_px) if self.draw_band_px is not None else 0.0

    @property
    def reserved_inner_radius_px(self) -> float:
        return float(self.reserved_band_px.inner_px) if self.reserved_band_px is not None else 0.0

    @property
    def reserved_outer_radius_px(self) -> float:
        return float(self.reserved_band_px.outer_px) if self.reserved_band_px is not None else 0.0

    @property
    def draw_width_px(self) -> float:
        return float(self.draw_band_px.width_px) if self.draw_band_px is not None else 0.0

    @property
    def reserved_width_px(self) -> float:
        return float(self.reserved_band_px.width_px) if self.reserved_band_px is not None else 0.0


@dataclass(frozen=True)
class CircularRadialLayout:
    axis: CircularAxisLayout
    slots: tuple[CircularResolvedSlot, ...]
    definition_reserved_band_px: RadialBand | None
    outer_content_radius_px: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "slots", tuple(self.slots))

    @property
    def features(self) -> CircularFeatureLayout | None:
        for slot in self.slots:
            if slot.renderer == "features" and isinstance(slot.payload, CircularFeatureLayout):
                return slot.payload
        return None

    @property
    def ticks(self) -> CircularTickLayout | None:
        for slot in self.slots:
            if slot.renderer == "ticks" and isinstance(slot.payload, CircularTickLayout):
                return slot.payload
        return None

    @property
    def tracks(self) -> tuple[CircularResolvedSlot, ...]:
        return tuple(slot for slot in self.slots if slot.renderer not in {"features", "ticks", "feature_labels"})


@dataclass(frozen=True)
class CircularRecordRenderContext:
    """Resolved layout state required to render one Circular record."""

    profile: CircularRenderProfile
    track_preset: str
    feature_lane_direction: CircularFeatureLaneDirection
    radial_layout: CircularRadialLayout

    @property
    def feature_layout(self) -> CircularFeatureLayout | None:
        return self.radial_layout.features


def band_union(bands: Sequence[RadialBand]) -> RadialBand | None:
    non_empty = [band for band in bands if band.width_px > LAYOUT_EPSILON]
    if not non_empty:
        return None
    return RadialBand(
        min(float(band.inner_px) for band in non_empty),
        max(float(band.outer_px) for band in non_empty),
    )


def bands_overlap(a: RadialBand, b: RadialBand) -> bool:
    return (
        float(a.inner_px) < (float(b.outer_px) - LAYOUT_EPSILON)
        and float(a.outer_px) > (float(b.inner_px) + LAYOUT_EPSILON)
    )


def feature_band_bounds_px(feature_layout: CircularFeatureLayout | None) -> tuple[float, float] | None:
    if feature_layout is None:
        return None
    return float(feature_layout.all_band_px.inner_px), float(feature_layout.all_band_px.outer_px)


def feature_radii_for_object(
    feature_object: Any,
    feature_layout: CircularFeatureLayout,
) -> tuple[float, float, float]:
    lane = feature_layout.lane_for_track_id(int(getattr(feature_object, "feature_track_id", 0)))
    return float(lane.inner_px), float(lane.center_px), float(lane.outer_px)


def feature_radii_for_coordinate(
    feature_object: Any,
    coordinate: Any,
    feature_layout: CircularFeatureLayout,
) -> tuple[float, float, float]:
    del coordinate
    return feature_radii_for_object(feature_object, feature_layout)


def feature_radius_intervals(
    feature_dict: Mapping[str, Any],
    total_length: int,
    feature_layout: CircularFeatureLayout,
) -> list[tuple[float, float, float, float]]:
    if total_length <= 0:
        return []

    intervals: list[tuple[float, float, float, float]] = []
    total_length_float = float(total_length)
    for feature_object in feature_dict.values():
        lane = feature_layout.lane_for_track_id(int(getattr(feature_object, "feature_track_id", 0)))
        feature_location_list = getattr(feature_object, "location", [])
        list_of_coordinates = getattr(feature_object, "coordinates", [])

        for coordinate_idx, coordinate in enumerate(list_of_coordinates):
            if (
                coordinate_idx < len(feature_location_list)
                and getattr(feature_location_list[coordinate_idx], "kind", None)
                == "line"
            ):
                continue

            coordinate_start_raw = int(coordinate.start)
            coordinate_end_raw = int(coordinate.end)
            if coordinate_start_raw == coordinate_end_raw:
                continue

            coordinate_start = float(coordinate_start_raw % total_length)
            coordinate_end = float(coordinate_end_raw % total_length)
            if coordinate_end_raw > coordinate_start_raw and coordinate_end_raw <= total_length:
                coordinate_end = float(coordinate_end_raw)

            if coordinate_start < coordinate_end:
                intervals.append((coordinate_start, coordinate_end, lane.inner_px, lane.outer_px))
            else:
                intervals.append((coordinate_start, total_length_float, lane.inner_px, lane.outer_px))
                intervals.append((0.0, coordinate_end, lane.inner_px, lane.outer_px))
    return intervals


def calculate_feature_position_factors_circular(
    total_length: int,
    strand: str,
    track_ratio: float,
    cds_ratio: float,
    offset: float,
    track_type: str = "tuckin",
    strandedness: bool = True,
    track_id: int = 0,
) -> list[float]:
    """
    Calculates position factors for a feature based on its strand orientation on a circular canvas.

    The factors determine the inner, middle, and outer radii of the feature arc relative to
    the base radius. When resolve_overlaps is enabled with strandedness=False, track_id is
    used to offset features to avoid visual overlap.

    Args:
        total_length: Total length of the genome (not used in calculation but kept for API consistency)
        strand: "positive", "negative", or "undefined"
        track_ratio: Base track ratio from config
        cds_ratio: Calculated CDS ratio for feature height
        offset: Base offset value
        track_type: "tuckin", "middle", or "spreadout"
        strandedness: Whether to separate strands
        track_id: Track number for overlap resolution. Only effective when strandedness=False.
                  - spreadout: displaced outward
                  - middle: positive/undefined displaced outward, negative displaced inward
                  - tuckin: displaced inward

    Returns:
        List of three floats [inner_factor, middle_factor, outer_factor] used to calculate
        the actual radii by multiplying with the base radius.
    """
    BASE: float = 1.0
    cds_ratio = float(cds_ratio)
    offset = float(offset)

    lane_width = cds_ratio
    lane_spacing = 0.01
    lane_step = lane_width + lane_spacing

    def factors_from_center(center: float) -> list[float]:
        return [
            center - (0.5 * lane_width),
            center,
            center + (0.5 * lane_width),
        ]

    if strandedness is True:
        if track_type == "middle":
            center_positive = BASE + (0.5 * lane_spacing) + (0.5 * lane_width)
            center_negative = BASE - (0.5 * lane_spacing) - (0.5 * lane_width)
        elif track_type == "spreadout":
            center_positive = BASE + lane_spacing + (0.5 * lane_width)
            center_negative = center_positive + lane_step
        elif track_type == "tuckin":
            center_positive = BASE - lane_spacing - (0.5 * lane_width)
            center_negative = center_positive - lane_step
        else:
            center_positive = BASE + (0.5 * lane_spacing) + (0.5 * lane_width)
            center_negative = BASE - (0.5 * lane_spacing) - (0.5 * lane_width)

        if strand == "positive":
            factors = factors_from_center(center_positive + offset)
        else:
            factors = factors_from_center(center_negative - offset)

        return factors

    else:
        if track_type == "middle":
            if track_id != 0:
                if int(track_id) < 0 or strand == "negative":
                    center = BASE - (abs(track_id) * lane_step)
                else:
                    center = BASE + (abs(track_id) * lane_step)
            else:
                center = BASE
        elif track_type == "spreadout":
            center = BASE + lane_spacing + (0.5 * lane_width) + (abs(track_id) * lane_step)
        elif track_type == "tuckin":
            center = BASE - lane_spacing - (0.5 * lane_width) - (abs(track_id) * lane_step)
        else:
            center = BASE

        return factors_from_center(center)


__all__ = [
    "CircularAxisLayout",
    "CircularFeatureLabelLayout",
    "CircularFeatureLaneDirection",
    "CircularFeatureLane",
    "CircularFeatureLayout",
    "CircularFeatureStackMetrics",
    "CircularRadialLayout",
    "CircularRecordRenderContext",
    "CircularResolvedSlot",
    "CircularSlotPayload",
    "CircularTickLayout",
    "RadialBand",
    "band_union",
    "bands_overlap",
    "calculate_feature_position_factors_circular",
    "feature_band_bounds_px",
    "feature_radii_for_coordinate",
    "feature_radii_for_object",
    "feature_radius_intervals",
]
