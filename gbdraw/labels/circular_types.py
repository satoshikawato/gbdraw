"""Small immutable contracts shared by circular label placement strategies."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal, Mapping, TypeAlias

from ..exceptions import ValidationError
from ..layout.text_geometry import AabbTuple, Point

CircularLabelPlacement: TypeAlias = Literal["horizontal", "radial"]


def normalize_circular_label_placement(value: object) -> CircularLabelPlacement:
    placement = str(value or "horizontal").strip().lower()
    if placement not in {"horizontal", "radial"}:
        raise ValidationError(
            "circular label placement must be 'horizontal' or 'radial' "
            f"(received {value!r})"
        )
    return placement  # type: ignore[return-value]


@dataclass(frozen=True)
class CircularLabelCandidate:
    """Placement-independent semantic data and measured single-line text metrics."""

    stable_id: str
    input_order: int
    text: str
    font_family: str
    font_size: float
    width_px: float
    height_px: float
    feature_type: str
    strand: str
    track_id: int
    directional: bool
    segment_start_bp: float
    segment_end_bp: float
    segment_middle_bp: float
    segment_span_bp: float
    feature_coordinates: tuple[tuple[float, float], ...]


@dataclass(frozen=True)
class CircularPlacedLabel:
    """Authoritative geometry for one external circular feature label."""

    preferred_angle_deg: float
    placed_angle_deg: float
    text_x: float
    text_y: float
    rotation_deg: float
    text_anchor: Literal["start", "middle", "end"]
    dominant_baseline: str
    width_px: float
    height_px: float
    corners: tuple[Point, Point, Point, Point]
    aabb: AabbTuple
    feature_anchor: Point
    elbow: Point
    leader_start: Point
    is_inner: bool
    metadata: Mapping[str, Any]


@dataclass(frozen=True)
class CircularLabelLayoutResult:
    """Complete placement result consumed by rendering and canvas orchestration."""

    labels: tuple[Mapping[str, Any], ...]
    content_bounds: AabbTuple | None
    required_radius_growth_px: float = 0.0
    text_collision_count: int = 0
    leader_text_collision_count: int = 0
    leader_crossing_count: int = 0
    order_violation_count: int = 0

    @property
    def collision_free(self) -> bool:
        return (
            self.text_collision_count == 0
            and self.leader_text_collision_count == 0
            and self.leader_crossing_count == 0
            and self.order_violation_count == 0
        )


__all__ = [
    "CircularLabelCandidate",
    "CircularLabelLayoutResult",
    "CircularLabelPlacement",
    "CircularPlacedLabel",
    "normalize_circular_label_placement",
]
