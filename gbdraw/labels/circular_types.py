"""Small immutable contracts shared by circular label placement strategies."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, TypeAlias, TypedDict, cast

from ..exceptions import ValidationError
from ..layout.text_geometry import AabbTuple, Point
from .label_style import LabelFontStyle

CircularLabelPlacement: TypeAlias = Literal["horizontal", "radial"]
# The margins of a horizontal label's box in the fit count, in px: the box is this much wider and taller than the
# label (half on each side). The capacity check reads the same box.
LABEL_FIT_X_MARGIN_PX = 1.0
LABEL_FIT_Y_MARGIN_PX = 0.1


def normalize_circular_label_placement(value: object) -> CircularLabelPlacement:
    if not isinstance(value, str):
        raise ValidationError(
            "circular label placement must be 'horizontal' or 'radial' "
            f"(received {value!r})"
        )
    placement = value.strip().lower()
    if placement not in {"horizontal", "radial"}:
        raise ValidationError(
            "circular label placement must be 'horizontal' or 'radial' "
            f"(received {value!r})"
        )
    return cast(CircularLabelPlacement, placement)


@dataclass(frozen=True)
class CircularLabelCandidate:
    """Placement-independent semantic data and measured single-line text metrics."""

    stable_id: str
    input_order: int
    text: str
    font_family: str
    font_size: float
    font_style: LabelFontStyle
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
    source_anchor_span: tuple[int, int] | None = None
    source_part_index: int | None = None
    source_middle_bp: float | None = None


# The measured candidates of one label font size ("candidates"), shared by the
# placements of one record.
CircularLabelCandidateCache: TypeAlias = dict[str, tuple[CircularLabelCandidate, ...]]


class _CircularLabelRequired(TypedDict):
    """Keys every Circular label has from ``prepare_label_list`` on."""

    track_id: int
    label_text: str
    feature_id: str | None
    stable_id: str
    input_order: int
    preferred_middle: float
    preferred_angle_deg: float
    middle: float
    start: float
    end: float
    middle_x: float
    middle_y: float
    feature_middle_x: float
    feature_middle_y: float
    feature_anchor_x: float
    feature_anchor_y: float
    feature_inner_radius_px: float
    feature_center_radius_px: float
    feature_outer_radius_px: float
    track_preset: str
    width_px: float
    height_px: float
    strand: str
    is_embedded: bool
    font_size: float
    font_family: str


class CircularLabel(_CircularLabelRequired, total=False):
    """One Circular feature label as placement, fit check, repair and render share it.

    Embedded and external labels, horizontal and radial, travel in one list and
    are told apart by ``is_embedded``, ``is_inner`` and ``placement``, so the
    keys below are optional for every kind:

    - side of an external label: ``is_inner``;
    - feature clearance of a displaced lane: ``min_outer_start_radius_px``,
      ``max_inner_start_radius_px``, ``max_inner_middle_radius_px``,
      ``min_outer_middle_radius_px``;
    - placed text: ``start_x``, ``start_y``, ``angle_unwrapped``,
      ``target_angle_unwrapped``;
    - placement stage only: ``_spacing_px``, ``_enforce_spacing_px``,
      ``_fixed_radius_px``, ``_target_angle_unwrapped_fixed``,
      ``_angle_unwrapped_fixed``, ``_opt_original_idx``;
    - leader: ``leader_start_x``, ``leader_start_y``, ``leader_contact`` (a
      contact as fractions of the ink box);
    - radial placement: ``placement``, ``placed_angle_deg``, ``text_x``,
      ``text_y``, ``rotation_deg``, ``text_anchor``, ``dominant_baseline``,
      ``oriented_corners``, ``aabb_local``;
    - radial layout summary: ``required_radius_growth_px``,
      ``layout_collision_counts``;
    - the text's style: ``font_style`` (``FeatureLabelStyle.font_style``; absent for a list made before the
      style existed, which draws upright, ``normal``).
    """

    font_style: LabelFontStyle

    is_inner: bool
    min_outer_start_radius_px: float
    max_inner_start_radius_px: float
    max_inner_middle_radius_px: float
    min_outer_middle_radius_px: float
    start_x: float
    start_y: float
    angle_unwrapped: float
    target_angle_unwrapped: float
    _spacing_px: float
    _enforce_spacing_px: bool
    _fixed_radius_px: float
    _target_angle_unwrapped_fixed: float
    _angle_unwrapped_fixed: float
    _opt_original_idx: int
    leader_start_x: float
    leader_start_y: float
    leader_contact: tuple[float, float]
    placement: Literal["radial"]
    placed_angle_deg: float
    text_x: float
    text_y: float
    rotation_deg: float
    text_anchor: str
    dominant_baseline: str
    oriented_corners: tuple[Point, ...]
    aabb_local: AabbTuple
    required_radius_growth_px: float
    layout_collision_counts: tuple[int, int, int, int]


@dataclass(frozen=True)
class CircularLabelLayoutResult:
    """Complete placement result consumed by rendering and canvas orchestration."""

    labels: tuple[CircularLabel, ...]
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
    "CircularLabel",
    "CircularLabelCandidate",
    "CircularLabelCandidateCache",
    "CircularLabelLayoutResult",
    "CircularLabelPlacement",
    "normalize_circular_label_placement",
]
