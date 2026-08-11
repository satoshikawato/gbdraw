"""Pure, mode-neutral composition geometry for complete diagrams.

The planner in this module operates only on authoritative axis-aligned bounds.
It does not inspect rendered SVG, configuration files, or mode-specific feature
geometry.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from itertools import product

from .spatial import Aabb, union_aabbs


@dataclass(frozen=True, slots=True)
class CompositionSpacing:
    """Shared spacing policy for final diagram composition."""

    edge_padding_px: float = 16.0
    dock_gap_px: float = 24.0
    title_gap_px: float = 20.0
    stack_gap_px: float = 20.0
    overlay_clearance_px: float = 8.0

    def __post_init__(self) -> None:
        for field_name in (
            "edge_padding_px",
            "dock_gap_px",
            "title_gap_px",
            "stack_gap_px",
            "overlay_clearance_px",
        ):
            value = float(getattr(self, field_name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{field_name} must be a finite non-negative value")
            object.__setattr__(self, field_name, value)


DEFAULT_COMPOSITION_SPACING = CompositionSpacing()

# The Web editor replays the same pure composition algorithm after legend or
# title edits.  These ordered metric names are the policy used by the Python
# planner and are emitted with every schema-1 composition payload so the Web
# interpreter does not carry a second, implicit ordering policy.
OVERLAY_CANDIDATE_SCORE_ORDER = (
    "totalAnchorDistance",
    "xAnchorDistance",
    "yAnchorDistance",
    "nearEdgeX",
    "nearEdgeY",
)
OVERLAY_CANVAS_GROWTH_CANDIDATE_ORDER = ("horizontal", "vertical")
OVERLAY_CANVAS_GROWTH_SCORE_ORDER = (
    "addedArea",
    "addedExtent",
    "candidateOrder",
)
OVERLAY_QUADRANT_BOUNDARY_RATIO = 0.5


class LegendPlacement(str, Enum):
    """Closed set of supported legend placement semantics."""

    LEFT = "left"
    RIGHT = "right"
    TOP = "top"
    BOTTOM = "bottom"
    UPPER_LEFT = "upper_left"
    UPPER_RIGHT = "upper_right"
    LOWER_LEFT = "lower_left"
    LOWER_RIGHT = "lower_right"
    NONE = "none"

    @property
    def is_dock(self) -> bool:
        return self in {
            LegendPlacement.LEFT,
            LegendPlacement.RIGHT,
            LegendPlacement.TOP,
            LegendPlacement.BOTTOM,
        }

    @property
    def is_overlay(self) -> bool:
        return self in {
            LegendPlacement.UPPER_LEFT,
            LegendPlacement.UPPER_RIGHT,
            LegendPlacement.LOWER_LEFT,
            LegendPlacement.LOWER_RIGHT,
        }


class TitlePlacement(str, Enum):
    """Closed set of supported plot-title placement semantics."""

    TOP = "top"
    BOTTOM = "bottom"
    CENTER = "center"
    NONE = "none"


class OverlayResolution(str, Enum):
    """How an overlay legend's automatic placement was resolved."""

    ANCHORED = "anchored"
    SHIFTED = "shifted"
    CANVAS_GROWTH = "canvas_growth"


@dataclass(frozen=True, slots=True)
class CompositionItem:
    """One item with authoritative bounds in its own local coordinates."""

    role: str
    local_bounds: Aabb

    def __post_init__(self) -> None:
        if not isinstance(self.role, str) or not self.role.strip():
            raise ValueError("composition item role must be a non-empty string")
        if not isinstance(self.local_bounds, Aabb):
            raise TypeError("composition item local_bounds must be an Aabb")

    @property
    def is_empty(self) -> bool:
        """Whether this item has no two-dimensional painted extent."""
        return self.local_bounds.width <= 0.0 or self.local_bounds.height <= 0.0


def _coerce_legend_placement(value: LegendPlacement | str) -> LegendPlacement:
    if isinstance(value, LegendPlacement):
        return value
    try:
        return LegendPlacement(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"unknown legend placement: {value!r}") from exc


def _coerce_title_placement(value: TitlePlacement | str) -> TitlePlacement:
    if isinstance(value, TitlePlacement):
        return value
    try:
        return TitlePlacement(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"unknown title placement: {value!r}") from exc


@dataclass(frozen=True, slots=True)
class CompositionRequest:
    """Complete immutable input to the composition planner.

    Overlay obstacles use the same plot-space coordinate system as the primary
    item. Their tuple order is retained so conflict indices can be mapped back
    to mode-owned geometry without putting that geometry in the shared planner.
    """

    primary: CompositionItem
    legend: CompositionItem | None = None
    title: CompositionItem | None = None
    legend_placement: LegendPlacement | str = LegendPlacement.NONE
    title_placement: TitlePlacement | str = TitlePlacement.NONE
    overlay_obstacles: tuple[Aabb, ...] = ()
    spacing: CompositionSpacing = DEFAULT_COMPOSITION_SPACING

    def __post_init__(self) -> None:
        if not isinstance(self.primary, CompositionItem):
            raise ValueError("composition request requires a primary item")
        if self.primary.is_empty:
            raise ValueError("primary bounds must have positive width and height")
        if self.legend is not None and not isinstance(self.legend, CompositionItem):
            raise TypeError("legend must be a CompositionItem or None")
        if self.title is not None and not isinstance(self.title, CompositionItem):
            raise TypeError("title must be a CompositionItem or None")
        if not isinstance(self.spacing, CompositionSpacing):
            raise TypeError("spacing must be CompositionSpacing")

        object.__setattr__(
            self,
            "legend_placement",
            _coerce_legend_placement(self.legend_placement),
        )
        object.__setattr__(
            self,
            "title_placement",
            _coerce_title_placement(self.title_placement),
        )

        obstacles = tuple(self.overlay_obstacles)
        if any(not isinstance(obstacle, Aabb) for obstacle in obstacles):
            raise TypeError("overlay obstacles must be Aabb instances")
        object.__setattr__(self, "overlay_obstacles", obstacles)

        items = tuple(
            item
            for item in (self.primary, self.legend, self.title)
            if item is not None
        )
        roles = tuple(item.role for item in items)
        if len(set(roles)) != len(roles):
            raise ValueError("composition item roles must be unique")


@dataclass(frozen=True, slots=True)
class CompositionPlacement:
    """Final automatic translation and bounds for one composition item."""

    role: str
    translation: tuple[float, float]
    final_bounds: Aabb

    def __post_init__(self) -> None:
        if not isinstance(self.role, str) or not self.role.strip():
            raise ValueError("composition placement role must be a non-empty string")
        if len(self.translation) != 2:
            raise ValueError("composition translation must contain two values")
        translation = (float(self.translation[0]), float(self.translation[1]))
        if not all(math.isfinite(value) for value in translation):
            raise ValueError("composition translation must be finite")
        if not isinstance(self.final_bounds, Aabb):
            raise TypeError("composition final_bounds must be an Aabb")
        object.__setattr__(self, "translation", translation)

    @property
    def dx(self) -> float:
        return self.translation[0]

    @property
    def dy(self) -> float:
        return self.translation[1]


@dataclass(frozen=True, slots=True)
class CompositionPlan:
    """Resolved final canvas and automatic item placements."""

    canvas_bounds: Aabb
    view_box: tuple[float, float, float, float]
    primary_bounds: Aabb
    placements: tuple[CompositionPlacement, ...]
    spacing: CompositionSpacing
    overlay_obstacles: tuple[Aabb, ...] = ()
    overlay_conflict_indices: tuple[int, ...] = ()
    overlay_resolution: OverlayResolution | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.canvas_bounds, Aabb):
            raise TypeError("canvas_bounds must be an Aabb")
        if not isinstance(self.primary_bounds, Aabb):
            raise TypeError("primary_bounds must be an Aabb")
        view_box = tuple(float(value) for value in self.view_box)
        if len(view_box) != 4 or not all(math.isfinite(value) for value in view_box):
            raise ValueError("view_box must contain four finite values")
        object.__setattr__(self, "view_box", view_box)
        object.__setattr__(self, "placements", tuple(self.placements))
        obstacles = tuple(self.overlay_obstacles)
        if any(not isinstance(obstacle, Aabb) for obstacle in obstacles):
            raise TypeError("overlay_obstacles must be Aabb instances")
        object.__setattr__(self, "overlay_obstacles", obstacles)
        object.__setattr__(
            self,
            "overlay_conflict_indices",
            tuple(int(index) for index in self.overlay_conflict_indices),
        )

    @property
    def width(self) -> float:
        return self.canvas_bounds.width

    @property
    def height(self) -> float:
        return self.canvas_bounds.height

    def placement_for(self, role: str) -> CompositionPlacement | None:
        """Return the placement for ``role``, if that role is present."""
        return next(
            (placement for placement in self.placements if placement.role == role),
            None,
        )


@dataclass(frozen=True, slots=True)
class _WorkingPlacement:
    item: CompositionItem
    dx: float
    dy: float

    @property
    def bounds(self) -> Aabb:
        return self.item.local_bounds.translated(self.dx, self.dy)


def _align_min(
    item: CompositionItem,
    *,
    min_x: float,
    min_y: float,
) -> _WorkingPlacement:
    return _WorkingPlacement(
        item=item,
        dx=float(min_x) - item.local_bounds.min_x,
        dy=float(min_y) - item.local_bounds.min_y,
    )


def _align_center(
    item: CompositionItem,
    *,
    center_x: float,
    center_y: float,
) -> _WorkingPlacement:
    bounds = item.local_bounds
    return _align_min(
        item,
        min_x=float(center_x) - (0.5 * bounds.width),
        min_y=float(center_y) - (0.5 * bounds.height),
    )


def _place_docked_legend(
    primary: Aabb,
    legend: CompositionItem,
    placement: LegendPlacement,
    spacing: CompositionSpacing,
) -> _WorkingPlacement:
    gap = spacing.dock_gap_px
    if placement is LegendPlacement.LEFT:
        return _align_min(
            legend,
            min_x=primary.min_x - gap - legend.local_bounds.width,
            min_y=0.5 * (primary.min_y + primary.max_y - legend.local_bounds.height),
        )
    if placement is LegendPlacement.RIGHT:
        return _align_min(
            legend,
            min_x=primary.max_x + gap,
            min_y=0.5 * (primary.min_y + primary.max_y - legend.local_bounds.height),
        )
    if placement is LegendPlacement.TOP:
        return _align_min(
            legend,
            min_x=0.5 * (primary.min_x + primary.max_x - legend.local_bounds.width),
            min_y=primary.min_y - gap - legend.local_bounds.height,
        )
    if placement is LegendPlacement.BOTTOM:
        return _align_min(
            legend,
            min_x=0.5 * (primary.min_x + primary.max_x - legend.local_bounds.width),
            min_y=primary.max_y + gap,
        )
    raise ValueError(f"legend placement is not docked: {placement.value}")


def _overlay_axis_range(
    *,
    primary_min: float,
    primary_max: float,
    item_size: float,
    near_min: bool,
    boundary_ratio: float,
) -> tuple[float, float] | None:
    minimum = primary_min
    maximum = primary_max - item_size
    midpoint_start = minimum + (boundary_ratio * (maximum - minimum))
    if near_min:
        maximum = min(maximum, midpoint_start)
    else:
        minimum = max(minimum, midpoint_start)
    if minimum > maximum:
        return None
    return minimum, maximum


def _overlay_candidate_axis_values(
    *,
    minimum: float,
    maximum: float,
    item_size: float,
    obstacles: tuple[Aabb, ...],
    clearance: float,
    axis: str,
) -> tuple[float, ...]:
    values = {minimum, maximum}
    for obstacle in obstacles:
        if axis == "x":
            values.add(obstacle.min_x - clearance - item_size)
            values.add(obstacle.max_x + clearance)
        else:
            values.add(obstacle.min_y - clearance - item_size)
            values.add(obstacle.max_y + clearance)
    return tuple(sorted(value for value in values if minimum <= value <= maximum))


def _conflict_indices(
    bounds: Aabb,
    obstacles: tuple[Aabb, ...],
    clearance: float,
) -> tuple[int, ...]:
    return tuple(
        index
        for index, obstacle in enumerate(obstacles)
        if bounds.intersects(obstacle, clearance=clearance)
    )


def _fallback_overlay_candidates(
    primary: Aabb,
    legend: CompositionItem,
    placement: LegendPlacement,
    clearance: float,
) -> tuple[_WorkingPlacement, _WorkingPlacement]:
    bounds = legend.local_bounds
    left = placement in {LegendPlacement.UPPER_LEFT, LegendPlacement.LOWER_LEFT}
    upper = placement in {LegendPlacement.UPPER_LEFT, LegendPlacement.UPPER_RIGHT}
    aligned_x = (
        primary.min_x
        if left
        else primary.max_x - bounds.width
    )
    aligned_y = (
        primary.min_y
        if upper
        else primary.max_y - bounds.height
    )

    horizontal_x = (
        primary.min_x - clearance - bounds.width
        if left
        else primary.max_x + clearance
    )
    vertical_y = (
        primary.min_y - clearance - bounds.height
        if upper
        else primary.max_y + clearance
    )
    return (
        _align_min(legend, min_x=horizontal_x, min_y=aligned_y),
        _align_min(legend, min_x=aligned_x, min_y=vertical_y),
    )


def _place_overlay_legend(
    primary: Aabb,
    legend: CompositionItem,
    placement: LegendPlacement,
    obstacles: tuple[Aabb, ...],
    spacing: CompositionSpacing,
) -> tuple[_WorkingPlacement, tuple[int, ...], OverlayResolution]:
    clearance = spacing.overlay_clearance_px
    left = placement in {LegendPlacement.UPPER_LEFT, LegendPlacement.LOWER_LEFT}
    upper = placement in {LegendPlacement.UPPER_LEFT, LegendPlacement.UPPER_RIGHT}
    anchor_x = (
        primary.min_x
        if left
        else primary.max_x - legend.local_bounds.width
    )
    anchor_y = (
        primary.min_y
        if upper
        else primary.max_y - legend.local_bounds.height
    )
    anchor = _align_min(legend, min_x=anchor_x, min_y=anchor_y)
    initial_conflicts = _conflict_indices(anchor.bounds, obstacles, clearance)
    x_range = _overlay_axis_range(
        primary_min=primary.min_x,
        primary_max=primary.max_x,
        item_size=legend.local_bounds.width,
        near_min=left,
        boundary_ratio=OVERLAY_QUADRANT_BOUNDARY_RATIO,
    )
    y_range = _overlay_axis_range(
        primary_min=primary.min_y,
        primary_max=primary.max_y,
        item_size=legend.local_bounds.height,
        near_min=upper,
        boundary_ratio=OVERLAY_QUADRANT_BOUNDARY_RATIO,
    )

    if x_range is not None and y_range is not None:
        if not initial_conflicts:
            return anchor, (), OverlayResolution.ANCHORED

        x_values = _overlay_candidate_axis_values(
            minimum=x_range[0],
            maximum=x_range[1],
            item_size=legend.local_bounds.width,
            obstacles=obstacles,
            clearance=clearance,
            axis="x",
        )
        y_values = _overlay_candidate_axis_values(
            minimum=y_range[0],
            maximum=y_range[1],
            item_size=legend.local_bounds.height,
            obstacles=obstacles,
            clearance=clearance,
            axis="y",
        )

        def candidate_score(coordinate: tuple[float, float]) -> tuple[float, ...]:
            metrics = {
                "totalAnchorDistance": (
                    abs(coordinate[0] - anchor_x)
                    + abs(coordinate[1] - anchor_y)
                ),
                "xAnchorDistance": abs(coordinate[0] - anchor_x),
                "yAnchorDistance": abs(coordinate[1] - anchor_y),
                "nearEdgeX": coordinate[0] if left else -coordinate[0],
                "nearEdgeY": coordinate[1] if upper else -coordinate[1],
            }
            return tuple(metrics[name] for name in OVERLAY_CANDIDATE_SCORE_ORDER)

        coordinates = sorted(product(x_values, y_values), key=candidate_score)
        for min_x, min_y in coordinates:
            candidate = _align_min(legend, min_x=min_x, min_y=min_y)
            if not _conflict_indices(candidate.bounds, obstacles, clearance):
                return candidate, initial_conflicts, OverlayResolution.SHIFTED

    horizontal, vertical = _fallback_overlay_candidates(
        primary,
        legend,
        placement,
        clearance,
    )
    fallback_by_name = {"horizontal": horizontal, "vertical": vertical}
    fallback_candidates = tuple(
        fallback_by_name[name]
        for name in OVERLAY_CANVAS_GROWTH_CANDIDATE_ORDER
    )

    def growth_key(candidate: _WorkingPlacement) -> tuple[float, ...]:
        union = union_aabbs((primary, candidate.bounds))
        if union is None:  # pragma: no cover - primary always makes this non-empty
            raise RuntimeError("overlay union unexpectedly empty")
        added_area = (union.width * union.height) - (primary.width * primary.height)
        added_extent = (union.width - primary.width) + (union.height - primary.height)
        metrics = {
            "addedArea": added_area,
            "addedExtent": added_extent,
            "candidateOrder": fallback_candidates.index(candidate),
        }
        return tuple(metrics[name] for name in OVERLAY_CANVAS_GROWTH_SCORE_ORDER)

    fallback = min(fallback_candidates, key=growth_key)
    return fallback, initial_conflicts, OverlayResolution.CANVAS_GROWTH


def _place_title(
    primary: Aabb,
    title: CompositionItem,
    placement: TitlePlacement,
    spacing: CompositionSpacing,
    legend: _WorkingPlacement | None,
    legend_placement: LegendPlacement,
) -> _WorkingPlacement:
    title_bounds = title.local_bounds
    primary_center_x = 0.5 * (primary.min_x + primary.max_x)
    primary_center_y = 0.5 * (primary.min_y + primary.max_y)
    if placement is TitlePlacement.CENTER:
        return _align_center(
            title,
            center_x=primary_center_x,
            center_y=primary_center_y,
        )

    same_side = (
        legend is not None
        and (
            (placement is TitlePlacement.TOP and legend_placement is LegendPlacement.TOP)
            or (
                placement is TitlePlacement.BOTTOM
                and legend_placement is LegendPlacement.BOTTOM
            )
        )
    )
    if placement is TitlePlacement.TOP:
        target_max_y = (
            legend.bounds.min_y - spacing.stack_gap_px
            if same_side and legend is not None
            else primary.min_y - spacing.title_gap_px
        )
        result = _align_min(
            title,
            min_x=primary_center_x - (0.5 * title_bounds.width),
            min_y=target_max_y - title_bounds.height,
        )
    elif placement is TitlePlacement.BOTTOM:
        target_min_y = (
            legend.bounds.max_y + spacing.stack_gap_px
            if same_side and legend is not None
            else primary.max_y + spacing.title_gap_px
        )
        result = _align_min(
            title,
            min_x=primary_center_x - (0.5 * title_bounds.width),
            min_y=target_min_y,
        )
    else:
        raise ValueError(f"unsupported title placement: {placement.value}")

    if (
        legend is not None
        and legend_placement in {LegendPlacement.LEFT, LegendPlacement.RIGHT}
        and result.bounds.intersects(legend.bounds)
    ):
        if placement is TitlePlacement.TOP:
            result = _align_min(
                title,
                min_x=primary_center_x - (0.5 * title_bounds.width),
                min_y=legend.bounds.min_y - spacing.stack_gap_px - title_bounds.height,
            )
        else:
            result = _align_min(
                title,
                min_x=primary_center_x - (0.5 * title_bounds.width),
                min_y=legend.bounds.max_y + spacing.stack_gap_px,
            )
    return result


def plan_composition(request: CompositionRequest) -> CompositionPlan:
    """Resolve a complete immutable composition plan for ``request``."""
    if not isinstance(request, CompositionRequest):
        raise TypeError("request must be a CompositionRequest")

    primary = _WorkingPlacement(request.primary, 0.0, 0.0)
    working: list[_WorkingPlacement] = [primary]
    legend: _WorkingPlacement | None = None
    overlay_conflict_indices: tuple[int, ...] = ()
    overlay_resolution: OverlayResolution | None = None

    if (
        request.legend is not None
        and not request.legend.is_empty
        and request.legend_placement is not LegendPlacement.NONE
    ):
        if request.legend_placement.is_dock:
            legend = _place_docked_legend(
                primary.bounds,
                request.legend,
                request.legend_placement,
                request.spacing,
            )
        elif request.legend_placement.is_overlay:
            legend, overlay_conflict_indices, overlay_resolution = _place_overlay_legend(
                primary.bounds,
                request.legend,
                request.legend_placement,
                request.overlay_obstacles,
                request.spacing,
            )
        else:  # pragma: no cover - closed enum exhaustiveness
            raise ValueError(f"unsupported legend placement: {request.legend_placement.value}")
        working.append(legend)

    if (
        request.title is not None
        and not request.title.is_empty
        and request.title_placement is not TitlePlacement.NONE
    ):
        title = _place_title(
            primary.bounds,
            request.title,
            request.title_placement,
            request.spacing,
            legend,
            request.legend_placement,
        )
        working.append(title)

    painted_union = union_aabbs(placement.bounds for placement in working)
    if painted_union is None:  # pragma: no cover - validated primary is always present
        raise RuntimeError("composition has no painted bounds")
    padded_union = painted_union.expanded(request.spacing.edge_padding_px)
    outer_dx = -padded_union.min_x
    outer_dy = -padded_union.min_y

    placements = tuple(
        CompositionPlacement(
            role=placement.item.role,
            translation=(placement.dx + outer_dx, placement.dy + outer_dy),
            final_bounds=placement.bounds.translated(outer_dx, outer_dy),
        )
        for placement in working
    )
    primary_placement = placements[0]
    canvas_bounds = Aabb(0.0, 0.0, padded_union.width, padded_union.height)
    return CompositionPlan(
        canvas_bounds=canvas_bounds,
        view_box=(0.0, 0.0, canvas_bounds.width, canvas_bounds.height),
        primary_bounds=primary_placement.final_bounds,
        placements=placements,
        spacing=request.spacing,
        overlay_obstacles=tuple(
            obstacle.translated(outer_dx, outer_dy)
            for obstacle in request.overlay_obstacles
        ),
        overlay_conflict_indices=overlay_conflict_indices,
        overlay_resolution=overlay_resolution,
    )


__all__ = [
    "CompositionItem",
    "CompositionPlacement",
    "CompositionPlan",
    "CompositionRequest",
    "CompositionSpacing",
    "DEFAULT_COMPOSITION_SPACING",
    "LegendPlacement",
    "OVERLAY_CANDIDATE_SCORE_ORDER",
    "OVERLAY_CANVAS_GROWTH_CANDIDATE_ORDER",
    "OVERLAY_CANVAS_GROWTH_SCORE_ORDER",
    "OVERLAY_QUADRANT_BOUNDARY_RATIO",
    "OverlayResolution",
    "TitlePlacement",
    "plan_composition",
]
