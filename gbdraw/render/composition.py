"""Apply pure composition plans to internal SVG drawing targets."""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from typing import Any

from svgwrite.container import Group
from svgwrite.drawing import Drawing
from svgwrite.params import Parameter

from gbdraw.layout.composition import (
    CompositionPlacement,
    CompositionPlan,
    LegendPlacement,
    OVERLAY_CANDIDATE_SCORE_ORDER,
    OVERLAY_CANVAS_GROWTH_CANDIDATE_ORDER,
    OVERLAY_CANVAS_GROWTH_SCORE_ORDER,
    OVERLAY_QUADRANT_BOUNDARY_RATIO,
    TitlePlacement,
)
from gbdraw.layout.spatial import Aabb

COMPOSITION_SCHEMA_VERSION = 1
COMPOSITION_SCHEMA_ATTRIBUTE = "data-gbdraw-composition-schema"
COMPOSITION_METADATA_ATTRIBUTE = "data-gbdraw-composition"
COMPOSITION_ROLE_ATTRIBUTE = "data-gbdraw-composition-role"

_PRIMARY_ROLE = "primary"
_LEGEND_ROLE = "legend"
_TITLE_ROLE = "title"
_ROLE_SELECTORS = {
    role: f'[{COMPOSITION_ROLE_ATTRIBUTE}="{role}"]'
    for role in (_PRIMARY_ROLE, _LEGEND_ROLE, _TITLE_ROLE)
}


def _wire_number(value: float) -> float:
    value = float(value)
    return 0.0 if value == 0.0 else value


def _bounds_payload(bounds: Aabb) -> dict[str, float]:
    return {
        "height": _wire_number(bounds.height),
        "width": _wire_number(bounds.width),
        "x": _wire_number(bounds.min_x),
        "y": _wire_number(bounds.min_y),
    }


def _translation_payload(placement: CompositionPlacement) -> list[float]:
    return [_wire_number(placement.dx), _wire_number(placement.dy)]


def _local_bounds(placement: CompositionPlacement) -> Aabb:
    return placement.final_bounds.translated(-placement.dx, -placement.dy)


def _target_payload(
    placement: CompositionPlacement,
    *,
    composition_role: str,
    include_final_bounds: bool,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "automaticTranslation": _translation_payload(placement),
        "role": composition_role,
        "selector": _ROLE_SELECTORS[composition_role],
    }
    bounds = (
        placement.final_bounds if include_final_bounds else _local_bounds(placement)
    )
    payload["finalBounds" if include_final_bounds else "localBounds"] = _bounds_payload(
        bounds
    )
    return payload


def _metadata_payload(
    plan: CompositionPlan,
    *,
    primary_placement: CompositionPlacement,
    legend_placement: CompositionPlacement | None,
    title_placement: CompositionPlacement | None,
    legend_side: LegendPlacement,
    title_side: TitlePlacement,
    legend_reflow_metrics: Mapping[str, float] | None,
) -> str:
    spacing = plan.spacing
    payload = {
        "legend": (
            _target_payload(
                legend_placement,
                composition_role=_LEGEND_ROLE,
                include_final_bounds=False,
            )
            if legend_placement is not None
            else None
        ),
        "legendReflow": (
            {
                key: _wire_number(value)
                for key, value in legend_reflow_metrics.items()
            }
            if legend_reflow_metrics is not None
            else None
        ),
        "legendSide": legend_side.value,
        "overlayObstacles": [
            _bounds_payload(bounds) for bounds in plan.overlay_obstacles
        ],
        "overlayPolicy": {
            "candidateScoreOrder": list(OVERLAY_CANDIDATE_SCORE_ORDER),
            "canvasGrowthCandidateOrder": list(
                OVERLAY_CANVAS_GROWTH_CANDIDATE_ORDER
            ),
            "canvasGrowthScoreOrder": list(OVERLAY_CANVAS_GROWTH_SCORE_ORDER),
            "quadrantBoundaryRatio": _wire_number(
                OVERLAY_QUADRANT_BOUNDARY_RATIO
            ),
        },
        "primary": _target_payload(
            primary_placement,
            composition_role=_PRIMARY_ROLE,
            include_final_bounds=True,
        ),
        "spacing": {
            "dockGapPx": _wire_number(spacing.dock_gap_px),
            "edgePaddingPx": _wire_number(spacing.edge_padding_px),
            "overlayClearancePx": _wire_number(spacing.overlay_clearance_px),
            "stackGapPx": _wire_number(spacing.stack_gap_px),
            "titleGapPx": _wire_number(spacing.title_gap_px),
        },
        "title": (
            _target_payload(
                title_placement,
                composition_role=_TITLE_ROLE,
                include_final_bounds=False,
            )
            if title_placement is not None
            else None
        ),
        "titleSide": title_side.value,
    }
    return json.dumps(
        payload,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    )


def _placement_for_role(
    plan: CompositionPlan, role: str
) -> CompositionPlacement | None:
    matches = [placement for placement in plan.placements if placement.role == role]
    if len(matches) > 1:
        raise ValueError(f"composition plan contains duplicate role {role!r}")
    return matches[0] if matches else None


def _require_group(target: object, *, name: str) -> Group:
    if not isinstance(target, Group):
        raise TypeError(f"{name} must be an svgwrite Group")
    return target


def _validate_targets(
    drawing: Drawing,
    *,
    primary_targets: Sequence[Group],
    legend_target: Group | None,
    title_target: Group | None,
) -> tuple[Group, ...]:
    if not isinstance(drawing, Drawing):
        raise TypeError("drawing must be an svgwrite Drawing")
    targets = tuple(
        _require_group(target, name="primary target") for target in primary_targets
    )
    if not targets:
        raise ValueError("at least one primary target is required")
    if legend_target is not None:
        targets += (_require_group(legend_target, name="legend target"),)
    if title_target is not None:
        targets += (_require_group(title_target, name="title target"),)
    if len({id(target) for target in targets}) != len(targets):
        raise ValueError("composition targets must be distinct")
    if drawing.debug:
        raise ValueError("composition metadata requires a drawing with debug=False")
    return targets


def _allow_internal_role_attribute(target: Group) -> None:
    """Disable SVG 1.1 validation only on a role-marked target group.

    Several existing render groups are constructed directly instead of through
    the ``Drawing`` factory and therefore retain svgwrite's default strict
    validator.  SVG 1.1 does not know the internal ``data-gbdraw-*`` role
    attribute, so those otherwise valid groups must use svgwrite's permissive
    serialization mode once the renderer marks them.  Child elements keep
    their own validation parameters.
    """
    if target.debug:
        target.set_parameter(Parameter(debug=False, profile=target.profile))


def _ensure_not_applied(drawing: Drawing, targets: Sequence[Group]) -> None:
    root_attributes = drawing.attribs
    if (
        COMPOSITION_SCHEMA_ATTRIBUTE in root_attributes
        or COMPOSITION_METADATA_ATTRIBUTE in root_attributes
    ):
        raise RuntimeError("composition plan has already been applied to this drawing")
    if any(COMPOSITION_ROLE_ATTRIBUTE in target.attribs for target in targets):
        raise RuntimeError("a composition target has already been marked")


def _compose_outer_translation(target: Group, placement: CompositionPlacement) -> None:
    translation = (
        f"translate({_wire_number(placement.dx)},{_wire_number(placement.dy)})"
    )
    existing = str(target.attribs.get("transform", "")).strip()
    target.attribs["transform"] = f"{translation} {existing}".strip()


def _apply_target(
    target: Group,
    placement: CompositionPlacement,
    *,
    composition_role: str,
) -> None:
    _allow_internal_role_attribute(target)
    _compose_outer_translation(target, placement)
    target.attribs[COMPOSITION_ROLE_ATTRIBUTE] = composition_role


def _coerce_legend_side(value: LegendPlacement | str) -> LegendPlacement:
    try:
        return value if isinstance(value, LegendPlacement) else LegendPlacement(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"unknown legend side: {value!r}") from exc


def _coerce_title_side(value: TitlePlacement | str) -> TitlePlacement:
    try:
        return value if isinstance(value, TitlePlacement) else TitlePlacement(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"unknown title side: {value!r}") from exc


def _validate_legend_reflow_metrics(
    value: Mapping[str, float] | None,
    *,
    required: bool,
) -> Mapping[str, float] | None:
    if value is None:
        if required:
            raise ValueError("legend target requires legend reflow metrics")
        return None
    if not required:
        raise ValueError("legend reflow metrics require a legend target")
    required_keys = {"colorRectSize", "lineHeight", "textXOffset"}
    if set(value) != required_keys:
        raise ValueError("legend reflow metrics have an invalid field set")
    if any(
        isinstance(raw, bool) or not isinstance(raw, (int, float))
        for raw in value.values()
    ):
        raise ValueError("legend reflow metrics must be finite and positive")
    normalized = {key: float(raw) for key, raw in value.items()}
    if any(not math.isfinite(number) or number <= 0.0 for number in normalized.values()):
        raise ValueError("legend reflow metrics must be finite and positive")
    return normalized


def apply_composition_plan(
    drawing: Drawing,
    plan: CompositionPlan,
    *,
    primary_targets: Sequence[Group],
    primary_role: str = _PRIMARY_ROLE,
    legend_target: Group | None = None,
    legend_role: str = _LEGEND_ROLE,
    legend_side: LegendPlacement | str = LegendPlacement.NONE,
    legend_reflow_metrics: Mapping[str, float] | None = None,
    title_target: Group | None = None,
    title_role: str = _TITLE_ROLE,
    title_side: TitlePlacement | str = TitlePlacement.NONE,
) -> None:
    """Apply ``plan`` once without changing the drawing's group topology.

    ``primary_targets`` share one planner placement, which lets callers retain
    several existing top-level primary groups instead of inserting a wrapper.
    The role arguments name placements in ``plan``; the SVG bookkeeping roles
    remain the stable internal values ``primary``, ``legend``, and ``title``.
    """
    if not isinstance(plan, CompositionPlan):
        raise TypeError("plan must be a CompositionPlan")
    all_targets = _validate_targets(
        drawing,
        primary_targets=primary_targets,
        legend_target=legend_target,
        title_target=title_target,
    )
    resolved_legend_side = _coerce_legend_side(legend_side)
    resolved_title_side = _coerce_title_side(title_side)
    resolved_legend_reflow = _validate_legend_reflow_metrics(
        legend_reflow_metrics,
        required=legend_target is not None,
    )

    roles = (primary_role, legend_role, title_role)
    if any(not isinstance(role, str) or not role.strip() for role in roles):
        raise ValueError("composition placement roles must be non-empty strings")
    if len(set(roles)) != len(roles):
        raise ValueError("primary, legend, and title placement roles must be distinct")

    primary_placement = _placement_for_role(plan, primary_role)
    if primary_placement is None:
        raise ValueError(f"composition plan has no primary role {primary_role!r}")
    if primary_placement.final_bounds != plan.primary_bounds:
        raise ValueError("primary placement does not match plan.primary_bounds")
    legend_placement = _placement_for_role(plan, legend_role)
    title_placement = _placement_for_role(plan, title_role)
    if (legend_target is None) != (legend_placement is None):
        raise ValueError(
            "legend target and plan placement must either both exist or be absent"
        )
    if (title_target is None) != (title_placement is None):
        raise ValueError(
            "title target and plan placement must either both exist or be absent"
        )
    expected_roles = {
        primary_role,
        *(() if legend_placement is None else (legend_role,)),
        *(() if title_placement is None else (title_role,)),
    }
    actual_roles = {placement.role for placement in plan.placements}
    if actual_roles != expected_roles:
        raise ValueError("composition plan contains unbound placement roles")

    _ensure_not_applied(drawing, all_targets)
    metadata = _metadata_payload(
        plan,
        primary_placement=primary_placement,
        legend_placement=legend_placement,
        title_placement=title_placement,
        legend_side=resolved_legend_side,
        title_side=resolved_title_side,
        legend_reflow_metrics=resolved_legend_reflow,
    )

    primary_count = len(primary_targets)
    for target in all_targets[:primary_count]:
        _apply_target(target, primary_placement, composition_role=_PRIMARY_ROLE)
    if legend_target is not None and legend_placement is not None:
        _apply_target(legend_target, legend_placement, composition_role=_LEGEND_ROLE)
    if title_target is not None and title_placement is not None:
        _apply_target(title_target, title_placement, composition_role=_TITLE_ROLE)

    drawing.attribs["width"] = f"{_wire_number(plan.width)}px"
    drawing.attribs["height"] = f"{_wire_number(plan.height)}px"
    drawing.attribs["viewBox"] = " ".join(
        str(_wire_number(value)) for value in plan.view_box
    )
    drawing.attribs[COMPOSITION_SCHEMA_ATTRIBUTE] = str(COMPOSITION_SCHEMA_VERSION)
    drawing.attribs[COMPOSITION_METADATA_ATTRIBUTE] = metadata


__all__ = [
    "COMPOSITION_METADATA_ATTRIBUTE",
    "COMPOSITION_ROLE_ATTRIBUTE",
    "COMPOSITION_SCHEMA_ATTRIBUTE",
    "COMPOSITION_SCHEMA_VERSION",
    "apply_composition_plan",
]
