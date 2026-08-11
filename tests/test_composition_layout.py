from __future__ import annotations

from dataclasses import FrozenInstanceError

import pytest

from gbdraw.layout.composition import (
    DEFAULT_COMPOSITION_SPACING,
    CompositionItem,
    CompositionRequest,
    CompositionSpacing,
    LegendPlacement,
    OverlayResolution,
    TitlePlacement,
    plan_composition,
)
from gbdraw.layout.spatial import Aabb, union_aabbs


PRIMARY_BOUNDS = Aabb(10.0, 20.0, 110.0, 100.0)
LEGEND_BOUNDS = Aabb(-2.0, -5.0, 28.0, 15.0)
TITLE_BOUNDS = Aabb(-5.0, -10.0, 55.0, 10.0)


def _item(role: str, bounds: Aabb) -> CompositionItem:
    return CompositionItem(role=role, local_bounds=bounds)


def _request(**changes: object) -> CompositionRequest:
    values: dict[str, object] = {
        "primary": _item("primary", PRIMARY_BOUNDS),
    }
    values.update(changes)
    return CompositionRequest(**values)  # type: ignore[arg-type]


def _bounds(plan, role: str) -> Aabb:
    placement = plan.placement_for(role)
    assert placement is not None
    return placement.final_bounds


def _painted_union(plan) -> Aabb:
    result = union_aabbs(placement.final_bounds for placement in plan.placements)
    assert result is not None
    return result


def test_aabb_pure_operations_and_empty_union() -> None:
    normalized = Aabb(8, 7, 2, 3)
    assert normalized == Aabb(2, 3, 8, 7)
    assert normalized.width == 6.0
    assert normalized.height == 4.0
    assert normalized.translated(-4, 5) == Aabb(-2, 8, 4, 12)
    assert normalized.expanded(2) == Aabb(0, 1, 10, 9)
    assert normalized.expanded(2, 3) == Aabb(0, 0, 10, 10)
    assert union_aabbs([]) is None
    assert union_aabbs(box for box in (normalized, Aabb(-1, 4, 4, 12))) == Aabb(
        -1,
        3,
        8,
        12,
    )


def test_aabb_intersection_allows_edge_contact_and_exact_clearance() -> None:
    left = Aabb(0, 0, 10, 10)
    touching = Aabb(10, 0, 20, 10)
    eight_away = Aabb(18, 0, 28, 10)
    assert not left.intersects(touching)
    assert not left.intersects(eight_away, clearance=8)
    assert left.intersects(eight_away, clearance=8.001)


@pytest.mark.parametrize(
    ("operation", "match"),
    [
        (lambda: Aabb(0, 0, 1, 1).translated(float("nan"), 0), "translation"),
        (lambda: Aabb(0, 0, 1, 1).expanded(-1), "non-negative"),
        (
            lambda: Aabb(0, 0, 1, 1).intersects(Aabb(2, 2, 3, 3), clearance=-1),
            "non-negative",
        ),
        (lambda: union_aabbs([Aabb(0, 0, 1, 1), object()]), "Aabb"),
    ],
)
def test_aabb_operations_validate_inputs(operation, match: str) -> None:
    with pytest.raises((TypeError, ValueError), match=match):
        operation()


def test_shared_spacing_defaults_are_frozen() -> None:
    assert DEFAULT_COMPOSITION_SPACING == CompositionSpacing(
        edge_padding_px=16,
        dock_gap_px=24,
        title_gap_px=20,
        stack_gap_px=20,
        overlay_clearance_px=8,
    )
    with pytest.raises(FrozenInstanceError):
        DEFAULT_COMPOSITION_SPACING.edge_padding_px = 99  # type: ignore[misc]


@pytest.mark.parametrize(
    ("placement", "expected_primary", "expected_legend", "expected_canvas"),
    [
        (
            LegendPlacement.RIGHT,
            Aabb(16, 16, 116, 96),
            Aabb(140, 46, 170, 66),
            Aabb(0, 0, 186, 112),
        ),
        (
            LegendPlacement.LEFT,
            Aabb(70, 16, 170, 96),
            Aabb(16, 46, 46, 66),
            Aabb(0, 0, 186, 112),
        ),
        (
            LegendPlacement.TOP,
            Aabb(16, 60, 116, 140),
            Aabb(51, 16, 81, 36),
            Aabb(0, 0, 132, 156),
        ),
        (
            LegendPlacement.BOTTOM,
            Aabb(16, 16, 116, 96),
            Aabb(51, 120, 81, 140),
            Aabb(0, 0, 132, 156),
        ),
    ],
)
def test_exact_dock_gap_and_outer_padding(
    placement: LegendPlacement,
    expected_primary: Aabb,
    expected_legend: Aabb,
    expected_canvas: Aabb,
) -> None:
    plan = plan_composition(
        _request(
            legend=_item("legend", LEGEND_BOUNDS),
            legend_placement=placement,
        )
    )
    assert plan.primary_bounds == expected_primary
    assert _bounds(plan, "legend") == expected_legend
    assert plan.canvas_bounds == expected_canvas
    assert plan.view_box == (0.0, 0.0, expected_canvas.width, expected_canvas.height)

    primary = plan.primary_bounds
    legend = _bounds(plan, "legend")
    if placement is LegendPlacement.RIGHT:
        assert legend.min_x - primary.max_x == 24.0
        assert (legend.min_y + legend.max_y) == (primary.min_y + primary.max_y)
    elif placement is LegendPlacement.LEFT:
        assert primary.min_x - legend.max_x == 24.0
        assert (legend.min_y + legend.max_y) == (primary.min_y + primary.max_y)
    elif placement is LegendPlacement.TOP:
        assert primary.min_y - legend.max_y == 24.0
        assert (legend.min_x + legend.max_x) == (primary.min_x + primary.max_x)
    else:
        assert legend.min_y - primary.max_y == 24.0
        assert (legend.min_x + legend.max_x) == (primary.min_x + primary.max_x)

    painted = _painted_union(plan)
    assert painted.min_x == 16.0
    assert painted.min_y == 16.0
    assert plan.width - painted.max_x == 16.0
    assert plan.height - painted.max_y == 16.0


@pytest.mark.parametrize(
    ("placement", "expected_primary", "expected_title", "expected_canvas"),
    [
        (
            TitlePlacement.TOP,
            Aabb(16, 56, 116, 136),
            Aabb(36, 16, 96, 36),
            Aabb(0, 0, 132, 152),
        ),
        (
            TitlePlacement.BOTTOM,
            Aabb(16, 16, 116, 96),
            Aabb(36, 116, 96, 136),
            Aabb(0, 0, 132, 152),
        ),
        (
            TitlePlacement.CENTER,
            Aabb(16, 16, 116, 96),
            Aabb(36, 46, 96, 66),
            Aabb(0, 0, 132, 112),
        ),
    ],
)
def test_title_is_centered_on_primary_content(
    placement: TitlePlacement,
    expected_primary: Aabb,
    expected_title: Aabb,
    expected_canvas: Aabb,
) -> None:
    plan = plan_composition(
        _request(
            title=_item("title", TITLE_BOUNDS),
            title_placement=placement,
        )
    )
    assert plan.primary_bounds == expected_primary
    assert _bounds(plan, "title") == expected_title
    assert plan.canvas_bounds == expected_canvas
    assert (
        _bounds(plan, "title").min_x + _bounds(plan, "title").max_x
        == plan.primary_bounds.min_x + plan.primary_bounds.max_x
    )


@pytest.mark.parametrize(
    ("legend_placement", "title_placement", "expected_primary", "expected_legend", "expected_title"),
    [
        (
            LegendPlacement.TOP,
            TitlePlacement.TOP,
            Aabb(16, 100, 116, 180),
            Aabb(51, 56, 81, 76),
            Aabb(36, 16, 96, 36),
        ),
        (
            LegendPlacement.BOTTOM,
            TitlePlacement.BOTTOM,
            Aabb(16, 16, 116, 96),
            Aabb(51, 120, 81, 140),
            Aabb(36, 160, 96, 180),
        ),
    ],
)
def test_same_side_title_is_outer_of_legend(
    legend_placement: LegendPlacement,
    title_placement: TitlePlacement,
    expected_primary: Aabb,
    expected_legend: Aabb,
    expected_title: Aabb,
) -> None:
    plan = plan_composition(
        _request(
            legend=_item("legend", LEGEND_BOUNDS),
            legend_placement=legend_placement,
            title=_item("title", TITLE_BOUNDS),
            title_placement=title_placement,
        )
    )
    assert plan.canvas_bounds == Aabb(0, 0, 132, 196)
    assert plan.primary_bounds == expected_primary
    assert _bounds(plan, "legend") == expected_legend
    assert _bounds(plan, "title") == expected_title
    if title_placement is TitlePlacement.TOP:
        assert expected_primary.min_y - expected_legend.max_y == 24.0
        assert expected_legend.min_y - expected_title.max_y == 20.0
    else:
        assert expected_legend.min_y - expected_primary.max_y == 24.0
        assert expected_title.min_y - expected_legend.max_y == 20.0


def test_long_top_title_moves_outward_from_tall_side_legend() -> None:
    plan = plan_composition(
        _request(
            legend=_item("legend", Aabb(0, 0, 30, 200)),
            legend_placement=LegendPlacement.RIGHT,
            title=_item("title", Aabb(0, 0, 200, 20)),
            title_placement=TitlePlacement.TOP,
        )
    )
    primary = plan.primary_bounds
    legend = _bounds(plan, "legend")
    title = _bounds(plan, "title")
    assert legend.min_x - primary.max_x == 24.0
    assert legend.min_y + legend.max_y == primary.min_y + primary.max_y
    assert legend.min_y - title.max_y == 20.0
    assert not title.intersects(legend)


@pytest.mark.parametrize(
    ("legend_bounds", "placement", "expected_canvas", "expected_primary", "expected_legend"),
    [
        (
            Aabb(0, 0, 30, 200),
            LegendPlacement.RIGHT,
            Aabb(0, 0, 186, 232),
            Aabb(16, 76, 116, 156),
            Aabb(140, 16, 170, 216),
        ),
        (
            Aabb(0, 0, 200, 20),
            LegendPlacement.TOP,
            Aabb(0, 0, 232, 156),
            Aabb(66, 60, 166, 140),
            Aabb(16, 16, 216, 36),
        ),
    ],
)
def test_tall_and_wide_decorations_drive_the_union(
    legend_bounds: Aabb,
    placement: LegendPlacement,
    expected_canvas: Aabb,
    expected_primary: Aabb,
    expected_legend: Aabb,
) -> None:
    plan = plan_composition(
        _request(
            legend=_item("legend", legend_bounds),
            legend_placement=placement,
        )
    )
    assert plan.canvas_bounds == expected_canvas
    assert plan.primary_bounds == expected_primary
    assert _bounds(plan, "legend") == expected_legend


def test_docked_placements_are_mirror_symmetric() -> None:
    legend = _item("legend", LEGEND_BOUNDS)
    left = plan_composition(_request(legend=legend, legend_placement="left"))
    right = plan_composition(_request(legend=legend, legend_placement="right"))
    assert left.canvas_bounds == right.canvas_bounds
    for role in ("primary", "legend"):
        left_bounds = _bounds(left, role)
        right_bounds = _bounds(right, role)
        assert left_bounds.min_x == right.width - right_bounds.max_x
        assert left_bounds.max_x == right.width - right_bounds.min_x
        assert (left_bounds.min_y, left_bounds.max_y) == (
            right_bounds.min_y,
            right_bounds.max_y,
        )

    top = plan_composition(_request(legend=legend, legend_placement="top"))
    bottom = plan_composition(_request(legend=legend, legend_placement="bottom"))
    assert top.canvas_bounds == bottom.canvas_bounds
    for role in ("primary", "legend"):
        top_bounds = _bounds(top, role)
        bottom_bounds = _bounds(bottom, role)
        assert top_bounds.min_y == bottom.height - bottom_bounds.max_y
        assert top_bounds.max_y == bottom.height - bottom_bounds.min_y
        assert (top_bounds.min_x, top_bounds.max_x) == (
            bottom_bounds.min_x,
            bottom_bounds.max_x,
        )


def test_negative_primary_and_decoration_coordinates_do_not_change_layout() -> None:
    positive = plan_composition(
        _request(
            legend=_item("legend", LEGEND_BOUNDS),
            legend_placement=LegendPlacement.RIGHT,
        )
    )
    negative = plan_composition(
        CompositionRequest(
            primary=_item("primary", Aabb(-40, -10, 60, 70)),
            legend=_item("legend", LEGEND_BOUNDS),
            legend_placement=LegendPlacement.RIGHT,
        )
    )
    assert negative.canvas_bounds == positive.canvas_bounds
    assert negative.primary_bounds == positive.primary_bounds
    assert _bounds(negative, "legend") == _bounds(positive, "legend")


@pytest.mark.parametrize(
    "changes",
    [
        {"legend": None, "legend_placement": LegendPlacement.RIGHT},
        {
            "legend": _item("legend", Aabb(0, 0, 0, 20)),
            "legend_placement": LegendPlacement.RIGHT,
        },
        {
            "legend": _item("legend", LEGEND_BOUNDS),
            "legend_placement": LegendPlacement.NONE,
        },
        {
            "title": _item("title", Aabb(0, 0, 20, 0)),
            "title_placement": TitlePlacement.TOP,
        },
        {
            "title": _item("title", TITLE_BOUNDS),
            "title_placement": TitlePlacement.NONE,
        },
    ],
)
def test_empty_or_absent_decorations_reserve_no_space(changes: dict[str, object]) -> None:
    primary_only = plan_composition(_request())
    plan = plan_composition(_request(**changes))
    assert plan == primary_only
    assert plan.canvas_bounds == Aabb(0, 0, 132, 112)
    assert plan.primary_bounds == Aabb(16, 16, 116, 96)


@pytest.mark.parametrize(
    ("placement", "expected"),
    [
        (LegendPlacement.UPPER_LEFT, Aabb(16, 16, 36, 26)),
        (LegendPlacement.UPPER_RIGHT, Aabb(96, 16, 116, 26)),
        (LegendPlacement.LOWER_LEFT, Aabb(16, 106, 36, 116)),
        (LegendPlacement.LOWER_RIGHT, Aabb(96, 106, 116, 116)),
    ],
)
def test_overlay_legends_anchor_to_requested_primary_corner(
    placement: LegendPlacement,
    expected: Aabb,
) -> None:
    plan = plan_composition(
        CompositionRequest(
            primary=_item("primary", Aabb(0, 0, 100, 100)),
            legend=_item("legend", Aabb(-3, -2, 17, 8)),
            legend_placement=placement,
        )
    )
    assert plan.canvas_bounds == Aabb(0, 0, 132, 132)
    assert _bounds(plan, "legend") == expected
    assert plan.overlay_conflict_indices == ()
    assert plan.overlay_resolution is OverlayResolution.ANCHORED


def test_overlay_reports_initial_conflicts_then_moves_within_quadrant() -> None:
    obstacles = (
        Aabb(70, 70, 90, 90),
        Aabb(0, 0, 30, 30),
    )
    plan = plan_composition(
        CompositionRequest(
            primary=_item("primary", Aabb(0, 0, 100, 100)),
            legend=_item("legend", Aabb(0, 0, 20, 20)),
            legend_placement=LegendPlacement.UPPER_LEFT,
            overlay_obstacles=obstacles,
        )
    )
    legend = _bounds(plan, "legend")
    assert plan.overlay_conflict_indices == (1,)
    assert plan.overlay_resolution is OverlayResolution.SHIFTED
    assert legend == Aabb(16, 54, 36, 74)
    assert 0.5 * (legend.min_x + legend.max_x) <= 0.5 * (
        plan.primary_bounds.min_x + plan.primary_bounds.max_x
    )
    assert 0.5 * (legend.min_y + legend.max_y) <= 0.5 * (
        plan.primary_bounds.min_y + plan.primary_bounds.max_y
    )
    translated_obstacles = tuple(obstacle.translated(16, 16) for obstacle in obstacles)
    assert plan.overlay_obstacles == translated_obstacles
    assert all(
        not legend.intersects(obstacle, clearance=8)
        for obstacle in translated_obstacles
    )


def test_overlay_canvas_growth_fallback_is_exact_and_deterministic() -> None:
    request = CompositionRequest(
        primary=_item("primary", Aabb(0, 0, 40, 40)),
        legend=_item("legend", Aabb(0, 0, 20, 20)),
        legend_placement=LegendPlacement.UPPER_LEFT,
        overlay_obstacles=(Aabb(0, 0, 40, 40),),
    )
    first = plan_composition(request)
    second = plan_composition(request)
    assert first == second
    assert first.overlay_conflict_indices == (0,)
    assert first.overlay_resolution is OverlayResolution.CANVAS_GROWTH
    assert first.canvas_bounds == Aabb(0, 0, 100, 72)
    assert _bounds(first, "legend") == Aabb(16, 16, 36, 36)
    assert first.primary_bounds == Aabb(44, 16, 84, 56)
    assert first.primary_bounds.min_x - _bounds(first, "legend").max_x == 8.0


def test_request_freezes_obstacles_and_planning_does_not_mutate_inputs() -> None:
    obstacle_list = [Aabb(0, 0, 10, 10)]
    primary = _item("primary", Aabb(-5, -7, 95, 73))
    legend = _item("legend", Aabb(-2, -1, 18, 9))
    request = CompositionRequest(
        primary=primary,
        legend=legend,
        legend_placement="upper_right",
        overlay_obstacles=obstacle_list,  # type: ignore[arg-type]
    )
    obstacle_list.append(Aabb(50, 50, 60, 60))
    before = (request.primary, request.legend, request.overlay_obstacles)
    first = plan_composition(request)
    second = plan_composition(request)
    assert first == second
    assert (request.primary, request.legend, request.overlay_obstacles) == before
    assert request.overlay_obstacles == (Aabb(0, 0, 10, 10),)
    assert isinstance(request.legend_placement, LegendPlacement)
    with pytest.raises(FrozenInstanceError):
        request.legend = None  # type: ignore[misc]


@pytest.mark.parametrize(
    "field_name",
    [
        "edge_padding_px",
        "dock_gap_px",
        "title_gap_px",
        "stack_gap_px",
        "overlay_clearance_px",
    ],
)
@pytest.mark.parametrize("invalid", [-1.0, float("nan"), float("inf")])
def test_spacing_rejects_negative_or_nonfinite_values(
    field_name: str,
    invalid: float,
) -> None:
    values = {
        "edge_padding_px": 16.0,
        "dock_gap_px": 24.0,
        "title_gap_px": 20.0,
        "stack_gap_px": 20.0,
        "overlay_clearance_px": 8.0,
    }
    values[field_name] = invalid
    with pytest.raises(ValueError, match=field_name):
        CompositionSpacing(**values)


@pytest.mark.parametrize(
    ("changes", "match"),
    [
        ({"primary": None}, "primary"),
        ({"primary": _item("primary", Aabb(0, 0, 0, 10))}, "positive"),
        ({"legend_placement": "beside"}, "unknown legend"),
        ({"title_placement": "middle"}, "unknown title"),
        ({"overlay_obstacles": (object(),)}, "Aabb"),
        ({"legend": _item("primary", LEGEND_BOUNDS)}, "unique"),
    ],
)
def test_request_validation_is_explicit(changes: dict[str, object], match: str) -> None:
    with pytest.raises((TypeError, ValueError), match=match):
        _request(**changes)


def test_unknown_placements_fail_even_when_decorations_are_absent() -> None:
    with pytest.raises(ValueError, match="unknown legend"):
        _request(legend=None, legend_placement="outside")
    with pytest.raises(ValueError, match="unknown title"):
        _request(title=None, title_placement="above")
