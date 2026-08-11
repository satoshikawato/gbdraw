from __future__ import annotations

from types import SimpleNamespace

import pytest
from svgwrite import Drawing
from svgwrite.container import Group

from gbdraw.diagrams.circular.assemble import (
    _CircularPlotAssembly,
    _compose_circular_plot,
)
from gbdraw.layout.composition import LegendPlacement, TitlePlacement
from gbdraw.layout.spatial import Aabb, union_aabbs


def _center_x(bounds: Aabb) -> float:
    return 0.5 * (bounds.min_x + bounds.max_x)


def _center_y(bounds: Aabb) -> float:
    return 0.5 * (bounds.min_y + bounds.max_y)


def _plot(
    legend_placement: LegendPlacement,
    *,
    legend_bounds: Aabb | None = Aabb(-3.0, -5.0, 117.0, 55.0),
) -> tuple[_CircularPlotAssembly, Group]:
    drawing = Drawing(
        size=("1000px", "1000px"),
        viewBox="0 0 1000 1000",
        debug=False,
    )
    primary_target = Group(id="primary-target")
    primary_target.translate(500.0, 500.0)
    drawing.add(primary_target)

    legend_target = None
    if legend_bounds is not None and legend_bounds.width > 0.0 and legend_bounds.height > 0.0:
        legend_target = Group(id="legend-target")
        drawing.add(legend_target)

    canvas_config = SimpleNamespace(
        legend_position=legend_placement.value,
        total_width=1000.0,
        total_height=1000.0,
        offset_x=500.0,
        offset_y=500.0,
        legend_offset_x=0.0,
        legend_offset_y=0.0,
    )
    return (
        _CircularPlotAssembly(
            drawing=drawing,
            canvas_config=canvas_config,
            legend_measurement=SimpleNamespace(
                local_bounds=legend_bounds or Aabb(0.0, 0.0, 0.0, 0.0),
                reflow_metrics={
                    "colorRectSize": 24.0,
                    "lineHeight": 24.0 * (24.0 / 14.0),
                    "textXOffset": 24.0 * (22.0 / 14.0),
                },
            ),
            legend_table={},
            source_content_bounds=Aabb(100.0, 200.0, 500.0, 600.0),
            source_overlay_obstacles=(),
            primary_targets=(primary_target,),
            legend_target=legend_target,
            track_slot_geometry={"slots": ()},
        ),
        primary_target,
    )


@pytest.mark.circular
@pytest.mark.parametrize(
    "side",
    (
        LegendPlacement.LEFT,
        LegendPlacement.RIGHT,
        LegendPlacement.TOP,
        LegendPlacement.BOTTOM,
    ),
)
def test_circular_composition_uses_exact_dock_gap_and_edge_padding(
    side: LegendPlacement,
) -> None:
    plot, primary_target = _plot(side)

    result = _compose_circular_plot(plot)

    primary = result.content_bounds
    legend_placement = result.composition_plan.placement_for("legend")
    assert legend_placement is not None
    legend = legend_placement.final_bounds
    if side is LegendPlacement.LEFT:
        assert primary.min_x - legend.max_x == pytest.approx(24.0)
        assert _center_y(legend) == pytest.approx(_center_y(primary))
    elif side is LegendPlacement.RIGHT:
        assert legend.min_x - primary.max_x == pytest.approx(24.0)
        assert _center_y(legend) == pytest.approx(_center_y(primary))
    elif side is LegendPlacement.TOP:
        assert primary.min_y - legend.max_y == pytest.approx(24.0)
        assert _center_x(legend) == pytest.approx(_center_x(primary))
    else:
        assert legend.min_y - primary.max_y == pytest.approx(24.0)
        assert _center_x(legend) == pytest.approx(_center_x(primary))

    painted_bounds = union_aabbs(
        placement.final_bounds for placement in result.composition_plan.placements
    )
    assert painted_bounds == Aabb(
        16.0,
        16.0,
        result.composition_plan.width - 16.0,
        result.composition_plan.height - 16.0,
    )
    assert result.content_bounds == result.composition_plan.primary_bounds
    assert result.source_content_bounds == Aabb(100.0, 200.0, 500.0, 600.0)
    assert result.drawing._gbdraw_circular_assembly_result is result

    primary_plan = result.composition_plan.placement_for("primary")
    assert primary_plan is not None
    assert str(primary_target.attribs["transform"]).endswith("translate(500.0,500.0)")
    assert str(primary_target.attribs["transform"]).startswith(
        f"translate({primary_plan.dx},{primary_plan.dy})"
    )


@pytest.mark.circular
def test_circular_composition_omits_empty_legend_and_places_title_once() -> None:
    plot, _primary_target = _plot(LegendPlacement.RIGHT, legend_bounds=None)
    title_target = Group(id="plot-title")
    title_bounds = Aabb(-80.0, -10.0, 80.0, 10.0)

    result = _compose_circular_plot(
        plot,
        title_target=title_target,
        title_bounds=title_bounds,
        title_placement=TitlePlacement.TOP,
    )

    assert result.legend_target is None
    assert result.composition_plan.placement_for("legend") is None
    title_placement = result.composition_plan.placement_for("title")
    assert title_placement is not None
    assert result.content_bounds.min_y - title_placement.final_bounds.max_y == pytest.approx(20.0)
    assert _center_x(title_placement.final_bounds) == pytest.approx(_center_x(result.content_bounds))
    assert sum(element is title_target for element in result.drawing.elements) == 1
    assert title_target.attribs["data-gbdraw-composition-role"] == "title"


@pytest.mark.circular
@pytest.mark.parametrize(
    "side",
    (
        LegendPlacement.UPPER_LEFT,
        LegendPlacement.UPPER_RIGHT,
        LegendPlacement.LOWER_LEFT,
        LegendPlacement.LOWER_RIGHT,
    ),
)
def test_circular_overlay_compatibility_offsets_copy_final_plan(
    side: LegendPlacement,
) -> None:
    plot, _primary_target = _plot(side)

    result = _compose_circular_plot(plot)

    legend = result.composition_plan.placement_for("legend")
    primary = result.composition_plan.placement_for("primary")
    assert legend is not None
    assert primary is not None
    assert plot.canvas_config.legend_offset_x == legend.dx
    assert plot.canvas_config.legend_offset_y == legend.dy
    assert plot.canvas_config.offset_x == 500.0 + primary.dx
    assert plot.canvas_config.offset_y == 500.0 + primary.dy
