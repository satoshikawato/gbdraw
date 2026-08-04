"""Emit Python composition plans and their schema-1 Web replay payloads."""

from __future__ import annotations

import json
import sys
from pathlib import Path

from svgwrite import Drawing

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from gbdraw.layout.composition import (  # noqa: E402
    CompositionItem,
    CompositionRequest,
    LegendPlacement,
    TitlePlacement,
    plan_composition,
)
from gbdraw.layout.spatial import Aabb  # noqa: E402
from gbdraw.render.composition import (  # noqa: E402
    COMPOSITION_METADATA_ATTRIBUTE,
    COMPOSITION_SCHEMA_ATTRIBUTE,
    apply_composition_plan,
)

LEGEND_REFLOW = {
    "colorRectSize": 14.0,
    "lineHeight": 24.0,
    "textXOffset": 22.0,
}


def _wire(value: float) -> float:
    value = float(value)
    return 0.0 if value == 0.0 else value


def _bounds(bounds: Aabb) -> dict[str, float]:
    return {
        "x": _wire(bounds.min_x),
        "y": _wire(bounds.min_y),
        "width": _wire(bounds.width),
        "height": _wire(bounds.height),
    }


def _case(name: str, request: CompositionRequest) -> dict[str, object]:
    plan = plan_composition(request)
    drawing = Drawing(debug=False)
    primary = drawing.g(id=f"{name}-primary")
    drawing.add(primary)
    legend = drawing.g(id=f"{name}-legend") if plan.placement_for("legend") else None
    title = drawing.g(id=f"{name}-title") if plan.placement_for("title") else None
    if legend is not None:
        drawing.add(legend)
    if title is not None:
        drawing.add(title)
    apply_composition_plan(
        drawing,
        plan,
        primary_targets=(primary,),
        legend_target=legend,
        legend_side=request.legend_placement,
        legend_reflow_metrics=LEGEND_REFLOW if legend is not None else None,
        title_target=title,
        title_side=request.title_placement,
    )
    metadata = json.loads(drawing.attribs[COMPOSITION_METADATA_ATTRIBUTE])
    placements = {
        placement.role: {
            "automaticTranslation": [_wire(placement.dx), _wire(placement.dy)],
            "finalBounds": _bounds(placement.final_bounds),
        }
        for placement in plan.placements
    }
    return {
        "name": name,
        "schema": drawing.attribs[COMPOSITION_SCHEMA_ATTRIBUTE],
        "metadata": metadata,
        "expected": {
            "width": _wire(plan.width),
            "height": _wire(plan.height),
            "placements": placements,
            "overlayObstacles": [_bounds(bounds) for bounds in plan.overlay_obstacles],
            "overlayPolicy": metadata["overlayPolicy"],
            "spacing": metadata["spacing"],
        },
    }


def main() -> None:
    primary = CompositionItem("primary", Aabb(10.0, 20.0, 110.0, 100.0))
    legend = CompositionItem("legend", Aabb(-2.0, -5.0, 28.0, 15.0))
    title = CompositionItem("title", Aabb(-5.0, -10.0, 55.0, 10.0))
    cases = [
        _case(
            "dock",
            CompositionRequest(
                primary=primary,
                legend=legend,
                legend_placement=LegendPlacement.RIGHT,
            ),
        ),
        _case(
            "overlay",
            CompositionRequest(
                primary=CompositionItem("primary", Aabb(0.0, 0.0, 100.0, 100.0)),
                legend=CompositionItem("legend", Aabb(0.0, 0.0, 20.0, 20.0)),
                legend_placement=LegendPlacement.UPPER_LEFT,
                overlay_obstacles=(
                    Aabb(70.0, 70.0, 90.0, 90.0),
                    Aabb(0.0, 0.0, 30.0, 30.0),
                ),
            ),
        ),
        _case(
            "title_stack",
            CompositionRequest(
                primary=primary,
                legend=legend,
                title=title,
                legend_placement=LegendPlacement.TOP,
                title_placement=TitlePlacement.TOP,
            ),
        ),
        _case("no_legend", CompositionRequest(primary=primary)),
        _case(
            "canvas_growth",
            CompositionRequest(
                primary=CompositionItem("primary", Aabb(0.0, 0.0, 40.0, 40.0)),
                legend=CompositionItem("legend", Aabb(0.0, 0.0, 20.0, 20.0)),
                legend_placement=LegendPlacement.UPPER_LEFT,
                overlay_obstacles=(Aabb(0.0, 0.0, 40.0, 40.0),),
            ),
        ),
    ]
    print(json.dumps(cases, allow_nan=False, separators=(",", ":")))


if __name__ == "__main__":
    main()
