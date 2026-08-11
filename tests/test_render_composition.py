from __future__ import annotations

import json
import math
import re
import xml.etree.ElementTree as ET

import pytest
from svgwrite import Drawing
from svgwrite.container import Group

from gbdraw.layout.composition import (
    CompositionItem,
    CompositionRequest,
    LegendPlacement,
    TitlePlacement,
    plan_composition,
)
from gbdraw.layout.spatial import Aabb
from gbdraw.render.composition import (
    COMPOSITION_METADATA_ATTRIBUTE,
    COMPOSITION_ROLE_ATTRIBUTE,
    COMPOSITION_SCHEMA_ATTRIBUTE,
    apply_composition_plan,
)

_URL_REFERENCE_RE = re.compile(r"url\(#([^)]+)\)")
_TRANSFORM_RE = re.compile(r"(matrix|rotate|scale|translate)\(([^)]*)\)")
_LEGEND_REFLOW_METRICS = {
    "colorRectSize": 14.0,
    "lineHeight": 24.0,
    "textXOffset": 22.0,
}


def _multiply_affine(
    left: tuple[float, float, float, float, float, float],
    right: tuple[float, float, float, float, float, float],
) -> tuple[float, float, float, float, float, float]:
    a1, b1, c1, d1, e1, f1 = left
    a2, b2, c2, d2, e2, f2 = right
    return (
        a1 * a2 + c1 * b2,
        b1 * a2 + d1 * b2,
        a1 * c2 + c1 * d2,
        b1 * c2 + d1 * d2,
        a1 * e2 + c1 * f2 + e1,
        b1 * e2 + d1 * f2 + f1,
    )


def _transform_point(transform: str, point: tuple[float, float]) -> tuple[float, float]:
    combined = (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)
    for name, raw_arguments in _TRANSFORM_RE.findall(transform):
        values = tuple(
            float(value) for value in re.split(r"[ ,]+", raw_arguments.strip()) if value
        )
        if name == "matrix":
            matrix = values
        elif name == "translate":
            matrix = (
                1.0,
                0.0,
                0.0,
                1.0,
                values[0],
                values[1] if len(values) > 1 else 0.0,
            )
        elif name == "scale":
            matrix = (
                values[0],
                0.0,
                0.0,
                values[1] if len(values) > 1 else values[0],
                0.0,
                0.0,
            )
        else:
            radians = math.radians(values[0])
            matrix = (
                math.cos(radians),
                math.sin(radians),
                -math.sin(radians),
                math.cos(radians),
                0.0,
                0.0,
            )
        combined = _multiply_affine(combined, matrix)  # type: ignore[arg-type]
    a, b, c, d, e, f = combined
    x, y = point
    return a * x + c * y + e, b * x + d * y + f


def _complete_plan():
    return plan_composition(
        CompositionRequest(
            primary=CompositionItem("primary", Aabb(10, 20, 110, 100)),
            legend=CompositionItem("legend", Aabb(-2, -5, 28, 15)),
            title=CompositionItem("title", Aabb(-5, -10, 55, 10)),
            legend_placement=LegendPlacement.RIGHT,
            title_placement=TitlePlacement.TOP,
        )
    )


def _drawing_with_targets():
    drawing = Drawing(
        size=("999px", "888px"),
        viewBox="4 5 999 888",
        debug=False,
    )
    drawing.attribs.update(
        {
            "class": "figure-root",
            "data-owner": "retained",
            "id": "diagram-root",
        }
    )
    first = drawing.g(
        id="axis",
        transform=("translate(1,2) scale(2) rotate(7) matrix(1 0 0 1 3 4)"),
    )
    first.attribs["data-gbdraw-role"] = "supported-axis-role"
    second = drawing.g(id="labels", transform="scale(0.5)")
    legend = drawing.g(id="legend", transform="rotate(3)")
    title = drawing.g(id="title", transform="matrix(1 0 0 1 0 0)")
    for target in (first, second, legend, title):
        drawing.add(target)
    return drawing, first, second, legend, title


def _ids_and_references(svg: str) -> tuple[set[str], set[str]]:
    root = ET.fromstring(svg)
    identifiers = {
        element.attrib["id"] for element in root.iter() if element.attrib.get("id")
    }
    references: set[str] = set()
    for element in root.iter():
        for name, value in element.attrib.items():
            references.update(_URL_REFERENCE_RE.findall(value))
            if name.rsplit("}", 1)[-1] == "href" and value.startswith("#"):
                references.add(value[1:])
    return identifiers, references


def test_apply_plan_composes_outer_transforms_and_syncs_root_once() -> None:
    plan = _complete_plan()
    drawing, first, second, legend, title = _drawing_with_targets()
    old_first_transform = first.attribs["transform"]

    apply_composition_plan(
        drawing,
        plan,
        primary_targets=(first, second),
        legend_target=legend,
        legend_side="right",
        legend_reflow_metrics=_LEGEND_REFLOW_METRICS,
        title_target=title,
        title_side="top",
    )

    assert first.attribs["transform"] == (
        "translate(6.0,36.0) translate(1,2) scale(2) rotate(7) matrix(1 0 0 1 3 4)"
    )
    assert second.attribs["transform"] == "translate(6.0,36.0) scale(0.5)"
    assert legend.attribs["transform"] == "translate(142.0,91.0) rotate(3)"
    assert title.attribs["transform"] == ("translate(41.0,26.0) matrix(1 0 0 1 0 0)")
    old_point = _transform_point(old_first_transform, (11.0, 13.0))
    new_point = _transform_point(first.attribs["transform"], (11.0, 13.0))
    assert new_point == pytest.approx((old_point[0] + 6.0, old_point[1] + 36.0))
    assert first.attribs["data-gbdraw-role"] == "supported-axis-role"
    assert [
        target.attribs[COMPOSITION_ROLE_ATTRIBUTE]
        for target in (first, second, legend, title)
    ] == ["primary", "primary", "legend", "title"]

    assert drawing.attribs["width"] == "186.0px"
    assert drawing.attribs["height"] == "152.0px"
    assert drawing.attribs["viewBox"] == "0.0 0.0 186.0 152.0"
    assert drawing.attribs["class"] == "figure-root"
    assert drawing.attribs["data-owner"] == "retained"
    assert drawing.attribs["id"] == "diagram-root"


def test_composition_metadata_schema_one_is_exact_and_deterministic() -> None:
    expected = (
        '{"legend":{"automaticTranslation":[142.0,91.0],"localBounds":'
        '{"height":20.0,"width":30.0,"x":-2.0,"y":-5.0},"role":"legend",'
        '"selector":"[data-gbdraw-composition-role=\\"legend\\"]"},'
        '"legendReflow":{"colorRectSize":14.0,"lineHeight":24.0,'
        '"textXOffset":22.0},"legendSide":"right","overlayObstacles":[],'
        '"overlayPolicy":{"candidateScoreOrder":["totalAnchorDistance",'
        '"xAnchorDistance","yAnchorDistance","nearEdgeX","nearEdgeY"],'
        '"canvasGrowthCandidateOrder":["horizontal","vertical"],'
        '"canvasGrowthScoreOrder":["addedArea","addedExtent",'
        '"candidateOrder"],"quadrantBoundaryRatio":0.5},'
        '"primary":'
        '{"automaticTranslation":[6.0,36.0],'
        '"finalBounds":{"height":80.0,"width":100.0,"x":16.0,"y":56.0},'
        '"role":"primary","selector":"[data-gbdraw-composition-role='
        '\\"primary\\"]"},"spacing":{"dockGapPx":24.0,"edgePaddingPx":16.0,'
        '"overlayClearancePx":8.0,"stackGapPx":20.0,"titleGapPx":20.0},'
        '"title":{"automaticTranslation":[41.0,26.0],"localBounds":'
        '{"height":20.0,"width":60.0,"x":-5.0,"y":-10.0},"role":"title",'
        '"selector":"[data-gbdraw-composition-role=\\"title\\"]"},'
        '"titleSide":"top"}'
    )
    emitted: list[str] = []
    for _ in range(2):
        drawing, first, second, legend, title = _drawing_with_targets()
        apply_composition_plan(
            drawing,
            _complete_plan(),
            primary_targets=(first, second),
            legend_target=legend,
            legend_side=LegendPlacement.RIGHT,
            legend_reflow_metrics=_LEGEND_REFLOW_METRICS,
            title_target=title,
            title_side=TitlePlacement.TOP,
        )
        assert drawing.attribs[COMPOSITION_SCHEMA_ATTRIBUTE] == "1"
        emitted.append(drawing.attribs[COMPOSITION_METADATA_ATTRIBUTE])
        parsed_root = ET.fromstring(drawing.tostring())
        assert parsed_root.attrib[COMPOSITION_METADATA_ATTRIBUTE] == expected

    assert emitted == [expected, expected]
    payload = json.loads(expected)
    assert payload["primary"]["finalBounds"] == {
        "height": 80.0,
        "width": 100.0,
        "x": 16.0,
        "y": 56.0,
    }
    assert payload["primary"]["selector"] == (
        '[data-gbdraw-composition-role="primary"]'
    )
    assert payload["overlayObstacles"] == []
    assert payload["overlayPolicy"] == {
        "candidateScoreOrder": [
            "totalAnchorDistance",
            "xAnchorDistance",
            "yAnchorDistance",
            "nearEdgeX",
            "nearEdgeY",
        ],
        "canvasGrowthCandidateOrder": ["horizontal", "vertical"],
        "canvasGrowthScoreOrder": [
            "addedArea",
            "addedExtent",
            "candidateOrder",
        ],
        "quadrantBoundaryRatio": 0.5,
    }


def test_apply_plan_preserves_defs_ids_and_local_references() -> None:
    drawing, first, second, legend, title = _drawing_with_targets()
    pattern = drawing.pattern(id="stripe-pattern")
    pattern.add(drawing.rect(insert=(0, 0), size=(4, 4), fill="#ffffff"))
    drawing.defs.add(pattern)
    source = drawing.circle(id="source-shape", center=(2, 2), r=2)
    drawing.defs.add(source)
    first.add(drawing.rect(id="painted", size=(10, 10), fill="url(#stripe-pattern)"))
    second.add(drawing.use("#source-shape", id="source-use"))

    before_svg = drawing.tostring()
    before_defs = drawing.defs.tostring()
    before_ids, before_references = _ids_and_references(before_svg)

    apply_composition_plan(
        drawing,
        _complete_plan(),
        primary_targets=(first, second),
        legend_target=legend,
        legend_side="right",
        legend_reflow_metrics=_LEGEND_REFLOW_METRICS,
        title_target=title,
        title_side="top",
    )

    after_svg = drawing.tostring()
    after_ids, after_references = _ids_and_references(after_svg)
    assert drawing.defs.tostring() == before_defs
    assert after_ids == before_ids
    assert (
        after_references
        == before_references
        == {
            "source-shape",
            "stripe-pattern",
        }
    )
    assert after_references <= after_ids


def test_absent_decorations_are_explicit_in_metadata() -> None:
    plan = plan_composition(
        CompositionRequest(
            primary=CompositionItem("primary", Aabb(10, 20, 110, 100)),
            legend_placement=LegendPlacement.RIGHT,
            title_placement=TitlePlacement.BOTTOM,
        )
    )
    drawing = Drawing(debug=False)
    primary = drawing.g(id="primary")
    drawing.add(primary)

    apply_composition_plan(
        drawing,
        plan,
        primary_targets=(primary,),
        legend_side="right",
        title_side="bottom",
    )

    payload = json.loads(drawing.attribs[COMPOSITION_METADATA_ATTRIBUTE])
    assert payload["legend"] is None
    assert payload["legendReflow"] is None
    assert payload["legendSide"] == "right"
    assert payload["title"] is None
    assert payload["titleSide"] == "bottom"
    assert primary.attribs[COMPOSITION_ROLE_ATTRIBUTE] == "primary"
    assert drawing.attribs["width"] == "132.0px"
    assert drawing.attribs["height"] == "112.0px"
    assert drawing.attribs["viewBox"] == "0.0 0.0 132.0 112.0"


def test_directly_constructed_target_accepts_internal_role_metadata() -> None:
    plan = plan_composition(
        CompositionRequest(
            primary=CompositionItem("primary", Aabb(0, 0, 10, 10)),
        )
    )
    drawing = Drawing(debug=False)
    primary = Group(id="direct-primary")
    child = Group(id="validated-child")
    primary.add(child)
    drawing.add(primary)

    assert primary.debug is True
    assert child.debug is True

    apply_composition_plan(drawing, plan, primary_targets=(primary,))

    assert primary.debug is False
    assert child.debug is True
    parsed = ET.fromstring(drawing.tostring())
    emitted = next(node for node in parsed.iter() if node.attrib.get("id") == "direct-primary")
    assert emitted.attrib[COMPOSITION_ROLE_ATTRIBUTE] == "primary"


def test_reapplication_fails_without_mutating_applied_state() -> None:
    drawing, first, second, legend, title = _drawing_with_targets()
    plan = _complete_plan()
    arguments = {
        "primary_targets": (first, second),
        "legend_target": legend,
        "legend_side": "right",
        "legend_reflow_metrics": _LEGEND_REFLOW_METRICS,
        "title_target": title,
        "title_side": "top",
    }
    apply_composition_plan(drawing, plan, **arguments)
    serialized = drawing.tostring()

    with pytest.raises(RuntimeError, match="already been applied"):
        apply_composition_plan(drawing, plan, **arguments)

    assert drawing.tostring() == serialized


def test_binding_error_is_reported_before_any_svg_mutation() -> None:
    drawing, first, second, _legend, _title = _drawing_with_targets()
    before = drawing.tostring()

    with pytest.raises(ValueError, match="legend target and plan placement"):
        apply_composition_plan(
            drawing,
            _complete_plan(),
            primary_targets=(first, second),
            legend_side="right",
            title_side="top",
        )

    assert drawing.tostring() == before


@pytest.mark.parametrize(
    ("metrics", "message"),
    [
        (None, "legend target requires legend reflow metrics"),
        ({"colorRectSize": 14.0}, "invalid field set"),
        (
            {"colorRectSize": 14.0, "lineHeight": float("nan"), "textXOffset": 22.0},
            "finite and positive",
        ),
        (
            {"colorRectSize": 14.0, "lineHeight": "invalid", "textXOffset": 22.0},
            "finite and positive",
        ),
        (
            {"colorRectSize": 14.0, "lineHeight": "24", "textXOffset": 22.0},
            "finite and positive",
        ),
        (
            {"colorRectSize": 14.0, "lineHeight": True, "textXOffset": 22.0},
            "finite and positive",
        ),
    ],
)
def test_legend_reflow_metrics_are_required_and_strict(
    metrics: dict[str, float] | None,
    message: str,
) -> None:
    drawing, first, second, legend, title = _drawing_with_targets()
    before = drawing.tostring()

    with pytest.raises(ValueError, match=message):
        apply_composition_plan(
            drawing,
            _complete_plan(),
            primary_targets=(first, second),
            legend_target=legend,
            legend_side="right",
            legend_reflow_metrics=metrics,
            title_target=title,
            title_side="top",
        )

    assert drawing.tostring() == before
