from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api import LinearMultiRecordOptions
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.layout.composition import LegendPlacement, TitlePlacement
from gbdraw.layout.spatial import Aabb, union_aabbs


_TRANSLATE_RE = re.compile(
    r"translate\(\s*([-+0-9.eE]+)[,\s]+([-+0-9.eE]+)\s*\)"
)


def _record(record_id: str = "record_a", *, length: int = 1_000) -> SeqRecord:
    record = SeqRecord(Seq("A" * length), id=record_id, description=record_id)
    record.features = [
        SeqFeature(FeatureLocation(100, min(700, length), strand=1), type="CDS")
    ]
    return record


def _config(*, track_layout: str = "middle"):
    return apply_config_overrides(
        None,
        {
            "canvas.show_gc": False,
            "canvas.show_skew": False,
            "canvas.linear.track_layout": track_layout,
            "labels.linear.scope": "none",
        },
    )


def _assemble(
    *,
    legend: str = "none",
    title: str | None = None,
    title_position: str = "bottom",
):
    return assemble_linear_diagram_from_records(
        [_record()],
        cfg=_config(),
        selected_features_set=["CDS"],
        legend=legend,
        plot_title=title,
        plot_title_position=title_position,
    )


def _center_x(bounds: Aabb) -> float:
    return 0.5 * (bounds.min_x + bounds.max_x)


def _center_y(bounds: Aabb) -> float:
    return 0.5 * (bounds.min_y + bounds.max_y)


def _translations(element: ET.Element) -> tuple[tuple[float, float], ...]:
    return tuple(
        (float(x_value), float(y_value))
        for x_value, y_value in _TRANSLATE_RE.findall(
            element.attrib.get("transform", "")
        )
    )


def _record_group(root: ET.Element, record_id: str) -> ET.Element:
    return next(
        element
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "g"
        and element.attrib.get("data-gbdraw-record-id") == record_id
    )


@pytest.mark.linear
@pytest.mark.parametrize(
    "side",
    (
        LegendPlacement.LEFT,
        LegendPlacement.RIGHT,
        LegendPlacement.TOP,
        LegendPlacement.BOTTOM,
    ),
)
def test_linear_composition_uses_exact_dock_gap_and_orthogonal_center(
    side: LegendPlacement,
) -> None:
    drawing = _assemble(legend=side.value)
    plan = drawing._gbdraw_linear_composition_plan
    primary = plan.primary_bounds
    legend = plan.placement_for("legend")

    assert legend is not None
    assert plan.spacing.dock_gap_px == pytest.approx(24.0)
    if side is LegendPlacement.LEFT:
        assert primary.min_x - legend.final_bounds.max_x == pytest.approx(24.0)
        assert _center_y(legend.final_bounds) == pytest.approx(_center_y(primary))
    elif side is LegendPlacement.RIGHT:
        assert legend.final_bounds.min_x - primary.max_x == pytest.approx(24.0)
        assert _center_y(legend.final_bounds) == pytest.approx(_center_y(primary))
    elif side is LegendPlacement.TOP:
        assert primary.min_y - legend.final_bounds.max_y == pytest.approx(24.0)
        assert _center_x(legend.final_bounds) == pytest.approx(_center_x(primary))
    else:
        assert legend.final_bounds.min_y - primary.max_y == pytest.approx(24.0)
        assert _center_x(legend.final_bounds) == pytest.approx(_center_x(primary))

    painted = union_aabbs(placement.final_bounds for placement in plan.placements)
    assert painted is not None
    assert (
        painted.min_x,
        painted.min_y,
        painted.max_x,
        painted.max_y,
    ) == pytest.approx(
        (16.0, 16.0, plan.width - 16.0, plan.height - 16.0),
        abs=1e-9,
    )


@pytest.mark.linear
def test_linear_composition_emits_schema_one_without_legacy_viewboxes() -> None:
    drawing = _assemble()
    root = ET.fromstring(drawing.tostring())
    metadata = json.loads(root.attrib["data-gbdraw-composition"])

    assert root.attrib["data-gbdraw-composition-schema"] == "1"
    assert "data-horizontal-viewbox" not in root.attrib
    assert "data-vertical-viewbox" not in root.attrib
    assert metadata["legend"] is None
    assert metadata["title"] is None
    assert metadata["legendSide"] == "none"
    assert metadata["titleSide"] == "none"
    assert metadata["spacing"] == {
        "dockGapPx": 24.0,
        "edgePaddingPx": 16.0,
        "overlayClearancePx": 8.0,
        "stackGapPx": 20.0,
        "titleGapPx": 20.0,
    }
    roles = [
        element.attrib.get("data-gbdraw-composition-role")
        for element in root
        if element.tag.rsplit("}", 1)[-1] == "g"
    ]
    assert roles.count("primary") >= 1
    assert "legend" not in roles
    assert "title" not in roles


@pytest.mark.linear
@pytest.mark.parametrize(
    "side",
    (TitlePlacement.TOP, TitlePlacement.BOTTOM, TitlePlacement.CENTER),
)
def test_linear_title_uses_shared_gap_or_center_overlay(
    side: TitlePlacement,
) -> None:
    drawing = _assemble(title="Linear composition", title_position=side.value)
    plan = drawing._gbdraw_linear_composition_plan
    primary = plan.primary_bounds
    title = plan.placement_for("title")

    assert title is not None
    assert _center_x(title.final_bounds) == pytest.approx(_center_x(primary))
    if side is TitlePlacement.TOP:
        assert primary.min_y - title.final_bounds.max_y == pytest.approx(20.0)
    elif side is TitlePlacement.BOTTOM:
        assert title.final_bounds.min_y - primary.max_y == pytest.approx(20.0)
    else:
        assert _center_y(title.final_bounds) == pytest.approx(_center_y(primary))


@pytest.mark.linear
@pytest.mark.parametrize("side", ("top", "bottom"))
def test_linear_same_side_legend_is_nearest_plot_and_title_is_outward(
    side: str,
) -> None:
    drawing = _assemble(
        legend=side,
        title="Linear composition",
        title_position=side,
    )
    plan = drawing._gbdraw_linear_composition_plan
    primary = plan.primary_bounds
    legend = plan.placement_for("legend")
    title = plan.placement_for("title")

    assert legend is not None
    assert title is not None
    if side == "top":
        assert primary.min_y - legend.final_bounds.max_y == pytest.approx(24.0)
        assert legend.final_bounds.min_y - title.final_bounds.max_y == pytest.approx(
            20.0
        )
    else:
        assert legend.final_bounds.min_y - primary.max_y == pytest.approx(24.0)
        assert title.final_bounds.min_y - legend.final_bounds.max_y == pytest.approx(
            20.0
        )


@pytest.mark.linear
def test_linear_legend_side_does_not_change_record_local_geometry() -> None:
    snapshots: list[tuple[Aabb, tuple[tuple[float, float], ...], tuple[str, ...]]] = []
    for side in ("none", "left", "right", "top", "bottom"):
        drawing = _assemble(legend=side)
        root = ET.fromstring(drawing.tostring())
        record_group = _record_group(root, "record_a")
        translations = _translations(record_group)
        paths = tuple(
            element.attrib["d"]
            for element in record_group.iter()
            if element.tag.rsplit("}", 1)[-1] == "path" and "d" in element.attrib
        )

        assert len(translations) >= 2
        snapshots.append(
            (
                drawing._gbdraw_linear_source_content_bounds,
                translations[1:],
                paths,
            )
        )

    assert all(snapshot == snapshots[0] for snapshot in snapshots[1:])


@pytest.mark.linear
def test_linear_multi_row_track_metadata_uses_final_primary_translation() -> None:
    records = [_record("record_a"), _record("record_b", length=800)]
    drawing = assemble_linear_diagram_from_records(
        records,
        cfg=_config(track_layout="above"),
        selected_features_set=["CDS"],
        layout=LinearMultiRecordOptions(
            multi_record_positions=("#1@1", "#2@2"),
        ),
        legend="top",
        plot_title="Two rows",
        plot_title_position="top",
    )
    plan = drawing._gbdraw_linear_composition_plan
    primary = plan.placement_for("primary")
    root = ET.fromstring(drawing.tostring())
    geometry = drawing._gbdraw_track_slot_geometry["records"]

    assert primary is not None
    assert plan.primary_bounds == drawing._gbdraw_linear_source_content_bounds.translated(
        primary.dx,
        primary.dy,
    )
    for record_payload in geometry:
        record_group = _record_group(root, record_payload["recordId"])
        translations = _translations(record_group)
        assert translations[0] == pytest.approx((primary.dx, primary.dy))
        assert sum(y_value for _x_value, y_value in translations) == pytest.approx(
            record_payload["axisYpx"]
        )
