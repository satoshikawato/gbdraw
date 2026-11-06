"""The Legend row of an annotation gradient: one bar from ``start`` to ``end``.

An annotation whose style has a gradient (a sequence-search hit drawn pale at
the query's start and dark at its end) gets one Legend row per caption and pair
of colours, drawn as the identity colour bar of Pairwise Comparisons with the
scale texts ``start`` and ``end``. Gradient rows that share scale texts share
one scale line; rows with other texts get a scale line of their own below.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

from gbdraw.annotations.legend import sync_annotation_legend_entries
from gbdraw.annotations.models import (
    AnnotationGradient,
    RegionAnnotationStyle,
    ResolvedAnnotationBundle,
    ResolvedRegionAnnotation,
)
from gbdraw.canvas import CircularCanvasConfigurator, LinearCanvasConfigurator
from gbdraw.config.models import CircularRenderProfile, GbdrawConfig, LinearRenderProfile
from gbdraw.config.toml import load_config_toml
from gbdraw.configurators import LegendDrawingConfigurator, LegendMeasurement
from gbdraw.configurators.legend import _circular_legend_local_bounds, measure_legend
from gbdraw.core.text import drawn_text_box
from gbdraw.layout.spatial import Aabb, union_aabbs
from gbdraw.legend.metrics import gradient_scale_texts
from gbdraw.legend.sources import LegendSlotSource
from gbdraw.legend.table import LegendTable
from gbdraw.render.groups.circular.legend import LegendGroup as CircularLegendGroup
from gbdraw.render.groups.linear.legend import LegendGroup as LinearLegendGroup

SVG = "{http://www.w3.org/2000/svg}"

PAIRWISE_ROW = {
    "type": "gradient",
    "min_color": "#ffffff",
    "max_color": "#000000",
    "stroke": "none",
    "width": 0,
    "min_value": 70,
    "max_value": 100,
}
ANNOTATION_ROW = {
    "type": "gradient",
    "min_color": "#dbeafe",
    "max_color": "#1e3a8a",
    "stroke": "none",
    "width": 0,
    "min_label": "start",
    "max_label": "end",
}
RING_ROW = {
    "type": "ring",
    "source_index": 0,
    "track_index": 1,
    "fill": None,
    "stroke": "none",
    "width": 0,
    "min_color": "#f0f1f5",
    "max_color": "#8b9cc1",
    "min_value": 70,
    "max_value": 100,
}


def _config() -> GbdrawConfig:
    return GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))


def _record() -> SeqRecord:
    record = SeqRecord(Seq("ATGC" * 250), id="legend-test")
    record.annotations["molecule_type"] = "DNA"
    return record


def _linear_measurement(table: LegendTable | dict, legend: str = "bottom") -> LegendMeasurement:
    profile = LinearRenderProfile(_config())
    canvas = LinearCanvasConfigurator(num_of_entries=1, longest_genome=1_000, profile=profile, legend=legend)
    return measure_legend(
        LegendDrawingConfigurator(profile=profile, canvas_config=canvas), table,
        placement=canvas.legend_position, wrap_width=canvas.total_width,
    )


def _circular_measurement(table: LegendTable | dict, legend: str = "right") -> LegendMeasurement:
    profile = CircularRenderProfile(_config())
    canvas = CircularCanvasConfigurator(output_prefix="out", profile=profile, legend=legend, gb_record=_record())
    return measure_legend(
        LegendDrawingConfigurator(profile=profile, canvas_config=canvas), table,
        placement=canvas.legend_position, wrap_width=canvas.total_width,
    )


def _svg(group: object) -> ET.Element:
    drawing = Drawing(debug=False)
    drawing.add(group)  # type: ignore[arg-type]
    return ET.fromstring(drawing.tostring())


def _texts(element: ET.Element) -> list[str]:
    return ["".join(text.itertext()) for text in element.iter(f"{SVG}text")]


def _translate(element: ET.Element) -> tuple[float, float]:
    transform = element.get("transform") or ""
    if not transform.startswith("translate("):
        return 0.0, 0.0
    x, _, y = transform[len("translate("):-1].partition(",")
    return float(x), float(y or 0)


def _child(parent: ET.Element, element_id: str) -> ET.Element:
    found = parent.find(f".//{SVG}g[@id='{element_id}']")
    assert found is not None, element_id
    return found


def test_scale_texts_come_from_one_owner() -> None:
    assert gradient_scale_texts(PAIRWISE_ROW) == ("70%", "100%")
    assert gradient_scale_texts({**PAIRWISE_ROW, "min_value": 72.5}) == ("72.5%", "100%")
    assert gradient_scale_texts(ANNOTATION_ROW) == ("start", "end")


@pytest.mark.parametrize("legend", ["bottom", "right"])
def test_linear_pairwise_and_annotation_gradient_rows_draw_two_scale_lines(legend: str) -> None:
    table = {"Pairwise match identity": PAIRWISE_ROW, "Query position": ANNOTATION_ROW}
    measurement = _linear_measurement(table, legend)
    layout = measurement.linear_layout
    assert layout is not None
    root = _svg(LinearLegendGroup(measurement, LegendTable(table)).get_group())
    rect = measurement.color_rect_size
    for orientation, suffix in ((layout.horizontal, "h"), (layout.vertical, "v")):
        gradient = orientation.gradient
        assert gradient is not None and gradient.compact
        # The first group's scale keeps the existing fields; the second one follows its own row.
        assert [entry.key for entry in gradient.entries] == ["Pairwise match identity", "Query position"]
        assert (gradient.min_label_text, gradient.max_label_text) == ("70%", "100%")
        assert [(scale.min_label_text, scale.max_label_text) for scale in gradient.more_scales] == [("start", "end")]
        second_row_top = gradient.entries[1].bar_y - (rect / 2)
        assert second_row_top > gradient.scale_label_y
        assert gradient.more_scales[0].scale_label_y > gradient.entries[1].bar_y
        assert gradient.height > gradient.more_scales[0].scale_label_y
        group = _child(root, f"pairwise_legend_{suffix}")
        assert _texts(group) == ["Pairwise match identity", "Query position", "70%", "100%", "start", "end"]


def _hit(annotation_id: str, direction: str) -> ResolvedRegionAnnotation:
    style = RegionAnnotationStyle(gradient=AnnotationGradient("#dbeafe", "#1e3a8a", direction))  # type: ignore[arg-type]
    return ResolvedRegionAnnotation(
        id=annotation_id, set_id="hits", record_index=0, segments=((10, 40),), midpoint_bp=25.0, span_bp=30,
        label="", mark="highlight", lane=None, style=style, legend_label="Query position",
    )


def test_linear_annotation_gradient_rows_of_one_caption_and_colours_are_one_legend_row() -> None:
    bundle = ResolvedAnnotationBundle(annotations=(_hit("hit1", "forward"), _hit("hit2", "reverse")))
    table = sync_annotation_legend_entries(
        LegendTable(), bundle, [LegendSlotSource("search", "annotations", {"set_id": "hits"})],
    )
    assert list(table) == ["Query position"]
    measurement = _linear_measurement(table)
    root = _svg(LinearLegendGroup(measurement, table).get_group())
    group = _child(root, "pairwise_legend_h")
    rows = group.findall(f"{SVG}g[@data-legend-key]")
    assert [row.get("data-legend-key") for row in rows] == ["Query position"]
    assert _texts(group) == ["Query position", "start", "end"]
    assert len(group.findall(f".//{SVG}linearGradient")) == 1


def test_circular_annotation_gradient_row_is_laid_out_and_drawn() -> None:
    table = {"CDS": {"type": "solid", "fill": "#54bcf8", "stroke": "none", "width": 0}, "Query position": ANNOTATION_ROW}
    measurement = _circular_measurement(table)
    layout = measurement.circular_layout
    assert layout is not None and layout.gradient is None and layout.gradient_rows is not None
    assert layout.has_gradient
    rows = layout.gradient_rows
    assert (rows.min_label_text, rows.max_label_text) == ("start", "end")
    root = _svg(CircularLegendGroup(measurement, LegendTable(table)).get_group())
    block = _child(root, "conservation_identity_legend")
    assert _texts(block) == ["Query position", "start", "end"]
    bars = [path for path in block.iter(f"{SVG}path") if (path.get("fill") or "").startswith("url(#")]
    assert len(bars) == 1
    # The drawn row names its scale's ends for the Web's Legend layout; one row has no trailing gap (D-44).
    entry_group = next(group for group in block.iter(f"{SVG}g") if group.get("data-legend-key") is not None)
    assert (entry_group.get("data-legend-min-label"), entry_group.get("data-legend-max-label")) == ("start", "end")
    label_height = drawn_text_box("100%", measurement.font_family, measurement.font_size)[1]
    assert rows.height == (rows.scale_label_y - rows.single_entries[0].title_y) + label_height
    # The bar lies inside the measured bounds.
    gx, gy = layout.gradient_x, layout.gradient_y
    entry = rows.single_entries[0]
    bar = Aabb(
        gx + rows.bar_x, gy + entry.bar_y - measurement.color_rect_size / 2,
        gx + rows.bar_x + rows.bar_width, gy + entry.bar_y + measurement.color_rect_size / 2,
    )
    assert union_aabbs([measurement.local_bounds, bar]) == measurement.local_bounds


@pytest.mark.parametrize("legend", ["right", "bottom"])
def test_circular_ring_block_and_annotation_gradient_rows_are_both_drawn(legend: str) -> None:
    table = {
        "CDS": {"type": "solid", "fill": "#54bcf8", "stroke": "none", "width": 0},
        "Reference A": RING_ROW,
        "Query position": ANNOTATION_ROW,
    }
    measurement = _circular_measurement(table, legend)
    layout = measurement.circular_layout
    assert layout is not None and layout.gradient is not None and layout.gradient_rows is not None
    ring, rows = layout.gradient, layout.gradient_rows
    rect = measurement.color_rect_size
    # The gradient rows sit below the ring block, half a legend line under it, in one side block.
    rows_top = rows.single_entries[0].title_y
    assert rows_top == pytest.approx(ring.height + layout.line_margin / 2)
    assert layout.width >= max(ring.width, rows.width)
    side_bottom = layout.gradient_y + rows_top + rows.height
    assert layout.height >= side_bottom
    root = _svg(CircularLegendGroup(measurement, LegendTable(table)).get_group())
    block = _child(root, "conservation_identity_legend")
    assert _texts(block) == ["outer", "Reference A", "inner", "identity", "70%", "100%", "Query position", "start", "end"]
    assert len(root.findall(f".//{SVG}g[@id='conservation_identity_legend']")) == 1
    # No overlap: the gradient row's title starts below the ring block's scale labels.
    assert rows_top > ring.scale_label_y
    bounds = measurement.local_bounds
    gx, gy = layout.gradient_x, layout.gradient_y
    band = Aabb(gx, gy + ring.band_y - rect / 2, gx + ring.bar_width, gy + ring.band_y + rect / 2)
    entry = rows.single_entries[0]
    bar = Aabb(
        gx + rows.bar_x, gy + entry.bar_y - rect / 2, gx + rows.bar_x + rows.bar_width, gy + entry.bar_y + rect / 2,
    )
    assert union_aabbs([bounds, band]) == bounds
    assert union_aabbs([bounds, bar]) == bounds
    assert bounds == _circular_legend_local_bounds(layout, color_rect_size=rect)


def test_circular_gradient_rows_of_two_scales_draw_two_scale_lines() -> None:
    table = {"Pairwise match identity": PAIRWISE_ROW, "Query position": ANNOTATION_ROW}
    measurement = _circular_measurement(table)
    rows = measurement.circular_layout.gradient_rows  # type: ignore[union-attr]
    assert rows is not None and rows.compact
    assert [entry.key for entry in rows.compact_entries] == ["Pairwise match identity", "Query position"]
    assert [(scale.min_label_text, scale.max_label_text) for scale in rows.more_scales] == [("start", "end")]
    assert rows.compact_entries[1].bar_y - measurement.color_rect_size / 2 > rows.scale_label_y
    root = _svg(CircularLegendGroup(measurement, LegendTable(table)).get_group())
    block = _child(root, "conservation_identity_legend")
    assert _texts(block) == ["Pairwise match identity", "Query position", "70%", "100%", "start", "end"]
