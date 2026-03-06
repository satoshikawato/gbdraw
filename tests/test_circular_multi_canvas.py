from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.circular as circular_cli_module
import gbdraw.api.diagram as diagram_api_module
import gbdraw.diagrams.circular.assemble as circular_assemble_module
from gbdraw.api.diagram import assemble_circular_diagram_from_records
from gbdraw.api.options import DiagramOptions, OutputOptions
from gbdraw.core.text import calculate_bbox_dimensions
from gbdraw.features.colors import compute_feature_hash


def _build_record(record_id: str, feature_start: int, length: int = 1200) -> SeqRecord:
    record = SeqRecord(Seq("A" * length), id=record_id)
    record.features = [
        SeqFeature(
            FeatureLocation(feature_start, feature_start + 120, strand=1),
            type="CDS",
            qualifiers={"product": [f"protein_{record_id}"]},
        )
    ]
    return record


def _build_record_with_source(
    record_id: str,
    *,
    organism: str,
    strain: str,
    feature_start: int = 30,
    length: int = 1200,
) -> SeqRecord:
    record = SeqRecord(Seq("A" * length), id=record_id)
    record.features = [
        SeqFeature(
            FeatureLocation(0, length, strand=1),
            type="source",
            qualifiers={
                "organism": [organism],
                "strain": [strain],
            },
        ),
        SeqFeature(
            FeatureLocation(feature_start, feature_start + 120, strand=1),
            type="CDS",
            qualifiers={"product": [f"protein_{record_id}"]},
        ),
    ]
    return record


def _build_dense_record_with_source(
    record_id: str,
    *,
    organism: str,
    strain: str,
    length: int = 2400,
) -> SeqRecord:
    record = SeqRecord(Seq("A" * length), id=record_id)
    record.features = [
        SeqFeature(
            FeatureLocation(0, length, strand=1),
            type="source",
            qualifiers={
                "organism": [organism],
                "strain": [strain],
            },
        ),
        SeqFeature(
            FeatureLocation(120, 320, strand=1),
            type="CDS",
            qualifiers={"product": [f"protein_{record_id}_a"]},
        ),
        SeqFeature(
            FeatureLocation(780, 930, strand=-1),
            type="CDS",
            qualifiers={"product": [f"protein_{record_id}_b"]},
        ),
        SeqFeature(
            FeatureLocation(1180, 1320, strand=1),
            type="tRNA",
            qualifiers={"product": [f"trna_{record_id}"]},
        ),
        SeqFeature(
            FeatureLocation(1710, 1960, strand=-1),
            type="rRNA",
            qualifiers={"product": [f"rrna_{record_id}"]},
        ),
    ]
    return record


def _build_multi_feature_record_with_source(
    record_id: str,
    *,
    organism: str,
    strain: str,
    length: int = 2600,
) -> SeqRecord:
    record = SeqRecord(Seq("A" * length), id=record_id)
    record.features = [
        SeqFeature(
            FeatureLocation(0, length, strand=1),
            type="source",
            qualifiers={
                "organism": [organism],
                "strain": [strain],
            },
        ),
        SeqFeature(
            FeatureLocation(80, 260, strand=1),
            type="CDS",
            qualifiers={"product": [f"cds_{record_id}"]},
        ),
        SeqFeature(
            FeatureLocation(340, 480, strand=1),
            type="tRNA",
            qualifiers={"product": [f"trna_{record_id}"]},
        ),
        SeqFeature(
            FeatureLocation(560, 760, strand=1),
            type="rRNA",
            qualifiers={"product": [f"rrna_{record_id}"]},
        ),
        SeqFeature(
            FeatureLocation(840, 980, strand=-1),
            type="tmRNA",
            qualifiers={"product": [f"tmrna_{record_id}"]},
        ),
        SeqFeature(
            FeatureLocation(1060, 1180, strand=1),
            type="ncRNA",
            qualifiers={"product": [f"ncrna_{record_id}"]},
        ),
        SeqFeature(
            FeatureLocation(1260, 1410, strand=-1),
            type="misc_RNA",
            qualifiers={"product": [f"misc_{record_id}"]},
        ),
        SeqFeature(
            FeatureLocation(1490, 1670, strand=1),
            type="repeat_region",
            qualifiers={"rpt_type": [f"repeat_{record_id}"]},
        ),
    ]
    return record


def _extract_group_texts(root: ET.Element, group_id: str) -> list[str]:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    group = root.find(f".//svg:g[@id='{group_id}']", ns)
    assert group is not None
    return [
        "".join(text.itertext()).strip()
        for text in group.findall("./svg:text", ns)
        if "".join(text.itertext()).strip()
    ]


def _extract_group_font_sizes(root: ET.Element, group_id: str) -> list[float]:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    group = root.find(f".//svg:g[@id='{group_id}']", ns)
    assert group is not None
    sizes: list[float] = []
    for text in group.findall("./svg:text", ns):
        if not "".join(text.itertext()).strip():
            continue
        font_size = text.attrib.get("font-size")
        assert font_size is not None
        sizes.append(float(font_size))
    return sizes


def _extract_group_font_weights(root: ET.Element, group_id: str) -> list[str]:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    group = root.find(f".//svg:g[@id='{group_id}']", ns)
    assert group is not None
    weights: list[str] = []
    for text in group.findall("./svg:text", ns):
        if not "".join(text.itertext()).strip():
            continue
        weights.append(text.attrib.get("font-weight", "normal"))
    return weights


def _extract_group_translate_y(root: ET.Element, group_id: str) -> float:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    group = root.find(f".//svg:g[@id='{group_id}']", ns)
    assert group is not None
    transform = group.attrib.get("transform", "")
    match = re.search(
        r"translate\(\s*[-+0-9.eE]+\s*,\s*([-+0-9.eE]+)\s*\)",
        transform,
    )
    assert match is not None
    return float(match.group(1))


def _extract_group_translate_xy(root: ET.Element, group_id: str) -> tuple[float, float]:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    group = root.find(f".//svg:g[@id='{group_id}']", ns)
    assert group is not None
    transform = group.attrib.get("transform", "")
    match = re.search(
        r"translate\(\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\)",
        transform,
    )
    assert match is not None
    return float(match.group(1)), float(match.group(2))


def _extract_viewbox_height(root: ET.Element) -> float:
    view_box = root.attrib.get("viewBox", "")
    parts = [float(part) for part in view_box.split()]
    assert len(parts) == 4
    return parts[3]


def _extract_viewbox_width(root: ET.Element) -> float:
    view_box = root.attrib.get("viewBox", "")
    parts = [float(part) for part in view_box.split()]
    assert len(parts) == 4
    return parts[2]


def _build_mock_circular_subcanvas(width: float, height: float, inset: float) -> Drawing:
    canvas = Drawing(
        filename="mock.svg",
        size=(f"{width}px", f"{height}px"),
        viewBox=f"0 0 {width} {height}",
        debug=False,
    )
    axis_group = canvas.g(id="Axis")
    axis_group.translate(width * 0.5, height * 0.5)
    radius = max(1.0, (width * 0.5) - float(inset))
    axis_group.add(canvas.circle(center=(0, 0), r=radius, fill="none", stroke="gray"))
    canvas.add(axis_group)
    return canvas


def _extract_record_axis_outer_edges(root: ET.Element, record_index: int) -> tuple[float, float]:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    record_group = root.find(f".//svg:g[@id='record_{record_index}']", ns)
    assert record_group is not None
    record_x, _record_y = _parse_translate(record_group.attrib.get("transform", ""))

    axis_group = record_group.find(f"./svg:g[@id='Axis_{record_index}']", ns)
    assert axis_group is not None
    axis_x, _axis_y = _parse_translate(axis_group.attrib.get("transform", ""))

    circle = axis_group.find("./svg:circle", ns)
    assert circle is not None
    radius = float(circle.attrib.get("r", "0"))

    center_x = record_x + axis_x
    return center_x - radius, center_x + radius


def _extract_legend_text_transforms(root: ET.Element) -> list[tuple[float, float, str]]:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    legend = root.find(".//svg:g[@id='legend']", ns)
    assert legend is not None
    text_positions: list[tuple[float, float, str]] = []
    for text in legend.findall(".//svg:text", ns):
        content = "".join(text.itertext()).strip()
        if not content:
            continue
        transform = text.attrib.get("transform", "")
        match = re.search(
            r"translate\(\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\)",
            transform,
        )
        assert match is not None
        text_positions.append((float(match.group(1)), float(match.group(2)), content))
    return text_positions


def _extract_horizontal_legend_row_centers(root: ET.Element) -> list[tuple[float, float]]:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    legend = root.find(".//svg:g[@id='legend']", ns)
    assert legend is not None
    legend_x, _legend_y = _parse_translate(legend.attrib.get("transform", ""))

    row_bounds: dict[float, tuple[float, float]] = {}
    for entry in legend.findall("./svg:g[@data-legend-key]", ns):
        text = entry.find("./svg:text", ns)
        if text is None:
            continue

        rect: ET.Element | None = None
        for candidate in entry.findall("./svg:path", ns):
            fill = candidate.attrib.get("fill", "")
            if fill and fill != "none" and not fill.startswith("url("):
                rect = candidate
                break
        if rect is None:
            continue

        text_x, text_y = _parse_translate(text.attrib.get("transform", ""))
        rect_x, _rect_y = _parse_translate(rect.attrib.get("transform", ""))
        x_margin = text_x - rect_x

        caption = "".join(text.itertext()).strip()
        font_size = float(text.attrib.get("font-size", "0").replace("px", ""))
        font_family = text.attrib.get("font-family", "")
        text_width, _ = calculate_bbox_dimensions(caption, font_family, font_size, 72)

        entry_left = legend_x + rect_x
        entry_right = entry_left + float(text_width) + 2.0 * float(x_margin)
        row_key = round(text_y, 3)
        if row_key not in row_bounds:
            row_bounds[row_key] = (entry_left, entry_right)
        else:
            row_min, row_max = row_bounds[row_key]
            row_bounds[row_key] = (min(row_min, entry_left), max(row_max, entry_right))

    centers: list[tuple[float, float]] = []
    for row_key, (row_min, row_max) in row_bounds.items():
        centers.append((row_key, (row_min + row_max) * 0.5))
    centers.sort(key=lambda item: item[0])
    return centers


def _parse_translate(transform: str) -> tuple[float, float]:
    match = re.search(
        r"translate\(\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\)",
        transform or "",
    )
    if not match:
        return 0.0, 0.0
    return float(match.group(1)), float(match.group(2))


def _extract_legend_vertical_bounds(root: ET.Element) -> tuple[float, float]:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    legend = root.find(".//svg:g[@id='legend']", ns)
    assert legend is not None
    _, legend_y = _parse_translate(legend.attrib.get("transform", ""))

    min_y = float("inf")
    max_y = float("-inf")
    for path in legend.findall(".//svg:path", ns):
        _, path_y = _parse_translate(path.attrib.get("transform", ""))
        d_attr = path.attrib.get("d", "")
        for match in re.finditer(
            r"[ML]\s*[-+0-9.eE]+\s*,\s*([-+0-9.eE]+)",
            d_attr,
        ):
            y_value = legend_y + path_y + float(match.group(1))
            min_y = min(min_y, y_value)
            max_y = max(max_y, y_value)

    assert min_y != float("inf")
    assert max_y != float("-inf")
    return min_y, max_y


def _extract_definition_top_y(root: ET.Element, group_id: str) -> float:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    group = root.find(f".//svg:g[@id='{group_id}']", ns)
    assert group is not None
    _, group_y = _parse_translate(group.attrib.get("transform", ""))

    min_top = float("inf")
    for text in group.findall("./svg:text", ns):
        y_raw = text.attrib.get("y", "0")
        font_size_raw = text.attrib.get("font-size", "0")
        y_val = float(str(y_raw).replace("px", ""))
        font_size = float(str(font_size_raw).replace("px", ""))
        top = group_y + y_val - (0.5 * font_size)
        min_top = min(min_top, top)

    assert min_top != float("inf")
    return min_top


def _extract_definition_bottom_y(root: ET.Element, group_id: str) -> float:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    group = root.find(f".//svg:g[@id='{group_id}']", ns)
    assert group is not None
    _, group_y = _parse_translate(group.attrib.get("transform", ""))

    max_bottom = float("-inf")
    for text in group.findall("./svg:text", ns):
        y_raw = text.attrib.get("y", "0")
        font_size_raw = text.attrib.get("font-size", "0")
        y_val = float(str(y_raw).replace("px", ""))
        font_size = float(str(font_size_raw).replace("px", ""))
        bottom = group_y + y_val + (0.5 * font_size)
        max_bottom = max(max_bottom, bottom)

    assert max_bottom != float("-inf")
    return max_bottom


@pytest.mark.circular
def test_assemble_circular_diagram_from_records_shared_legend_and_unique_ids() -> None:
    records = [
        _build_record("rec_a", 30),
        _build_record("rec_b", 260),
    ]

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="right",
        config_overrides={"show_labels": True},
    )

    root = ET.fromstring(canvas.tostring())
    ns = {"svg": "http://www.w3.org/2000/svg"}

    legends = root.findall(".//svg:g[@id='legend']", ns)
    assert len(legends) == 1

    expected_feature_ids = {
        compute_feature_hash(record.features[0], record_id=record.id) for record in records
    }
    rendered_feature_ids = {
        path.attrib["id"]
        for path in root.findall(".//svg:path[@id]", ns)
        if path.attrib.get("id", "").startswith("f")
    }
    assert expected_feature_ids.issubset(rendered_feature_ids)

    all_group_ids = [group.attrib.get("id", "") for group in root.findall(".//svg:g[@id]", ns)]

    for base_id in ("Axis", "tick", "labels", "gc_content"):
        ids = [gid for gid in all_group_ids if gid == base_id or gid.startswith(f"{base_id}_")]
        assert len(ids) == len(records)
        assert len(ids) == len(set(ids))
        assert all(gid != base_id for gid in ids)

    skew_ids = [
        gid
        for gid in all_group_ids
        if gid in {"skew", "gc_skew"}
        or gid.startswith("skew_")
        or gid.startswith("gc_skew_")
    ]
    assert len(skew_ids) == len(records)
    assert len(skew_ids) == len(set(skew_ids))
    assert all(gid not in {"skew", "gc_skew"} for gid in skew_ids)


@pytest.mark.circular
def test_assemble_circular_diagram_from_records_default_sqrt_scaling(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("max_len", 30, length=1000),
        _build_record("mid_len", 60, length=640),
    ]
    captured_radii: list[float] = []

    def fake_single(gb_record: SeqRecord, **kwargs: Any) -> Drawing:
        cfg = kwargs["cfg"]
        radius = float(cfg.canvas.circular.radius)
        width = float(cfg.canvas.circular.width.without_labels)
        height = float(cfg.canvas.circular.height)
        captured_radii.append(radius)
        return Drawing(
            filename=f"{gb_record.id}.svg",
            size=(f"{width}px", f"{height}px"),
            viewBox=f"0 0 {width} {height}",
            debug=False,
        )

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
    )

    assert len(captured_radii) == 2
    assert captured_radii[0] > captured_radii[1]
    assert captured_radii[1] / captured_radii[0] == pytest.approx(0.8, rel=1e-6)


@pytest.mark.circular
@pytest.mark.parametrize(
    ("size_mode", "min_ratio", "expected_ratio"),
    [
        ("linear", 0.55, 0.64),
        ("equal", 0.55, 1.0),
    ],
)
def test_assemble_circular_diagram_from_records_scaling_mode_selection(
    monkeypatch: pytest.MonkeyPatch,
    size_mode: str,
    min_ratio: float,
    expected_ratio: float,
) -> None:
    records = [
        _build_record("max_len", 30, length=1000),
        _build_record("mid_len", 60, length=640),
    ]
    captured_radii: list[float] = []

    def fake_single(gb_record: SeqRecord, **kwargs: Any) -> Drawing:
        cfg = kwargs["cfg"]
        radius = float(cfg.canvas.circular.radius)
        width = float(cfg.canvas.circular.width.without_labels)
        height = float(cfg.canvas.circular.height)
        captured_radii.append(radius)
        return Drawing(
            filename=f"{gb_record.id}.svg",
            size=(f"{width}px", f"{height}px"),
            viewBox=f"0 0 {width} {height}",
            debug=False,
        )

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        multi_record_size_mode=size_mode,
        multi_record_min_radius_ratio=min_ratio,
    )

    assert len(captured_radii) == 2
    assert captured_radii[1] / captured_radii[0] == pytest.approx(expected_ratio, rel=1e-6)


@pytest.mark.circular
def test_assemble_circular_diagram_from_records_min_ratio_clamp(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("max_len", 30, length=1000),
        _build_record("tiny_len", 60, length=10),
    ]
    captured_radii: list[float] = []

    def fake_single(gb_record: SeqRecord, **kwargs: Any) -> Drawing:
        cfg = kwargs["cfg"]
        radius = float(cfg.canvas.circular.radius)
        width = float(cfg.canvas.circular.width.without_labels)
        height = float(cfg.canvas.circular.height)
        captured_radii.append(radius)
        return Drawing(
            filename=f"{gb_record.id}.svg",
            size=(f"{width}px", f"{height}px"),
            viewBox=f"0 0 {width} {height}",
            debug=False,
        )

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        multi_record_size_mode="linear",
        multi_record_min_radius_ratio=0.55,
    )

    assert len(captured_radii) == 2
    assert captured_radii[1] / captured_radii[0] == pytest.approx(0.55, rel=1e-6)


@pytest.mark.circular
def test_multi_record_variable_grid_width_is_tighter_than_fixed_cell_layout(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("grid_a", 20, length=1000),
        _build_record("grid_b", 220, length=900),
        _build_record("grid_c", 420, length=800),
    ]
    planned_sizes = [
        (1000.0, 900.0),
        (700.0, 900.0),
        (520.0, 900.0),
    ]
    captured_radii: list[float] = []

    def fake_single(gb_record: SeqRecord, **kwargs: Any) -> Drawing:
        index_map = {record.id: idx for idx, record in enumerate(records)}
        width, height = planned_sizes[index_map[gb_record.id]]
        captured_radii.append(float(kwargs["cfg"].canvas.circular.radius))
        return Drawing(
            filename=f"{gb_record.id}.svg",
            size=(f"{width}px", f"{height}px"),
            viewBox=f"0 0 {width} {height}",
            debug=False,
        )

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
    )

    root = ET.fromstring(canvas.tostring())
    actual_width = _extract_viewbox_width(root)
    expected_gap = max(captured_radii) * 0.1
    fixed_cell_width = 2.0 * 1000.0
    expected_variable_width = 1000.0 + 700.0 + expected_gap

    assert actual_width == pytest.approx(expected_variable_width, rel=1e-6)
    assert actual_width < fixed_cell_width


@pytest.mark.circular
def test_multi_record_row_gap_uses_ten_percent_of_max_radius_and_keeps_row_centered(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("gap_a", 20, length=1000),
        _build_record("gap_b", 220, length=900),
        _build_record("gap_c", 420, length=800),
        _build_record("gap_d", 620, length=700),
        _build_record("gap_e", 820, length=600),
        _build_record("gap_f", 1020, length=500),
    ]
    planned_sizes = [
        (1000.0, 900.0),
        (700.0, 900.0),
        (520.0, 900.0),
        (550.0, 900.0),
        (420.0, 900.0),
        (390.0, 900.0),
    ]
    captured_radii: list[float] = []

    def fake_single(gb_record: SeqRecord, **kwargs: Any) -> Drawing:
        index_map = {record.id: idx for idx, record in enumerate(records)}
        width, height = planned_sizes[index_map[gb_record.id]]
        captured_radii.append(float(kwargs["cfg"].canvas.circular.radius))
        return Drawing(
            filename=f"{gb_record.id}.svg",
            size=(f"{width}px", f"{height}px"),
            viewBox=f"0 0 {width} {height}",
            debug=False,
        )

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
    )
    root = ET.fromstring(canvas.tostring())

    expected_gap = max(captured_radii) * 0.1
    record_x = [
        _extract_group_translate_xy(root, f"record_{index}")[0]
        for index in range(len(records))
    ]

    assert record_x[1] - (record_x[0] + planned_sizes[0][0]) == pytest.approx(expected_gap, abs=1e-6)
    assert record_x[2] - (record_x[1] + planned_sizes[1][0]) == pytest.approx(expected_gap, abs=1e-6)
    assert record_x[4] - (record_x[3] + planned_sizes[3][0]) == pytest.approx(expected_gap, abs=1e-6)
    assert record_x[5] - (record_x[4] + planned_sizes[4][0]) == pytest.approx(expected_gap, abs=1e-6)

    row0_width = planned_sizes[0][0] + planned_sizes[1][0] + planned_sizes[2][0] + 2.0 * expected_gap
    row1_width = planned_sizes[3][0] + planned_sizes[4][0] + planned_sizes[5][0] + 2.0 * expected_gap
    expected_row1_start = (row0_width - row1_width) * 0.5
    assert record_x[3] == pytest.approx(expected_row1_start, abs=1e-6)


@pytest.mark.circular
@pytest.mark.parametrize("legend_position", ["none", "top", "bottom"])
def test_multi_record_row_outer_margins_equalize_to_larger_side_for_none_top_bottom(
    monkeypatch: pytest.MonkeyPatch,
    legend_position: str,
) -> None:
    records = [
        _build_record("sym_a", 20, length=1000),
        _build_record("sym_b", 220, length=900),
        _build_record("sym_c", 420, length=800),
        _build_record("sym_d", 620, length=700),
        _build_record("sym_e", 820, length=600),
        _build_record("sym_f", 1020, length=500),
    ]
    planned_sizes = [
        (1000.0, 900.0),
        (700.0, 900.0),
        (520.0, 900.0),
        (550.0, 900.0),
        (420.0, 900.0),
        (390.0, 900.0),
    ]
    planned_insets = [120.0, 85.0, 60.0, 70.0, 58.0, 45.0]
    captured_radii: list[float] = []

    def fake_single(gb_record: SeqRecord, **kwargs: Any) -> Drawing:
        index_map = {record.id: idx for idx, record in enumerate(records)}
        index = index_map[gb_record.id]
        width, height = planned_sizes[index]
        captured_radii.append(float(kwargs["cfg"].canvas.circular.radius))
        return _build_mock_circular_subcanvas(width, height, planned_insets[index])

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend=legend_position,
    )
    root = ET.fromstring(canvas.tostring())

    expected_gap = max(captured_radii) * 0.1
    record_x = [
        _extract_group_translate_xy(root, f"record_{index}")[0]
        for index in range(len(records))
    ]

    assert record_x[1] - (record_x[0] + planned_sizes[0][0]) == pytest.approx(expected_gap, abs=1e-6)
    assert record_x[2] - (record_x[1] + planned_sizes[1][0]) == pytest.approx(expected_gap, abs=1e-6)
    assert record_x[4] - (record_x[3] + planned_sizes[3][0]) == pytest.approx(expected_gap, abs=1e-6)
    assert record_x[5] - (record_x[4] + planned_sizes[4][0]) == pytest.approx(expected_gap, abs=1e-6)

    row0_inner_width = planned_sizes[0][0] + planned_sizes[1][0] + planned_sizes[2][0] + 2.0 * expected_gap
    row1_inner_width = planned_sizes[3][0] + planned_sizes[4][0] + planned_sizes[5][0] + 2.0 * expected_gap
    row0_target = max(planned_insets[0], planned_insets[2])
    row1_target = max(planned_insets[3], planned_insets[5])
    row0_extra_left = row0_target - planned_insets[0]
    row0_extra_right = row0_target - planned_insets[2]
    row1_extra_left = row1_target - planned_insets[3]
    row1_extra_right = row1_target - planned_insets[5]
    row0_total_width = row0_inner_width + row0_extra_left + row0_extra_right
    row1_total_width = row1_inner_width + row1_extra_left + row1_extra_right
    grid_width = max(row0_total_width, row1_total_width)
    expected_row1_start = ((grid_width - row1_total_width) * 0.5) + row1_extra_left
    assert record_x[3] == pytest.approx(expected_row1_start, abs=1e-6)

    canvas_width = _extract_viewbox_width(root)
    row0_left_edge = _extract_record_axis_outer_edges(root, 0)[0]
    row0_right_edge = _extract_record_axis_outer_edges(root, 2)[1]
    row0_left_margin = row0_left_edge
    row0_right_margin = canvas_width - row0_right_edge
    assert row0_left_margin == pytest.approx(row0_right_margin, abs=1e-6)
    assert row0_left_margin == pytest.approx(row0_target, abs=1e-6)

    row1_left_edge = _extract_record_axis_outer_edges(root, 3)[0]
    row1_right_edge = _extract_record_axis_outer_edges(root, 5)[1]
    row1_left_margin = row1_left_edge
    row1_right_margin = canvas_width - row1_right_edge
    assert row1_left_margin == pytest.approx(row1_right_margin, abs=1e-6)


@pytest.mark.circular
def test_multi_record_row_margin_symmetry_not_applied_for_right_legend(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("scope_a", 20, length=1000),
        _build_record("scope_b", 220, length=900),
        _build_record("scope_c", 420, length=800),
        _build_record("scope_d", 620, length=700),
        _build_record("scope_e", 820, length=600),
        _build_record("scope_f", 1020, length=500),
    ]
    planned_sizes = [
        (1000.0, 900.0),
        (700.0, 900.0),
        (520.0, 900.0),
        (550.0, 900.0),
        (420.0, 900.0),
        (390.0, 900.0),
    ]
    planned_insets = [120.0, 85.0, 60.0, 70.0, 58.0, 45.0]
    captured_radii: list[float] = []

    def fake_single(gb_record: SeqRecord, **kwargs: Any) -> Drawing:
        index_map = {record.id: idx for idx, record in enumerate(records)}
        index = index_map[gb_record.id]
        width, height = planned_sizes[index]
        captured_radii.append(float(kwargs["cfg"].canvas.circular.radius))
        return _build_mock_circular_subcanvas(width, height, planned_insets[index])

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="right",
    )
    root = ET.fromstring(canvas.tostring())

    expected_gap = max(captured_radii) * 0.1
    record_x = [
        _extract_group_translate_xy(root, f"record_{index}")[0]
        for index in range(len(records))
    ]

    row0_inner_width = planned_sizes[0][0] + planned_sizes[1][0] + planned_sizes[2][0] + 2.0 * expected_gap
    row1_inner_width = planned_sizes[3][0] + planned_sizes[4][0] + planned_sizes[5][0] + 2.0 * expected_gap
    expected_row1_start_without_symmetry = (row0_inner_width - row1_inner_width) * 0.5
    assert record_x[3] == pytest.approx(expected_row1_start_without_symmetry, abs=1e-6)

@pytest.mark.circular
def test_circular_cli_multi_record_canvas_opt_in_saves_once(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    records = [_build_record("cli_a", 20), _build_record("cli_b", 220)]
    calls: dict[str, int] = {"single": 0, "multi": 0, "save": 0}
    captured_kwargs: dict[str, Any] = {}

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda *_args, **_kwargs: records)
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "read_feature_visibility_file", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)

    def fake_single(*_args: Any, **_kwargs: Any) -> Drawing:
        calls["single"] += 1
        return Drawing(filename=str(tmp_path / "single.svg"))

    def fake_multi(*_args: Any, **_kwargs: Any) -> Drawing:
        calls["multi"] += 1
        captured_kwargs.update(_kwargs)
        return Drawing(filename=str(tmp_path / "multi.svg"))

    def fake_save(*_args: Any, **_kwargs: Any) -> None:
        calls["save"] += 1

    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_record", fake_single)
    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_records", fake_multi)
    monkeypatch.setattr(circular_cli_module, "save_figure", fake_save)

    circular_cli_module.circular_main(
        [
            "--gbk",
            "dummy.gb",
            "--format",
            "svg",
            "--multi_record_canvas",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert calls["single"] == 0
    assert calls["multi"] == 1
    assert calls["save"] == 1
    assert captured_kwargs["multi_record_size_mode"] == "sqrt"
    assert captured_kwargs["multi_record_min_radius_ratio"] == pytest.approx(0.55)
    assert captured_kwargs["definition_position"] == "center"
    assert captured_kwargs["multi_record_definition_mode"] == "shared"
    assert captured_kwargs["shared_definition_position"] == "bottom"


@pytest.mark.circular
def test_circular_cli_multi_record_canvas_passes_size_scaling_options(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    records = [_build_record("cli_a", 20), _build_record("cli_b", 220)]
    calls: dict[str, int] = {"single": 0, "multi": 0, "save": 0}
    captured_kwargs: dict[str, Any] = {}

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda *_args, **_kwargs: records)
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "read_feature_visibility_file", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)

    def fake_single(*_args: Any, **_kwargs: Any) -> Drawing:
        calls["single"] += 1
        return Drawing(filename=str(tmp_path / "single.svg"))

    def fake_multi(*_args: Any, **_kwargs: Any) -> Drawing:
        calls["multi"] += 1
        captured_kwargs.update(_kwargs)
        return Drawing(filename=str(tmp_path / "multi.svg"))

    def fake_save(*_args: Any, **_kwargs: Any) -> None:
        calls["save"] += 1

    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_record", fake_single)
    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_records", fake_multi)
    monkeypatch.setattr(circular_cli_module, "save_figure", fake_save)

    circular_cli_module.circular_main(
        [
            "--gbk",
            "dummy.gb",
            "--format",
            "svg",
            "--multi_record_canvas",
            "--multi_record_size_mode",
            "linear",
            "--multi_record_min_radius_ratio",
            "0.4",
            "--definition_position",
            "top",
            "--multi_record_definition_mode",
            "legacy",
            "--shared_definition_position",
            "top",
            "--shared_definition_font_size",
            "30",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert calls["single"] == 0
    assert calls["multi"] == 1
    assert calls["save"] == 1
    assert captured_kwargs["multi_record_size_mode"] == "linear"
    assert captured_kwargs["multi_record_min_radius_ratio"] == pytest.approx(0.4)
    assert captured_kwargs["definition_position"] == "top"
    assert captured_kwargs["multi_record_definition_mode"] == "legacy"
    assert captured_kwargs["shared_definition_position"] == "top"
    assert captured_kwargs["cfg"].objects.definition.circular.shared_font_size == pytest.approx(30.0)


@pytest.mark.circular
def test_circular_cli_without_multi_record_canvas_keeps_per_record_saves(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    records = [_build_record("cli_c", 20), _build_record("cli_d", 220)]
    calls: dict[str, int] = {"single": 0, "multi": 0, "save": 0}
    single_kwargs: list[dict[str, Any]] = []

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda *_args, **_kwargs: records)
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "read_feature_visibility_file", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)

    def fake_single(*_args: Any, **_kwargs: Any) -> Drawing:
        calls["single"] += 1
        single_kwargs.append(dict(_kwargs))
        return Drawing(filename=str(tmp_path / "single.svg"))

    def fake_multi(*_args: Any, **_kwargs: Any) -> Drawing:
        calls["multi"] += 1
        return Drawing(filename=str(tmp_path / "multi.svg"))

    def fake_save(*_args: Any, **_kwargs: Any) -> None:
        calls["save"] += 1

    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_record", fake_single)
    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_records", fake_multi)
    monkeypatch.setattr(circular_cli_module, "save_figure", fake_save)

    circular_cli_module.circular_main(
        [
            "--gbk",
            "dummy.gb",
            "--format",
            "svg",
            "--definition_position",
            "top",
            "--multi_record_definition_mode",
            "legacy",
            "--shared_definition_position",
            "top",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert calls["single"] == len(records)
    assert calls["multi"] == 0
    assert calls["save"] == len(records)
    assert single_kwargs
    assert all(kwargs.get("definition_position") == "top" for kwargs in single_kwargs)
    assert all("multi_record_definition_mode" not in kwargs for kwargs in single_kwargs)
    assert all("shared_definition_position" not in kwargs for kwargs in single_kwargs)


@pytest.mark.circular
@pytest.mark.parametrize("ratio", ["0", "-0.1", "1.5"])
def test_circular_cli_rejects_invalid_multi_record_min_radius_ratio(ratio: str) -> None:
    with pytest.raises(SystemExit):
        circular_cli_module._get_args(
            [
                "--gbk",
                "dummy.gb",
                "--multi_record_min_radius_ratio",
                ratio,
            ]
        )


@pytest.mark.circular
def test_circular_cli_definition_layout_defaults() -> None:
    args = circular_cli_module._get_args(["--gbk", "dummy.gb"])
    assert args.definition_position == "center"
    assert args.multi_record_definition_mode == "shared"
    assert args.shared_definition_position == "bottom"
    assert args.shared_definition_font_size is None


@pytest.mark.circular
@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("--definition_position", "left"),
        ("--multi_record_definition_mode", "invalid"),
        ("--shared_definition_position", "left"),
    ],
)
def test_circular_cli_rejects_invalid_definition_layout_options(
    option: str,
    value: str,
) -> None:
    with pytest.raises(SystemExit):
        circular_cli_module._get_args(
            [
                "--gbk",
                "dummy.gb",
                option,
                value,
            ]
        )


@pytest.mark.circular
def test_multi_record_default_shared_definition_and_record_summary_content() -> None:
    records = [
        _build_record_with_source(
            "rec_alpha",
            organism="Organism alpha",
            strain="Strain A",
            feature_start=20,
            length=1500,
        ),
        _build_record_with_source(
            "rec_beta",
            organism="Organism beta",
            strain="Strain B",
            feature_start=240,
            length=900,
        ),
    ]

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
    )
    root = ET.fromstring(canvas.tostring())
    ns = {"svg": "http://www.w3.org/2000/svg"}

    shared_groups = root.findall(".//svg:g[@id='shared_definition']", ns)
    assert len(shared_groups) == 1

    shared_texts = _extract_group_texts(root, "shared_definition")
    assert shared_texts == ["Organism alpha Strain A"]
    shared_font_sizes = _extract_group_font_sizes(root, "shared_definition")
    assert shared_font_sizes == [32.0]
    shared_font_weights = _extract_group_font_weights(root, "shared_definition")
    assert shared_font_weights == ["normal"]

    for record in records:
        record_texts = _extract_group_texts(root, f"{record.id}_definition")
        assert record_texts
        assert record_texts[0] == record.id
        assert record.id in record_texts
        assert any(text.endswith("bp") for text in record_texts)
        assert any(text.endswith("% GC") for text in record_texts)
        assert record.features[0].qualifiers["organism"][0] not in record_texts
        assert record.features[0].qualifiers["strain"][0] not in record_texts
        record_font_sizes = _extract_group_font_sizes(root, f"{record.id}_definition")
        assert all(size == pytest.approx(18.0) for size in record_font_sizes)


@pytest.mark.circular
def test_shared_definition_font_size_override_only_changes_shared_definition() -> None:
    records = [
        _build_record_with_source(
            "rec_font_a",
            organism="Organism font A",
            strain="Strain A",
            feature_start=20,
            length=1500,
        ),
        _build_record_with_source(
            "rec_font_b",
            organism="Organism font B",
            strain="Strain B",
            feature_start=240,
            length=900,
        ),
    ]

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        config_overrides={"shared_definition_font_size": 30},
    )
    root = ET.fromstring(canvas.tostring())

    shared_font_sizes = _extract_group_font_sizes(root, "shared_definition")
    assert shared_font_sizes == [30.0]
    shared_font_weights = _extract_group_font_weights(root, "shared_definition")
    assert shared_font_weights == ["normal"]

    for record in records:
        record_font_sizes = _extract_group_font_sizes(root, f"{record.id}_definition")
        assert all(size == pytest.approx(18.0) for size in record_font_sizes)


@pytest.mark.circular
@pytest.mark.parametrize(
    ("organism", "strain", "expected_shared"),
    [
        ("Species only", "", "Species only"),
        ("", "Strain only", "Strain only"),
    ],
)
def test_shared_definition_single_line_with_missing_species_or_strain(
    organism: str,
    strain: str,
    expected_shared: str,
) -> None:
    records = [
        _build_record_with_source(
            "rec_partial_a",
            organism=organism,
            strain=strain,
            length=1000,
        ),
        _build_record_with_source(
            "rec_partial_b",
            organism="Other species",
            strain="Other strain",
            length=800,
        ),
    ]

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
    )
    root = ET.fromstring(canvas.tostring())

    shared_texts = _extract_group_texts(root, "shared_definition")
    assert shared_texts == [expected_shared]


@pytest.mark.circular
def test_shared_definition_position_moves_shared_group_vertically() -> None:
    records = [
        _build_record_with_source(
            "rec_top_bottom_a",
            organism="Organism top bottom A",
            strain="Strain A",
            length=1300,
        ),
        _build_record_with_source(
            "rec_top_bottom_b",
            organism="Organism top bottom B",
            strain="Strain B",
            length=700,
        ),
    ]

    top_canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        shared_definition_position="top",
    )
    center_canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        shared_definition_position="center",
    )
    bottom_canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        shared_definition_position="bottom",
    )

    y_top = _extract_group_translate_y(ET.fromstring(top_canvas.tostring()), "shared_definition")
    y_center = _extract_group_translate_y(
        ET.fromstring(center_canvas.tostring()),
        "shared_definition",
    )
    y_bottom = _extract_group_translate_y(
        ET.fromstring(bottom_canvas.tostring()),
        "shared_definition",
    )

    assert y_top < y_center < y_bottom


@pytest.mark.circular
def test_single_record_definition_position_moves_group_vertically() -> None:
    record = _build_record_with_source(
        "single_layout",
        organism="Single organism",
        strain="Single strain",
        length=1200,
    )

    top_canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS"],
        legend="none",
        definition_position="top",
    )
    center_canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS"],
        legend="none",
        definition_position="center",
    )
    bottom_canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS"],
        legend="none",
        definition_position="bottom",
    )

    group_id = f"{record.id}_definition"
    y_top = _extract_group_translate_y(ET.fromstring(top_canvas.tostring()), group_id)
    y_center = _extract_group_translate_y(ET.fromstring(center_canvas.tostring()), group_id)
    y_bottom = _extract_group_translate_y(ET.fromstring(bottom_canvas.tostring()), group_id)

    assert y_top < y_center < y_bottom


@pytest.mark.circular
def test_single_record_legend_top_bottom_positions_expand_and_move() -> None:
    record = _build_record_with_source(
        "legend_top_bottom_single",
        organism="Single legend organism",
        strain="Single legend strain",
        length=1400,
    )

    top_canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS"],
        legend="top",
    )
    left_canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS"],
        legend="left",
    )
    bottom_canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS"],
        legend="bottom",
    )
    none_canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS"],
        legend="none",
    )

    top_root = ET.fromstring(top_canvas.tostring())
    left_root = ET.fromstring(left_canvas.tostring())
    bottom_root = ET.fromstring(bottom_canvas.tostring())
    none_root = ET.fromstring(none_canvas.tostring())

    top_x, top_y = _extract_group_translate_xy(top_root, "legend")
    bottom_x, bottom_y = _extract_group_translate_xy(bottom_root, "legend")
    left_y = _extract_group_translate_y(left_root, "legend")

    assert top_x != 0.0 or top_y != 0.0
    assert bottom_x != 0.0 or bottom_y != 0.0
    assert top_y < left_y < bottom_y

    none_height = _extract_viewbox_height(none_root)
    top_height = _extract_viewbox_height(top_root)
    bottom_height = _extract_viewbox_height(bottom_root)
    assert top_height > none_height
    assert bottom_height > none_height


@pytest.mark.circular
@pytest.mark.parametrize(
    ("legend_position", "show_labels", "track_type"),
    [
        ("top", False, "tuckin"),
        ("bottom", False, "tuckin"),
        ("top", True, "spreadout"),
        ("bottom", True, "spreadout"),
    ],
)
def test_single_record_top_bottom_legend_centers_between_edge_and_content(
    monkeypatch: pytest.MonkeyPatch,
    legend_position: str,
    show_labels: bool,
    track_type: str,
) -> None:
    record = _build_dense_record_with_source(
        f"legend_midpoint_{legend_position}_{track_type}_{'labels' if show_labels else 'nolabels'}",
        organism="Legend midpoint organism",
        strain="Legend midpoint strain",
        length=2600,
    )
    captured: dict[str, float] = {}
    original = circular_assemble_module._position_single_top_bottom_legend_between_edge_and_content

    def wrapped(
        canvas: Drawing,
        canvas_config: Any,
        legend_config: Any,
        *,
        position: str,
        content_top: float,
        content_bottom: float,
        edge_min_px: float = 16.0,
        content_gap_px: float = 12.0,
    ) -> None:
        captured["position"] = 0.0 if position == "top" else 1.0
        captured["content_top"] = float(content_top)
        captured["content_bottom"] = float(content_bottom)
        captured["edge_min_px"] = float(edge_min_px)
        captured["content_gap_px"] = float(content_gap_px)
        captured["legend_height"] = float(legend_config.legend_height)
        captured["pre_total_height"] = float(canvas_config.total_height)
        original(
            canvas,
            canvas_config,
            legend_config,
            position=position,
            content_top=content_top,
            content_bottom=content_bottom,
            edge_min_px=edge_min_px,
            content_gap_px=content_gap_px,
        )
        captured["post_total_height"] = float(canvas_config.total_height)

    monkeypatch.setattr(
        circular_assemble_module,
        "_position_single_top_bottom_legend_between_edge_and_content",
        wrapped,
    )

    canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS", "rRNA", "tRNA"],
        legend=legend_position,
        config_overrides={
            "show_labels": show_labels,
            "track_type": track_type,
        },
    )
    root = ET.fromstring(canvas.tostring())
    legend_top, legend_bottom = _extract_legend_vertical_bounds(root)
    viewbox_height = _extract_viewbox_height(root)

    assert "legend_height" in captured
    edge_min_px = captured["edge_min_px"]
    content_gap_px = captured["content_gap_px"]
    legend_height = captured["legend_height"]

    assert legend_top >= edge_min_px - 1e-6
    assert legend_bottom <= viewbox_height - edge_min_px + 1e-6

    actual_center = (legend_top + legend_bottom) / 2.0
    if legend_position == "top":
        lane_top = edge_min_px
        lane_bottom_before = captured["content_top"] - content_gap_px
        free_height = lane_bottom_before - lane_top
        missing = max(0.0, legend_height - free_height)
        final_content_top = captured["content_top"] + missing
        expected_center = (lane_top + (final_content_top - content_gap_px)) / 2.0
        assert captured["post_total_height"] == pytest.approx(
            captured["pre_total_height"] + missing, rel=1e-6
        )
    else:
        lane_top = captured["content_bottom"] + content_gap_px
        lane_bottom_before = captured["pre_total_height"] - edge_min_px
        free_height = lane_bottom_before - lane_top
        missing = max(0.0, legend_height - free_height)
        final_lane_bottom = captured["pre_total_height"] + missing - edge_min_px
        expected_center = (lane_top + final_lane_bottom) / 2.0
        assert captured["post_total_height"] == pytest.approx(
            captured["pre_total_height"] + missing, rel=1e-6
        )

    assert actual_center == pytest.approx(expected_center, abs=2.0)


@pytest.mark.circular
def test_multi_record_legend_top_bottom_positions_expand_and_move() -> None:
    records = [
        _build_record_with_source(
            "legend_grid_top_bottom_a",
            organism="Legend grid A",
            strain="Strain A",
            length=1300,
        ),
        _build_record_with_source(
            "legend_grid_top_bottom_b",
            organism="Legend grid B",
            strain="Strain B",
            length=900,
        ),
    ]

    top_canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="top",
    )
    left_canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="left",
    )
    bottom_canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="bottom",
    )
    none_canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
    )

    top_root = ET.fromstring(top_canvas.tostring())
    left_root = ET.fromstring(left_canvas.tostring())
    bottom_root = ET.fromstring(bottom_canvas.tostring())
    none_root = ET.fromstring(none_canvas.tostring())

    _, top_y = _extract_group_translate_xy(top_root, "legend")
    left_y = _extract_group_translate_y(left_root, "legend")
    _, bottom_y = _extract_group_translate_xy(bottom_root, "legend")
    assert top_y < left_y < bottom_y

    none_height = _extract_viewbox_height(none_root)
    top_height = _extract_viewbox_height(top_root)
    bottom_height = _extract_viewbox_height(bottom_root)
    assert top_height > none_height
    assert bottom_height > none_height


@pytest.mark.circular
def test_multi_record_top_legend_keeps_minimum_top_padding() -> None:
    records = [
        _build_record_with_source(
            "legend_top_padding_a",
            organism="Legend top padding A",
            strain="Strain A",
            length=1300,
        ),
        _build_record_with_source(
            "legend_top_padding_b",
            organism="Legend top padding B",
            strain="Strain B",
            length=900,
        ),
    ]

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="top",
    )
    root = ET.fromstring(canvas.tostring())
    legend_top, _legend_bottom = _extract_legend_vertical_bounds(root)
    assert legend_top >= 32.0 - 1e-6


@pytest.mark.circular
def test_multi_record_bottom_shared_definition_stays_below_legend() -> None:
    records = [
        _build_record_with_source(
            "legend_bottom_shared_a",
            organism="Legend bottom shared A",
            strain="Strain A",
            length=1300,
        ),
        _build_record_with_source(
            "legend_bottom_shared_b",
            organism="Legend bottom shared B",
            strain="Strain B",
            length=900,
        ),
    ]

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="bottom",
        multi_record_definition_mode="shared",
        shared_definition_position="bottom",
    )
    root = ET.fromstring(canvas.tostring())
    _legend_top, legend_bottom = _extract_legend_vertical_bounds(root)
    shared_top = _extract_definition_top_y(root, "shared_definition")
    assert shared_top >= legend_bottom + 20.0 - 1e-6


@pytest.mark.circular
@pytest.mark.parametrize("legend_position", ["left", "right"])
def test_multi_record_left_right_shared_definition_bottom_keeps_margins(
    legend_position: str,
) -> None:
    records = [
        _build_record_with_source(
            "legend_side_shared_a",
            organism="Legend side shared A",
            strain="Strain A",
            length=1300,
        ),
        _build_record_with_source(
            "legend_side_shared_b",
            organism="Legend side shared B",
            strain="Strain B",
            length=900,
        ),
    ]

    baseline_canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend=legend_position,
        multi_record_definition_mode="shared",
        shared_definition_position="center",
    )
    baseline_root = ET.fromstring(baseline_canvas.tostring())
    baseline_height = _extract_viewbox_height(baseline_root)

    bottom_canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend=legend_position,
        multi_record_definition_mode="shared",
        shared_definition_position="bottom",
    )
    bottom_root = ET.fromstring(bottom_canvas.tostring())
    bottom_height = _extract_viewbox_height(bottom_root)
    shared_top = _extract_definition_top_y(bottom_root, "shared_definition")
    shared_bottom = _extract_definition_bottom_y(bottom_root, "shared_definition")

    assert shared_top >= baseline_height + 20.0 - 1e-6
    assert bottom_height - shared_bottom >= 24.0 - 1e-6


@pytest.mark.circular
def test_circular_top_legend_uses_horizontal_entry_layout() -> None:
    record = SeqRecord(Seq("A" * 1600), id="legend_layout_top")
    record.features = [
        SeqFeature(
            FeatureLocation(0, 1600, strand=1),
            type="source",
            qualifiers={"organism": ["Legend layout"], "strain": ["Layout"]},
        ),
        SeqFeature(
            FeatureLocation(50, 200, strand=1),
            type="CDS",
            qualifiers={"product": ["protein_cds"]},
        ),
        SeqFeature(
            FeatureLocation(350, 460, strand=1),
            type="tRNA",
            qualifiers={"product": ["protein_trna"]},
        ),
        SeqFeature(
            FeatureLocation(700, 920, strand=1),
            type="rRNA",
            qualifiers={"product": ["protein_rrna"]},
        ),
    ]

    top_canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS", "tRNA", "rRNA"],
        legend="top",
    )
    right_canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS", "tRNA", "rRNA"],
        legend="right",
    )

    top_texts = _extract_legend_text_transforms(ET.fromstring(top_canvas.tostring()))
    right_texts = _extract_legend_text_transforms(ET.fromstring(right_canvas.tostring()))
    assert len(top_texts) >= 2
    assert len(right_texts) >= 2

    top_unique_y = {round(y, 3) for _, y, _ in top_texts}
    right_unique_y = {round(y, 3) for _, y, _ in right_texts}
    assert len(top_unique_y) < len(top_texts)
    assert len(right_unique_y) == len(right_texts)


@pytest.mark.circular
@pytest.mark.parametrize("legend_position", ["top", "bottom"])
def test_single_record_top_bottom_legend_centers_each_wrapped_row(
    legend_position: str,
) -> None:
    record = _build_multi_feature_record_with_source(
        f"single_legend_center_{legend_position}",
        organism="Legend row center organism",
        strain="Legend row center strain",
        length=2600,
    )
    selected_features = ["CDS", "tRNA", "rRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]

    canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=selected_features,
        legend=legend_position,
    )
    root = ET.fromstring(canvas.tostring())
    row_centers = _extract_horizontal_legend_row_centers(root)

    assert len(row_centers) >= 2
    expected_center = row_centers[0][1]
    for _row_y, center_x in row_centers[1:]:
        assert center_x == pytest.approx(expected_center, abs=2.0)


@pytest.mark.circular
@pytest.mark.parametrize("legend_position", ["top", "bottom"])
def test_multi_record_top_bottom_legend_centers_each_wrapped_row(
    legend_position: str,
) -> None:
    records = [
        _build_multi_feature_record_with_source(
            "multi_legend_center_a",
            organism="Multi legend center A",
            strain="Strain A",
            length=2600,
        ),
        _build_multi_feature_record_with_source(
            "multi_legend_center_b",
            organism="Multi legend center B",
            strain="Strain B",
            length=2100,
        ),
    ]
    selected_features = ["CDS", "tRNA", "rRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=selected_features,
        legend=legend_position,
    )
    root = ET.fromstring(canvas.tostring())
    row_centers = _extract_horizontal_legend_row_centers(root)

    assert len(row_centers) >= 2
    expected_center = row_centers[0][1]
    for _row_y, center_x in row_centers[1:]:
        assert center_x == pytest.approx(expected_center, abs=2.0)


@pytest.mark.circular
def test_multi_record_legacy_keeps_full_per_record_definition() -> None:
    records = [
        _build_record_with_source(
            "legacy_a",
            organism="Legacy organism A",
            strain="Legacy strain A",
            length=1400,
        ),
        _build_record_with_source(
            "legacy_b",
            organism="Legacy organism B",
            strain="Legacy strain B",
            length=800,
        ),
    ]

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        multi_record_definition_mode="legacy",
    )
    root = ET.fromstring(canvas.tostring())
    ns = {"svg": "http://www.w3.org/2000/svg"}

    shared_groups = root.findall(".//svg:g[@id='shared_definition']", ns)
    assert len(shared_groups) == 0

    for record in records:
        record_texts = _extract_group_texts(root, f"{record.id}_definition")
        assert record.features[0].qualifiers["organism"][0] in record_texts
        assert record.features[0].qualifiers["strain"][0] in record_texts


@pytest.mark.circular
def test_build_circular_diagram_passes_definition_position_option(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    record = _build_record("build_options_record", 30)
    captured_kwargs: dict[str, Any] = {}

    def fake_single(_record: SeqRecord, **kwargs: Any) -> Drawing:
        captured_kwargs.update(kwargs)
        return Drawing(filename="dummy.svg")

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    diagram_api_module.build_circular_diagram(
        record,
        options=DiagramOptions(
            output=OutputOptions(definition_position="bottom"),
        ),
    )

    assert captured_kwargs["definition_position"] == "bottom"
