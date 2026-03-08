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
import gbdraw.render.groups.circular.ticks as circular_ticks_group_module
from gbdraw.api.diagram import assemble_circular_diagram_from_records
from gbdraw.api.options import DiagramOptions, OutputOptions
from gbdraw.core.text import calculate_bbox_dimensions
from gbdraw.exceptions import ValidationError
from gbdraw.features.colors import compute_feature_hash
from gbdraw.svg.circular_ticks import get_circular_tick_path_ratio_bounds


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


def _build_distinct_circular_style_config_dict() -> dict[str, Any]:
    config_dict = diagram_api_module.load_config_toml("gbdraw.data", "config.toml")
    config_dict["canvas"]["circular"]["track_ratio_factors"]["short"][0] = 0.77
    config_dict["canvas"]["circular"]["track_ratio_factors"]["long"][0] = 0.33
    config_dict["canvas"]["circular"]["track_ratio_factors"]["short"][1] = 0.91
    config_dict["canvas"]["circular"]["track_ratio_factors"]["long"][1] = 0.41
    config_dict["canvas"]["circular"]["track_ratio_factors"]["short"][2] = 0.81
    config_dict["canvas"]["circular"]["track_ratio_factors"]["long"][2] = 0.31
    for track_type in ("tuckin", "middle", "spreadout"):
        config_dict["canvas"]["circular"]["track_dict"]["short"][track_type]["2"] = 0.24
        config_dict["canvas"]["circular"]["track_dict"]["short"][track_type]["3"] = 0.19
        config_dict["canvas"]["circular"]["track_dict"]["long"][track_type]["2"] = 0.74
        config_dict["canvas"]["circular"]["track_dict"]["long"][track_type]["3"] = 0.59
    config_dict["objects"]["features"]["block_stroke_width"]["short"] = 9.0
    config_dict["objects"]["features"]["block_stroke_width"]["long"] = 1.0
    config_dict["objects"]["features"]["line_stroke_width"]["short"] = 8.0
    config_dict["objects"]["features"]["line_stroke_width"]["long"] = 2.0
    config_dict["objects"]["axis"]["circular"]["stroke_width"]["short"] = 7.0
    config_dict["objects"]["axis"]["circular"]["stroke_width"]["long"] = 3.0
    return config_dict


def _build_record_with_source(
    record_id: str,
    *,
    organism: str,
    strain: str,
    feature_start: int = 30,
    length: int = 1200,
    chromosome: str | None = None,
    plasmid: str | None = None,
) -> SeqRecord:
    source_qualifiers: dict[str, list[str]] = {
        "organism": [organism],
        "strain": [strain],
    }
    if chromosome:
        source_qualifiers["chromosome"] = [chromosome]
    if plasmid:
        source_qualifiers["plasmid"] = [plasmid]

    record = SeqRecord(Seq("A" * length), id=record_id)
    record.features = [
        SeqFeature(
            FeatureLocation(0, length, strand=1),
            type="source",
            qualifiers=source_qualifiers,
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


def _extract_group_italic_tspan_texts(root: ET.Element, group_id: str) -> list[str]:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    group = root.find(f".//svg:g[@id='{group_id}']", ns)
    assert group is not None
    return [
        "".join(tspan.itertext()).strip()
        for tspan in group.findall(".//svg:tspan[@font-style='italic']", ns)
        if "".join(tspan.itertext()).strip()
    ]


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


def _extract_record_axis_vertical_edges(root: ET.Element, record_index: int) -> tuple[float, float]:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    record_group = root.find(f".//svg:g[@id='record_{record_index}']", ns)
    assert record_group is not None
    _record_x, record_y = _parse_translate(record_group.attrib.get("transform", ""))

    axis_group = record_group.find(f"./svg:g[@id='Axis_{record_index}']", ns)
    assert axis_group is not None
    _axis_x, axis_y = _parse_translate(axis_group.attrib.get("transform", ""))

    circle = axis_group.find("./svg:circle", ns)
    assert circle is not None
    radius = float(circle.attrib.get("r", "0"))

    center_y = record_y + axis_y
    return center_y - radius, center_y + radius


def _extract_record_axis_and_definition_center_y(
    root: ET.Element,
    *,
    record_index: int,
    record_id: str,
) -> tuple[float, float]:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    record_group = root.find(f".//svg:g[@id='record_{record_index}']", ns)
    assert record_group is not None
    _record_x, record_y = _parse_translate(record_group.attrib.get("transform", ""))

    axis_group = record_group.find(f"./svg:g[@id='Axis_{record_index}']", ns)
    assert axis_group is not None
    _axis_x, axis_y = _parse_translate(axis_group.attrib.get("transform", ""))
    axis_center_y = record_y + axis_y

    track_id = str(record_id).replace(" ", "_")
    definition_group = record_group.find(f"./svg:g[@id='{track_id}_definition']", ns)
    assert definition_group is not None
    _definition_x, definition_y = _parse_translate(definition_group.attrib.get("transform", ""))

    min_top = float("inf")
    max_bottom = float("-inf")
    for text in definition_group.findall("./svg:text", ns):
        y_raw = text.attrib.get("y", "0")
        font_size_raw = text.attrib.get("font-size", "0")
        y_val = float(str(y_raw).replace("px", ""))
        font_size = float(str(font_size_raw).replace("px", ""))
        min_top = min(min_top, y_val - (0.5 * font_size))
        max_bottom = max(max_bottom, y_val + (0.5 * font_size))

    assert min_top != float("inf")
    assert max_bottom != float("-inf")
    definition_center_y = record_y + definition_y + ((min_top + max_bottom) * 0.5)
    return axis_center_y, definition_center_y


def _extract_single_axis_and_definition_center_y(
    root: ET.Element,
    *,
    record_id: str,
) -> tuple[float, float]:
    ns = {"svg": "http://www.w3.org/2000/svg"}
    axis_group = root.find(".//svg:g[@id='Axis']", ns)
    assert axis_group is not None
    _axis_x, axis_y = _parse_translate(axis_group.attrib.get("transform", ""))

    circle = axis_group.find("./svg:circle", ns)
    assert circle is not None
    _radius = float(circle.attrib.get("r", "0"))
    axis_center_y = axis_y

    track_id = str(record_id).replace(" ", "_")
    definition_group = root.find(f".//svg:g[@id='{track_id}_definition']", ns)
    assert definition_group is not None
    _definition_x, definition_y = _parse_translate(definition_group.attrib.get("transform", ""))

    min_top = float("inf")
    max_bottom = float("-inf")
    for text in definition_group.findall("./svg:text", ns):
        y_raw = text.attrib.get("y", "0")
        font_size_raw = text.attrib.get("font-size", "0")
        y_val = float(str(y_raw).replace("px", ""))
        font_size = float(str(font_size_raw).replace("px", ""))
        min_top = min(min_top, y_val - (0.5 * font_size))
        max_bottom = max(max_bottom, y_val + (0.5 * font_size))

    assert min_top != float("inf")
    assert max_bottom != float("-inf")
    definition_center_y = definition_y + ((min_top + max_bottom) * 0.5)
    return axis_center_y, definition_center_y


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

    skew_clip_ids: list[str] = []
    for skew_id in skew_ids:
        skew_group = root.find(f".//svg:g[@id='{skew_id}']", ns)
        assert skew_group is not None

        local_clip_paths = skew_group.findall("./svg:clipPath[@id]", ns)
        assert len(local_clip_paths) == 1
        clip_id = local_clip_paths[0].attrib["id"]
        skew_clip_ids.append(clip_id)

        low_paths = [
            path
            for path in skew_group.findall("./svg:path", ns)
            if path.attrib.get("clip-path")
        ]
        assert len(low_paths) == 1
        assert low_paths[0].attrib.get("clip-path") == f"url(#{clip_id})"

    assert len(skew_clip_ids) == len(set(skew_clip_ids))


@pytest.mark.circular
def test_assemble_circular_diagram_from_records_default_auto_scaling(
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
def test_multi_record_mixed_lengths_harmonize_short_feature_axis_style_to_long(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("long_len", 30, length=120_000),
        _build_record("short_len", 60, length=20_000),
    ]
    config_dict = _build_distinct_circular_style_config_dict()
    captured_styles: dict[str, dict[str, float]] = {}
    captured_tick_channels: dict[str, str | None] = {}

    def fake_single(gb_record: SeqRecord, **kwargs: Any) -> Drawing:
        cfg = kwargs["cfg"]
        threshold = int(cfg.labels.length_threshold.circular)
        length_param = "short" if len(gb_record.seq) < threshold else "long"
        track_type = str(cfg.canvas.circular.track_type)
        captured_tick_channels[gb_record.id] = kwargs.get("_tick_track_channel_override")
        captured_styles[gb_record.id] = {
            "feature_ratio_factor": float(cfg.canvas.circular.track_ratio_factors[length_param][0]),
            "gc_width_ratio_factor": float(cfg.canvas.circular.track_ratio_factors[length_param][1]),
            "skew_width_ratio_factor": float(cfg.canvas.circular.track_ratio_factors[length_param][2]),
            "block_stroke_width": float(cfg.objects.features.block_stroke_width.for_length_param(length_param)),
            "line_stroke_width": float(cfg.objects.features.line_stroke_width.for_length_param(length_param)),
            "axis_stroke_width": float(cfg.objects.axis.circular.stroke_width.for_length_param(length_param)),
            "gc_center_ratio": float(cfg.canvas.circular.track_dict[length_param][track_type]["2"]),
            "skew_center_ratio": float(cfg.canvas.circular.track_dict[length_param][track_type]["3"]),
        }
        width = float(cfg.canvas.circular.width.without_labels)
        height = float(cfg.canvas.circular.height)
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
        config_dict=config_dict,
        selected_features_set=["CDS"],
        legend="none",
    )

    assert captured_styles["short_len"] == pytest.approx(captured_styles["long_len"], rel=1e-9)
    assert captured_tick_channels == {"long_len": "long", "short_len": "long"}
    assert captured_styles["short_len"]["feature_ratio_factor"] == pytest.approx(0.33, rel=1e-9)
    assert captured_styles["short_len"]["gc_width_ratio_factor"] == pytest.approx(0.41, rel=1e-9)
    assert captured_styles["short_len"]["skew_width_ratio_factor"] == pytest.approx(0.31, rel=1e-9)
    assert captured_styles["short_len"]["block_stroke_width"] == pytest.approx(1.0, rel=1e-9)
    assert captured_styles["short_len"]["line_stroke_width"] == pytest.approx(2.0, rel=1e-9)
    assert captured_styles["short_len"]["axis_stroke_width"] == pytest.approx(3.0, rel=1e-9)
    assert captured_styles["short_len"]["gc_center_ratio"] == pytest.approx(0.74, rel=1e-9)
    assert captured_styles["short_len"]["skew_center_ratio"] == pytest.approx(0.59, rel=1e-9)


@pytest.mark.circular
def test_multi_record_all_short_keeps_short_feature_axis_style(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("short_a", 30, length=20_000),
        _build_record("short_b", 60, length=30_000),
    ]
    config_dict = _build_distinct_circular_style_config_dict()
    captured_styles: dict[str, dict[str, float]] = {}
    captured_tick_channels: dict[str, str | None] = {}

    def fake_single(gb_record: SeqRecord, **kwargs: Any) -> Drawing:
        cfg = kwargs["cfg"]
        threshold = int(cfg.labels.length_threshold.circular)
        length_param = "short" if len(gb_record.seq) < threshold else "long"
        track_type = str(cfg.canvas.circular.track_type)
        captured_tick_channels[gb_record.id] = kwargs.get("_tick_track_channel_override")
        captured_styles[gb_record.id] = {
            "feature_ratio_factor": float(cfg.canvas.circular.track_ratio_factors[length_param][0]),
            "gc_width_ratio_factor": float(cfg.canvas.circular.track_ratio_factors[length_param][1]),
            "skew_width_ratio_factor": float(cfg.canvas.circular.track_ratio_factors[length_param][2]),
            "block_stroke_width": float(cfg.objects.features.block_stroke_width.for_length_param(length_param)),
            "line_stroke_width": float(cfg.objects.features.line_stroke_width.for_length_param(length_param)),
            "axis_stroke_width": float(cfg.objects.axis.circular.stroke_width.for_length_param(length_param)),
            "gc_center_ratio": float(cfg.canvas.circular.track_dict[length_param][track_type]["2"]),
            "skew_center_ratio": float(cfg.canvas.circular.track_dict[length_param][track_type]["3"]),
        }
        width = float(cfg.canvas.circular.width.without_labels)
        height = float(cfg.canvas.circular.height)
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
        config_dict=config_dict,
        selected_features_set=["CDS"],
        legend="none",
    )

    for style in captured_styles.values():
        assert style["gc_width_ratio_factor"] == pytest.approx(0.91, rel=1e-9)
        assert style["skew_width_ratio_factor"] == pytest.approx(0.81, rel=1e-9)
        assert style["feature_ratio_factor"] == pytest.approx(0.77, rel=1e-9)
        assert style["block_stroke_width"] == pytest.approx(9.0, rel=1e-9)
        assert style["line_stroke_width"] == pytest.approx(8.0, rel=1e-9)
        assert style["axis_stroke_width"] == pytest.approx(7.0, rel=1e-9)
        assert style["gc_center_ratio"] == pytest.approx(0.24, rel=1e-9)
        assert style["skew_center_ratio"] == pytest.approx(0.19, rel=1e-9)
    assert captured_tick_channels == {"short_a": None, "short_b": None}


@pytest.mark.circular
def test_multi_record_mixed_lengths_force_long_tick_channel_for_short_record(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("long_len", 30, length=120_000),
        _build_record("short_len", 60, length=20_000),
    ]
    captured_tick_path_channels: dict[int, str | None] = {}
    captured_tick_label_channels: dict[int, str | None] = {}

    def fake_generate_tick_paths(
        radius: float,
        total_len: int,
        size: str,
        ticks: list[int],
        tick_width: float,
        track_type: str,
        strandedness: bool,
        tick_track_channel_override: str | None = None,
    ) -> list[Any]:
        captured_tick_path_channels[int(total_len)] = tick_track_channel_override
        return []

    def fake_generate_tick_labels(
        radius: float,
        total_len: int,
        size: str,
        ticks: list[int],
        stroke: str,
        fill: str,
        font_size: float,
        font_weight: str,
        font_family: str,
        track_type: str,
        strandedness: bool,
        dpi: int,
        tick_track_channel_override: str | None = None,
    ) -> list[Any]:
        captured_tick_label_channels[int(total_len)] = tick_track_channel_override
        return []

    monkeypatch.setattr(
        circular_ticks_group_module,
        "generate_circular_tick_paths",
        fake_generate_tick_paths,
    )
    monkeypatch.setattr(
        circular_ticks_group_module,
        "generate_circular_tick_labels",
        fake_generate_tick_labels,
    )

    assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        config_overrides={"show_gc": False, "show_skew": False},
    )

    assert captured_tick_path_channels[20_000] == "long"
    assert captured_tick_label_channels[20_000] == "long"
    assert captured_tick_path_channels[120_000] == "long"
    assert captured_tick_label_channels[120_000] == "long"

    default_short_tick_bounds = get_circular_tick_path_ratio_bounds(
        20_000,
        "tuckin",
        False,
    )
    forced_long_tick_bounds = get_circular_tick_path_ratio_bounds(
        20_000,
        "tuckin",
        False,
        tick_track_channel_override="long",
    )
    assert default_short_tick_bounds[1] >= 0.98
    assert forced_long_tick_bounds[1] < 1.0


@pytest.mark.circular
def test_multi_record_mixed_lengths_keep_gc_window_step_per_record_defaults(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("long_len", 30, length=1_200_000),
        _build_record("short_len", 60, length=20_000),
    ]
    captured_gc_window_steps: dict[str, tuple[int, int]] = {}

    def fake_assemble(**kwargs: Any) -> Drawing:
        gb_record = kwargs["gb_record"]
        canvas_config = kwargs["canvas_config"]
        gc_config = kwargs["gc_config"]
        captured_gc_window_steps[gb_record.id] = (
            int(gc_config.window),
            int(gc_config.step),
        )
        width = float(canvas_config.total_width)
        height = float(canvas_config.total_height)
        return Drawing(
            filename=f"{gb_record.id}.svg",
            size=(f"{width}px", f"{height}px"),
            viewBox=f"0 0 {width} {height}",
            debug=False,
        )

    monkeypatch.setattr(diagram_api_module, "assemble_circular_diagram", fake_assemble)

    expected_cfg = diagram_api_module.GbdrawConfig.from_dict(
        diagram_api_module.load_config_toml("gbdraw.data", "config.toml")
    )
    expected_short = tuple(expected_cfg.objects.sliding_window.default)
    expected_long = tuple(expected_cfg.objects.sliding_window.up1m)

    assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        config_overrides={"show_gc": False, "show_skew": False},
    )

    assert captured_gc_window_steps["short_len"] == expected_short
    assert captured_gc_window_steps["long_len"] == expected_long


@pytest.mark.circular
@pytest.mark.parametrize(
    ("size_mode", "min_ratio", "expected_ratio"),
    [
        ("auto", 0.55, 0.8),
        ("sqrt", 0.55, 0.8),
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
def test_assemble_circular_diagram_from_records_auto_renormalizes_when_multiple_would_clamp(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("max_len", 30, length=1000),
        _build_record("mid_len", 60, length=300),
        _build_record("small_len", 90, length=200),
        _build_record("tiny_len", 120, length=100),
    ]
    captured_radii_by_id: dict[str, float] = {}

    def fake_single(gb_record: SeqRecord, **kwargs: Any) -> Drawing:
        cfg = kwargs["cfg"]
        radius = float(cfg.canvas.circular.radius)
        width = float(cfg.canvas.circular.width.without_labels)
        height = float(cfg.canvas.circular.height)
        captured_radii_by_id[gb_record.id] = radius
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
        multi_record_size_mode="auto",
        multi_record_min_radius_ratio=0.55,
    )

    assert len(captured_radii_by_id) == 4
    max_radius = captured_radii_by_id["max_len"]
    tiny_ratio = captured_radii_by_id["tiny_len"] / max_radius
    small_ratio = captured_radii_by_id["small_len"] / max_radius
    mid_ratio = captured_radii_by_id["mid_len"] / max_radius

    assert tiny_ratio == pytest.approx(0.55, rel=1e-6)
    assert small_ratio > 0.55
    assert mid_ratio > small_ratio


@pytest.mark.circular
def test_assemble_circular_diagram_from_records_auto_keeps_single_clamp_behavior(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("max_len", 30, length=1000),
        _build_record("mid_len", 60, length=640),
        _build_record("tiny_len", 90, length=100),
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
        multi_record_size_mode="auto",
        multi_record_min_radius_ratio=0.55,
    )

    assert len(captured_radii) == 3
    assert captured_radii[1] / captured_radii[0] == pytest.approx(0.8, rel=1e-6)
    assert captured_radii[2] / captured_radii[0] == pytest.approx(0.55, rel=1e-6)


@pytest.mark.circular
def test_assemble_circular_diagram_from_records_sqrt_alias_matches_auto(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("max_len", 30, length=1000),
        _build_record("mid_len", 60, length=300),
        _build_record("small_len", 90, length=200),
        _build_record("tiny_len", 120, length=100),
    ]
    captured_radii: dict[str, list[float]] = {"auto": [], "sqrt": []}
    active_mode = {"value": "auto"}

    def fake_single(gb_record: SeqRecord, **kwargs: Any) -> Drawing:
        cfg = kwargs["cfg"]
        radius = float(cfg.canvas.circular.radius)
        width = float(cfg.canvas.circular.width.without_labels)
        height = float(cfg.canvas.circular.height)
        captured_radii[active_mode["value"]].append(radius)
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

    active_mode["value"] = "auto"
    assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        multi_record_size_mode="auto",
        multi_record_min_radius_ratio=0.55,
    )
    active_mode["value"] = "sqrt"
    assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        multi_record_size_mode="sqrt",
        multi_record_min_radius_ratio=0.55,
    )

    assert captured_radii["auto"] == pytest.approx(captured_radii["sqrt"], rel=1e-6)


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
def test_multi_record_default_column_and_row_gap_ratios_are_ten_and_five_percent(
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

    expected_column_gap = max(captured_radii) * 0.1
    expected_row_gap = max(captured_radii) * 0.05
    record_x = [
        _extract_group_translate_xy(root, f"record_{index}")[0]
        for index in range(len(records))
    ]
    record_y = [
        _extract_group_translate_xy(root, f"record_{index}")[1]
        for index in range(len(records))
    ]

    assert record_x[1] - (record_x[0] + planned_sizes[0][0]) == pytest.approx(expected_column_gap, abs=1e-6)
    assert record_x[2] - (record_x[1] + planned_sizes[1][0]) == pytest.approx(expected_column_gap, abs=1e-6)
    assert record_x[4] - (record_x[3] + planned_sizes[3][0]) == pytest.approx(expected_column_gap, abs=1e-6)
    assert record_x[5] - (record_x[4] + planned_sizes[4][0]) == pytest.approx(expected_column_gap, abs=1e-6)
    assert record_y[3] - record_y[0] == pytest.approx(planned_sizes[0][1] + expected_row_gap, abs=1e-6)

    row0_width = planned_sizes[0][0] + planned_sizes[1][0] + planned_sizes[2][0] + 2.0 * expected_column_gap
    row1_width = planned_sizes[3][0] + planned_sizes[4][0] + planned_sizes[5][0] + 2.0 * expected_column_gap
    expected_row1_start = (row0_width - row1_width) * 0.5
    assert record_x[3] == pytest.approx(expected_row1_start, abs=1e-6)


@pytest.mark.circular
def test_multi_record_column_gap_ratio_override_controls_horizontal_spacing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("column_gap_a", 20, length=1000),
        _build_record("column_gap_b", 220, length=900),
        _build_record("column_gap_c", 420, length=800),
        _build_record("column_gap_d", 620, length=700),
    ]
    planned_sizes = [
        (1000.0, 900.0),
        (700.0, 900.0),
        (520.0, 900.0),
        (550.0, 900.0),
    ]
    planned_insets = [120.0, 85.0, 60.0, 70.0]
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
        legend="none",
        multi_record_column_gap_ratio=0.2,
    )
    root = ET.fromstring(canvas.tostring())

    expected_column_gap = max(captured_radii) * 0.2
    gap_row0 = _extract_record_axis_outer_edges(root, 1)[0] - _extract_record_axis_outer_edges(root, 0)[1]
    gap_row1 = _extract_record_axis_outer_edges(root, 3)[0] - _extract_record_axis_outer_edges(root, 2)[1]

    assert gap_row0 == pytest.approx(expected_column_gap, abs=1e-6)
    assert gap_row1 == pytest.approx(expected_column_gap, abs=1e-6)


@pytest.mark.circular
def test_multi_record_column_gap_ratio_zero_removes_gap_between_records_in_same_row(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("column_gap_zero_a", 20, length=1000),
        _build_record("column_gap_zero_b", 220, length=900),
        _build_record("column_gap_zero_c", 420, length=800),
        _build_record("column_gap_zero_d", 620, length=700),
    ]
    planned_sizes = [
        (1000.0, 900.0),
        (700.0, 900.0),
        (520.0, 900.0),
        (550.0, 900.0),
    ]
    planned_insets = [120.0, 85.0, 60.0, 70.0]

    def fake_single(gb_record: SeqRecord, **_kwargs: Any) -> Drawing:
        index_map = {record.id: idx for idx, record in enumerate(records)}
        index = index_map[gb_record.id]
        width, height = planned_sizes[index]
        return _build_mock_circular_subcanvas(width, height, planned_insets[index])

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        multi_record_column_gap_ratio=0.0,
    )
    root = ET.fromstring(canvas.tostring())

    gap_row0 = _extract_record_axis_outer_edges(root, 1)[0] - _extract_record_axis_outer_edges(root, 0)[1]
    gap_row1 = _extract_record_axis_outer_edges(root, 3)[0] - _extract_record_axis_outer_edges(root, 2)[1]

    assert gap_row0 == pytest.approx(0.0, abs=1e-6)
    assert gap_row1 == pytest.approx(0.0, abs=1e-6)


@pytest.mark.circular
def test_multi_record_row_gap_ratio_override_accepts_legacy_ten_percent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("row_gap_a", 20, length=1000),
        _build_record("row_gap_b", 220, length=900),
        _build_record("row_gap_c", 420, length=800),
        _build_record("row_gap_d", 620, length=700),
    ]
    planned_sizes = [
        (1000.0, 900.0),
        (700.0, 900.0),
        (520.0, 900.0),
        (550.0, 900.0),
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
        multi_record_row_gap_ratio=0.1,
    )
    root = ET.fromstring(canvas.tostring())

    expected_row_gap = max(captured_radii) * 0.1
    record_y = [
        _extract_group_translate_xy(root, f"record_{index}")[1]
        for index in range(len(records))
    ]

    assert record_y[2] - record_y[0] == pytest.approx(planned_sizes[0][1] + expected_row_gap, abs=1e-6)


@pytest.mark.circular
def test_multi_record_row_gap_ratio_zero_removes_visible_gap_between_row_content(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("visible_gap_a", 20, length=1000),
        _build_record("visible_gap_b", 220, length=900),
        _build_record("visible_gap_c", 420, length=800),
        _build_record("visible_gap_d", 620, length=700),
    ]

    def fake_single(gb_record: SeqRecord, **_kwargs: Any) -> Drawing:
        _ = gb_record
        return _build_mock_circular_subcanvas(600.0, 1000.0, 0.0)

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        multi_record_row_gap_ratio=0.0,
    )
    root = ET.fromstring(canvas.tostring())

    row0_bottom = max(
        _extract_record_axis_vertical_edges(root, 0)[1],
        _extract_record_axis_vertical_edges(root, 1)[1],
    )
    row1_top = min(
        _extract_record_axis_vertical_edges(root, 2)[0],
        _extract_record_axis_vertical_edges(root, 3)[0],
    )

    assert row1_top - row0_bottom == pytest.approx(0.0, abs=1e-6)


@pytest.mark.circular
def test_multi_record_row_gap_ratio_adds_visible_gap_between_row_content(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("visible_gap_ratio_a", 20, length=1000),
        _build_record("visible_gap_ratio_b", 220, length=900),
        _build_record("visible_gap_ratio_c", 420, length=800),
        _build_record("visible_gap_ratio_d", 620, length=700),
    ]
    captured_radii: list[float] = []

    def fake_single(gb_record: SeqRecord, **kwargs: Any) -> Drawing:
        _ = gb_record
        captured_radii.append(float(kwargs["cfg"].canvas.circular.radius))
        return _build_mock_circular_subcanvas(600.0, 1000.0, 0.0)

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        multi_record_row_gap_ratio=0.1,
    )
    root = ET.fromstring(canvas.tostring())

    row0_bottom = max(
        _extract_record_axis_vertical_edges(root, 0)[1],
        _extract_record_axis_vertical_edges(root, 1)[1],
    )
    row1_top = min(
        _extract_record_axis_vertical_edges(root, 2)[0],
        _extract_record_axis_vertical_edges(root, 3)[0],
    )
    expected_visible_gap = max(captured_radii) * 0.1

    assert row1_top - row0_bottom == pytest.approx(expected_visible_gap, abs=1e-6)


@pytest.mark.circular
def test_multi_record_positions_group_rows_and_preserve_within_row_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("row_pattern_a", 20, length=1000),
        _build_record("row_pattern_b", 220, length=900),
        _build_record("row_pattern_c", 420, length=800),
        _build_record("row_pattern_d", 620, length=700),
        _build_record("row_pattern_e", 820, length=600),
        _build_record("row_pattern_f", 1020, length=500),
    ]

    def fake_single(gb_record: SeqRecord, **_kwargs: Any) -> Drawing:
        canvas = Drawing(
            filename=f"{gb_record.id}.svg",
            size=("400px", "300px"),
            viewBox="0 0 400 300",
            debug=False,
        )
        marker = canvas.g(id=f"marker_{gb_record.id}")
        marker.add(canvas.circle(center=(10, 10), r=1))
        canvas.add(marker)
        return canvas

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        multi_record_positions=["#2@2", "#5@1", "#1@1", "#6@2", "#3@1", "#4@2"],
    )
    root = ET.fromstring(canvas.tostring())
    ns = {"svg": "http://www.w3.org/2000/svg"}
    record_y = [
        _extract_group_translate_xy(root, f"record_{index}")[1]
        for index in range(len(records))
    ]
    expected_ids = [
        "row_pattern_e",
        "row_pattern_a",
        "row_pattern_c",
        "row_pattern_b",
        "row_pattern_f",
        "row_pattern_d",
    ]

    for index, expected in enumerate(expected_ids):
        record_group = root.find(f".//svg:g[@id='record_{index}']", ns)
        assert record_group is not None
        marker_group = record_group.find(f"./svg:g[@id='marker_{expected}']", ns)
        assert marker_group is not None
    assert record_y[0] == pytest.approx(record_y[1], abs=1e-6)
    assert record_y[0] == pytest.approx(record_y[2], abs=1e-6)
    assert record_y[3] == pytest.approx(record_y[4], abs=1e-6)
    assert record_y[3] == pytest.approx(record_y[5], abs=1e-6)
    assert record_y[3] > record_y[0]


@pytest.mark.circular
def test_multi_record_positions_compress_row_gaps(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("row_gap_a", 20, length=1000),
        _build_record("row_gap_b", 220, length=900),
        _build_record("row_gap_c", 420, length=800),
        _build_record("row_gap_d", 620, length=700),
        _build_record("row_gap_e", 820, length=600),
    ]

    def fake_single(gb_record: SeqRecord, **_kwargs: Any) -> Drawing:
        _ = gb_record
        return Drawing(
            filename="row_gap.svg",
            size=("400px", "300px"),
            viewBox="0 0 400 300",
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
        multi_record_positions=["#1@1", "#2@1", "#3@3", "#4@3", "#5@3"],
    )
    root = ET.fromstring(canvas.tostring())
    record_y = [
        _extract_group_translate_xy(root, f"record_{index}")[1]
        for index in range(len(records))
    ]

    assert record_y[0] == pytest.approx(record_y[1], abs=1e-6)
    assert record_y[2] == pytest.approx(record_y[3], abs=1e-6)
    assert record_y[2] == pytest.approx(record_y[4], abs=1e-6)
    assert record_y[2] > record_y[0]


@pytest.mark.circular
def test_multi_record_positions_default_layout_keeps_auto_square_when_unspecified(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("row_auto_a", 20, length=1000),
        _build_record("row_auto_b", 220, length=900),
        _build_record("row_auto_c", 420, length=800),
        _build_record("row_auto_d", 620, length=700),
        _build_record("row_auto_e", 820, length=600),
    ]

    def fake_single(gb_record: SeqRecord, **_kwargs: Any) -> Drawing:
        _ = gb_record
        return Drawing(
            filename="row_auto.svg",
            size=("400px", "300px"),
            viewBox="0 0 400 300",
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
    record_y = [
        _extract_group_translate_xy(root, f"record_{index}")[1]
        for index in range(len(records))
    ]

    assert record_y[0] == pytest.approx(record_y[1], abs=1e-6)
    assert record_y[0] == pytest.approx(record_y[2], abs=1e-6)
    assert record_y[3] == pytest.approx(record_y[4], abs=1e-6)
    assert record_y[3] > record_y[0]


@pytest.mark.circular
def test_multi_record_positions_accept_record_id_selectors(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("order_a", 20, length=1000),
        _build_record("order_b", 220, length=900),
        _build_record("order_c", 420, length=800),
        _build_record("order_d", 620, length=700),
    ]

    def fake_single(gb_record: SeqRecord, **_kwargs: Any) -> Drawing:
        canvas = Drawing(
            filename=f"{gb_record.id}.svg",
            size=("400px", "300px"),
            viewBox="0 0 400 300",
            debug=False,
        )
        marker = canvas.g(id=f"marker_{gb_record.id}")
        marker.add(canvas.circle(center=(10, 10), r=1))
        canvas.add(marker)
        return canvas

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        multi_record_positions=["order_c@1", "order_a@2", "order_b@2", "order_d@1"],
    )
    root = ET.fromstring(canvas.tostring())
    ns = {"svg": "http://www.w3.org/2000/svg"}
    expected_ids = ["order_c", "order_d", "order_a", "order_b"]

    for index, expected in enumerate(expected_ids):
        record_group = root.find(f".//svg:g[@id='record_{index}']", ns)
        assert record_group is not None
        marker_group = record_group.find(f"./svg:g[@id='marker_{expected}']", ns)
        assert marker_group is not None


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
    row0_gap_01 = _extract_record_axis_outer_edges(root, 1)[0] - _extract_record_axis_outer_edges(root, 0)[1]
    row0_gap_12 = _extract_record_axis_outer_edges(root, 2)[0] - _extract_record_axis_outer_edges(root, 1)[1]
    row1_gap_34 = _extract_record_axis_outer_edges(root, 4)[0] - _extract_record_axis_outer_edges(root, 3)[1]
    row1_gap_45 = _extract_record_axis_outer_edges(root, 5)[0] - _extract_record_axis_outer_edges(root, 4)[1]

    assert row0_gap_01 == pytest.approx(expected_gap, abs=1e-6)
    assert row0_gap_12 == pytest.approx(expected_gap, abs=1e-6)
    assert row1_gap_34 == pytest.approx(expected_gap, abs=1e-6)
    assert row1_gap_45 == pytest.approx(expected_gap, abs=1e-6)

    content_widths = [
        planned_sizes[index][0] - (2.0 * planned_insets[index])
        for index in range(len(planned_sizes))
    ]
    row0_content_width = (
        content_widths[0] + content_widths[1] + content_widths[2] + 2.0 * expected_gap
    )
    row1_content_width = (
        content_widths[3] + content_widths[4] + content_widths[5] + 2.0 * expected_gap
    )
    row0_left = planned_insets[0]
    row0_right = planned_insets[2]
    row1_left = planned_insets[3]
    row1_right = planned_insets[5]
    row0_physical_width = row0_left + row0_content_width + row0_right
    row1_physical_width = row1_left + row1_content_width + row1_right
    row0_target = max(row0_left, row0_right)
    row1_target = max(row1_left, row1_right)
    row0_extra_left = row0_target - row0_left
    row0_extra_right = row0_target - row0_right
    row1_extra_left = row1_target - row1_left
    row1_extra_right = row1_target - row1_right
    row0_total_width = row0_physical_width + row0_extra_left + row0_extra_right
    row1_total_width = row1_physical_width + row1_extra_left + row1_extra_right
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
    row0_gap_01 = _extract_record_axis_outer_edges(root, 1)[0] - _extract_record_axis_outer_edges(root, 0)[1]
    row0_gap_12 = _extract_record_axis_outer_edges(root, 2)[0] - _extract_record_axis_outer_edges(root, 1)[1]
    row1_gap_34 = _extract_record_axis_outer_edges(root, 4)[0] - _extract_record_axis_outer_edges(root, 3)[1]
    row1_gap_45 = _extract_record_axis_outer_edges(root, 5)[0] - _extract_record_axis_outer_edges(root, 4)[1]

    assert row0_gap_01 == pytest.approx(expected_gap, abs=1e-6)
    assert row0_gap_12 == pytest.approx(expected_gap, abs=1e-6)
    assert row1_gap_34 == pytest.approx(expected_gap, abs=1e-6)
    assert row1_gap_45 == pytest.approx(expected_gap, abs=1e-6)

    content_widths = [
        planned_sizes[index][0] - (2.0 * planned_insets[index])
        for index in range(len(planned_sizes))
    ]
    row0_content_width = (
        content_widths[0] + content_widths[1] + content_widths[2] + 2.0 * expected_gap
    )
    row1_content_width = (
        content_widths[3] + content_widths[4] + content_widths[5] + 2.0 * expected_gap
    )
    row0_physical_width = planned_insets[0] + row0_content_width + planned_insets[2]
    row1_physical_width = planned_insets[3] + row1_content_width + planned_insets[5]
    expected_row1_start_without_symmetry = (row0_physical_width - row1_physical_width) * 0.5
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
    assert captured_kwargs["multi_record_size_mode"] == "auto"
    assert captured_kwargs["multi_record_min_radius_ratio"] == pytest.approx(0.55)
    assert captured_kwargs["multi_record_column_gap_ratio"] == pytest.approx(0.10)
    assert captured_kwargs["multi_record_row_gap_ratio"] == pytest.approx(0.05)
    assert captured_kwargs["plot_title_position"] == "none"


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
            "--multi_record_column_gap_ratio",
            "0.2",
            "--multi_record_row_gap_ratio",
            "0.12",
            "--plot_title_position",
            "top",
            "--plot_title_font_size",
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
    assert captured_kwargs["multi_record_column_gap_ratio"] == pytest.approx(0.2)
    assert captured_kwargs["multi_record_row_gap_ratio"] == pytest.approx(0.12)
    assert captured_kwargs["plot_title_position"] == "top"
    assert captured_kwargs["cfg"].objects.definition.circular.plot_title_font_size == pytest.approx(30.0)


@pytest.mark.circular
def test_circular_cli_multi_record_canvas_accepts_sqrt_alias(
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
            "sqrt",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert calls["single"] == 0
    assert calls["multi"] == 1
    assert calls["save"] == 1
    assert captured_kwargs["multi_record_size_mode"] == "sqrt"


@pytest.mark.circular
def test_circular_cli_multi_record_canvas_passes_positions(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    records = [
        _build_record("cli_a", 20),
        _build_record("cli_b", 220),
        _build_record("cli_c", 420),
    ]
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
            "--multi_record_position",
            "#3@1",
            "--multi_record_position",
            "cli_a@2",
            "--multi_record_position",
            "#2@2",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert calls["single"] == 0
    assert calls["multi"] == 1
    assert calls["save"] == 1
    assert captured_kwargs["multi_record_positions"] == ["#3@1", "cli_a@2", "#2@2"]


@pytest.mark.circular
@pytest.mark.parametrize(
    "position",
    ["#1@", "@1", "#1@0", "#1@-1", "#1@abc", "#0@1", "none@1", "selector_only"],
)
def test_circular_cli_rejects_invalid_multi_record_position(
    monkeypatch: pytest.MonkeyPatch,
    position: str,
) -> None:
    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda *_args, **_kwargs: [_build_record("cli_a", 20)])
    with pytest.raises(SystemExit):
        circular_cli_module.circular_main(
            [
                "--gbk",
                "dummy.gb",
                "--format",
                "svg",
                "--multi_record_position",
                position,
            ]
        )


@pytest.mark.circular
@pytest.mark.parametrize(
    ("removed_option", "value"),
    [("--multi_record_row_pattern", "2,4"), ("--multi_record_order", "#2")],
)
def test_circular_cli_rejects_removed_multi_record_layout_options(
    monkeypatch: pytest.MonkeyPatch,
    removed_option: str,
    value: str,
) -> None:
    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda *_args, **_kwargs: [_build_record("cli_a", 20)])
    with pytest.raises(SystemExit):
        circular_cli_module.circular_main(
            [
                "--gbk",
                "dummy.gb",
                "--format",
                "svg",
                removed_option,
                value,
            ]
        )


@pytest.mark.circular
def test_multi_record_positions_invalid_selector_raises_validation_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [_build_record("order_err_a", 20, length=1000), _build_record("order_err_b", 220, length=900)]

    def fake_single(gb_record: SeqRecord, **_kwargs: Any) -> Drawing:
        _ = gb_record
        return Drawing(
            filename="order_err.svg",
            size=("400px", "300px"),
            viewBox="0 0 400 300",
            debug=False,
        )

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    with pytest.raises(ValidationError):
        assemble_circular_diagram_from_records(
            records,
            selected_features_set=["CDS"],
            legend="none",
            multi_record_positions=["none@1", "#2@2"],
        )


@pytest.mark.circular
def test_multi_record_positions_with_duplicate_selector_raises_validation_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("row_err_a", 20, length=1000),
        _build_record("row_err_b", 220, length=900),
        _build_record("row_err_c", 420, length=800),
    ]

    def fake_single(gb_record: SeqRecord, **_kwargs: Any) -> Drawing:
        _ = gb_record
        return Drawing(
            filename="row_err.svg",
            size=("400px", "300px"),
            viewBox="0 0 400 300",
            debug=False,
        )

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    with pytest.raises(ValidationError, match="specified more than once"):
        assemble_circular_diagram_from_records(
            records,
            selected_features_set=["CDS"],
            legend="none",
            multi_record_positions=["#1@1", "#1@2", "#3@2"],
        )


@pytest.mark.circular
def test_multi_record_positions_with_missing_record_raises_validation_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("missing_row_a", 20, length=1000),
        _build_record("missing_row_b", 220, length=900),
        _build_record("missing_row_c", 420, length=800),
    ]

    def fake_single(gb_record: SeqRecord, **_kwargs: Any) -> Drawing:
        _ = gb_record
        return Drawing(
            filename="missing_row.svg",
            size=("400px", "300px"),
            viewBox="0 0 400 300",
            debug=False,
        )

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    with pytest.raises(ValidationError, match="must include each loaded record exactly once"):
        assemble_circular_diagram_from_records(
            records,
            selected_features_set=["CDS"],
            legend="none",
            multi_record_positions=["#1@1", "#2@2"],
        )


@pytest.mark.circular
def test_multi_record_positions_with_non_positive_row_raises_validation_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _build_record("row_err_a", 20, length=1000),
        _build_record("row_err_b", 220, length=900),
    ]

    def fake_single(gb_record: SeqRecord, **_kwargs: Any) -> Drawing:
        _ = gb_record
        return Drawing(
            filename="row_err.svg",
            size=("400px", "300px"),
            viewBox="0 0 400 300",
            debug=False,
        )

    monkeypatch.setattr(
        diagram_api_module,
        "assemble_circular_diagram_from_record",
        fake_single,
    )

    with pytest.raises(ValidationError, match="positive integer row"):
        assemble_circular_diagram_from_records(
            records,
            selected_features_set=["CDS"],
            legend="none",
            multi_record_positions=["#1@0", "#2@1"],
        )


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
            "--plot_title",
            "CLI Shared Title",
            "--plot_title_position",
            "top",
            "--plot_title_font_size",
            "30",
            "--keep_full_definition_with_plot_title",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert calls["single"] == len(records)
    assert calls["multi"] == 0
    assert calls["save"] == len(records)
    assert single_kwargs
    assert all(kwargs.get("plot_title") == "CLI Shared Title" for kwargs in single_kwargs)
    assert all(kwargs.get("plot_title_position") == "top" for kwargs in single_kwargs)
    assert all(kwargs.get("plot_title_font_size") == pytest.approx(30.0) for kwargs in single_kwargs)
    assert all(kwargs.get("keep_full_definition_with_plot_title") is True for kwargs in single_kwargs)


@pytest.mark.circular
def test_circular_cli_multi_record_canvas_passes_keep_full_definition_option(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    records = [_build_record("cli_a", 20), _build_record("cli_b", 220)]
    captured_kwargs: dict[str, Any] = {}

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda *_args, **_kwargs: records)
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "read_feature_visibility_file", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(
        circular_cli_module,
        "assemble_circular_diagram_from_record",
        lambda *_args, **_kwargs: Drawing(filename=str(tmp_path / "single.svg")),
    )

    def fake_multi(*_args: Any, **_kwargs: Any) -> Drawing:
        captured_kwargs.update(_kwargs)
        return Drawing(filename=str(tmp_path / "multi.svg"))

    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_records", fake_multi)
    monkeypatch.setattr(circular_cli_module, "save_figure", lambda *_args, **_kwargs: None)

    circular_cli_module.circular_main(
        [
            "--gbk",
            "dummy.gb",
            "--format",
            "svg",
            "--multi_record_canvas",
            "--keep_full_definition_with_plot_title",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured_kwargs["keep_full_definition_with_plot_title"] is True


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
@pytest.mark.parametrize("ratio", ["-0.1", "nan", "inf"])
def test_circular_cli_rejects_invalid_multi_record_column_gap_ratio(ratio: str) -> None:
    with pytest.raises(SystemExit):
        circular_cli_module._get_args(
            [
                "--gbk",
                "dummy.gb",
                "--multi_record_column_gap_ratio",
                ratio,
            ]
        )


@pytest.mark.circular
@pytest.mark.parametrize("ratio", ["-0.1", "nan", "inf"])
def test_circular_cli_rejects_invalid_multi_record_row_gap_ratio(ratio: str) -> None:
    with pytest.raises(SystemExit):
        circular_cli_module._get_args(
            [
                "--gbk",
                "dummy.gb",
                "--multi_record_row_gap_ratio",
                ratio,
            ]
        )


@pytest.mark.circular
def test_circular_cli_definition_layout_defaults() -> None:
    args = circular_cli_module._get_args(["--gbk", "dummy.gb"])
    assert args.plot_title is None
    assert args.multi_record_column_gap_ratio == pytest.approx(0.10)
    assert args.multi_record_row_gap_ratio == pytest.approx(0.05)
    assert args.plot_title_position == "none"
    assert args.plot_title_font_size is None
    assert args.keep_full_definition_with_plot_title is False


@pytest.mark.circular
def test_circular_cli_parses_keep_full_definition_option() -> None:
    args = circular_cli_module._get_args(
        ["--gbk", "dummy.gb", "--keep_full_definition_with_plot_title"]
    )
    assert args.keep_full_definition_with_plot_title is True


@pytest.mark.circular
@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("--plot_title_position", "center"),
        ("--plot_title_position", "left"),
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
@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("--definition_position", "bottom"),
        ("--multi_record_title_mode", "shared"),
        ("--multi_record_definition_mode", "shared"),
        ("--shared_definition_position", "bottom"),
        ("--shared_definition_font_size", "30"),
    ],
)
def test_circular_cli_rejects_legacy_shared_definition_options(
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
def test_multi_record_default_hides_plot_title_and_keeps_record_summary_content() -> None:
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

    shared_groups = root.findall(".//svg:g[@id='plot_title']", ns)
    assert len(shared_groups) == 0

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
def test_multi_record_shared_record_summary_includes_replicon_when_available() -> None:
    records = [
        _build_record_with_source(
            "rec_chr",
            organism="Organism chr",
            strain="Strain chr",
            feature_start=20,
            length=1500,
            chromosome="1",
        ),
        _build_record_with_source(
            "rec_plasmid",
            organism="Organism plasmid",
            strain="Strain plasmid",
            feature_start=240,
            length=900,
            plasmid="pVNTG2",
        ),
        _build_record_with_source(
            "rec_no_replicon",
            organism="Organism none",
            strain="Strain none",
            feature_start=120,
            length=1000,
        ),
    ]

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
    )
    root = ET.fromstring(canvas.tostring())

    expected_replicon_by_record = {
        "rec_chr": "Chromosome 1",
        "rec_plasmid": "pVNTG2",
    }
    for record in records:
        group_id = f"{record.id}_definition"
        record_texts = _extract_group_texts(root, group_id)
        assert record_texts
        assert any(text.endswith("bp") for text in record_texts)
        assert any(text.endswith("% GC") for text in record_texts)

        expected_replicon = expected_replicon_by_record.get(record.id)
        if expected_replicon is not None:
            assert record_texts[0] == expected_replicon
            assert record_texts[1] == record.id
            record_font_weights = _extract_group_font_weights(root, group_id)
            assert record_font_weights[0] == "bold"
        else:
            # No replicon qualifier: keep current summary fallback (accession first).
            assert record_texts[0] == record.id


@pytest.mark.circular
def test_plot_title_font_size_override_only_changes_plot_title() -> None:
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
        plot_title_position="top",
        config_overrides={"plot_title_font_size": 30},
    )
    root = ET.fromstring(canvas.tostring())

    plot_title_font_sizes = _extract_group_font_sizes(root, "plot_title")
    assert plot_title_font_sizes == [30.0]
    shared_font_weights = _extract_group_font_weights(root, "plot_title")
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
def test_plot_title_single_line_with_missing_species_or_strain(
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
        plot_title_position="top",
    )
    root = ET.fromstring(canvas.tostring())

    shared_texts = _extract_group_texts(root, "plot_title")
    assert shared_texts == [expected_shared]


@pytest.mark.circular
def test_plot_title_position_moves_shared_group_vertically() -> None:
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

    none_canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        plot_title_position="none",
    )
    top_canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        plot_title_position="top",
    )
    bottom_canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        plot_title_position="bottom",
    )

    none_root = ET.fromstring(none_canvas.tostring())
    ns = {"svg": "http://www.w3.org/2000/svg"}
    assert none_root.findall(".//svg:g[@id='plot_title']", ns) == []

    y_top = _extract_group_translate_y(ET.fromstring(top_canvas.tostring()), "plot_title")
    y_bottom = _extract_group_translate_y(
        ET.fromstring(bottom_canvas.tostring()),
        "plot_title",
    )

    assert y_top < y_bottom


@pytest.mark.circular
def test_single_record_bottom_plot_title_uses_summary_center_definition() -> None:
    record = _build_record_with_source(
        "single_shared_bottom",
        organism="Single shared organism",
        strain="Single shared strain",
        length=1200,
        chromosome="1",
    )

    canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS"],
        legend="none",
        plot_title_position="bottom",
    )
    root = ET.fromstring(canvas.tostring())

    assert _extract_group_texts(root, "plot_title") == ["Single shared organism Single shared strain"]
    record_texts = _extract_group_texts(root, f"{record.id}_definition")
    assert record_texts[:2] == ["Chromosome 1", record.id]
    assert any(text.endswith("bp") for text in record_texts)
    assert any(text.endswith("% GC") for text in record_texts)
    assert "Single shared organism" not in record_texts
    assert "Single shared strain" not in record_texts

    axis_center_y, definition_center_y = _extract_single_axis_and_definition_center_y(
        root,
        record_id=record.id,
    )
    assert definition_center_y == pytest.approx(axis_center_y, abs=1e-6)


@pytest.mark.circular
def test_single_record_plot_title_can_keep_full_center_definition() -> None:
    record = _build_record_with_source(
        "single_shared_full_definition",
        organism="Single shared organism",
        strain="Single shared strain",
        length=1200,
        chromosome="1",
    )

    canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS"],
        legend="none",
        plot_title_position="top",
        keep_full_definition_with_plot_title=True,
    )
    root = ET.fromstring(canvas.tostring())

    assert _extract_group_texts(root, "plot_title") == ["Single shared organism Single shared strain"]
    record_texts = _extract_group_texts(root, f"{record.id}_definition")
    assert "Single shared organism" in record_texts
    assert "Single shared strain" in record_texts
    assert any(text.endswith("bp") for text in record_texts)
    assert any(text.endswith("% GC") for text in record_texts)

    axis_center_y, definition_center_y = _extract_single_axis_and_definition_center_y(
        root,
        record_id=record.id,
    )
    assert definition_center_y == pytest.approx(axis_center_y, abs=1e-6)


@pytest.mark.circular
def test_single_record_custom_plot_title_overrides_default_when_visible() -> None:
    record = _build_record_with_source(
        "single_custom_title",
        organism="Default single organism",
        strain="Default single strain",
        length=1200,
    )

    canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS"],
        legend="none",
        plot_title="Custom Circular Title",
        plot_title_position="top",
    )
    root = ET.fromstring(canvas.tostring())

    assert _extract_group_texts(root, "plot_title") == ["Custom Circular Title"]
    record_texts = _extract_group_texts(root, f"{record.id}_definition")
    assert record_texts[0] == record.id
    assert any(text.endswith("bp") for text in record_texts)
    assert any(text.endswith("% GC") for text in record_texts)
    assert "Default single organism" not in record_texts
    assert "Default single strain" not in record_texts


@pytest.mark.circular
def test_single_record_plot_title_keeps_plain_text_around_inline_italics() -> None:
    record = _build_record_with_source(
        "single_mixed_content_title",
        organism="Default single organism",
        strain="Default single strain",
        length=1200,
    )
    plot_title = (
        "Erythromycin A biosynthetic gene cluster from "
        "<i>Saccharopolyspora erythraea</i>"
    )

    canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS"],
        legend="none",
        plot_title=plot_title,
        plot_title_position="top",
    )
    root = ET.fromstring(canvas.tostring())

    assert _extract_group_texts(root, "plot_title") == [
        "Erythromycin A biosynthetic gene cluster from Saccharopolyspora erythraea"
    ]
    assert _extract_group_italic_tspan_texts(root, "plot_title") == ["Saccharopolyspora erythraea"]


@pytest.mark.circular
def test_single_record_hidden_plot_title_ignores_custom_title_and_keeps_full_definition() -> None:
    record = _build_record_with_source(
        "single_layout",
        organism="Single organism",
        strain="Single strain",
        length=1200,
    )

    canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS"],
        legend="none",
        plot_title="Hidden Circular Title",
    )
    root = ET.fromstring(canvas.tostring())
    ns = {"svg": "http://www.w3.org/2000/svg"}

    assert root.findall(".//svg:g[@id='plot_title']", ns) == []
    record_texts = _extract_group_texts(root, f"{record.id}_definition")
    assert "Single organism" in record_texts
    assert "Single strain" in record_texts


@pytest.mark.circular
@pytest.mark.parametrize(
    ("plot_title_position", "keep_full_definition_with_plot_title"),
    [
        ("none", False),
        ("top", False),
        ("bottom", False),
        ("top", True),
        ("bottom", True),
    ],
)
def test_multi_record_center_definition_aligns_with_record_axis(
    plot_title_position: str,
    keep_full_definition_with_plot_title: bool,
) -> None:
    records = [
        _build_record_with_source(
            "center_align_a",
            organism="Center align A",
            strain="Strain A",
            length=1400,
        ),
        _build_record_with_source(
            "center_align_b",
            organism="Center align B",
            strain="Strain B",
            length=900,
        ),
        _build_record_with_source(
            "center_align_c",
            organism="Center align C",
            strain="Strain C",
            length=700,
        ),
    ]

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        plot_title_position=plot_title_position,
        keep_full_definition_with_plot_title=keep_full_definition_with_plot_title,
    )
    root = ET.fromstring(canvas.tostring())

    for index, record in enumerate(records):
        axis_center_y, definition_center_y = _extract_record_axis_and_definition_center_y(
            root,
            record_index=index,
            record_id=record.id,
        )
        assert definition_center_y == pytest.approx(axis_center_y, abs=1e-6)


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
def test_multi_record_bottom_plot_title_stays_below_legend() -> None:
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
        plot_title_position="bottom",
        keep_full_definition_with_plot_title=True,
    )
    root = ET.fromstring(canvas.tostring())
    _legend_top, legend_bottom = _extract_legend_vertical_bounds(root)
    shared_top = _extract_definition_top_y(root, "plot_title")
    assert shared_top >= legend_bottom + 20.0 - 1e-6


@pytest.mark.circular
def test_single_record_bottom_plot_title_stays_below_legend() -> None:
    record = _build_record_with_source(
        "single_bottom_shared",
        organism="Single bottom shared",
        strain="Strain shared",
        length=1300,
    )

    canvas = diagram_api_module.assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS"],
        legend="bottom",
        plot_title_position="bottom",
    )
    root = ET.fromstring(canvas.tostring())
    _legend_top, legend_bottom = _extract_legend_vertical_bounds(root)
    shared_top = _extract_definition_top_y(root, "plot_title")
    assert shared_top >= legend_bottom + 20.0 - 1e-6


@pytest.mark.circular
def test_multi_record_custom_plot_title_overrides_default_shared_title() -> None:
    records = [
        _build_record_with_source(
            "custom_shared_a",
            organism="Custom shared A",
            strain="Strain A",
            length=1300,
        ),
        _build_record_with_source(
            "custom_shared_b",
            organism="Custom shared B",
            strain="Strain B",
            length=900,
        ),
    ]

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        plot_title="Custom Shared Plot Title",
        plot_title_position="top",
    )
    root = ET.fromstring(canvas.tostring())

    assert _extract_group_texts(root, "plot_title") == ["Custom Shared Plot Title"]


@pytest.mark.circular
@pytest.mark.parametrize("legend_position", ["left", "right"])
def test_multi_record_left_right_plot_title_bottom_keeps_margins(
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
        plot_title_position="none",
    )
    baseline_root = ET.fromstring(baseline_canvas.tostring())
    baseline_height = _extract_viewbox_height(baseline_root)

    bottom_canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend=legend_position,
        plot_title_position="bottom",
    )
    bottom_root = ET.fromstring(bottom_canvas.tostring())
    bottom_height = _extract_viewbox_height(bottom_root)
    shared_top = _extract_definition_top_y(bottom_root, "plot_title")
    shared_bottom = _extract_definition_bottom_y(bottom_root, "plot_title")

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
def test_multi_record_visible_plot_title_uses_default_shared_title_when_blank() -> None:
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
        plot_title_position="top",
    )
    root = ET.fromstring(canvas.tostring())
    assert _extract_group_texts(root, "plot_title") == ["Legacy organism A Legacy strain A"]

    for record in records:
        record_texts = _extract_group_texts(root, f"{record.id}_definition")
        assert record.features[0].qualifiers["organism"][0] not in record_texts
        assert record.features[0].qualifiers["strain"][0] not in record_texts


@pytest.mark.circular
def test_multi_record_plot_title_can_keep_full_definitions() -> None:
    records = [
        _build_record_with_source(
            "legacy_full_a",
            organism="Legacy organism A",
            strain="Legacy strain A",
            length=1400,
        ),
        _build_record_with_source(
            "legacy_full_b",
            organism="Legacy organism B",
            strain="Legacy strain B",
            length=800,
        ),
    ]

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        plot_title_position="top",
        keep_full_definition_with_plot_title=True,
    )
    root = ET.fromstring(canvas.tostring())

    assert _extract_group_texts(root, "plot_title") == ["Legacy organism A Legacy strain A"]
    for record in records:
        record_texts = _extract_group_texts(root, f"{record.id}_definition")
        assert record.features[0].qualifiers["organism"][0] in record_texts
        assert record.features[0].qualifiers["strain"][0] in record_texts


@pytest.mark.circular
def test_multi_record_hidden_plot_title_keeps_summary_when_keep_full_enabled() -> None:
    records = [
        _build_record_with_source(
            "legacy_none_a",
            organism="Legacy organism A",
            strain="Legacy strain A",
            length=1400,
        ),
        _build_record_with_source(
            "legacy_none_b",
            organism="Legacy organism B",
            strain="Legacy strain B",
            length=800,
        ),
    ]

    canvas = assemble_circular_diagram_from_records(
        records,
        selected_features_set=["CDS"],
        legend="none",
        plot_title_position="none",
        keep_full_definition_with_plot_title=True,
    )
    root = ET.fromstring(canvas.tostring())

    for record in records:
        record_texts = _extract_group_texts(root, f"{record.id}_definition")
        assert record.features[0].qualifiers["organism"][0] not in record_texts
        assert record.features[0].qualifiers["strain"][0] not in record_texts


@pytest.mark.circular
def test_build_circular_diagram_passes_plot_title_position_option(
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
            plot_title="Build Shared Title",
            plot_title_font_size=28,
            keep_full_definition_with_plot_title=True,
            output=OutputOptions(plot_title_position="bottom"),
        ),
    )

    assert captured_kwargs["plot_title"] == "Build Shared Title"
    assert captured_kwargs["plot_title_font_size"] == pytest.approx(28.0)
    assert captured_kwargs["plot_title_position"] == "bottom"
    assert captured_kwargs["keep_full_definition_with_plot_title"] is True
