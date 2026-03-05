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
from gbdraw.api.diagram import assemble_circular_diagram_from_records
from gbdraw.api.options import DiagramOptions, OutputOptions
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
    assert shared_font_sizes == [64.0]
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
