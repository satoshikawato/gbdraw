from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.circular as circular_cli_module
import gbdraw.linear as linear_cli_module
from gbdraw.api.diagram import assemble_circular_diagram_from_record
from gbdraw.features.colors import compute_feature_hash
from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.factory import create_feature_dict
from gbdraw.features.objects import FeatureLocationPart, FeatureObject
from gbdraw.features.shapes import (
    DEFAULT_DIRECTIONAL_FEATURE_TYPES,
    parse_feature_shape_assignment,
    parse_feature_shape_overrides,
    resolve_directional_feature_types,
)
from gbdraw.io.colors import load_default_colors
from gbdraw.io.genome import load_gbks
from gbdraw.labels.filtering import preprocess_label_filtering
from gbdraw.render.drawers.circular.features import FeaturePathGenerator as CircularFeaturePathGenerator
from gbdraw.render.drawers.linear.features import FeaturePathGenerator as LinearFeaturePathGenerator


def _build_test_record() -> SeqRecord:
    record = SeqRecord(Seq("A" * 400), id="record_1")
    record.features = [
        SeqFeature(
            FeatureLocation(10, 90, strand=1),
            type="CDS",
            qualifiers={"product": ["dummy cds"]},
        ),
        SeqFeature(
            FeatureLocation(120, 180, strand=1),
            type="repeat_region",
            qualifiers={"note": ["dummy repeat"]},
        ),
    ]
    return record


def _build_test_feature(*, is_directional: bool) -> FeatureObject:
    location = [FeatureLocationPart("block", "001", "positive", 100, 200, True)]
    return FeatureObject(
        feature_id="feature_000000001",
        location=location,
        is_directional=is_directional,
        color="#cccccc",
        note="",
        label_text="",
        coordinates=location,
        type="CDS",
        qualifiers={},
    )


def test_parse_feature_shape_assignment_valid() -> None:
    assert parse_feature_shape_assignment("CDS=arrow") == ("CDS", "arrow")
    assert parse_feature_shape_assignment("repeat_region=RECTANGLE") == (
        "repeat_region",
        "rectangle",
    )


@pytest.mark.parametrize(
    "raw",
    ["CDSarrow", "=arrow", "CDS=", "CDS=triangle"],
)
def test_parse_feature_shape_assignment_invalid(raw: str) -> None:
    with pytest.raises(ValueError):
        parse_feature_shape_assignment(raw)


def test_parse_feature_shape_overrides_last_wins() -> None:
    overrides = parse_feature_shape_overrides(
        ["CDS=arrow", "repeat_region=rectangle", "CDS=rectangle"]
    )
    assert overrides == {"CDS": "rectangle", "repeat_region": "rectangle"}


def test_resolve_directional_feature_types_with_overrides() -> None:
    directional_types = resolve_directional_feature_types(
        {"CDS": "rectangle", "repeat_region": "arrow"}
    )
    assert "CDS" not in directional_types
    assert "repeat_region" in directional_types
    assert "rRNA" in directional_types


def test_create_feature_dict_applies_directional_feature_types() -> None:
    default_colors = load_default_colors("", "default")
    color_table, default_colors = preprocess_color_tables(None, default_colors)
    label_filtering = preprocess_label_filtering(
        {"blacklist_keywords": [], "whitelist_df": None, "qualifier_priority_df": None}
    )
    record = _build_test_record()
    directional_types = set(DEFAULT_DIRECTIONAL_FEATURE_TYPES)
    directional_types.discard("CDS")
    directional_types.add("repeat_region")

    feature_dict, _ = create_feature_dict(
        record,
        color_table,
        ["CDS", "repeat_region"],
        default_colors,
        separate_strands=False,
        resolve_overlaps=False,
        label_filtering=label_filtering,
        directional_feature_types=directional_types,
    )

    cds_objects = [feature for feature in feature_dict.values() if feature.feature_type == "CDS"]
    repeat_objects = [
        feature for feature in feature_dict.values() if feature.feature_type == "repeat_region"
    ]
    assert cds_objects and repeat_objects
    assert cds_objects[0].is_directional is False
    assert repeat_objects[0].is_directional is True


def test_linear_feature_paths_change_with_directionality() -> None:
    generator = LinearFeaturePathGenerator(
        genome_length=1000,
        alignment_width=800,
        cds_height=20,
        genome_size_normalization_factor=1.0,
        feature_strand="positive",
        separate_strands=False,
        arrow_length=50,
        track_layout="middle",
        track_axis_gap=None,
    )
    directional_paths = generator.generate_linear_gene_path(
        _build_test_feature(is_directional=True)
    )
    rectangle_paths = generator.generate_linear_gene_path(
        _build_test_feature(is_directional=False)
    )
    assert directional_paths[0][0] == "block"
    assert rectangle_paths[0][0] == "block"
    assert directional_paths[0][1] != rectangle_paths[0][1]


def test_circular_feature_paths_change_with_directionality() -> None:
    generator = CircularFeaturePathGenerator(
        radius=250.0,
        total_length=1000,
        track_ratio=0.19,
        cds_ratio=0.03,
        offset=0.0,
        track_type="tuckin",
        strandedness=True,
        track_id=0,
    )
    directional_paths = generator.generate_circular_gene_path(
        _build_test_feature(is_directional=True)
    )
    rectangle_paths = generator.generate_circular_gene_path(
        _build_test_feature(is_directional=False)
    )
    assert directional_paths[0][0] == "block"
    assert rectangle_paths[0][0] == "block"
    assert directional_paths[0][1] != rectangle_paths[0][1]


def test_hmmtdna_origin_spanning_d_loop_arrow_renders_as_single_block() -> None:
    record = load_gbks([str(Path("tests/test_inputs/HmmtDNA.gbk"))], "circular")[0]
    d_loop_feature = next(feature for feature in record.features if feature.type == "D-loop")
    d_loop_id = compute_feature_hash(d_loop_feature, record_id=record.id)

    canvas = assemble_circular_diagram_from_record(
        record,
        selected_features_set=["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "repeat_region", "D-loop"],
        feature_shapes={"D-loop": "arrow"},
        legend="none",
    )

    root = ET.fromstring(canvas.tostring())
    ns = {"svg": "http://www.w3.org/2000/svg"}
    d_loop_paths = [
        path for path in root.findall(".//svg:path", ns) if path.attrib.get("id") == d_loop_id
    ]

    assert len(d_loop_paths) == 1


def test_circular_cli_feature_shape_forwards(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    record = _build_test_record()
    captured: dict[str, Any] = {}

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda *_args, **_kwargs: [record])
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda _path, _palette: None)
    monkeypatch.setattr(circular_cli_module, "save_figure", lambda _canvas, _formats: None)

    def fake_assemble(*_args, **kwargs):
        captured["feature_shapes"] = kwargs.get("feature_shapes")
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_record", fake_assemble)

    circular_cli_module.circular_main(
        [
            "--gbk",
            "dummy.gb",
            "--feature_shape",
            "CDS=rectangle",
            "--feature_shape",
            "repeat_region=arrow",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )
    assert captured["feature_shapes"] == {"CDS": "rectangle", "repeat_region": "arrow"}


def test_linear_cli_feature_shape_forwards(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    record = _build_test_record()
    captured: dict[str, Any] = {}

    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *_args, **_kwargs: [record])
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "save_figure", lambda _canvas, _formats: None)

    def fake_assemble(*_args, **kwargs):
        captured["feature_shapes"] = kwargs.get("feature_shapes")
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(linear_cli_module, "assemble_linear_diagram_from_records", fake_assemble)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "dummy.gb",
            "--feature_shape",
            "CDS=rectangle",
            "--feature_shape",
            "repeat_region=arrow",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )
    assert captured["feature_shapes"] == {"CDS": "rectangle", "repeat_region": "arrow"}


@pytest.mark.parametrize(
    "cmd_args",
    [
        ["--gbk", "dummy.gb", "--feature_shape", "CDSarrow"],
        ["--gbk", "dummy.gb", "--feature_shape", "=arrow"],
        ["--gbk", "dummy.gb", "--feature_shape", "CDS=triangle"],
    ],
)
def test_circular_cli_feature_shape_validation(cmd_args: list[str]) -> None:
    with pytest.raises(SystemExit):
        circular_cli_module._get_args(cmd_args)


@pytest.mark.parametrize(
    "cmd_args",
    [
        ["--gbk", "dummy.gb", "--feature_shape", "CDSarrow"],
        ["--gbk", "dummy.gb", "--feature_shape", "=arrow"],
        ["--gbk", "dummy.gb", "--feature_shape", "CDS=triangle"],
    ],
)
def test_linear_cli_feature_shape_validation(cmd_args: list[str]) -> None:
    with pytest.raises(SystemExit):
        linear_cli_module._get_args(cmd_args)
