from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.api.diagram as api_diagram_module
import gbdraw.linear as linear_cli_module
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.api.options import DiagramOptions
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import load_comparisons
from gbdraw.render.groups.linear.pairwise_match import PairWiseMatchGroup


def _build_record() -> SeqRecord:
    return SeqRecord(Seq("A" * 400), id="record_1")


def _write_blast_table(path: Path, rows: list[tuple[object, ...]]) -> None:
    contents = [
        "# Fields: query acc., subject acc., % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score"
    ]
    contents.extend("\t".join(map(str, row)) for row in rows)
    path.write_text("\n".join(contents) + "\n", encoding="utf-8")


def _build_pairwise_group(match_style: str = "ribbon") -> PairWiseMatchGroup:
    group = PairWiseMatchGroup.__new__(PairWiseMatchGroup)
    group.canvas_config = SimpleNamespace(
        normalize_length=False,
        alignment_width=1000,
        longest_genome=400,
    )
    group.records = [_build_record(), _build_record()]
    group.comparison_count = 1
    group.comparison_height = 40
    group.query_offset_x = 0
    group.subject_offset_x = 0
    group.query_alignment_offset_x = 0
    group.subject_alignment_offset_x = 0
    group.min_identity = 0
    group.match_min_color = "#ffffff"
    group.match_max_color = "#000000"
    group.match_fill_opacity = 0.75
    group.match_stroke_color = "none"
    group.match_stroke_width = 0
    group.match_style = match_style
    group.curve_tension = 0.5
    return group


def _build_match_row() -> SimpleNamespace:
    return SimpleNamespace(
        identity=95,
        qstart=20,
        qend=80,
        sstart=40,
        send=120,
    )


@pytest.mark.linear
def test_blast_match_config_defaults_pairwise_style_for_legacy_config() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    legacy_config = deepcopy(config_dict)
    legacy_config["objects"]["blast_match"].pop("style", None)
    legacy_config["objects"]["blast_match"].pop("curve_tension", None)

    cfg = GbdrawConfig.from_dict(legacy_config)

    assert cfg.objects.blast_match.style == "ribbon"
    assert cfg.objects.blast_match.curve_tension == 0.5


@pytest.mark.linear
def test_blast_match_config_rejects_invalid_pairwise_style() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict["objects"]["blast_match"]["style"] = "line"

    with pytest.raises(ValidationError, match="pairwise_match_style must be one of: ribbon, curve"):
        GbdrawConfig.from_dict(config_dict)


@pytest.mark.linear
def test_blast_match_config_rejects_invalid_curve_tension() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict["objects"]["blast_match"]["curve_tension"] = 1.5

    with pytest.raises(ValidationError, match="curve_tension"):
        GbdrawConfig.from_dict(config_dict)


@pytest.mark.linear
def test_modify_config_dict_maps_pairwise_match_style() -> None:
    config_dict = load_config_toml("gbdraw.data", "config.toml")

    modified = modify_config_dict(config_dict, pairwise_match_style="curve")

    assert GbdrawConfig.from_dict(modified).objects.blast_match.style == "curve"


@pytest.mark.linear
def test_load_comparisons_filters_hits_below_alignment_length(tmp_path: Path) -> None:
    blast_path = tmp_path / "comparisons.tsv"
    _write_blast_table(
        blast_path,
        [
            ("q1", "s1", 95, 80, 0, 0, 1, 80, 1, 80, 1e-40, 500),
            ("q1", "s1", 95, 150, 0, 0, 1, 150, 1, 150, 1e-40, 500),
        ],
    )

    config = SimpleNamespace(evalue=1e-2, bitscore=50, identity=0, alignment_length=100)

    comparisons = load_comparisons([str(blast_path)], config)

    assert len(comparisons) == 1
    assert comparisons[0]["alignment_length"].tolist() == [150]


@pytest.mark.linear
def test_load_comparisons_applies_alignment_length_with_other_thresholds(tmp_path: Path) -> None:
    blast_path = tmp_path / "comparisons.tsv"
    _write_blast_table(
        blast_path,
        [
            ("q1", "s1", 95, 90, 0, 0, 1, 90, 1, 90, 1e-40, 500),
            ("q1", "s1", 60, 180, 0, 0, 1, 180, 1, 180, 1e-40, 500),
            ("q1", "s1", 95, 180, 0, 0, 1, 180, 1, 180, 1e-40, 500),
        ],
    )

    config = SimpleNamespace(evalue=1e-10, bitscore=50, identity=90, alignment_length=100)

    comparisons = load_comparisons([str(blast_path)], config)

    assert len(comparisons) == 1
    assert comparisons[0].shape[0] == 1
    assert comparisons[0]["alignment_length"].tolist() == [180]
    assert comparisons[0]["identity"].tolist() == [95]


@pytest.mark.linear
def test_linear_cli_alignment_length_is_forwarded(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    record = _build_record()
    captured: dict[str, object] = {}

    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *_args, **_kwargs: [record])
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "read_feature_visibility_file", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "save_figure", lambda _canvas, _formats: None)

    def fake_assemble(*_args, **kwargs):
        captured["alignment_length"] = kwargs.get("alignment_length")
        captured["pairwise_match_style"] = kwargs.get("pairwise_match_style")
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(linear_cli_module, "assemble_linear_diagram_from_records", fake_assemble)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "dummy.gb",
            "--alignment_length",
            "123",
            "--pairwise_match_style",
            "curve",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["alignment_length"] == 123
    assert captured["pairwise_match_style"] == "curve"


@pytest.mark.linear
def test_linear_cli_rejects_invalid_pairwise_match_style(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit):
        linear_cli_module._get_args(["--gbk", "dummy.gb", "--pairwise_match_style", "line"])

    assert "pairwise_match_style must be one of: ribbon, curve" in capsys.readouterr().err


@pytest.mark.linear
def test_build_linear_diagram_forwards_alignment_length(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(*_args, **kwargs):
        captured["alignment_length"] = kwargs.get("alignment_length")
        captured["pairwise_match_style"] = kwargs.get("pairwise_match_style")
        return Drawing(filename="dummy.svg")

    monkeypatch.setattr(api_diagram_module, "assemble_linear_diagram_from_records", fake_assemble)

    canvas = api_diagram_module.build_linear_diagram(
        [_build_record()],
        options=DiagramOptions(alignment_length=321, pairwise_match_style="curve"),
    )

    assert isinstance(canvas, Drawing)
    assert captured["alignment_length"] == 321
    assert captured["pairwise_match_style"] == "curve"


@pytest.mark.linear
def test_assemble_linear_diagram_rejects_negative_alignment_length() -> None:
    with pytest.raises(ValidationError, match="alignment_length must be >= 0"):
        assemble_linear_diagram_from_records(
            [_build_record()],
            alignment_length=-1,
        )


@pytest.mark.linear
def test_assemble_linear_diagram_rejects_invalid_pairwise_match_style() -> None:
    with pytest.raises(ValidationError, match="pairwise_match_style must be one of: ribbon, curve"):
        assemble_linear_diagram_from_records(
            [_build_record()],
            pairwise_match_style="line",
        )


@pytest.mark.linear
def test_pairwise_match_ribbon_path_uses_straight_segments() -> None:
    path = _build_pairwise_group("ribbon").generate_linear_match_path(_build_match_row())
    path_desc = str(path.get_xml().attrib["d"])

    assert "L" in path_desc
    assert "C" not in path_desc
    assert path.attribs["data-pairwise-match-style"] == "ribbon"
    assert float(path.attribs["data-identity-factor"]) == pytest.approx(0.95)


@pytest.mark.linear
def test_pairwise_match_path_serializes_orthogroup_metadata_without_collinearity_block() -> None:
    row = SimpleNamespace(
        **_build_match_row().__dict__,
        orthogroup_id="og_1",
        query_protein_id="query_protein",
        subject_protein_id="subject_protein",
        query_feature_svg_id="fquery",
        subject_feature_svg_id="fsubject",
    )

    path = _build_pairwise_group("ribbon").generate_linear_match_path(row)

    assert path.attribs["data-orthogroup-id"] == "og_1"
    assert path.attribs["data-query-protein-id"] == "query_protein"
    assert path.attribs["data-subject-protein-id"] == "subject_protein"
    assert path.attribs["data-query-feature-svg-id"] == "fquery_record_1"
    assert path.attribs["data-subject-feature-svg-id"] == "fsubject_record_2"
    assert path.attribs["data-query-stable-feature-svg-id"] == "fquery"
    assert path.attribs["data-subject-stable-feature-svg-id"] == "fsubject"
    assert "data-collinearity-block-id" not in path.attribs


@pytest.mark.linear
def test_pairwise_match_curve_path_preserves_endpoint_spans() -> None:
    path = _build_pairwise_group("curve").generate_linear_match_path(_build_match_row())
    path_desc = str(path.get_xml().attrib["d"])

    assert "C" in path_desc
    for endpoint_x in ("50.0", "200.0", "100.0", "300.0"):
        assert endpoint_x in path_desc
    assert path_desc.strip().lower().endswith("z")
    assert path.attribs["data-pairwise-match-style"] == "curve"
    assert float(path.attribs["data-identity-factor"]) == pytest.approx(0.95)


@pytest.mark.linear
def test_pairwise_match_group_draws_smaller_high_identity_hits_above_broad_hits() -> None:
    rows = [
        {
            "identity": 95,
            "qstart": 60,
            "qend": 80,
            "sstart": 70,
            "send": 90,
        },
        {
            "identity": 10,
            "qstart": 1,
            "qend": 220,
            "sstart": 1,
            "send": 220,
        },
    ]
    comparison_df = pd.DataFrame.from_records(rows)
    canvas_config = SimpleNamespace(
        normalize_length=False,
        alignment_width=1000,
        longest_genome=400,
        align_center=False,
    )
    blast_config = SimpleNamespace(
        fill_color="#d3d3d3",
        min_color="#ffffff",
        max_color="#000000",
        fill_opacity=1.0,
        stroke_color="none",
        stroke_width=0,
        identity=0,
        sequence_length_dict={},
    )
    match_group = PairWiseMatchGroup(
        canvas_config,
        blast_config.sequence_length_dict,
        comparison_df,
        40,
        1,
        blast_config,
        [_build_record(), _build_record()],
    ).get_group()

    drawing = Drawing(debug=False)
    drawing.add(match_group)
    svg = drawing.tostring()

    assert svg.index('data-identity-factor="0.1"') < svg.index('data-identity-factor="0.95"')
