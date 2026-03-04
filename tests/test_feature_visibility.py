from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.circular as circular_cli_module
import gbdraw.linear as linear_cli_module
from gbdraw.core.sequence import check_feature_presence
from gbdraw.exceptions import InputFileError, ParseError, ValidationError
from gbdraw.features.colors import compute_feature_hash, precompute_used_color_rules, preprocess_color_tables
from gbdraw.features.factory import create_feature_dict
from gbdraw.features.visibility import (
    compile_feature_visibility_rules,
    read_feature_visibility_file,
    should_render_feature,
)
from gbdraw.io.colors import load_default_colors
from gbdraw.labels.filtering import preprocess_label_filtering


def _visibility_df(rows: list[list[str]]) -> pd.DataFrame:
    return pd.DataFrame(
        rows,
        columns=["record_id", "feature_type", "qualifier", "value", "action"],
    )


def _make_record() -> SeqRecord:
    record = SeqRecord(Seq("A" * 500), id="rec1")
    record.features = [
        SeqFeature(
            FeatureLocation(10, 90, strand=1),
            type="CDS",
            qualifiers={
                "gene": ["geneA"],
                "product": ["enzyme alpha"],
                "note": ["core enzyme"],
            },
        ),
        SeqFeature(
            FeatureLocation(150, 210, strand=1),
            type="misc_feature",
            qualifiers={
                "gene": ["insA"],
                "note": ["transposase island"],
            },
        ),
    ]
    return record


def _base_label_filtering() -> dict[str, Any]:
    return {
        "blacklist_keywords": [],
        "whitelist_df": None,
        "qualifier_priority_df": None,
    }


def test_read_feature_visibility_file_ok(tmp_path: Path) -> None:
    table = tmp_path / "feature_visibility.tsv"
    table.write_text(
        "# record_id\tfeature_type\tqualifier\tvalue\taction\n"
        "rec1\tCDS\tgene\t^geneA$\thide\n"
        "*\tmisc_feature\tnote\ttransposase\tshow\n",
        encoding="utf-8",
    )

    df = read_feature_visibility_file(str(table))
    assert df is not None
    assert list(df.columns) == ["record_id", "feature_type", "qualifier", "value", "action"]
    assert len(df) == 2
    assert df.iloc[0]["record_id"] == "rec1"
    assert df.iloc[1]["action"] == "show"


def test_read_feature_visibility_file_missing_columns_raises_validation_error(tmp_path: Path) -> None:
    table = tmp_path / "feature_visibility_missing.tsv"
    table.write_text("rec1\tCDS\tgene\t^geneA$\n", encoding="utf-8")
    with pytest.raises(ValidationError):
        read_feature_visibility_file(str(table))


def test_read_feature_visibility_file_extra_columns_raises_parse_error(tmp_path: Path) -> None:
    table = tmp_path / "feature_visibility_extra.tsv"
    table.write_text("rec1\tCDS\tgene\t^geneA$\thide\textra\n", encoding="utf-8")
    with pytest.raises(ParseError):
        read_feature_visibility_file(str(table))


def test_read_feature_visibility_file_not_found_raises_input_file_error() -> None:
    with pytest.raises(InputFileError):
        read_feature_visibility_file("tests/test_inputs/does_not_exist.feature_table.tsv")


def test_compile_feature_visibility_rules_skips_header_and_normalizes_actions() -> None:
    df = _visibility_df(
        [
            ["record_id", "feature_type", "qualifier", "value", "action"],
            ["*", "*", "gene", "^geneA$", "ON"],
            ["*", "*", "note", "transposase", "suppress"],
        ]
    )
    rules = compile_feature_visibility_rules(df)
    assert rules is not None
    assert len(rules) == 2
    assert rules[0]["action"] == "show"
    assert rules[1]["action"] == "hide"


def test_compile_feature_visibility_rules_invalid_regex_raises_parse_error() -> None:
    df = _visibility_df([["*", "*", "gene", "([", "show"]])
    with pytest.raises(ParseError):
        compile_feature_visibility_rules(df)


def test_compile_feature_visibility_rules_invalid_action_raises_validation_error() -> None:
    df = _visibility_df([["*", "*", "gene", "^geneA$", "maybe"]])
    with pytest.raises(ValidationError):
        compile_feature_visibility_rules(df)


@pytest.mark.parametrize(
    ("qualifier", "value_pattern"),
    [
        ("hash", "^f[0-9a-f]{8}$"),
        ("location", "^10\\.\\.90$"),
        ("record_location", "^rec1:10\\.\\.90:\\+$"),
        ("gene", "^geneA$"),
    ],
)
def test_should_render_feature_matches_supported_qualifiers(
    qualifier: str, value_pattern: str
) -> None:
    record = _make_record()
    feature = record.features[0]
    if qualifier == "hash":
        feature_hash = compute_feature_hash(feature, record_id=record.id)
        value_pattern = f"^{feature_hash}$"
    rules = compile_feature_visibility_rules(
        _visibility_df([["*", "*", qualifier, value_pattern, "hide"]])
    )

    assert (
        should_render_feature(
            feature,
            selected_features_set=["CDS"],
            feature_visibility_rules=rules,
            record_id=record.id,
        )
        is False
    )


def test_should_render_feature_first_match_wins() -> None:
    record = _make_record()
    feature = record.features[0]
    rules_hide_first = compile_feature_visibility_rules(
        _visibility_df(
            [
                ["*", "CDS", "gene", "^geneA$", "hide"],
                ["*", "CDS", "gene", "^geneA$", "show"],
            ]
        )
    )
    rules_show_first = compile_feature_visibility_rules(
        _visibility_df(
            [
                ["*", "CDS", "gene", "^geneA$", "show"],
                ["*", "CDS", "gene", "^geneA$", "hide"],
            ]
        )
    )

    assert (
        should_render_feature(
            feature,
            selected_features_set=["CDS"],
            feature_visibility_rules=rules_hide_first,
            record_id=record.id,
        )
        is False
    )
    assert (
        should_render_feature(
            feature,
            selected_features_set=["CDS"],
            feature_visibility_rules=rules_show_first,
            record_id=record.id,
        )
        is True
    )


def test_create_feature_dict_applies_show_and_hide_overrides() -> None:
    record = _make_record()
    default_colors_df = load_default_colors("", "default")
    color_table, default_colors = preprocess_color_tables(None, default_colors_df)
    label_filtering = preprocess_label_filtering(_base_label_filtering())
    rules = compile_feature_visibility_rules(
        _visibility_df(
            [
                ["*", "CDS", "gene", "^geneA$", "hide"],
                ["*", "misc_feature", "note", "transposase", "show"],
            ]
        )
    )

    feature_dict, _ = create_feature_dict(
        record,
        color_table,
        ["CDS"],
        default_colors,
        separate_strands=False,
        resolve_overlaps=False,
        label_filtering=label_filtering,
        feature_visibility_rules=rules,
    )

    rendered_types = {obj.feature_type for obj in feature_dict.values()}
    assert rendered_types == {"misc_feature"}


def test_legend_presence_and_color_usage_follow_visibility_rules() -> None:
    record = _make_record()
    default_colors_df = load_default_colors("", "default")
    color_map, default_color_map = preprocess_color_tables(None, default_colors_df)
    rules = compile_feature_visibility_rules(
        _visibility_df(
            [
                ["*", "CDS", "gene", "^geneA$", "hide"],
                ["*", "misc_feature", "note", "transposase", "show"],
            ]
        )
    )

    features_present = check_feature_presence(
        record,
        ["CDS"],
        feature_visibility_rules=rules,
    )
    used_rules, default_used_features = precompute_used_color_rules(
        record,
        color_map,
        default_color_map,
        {"CDS"},
        feature_visibility_rules=rules,
    )

    assert "CDS" not in features_present
    assert "misc_feature" in features_present
    assert "CDS" not in default_used_features
    assert "misc_feature" in default_used_features
    assert used_rules == set()


def test_circular_cli_feature_table_is_forwarded(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    record = _make_record()
    feature_table_df = _visibility_df([["*", "*", "gene", "^geneA$", "hide"]])
    captured: dict[str, Any] = {}

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda *_args, **_kwargs: [record])
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda _path, _palette: None)
    monkeypatch.setattr(circular_cli_module, "save_figure", lambda _canvas, _formats: None)
    monkeypatch.setattr(
        circular_cli_module,
        "read_feature_visibility_file",
        lambda _path: feature_table_df,
    )

    def fake_assemble(*_args, **kwargs):
        captured["feature_table"] = kwargs.get("feature_table")
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_record", fake_assemble)

    circular_cli_module.circular_main(
        [
            "--gbk",
            "dummy.gb",
            "--feature_table",
            "feature_table.tsv",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["feature_table"] is feature_table_df


def test_linear_cli_feature_table_is_forwarded(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    record = _make_record()
    feature_table_df = _visibility_df([["*", "*", "gene", "^geneA$", "hide"]])
    captured: dict[str, Any] = {}

    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *_args, **_kwargs: [record])
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "save_figure", lambda _canvas, _formats: None)
    monkeypatch.setattr(
        linear_cli_module,
        "read_feature_visibility_file",
        lambda _path: feature_table_df,
    )

    def fake_assemble(*_args, **kwargs):
        captured["feature_table"] = kwargs.get("feature_table")
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(linear_cli_module, "assemble_linear_diagram_from_records", fake_assemble)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "dummy.gb",
            "--feature_table",
            "feature_table.tsv",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["feature_table"] is feature_table_df


def test_circular_gff_loader_uses_keep_all_features_when_feature_table_given(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    record = _make_record()
    captured: dict[str, Any] = {}

    def fake_load_gff_fasta(*_args, **kwargs):
        captured["keep_all_features"] = kwargs.get("keep_all_features")
        return [record]

    monkeypatch.setattr(circular_cli_module, "load_gff_fasta", fake_load_gff_fasta)
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda _path, _palette: None)
    monkeypatch.setattr(circular_cli_module, "save_figure", lambda _canvas, _formats: None)
    monkeypatch.setattr(
        circular_cli_module,
        "read_feature_visibility_file",
        lambda _path: _visibility_df([["*", "*", "gene", "^geneA$", "hide"]]),
    )
    monkeypatch.setattr(
        circular_cli_module,
        "assemble_circular_diagram_from_record",
        lambda *_args, **_kwargs: Drawing(filename=str(tmp_path / "dummy.svg")),
    )

    circular_cli_module.circular_main(
        [
            "--gff",
            "dummy.gff",
            "--fasta",
            "dummy.fasta",
            "--feature_table",
            "feature_table.tsv",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["keep_all_features"] is True


def test_linear_gff_loader_uses_keep_all_features_when_feature_table_given(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    record = _make_record()
    captured: dict[str, Any] = {}

    def fake_load_gff_fasta(*_args, **kwargs):
        captured["keep_all_features"] = kwargs.get("keep_all_features")
        return [record]

    monkeypatch.setattr(linear_cli_module, "load_gff_fasta", fake_load_gff_fasta)
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "save_figure", lambda _canvas, _formats: None)
    monkeypatch.setattr(
        linear_cli_module,
        "read_feature_visibility_file",
        lambda _path: _visibility_df([["*", "*", "gene", "^geneA$", "hide"]]),
    )
    monkeypatch.setattr(
        linear_cli_module,
        "assemble_linear_diagram_from_records",
        lambda *_args, **_kwargs: Drawing(filename=str(tmp_path / "dummy.svg")),
    )

    linear_cli_module.linear_main(
        [
            "--gff",
            "dummy.gff",
            "--fasta",
            "dummy.fasta",
            "--feature_table",
            "feature_table.tsv",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["keep_all_features"] is True
