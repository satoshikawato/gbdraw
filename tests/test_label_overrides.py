from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.circular as circular_cli_module
import gbdraw.linear as linear_cli_module
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import InputFileError, ParseError, ValidationError
from gbdraw.features.colors import compute_feature_hash, preprocess_color_tables
from gbdraw.features.factory import create_feature_dict
from gbdraw.features.visibility import compile_feature_visibility_rules
from gbdraw.features.objects import FeatureLocationPart, FeatureObject
from gbdraw.io.colors import load_default_colors
from gbdraw.labels.circular import prepare_label_list
from gbdraw.labels.filtering import (
    get_label_text,
    preprocess_label_filtering,
    read_label_override_file,
)


def _rules_df(rows: list[list[str]]) -> pd.DataFrame:
    return pd.DataFrame(
        rows,
        columns=["record_id", "feature_type", "qualifier", "value", "label_text"],
    )


def _whitelist_df(rows: list[list[str]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["feature_type", "qualifier", "keyword"])


def _base_filtering(
    *,
    blacklist_keywords: list[str] | None = None,
    whitelist_df: pd.DataFrame | None = None,
    qualifier_priority_df: pd.DataFrame | None = None,
    label_override_df: pd.DataFrame | None = None,
) -> dict:
    return {
        "blacklist_keywords": blacklist_keywords or [],
        "whitelist_df": whitelist_df,
        "qualifier_priority_df": qualifier_priority_df,
        "label_override_df": label_override_df,
    }


def _make_seq_feature() -> SeqFeature:
    return SeqFeature(
        FeatureLocation(10, 90, strand=1),
        type="CDS",
        qualifiers={
            "product": ["enzyme alpha"],
            "gene": ["geneA"],
            "locus_tag": ["LT0001"],
            "protein_id": ["WP_000001"],
        },
    )


def _make_record() -> SeqRecord:
    record = SeqRecord(Seq("A" * 400), id="rec1")
    record.features = [_make_seq_feature()]
    return record


def _make_origin_spanning_seq_feature() -> SeqFeature:
    return SeqFeature(
        CompoundLocation(
            [
                SimpleLocation(0, 576, strand=-1),
                SimpleLocation(16023, 16569, strand=-1),
            ],
            operator="join",
        ),
        type="D-loop",
        qualifiers={},
    )


def _make_origin_spanning_feature_object(record_id: str = "rec1") -> FeatureObject:
    return FeatureObject(
        feature_id="feature_000000099",
        location=[
            FeatureLocationPart("block", "001", "negative", 1, 576, False),
            FeatureLocationPart("block", "002", "negative", 16023, 16569, True),
        ],
        is_directional=False,
        color="#cccccc",
        note="",
        label_text="",
        coordinates=[
            SimpleLocation(0, 576, strand=-1),
            SimpleLocation(16023, 16569, strand=-1),
        ],
        type="D-loop",
        qualifiers={},
        record_id=record_id,
    )


def test_read_label_override_file_ok(tmp_path: Path) -> None:
    table = tmp_path / "label_override.tsv"
    table.write_text(
        "# record_id\tfeature_type\tqualifier\tvalue\tlabel_text\n"
        "rec1\tCDS\tgene\t^geneA$\tGene A\n"
        "*\t*\tlabel\t^enzyme alpha$\tAnnotated enzyme\n",
        encoding="utf-8",
    )

    df = read_label_override_file(str(table))
    assert df is not None
    assert list(df.columns) == ["record_id", "feature_type", "qualifier", "value", "label_text"]
    assert len(df) == 2
    assert df.iloc[0]["record_id"] == "rec1"
    assert df.iloc[1]["qualifier"] == "label"


def test_read_label_override_file_missing_columns_raises_validation_error(tmp_path: Path) -> None:
    table = tmp_path / "label_override_missing.tsv"
    table.write_text("CDS\tlabel\t^geneA$\tGene A\n", encoding="utf-8")
    with pytest.raises(ValidationError):
        read_label_override_file(str(table))


def test_read_label_override_file_extra_columns_raises_parse_error(tmp_path: Path) -> None:
    table = tmp_path / "label_override_malformed.tsv"
    table.write_text("*\tCDS\tlabel\t^geneA$\tGene A\textra\n", encoding="utf-8")
    with pytest.raises(ParseError):
        read_label_override_file(str(table))


def test_read_label_override_file_not_found_raises_input_file_error() -> None:
    with pytest.raises(InputFileError):
        read_label_override_file("tests/test_inputs/does_not_exist.label_override.tsv")


def test_read_label_override_file_empty_label_text_allowed(tmp_path: Path) -> None:
    table = tmp_path / "label_override_empty_text.tsv"
    table.write_text("rec1\tCDS\thash\t^f123$\t\n", encoding="utf-8")

    df = read_label_override_file(str(table))
    assert df is not None
    assert str(df.iloc[0]["label_text"]) == ""


def test_get_label_text_row_order_first_wins() -> None:
    feature = _make_seq_feature()
    filtering = preprocess_label_filtering(
        _base_filtering(
            label_override_df=_rules_df(
                [
                    ["*", "*", "label", "^enzyme", "First match"],
                    ["*", "*", "label", "^enzyme alpha$", "Second match"],
                ]
            )
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "First match"


def test_get_label_text_record_constraint_works() -> None:
    feature = _make_seq_feature()
    filtering = preprocess_label_filtering(
        _base_filtering(
            label_override_df=_rules_df(
                [
                    ["rec2", "CDS", "label", "^enzyme alpha$", "Wrong record"],
                    ["rec1", "CDS", "label", "^enzyme alpha$", "Exact match"],
                ]
            )
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "Exact match"


def test_get_label_text_hash_override_works() -> None:
    feature = _make_seq_feature()
    feature_hash = compute_feature_hash(feature, record_id="rec1")
    filtering = preprocess_label_filtering(
        _base_filtering(
            label_override_df=_rules_df(
                [
                    ["*", "*", "hash", f"^{re.escape(feature_hash)}$", "From hash key"],
                    ["*", "*", "label", "^enzyme alpha$", "From base label"],
                ]
            )
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "From hash key"


def test_get_label_text_record_location_override_works() -> None:
    feature = _make_seq_feature()
    filtering = preprocess_label_filtering(
        _base_filtering(
            label_override_df=_rules_df(
                [
                    ["*", "*", "record_location", "^rec1:10..90:\\+$", "From record_location key"],
                ]
            )
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "From record_location key"


def test_get_label_text_feature_type_constraint_works() -> None:
    feature = _make_seq_feature()
    filtering = preprocess_label_filtering(
        _base_filtering(
            label_override_df=_rules_df(
                [
                    ["*", "tRNA", "label", "^enzyme alpha$", "Wrong feature type"],
                    ["*", "CDS", "label", "^enzyme alpha$", "Right feature type"],
                ]
            )
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "Right feature type"


def test_get_label_text_does_not_override_hidden_label() -> None:
    feature = _make_seq_feature()
    filtering = preprocess_label_filtering(
        _base_filtering(
            blacklist_keywords=["enzyme"],
            label_override_df=_rules_df(
                [
                    ["*", "*", "label", "^enzyme alpha$", "Should not re-enable"],
                ]
            ),
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == ""


def test_get_label_text_hash_override_bypasses_blacklist() -> None:
    feature = _make_seq_feature()
    feature_hash = compute_feature_hash(feature, record_id="rec1")
    filtering = preprocess_label_filtering(
        _base_filtering(
            blacklist_keywords=["enzyme"],
            label_override_df=_rules_df(
                [
                    ["*", "*", "hash", f"^{re.escape(feature_hash)}$", "Forced visible label"],
                ]
            ),
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "Forced visible label"


def test_get_label_text_hash_override_bypasses_whitelist() -> None:
    feature = _make_seq_feature()
    feature_hash = compute_feature_hash(feature, record_id="rec1")
    filtering = preprocess_label_filtering(
        _base_filtering(
            whitelist_df=_whitelist_df([["CDS", "product", "non-matching-keyword"]]),
            label_override_df=_rules_df(
                [
                    ["*", "*", "hash", f"^{re.escape(feature_hash)}$", "Whitelist bypass label"],
                ]
            ),
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "Whitelist bypass label"


def test_get_label_text_whitelist_hash_rule_matches_feature_hash() -> None:
    feature = _make_seq_feature()
    feature_hash = compute_feature_hash(feature, record_id="rec1")
    filtering = preprocess_label_filtering(
        _base_filtering(
            whitelist_df=_whitelist_df([["CDS", "hash", feature_hash]]),
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "enzyme alpha"


def test_get_label_text_whitelist_record_location_rule_matches_feature() -> None:
    feature = _make_seq_feature()
    filtering = preprocess_label_filtering(
        _base_filtering(
            whitelist_df=_whitelist_df([["CDS", "record_location", "rec1:10..90:+"]]),
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "enzyme alpha"


def test_get_label_text_whitelist_record_location_rule_non_matching_hides_label() -> None:
    feature = _make_seq_feature()
    filtering = preprocess_label_filtering(
        _base_filtering(
            whitelist_df=_whitelist_df([["CDS", "record_location", "rec1:10..91:+"]]),
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == ""


def test_get_label_text_hash_override_empty_text_hides_label() -> None:
    feature = _make_seq_feature()
    feature_hash = compute_feature_hash(feature, record_id="rec1")
    filtering = preprocess_label_filtering(
        _base_filtering(
            whitelist_df=_whitelist_df([["CDS", "product", "enzyme alpha"]]),
            label_override_df=_rules_df(
                [
                    ["*", "*", "hash", f"^{re.escape(feature_hash)}$", ""],
                ]
            ),
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == ""


def test_get_label_text_feature_object_hash_override_works() -> None:
    seq_feature = _make_seq_feature()
    feature_hash = compute_feature_hash(seq_feature, record_id="rec1")
    feature_object = FeatureObject(
        feature_id="feature_000000001",
        location=[FeatureLocationPart("block", "001", "positive", 10, 90, True)],
        is_directional=True,
        color="#cccccc",
        note="",
        label_text="",
        coordinates=[],
        type="CDS",
        qualifiers={"product": ["enzyme alpha"], "gene": ["geneA"]},
        record_id="rec1",
    )
    filtering = preprocess_label_filtering(
        _base_filtering(
            label_override_df=_rules_df(
                [
                    ["*", "*", "hash", f"^{re.escape(feature_hash)}$", "Object hash match"],
                ]
            )
        )
    )

    assert get_label_text(feature_object, filtering) == "Object hash match"


def test_get_label_text_feature_object_record_location_override_works() -> None:
    feature_object = FeatureObject(
        feature_id="feature_000000001",
        location=[FeatureLocationPart("block", "001", "positive", 10, 90, True)],
        is_directional=True,
        color="#cccccc",
        note="",
        label_text="",
        coordinates=[],
        type="CDS",
        qualifiers={"product": ["enzyme alpha"], "gene": ["geneA"]},
        record_id="rec1",
    )
    filtering = preprocess_label_filtering(
        _base_filtering(
            label_override_df=_rules_df(
                [
                    ["*", "*", "record_location", "^rec1:10..90:\\+$", "Object record location match"],
                ]
            )
        )
    )

    assert get_label_text(feature_object, filtering) == "Object record location match"


def test_get_label_text_feature_object_origin_spanning_hash_override_uses_coordinates() -> None:
    seq_feature = _make_origin_spanning_seq_feature()
    feature_hash = compute_feature_hash(seq_feature, record_id="rec1")
    feature_object = _make_origin_spanning_feature_object(record_id="rec1")
    filtering = preprocess_label_filtering(
        _base_filtering(
            label_override_df=_rules_df(
                [
                    ["*", "*", "hash", f"^{re.escape(feature_hash)}$", "D-loop"],
                ]
            )
        )
    )

    assert get_label_text(feature_object, filtering) == "D-loop"


def test_get_label_text_feature_object_origin_spanning_record_location_override_uses_coordinates() -> None:
    feature_object = _make_origin_spanning_feature_object(record_id="rec1")
    filtering = preprocess_label_filtering(
        _base_filtering(
            label_override_df=_rules_df(
                [
                    ["*", "*", "record_location", "^rec1:0..16569:-$", "D-loop"],
                ]
            )
        )
    )

    assert get_label_text(feature_object, filtering) == "D-loop"


def test_get_label_text_hmmtdna_d_loop_hash_override_matches_origin_spanning_feature_object() -> None:
    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    record = SeqIO.read(str(input_path), "genbank")
    d_loop_feature = next(feature for feature in record.features if feature.type == "D-loop")
    d_loop_hash = compute_feature_hash(d_loop_feature, record_id=record.id)
    feature_object = _make_origin_spanning_feature_object(record_id=record.id)
    filtering = preprocess_label_filtering(
        _base_filtering(
            label_override_df=_rules_df(
                [
                    [record.id, "D-loop", "hash", f"^{re.escape(d_loop_hash)}$", "D-loop"],
                ]
            )
        )
    )

    assert get_label_text(feature_object, filtering) == "D-loop"


def test_prepare_label_list_hmmtdna_d_loop_hash_override_survives_feature_object_recheck() -> None:
    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    record = SeqIO.read(str(input_path), "genbank")
    d_loop_feature = next(feature for feature in record.features if feature.type == "D-loop")
    d_loop_hash = compute_feature_hash(d_loop_feature, record_id=record.id)

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict = modify_config_dict(
        config_dict,
        show_labels=True,
        strandedness=True,
        track_type="tuckin",
        resolve_overlaps=False,
        allow_inner_labels=False,
    )
    override_df = _rules_df(
        [[record.id, "D-loop", "hash", f"^{re.escape(d_loop_hash)}$", "D-loop"]]
    )
    config_dict["labels"]["filtering"]["label_override_df"] = override_df
    cfg = GbdrawConfig.from_dict(config_dict)

    label_filtering = preprocess_label_filtering(cfg.labels.filtering.as_dict())
    default_colors = load_default_colors("", "default")
    color_table, default_colors = preprocess_color_tables(None, default_colors)
    visibility_rules = compile_feature_visibility_rules(
        pd.DataFrame(
            [[record.id, "D-loop", "hash", f"^{re.escape(d_loop_hash)}$", "show"]],
            columns=["record_id", "feature_type", "qualifier", "value", "action"],
        )
    )
    selected_features = ["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "repeat_region"]
    feature_dict, _ = create_feature_dict(
        record,
        color_table,
        selected_features,
        default_colors,
        cfg.canvas.strandedness,
        cfg.canvas.resolve_overlaps,
        label_filtering,
        feature_visibility_rules=visibility_rules,
    )

    labels = prepare_label_list(
        feature_dict,
        len(record.seq),
        cfg.canvas.circular.radius,
        cfg.canvas.circular.track_ratio,
        config_dict,
        cfg=cfg,
    )

    d_loop_label = next((label for label in labels if label.get("label_text") == "D-loop"), None)
    assert d_loop_label is not None

    d_loop_parts = list(d_loop_feature.location.parts)
    merged_start = max(int(part.start) for part in d_loop_parts)
    merged_end = min(int(part.end) for part in d_loop_parts)
    wrapped_span = (len(record.seq) - merged_start) + merged_end
    expected_midpoint = (merged_start + (wrapped_span / 2.0)) % len(record.seq)
    assert abs(float(d_loop_label["middle"]) - float(expected_midpoint)) <= 1e-6
    assert bool(d_loop_label["is_embedded"])


def test_get_label_text_invalid_override_regex_raises_parse_error() -> None:
    feature = _make_seq_feature()
    with pytest.raises(ParseError):
        preprocess_label_filtering(
            _base_filtering(
                label_override_df=_rules_df(
                    [
                        ["*", "CDS", "label", "[", "broken"],
                    ]
                )
            )
        )

    baseline = preprocess_label_filtering(_base_filtering())
    assert get_label_text(feature, baseline, record_id="rec1") == "enzyme alpha"


def test_get_label_text_location_override_works() -> None:
    feature = _make_seq_feature()
    filtering = preprocess_label_filtering(
        _base_filtering(
            label_override_df=_rules_df(
                [
                    ["*", "*", "location", "^10..90$", "Location match"],
                ]
            )
        )
    )
    assert get_label_text(feature, filtering, record_id="rec1") == "Location match"


def test_get_label_text_normal_qualifier_override_works() -> None:
    feature = _make_seq_feature()
    filtering = preprocess_label_filtering(
        _base_filtering(
            label_override_df=_rules_df(
                [
                    ["rec1", "CDS", "protein_id", "^WP_000001$", "Protein ID match"],
                ]
            )
        )
    )
    assert get_label_text(feature, filtering, record_id="rec1") == "Protein ID match"


def test_circular_cli_label_table_injects_override_df(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    record = _make_record()
    override_df = _rules_df([["*", "*", "label", "^enzyme alpha$", "CLI label"]])
    captured: dict[str, Any] = {}
    real_modify_config_dict = modify_config_dict

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda *_args, **_kwargs: [record])
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(circular_cli_module, "save_figure", lambda _canvas, _formats: None)
    monkeypatch.setattr(circular_cli_module, "read_label_override_file", lambda _path: override_df)

    def fake_modify(config_dict: dict, **kwargs: Any) -> dict:
        captured["label_table_arg"] = kwargs.get("label_table")
        captured["label_override_df"] = config_dict["labels"]["filtering"].get("label_override_df")
        return real_modify_config_dict(config_dict, **kwargs)

    def fake_assemble(*_args: Any, **_kwargs: Any) -> Drawing:
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(circular_cli_module, "modify_config_dict", fake_modify)
    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_record", fake_assemble)

    circular_cli_module.circular_main(
        [
            "--gbk",
            "dummy.gb",
            "--label_table",
            "table.tsv",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["label_table_arg"] == "table.tsv"
    assert captured["label_override_df"] is override_df


def test_linear_cli_label_table_injects_override_df(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    record = _make_record()
    override_df = _rules_df([["*", "*", "label", "^enzyme alpha$", "CLI label"]])
    captured: dict[str, Any] = {}
    real_modify_config_dict = modify_config_dict

    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *_args, **_kwargs: [record])
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "save_figure", lambda _canvas, _formats: None)
    monkeypatch.setattr(linear_cli_module, "read_label_override_file", lambda _path: override_df)

    def fake_modify(config_dict: dict, **kwargs: Any) -> dict:
        captured["label_table_arg"] = kwargs.get("label_table")
        captured["label_override_df"] = config_dict["labels"]["filtering"].get("label_override_df")
        return real_modify_config_dict(config_dict, **kwargs)

    def fake_assemble(*_args: Any, **_kwargs: Any) -> Drawing:
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(linear_cli_module, "modify_config_dict", fake_modify)
    monkeypatch.setattr(linear_cli_module, "assemble_linear_diagram_from_records", fake_assemble)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "dummy.gb",
            "--label_table",
            "table.tsv",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["label_table_arg"] == "table.tsv"
    assert captured["label_override_df"] is override_df
