from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw
import gbdraw.circular as circular_cli_module
import gbdraw.linear as linear_cli_module
import gbdraw.api.request_render as request_render_module
from gbdraw.config.models import CircularRenderProfile, GbdrawConfig
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
    has_forced_label_overrides,
    preprocess_label_filtering,
    read_label_override_file,
)
from tests.utils.feature_fixtures import (
    make_origin_spanning_feature_object as _make_origin_spanning_feature_object,
    make_origin_spanning_seq_feature as _make_origin_spanning_seq_feature,
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


def _make_seq_feature(
    *,
    product: str = "enzyme alpha",
    gene: str = "geneA",
    locus_tag: str = "LT0001",
    protein_id: str = "WP_000001",
) -> SeqFeature:
    return SeqFeature(
        FeatureLocation(10, 90, strand=1),
        type="CDS",
        qualifiers={
            "product": [product],
            "gene": [gene],
            "locus_tag": [locus_tag],
            "protein_id": [protein_id],
        },
    )


def _make_record() -> SeqRecord:
    record = SeqRecord(Seq("A" * 400), id="rec1")
    record.features = [_make_seq_feature()]
    return record


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


def test_get_label_text_blacklist_matching_is_case_insensitive() -> None:
    feature = _make_seq_feature(product="dUTPase")
    filtering = preprocess_label_filtering(
        _base_filtering(
            blacklist_keywords=["dUTPase"],
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


def test_get_label_text_whitelist_regex_matches_product() -> None:
    feature = _make_seq_feature(product="wsv360-like protein")
    filtering = preprocess_label_filtering(
        _base_filtering(
            whitelist_df=_whitelist_df([["CDS", "product", "wsv.+-like protein"]]),
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "wsv360-like protein"


def test_get_label_text_whitelist_regex_is_case_insensitive() -> None:
    feature = _make_seq_feature(product="WSV514-like protein")
    filtering = preprocess_label_filtering(
        _base_filtering(
            whitelist_df=_whitelist_df([["CDS", "product", "wsv.+-like protein"]]),
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "WSV514-like protein"


def test_get_label_text_whitelist_anchored_regex_matches_exact_label() -> None:
    feature = _make_seq_feature()
    filtering = preprocess_label_filtering(
        _base_filtering(
            whitelist_df=_whitelist_df([["CDS", "product", "^enzyme alpha$"]]),
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "enzyme alpha"


def test_get_label_text_whitelist_literal_pattern_still_matches() -> None:
    feature = _make_seq_feature()
    filtering = preprocess_label_filtering(
        _base_filtering(
            whitelist_df=_whitelist_df([["CDS", "product", "enzyme alpha"]]),
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "enzyme alpha"


def test_get_label_text_whitelist_hash_rule_matches_feature_hash_regex() -> None:
    feature = _make_seq_feature()
    filtering = preprocess_label_filtering(
        _base_filtering(
            whitelist_df=_whitelist_df([["CDS", "hash", "^f[0-9a-f]{8}$"]]),
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "enzyme alpha"


def test_get_label_text_whitelist_location_rule_matches_feature_regex() -> None:
    feature = _make_seq_feature()
    filtering = preprocess_label_filtering(
        _base_filtering(
            whitelist_df=_whitelist_df([["CDS", "location", "^10\\.\\.9\\d$"]]),
        )
    )

    assert get_label_text(feature, filtering, record_id="rec1") == "enzyme alpha"


def test_get_label_text_whitelist_record_location_rule_matches_feature_regex() -> None:
    feature = _make_seq_feature()
    filtering = preprocess_label_filtering(
        _base_filtering(
            whitelist_df=_whitelist_df([["CDS", "record_location", "^rec1:10\\.\\.9\\d:\\+$"]]),
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
    # The factory stores each FeatureObject's source hash (OV-401).
    feature_object.feature_hash = feature_hash
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


def test_get_label_text_feature_object_origin_spanning_hash_override_uses_feature_hash() -> None:
    seq_feature = _make_origin_spanning_seq_feature()
    feature_hash = compute_feature_hash(seq_feature, record_id="rec1")
    feature_object = _make_origin_spanning_feature_object(record_id="rec1")
    feature_object.feature_hash = feature_hash
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
    feature_object.feature_hash = d_loop_hash
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
        {
            "labels.circular.scope": "outer",
            "canvas.strandedness": True,
            "canvas.circular.track_type": 'tuckin',
            "canvas.resolve_overlaps": False,
        },
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
        CircularRenderProfile(cfg),
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


def test_get_label_text_invalid_whitelist_regex_raises_parse_error() -> None:
    with pytest.raises(ParseError):
        preprocess_label_filtering(
            _base_filtering(
                whitelist_df=_whitelist_df(
                    [
                        ["CDS", "product", "["],
                    ]
                )
            )
        )


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
    stub_typed_request_export,
) -> None:
    record = _make_record()
    override_df = _rules_df([["*", "*", "label", "^enzyme alpha$", "CLI label"]])
    captured: dict[str, Any] = {}

    monkeypatch.setattr(request_render_module, "load_gbks", lambda *_args, **_kwargs: [record])
    monkeypatch.setattr(request_render_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(request_render_module, "read_label_override_file", lambda _path: override_df)

    def fake_assemble(*_args: Any, **kwargs: Any) -> Drawing:
        captured["label_override_df"] = kwargs[
            "options"
        ].label_override_table
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(request_render_module, "build_circular_diagram", fake_assemble)

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

    assert captured["label_override_df"] is override_df


def test_linear_cli_label_table_injects_override_df(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    record = _make_record()
    override_df = _rules_df([["*", "*", "label", "^enzyme alpha$", "CLI label"]])
    captured: dict[str, Any] = {}

    monkeypatch.setattr(request_render_module, "load_gbks", lambda *_args, **_kwargs: [record])
    monkeypatch.setattr(request_render_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(request_render_module, "read_label_override_file", lambda _path: override_df)

    def fake_render_request(request, **_kwargs):
        resolved = request_render_module.resolve_request(request)
        captured["canonical_request"] = resolved
        return SimpleNamespace(
            drawing=Drawing(filename=str(tmp_path / "dummy.svg")),
            interactive_context=None,
            records=tuple(item.source.record for item in resolved.records),
            losat_cache_entries=(),
            losat_derived_cache_entries=(),
            protein_identity_manifest=None,
            request=resolved,
            annotation_warnings=(),
            feature_identity_notices=(),
        )

    monkeypatch.setattr(linear_cli_module, "render_request", fake_render_request)

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

    label_override_table = captured["canonical_request"].options.label_override_table
    assert label_override_table is not None
    pd.testing.assert_frame_equal(label_override_table, override_df)


# A per-feature override (a `hash` row) decides one feature's label whatever
# the label display scope selects: non-empty text shows it, empty text hides it.
_BGC_INPUTS = Path(__file__).parent / "test_inputs"


def _rendered_feature_labels(diagram: Any) -> list[tuple[str, str]]:
    svg = diagram.to_svg()
    root = ET.fromstring(svg if isinstance(svg, str) else svg.decode("utf-8"))
    return [
        (str(node.get("data-label-feature-id")), "".join(node.itertext()))
        for node in root.iter("{http://www.w3.org/2000/svg}text")
        if node.get("data-label-feature-id")
    ]


def _bgc_records_and_neor_hash() -> tuple[list[SeqRecord], str]:
    records = list(gbdraw.read_genbank([
        str(_BGC_INPUTS / "BGC0000708.gbk"),
        str(_BGC_INPUTS / "BGC0000709.gbk"),
    ]))
    neor = next(
        feature for feature in records[1].features
        if feature.type == "CDS" and feature.qualifiers.get("product") == ["putative regulator, NeoR"]
    )
    return records, compute_feature_hash(neor, record_id=records[1].id)


def _draw_linear_labels(records: list[SeqRecord], scope: str, rows: list[list[str]] | None):
    options = gbdraw.LinearOptions(
        labels=gbdraw.LabelOptions(overrides=_rules_df(rows) if rows else None),
        config_overrides={"labels.linear.scope": scope},
    )
    return _rendered_feature_labels(gbdraw.draw_linear(records, options=options))


def test_has_forced_label_overrides_counts_only_hash_rows_with_text() -> None:
    def forced(rows: list[list[str]]) -> bool:
        return has_forced_label_overrides(_base_filtering(label_override_df=_rules_df(rows)))

    assert forced([["*", "*", "hash", "^f1$", "Shown"]]) is True
    assert forced([["*", "*", "hash", "^f1$", ""]]) is False
    assert forced([["*", "CDS", "gene", "^geneA$", "Renamed"]]) is False
    assert has_forced_label_overrides(_base_filtering()) is False


def test_get_label_text_overrides_only_keeps_per_feature_decisions() -> None:
    feature = _make_seq_feature()
    feature_hash = compute_feature_hash(feature, record_id="rec1")
    filtering = preprocess_label_filtering(
        _base_filtering(
            label_override_df=_rules_df([
                ["*", "CDS", "gene", "^geneA$", "Ordinary override"],
            ])
        )
    )
    assert get_label_text(feature, filtering, record_id="rec1") == "Ordinary override"
    assert get_label_text(feature, filtering, record_id="rec1", overrides_only=True) == ""

    forced = preprocess_label_filtering(
        _base_filtering(
            label_override_df=_rules_df([
                ["*", "*", "hash", f"^{re.escape(feature_hash)}$", "Forced"],
            ])
        )
    )
    assert get_label_text(feature, forced, record_id="rec1", overrides_only=True) == "Forced"


@pytest.mark.linear
@pytest.mark.parametrize("scope", ["first", "none"])
def test_linear_per_feature_override_shows_label_outside_scope(scope: str) -> None:
    records, neor_hash = _bgc_records_and_neor_hash()
    base = _draw_linear_labels(records, scope, None)
    forced = _draw_linear_labels(
        records, scope, [[records[1].id, "CDS", "hash", f"^{re.escape(neor_hash)}$", "NeoR"]]
    )

    assert not any(feature_id.endswith("_record_2") for feature_id, _ in base)
    assert [label for label in forced if label not in base] == [(f"{neor_hash}_record_2", "NeoR")]
    assert len(forced) == len(base) + 1


@pytest.mark.linear
def test_linear_override_label_without_scoped_records_uses_all_record_label_size() -> None:
    records, neor_hash = _bgc_records_and_neor_hash()
    rows = [[records[1].id, "CDS", "hash", f"^{re.escape(neor_hash)}$", "NeoR"]]

    def neor_font_size(scope: str) -> str:
        options = gbdraw.LinearOptions(
            labels=gbdraw.LabelOptions(overrides=_rules_df(rows)),
            config_overrides={"labels.linear.scope": scope},
        )
        svg = gbdraw.draw_linear(records, options=options).to_svg()
        root = ET.fromstring(svg if isinstance(svg, str) else svg.decode("utf-8"))
        node = next(
            node for node in root.iter("{http://www.w3.org/2000/svg}text")
            if "".join(node.itertext()) == "NeoR"
        )
        return str(node.get("font-size"))

    assert neor_font_size("none") == neor_font_size("all")


@pytest.mark.linear
def test_linear_ordinary_override_does_not_extend_label_scope() -> None:
    records, _neor_hash = _bgc_records_and_neor_hash()
    base = _draw_linear_labels(records, "first", None)
    renamed = _draw_linear_labels(records, "first", [[records[1].id, "CDS", "gene", "^neoR$", "NeoR"]])

    assert renamed == base


@pytest.mark.linear
def test_linear_per_feature_override_hides_label_in_scope() -> None:
    records, _neor_hash = _bgc_records_and_neor_hash()
    base = _draw_linear_labels(records, "all", None)
    hidden_id, _text = base[0]
    feature_hash = hidden_id.rsplit("_record_", 1)[0]
    record_id = records[int(hidden_id.rsplit("_record_", 1)[1]) - 1].id
    hidden = _draw_linear_labels(
        records, "all", [[record_id, "*", "hash", f"^{re.escape(feature_hash)}$", ""]]
    )

    assert hidden == [label for label in base if label[0] != hidden_id]


@pytest.mark.circular
def test_circular_per_feature_override_shows_only_its_label_when_scope_is_none() -> None:
    records, neor_hash = _bgc_records_and_neor_hash()
    record = records[1]

    def labels(rows: list[list[str]] | None) -> list[tuple[str, str]]:
        options = gbdraw.CircularOptions(
            labels=gbdraw.LabelOptions(overrides=_rules_df(rows) if rows else None),
            config_overrides={"labels.circular.scope": "none"},
        )
        return _rendered_feature_labels(gbdraw.draw_circular([record], options=options))

    assert labels(None) == []
    assert labels([[record.id, "CDS", "hash", f"^{re.escape(neor_hash)}$", "NeoR"]]) == [(neor_hash, "NeoR")]
    assert labels([[record.id, "CDS", "hash", f"^{re.escape(neor_hash)}$", ""]]) == []
