from __future__ import annotations

import re

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.features.selector_values import build_feature_selector_values
from gbdraw.features.visibility import compile_feature_visibility_rules, should_render_feature


def _visibility_df(rows: list[list[str]]) -> pd.DataFrame:
    return pd.DataFrame(
        rows,
        columns=["record_id", "feature_type", "qualifier", "value", "action"],
    )


def test_selector_values_match_visibility_supported_selectors() -> None:
    record = SeqRecord(Seq("A" * 200), id="rec1")
    feature = SeqFeature(
        FeatureLocation(10, 90, strand=1),
        type="CDS",
        qualifiers={
            "protein_id": ["ABC123"],
            "locus_tag": ["gene_001"],
            "product": ["example protein"],
        },
    )
    record.features = [feature]

    selector = build_feature_selector_values(feature, record_id=record.id)

    assert selector["location"] == "10..90"
    assert selector["record_location"] == "rec1:10..90:+"
    assert selector["qualifiers"]["protein_id"] == ["ABC123"]
    assert selector["qualifiers"]["locus_tag"] == ["gene_001"]

    for qualifier, value in [
        ("hash", selector["hash"]),
        ("location", selector["location"]),
        ("record_location", selector["record_location"]),
        ("protein_id", "ABC123"),
    ]:
        rules = compile_feature_visibility_rules(
            _visibility_df([["rec1", "CDS", qualifier, f"^{re.escape(str(value))}$", "off"]])
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


def test_selector_values_match_compound_location_semantics() -> None:
    feature = SeqFeature(
        CompoundLocation(
            [
                FeatureLocation(0, 3, strand=1),
                FeatureLocation(6, 9, strand=1),
            ]
        ),
        type="misc_feature",
        qualifiers={"note": ["joined"]},
    )

    selector = build_feature_selector_values(feature, record_id="rec1")

    assert selector["location"] == "0..9"
    assert selector["record_location"] == "rec1:0..9:+"


def test_selector_values_normalize_strand_tokens() -> None:
    assert (
        build_feature_selector_values(
            SeqFeature(FeatureLocation(0, 3, strand=1), type="misc_feature"),
            record_id="rec1",
        )["record_location"]
        == "rec1:0..3:+"
    )
    assert (
        build_feature_selector_values(
            SeqFeature(FeatureLocation(0, 3, strand=-1), type="misc_feature"),
            record_id="rec1",
        )["record_location"]
        == "rec1:0..3:-"
    )
    assert (
        build_feature_selector_values(
            SeqFeature(FeatureLocation(0, 3, strand=None), type="misc_feature"),
            record_id="rec1",
        )["record_location"]
        == "rec1:0..3:undefined"
    )


def test_selector_values_do_not_invent_display_id_or_name_qualifiers() -> None:
    feature = SeqFeature(
        FeatureLocation(0, 9, strand=1),
        type="CDS",
        qualifiers={"locus_tag": ["LT001"]},
    )
    feature.id = "display-id"
    feature.name = "display-name"  # type: ignore[attr-defined]

    selector = build_feature_selector_values(feature, record_id="rec1")

    assert selector["qualifiers"] == {"locus_tag": ["LT001"]}
