from __future__ import annotations

import json
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw.analysis.protein_colinearity import OrthogroupMember, OrthogroupResult
from gbdraw.api.options import ColorOptions, LinearDiagramOptions, LinearOutputOptions
from gbdraw.api.request_render import build_request_diagram
from gbdraw.api.requests import InMemoryRecordSource, LinearDiagramRequest, RecordInput
from gbdraw.features.colors import get_color_with_info, preprocess_color_tables
from gbdraw.features.semantic_selectors import (
    FEATURE_SEMANTIC_SCOPE_QUALIFIER,
    build_feature_semantic_selector_context,
    encode_feature_semantic_selector,
    parse_feature_semantic_selector,
)
from gbdraw.io.colors import load_default_colors
from gbdraw.labels.filtering import preprocess_label_filtering


_VECTORS = json.loads(
    (Path(__file__).parent / "fixtures" / "feature_semantic_selector_vectors.json")
    .read_text(encoding="utf-8")
)


def _member(group_id: str) -> OrthogroupMember:
    return OrthogroupMember(
        orthogroup_id=group_id,
        protein_id="protein-1",
        record_index=0,
        feature_index=0,
        record_id="record",
        label="Shared",
        start=5,
        end=35,
        strand=1,
        feature_svg_id=None,
        source_protein_id=None,
    )


def _semantic_color_map(
    kind: str,
    value: str,
    *,
    feature_type: str = "*",
) -> tuple[dict, dict]:
    table = DataFrame(
        [[
            feature_type,
            FEATURE_SEMANTIC_SCOPE_QUALIFIER,
            encode_feature_semantic_selector(kind, value),
            "#abcdef",
            kind,
        ]],
        columns=["feature_type", "qualifier_key", "value", "color", "caption"],
    )
    table.attrs["gbdraw_specific_rule_schema"] = 6
    return preprocess_color_tables(table, load_default_colors("", "default"))


def test_semantic_selector_literals_match_shared_python_web_vectors() -> None:
    for vector in _VECTORS:
        assert encode_feature_semantic_selector(
            vector["kind"], vector["value"]
        ) == vector["literal"]
        assert parse_feature_semantic_selector(vector["literal"]) == (
            vector["kind"],
            vector["value"],
        )


def test_rendered_source_and_similarity_semantics_share_one_context() -> None:
    feature = SeqFeature(
        SimpleLocation(5, 35, strand=1),
        type="CDS",
        qualifiers={"product": ["Shared annotation"], "translation": ["M" * 10]},
    )
    record = SeqRecord(
        Seq("ATG" * 30),
        id="record",
        annotations={"molecule_type": "DNA", "gbdraw_record_key": "record-key"},
    )
    record.features = [feature]
    groups = OrthogroupResult(
        orthogroups={"og-alpha": [_member("og-alpha")], "og-beta": [_member("og-beta")]},
        member_by_protein_id={},
    )
    context = build_feature_semantic_selector_context(
        [record],
        label_filtering=preprocess_label_filtering({}),
        orthogroups=groups,
    )
    assert context.values_for_feature(feature).similarity_group_ids == (
        "og-alpha",
        "og-beta",
    )

    for kind, value in (
        ("source-annotation-label", "Shared annotation"),
        ("rendered-label", "Shared annotation"),
        ("similarity-group", "og-alpha"),
        ("similarity-group", "og-beta"),
    ):
        rules, defaults = _semantic_color_map(kind, value)
        assert get_color_with_info(
            feature,
            rules,
            defaults,
            record_id=record.id,
            feature_semantic_selector_context=context,
        ) == ("#abcdef", kind)


def test_linear_render_keeps_wildcard_semantic_color_in_legend() -> None:
    feature = SeqFeature(
        SimpleLocation(5, 35, strand=1),
        type="CDS",
        qualifiers={"product": ["Shared annotation"], "translation": ["M" * 10]},
    )
    record = SeqRecord(
        Seq("ATG" * 30),
        id="record",
        annotations={"molecule_type": "DNA"},
    )
    record.features = [feature]
    color_table = DataFrame(
        [[
            "*",
            FEATURE_SEMANTIC_SCOPE_QUALIFIER,
            encode_feature_semantic_selector(
                "source-annotation-label",
                "Shared annotation",
            ),
            "#abcdef",
            "Shared annotation",
        ]],
        columns=["feature_type", "qualifier_key", "value", "color", "caption"],
    )
    color_table.attrs["gbdraw_specific_rule_schema"] = 6

    drawing = build_request_diagram(
        LinearDiagramRequest(
            records=(
                RecordInput(
                    source=InMemoryRecordSource(record),
                    record_key="record-key",
                ),
            ),
            options=LinearDiagramOptions(
                selected_features_set=("CDS",),
                colors=ColorOptions(color_table=color_table),
                output=LinearOutputOptions(legend="right"),
            ),
        )
    ).drawing
    svg = drawing.tostring()

    assert 'data-legend-key="Shared annotation"' in svg
    assert 'fill="#abcdef"' in svg
    assert 'data-legend-key="CDS"' not in svg
