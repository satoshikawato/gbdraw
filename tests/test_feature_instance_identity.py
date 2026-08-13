from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import re

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw.api.options import (
    CircularDiagramOptions,
    CircularOutputOptions,
    ColorOptions,
    LinearDiagramOptions,
    LinearOutputOptions,
)
from gbdraw.api.request_render import (
    build_prepared_interactive_context,
    build_request_diagram,
)
from gbdraw.api.requests import (
    CircularDiagramRequest,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
)
from gbdraw.exceptions import ValidationError
from gbdraw.features.colors import (
    find_specific_color_rule,
    get_color_with_info,
    preprocess_color_tables,
)
from gbdraw.features.instance_identity import (
    FEATURE_INSTANCE_HASH_PATTERN,
    FEATURE_INSTANCE_HASH_QUALIFIER,
    build_feature_instance_identity_plan,
    compute_feature_instance_hash,
)
from gbdraw.features.semantic_selectors import (
    FEATURE_SEMANTIC_SCOPE_QUALIFIER,
    parse_feature_semantic_selector,
)
from gbdraw.io.colors import load_default_colors
from gbdraw.web_support.feature_catalog import build_feature_catalog_item


_PRECEDENCE_ORACLE = json.loads(
    (Path(__file__).parent / "fixtures" / "specific_color_precedence_oracle.json").read_text(
        encoding="utf-8"
    )
)


def _feature(label: str, *, start: int = 10, end: int = 40) -> SeqFeature:
    return SeqFeature(
        SimpleLocation(start, end, strand=1),
        type="CDS",
        qualifiers={"label": [label], "translation": ["M" * 10]},
    )


def _record(
    public_id: str,
    record_key: str,
    features: list[SeqFeature],
) -> SeqRecord:
    record = SeqRecord(
        Seq("ATG" * 80),
        id=public_id,
        annotations={
            "molecule_type": "DNA",
            "gbdraw_record_key": record_key,
        },
    )
    record.features = features
    return record


def _color_table(rows: list[tuple[str, str, str, str, str]]) -> DataFrame:
    return DataFrame(
        rows,
        columns=["feature_type", "qualifier_key", "value", "color", "caption"],
    )


def test_feature_instance_hash_has_frozen_cross_runtime_vectors() -> None:
    assert (
        compute_feature_instance_hash("record-1", "f12345678")
        == "fi1_f6lpj7ii45v7lrrthhg4pci24i"
    )
    assert (
        compute_feature_instance_hash("rec-α", "fdeadbeef~7")
        == "fi1_3rqvxium7ztier2dwg2nzdkp4q"
    )
    assert FEATURE_INSTANCE_HASH_PATTERN.fullmatch(
        compute_feature_instance_hash("record-1", "f12345678")
    )


def test_plan_scopes_duplicate_public_ids_by_record_key() -> None:
    first = _record("duplicate", "record-a", [_feature("a")])
    second = _record("duplicate", "record-b", [_feature("b")])

    plan = build_feature_instance_identity_plan([first, second])
    first_identity = plan.identity_for_feature(first.features[0])
    second_identity = plan.identity_for_feature(second.features[0])

    assert first_identity.stable_feature_id == second_identity.stable_feature_id
    assert first_identity.biological_feature_id == second_identity.biological_feature_id
    assert first_identity.instance_hash != second_identity.instance_hash


def test_plan_disambiguates_same_record_stable_id_collisions() -> None:
    record = _record(
        "collision",
        "record-a",
        [_feature("first"), _feature("second")],
    )

    plan = build_feature_instance_identity_plan(record)
    first = plan.identity_for_feature(record.features[0])
    second = plan.identity_for_feature(record.features[1])

    assert first.stable_feature_id == second.stable_feature_id
    assert first.biological_feature_id == f"{first.stable_feature_id}~0"
    assert second.biological_feature_id == f"{second.stable_feature_id}~1"
    assert first.instance_hash != second.instance_hash


def test_plan_includes_nested_features_in_authoritative_source_order() -> None:
    parent = SeqFeature(SimpleLocation(1, 50, strand=1), type="gene")
    nested = _feature("nested", start=5, end=20)
    parent.sub_features = [nested]
    record = _record("nested-record", "record-a", [parent])

    plan = build_feature_instance_identity_plan(record)

    assert [identity.source_feature_index for identity in plan.identities] == [0, 1]
    assert plan.identity_for_feature(nested) is plan.identities[1]


def test_reserved_instance_literal_precedes_legacy_selectors_within_rule_set() -> None:
    record = _record("record", "record-a", [_feature("target")])
    identity = build_feature_instance_identity_plan(record).identity_for_feature(
        record.features[0]
    )
    default_colors = load_default_colors("", "default")
    rules, defaults = preprocess_color_tables(
        _color_table(
            [
                ("CDS", "hash", identity.stable_feature_id, "#0000ff", "hash"),
                (
                    "CDS",
                    FEATURE_INSTANCE_HASH_QUALIFIER,
                    identity.instance_hash,
                    "#ff0000",
                    "instance",
                ),
            ]
        ),
        default_colors,
    )

    assert get_color_with_info(
        record.features[0],
        rules,
        defaults,
        record_id=record.id,
        feature_instance_identity=identity,
    ) == ("#ff0000", "instance")


def test_exact_type_hash_precedes_wildcard_instance_literal() -> None:
    record = _record("record", "record-a", [_feature("target")])
    identity = build_feature_instance_identity_plan(record).identity_for_feature(
        record.features[0]
    )
    rules, defaults = preprocess_color_tables(
        _color_table(
            [
                (
                    "*",
                    FEATURE_INSTANCE_HASH_QUALIFIER,
                    identity.instance_hash,
                    "#ff0000",
                    "instance",
                ),
                ("CDS", "hash", identity.stable_feature_id, "#0000ff", "hash"),
            ]
        ),
        load_default_colors("", "default"),
    )

    assert get_color_with_info(
        record.features[0],
        rules,
        defaults,
        record_id=record.id,
        feature_instance_identity=identity,
    ) == ("#0000ff", "hash")


def test_biological_instance_hash_qualifier_keeps_legacy_regex_semantics() -> None:
    feature = _feature("legacy")
    feature.qualifiers["instance_hash"] = ["Biological-Value-42"]
    record = _record("record", "record-a", [feature])
    rules, defaults = preprocess_color_tables(
        _color_table(
            [
                (
                    "CDS",
                    "instance_hash",
                    "biological-value-[0-9]+",
                    "#123456",
                    "legacy",
                )
            ]
        ),
        load_default_colors("", "default"),
    )

    assert get_color_with_info(
        feature,
        rules,
        defaults,
        record_id=record.id,
    ) == ("#123456", "legacy")


def test_python_specific_color_resolution_matches_shared_web_precedence_oracle() -> None:
    payload = _PRECEDENCE_ORACLE["feature"]
    feature = _feature(
        "oracle",
        start=payload["start"],
        end=payload["end"],
    )
    feature.qualifiers = deepcopy(payload["qualifiers"])
    feature.instance_hash = payload["instanceHash"]

    for case in _PRECEDENCE_ORACLE["cases"]:
        color_map: dict[str, dict[str, list[tuple[object, str, str]]]] = {}
        for feature_type, qualifier, value, color, caption, match in case["rules"]:
            if match == "literal" and qualifier == FEATURE_SEMANTIC_SCOPE_QUALIFIER:
                pattern: object = parse_feature_semantic_selector(value)
            else:
                pattern = value if match == "literal" else re.compile(value, re.IGNORECASE)
            color_map.setdefault(feature_type, {}).setdefault(qualifier, []).append(
                (pattern, color, caption)
            )
        assert find_specific_color_rule(
            feature,
            color_map,
            record_id=payload["recordId"],
        ) == next(
            (rule[3], rule[4])
            for rule in case["rules"]
            if rule[4] == case["expectedCaption"]
        ), case["id"]


def test_schema5_reserved_spelling_keeps_legacy_source_qualifier_semantics() -> None:
    feature = _feature("legacy")
    feature.qualifiers[FEATURE_INSTANCE_HASH_QUALIFIER] = ["Legacy-Value-42"]
    record = _record("record", "record-a", [feature])
    table = _color_table(
        [
            (
                "CDS",
                FEATURE_INSTANCE_HASH_QUALIFIER,
                "legacy-value-[0-9]+",
                "#654321",
                "legacy schema",
            )
        ]
    )
    table.attrs["gbdraw_specific_rule_schema"] = 5
    rules, defaults = preprocess_color_tables(
        table,
        load_default_colors("", "default"),
    )

    assert get_color_with_info(
        feature,
        rules,
        defaults,
        record_id=record.id,
    ) == ("#654321", "legacy schema")


def test_reserved_instance_rule_requires_context_and_valid_literal() -> None:
    record = _record("record", "record-a", [_feature("target")])
    identity = build_feature_instance_identity_plan(record).identity_for_feature(
        record.features[0]
    )
    rules, defaults = preprocess_color_tables(
        _color_table(
            [
                (
                    "CDS",
                    FEATURE_INSTANCE_HASH_QUALIFIER,
                    identity.instance_hash,
                    "#ff0000",
                    "instance",
                )
            ]
        ),
        load_default_colors("", "default"),
    )
    with pytest.raises(ValidationError, match="shared feature instance identity plan"):
        get_color_with_info(record.features[0], rules, defaults, record_id=record.id)

    with pytest.raises(ValidationError, match="26 lower-case base32"):
        preprocess_color_tables(
            _color_table(
                [
                    (
                        "CDS",
                        FEATURE_INSTANCE_HASH_QUALIFIER,
                        identity.instance_hash.upper(),
                        "#ff0000",
                        "instance",
                    )
                ]
            ),
            load_default_colors("", "default"),
        )


def _catalog_for_exact_request(
    request: CircularDiagramRequest | LinearDiagramRequest,
) -> tuple[dict, object]:
    prepared = build_request_diagram(request)
    assert prepared.feature_instance_identity_plan is not None
    context = build_prepared_interactive_context(prepared)
    item = build_feature_catalog_item(
        prepared.drawing.tostring(),
        context,
        result_index=0,
        result_name="result.svg",
    )
    return item, prepared


def _assert_exact_catalog_fill(
    item: dict,
    *,
    target_instance_hash: str,
    target_color: str,
) -> None:
    biological_by_reference = {
        (feature["recordKey"], feature["biologicalFeatureId"]): feature
        for feature in item["biologicalFeatures"]
    }
    target_references = {
        reference
        for reference, feature in biological_by_reference.items()
        if feature["instanceHash"] == target_instance_hash
    }
    assert len(target_references) == 1
    assert all(
        FEATURE_INSTANCE_HASH_PATTERN.fullmatch(feature["instanceHash"])
        for feature in item["biologicalFeatures"]
    )
    target_rendered = [
        feature
        for feature in item["features"]
        if (feature["recordKey"], feature["biologicalFeatureId"]) in target_references
    ]
    assert target_rendered
    assert {feature["fillColor"] for feature in target_rendered} == {target_color}


def test_circular_render_and_catalog_share_same_record_collision_identity() -> None:
    record = _record(
        "collision",
        "record-a",
        [_feature("first"), _feature("second")],
    )
    target_hash = (
        build_feature_instance_identity_plan(record)
        .identity_for_feature(record.features[1])
        .instance_hash
    )
    exact_color = "#ed1c24"
    options = CircularDiagramOptions(
        selected_features_set=("CDS",),
        colors=ColorOptions(
            color_table=_color_table(
                [
                    (
                        "CDS",
                        FEATURE_INSTANCE_HASH_QUALIFIER,
                        target_hash,
                        exact_color,
                        "target",
                    )
                ]
            )
        ),
        output=CircularOutputOptions(legend="none"),
    )
    item, _prepared = _catalog_for_exact_request(
        CircularDiagramRequest(
            records=(
                RecordInput(
                    source=InMemoryRecordSource(deepcopy(record)),
                    record_key="record-a",
                ),
            ),
            options=options,
        )
    )

    _assert_exact_catalog_fill(
        item,
        target_instance_hash=target_hash,
        target_color=exact_color,
    )


def test_linear_render_and_catalog_scope_duplicate_public_ids_by_record_key() -> None:
    first = _record("duplicate", "record-a", [_feature("first")])
    second = _record("duplicate", "record-b", [_feature("second")])
    target_hash = (
        build_feature_instance_identity_plan([first, second])
        .identity_for_feature(second.features[0])
        .instance_hash
    )
    exact_color = "#ed1c24"
    options = LinearDiagramOptions(
        selected_features_set=("CDS",),
        colors=ColorOptions(
            color_table=_color_table(
                [
                    (
                        "CDS",
                        FEATURE_INSTANCE_HASH_QUALIFIER,
                        target_hash,
                        exact_color,
                        "target",
                    )
                ]
            )
        ),
        output=LinearOutputOptions(legend="none"),
    )
    item, _prepared = _catalog_for_exact_request(
        LinearDiagramRequest(
            records=(
                RecordInput(
                    source=InMemoryRecordSource(deepcopy(first)),
                    record_key="record-a",
                ),
                RecordInput(
                    source=InMemoryRecordSource(deepcopy(second)),
                    record_key="record-b",
                ),
            ),
            options=options,
        )
    )

    _assert_exact_catalog_fill(
        item,
        target_instance_hash=target_hash,
        target_color=exact_color,
    )
