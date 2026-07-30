from __future__ import annotations

import json

import pytest

import gbdraw.web_support.feature_catalog as feature_catalog_module
from gbdraw.exceptions import GbdrawError
from gbdraw.render.interactive_svg import InteractiveSvgContext
from gbdraw.web_support.feature_catalog import (
    build_feature_catalog_item,
    canonical_catalog_sequence_sources,
    materialize_catalog_nucleotide_sequence,
    select_feature_catalog_item,
)


def test_feature_catalog_normalizes_references_and_sequence_ownership() -> None:
    full_note = ("x" * 49) + "😀tail"
    compact_note = ("x" * 49) + "😀"
    context = InteractiveSvgContext(
        record_keys=("stable-record-key",),
        features=(
            {
                "svg_id": "stable-feature",
                "stable_feature_id": "stable-feature",
                "rendered_feature_svg_id": "rendered-feature",
                "record_id": "public-record-id",
                "record_idx": 0,
                "type": "CDS",
                "qualifiers": {
                    "locus_tag": ["LOCUS_001"],
                    "translation": ["MK"],
                },
                "nucleotide_sequence": "ATGAAATAA",
                "amino_acid_sequence": "MK",
            },
        ),
        biological_features=(
            {
                "id": "f0",
                "stable_feature_id": "source-feature",
                "record_id": "public-record-id",
                "record_idx": 0,
                "feature_index": 0,
                "type": "source",
                "qualifiers": {"organism": ["Example organism"]},
                "nucleotide_sequence": "ATGAAATAAGGG",
            },
            {
                "id": "f1",
                "svg_id": "stable-feature",
                "stable_svg_id": "stable-feature",
                "stable_feature_id": "stable-feature",
                "record_id": "public-record-id",
                "record_idx": 0,
                "feature_index": 1,
                "type": "CDS",
                "start": 0,
                "end": 9,
                "strand": "+",
                "location_parts": [
                    {
                        "start": 0,
                        "end": 9,
                        "strand": "+",
                        "display": "1..9",
                    }
                ],
                "protein_id": "PROTEIN_001",
                "locus_tag": "LOCUS_001",
                "old_locus_tag": "OLD_LOCUS_001",
                "gene": "example_gene",
                "product": "Example enzyme",
                "note": compact_note,
                "source_protein_id": "PROTEIN_001",
                "qualifiers": {
                    "protein_id": ["PROTEIN_001"],
                    "locus_tag": ["LOCUS_001"],
                    "old_locus_tag": ["OLD_LOCUS_001"],
                    "gene": ["example_gene"],
                    "product": ["Example enzyme"],
                    "note": [full_note],
                    "translation": ["MK"],
                },
                "selector": {
                    "hash": "stable-feature",
                    "qualifiers": {
                        "locus_tag": ["LOCUS_001"],
                        "translation": ["MK"],
                    },
                },
                "nucleotide_sequence": "ATGAAATAA",
                "amino_acid_sequence": "MK",
            },
        ),
        orthogroups=(
            {
                "id": "og-1",
                "member_count": 1,
                "record_coverage_count": 1,
                "orthologEdges": [
                    {
                        "orthogroupId": "og-1",
                        "queryProteinId": "PROTEIN_001",
                        "subjectProteinId": "PROTEIN_001",
                    }
                ],
                "orthologPaths": [
                    {
                        "pathId": "og-1.path-1",
                        "proteinIds": ["PROTEIN_001"],
                    }
                ],
                "relatedEdges": [
                    {
                        "edgeKind": "weak_bridge",
                        "queryProteinId": "PROTEIN_001",
                    }
                ],
                "members": [
                    {
                        "orthogroupId": "og-1",
                        "recordIndex": 0,
                        "featureIndex": 1,
                        "stableFeatureSvgId": "stable-feature",
                        "proteinId": "PROTEIN_001",
                        "representative": True,
                        "role": "anchor",
                        "confidence": "high",
                        "assignment_reason": "reciprocal best hit",
                        "supportingEdges": ["internal-edge-id"],
                        "bestCoreSupport": 0.0,
                        "secondBestCoreSupport": 0.0,
                        "start": 0,
                        "end": 9,
                        "gene": "duplicated gene",
                        "product": "duplicated product",
                    }
                ],
            },
        ),
        sequence_sources=(
            {
                "key": "linear:record:0",
                "recordId": "public-record-id",
                "sequence": "ATGAAATAAGGG",
                "origin": "linear-record",
                "recordIndex": 0,
            },
            {
                "key": "linear:record:1",
                "recordId": "unused-record",
                "sequence": "CCCC",
                "origin": "linear-record",
                "recordIndex": 1,
            },
        ),
    )
    svg = """\
<svg xmlns="http://www.w3.org/2000/svg">
  <path id="rendered-feature"
        data-gbdraw-feature-id="rendered-feature"
        data-gbdraw-stable-feature-id="stable-feature"
        data-gbdraw-record-index="0"
        fill="#123456" />
  <path id="comparison-match"
        data-gbdraw-match-id="match-1"
        data-match-kind="pairwise"
        data-query-record-index="0"
        data-subject-record-index="0"
        data-query-feature-svg-id="rendered-feature"
        data-subject-feature-svg-id="rendered-feature"
        data-query-protein-id="PROTEIN_001"
        data-subject-protein-id="PROTEIN_001" />
</svg>
"""

    item = build_feature_catalog_item(
        svg,
        context,
        result_index=2,
        result_name="result.svg",
    )

    assert item["resultIndex"] == 2
    assert item["resultName"] == "result.svg"
    assert item["recordKeys"] == ["stable-record-key"]
    assert [feature["type"] for feature in item["biologicalFeatures"]] == ["CDS"]

    biological = item["biologicalFeatures"][0]
    biological_id = biological["biologicalFeatureId"]
    assert biological["recordKey"] == "stable-record-key"
    assert biological_id == "stable-feature"
    assert "stableFeatureId" not in biological
    assert "record_idx" not in biological
    assert "translation" not in biological["qualifiers"]
    assert "selector" not in biological
    assert "source_protein_id" not in biological
    assert "protein_id" not in biological
    assert "locus_tag" not in biological
    assert "old_locus_tag" not in biological
    assert "gene" not in biological
    assert "product" not in biological
    assert "note" not in biological
    assert "location_parts" not in biological
    assert biological["translationFromAminoAcidSequence"] is True
    assert biological["sequenceSourceIndex"] == 0
    assert "nucleotide_sequence" not in biological

    assert item["features"] == [
        {
            "svgId": "rendered-feature",
            "recordKey": "stable-record-key",
            "biologicalFeatureId": biological_id,
            "fillColor": "#123456",
        }
    ]
    group = item["orthogroups"][0]
    assert group["orthologEdgeCount"] == 1
    assert group["orthologPathCount"] == 1
    assert group["relatedEdgeCount"] == 1
    assert "orthologEdges" not in group
    assert "orthologPaths" not in group
    assert "relatedEdges" not in group
    member = group["members"][0]
    assert member["recordKey"] == "stable-record-key"
    assert member["biologicalFeatureId"] == biological_id
    assert member["representative"] is True
    assert member["assignmentReason"] == "reciprocal best hit"
    assert "assignment_reason" not in member
    assert "supportingEdges" not in member
    assert "bestCoreSupport" not in member
    assert "secondBestCoreSupport" not in member
    assert "start" not in member
    assert "end" not in member
    assert "gene" not in member
    assert "product" not in member

    match = item["comparisonMatches"][0]
    assert match["queryBiologicalFeatureId"] == biological_id
    assert match["subjectBiologicalFeatureId"] == biological_id
    assert match["queryRecordKey"] == "stable-record-key"
    assert match["subjectRecordKey"] == "stable-record-key"
    assert "query_feature_svg_id" not in match
    assert "subject_feature_svg_id" not in match
    assert "query_protein_id" not in match
    assert "subject_protein_id" not in match
    assert [source["key"] for source in item["sequenceSources"]] == [
        "linear:record:0"
    ]

    encoded = json.dumps(item)
    assert encoded.count('"nucleotide_sequence": "ATGAAATAA"') == 0
    assert encoded.count('"sequence": "ATGAAATAAGGG"') == 1
    assert encoded.count('"MK"') == 1
    assert (
        materialize_catalog_nucleotide_sequence(
            biological,
            item["sequenceSources"],
        )
        == "ATGAAATAA"
    )
    invalid_item = json.loads(json.dumps(item))
    invalid_item["biologicalFeatures"][0]["end"] = 99
    with pytest.raises(
        GbdrawError,
        match="invalid nucleotide sequence source reference",
    ):
        select_feature_catalog_item(
            {"schema": 3, "items": [invalid_item]},
            result_index=2,
            result_name="result.svg",
        )
    for conflict in ("DIFFERENT", "", " ", None, [], [None], 0, False):
        invalid_translation = json.loads(json.dumps(item))
        invalid_translation["biologicalFeatures"][0]["qualifiers"][
            "translation"
        ] = conflict
        with pytest.raises(
            GbdrawError,
            match="invalid translation source reference",
        ):
            select_feature_catalog_item(
                {"schema": 3, "items": [invalid_translation]},
                result_index=2,
                result_name="result.svg",
            )
    for invalid_amino_acid in (0, False, {}, []):
        invalid_translation = json.loads(json.dumps(item))
        invalid_translation["biologicalFeatures"][0][
            "amino_acid_sequence"
        ] = invalid_amino_acid
        with pytest.raises(
            GbdrawError,
            match="invalid translation source reference",
        ):
            select_feature_catalog_item(
                {"schema": 3, "items": [invalid_translation]},
                result_index=2,
                result_name="result.svg",
            )
    camel_amino_acid = json.loads(json.dumps(item))
    camel_feature = camel_amino_acid["biologicalFeatures"][0]
    camel_feature["aminoAcidSequence"] = camel_feature.pop(
        "amino_acid_sequence"
    )
    select_feature_catalog_item(
        {"schema": 3, "items": [camel_amino_acid]},
        result_index=2,
        result_name="result.svg",
    )
    for shadowing_value in (None, ""):
        shadowed_amino_acid = json.loads(json.dumps(camel_amino_acid))
        shadowed_amino_acid["biologicalFeatures"][0][
            "amino_acid_sequence"
        ] = shadowing_value
        with pytest.raises(
            GbdrawError,
            match="invalid translation source reference",
        ):
            select_feature_catalog_item(
                {"schema": 3, "items": [shadowed_amino_acid]},
                result_index=2,
                result_name="result.svg",
            )
    for invalid_sequence in (123, "ATG AAA"):
        invalid_source = json.loads(json.dumps(item))
        invalid_source["sequenceSources"][0]["sequence"] = invalid_sequence
        with pytest.raises(
            GbdrawError,
            match="invalid nucleotide sequence source reference",
        ):
            select_feature_catalog_item(
                {"schema": 3, "items": [invalid_source]},
                result_index=2,
                result_name="result.svg",
            )
    unreferenced_invalid_source = json.loads(json.dumps(item))
    unreferenced_invalid_source["sequenceSources"].append(
        {
            "origin": "linear-record",
            "recordIndex": 0,
            "sequence": "ATG AAA",
        }
    )
    with pytest.raises(
        GbdrawError,
        match="invalid nucleotide sequence source reference",
    ):
        select_feature_catalog_item(
            {"schema": 3, "items": [unreferenced_invalid_source]},
            result_index=2,
            result_name="result.svg",
        )
    coexisting_sequence = json.loads(json.dumps(item))
    coexisting_sequence["biologicalFeatures"][0][
        "nucleotide_sequence"
    ] = "ATGAAATAA"
    with pytest.raises(
        GbdrawError,
        match="invalid nucleotide sequence source reference",
    ):
        select_feature_catalog_item(
            {"schema": 3, "items": [coexisting_sequence]},
            result_index=2,
            result_name="result.svg",
        )


def test_feature_catalog_references_only_exactly_reconstructible_sequences() -> None:
    source_sequence = "AAAACCCCGGGGTTTT"
    biological_features = (
        {
            "stable_feature_id": "plus",
            "record_idx": 0,
            "feature_index": 0,
            "type": "CDS",
            "start": 0,
            "end": 4,
            "strand": "+",
            "nucleotide_sequence": "AAAA",
        },
        {
            "stable_feature_id": "minus",
            "record_idx": 0,
            "feature_index": 1,
            "type": "CDS",
            "start": 4,
            "end": 8,
            "strand": "-",
            "nucleotide_sequence": "GGGG",
        },
        {
            "stable_feature_id": "multipart",
            "record_idx": 0,
            "feature_index": 2,
            "type": "CDS",
            "start": 0,
            "end": 16,
            "strand": "+",
            "location_parts": [
                {"start": 0, "end": 2, "strand": "+"},
                {"start": 12, "end": 16, "strand": "-"},
            ],
            "nucleotide_sequence": "AAAAAA",
        },
        {
            "stable_feature_id": "mismatch",
            "record_idx": 0,
            "feature_index": 3,
            "type": "CDS",
            "start": 8,
            "end": 12,
            "strand": "+",
            "nucleotide_sequence": "TTTT",
            "amino_acid_sequence": "MK",
            "qualifiers": {"translation": ["DIFFERENT"]},
        },
    )
    context = InteractiveSvgContext(
        record_keys=("record-key",),
        features=(biological_features[0],),
        biological_features=biological_features,
        sequence_sources=(
            {
                "key": "linear:record:0",
                "recordId": "record-id",
                "sequence": source_sequence,
                "origin": "linear-record",
                "recordIndex": 0,
            },
        ),
    )
    svg = """\
<svg xmlns="http://www.w3.org/2000/svg">
  <path id="plus" data-gbdraw-feature-id="plus"
        data-gbdraw-stable-feature-id="plus" data-gbdraw-record-index="0" />
  <path id="match" data-gbdraw-match-id="match-1"
        data-match-kind="pairwise" data-query-record-index="0"
        data-subject-record-index="0" data-query-feature-svg-id="plus"
        data-subject-feature-svg-id="plus" />
</svg>
"""

    item = build_feature_catalog_item(
        svg,
        context,
        result_index=0,
        result_name="diagram.svg",
    )

    by_id = {
        feature["biologicalFeatureId"]: feature
        for feature in item["biologicalFeatures"]
    }
    for feature_id, expected in (
        ("plus", "AAAA"),
        ("minus", "GGGG"),
        ("multipart", "AAAAAA"),
    ):
        feature = by_id[feature_id]
        assert feature["sequenceSourceIndex"] == 0
        assert "nucleotide_sequence" not in feature
        assert (
            materialize_catalog_nucleotide_sequence(
                feature,
                item["sequenceSources"],
            )
            == expected
        )

    mismatch = by_id["mismatch"]
    assert "sequenceSourceIndex" not in mismatch
    assert mismatch["nucleotide_sequence"] == "TTTT"
    assert "translationFromAminoAcidSequence" not in mismatch
    assert mismatch["qualifiers"]["translation"] == ["DIFFERENT"]


def test_feature_catalog_disambiguates_identical_biological_feature_ids() -> None:
    context = InteractiveSvgContext(
        record_keys=("record-key",),
        biological_features=(
            {
                "stable_feature_id": "duplicate",
                "record_idx": 0,
                "feature_index": 4,
                "type": "CDS",
            },
            {
                "stable_feature_id": "duplicate",
                "record_idx": 0,
                "feature_index": 7,
                "type": "CDS",
            },
        ),
    )

    item = build_feature_catalog_item(
        '<svg xmlns="http://www.w3.org/2000/svg" />',
        context,
        result_index=0,
        result_name="empty.svg",
    )

    identities = [
        feature["biologicalFeatureId"]
        for feature in item["biologicalFeatures"]
    ]
    assert len(set(identities)) == 2
    assert all(identity.startswith("duplicate~") for identity in identities)


def test_feature_catalog_references_one_biological_payload_from_many_renderings() -> None:
    rendered_count = 128
    stable_id = "stable-repeated-feature"
    sequence = "ATG" + ("ACGT" * 128) + "TAA"
    amino_acid = "M" + ("ACDEFGHIKLMNPQRSTVWY" * 8)
    context = InteractiveSvgContext(
        record_keys=("record-key",),
        features=tuple(
            {
                "stable_feature_id": stable_id,
                "rendered_feature_svg_id": f"rendered-{index}",
                "record_idx": 0,
                "feature_index": 17,
                "type": "CDS",
                "qualifiers": {
                    "product": ["Repeated protein"],
                    "translation": [amino_acid],
                },
                "nucleotide_sequence": sequence,
                "amino_acid_sequence": amino_acid,
            }
            for index in range(rendered_count)
        ),
        biological_features=(
            {
                "stable_feature_id": stable_id,
                "record_idx": 0,
                "feature_index": 17,
                "type": "CDS",
                "qualifiers": {
                    "product": ["Repeated protein"],
                    "translation": [amino_acid],
                },
                "nucleotide_sequence": sequence,
                "amino_acid_sequence": amino_acid,
            },
        ),
    )
    paths = "\n".join(
        (
            f'<path id="rendered-{index}" '
            f'data-gbdraw-feature-id="rendered-{index}" '
            f'data-gbdraw-stable-feature-id="{stable_id}" '
            'data-gbdraw-record-index="0" fill="#123456" />'
        )
        for index in range(rendered_count)
    )

    item = build_feature_catalog_item(
        f'<svg xmlns="http://www.w3.org/2000/svg">{paths}</svg>',
        context,
        result_index=0,
        result_name="repeated.svg",
    )

    assert len(item["features"]) == rendered_count
    assert len(item["biologicalFeatures"]) == 1
    references = {
        (feature["recordKey"], feature["biologicalFeatureId"])
        for feature in item["features"]
    }
    biological = item["biologicalFeatures"][0]
    assert references == {
        (biological["recordKey"], biological["biologicalFeatureId"])
    }
    assert all(
        set(feature)
        <= {"svgId", "recordKey", "biologicalFeatureId", "fillColor"}
        for feature in item["features"]
    )
    encoded = json.dumps(item)
    assert encoded.count(sequence) == 1
    assert encoded.count(amino_acid) == 1
    assert encoded.count('"translation"') == 0


def test_sequence_source_validation_is_bounded_by_source_count(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_sequence = "ATG" + ("A" * 1_000_000)
    source = {
        "origin": "linear-record",
        "recordIndex": 0,
        "sequence": source_sequence,
    }
    features = [
        {
            "recordKey": "record-key",
            "biologicalFeatureId": f"feature-{index}",
            "start": 0,
            "end": 3,
            "strand": "+",
            "nucleotide_sequence": "ATG",
        }
        for index in range(512)
    ]
    validation_calls = 0
    original = feature_catalog_module._canonical_sequence_source

    def counted(source_payload: object) -> str | None:
        nonlocal validation_calls
        validation_calls += 1
        return original(source_payload)  # type: ignore[arg-type]

    monkeypatch.setattr(
        feature_catalog_module,
        "_canonical_sequence_source",
        counted,
    )

    feature_catalog_module._reference_reconstructible_nucleotide_sequences(
        features,
        ["record-key"],
        [source],
    )
    assert validation_calls == 1
    assert all(feature["sequenceSourceIndex"] == 0 for feature in features)

    validation_calls = 0
    catalog = {
        "schema": 3,
        "items": [
            {
                "resultIndex": 0,
                "resultName": "diagram.svg",
                "recordKeys": ["record-key"],
                "features": [],
                "biologicalFeatures": features,
                "orthogroups": [],
                "annotations": [],
                "comparisonMatches": [],
                "sequenceSources": [source],
            }
        ],
    }
    select_feature_catalog_item(
        catalog,
        result_index=0,
        result_name="diagram.svg",
    )
    assert validation_calls == 1

    validation_calls = 0
    canonical_sources = canonical_catalog_sequence_sources([source])
    assert validation_calls == 1
    for feature in features:
        assert materialize_catalog_nucleotide_sequence(
            feature,
            [source],
            canonical_sources=canonical_sources,
        ) == "ATG"
    assert validation_calls == 1

    with pytest.raises(
        GbdrawError,
        match="invalid nucleotide sequence source",
    ):
        feature_catalog_module._reference_reconstructible_nucleotide_sequences(
            features,
            ["record-key"],
            [
                source,
                {
                    "origin": "linear-record",
                    "recordIndex": 0,
                    "sequence": "ATG AAA",
                },
            ],
        )
