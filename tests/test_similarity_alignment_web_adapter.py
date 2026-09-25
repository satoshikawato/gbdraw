from __future__ import annotations

import json

import pytest

from gbdraw.exceptions import ValidationError
from gbdraw.web_support.similarity_alignment import (
    resolve_similarity_alignment_json,
    resolve_similarity_alignment_payload,
)


def _anchor(record_key: str, feature_id: str, index: int) -> dict[str, object]:
    return {
        "recordKey": record_key,
        "biologicalFeatureId": feature_id,
        "sourceFeatureIndex": index,
        "stableFeatureSvgId": f"stable-{feature_id}",
    }


def _member(
    record_key: str,
    feature_id: str,
    index: int,
    *,
    start: int,
    end: int,
    strand: int = 1,
) -> dict[str, object]:
    return {
        "groupId": "og-1",
        "anchor": _anchor(record_key, feature_id, index),
        "sourceStart": start,
        "sourceEnd": end,
        "sourceStrand": strand,
        "identityIsUnique": True,
        "hidden": False,
        "representative": False,
        "role": "inparalog",
    }


def _request() -> dict[str, object]:
    reference = _anchor("record-a", "reference", 0)
    return {
        "schema": 2,
        "groupId": "og-1",
        "records": [
            {
                "recordKey": "record-a",
                "recordLength": 1000,
                "region": None,
                "presentation": {"reverseComplement": False},
            },
            {
                "recordKey": "record-b",
                "recordLength": 400,
                "region": {"start": 101, "end": 300, "reverseComplement": True},
                "presentation": {"reverseComplement": False},
            },
        ],
        "reference": reference,
        "members": [
            _member("record-a", "reference", 0, start=10, end=40),
            _member("record-b", "target", 1, start=140, end=170, strand=-1),
        ],
        "directEdges": [
            {
                "groupId": "og-1",
                "query": reference,
                "subject": _anchor("record-b", "target", 1),
                "edgeKind": "rbh",
            }
        ],
        "choices": [],
    }


def test_web_adapter_projects_crop_facts_and_returns_canonical_plan() -> None:
    result = resolve_similarity_alignment_payload(_request())

    assert result["status"] == "resolved"
    assert result["plan"] == {
        "schema": 2,
        "groupId": "og-1",
        "reference": _anchor("record-a", "reference", 0),
        "records": [
            {
                "recordKey": "record-a",
                "status": "reference",
                "rationale": "reference",
                "anchor": _anchor("record-a", "reference", 0),
                "orientationPolicy": "preserve",
                "effectiveReverseComplement": None,
            },
            {
                "recordKey": "record-b",
                "status": "aligned",
                "rationale": "only_usable_candidate",
                "anchor": _anchor("record-b", "target", 1),
                "orientationPolicy": "preserve",
                "effectiveReverseComplement": True,
            },
        ],
    }
    assert json.loads(resolve_similarity_alignment_json(json.dumps(_request()))) == result


def test_web_adapter_returns_ambiguity_then_uses_explicit_skip() -> None:
    request = _request()
    request["directEdges"] = []
    request["members"].append(  # type: ignore[union-attr]
        _member("record-b", "other", 2, start=190, end=220)
    )

    ambiguous = resolve_similarity_alignment_payload(request)
    assert ambiguous["status"] == "ambiguous"
    assert ambiguous["plan"] is None
    assert ambiguous["records"][1]["kind"] == "ambiguous"  # type: ignore[index]

    request["choices"] = [
        {"recordKey": "record-b", "kind": "skip", "anchor": None,
         "orientationPolicy": "preserve"}
    ]
    resolved = resolve_similarity_alignment_payload(request)
    assert resolved["status"] == "resolved"
    assert resolved["plan"]["records"][1]["rationale"] == "skipped_by_user"  # type: ignore[index]


def test_web_adapter_rejects_unknown_fields_and_unmappable_reference() -> None:
    unknown = _request()
    unknown["unexpected"] = True
    with pytest.raises(ValidationError, match="invalid fields"):
        resolve_similarity_alignment_payload(unknown)

    cropped = _request()
    cropped["records"][0]["region"] = {  # type: ignore[index]
        "start": 500,
        "end": 600,
        "reverseComplement": False,
    }
    with pytest.raises(ValidationError, match="exact reference.*current crop"):
        resolve_similarity_alignment_payload(cropped)


def test_web_projection_discloses_candidate_evidence_and_orientation_facts() -> None:
    request = _request()
    request["members"][1]["displayName"] = "target gene"  # type: ignore[index]
    result = resolve_similarity_alignment_payload(request)
    row = result["records"][1]
    assert row["reviewReason"] == "only_usable_candidate"
    assert row["anchor"] == _anchor("record-b", "target", 1)
    assert row["orientationPolicy"] == "preserve"
    candidate = row["candidates"][0]
    assert candidate["displayName"] == "target gene"
    assert (candidate["sourceStart"], candidate["sourceEnd"]) == (140, 170)
    assert candidate["displayCenter"] == 145.0
    assert candidate["displayedStrand"] == 1
    assert candidate["representative"] is False
    assert candidate["directEvidence"] == ["rbh"]
    assert candidate["orientation"] == {
        "preserve": {"effect": "preserve", "effectiveReverseComplement": True},
        "match_reference": {"effect": "preserve", "effectiveReverseComplement": True},
    }
    assert result["referenceDisplayedStrand"] == 1
    assert result["referenceDisplayCenter"] == 25.0


def test_web_recommendation_is_transient_and_explicit_match_is_validated() -> None:
    request = _request()
    request["directEdges"] = []
    request["members"][1]["sourceStrand"] = 1  # type: ignore[index]
    replacement = _member(
        "record-b", "other", 2, start=190, end=220, strand=1
    )
    replacement["representative"] = True
    request["members"].append(replacement)  # type: ignore[union-attr]
    unresolved = resolve_similarity_alignment_payload(request)
    row = unresolved["records"][1]
    assert row["recommendationReason"] == "unique_representative"
    assert row["recommendedAnchor"] == _anchor("record-b", "other", 2)
    assert row["candidates"][0]["orientation"]["match_reference"] == {
        "effect": "reverse_whole_record",
        "effectiveReverseComplement": False,
    }
    assert unresolved["plan"] is None

    request["choices"] = [{
        "recordKey": "record-b", "kind": "select",
        "anchor": _anchor("record-b", "target", 1),
        "orientationPolicy": "match_reference",
    }]
    selected = resolve_similarity_alignment_payload(request)
    assert selected["plan"]["schema"] == 2
    assert selected["plan"]["records"][1]["anchor"] == _anchor("record-b", "target", 1)
    assert selected["plan"]["records"][1]["orientationPolicy"] == "match_reference"
    assert selected["plan"]["records"][1]["effectiveReverseComplement"] is False


def test_web_rejects_incomplete_choice_and_invalid_crop() -> None:
    request = _request()
    request["choices"] = [{
        "recordKey": "record-b", "kind": "skip", "anchor": None,
    }]
    with pytest.raises(ValidationError, match="invalid fields"):
        resolve_similarity_alignment_payload(request)
    request["choices"][0]["orientationPolicy"] = "match_reference"  # type: ignore[index]
    with pytest.raises(ValidationError, match="Skip choice must preserve"):
        resolve_similarity_alignment_payload(request)

    cropped = _request()
    cropped["records"][1]["recordLength"] = 200  # type: ignore[index]
    with pytest.raises(ValidationError, match="exceeds recordLength"):
        resolve_similarity_alignment_payload(cropped)


def test_web_crop_after_presentation_reverse_uses_source_coordinates() -> None:
    request = _request()
    target_record = request["records"][1]
    target_record["presentation"]["reverseComplement"] = True
    target_record["region"]["reverseComplement"] = False
    request["members"][1]["sourceStrand"] = 1
    result = resolve_similarity_alignment_payload(request)
    candidate = result["records"][1]["candidates"][0]
    # Source center 155 maps to 245 after the 400-bp presentation reversal,
    # then to 145 in the 101..300 crop.
    assert candidate["displayCenter"] == 145.0
    assert candidate["displayedStrand"] == -1
    assert candidate["orientation"]["match_reference"] == {
        "effect": "reverse_whole_record",
        "effectiveReverseComplement": False,
    }

    target_record["region"]["reverseComplement"] = True
    with pytest.raises(ValidationError, match="both region and presentation"):
        resolve_similarity_alignment_payload(request)
