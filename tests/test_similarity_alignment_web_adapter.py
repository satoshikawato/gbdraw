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
        "schema": 1,
        "mode": "position",
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
                "recordLength": 200,
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
        "schema": 1,
        "mode": "position",
        "groupId": "og-1",
        "reference": _anchor("record-a", "reference", 0),
        "records": [
            {
                "recordKey": "record-a",
                "status": "reference",
                "rationale": "reference",
                "anchor": _anchor("record-a", "reference", 0),
                "effectiveReverseComplement": None,
            },
            {
                "recordKey": "record-b",
                "status": "aligned",
                "rationale": "only_usable_candidate",
                "anchor": _anchor("record-b", "target", 1),
                "effectiveReverseComplement": None,
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
        {"recordKey": "record-b", "kind": "skip", "anchor": None}
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
