"""The structural checker gives each shared vector's violations (tests/fixtures/contract_shape_vectors.json).

The Web twin runs the same file: tests/web/contract-shape-vectors.test.mjs.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, get_args

import pytest

from gbdraw.contract_shapes import ViolationReason, shape_violations, violations

VECTORS: dict[str, Any] = json.loads(
    (Path(__file__).resolve().parent / "fixtures" / "contract_shape_vectors.json").read_text(encoding="utf-8")
)


def _rows(value: object, shape: dict[str, Any]) -> list[dict[str, str]]:
    return [
        {"path": problem.path, "reason": problem.reason}
        for problem in shape_violations(value, shape, VECTORS["types"], "value")
    ]


def test_vector_ids_are_unique() -> None:
    ids = [vector["id"] for vector in [*VECTORS["vectors"], *VECTORS["jsonDifferences"]]]
    assert len(ids) == len(set(ids))


def test_the_vectors_list_every_reason_and_report_each() -> None:
    assert VECTORS["reasons"] == list(get_args(ViolationReason))
    reported = {row["reason"] for vector in VECTORS["vectors"] for row in vector["violations"]}
    assert reported == set(VECTORS["reasons"])


@pytest.mark.parametrize("vector", VECTORS["vectors"], ids=lambda vector: vector["id"])
def test_the_checker_reports_the_vector_violations(vector: dict[str, Any]) -> None:
    assert _rows(vector["value"], vector["shape"]) == vector["violations"]


@pytest.mark.parametrize("difference", VECTORS["jsonDifferences"], ids=lambda difference: difference["id"])
def test_the_checker_reports_python_s_reading_of_a_json_difference(difference: dict[str, Any]) -> None:
    assert _rows(json.loads(difference["text"]), difference["shape"]) == difference["python"]


def test_a_message_names_the_expectation_and_never_the_value() -> None:
    def messages(value: object, shape: dict[str, Any]) -> list[str]:
        return violations(value, shape, VECTORS["types"], "value")

    secret = "private-genome-name"
    closed = {"kind": "object", "fields": {"a": {"kind": "integer"}}, "optional": [], "open": False}
    assert messages({"extra": secret}, closed) == ["value.extra: unknown field", "value.a: required field is absent"]
    assert messages({"a": secret}, closed) == ["value.a: expected integer"]
    assert messages(secret, {"kind": "literal", "values": [1, "a"]}) == ['value: expected one of [1,"a"]']
    assert messages([secret], {"kind": "tuple", "items": []}) == ["value: expected 0 items"]
    assert messages(None, {"kind": "array", "items": {"kind": "string"}}) == ["value: expected an array"]
    assert messages(0.5, {"kind": "integer"}) == ["value: expected integer"]
