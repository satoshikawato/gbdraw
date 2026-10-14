"""The structural checker gives each shared vector's violations (tests/fixtures/contract_shape_vectors.json)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from gbdraw.contract_shapes import violations

VECTORS: dict[str, Any] = json.loads(
    (Path(__file__).resolve().parent / "fixtures" / "contract_shape_vectors.json").read_text(encoding="utf-8")
)


def test_vector_ids_are_unique() -> None:
    ids = [vector["id"] for vector in VECTORS["vectors"]]
    assert len(ids) == len(set(ids))


@pytest.mark.parametrize("vector", VECTORS["vectors"], ids=lambda vector: vector["id"])
def test_the_checker_reports_the_vector_violations(vector: dict[str, Any]) -> None:
    assert violations(vector["value"], vector["shape"], VECTORS["types"], "value") == vector["violations"]
