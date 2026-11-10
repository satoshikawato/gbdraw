"""The Python admission of a Linear comparison-plan edge decides every shared vector case."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gbdraw.exceptions import ValidationError
from gbdraw.session_io import _validate_comparison_plan_edge

VECTORS = json.loads(
    (Path(__file__).parent / "fixtures" / "linear_comparison_plan_edge_cases.json").read_text(encoding="utf-8")
)
ALLOWED_FIELDS = {
    "id",
    "queryUid",
    "subjectUid",
    "included",
    "fileActive",
    "losatFilenameActive",
    "source",
    "losatFilename",
    "windowLosatFilenames",
}


@pytest.mark.parametrize("case", VECTORS["cases"], ids=lambda case: case["name"])
def test_python_admission_matches_the_shared_vector(case: dict) -> None:
    edge_ids: set[str] = set()

    def admit() -> None:
        for edge in case["edges"]:
            _validate_comparison_plan_edge(edge, ALLOWED_FIELDS, edge_ids)

    if case["accepted"]:
        admit()
    else:
        with pytest.raises(ValidationError):
            admit()
