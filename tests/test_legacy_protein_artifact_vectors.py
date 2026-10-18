"""Legacy protein artifacts a written Session keeps, shared with the Web (OV-1132, OV-1143).

``tests/fixtures/legacy_protein_artifact_vectors.json`` is read by this file
(``written_legacy_protein_artifacts``) and by
``tests/web/legacy-protein-artifact-vectors.test.mjs``
(``writtenLegacyProteinArtifacts``).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from gbdraw.api.session_compat import written_legacy_protein_artifacts

REPO_ROOT = Path(__file__).resolve().parents[1]
CASES = json.loads(
    (REPO_ROOT / "tests" / "fixtures" / "legacy_protein_artifact_vectors.json").read_text(encoding="utf-8")
)["cases"]


def _candidate(index: int, state: str) -> dict[str, Any]:
    return {
        "schema": 1,
        "state": state,
        "rejectionReason": "no match" if state == "rejected" else None,
        "originalEntry": {"schema": 2, "kind": "raw-losat", "key": f"legacy-{index}", "program": "blastp", "text": ""},
    }


def _evidence(index: int) -> dict[str, Any]:
    return {"schema": 1, "kind": "derived-losatp-payload", "key": f"evidence-{index}", "payload": {}}


@pytest.mark.parametrize("case", CASES, ids=[case["id"] for case in CASES])
def test_written_legacy_protein_artifacts_follow_the_vector(case: dict[str, Any]) -> None:
    candidates = tuple(_candidate(index, state) for index, state in enumerate(case["candidates"]))
    evidence = tuple(_evidence(index) for index in range(case["evidence"]))
    written_candidates, written_evidence = written_legacy_protein_artifacts(
        candidates,
        evidence,
        promoted=frozenset(case["promoted"]),
        has_current_derived=case["currentDerived"],
    )
    kept = [index for index in range(len(candidates)) if index not in case["promoted"]]
    assert [candidate["state"] for candidate in written_candidates] == case["written"]["candidates"]
    assert list(written_candidates) == [candidates[index] for index in kept]
    assert list(written_evidence) == list(evidence[: case["written"]["evidence"]])
