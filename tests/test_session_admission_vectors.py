"""Session admission vectors shared with the Web: what Python's whole Load decides.

``tests/fixtures/session_admission_vectors.json`` is read by this file and by
``tests/web/session-admission-vectors.test.mjs`` (the Web's ``importSession``).
Each vector edits a copy of its base Session and states ``accept`` or
``refuse`` per side. A vector whose sides disagree names its open finding in
``difference``: it records today's behaviour, not a decision.
"""

from __future__ import annotations

import copy
import gzip
import json
from pathlib import Path
from typing import Any

import pytest

from gbdraw.exceptions import GbdrawError
from gbdraw.session_io import load_session

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE = json.loads(
    (REPO_ROOT / "tests" / "fixtures" / "session_admission_vectors.json").read_text(encoding="utf-8")
)
VECTORS = FIXTURE["vectors"]


def _base(name: str) -> dict[str, Any]:
    base = FIXTURE["bases"][name]
    if "document" in base:
        return copy.deepcopy(base["document"])
    data = (REPO_ROOT / base["fixture"]).read_bytes()
    return json.loads(gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data)


def _edited(vector: dict[str, Any]) -> dict[str, Any]:
    document = _base(vector["base"])
    for edit in vector["edits"]:
        *parents, last = (
            token.replace("~1", "/").replace("~0", "~") for token in edit["path"].split("/")[1:]
        )
        parent: Any = document
        for token in parents:
            parent = parent[int(token)] if isinstance(parent, list) else parent[token]
        if edit["op"] == "set":
            parent[last] = edit["value"]
        else:
            assert edit["op"] == "remove", edit
            del parent[last]
    return document


def test_vector_ids_are_unique_and_a_disagreement_names_its_finding() -> None:
    assert len({vector["id"] for vector in VECTORS}) == len(VECTORS)
    for vector in VECTORS:
        assert {vector["python"], vector["web"]} <= {"accept", "refuse"}, vector["id"]
        assert (vector["python"] != vector["web"]) == ("difference" in vector), vector["id"]
        # A refusal names a fragment of its message; without one any error would pass.
        assert (vector["python"] == "refuse") == bool(vector.get("pythonMessage")), vector["id"]


@pytest.mark.parametrize("vector", VECTORS, ids=[vector["id"] for vector in VECTORS])
def test_python_load_gives_the_vector_verdict(vector: dict[str, Any], tmp_path: Path) -> None:
    path = tmp_path / "session.json"
    path.write_text(json.dumps(_edited(vector)), encoding="utf-8")
    if vector["python"] == "accept":
        load_session(path)
        return
    with pytest.raises(GbdrawError) as refused:
        load_session(path)
    assert vector["pythonMessage"] in str(refused.value)
