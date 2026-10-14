"""Session admission vectors shared with the Web: what Python's whole Load decides.

``tests/fixtures/session_admission_vectors.json`` is read by this file and by
``tests/web/session-admission-vectors.test.mjs`` (the Web's ``importSession``).
Each vector edits a copy of its base Session and states Python's verdict
(``accept`` or ``refuse``) and the Web's (``loaded``, ``repaired`` or
``refused``). A vector Python accepts and the Web refuses names its open
finding in ``difference``: it records today's behaviour, not a decision.
"""

from __future__ import annotations

import copy
import gzip
import json
from pathlib import Path
from typing import Any

import pytest

from gbdraw.contract_shapes import shape_violations
from gbdraw.exceptions import GbdrawError
from gbdraw.session_contract import describe_session_admission
from gbdraw.session_io import FLAT_DRAFT_TOP_LEVEL_FIELDS, load_session

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
            # A copy: a later edit of the same vector may write inside this value.
            parent[last] = copy.deepcopy(edit["value"])
        else:
            assert edit["op"] == "remove", edit
            del parent[last]
    return document


# R1-R4: the only refusals the Load repair receipt keeps (BA-D25).
WEB_CLASSES = {"R1", "R2", "R3", "R4"}


def test_vector_ids_are_unique_and_each_verdict_is_complete() -> None:
    assert len({vector["id"] for vector in VECTORS}) == len(VECTORS)
    for vector in VECTORS:
        assert vector["python"] in {"accept", "refuse"}, vector["id"]
        assert vector["web"] in {"loaded", "repaired", "refused"}, vector["id"]
        # A refusal names a fragment of its message; without one any error would pass.
        assert (vector["python"] == "refuse") == bool(vector.get("pythonMessage")), vector["id"]
        assert ("webClass" in vector) == (vector["web"] == "refused"), vector["id"]
        assert ("webRepairs" in vector) == (vector["web"] == "repaired"), vector["id"]
        assert vector.get("webRepairs", ["code@path"]) and all(
            "@" in repair for repair in vector.get("webRepairs", [])
        ), vector["id"]


def test_the_web_refuses_only_with_a_class_and_repairs_only_what_python_does_not_read() -> None:
    web_only = {f"{row['code']}@{row['path']}" for row in FIXTURE["webOnlyRepairs"]}
    assert all(set(row) == {"code", "path", "finding"} for row in FIXTURE["webOnlyRepairs"])
    for vector in VECTORS:
        # S1: what Python accepts, the Web does not refuse (or the vector names the finding).
        assert (vector["python"] == "accept" and vector["web"] == "refused") == ("difference" in vector), vector["id"]
        # S2
        assert vector["web"] != "refused" or vector["webClass"] in WEB_CLASSES, vector["id"]
        # S3
        if vector["python"] == "accept" and vector["web"] == "repaired":
            assert set(vector["webRepairs"]) <= web_only, vector["id"]


ADMISSION = describe_session_admission()
ADMISSION_TYPES: dict[str, Any] = ADMISSION["types"]
BAND_ROOT = {version: band["root"] for band in ADMISSION["bands"] for version in band["versions"]}
Row = tuple[str, str]  # (object type, field), or (object type, "*") for another key of a closed object
# Rows inside `otherModeResult`: no stored Session holds one, so no vector base
# does; the soundness test of tests/test_session_contract.py covers them.
ROWS_WITHOUT_A_VECTOR = {
    "CheckedSessionOtherModeResult", "CheckedSessionOtherModeEditorState", "CheckedSessionOtherModeLegend",
    "CheckedSessionOtherModeUi", "CheckedSessionRunMetadata",
}


def _description_rows() -> set[Row]:
    """Each check of the admission description: a field whose kind is declared
    or that is required, and another key of a closed object."""

    rows: set[Row] = set()
    for name, shape in ADMISSION_TYPES.items():
        if shape["kind"] != "object":
            continue
        rows |= {
            (name, key) for key, field in shape["fields"].items()
            if field["kind"] != "unknown" or key not in shape["optional"]
        }
        if not shape["open"]:
            rows.add((name, "*"))
    return rows


def _rows_left(value: object, shape: dict[str, Any], owner: str, field: str) -> set[Row]:
    """The rows of the description that ``value`` leaves, as ``field`` of ``owner``."""

    kind = shape["kind"]
    if kind == "ref":
        target = ADMISSION_TYPES[shape["name"]]
        if target["kind"] != "object":
            return _rows_left(value, target, owner, field)
        return _rows_left(value, target, shape["name"], "") if isinstance(value, dict) else {(owner, field)}
    if kind == "union":
        found = [shape_violations(value, option, ADMISSION_TYPES, "") for option in shape["options"]]
        if not all(found):
            return set()
        option = shape["options"][found.index(min(found, key=len))]
        nested = option["kind"] in {"ref", "array", "record"} and isinstance(value, (dict, list))
        return _rows_left(value, option, owner, field) if nested else {(owner, field)}
    if kind in {"array", "record"}:
        if not isinstance(value, list if kind == "array" else dict):
            return {(owner, field)}
        items = value if isinstance(value, list) else list(value.values())
        item_shape = shape["items"] if kind == "array" else shape["values"]
        return set().union(*(_rows_left(item, item_shape, owner, field) for item in items))
    if kind == "object":
        rows = {(owner, "*") for key in value if key not in shape["fields"]} if not shape["open"] else set()
        rows |= {(owner, key) for key in shape["fields"] if key not in value and key not in shape["optional"]}
        return rows.union(*(
            _rows_left(item, shape["fields"][key], owner, key) for key, item in value.items() if key in shape["fields"]
        ))
    return {(owner, field)} if shape_violations(value, shape, ADMISSION_TYPES, "") else set()


def _holds_nothing(value: object) -> bool:
    return value is None or value == [] or value == {}


def test_every_row_of_the_admission_description_has_a_vector() -> None:
    """PLAN 9.3: each check the description makes is one Python refuses and
    the Web refuses or repairs in at least one vector (K13: none refuses a
    file the Web loads unchanged)."""

    left: dict[Row, list[str]] = {}
    for vector in VECTORS:
        document = _edited(vector)
        version = document.get("version")
        root = BAND_ROOT.get(version) if isinstance(version, int) else None
        root = root or BAND_ROOT[_base(vector["base"])["version"]]
        for row in _rows_left(document, {"kind": "ref", "name": root}, root, ""):
            assert vector["python"] == "refuse", (row, vector["id"])
            if vector["web"] == "loaded" and row[1] == "*":
                # E-OD03: a retired home whose place a slice fills is dropped with no row.
                assert all(edit["path"][1:] in FLAT_DRAFT_TOP_LEVEL_FIELDS for edit in vector["edits"]), (
                    row,
                    vector["id"],
                )
                continue
            if vector["web"] == "loaded":
                # Rule 1 (E-OD06): a value that holds nothing at an optional field is dropped with no row.
                assert row[1] in ADMISSION_TYPES[row[0]]["optional"] and all(
                    edit["op"] == "set" and _holds_nothing(edit["value"]) for edit in vector["edits"]
                ), (row, vector["id"])
                continue
            left.setdefault(row, []).append(vector["id"])
    rows = _description_rows()
    assert set(left) <= rows
    assert {row for row in rows - set(left) if row[0] not in ROWS_WITHOUT_A_VECTOR} == set()
    assert {row[0] for row in rows - set(left)} == ROWS_WITHOUT_A_VECTOR


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
