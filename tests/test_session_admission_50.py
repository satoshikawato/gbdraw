"""Python admits Session 50 (DESIGN §1) before the Web does (D-15-81).

``validate_session`` checks the closed project and its drawings, then checks
each drawing's one-drawing view as it checks a Session 46-49. The v49
fixtures were written by the 49 writer (``write_v49_fixtures.py``, named in
each ``.provenance.json``).
"""

from __future__ import annotations

import copy
import gzip
import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from gbdraw.exceptions import ValidationError
from gbdraw.session import (
    _build_session_document_from_drawings,
    _DrawingBuild,
    load_session_document,
    upgrade_session_document,
)
from gbdraw.session_drawings import read_session_drawings, write_session_drawings
from gbdraw.session_io import (
    CURRENT_SESSION_VERSION,
    DRAWINGS_SESSION_MIN_VERSION,
    SUPPORTED_SESSION_VERSIONS,
    validate_session,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "sessions"
V49 = sorted(FIXTURES.glob("drawings-*.v50.gbdraw-session.json.gz"))
Document = dict[str, Any]


def _load(name: str) -> Document:
    with gzip.open(FIXTURES / f"{name}.v50.gbdraw-session.json.gz", "rt", encoding="utf-8") as handle:
        return json.load(handle)


def test_python_reads_and_writes_49() -> None:
    assert CURRENT_SESSION_VERSION == DRAWINGS_SESSION_MIN_VERSION == max(SUPPORTED_SESSION_VERSIONS)
    assert len(V49) == 4


@pytest.mark.parametrize("path", V49, ids=lambda path: path.name.split(".")[0])
def test_each_v49_fixture_is_admitted_and_rewritten_unchanged(path: Path) -> None:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        document = json.load(handle)
    validate_session(document)
    loaded = load_session_document(path)
    drawings = read_session_drawings(document)
    active = next(index for index, drawing in enumerate(drawings) if drawing.id == document["ui"]["activeDrawingId"])

    rewritten = write_session_drawings(
        [drawing.artifacts.fields for drawing in drawings],
        ids=[drawing.id for drawing in drawings],
        names=[drawing.name for drawing in drawings],
        active=active,
    )

    assert {key: value for key, value in rewritten.items() if key != "createdAt"} == {
        key: value for key, value in document.items() if key != "createdAt"
    }
    assert [drawing.id for drawing in loaded.drawings] == [drawing["id"] for drawing in document["drawings"]]
    assert upgrade_session_document(loaded).document is loaded


def test_each_drawing_keeps_its_own_cards_and_losat_entries() -> None:
    cards = {drawing["id"]: drawing["webFiles"]["bindings"] for drawing in _load("drawings-per-drawing-cards")["drawings"]}
    assert cards["circular"]["d_color"] is not None
    assert cards["circular-2"]["d_color"] is None
    pair = {drawing["id"]: drawing["losatCache"]["entries"] for drawing in _load("drawings-linear-pair-losat")["drawings"]}
    assert (len(pair["linear"]), len(pair["linear-2"])) == (1, 0)


def _drawing(document: Document, drawing_id: str) -> Document:
    return next(drawing for drawing in document["drawings"] if drawing["id"] == drawing_id)


def _set(path: str, value: object) -> Callable[[Document], None]:
    def change(document: Document) -> None:
        *parents, key = path.split(".")
        target = document
        for part in parents:
            target = _drawing(document, part[1:]) if part.startswith("@") else target[part]
        target[key] = value

    return change


def _first_drawing_id(value: str) -> Callable[[Document], None]:
    def change(document: Document) -> None:
        document["drawings"][0]["id"] = value

    return change


def _duplicate_id(document: Document) -> None:
    document["drawings"][1]["id"] = document["drawings"][0]["id"]


def _unreferenced_resource(document: Document) -> None:
    document["resources"]["unused-gb"] = copy.deepcopy(next(iter(document["resources"].values())))


def _undrawn_with_results(document: Document) -> None:
    document["drawings"][0]["renderRequest"] = None


def _request_of_other_mode(document: Document) -> None:
    linear = _drawing(document, "linear")
    linear["renderRequest"] = copy.deepcopy(_drawing(document, "circular")["renderRequest"])


def _without_own_input(document: Document) -> None:
    del _drawing(document, "circular")["webFiles"]["bindings"]["c_gb"]


def _other_mode_binding(document: Document) -> None:
    _drawing(document, "circular-2")["webFiles"]["bindings"]["linearSeqs"] = []


@pytest.mark.parametrize(
    "fixture, change, message",
    [
        ("drawings-circular-linear", _set("modes", {}), "unclassified field.*modes"),
        ("drawings-circular-linear", _set("otherModeResult", {}), "unclassified field.*otherModeResult"),
        ("drawings-circular-linear", _set("renderRequest", None), "unclassified field.*renderRequest"),
        ("drawings-circular-linear", _set("results", []), "unclassified field.*results"),
        ("drawings-circular-linear", _set("config", {}), "unclassified field.*config"),
        ("drawings-circular-linear", _set("ui.mode", "linear"), "ui contains unclassified field.*mode"),
        ("drawings-circular-linear", _set("@linear.ui.mode", "linear"), "drawing linear ui"),
        ("drawings-circular-linear", _duplicate_id, "two drawings share an ID"),
        ("drawings-circular-linear", _set("ui.activeDrawingId", "linear-2"), "activeDrawingId"),
        ("drawings-circular-linear", _unreferenced_resource, "unused-gb"),
        ("drawings-circular-linear", _first_drawing_id("Circular"), "drawing ID"),
        ("drawings-circular-linear", _first_drawing_id("c" * 41), "drawing ID"),
        ("drawings-circular-linear", _set("drawings", []), "non-empty drawings"),
        ("drawings-circular-linear", _undrawn_with_results, "Results without a request"),
        ("drawings-circular-linear", _request_of_other_mode, "another mode"),
        ("drawings-circular-linear", _set("@linear.webFiles", {"circularInputOriginalName": "x"}), "webFiles"),
        ("drawings-per-drawing-cards", _other_mode_binding, "no other mode's input"),
        ("drawings-per-drawing-cards", _without_own_input, "requires c_gb"),
    ],
)
def test_python_refuses_a_49_document_outside_its_layout(
    fixture: str, change: Callable[[Document], None], message: str
) -> None:
    document = _load(fixture)
    validate_session(document)
    change(document)

    with pytest.raises(ValidationError, match=message):
        validate_session(document)


def test_the_conservation_web_files_are_the_circular_drawings() -> None:
    document = _load("drawings-per-drawing-cards")
    view = copy.deepcopy(dict(read_session_drawings(document)[0].artifacts.fields))
    resource_id = next(iter(view["resources"]))
    view["webFiles"]["conservationLosatFastaSources"] = [resource_id]
    view["webFiles"]["conservationBlastSource"] = "losat-cache"

    drawing = write_session_drawings([view])["drawings"][0]

    assert drawing["webFiles"]["conservationLosatFastaSources"] == [resource_id]
    assert drawing["webFiles"]["conservationBlastSource"] == "losat-cache"


def test_saving_a_49_again_drops_the_resources_no_drawing_names() -> None:
    base = load_session_document(FIXTURES / "drawings-circular-linear.v50.gbdraw-session.json.gz")
    linear = base.drawing("linear")
    state = {
        key: value
        for key, value in copy.deepcopy(dict(read_session_drawings(base.to_dict())[1].artifacts.fields)).items()
        if key not in {"format", "version", "createdAt", "renderRequest", "resources"}
    }
    state.update(results=[], runMetadata={}, editorState={**state["editorState"], "featureCatalog": None})

    saved = _build_session_document_from_drawings(
        (_DrawingBuild(mode="linear", request=None, state=state, id=linear.id, name=linear.name),),
        base=base,
        resources=base.to_dict()["resources"],
    ).to_dict()

    assert saved["version"] == CURRENT_SESSION_VERSION
    assert sorted(saved["resources"]) == ["colors-default-colors", "record-1-genbank"]
    assert _drawing(saved, "linear")["renderRequest"] is None
    assert _drawing(saved, "circular") == _drawing(base.to_dict(), "circular")
    assert saved["ui"]["activeDrawingId"] == "linear"
