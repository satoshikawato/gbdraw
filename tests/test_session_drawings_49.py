"""The Session 49 drawing layout, tested on documents directly.

``read_session_drawings`` sees each drawing of a Session 49 as the one-drawing
view that Sessions 46-48 give, and ``write_session_drawings(version=49)`` is its
inverse (DESIGN §1, D-15-79). These tests build and compare documents; the
admission of 49 is tested in ``tests/test_session_admission_49.py``.
"""

from __future__ import annotations

import copy
import gzip
import json
from pathlib import Path
from typing import Any

import pytest

from gbdraw.exceptions import ValidationError
from gbdraw.session_drawings import (
    active_session_drawing,
    read_session_drawings,
    replace_session_drawings,
    write_session_drawings,
)
from tests.utils.two_mode_session import two_mode_session

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "sessions"
CARD_FIXTURE = FIXTURES / "file-card-imports-circular.v48.gbdraw-session.json.gz"
CARDS = {"d_color": None, "t_color": None, "blacklist": None, "whitelist": None, "qualifier_priority": None}


def _linear_drawing(drawing_id: str, name: str, resource_id: str, entry: str) -> dict[str, Any]:
    source = {"kind": "genbank", "resourceId": resource_id}
    return {
        "id": drawing_id,
        "name": name,
        "mode": "linear",
        "renderRequest": {"schema": 9, "mode": "linear", "records": [{"source": source}]},
        "results": [{"svg": f"<svg>{drawing_id}</svg>"}],
        "editorState": {
            "featureCatalog": None,
            "legend": {"originalOrder": ["CDS"], "entries": [{"caption": "CDS"}]},
            "featureStrokes": {"f1": {"width": 2}},
        },
        "config": {"form": {"title": name}},
        "features": {"featureOverrides": {}},
        "ui": {"inputType": "gb", "selectedResultIndex": 0, "canvasPadding": {"top": 1}, "zoom": 2},
        "webFiles": {
            "bindings": {"schema": 2, "linearSeqs": [{"gb": resource_id}], "linearComparisons": [], **CARDS},
            "linearRecordMetadata": [{"definition": name}],
        },
        "losatCache": {"entries": [{"key": entry}]},
        "losatDerivedCache": {"entries": []},
    }


def _minimal_49() -> dict[str, Any]:
    """Two Linear drawings and a settings-only Circular one; the second Linear is shown."""

    return {
        "format": "gbdraw-session",
        "version": 49,
        "createdAt": "2026-10-23T00:00:00+00:00",
        "title": "three drawings",
        "resources": {"lambda-gb": {"name": "lambda.gb"}, "phage-gb": {"name": "phage.gb"}},
        "ui": {"activeDrawingId": "linear-2", "downloadDpi": 300, "featurePanelTab": "features"},
        "legacyArtifacts": {},
        "webFiles": {"resourceOriginalNames": {"lambda-gb": "lambda.gb"}},
        "drawings": [
            _linear_drawing("linear", "Linear", "lambda-gb", "shared"),
            _linear_drawing("linear-2", "Phage map", "phage-gb", "own"),
            {
                "id": "circular",
                "name": "Circular",
                "mode": "circular",
                "renderRequest": None,
                "results": [],
                "config": {"form": {"title": "draft"}},
                "ui": {"inputType": "gb"},
                "webFiles": {
                    "bindings": {"schema": 2, "c_gb": None, "c_conservation_blasts": [], **CARDS},
                    "circularOutputPrefixExplicit": False,
                },
            },
        ],
    }


def _identity(document: dict[str, Any]) -> list[tuple[str, str, str]]:
    return [(part.id, part.name, part.mode) for part in read_session_drawings(document)]


def _rewritten(document: dict[str, Any]) -> dict[str, Any]:
    parts = read_session_drawings(document)
    active = active_session_drawing(document, parts)
    assert active is not None
    return write_session_drawings(
        [part.artifacts.fields for part in parts],
        version=49,
        ids=[part.id for part in parts],
        names=[part.name for part in parts],
        active=parts.index(active),
    )


def test_a_49_document_reads_its_drawings_by_id() -> None:
    document = _minimal_49()
    parts = read_session_drawings(document)

    assert _identity(document) == [
        ("linear", "Linear", "linear"),
        ("linear-2", "Phage map", "linear"),
        ("circular", "Circular", "circular"),
    ]
    assert active_session_drawing(document, parts) is parts[1]
    assert parts[1].resource_ids == {"phage-gb"}
    assert parts[2].request is None and parts[2].results == ()
    second = parts[1].artifacts.fields
    assert second["ui"]["mode"] == "linear"
    assert second["ui"]["lInputType"] == "gb"
    assert second["ui"]["downloadDpi"] == 300
    assert second["modes"]["linear"]["ui"] == {"canvasPadding": {"top": 1}}
    assert second["modes"]["linear"]["editorState"] == {
        "featureStrokes": {"f1": {"width": 2}},
        "legend": {"entries": [{"caption": "CDS"}]},
    }
    assert second["editorState"] == {"featureCatalog": None, "legend": {"originalOrder": ["CDS"]}}
    assert second["webFiles"]["resourceOriginalNames"] == {"lambda-gb": "lambda.gb"}
    assert second["losatCache"] == {"entries": [{"key": "own"}]}


def test_a_49_document_rewritten_from_its_drawings_is_unchanged() -> None:
    document = _minimal_49()

    assert _rewritten(document) == document


def test_written_drawings_read_back_with_their_bindings_and_losat_entries() -> None:
    session = two_mode_session()
    session["webFiles"]["bindings"] = {
        "schema": 2,
        "c_gb": {"resourceId": "record-1-genbank", "name": "HmmtDNA.gbk"},
        "c_conservation_blasts": [],
        "linearSeqs": [{"gb": {"resourceId": "linear-record-1-genbank", "name": "NC_001416.gb"}}],
        "linearComparisons": [],
        **CARDS,
    }
    circular, linear = (part.artifacts.fields for part in read_session_drawings(session))
    linear = {**linear, "losatCache": {"entries": [{"key": "linear"}]}}
    second = {**copy.deepcopy(linear), "losatCache": {"entries": [{"key": "linear"}, {"key": "second"}]}}
    second["webFiles"] = {**second.get("webFiles", {}), "linearRecordMetadata": [{"definition": "second"}]}

    document = write_session_drawings(
        [circular, linear, second],
        version=49,
        ids=["circular", "linear", "linear-2"],
        names=["Circular", "Linear", "Linear 2"],
        active=2,
    )

    assert set(document) <= {"format", "version", "createdAt", "title", "resources", "ui", "legacyArtifacts",
                             "webFiles", "drawings"}
    assert document["version"] == 49
    assert "mode" not in document["ui"]
    assert set(document.get("webFiles", {})) <= {"resourceOriginalNames"}
    parts = read_session_drawings(document)
    assert _identity(document) == [
        ("circular", "Circular", "circular"),
        ("linear", "Linear", "linear"),
        ("linear-2", "Linear 2", "linear"),
    ]
    assert active_session_drawing(document, parts) is parts[2]
    circular_bindings = parts[0].artifacts.fields["webFiles"]["bindings"]
    assert "c_gb" in circular_bindings and "linearSeqs" not in circular_bindings
    assert circular_bindings == {
        key: value for key, value in circular["webFiles"]["bindings"].items()
        if key not in {"linearSeqs", "linearComparisons"}
    }
    for part in parts[1:]:
        assert part.artifacts.fields["webFiles"]["bindings"]["linearSeqs"] == circular["webFiles"]["bindings"]["linearSeqs"]
        assert "c_gb" not in part.artifacts.fields["webFiles"]["bindings"]
    assert parts[1].artifacts.fields["losatCache"] == {"entries": [{"key": "linear"}]}
    assert parts[2].artifacts.fields["losatCache"] == {"entries": [{"key": "linear"}, {"key": "second"}]}
    assert parts[2].artifacts.fields["webFiles"]["linearRecordMetadata"] == [{"definition": "second"}]
    for written, part in zip((circular, linear, second), parts):
        assert part.request == written["renderRequest"]
        assert list(part.results) == written["results"]
        assert part.artifacts.fields["modes"][part.mode] == written.get("modes", {}).get(part.mode, {})
    assert _rewritten(document) == document


def test_replacing_a_49_drawing_keeps_the_others_in_place() -> None:
    document = _minimal_49()
    replaced = copy.deepcopy(dict(read_session_drawings(document)[0].artifacts.fields))
    replaced["results"] = [{"svg": "<svg>new</svg>"}]
    replaced["title"] = "renamed"

    rewritten = replace_session_drawings(document, {"linear": replaced})

    assert rewritten["title"] == "renamed"
    assert rewritten["ui"]["activeDrawingId"] == "linear-2"
    assert rewritten["drawings"][0]["results"] == [{"svg": "<svg>new</svg>"}]
    assert rewritten["drawings"][1:] == document["drawings"][1:]


@pytest.mark.parametrize(
    "change, message",
    [
        (lambda view: view["webFiles"].update(unplacedWebFile=[]), "webFiles.unplacedWebFile"),
        (lambda view: view["resources"].update({"unused-gb": {"name": "unused.gb"}}), "unused-gb"),
        (lambda view: view.update(config={}), "config"),
    ],
)
def test_a_49_writer_refuses_what_it_cannot_place(change: Any, message: str) -> None:
    view = copy.deepcopy(dict(read_session_drawings(_minimal_49())[0].artifacts.fields))
    view["resources"] = {"lambda-gb": {"name": "lambda.gb"}}
    assert write_session_drawings([view], version=49)["drawings"][0]["id"] == "linear"
    change(view)

    with pytest.raises(ValidationError, match=message):
        write_session_drawings([view], version=49)


def test_49_drawing_names_are_unique_ignoring_case() -> None:
    view = read_session_drawings(_minimal_49())[0].artifacts.fields

    with pytest.raises(ValidationError, match="name of its own"):
        write_session_drawings([view, view], version=49, ids=["linear", "linear-2"], names=["Map", "map"])


def test_the_48_writer_is_unchanged_for_an_existing_48_session() -> None:
    with gzip.open(CARD_FIXTURE, "rt", encoding="utf-8") as handle:
        document = json.load(handle)
    assert document["version"] == 48

    views = [part.artifacts.fields for part in read_session_drawings(document)]

    assert write_session_drawings(views) == document
    with pytest.raises(ValidationError, match="takes no drawing IDs"):
        write_session_drawings(views, ids=["circular"], names=["Circular"])
