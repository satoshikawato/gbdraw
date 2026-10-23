"""Web Load of a Session 27-44: one Worker call returns its Session 49."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from gbdraw.session import load_session_document, materialize_session, render_session
from gbdraw.session_io import DRAWINGS_SESSION_MIN_VERSION, SessionDraftMigration
from gbdraw.session_migration import session_draft_migration_repairs
from gbdraw.web_support.legacy_session import legacy_session_to_drawings_json

SESSIONS = Path(__file__).parent / "fixtures" / "sessions"


def _text(name: str) -> str:
    return gzip.decompress((SESSIONS / name).read_bytes()).decode("utf-8")


@pytest.mark.parametrize(
    ("name", "mode"),
    [
        ("scale-interval-zero-circular-cli.v44.gbdraw-session.json.gz", "circular"),
        # Session 42: the stored comparison table of a reversed record moves to the search frame.
        ("q-frame-main-linear-reverse.v42.gbdraw-session.json.gz", "linear"),
    ],
)
def test_a_legacy_session_becomes_session_49_with_its_results_and_replays_as_before(
    name: str, mode: str, tmp_path: Path
) -> None:
    source = json.loads(_text(name))
    payload = json.loads(legacy_session_to_drawings_json(_text(name)))
    assert set(payload) == {"document", "storedResults", "warnings"}
    assert payload["warnings"] == []
    # A Session 40-44 keeps its Results with their catalog in the document.
    assert payload["storedResults"] == {}
    document = payload["document"]
    assert document["version"] == DRAWINGS_SESSION_MIN_VERSION
    first = document["drawings"][0]
    assert first["mode"] == mode
    assert document["ui"]["activeDrawingId"] == first["id"]
    assert first["results"] == source["results"]
    # The catalog is promoted to the current schema, one item per kept Result.
    assert len(first["editorState"]["featureCatalog"]["items"]) == len(source["results"])
    with materialize_session(load_session_document(source), output_directory=tmp_path / "old") as materialized:
        old = render_session(materialized)
    with materialize_session(load_session_document(document), output_directory=tmp_path / "new") as materialized:
        new = render_session(materialized, drawing=first["id"])
    assert len(old.output_paths) == len(new.output_paths) == 1
    replayed = new.output_paths[0].read_text(encoding="utf-8")
    assert "<svg" in replayed
    assert replayed == old.output_paths[0].read_text(encoding="utf-8")


@pytest.mark.parametrize(
    ("name", "mode"),
    [
        # Results without a request (27-30) or a feature catalog (31-39): a Session 49 cannot hold
        # them, so they travel beside the document, keyed by their drawing (OV-1659).
        ("scale-interval-zero-circular-cli.v30.gbdraw-session.json.gz", "circular"),
        ("mode-switch-circular-result-linear-shown.v30.gbdraw-session.json.gz", "circular"),
        ("conservation-fasta.v39.gbdraw-session.json.gz", "circular"),
    ],
)
def test_stored_results_session_49_cannot_hold_travel_beside_the_document(name: str, mode: str) -> None:
    source = json.loads(_text(name))
    payload = json.loads(legacy_session_to_drawings_json(_text(name)))
    assert set(payload) == {"document", "storedResults", "warnings"}, payload.get("error")
    document = payload["document"]
    assert load_session_document(document).version == DRAWINGS_SESSION_MIN_VERSION
    (drawing_id, stored), = payload["storedResults"].items()
    (drawing,) = (drawing for drawing in document["drawings"] if drawing["id"] == drawing_id)
    assert drawing["mode"] == mode
    assert drawing["results"] == []
    assert drawing["editorState"]["featureCatalog"] is None
    if source["version"] < 31:
        # A Session 27-30 has no request: an undrawn drawing that keeps its bindings.
        assert drawing["renderRequest"] is None
        assert drawing["webFiles"]["bindings"]["c_gb"]["resourceId"] in document["resources"]
    # Every stored Result, byte for byte.
    assert json.dumps(stored, sort_keys=True) == json.dumps(source["results"], sort_keys=True)


def test_a_document_that_is_not_json_is_refused_in_the_error_shape() -> None:
    payload = json.loads(legacy_session_to_drawings_json("{not json"))
    assert payload["error"]["code"] == "HELPER_PROTOCOL"


def test_a_session_49_is_not_moved_again() -> None:
    document = json.loads(legacy_session_to_drawings_json(_text("scale-interval-zero-circular-cli.v44.gbdraw-session.json.gz")))
    payload = json.loads(legacy_session_to_drawings_json(json.dumps(document["document"])))
    assert set(payload) == {"error"}


def test_draft_migration_repairs_are_the_web_load_repairs() -> None:
    migration = SessionDraftMigration(
        session={}, dropped_feature_edit_count=2, narrowed_visibility_count=1, unmapped_hash_selector_count=3,
    )
    assert session_draft_migration_repairs(migration) == (
        {"code": "feature-edit-unmatched", "path": "features", "effect": "dropped", "count": 2},
        {"code": "feature-visibility-narrowed", "path": "features.featureVisibilityOverrides",
         "effect": "changed", "count": 1},
        {"code": "hash-selector-unmapped", "path": "config.rules", "effect": "changed", "count": 3},
    )
    assert session_draft_migration_repairs(SessionDraftMigration(session={})) == ()
