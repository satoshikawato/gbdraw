"""Web Load of a Session 27-44: one Worker call returns its Session 49."""

from __future__ import annotations

import base64
import gzip
import json
from pathlib import Path
from typing import Any

import pytest

from gbdraw.exceptions import GbdrawError

from gbdraw.session import load_session_document, materialize_session, render_session
from gbdraw import session_drawings
from gbdraw.session_io import (
    CLI_SESSION_ADMISSION,
    DRAWINGS_SESSION_MIN_VERSION,
    WEB_READ_SESSION_ADMISSION,
    SessionAdmission,
    SessionDraftMigration,
)
from gbdraw.session_migration import session_draft_migration_repairs
from gbdraw.web_support.feature_catalog import promote_legacy_feature_catalog
from gbdraw.web_support.legacy_session import legacy_session_to_drawings_json

SESSIONS = Path(__file__).parent / "fixtures" / "sessions"
REPLY_FIELDS = {"document", "storedResults", "resultSetSources", "warnings"}


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
    assert set(payload) == REPLY_FIELDS
    assert payload["warnings"] == []
    # A Session 40-44 keeps its Results with their catalog in the document.
    assert payload["storedResults"] == {}
    document = payload["document"]
    assert document["version"] == DRAWINGS_SESSION_MIN_VERSION
    first = document["drawings"][0]
    assert first["mode"] == mode
    assert document["ui"]["activeDrawingId"] == first["id"]
    # The file's one Result set (top level) is the drawing's that holds its request.
    assert payload["resultSetSources"] == {first["id"]: ""}
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
    assert set(payload) == REPLY_FIELDS, payload.get("error")
    document = payload["document"]
    assert load_session_document(document).version == DRAWINGS_SESSION_MIN_VERSION
    (drawing_id, stored), = payload["storedResults"].items()
    assert payload["resultSetSources"] == {drawing_id: ""}
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


@pytest.mark.parametrize(
    ("name", "result_mode", "shown_mode"),
    [
        ("mode-switch-circular-result-linear-shown.v30.gbdraw-session.json.gz", "circular", "linear"),
        ("mode-switch-linear-result-circular-shown.v30.gbdraw-session.json.gz", "linear", "circular"),
    ],
)
@pytest.mark.parametrize("settings_used", [True, False], ids=["settings-used", "settings-at-default"])
def test_a_legacy_session_keeps_its_shown_mode_and_opens_its_stored_result(
    name: str, result_mode: str, shown_mode: str, settings_used: bool, monkeypatch: pytest.MonkeyPatch
) -> None:
    # The mode shown at Save gets a drawing; with no Result of its own, the
    # drawing with the stored Result opens and comes first.
    # Being shown is a use on its own: rule d (other-mode settings) need not
    # hold for it (OV-1666).
    if not settings_used:
        monkeypatch.setattr(session_drawings, "_used_by_settings", lambda *_args: False)
    payload = json.loads(legacy_session_to_drawings_json(_text(name)))
    assert set(payload) == REPLY_FIELDS, payload.get("error")
    document = payload["document"]
    # The shown mode's drawing took no Result set of the file.
    assert payload["resultSetSources"] == {result_mode: ""}
    assert [(drawing["id"], drawing["mode"]) for drawing in document["drawings"]] == [
        (result_mode, result_mode),
        (shown_mode, shown_mode),
    ]
    assert document["ui"]["activeDrawingId"] == result_mode
    assert set(payload["storedResults"]) == {result_mode}
    shown = document["drawings"][1]
    assert shown["renderRequest"] is None
    assert shown["results"] == []


def _web_read_44_without_catalog() -> dict[str, Any]:
    # Web Load's readers dropped the catalog of a Session 44's Result (vector a7b).
    document = json.loads(_text("HmmtDNA_basic_circular.v44-schema8.gbdraw-session.json.gz"))
    document["editorState"] = {"featureCatalog": None}
    return document


def _web_read_39_with_unresolved_protein_cache() -> dict[str, Any]:
    # Web Load's LOSAT reader kept the first of two entries with one key (vector a13);
    # Python's manifest does not resolve it under that key.
    document = json.loads(_text("BGC0000708-BGC0000713.v39.gbdraw-session.json.gz"))
    entries = document["losatCache"]["entries"]
    entries[0]["key"] = "DUP"
    del entries[1]
    return document


def test_the_web_entry_does_not_check_again_what_the_web_readers_read() -> None:
    # D-15-131: the CLI refuses these Sessions; the Web entry converts what the
    # Web's readers kept, here and in the check of the Session 49.
    for document, refusal in (
        (_web_read_44_without_catalog(), "Results require editorState.featureCatalog"),
        (_web_read_39_with_unresolved_protein_cache(), "does not resolve through the manifest"),
    ):
        with pytest.raises(GbdrawError, match=refusal):
            load_session_document(document)
        payload = json.loads(legacy_session_to_drawings_json(json.dumps(document)))
        assert set(payload) == REPLY_FIELDS, payload.get("error")
        assert payload["document"]["version"] == DRAWINGS_SESSION_MIN_VERSION
        (drawing_id, stored), = payload["storedResults"].items()
        assert payload["resultSetSources"] == {drawing_id: ""}
        assert json.dumps(stored, sort_keys=True) == json.dumps(document["results"], sort_keys=True)


def test_the_web_admission_skips_only_the_catalog_and_manifest_checks() -> None:
    assert WEB_READ_SESSION_ADMISSION == SessionAdmission(feature_catalog=False, losat_manifest=False)
    assert CLI_SESSION_ADMISSION == SessionAdmission(feature_catalog=True, losat_manifest=True)
    # Another refusal of the CLI's admission stays a refusal of the Web entry.
    document = _web_read_44_without_catalog()
    document["results"] = {"not": "a list"}
    payload = json.loads(legacy_session_to_drawings_json(json.dumps(document)))
    assert set(payload) == {"error"}


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


def test_a_session_27_30_inline_file_takes_its_decoded_size() -> None:
    # OV-1993 (D-15-144 A): the Web reader of a Session 27-30 ignored an inline
    # file's stored `size`; Python's projection of that `files` takes the
    # decoded bytes' size, so the file converts (OV-712's shape: bytes
    # replaced, `size` unchanged). The bytes become the drawing's resource.
    document = json.loads(_text("scale-interval-negative-linear-cli.v30.gbdraw-session.json.gz"))
    replaced = b"LOCUS broken\nORIGIN\n 1 zzzz\n"
    entry = document["files"]["linearSeqs"][0]["gb"]
    assert entry["size"] != len(replaced)
    entry["data"] = base64.b64encode(replaced).decode("ascii")
    payload = json.loads(legacy_session_to_drawings_json(json.dumps(document)))
    assert set(payload) == REPLY_FIELDS, payload.get("error")
    (drawing,) = payload["document"]["drawings"]
    (binding,) = drawing["webFiles"]["bindings"]["linearSeqs"]
    resource = payload["document"]["resources"][binding["gb"]["resourceId"]]
    assert (resource["size"], base64.b64decode(resource["data"])) == (len(replaced), replaced)


def test_a_session_40_legacy_alignment_reads_the_catalog_the_web_reader_promoted() -> None:
    # OV-1995: the Web reader sends Python a Session 40's feature catalog
    # promoted to the current schema; the legacy alignment target (a similarity
    # group of that catalog) is found as in the stored schema-3 catalog.
    document = json.loads(
        (SESSIONS / "BGC0000708-BGC0000713.v40-schema5.json").read_text(encoding="utf-8")
    )
    stored = json.loads(legacy_session_to_drawings_json(json.dumps(document)))
    editor_state = document["editorState"]
    editor_state["featureCatalog"] = promote_legacy_feature_catalog(editor_state["featureCatalog"])
    promoted = json.loads(legacy_session_to_drawings_json(json.dumps(document)))
    assert set(promoted) == REPLY_FIELDS, promoted.get("error")
    (drawing, _) = promoted["document"]["drawings"]
    alignment = drawing["renderRequest"]["layout"]["similarityAlignment"]
    assert alignment["groupId"] == "og_1"
    assert promoted["document"]["drawings"] == stored["document"]["drawings"]
    # Each anchor carries the stable feature ID the Web's legacy promotion
    # writes (the catalog's stableFeatureId, else the biological ID), so the
    # promoted request the Web sends agrees with Python's plan.
    (item,) = editor_state["featureCatalog"]["items"]
    stable = {
        (feature["recordKey"], feature["biologicalFeatureId"]): feature.get("stableFeatureId")
        or feature["biologicalFeatureId"]
        for feature in item["biologicalFeatures"]
    }
    anchors = [alignment["reference"], *(r["anchor"] for r in alignment["records"] if r["anchor"])]
    assert [a["stableFeatureSvgId"] for a in anchors] == [
        stable[(a["recordKey"], a["biologicalFeatureId"])] for a in anchors
    ]
