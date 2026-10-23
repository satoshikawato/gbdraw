"""The drawing vectors that Python and the Web share (B2-D4, DESIGN §4.2).

``tests/fixtures/sessions/drawings-vectors.json`` pins the drawing ID rule:
``gbdraw.session_drawings.allocate_drawing_id`` here and
``allocateDrawingId`` in ``gbdraw/web/js/services/drawing-view.js``
(``tests/web/drawing-ids.test.mjs``) read the same ``ids`` cases. Its
``upgrades`` cases pin the move of each Session 46-48 fixture to drawings
(``gbdraw.session_drawings.drawings_document_from_mode_scoped``); the Web's
``upgradeModeScopedSession`` reads them from B2d on.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from gbdraw.exceptions import ValidationError
from gbdraw.session import SessionDocument
from gbdraw.session_drawings import (
    allocate_drawing_id,
    allocate_session_drawing,
    drawings_document_from_mode_scoped,
)

REPO_ROOT = Path(__file__).resolve().parent.parent
FIXTURES = REPO_ROOT / "tests" / "fixtures" / "sessions"
VECTORS = json.loads((FIXTURES / "drawings-vectors.json").read_text(encoding="utf-8"))


def _load(path: Path) -> dict[str, object]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        data = json.load(handle)
    assert isinstance(data, dict)
    return data


def test_vectors_declare_their_layout() -> None:
    assert VECTORS["schema"] == 1
    assert set(VECTORS) == {"schema", "ids", "upgrades"}


@pytest.mark.parametrize("case", VECTORS["ids"], ids=lambda case: f"{case['mode']}+{'.'.join(case['taken'])}")
def test_drawing_id_follows_the_shared_vectors(case: dict[str, object]) -> None:
    taken = case["taken"]
    assert isinstance(taken, list)
    assert allocate_drawing_id(str(case["mode"]), taken) == case["id"]


def test_a_second_drawing_of_a_mode_is_numbered() -> None:
    assert allocate_session_drawing("linear", requested_id=None, requested_name=None, taken_ids=["linear"]) == (
        "linear-2",
        "Linear 2",
    )
    assert allocate_session_drawing(
        "linear", requested_id="linear-2", requested_name="Linear 2", taken_ids=["circular", "linear"]
    ) == ("linear-2", "Linear 2")


def test_upgrades_cover_every_session_46_to_48_fixture() -> None:
    fixtures = {
        str(path.relative_to(REPO_ROOT))
        for path in FIXTURES.glob("*.gbdraw-session.json.gz")
        if _load(path).get("version") in (46, 47, 48)
    }
    assert {case["input"] for case in VECTORS["upgrades"]} == fixtures


@pytest.mark.parametrize("case", VECTORS["upgrades"], ids=lambda case: str(case["name"]))
def test_session_46_to_48_moves_to_the_vector_drawings(case: dict[str, object]) -> None:
    document, facts = drawings_document_from_mode_scoped(_load(REPO_ROOT / str(case["input"])))
    SessionDocument(document)
    assert [
        {"id": drawing["id"], "mode": drawing["mode"], "keys": sorted(drawing)}
        for drawing in document["drawings"]
    ] == case["drawings"]
    assert document["ui"]["activeDrawingId"] == case["activeDrawingId"]
    assert facts == case["sourceFacts"]


def _settings_only_v48() -> dict[str, object]:
    return _load(FIXTURES / "file-card-imports-circular.v48.gbdraw-session.json.gz")


def _slices(data: dict[str, object]) -> dict[str, dict[str, object]]:
    modes = data["modes"]
    assert isinstance(modes, dict)
    return modes


def _ids(data: dict[str, object]) -> list[str]:
    document, _facts = drawings_document_from_mode_scoped(data)
    SessionDocument(document)
    return [drawing["id"] for drawing in document["drawings"]]


def test_other_mode_settings_count_only_when_they_differ_from_its_default_and_the_shown_mode() -> None:
    data = _settings_only_v48()
    linear = _slices(data)["linear"]
    assert _ids(data) == ["circular"]
    # A value of a row both modes hold, equal to the shown mode's: not a use.
    linear["config"]["colors"] = _slices(data)["circular"]["config"]["colors"]
    assert _ids(data) == ["circular"]
    # A Linear-only row that differs from its default is a use.
    linear["config"]["form"]["show_labels_linear"] = "all"
    assert _ids(data) == ["circular", "linear"]


def test_a_linear_record_layout_counts_only_when_it_arranges_records() -> None:
    data = _settings_only_v48()
    layout = _slices(data)["linear"]["config"]["linearRecordLayout"]
    # The rows the Web fills for its sequence rows, with their run-time uids, arrange nothing.
    layout["rows"] = [{"row": 1, "uid": "linear-seq-a"}, {"row": 2, "uid": "linear-seq-b"}]
    assert _ids(data) == ["circular"]
    # Two records in one row are an arrangement.
    layout["rows"][1]["row"] = 1
    assert _ids(data) == ["circular", "linear"]


def test_other_mode_keyed_rows_count_entry_by_entry() -> None:
    data = _settings_only_v48()
    circular, linear = _slices(data)["circular"], _slices(data)["linear"]
    edits = _slices(_load(FIXTURES / "rule-captions-with-color-id.v46.gbdraw-session.json.gz"))["circular"]
    (key, entry), *_ = edits["features"]["featureColorOverrides"].items()
    circular["features"]["featureColorOverrides"] = {key: entry}
    linear["features"]["featureColorOverrides"] = {key: entry}
    assert _ids(data) == ["circular"]
    linear["features"]["featureColorOverrides"] = {key: entry, f"{key}-2": entry}
    assert _ids(data) == ["circular", "linear"]


def test_other_mode_input_binding_counts_and_an_empty_sequence_row_does_not() -> None:
    data = _settings_only_v48()
    bindings = data["webFiles"]["bindings"]
    assert _ids(data) == ["circular"]
    # A Depth file alone counts.
    bindings["linearSeqs"][0]["depth"] = dict(bindings["d_color"])
    assert _ids(data) == ["circular", "linear"]


def test_the_shown_mode_comes_first_and_is_active() -> None:
    data = _load(FIXTURES / "legend-renamed-rows.v46.gbdraw-session.json.gz")
    data["ui"]["mode"] = "linear"
    document, _facts = drawings_document_from_mode_scoped(data)
    SessionDocument(document)
    assert [(drawing["id"], drawing["renderRequest"] is None) for drawing in document["drawings"]] == [
        ("linear", True),
        ("circular", False),
    ]
    assert document["ui"]["activeDrawingId"] == "linear"


def test_session_46_legend_edits_are_keyed_by_row_key() -> None:
    from gbdraw.session_io import session_legend_edits_by_row_key

    data = _load(FIXTURES / "legend-renamed-rows.v46.gbdraw-session.json.gz")
    document, _facts = drawings_document_from_mode_scoped(data)
    keyed, _outcomes = session_legend_edits_by_row_key(data)
    expected = keyed["modes"]["circular"]["editorState"]["legend"]
    assert expected != _slices(data)["circular"]["editorState"]["legend"]
    legend = document["drawings"][0]["editorState"]["legend"]
    assert {key: legend[key] for key in expected} == expected


def test_source_facts_name_a_session_48_without_a_web_draft() -> None:
    data = _settings_only_v48()
    del data["modes"]
    _document, facts = drawings_document_from_mode_scoped(data)
    assert facts == {"cardRowsAreFileRows": True, "proteinResultsInDrawnFrame": False}


@pytest.mark.parametrize("version", [44, 49])
def test_only_sessions_46_to_48_move_to_drawings(version: int) -> None:
    data = _settings_only_v48()
    data["version"] = version
    with pytest.raises(ValidationError, match="moves to drawings"):
        drawings_document_from_mode_scoped(data)


def test_a_circular_only_session_48_upgrades_to_one_drawing_with_the_packaged_defaults() -> None:
    from gbdraw.session import upgrade_session_document

    upgraded = upgrade_session_document(_settings_only_v48(), version=49).document.to_dict()
    assert [drawing["id"] for drawing in upgraded["drawings"]] == ["circular"]


def test_the_packaged_defaults_hold_every_registry_row_of_each_mode() -> None:
    from importlib import resources

    from gbdraw.web_support.mode_scoped_settings import MODE_SCOPED_SETTINGS

    packaged = json.loads(resources.files("gbdraw.data").joinpath("mode-slice-defaults.json").read_text("utf-8"))
    for mode in ("circular", "linear"):
        assert [(row["domain"], row["path"]) for row in packaged["defaults"][mode]] == [
            (row.domain, row.path) for row in MODE_SCOPED_SETTINGS if row.modes in ("both", mode)
        ]
