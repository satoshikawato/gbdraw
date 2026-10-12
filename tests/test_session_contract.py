"""`describe_session_document()` is true for the Sessions the writers write."""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Iterator

import pytest

from gbdraw import session_contract, session_drawings, session_io
from gbdraw.contract_shapes import typed_dict_shape
from gbdraw.exceptions import ValidationError
from gbdraw.session import build_session_document, upgrade_session_document
from gbdraw.session_contract import describe_session_document, keys_of
from gbdraw.session_io import CURRENT_SESSION_VERSION, SESSION_FORMAT, validate_session
from gbdraw.session_request_codec import describe_canonical_request
from gbdraw.web_support.mode_scoped_settings import DIAGRAM_MODES, MODE_SCOPED_SETTINGS
from tests.test_canonical_request_shape import _request
from tests.utils.contract_shape_check import REPO_ROOT, load_json, stored_json_files, violations
from tests.utils.two_mode_session import two_mode_session

DESCRIPTION = describe_session_document()
TYPES: dict[str, Any] = {**describe_canonical_request()["types"], **DESCRIPTION["types"]}
ROOT: dict[str, Any] = TYPES[DESCRIPTION["root"]]


def _problems(session: object) -> list[str]:
    return violations(session, ROOT, TYPES, "session")


def _stored_sessions() -> Iterator[tuple[str, dict[str, Any]]]:
    for path in stored_json_files():
        value = load_json(path)
        if isinstance(value, dict) and value.get("format") == SESSION_FORMAT:
            yield str(path.relative_to(REPO_ROOT)), value


def test_stored_sessions_of_the_current_version_match_the_description() -> None:
    """The Web wrote the Gallery Sessions; the CLI and the Python API wrote the documentation Sessions."""

    current = {
        name: _problems(session)
        for name, session in _stored_sessions()
        if session.get("version") == CURRENT_SESSION_VERSION
    }
    assert {name: found for name, found in current.items() if found} == {}
    assert len(current) >= 12


def test_sessions_that_python_upgrades_match_the_description() -> None:
    """An older Session saved again holds the migrated Web draft in `modes`."""

    upgraded = 0
    problems: dict[str, list[str]] = {}
    for name, session in _stored_sessions():
        if session.get("version") == CURRENT_SESSION_VERSION:
            continue
        try:
            document = upgrade_session_document(session).document
        except ValidationError:
            continue  # a fixture of a Session the readers refuse, or one that needs its sources
        upgraded += 1
        found = _problems(document.to_dict())
        if found:
            problems[name] = found
    assert problems == {}
    assert upgraded >= 20


@pytest.mark.parametrize("mode", DIAGRAM_MODES)
def test_a_document_built_from_a_request_matches_the_description(mode: str) -> None:
    session = build_session_document(_request(mode)).to_dict()

    assert _problems(session) == []
    assert session["renderRequest"]["mode"] == mode


def test_a_session_with_a_result_set_of_each_mode_matches_the_description() -> None:
    session = two_mode_session()
    validate_session(session)

    assert _problems(session) == []


def _closed_objects(value: object, shape: dict[str, Any], path: tuple[str | int, ...] = ()) -> Iterator[
    tuple[tuple[str | int, ...], str]
]:
    """Each object of a Session whose description takes no other key, outside the request."""

    kind = shape["kind"]
    if kind == "ref" and shape["name"] in DESCRIPTION["types"]:
        for found, _name in _closed_objects(value, TYPES[shape["name"]], path):
            yield found, _name or shape["name"]
    elif kind == "union":
        for option in shape["options"]:
            if not violations(value, option, TYPES, "session"):
                yield from _closed_objects(value, option, path)
                break
    elif kind == "array" and isinstance(value, list):
        for index, item in enumerate(value[:1]):
            yield from _closed_objects(item, shape["items"], (*path, index))
    elif kind == "record" and isinstance(value, dict):
        for key, item in list(value.items())[:1]:
            yield from _closed_objects(item, shape["values"], (*path, key))
    elif kind == "object" and isinstance(value, dict):
        if not shape["open"]:
            yield path, ""
        for key, item in value.items():
            if key in shape["fields"]:
                yield from _closed_objects(item, shape["fields"][key], (*path, key))


def test_an_unknown_key_in_each_closed_object_is_a_violation() -> None:
    edited = two_mode_session()
    edited["legacyArtifacts"] = {}
    edited["modes"]["circular"].setdefault("features", {})["featureOverrides"] = {
        '["record-1","feature-1"]': {
            "recordKey": "record-1", "biologicalFeatureId": "feature-1", "featureVisibility": "off",
            "labelVisibility": None, "labelText": None, "labelSourceText": None,
        }
    }
    validate_session(edited)
    sessions = [edited, *(
        session for _name, session in _stored_sessions() if session.get("version") == CURRENT_SESSION_VERSION
    )]
    reached: set[str] = set()
    for session in sessions:
        for path, name in _closed_objects(session, {"kind": "ref", "name": DESCRIPTION["root"]}):
            reached.add(name)
            changed = deepcopy(session)
            holder: Any = changed
            for step in path:
                holder = holder[step]
            holder["keyTheDescriptionDoesNotName"] = None
            where = ".".join(["session", *map(str, path)]).replace(".0", "[0]")
            assert _problems(changed) == [f"{where}.keyTheDescriptionDoesNotName: unknown field"], name
    closed = {
        name for name, shape in DESCRIPTION["types"].items()
        if shape["kind"] == "object" and not shape["open"] and name != "SessionInput"
    }
    assert closed - reached == set()


def test_an_open_object_takes_other_keys() -> None:
    session = two_mode_session()
    for field in ("ui", "editorState", "orthogroupState", "webFiles"):
        assert TYPES[ROOT["fields"][field]["name"]]["open"], field
        session[field] = {**session.get(field, {}), "keyOnlyTheWebReads": 1}

    assert _problems(session) == []
    validate_session(session)


def test_the_validators_read_their_field_sets_from_the_declaration() -> None:
    top_level = set(ROOT["fields"])

    assert session_io.CURRENT_SESSION_TOP_LEVEL_FIELDS == top_level | keys_of(session_contract.LegacySessionFields)
    assert set(TYPES["SessionInput"]["fields"]) == session_io.CURRENT_SESSION_TOP_LEVEL_FIELDS
    assert session_io.FLAT_DRAFT_TOP_LEVEL_FIELDS < keys_of(session_contract.LegacySessionFields)
    assert top_level.isdisjoint(keys_of(session_contract.LegacySessionFields))
    # One committed Result set is named the same way at the top level and in `otherModeResult`.
    assert session_io.OTHER_MODE_RESULT_FIELDS == set(TYPES["SessionOtherModeResult"]["fields"]) <= top_level
    assert session_drawings._SET_FIELDS == keys_of(session_contract.SessionResultSet) | {"renderRequest"}
    for part, whole in (("SessionOtherModeUi", "SessionUi"), ("SessionOtherModeEditorState", "SessionEditorState")):
        assert set(TYPES[part]["fields"]) <= set(TYPES[whole]["fields"])
        assert not TYPES[part]["open"] and TYPES[whole]["open"]


def test_a_mode_slice_is_described_from_the_registry() -> None:
    described = {
        f"{container}.{key}"
        for container, name in (
            ("config.form", "SessionModeSliceConfigForm"), ("config.adv", "SessionModeSliceConfigAdv"),
            ("config.losat", "SessionModeSliceConfigLosat"), ("config", "SessionModeSliceConfig"),
            ("features", "SessionModeSliceFeatures"), ("editorState", "SessionModeSliceEditorState"),
            ("editorState.legend", "SessionModeSliceEditorStateLegend"),
            ("editorState.featureStrokes", "SessionModeSliceEditorStateFeatureStrokes"),
            ("ui", "SessionModeSliceUi"),
        )
        for key in TYPES[name]["fields"]
        if f"{container}.{key}" not in {"config.form", "config.adv", "config.losat", "editorState.legend", "editorState.featureStrokes"}
    }

    assert described == {f"{row.domain}.{row.path}" for row in MODE_SCOPED_SETTINGS}
    assert set(TYPES["SessionModes"]["fields"]) == set(DIAGRAM_MODES)
    assert TYPES["SessionModeSliceFeatures"]["fields"]["featureOverrides"] == {
        "kind": "record", "values": {"kind": "ref", "name": "SessionFeatureOverrideDraftRow"}
    }


def test_a_typed_dict_is_described_with_its_optional_keys() -> None:
    shape = typed_dict_shape(session_contract.SessionCliInvocation)

    assert shape["optional"] == ["sessionExportable"] and not shape["open"]
    assert shape["fields"]["args"] == {"kind": "array", "items": {"kind": "string"}}
    assert typed_dict_shape(session_contract.SessionWebFiles, open=True)["open"]
