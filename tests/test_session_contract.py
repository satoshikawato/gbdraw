"""`describe_session_document()` is true for the Sessions the writers write."""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Iterator

import pytest

from gbdraw import session_contract, session_drawings, session_io
from gbdraw.contract_shapes import shape_violations, typed_dict_shape, violations
from gbdraw.exceptions import ValidationError
from gbdraw.session import build_session_document, upgrade_session_document
from gbdraw.session_contract import describe_session_admission, describe_session_document, keys_of
from gbdraw.session_io import CURRENT_SESSION_VERSION, SESSION_FORMAT, validate_session
from gbdraw.session_request_codec import describe_canonical_request
from gbdraw.web_support.mode_scoped_settings import DIAGRAM_MODES, MODE_SCOPED_SETTINGS
from tests.test_canonical_request_shape import _request
from tests.test_session_admission_vectors import VECTORS, _edited
from tests.utils.contract_shape_check import REPO_ROOT, load_json, stored_json_files
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
        "kind": "record", "values": {"kind": "ref", "name": "FeatureOverrideDraftRow"}
    }


def test_a_typed_dict_is_described_with_its_optional_keys() -> None:
    shape = typed_dict_shape(session_contract.SessionCliInvocation)

    assert shape["optional"] == ["sessionExportable"] and not shape["open"]
    assert shape["fields"]["args"] == {"kind": "array", "items": {"kind": "string"}}
    assert typed_dict_shape(session_contract.SessionWebFiles, open=True)["open"]


# --- The admission description: what validate_session proves per version band ---

ADMISSION = describe_session_admission()
ADMISSION_TYPES: dict[str, Any] = ADMISSION["types"]
BAND_ROOT = {version: band["root"] for band in ADMISSION["bands"] for version in band["versions"]}
# Supported versions with no stored Session that Python accepts: no tracked file
# was written at 28; the tracked files at 27 are test inputs of 15 MB and more;
# the Gallery Sessions of `main` at 31 and 32 are refused by the Web's Load
# (its config override check reads their stored config unmigrated).
VERSIONS_WITHOUT_A_BASE = {27, 28, 31, 32}
# One value of each JSON kind, and an array and an object that hold a value.
_CANDIDATES: tuple[object, ...] = (None, True, 1, 0.5, "x", [], {}, [1], {"keyNoWriterWrites": 1})
_ABSENT = object()


def _admission_problems(session: dict[str, Any], root: str | None = None) -> list[str]:
    """Where ``session`` leaves the band of its version, or of the band ``root``."""
    shape = ADMISSION_TYPES[root or BAND_ROOT[session["version"]]]
    return violations(session, shape, ADMISSION_TYPES, "")


def _admission_bases() -> list[tuple[str, dict[str, Any]]]:
    """Per supported version the smallest stored Session that Python accepts, a
    Session with a Result set of each mode and drafts, and every vector Python accepts."""

    smallest: dict[int, tuple[int, str, dict[str, Any]]] = {}
    for path in stored_json_files():
        if path.stat().st_size > 1_000_000:
            continue  # a base is small; this keeps the Gallery's large Sessions out
        session = load_json(path)
        if not isinstance(session, dict) or session.get("format") != SESSION_FORMAT:
            continue
        try:
            validate_session(session)
        except ValidationError:
            continue
        size = path.stat().st_size
        if session["version"] not in smallest or size < smallest[session["version"]][0]:
            smallest[session["version"]] = (size, str(path.relative_to(REPO_ROOT)), session)
    two_mode = two_mode_session()
    two_mode["legacyArtifacts"] = {}
    two_mode["modes"]["circular"].setdefault("features", {})["featureOverrides"] = {
        '["record-1","feature-1"]': {
            "recordKey": "record-1", "biologicalFeatureId": "feature-1", "featureVisibility": "off",
            "labelVisibility": None, "labelText": None, "labelSourceText": None,
        }
    }
    return [
        *((name, session) for _size, name, session in smallest.values()),
        ("two-mode Session with drafts", two_mode),
        *((f"vector {vector['id']}", _edited(vector)) for vector in VECTORS if vector["python"] == "accept"),
    ]


def _leaves_shape(value: object, shape: dict[str, Any]) -> bool:
    return bool(shape_violations(value, shape, ADMISSION_TYPES, ""))


def _admission_mutations(
    value: object, shape: dict[str, Any], path: tuple[str | int, ...] = ()
) -> Iterator[tuple[tuple[str | int, ...], str | int, object]]:
    """Each edit that takes ``value`` out of ``shape``, found by walking the
    description: (path of the container, key, new value or ``_ABSENT``).
    A value of another kind or ``null``, an absent required field, and an
    unknown key of a closed object; two items of an array or a record."""

    kind = shape["kind"]
    if kind == "ref":
        yield from _admission_mutations(value, ADMISSION_TYPES[shape["name"]], path)
        return
    if kind == "union":
        option = next((option for option in shape["options"] if not _leaves_shape(value, option)), None)
        if option is not None:
            yield from _admission_mutations(value, option, path)
        return

    def field(key: str | int, item: object, item_shape: dict[str, Any]) -> Iterator[
        tuple[tuple[str | int, ...], str | int, object]
    ]:
        yield from ((path, key, candidate) for candidate in _CANDIDATES if _leaves_shape(candidate, item_shape))
        yield from _admission_mutations(item, item_shape, (*path, key))

    if kind == "array" and isinstance(value, list):
        for index, item in enumerate(value[:2]):
            yield from field(index, item, shape["items"])
    elif kind == "record" and isinstance(value, dict):
        for key, item in list(value.items())[:2]:
            yield from field(key, item, shape["values"])
    elif kind == "object" and isinstance(value, dict):
        for key, item_shape in shape["fields"].items():
            if key in value:
                yield from field(key, value[key], item_shape)
                if key not in shape["optional"]:
                    yield path, key, _ABSENT
            else:
                yield from ((path, key, candidate) for candidate in _CANDIDATES if _leaves_shape(candidate, item_shape))
        if not shape["open"]:
            yield path, "keyNoWriterWrites", 1


def test_every_supported_version_has_an_admission_band_and_a_base() -> None:
    assert sorted(BAND_ROOT) == sorted(session_io.SUPPORTED_SESSION_VERSIONS)
    assert [band["minVersion"] for band in ADMISSION["bands"]] == [
        min(session_io.SUPPORTED_SESSION_VERSIONS), session_io.CANONICAL_SESSION_MIN_VERSION,
        session_io.CURRENT_ARTIFACT_SESSION_MIN_VERSION, session_io.CURRENT_AUTHORITY_SESSION_MIN_VERSION,
        session_io.MODE_SCOPED_SESSION_MIN_VERSION,
    ]
    based = {session["version"] for _name, session in _admission_bases()}
    assert session_io.SUPPORTED_SESSION_VERSIONS - based <= VERSIONS_WITHOUT_A_BASE
    assert [band["root"] for band in ADMISSION["bands"] if not based & set(band["versions"])] == []


def test_stored_sessions_python_accepts_are_inside_their_band() -> None:
    problems = {}
    for name, session in _stored_sessions():
        try:
            validate_session(session)
        except ValidationError:
            continue
        if found := _admission_problems(session):
            problems[name] = found
    assert problems == {}


def test_a_session_outside_its_band_is_one_python_refuses() -> None:
    """Soundness (PLAN 9.3): for every base and every edit the description
    makes, a violation implies that validate_session refuses the document."""

    accepted: list[str] = []
    crashed: list[str] = []
    edits = 0
    for name, session in _admission_bases():
        validate_session(session)
        assert _admission_problems(session) == [], name
        band = BAND_ROOT[session["version"]]
        for path, key, value in list(_admission_mutations(session, {"kind": "ref", "name": band})):
            container: Any = session
            for step in path:
                container = container[step]
            before = container[key] if (isinstance(container, list) or key in container) else _ABSENT
            if value is _ABSENT:
                del container[key]
            else:
                container[key] = deepcopy(value)
            where = f"{name}: {'/'.join(map(str, (*path, key)))} = {'(absent)' if value is _ABSENT else value!r}"
            try:
                assert _admission_problems(session, band), where
                edits += 1
                validate_session(session)
                accepted.append(where)
            except ValidationError:
                pass
            except (AttributeError, TypeError, KeyError, ValueError) as error:
                crashed.append(f"{where}: {type(error).__name__}")
            finally:
                if before is _ABSENT:
                    del container[key]
                else:
                    container[key] = before
    assert accepted == []
    assert crashed == []
    assert edits > 500


def test_each_default_of_the_admission_description_fits_its_field() -> None:
    defaulted = []
    for name, shape in ADMISSION_TYPES.items():
        for field, value in shape.get("defaults", {}).items():
            defaulted.append(f"{name}.{field}")
            assert not shape_violations(value, shape["fields"][field], ADMISSION_TYPES, field), field
    assert sorted(defaulted) == sorted(
        [f"SessionChecked{minimum}.resources" for minimum in (31, 39, 40, 46)]
        + [f"SessionChecked{minimum}.{field}" for minimum in (40, 46) for field in ("results", "editorState")]
    )
