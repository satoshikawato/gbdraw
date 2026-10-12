"""`describe_canonical_request()` is true for what the codec writes and reads."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import fields
import gzip
import inspect
import json
from pathlib import Path
from typing import Any, Iterator

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.requests import (
    CircularDiagramRequest,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
)
from gbdraw.session_request_codec import (
    CANONICAL_REQUEST_SCHEMA,
    CanonicalRequestCodecError,
    CanonicalRequestEncodingError,
    _SHARED_OPTION_WRONG_MODE_DEFAULTS,
    _option_wire_key,
    decode_canonical_request,
    describe_canonical_request,
    encode_canonical_request,
)
from tests import test_session_request_codec as codec_tests

REPO_ROOT = Path(__file__).resolve().parents[1]
DESCRIPTION = describe_canonical_request()
TYPES: dict[str, Any] = DESCRIPTION["types"]
ROOT: dict[str, Any] = {"kind": "ref", "name": DESCRIPTION["root"]}

_SCALAR_CHECKS = {
    "string": lambda value: isinstance(value, str),
    "boolean": lambda value: isinstance(value, bool),
    "integer": lambda value: isinstance(value, int) and not isinstance(value, bool),
    "number": lambda value: isinstance(value, (int, float)) and not isinstance(value, bool),
    "null": lambda value: value is None,
}


def violations(value: object, shape: dict[str, Any], path: str = "renderRequest") -> list[str]:
    """Where ``value`` leaves ``shape``: an unknown field, a missing one, or a wrong value."""

    kind = shape["kind"]
    if kind == "unknown":
        return []
    if kind == "ref":
        return violations(value, TYPES[shape["name"]], path)
    if kind in _SCALAR_CHECKS:
        return [] if _SCALAR_CHECKS[kind](value) else [f"{path}: {value!r} is not {kind}"]
    if kind == "literal":
        allowed = any(value == item and type(value) is type(item) for item in shape["values"])
        return [] if allowed else [f"{path}: {value!r} is not one of {shape['values']!r}"]
    if kind == "union":
        found = [violations(value, option, path) for option in shape["options"]]
        return min(found, key=len) if all(found) else []
    if kind in {"array", "tuple"}:
        if not isinstance(value, list):
            return [f"{path}: {value!r} is not an array"]
        if kind == "tuple" and len(value) != len(shape["items"]):
            return [f"{path}: expected {len(shape['items'])} items"]
        items = shape["items"] if kind == "tuple" else [shape["items"]] * len(value)
        return [
            problem
            for index, (item, item_shape) in enumerate(zip(value, items))
            for problem in violations(item, item_shape, f"{path}[{index}]")
        ]
    if not isinstance(value, dict):
        return [f"{path}: {value!r} is not an object"]
    if kind == "record":
        return [
            problem
            for key, item in value.items()
            for problem in violations(item, shape["values"], f"{path}.{key}")
        ]
    assert kind == "object", kind
    fields = shape["fields"]
    problems = [f"{path}.{key}: unknown field" for key in value if key not in fields]
    problems += [
        f"{path}.{key}: required field is absent"
        for key in fields
        if key not in value and key not in shape["optional"]
    ]
    for key, item in value.items():
        if key in fields:
            problems += violations(item, fields[key], f"{path}.{key}")
    return problems


def _resolved(shape: dict[str, Any]) -> dict[str, Any]:
    return _resolved(TYPES[shape["name"]]) if shape["kind"] == "ref" else shape


def _children(shape: dict[str, Any]) -> list[dict[str, Any]]:
    kind = shape["kind"]
    if kind == "union":
        return shape["options"]
    if kind == "tuple":
        return shape["items"]
    if kind == "array":
        return [shape["items"]]
    return [shape["values"]] if kind == "record" else []


def _described_objects() -> list[tuple[str, dict[str, Any]]]:
    """Each object of the description with a name: its typedef, then the path inside it."""

    found: list[tuple[str, dict[str, Any]]] = []

    def walk(shape: dict[str, Any], label: str) -> None:
        if shape["kind"] == "object":
            tag = shape["fields"].get("kind", {}).get("values")
            label = f"{label}[{tag[0]}]" if tag and "." in label else label
            found.append((label, shape))
            for key, item in shape["fields"].items():
                walk(item, f"{label}.{key}")
        for item in _children(shape):
            walk(item, label)

    for name, shape in TYPES.items():
        walk(shape, name)
    return found


OBJECTS = _described_objects()
LABELS = {id(shape): label for label, shape in OBJECTS}


def _objects(
    value: object, shape: dict[str, Any], path: tuple[str | int, ...] = ()
) -> Iterator[tuple[tuple[str | int, ...], dict[str, Any]]]:
    """Each object of a valid ``value`` with the object shape it matches."""

    kind = shape["kind"]
    if kind == "ref":
        yield from _objects(value, TYPES[shape["name"]], path)
    elif kind == "union":
        match = next(option for option in shape["options"] if not violations(value, option))
        yield from _objects(value, match, path)
    elif kind in {"array", "tuple"}:
        assert isinstance(value, list)
        items = shape["items"] if kind == "tuple" else [shape["items"]] * len(value)
        for index, (item, item_shape) in enumerate(zip(value, items)):
            yield from _objects(item, item_shape, (*path, index))
    elif kind == "record":
        assert isinstance(value, dict)
        for key, item in value.items():
            yield from _objects(item, shape["values"], (*path, key))
    elif kind == "object":
        assert isinstance(value, dict)
        yield path, shape
        for key, item in value.items():
            yield from _objects(item, shape["fields"][key], (*path, key))


def _literal_values(shape: dict[str, Any]) -> list[object]:
    """The values a field names one by one; none when any string is allowed beside them."""

    options = [_resolved(item) for item in (_children(shape) if shape["kind"] == "union" else [shape])]
    if any(item["kind"] == "string" for item in options):
        return []
    return [value for item in options if item["kind"] == "literal" for value in item["values"]]


def _with(payload: dict[str, Any], path: tuple[str | int, ...], key: str, *value: object) -> dict[str, Any]:
    """``payload`` with ``key`` of the object at ``path`` set to ``value``, or removed."""

    result = deepcopy(payload)
    holder: Any = result
    for step in path:
        holder = holder[step]
    if value:
        holder[key] = value[0]
    else:
        del holder[key]
    return result


def _decodes(payload: dict[str, Any], resource_paths: dict[str, Path], output: Path) -> bool:
    try:
        decode_canonical_request(payload, resource_paths=resource_paths, output_directory=output)
    except CanonicalRequestCodecError:
        return False
    return True


Encoded = tuple[dict[str, Any], dict[str, Path]]


@pytest.fixture(scope="module")
def encoded_requests(tmp_path_factory: pytest.TempPathFactory) -> dict[str, list[Encoded]]:
    """What each test of the codec that takes only a temporary directory encodes.

    Those tests build the representative typed requests of both modes; the
    payload is copied as encoded, before the test changes it.
    """

    result: dict[str, list[Encoded]] = {}
    with pytest.MonkeyPatch.context() as patch:
        for name in _CODEC_TESTS:
            folder = tmp_path_factory.mktemp("codec")
            seen: list[Encoded] = result.setdefault(name, [])

            def encode(*args: Any, _seen: list[Encoded] = seen, _folder: Path = folder, **kwargs: Any):
                encoded = encode_canonical_request(*args, **kwargs)
                resources = codec_tests._materialize_resources(encoded, _folder / f"shape-{len(_seen)}")
                _seen.append((deepcopy(encoded.payload), resources))
                return encoded

            patch.setattr(codec_tests, "encode_canonical_request", encode)
            getattr(codec_tests, name)(folder)
    return result


_CODEC_TESTS = sorted(
    name
    for name, function in vars(codec_tests).items()
    if name.startswith("test_")
    and inspect.isfunction(function)
    and list(inspect.signature(function).parameters) == ["tmp_path"]
)


@pytest.mark.parametrize("name", _CODEC_TESTS)
def test_every_encoded_request_matches_the_description(
    name: str, encoded_requests: dict[str, list[Encoded]]
) -> None:
    if not encoded_requests[name]:
        pytest.skip("This test of the codec encodes no request.")
    assert [violations(payload, ROOT) for payload, _resources in encoded_requests[name]] == [
        [] for _item in encoded_requests[name]
    ]


def _record() -> SeqRecord:
    return SeqRecord(Seq("ATGC" * 10), id="record", annotations={"molecule_type": "DNA"})


_REQUEST_TYPES: dict[str, tuple[Any, Any]] = {
    "circular": (CircularDiagramRequest, CircularDiagramOptions),
    "linear": (LinearDiagramRequest, LinearDiagramOptions),
}


def _request(mode: str, **options: object) -> Any:
    request_type, options_type = _REQUEST_TYPES[mode]
    return request_type(
        records=(RecordInput(source=InMemoryRecordSource(_record())),),
        options=options_type(**options),
    )


def _smallest_request(mode: str, tmp_path: Path) -> Encoded:
    encoded = encode_canonical_request(_request(mode))
    return encoded.payload, codec_tests._materialize_resources(encoded, tmp_path / "resources")


# Where one change cannot show that the decoder and the description agree, and why.
_COERCED = "accepted", "the typed constructor reads null as a value (str, bool, a default); no writer writes null"
_TOGETHER = "refused", "a decision's status, rationale and anchor are checked together; none encoded is skipped"
_KNOWN_DIFFERENCES = {
    "CanonicalAnnotationSet.id = null": _COERCED,
    "CanonicalRegionAnnotation.id = null": _COERCED,
    "CanonicalRegionAnnotation.label = null": _COERCED,
    "CanonicalCoordinateSpan.wrapsOrigin = null": _COERCED,
    "CanonicalFeatureSelector.value = null": _COERCED,
    "CanonicalHatchStyle.cross = null": _COERCED,
    "CanonicalGeneratedProteinComparison.mode = null": _COERCED,
    "CanonicalCircularDiagramOptions.conservationReference = null": _COERCED,
    "CanonicalLinearDiagramOptions.pairwiseMatchStyle = null": _COERCED,
    "CanonicalAlignmentRecordDecision.anchor = null": _TOGETHER,
    "CanonicalAlignmentRecordDecision.status = 'skipped'": _TOGETHER,
    "CanonicalAlignmentRecordDecision.rationale = 'skipped_by_user'": _TOGETHER,
    "CanonicalAlignmentRecordDecision.rationale = 'skipped_no_candidate'": _TOGETHER,
    "CanonicalAlignmentRecordDecision.rationale = 'skipped_unmappable'": _TOGETHER,
    "CanonicalCollinearityResultComparison.valueKind = 'blocks'": (
        "refused", "the resource of the encoded comparison holds a result, not blocks"
    ),
    "CanonicalLinearTrackOptions.linearTrackSlots = null": (
        "refused", "the encoded request has an axis index, which requires the slots"
    ),
}


def test_the_decoder_agrees_with_the_description_on_every_object(
    encoded_requests: dict[str, list[Encoded]], tmp_path: Path
) -> None:
    """Remove each field, set it to null and to each named value, add an unknown one.

    A cross-field rule can refuse one request and accept another, so a change
    the description allows must decode in some request, and one it refuses in none.
    """

    accepted: dict[tuple[str, bool], set[bool]] = {}
    reached: set[str] = set()
    decodable = [
        item
        for found in encoded_requests.values()
        for item in found
        if not violations(item[0], ROOT) and _decodes(*item, tmp_path)
    ]
    for mode in ("circular", "linear"):
        decodable.append(_smallest_request(mode, tmp_path / mode))
    for payload, resources in decodable:
        for path, shape in _objects(payload, ROOT):
            label = LABELS[id(shape)]
            reached.add(label)
            holder: Any = payload
            for step in path:
                holder = holder[step]
            changes: list[tuple[str, bool, tuple[object, ...]]] = [
                (f"{label}: an unknown field", False, ("fieldTheDescriptionDoesNotName", None)),
            ]
            for key, field in shape["fields"].items():
                if key in holder:
                    changes.append((f"{label}.{key} removed", key in shape["optional"], (key,)))
                changes.append((f"{label}.{key} = null", violations(None, field) == [], (key, None)))
                changes += [
                    (f"{label}.{key} = {value!r}", True, (key, value))
                    for value in _literal_values(field)
                    if value is not None
                ]
            for change, allowed, arguments in changes:
                unchanged = len(arguments) == 2 and holder.get(arguments[0], 0) == arguments[1]
                accepted.setdefault((change, allowed), set()).add(
                    unchanged or _decodes(_with(payload, path, *arguments), resources, tmp_path)
                )
    differences = {
        change: "refused" if allowed else "accepted"
        for (change, allowed), found in accepted.items()
        if (True in found) != allowed
    }
    assert differences == {key: value for key, (value, _reason) in _KNOWN_DIFFERENCES.items()}
    assert {label for label, _shape in OBJECTS} - reached == set()


# The decoder reads these keys and no writer writes them: the depth inputs of
# older writers (the writer stores every depth input as depthTracks), and a
# comparison table path, which the encoder refuses (a file is a resource).
_DEPTH_INPUTS_OF_OLDER_WRITERS = {
    "depthFile", "depthFiles", "depthTable", "depthTables", "depthTrackColors", "depthTrackFiles",
    "depthTrackLabels", "depthTrackLargeTickIntervals", "depthTrackSmallTickIntervals",
    "depthTrackTables", "depthTrackTickFontSizes",
}
_DECODER_ONLY_OPTIONS = {
    "circular": ("conservation_table_file", _DEPTH_INPUTS_OF_OLDER_WRITERS),
    "linear": ("comparison_table_file", _DEPTH_INPUTS_OF_OLDER_WRITERS | {"depthTrackHeights"}),
}


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_the_description_names_every_option_key_a_writer_can_write(
    mode: str, tmp_path: Path
) -> None:
    table_path, older_depth_inputs = _DECODER_ONLY_OPTIONS[mode]
    described = TYPES[f"Canonical{mode.capitalize()}DiagramOptions"]["fields"]
    payload, resources = _smallest_request(mode, tmp_path)

    # A shared option of the other mode is a key only at its default.
    defaults = {
        _option_wire_key(name): value
        for name, value in _SHARED_OPTION_WRONG_MODE_DEFAULTS[mode].items()
    }

    def decoder_names(key: str) -> bool:
        for value in (None, defaults.get(key)):
            try:
                decode_canonical_request(
                    _with(payload, ("diagramOptions",), key, value),
                    resource_paths=resources,
                    output_directory=tmp_path / "output",
                )
            except CanonicalRequestCodecError as error:
                if "Unknown field" in str(error):
                    continue
            return True
        return False

    candidates = {_option_wire_key(item.name) for item in fields(_REQUEST_TYPES[mode][1])}
    decoder_keys = {key for key in candidates | set(defaults) | set(described) if decoder_names(key)}
    assert set(described) <= decoder_keys
    assert decoder_keys - set(described) == older_depth_inputs | {_option_wire_key(table_path)}
    with pytest.raises(CanonicalRequestEncodingError, match="comparison table"):
        encode_canonical_request(_request(mode, **{table_path: "table.tsv"}))


def test_every_described_field_is_held_by_a_request_or_named_by_the_decoder(
    encoded_requests: dict[str, list[Encoded]],
) -> None:
    """Removing a field from the description must make some test fail.

    A request that holds the field would then hold an unknown field; an option
    key that no request holds is compared with the decoder's keys above.
    """

    held: set[str] = set()
    requests = [payload for found in encoded_requests.values() for payload, _resources in found]
    requests += [request for path in _session_files() for request in _stored_requests(_load(path))]
    for payload in requests:
        if violations(payload, ROOT):
            continue
        for path, shape in _objects(payload, ROOT):
            holder: Any = payload
            for step in path:
                holder = holder[step]
            held.update(f"{LABELS[id(shape)]}.{key}" for key in holder)
    described = {f"{label}.{key}" for label, shape in OBJECTS for key in shape["fields"]}
    assert {name for name in described - held if "DiagramOptions." not in name} == set()


def _stored_requests(value: object) -> Iterator[dict[str, Any]]:
    if isinstance(value, dict):
        if value.get("schema") == CANONICAL_REQUEST_SCHEMA and "diagramOptions" in value:
            yield value
        for item in value.values():
            yield from _stored_requests(item)
    elif isinstance(value, list):
        for item in value:
            yield from _stored_requests(item)


def _session_files() -> list[Path]:
    return sorted(
        path
        for folder in (
            "tests/fixtures", "tests/web/fixtures", "tests/test_inputs", "docs/images",
            "gbdraw/web/gallery/sessions",
        )
        for path in (REPO_ROOT / folder).rglob("*.json*")
        if path.suffix in {".json", ".gz"}
    )


def _load(path: Path) -> object:
    raw = path.read_bytes()
    return json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)


def test_requests_stored_in_current_sessions_match_the_description() -> None:
    """The Web wrote the Gallery Sessions; the CLI wrote the documentation Sessions."""

    modes: list[str] = []
    problems: dict[str, list[str]] = {}
    for path in _session_files():
        for request in _stored_requests(_load(path)):
            modes.append(request["mode"])
            found = violations(request, ROOT)
            if found:
                problems[str(path.relative_to(REPO_ROOT))] = found
    assert problems == {}
    assert modes.count("circular") >= 4 and modes.count("linear") >= 4


@pytest.mark.parametrize(
    ("value", "shape", "expected"),
    [
        ({"a": 1}, {"kind": "object", "fields": {}, "optional": []}, ["renderRequest.a: unknown field"]),
        (
            {},
            {"kind": "object", "fields": {"a": {"kind": "integer"}}, "optional": []},
            ["renderRequest.a: required field is absent"],
        ),
        ({}, {"kind": "object", "fields": {"a": {"kind": "integer"}}, "optional": ["a"]}, []),
        (True, {"kind": "integer"}, ["renderRequest: True is not integer"]),
        (1, {"kind": "literal", "values": [True]}, ["renderRequest: 1 is not one of [True]"]),
        ([1, "a"], {"kind": "array", "items": {"kind": "integer"}}, ["renderRequest[1]: 'a' is not integer"]),
        (None, {"kind": "union", "options": [{"kind": "string"}, {"kind": "null"}]}, []),
        ({"k": 1}, {"kind": "record", "values": {"kind": "string"}}, ["renderRequest.k: 1 is not string"]),
    ],
)
def test_the_structural_checker_reports_each_kind_of_violation(
    value: object, shape: dict[str, Any], expected: list[str]
) -> None:
    assert violations(value, shape) == expected
