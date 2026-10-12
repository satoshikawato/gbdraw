"""`describe_canonical_request()` is true for what the codec writes and reads."""

from __future__ import annotations

from copy import deepcopy
import gzip
import inspect
import json
from pathlib import Path
from typing import Any, Iterator

import pytest

from gbdraw.session_request_codec import (
    CANONICAL_REQUEST_SCHEMA,
    CanonicalRequestCodecError,
    EncodedCanonicalRequest,
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


def _mutated(payload: dict[str, Any], path: tuple[str | int, ...], key: str, *, add: bool) -> dict[str, Any]:
    result = deepcopy(payload)
    holder: Any = result
    for step in path:
        holder = holder[step]
    if add:
        holder[key] = None
    else:
        del holder[key]
    return result


def _checked_encoder(seen: list[EncodedCanonicalRequest]):
    def encode(*args: Any, **kwargs: Any) -> EncodedCanonicalRequest:
        encoded = encode_canonical_request(*args, **kwargs)
        assert violations(encoded.payload, ROOT) == []
        seen.append(encoded)
        return encoded

    return encode


# Every test of the codec that takes only a temporary directory builds its own
# typed request; they are the representative requests of both modes.
_CODEC_TESTS = sorted(
    name
    for name, function in vars(codec_tests).items()
    if name.startswith("test_")
    and inspect.isfunction(function)
    and list(inspect.signature(function).parameters) == ["tmp_path"]
)
# Requests with records, layout, comparisons, annotations, track slots, depth
# tracks, resources and a batch: the decoder is checked object by object.
_DECODED_CODEC_TESTS = (
    "test_schema6_round_trips_circular_batch_grouping_and_outputs",
    "test_schema8_round_trips_typed_similarity_alignment_and_translations",
    "test_circular_request_payload_round_trip_uses_stable_resources",
    "test_linear_comparison_kinds_and_payload_round_trip",
    "test_file_backed_options_and_typed_config_round_trip",
    "test_canonical_depth_tracks_round_trip_mixed_sources",
    "test_annotation_targets_and_styles_round_trip",
)


@pytest.mark.parametrize("name", _CODEC_TESTS)
def test_every_encoded_request_matches_the_description(
    name: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[EncodedCanonicalRequest] = []
    monkeypatch.setattr(codec_tests, "encode_canonical_request", _checked_encoder(seen))

    getattr(codec_tests, name)(tmp_path)


def test_the_decoded_codec_tests_exist() -> None:
    assert set(_DECODED_CODEC_TESTS) <= set(_CODEC_TESTS)


@pytest.mark.parametrize("name", _DECODED_CODEC_TESTS)
def test_the_decoder_refuses_what_the_description_refuses(
    name: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[EncodedCanonicalRequest] = []
    monkeypatch.setattr(codec_tests, "encode_canonical_request", _checked_encoder(seen))
    getattr(codec_tests, name)(tmp_path)
    assert seen

    accepted: list[str] = []
    for number, encoded in enumerate(seen):
        resource_paths = codec_tests._materialize_resources(encoded, tmp_path / f"shape-{number}")
        mutations = [
            (path, key, add)
            for path, shape in _objects(encoded.payload, ROOT)
            for key, add in [
                *((key, False) for key in shape["fields"] if key not in shape["optional"]),
                ("fieldTheDescriptionDoesNotName", True),
            ]
        ]
        for path, key, add in mutations:
            try:
                decode_canonical_request(
                    _mutated(encoded.payload, path, key, add=add),
                    resource_paths=resource_paths,
                    output_directory=tmp_path / "output",
                )
            except CanonicalRequestCodecError:
                continue
            accepted.append(f"{'added' if add else 'removed'} {'.'.join(map(str, (*path, key)))}")
    assert accepted == []


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
