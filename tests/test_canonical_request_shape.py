"""`describe_canonical_request()` is true for what the codec writes and reads."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import fields
import inspect
from pathlib import Path
from typing import Any, Iterator

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.contract_shapes import INTEGER, STRING, array, literal, nullable, obj, record
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
from tests.utils.contract_shape_check import load_json, stored_json_files, violations

REPO_ROOT = Path(__file__).resolve().parents[1]
DESCRIPTION = describe_canonical_request()
TYPES: dict[str, Any] = DESCRIPTION["types"]
ROOT: dict[str, Any] = {"kind": "ref", "name": DESCRIPTION["root"]}



def _problems(value: object, shape: dict[str, Any]) -> list[str]:
    return violations(value, shape, TYPES, "renderRequest")


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
        match = next(option for option in shape["options"] if not _problems(value, option))
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
    assert [_problems(payload, ROOT) for payload, _resources in encoded_requests[name]] == [
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
# An entry with positions differs only there and agrees at every other position.
_COERCED = "accepted", "the typed constructor reads null as a value (str, bool, a default); no writer writes null"
_TOGETHER = "refused", "a decision's status, rationale and anchor are checked together; none encoded is skipped"
_REFERENCE = (
    "refused", "the plan's reference must equal the anchor of its reference record",
    ("linear:layout.similarityAlignment.reference",),
)
_CIRCULAR_PLACEMENT = "circular:diagramOptions.featurePlacements.[].placement"
_LINEAR_PLACEMENT = "linear:diagramOptions.featurePlacements.[].placement"
_LANE_SIDE = "a lane side belongs to one mode (FeaturePlacementOverride.from_mapping); the row is one typedef"
_KNOWN_DIFFERENCES: dict[str, tuple[str, str] | tuple[str, str, tuple[str, ...]]] = {
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
    "CanonicalAlignmentAnchorIdentity.biologicalFeatureId = another string": _REFERENCE,
    "CanonicalAlignmentAnchorIdentity.sourceFeatureIndex = another integer": _REFERENCE,
    "CanonicalAlignmentAnchorIdentity.sourceFeatureIndex = null": _REFERENCE,
    "CanonicalAlignmentAnchorIdentity.stableFeatureSvgId = another string": _REFERENCE,
    "CanonicalAlignmentAnchorIdentity.stableFeatureSvgId = null": _REFERENCE,
    "CanonicalFeaturePlacement.placement[lane].side = 'above'": ("refused", _LANE_SIDE, (_CIRCULAR_PLACEMENT,)),
    "CanonicalFeaturePlacement.placement[lane].side = 'below'": ("refused", _LANE_SIDE, (_CIRCULAR_PLACEMENT,)),
    "CanonicalFeaturePlacement.placement[lane].side = 'inward'": ("refused", _LANE_SIDE, (_LINEAR_PLACEMENT,)),
    "CanonicalFeaturePlacement.placement[lane].side = 'outward'": ("refused", _LANE_SIDE, (_LINEAR_PLACEMENT,)),
    "CanonicalRecord.cardinality = 'all'": (
        "refused", "a Circular record of every cardinality needs a grouping, a selector and a region that agree",
        ("circular:records.[]",),
    ),
    "CanonicalRequestOutput.interactiveMetadataPolicy = 'required'": (
        "refused", "requires interactive_svg among the formats; no request encoded there has it",
        ("circular:output.[]", "linear:output"),
    ),
    "CanonicalDepthTrack.height = a fraction": (
        "refused", "a depth track has a height only in a Linear diagram", ("circular:diagramOptions.depthTracks.[]",)
    ),
    "CanonicalDepthTrack.height = another integer": (
        "refused", "a depth track has a height only in a Linear diagram", ("circular:diagramOptions.depthTracks.[]",)
    ),
    "CanonicalRecordPresentation.gridColumn = another integer": (
        "refused", "grid_column requires grid_row; no Linear request encoded has a grid row",
        ("linear:records.[].presentation",),
    ),
}
# The decoder refuses another string, at every position: Python names the
# values the field takes and the description types it ``string``. The type is
# wider than Python there, as the Linear slot scalar was (OV-530, OV-597); a
# field joins this list only with a reason it cannot be a literal union.
_NAMES_TYPED_STRING: tuple[str, ...] = (
    "CanonicalCircularOutputOptions.legend",
    "CanonicalCircularTrackSlot.renderer",
    "CanonicalCircularTrackSlot.side",
    "CanonicalLinearOutputOptions.legend",
    "CanonicalLinearTrackSlot.renderer",
    "CanonicalLinearTrackSlot.side",
    "CanonicalProteinComparisonSettings.collinearityAnchorMode",
    "CanonicalProteinComparisonSettings.collinearityColorMode",
    "CanonicalProteinComparisonSettings.collinearitySearchScope",
    "CanonicalProteinComparisonSettings.collinearityUnitMode",
    "CanonicalProteinComparisonSettings.orthogroupMembershipMode",
)
# The decoder refuses a value of a kind the description allows, at every
# position: Python checks the value itself (a resource, record or feature that
# must exist, a range, a colour, an array that must hold an item), which no
# type can say.
_VALUE_CHECKED: dict[str, tuple[str, ...]] = {
    "CanonicalAlignmentAnchorIdentity.recordKey": ("another string",),
    "CanonicalAlignmentRecordDecision.recordKey": ("another string",),
    "CanonicalCircularDiagramOptions.depthTracks": ("an array",),
    "CanonicalCircularDiagramOptions.dinucleotide": ("another string",),
    "CanonicalCircularLayout.multiRecordMinRadiusRatio": ("another integer",),
    "CanonicalCircularRenderRequest.output": ("an array",),
    "CanonicalCircularRenderRequest.records": ("an array",),
    "CanonicalCollinearityResultComparison.resourceId": ("another string",),
    "CanonicalDepthTrack.source": ("an array",),
    "CanonicalFeaturePlacement.placement[lane].level": ("another integer",),
    "CanonicalFeatureSpan.selectors": ("an array",),
    "CanonicalFileRef.resourceId": ("another string",),
    "CanonicalGenbankRecordSource.resourceId": ("another string",),
    "CanonicalGffFastaRecordSource.fastaResourceId": ("another string",),
    "CanonicalGffFastaRecordSource.gffResourceId": ("another string",),
    "CanonicalHatchStyle.color": ("another string",),
    "CanonicalLinearDiagramOptions.depthTracks": ("an array",),
    "CanonicalLinearDiagramOptions.dinucleotide": ("another string",),
    "CanonicalLinearRecordTranslation.recordKey": ("another string",),
    "CanonicalLinearRenderRequest.records": ("an array",),
    "CanonicalLinearTrackOptions.linearTrackAxisIndex": ("another integer",),
    "CanonicalLinearTrackOptions.linearTrackSlots": ("an array",),
    "CanonicalNucleotideBlastComparison.queryRecordIndex": ("another integer",),
    "CanonicalNucleotideBlastComparison.resourceId": ("another string",),
    "CanonicalOrthogroupResultComparison.resourceId": ("another string",),
    "CanonicalPrecomputedProteinComparison.queryRecordIndex": ("another integer",),
    "CanonicalPrecomputedProteinComparison.resourceId": ("another string",),
    "CanonicalRegionAnnotationStyle.fill": ("another string",),
    "CanonicalRegionAnnotationStyle.fillOpacity": ("another integer",),
    "CanonicalRegionAnnotationStyle.labelColor": ("another string",),
    "CanonicalRegionAnnotationStyle.stroke": ("another string",),
    "CanonicalRequestOutput.formats": ("an array",),
    "CanonicalSimilarityAlignmentPlan.records": ("an array",),
    "CanonicalTableRef.resourceId": ("another string",),
}
# The decoder takes a value of a kind the description refuses, at every
# position: the typed constructor converts it (str, int, bool) or reads an
# empty array as no value. No writer writes it; the type is narrower than
# Python there (OV-528, OV-596).
_CONVERTED: dict[str, tuple[str, ...]] = {
    "CanonicalAnnotationSet.id": ("a fraction", "a null item", "an array", "another integer"),
    "CanonicalAnnotationSet.legendLabel": ("a fraction", "a null item", "an array", "another integer"),
    "CanonicalCircularDiagramOptions.conservationReference": ("an array",),
    "CanonicalCoordinateSpan.wrapsOrigin": ("an array",),
    "CanonicalFeatureSelector.key": ("a fraction", "a null item", "an array", "another integer"),
    "CanonicalFeatureSelector.value": ("a fraction", "a null item", "an array", "another integer"),
    "CanonicalGeneratedProteinComparison.mode": ("an array",),
    "CanonicalHatchStyle.cross": ("a fraction", "a null item", "an array", "another integer", "another string"),
    "CanonicalRegionAnnotation.id": ("a fraction", "a null item", "an array", "another integer"),
    "CanonicalRegionAnnotation.label": ("a fraction", "a null item", "an array", "another integer"),
    "CanonicalRegionAnnotation.lane": ("a fraction",),
    "CanonicalRegionAnnotation.legendLabel": ("a fraction", "a null item", "an array", "another integer"),
    "CanonicalRegionAnnotationStyle.fill": ("an array",),
}
_EXPECTED_DIFFERENCES: dict[str, tuple[str, list[str] | None]] = {
    **{
        change: (entry[0], list(entry[2]) if len(entry) == 3 else None)
        for change, entry in _KNOWN_DIFFERENCES.items()
    },
    **{
        f"{field} = {value}": (verdict, None)
        for verdict, table in (
            ("refused", {field: ("another string",) for field in _NAMES_TYPED_STRING}),
            ("refused", _VALUE_CHECKED),
            ("accepted", _CONVERTED),
        )
        for field, values in table.items()
        for value in values
    },
}


_UNNAMED = "valueTheDescriptionDoesNotName"


def _other_values(held: object) -> list[tuple[str, object]]:
    """A value of each JSON kind beside the one a field holds.

    The values do not depend on the description, so a description made wider
    (a literal as its scalar, an integer as a number, a string as a string or
    an integer, an array that takes null, a value beside null) allows one of
    them where the decoder refuses it.
    """

    integer = isinstance(held, int) and not isinstance(held, bool)
    return [
        ("another integer", held + 1 if integer else 7),
        ("a fraction", 0.5),
        ("another string", _UNNAMED),
        ("an array", []),
        ("a null item", [*held, None] if isinstance(held, list) else [None]),
    ]


def _position(mode: str, path: tuple[str | int, ...]) -> str:
    """Where an object is in a request of a mode; every item of an array is one position."""

    return f"{mode}:" + ".".join("[]" if isinstance(step, int) else step for step in path)


def _differences(
    encoded_requests: dict[str, list[Encoded]], tmp_path: Path
) -> tuple[dict[str, tuple[str, list[str] | None]], set[str]]:
    """Each change the decoder and the description judge differently, and the objects reached.

    Returns ``{change: (what the decoder does, positions)}``; the positions
    are ``None`` when the two differ wherever the change was tried. A
    cross-field rule can refuse one request and accept another, so at one
    position a change the description allows must decode in some request, and
    one it refuses in none.
    """

    accepted: dict[tuple[str, str, bool], set[bool]] = {}
    reached: set[str] = set()
    decodable = [
        item
        for found in encoded_requests.values()
        for item in found
        if not _problems(item[0], ROOT) and _decodes(*item, tmp_path)
    ]
    for mode in ("circular", "linear"):
        decodable.append(_smallest_request(mode, tmp_path / mode))
    for payload, resources in decodable:
        for path, shape in _objects(payload, ROOT):
            label = LABELS[id(shape)]
            reached.add(label)
            where = _position(payload["mode"], path)
            holder: Any = payload
            for step in path:
                holder = holder[step]
            changes: list[tuple[str, bool, tuple[object, ...]]] = [
                (f"{label}: an unknown field", shape["open"], ("fieldTheDescriptionDoesNotName", None)),
            ]
            for key, field in shape["fields"].items():
                if key in holder:
                    changes.append((f"{label}.{key} removed", key in shape["optional"], (key,)))
                changes.append((f"{label}.{key} = null", _problems(None, field) == [], (key, None)))
                changes += [
                    (f"{label}.{key} = {value!r}", True, (key, value))
                    for value in _literal_values(field)
                    if value is not None
                ]
                changes += [
                    (f"{label}.{key} = {name}", _problems(value, field) == [], (key, value))
                    for name, value in _other_values(holder.get(key))
                ]
            for change, allowed, arguments in changes:
                unchanged = len(arguments) == 2 and arguments[0] in holder and holder[arguments[0]] == arguments[1]
                accepted.setdefault((change, where, allowed), set()).add(
                    unchanged or _decodes(_with(payload, path, *arguments), resources, tmp_path)
                )
    found_by_change: dict[str, tuple[str, list[str], list[str]]] = {}
    for (change, where, allowed), found in sorted(accepted.items()):
        entry = found_by_change.setdefault(change, ("refused" if allowed else "accepted", [], []))
        entry[1 if (True in found) != allowed else 2].append(where)
    return {
        change: (verdict, differs if agrees else None)
        for change, (verdict, differs, agrees) in found_by_change.items()
        if differs
    }, reached


def test_the_decoder_agrees_with_the_description_at_every_position(
    encoded_requests: dict[str, list[Encoded]], tmp_path: Path
) -> None:
    """Remove each field, set it to null, to each named value and to a value just
    outside its shape, and add an unknown field: at each position of each mode.

    A typedef used at several positions is judged at each one, so a value one
    position accepts cannot hide that another always refuses it (OV-530).
    """

    differences, reached = _differences(encoded_requests, tmp_path)
    assert differences == _EXPECTED_DIFFERENCES
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
    requests += [request for path in stored_json_files() for request in _stored_requests(load_json(path))]
    for payload in requests:
        if _problems(payload, ROOT):
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


def test_requests_stored_in_current_sessions_match_the_description() -> None:
    """The Web wrote the Gallery Sessions; the CLI wrote the documentation Sessions."""

    modes: list[str] = []
    problems: dict[str, list[str]] = {}
    for path in stored_json_files():
        for request in _stored_requests(load_json(path)):
            modes.append(request["mode"])
            found = _problems(request, ROOT)
            if found:
                problems[str(path.relative_to(REPO_ROOT))] = found
    assert problems == {}
    assert modes.count("circular") >= 4 and modes.count("linear") >= 4


@pytest.mark.parametrize(
    ("value", "shape", "expected"),
    [
        ({"a": 1}, obj({}), ["renderRequest.a: unknown field"]),
        ({"a": 1}, obj({}, open=True), []),
        ({}, obj({"a": INTEGER}), ["renderRequest.a: required field is absent"]),
        ({"a": "x", "b": 1}, obj({"a": INTEGER}, open=True), ["renderRequest.a: 'x' is not integer"]),
        ({}, obj({"a": INTEGER}, optional=["a"]), []),
        (True, INTEGER, ["renderRequest: True is not integer"]),
        (1, literal(True), ["renderRequest: 1 is not one of [True]"]),
        ([1, "a"], array(INTEGER), ["renderRequest[1]: 'a' is not integer"]),
        (None, nullable(STRING), []),
        ({"k": 1}, record(STRING), ["renderRequest.k: 1 is not string"]),
    ],
)
def test_the_structural_checker_reports_each_kind_of_violation(
    value: object, shape: dict[str, Any], expected: list[str]
) -> None:
    assert _problems(value, shape) == expected
