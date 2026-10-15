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

from gbdraw.api.options import (
    CircularDiagramOptions,
    CircularRequestTrackOptions,
    LinearDiagramOptions,
    LinearRequestTrackOptions,
    LosatSearchOptions,
)
from gbdraw.contract_shapes import violations
from gbdraw.api.requests import (
    CircularBatchRequest,
    CircularDiagramRequest,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.session_request_codec import (
    CANONICAL_REQUEST_SCHEMA,
    SUPPORTED_CANONICAL_REQUEST_SCHEMAS,
    CanonicalRequestCodecError,
    CanonicalRequestDecodingError,
    CanonicalRequestEncodingError,
    _SHARED_OPTION_WRONG_MODE_DEFAULTS,
    _option_wire_key,
    decode_canonical_request,
    describe_canonical_request,
    encode_canonical_request,
)
from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot
from gbdraw.web_support.error_adapter import serialize_web_error
from tests import test_session_request_codec as codec_tests
from tests.utils.contract_shape_check import load_json, stored_json_files

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


def _encoded(request: Any, tmp_path: Path) -> Encoded:
    encoded = encode_canonical_request(request)
    return encoded.payload, codec_tests._materialize_resources(encoded, tmp_path / "resources")


def _smallest_request(mode: str, tmp_path: Path) -> Encoded:
    return _encoded(_request(mode), tmp_path)


# Where one change cannot show that the decoder and the description agree, and why.
# An entry with positions differs only there and agrees at every other position.
_TOGETHER = "refused", "a decision's status, rationale and anchor are checked together; none encoded is skipped"
_REFERENCE = (
    "refused", "the plan's reference must equal the anchor of its reference record",
    ("linear:layout.similarityAlignment.reference",),
)
_KNOWN_DIFFERENCES: dict[str, tuple[str, str] | tuple[str, str, tuple[str, ...]]] = {
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
    "CanonicalRecord.cardinality = 'all'": (
        "refused", "a Circular record of every cardinality needs a grouping, a selector and a region that agree",
        ("circular:records.[]",),
    ),
    "CanonicalRequestOutput.interactiveMetadataPolicy = 'required'": (
        "refused", "requires interactive_svg among the formats; no request encoded there has it",
        ("circular:output.[]", "linear:output"),
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
# - A legend position: the Web writer writes the text its state holds, which a
#   stored Session can set to any text (OV-652); it narrows once the Web writes
#   only a position of the mode's set.
# - A slot's renderer and side: the track planner reads other spellings of a
#   name, and the Session writer of request schemas 5 to 9 stored a caller's
#   spelling until the encoder wrote the name (_SLOT_SPELLINGS, OV-653).
_NAMES_TYPED_STRING: tuple[str, ...] = (
    "CanonicalCircularOutputOptions.legend",
    "CanonicalCircularTrackSlot.renderer",
    "CanonicalCircularTrackSlot.side",
    "CanonicalLinearOutputOptions.legend",
    "CanonicalLinearTrackSlot.renderer",
    "CanonicalLinearTrackSlot.side",
)
# The decoder refuses a value of a kind the description allows, at every
# position: Python checks that a resource, a record or a slot the value names
# exists, which the description cannot say (a reader of a stored request
# checks it against the file). The value domains Python checks, colours
# included, are bounds of the description.
_VALUE_CHECKED: dict[str, tuple[str, ...]] = {
    "CanonicalAlignmentAnchorIdentity.recordKey": ("another string",),
    "CanonicalAlignmentRecordDecision.recordKey": ("another string",),
    "CanonicalCollinearityResultComparison.resourceId": ("another string",),
    "CanonicalFileRef.resourceId": ("another string",),
    "CanonicalGenbankRecordSource.resourceId": ("another string",),
    "CanonicalGffFastaRecordSource.fastaResourceId": ("another string",),
    "CanonicalGffFastaRecordSource.gffResourceId": ("another string",),
    "CanonicalLinearRecordTranslation.recordKey": ("another string",),
    "CanonicalLinearTrackOptions.linearTrackAxisIndex": ("another integer",),
    "CanonicalNucleotideBlastComparison.queryRecordIndex": ("another integer",),
    "CanonicalNucleotideBlastComparison.resourceId": ("another string",),
    "CanonicalOrthogroupResultComparison.resourceId": ("another string",),
    "CanonicalPrecomputedProteinComparison.queryRecordIndex": ("another integer",),
    "CanonicalPrecomputedProteinComparison.resourceId": ("another string",),
    "CanonicalTableRef.resourceId": ("another string",),
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
        )
        for field, values in table.items()
        for value in values
    },
}


_UNNAMED = "valueTheDescriptionDoesNotName"
# Numbers of the right kind outside the range the decoder checks: a depth
# track's height (Linear only; a Circular one holds null), tick intervals and
# tick font size are greater than 0 (DEPTH_TRACK_POSITIVE_FIELDS), which the
# description bounds (U7, BA-OD73).
_OUT_OF_RANGE: dict[str, tuple[tuple[str, object], ...]] = {
    f"Canonical{mode}DepthTrack.{key}": (("zero", 0), ("a negative number", -1))
    for mode, keys in (("Circular", ()), ("Linear", ("height",)))
    for key in (*keys, "largeTickInterval", "smallTickInterval", "tickFontSize")
}


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
                    for name, value in [*_other_values(holder.get(key)), *_OUT_OF_RANGE.get(f"{label}.{key}", ())]
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


# A slot built with another spelling of a renderer or a side, and the names
# the track planner reads it as (a renderer alias is the name it stands for).
_SLOT_SPELLINGS: dict[str, tuple[Any, Any, tuple[tuple[str, tuple[str, str], tuple[str, str]], ...]]] = {
    "circular": (CircularTrackSlot, CircularRequestTrackOptions, (
        ("features", ("Features", "features"), (" OUTSIDE ", "outside")),
        ("gc_content", ("gc_content", "dinucleotide_content"), ("INSIDE", "inside")),
    )),
    "linear": (LinearTrackSlot, LinearRequestTrackOptions, (
        ("features", (" FEATURES ", "features"), ("Above", "above")),
        ("gc_content", ("GC_CONTENT", "dinucleotide_content"), ("BELOW", "below")),
    )),
}


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_the_encoder_writes_the_name_of_a_slot_renderer_and_side(mode: str, tmp_path: Path) -> None:
    """A Python caller's spelling is written as the name (OV-653, BA-D23)."""

    slot_type, tracks_type, spellings = _SLOT_SPELLINGS[mode]

    def tracks(index: int) -> Any:
        return tracks_type(**{
            f"{mode}_track_slots": tuple(
                slot_type(id=slot_id, renderer=renderer[index], side=side[index])
                for slot_id, renderer, side in spellings
            ),
            f"{mode}_track_axis_index": 1,
        })

    payload, resources = _encoded(_request(mode, tracks=tracks(0)), tmp_path)
    written = payload["diagramOptions"]["tracks"][f"{mode}TrackSlots"]
    assert [(slot["id"], slot["renderer"], slot["side"]) for slot in written] == [
        (slot_id, renderer[1], side[1]) for slot_id, renderer, side in spellings
    ]
    assert _problems(payload, ROOT) == []
    decoded = decode_canonical_request(payload, resource_paths=resources, output_directory=tmp_path / "output")
    assert decoded.options.tracks == tracks(1)


# The decoder reads these keys and no writer writes them: the depth inputs of
# older writers (the writer stores every depth input as depthTracks).
_DEPTH_INPUTS_OF_OLDER_WRITERS = {
    "depthFile", "depthFiles", "depthTable", "depthTables", "depthTrackColors", "depthTrackFiles",
    "depthTrackLabels", "depthTrackLargeTickIntervals", "depthTrackSmallTickIntervals",
    "depthTrackTables", "depthTrackTickFontSizes",
}
_DECODER_ONLY_OPTIONS = {
    "circular": ("conservation_table_file", _DEPTH_INPUTS_OF_OLDER_WRITERS),
    "linear": ("comparison_table_file", _DEPTH_INPUTS_OF_OLDER_WRITERS | {"depthTrackHeights"}),
}


@pytest.mark.parametrize("schema", [8, CANONICAL_REQUEST_SCHEMA])
@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_the_description_names_every_option_key_a_writer_can_write(
    mode: str, schema: int, tmp_path: Path
) -> None:
    """At the current schema the decoder reads the described option keys only;
    an older schema also reads the depth inputs of older writers."""

    table_path, older_depth_inputs = _DECODER_ONLY_OPTIONS[mode]
    described = TYPES[f"Canonical{mode.capitalize()}DiagramOptions"]["fields"]
    encoded = encode_canonical_request(_request(mode))
    payload = codec_tests._payload_for_schema(encoded, schema)
    resources = codec_tests._materialize_resources(encoded, tmp_path / "resources")

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
                if "unknown field" in str(error).lower() and key in str(error):
                    continue
            return True
        return False

    candidates = {_option_wire_key(item.name) for item in fields(_REQUEST_TYPES[mode][1])}
    if schema == CANONICAL_REQUEST_SCHEMA:
        assert {key for key in described if not decoder_names(key)} == set()
        # Any other key is refused, by the decoder or as outside the description.
        for key in (candidates | set(defaults)) - set(described):
            with pytest.raises(CanonicalRequestCodecError):
                decode_canonical_request(
                    _with(payload, ("diagramOptions",), key, None),
                    resource_paths=resources,
                    output_directory=tmp_path / "output",
                )
    else:
        decoder_keys = {key for key in candidates | set(defaults) | set(described) if decoder_names(key)}
        assert set(described) - decoder_keys == {"featureOverrides"}
        assert decoder_keys - set(described) == older_depth_inputs
    with pytest.raises(CanonicalRequestEncodingError, match="comparison table"):
        encode_canonical_request(_request(mode, **{table_path: "table.tsv"}))


def test_the_encoder_refuses_a_comparison_table_path_in_a_circular_batch() -> None:
    request = CircularBatchRequest(
        records=(RecordInput(source=InMemoryRecordSource(_record())),),
        options=CircularDiagramOptions(conservation_table_file="table.tsv"),
        outputs=(RenderOutputRequest(),),
    )
    with pytest.raises(CanonicalRequestEncodingError, match="Circular comparison table"):
        encode_canonical_request(request)


SCHEMAS = sorted(SUPPORTED_CANONICAL_REQUEST_SCHEMAS)


def _decode_at(schema: int, request: Any, tmp_path: Path, change: Any = None) -> Any:
    """Decode ``request`` as a payload of ``schema``, after ``change`` edits the payload."""

    encoded = encode_canonical_request(request)
    payload = codec_tests._payload_for_schema(encoded, schema)
    if schema == 1:
        # Schema 1 stored no record pairs of a protein comparison.
        for comparison in payload["comparisons"]:
            comparison.pop("pairs", None)
    if change is not None:
        change(payload)
    return decode_canonical_request(
        payload,
        resource_paths=codec_tests._materialize_resources(encoded, tmp_path / f"resources-{schema}"),
        output_directory=tmp_path / "output",
    )


@pytest.mark.parametrize("schema", SCHEMAS)
@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_the_decoder_refuses_a_comparison_table_path_at_every_schema(
    mode: str, schema: int, tmp_path: Path
) -> None:
    """No writer wrote a table path: the Web never, the encoder refuses it (OV-324)."""

    key = _option_wire_key(_DECODER_ONLY_OPTIONS[mode][0])
    _decode_at(schema, _request(mode), tmp_path)

    def add_path(payload: dict[str, Any]) -> None:
        payload["diagramOptions"][key] = "table.tsv"

    with pytest.raises(CanonicalRequestDecodingError, match=key):
        _decode_at(schema, _request(mode), tmp_path, add_path)


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


# Defaults that a rule across records or the grouping can refuse: one record's
# grid cell (every Circular record has one, or none), so also its whole
# presentation, a Linear layout's translations (an alignment plan needs one
# per record), and a batch's output (one per record).
_DEFAULTS_DECIDED_BY_OTHER_PARTS = frozenset({
    "CanonicalRecord.presentation",
    "CanonicalLinearLayout.recordTranslations",
    "CanonicalRecordPresentation.gridColumn",
    "CanonicalRecordPresentation.gridRow",
    "CanonicalCircularRenderRequest.output",
})


def test_a_described_default_decodes_as_the_absent_field(
    encoded_requests: dict[str, list[Encoded]], tmp_path: Path
) -> None:
    """The Web writes a description's default where a stored field is absent
    or unreadable (stage 1d PLAN 12.4 rule 3): the decoder must read an
    optional field's default as the absent field, and a required field's
    default as a value it takes."""

    exercised: set[str] = set()
    for found in encoded_requests.values():
        for payload, resources in found:
            for path, shape in _objects(payload, ROOT):
                for key, default in shape.get("defaults", {}).items():
                    label = f"{LABELS[id(shape)]}.{key}"
                    # A required field's default is a value it takes in every
                    # request (OV-777), except where another record or the
                    # grouping decides it (_DEFAULTS_DECIDED_BY_OTHER_PARTS).
                    if label in exercised and (
                        key in shape["optional"] or label in _DEFAULTS_DECIDED_BY_OTHER_PARTS
                    ):
                        continue
                    absent = deepcopy(payload)
                    target: Any = absent
                    for step in path:
                        target = target[step]
                    target.pop(key, None)
                    output = tmp_path / "output"
                    # The encoder writes a decoded request in its one form.
                    decode = lambda value: encode_canonical_request(decode_canonical_request(  # noqa: E731
                        value, resource_paths=resources, output_directory=output
                    )).payload
                    if key in shape["optional"]:
                        assert decode(_with(absent, path, key, default)) == decode(absent), label
                    else:
                        decode(_with(absent, path, key, default))
                    exercised.add(label)
    described = {f"{label}.{key}" for label, shape in OBJECTS for key in shape.get("defaults", {})}
    assert exercised == described


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


def _annotated_request(mode: str) -> Any:
    """A request whose annotations hold each annotation field the decoder converted (OV-528, OV-596)."""

    from gbdraw.annotations import (
        AnnotationOptions, AnnotationSet, CoordinateSpan, FeatureSelector, FeatureSpan, HatchStyle,
        RegionAnnotation, RegionAnnotationStyle,
    )
    from gbdraw.io.record_select import parse_record_selector

    annotations = AnnotationOptions(sets=(AnnotationSet("review", (
        RegionAnnotation(
            "coordinates", CoordinateSpan(parse_record_selector("#1"), 2, 20), label="Review", lane=1,
            style=RegionAnnotationStyle(fill="#88aaff", hatch=HatchStyle(cross=True)), legend_label="Reviewed",
        ),
        RegionAnnotation(
            "features", FeatureSpan(parse_record_selector("#1"), (FeatureSelector("dnaA", "gene"),)),
        ),
    ), legend_label="Review regions"),))
    options: dict[str, object] = {"annotations": annotations}
    if mode == "linear":
        options["losat_search"] = LosatSearchOptions(program="losatp", losatp_mode="collinear")
    return _request(mode, **options)


_SET = ("diagramOptions", "annotations", "sets", 0)
_COORDINATES = (*_SET, "annotations", 0)
_FEATURES = (*_SET, "annotations", 1)
_TEXT_FORMS = (("null", None), ("a fraction", 0.5), ("an integer", 7), ("an array", []), ("a null item", [None]))
_OPTIONAL_TEXT_FORMS = _TEXT_FORMS[1:]
# Each field and the forms of another kind the decoder converted at every
# schema (the test tables _COERCED and _CONVERTED before stage 1d PR-B); no
# writer wrote any of them (stage1d/audit.md section 3).
_UNWRITTEN_FORMS: dict[tuple[str, tuple[str | int, ...], str], tuple[tuple[str, object], ...]] = {
    ("linear", _SET, "id"): _TEXT_FORMS,
    ("linear", _SET, "legendLabel"): _OPTIONAL_TEXT_FORMS,
    ("linear", _COORDINATES, "id"): _TEXT_FORMS,
    ("linear", _COORDINATES, "label"): _TEXT_FORMS,
    ("linear", _COORDINATES, "legendLabel"): _OPTIONAL_TEXT_FORMS,
    ("linear", (*_COORDINATES, "target"), "wrapsOrigin"): (("null", None), ("an array", []), ("an integer", 0)),
    ("linear", (*_FEATURES, "target", "selectors", 0), "key"): _OPTIONAL_TEXT_FORMS,
    ("linear", (*_FEATURES, "target", "selectors", 0), "value"): _TEXT_FORMS,
    ("linear", (*_COORDINATES, "style", "hatch"), "cross"): (*_TEXT_FORMS, ("a string", "true")),
    ("linear", (*_COORDINATES, "style"), "fill"): (("an array", []),),
    ("linear", ("comparisons", 0), "mode"): (("null", None), ("an array", [])),
    ("linear", ("diagramOptions",), "pairwiseMatchStyle"): (("null", None),),
    ("circular", ("diagramOptions",), "conservationReference"): (("null", None), ("an array", [])),
}


@pytest.mark.parametrize("schema", SCHEMAS)
@pytest.mark.parametrize(
    ("mode", "path", "key", "name", "value"),
    [
        (mode, path, key, name, value)
        for (mode, path, key), forms in _UNWRITTEN_FORMS.items()
        for name, value in forms
    ],
)
def test_the_decoder_refuses_a_value_of_another_kind_at_every_schema(
    mode: str, path: tuple[str | int, ...], key: str, name: str, value: object, schema: int, tmp_path: Path
) -> None:
    request = _annotated_request(mode)
    _decode_at(schema, request, tmp_path)

    def change(payload: dict[str, Any]) -> None:
        holder: Any = payload
        for step in path:
            holder = holder[step]
        holder[key] = value

    with pytest.raises(CanonicalRequestDecodingError, match=key):
        _decode_at(schema, request, tmp_path, change)


def _set_at(path: tuple[str | int, ...], key: str, value: object) -> Any:
    def change(payload: dict[str, Any]) -> None:
        holder: Any = payload
        for step in path:
            holder = holder[step]
        holder[key] = value

    return change


@pytest.mark.parametrize("schema", SCHEMAS)
@pytest.mark.parametrize(
    ("path", "key", "value"),
    [
        ((*_COORDINATES, "target"), "wrapsOrigin", 0),
        (_COORDINATES, "legendLabel", 7),
        (("diagramOptions",), "comparisonTableFile", "table.tsv"),
    ],
)
def test_the_web_reports_a_refused_request_value_as_an_invalid_input(
    path: tuple[str | int, ...], key: str, value: object, schema: int, tmp_path: Path
) -> None:
    """The Web shows each refusal as INPUT_INVALID, not as an unclassified error (R6)."""

    with pytest.raises(CanonicalRequestDecodingError, match=key) as refused:
        _decode_at(schema, _annotated_request("linear"), tmp_path, _set_at(path, key, value))
    shown = serialize_web_error(refused.value, operation="generate", stage="request-validation")
    assert (shown["code"], shown["context"]) == ("INPUT_INVALID", {"field": "schema", "reason": "FIELDS"})


# The Web's Lane, Start and End number inputs write fractions at every schema
# (OV-630, BA-OD74); Python reads the integer part.
_FRACTIONS: tuple[tuple[tuple[str | int, ...], str, float, int], ...] = (
    (_COORDINATES, "lane", 0.5, 0),
    ((*_COORDINATES, "target"), "start", 2.5, 2),
    ((*_COORDINATES, "target"), "end", 20.5, 20),
)


@pytest.mark.parametrize("schema", SCHEMAS)
@pytest.mark.parametrize(("path", "key", "fraction", "integer"), _FRACTIONS)
def test_a_fraction_at_an_annotation_number_reads_as_its_integer_part_at_every_schema(
    path: tuple[str | int, ...], key: str, fraction: float, integer: int, schema: int, tmp_path: Path
) -> None:
    def change(payload: dict[str, Any]) -> None:
        holder: Any = payload
        for step in path:
            holder = holder[step]
        holder[key] = fraction

    annotation = _decode_at(schema, _annotated_request("linear"), tmp_path, change).options.annotations.sets[0].annotations[0]
    assert {"lane": annotation.lane, **vars(annotation.target)}[key] == integer


def test_the_encoder_writes_a_boolean_where_a_coordinate_span_wraps_the_origin(tmp_path: Path) -> None:
    """A caller that gives wraps_origin another kind still gets a boolean (OV-632)."""

    from gbdraw.annotations import AnnotationOptions, AnnotationSet, CoordinateSpan, RegionAnnotation

    written = []
    for wraps, start, end in ((None, 10, 30), (1, 30, 10), ([], 10, 30), ("yes", 30, 10)):
        span = CoordinateSpan(None, start, end, wraps_origin=wraps)
        annotations = AnnotationOptions(sets=(AnnotationSet("s", (RegionAnnotation("a", span),)),))
        payload, resources = _encoded(_request("linear", annotations=annotations), tmp_path / str(len(written)))
        written.append(payload["diagramOptions"]["annotations"]["sets"][0]["annotations"][0]["target"]["wrapsOrigin"])
        decode_canonical_request(payload, resource_paths=resources, output_directory=tmp_path / "output")
    assert written == [False, True, False, True]


# Another spelling of a comparison setting and the name the normalizers read
# it as. Python stores a setting normalized, and the Web writes the name from
# schema 6 (OV-651).
_SETTING_SPELLINGS: dict[str, tuple[tuple[str, str], ...]] = {
    "collinearityUnitMode": ((" LOCUS ", "locus"), ("", "auto")),
    "collinearityAnchorMode": (("Reciprocal-Best", "rbh"), ("top_n", "all"), ("", "rbh")),
    "collinearitySearchScope": (("ALL", "all"), ("", "adjacent")),
    "collinearityColorMode": (("identity", "average_identity"), ("Orientation-Identity", "orientation_identity")),
    "orthogroupMembershipMode": (("rbh", "anchor_core_v1"), ("family_merge", "anchor_core_v1"), ("", "anchor_core_v1")),
}


@pytest.mark.parametrize("schema", SCHEMAS)
@pytest.mark.parametrize(
    ("key", "spelling", "name"),
    [(key, spelling, name) for key, found in _SETTING_SPELLINGS.items() for spelling, name in found],
)
def test_a_comparison_setting_holds_its_name_from_schema_9(
    key: str, spelling: str, name: str, schema: int, tmp_path: Path
) -> None:
    """Older schemas keep reading a spelling; schema 9 reads the name only."""

    request = _request("linear", losat_search=LosatSearchOptions(program="losatp", losatp_mode="collinear"))

    def spell(payload: dict[str, Any]) -> None:
        payload["comparisons"][0]["settings"][key] = spelling

    if schema >= 9:
        with pytest.raises(CanonicalRequestDecodingError, match=key):
            _decode_at(schema, request, tmp_path, spell)
        return
    decoded = _decode_at(schema, request, tmp_path, spell)
    assert encode_canonical_request(decoded).payload["comparisons"][0]["settings"][key] == name


@pytest.mark.parametrize("schema", SCHEMAS)
@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_a_current_request_outside_its_description_is_refused(mode: str, schema: int, tmp_path: Path) -> None:
    """The decoder reads a depth input of older writers, which schema 9 writers never wrote (K11)."""

    def older_depth_input(payload: dict[str, Any]) -> None:
        payload["diagramOptions"]["depthTrackColors"] = None

    if schema < CANONICAL_REQUEST_SCHEMA:
        _decode_at(schema, _request(mode), tmp_path, older_depth_input)
        return
    with pytest.raises(
        CanonicalRequestDecodingError,
        match=r"does not match its description: renderRequest\.diagramOptions\.depthTrackColors: unknown field",
    ):
        _decode_at(schema, _request(mode), tmp_path, older_depth_input)


@pytest.mark.parametrize("schema", [schema for schema in SCHEMAS if 2 <= schema < 6])
@pytest.mark.parametrize(("row", "grid_row"), [("2", 2), ("١", None), ("²", None), ("x", None)])
def test_older_linear_card_rows_are_read_as_ascii_integers_only(
    schema: int, row: str, grid_row: int | None, tmp_path: Path
) -> None:
    """Schemas 2 to 5 wrote ``#<card>@<row>``; a row with other digits is not a row."""

    def positions(payload: dict[str, Any]) -> None:
        payload["layout"] = {"recordGapPx": 10, "multiRecordPositions": [f"#1@{row}"]}

    decoded = _decode_at(schema, _request("linear"), tmp_path, positions)
    assert [record.presentation.grid_row for record in decoded.records] == [grid_row]
