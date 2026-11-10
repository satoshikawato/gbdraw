"""Worker helpers that plan the drawn units of a Linear request for the Web: its comparison searches (D-01, D-37)
and the numbers and omitted lengths a record card shows beside its windows (PD-OI-124).

Before its comparison searches run, a Generate sends the records its request will draw (each record key with its
windows and grid row, as the request writes them) and the record pairs its comparison plan compares. Python answers
with the drawn units, in drawing order, with the unit pairs of each record pair: every window of one record against
every window of the other, the pairs a CLI search runs (``comparison_unit_pairs``), and with the unit pairs the
drawing displays when its units stand in a grid (``grid_display_pairs``, the rule the CLI drawing and the derived
LOSATP identity read; OV-3248). The page runs its searches over these
units and writes their indexes into the request, which Python reads as drawn-unit indexes. A comparison file
(an uploaded BLAST TSV) gives the coordinates of one drawn unit, so a pair of a record drawn in several windows must
be searched.

A record card that lists its windows shows each window's circled number in the drawing and the bases left out
between two windows of the record that stand next to each other in their row (``plan_card_windows``). Python answers
from the owners the drawing uses: the row order of the drawn units (``grid_placement_order``) and the row's captions
(``plan_row_windows``), so the card shows the drawing's own values.

Before its LOSAT searches run, a Generate also sends its drawing's request, as Save writes it, and Python answers
with the E-value scale of every search between two drawn units (``plan_genome_scales``): the search space of the
genome pair (LOSATN, TLOSATX) or the proteins of the subject genome (LOSATP), from the owners the command line's
searches read (``linear_genome_scales``), so the page passes the command line's values and never sums a genome.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Literal, TypedDict, cast, is_typeddict

from gbdraw.api.record_planning import (
    ResolvedWindow,
    comparison_unit_pairs,
    numbered_places,
    record_drawn_windows,
    refuse_file_comparison_of_windows,
)
from gbdraw.api.request_render import (
    GENOME_SCALE_PROGRAMS,
    GenomeScaleProgram,
    grid_placement_order,
    linear_genome_scales,
)
from gbdraw.api.requests import InputPlace, LinearDiagramRequest, RecordPresentation, expanded_record_key
from gbdraw.contract_shapes import Shape, ref, typed_dict_shape
from gbdraw.exceptions import ValidationError
from gbdraw.io.regions import window_drawn_reverse
from gbdraw.layout.linear_multi_record import drawn_unit_rows, grid_display_pairs
from gbdraw.layout.linear_row_windows import plan_row_windows
from gbdraw.layout.record_units import RecordUnit
from gbdraw.session_request_codec import (
    CanonicalRequestDecodingError,
    _decode_windows,
    decode_canonical_request,
)
from gbdraw.web_support.request_render import _staged_source_names


class LinearUnitRecordInput(TypedDict):
    """A record the request draws: its key, its windows (``renderRequest.records[].windows``) and its grid row
    (``renderRequest.records[].presentation.gridRow``: every record has one or none has)."""

    recordKey: str
    windows: list[object]
    gridRow: int | None


class LinearUnitPairInput(TypedDict):
    """A record pair of the comparison plan, by record keys, and how it is compared."""

    queryRecordKey: str
    subjectRecordKey: str
    source: Literal["losat", "upload"]


class LinearUnitWindow(TypedDict):
    """The span of its record a unit draws, in source coordinates, 1-based and inclusive (start after end: across
    the origin), and ``searchRegion``, the region string of the span a comparison search reads: the crop before its
    reverse complement, in the form the Worker's extraction helpers read (``region_text``)."""

    start: int
    end: int
    reverseComplement: bool
    searchRegion: str


class LinearDrawnUnit(TypedDict):
    """One drawn unit: its key, the key of its record, ``windowOf`` (the record key when the unit is one of several
    windows of it, else null), and its window (null for a record drawn whole)."""

    unitKey: str
    recordKey: str
    windowOf: str | None
    window: LinearUnitWindow | None


class LinearUnitPair(TypedDict):
    """One unit pair: the index of the record pair it belongs to and its query and subject units."""

    pair: int
    queryUnit: int
    subjectUnit: int


class LinearUnitGridPair(TypedDict):
    """A pair of units in adjacent rows that the drawing displays, the unit drawn first in request order first."""

    queryUnit: int
    subjectUnit: int


class LinearUnitPlan(TypedDict):
    """The drawn units in drawing order, the unit pairs, in the order of the record pairs, and ``gridPairs``: the unit
    pairs the drawing displays when its units stand in a grid (every pair of units in adjacent rows), else null
    (each row holds one unit, in unit order: its comparisons join consecutive units)."""

    units: list[LinearDrawnUnit]
    pairs: list[LinearUnitPair]
    gridPairs: list[LinearUnitGridPair] | None


def _protocol_error(message: str) -> ValidationError:
    return ValidationError(message, diagnostic={"code": "HELPER_PROTOCOL"})


def _key(value: object, what: str) -> str:
    if not isinstance(value, str) or not value:
        raise _protocol_error(f"A {what} needs a record key.")
    return value


def plan_linear_units(records: Sequence[object], pairs: Sequence[object]) -> LinearUnitPlan:
    """The drawn units of ``records`` and the unit pairs of the record ``pairs``."""

    units: list[LinearDrawnUnit] = []
    units_by_record: dict[str, list[int]] = {}
    unit_grid_rows: list[int | None] = []
    for position, record in enumerate(records):
        if not isinstance(record, Mapping):
            raise _protocol_error("A record of a unit plan must be an object.")
        record_key = _key(record.get("recordKey"), "record")
        if record_key in units_by_record:
            raise _protocol_error(f"Record {record_key!r} is named twice.")
        if "gridRow" not in record:
            raise _protocol_error(f"Record {record_key!r} has no gridRow: give a positive whole number or null.")
        grid_row = record["gridRow"]
        if grid_row is not None and not _positive_whole_number(grid_row):
            raise _protocol_error(f"Record {record_key!r} has an invalid gridRow: give a positive whole number or null.")
        try:
            windows = _decode_windows(record.get("windows"), path=f"records[{position}].windows")
        except (CanonicalRequestDecodingError, ValueError) as exc:
            raise _protocol_error(str(exc)) from exc
        indexes = units_by_record[record_key] = []
        for window in record_drawn_windows(record_key, windows):
            indexes.append(len(units))
            unit_grid_rows.append(cast(int | None, grid_row))
            units.append({
                "unitKey": record_key if window is None else window.unit_key,
                "recordKey": record_key,
                "windowOf": None if window is None else window.window_of,
                "window": None if window is None else {
                    "start": window.span.start,
                    "end": window.span.end,
                    "reverseComplement": window.span.reverse_complement,
                    "searchRegion": window.search_region,
                },
            })
    planned: list[LinearUnitPair] = []
    for number, pair in enumerate(pairs):
        if not isinstance(pair, Mapping) or pair.get("source") not in ("losat", "upload"):
            raise _protocol_error("A record pair needs record keys and a source, losat or upload.")
        query, subject = (
            units_by_record.get(_key(pair.get(side), "record pair")) for side in ("queryRecordKey", "subjectRecordKey")
        )
        if query is None or subject is None:
            raise _protocol_error("A record pair names a record the request does not draw.")
        if pair["source"] == "upload":
            refuse_file_comparison_of_windows(
                query, subject, records=(pair["queryRecordKey"], pair["subjectRecordKey"]), where="A comparison"
            )
        planned.extend(
            {"pair": number, "queryUnit": query_unit, "subjectUnit": subject_unit}
            for query_unit, subject_unit in comparison_unit_pairs(query, subject)
        )
    try:
        rows = drawn_unit_rows(unit_grid_rows, [unit["recordKey"] for unit in units])
    except ValidationError as exc:
        raise _protocol_error(str(exc)) from exc
    grid_pairs = grid_display_pairs(rows)
    return {
        "units": units,
        "pairs": planned,
        "gridPairs": None if grid_pairs is None else [
            {"queryUnit": query_unit, "subjectUnit": subject_unit} for query_unit, subject_unit in grid_pairs
        ],
    }


def plan_linear_units_json(records_json: str, pairs_json: str) -> str:
    """JSON boundary of the diagram Worker."""

    return json.dumps(plan_linear_units(json.loads(str(records_json)), json.loads(str(pairs_json))))


class GenomeScaleWarning(TypedDict):
    """A unit pair searched per file, its genome pair's search space past what LOSAT reads, and the warning the
    command line logs for it (Owner decision 2026-11-10)."""

    queryUnit: int
    subjectUnit: int
    message: str


class LinearGenomeScales(TypedDict):
    """The E-value scale of each search between two drawn units of a Linear drawing's LOSAT search: ``recordKeys``,
    the key of the record each drawn unit draws, in drawing order, and ``scales[q][s]``, the ``-searchsp`` (LOSATN,
    TLOSATX) or ``-dbsize`` (LOSATP) a search of unit ``q`` against unit ``s`` passes, null for LOSAT's own
    (PD-OI-018 revision 6). A scale is a decimal string, as the page passes it on unchanged: the product of two
    genomes of about 95 Mb is past what a JavaScript number holds exactly. ``warnings`` lists the unit pairs searched
    per file (:class:`GenomeScaleWarning`); the page shows the warnings of the pairs it searches."""

    recordKeys: list[str]
    scales: list[list[str | None]]
    warnings: list[GenomeScaleWarning]


def plan_genome_scales(
    request_payload: Mapping[str, object],
    resource_paths: Mapping[str, str],
    output_directory: str,
    program: object,
    resource_names: Mapping[str, object] | None = None,
) -> LinearGenomeScales:
    """The E-value scales of the ``program`` searches (``losatn``, ``tlosatx`` or ``losatp``) of a Linear request
    (``linear_genome_scales``). A request holds a LOSATP search but no nucleotide one, which the page resolves before
    it renders, so the page names the program. ``resource_names`` (``resourceManifest[].name``) names the files a
    warning names, as the input summary does."""

    request = decode_canonical_request(
        request_payload, resource_paths=resource_paths, output_directory=output_directory,
    )
    if not isinstance(request, LinearDiagramRequest) or program not in GENOME_SCALE_PROGRAMS:
        raise _protocol_error("Genome scales need a Linear request and a LOSAT program.")
    keys, scale, warning = linear_genome_scales(
        request,
        cast(GenomeScaleProgram, program),
        _staged_source_names(resource_names, {key: Path(path) for key, path in resource_paths.items()}),
    )
    units = range(len(keys))

    def decimal(query: int, subject: int) -> str | None:
        value = scale(query, subject)
        return None if value is None else str(value)

    return {
        "recordKeys": list(keys),
        "scales": [[decimal(query, subject) for subject in units] for query in units],
        "warnings": [
            {"queryUnit": query, "subjectUnit": subject, "message": message}
            for query in units for subject in units if (message := warning(query, subject)) is not None
        ],
    }


def plan_genome_scales_json(
    request_json: str, resource_paths_json: str, output_directory: str, program: str, resource_names_json: str
) -> str:
    """JSON boundary of the diagram Worker."""

    return json.dumps(plan_genome_scales(
        json.loads(str(request_json)), json.loads(str(resource_paths_json)), str(output_directory), str(program),
        json.loads(str(resource_names_json)),
    ))


class CardWindowInput(TypedDict):
    """A window as the request writes it (``renderRequest.records[].windows[]``)."""

    start: int
    end: int
    reverseComplement: bool
    column: int | None


class CardWindowRecordInput(TypedDict):
    """A record the request draws, as its card's window plan reads it: its key, its windows and its place
    (``gridRow``, ``gridColumn``, ``reverseComplement``) as the request writes them, the length of the record it
    reads (needed when it lists windows; a record drawn whole may leave it out, as no number or omitted length
    reads its span), how many records it draws (``drawnRecords``: a file card drawn whole draws each record of its
    file, one after another in its row, and one with windows draws one record), and whether it is drawn
    circular."""

    recordKey: str
    windows: list[CardWindowInput]
    gridRow: int | None
    gridColumn: int | None
    reverseComplement: bool
    sourceLength: int | None
    drawnRecords: int
    circular: bool


class CardWindowNumber(TypedDict):
    """A window of a record drawn in several and its circled number in the drawing (its place in its row)."""

    recordKey: str
    window: int
    number: int


class CardWindowGap(TypedDict):
    """The bases of a record left out between two of its windows that stand next to each other in their row and
    run in the same direction (``leftWindow`` drawn before ``rightWindow``): the length the drawing writes above the
    "//" between them."""

    recordKey: str
    leftWindow: int
    rightWindow: int
    omittedBp: int


class CardWindowPlan(TypedDict):
    """The numbers of the windows of the records drawn in several and the bases left out between them."""

    windows: list[CardWindowNumber]
    gaps: list[CardWindowGap]


def _card_record(record: object, position: int) -> tuple[CardWindowRecordInput, tuple[ResolvedWindow | None, ...]]:
    if not isinstance(record, Mapping):
        raise _protocol_error("A record of a card window plan must be an object.")
    record_key = _key(record.get("recordKey"), "record")
    length, row, column = record.get("sourceLength"), record.get("gridRow"), record.get("gridColumn")
    drawn = record.get("drawnRecords")
    if length is not None and (not isinstance(length, int) or isinstance(length, bool) or length < 1):
        raise _protocol_error(f"Record {record_key!r} has a length that is not a positive whole number.")
    if not isinstance(drawn, int) or isinstance(drawn, bool) or drawn < 1:
        raise _protocol_error(f"Record {record_key!r} needs how many records it draws.")
    if any(value is not None and (not isinstance(value, int) or isinstance(value, bool)) for value in (row, column)):
        raise _protocol_error(f"Record {record_key!r} has a grid row or column that is not a whole number.")
    if not isinstance(record.get("reverseComplement"), bool) or not isinstance(record.get("circular"), bool):
        raise _protocol_error(f"Record {record_key!r} needs reverseComplement and circular.")
    try:
        windows = _decode_windows(record.get("windows"), path=f"records[{position}].windows")
    except (CanonicalRequestDecodingError, ValueError) as exc:
        raise _protocol_error(str(exc)) from exc
    if windows and (length is None or drawn != 1):
        raise _protocol_error(f"Record {record_key!r} lists windows: it needs its length and draws one record.")
    return cast(CardWindowRecordInput, record), tuple(record_drawn_windows(record_key, windows))


def plan_card_windows(records: Sequence[object]) -> CardWindowPlan:
    """The circled numbers and omitted lengths of the windows of ``records`` (in request order) as the drawing writes
    them: the drawn units in each row in drawing order (``grid_placement_order``; without grid rows each record is a
    row of its own), read by the row's captions (``plan_row_windows``)."""

    units: list[tuple[CardWindowRecordInput, int | None, RecordUnit, RecordPresentation]] = []
    cards = [_card_record(record, position) for position, record in enumerate(records)]
    # Each card at the place the drawing gives it: a placed file card drawn whole takes a column per record.
    places = numbered_places([
        InputPlace.at(
            card["gridRow"], card["gridColumn"], [window.span.column for window in windows if window is not None],
            card["drawnRecords"],
        )
        for card, windows in cards
    ])
    for (card, windows), place in zip(cards, places, strict=True):
        # A file card drawn whole draws each of its records as a unit of its own; a record drawn whole without a known
        # length takes a one-base span, which no number or omitted length reads.
        many = card["drawnRecords"] > 1
        drawn = (
            [
                (expanded_record_key(card["recordKey"], index + 1) if many else card["recordKey"], None,
                 place.record_column(index))
                for index in range(card["drawnRecords"])
            ]
            if windows == (None,)
            else [
                (card["recordKey"], window, place.window_columns[at] if place.window_columns else place.column)
                for at, window in enumerate(windows)
            ]
        )
        for identity, window, column in drawn:
            length = card["sourceLength"] or 1
            span = (1, length) if window is None else (window.span.start, window.span.end)
            reverse = window_drawn_reverse(window is not None and window.span.reverse_complement, card["reverseComplement"])
            unit = RecordUnit(
                record_index=len(units), record_key=identity, record_id=identity,
                span=span, reverse=reverse, source_length=length, circular=card["circular"],
                is_window=window is not None, window_index=None if window is None or len(windows) < 2 else window.index,
            )
            placement = RecordPresentation(grid_row=card["gridRow"], grid_column=column)
            units.append((card, None if window is None else window.index, unit, placement))
    positioned = [(index, placement) for index, (*_, placement) in enumerate(units) if placement.grid_row is not None]
    if positioned:
        ordered = grid_placement_order(positioned)
        rows: dict[object, list[int]] = {}
        for index, placement in ordered:
            rows.setdefault(placement.grid_row, []).append(index)
        unplaced = [[index] for index, (*_, placement) in enumerate(units) if placement.grid_row is None]
        row_units = [*rows.values(), *unplaced]
    else:
        by_record: dict[str, list[int]] = {}
        for index, (_card, _window, unit, _placement) in enumerate(units):
            by_record.setdefault(unit.record_key, []).append(index)
        row_units = list(by_record.values())
    numbers: list[CardWindowNumber] = []
    gaps: list[CardWindowGap] = []
    for row in row_units:
        plan = plan_row_windows([units[index][2] for index in row], [""] * len(row))
        if plan is None:
            continue
        for index in row:
            card, window_index, _unit, _placement = units[index]
            if window_index is not None and len(card["windows"]) > 1:
                numbers.append({"recordKey": card["recordKey"], "window": window_index, "number": plan.number_of(index)})
        for gap in plan.gaps:
            left, right = units[gap.left_record_index], units[gap.right_record_index]
            if gap.omitted_bp is not None and left[1] is not None and right[1] is not None:
                gaps.append({
                    "recordKey": left[0]["recordKey"], "leftWindow": left[1], "rightWindow": right[1],
                    "omittedBp": gap.omitted_bp,
                })
    return {"windows": numbers, "gaps": gaps}


def plan_card_windows_json(records_json: str) -> str:
    """JSON boundary of the diagram Worker."""

    return json.dumps(plan_card_windows(json.loads(str(records_json))))


class RecordPlaceInput(TypedDict):
    """Where a request record stands before its records expand (``renderRequest.records[]``): its grid row and
    column, the column of each of its windows, and how many records it draws (a file drawn whole: each record of the
    file; any other record: one)."""

    gridRow: int | None
    gridColumn: int | None
    windowColumns: list[int | None]
    drawnRecords: int


class RecordPlace(TypedDict):
    """Where the drawing places a request record: the column of each record or window it draws, in drawing order
    (``InputPlace.drawn_columns`` after ``numbered_places``); its row stays the one the request writes."""

    columns: list[int | None]


def _positive_whole_number(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 1


def plan_record_places(records: Sequence[object]) -> list[RecordPlace]:
    """The columns the drawing gives the records of ``records`` (in request order): a placed file drawn whole takes
    one column per record and moves the later columns of its row (``numbered_places``). The Source recipe writes a
    records table row per drawn record at these columns."""

    places: list[InputPlace] = []
    for record in records:
        if not isinstance(record, Mapping):
            raise _protocol_error("A record of a place plan must be an object.")
        row, column, windows, drawn = (
            record.get("gridRow"), record.get("gridColumn"), record.get("windowColumns"), record.get("drawnRecords")
        )
        if (
            any(value is not None and not _positive_whole_number(value) for value in (row, column))
            or not isinstance(windows, list)
            or any(value is not None and not _positive_whole_number(value) for value in windows)
            or not _positive_whole_number(drawn)
            or (windows and drawn != 1)
        ):
            raise _protocol_error("A record of a place plan needs whole-number places and how many records it draws.")
        places.append(InputPlace.at(cast(int | None, row), cast(int | None, column), windows, cast(int, drawn)))
    return [{"columns": list(place.drawn_columns())} for place in numbered_places(places)]


def plan_record_places_json(records_json: str) -> str:
    """JSON boundary of the diagram Worker."""

    return json.dumps(plan_record_places(json.loads(str(records_json))))


def describe_card_window_plan() -> dict[str, dict[str, Shape]]:
    """The requests and the answers of :func:`plan_card_windows_json` and :func:`plan_record_places_json` as named
    shapes of :mod:`gbdraw.contract_shapes` (tools/generate_web_contracts.py)."""

    types_: dict[str, Shape] = {}

    def special(hint: object) -> Shape | None:
        if isinstance(hint, type) and is_typeddict(hint):
            if hint.__name__ not in types_:
                types_[hint.__name__] = {"kind": "unknown"}  # reserves its place before its members'
                types_[hint.__name__] = typed_dict_shape(hint, special)
            return ref(hint.__name__)
        return None

    for root in (CardWindowRecordInput, CardWindowPlan, RecordPlaceInput, RecordPlace):
        special(root)
    return {"types": types_}


__all__ = [
    "CardWindowGap",
    "CardWindowInput",
    "CardWindowNumber",
    "CardWindowPlan",
    "CardWindowRecordInput",
    "describe_card_window_plan",
    "GenomeScaleWarning",
    "LinearGenomeScales",
    "plan_genome_scales",
    "plan_genome_scales_json",
    "plan_card_windows",
    "plan_card_windows_json",
    "LinearDrawnUnit",
    "LinearUnitGridPair",
    "LinearUnitPair",
    "LinearUnitPairInput",
    "LinearUnitPlan",
    "LinearUnitRecordInput",
    "LinearUnitWindow",
    "plan_linear_units",
    "plan_linear_units_json",
    "RecordPlace",
    "RecordPlaceInput",
    "plan_record_places",
    "plan_record_places_json",
]
