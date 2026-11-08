"""Worker helper that plans the drawn units of a Linear request for the Web's comparison searches (D-01, D-37).

Before its comparison searches run, a Generate sends the records its request will draw (each record key with its
windows, as the request writes them) and the record pairs its comparison plan compares. Python answers with the
drawn units, in drawing order, and with the unit pairs of each record pair: every window of one record against every
window of the other, the pairs a CLI search runs (``comparison_unit_pairs``). The page runs its searches over these
units and writes their indexes into the request, which Python reads as drawn-unit indexes. A comparison file
(an uploaded BLAST TSV) gives the coordinates of one drawn unit, so a pair of a record drawn in several windows must
be searched.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from typing import Literal, TypedDict

from gbdraw.api.record_planning import comparison_unit_pairs, record_drawn_windows, refuse_file_comparison_of_windows
from gbdraw.exceptions import ValidationError
from gbdraw.session_request_codec import CanonicalRequestDecodingError, _decode_windows


class LinearUnitRecordInput(TypedDict):
    """A record the request draws: its key and its windows (``renderRequest.records[].windows``)."""

    recordKey: str
    windows: list[object]


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


class LinearUnitPlan(TypedDict):
    """The drawn units in drawing order and the unit pairs, in the order of the record pairs."""

    units: list[LinearDrawnUnit]
    pairs: list[LinearUnitPair]


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
    for position, record in enumerate(records):
        if not isinstance(record, Mapping):
            raise _protocol_error("A record of a unit plan must be an object.")
        record_key = _key(record.get("recordKey"), "record")
        if record_key in units_by_record:
            raise _protocol_error(f"Record {record_key!r} is named twice.")
        try:
            windows = _decode_windows(record.get("windows"), path=f"records[{position}].windows")
        except (CanonicalRequestDecodingError, ValueError) as exc:
            raise _protocol_error(str(exc)) from exc
        indexes = units_by_record[record_key] = []
        for window in record_drawn_windows(record_key, windows):
            indexes.append(len(units))
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
    return {"units": units, "pairs": planned}


def plan_linear_units_json(records_json: str, pairs_json: str) -> str:
    """JSON boundary of the diagram Worker."""

    return json.dumps(plan_linear_units(json.loads(str(records_json)), json.loads(str(pairs_json))))


__all__ = [
    "LinearDrawnUnit",
    "LinearUnitPair",
    "LinearUnitPairInput",
    "LinearUnitPlan",
    "LinearUnitRecordInput",
    "LinearUnitWindow",
    "plan_linear_units",
    "plan_linear_units_json",
]
