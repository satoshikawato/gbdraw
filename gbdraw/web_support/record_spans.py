"""Worker helper that reads spans of the records a request draws (Copy FASTA of Select region, PD-OI-115).

The page sends the request the shown Result was made from (the committed
request, with its staged resources) and spans of its records by record key, in
source coordinates, 1-based and inclusive. Python resolves the records as the
render does and slices each record's source sequence, on the source strand
(a reverse-complemented display gives the same bases as the forward one); a
span whose start is after its end crosses the origin of a circular record and
reads the end of the record, then its start. The reply holds each span's
sequence and header and the FASTA text of all of them. It renders nothing.

A sequence search reads its subjects here as well, naming the drawn records by
key alone (``DrawnRecordInput``): Python reads the bases each record draws,
from the transform it is drawn by, and the origin window of each record a span
may cross the origin of.
The headers of one read are its LOSAT subject IDs, so they are unique: spans
that would share a header (two records of one ID and length) name their record
keys in it.
"""

from __future__ import annotations

import json
from collections import Counter
from collections.abc import Mapping, Sequence
from typing import TypedDict, final

from gbdraw.annotations.search_hits import origin_window_span
from gbdraw.api.request_render import load_render_source_records, resolve_request_records
from gbdraw.comparisons.losat_jobs import losat_fasta_text
from gbdraw.exceptions import ValidationError
from gbdraw.core.record_metadata import _record_source_start, _source_period, _wrap_source_base
from gbdraw.layout.record_transform import spans_may_cross_origin
from gbdraw.session_request_codec import decode_canonical_request


@final  # closed, so a key test tells the two inputs apart
class RecordSpanInput(TypedDict):
    """One span: a record key of the request and source coordinates, 1-based and inclusive."""

    recordKey: str
    start: int
    end: int


@final
class DrawnRecordInput(TypedDict):
    """A record the request draws, read over the bases it draws (from ``_record_source_start``, the source base of
    its first base on the source strand that a crop records, for as many bases as it has); and, when a span of it may
    cross its origin (``spans_may_cross_origin``, the axis source map's ``isCircular``), over its origin window
    (``origin_window_span``): its last and first bases, as many as ``originWindow`` (a sequence search's longest
    query) and at most half the record each. A record asked twice is read once; the origin windows follow every
    drawn span."""

    recordKey: str
    originWindow: int


class RecordSpanSequence(TypedDict):
    """One span's bases on the source strand, with the FASTA header that names them; ``originWindow`` when it was
    read as an origin window (``DrawnRecordInput``)."""

    recordKey: str
    recordId: str
    start: int
    end: int
    length: int
    wrapsOrigin: bool
    originWindow: bool
    recordLength: int
    header: str
    sequence: str


class RecordSpanSequencesReply(TypedDict):
    """The spans in the order asked (the origin windows of drawn records last), and their FASTA text."""

    spans: list[RecordSpanSequence]
    fasta: str


def _protocol_error(message: str) -> ValidationError:
    return ValidationError(message, diagnostic={"code": "HELPER_PROTOCOL"})


def _span_inputs(value: object) -> list[RecordSpanInput | DrawnRecordInput]:
    if not isinstance(value, list) or not value:
        raise _protocol_error("A record span read needs a list of spans.")
    spans: list[RecordSpanInput | DrawnRecordInput] = []
    for item in value:
        if not (isinstance(item, Mapping) and isinstance(item.get("recordKey"), str)):
            raise _protocol_error("A span needs a record key.")
        if set(item) == {"recordKey", "originWindow"} and type(item["originWindow"]) is int and item["originWindow"] >= 1:
            spans.append({"recordKey": item["recordKey"], "originWindow": item["originWindow"]})
        elif all(type(item.get(key)) is int for key in ("start", "end")):
            spans.append({"recordKey": item["recordKey"], "start": item["start"], "end": item["end"]})
        else:
            raise _protocol_error("A span needs a start and an end, or, for a drawn record, the reach of its origin window.")
    # One read is spans by coordinates (Select region) or drawn records (a search), so the page can check its reply.
    if len({"originWindow" in span for span in spans}) > 1:
        raise _protocol_error("A record span read asks for spans by coordinates or for drawn records, not both.")
    return spans


def span_header(record_id: str, start: int, end: int) -> str:
    """The FASTA header of a span: the record ID and its source coordinates (``start`` after ``end`` across the origin)."""

    return f"{record_id}:{start}-{end}"


def read_record_span_sequences(
    request_payload: Mapping[str, object],
    resource_paths: Mapping[str, str],
    output_directory: str,
    spans: Sequence[RecordSpanInput | DrawnRecordInput],
) -> RecordSpanSequencesReply:
    """The bases of ``spans`` of the records that the canonical request ``request_payload`` draws."""

    asked = _span_inputs(list(spans))
    request = decode_canonical_request(request_payload, resource_paths=resource_paths, output_directory=output_directory)
    collection = resolve_request_records(request)
    # A record drawn through windows (D-01) is named by its own key and draws one unit per window.
    by_key: dict[str, list[int]] = {}
    for index, item in enumerate(collection.provenance):
        by_key.setdefault(item.record_identity_key, []).append(index)
    # (record index, start, end, origin window) in reading order: the drawn records once each, their windows last.
    reads: list[tuple[int, int, int, bool]] = []
    windows: list[tuple[int, int]] = []
    drawn: set[str] = set()
    for span in asked:
        units = by_key.get(span["recordKey"])
        if units is None:
            raise ValidationError(
                f"The drawing has no record {span['recordKey']!r}.",
                diagnostic={"code": "INPUT_INVALID", "field": "record_selector"},
            )
        if "originWindow" in span:  # the key test narrows the TypedDict
            if span["recordKey"] in drawn:
                continue
            drawn.add(span["recordKey"])
            for index in units:
                record = collection.records[index]
                start = _record_source_start(record)
                # A window across the origin ends at a base of the record's start (D-28).
                reads.append((index, start, _wrap_source_base(start + len(record.seq) - 1, _source_period(record)), False))
            if any(spans_may_cross_origin(collection.records[index], collection.transforms[index]) for index in units):
                windows.append((units[0], span["originWindow"]))
        else:
            # Spans are in source coordinates, which every window of a record shares.
            reads.append((units[0], span["start"], span["end"], False))
    sources: dict[int, str] = {}

    def source_bases(index: int) -> str:
        if index not in sources:
            item = collection.provenance[index]
            source = load_render_source_records(request.records[item.input_index].source)
            sources[index] = str(source.records[item.source_record_index].seq)
        return sources[index]

    reads.extend((index, *origin_window_span(len(source_bases(index)), reach), True) for index, reach in windows)
    read: list[RecordSpanSequence] = []
    for index, start, end, window in reads:
        item = collection.provenance[index]
        bases = source_bases(index)
        length = len(bases)
        if not (1 <= start <= length and 1 <= end <= length) or (start > end and not item.is_circular):
            raise ValidationError(
                f"Span {start}..{end} of {item.record_key!r} is outside 1..{length}.",
                diagnostic={"code": "REGION_INVALID", "field": "region", "reason": "RECORD_BOUNDS"},
            )
        wraps = start > end
        sequence = bases[start - 1:] + bases[:end] if wraps else bases[start - 1:end]
        read.append({
            "recordKey": item.record_identity_key,
            "recordId": item.source_record_id,
            "start": start,
            "end": end,
            "length": len(sequence),
            "wrapsOrigin": wraps,
            "originWindow": window,
            "recordLength": length,
            "header": span_header(item.source_record_id, start, end),
            "sequence": sequence,
        })
    shared = {header for header, count in Counter(span["header"] for span in read).items() if count > 1}
    for span in read:
        if span["header"] in shared:
            span["header"] = f"{span['header']}@{span['recordKey']}"
    return {"spans": read, "fasta": "".join(losat_fasta_text(span["header"], span["sequence"]) for span in read)}


def read_record_span_sequences_json(
    request_json: str,
    resource_paths_json: str,
    output_directory: str,
    spans_json: str,
) -> str:
    """JSON boundary of the diagram Worker."""

    return json.dumps(
        read_record_span_sequences(
            json.loads(str(request_json)),
            json.loads(str(resource_paths_json)),
            str(output_directory),
            json.loads(str(spans_json)),
        )
    )


__all__ = [
    "DrawnRecordInput",
    "RecordSpanInput",
    "RecordSpanSequence",
    "RecordSpanSequencesReply",
    "read_record_span_sequences",
    "read_record_span_sequences_json",
    "span_header",
]
