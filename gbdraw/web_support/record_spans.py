"""Worker helper that reads spans of the records a request draws (Copy FASTA of Select region, PD-OI-115).

The page sends the request the shown Result was made from (the committed
request, with its staged resources) and spans of its records by record key, in
source coordinates, 1-based and inclusive. Python resolves the records as the
render does and slices each record's source sequence, on the source strand
(a reverse-complemented display gives the same bases as the forward one); a
span whose start is after its end crosses the origin of a circular record and
reads the end of the record, then its start. The reply holds each span's
sequence and header and the FASTA text of all of them. It renders nothing.

A sequence search over drawn records can read its subjects here as well.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from typing import TypedDict

from gbdraw.api.request_render import load_render_source_records, resolve_request_records
from gbdraw.comparisons.losat_jobs import losat_fasta_text
from gbdraw.exceptions import ValidationError
from gbdraw.session_request_codec import decode_canonical_request


class RecordSpanInput(TypedDict):
    """One span: a record key of the request and source coordinates, 1-based and inclusive."""

    recordKey: str
    start: int
    end: int


class RecordSpanSequence(TypedDict):
    """One span's bases on the source strand, with the FASTA header that names them."""

    recordKey: str
    recordId: str
    start: int
    end: int
    length: int
    wrapsOrigin: bool
    header: str
    sequence: str


class RecordSpanSequencesReply(TypedDict):
    """The spans in the order asked, and their FASTA text."""

    spans: list[RecordSpanSequence]
    fasta: str


def _protocol_error(message: str) -> ValidationError:
    return ValidationError(message, diagnostic={"code": "HELPER_PROTOCOL"})


def _span_inputs(value: object) -> list[RecordSpanInput]:
    if not isinstance(value, list) or not value:
        raise _protocol_error("A record span read needs a list of spans.")
    spans: list[RecordSpanInput] = []
    for item in value:
        if not (
            isinstance(item, Mapping) and isinstance(item.get("recordKey"), str)
            and all(type(item.get(key)) is int for key in ("start", "end"))
        ):
            raise _protocol_error("A span needs a record key, a start and an end.")
        spans.append({"recordKey": item["recordKey"], "start": item["start"], "end": item["end"]})
    return spans


def span_header(record_id: str, start: int, end: int) -> str:
    """The FASTA header of a span: the record ID and its source coordinates (``start`` after ``end`` across the origin)."""

    return f"{record_id}:{start}-{end}"


def read_record_span_sequences(
    request_payload: Mapping[str, object],
    resource_paths: Mapping[str, str],
    output_directory: str,
    spans: Sequence[RecordSpanInput],
) -> RecordSpanSequencesReply:
    """The bases of ``spans`` of the records that the canonical request ``request_payload`` draws."""

    asked = _span_inputs(list(spans))
    request = decode_canonical_request(request_payload, resource_paths=resource_paths, output_directory=output_directory)
    collection = resolve_request_records(request)
    by_key = {item.record_key: item for item in collection.provenance}
    read: list[RecordSpanSequence] = []
    for span in asked:
        item = by_key.get(span["recordKey"])
        if item is None:
            raise ValidationError(
                f"The drawing has no record {span['recordKey']!r}.",
                diagnostic={"code": "INPUT_INVALID", "field": "record_selector"},
            )
        record = load_render_source_records(request.records[item.input_index].source).records[item.source_record_index]
        bases = str(record.seq)
        length = len(bases)
        start, end = span["start"], span["end"]
        if not (1 <= start <= length and 1 <= end <= length) or (start > end and not item.is_circular):
            raise ValidationError(
                f"Span {start}..{end} of {span['recordKey']!r} is outside 1..{length}.",
                diagnostic={"code": "REGION_INVALID", "field": "region", "reason": "RECORD_BOUNDS"},
            )
        wraps = start > end
        sequence = bases[start - 1:] + bases[:end] if wraps else bases[start - 1:end]
        read.append({
            "recordKey": item.record_key,
            "recordId": item.source_record_id,
            "start": start,
            "end": end,
            "length": len(sequence),
            "wrapsOrigin": wraps,
            "header": span_header(item.source_record_id, start, end),
            "sequence": sequence,
        })
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
    "RecordSpanInput",
    "RecordSpanSequence",
    "RecordSpanSequencesReply",
    "read_record_span_sequences",
    "read_record_span_sequences_json",
    "span_header",
]
