"""Worker helper that turns one Sequence search's LOSAT rows into an annotation set (PD-OI-117).

The page sends the raw outfmt 6 text LOSAT wrote, the sequences it searched
(each with the ID LOSAT names it by, the drawing's record selector, the
request's record key, and the source coordinate of its first base), the
panel's thresholds, the queries in order with their lengths and programs, and the colour mode.
Python reads the table with ``read_comparison_table`` and builds the set with
``search_hit_annotation_set``; the reply is data (the set as the canonical
request writes it, the hit table, and the threshold tally). It renders nothing.

The queries themselves are read here too (``read_search_queries``), with the
FASTA reader of search queries (``load_query_fasta``), which tells a nucleotide
query from a protein query; the program that searches each is decided here
(``search_program``), so the Web has no sequence reader or program rule of its own.
"""

from __future__ import annotations

import io
import json
from collections import Counter
from collections.abc import Mapping
from typing import Literal, TypedDict, cast

from gbdraw.annotations.search_hits import (
    SearchColorMode,
    SearchProgram,
    SearchQuery,
    SearchSubject,
    search_hit_annotation_set,
    search_program,
)
from gbdraw.comparisons.losat_jobs import losat_fasta_text
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import read_comparison_table
from gbdraw.io.genome import load_query_fasta
from gbdraw.mode_profiles import ComparisonThresholds
from gbdraw.session_request_codec import decode_record_selector, encode_annotation_set


class SearchSubjectInput(TypedDict):
    """One searched sequence: its LOSAT ID, the record selector (``recordId``/``recordIndex``), its record key,
    the source coordinate of its first base, and, as the span reader read it (``RecordSpanSequence``), whether it
    is an origin window (``origin_window_span``) and its record's length."""

    subjectId: str
    record: dict[str, object]
    recordKey: str
    sourceStart: int
    originWindow: bool
    recordLength: int


class SearchThresholdsInput(TypedDict):
    evalue: float
    bitscore: float
    identity: float
    alignmentLength: int


class SearchQueryInput(TypedDict):
    """One query: its FASTA ID, its length (bases, or residues for a protein query) and the program that searched it."""

    name: str
    length: int
    program: SearchProgram


class SequenceSearchInput(TypedDict):
    table: str
    setId: str
    subjects: list[SearchSubjectInput]
    thresholds: SearchThresholdsInput
    queries: list[SearchQueryInput]
    colorMode: SearchColorMode
    queryColors: dict[str, str]


class SearchHitReply(TypedDict):
    annotationId: str
    query: str
    recordKey: str
    start: int
    end: int
    strand: Literal["+", "-"]
    identity: float
    alignmentLength: int
    evalue: float
    bitscore: float
    qstart: int
    qend: int


class ThresholdTallyReply(TypedDict):
    name: str
    threshold: float
    rowsLeftOut: int
    bestValue: float


class SequenceSearchReply(TypedDict):
    """The set (null when no hit passed the thresholds), its hits, and what the thresholds and IDs left out."""

    annotationSet: dict[str, object] | None
    hits: list[SearchHitReply]
    rowsRead: int
    rowsKept: int
    thresholds: list[ThresholdTallyReply]
    unmatchedRows: int
    unmatchedSubjectIds: list[str]


class SearchQueryRead(TypedDict):
    """One query as read: its FASTA ID, its length (bases, or residues for a protein query), the program that
    searches it, the most subject bases one of its hits covers (``SearchQuery.reach``, which places a record's origin
    window), and the FASTA text LOSAT reads."""

    name: str
    length: int
    program: SearchProgram
    reach: int
    fasta: str


class SearchQueriesReply(TypedDict):
    queries: list[SearchQueryRead]


# The queries name their hits and their colours, so two queries may not share an ID.
_DUPLICATE_QUERY_IDS = {"code": "INPUT_INVALID", "field": "fasta", "reason": "UNIQUE_IDS"}


def read_search_queries(path: str) -> SearchQueriesReply:
    """The queries of a FASTA file, each with its own ID: a nucleotide query is searched with blastn, a protein
    query with tblastn.

    The query reader (``load_query_fasta``) refuses a file without an entry, a
    record without sequence and a character that is neither a nucleotide nor an
    amino acid code, each with its diagnostic.
    """

    records = load_query_fasta(path)
    repeated = sorted(name for name, count in Counter(str(record.id) for record, _ in records).items() if count > 1)
    if repeated:
        raise ValidationError(
            f"Query IDs repeat: {', '.join(repeated)}. Give each query its own ID.", diagnostic=_DUPLICATE_QUERY_IDS
        )
    queries = [(record, SearchQuery(str(record.id), len(record.seq), search_program(molecule))) for record, molecule in records]
    return {
        "queries": [
            {
                "name": query.name, "length": query.length, "program": query.program, "reach": query.reach,
                "fasta": losat_fasta_text(query.name, str(record.seq)),
            }
            for record, query in queries
        ]
    }


def read_search_queries_json(path: str) -> str:
    """JSON boundary of the diagram Worker."""

    return json.dumps(read_search_queries(path))


_FIELDS = ("table", "setId", "subjects", "thresholds", "queries", "colorMode", "queryColors")


def _protocol_error(message: str) -> ValidationError:
    return ValidationError(message, diagnostic={"code": "HELPER_PROTOCOL", "reason": "FIELDS"})


def _subject(value: object, index: int) -> SearchSubject:
    if not (
        isinstance(value, Mapping) and isinstance(value.get("subjectId"), str)
        and isinstance(value.get("recordKey"), str) and type(value.get("sourceStart")) is int
        and type(value.get("originWindow")) is bool and type(value.get("recordLength")) is int
    ):
        raise _protocol_error(
            "A search subject needs an ID, a record, a record key, a source start, its origin window flag and its record's length."
        )
    record = decode_record_selector(value.get("record"), path=f"subjects[{index}].record")
    if record is None:
        raise _protocol_error("A search subject needs a record.")
    return SearchSubject(
        value["subjectId"], record, value["recordKey"], value["sourceStart"],
        record_length=value["recordLength"], origin_window=value["originWindow"],
    )


def _query(value: object) -> SearchQuery:
    if not (isinstance(value, Mapping) and set(value) == {"name", "length", "program"} and isinstance(value["name"], str)):
        raise _protocol_error("A search query has a name, a length and a program.")
    # SearchQuery checks the length and the program.
    return SearchQuery(value["name"], cast(int, value["length"]), cast(SearchProgram, value["program"]))


def sequence_search_annotation_set(payload: Mapping[str, object]) -> SequenceSearchReply:
    """The annotation set of one search, from the page's ``payload`` (``SequenceSearchInput``)."""

    if not isinstance(payload, Mapping) or set(payload) != set(_FIELDS):
        raise _protocol_error("A sequence search payload has the fields " + ", ".join(_FIELDS) + ".")
    subjects, thresholds, queries, colors = (
        payload["subjects"], payload["thresholds"], payload["queries"], payload["queryColors"]
    )
    if not (
        isinstance(payload["table"], str) and isinstance(payload["setId"], str) and isinstance(subjects, list)
        and isinstance(thresholds, Mapping) and isinstance(queries, list) and isinstance(colors, Mapping)
        and all(isinstance(key, str) and isinstance(color, str) for key, color in colors.items())
    ):
        raise _protocol_error("A sequence search payload has a value of the wrong type.")
    result = search_hit_annotation_set(
        read_comparison_table(io.StringIO(payload["table"]), label="sequence search"),
        set_id=payload["setId"],
        subjects=[_subject(item, index) for index, item in enumerate(subjects)],
        # The casts only satisfy the checker: ComparisonThresholds and search_hit_annotation_set check each value.
        thresholds=ComparisonThresholds(
            evalue=cast(float, thresholds.get("evalue")),
            bitscore=cast(float, thresholds.get("bitscore")),
            identity=cast(float, thresholds.get("identity")),
            alignment_length=cast(int, thresholds.get("alignmentLength")),
        ),
        queries=[_query(item) for item in queries],
        color_mode=cast(SearchColorMode, payload["colorMode"]),
        query_colors=dict(colors),
    )
    return {
        "annotationSet": None if result.annotation_set is None else encode_annotation_set(result.annotation_set),
        "hits": [
            {
                "annotationId": hit.annotation_id, "query": hit.query, "recordKey": hit.record_key,
                "start": hit.start, "end": hit.end, "strand": hit.strand, "identity": hit.identity,
                "alignmentLength": hit.alignment_length, "evalue": hit.evalue, "bitscore": hit.bitscore,
                "qstart": hit.qstart, "qend": hit.qend,
            }
            for hit in result.hits
        ],
        "rowsRead": result.tally.rows_read,
        "rowsKept": result.tally.rows_kept,
        "thresholds": [
            {"name": item.name, "threshold": item.threshold, "rowsLeftOut": item.rows_left_out, "bestValue": item.best_value}
            for item in result.tally.thresholds
        ],
        "unmatchedRows": result.unmatched_rows,
        "unmatchedSubjectIds": list(result.unmatched_subject_ids),
    }


def sequence_search_annotation_set_json(payload_json: str) -> str:
    """JSON boundary of the diagram Worker."""

    return json.dumps(sequence_search_annotation_set(json.loads(str(payload_json))))


__all__ = [
    "SearchHitReply",
    "SearchQueriesReply",
    "SearchQueryRead",
    "SearchQueryInput",
    "SearchSubjectInput",
    "SearchThresholdsInput",
    "SequenceSearchInput",
    "SequenceSearchReply",
    "ThresholdTallyReply",
    "read_search_queries",
    "read_search_queries_json",
    "sequence_search_annotation_set",
    "sequence_search_annotation_set_json",
]
