"""Sequence-search hits as one annotation set (PD-OI-117, the Web's Sequence search).

A search runs its queries with LOSAT (blastn for a nucleotide query, tblastn
for a protein query, ``search_program``) against sequences of the records a
drawing draws. Its outfmt 6 rows, read by
``read_comparison_table``, are filtered with the comparison thresholds
(``filter_comparison_dataframe``, the rule of Pairwise Comparisons) and become
one ``highlight`` annotation per hit, in source coordinates of the record the
row's subject names. The hit colour is one colour per query, or a gradient by
query position whose Legend row is ``QUERY_POSITION_LEGEND_LABEL``.

A circular record drawn whole is searched as its sequence from base 1 and over
its origin window (``origin_window_span``), so a hit across the origin is one
row with ``wraps_origin`` and is not split at the seam (OV-2571).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

from pandas import DataFrame

from gbdraw.exceptions import ValidationError
from gbdraw.io.genome import SequenceMolecule
from gbdraw.io.comparisons import (
    COMPARISON_COLUMNS,
    ComparisonFilter,
    ComparisonRowTally,
    comparison_row_tally,
    filter_comparison_dataframe,
    normalize_comparison_dataframe,
)
from gbdraw.io.record_select import RecordSelector
from gbdraw.linear_comparison import table_id_index, table_id_matches

from .models import AnnotationGradient, AnnotationSet, CoordinateSpan, RegionAnnotation, RegionAnnotationStyle

SearchProgram = Literal["blastn", "tblastn"]
# The program that searches a query of each molecule (PD-OI-117), and the subject bases one query letter covers.
_SEARCH_PROGRAMS: Mapping[SequenceMolecule, SearchProgram] = {"nucleotide": "blastn", "protein": "tblastn"}
_BASES_PER_QUERY_LETTER: Mapping[SearchProgram, int] = {"blastn": 1, "tblastn": 3}
SearchColorMode = Literal["query", "gradient"]

# The default hit colours, one per query in query order (repeated past the last), and the
# gradient by query position from its pale start to its dark end. Moderate chroma: a hit is
# a large area on a short record. The Web reads them from utils/sequence-search-colors.generated.js.
SEARCH_QUERY_COLORS: tuple[str, ...] = (
    "#c8607e", "#8670c4", "#3f9688", "#c08a3e", "#4f7fb8", "#7f9a4a", "#a0705a", "#6b7a8f",
)
SEARCH_GRADIENT_START_COLOR = "#f6c9d8"
SEARCH_GRADIENT_END_COLOR = "#5b1a4a"
# The Legend row of every gradient hit (one bar from "start" to "end").
QUERY_POSITION_LEGEND_LABEL = "Query position"
# Every hit covers the features under it at this fill opacity, a gradient hit too (receipt web.sequence-search).
SEARCH_HIT_FILL_OPACITY = 0.85
SEARCH_GRADIENT_STROKE_WIDTH = 0.5
# The keys of each hit's annotation metadata (D-02), as the comparison columns name them.
SEARCH_HIT_METADATA_KEYS = ("query", "program", "identity", "alignment_length", "evalue", "bitscore")


@dataclass(frozen=True)
class SearchSubject:
    """One sequence a search ran against.

    ``subject_id`` is the ID the LOSAT table names it by (its FASTA header),
    ``record`` the drawing's record the annotation names, ``record_key`` the
    request's identity of that record, and ``source_start`` the source
    coordinate of the sequence's first base (1, or the start of a crop).
    ``record_length`` is the record's length, after which a sequence read
    across the origin goes on at 1. ``origin_window`` marks the origin window
    of a record (``origin_window_span``: its last bases, then its first), whose
    rows that do not cross the origin are rows of the record's own sequence.
    """

    subject_id: str
    record: RecordSelector
    record_key: str
    source_start: int = 1
    record_length: int | None = None
    origin_window: bool = False

    def __post_init__(self) -> None:
        if not str(self.subject_id).strip() or not str(self.record_key).strip():
            raise ValidationError("A search subject needs an ID and a record key.", diagnostic={"code": "INPUT_INVALID", "reason": "REQUIRED"})
        if not isinstance(self.record, RecordSelector):
            raise ValidationError("A search subject's record must be a RecordSelector.", diagnostic={"code": "INPUT_INVALID", "reason": "FIELDS"})
        if isinstance(self.source_start, bool) or not isinstance(self.source_start, int) or self.source_start < 1:
            raise ValidationError(
                "A search subject's source start must be a positive integer.",
                diagnostic={"code": "INPUT_INVALID", "reason": "POSITIVE_INTEGER"},
            )
        length = self.record_length
        if length is not None and (
            isinstance(length, bool) or not isinstance(length, int) or not self.source_start <= length
        ):
            raise ValidationError(
                "A search subject starts within its record's length.",
                diagnostic={"code": "INPUT_INVALID", "reason": "FIELDS"},
            )
        if type(self.origin_window) is not bool or (self.origin_window and (length is None or self.source_start == 1)):
            raise ValidationError(
                "An origin window has its record's length and starts after base 1.",
                diagnostic={"code": "INPUT_INVALID", "reason": "FIELDS"},
            )

    def source_span(self, low: int, high: int) -> tuple[int, int, bool]:
        """The source span of the subject's bases ``low``..``high`` (1-based), and whether it crosses the origin
        (then the start is after the end)."""

        start, end = low + self.source_start - 1, high + self.source_start - 1
        length = self.record_length
        if length is None:
            return start, end, False
        return (start - 1) % length + 1, (end - 1) % length + 1, start <= length < end


def origin_window_span(record_length: int, reach: int) -> tuple[int, int]:
    """The origin window of a circular record, as a source span across its origin: its last ``w`` bases, then its
    first ``w``, with ``w`` the bases the longest query can cover (``reach``, ``SearchQuery.reach``) and at most half
    the record. A hit across the origin lies in it whole, when it is no longer than the query."""

    width = min(int(reach), int(record_length) // 2)
    if width < 1:
        raise ValidationError(
            f"A record of {record_length} bases has no origin window.",
            diagnostic={"code": "REGION_INVALID", "field": "region", "reason": "RECORD_BOUNDS"},
        )
    return record_length - width + 1, width


def search_program(molecule: SequenceMolecule) -> SearchProgram:
    """The program that searches a query of ``molecule``: blastn for a nucleotide query, tblastn for a protein query."""

    return _SEARCH_PROGRAMS[molecule]


@dataclass(frozen=True)
class SearchQuery:
    """One query of a search: the name the LOSAT table gives it (its FASTA ID), its length (bases, or residues
    for a protein query), which places each hit in the query, and the program that searched it
    (``search_program``), which each hit's metadata names."""

    name: str
    length: int
    program: SearchProgram

    def __post_init__(self) -> None:
        if not str(self.name).strip():
            raise ValidationError("A search query needs a name.", diagnostic={"code": "INPUT_INVALID", "reason": "REQUIRED"})
        if isinstance(self.length, bool) or not isinstance(self.length, int) or self.length < 1:
            raise ValidationError(
                "A search query's length must be a positive integer.",
                diagnostic={"code": "INPUT_INVALID", "reason": "POSITIVE_INTEGER"},
            )
        if self.program not in _BASES_PER_QUERY_LETTER:
            raise ValidationError("A search program is blastn or tblastn.", diagnostic={"code": "INPUT_INVALID", "reason": "LOSAT_OPTION_PROGRAM"})

    @property
    def reach(self) -> int:
        """The most subject bases one hit of the query covers: its length, three bases per residue for tblastn."""

        return self.length * _BASES_PER_QUERY_LETTER[self.program]


@dataclass(frozen=True)
class SearchHit:
    """One hit as the hit table lists it, and the annotation that draws it."""

    annotation_id: str
    query: str
    record_key: str
    start: int
    end: int
    strand: Literal["+", "-"]
    identity: float
    alignment_length: int
    evalue: float
    bitscore: float
    qstart: int
    qend: int


@dataclass(frozen=True)
class SearchAnnotationResult:
    """The annotation set of one search, or ``None`` when no hit passed the thresholds.

    ``tally`` counts the rows read and those the thresholds removed (the
    panel's "No hits" text); ``unmatched_rows`` are rows whose subject names no
    searched sequence, left out (PD-OI-074), with up to three of their IDs.
    """

    annotation_set: AnnotationSet | None
    hits: tuple[SearchHit, ...]
    tally: ComparisonRowTally
    unmatched_rows: int = 0
    unmatched_subject_ids: tuple[str, ...] = ()

    @property
    def rows_removed_by_thresholds(self) -> int:
        return self.tally.rows_read - self.tally.rows_kept


def search_query_color(index: int) -> str:
    """The default hit colour of the ``index``-th query (0-based)."""

    return SEARCH_QUERY_COLORS[index % len(SEARCH_QUERY_COLORS)]


def _number_text(value: float) -> str:
    """A table value as the shortest text that reads back as it: ``1083``, ``70.968``, ``2.89e-15``."""

    text = repr(float(value))
    return text[:-2] if text.endswith(".0") else text


def _strand(sstart: int, send: int, qstart: int, qend: int) -> Literal["+", "-"]:
    """The hit's strand on the subject, the source strand of its record: the query runs along it forward."""

    return "+" if (sstart <= send) == (qstart <= qend) else "-"


@dataclass(frozen=True)
class _Placed:
    """One row where it lies: its query, record, strand, source span and query span."""

    index: object
    query: str
    record_key: str
    strand: Literal["+", "-"]
    start: int
    end: int
    qlow: int
    qhigh: int

    def holds(self, other: _Placed) -> bool:
        """Whether ``other``, a row of the record's own sequence, is a part of this alignment across the origin
        (its start after its end): the same query and strand, inside it on the record and in the query."""

        return (
            (other.query, other.record_key, other.strand) == (self.query, self.record_key, self.strand)
            and (other.start >= self.start or other.end <= self.end)
            and self.qlow <= other.qlow and other.qhigh <= self.qhigh
        )


def _one_row_per_alignment(
    frame: DataFrame, subject_of: Mapping[str, SearchSubject], thresholds: ComparisonFilter
) -> DataFrame:
    """The rows of ``frame`` with each alignment once (OV-2571).

    A row of an origin window (``SearchSubject.origin_window``) is kept only
    when it crosses the origin: one that
    does not is the same alignment as a row of the record's own sequence. A
    crossing row stands for its halves, the rows of the record's own sequence
    inside it (same query and strand); when it does not pass ``thresholds`` but a
    half does, the halves that pass stand for it instead.
    """

    passing = set(filter_comparison_dataframe(frame, thresholds).index)
    crossing: list[_Placed] = []
    record_rows: list[_Placed] = []
    dropped: set[object] = set()
    for index, row in zip(frame.index, frame.itertuples(index=False)):
        values = dict(zip(COMPARISON_COLUMNS, row))
        subject = subject_of[str(values["subject"]).strip()]
        sstart, send, qstart, qend = (int(values[key]) for key in ("sstart", "send", "qstart", "qend"))
        start, end, wraps = subject.source_span(min(sstart, send), max(sstart, send))
        placed = _Placed(
            index, str(values["query"]).strip(), subject.record_key, _strand(sstart, send, qstart, qend),
            start, end, min(qstart, qend), max(qstart, qend),
        )
        if not subject.origin_window:
            record_rows.append(placed)
        elif wraps:
            crossing.append(placed)
        else:
            dropped.add(index)
    for whole in crossing:
        halves = [half for half in record_rows if whole.holds(half)]
        stands = whole.index in passing or not any(half.index in passing for half in halves)
        if not stands:
            dropped.add(whole.index)
        dropped.update(half.index for half in halves if stands or half.index not in passing)
    return frame.drop(index=list(dropped))


def search_hit_annotation_set(
    table: DataFrame,
    *,
    set_id: str,
    subjects: Sequence[SearchSubject],
    thresholds: ComparisonFilter,
    queries: Sequence[SearchQuery],
    color_mode: SearchColorMode = "query",
    query_colors: Mapping[str, str] | None = None,
) -> SearchAnnotationResult:
    """The annotation set of one search: one ``highlight`` row per hit that passes ``thresholds``.

    ``table`` holds outfmt 6 rows (``read_comparison_table``); ``queries`` are
    the queries in the panel's order, which gives each query its default
    colour, with their lengths and programs (a search of nucleotide and protein
    queries holds the rows of both programs); ``query_colors`` replaces a query's colour. A hit's target is the
    span ``sstart``..``send`` of the subject it names, in source coordinates,
    and its ``strand`` the strand the query runs forward on there;
    its metadata holds ``SEARCH_HIT_METADATA_KEYS``. A row of a subject read
    across the origin (``SearchSubject.record_length``) gives a target with ``wraps_origin`` when
    it crosses the origin; each alignment is one row and counts once in the
    tally (``_one_row_per_alignment``). Every hit is drawn at
    ``SEARCH_HIT_FILL_OPACITY``. With ``color_mode``
    ``gradient`` the scale runs from ``SEARCH_GRADIENT_START_COLOR`` at the
    query's start to ``SEARCH_GRADIENT_END_COLOR`` at its end; each hit draws
    its part of the scale (``qstart``..``qend`` of the query's length),
    following the hit's strand, and its Legend row is
    ``QUERY_POSITION_LEGEND_LABEL``. A row whose query is not in ``queries``,
    or that runs past its query's length, is refused.
    """

    if color_mode not in {"query", "gradient"}:
        raise ValidationError("A search colour mode is query or gradient.", diagnostic={"code": "INPUT_INVALID", "reason": "FIELDS"})
    if not str(set_id).strip():
        raise ValidationError("A search annotation set needs an id.", diagnostic={"code": "INPUT_INVALID", "reason": "REQUIRED"})
    frame = normalize_comparison_dataframe(table)
    index = table_id_index([[subject.subject_id] for subject in subjects])
    named = {
        value: table_id_matches(value, index) for value in dict.fromkeys(frame["subject"].astype(str).str.strip())
    }
    ambiguous = next((value for value, items in named.items() if len(items) > 1), None)
    if ambiguous is not None:
        raise ValidationError(
            f"The search subject ID {ambiguous!r} names more than one searched sequence.",
            diagnostic={"code": "COMPARISON_INPUT", "reason": "RECORD_ID", "column": 2},
        )
    unknown = [value for value, items in named.items() if not items]
    subject_of = {value: subjects[items[0]] for value, items in named.items() if items}
    known = _one_row_per_alignment(
        frame[~frame["subject"].astype(str).str.strip().isin(unknown)], subject_of, thresholds
    )
    tally = comparison_row_tally(known, thresholds)
    kept = filter_comparison_dataframe(known, thresholds)
    lengths = {query.name: query.length for query in queries}
    programs = {query.name: query.program for query in queries}
    colors = {
        query.name: (query_colors or {}).get(query.name) or search_query_color(order)
        for order, query in enumerate(queries)
    }
    annotations: list[RegionAnnotation] = []
    hits: list[SearchHit] = []
    taken: set[str] = set()
    counts: dict[str, int] = {}
    for row in kept.itertuples(index=False):
        values = dict(zip(COMPARISON_COLUMNS, row))
        query = str(values["query"]).strip()
        subject = subject_of[str(values["subject"]).strip()]
        sstart, send = int(values["sstart"]), int(values["send"])
        qstart, qend = int(values["qstart"]), int(values["qend"])
        if query not in lengths:
            raise ValidationError(
                f"The search query {query!r} is not one of the queries searched.",
                diagnostic={"code": "COMPARISON_INPUT", "column": 1},
            )
        if max(qstart, qend) > lengths[query]:
            raise ValidationError(
                f"A hit of {query!r} ends past the query's length {lengths[query]}.",
                diagnostic={"code": "COMPARISON_INPUT", "column": 8},
            )
        strand = _strand(sstart, send, qstart, qend)
        forward = strand == "+"
        counts[query] = counts.get(query, 0)
        annotation_id = ""
        while not annotation_id or annotation_id in taken:
            counts[query] += 1
            annotation_id = f"{query}#{counts[query]}"
        taken.add(annotation_id)
        if color_mode == "gradient":
            style = RegionAnnotationStyle(
                stroke=SEARCH_GRADIENT_END_COLOR,
                stroke_width=SEARCH_GRADIENT_STROKE_WIDTH,
                fill=SEARCH_GRADIENT_END_COLOR,
                fill_opacity=SEARCH_HIT_FILL_OPACITY,
                gradient=AnnotationGradient(
                    SEARCH_GRADIENT_START_COLOR, SEARCH_GRADIENT_END_COLOR, "forward" if forward else "reverse",
                    query_start=(min(qstart, qend) - 1) / lengths[query], query_end=max(qstart, qend) / lengths[query],
                ),
            )
            legend_label = QUERY_POSITION_LEGEND_LABEL
        else:
            color = colors[query]
            style = RegionAnnotationStyle(stroke=color, fill=color, fill_opacity=SEARCH_HIT_FILL_OPACITY)
            legend_label = query
        start, end, wraps = subject.source_span(min(sstart, send), max(sstart, send))
        metadata = {
            "query": query,
            "program": programs[query],
            "identity": _number_text(values["identity"]),
            "alignment_length": str(int(values["alignment_length"])),
            "evalue": _number_text(values["evalue"]),
            "bitscore": _number_text(values["bitscore"]),
        }
        annotations.append(RegionAnnotation(
            id=annotation_id,
            target=CoordinateSpan(record=subject.record, start=start, end=end, wraps_origin=wraps),
            label="",
            mark="highlight",
            style=style,
            legend_label=legend_label,
            metadata=metadata,
            strand=strand,
        ))
        hits.append(SearchHit(
            annotation_id=annotation_id,
            query=query,
            record_key=subject.record_key,
            start=start,
            end=end,
            strand=strand,
            identity=float(values["identity"]),
            alignment_length=int(values["alignment_length"]),
            evalue=float(values["evalue"]),
            bitscore=float(values["bitscore"]),
            qstart=qstart,
            qend=qend,
        ))
    return SearchAnnotationResult(
        annotation_set=(
            AnnotationSet(id=set_id, annotations=tuple(annotations), keep_highlight_overlay=True)
            if annotations else None
        ),
        hits=tuple(hits),
        tally=tally,
        unmatched_rows=len(frame) - len(known),
        unmatched_subject_ids=tuple(unknown[:3]),
    )


__all__ = [
    "QUERY_POSITION_LEGEND_LABEL",
    "SEARCH_GRADIENT_END_COLOR",
    "SEARCH_GRADIENT_START_COLOR",
    "SEARCH_HIT_METADATA_KEYS",
    "SEARCH_QUERY_COLORS",
    "SearchAnnotationResult",
    "SearchColorMode",
    "SearchHit",
    "SearchProgram",
    "SearchQuery",
    "SearchSubject",
    "origin_window_span",
    "search_hit_annotation_set",
    "search_program",
    "search_query_color",
]
