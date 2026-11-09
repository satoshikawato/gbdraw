#!/usr/bin/env python
# coding: utf-8

"""The locate owner: where a reference span lies in other genomes.

A reference record's windows (the user's pick) give one query per gene: each
CDS wholly inside a window, as protein (``tlosatn``, tblastn) or as nucleotide
(``losatn``, blastn). One search per target File, planned by the shared LOSAT
job owner (:func:`plan_locate_search`), returns raw outfmt 6 rows, and
:func:`plan_locate` turns the reference, the target records and those rows
into one :class:`LocatePlan`: per File the windows to draw (``RecordWindow``
values, the windows of PR-W's records), the candidates not selected, gene
coverage and warnings. The command line and the Web run the same search plan
and hand the rows back to the same function; the records table, the
comparisons table, the new drawing and the dialog read the plan.

Rules (every step is a pure function of its inputs, and every iteration runs
over a sorted sequence):

1. Rows pass the comparison thresholds (:func:`locate_default_thresholds`).
2. One locus per gene: the hits of one gene on one record and strand chain
   when they lie within the gene's length; a locus keeps the best hits that
   each add at least half new query letters. Anchors are a gene's loci that
   score at least ``best_ratio`` of its best locus in the File.
3. An anchor takes the whole target CDS on its strand that holds at least half
   of it.
4. Anchors of one record cluster while the bases between them stay within
   ``merge_distance``; a circular record joins its last and first cluster
   across the origin.
5. Clusters, ranked by evidence, are selected while they add enough reference
   genes not covered yet (:func:`_genes_needed`), at most ``max_windows``
   windows; the others stay candidates.
6. A cluster reads in its orientation (reversed when most of its genes lie on
   the opposite strand); it splits where the reference order falls back by a jump
   (:func:`_jump`, :func:`_split_runs`). On a linear record a stored rotation
   of the reference gives abutting windows in reference order. A circular
   record has no ends: a cluster around the whole record is read from where
   the reference order starts (:func:`_circular_frame`), and a window across
   the origin is one window with ``start > end`` (Owner D-28, 2026-11-08).
7. A window is its genes plus ``flank`` on each side, or to the end of a
   linear record when less than ``flank`` remains (PD-OI-123). It is
   reverse-complemented when most of its genes lie on the opposite strand of
   the reference as drawn; a tie keeps it as drawn (PD-OI-123). Windows are
   ordered by the reference order of their genes.
"""

from __future__ import annotations

from dataclasses import dataclass, fields, is_dataclass, replace
import io
import math
import types
from typing import (
    Literal, Mapping, NamedTuple, Sequence, TypedDict, Union, cast, get_args, get_origin, get_type_hints, is_typeddict,
)

from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.protein_colinearity import extract_web_stable_cds_proteins, proteins_to_fasta
from gbdraw.comparisons.losat_jobs import (
    LosatBatch,
    losat_fasta_text,
    nucleotide_losat_cache_key,
    plan_losat_jobs,
    prepare_losat_batches,
    split_losat_batch_result,
)
from gbdraw.comparisons.losat_runtime import (
    DEFAULT_LOSATN_TASK,
    LOSAT_PROGRAMS,
    LosatRuntimeRecord,
    LosatSearchArgs,
    losat_cache_args,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import (
    ComparisonRowTally,
    comparison_row_tally,
    filter_comparison_dataframe,
    read_comparison_table,
)
from gbdraw.io.genome import filter_features_by_type, load_fasta, load_gbks, parse_gff_fasta
from gbdraw.io.record_select import parse_record_selector, select_record
from gbdraw.io.regions import RecordWindow
from gbdraw.mode_profiles import LINEAR_MODE_PROFILE, ComparisonThresholds

Strand = Literal["+", "-"]
LocateProgramName = Literal["tlosatn", "losatn"]
# The programs a locate search runs, the default first (D-03: tblastn).
LOCATE_PROGRAMS: tuple[LocateProgramName, ...] = get_args(LocateProgramName)
GeneStatus = Literal["complete", "partial"]
CandidateReason = Literal["TOO_FEW_GENES", "COVERED_BY_OTHER_WINDOW", "MAX_WINDOWS"]
LocateWarningCode = Literal[
    "FILE_NO_WINDOWS",
    "FILE_SEARCH_FAILED",
    "GENES_ABSENT",
    "GENES_PARTIAL",
    "MAX_WINDOWS",
    "REFERENCE_GENE_CUT",
    "TARGET_NOT_ANNOTATED",
]
LocateArg = Union[str, int, float]
# Why a window cannot be a reference window or a window of the drawing (:class:`LocateWindowRefusal`).
LocateWindowRefusalReason = Literal[
    "WINDOW_PAST_END",
    "WINDOW_ACROSS_LINEAR_ORIGIN",
    "WINDOWS_OVERLAP",
    "WINDOW_WITHOUT_CDS",
]

# A gene found with at least this share of its query aligned, away from a record end, is complete.
COMPLETE_QUERY_COVERAGE = 0.9
_INVALID = {"code": "INPUT_INVALID", "reason": "FIELDS"}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValidationError(message, diagnostic=_INVALID)


def _positive_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value > 0


class LocateWindowRefusal(ValidationError):
    """A window a locate plan cannot take, by the rule it breaks (``reason``), with the window's span, and for a
    window of a File the File key and column that name it (``place``)."""

    def __init__(
        self, message: str, reason: LocateWindowRefusalReason, start: int, end: int,
        place: tuple[str, int] | None = None,
    ) -> None:
        super().__init__(message, diagnostic=_INVALID)
        self.reason: LocateWindowRefusalReason = reason
        self.start = start
        self.end = end
        self.place = place


def _refuse_window(condition: bool, message: str, reason: LocateWindowRefusalReason, window: RecordWindow) -> None:
    if not condition:
        raise LocateWindowRefusal(message, reason, window.start, window.end)


def check_record_window(window: RecordWindow, record_length: int, circular: bool) -> None:
    """Refuse a window that runs past the record end, or that crosses the origin of a linear record."""

    _refuse_window(
        max(window.start, window.end) <= record_length,
        f"Reference window {window.start}-{window.end} runs past the record end ({record_length}).",
        "WINDOW_PAST_END",
        window,
    )
    _refuse_window(
        circular or window.start <= window.end,
        f"Reference window {window.start}-{window.end} crosses the origin of a linear record.",
        "WINDOW_ACROSS_LINEAR_ORIGIN",
        window,
    )


def locate_default_thresholds(program: LocateProgramName) -> ComparisonThresholds:
    """The row thresholds of a locate search with ``program``.

    Protein queries (tblastn) use the Linear mode's comparison thresholds
    too, with no separate protein set (Owner, 2026-11-08).
    """

    _require(program in LOCATE_PROGRAMS, f"Unknown locate program {program!r}.")
    return LINEAR_MODE_PROFILE.comparison


# Inputs.


@dataclass(frozen=True)
class ReferenceCds:
    """One CDS of the reference record: 1-based ``start``..``end`` on ``strand``.

    ``start > end`` is a CDS across the origin of a circular record.
    ``protein_id`` names its query in the search rows; ``length_aa`` and
    ``length_nt`` are its protein and coding lengths.
    """

    protein_id: str
    locus_tag: str | None
    gene: str | None
    product: str | None
    start: int
    end: int
    strand: Strand
    length_aa: int
    length_nt: int

    @property
    def name(self) -> str:
        """The name a report gives the CDS: its gene name, else its locus tag, else its protein ID."""

        return self.gene or self.locus_tag or self.protein_id

    def __post_init__(self) -> None:
        _require(bool(str(self.protein_id).strip()), "A reference CDS needs a protein ID.")
        _require(
            _positive_int(self.start) and _positive_int(self.end),
            f"Reference CDS {self.protein_id} needs 1-based coordinates.",
        )
        _require(self.strand in ("+", "-"), f"Reference CDS {self.protein_id} needs strand + or -.")
        _require(
            _positive_int(self.length_aa) and _positive_int(self.length_nt),
            f"Reference CDS {self.protein_id} needs positive lengths.",
        )


@dataclass(frozen=True)
class ReferenceGene:
    """A query of the search: a reference CDS wholly inside window ``window_index``.

    ``order`` is its position along the drawn reference: windows in column
    order, the genes of a reversed window in reverse.
    """

    cds: ReferenceCds
    window_index: int
    order: int

    @property
    def protein_id(self) -> str:
        return self.cds.protein_id


@dataclass(frozen=True)
class LocateReference:
    """The reference windows and their genes (:func:`locate_reference`).

    ``windows`` are in column order; ``cut_genes`` are the CDS a window edge
    cuts, by protein ID, which are not searched.
    """

    file_key: str
    record_key: str
    record_id: str
    record_length: int
    circular: bool
    windows: tuple[RecordWindow, ...]
    genes: tuple[ReferenceGene, ...]
    cut_genes: tuple[ReferenceCds, ...]


def _gene_strand_as_drawn(gene: ReferenceGene, windows: Sequence[RecordWindow]) -> Strand:
    if not windows[gene.window_index].reverse_complement:
        return gene.cds.strand
    return "-" if gene.cds.strand == "+" else "+"


def _arc(start: int, end: int, length: int) -> tuple[int, int]:
    """A span as an interval that may run past the record end (``start > end`` crosses the origin)."""

    return (start, end) if start <= end else (start, end + length)


def _shifted(span: tuple[int, int], length: int) -> tuple[tuple[int, int], ...]:
    return tuple((span[0] + shift, span[1] + shift) for shift in (-length, 0, length))


def _contains(outer: tuple[int, int], inner: tuple[int, int], length: int) -> bool:
    return any(outer[0] <= low and high <= outer[1] for low, high in _shifted(inner, length))


def _overlaps(first: tuple[int, int], second: tuple[int, int], length: int) -> bool:
    return any(low <= first[1] and first[0] <= high for low, high in _shifted(second, length))


def locate_reference(
    *,
    file_key: str,
    record_key: str,
    record_id: str,
    record_length: int,
    circular: bool,
    windows: Sequence[RecordWindow],
    cds: Sequence[ReferenceCds],
) -> LocateReference:
    """The reference of a locate search: the genes of each window in drawn order.

    Windows are read in column order (a window without a column keeps its
    place). A window must lie in the record, share no base with another, and
    hold a whole CDS; on a circular record it may cross the origin.
    """

    _require(_positive_int(record_length), "The reference record needs a positive length.")
    _require(bool(windows), "A locate reference needs at least one window.")
    ordered = tuple(
        window
        for _key, window in sorted(
            ((window.column if window.column is not None else index + 1, index), window)
            for index, window in enumerate(windows)
        )
    )
    arcs = [_arc(window.start, window.end, record_length) for window in ordered]
    for window in ordered:
        check_record_window(window, record_length, circular)
    for index, arc in enumerate(arcs):
        for other_index, other in enumerate(arcs[index + 1:], start=index + 1):
            _refuse_window(
                not _overlaps(arc, other, record_length), "Reference windows overlap.", "WINDOWS_OVERLAP",
                ordered[other_index],
            )
    genes: list[ReferenceGene] = []
    cut: dict[str, ReferenceCds] = {}
    for window_index, (window, arc) in enumerate(zip(ordered, arcs)):
        inside = []
        for item in cds:
            span = _arc(item.start, item.end, record_length)
            if _contains(arc, span, record_length):
                inside.append(item)
            elif _overlaps(arc, span, record_length):
                cut[item.protein_id] = item
        _refuse_window(
            bool(inside), f"Reference window {window.start}-{window.end} holds no whole CDS.", "WINDOW_WITHOUT_CDS", window
        )
        # Positions along the window: a CDS after the origin follows the record end.
        inside.sort(key=lambda item: (item.start + (record_length if item.start < window.start else 0), item.protein_id))
        if window.reverse_complement:
            inside.reverse()
        first = len(genes)
        genes.extend(ReferenceGene(item, window_index, first + rank) for rank, item in enumerate(inside))
    ids = [gene.protein_id for gene in genes]
    _require(len(set(ids)) == len(ids), "Reference CDS protein IDs must be unique.")
    return LocateReference(
        file_key=str(file_key),
        record_key=str(record_key),
        record_id=str(record_id),
        record_length=int(record_length),
        circular=bool(circular),
        windows=ordered,
        genes=tuple(genes),
        cut_genes=tuple(cut[key] for key in sorted(cut)),
    )


@dataclass(frozen=True)
class CdsExtent:
    """One location part of a target CDS (a CDS across the origin is two parts)."""

    start: int
    end: int
    strand: Strand


@dataclass(frozen=True)
class LocateRecord:
    """One record of a target File; ``index`` is its order in the File."""

    record_key: str
    record_id: str
    index: int
    length: int
    circular: bool
    cds: tuple[CdsExtent, ...]


@dataclass(frozen=True)
class LocateTarget:
    """One target File: one genome, searched whole."""

    file_key: str
    records: tuple[LocateRecord, ...]


@dataclass(frozen=True)
class LocateOptions:
    """The program, row thresholds and plan rules of a locate search (D-L3)."""

    program: LocateProgramName
    thresholds: ComparisonThresholds
    best_ratio: float = 0.8
    merge_distance: int = 10_000
    flank: int = 1_000
    min_genes: int = 2
    min_share: float = 0.05
    jump_min: int = 5
    jump_share: float = 0.25
    max_windows: int = 8

    def __post_init__(self) -> None:
        _require(self.program in LOCATE_PROGRAMS, f"Unknown locate program {self.program!r}.")
        _require(isinstance(self.thresholds, ComparisonThresholds), "Locate thresholds must be ComparisonThresholds.")
        for name in ("best_ratio", "min_share", "jump_share"):
            value = getattr(self, name)
            _require(
                isinstance(value, (int, float)) and not isinstance(value, bool) and 0 <= value <= 1,
                f"Locate option {name} must be between 0 and 1.",
            )
        for name in ("merge_distance", "flank"):
            value = getattr(self, name)
            _require(
                isinstance(value, int) and not isinstance(value, bool) and value >= 0,
                f"Locate option {name} must be a non-negative integer.",
            )
        for name in ("min_genes", "jump_min", "max_windows"):
            _require(_positive_int(getattr(self, name)), f"Locate option {name} must be a positive integer.")

    @classmethod
    def for_program(cls, program: LocateProgramName = "tlosatn") -> LocateOptions:
        """The default options of ``program``."""

        return cls(program, locate_default_thresholds(program))


@dataclass(frozen=True)
class LocateSearchResult:
    """The raw outfmt 6 rows of one File's search, with original record IDs.

    ``error`` holds why the search failed; its File gets ``FILE_SEARCH_FAILED``.
    """

    file_key: str
    text: str
    error: str | None = None


# The plan.


@dataclass(frozen=True)
class WindowGene:
    """One reference gene found in a window: its locus on the target record."""

    protein_id: str
    order: int
    window_index: int
    start: int
    end: int
    strand: Strand
    query_coverage: float
    bits: float
    status: GeneStatus


@dataclass(frozen=True)
class PlannedWindow:
    """One window of a target record, selected or a candidate.

    ``window.column`` is the window's place among every window of its File in
    reference order, candidates included, so selecting a candidate changes
    ``selected`` only. Windows that share ``cluster`` come from one cluster: a
    rotation of the reference or an arc across the origin.
    ``gap_to_next_same_record`` is the number of bases left out between this
    window and the next selected window when that one is of the same record
    (0 when they abut), else None. ``ref_span`` is the first and last
    reference coordinate of its genes. ``touches_start`` and ``touches_end``
    say the window reaches the start or end of a linear record (a contig
    end); a circular record (``record_circular``) has none.
    """

    record_key: str
    record_id: str
    record_index: int
    record_length: int
    record_circular: bool
    window: RecordWindow
    selected: bool
    reason: CandidateReason | None
    cluster: int
    touches_start: bool
    touches_end: bool
    across_origin: bool
    gap_to_next_same_record: int | None
    evidence_bits: float
    genes: tuple[WindowGene, ...]
    ref_span: tuple[int, int]


@dataclass(frozen=True)
class GeneCoverage:
    """The reference genes by what the selected windows hold, as protein IDs in reference order.

    ``reference_not_found`` are the runs of reference genes that are not
    complete, as reference coordinates.
    """

    complete: tuple[str, ...]
    partial: tuple[str, ...]
    absent: tuple[str, ...]
    reference_not_found: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class LocateWarning:
    code: LocateWarningCode
    args: tuple[tuple[str, LocateArg], ...] = ()


@dataclass(frozen=True)
class FilePlan:
    file_key: str
    windows: tuple[PlannedWindow, ...]
    candidates: tuple[PlannedWindow, ...]
    coverage: GeneCoverage
    tally: ComparisonRowTally
    warnings: tuple[LocateWarning, ...]


@dataclass(frozen=True)
class LocatePlan:
    """The plan of one locate search; it holds no file paths and no Session fields."""

    reference: LocateReference
    options: LocateOptions
    runtime: LosatRuntimeRecord | None
    files: tuple[FilePlan, ...]
    warnings: tuple[LocateWarning, ...]

    def to_json(self) -> LocatePlanJson:
        """The plan as JSON data with stable keys: each dataclass as an object of its fields (``LocatePlanJson``),
        a warning's ``args`` as an object."""

        value = _json_value(self)
        assert isinstance(value, dict)
        return cast("LocatePlanJson", value)

    @classmethod
    def from_json(cls, value: LocatePlanJson) -> LocatePlan:
        """The plan that :meth:`to_json` wrote as ``value``; a value of another shape is refused."""

        plan = _from_json(cls, value, "plan")
        assert isinstance(plan, LocatePlan)
        return plan


# The plan as JSON (``LocatePlan.to_json``): one object type per dataclass, with its field names.


class RecordWindowJson(TypedDict):
    start: int
    end: int
    reverse_complement: bool
    column: int | None


class ReferenceCdsJson(TypedDict):
    protein_id: str
    locus_tag: str | None
    gene: str | None
    product: str | None
    start: int
    end: int
    strand: Strand
    length_aa: int
    length_nt: int


class ReferenceGeneJson(TypedDict):
    cds: ReferenceCdsJson
    window_index: int
    order: int


class LocateReferenceJson(TypedDict):
    file_key: str
    record_key: str
    record_id: str
    record_length: int
    circular: bool
    windows: list[RecordWindowJson]
    genes: list[ReferenceGeneJson]
    cut_genes: list[ReferenceCdsJson]


class ComparisonThresholdsJson(TypedDict):
    evalue: float
    bitscore: float
    identity: float
    alignment_length: int


class LocateOptionsJson(TypedDict):
    program: LocateProgramName
    thresholds: ComparisonThresholdsJson
    best_ratio: float
    merge_distance: int
    flank: int
    min_genes: int
    min_share: float
    jump_min: int
    jump_share: float
    max_windows: int


class WindowGeneJson(TypedDict):
    protein_id: str
    order: int
    window_index: int
    start: int
    end: int
    strand: Strand
    query_coverage: float
    bits: float
    status: GeneStatus


class PlannedWindowJson(TypedDict):
    record_key: str
    record_id: str
    record_index: int
    record_length: int
    record_circular: bool
    window: RecordWindowJson
    selected: bool
    reason: CandidateReason | None
    cluster: int
    touches_start: bool
    touches_end: bool
    across_origin: bool
    gap_to_next_same_record: int | None
    evidence_bits: float
    genes: list[WindowGeneJson]
    ref_span: list[int]


class GeneCoverageJson(TypedDict):
    complete: list[str]
    partial: list[str]
    absent: list[str]
    reference_not_found: list[list[int]]


class ComparisonThresholdTallyJson(TypedDict):
    name: str
    threshold: float
    rows_left_out: int
    best_value: float


class ComparisonRowTallyJson(TypedDict):
    rows_read: int
    rows_kept: int
    thresholds: list[ComparisonThresholdTallyJson]


class LocateWarningJson(TypedDict):
    code: LocateWarningCode
    args: dict[str, LocateArg]


class FilePlanJson(TypedDict):
    file_key: str
    windows: list[PlannedWindowJson]
    candidates: list[PlannedWindowJson]
    coverage: GeneCoverageJson
    tally: ComparisonRowTallyJson
    warnings: list[LocateWarningJson]


class LocatePlanJson(TypedDict):
    reference: LocateReferenceJson
    options: LocateOptionsJson
    runtime: LosatRuntimeRecord | None
    files: list[FilePlanJson]
    warnings: list[LocateWarningJson]


def _json_value(value: object) -> object:
    if isinstance(value, LocateWarning):
        return {"code": value.code, "args": {key: item for key, item in value.args}}
    if is_dataclass(value) and not isinstance(value, type):
        return {item.name: _json_value(getattr(value, item.name)) for item in fields(value)}
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


def _from_json(hint: object, value: object, path: str) -> object:
    """A value of type ``hint`` from the JSON :func:`_json_value` wrote for it."""

    if hint is LocateWarning:
        _require(isinstance(value, Mapping) and isinstance(value.get("args"), Mapping), f"{path} is no locate warning.")
        assert isinstance(value, Mapping)
        return LocateWarning(
            cast(LocateWarningCode, _from_json(LocateWarningCode, value.get("code"), f"{path}.code")),
            tuple((str(key), cast(LocateArg, item)) for key, item in value["args"].items()),
        )
    if isinstance(hint, type) and is_dataclass(hint):
        _require(isinstance(value, Mapping), f"{path} must be an object.")
        assert isinstance(value, Mapping)
        hints = get_type_hints(hint)
        names = [item.name for item in fields(hint) if item.init]
        _require(sorted(value) == sorted(names), f"{path} must hold the fields {sorted(names)}.")
        return hint(**{name: _from_json(hints[name], value[name], f"{path}.{name}") for name in names})
    if isinstance(hint, type) and is_typeddict(hint):
        _require(isinstance(value, Mapping), f"{path} must be an object.")
        return dict(cast(Mapping[str, object], value))
    origin, args = get_origin(hint), get_args(hint)
    if origin in (Union, types.UnionType):
        if value is None and type(None) in args:
            return None
        [branch] = [item for item in args if item is not type(None)]
        return _from_json(branch, value, path)
    if origin is Literal:
        _require(value in args, f"{path} must be one of {list(args)}.")
        return value
    if origin is tuple:
        _require(isinstance(value, list), f"{path} must be a list.")
        assert isinstance(value, list)
        if len(args) == 2 and args[1] is Ellipsis:
            return tuple(_from_json(args[0], item, f"{path}[{index}]") for index, item in enumerate(value))
        _require(len(value) == len(args), f"{path} must hold {len(args)} items.")
        return tuple(_from_json(kind, item, f"{path}[{index}]") for index, (kind, item) in enumerate(zip(args, value)))
    if hint is float:
        _require(isinstance(value, (int, float)) and not isinstance(value, bool), f"{path} must be a number.")
        return float(cast(float, value))
    if hint in (int, str, bool):
        _require(type(value) is hint, f"{path} must be of type {cast(type, hint).__name__}.")
        return value
    raise TypeError(f"No JSON reading is defined for {hint!r}.")


# Choices on a plan: the windows the user keeps and their edits.


@dataclass(frozen=True)
class LocateWindowChoice:
    """A window the user keeps for the drawing: its File, its ``column`` among the File's windows and candidates
    (its identity in the plan), and its span and orientation (the search's, or the user's edit).

    ``start > end`` is a window across the origin of a circular record.
    """

    file_key: str
    column: int
    start: int
    end: int
    reverse_complement: bool


def _chosen_window(planned: PlannedWindow, choice: LocateWindowChoice | None) -> PlannedWindow:
    if choice is None:
        return replace(planned, selected=False, gap_to_next_same_record=None)
    window = RecordWindow(choice.start, choice.end, choice.reverse_complement, planned.window.column)
    try:
        check_record_window(window, planned.record_length, planned.record_circular)
    except LocateWindowRefusal as refusal:
        raise LocateWindowRefusal(
            str(refusal), refusal.reason, refusal.start, refusal.end, (choice.file_key, choice.column)
        ) from None
    touches_start, touches_end = _record_ends(window, planned.record_length, planned.record_circular)
    return replace(
        planned,
        window=window,
        selected=True,
        gap_to_next_same_record=None,
        touches_start=touches_start,
        touches_end=touches_end,
        across_origin=window.crosses_origin,
    )


def _genes_warnings(warnings: Sequence[LocateWarning], coverage: GeneCoverage, has_windows: bool) -> tuple[LocateWarning, ...]:
    """The search's warnings with the gene counts of ``coverage`` (no gene warnings for a File without windows)."""

    kept = [warning for warning in warnings if warning.code not in ("GENES_PARTIAL", "GENES_ABSENT")]
    if has_windows and coverage.partial:
        kept.append(LocateWarning("GENES_PARTIAL", (("count", len(coverage.partial)),)))
    if has_windows and coverage.absent:
        kept.append(LocateWarning("GENES_ABSENT", (("count", len(coverage.absent)),)))
    return tuple(kept)


def choose_locate_windows(plan: LocatePlan, choices: Sequence[LocateWindowChoice]) -> LocatePlan:
    """``plan`` with the windows the user keeps: each window of ``choices`` is selected with its span and
    orientation; every other window of the plan is not. A window across the origin of a linear record, or past
    its end, is refused (:class:`LocateWindowRefusal`).

    A window keeps the genes the search found in it (an edit does not count them again, PD-OI-123); ``reason``
    stays the search's. Coverage and the gene warnings follow the windows kept.
    """

    by_place = {(choice.file_key, choice.column): choice for choice in choices}
    _require(len(by_place) == len(choices), "A locate window is chosen twice.")
    columns = {(file.file_key, item.window.column) for file in plan.files for item in (*file.windows, *file.candidates)}
    _require(set(by_place) <= columns, "A chosen locate window is not in the plan.")
    files = []
    for file in plan.files:
        every = sorted((*file.windows, *file.candidates), key=lambda item: item.window.column or 0)
        chosen = [_chosen_window(item, by_place.get((file.file_key, item.window.column or 0))) for item in every]
        selected = [item for item in chosen if item.selected]
        gaps = {
            id(item): _gap(item.window, following.window)
            for item, following in zip(selected, selected[1:])
            if item.record_index == following.record_index
        }
        chosen = [replace(item, gap_to_next_same_record=gaps.get(id(item))) for item in chosen]
        windows = tuple(item for item in chosen if item.selected)
        coverage = _coverage(plan.reference, windows)
        files.append(replace(
            file,
            windows=windows,
            candidates=tuple(item for item in chosen if not item.selected),
            coverage=coverage,
            warnings=_genes_warnings(file.warnings, coverage, bool(windows)),
        ))
    return replace(plan, files=tuple(files))


# Algorithm.


@dataclass(frozen=True)
class _Locus:
    gene: ReferenceGene
    record: LocateRecord
    strand: Strand
    start: int
    end: int
    bits: float
    query_coverage: float


@dataclass(frozen=True)
class _Anchor:
    locus: _Locus
    # The snapped span: whole CDS; ``offset`` moves it into an origin-crossing cluster's frame.
    start: int
    end: int
    offset: int
    opposite: bool

    @property
    def low(self) -> int:
        return self.start + self.offset

    @property
    def high(self) -> int:
        return self.end + self.offset

    @property
    def order(self) -> int:
        return self.locus.gene.order


@dataclass(frozen=True)
class _Cluster:
    record: LocateRecord
    anchors: tuple[_Anchor, ...]

    @property
    def bits(self) -> float:
        return math.fsum(anchor.locus.bits for anchor in self.anchors)

    @property
    def genes(self) -> frozenset[int]:
        return frozenset(anchor.order for anchor in self.anchors)

    @property
    def low(self) -> int:
        return min(anchor.low for anchor in self.anchors)


class _Hsp(NamedTuple):
    gene: ReferenceGene
    record: LocateRecord
    strand: Strand
    low: int
    high: int
    qstart: int
    qend: int
    bits: float

    @property
    def chain_key(self) -> tuple[int, int, Strand]:
        return self.gene.order, self.record.index, self.strand


def _hsps(
    result: LocateSearchResult,
    reference: LocateReference,
    target: LocateTarget,
    thresholds: ComparisonThresholds,
) -> tuple[ComparisonRowTally, list[_Hsp]]:
    """The rows that pass ``thresholds``, in a total order that does not depend on the row order (rule 1)."""

    frame = read_comparison_table(io.StringIO(result.text), label=f"locate rows of {result.file_key}")
    tally = comparison_row_tally(frame, thresholds)
    kept = filter_comparison_dataframe(frame, thresholds)
    genes = {gene.protein_id: gene for gene in reference.genes}
    records = {record.record_id: record for record in target.records}
    hsps = []
    for query, subject, qstart, qend, sstart, send, bits in zip(
        kept["query"], kept["subject"], kept["qstart"], kept["qend"], kept["sstart"], kept["send"], kept["bitscore"]
    ):
        gene = genes.get(str(query))
        record = records.get(str(subject))
        if gene is None or record is None:
            raise ValidationError(
                f"Locate rows of {result.file_key} name a query or record the search did not have.",
                diagnostic={"code": "COMPARISON_INPUT", "reason": "FIELDS"},
            )
        low, high = sorted((int(sstart), int(send)))
        hsps.append(_Hsp(
            gene, record, "+" if int(send) >= int(sstart) else "-", low, high,
            min(int(qstart), int(qend)), max(int(qstart), int(qend)), float(bits),
        ))
    hsps.sort(key=lambda hsp: (*hsp.chain_key, hsp.low, hsp.high, hsp.qstart, hsp.qend, -hsp.bits))
    return tally, hsps


def _query_length(gene: ReferenceGene, program: LocateProgramName) -> int:
    query_molecule = LOSAT_PROGRAMS[program].molecules[0]
    return gene.cds.length_aa if query_molecule == "protein" else gene.cds.length_nt


def _loci(hsps: Sequence[_Hsp], program: LocateProgramName) -> list[_Locus]:
    """One locus per chain of a gene's hits on one record and strand within the gene's length (rule 2)."""

    loci: list[_Locus] = []
    index = 0
    while index < len(hsps):
        chain = [hsps[index]]
        high = chain[0].high
        index += 1
        while (
            index < len(hsps)
            and hsps[index].chain_key == chain[0].chain_key
            and hsps[index].low - high <= chain[0].gene.cds.length_nt
        ):
            chain.append(hsps[index])
            high = max(high, hsps[index].high)
            index += 1
        loci.append(_locus(chain, _query_length(chain[0].gene, program)))
    return loci


def _locus(chain: Sequence[_Hsp], query_length: int) -> _Locus:
    """The best hits of a chain that each add at least half new query letters."""

    covered = [False] * (query_length + 1)
    kept = []
    for hsp in sorted(chain, key=lambda item: (-item.bits, item.low, item.high, item.qstart, item.qend)):
        letters = range(hsp.qstart, min(hsp.qend, query_length) + 1)
        new = sum(1 for letter in letters if not covered[letter])
        if letters and new >= 0.5 * len(letters):
            kept.append(hsp)
            for letter in letters:
                covered[letter] = True
    if not kept:
        kept = [chain[0]]
    return _Locus(
        gene=chain[0].gene,
        record=chain[0].record,
        strand=chain[0].strand,
        start=min(hsp.low for hsp in kept),
        end=max(hsp.high for hsp in kept),
        bits=math.fsum(hsp.bits for hsp in kept),
        query_coverage=sum(covered[1:]) / query_length,
    )


def _snap(locus: _Locus) -> tuple[int, int]:
    """The locus widened to the same-strand CDS part that holds at least half of it (rule 3)."""

    width = locus.end - locus.start + 1
    for extent in sorted(locus.record.cds, key=lambda item: (item.start, item.end, item.strand)):
        overlap = min(extent.end, locus.end) - max(extent.start, locus.start) + 1
        if extent.strand == locus.strand and overlap >= 0.5 * width:
            return min(locus.start, extent.start), max(locus.end, extent.end)
    return locus.start, locus.end


def _anchors(loci: Sequence[_Locus], reference: LocateReference, options: LocateOptions) -> list[_Anchor]:
    best: dict[str, float] = {}
    for locus in loci:
        best[locus.gene.protein_id] = max(best.get(locus.gene.protein_id, 0.0), locus.bits)
    anchors = []
    for locus in loci:
        if locus.bits < options.best_ratio * best[locus.gene.protein_id]:
            continue
        start, end = _snap(locus)
        opposite = locus.strand != _gene_strand_as_drawn(locus.gene, reference.windows)
        anchors.append(_Anchor(locus, start, end, 0, opposite))
    return anchors


def _clusters(anchors: Sequence[_Anchor], merge_distance: int) -> list[_Cluster]:
    """Anchors of one record, clustered by position; a circular record joins across its origin (rule 4)."""

    by_record: dict[int, list[_Anchor]] = {}
    records: dict[int, LocateRecord] = {}
    for anchor in anchors:
        by_record.setdefault(anchor.locus.record.index, []).append(anchor)
        records[anchor.locus.record.index] = anchor.locus.record
    clusters: list[_Cluster] = []
    for record_index in sorted(by_record):
        record = records[record_index]
        ordered = sorted(by_record[record_index], key=lambda item: (item.start, item.end, item.order, item.locus.strand))
        groups: list[list[_Anchor]] = [[ordered[0]]]
        high = ordered[0].end
        for anchor in ordered[1:]:
            if anchor.start - high - 1 > merge_distance:
                groups.append([anchor])
                high = anchor.end
            else:
                groups[-1].append(anchor)
                high = max(high, anchor.end)
        if record.circular and len(groups) > 1:
            gap = (groups[0][0].start - 1) + (record.length - max(anchor.end for anchor in groups[-1]))
            if gap <= merge_distance:
                groups[0] = [_moved(item, -record.length) for item in groups.pop()] + groups[0]
        clusters.extend(_Cluster(record, tuple(group)) for group in groups)
    return clusters


def _moved(anchor: _Anchor, shift: int) -> _Anchor:
    return _Anchor(anchor.locus, anchor.start, anchor.end, anchor.offset + shift, anchor.opposite)


def _majority_window(cluster: _Cluster) -> int:
    weights: dict[int, float] = {}
    for anchor in cluster.anchors:
        weights[anchor.locus.gene.window_index] = weights.get(anchor.locus.gene.window_index, 0.0) + anchor.locus.bits
    return min(weights, key=lambda window_index: (-weights[window_index], window_index))


def _window_gene_counts(reference: LocateReference) -> dict[int, int]:
    counts: dict[int, int] = {}
    for gene in reference.genes:
        counts[gene.window_index] = counts.get(gene.window_index, 0) + 1
    return counts


def _genes_needed(window_genes: int, options: LocateOptions) -> int:
    """The new reference genes a cluster must add: all genes of a tiny window, else a share (rule 5)."""

    return min(window_genes, max(options.min_genes, math.ceil(options.min_share * window_genes)))


def _jump(window_genes: int, options: LocateOptions) -> float:
    """The fall in reference order that breaks a cluster (rule 6)."""

    return max(options.jump_min, options.jump_share * window_genes)


def _reversed(anchors: Sequence[_Anchor]) -> bool:
    """Whether most of the genes of ``anchors`` lie on the opposite strand of the reference as drawn (PD-OI-123,
    rules 6 and 7); a tie keeps them as drawn."""

    opposite = {anchor.locus.gene.protein_id for anchor in anchors if anchor.opposite}
    same = {anchor.locus.gene.protein_id for anchor in anchors if not anchor.opposite}
    return len(opposite) > len(same)


def _split_runs(cluster: _Cluster, jump: float, reverse: bool) -> list[list[_Anchor]]:
    """The parts of a cluster in position order: it breaks where the reference order falls by ``jump``.

    The cluster reads in its orientation; one rule covers a stored rotation of
    the reference. An inversion inside a cluster does not break it.
    """

    reading = sorted(cluster.anchors, key=lambda item: (item.low, item.high, item.order))
    if reverse:
        reading.reverse()
    runs = [[reading[0]]]
    for previous, anchor in zip(reading, reading[1:]):
        if previous.order - anchor.order >= jump:
            runs.append([])
        runs[-1].append(anchor)
    if reverse:
        runs = [list(reversed(run)) for run in reversed(runs)]
    return runs


def _circular_frame(
    cluster: _Cluster, runs: list[list[_Anchor]], jump: float, reverse: bool, merge_distance: int
) -> tuple[_Cluster, list[list[_Anchor]], bool]:
    """A cluster around the whole of a circular record, read from where the reference order starts.

    When a circular record holds the cluster all around (what lies outside
    it is within ``merge_distance``) and the cluster breaks, the record is
    read from the start of the part that holds the lowest reference order, so
    a stored rotation is one window across the origin. Returns the cluster,
    its parts, and whether the parts close the circle.
    """

    record = cluster.record
    span = max(anchor.high for anchor in cluster.anchors) - min(anchor.low for anchor in cluster.anchors) + 1
    if not record.circular or len(runs) < 2 or record.length - span > merge_distance:
        return cluster, runs, False
    first = min(range(len(runs)), key=lambda index: (min(anchor.order for anchor in runs[index]), index))
    if reverse:
        moved = [anchor for run in runs[:first + 1] for anchor in run]
        moved += [_moved(anchor, -record.length) for run in runs[first + 1:] for anchor in run]
    else:
        moved = [anchor for run in runs[first:] for anchor in run]
        moved += [_moved(anchor, record.length) for run in runs[:first] for anchor in run]
    framed = _Cluster(record, tuple(moved))
    return framed, _split_runs(framed, jump, reverse), True


@dataclass(frozen=True)
class _Piece:
    """A window before its column: a span of one record (``start > end`` across the origin) and its anchors."""

    cluster_rank: int
    cluster: _Cluster
    start: int
    end: int
    anchors: tuple[_Anchor, ...]
    reverse: bool

    @property
    def across_origin(self) -> bool:
        return self.start > self.end


def _run_spans(runs: Sequence[Sequence[_Anchor]], flank: int, closed_length: int | None) -> list[tuple[int, int]]:
    """Each run's span in the cluster frame.

    Runs of one cluster abut at the middle of the gap between them, with no
    flank there; the runs of a closed circle (``closed_length``) also abut
    across the frame's seam.
    """

    lows = [min(anchor.low for anchor in run) for run in runs]
    highs = [max(anchor.high for anchor in run) for run in runs]
    cuts = [(highs[index] + lows[index + 1]) // 2 for index in range(len(runs) - 1)]
    start, end = lows[0] - flank, highs[-1] + flank
    if closed_length is not None:
        seam = (highs[-1] - closed_length + lows[0]) // 2
        start, end = seam + 1, seam + closed_length
    edges = [start - 1, *cuts, end]
    return [(edges[index] + 1, edges[index + 1]) for index in range(len(runs))]


def _record_span(start: int, end: int, record: LocateRecord) -> tuple[int, int]:
    """A span of the cluster frame as record coordinates (rule 7).

    A circular record has no ends: the span wraps, crossing the origin when
    it must. A linear record clamps it at its ends.
    """

    length = record.length
    if record.circular:
        if end - start + 1 >= length:
            low = (start - 1) % length + 1
            return low, low - 1 if low > 1 else length
        return (start - 1) % length + 1, (end - 1) % length + 1
    return max(1, start), min(length, end)


def _pieces(rank: int, cluster: _Cluster, runs: Sequence[Sequence[_Anchor]], flank: int, closed: bool) -> list[_Piece]:
    """The windows of one cluster: one per part, each oriented by its own genes (rule 7)."""

    record = cluster.record
    pieces = []
    for run, (start, end) in zip(runs, _run_spans(runs, flank, record.length if closed else None)):
        # A cut between parts is not a record end: only the outer edges take one.
        low, high = _record_span(start, end, record)
        pieces.append(_Piece(rank, cluster, low, high, tuple(run), _reversed(run)))
    return pieces


def _gene_status(anchor: _Anchor) -> GeneStatus:
    record = anchor.locus.record
    at_end = not record.circular and (anchor.start <= 1 or anchor.end >= record.length)
    complete = anchor.locus.query_coverage >= COMPLETE_QUERY_COVERAGE and not at_end
    return "complete" if complete else "partial"


def _window_genes(anchors: Sequence[_Anchor]) -> tuple[WindowGene, ...]:
    return tuple(
        WindowGene(
            protein_id=anchor.locus.gene.protein_id,
            order=anchor.order,
            window_index=anchor.locus.gene.window_index,
            start=anchor.locus.start,
            end=anchor.locus.end,
            strand=anchor.locus.strand,
            query_coverage=anchor.locus.query_coverage,
            bits=anchor.locus.bits,
            status=_gene_status(anchor),
        )
        for anchor in sorted(anchors, key=lambda item: (item.order, item.locus.start, item.locus.end))
    )


def _median_order(piece: _Piece) -> tuple[int, int]:
    orders = sorted(anchor.order for anchor in piece.anchors)
    return orders[len(orders) // 2], orders[0]


def _select(
    ranked: Sequence[tuple[_Cluster, list[_Piece]]],
    reference: LocateReference,
    options: LocateOptions,
) -> tuple[dict[int, CandidateReason | None], int]:
    """Which clusters become windows (rule 5), and how many windows ``max_windows`` left out."""

    counts = _window_gene_counts(reference)
    covered: set[int] = set()
    chosen: dict[int, CandidateReason | None] = {}
    windows = 0
    capped = 0
    for rank, (cluster, pieces) in enumerate(ranked):
        needed = _genes_needed(counts[_majority_window(cluster)], options)
        if len(cluster.genes) < needed:
            chosen[rank] = "TOO_FEW_GENES"
        elif len(cluster.genes - covered) < needed:
            chosen[rank] = "COVERED_BY_OTHER_WINDOW"
        elif windows + len(pieces) > options.max_windows:
            chosen[rank] = "MAX_WINDOWS"
            capped += len(pieces)
        else:
            chosen[rank] = None
            covered |= cluster.genes
            windows += len(pieces)
    return chosen, capped


def _gap(first: RecordWindow, second: RecordWindow) -> int | None:
    """The bases left out between two windows of one record (0 when they abut), or None."""

    if first.crosses_origin or second.crosses_origin:
        return None
    if first.end < second.start:
        return second.start - first.end - 1
    if second.end < first.start:
        return first.start - second.end - 1
    return None


def _record_ends(window: RecordWindow, length: int, circular: bool) -> tuple[bool, bool]:
    """Whether ``window`` reaches the start and the end of its record; a circular record has none (D-28)."""

    return (not circular and window.start == 1, not circular and window.end == length)


def _planned_windows(
    pieces: Sequence[_Piece], reasons: Mapping[int, CandidateReason | None], reference: LocateReference
) -> list[PlannedWindow]:
    ordered = sorted(
        pieces,
        key=lambda piece: (*_median_order(piece), piece.cluster.record.index, piece.start, piece.end),
    )
    selected = [piece for piece in ordered if reasons[piece.cluster_rank] is None]
    gaps: dict[int, int | None] = {}
    for piece, following in zip(selected, selected[1:]):
        if piece.cluster.record.index == following.cluster.record.index:
            gaps[id(piece)] = _gap(RecordWindow(piece.start, piece.end), RecordWindow(following.start, following.end))
    by_order = {gene.order: gene for gene in reference.genes}
    planned = []
    for column, piece in enumerate(ordered, start=1):
        record = piece.cluster.record
        genes = [by_order[anchor.order].cds for anchor in piece.anchors]
        window = RecordWindow(piece.start, piece.end, piece.reverse, column)
        touches_start, touches_end = _record_ends(window, record.length, record.circular)
        planned.append(PlannedWindow(
            record_key=record.record_key,
            record_id=record.record_id,
            record_index=record.index,
            record_length=record.length,
            record_circular=record.circular,
            window=window,
            selected=reasons[piece.cluster_rank] is None,
            reason=reasons[piece.cluster_rank],
            cluster=piece.cluster_rank + 1,
            touches_start=touches_start,
            touches_end=touches_end,
            across_origin=piece.across_origin,
            gap_to_next_same_record=gaps.get(id(piece)),
            evidence_bits=math.fsum(anchor.locus.bits for anchor in piece.anchors),
            genes=_window_genes(piece.anchors),
            ref_span=(min(cds.start for cds in genes), max(cds.end for cds in genes)),
        ))
    return planned


def _coverage(reference: LocateReference, windows: Sequence[PlannedWindow]) -> GeneCoverage:
    status: dict[str, GeneStatus] = {}
    for window in windows:
        for found in window.genes:
            if status.get(found.protein_id) != "complete":
                status[found.protein_id] = found.status
    groups: list[list[ReferenceGene]] = []
    previous: ReferenceGene | None = None
    for gene in reference.genes:
        if status.get(gene.protein_id) == "complete":
            previous = None
            continue
        if previous is not None and previous.window_index == gene.window_index:
            groups[-1].append(gene)
        else:
            groups.append([gene])
        previous = gene
    runs = []
    for group in groups:
        # Genes in drawn order: a reversed window draws its last source gene first.
        first, last = (group[-1], group[0]) if reference.windows[group[0].window_index].reverse_complement else (group[0], group[-1])
        runs.append((first.cds.start, last.cds.end))
    ids = [gene.protein_id for gene in reference.genes]
    return GeneCoverage(
        complete=tuple(item for item in ids if status.get(item) == "complete"),
        partial=tuple(item for item in ids if status.get(item) == "partial"),
        absent=tuple(item for item in ids if item not in status),
        reference_not_found=tuple(runs),
    )


def _file_warnings(
    target: LocateTarget,
    windows: Sequence[PlannedWindow],
    coverage: GeneCoverage,
    tally: ComparisonRowTally,
    capped: int,
    options: LocateOptions,
) -> tuple[LocateWarning, ...]:
    warnings: list[LocateWarning] = []
    if not any(record.cds for record in target.records):
        warnings.append(LocateWarning("TARGET_NOT_ANNOTATED"))
    if not windows:
        warnings.append(LocateWarning(
            "FILE_NO_WINDOWS", (("rows_read", tally.rows_read), ("rows_kept", tally.rows_kept))
        ))
        return tuple(warnings)
    if capped:
        warnings.append(LocateWarning("MAX_WINDOWS", (("max_windows", options.max_windows), ("windows", capped))))
    if coverage.partial:
        warnings.append(LocateWarning("GENES_PARTIAL", (("count", len(coverage.partial)),)))
    if coverage.absent:
        warnings.append(LocateWarning("GENES_ABSENT", (("count", len(coverage.absent)),)))
    return tuple(warnings)


def _plan_file(
    reference: LocateReference, target: LocateTarget, result: LocateSearchResult, options: LocateOptions
) -> FilePlan:
    tally, hsps = _hsps(result, reference, target, options.thresholds)
    clusters = _clusters(_anchors(_loci(hsps, options.program), reference, options), options.merge_distance)
    clusters.sort(key=lambda cluster: (-cluster.bits, cluster.record.index, cluster.low))
    counts = _window_gene_counts(reference)
    ranked: list[tuple[_Cluster, list[_Piece]]] = []
    for rank, cluster in enumerate(clusters):
        reverse = _reversed(cluster.anchors)
        jump = _jump(counts[_majority_window(cluster)], options)
        framed, runs, closed = _circular_frame(
            cluster, _split_runs(cluster, jump, reverse), jump, reverse, options.merge_distance
        )
        ranked.append((cluster, _pieces(rank, framed, runs, options.flank, closed)))
    reasons, capped = _select(ranked, reference, options)
    planned = _planned_windows([piece for _cluster, pieces in ranked for piece in pieces], reasons, reference)
    windows = tuple(window for window in planned if window.selected)
    coverage = _coverage(reference, windows)
    return FilePlan(
        file_key=target.file_key,
        windows=windows,
        candidates=tuple(window for window in planned if not window.selected),
        coverage=coverage,
        tally=tally,
        warnings=_file_warnings(target, windows, coverage, tally, capped, options),
    )


def _failed_file(reference: LocateReference, result: LocateSearchResult) -> FilePlan:
    ids = tuple(gene.protein_id for gene in reference.genes)
    return FilePlan(
        file_key=result.file_key,
        windows=(),
        candidates=(),
        coverage=GeneCoverage((), (), ids, ()),
        tally=ComparisonRowTally(0, 0, ()),
        warnings=(LocateWarning("FILE_SEARCH_FAILED", (("error", str(result.error)),)),),
    )


def plan_locate(
    reference: LocateReference,
    targets: Sequence[LocateTarget],
    results: Sequence[LocateSearchResult],
    options: LocateOptions,
    *,
    runtime: LosatRuntimeRecord | None = None,
) -> LocatePlan:
    """The windows of each target File where the reference windows' genes lie.

    ``results`` holds one search result per target, by ``file_key``. The plan
    keeps the targets' order. ``runtime`` records the runtime that searched.
    """

    by_key = {result.file_key: result for result in results}
    _require(
        len(by_key) == len(results) and sorted(by_key) == sorted(target.file_key for target in targets),
        "Locate needs one search result per target File.",
    )
    files = tuple(
        _failed_file(reference, by_key[target.file_key])
        if by_key[target.file_key].error is not None
        else _plan_file(reference, target, by_key[target.file_key], options)
        for target in targets
    )
    warnings: tuple[LocateWarning, ...] = ()
    if reference.cut_genes:
        warnings = (LocateWarning("REFERENCE_GENE_CUT", (("count", len(reference.cut_genes)),)),)
    return LocatePlan(reference, options, runtime, files, warnings)


# The drawing a plan gives: what the records table, the comparisons table and the new drawing read.


@dataclass(frozen=True)
class LocateDrawnRecord:
    """One record of a row of the drawing: its windows in column order, numbered 1..k over the row's windows.

    ``window_genes`` are the protein IDs of the reference genes each window holds, in the order of ``windows``.
    """

    file_key: str
    record_key: str
    record_id: str
    record_length: int
    windows: tuple[RecordWindow, ...]
    window_genes: tuple[frozenset[str], ...]

    @property
    def whole(self) -> bool:
        """Whether its one window is the whole record as stored, drawn as the record itself."""

        return len(self.windows) == 1 and (self.windows[0].start, self.windows[0].end) == (1, self.record_length)


@dataclass(frozen=True)
class LocateDrawnRow:
    """One row of the drawing: the reference windows, or the selected windows of one target File."""

    file_key: str
    records: tuple[LocateDrawnRecord, ...]


# A window of the drawing as (row, record, window) indexes, 0-based.
LocateWindowPlace = tuple[int, int, int]


@dataclass(frozen=True)
class LocateDrawing:
    """The rows of the drawing a plan gives, and the window pairs it compares.

    The first row is the reference; each target File with a selected window
    follows in plan order. ``pairs`` are the windows of adjacent rows that
    share a reference gene (PD-OI-123), upper row first, in row, record and
    window order.
    """

    rows: tuple[LocateDrawnRow, ...]
    pairs: tuple[tuple[LocateWindowPlace, LocateWindowPlace], ...]


_DrawnWindow = tuple[str, str, int, RecordWindow, frozenset[str]]


def _drawn_row(file_key: str, windows: Sequence[_DrawnWindow]) -> LocateDrawnRow:
    """A row from its windows in column order, each as (record key, record ID, length, window, genes)."""

    grouped: dict[str, list[tuple[int, _DrawnWindow]]] = {}
    for column, item in enumerate(windows, start=1):
        grouped.setdefault(item[0], []).append((column, item))
    return LocateDrawnRow(
        file_key=file_key,
        records=tuple(
            LocateDrawnRecord(
                file_key=file_key,
                record_key=key,
                record_id=items[0][1][1],
                record_length=items[0][1][2],
                windows=tuple(
                    RecordWindow(item[3].start, item[3].end, item[3].reverse_complement, column)
                    for column, item in items
                ),
                window_genes=tuple(item[4] for _column, item in items),
            )
            for key, items in grouped.items()
        ),
    )


def locate_drawing(plan: LocatePlan) -> LocateDrawing:
    """The drawing of ``plan``: the reference windows as picked, then each File's selected windows."""

    reference = plan.reference
    genes_by_window: dict[int, frozenset[str]] = {}
    for gene in reference.genes:
        genes_by_window[gene.window_index] = genes_by_window.get(gene.window_index, frozenset()) | {gene.protein_id}
    rows = [
        _drawn_row(
            reference.file_key,
            [
                (reference.record_key, reference.record_id, reference.record_length, window, genes_by_window[index])
                for index, window in enumerate(reference.windows)
            ],
        )
    ]
    rows.extend(
        _drawn_row(
            file_plan.file_key,
            [
                (
                    planned.record_key,
                    planned.record_id,
                    planned.record_length,
                    planned.window,
                    frozenset(gene.protein_id for gene in planned.genes),
                )
                for planned in file_plan.windows
            ],
        )
        for file_plan in plan.files
        if file_plan.windows
    )
    places = [
        [((row, record_index, window_index), genes)
         for record_index, record in enumerate(drawn.records)
         for window_index, genes in enumerate(record.window_genes)]
        for row, drawn in enumerate(rows)
    ]
    pairs = tuple(
        (upper, lower)
        for row in range(len(rows) - 1)
        for upper, upper_genes in places[row]
        for lower, lower_genes in places[row + 1]
        if upper_genes & lower_genes
    )
    return LocateDrawing(tuple(rows), pairs)


# Inputs from records and the search plan.


def _strand(value: int | None) -> Strand:
    return "-" if value == -1 else "+"


def locate_target(file_key: str, records: Sequence[SeqRecord], record_keys: Sequence[str]) -> LocateTarget:
    """A target File from its records: lengths, topology and CDS location parts."""

    _require(len(records) == len(record_keys), "Locate needs one record key per target record.")
    located = []
    for index, (record, key) in enumerate(zip(records, record_keys)):
        extents = sorted(
            {
                CdsExtent(int(part.start) + 1, int(part.end), _strand(part.strand))
                for feature in record.features
                if feature.type == "CDS"
                for part in feature.location.parts
            },
            key=lambda extent: (extent.start, extent.end, extent.strand),
        )
        located.append(LocateRecord(
            record_key=str(key),
            record_id=str(record.id),
            index=index,
            length=len(record.seq),
            circular=str(record.annotations.get("topology", "")).lower() == "circular",
            cds=tuple(extents),
        ))
    return LocateTarget(str(file_key), tuple(located))


@dataclass(frozen=True)
class LocateQuery:
    """The reference of a search and the FASTA text of its queries."""

    reference: LocateReference
    fasta: str


def locate_query(
    record: SeqRecord,
    *,
    file_key: str,
    record_key: str,
    windows: Sequence[RecordWindow],
    program: LocateProgramName,
) -> LocateQuery:
    """The reference genes of ``record`` in ``windows`` and their query FASTA for ``program``.

    The proteins and their IDs are those LOSATP uses
    (``extract_web_stable_cds_proteins``); a nucleotide query is the CDS's
    coding sequence.
    """

    proteins = extract_web_stable_cds_proteins([record], record_instance_keys=[str(record_key)]).proteins_by_record[0]
    cds = [
        ReferenceCds(
            protein_id=protein.protein_id,
            locus_tag=protein.locus_tag,
            gene=protein.gene,
            product=protein.product,
            start=int(protein.arc[0]) + 1,
            end=int(protein.arc[1]),
            strand=_strand(protein.strand),
            length_aa=int(protein.protein_length),
            length_nt=sum(int(end) - int(start) for start, end, _strand_value in protein.feature_hash_parts),
        )
        for protein in proteins
    ]
    reference = locate_reference(
        file_key=file_key,
        record_key=record_key,
        record_id=str(record.id),
        record_length=len(record.seq),
        circular=str(record.annotations.get("topology", "")).lower() == "circular",
        windows=windows,
        cds=cds,
    )
    queried = {gene.protein_id for gene in reference.genes}
    chosen = [protein for protein in proteins if protein.protein_id in queried]
    if LOSAT_PROGRAMS[program].molecules[0] == "protein":
        return LocateQuery(reference, proteins_to_fasta(chosen))
    texts = []
    for protein in chosen:
        if protein.origin_arc is not None:
            coding = record.seq[protein.origin_arc[0]:] + record.seq[:protein.origin_arc[1]]
        else:
            coding = record.seq[0:0]
            for start, end, _strand_value in sorted(protein.feature_hash_parts):
                coding += record.seq[int(start):int(end)]
        sequence = coding.reverse_complement() if protein.strand == -1 else coding
        texts.append(losat_fasta_text(protein.protein_id, str(sequence)))
    return LocateQuery(reference, "".join(texts))


@dataclass(frozen=True)
class LocateSearchJob:
    """One search: the queries against every record of one target File.

    ``raw_key`` is the job's raw-cache key; ``batch`` maps the search IDs back
    (:func:`locate_search_rows`).
    """

    file_key: str
    query_fasta: str
    subject_fasta: str
    raw_key: str
    batch: LosatBatch


@dataclass(frozen=True)
class LocateSearchPlan:
    program: LocateProgramName
    # The search options every job runs with, and the same in the raw-cache form.
    search_args: LosatSearchArgs
    cache_args: tuple[str, ...]
    jobs: tuple[LocateSearchJob, ...]


def locate_search_args(program: LocateProgramName) -> LosatSearchArgs:
    """The search options of a locate search with ``program``: the runtime defaults, and the default
    ``losatn`` task (``losatn`` needs one)."""

    _require(program in LOCATE_PROGRAMS, f"Unknown locate program {program!r}.")
    return LosatSearchArgs(task=DEFAULT_LOSATN_TASK if "task" in LOSAT_PROGRAMS[program].required_options else None)


def plan_locate_search(
    query: LocateQuery,
    targets: Sequence[tuple[str, Sequence[SeqRecord], Sequence[str]]],
    *,
    program: LocateProgramName,
) -> LocateSearchPlan:
    """One job per target File through the shared job owner (``plan_losat_jobs``).

    ``targets`` holds ``(file_key, records, record_keys)`` per File. The
    queries are one source; each File is one source searched whole, with
    :func:`locate_search_args`.
    """

    search_args = locate_search_args(program)
    spec = LOSAT_PROGRAMS[program]
    source_ids: list[str] = [f"reference:{query.reference.file_key}"]
    uids: list[str] = [query.reference.record_key]
    fastas: list[str] = [query.fasta]
    file_keys: list[str] = [query.reference.file_key]
    for file_key, records, record_keys in targets:
        _require(len(records) == len(record_keys), "Locate needs one record key per target record.")
        for record, key in zip(records, record_keys):
            source_ids.append(f"target:{file_key}")
            uids.append(str(key))
            fastas.append(losat_fasta_text(str(record.id), str(record.seq)))
            file_keys.append(str(file_key))
    cache_args = tuple(losat_cache_args(spec, search_args))
    jobs = plan_losat_jobs(
        source_ids=source_ids,
        specs=[(0, index) for index in range(1, len(source_ids))],
        build_args=lambda _query, _subject: cache_args,
    )
    batches = prepare_losat_batches(jobs, uids=uids, record_fasta=lambda index: fastas[index])
    return LocateSearchPlan(
        program=program,
        search_args=search_args,
        cache_args=cache_args,
        jobs=tuple(
            LocateSearchJob(
                file_key=file_keys[batch.subject.indexes[0]],
                query_fasta=batch.query.fasta,
                subject_fasta=batch.subject.fasta,
                raw_key=nucleotide_losat_cache_key(
                    program=spec.search,
                    args=cache_args,
                    query_hash=batch.query.hash,
                    subject_hash=batch.subject.hash,
                    search_context=batch.search_context,
                ),
                batch=batch,
            )
            for batch in batches
        ),
    )


def locate_search_rows(job: LocateSearchJob, text: str) -> LocateSearchResult:
    """A job's raw rows with the original record IDs, in runtime order per record."""

    split = split_losat_batch_result(text, job.batch, job.batch.job.specs)
    return LocateSearchResult(job.file_key, "".join(split[spec] for spec in job.batch.job.specs))


# One locate run from input files: what the command line and the Web's Worker share; each runs the jobs itself.


class LocateInputFile(NamedTuple):
    """One input File: its key in the plan (the command line's path, the Web's File key) and its input paths:
    GenBank, GFF3 with FASTA, or FASTA alone."""

    key: str
    gbk: str | None
    gff: str | None
    fasta: str | None


def read_locate_file(file: LocateInputFile) -> list[SeqRecord]:
    """The records of ``file``, nested GFF3 features (a CDS under an mRNA) read as the drawing reads them."""

    if file.gbk:
        return load_gbks([file.gbk])
    if file.gff and file.fasta:
        return [filter_features_by_type(record, None) for record in parse_gff_fasta(file.gff, file.fasta)]
    _require(bool(file.fasta), f"Locate File {file.key!r} names no input.")
    return load_fasta([str(file.fasta)])


@dataclass(frozen=True)
class LocateSearch:
    """The queries, the targets and the search plan of one locate run (:func:`prepare_locate_search`)."""

    query: LocateQuery
    targets: tuple[LocateTarget, ...]
    plan: LocateSearchPlan


def _reference_record(records: Sequence[SeqRecord], selector: str | None, file_key: str) -> SeqRecord:
    if selector is None:
        _require(
            len(records) == 1,
            f"The reference file {file_key} holds {len(records)} records; choose one with --reference_record.",
        )
        return records[0]
    try:
        [record] = select_record(list(records), parse_record_selector(selector))
    except ValueError as exc:
        raise ValidationError(f"The reference file {file_key} has no record {selector!r}: {exc}", diagnostic=_INVALID) from exc
    return record


def prepare_locate_search(
    reference: LocateInputFile,
    targets: Sequence[LocateInputFile],
    *,
    record: str | None,
    record_key: str,
    windows: Sequence[RecordWindow],
    program: LocateProgramName,
    target_record_keys: Mapping[str, Mapping[int, str]] | None = None,
) -> LocateSearch:
    """Read the reference record (``record``: an ID or ``#index``, required when the file holds several) and its
    windows (none: the whole record), and plan one search per target File.

    ``record_key`` is the reference record's key in the drawing it was chosen from, so the new drawing carries
    that record's edits (:func:`gbdraw.api.locate_drawing.derive_locate_request`). ``target_record_keys`` gives,
    by File key and 1-based record index, the keys that a target record has in that drawing; any other record of a
    target File is ``FILE#index``.
    """

    chosen = _reference_record(read_locate_file(reference), record, reference.key)
    query = locate_query(
        chosen, file_key=reference.key, record_key=record_key,
        windows=windows or [RecordWindow(1, len(chosen.seq))], program=program,
    )
    keys = target_record_keys or {}
    files = [
        (
            file.key,
            records,
            [keys.get(file.key, {}).get(index, f"{file.key}#{index}") for index in range(1, len(records) + 1)],
        )
        for file in targets
        for records in [read_locate_file(file)]
    ]
    return LocateSearch(
        query=query,
        targets=tuple(locate_target(key, records, keys) for key, records, keys in files),
        plan=plan_locate_search(query, files, program=program),
    )


def locate_results(search: LocateSearch, outcomes: Sequence[LocateSearchResult]) -> tuple[LocateSearchResult, ...]:
    """Each job's runtime output (the raw text as the runtime wrote it, or the error that stopped it) with the
    original record IDs."""

    jobs = {job.file_key: job for job in search.plan.jobs}
    return tuple(
        outcome if outcome.error is not None else locate_search_rows(jobs[outcome.file_key], outcome.text)
        for outcome in outcomes
    )


def finish_locate(
    search: LocateSearch,
    results: Sequence[LocateSearchResult],
    options: LocateOptions,
    *,
    runtime: LosatRuntimeRecord | None = None,
) -> LocatePlan:
    """The plan of ``search`` from its results with original record IDs (:func:`locate_results`)."""

    return plan_locate(search.query.reference, search.targets, results, options, runtime=runtime)


__all__ = [
    "COMPLETE_QUERY_COVERAGE",
    "CandidateReason",
    "CdsExtent",
    "FilePlan",
    "GeneCoverage",
    "GeneStatus",
    "LOCATE_PROGRAMS",
    "LocateWindowChoice",
    "LocateWindowRefusal",
    "LocateWindowRefusalReason",
    "check_record_window",
    "choose_locate_windows",
    "LocateDrawing",
    "LocateInputFile",
    "LocateSearch",
    "finish_locate",
    "locate_results",
    "prepare_locate_search",
    "read_locate_file",
    "LocateDrawnRecord",
    "LocateDrawnRow",
    "LocateOptions",
    "LocatePlan",
    "LocatePlanJson",
    "LocateProgramName",
    "LocateQuery",
    "LocateRecord",
    "LocateWindowPlace",
    "LocateReference",
    "LocateSearchJob",
    "LocateSearchPlan",
    "LocateSearchResult",
    "LocateTarget",
    "LocateWarning",
    "LocateWarningCode",
    "PlannedWindow",
    "ReferenceCds",
    "ReferenceGene",
    "WindowGene",
    "locate_default_thresholds",
    "locate_drawing",
    "locate_query",
    "locate_reference",
    "locate_search_args",
    "locate_search_rows",
    "locate_target",
    "plan_locate",
    "plan_locate_search",
]
