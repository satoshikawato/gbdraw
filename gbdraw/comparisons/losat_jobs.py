#!/usr/bin/env python
# coding: utf-8

"""The LOSAT job plan: which searches run, against which database, and under
which raw-cache key.

This is the Python owner of the rules that the Web keeps in
``planLosatSourceJobs``, ``prepareLosatSourceBatches`` and
``splitLosatSourceResult`` (``gbdraw/web/js/app/linear-sources.js``), the
LOSATP record-pair searches of ``buildLosatJobSpecs``
(``gbdraw/web/js/app/linear-comparisons.js``), the nucleotide FASTA text of
``extractLosatFastaFast`` and the nucleotide raw key of
``buildLosatCachePayload`` (``gbdraw/web/js/app/run-analysis.js``). The Web
needs its own copy for the Settings job estimate without Python (CW-01, design
D15); the shared vectors ``tests/fixtures/losat_job_plan_cases.json`` and
``tests/fixtures/losat_fasta_extraction_cases.json`` hold both copies to the
same result, so a CLI Session reuses its raw cache in the Web and back.

Rules: a source file is one genome. Between two sources the query source
searches the whole subject source; within one source a query searches the
source without itself; a record searches itself only when that self search is
requested. Records whose explicit search args differ never share a batch.
LOSATN, TLOSATX and LOSATP all search these source batches (design D7).
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import re
from typing import Callable, Hashable, Literal, Mapping, Sequence

from Bio.SeqRecord import SeqRecord

from gbdraw.core.record_metadata import _read_coord_map
from gbdraw.exceptions import ValidationError
from gbdraw.io.filenames import unique_filenames

NUCLEOTIDE_LOSAT_CACHE_SCHEMA = 2
LOSAT_OUTFMT = "6"
_FASTA_LINE_WIDTH = 60
_OUTFMT6_COLUMNS = 12
_FASTA_HEADER = re.compile(r"^>(\S+)([^\r\n]*)", re.MULTILINE)

LosatJobScope = Literal["between-sources", "within-source", "self"]
RecordPair = tuple[int, int]


def sha256_text(text: str) -> str:
    """Hex SHA-256 of UTF-8 text (the Web ``hashText``)."""

    return hashlib.sha256(str(text).encode("utf-8")).hexdigest()


def unique_losat_filenames(
    names: Sequence[str], *, reserved: Sequence[str] = ()
) -> tuple[str, ...]:
    """``--losat_output_dir`` TSV names of LOSAT edges or rings, unique in order.

    A repeated ``<stem>.tsv`` becomes ``<stem>.2.tsv``, ``<stem>.3.tsv``, ...
    (:func:`gbdraw.io.filenames.unique_filenames`, the Session resource rule).
    """

    return unique_filenames(
        [name if name.endswith(".tsv") else f"{name or 'losat'}.tsv" for name in names],
        reserved=reserved,
    )


def _json_text(value: object) -> str:
    """``JSON.stringify`` of plain strings, numbers, lists and dicts."""

    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


# Nucleotide FASTA text of one record.


def losat_fasta_text(record_id: str, sequence: str) -> str:
    """One FASTA record as the Web writes LOSAT input: 60 columns, uppercase."""

    text = str(sequence).upper()
    lines = [
        text[start:start + _FASTA_LINE_WIDTH]
        for start in range(0, len(text), _FASTA_LINE_WIDTH)
    ]
    return f">{record_id}\n" + "\n".join(lines) + "\n"


def losat_search_frame_fasta(record: SeqRecord) -> str:
    """FASTA of a displayed record in its search frame.

    The search frame is the frame of ``-b`` tables: the selected (and cropped)
    source record before a display reverse complement (PD-OI-073).
    """

    _base, step = _read_coord_map(record)
    sequence = record.seq if step == 1 else record.seq.reverse_complement()
    return losat_fasta_text(str(record.id), str(sequence))


# Record sources and the LOSATP record-pair searches.


def record_source_paths(record: SeqRecord) -> tuple[str, ...]:
    """The source file paths the request planner recorded; empty for any other record."""

    paths = (getattr(record, "annotations", None) or {}).get("gbdraw_source_paths")
    return tuple(str(path) for path in paths) if paths else ()


def losat_source_ids(records: Sequence[SeqRecord]) -> tuple[Hashable, ...]:
    """The source file of each record: one file is one genome.

    Records read by the request planner carry their source paths; any other
    record is its own source.
    """

    return tuple(
        record_source_paths(record)
        or ("memory", (getattr(record, "annotations", None) or {}).get("gbdraw_input_index", index))
        for index, record in enumerate(records)
    )


def losat_record_uids(records: Sequence[SeqRecord]) -> tuple[str, ...]:
    """The record keys that order multi-record batch sides (Web sequence uids)."""

    return tuple(
        str((getattr(record, "annotations", None) or {}).get("gbdraw_record_key") or f"record-{index + 1}")
        for index, record in enumerate(records)
    )


LosatpSpecMode = Literal["pairwise", "orthogroup", "collinear"]


def losatp_job_specs(
    mode: LosatpSpecMode | str,
    *,
    record_count: int,
    pairs: Sequence[RecordPair] | None = None,
    infer_orthogroups: bool = True,
    search_scope: str = "adjacent",
) -> tuple[RecordPair, ...]:
    """The directed record-pair searches of one LOSATP run (Web ``buildLosatJobSpecs``).

    ``pairs`` are the displayed comparison edges; omitted, they are the
    consecutive records. Pairwise searches each edge. Similarity groups search
    every record against every record, itself included. Collinear searches
    each edge in both directions (every record pair with ``search_scope="all"``)
    and adds the within-record searches when inference is on.
    """

    count = max(0, int(record_count))
    edges = (
        tuple((int(query), int(subject)) for query, subject in pairs)
        if pairs is not None
        else tuple((index, index + 1) for index in range(max(0, count - 1)))
    )
    if mode == "pairwise":
        return edges
    specs: list[RecordPair] = []

    def every_pair(with_self: bool) -> None:
        for query in range(count):
            if with_self:
                specs.append((query, query))
            for subject in range(query + 1, count):
                specs.extend(((query, subject), (subject, query)))

    if mode == "orthogroup":
        every_pair(True)
        return tuple(specs)
    if mode != "collinear":
        raise ValidationError(
            f"Unsupported LOSATP mode: {mode!r}.",
            diagnostic={"code": "COMPARISON_INPUT", "reason": "LOSAT_PLAN"},
        )
    if infer_orthogroups:
        specs.extend((index, index) for index in range(count))
    if str(search_scope) == "all":
        every_pair(False)
    else:
        for query, subject in edges:
            specs.extend(((query, subject), (subject, query)))
    return tuple(specs)


# Job plan.


@dataclass(frozen=True)
class LosatJob:
    """One raw search: query records against a subject database."""

    query_indexes: tuple[int, ...]
    subject_indexes: tuple[int, ...]
    args: tuple[str, ...]
    scope: LosatJobScope
    specs: tuple[RecordPair, ...]


def plan_losat_jobs(
    *,
    source_ids: Sequence[Hashable],
    specs: Sequence[RecordPair],
    build_args: Callable[[int, int], Sequence[str]],
) -> tuple[LosatJob, ...]:
    """Return the source jobs that answer ``specs`` (Web ``planLosatSourceJobs``).

    ``source_ids[i]`` identifies the source file of record ``i``;
    ``build_args(query, subject)`` returns the search args of one record pair.
    """

    records_by_source: dict[Hashable, list[int]] = {}
    for index, source_id in enumerate(source_ids):
        records_by_source.setdefault(source_id, []).append(index)

    def args_of(query: int, subject: int) -> tuple[str, ...]:
        return tuple(str(arg) for arg in build_args(query, subject))

    order: list[tuple[object, ...]] = []
    drafts: dict[tuple[object, ...], dict[str, object]] = {}
    for query_index, subject_index in specs:
        args = args_of(query_index, subject_index)
        query_source = source_ids[query_index]
        subject_source = source_ids[subject_index]
        scope: LosatJobScope = (
            "self"
            if query_index == subject_index
            else "within-source" if query_source == subject_source else "between-sources"
        )
        key = (
            query_source,
            subject_source,
            args,
            scope,
            None if scope == "between-sources" else query_index,
        )
        draft = drafts.get(key)
        if draft is None:
            query_indexes = (
                tuple(
                    index
                    for index in records_by_source[query_source]
                    if args_of(index, subject_index) == args
                )
                if scope == "between-sources"
                else (query_index,)
            )
            subject_indexes = (
                (subject_index,)
                if scope == "self"
                else tuple(
                    index
                    for index in records_by_source[subject_source]
                    if index != query_index and args_of(query_index, index) == args
                )
            )
            draft = {
                "query": query_indexes,
                "subject": subject_indexes,
                "args": args,
                "scope": scope,
                "specs": [],
            }
            drafts[key] = draft
            order.append(key)
        draft["specs"].append((int(query_index), int(subject_index)))
    return tuple(
        LosatJob(
            query_indexes=drafts[key]["query"],  # type: ignore[arg-type]
            subject_indexes=drafts[key]["subject"],  # type: ignore[arg-type]
            args=drafts[key]["args"],  # type: ignore[arg-type]
            scope=drafts[key]["scope"],  # type: ignore[arg-type]
            specs=tuple(drafts[key]["specs"]),  # type: ignore[arg-type]
        )
        for key in order
    )


# Source batches.


@dataclass(frozen=True)
class LosatBatchSide:
    """The FASTA of one batch side; ``ids`` maps a search ID to (record, original ID)."""

    indexes: tuple[int, ...]
    ids: Mapping[str, tuple[int, str]]
    fasta: str
    hash: str


@dataclass(frozen=True)
class LosatBatch:
    job: LosatJob
    query: LosatBatchSide
    subject: LosatBatchSide
    # The searched database is part of raw-cache identity when a side has
    # several records.
    search_context: str | None


def prepare_losat_batches(
    jobs: Sequence[LosatJob],
    *,
    uids: Sequence[str],
    record_fasta: Callable[[int], str],
    protein: bool = False,
) -> tuple[LosatBatch, ...]:
    """Build the FASTA of each job side (Web ``prepareLosatSourceBatches``).

    ``uids`` are the record keys the Web uses as sequence UIDs. Nucleotide
    records of a multi-record side get collision-free search IDs.
    """

    sides: dict[tuple[int, ...], LosatBatchSide] = {}

    def prepare_side(indexes: Sequence[int]) -> LosatBatchSide:
        ordered = tuple(sorted(indexes, key=lambda index: str(uids[index])))
        cached = sides.get(ordered)
        if cached is not None:
            return cached
        ids: dict[str, tuple[int, str]] = {}
        parts: list[str] = []
        for index in ordered:
            fasta = record_fasta(index)
            prefix = (
                f"n_{sha256_text(str(uids[index]))}"
                if not protein and len(ordered) > 1
                else ""
            )
            ordinal = 0

            def rename(match: re.Match[str]) -> str:
                nonlocal ordinal
                original_id, description = match.group(1), match.group(2)
                search_id = f"{prefix}_{ordinal}" if prefix else original_id
                ordinal += 1
                if search_id in ids:
                    raise ValidationError(
                        f"Duplicate LOSAT source identifier: {search_id}",
                        diagnostic={"code": "COMPARISON_INPUT", "reason": "LOSAT_PLAN"},
                    )
                ids[search_id] = (index, original_id)
                return f">{search_id}{description}"

            fasta = _FASTA_HEADER.sub(rename, fasta)
            parts.append(fasta if fasta.endswith("\n") else f"{fasta}\n")
        text = "".join(parts)
        side = LosatBatchSide(ordered, ids, text, sha256_text(text))
        sides[ordered] = side
        return side

    batches: list[LosatBatch] = []
    for job in jobs:
        query = prepare_side(job.query_indexes)
        subject = prepare_side(job.subject_indexes)
        search_context = (
            sha256_text(_json_text([query.hash, subject.hash]))
            if len(query.indexes) > 1 or len(subject.indexes) > 1
            else None
        )
        batches.append(LosatBatch(job, query, subject, search_context))
    return tuple(batches)


def split_losat_batch_result(
    text: str,
    batch: LosatBatch,
    specs: Sequence[RecordPair],
) -> dict[RecordPair, str]:
    """Split one batch result into record-pair TSVs with original record IDs.

    Rows keep the runtime order (Web ``splitLosatSourceResult``).
    """

    rows: dict[RecordPair, list[str]] = {tuple(spec): [] for spec in specs}  # type: ignore[misc]
    for line_number, line in enumerate(re.split(r"\r?\n", str(text)), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        columns = line.split("\t")
        query = batch.query.ids.get(columns[0])
        subject = batch.subject.ids.get(columns[1]) if len(columns) > 1 else None
        if len(columns) != _OUTFMT6_COLUMNS or query is None or subject is None:
            raise ValidationError(
                "LOSAT result contains malformed or unrecognized record endpoints "
                f"(output line {line_number}).",
                diagnostic={"code": "LOSAT_RUNTIME", "reason": "OUTPUT", "row": line_number},
            )
        target = rows.get((query[0], subject[0]))
        if target is not None:
            columns[0] = query[1]
            columns[1] = subject[1]
            target.append("\t".join(columns))
    return {pair: "".join(f"{line}\n" for line in lines) for pair, lines in rows.items()}


# Raw-cache identity.


def nucleotide_losat_cache_key(
    *,
    program: str,
    args: Sequence[str],
    query_hash: str,
    subject_hash: str,
    flow: str | None = None,
    search_context: str | None = None,
) -> str:
    """Raw key of one nucleotide record pair (Web ``buildLosatCachePayload``).

    ``program`` is the search program (``blastn`` or ``tblastx``); the hashes
    are of each record's search-frame FASTA.
    """

    payload: dict[str, object] = {
        "cacheSchema": NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
        "program": str(program),
        "outfmt": LOSAT_OUTFMT,
        "args": [str(arg) for arg in args],
        "queryCanonicalHash": str(query_hash),
        "subjectCanonicalHash": str(subject_hash),
    }
    if flow:
        payload["flow"] = str(flow)
    if search_context:
        payload["searchContext"] = str(search_context)
    return sha256_text(_json_text(payload))


__all__ = [
    "LOSAT_OUTFMT",
    "LosatBatch",
    "LosatBatchSide",
    "LosatJob",
    "LosatpSpecMode",
    "NUCLEOTIDE_LOSAT_CACHE_SCHEMA",
    "losat_fasta_text",
    "losat_record_uids",
    "losat_search_frame_fasta",
    "losat_source_ids",
    "losatp_job_specs",
    "nucleotide_losat_cache_key",
    "plan_losat_jobs",
    "prepare_losat_batches",
    "sha256_text",
    "split_losat_batch_result",
    "unique_losat_filenames",
]
