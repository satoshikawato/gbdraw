#!/usr/bin/env python
# coding: utf-8

"""Resolve Linear LOSATN / TLOSATX intent into comparisons and raw cache entries.

The request planner calls :func:`resolve_linear_nucleotide_losat` once the
records are loaded. It plans the searches with :mod:`gbdraw.comparisons.losat_jobs`,
runs each source job once through :mod:`gbdraw.comparisons.losat_runtime`, and
returns the options with ``LinearComparison`` values (search-frame rows) in
place of the search intent, plus Web-identical raw cache entries. The resolved
request carries no search intent, so Session replay needs no LOSAT.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import io
import logging
import re
from typing import Any, Callable, Hashable, Mapping, Sequence, cast

from Bio.SeqRecord import SeqRecord

from gbdraw.api.options import LinearDiagramOptions, LosatSearchOptions
from gbdraw.comparisons.losat_jobs import (
    LOSAT_OUTFMT,
    NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
    GenomeSearchSpaces,
    genome_file_names,
    genome_search_spaces,
    losat_search_frame_fasta,
    nucleotide_losat_cache_key,
    plan_losat_jobs,
    prepare_losat_batches,
    sha256_text,
    source_file_name,
    split_losat_batch_result,
)
from gbdraw.comparisons.losat_runtime import (
    LOSAT_PROGRAMS,
    NUCLEOTIDE_LOSAT_PROGRAMS,
    LosatRawCache,
    LosatRuntimeRecord,
    LosatSearchArgs,
    losat_cache_args,
    losat_takes_genome_scale,
    require_losat_task_support,
    run_losat_search,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import read_comparison_table
from gbdraw.linear_comparison import LinearComparison

logger = logging.getLogger(__name__)

def _plan_error(message: str) -> ValidationError:
    return ValidationError(
        message, diagnostic={"code": "COMPARISON_INPUT", "reason": "LOSAT_PLAN"}
    )


def _adjacent_row_pairs(rows_by_record: Sequence[int]) -> list[tuple[int, int]]:
    """Every record of a row against every record of the next row (Web adjacent)."""

    rows = sorted(set(rows_by_record))
    by_row = {row: [i for i, value in enumerate(rows_by_record) if value == row] for row in rows}
    pairs: list[tuple[int, int]] = []
    for upper, lower in zip(rows, rows[1:]):
        pairs.extend((query, subject) for query in by_row[upper] for subject in by_row[lower])
    return pairs


def _validated_pairs(
    pairs: Sequence[tuple[int, int]],
    rows_by_record: Sequence[int],
) -> list[tuple[int, int]]:
    count = len(rows_by_record)
    result: list[tuple[int, int]] = []
    for query, subject in pairs:
        if not (0 <= query < count and 0 <= subject < count):
            raise _plan_error(
                f"LOSAT pair ({query}, {subject}) names a record outside the "
                f"{count} loaded record(s)."
            )
        if query == subject or abs(rows_by_record[query] - rows_by_record[subject]) != 1:
            raise _plan_error(
                f"LOSAT pair ({query}, {subject}) must connect records in adjacent rows."
            )
        if (query, subject) not in result:
            result.append((query, subject))
    return result


def _record_gencodes(
    search: LosatSearchOptions,
    input_indexes: Sequence[int],
) -> list[int | None]:
    gencodes = tuple(search.record_gencodes)
    if not gencodes:
        return [None] * len(input_indexes)
    if len(gencodes) == 1:
        return [gencodes[0]] * len(input_indexes)
    input_count = max(input_indexes) + 1 if input_indexes else 0
    if len(gencodes) != input_count:
        raise ValidationError(
            f"Give one TLOSATX translation table for all records or one per record "
            f"input ({input_count}); got {len(gencodes)}.",
            diagnostic={"code": "COMPARISON_INPUT", "field": "record_gencodes"},
        )
    return [gencodes[index] for index in input_indexes]


def _filename_label(text: str) -> str:
    dotted = re.sub(r"\.+", ".", re.sub(r"[\s/]+", ".", text)).strip(".")
    return re.sub(r"[^\w.-]+", "_", dotted).strip("_")


def losat_edge_filename(
    query_label: str,
    subject_label: str,
    program: str,
    query_window: int | None = None,
    subject_window: int | None = None,
) -> str:
    """Raw TSV name of one edge, as the Web names it (``*.losatn.tsv``).

    A side drawn in several windows adds its window's number (``w<n>``, joined by ``-``) before ``.tsv``, so the
    window pairs of one record pair get names of their own (OV-2742).
    """

    left = _filename_label(query_label) or "seq_1"
    right = _filename_label(subject_label) or "seq_2"
    windows = "-".join(f"w{window}" for window in (query_window, subject_window) if window is not None)
    return f"{left}.{right}.{program}{'.' + windows if windows else ''}.tsv"


@dataclass(frozen=True)
class NucleotideLosatRecord:
    """One drawn record a Linear nucleotide search reads, with what names and groups it.

    ``source_id`` identifies its source file (one genome) and ``source_length`` is the length of every record of that
    file (``source_genome_length``); ``record_key`` is the key the Web uses as its sequence UID;
    ``label`` names its raw TSVs; ``window`` is its window number when it is one of several windows of a record, else
    ``None``; ``input_index`` is the record input it came from (its translation table).
    """

    record: SeqRecord
    source_id: object
    source_length: int
    record_key: str
    label: str
    window: int | None
    input_index: int


def _pair_search_args(
    search: LosatSearchOptions,
    *,
    records: Sequence[SeqRecord],
    rows_by_record: Sequence[int],
    source_ids: Sequence[object],
    source_lengths: Sequence[int],
    input_indexes: Sequence[int],
    pairs: Sequence[tuple[int, int]],
) -> Callable[[int, int], LosatSearchArgs]:
    """The search options of each record pair, by record indexes.

    The files of a row are one genome, the length of every record of them
    counted (``source_lengths``), and a search between two rows passes the
    genome pair's search space as ``-searchsp``
    (:func:`gbdraw.comparisons.losat_jobs.genome_search_spaces`). A runtime
    without ``-searchsp`` searches as before, with a warning.
    """

    gencodes = _record_gencodes(search, input_indexes)
    space = _genome_search_space(
        search, records=records, rows_by_record=rows_by_record, source_ids=source_ids,
        source_lengths=source_lengths, pairs=pairs,
    )

    def search_args(query: int, subject: int) -> LosatSearchArgs:
        searchsp = space(query, subject)
        if search.program == "losatn":
            return LosatSearchArgs(task=search.losatn_task, searchsp=searchsp)
        return LosatSearchArgs(query_gencode=gencodes[query], db_gencode=gencodes[subject], searchsp=searchsp)

    return search_args


def linear_search_spaces(
    program: str,
    *,
    records: Sequence[SeqRecord],
    rows_by_record: Sequence[Hashable],
    source_ids: Sequence[Hashable],
    source_lengths: Sequence[int],
) -> GenomeSearchSpaces:
    """The search space each record-pair search of a LOSATN or TLOSATX ``program`` passes LOSAT, by record indexes:
    the genome pair's (:func:`gbdraw.comparisons.losat_jobs.genome_search_spaces`), each record searched whole as
    drawn (a crop or window is its own record), in the program's letters."""

    return genome_search_spaces(
        rows=rows_by_record,
        source_ids=source_ids,
        source_lengths=source_lengths,
        searched_lengths=[len(record) for record in records],
        letters=LOSAT_PROGRAMS[program].space_letters,
    )


def _genome_search_space(
    search: LosatSearchOptions,
    *,
    records: Sequence[SeqRecord],
    rows_by_record: Sequence[int],
    source_ids: Sequence[object],
    source_lengths: Sequence[int],
    pairs: Sequence[tuple[int, int]],
) -> Callable[[int, int], int | None]:

    spec = LOSAT_PROGRAMS[search.program]
    space = linear_search_spaces(
        search.program, records=records, rows_by_record=rows_by_record, source_ids=source_ids,
        source_lengths=source_lengths,
    )
    names = genome_file_names(rows_by_record, source_ids, source_file_name)
    for message in dict.fromkeys(space.per_file_warning(query, subject, names) for query, subject in pairs):
        if message is not None:
            logger.warning("%s", message)
    if all(space(query, subject) is None for query, subject in pairs) or losat_takes_genome_scale(
        spec, losat_bin=search.runtime.losat_executable, ncbi_blast_bin=search.runtime.ncbi_blast_executable
    ):
        return space
    return lambda _query, _subject: None


def resolve_linear_nucleotide_losat(
    options: LinearDiagramOptions,
    *,
    records: Sequence[NucleotideLosatRecord],
    rows_by_record: Sequence[int],
) -> tuple[LinearDiagramOptions, tuple[Mapping[str, Any], ...]]:
    """Run the requested nucleotide searches over the drawn ``records`` and return resolved options."""

    search = options.losat_search
    if search is None or search.program not in NUCLEOTIDE_LOSAT_PROGRAMS:
        return options, ()
    spec = LOSAT_PROGRAMS[search.program]
    if options.blast_files:
        raise ValidationError(
            f"{search.program.upper()} cannot be combined with blast_files (-b/--blast); "
            "use a comparisons table with source=losat and source=table rows.",
            diagnostic={"code": "COMPARISON_INPUT", "reason": "LOSAT_PLAN", "field": "blast"},
        )
    if len(records) < 2:
        raise _plan_error(f"{search.program.upper()} needs at least two records.")
    pairs = _validated_pairs(
        search.pairs if search.pairs is not None else _adjacent_row_pairs(rows_by_record),
        rows_by_record,
    )
    if not pairs:
        raise _plan_error(
            f"{search.program.upper()} found no record pair to compare: place records "
            "in two or more rows, or list pairs."
        )
    search_args = _pair_search_args(
        search,
        records=[item.record for item in records],
        rows_by_record=rows_by_record,
        source_ids=[item.source_id for item in records],
        source_lengths=[item.source_length for item in records],
        input_indexes=[item.input_index for item in records],
        pairs=pairs,
    )
    runtime = search.runtime
    if search.program == "losatn" and search.losatn_task != "megablast":
        require_losat_task_support(
            str(search.losatn_task),
            losat_bin=runtime.losat_executable,
            ncbi_blast_bin=runtime.ncbi_blast_executable,
        )

    # FASTA extraction once per record (CW-02).
    fasta_by_record: dict[int, str] = {}

    def record_fasta(index: int) -> str:
        if index not in fasta_by_record:
            fasta_by_record[index] = losat_search_frame_fasta(records[index].record)
        return fasta_by_record[index]

    jobs = plan_losat_jobs(
        source_ids=[item.source_id for item in records],
        specs=pairs,
        build_args=lambda query, subject: losat_cache_args(spec, search_args(query, subject)),
    )
    batches = prepare_losat_batches(jobs, uids=[item.record_key for item in records], record_fasta=record_fasta)
    cache = LosatRawCache()
    texts: dict[tuple[int, int], tuple[str, str]] = {}
    searched: dict[tuple[int, int], tuple[dict[str, object], Mapping[str, object] | None]] = {}
    for job_number, batch in enumerate(batches, start=1):
        keys = {
            pair: nucleotide_losat_cache_key(
                program=spec.search,
                args=batch.job.args,
                query_hash=sha256_text(record_fasta(pair[0])),
                subject_hash=sha256_text(record_fasta(pair[1])),
                search_context=batch.search_context,
            )
            for pair in batch.job.specs
        }
        cached = {pair: cache.cached_entry(key) for pair, key in keys.items()}
        if all(entry is not None for entry in cached.values()):
            # The all() above guarantees that no entry is None.
            split = {pair: str(cast("dict[str, object]", entry)["text"]) for pair, entry in cached.items()}
        else:
            logger.info(
                "%s source job %d/%d (%d query x %d subject record(s)).",
                search.program.upper(), job_number, len(batches),
                len(batch.query.indexes), len(batch.subject.indexes),
            )
            runtime_records: list[LosatRuntimeRecord] = []
            first_query, first_subject = batch.job.specs[0]
            raw_text = run_losat_search(
                spec,
                batch.query.fasta,
                batch.subject.fasta,
                options=search_args(first_query, first_subject),
                losat_bin=runtime.losat_executable,
                ncbi_blast_bin=runtime.ncbi_blast_executable,
                threads=runtime.threads,
                runtime_callback=runtime_records.append,
            )
            split = split_losat_batch_result(raw_text, batch, batch.job.specs)
            for pair, text in split.items():
                query, subject = pair
                entry: dict[str, object] = {
                    "schema": NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
                    "kind": "raw-losat",
                    "identityKind": "nucleotide",
                    "key": keys[pair],
                    "text": text,
                    "program": spec.search,
                    "outfmt": LOSAT_OUTFMT,
                    "args": list(batch.job.args),
                    **({"searchContext": batch.search_context} if batch.search_context else {}),
                    "queryCanonicalHash": sha256_text(record_fasta(query)),
                    "subjectCanonicalHash": sha256_text(record_fasta(subject)),
                }
                searched[pair] = (entry, runtime_records[0] if runtime_records else None)
        for pair, text in split.items():
            texts[pair] = (keys[pair], text)

    # Session display order is edge order, as in the Web.
    for query, subject in pairs:
        if (query, subject) not in searched:
            continue
        entry, runtime_record = searched[(query, subject)]
        cache.store_search_entry(
            str(entry["key"]),
            entry,
            runtime=runtime_record,
            filename=losat_edge_filename(
                records[query].label, records[subject].label, search.program,
                records[query].window, records[subject].window,
            ),
        )

    comparisons = []
    for query, subject in pairs:
        _key, text = texts[(query, subject)]
        matches = read_comparison_table(
            io.StringIO(text), label=f"{search.program.upper()} result {query + 1}->{subject + 1}"
        )
        comparisons.append(
            LinearComparison(query, subject, matches, search_frame_text=text)
        )
    resolved = replace(
        options,
        losat_search=None,
        linear_comparisons=(*(options.linear_comparisons or ()), *comparisons),
    )
    return resolved, cache.session_entries()


__all__ = [
    "NUCLEOTIDE_LOSAT_PROGRAMS",
    "NucleotideLosatRecord",
    "losat_edge_filename",
    "linear_search_spaces",
    "resolve_linear_nucleotide_losat",
]
