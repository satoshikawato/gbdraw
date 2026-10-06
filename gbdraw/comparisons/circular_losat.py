#!/usr/bin/env python
# coding: utf-8

"""Resolve Circular LOSATN / TLOSATX ring intent into ring rows and raw cache entries.

The Circular planners call :func:`resolve_circular_conservation_losat` once the
records are loaded. Each comparison genome (all records of one file, the
query) is searched against all displayed records (the subject database, design
3.3), so every ring has the reference genome as its E-value database. The
result replaces the search intent with :class:`ConservationSearchResult` rows
and Web-identical ``flow: 'circular-conservation'`` raw cache entries, so
Session replay needs no LOSAT. Rendering, thresholds and track slots are those
of precomputed rings.
"""

from __future__ import annotations

from dataclasses import replace
import logging
import os
import re
from typing import Any, Callable, Mapping, Sequence, cast

from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.conservation import ConservationSearchResult
from gbdraw.api.options import CircularDiagramOptions
from gbdraw.comparisons.losat_jobs import (
    LOSAT_OUTFMT,
    NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
    losat_fasta_text,
    losat_search_frame_fasta,
    nucleotide_losat_cache_key,
    sha256_text,
)
from gbdraw.comparisons.losat_runtime import (
    LOSAT_PROGRAMS,
    LosatRawCache,
    LosatSearchArgs,
    losat_cache_args,
    require_losat_task_support,
    run_losat_search,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparison_sequences import ComparisonSequenceFile

logger = logging.getLogger(__name__)

CIRCULAR_CONSERVATION_FLOW = "circular-conservation"
_OUTFMT6_COLUMNS = 12
# Effective translation table when none is set (design D17; never /transl_table).
_DEFAULT_GENCODE = 1


def ring_losat_filename(source_path: str, program: str) -> str:
    """Raw TSV name of one ring, as the Web names it."""

    stem = re.sub(r"\.[^.]+$", "", os.path.basename(str(source_path))) or "comparison"
    cleaned = re.sub(r"[^\w.-]+", "_", f"{stem}.circular_conservation.{program}.tsv")
    return cleaned.strip("_") or "gbdraw_session"


def _fasta_ids(fasta: str) -> set[str]:
    return {
        line[1:].split()[0]
        for line in fasta.splitlines()
        if line.startswith(">") and line[1:].split()
    }


def _validated_rows(text: str, *, query_ids: set[str], subject_ids: set[str]) -> str:
    for line_number, line in enumerate(re.split(r"\r?\n", str(text)), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        columns = line.split("\t")
        if (
            len(columns) != _OUTFMT6_COLUMNS
            or columns[0] not in query_ids
            or columns[1] not in subject_ids
        ):
            raise ValidationError(
                "LOSAT result contains malformed or unrecognized record endpoints "
                f"(output line {line_number}).",
                diagnostic={"code": "LOSAT_RUNTIME", "reason": "OUTPUT", "row": line_number},
            )
    return text


def comparison_query_fasta(comparison: ComparisonSequenceFile, *, ordinal: int) -> str:
    """LOSAT query FASTA of one comparison genome: every record, Web FASTA layout.

    The CLI ring search and the Web worker helper
    (``gbdraw.web_support.comparison_sequences``) hash this text, so one
    sequence has one raw key whatever its file format or FASTA layout.
    """

    if not comparison.records:
        raise ValidationError(
            f"Comparison sequence file #{ordinal} "
            f"({os.path.basename(comparison.path)}) has no sequence record.",
            diagnostic={
                "code": "INPUT_UNREADABLE",
                "reason": "SEQUENCE_MISSING",
                "field": "comparison_sequence",
            },
        )
    return "".join(
        losat_fasta_text(record.id, str(record.seq)) for record in comparison.records
    )


def _ring_gencodes(options: CircularDiagramOptions, count: int) -> list[int]:
    gencodes = tuple(options.conservation_losat_gencodes or ())
    if not gencodes:
        return [_DEFAULT_GENCODE] * count
    if len(gencodes) == 1:
        return [gencodes[0]] * count
    return list(gencodes)


def resolve_circular_conservation_losat(
    options: CircularDiagramOptions,
    *,
    records: Sequence[SeqRecord],
    load_sequences: Callable[[], Sequence[ComparisonSequenceFile]],
) -> tuple[CircularDiagramOptions, tuple[Mapping[str, Any], ...]]:
    """Run the ring searches and return resolved options and raw cache entries.

    ``load_sequences`` returns the comparison genomes in ring order (the
    request's memoized reader result).
    """

    search = options.losat_search
    if search is None:
        return options, ()
    spec = LOSAT_PROGRAMS[search.program]
    runtime = search.runtime
    comparison_files = tuple(load_sequences())
    paths = tuple(options.conservation_sequence_files or ())
    query_fastas = [
        comparison_query_fasta(comparison, ordinal=index)
        for index, comparison in enumerate(comparison_files, start=1)
    ]
    if search.program == "losatn" and search.losatn_task != "megablast":
        require_losat_task_support(
            str(search.losatn_task),
            losat_bin=runtime.losat_executable,
            ncbi_blast_bin=runtime.ncbi_blast_executable,
        )

    # The subject database: every displayed record, in its search frame.
    subject_fasta = "".join(losat_search_frame_fasta(record) for record in records)
    subject_hash = sha256_text(subject_fasta)
    subject_ids = _fasta_ids(subject_fasta)
    reference_gencodes = tuple(search.record_gencodes)
    reference_gencode = (
        reference_gencodes[0]
        if reference_gencodes and reference_gencodes[0] is not None
        else _DEFAULT_GENCODE
    )
    ring_gencodes = _ring_gencodes(options, len(comparison_files))

    def search_args(index: int) -> LosatSearchArgs:
        if search.program == "losatn":
            return LosatSearchArgs(task=search.losatn_task)
        # The Web ring flow always passes both tables (default 1).
        return LosatSearchArgs(query_gencode=ring_gencodes[index], db_gencode=reference_gencode)

    cache = LosatRawCache()
    results: list[ConservationSearchResult] = []
    for index, comparison in enumerate(comparison_files):
        query_fasta = query_fastas[index]
        query_hash = sha256_text(query_fasta)
        args = search_args(index)
        cache_args = losat_cache_args(spec, args)
        key = nucleotide_losat_cache_key(
            program=spec.search,
            args=cache_args,
            query_hash=query_hash,
            subject_hash=subject_hash,
            flow=CIRCULAR_CONSERVATION_FLOW,
        )
        filename = ring_losat_filename(
            # Ring LOSAT validation rejects an empty comparison sequence entry.
            cast(str, paths[index]) if index < len(paths) else comparison.path,
            search.program,
        )
        cached = cache.cached_entry(key)
        if cached is not None:
            text = str(cached["text"])
        else:
            logger.info(
                "INFO: %s ring %d/%d (%s against %d displayed record(s)).",
                search.program.upper(), index + 1, len(comparison_files),
                os.path.basename(comparison.path), len(records),
            )
            runtime_records: list[dict[str, object]] = []
            text = _validated_rows(
                run_losat_search(
                    spec,
                    query_fasta,
                    subject_fasta,
                    options=args,
                    losat_bin=runtime.losat_executable,
                    ncbi_blast_bin=runtime.ncbi_blast_executable,
                    threads=runtime.threads,
                    runtime_callback=runtime_records.append,
                ),
                query_ids=_fasta_ids(query_fasta),
                subject_ids=subject_ids,
            )
            entry: dict[str, object] = {
                "schema": NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
                "kind": "raw-losat",
                "identityKind": "nucleotide",
                "key": key,
                "text": text,
                "program": spec.search,
                "flow": CIRCULAR_CONSERVATION_FLOW,
                "outfmt": LOSAT_OUTFMT,
                "args": list(cache_args),
                "queryCanonicalHash": query_hash,
                "subjectCanonicalHash": subject_hash,
            }
            cache.store_search_entry(
                key,
                entry,
                runtime=runtime_records[0] if runtime_records else None,
                filename=filename,
            )
        results.append(ConservationSearchResult(name=filename, text=text))

    labels = (
        tuple(options.conservation_labels)
        if options.conservation_labels is not None
        else tuple(comparison.label for comparison in comparison_files)
    )
    resolved = replace(
        options,
        losat_search=None,
        conservation_losat_gencodes=None,
        conservation_search_results=tuple(results),
        conservation_reference="subject",
        conservation_labels=labels,
    )
    return resolved, cache.session_entries()


__all__ = [
    "CIRCULAR_CONSERVATION_FLOW",
    "comparison_query_fasta",
    "resolve_circular_conservation_losat",
    "ring_losat_filename",
]
