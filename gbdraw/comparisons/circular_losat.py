#!/usr/bin/env python
# coding: utf-8

"""Resolve Circular LOSATN / TLOSATX ring intent into ring rows and raw cache entries.

The Circular planners call :func:`resolve_circular_conservation_losat` once the
records are loaded. Each comparison genome (all records of its files, each
file whole in the listed order: the query) is searched against the subject database: every record of the
subject-side source files, whole (:func:`ring_subject_sources`, Owner D-25), as
the Web searches it. Records that are not drawn, a single presentation and a
crop do not change it, and the raw cache keys equal the Web's. A ring search
passes the search space of its genome pair, the reference genome and the
ring's (:func:`ring_search_spaces`, PD-OI-018 revision 6), so its E-values do
not depend on how either genome is split into records or files. The raw cache entries keep
the search output; :func:`ring_rows_in_search_frame` turns it into the ring
rows of the drawn records, in their search frame like every ring table
(PD-OI-073, Owner D-28), which replace the search intent as
:class:`ConservationSearchResult` rows. Session replay needs no LOSAT.
Rendering, thresholds and track slots are those of precomputed rings; the Web
converts its search output through the same owner
(``gbdraw.web_support.ring_search_frame``).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import logging
import os
import re
from typing import Any, Callable, Hashable, Iterable, Mapping, Sequence

from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.conservation import ConservationSearchResult
from gbdraw.api.options import CircularDiagramOptions, LosatSearchOptions
from gbdraw.api.requests import FastaInputSource, GenBankInputSource, GffFastaInputSource, RecordInputSource
from gbdraw.comparisons.losat_jobs import (
    LOSAT_OUTFMT,
    NUCLEOTIDE_LOSAT_CACHE_SCHEMA,
    PerFileWarning,
    genome_file_names,
    genome_search_spaces,
    losat_fasta_text,
    nucleotide_losat_cache_key,
    sha256_text,
)
from gbdraw.comparisons.losat_runtime import (
    LOSAT_PROGRAMS,
    LosatRawCache,
    LosatRuntimeRecord,
    LosatSearchArgs,
    losat_cache_args,
    losat_takes_genome_scale,
    require_losat_task_support,
    run_losat_search,
)
from gbdraw.core.record_metadata import _read_coord_map, _record_source_start
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparison_sequences import ComparisonGenome, comparison_file_stem
from gbdraw.io.genome import read_sequence_file_records, record_sequence_text, source_genome_length

logger = logging.getLogger(__name__)

CIRCULAR_CONSERVATION_FLOW = "circular-conservation"
_OUTFMT6_COLUMNS = 12
# Effective translation table when none is set (design D17; never /transl_table).
_DEFAULT_GENCODE = 1


def ring_losat_filename(source_path: str, program: str) -> str:
    """Raw TSV name of one ring, as the Web names it."""

    stem = comparison_file_stem(str(source_path)) or "comparison"
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


def comparison_query_fasta(comparison: ComparisonGenome, *, ordinal: int) -> str:
    """LOSAT query FASTA of one comparison genome: every record of every file, Web FASTA layout.

    The CLI ring search and the Web worker helper
    (``gbdraw.web_support.comparison_sequences``) hash this text, so one
    sequence has one raw key whatever its file format or FASTA layout.
    """

    if not comparison.records:
        raise ValidationError(
            f"Comparison genome #{ordinal} "
            f"({', '.join(os.path.basename(file.path) for file in comparison.files)}) has no sequence record.",
            diagnostic={
                "code": "INPUT_UNREADABLE",
                "reason": "SEQUENCE_MISSING",
                "field": "comparison_sequence",
            },
        )
    return "".join(
        losat_fasta_text(record.id, str(record.seq)) for record in comparison.records
    )


@dataclass(frozen=True)
class RingGenomeFile:
    """One file of a ring search's genome, the reference or a ring: the name a warning gives it, its length
    (every record, uncropped: :func:`gbdraw.io.genome.source_genome_length`) and the length of each record the
    search reads of it, every record whole. A frozen value."""

    name: str
    length: int
    record_lengths: tuple[int, ...]


def ring_genome_file(name: str, kind: str, paths: Sequence[str], records: Sequence[SeqRecord]) -> RingGenomeFile:
    """The genome file of ``records``, every record of one source (``kind`` and ``paths`` as
    :func:`gbdraw.io.genome.source_genome_length` reads them), as the CLI and the Web read it."""

    lengths = tuple(len(record) for record in records)
    return RingGenomeFile(name, source_genome_length(kind, paths, sum(lengths)), lengths)


@dataclass(frozen=True)
class RingSubjectSource:
    """One source of the ring subject database: its records, whole, in file order, and its genome file."""

    records: tuple[SeqRecord, ...]
    genome: RingGenomeFile


def ring_subject_sources(
    sources: Iterable[RecordInputSource],
    *,
    file_records: Callable[[GenBankInputSource | FastaInputSource], Sequence[SeqRecord]],
) -> tuple[RingSubjectSource, ...]:
    """The sources of the ring subject database, in request order (Owner D-25): the reference genome.

    Every distinct subject-side source the request reads gives all of its
    records, uncropped, in file order, on the file's strand, also records that
    are not drawn. A source whose file holds whole records (a GenBank file, a
    FASTA given without GFF3) gives the records the render's loader read
    (``file_records``: the planner's parse, so no source is parsed twice). A
    GFF3/FASTA pair gives every entry of its FASTA file, also one the GFF3 file
    does not name: the GFF3 parse keeps only named entries, so the FASTA file is
    read again, and that read is its genome length too (the FASTA file is the
    source here, as its key says). An in-memory record (Python API) is its own
    source, as given.
    """

    seen: set[Hashable] = set()
    found: list[RingSubjectSource] = []
    for source in sources:
        key: Hashable = (
            ("genbank", str(source.path)) if isinstance(source, GenBankInputSource)
            else ("fasta", str(source.fasta_path)) if isinstance(source, (GffFastaInputSource, FastaInputSource))
            else ("memory", id(source.record))
        )
        if key in seen:
            continue
        seen.add(key)
        if isinstance(source, GenBankInputSource):
            records = tuple(file_records(source))
            genome = ring_genome_file(os.path.basename(str(source.path)), "genbank", (str(source.path),), records)
        elif isinstance(source, FastaInputSource):
            records = tuple(file_records(source))
            genome = ring_genome_file(
                os.path.basename(str(source.fasta_path)), "fasta", (str(source.fasta_path),), records
            )
        elif isinstance(source, GffFastaInputSource):
            records = tuple(read_sequence_file_records(str(source.fasta_path), "fasta"))
            genome = ring_genome_file(
                os.path.basename(str(source.gff_path)), "fasta", (str(source.fasta_path),), records
            )
        else:
            records = (source.record,)
            genome = ring_genome_file("in-memory record", "memory", (), records)
        found.append(RingSubjectSource(records, genome))
    return tuple(found)


def comparison_genome_files(
    comparison: ComparisonGenome, names: Sequence[str] | None = None
) -> tuple[RingGenomeFile, ...]:
    """The genome files of one ring's comparison genome, in the listed order; ``names`` names them (the Web stages
    files under fixed names), else each file's base name."""

    return tuple(
        ring_genome_file(
            names[index] if names is not None else os.path.basename(file.path), file.format, (file.path,), file.records
        )
        for index, file in enumerate(comparison.files)
    )


@dataclass(frozen=True)
class RingSearchSpaces:
    """The ``-searchsp`` of each ring search, in ring order (``None``: LOSAT's own), and the warning of a ring whose
    genome pair is searched per file (Owner decision 2026-11-10), else ``None``. A frozen value."""

    spaces: tuple[int | None, ...]
    warnings: tuple[PerFileWarning | None, ...]


def ring_search_spaces(
    program: str, reference: Sequence[RingGenomeFile], rings: Sequence[Sequence[RingGenomeFile]]
) -> RingSearchSpaces:
    """The search space each ring search of ``program`` (``losatn`` or ``tlosatx``) passes LOSAT (PD-OI-018 revision 6).

    The reference genome is every file of ``reference`` (the subject database), ring ``k``'s genome every file of
    ``rings[k]``, each file counting every record (:class:`RingGenomeFile`). A ring search passes the two genome
    lengths multiplied, in the program's letters (:func:`gbdraw.comparisons.losat_jobs.genome_search_spaces`), so a
    hit has one E-value however either genome is split into records or files. A pair of two genomes that are each
    one file of one record passes none: LOSAT's own search space is then the genome pair's, and the figure, raw
    cache key and Session stay as they were. The one owner of the ring scales: the CLI ring search and the Web's
    (``gbdraw.web_support.comparison_sequences.plan_ring_scales_json``) read it.
    """

    rows: list[Hashable] = []
    source_ids: list[Hashable] = []
    source_lengths: list[int] = []
    searched_lengths: list[int] = []
    names: dict[Hashable, str] = {}
    firsts: list[int] = []
    genomes: list[tuple[Hashable, Sequence[RingGenomeFile]]] = [("reference", reference)]
    genomes.extend((("ring", index), ring) for index, ring in enumerate(rings))
    for genome, files in genomes:
        firsts.append(len(rows))
        for file_index, file in enumerate(files):
            source = (genome, file_index)
            names[source] = file.name
            for record_length in file.record_lengths:
                rows.append(genome)
                source_ids.append(source)
                source_lengths.append(file.length)
                searched_lengths.append(record_length)
        if firsts[-1] == len(rows):
            raise ValueError(f"The {genome} genome of the ring searches has no record.")
    spaces = genome_search_spaces(
        rows=rows, source_ids=source_ids, source_lengths=source_lengths, searched_lengths=searched_lengths,
        letters=LOSAT_PROGRAMS[program].space_letters,
    )
    genome_names = genome_file_names(rows, source_ids, names.__getitem__)
    reference_index, ring_indexes = firsts[0], firsts[1:]
    return RingSearchSpaces(
        spaces=tuple(spaces(ring, reference_index) for ring in ring_indexes),
        warnings=tuple(spaces.per_file_warning(ring, reference_index, genome_names) for ring in ring_indexes),
    )


def ring_subject_fasta(records: Iterable[SeqRecord]) -> str:
    """LOSAT subject FASTA of the rings: each record that has a sequence, Web FASTA layout.

    A record without a sequence (a CON record) is left out. The Web's
    ``readAllLosatFasta`` writes the same text for the same file, so the
    subject hash and the raw cache keys are the same on both surfaces.
    """

    fasta = "".join(
        losat_fasta_text(str(record.id), sequence)
        for record in records
        if (sequence := record_sequence_text(record)) is not None
    )
    if not fasta:
        raise ValidationError(
            "The reference of the similarity rings has no record with a sequence.",
            diagnostic={
                "code": "INPUT_UNREADABLE",
                "reason": "SEQUENCE_MISSING",
                "field": "records",
            },
        )
    return fasta


def _fasta_id(record: SeqRecord) -> str:
    """The ID a search reports for a record: the first word of its FASTA header."""

    words = str(record.id).split()
    return words[0] if words else ""


def _one_frame_per_id(frames: Iterable[tuple[str, tuple[int, int]]]) -> dict[str, tuple[int, int]]:
    by_id: dict[str, tuple[int, int]] = {}
    for record_id, frame in frames:
        if by_id.setdefault(record_id, frame) != frame:
            raise ValidationError(
                f"Similarity rings cannot place their hits on record {record_id}: the request has two records "
                "with this ID that are cropped or oriented differently. Draw them in separate requests, or "
                "give the records distinct IDs.",
                diagnostic={"code": "COMPARISON_INPUT", "reason": "UNIQUE_IDS", "field": "records"},
            )
    return by_id


def ring_rows_in_search_frame(
    text: str,
    *,
    drawn_records: Sequence[SeqRecord],
    subject_records: Sequence[SeqRecord] = (),
) -> str:
    """Ring search output as the ring rows of the drawn records, in their search frame.

    Every ring table is in the search frame of its record (PD-OI-073, Owner
    D-28): the drawn record after its crop, on its source strand, 1-based. The
    search reports subject positions on the database records (the whole
    reference files, :func:`ring_subject_records`), so each row is moved to
    the search frame of the drawn record with its subject ID. Rows of records
    that are not drawn, and rows wholly outside a cropped record, are left
    out. A row that crosses the crop edge keeps its HSP with the subject span
    clipped to the record (Owner-delegated); the query columns and the other
    columns stay as the search reported them.

    A database record read from a file is the whole source record. One passed
    in memory (``subject_records``, Python API) is searched as given: its
    coordinate map turns its positions into source positions, so an already
    cropped or reverse complemented record converts correctly. Records not
    named in ``subject_records`` are whole source records (the Web).
    """

    frames = _one_frame_per_id(
        (_fasta_id(record), (_record_source_start(record), len(record.seq))) for record in drawn_records
    )
    database_maps = _one_frame_per_id((_fasta_id(record), _read_coord_map(record)) for record in subject_records)
    rows: list[str] = []
    for line_number, line in enumerate(re.split(r"\r?\n", str(text)), start=1):
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        columns = line.split("\t")
        try:
            positions = (int(columns[8]), int(columns[9])) if len(columns) == _OUTFMT6_COLUMNS else None
        except ValueError:
            positions = None
        if positions is None:
            raise ValidationError(
                f"LOSAT result is not BLAST outfmt 6 rows (output line {line_number}).",
                diagnostic={"code": "LOSAT_RUNTIME", "reason": "OUTPUT", "row": line_number},
            )
        frame = frames.get(columns[1])
        if frame is None:
            continue
        frame_start, length = frame
        base, step = database_maps.get(columns[1], (1, 1))
        ends = [base + step * (position - 1) - frame_start + 1 for position in positions]
        if max(ends) < 1 or min(ends) > length:
            continue
        columns[8], columns[9] = (str(min(max(end, 1), length)) for end in ends)
        rows.append("\t".join(columns) + "\n")
    return "".join(rows)


def _ring_gencodes(options: CircularDiagramOptions, count: int) -> list[int]:
    gencodes = tuple(options.conservation_losat_gencodes or ())
    if not gencodes:
        return [_DEFAULT_GENCODE] * count
    if len(gencodes) == 1:
        return [gencodes[0]] * count
    return list(gencodes)


def _ring_search_spaces(
    search: LosatSearchOptions,
    subject_sources: Sequence[RingSubjectSource],
    comparison_files: Sequence[ComparisonGenome],
) -> tuple[int | None, ...]:
    """Each ring search's ``-searchsp`` (:func:`ring_search_spaces`), with the warning of each ring searched per file
    logged once. A runtime without ``-searchsp`` (LOSAT v0.1.0) searches as before, with a warning."""

    scales = ring_search_spaces(
        search.program,
        [source.genome for source in subject_sources],
        [comparison_genome_files(comparison) for comparison in comparison_files],
    )
    for warning in dict.fromkeys(scales.warnings):
        if warning is not None:
            logger.warning("%s", warning.message())
    if all(space is None for space in scales.spaces) or losat_takes_genome_scale(
        LOSAT_PROGRAMS[search.program],
        losat_bin=search.runtime.losat_executable,
        ncbi_blast_bin=search.runtime.ncbi_blast_executable,
    ):
        return scales.spaces
    return tuple(None for _space in scales.spaces)


def resolve_circular_conservation_losat(
    options: CircularDiagramOptions,
    *,
    subject_sources: Sequence[RingSubjectSource],
    drawn_records: Sequence[SeqRecord],
    load_sequences: Callable[[], Sequence[ComparisonGenome]],
) -> tuple[CircularDiagramOptions, tuple[Mapping[str, Any], ...]]:
    """Run the ring searches and return resolved options and raw cache entries.

    ``subject_sources`` are the sources of the subject database, the reference
    genome (:func:`ring_subject_sources`); ``drawn_records`` the records the request
    draws, whose search frame the ring rows are in
    (:func:`ring_rows_in_search_frame`); ``load_sequences`` returns the
    comparison genomes in ring order (the request's memoized reader result).
    """

    search = options.losat_search
    if search is None:
        return options, ()
    spec = LOSAT_PROGRAMS[search.program]
    runtime = search.runtime
    comparison_files = tuple(load_sequences())
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

    # The subject database: the whole subject-side source files (D-25).
    subject_records = tuple(record for source in subject_sources for record in source.records)
    subject_fasta = ring_subject_fasta(subject_records)
    subject_hash = sha256_text(subject_fasta)
    subject_ids = _fasta_ids(subject_fasta)
    reference_gencodes = tuple(search.record_gencodes)
    reference_gencode = (
        reference_gencodes[0]
        if reference_gencodes and reference_gencodes[0] is not None
        else _DEFAULT_GENCODE
    )
    ring_gencodes = _ring_gencodes(options, len(comparison_files))
    searchsp = _ring_search_spaces(search, subject_sources, comparison_files)

    def search_args(index: int) -> LosatSearchArgs:
        if search.program == "losatn":
            return LosatSearchArgs(task=search.losatn_task, searchsp=searchsp[index])
        # The Web ring flow always passes both tables (default 1).
        return LosatSearchArgs(
            query_gencode=ring_gencodes[index], db_gencode=reference_gencode, searchsp=searchsp[index]
        )

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
        # The raw TSV is named after the genome's first file.
        filename = ring_losat_filename(comparison.path, search.program)
        cached = cache.cached_entry(key)
        if cached is not None:
            text = str(cached["text"])
        else:
            logger.info(
                "%s ring %d/%d (%s against %d record(s) of the reference file(s)).",
                search.program.upper(), index + 1, len(comparison_files),
                ", ".join(os.path.basename(file.path) for file in comparison.files), len(subject_ids),
            )
            runtime_records: list[LosatRuntimeRecord] = []
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
        results.append(ConservationSearchResult(
            name=filename,
            text=ring_rows_in_search_frame(text, drawn_records=drawn_records, subject_records=subject_records),
        ))

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
    "ring_rows_in_search_frame",
    "comparison_genome_files",
    "ring_genome_file",
    "ring_search_spaces",
    "ring_subject_fasta",
    "ring_subject_sources",
    "RingGenomeFile",
    "RingSearchSpaces",
    "RingSubjectSource",
]
