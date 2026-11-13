#!/usr/bin/env python
# coding: utf-8

"""The LOSAT job plan: which searches run, against which database, and under
which raw-cache key.

This is the Python owner of the rules that the Web keeps in
``planLosatSourceJobs``, ``prepareLosatSourceBatches`` and
``splitLosatSourceResult`` (``gbdraw/web/js/services/linear-sources.js``), the
LOSATP record-pair searches of ``buildLosatJobSpecs``
(``gbdraw/web/js/services/linear-comparisons.js``), the nucleotide FASTA text of
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

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Hashable, Literal, Mapping, Sequence, TypedDict, cast

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


class _LosatJobDraft(TypedDict):
    query: tuple[int, ...]
    subject: tuple[int, ...]
    args: tuple[str, ...]
    scope: LosatJobScope
    specs: list[RecordPair]


def sha256_text(text: str) -> str:
    """Hex SHA-256 of UTF-8 text (the Web ``hashText``)."""

    return hashlib.sha256(str(text).encode("utf-8")).hexdigest()


def losat_entry_filename(entry: Mapping[str, object]) -> str:
    """The file name a raw LOSAT cache entry saved, or ``""`` when it saved no text.

    A Session is unvalidated data: a ``filename`` that is not text (``7``,
    ``["a.tsv"]``, ``{"name": ...}``) is read as none, as the Web reads it, so
    the caller's default name applies.
    """

    filename = entry.get("filename")
    return filename if isinstance(filename, str) else ""


def losat_pair_filename(entry: Mapping[str, object], ordinal: int) -> str:
    """The name of a shown raw LOSAT entry, ``ordinal`` counting entries from 1.

    Its saved text name, else ``losat_pair_<ordinal>.tsv``, the name the Web
    gives an entry without one (OV-1007).
    """

    return losat_entry_filename(entry) or f"losat_pair_{ordinal}.tsv"


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


def effective_collinear_search_scope(
    mode: LosatpSpecMode | str,
    *,
    infer_orthogroups: bool,
    search_scope: str,
) -> str:
    """The Search between scope a LOSATP run uses (Web ``effectiveCollinearSearchScope``).

    Search between widens only the similarity-group evidence: Collinear blocks
    without inferred groups read the displayed edges alone, so Collinear with
    inference off searches the edges whatever the stored scope holds.
    """

    if mode == "collinear" and not infer_orthogroups:
        return "adjacent"
    return str(search_scope)


def losatp_job_specs(
    mode: LosatpSpecMode | str,
    *,
    record_count: int,
    pairs: Sequence[RecordPair] | None = None,
    infer_orthogroups: bool,
    search_scope: str = "adjacent",
) -> tuple[RecordPair, ...]:
    """The directed record-pair searches of one LOSATP run (Web ``buildLosatJobSpecs``).

    ``pairs`` are the displayed comparison edges; omitted, they are the
    consecutive records. Pairwise searches each edge. Similarity groups search
    every record against every record, itself included. Collinear searches
    each edge in both directions and, with inference on, adds the
    within-record searches and searches every record pair when
    ``search_scope="all"`` (:func:`effective_collinear_search_scope`).
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
    scope = effective_collinear_search_scope(
        mode, infer_orthogroups=infer_orthogroups, search_scope=search_scope
    )
    if scope == "all":
        every_pair(False)
    else:
        for query, subject in edges:
            specs.extend(((query, subject), (subject, query)))
    return tuple(specs)


# One E-value scale per genome pair (OV-3143).


def linked_genomes(rows: Sequence[Hashable], source_ids: Sequence[Hashable]) -> list[int]:
    """The genome of each record, by record indexes: the source files linked through the rows they share.

    Two files with records in one row are one genome, and so is every file linked to them through other shared
    rows (``rows[i]`` and ``source_ids[i]`` are record ``i``'s row and source file); a file that shares no row
    with another file is a genome of its own, wherever its records stand. The one builder of genome keys:
    :func:`genome_search_spaces` reads them.
    """

    parent: dict[Hashable, Hashable] = {}

    def root(source: Hashable) -> Hashable:
        while parent.setdefault(source, source) != source:
            parent[source] = parent[parent[source]]
            source = parent[source]
        return source

    first_of_row: dict[Hashable, Hashable] = {}
    for row, source_id in zip(rows, source_ids, strict=True):
        parent[root(source_id)] = root(first_of_row.setdefault(row, source_id))
    roots: dict[Hashable, int] = {}
    return [roots.setdefault(root(source_id), len(roots)) for source_id in source_ids]


# The largest search space LOSAT reads (``-searchsp``): a signed 64-bit integer, as NCBI BLAST+'s.
LOSAT_SEARCH_SPACE_LIMIT = 2**63 - 1

# The warning of a genome pair searched per file (Owner decision 2026-11-10), with its parameters: the command line
# fills them with ``str.format``, the Web app with ``tm()``, which shows it translated (a catalog key, generated into
# gbdraw/web/js/contracts/card-windows.generated.js by tools/generate_web_contracts.py).
PER_FILE_WARNING_TEXT = (
    "The search space of genomes {first} and {second} ({space}) is past the largest LOSAT reads ({limit}), "
    "so their E-values are computed per file, not per genome."
)


@dataclass(frozen=True)
class PerFileWarning:
    """The warning of one genome pair searched per file (:data:`PER_FILE_WARNING_TEXT`): the two genomes' names,
    in sorted order, and the pair's search space and LOSAT's limit, each written with thousands separators."""

    first: str
    second: str
    space: str
    limit: str

    def params(self) -> dict[str, str]:
        return {"first": self.first, "second": self.second, "space": self.space, "limit": self.limit}

    def message(self) -> str:
        """The English text the command line logs and the Web app shows untranslated."""

        return PER_FILE_WARNING_TEXT.format(**self.params())


@dataclass(frozen=True)
class GenomeSearchSpaces:
    """The ``-searchsp`` of each record-pair search (:func:`genome_search_spaces`), called with record indexes.

    ``query_letters[i]`` and ``subject_letters[i]`` are record ``i``'s genome length in the program's letters on
    each side, ``genomes[i]`` its genome, ``whole[i]`` whether its genome is searched whole (:func:`_whole_genomes`).
    A pair of two whole genomes passes none (LOSAT's own); so does a genome pair whose search space is past
    :data:`LOSAT_SEARCH_SPACE_LIMIT`, whose E-values are then computed per file (Owner decision 2026-11-10;
    :meth:`per_file_warning` says so). A frozen value, so a render stage may take it as an argument.
    """

    query_letters: tuple[int, ...]
    subject_letters: tuple[int, ...]
    genomes: tuple[int, ...]
    whole: tuple[bool, ...]

    def _genome_space(self, query: int, subject: int) -> int | None:
        if self.whole[query] and self.whole[subject]:
            return None
        return self.query_letters[query] * self.subject_letters[subject]

    def __call__(self, query: int, subject: int) -> int | None:
        space = self._genome_space(query, subject)
        return None if space is None or space > LOSAT_SEARCH_SPACE_LIMIT else space

    def per_file_warning(self, query: int, subject: int, genome_names: Sequence[str]) -> PerFileWarning | None:
        """The warning of a pair searched per file because its genome pair's search space is past what LOSAT
        reads, or ``None``. ``genome_names[i]`` names record ``i``'s genome (:func:`genome_file_names`); the
        warning is one for the genome pair, whichever records and direction the search has."""

        space = self._genome_space(query, subject)
        if space is None or space <= LOSAT_SEARCH_SPACE_LIMIT:
            return None
        first, second = sorted((genome_names[query], genome_names[subject]))
        return PerFileWarning(first, second, f"{space:,}", f"{LOSAT_SEARCH_SPACE_LIMIT:,}")


def genome_search_spaces(
    *,
    rows: Sequence[Hashable],
    source_ids: Sequence[Hashable],
    source_lengths: Sequence[int],
    searched_lengths: Sequence[int],
    letters: tuple[int, int],
) -> GenomeSearchSpaces:
    """The search space each record-pair search passes LOSAT (``-searchsp``), by record indexes.

    A genome is the source files linked through the rows they share (:func:`linked_genomes`; ``rows[i]`` is
    record ``i``'s row); its length is the length of every record of those files, drawn or not
    (``source_lengths[i]`` is the length of record ``i``'s whole source, :func:`gbdraw.io.genome.source_genome_length`). A search
    between two genomes counts their lengths multiplied, each in the program's letters (``letters``: bases per
    letter on the query and subject sides), so every record, crop or window of the two genomes shares one
    E-value scale, whichever records a search reads. ``searched_lengths[i]`` is the length record ``i``'s search
    reads: a pair of two whole genomes (:func:`_whole_genomes`) needs no search space (``None``), since LOSAT's own
    is then the genome pair's.
    """

    lengths = _genome_totals(rows, source_ids, source_lengths, "lengths")
    query_letters, subject_letters = letters
    return GenomeSearchSpaces(
        query_letters=tuple(length // query_letters for length in lengths),
        subject_letters=tuple(length // subject_letters for length in lengths),
        genomes=tuple(linked_genomes(rows, source_ids)),
        whole=tuple(_whole_genomes(rows, source_ids, lengths, searched_lengths)),
    )


def source_file_name(source_id: Hashable) -> str:
    """The file name of a source (:func:`losat_source_ids`): its first path's name, as the input summary names it."""

    if isinstance(source_id, tuple) and source_id and source_id[0] != "memory":
        return Path(str(source_id[0])).name
    return "in-memory record"


def genome_file_names(
    rows: Sequence[Hashable], source_ids: Sequence[Hashable], name_of: Callable[[Hashable], str]
) -> list[str]:
    """The name of each record's genome: the names of its files (``name_of`` of each source), in record order."""

    genomes = linked_genomes(rows, source_ids)
    files: dict[int, list[str]] = {}
    for genome, source_id in zip(genomes, source_ids, strict=True):
        name = name_of(source_id)
        if name not in files.setdefault(genome, []):
            files[genome].append(name)
    return [" + ".join(files[genome]) for genome in genomes]


@dataclass(frozen=True)
class GenomeDatabaseSizes:
    """The ``-dbsize`` of each LOSATP record-pair search (:func:`genome_database_sizes`), called with record indexes.

    ``proteins[i]`` is record ``i``'s genome's protein letters, ``whole[i]`` whether record ``i``'s genome is
    searched whole (:func:`_whole_genomes`). A frozen value, so a render stage may take it as an argument.
    """

    proteins: tuple[int, ...]
    whole: tuple[bool, ...]

    def __call__(self, query: int, subject: int) -> int | None:
        return None if self.whole[query] and self.whole[subject] else self.proteins[subject]


def genome_database_sizes(
    *,
    rows: Sequence[Hashable],
    source_ids: Sequence[Hashable],
    source_lengths: Sequence[int],
    searched_lengths: Sequence[int],
    source_proteins: Sequence[int],
) -> GenomeDatabaseSizes:
    """The database size each LOSATP record-pair search passes LOSAT (``-dbsize``), by record indexes.

    The genomes are those of :func:`genome_search_spaces`; ``source_proteins[i]`` is the length of every protein
    of record ``i``'s whole source, extracted by the rules LOSATP searches with
    (``gbdraw.analysis.protein_colinearity.source_protein_letters``). A search counts its subject genome's
    proteins as its database, so every record, crop or window of the genome shares one E-value scale; the query
    side counts each protein as LOSAT does. A pair of two whole genomes needs no database size (``None``), as for
    the search space.
    """

    lengths = _genome_totals(rows, source_ids, source_lengths, "lengths")
    return GenomeDatabaseSizes(
        proteins=tuple(_genome_totals(rows, source_ids, source_proteins, "protein counts")),
        whole=tuple(_whole_genomes(rows, source_ids, lengths, searched_lengths)),
    )


def _whole_genomes(
    rows: Sequence[Hashable], source_ids: Sequence[Hashable], lengths: Sequence[int], searched_lengths: Sequence[int]
) -> list[bool]:
    """Whether each record's genome is searched whole: every drawn record of the genome reads all of it.

    ``lengths[i]`` is record ``i``'s genome length, ``searched_lengths[i]`` the length its search reads. One record
    that reads less (a crop, a window, one record of several) makes every search of its genome pass the genome's
    scale, so a source batch is never split by the scale (a record drawn whole and as a window is one batch).
    """

    genomes = linked_genomes(rows, source_ids)
    partial = {
        genome for genome, searched, length in zip(genomes, searched_lengths, lengths, strict=True) if searched != length
    }
    return [genome not in partial for genome in genomes]


def _genome_totals(
    rows: Sequence[Hashable], source_ids: Sequence[Hashable], values: Sequence[int], name: str
) -> list[int]:
    """Each record's genome total of a per-source value (``values[i]`` is record ``i``'s source's)."""

    by_source: dict[Hashable, int] = {}
    for source_id, value in zip(source_ids, values, strict=True):
        if by_source.setdefault(source_id, int(value)) != int(value):
            raise ValueError(f"Source {source_id!r} is given two {name}.")
    genomes = linked_genomes(rows, source_ids)
    totals: dict[int, int] = {}
    for source_id, genome in dict(zip(source_ids, genomes)).items():
        totals[genome] = totals.get(genome, 0) + by_source[source_id]
    return [totals[genome] for genome in genomes]


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
    drafts: dict[tuple[object, ...], _LosatJobDraft] = {}
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
            query_indexes=drafts[key]["query"],
            subject_indexes=drafts[key]["subject"],
            args=drafts[key]["args"],
            scope=drafts[key]["scope"],
            specs=tuple(drafts[key]["specs"]),
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

    @property
    def searchable(self) -> bool:
        """Whether both sides hold a sequence. A side of records without a CDS
        protein holds none: its pairs find no hit and LOSAT, which refuses an
        empty database, does not run (OV-3493; Web ``batch.searchable``)."""

        return bool(self.query.ids and self.subject.ids)


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

    rows: dict[RecordPair, list[str]] = {cast(RecordPair, tuple(spec)): [] for spec in specs}
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
    """Raw key of one search against nucleotide records (Web ``buildLosatCachePayload``).

    ``program`` is the search program (``blastn``, ``tblastx``, or ``tblastn``,
    whose query is protein FASTA); the hashes are of each side's FASTA: a
    record's search-frame FASTA, or a batch side's text.
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
    "LOSAT_SEARCH_SPACE_LIMIT",
    "PER_FILE_WARNING_TEXT",
    "PerFileWarning",
    "GenomeSearchSpaces",
    "genome_file_names",
    "source_file_name",
    "LosatBatch",
    "LosatBatchSide",
    "LosatJob",
    "LosatpSpecMode",
    "NUCLEOTIDE_LOSAT_CACHE_SCHEMA",
    "effective_collinear_search_scope",
    "GenomeDatabaseSizes",
    "genome_database_sizes",
    "genome_search_spaces",
    "linked_genomes",
    "losat_entry_filename",
    "losat_fasta_text",
    "losat_pair_filename",
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
