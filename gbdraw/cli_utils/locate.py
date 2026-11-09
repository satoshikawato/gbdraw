"""The ``gbdraw locate`` command: where a reference span lies in other genomes.

It searches the genes of the reference windows in every record of each target
File (tblastn by default, or blastn), plans the windows of each File with the
locate owner (:mod:`gbdraw.comparisons.locate`), and writes the plan as a
records table and a comparisons table that ``gbdraw linear`` draws, with the
full plan as JSON and a report per File on standard error.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
from dataclasses import replace
from pathlib import Path
from typing import Mapping, Sequence

from gbdraw.api.input_summary import _threshold_text
from gbdraw.api.record_planning import CLI_RECORD_KEY_PREFIX, DrawnUnitName, comparison_selector_units
from gbdraw.cli_utils.common import setup_logging
from gbdraw.cli_utils.losat_output import parse_positive_int, write_losat_output_files
from gbdraw.comparisons.locate import (
    LOCATE_PROGRAMS,
    FilePlan,
    LocateDrawing,
    LocateInputFile,
    LocateOptions,
    LocatePlan,
    LocateProgramName,
    LocateSearchResult,
    LocateWarning,
    LocateWindowPlace,
    PlannedWindow,
    finish_locate,
    locate_default_thresholds,
    locate_drawing,
    locate_file_counts,
    locate_results,
    prepare_locate_search,
)
from gbdraw.comparisons.losat_runtime import (
    LosatRuntimeRecord,
    require_losat_program,
    run_losat_search,
)
from gbdraw.core.numeric import argparse_float, argparse_integer
from gbdraw.exceptions import ValidationError
from gbdraw.io.cli_tables import (
    ComparisonsTableEntry,
    RecordsTableEntry,
    comparisons_table_text,
    records_table_text,
)
from gbdraw.io.regions import RecordWindow, parse_region_spec, region_text
from gbdraw.mode_profiles import ComparisonThresholds
from gbdraw.render.output_paths import commit_staged_outputs, preflight_output_paths, staged_output_files

logger = logging.getLogger(__name__)

DEFAULT_PREFIX = "locate"
OUTPUT_SUFFIXES = (".records.tsv", ".comparisons.tsv", ".plan.json")


def _parser() -> argparse.ArgumentParser:
    defaults = locate_default_thresholds(LOCATE_PROGRAMS[0])
    parser = argparse.ArgumentParser(
        prog="gbdraw locate",
        description=(
            "Find where the genes of a reference span lie in other genomes and write the windows to draw: "
            "PREFIX.records.tsv and PREFIX.comparisons.tsv for gbdraw linear, and PREFIX.plan.json."
        ),
    )
    reference = parser.add_argument_group("reference")
    reference.add_argument("--reference_gbk", metavar="GBK", help="Reference GenBank file.")
    reference.add_argument("--reference_gff", metavar="GFF", help="Reference GFF3 file (requires --reference_fasta).")
    reference.add_argument("--reference_fasta", metavar="FASTA", help="Reference FASTA file.")
    reference.add_argument(
        "--reference_record", metavar="ID",
        help="The reference record (record ID or #index); required when the reference file holds several records.",
    )
    reference.add_argument(
        "--reference_region", metavar="START-END[:origin][:rc]", action="append", default=[],
        help=(
            "One reference window, 1-based and inclusive; repeat it once per window, in column order. "
            ":origin marks a window across the origin of a circular record, :rc draws it reversed "
            "(default: the whole record)."
        ),
    )
    targets = parser.add_argument_group("targets")
    targets.add_argument("--targets", metavar="GBK", nargs="+", default=[], help="Target GenBank files, one genome each.")
    targets.add_argument("--target_gffs", metavar="GFF", nargs="+", default=[], help="Target GFF3 files (with --target_fastas).")
    targets.add_argument("--target_fastas", metavar="FASTA", nargs="+", default=[], help="Target FASTA files.")
    parser.add_argument("-o", "--output", metavar="PREFIX", default=DEFAULT_PREFIX, help=f"Output prefix (default: {DEFAULT_PREFIX}).")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing output files (default: refuse).")
    search = parser.add_argument_group("search")
    search.add_argument(
        "--losat", choices=LOCATE_PROGRAMS, default=LOCATE_PROGRAMS[0],
        help="tlosatn searches each gene's protein (tblastn); losatn its nucleotide sequence (blastn) (default: tlosatn).",
    )
    search.add_argument("--losat_bin", metavar="PATH", default=None, help="Native LOSAT executable (default: automatic runtime resolution).")
    search.add_argument(
        "--ncbi_blast_bin", metavar="PATH", default=None,
        help="NCBI BLAST+ executable of the --losat program (tblastn or blastn) (default: automatic runtime resolution).",
    )
    search.add_argument("--losat_threads", metavar="N", type=parse_positive_int, default=None, help="Threads passed to each search (default: runtime default).")
    search.add_argument("--evalue", type=argparse_float, default=defaults.evalue, help=f"Maximum e-value (default: {defaults.evalue:g}).")
    search.add_argument("--bitscore", type=argparse_float, default=defaults.bitscore, help=f"Minimum bitscore (default: {defaults.bitscore:g}).")
    search.add_argument("--identity", type=argparse_float, default=defaults.identity, help=f"Minimum identity percentage (default: {defaults.identity:g}).")
    search.add_argument(
        "--alignment_length", type=argparse_integer, default=defaults.alignment_length,
        help=f"Minimum alignment length (default: {defaults.alignment_length}).",
    )
    search.add_argument("--losat_output_dir", metavar="DIR", default=None, help="Write the raw rows of each target File to DIR.")
    plan = parser.add_argument_group("windows")
    options = LocateOptions.for_program()
    plan.add_argument("--flank", metavar="BP", type=argparse_integer, default=options.flank, help=f"Bases kept on each side of the genes (default: {options.flank}).")
    plan.add_argument(
        "--merge_distance", metavar="BP", type=argparse_integer, default=options.merge_distance,
        help=f"Hits of one record at most this many bases apart are one window (default: {options.merge_distance}).",
    )
    plan.add_argument("--max_windows", metavar="N", type=parse_positive_int, default=options.max_windows, help=f"Most windows selected per File (default: {options.max_windows}).")
    return parser


def _input_kind(parser: argparse.ArgumentParser, gbk: Sequence[str], gff: Sequence[str], fasta: Sequence[str], name: str) -> None:
    if gbk and (gff or fasta):
        parser.error(f"Give the {name} as GenBank, or as GFF3 and FASTA, not both.")
    if gff and len(gff) != len(fasta):
        parser.error(f"Give one FASTA file per GFF3 file of the {name}.")
    if not (gbk or fasta):
        parser.error(f"Give the {name}.")


def _reference_windows(regions: Sequence[str]) -> list[RecordWindow]:
    windows = []
    for column, text in enumerate(regions, start=1):
        try:
            spec = parse_region_spec(text)
        except ValueError as exc:
            raise ValidationError(f"--reference_region {text!r}: {exc}", diagnostic={"code": "REGION_INVALID"}) from exc
        if spec.record_id is not None or spec.record_index is not None or spec.file_selector is not None:
            raise ValidationError(
                f"--reference_region {text!r} names a record; choose the record with --reference_record.",
                diagnostic={"code": "REGION_INVALID"},
            )
        windows.append(RecordWindow.from_region(spec, column=column))
    return windows


def _same_file(first: str, second: str) -> bool:
    try:
        return os.path.samefile(first, second)
    except OSError:
        return os.path.abspath(first) == os.path.abspath(second)


def _table_path(path: str, table_dir: Path) -> str:
    """A records-table cell for ``path``: relative to the table, as the reader resolves it."""

    absolute = os.path.abspath(path)
    try:
        return os.path.relpath(absolute, table_dir)
    except ValueError:  # Another drive on Windows.
        return absolute


def _records_table(drawing: LocateDrawing, files: Mapping[str, LocateInputFile], table_dir: Path) -> str:
    """The records table of ``drawing``: one table row per window, row 1 the reference, each File's row
    labelled with its file name on its first window."""

    entries = []
    for row_number, row in enumerate(drawing.rows, start=1):
        file = files[row.file_key]
        cells = {
            name: _table_path(path, table_dir) if path else ""
            for name, path in (("gbk", file.gbk), ("gff", file.gff), ("fasta", file.fasta))
        }
        for record_index, record in enumerate(row.records):
            for window in record.windows:
                entries.append(RecordsTableEntry(
                    **cells,
                    record_label=Path(row.file_key).name if record_index == 0 and window.column == 1 else "",
                    record_id=record.record_id,
                    region="" if record.whole else region_text(window.start, window.end, False),
                    reverse_complement=window.reverse_complement,
                    row=row_number,
                    column=window.column,
                    record_group="" if record.whole else f"{row_number}:{record.record_id}",
                ))
    return records_table_text(entries)


def _comparisons_table(drawing: LocateDrawing) -> str:
    """The window pairs of adjacent rows that share a reference gene. A record is named by its record ID, and a
    window of a record drawn in several windows by the record and its span as the records table writes it
    (``NC_009725.2:1911335-1989004``); ``#n`` takes the place of the ID wherever the comparisons table reader
    would not read the ID back as that record (a repeated ID, an ID that reads as a region or a file selector)."""

    records = [record for row in drawing.rows for record in row.records]
    units = [
        DrawnUnitName(record.record_id, record.record_key, None if record.whole else (window.start, window.end))
        for record in records
        for window in record.windows
    ]
    first_unit: dict[str, int] = {}
    for index, unit in enumerate(units):
        first_unit.setdefault(unit.record_key, index)
    numbers = {key: number for number, key in enumerate(first_unit, start=1)}

    def selector(place: LocateWindowPlace) -> str:
        row, index, window_index = place
        record = drawing.rows[row].records[index]
        window = record.windows[window_index]
        wanted = (first_unit[record.record_key] + window_index,)
        span = "" if len(record.windows) == 1 else ":" + region_text(window.start, window.end, False)
        for name in (record.record_id + span, f"#{numbers[record.record_key]}{span}"):
            if comparison_selector_units(units, name) == wanted:
                return name
        raise ValidationError(
            f"Window {span[1:] or 'of record'} {record.record_id} cannot be named in a comparisons table.",
            diagnostic={"code": "INPUT_INVALID"},
        )

    return comparisons_table_text([ComparisonsTableEntry(selector(upper), selector(lower)) for upper, lower in drawing.pairs])


def locate_tables(drawing: LocateDrawing, files: Mapping[str, LocateInputFile], table_dir: Path) -> tuple[str, str]:
    """The records table and comparisons table of ``drawing``, with paths relative to ``table_dir``."""

    return _records_table(drawing, files, table_dir), _comparisons_table(drawing)


def _span(window: PlannedWindow) -> str:
    text = f"{window.window.start:,}-{window.window.end:,}"
    return text + (" (−)" if window.window.reverse_complement else "")


def _warning_text(file_key: str, warning: LocateWarning, file_plan: FilePlan) -> str:
    args = dict(warning.args)
    if warning.code == "FILE_NO_WINDOWS":
        tally = file_plan.tally
        text = f"{file_key}: no window; the thresholds kept {tally.rows_kept:,} of {tally.rows_read:,} rows"
        return text + "".join(f"; {_threshold_text(item)}" for item in tally.thresholds) + "."
    return {
        "FILE_SEARCH_FAILED": lambda: f"{file_key}: the search failed: {args.get('error')}",
        "GENES_ABSENT": lambda: f"{file_key}: {args.get('count')} gene(s) not found.",
        "GENES_PARTIAL": lambda: f"{file_key}: {args.get('count')} gene(s) found in part.",
        "MAX_WINDOWS": lambda: (
            f"{file_key}: {args.get('windows')} window(s) over --max_windows {args.get('max_windows')} are not selected."
        ),
        "TARGET_NOT_ANNOTATED": lambda: f"{file_key} has no CDS annotation; a window may cut a gene.",
    }[warning.code]()


def _report(plan: LocatePlan) -> None:
    """One block per File on standard error: windows, orientation, genes, genes found in part or not at all."""

    reference = plan.reference
    names = {gene.protein_id: gene.cds.name for gene in reference.genes}
    if reference.cut_genes:
        logger.warning(
            "%d reference CDS cut by a window edge are not searched: %s.",
            len(reference.cut_genes), ", ".join(cds.name for cds in reference.cut_genes),
        )
    for file_plan, counts in zip(plan.files, locate_file_counts(plan)):
        coverage = file_plan.coverage
        logger.info(
            "%s: %d window(s); %d of %d genes found, %d in full.",
            file_plan.file_key, len(file_plan.windows), counts.genes_found, counts.genes_total, len(coverage.complete),
        )
        # Each window by the number the drawing gives it (``locate_file_counts``), in that order.
        numbers = dict(counts.window_numbers)
        for number, window in sorted(
            ((numbers[window.window.column or 0], window) for window in file_plan.windows), key=lambda item: item[0]
        ):
            genes = sorted(window.genes, key=lambda gene: gene.order)
            first, last = names[genes[0].protein_id], names[genes[-1].protein_id]
            ends = [label for flag, label in ((window.touches_start, "contig start"), (window.touches_end, "contig end")) if flag]
            logger.info(
                "  %d %s %s: %d gene(s), %s%s.", number, window.record_id, _span(window), len(genes),
                first if first == last else f"{first} to {last}", "".join(f", {end}" for end in ends),
            )
        if coverage.partial:
            logger.info("  Found in part: %s.", ", ".join(names[item] for item in coverage.partial))
        if coverage.absent:
            logger.info("  Not found: %s.", ", ".join(names[item] for item in coverage.absent))
        for warning in file_plan.warnings:
            logger.warning("%s", _warning_text(file_plan.file_key, warning, file_plan))


def locate_main(cmd_args: Sequence[str]) -> None:
    setup_logging()
    parser = _parser()
    args = parser.parse_args(list(cmd_args))
    _input_kind(parser, [args.reference_gbk] if args.reference_gbk else [],
                [args.reference_gff] if args.reference_gff else [],
                [args.reference_fasta] if args.reference_fasta else [], "reference")
    _input_kind(parser, args.targets, args.target_gffs, args.target_fastas, "targets")
    program: LocateProgramName = args.losat
    reference = LocateInputFile(
        str(args.reference_gbk or args.reference_fasta), args.reference_gbk, args.reference_gff, args.reference_fasta
    )
    targets = (
        [LocateInputFile(path, path, None, None) for path in args.targets]
        or [LocateInputFile(fasta, None, gff, fasta) for gff, fasta in zip(args.target_gffs, args.target_fastas)]
        or [LocateInputFile(fasta, None, None, fasta) for fasta in args.target_fastas]
    )
    keys = [target.key for target in targets]
    if len(set(keys)) != len(keys):
        parser.error("Give each target file once.")
    for key in keys:
        if _same_file(key, reference.key):
            raise ValidationError(
                f"Target {key} is the reference file; the reference file is not searched.",
                diagnostic={"code": "INPUT_INVALID", "field": "files"},
            )
    outputs = [Path(str(args.output) + suffix) for suffix in OUTPUT_SUFFIXES]
    raw_names = [f"{index}.{Path(key).name}.{program}.tsv" for index, key in enumerate(keys, start=1)]
    preflight_output_paths(
        [*outputs, *(Path(args.losat_output_dir) / name for name in raw_names if args.losat_output_dir)],
        overwrite=bool(args.overwrite),
    )
    search = prepare_locate_search(
        reference, targets, record=args.reference_record,
        # The key `gbdraw linear` gives the record of its first input.
        record_key=f"{CLI_RECORD_KEY_PREFIX}1",
        windows=_reference_windows(args.reference_region), program=program,
    )
    require_losat_program(
        program, losat_bin=args.losat_bin, ncbi_blast_bin=args.ncbi_blast_bin,
        alternatives=[name for name in LOCATE_PROGRAMS if name != program],
    )
    runtimes: list[LosatRuntimeRecord] = []
    outcomes = []
    for job in search.plan.jobs:
        try:
            text = run_losat_search(
                program, job.query_fasta, job.subject_fasta, options=search.plan.search_args,
                losat_bin=args.losat_bin, ncbi_blast_bin=args.ncbi_blast_bin, threads=args.losat_threads,
                runtime_callback=runtimes.append,
            )
        except ValidationError as exc:
            outcomes.append(LocateSearchResult(job.file_key, "", error=str(exc)))
            continue
        outcomes.append(LocateSearchResult(job.file_key, text))
    options = replace(
        LocateOptions.for_program(program),
        thresholds=ComparisonThresholds(args.evalue, args.bitscore, args.identity, args.alignment_length),
        flank=args.flank, merge_distance=args.merge_distance, max_windows=args.max_windows,
    )
    results = locate_results(search, outcomes)
    plan = finish_locate(search, results, options, runtime=runtimes[0] if runtimes else None)
    _report(plan)
    if not any(file_plan.windows for file_plan in plan.files):
        failed = {fp.file_key: dict(w.args)["error"] for fp in plan.files for w in fp.warnings if w.code == "FILE_SEARCH_FAILED"}
        raise ValidationError(
            ("Every search failed: " if len(failed) == len(plan.files) else "No target file has a window: ")
            + "; ".join(
                f"{fp.file_key}: the search failed: {failed[fp.file_key]}" if fp.file_key in failed
                else f"{fp.file_key} kept {fp.tally.rows_kept:,} of {fp.tally.rows_read:,} rows"
                for fp in plan.files
            )
            + ("." if len(failed) == len(plan.files) else ". Lower the thresholds or choose another reference span."),
            diagnostic={"code": "INPUT_INVALID"},
        )
    files = {item.key: item for item in (reference, *targets)}
    texts = [
        *locate_tables(locate_drawing(plan), files, outputs[0].resolve().parent),
        json.dumps(plan.to_json(), indent=2) + "\n",
    ]
    with staged_output_files(outputs) as staged:
        for path, text in zip(outputs, texts):
            with staged[path].open("w", encoding="utf-8", newline="\n") as handle:
                handle.write(text)
        commit_staged_outputs(staged, overwrite=bool(args.overwrite))
    if args.losat_output_dir:
        rows = {result.file_key: result.text for result in results}
        write_losat_output_files(
            Path(args.losat_output_dir), [(name, rows[key]) for name, key in zip(raw_names, keys)],
            overwrite=bool(args.overwrite),
        )
    for path in outputs:
        logger.info("Wrote %s", path)


__all__ = ["DEFAULT_PREFIX", "OUTPUT_SUFFIXES", "locate_main", "locate_tables"]
