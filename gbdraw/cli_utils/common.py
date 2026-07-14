#!/usr/bin/env python
# coding: utf-8

"""Shared CLI argument definitions and utilities.

This module provides common argument parsers and validation logic
used by both circular and linear diagram commands.
"""

import argparse
import logging
import sys
from typing import Optional

from gbdraw.exceptions import ValidationError
from gbdraw.features.shapes import parse_feature_shape_assignment
from gbdraw.features.visibility import resolve_candidate_feature_types
from gbdraw.io.cli_tables import RecordsTable
from gbdraw.io.genome import load_gbks, load_gff_fasta
from gbdraw.render.export import CAIROSVG_AVAILABLE, has_cairosvg
from gbdraw.render.formats import CAIROSVG_FORMATS, SVG_FORMAT, is_cairosvg_format

logger = logging.getLogger(__name__)

_OUTPUT_FORMAT_HELP = (
    "Comma-separated list of output file formats "
    "(svg, interactive-svg, png, pdf, eps, ps; default: svg; "
    "png/pdf/eps/ps require CairoSVG)."
)


def setup_logging() -> None:
    """Setup logging with INFO level and stdout handler."""
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    if not root_logger.handlers:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter('%(message)s'))
        root_logger.addHandler(handler)


def add_input_args(parser: argparse.ArgumentParser) -> None:
    """Add input file arguments (--gbk, --gff, --fasta)."""
    parser.add_argument(
        "--gbk",
        metavar="GBK_FILE",
        help='GenBank/DDBJ flat file',
        type=str,
        nargs='*')
    parser.add_argument(
        "--gff",
        metavar="GFF3_FILE",
        help="GFF3 file (instead of --gbk; --fasta is required)",
        type=str,
        nargs='*')
    parser.add_argument(
        "--fasta",
        metavar="FASTA_FILE",
        help="FASTA file (required with --gff)",
        type=str,
        nargs='*')


def _add_window_step_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '-w',
        '--window',
        help='window size (optional; default: 1kb for genomes < 1Mb, 10kb for genomes <10Mb, 100kb for genomes >=10Mb)',
        type=int)
    parser.add_argument(
        '-s',
        '--step',
        help='step size (optional; default: 100 bp for genomes < 1Mb, 1kb for genomes <10Mb, 10kb for genomes >=10Mb)',
        type=int)


def _add_feature_shape_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--feature_shape',
        help='Feature shape override (repeatable): TYPE=SHAPE where SHAPE is arrow or rectangle.',
        type=parse_feature_shape_assignment_arg,
        action='append',
        default=[],
        metavar='TYPE=SHAPE',
    )


def _add_block_stroke_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--block_stroke_color',
        help='Block stroke color (str; default: "gray")',
        type=str)
    parser.add_argument(
        '--block_stroke_width',
        help='Block stroke width (optional; float; default: 2 pt for genomes <= 50 kb, 0 pt for genomes >= 50 kb)',
        type=float)


def _add_format_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '-f',
        '--format',
        help=_OUTPUT_FORMAT_HELP,
        type=str,
        default="svg")


def _add_depth_track_label_color_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--depth_track_label',
        metavar='LABEL',
        help='Depth track label(s). Provide one label or one per --depth_track.',
        type=str,
        nargs='+')
    parser.add_argument(
        '--depth_track_color',
        metavar='COLOR',
        help='Depth track fill color(s). Provide one color or one per --depth_track.',
        type=str,
        nargs='+')


def _add_depth_track_tick_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--depth_track_large_tick_interval',
        metavar='VALUE',
        help='Depth track large tick interval(s). Provide one value or one per --depth_track.',
        type=str,
        nargs='+')
    parser.add_argument(
        '--depth_track_small_tick_interval',
        metavar='VALUE',
        help='Depth track small tick interval(s). Provide one value or one per --depth_track.',
        type=str,
        nargs='+')
    parser.add_argument(
        '--depth_track_tick_font_size',
        metavar='VALUE',
        help='Depth track tick font size(s). Provide one value or one per --depth_track.',
        type=str,
        nargs='+')
    parser.add_argument(
        '--show_depth',
        help='Show depth coverage track. Implied when --depth is supplied.',
        action='store_true')


def _add_depth_axis_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--depth_window',
        help='Depth aggregation window size. Defaults to one tenth of the GC/skew window, with a 100 bp minimum.',
        type=int)
    parser.add_argument(
        '--depth_step',
        help='Depth aggregation step size. Defaults to one tenth of the GC/skew step.',
        type=int)
    parser.add_argument(
        '--share_depth_axis',
        help='Use one depth y-axis scale across records.',
        action='store_true')
    parser.add_argument(
        '--depth_min',
        help='Minimum depth for clipping/normalization (optional; must be >= 0).',
        type=float)
    parser.add_argument(
        '--depth_max',
        help='Maximum depth for clipping/normalization (optional; must be >= 0).',
        type=float)
    parser.add_argument(
        '--depth_log_scale',
        dest='depth_normalize',
        help='Render depth coverage on a log10 scale (IGV-style).',
        action='store_true',
        default=None)
    parser.add_argument(
        '--no_depth_log_scale',
        dest='depth_normalize',
        help='Render depth coverage on a linear scale.',
        action='store_false')
    parser.add_argument(
        '--show_depth_axis',
        dest='depth_show_axis',
        help='Show depth coverage axis line, ticks, and labels.',
        action='store_true',
        default=None)
    parser.add_argument(
        '--hide_depth_axis',
        dest='depth_show_axis',
        help='Hide depth coverage axis line, ticks, and labels.',
        action='store_false')
    parser.add_argument(
        '--show_depth_ticks',
        dest='depth_show_ticks',
        help='Show depth coverage axis ticks and labels.',
        action='store_true',
        default=None)
    parser.add_argument(
        '--hide_depth_ticks',
        dest='depth_show_ticks',
        help='Hide depth coverage axis ticks and labels.',
        action='store_false')
    parser.add_argument(
        '--depth_tick_interval',
        help='Depth coverage large tick interval in x coverage units (optional; must be > 0; legacy alias for --depth_large_tick_interval).',
        type=float)
    parser.add_argument(
        '--depth_large_tick_interval',
        help='Depth coverage large tick interval in x coverage units (optional; must be > 0).',
        type=float)
    parser.add_argument(
        '--depth_small_tick_interval',
        help='Depth coverage small tick interval in x coverage units (optional; must be > 0; hidden by default).',
        type=float)
    parser.add_argument(
        '--depth_tick_font_size',
        help='Depth coverage tick label font size (optional; must be > 0).',
        type=float)


def _add_comparison_filter_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--evalue',
        help='Maximum BLAST e-value retained for similarity rings (default: 1e-5).',
        type=float,
        default=1e-5)
    parser.add_argument(
        '--bitscore',
        help='Minimum BLAST bitscore retained for similarity rings (default: 50).',
        type=float,
        default=50.0)
    parser.add_argument(
        '--identity',
        help='Minimum BLAST identity percentage retained for similarity rings (default: 70).',
        type=float,
        default=70.0)
    parser.add_argument(
        '--alignment_length',
        help='Minimum BLAST alignment length retained for similarity rings (default: 0).',
        type=int,
        default=0)


def _add_gc_content_axis_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--gc_content_mode',
        help='GC content display mode: deviation (current centered behavior) or percent (absolute 0-100%% area track).',
        choices=['deviation', 'percent'],
        default=None)
    parser.add_argument(
        '--gc_content_min_percent',
        help='Minimum GC percent for percent-mode clipping/axis (optional; finite number).',
        type=float)
    parser.add_argument(
        '--gc_content_max_percent',
        help='Maximum GC percent for percent-mode clipping/axis (optional; finite number).',
        type=float)
    parser.add_argument(
        '--gc_content_tick_interval',
        help='GC content percent-mode large tick interval (optional; must be > 0; alias for --gc_content_large_tick_interval).',
        type=float)
    parser.add_argument(
        '--gc_content_large_tick_interval',
        help='GC content percent-mode large tick interval (optional; must be > 0).',
        type=float)
    parser.add_argument(
        '--gc_content_small_tick_interval',
        help='GC content percent-mode small tick interval (optional; must be > 0; hidden by default).',
        type=float)
    parser.add_argument(
        '--gc_content_tick_font_size',
        help='GC content percent-mode tick label font size (optional; must be > 0).',
        type=float)
    parser.add_argument(
        '--show_gc_content_axis',
        dest='gc_content_show_axis',
        help='Show GC content percent-mode axis line, ticks, and labels.',
        action='store_true',
        default=None)
    parser.add_argument(
        '--hide_gc_content_axis',
        dest='gc_content_show_axis',
        help='Hide GC content percent-mode axis line, ticks, and labels.',
        action='store_false')
    parser.add_argument(
        '--show_gc_content_ticks',
        dest='gc_content_show_ticks',
        help='Show GC content percent-mode axis ticks and labels.',
        action='store_true',
        default=None)
    parser.add_argument(
        '--hide_gc_content_ticks',
        dest='gc_content_show_ticks',
        help='Hide GC content percent-mode axis ticks and labels.',
        action='store_false')


def _add_legend_size_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--legend_box_size',
        help='Legend box size (optional; float; default: 24 (pixels, 96 dpi) for genomes <= 50 kb, 20 for genomes >= 50 kb).',
        type=float)
    parser.add_argument(
        '--legend_font_size',
        help='Legend font size (optional; float; default: 20 (pt) for genomes <= 50 kb, 16 for genomes >= 50 kb).',
        type=float)


def add_output_args(parser: argparse.ArgumentParser, default_output: Optional[str] = None) -> None:
    """Add output-related arguments (-o, -f)."""
    parser.add_argument(
        '-o',
        '--output',
        help='output file prefix (default: accession number of the sequence)' if default_output is None else f'output file prefix (default: {default_output})',
        type=str,
        default=default_output)

    parser.add_argument(
        '-f',
        '--format',
        help=_OUTPUT_FORMAT_HELP,
        type=str,
        default="svg")


def add_color_args(parser: argparse.ArgumentParser) -> None:
    """Add color-related arguments (-p, -t, -d)."""
    parser.add_argument(
        "-p", "--palette",
        metavar="PALETTE",
        default="default",
        help="Palette name (default: default)",
        type=str)
    parser.add_argument(
        '-t',
        '--table',
        help='color table (optional)',
        type=str,
        default="")
    parser.add_argument(
        '-d',
        '--default_colors',
        help='TSV file that overrides the color palette (optional)',
        type=str,
        default="")


def add_analysis_args(parser: argparse.ArgumentParser) -> None:
    """Add GC/skew analysis arguments (-n, -w, -s)."""
    parser.add_argument(
        '-n',
        '--nt',
        help='dinucleotide (default: GC). ',
        type=str,
        default="GC")
    parser.add_argument(
        '-w',
        '--window',
        help='window size (optional; default: 1kb for genomes < 1Mb, 10kb for genomes <10Mb, 100kb for genomes >=10Mb)',
        type=int)
    parser.add_argument(
        '-s',
        '--step',
        help='step size (optional; default: 100 bp for genomes < 1Mb, 1kb for genomes <10Mb, 10kb for genomes >=10Mb)',
        type=int)


def add_feature_args(parser: argparse.ArgumentParser) -> None:
    """Add feature selection arguments (-k)."""
    parser.add_argument(
        '-k',
        '--features',
        help='Comma-separated list of feature keys to draw (default: CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,repeat_region)',
        type=str,
        default="CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,repeat_region")


def add_stroke_args(parser: argparse.ArgumentParser) -> None:
    """Add stroke styling arguments."""
    parser.add_argument(
        '--block_stroke_color',
        help='Block stroke color (str; default: "gray")',
        type=str)
    parser.add_argument(
        '--block_stroke_width',
        help='Block stroke width (optional; float; default: 2 pt for genomes <= 50 kb, 0 pt for genomes >= 50 kb)',
        type=float)
    parser.add_argument(
        '--line_stroke_color',
        help='Line stroke color (str; default: "gray")',
        type=str)
    parser.add_argument(
        '--line_stroke_width',
        help='Line stroke width (optional; float; default: 5 pt for genomes <= 50 kb, 1 pt for genomes >= 50 kb)',
        type=float)


def add_label_args(parser: argparse.ArgumentParser) -> None:
    """Add label-related arguments."""
    label_list_group = parser.add_mutually_exclusive_group()
    label_list_group.add_argument(
        '--label_whitelist',
        help='path to a file for label whitelisting (optional); mutually exclusive with --label_blacklist',
        type=str,
        default="")
    label_list_group.add_argument(
        '--label_blacklist',
        help='Comma-separated keywords for label blacklisting (optional); mutually exclusive with --label_whitelist',
        type=str,
        default="")
    parser.add_argument(
        '--qualifier_priority',
        help='Path to a TSV file defining qualifier priority for labels (optional)',
        type=str,
        default="")
    parser.add_argument(
        '--label_table',
        help='Path to a TSV file defining post-filter label text overrides (optional)',
        type=str,
        default="")
    parser.add_argument(
        '--feature_visibility_table',
        dest='feature_table',
        help='Path to a TSV file defining per-feature visibility overrides (optional)',
        type=str,
        default="")
    parser.add_argument(
        '--feature_table',
        dest='feature_table',
        help=argparse.SUPPRESS,
        type=str,
        default=argparse.SUPPRESS)


def add_legend_args(parser: argparse.ArgumentParser, choices_help: str = '"right", "left", "none"') -> None:
    """Add legend-related arguments."""
    parser.add_argument(
        '-l',
        '--legend',
        help=f'Legend position (default: "right"; {choices_help})',
        type=str,
        default="right")
    parser.add_argument(
        '--legend_box_size',
        help='Legend box size (optional; float; default: 24 (pixels, 96 dpi) for genomes <= 50 kb, 20 for genomes >= 50 kb).',
        type=float)
    parser.add_argument(
        '--legend_font_size',
        help='Legend font size (optional; float; default: 20 (pt) for genomes <= 50 kb, 16 for genomes >= 50 kb).',
        type=float)


def validate_input_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Validate input file argument combinations."""
    records_table = getattr(args, "records_table", None)
    if records_table and (args.gbk or args.gff or args.fasta):
        parser.error("Error: --records_table cannot be used with --gbk, --gff, or --fasta.")
    if args.gbk and (args.gff or args.fasta):
        parser.error("Error: --gbk cannot be used with --gff or --fasta.")
    if args.gff and not args.fasta:
        parser.error("Error: --gff requires --fasta.")
    if args.fasta and not args.gff:
        parser.error("Error: --fasta requires --gff.")
    if not records_table and not args.gbk and not (args.gff and args.fasta):
        parser.error("Error: Either --records_table, --gbk, or both --gff and --fasta must be provided.")


def validate_label_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Validate label argument combinations."""
    if args.label_whitelist and args.label_blacklist:
        parser.error("Error: --label_whitelist and --label_blacklist are mutually exclusive.")


def parse_feature_shape_assignment_arg(value: str) -> str:
    try:
        parse_feature_shape_assignment(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return value


def handle_output_formats(out_formats: list[str]) -> list[str]:
    """Handle WebAssembly and CairoSVG availability for output formats."""
    cairo_formats = [f for f in out_formats if is_cairosvg_format(f)]
    if "pyodide" in sys.modules:
        if cairo_formats:
            logger.info("Running in WebAssembly mode: Image conversion is handled by the browser.")
            return [f for f in out_formats if not is_cairosvg_format(f)] or [SVG_FORMAT]
    else:
        if cairo_formats and not has_cairosvg():
            logger.warning(
                f"CairoSVG is not installed. Cannot generate: {', '.join(cairo_formats).upper()}\n"
                f"   Output restricted to SVG-compatible formats only.\n"
                f"   (To enable PNG/PDF, run: pip install gbdraw[export])"
            )
            return [f for f in out_formats if f not in CAIROSVG_FORMATS] or [SVG_FORMAT]
    return out_formats


def calculate_window_step(seq_length: int, cfg, manual_window: Optional[int], manual_step: Optional[int]) -> tuple[int, int]:
    """Calculate window and step sizes based on genome length."""
    if not manual_window:
        if seq_length < 1_000_000:
            window = cfg.objects.sliding_window.default[0]
        elif seq_length < 10_000_000:
            window = cfg.objects.sliding_window.up1m[0]
        else:
            window = cfg.objects.sliding_window.up10m[0]
    else:
        window = manual_window

    if not manual_step:
        if seq_length < 1_000_000:
            step = cfg.objects.sliding_window.default[1]
        elif seq_length < 10_000_000:
            step = cfg.objects.sliding_window.up1m[1]
        else:
            step = cfg.objects.sliding_window.up10m[1]
    else:
        step = manual_step

    return window, step


def load_records_table_records(
    records_table: RecordsTable,
    *,
    mode: str,
    selected_features_set,
    color_table,
    feature_table,
    gbk_loader=load_gbks,
    gff_loader=load_gff_fasta,
) -> list:
    records: list = []
    if records_table.input_kind == "gbk":
        for row in records_table.rows:
            loaded = gbk_loader(
                [row.gbk],
                mode,
                False,
                record_selectors=[row.record_id],
                reverse_flags=[row.reverse_complement],
            )
            records.append(_require_one_records_table_record(records_table, row, loaded))
        return records

    candidate_feature_types, keep_all_features = resolve_candidate_feature_types(
        selected_features_set,
        color_table=color_table,
        feature_visibility_table=feature_table,
    )
    for row in records_table.rows:
        loaded = gff_loader(
            [row.gff],
            [row.fasta],
            mode,
            candidate_feature_types,
            keep_all_features=keep_all_features,
            load_comparison=False,
            record_selectors=[row.record_id],
            reverse_flags=[row.reverse_complement],
        )
        records.append(_require_one_records_table_record(records_table, row, loaded))
    return records


def _require_one_records_table_record(records_table: RecordsTable, row, loaded):
    if len(loaded) == 1:
        return loaded[0]
    if len(loaded) > 1 and not row.record_id:
        raise ValidationError(
            f"{records_table.table_path}: row {row.row_number} loaded {len(loaded)} records; "
            "add record_id so the row selects exactly one displayed record."
        )
    raise ValidationError(
        f"{records_table.table_path}: row {row.row_number} must load exactly one displayed record; "
        f"loaded {len(loaded)}."
    )


def record_major_depth_track_files_from_cli(
    depth_track_groups: list[list[str]] | None,
    *,
    record_count: int,
) -> list[list[str | None]] | None:
    if not depth_track_groups:
        return None
    rows: list[list[str | None]] = [[] for _ in range(record_count)]
    for track_number, group in enumerate(depth_track_groups, start=1):
        values = [
            None
            if str(path).strip().lower() in {"", "-", "none", "null"}
            else str(path)
            for path in (group or [])
        ]
        if not values or all(value is None for value in values):
            raise ValidationError(f"--depth_track #{track_number} must include at least one file.")
        if len(values) == 1:
            expanded = values * record_count
        elif len(values) == record_count:
            expanded = values
        else:
            raise ValidationError(
                f"--depth_track #{track_number} must contain one file or one per record ({record_count}); got {len(values)}."
            )
        for record_index, path in enumerate(expanded):
            rows[record_index].append(path)
    return rows


__all__ = [
    "CAIROSVG_AVAILABLE",
    "add_analysis_args",
    "add_color_args",
    "add_feature_args",
    "add_input_args",
    "add_label_args",
    "add_legend_args",
    "add_output_args",
    "add_stroke_args",
    "calculate_window_step",
    "handle_output_formats",
    "setup_logging",
    "validate_input_args",
    "validate_label_args",
]
