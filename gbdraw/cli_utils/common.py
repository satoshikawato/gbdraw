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

try:
    import cairosvg  # type: ignore[reportMissingImports]
    CAIROSVG_AVAILABLE = True
except (ImportError, OSError):
    CAIROSVG_AVAILABLE = False


logger = logging.getLogger(__name__)


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
        help='Genbank/DDBJ flatfile',
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


def add_output_args(parser: argparse.ArgumentParser, default_output: Optional[str] = None) -> None:
    """Add output-related arguments (-o, -f)."""
    parser.add_argument(
        '-o',
        '--output',
        help='output file prefix (default: accession number of the sequence)' if default_output is None else f'output file prefix (default: {default_output})',
        type=str,
        default=default_output)

    if CAIROSVG_AVAILABLE:
        parser.add_argument(
            '-f',
            '--format',
            help='Comma-separated list of output file formats (svg, png, pdf, eps, ps; default: svg).',
            type=str,
            default="svg")
    else:
        parser.add_argument(
            '-f',
            '--format',
            help='Comma-separated list of output file formats (svg; install CairoSVG to enable png, pdf, eps, ps output).',
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
        help='Comma-separated keywords or path to a file for label blacklisting (optional); mutually exclusive with --label_whitelist',
        type=str,
        default="")
    parser.add_argument(
        '--qualifier_priority',
        help='Path to a TSV file defining qualifier priority for labels (optional)',
        type=str,
        default="")


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
    if args.gbk and (args.gff or args.fasta):
        parser.error("Error: --gbk cannot be used with --gff or --fasta.")
    if args.gff and not args.fasta:
        parser.error("Error: --gff requires --fasta.")
    if args.fasta and not args.gff:
        parser.error("Error: --fasta requires --gff.")
    if not args.gbk and not (args.gff and args.fasta):
        parser.error("Error: Either --gbk or both --gff and --fasta must be provided.")


def validate_label_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Validate label argument combinations."""
    if args.label_whitelist and args.label_blacklist:
        parser.error("Error: --label_whitelist and --label_blacklist are mutually exclusive.")


def handle_output_formats(out_formats: list[str]) -> list[str]:
    """Handle WebAssembly and CairoSVG availability for output formats."""
    if "pyodide" in sys.modules:
        if any(f != 'svg' for f in out_formats):
            logger.info("Running in WebAssembly mode: Output format constrained to SVG. (Image conversion is handled by the browser)")
            return ['svg']
    elif not CAIROSVG_AVAILABLE:
        non_svg_formats = [f for f in out_formats if f != 'svg']
        if non_svg_formats:
            logger.warning(
                f"CairoSVG is not installed. Cannot generate: {', '.join(non_svg_formats).upper()}\n"
                f"   Output restricted to SVG only.\n"
                f"   (To enable PNG/PDF, run: pip install gbdraw[export])"
            )
            return ['svg']
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
