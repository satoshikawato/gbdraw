#!/usr/bin/env python
# coding: utf-8


import argparse
import logging
import math
import sys
from pathlib import Path
from typing import Optional
from pandas import DataFrame  # type: ignore[reportMissingImports]
from .io.colors import load_default_colors, read_color_table
from .io.genome import load_gbks, load_gff_fasta
from .io.regions import apply_region_specs, parse_region_specs
from .config.toml import load_config_toml
from .render.export import parse_formats, save_figure
from .api.diagram import assemble_linear_diagram_from_records  # type: ignore[reportMissingImports]
from .analysis.collinearity import (
    LosslessCollinearityParameters,
    build_orthogroup_collinearity_blocks,
    convert_collinearity_blocks_to_comparisons,
    normalize_collinearity_anchor_mode,
    normalize_collinearity_search_scope,
)
from .analysis.protein_colinearity import (
    ORTHOGROUP_INFERENCE_VERSION,
    ORTHOGROUP_MEMBERSHIP_MODES,
    PROTEIN_BLASTP_MODES,
    normalize_orthogroup_membership_mode,
)
from .io.collinearity import parse_native_collinearity_tsv, write_native_collinearity_tsv
from .config.modify import modify_config_dict  # type: ignore[reportMissingImports]
from .config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from .labels.filtering import (
    read_filter_list_file,
    read_label_override_file,
    read_qualifier_priority_file,
)  # type: ignore[reportMissingImports]
from .features.shapes import parse_feature_shape_assignment, parse_feature_shape_overrides
from .features.visibility import read_feature_visibility_file
from .exceptions import ValidationError
from .tracks import (
    linear_track_slots_from_order,
    normalize_linear_track_slots_with_axis,
    parse_linear_track_slots,
)


from .cli_utils.common import (
    setup_logging,
    validate_input_args,
    validate_label_args,
    handle_output_formats,
    calculate_window_step,
)


def _parse_optional_positive_int(value: str) -> int | None:
    normalized = str(value or "").strip().lower()
    if normalized in {"", "none", "auto", "null"}:
        return None
    try:
        parsed = int(normalized)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a positive integer or 'none'") from exc
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer or 'none'")
    return parsed


def _parse_optional_nonnegative_float(value: str) -> float | None:
    normalized = str(value or "").strip().lower()
    if normalized in {"", "none", "null"}:
        return None
    try:
        parsed = float(normalized)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a non-negative number or 'none'") from exc
    if not math.isfinite(parsed) or parsed < 0:
        raise argparse.ArgumentTypeError("must be a non-negative number or 'none'")
    return parsed


# Setup for the logging system
logger = logging.getLogger()
setup_logging()

def _parse_linear_label_placement(value: str) -> str:
    normalized = str(value).strip().lower()
    if normalized == "on_feature":
        return "above_feature"
    if normalized in {"auto", "above_feature"}:
        return normalized
    raise argparse.ArgumentTypeError(
        "label placement must be 'auto' or 'above_feature' (legacy alias: 'on_feature')"
    )


def _parse_linear_track_layout(value: str) -> str:
    normalized = str(value).strip().lower()
    if normalized in {"above", "spreadout"}:
        return "above"
    if normalized in {"below", "tuckin"}:
        return "below"
    if normalized == "middle":
        return "middle"
    raise argparse.ArgumentTypeError(
        "track layout must be 'above', 'middle', or 'below' "
        "(legacy aliases: 'spreadout' -> 'above', 'tuckin' -> 'below')"
    )


def _parse_linear_track_axis_gap(value: str) -> float | None:
    normalized = str(value).strip().lower()
    if normalized in {"", "auto"}:
        return None
    try:
        axis_gap = float(normalized)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "track axis gap must be a non-negative number or 'auto'"
        ) from exc
    if axis_gap < 0:
        raise argparse.ArgumentTypeError("track axis gap must be >= 0")
    return axis_gap


def _parse_pairwise_match_style(value: str) -> str:
    normalized = str(value).strip().lower()
    if normalized not in {"ribbon", "curve"}:
        raise argparse.ArgumentTypeError("pairwise_match_style must be one of: ribbon, curve")
    return normalized


def _parse_collinear_color_mode(value: str) -> str:
    normalized = str(value).strip().lower().replace("-", "_")
    if normalized == "identity":
        normalized = "average_identity"
    if normalized not in {"average_identity", "orientation", "orientation_identity"}:
        raise argparse.ArgumentTypeError(
            "collinear_color_mode must be one of: average_identity, orientation, orientation_identity"
        )
    return normalized


def _parse_collinear_anchor_mode(value: str) -> str:
    try:
        normalize_collinearity_anchor_mode(value)
    except ValidationError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return "rbh"


def _parse_collinear_search_scope(value: str) -> str:
    try:
        return normalize_collinearity_search_scope(value)
    except ValidationError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def _parse_orthogroup_membership_mode(value: str) -> str:
    try:
        return normalize_orthogroup_membership_mode(value)
    except ValidationError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def _parse_feature_shape_assignment_arg(value: str) -> str:
    try:
        parse_feature_shape_assignment(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return value


def _get_args(args) -> argparse.Namespace:
    """
    Parses command-line arguments for generating linear genome diagrams.

    This internal function defines and parses command-line arguments using argparse.
    It sets up the necessary parameters required for the linear genome diagram creation,
    including input files, output configurations, and various visualization options.

    Args:
        args (list of str): The command-line arguments passed to the script.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments.

    The function supports a variety of arguments for input files (GenBank and BLAST),
    output Configurator (file format and naming), and visualization preferences (color
    tables, window size, and feature selection).
    """
    parser = argparse.ArgumentParser(
        description='Generate  plot in PNG/PDF/SVG/PS/EPS.')
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
    parser.add_argument(
        '-b',
        '--blast',
        help="input BLAST result file in tab-separated format (-outfmt 6 or 7) (optional)",
        type=str,
        nargs='*')
    parser.add_argument(
        '--losatp_bin',
        '--losatp-bin',
        dest='losatp_bin',
        help='LOSATP executable for --protein_blastp_mode pairwise/orthogroup/collinear (default: losat).',
        type=str,
        default='losat')
    parser.add_argument(
        '--losatp_threads',
        '--losatp-threads',
        dest='losatp_threads',
        help='Threads passed to LOSATP via --num-threads for --protein_blastp_mode pairwise/orthogroup/collinear (default: LOSAT default).',
        type=int,
        default=None)
    parser.add_argument(
        '--protein_blastp_mode',
        '--protein-blastp-mode',
        dest='protein_blastp_mode',
        help='LOSATP blastp mode: none, pairwise adjacent ribbons, all-record Orthogroups, or Collinear blocks (default: none).',
        choices=PROTEIN_BLASTP_MODES,
        default='none')
    parser.add_argument(
        '--protein_blastp_max_hits',
        '--protein-blastp-max-hits',
        dest='protein_blastp_max_hits',
        help='Maximum distinct subject protein hits per query protein for pairwise LOSATP blastp display links (default: 5).',
        type=int,
        default=5)
    parser.add_argument(
        '--protein_blastp_candidate_limit',
        '--protein-blastp-candidate-limit',
        dest='protein_blastp_candidate_limit',
        help="Optional LOSATP blastp candidate cap per query; use 'none' for no cap (default: none).",
        type=_parse_optional_positive_int,
        default=None)
    parser.add_argument(
        '--orthogroup_membership_mode',
        '--orthogroup-membership-mode',
        dest='orthogroup_membership_mode',
        help='Orthogroup inference model for LOSATP Orthogroup/Collinear modes. Legacy values rbh, family_merge, and distribution_split are accepted as aliases for anchor_core_v1 (default: anchor_core_v1).',
        type=_parse_orthogroup_membership_mode,
        choices=ORTHOGROUP_MEMBERSHIP_MODES,
        default=ORTHOGROUP_INFERENCE_VERSION)
    parser.add_argument(
        '--orthogroup_member_max_hits',
        '--orthogroup-member-max-hits',
        dest='orthogroup_member_max_hits',
        help='Deprecated compatibility option; anchor_core_v1 inference no longer caps membership evidence with this value (default: 5).',
        type=int,
        default=5)
    parser.add_argument(
        '--align_orthogroup_feature',
        '--align-orthogroup-feature',
        dest='align_orthogroup_feature',
        help='Align linear records by the LOSATP blastp orthogroup containing this feature SVG hash or protein ID.',
        type=str,
        default="")
    parser.add_argument(
        '--collinear_unit_mode',
        '--collinear-unit-mode',
        dest='collinear_unit_mode',
        help=argparse.SUPPRESS,
        choices=["auto", "cds", "locus"],
        default='auto')
    parser.add_argument(
        '--collinear_anchor_mode',
        '--collinear-anchor-mode',
        '--collinear_orthogroup_edge_mode',
        '--collinear-orthogroup-edge-mode',
        dest='collinear_anchor_mode',
        help='Deprecated compatibility option. Collinear blocks always use RBH anchors (default: rbh).',
        type=_parse_collinear_anchor_mode,
        choices=["all", "one_to_one", "rbh"],
        default='rbh')
    parser.add_argument(
        '--collinear_search_scope',
        '--collinear-search-scope',
        dest='collinear_search_scope',
        help='Collinear LOSATP evidence search scope: adjacent record pairs or all record pairs (default: adjacent).',
        type=_parse_collinear_search_scope,
        choices=["adjacent", "all"],
        default='adjacent')
    parser.add_argument(
        '--collinear_min_anchors',
        '--collinear-min-anchors',
        dest='collinear_min_anchors',
        help='Minimum anchors/genes required for a rendered Collinear block; 1 allows singleton links (default: 1).',
        type=int,
        default=1)
    parser.add_argument(
        '--collinear_max_unit_gap',
        '--collinear-max-unit-gap',
        '--collinear_max_gene_gap',
        '--collinear-max-gene-gap',
        dest='collinear_max_unit_gap',
        help='Maximum unit gap between neighboring collinear anchors (default: 0).',
        type=int,
        default=0)
    parser.add_argument(
        '--collinear_block_merge_gap',
        '--collinear-block-merge-gap',
        dest='collinear_block_merge_gap',
        help=argparse.SUPPRESS,
        type=int,
        default=50)
    parser.add_argument(
        '--collinear_singleton_merge_gap',
        '--collinear-singleton-merge-gap',
        dest='collinear_singleton_merge_gap',
        help=argparse.SUPPRESS,
        type=int,
        default=25)
    parser.add_argument(
        '--collinear_max_diagonal_drift',
        '--collinear-max-diagonal-drift',
        dest='collinear_max_diagonal_drift',
        help='Maximum order-space diagonal drift allowed within or between collinear runs (default: 0).',
        type=int,
        default=0)
    parser.add_argument(
        '--collinear_max_conflicts_in_merge_gap',
        '--collinear-max-conflicts-in-merge-gap',
        dest='collinear_max_conflicts_in_merge_gap',
        help=argparse.SUPPRESS,
        type=int,
        default=1)
    parser.add_argument(
        '--collinear_max_paralog_links_per_orthogroup',
        '--collinear-max-paralog-links-per-orthogroup',
        dest='collinear_max_paralog_links_per_orthogroup',
        help=argparse.SUPPRESS,
        type=int,
        default=2)
    parser.add_argument(
        '--collinear_gap_penalty',
        '--collinear-gap-penalty',
        dest='collinear_gap_penalty',
        help=argparse.SUPPRESS,
        type=float,
        default=1.0)
    parser.add_argument(
        '--collinear_nearby_duplicate_window',
        '--collinear-nearby-duplicate-window',
        dest='collinear_nearby_duplicate_window',
        help=argparse.SUPPRESS,
        type=int,
        default=0)
    parser.add_argument(
        '--collinear_score_mode',
        '--collinear-score-mode',
        dest='collinear_score_mode',
        help=argparse.SUPPRESS,
        choices=["constant", "bitscore"],
        default='constant')
    parser.add_argument(
        '--collinear_constant_anchor_score',
        '--collinear-constant-anchor-score',
        dest='collinear_constant_anchor_score',
        help=argparse.SUPPRESS,
        type=float,
        default=50.0)
    parser.add_argument(
        '--collinear_min_block_score',
        '--collinear-min-block-score',
        dest='collinear_min_block_score',
        help=argparse.SUPPRESS,
        type=float,
        default=None)
    parser.add_argument(
        '--collinear_block_evalue',
        '--collinear-block-evalue',
        dest='collinear_block_evalue',
        help=argparse.SUPPRESS,
        type=_parse_optional_nonnegative_float,
        default=None)
    parser.add_argument(
        '--collinear_color_mode',
        '--collinear-color-mode',
        dest='collinear_color_mode',
        help='Collinear ribbon color mode: average_identity, orientation, or orientation_identity (default: orientation).',
        type=_parse_collinear_color_mode,
        choices=["average_identity", "orientation", "orientation_identity"],
        default='orientation')
    parser.add_argument(
        '--collinear_blocks',
        '--collinear-blocks',
        dest='collinear_blocks',
        help='Headered native .collinear.tsv file to import instead of running LOSATP.',
        type=str,
        default="")
    parser.add_argument(
        '--save_collinear_blocks',
        '--save-collinear-blocks',
        dest='save_collinear_blocks',
        help='Write accepted or validated native collinear blocks to this TSV path.',
        type=str,
        default="")
    parser.add_argument(
        '-t',
        '--table',
        help='color table (optional)',
        type=str,
        default="")
    parser.add_argument(
        "-p", "--palette",
        metavar="PALETTE",
        default="default",
        help="Palette name (default: default)",
        type=str
    )
    parser.add_argument(
        '-d',
        '--default_colors',
        help='TSV file that overrides the color palette (optional)',
        type=str,
        default="")
    parser.add_argument(
        '-o',
        '--output',
        help='output file prefix (default: out)',
        type=str,
        default="out")
    parser.add_argument(
        '-n',
        '--nt',
        help='dinucleotide skew (default: GC). ',
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
    parser.add_argument(
        '--separate_strands',
        help='separate forward and reverse strands (default: False). Features of undefined strands are shown on the forward strand.',
        action='store_true')
    parser.add_argument(
        '--show_gc',
        help='plot GC content below genome (default: False). ',
        action='store_true')
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
    parser.add_argument(
        '--show_skew',
        help='plot GC skew below genome (default: False). ',
        action='store_true')
    parser.add_argument(
        '--depth',
        help='Depth TSV file(s) in samtools depth format. Provide one for all records or one per input record.',
        type=str,
        nargs='+')
    parser.add_argument(
        '--depth_track',
        metavar='DEPTH',
        help='Repeatable logical depth track. Provide one file for all records or one file per input record.',
        type=str,
        nargs='+',
        action='append')
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
    parser.add_argument(
        '--depth_track_height',
        metavar='PX',
        help='Linear depth track height(s) in px. Provide one value or one per --depth_track.',
        type=str,
        nargs='+')
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
    parser.add_argument(
        '--depth_color',
        help='Depth track fill color (optional; default: #4A90E2).',
        type=str)
    parser.add_argument(
        '--depth_height',
        help='Depth track height for linear mode (in px; must be > 0).',
        type=float)
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
    parser.add_argument(
        '--align_center',
        help='Align genomes to the center (default: False). ',
        action='store_true')
    parser.add_argument(
        '--keep_definition_left_aligned',
        '--keep-definition-left-aligned',
        dest='keep_definition_left_aligned',
        help='Keep linear definition labels in the left column when records are center-aligned or aligned by orthogroup (default: False).',
        action='store_true')
    parser.add_argument(
        '--evalue',
        help='evalue threshold (default=1e-2)',
        type=float,
        default="1e-2")
    parser.add_argument(
        '--bitscore',
        help='bitscore threshold (default=50)',
        type=float,
        default="50")
    parser.add_argument(
        '--identity',
        help='identity threshold (default=0)',
        type=float,
        default="0")
    parser.add_argument(
        '--alignment_length',
        help='minimum BLAST alignment length threshold (default=0)',
        type=int,
        default=0)
    parser.add_argument(
        '--pairwise_match_style',
        '--pairwise-match-style',
        dest='pairwise_match_style',
        help=(
            'Pairwise comparison link style: ribbon keeps straight filled ribbons; '
            'curve draws curved filled ribbons that preserve alignment spans.'
        ),
        type=_parse_pairwise_match_style,
        choices=["ribbon", "curve"],
        default="ribbon")
    parser.add_argument(
        '-k',
        '--features',
        help='Comma-separated list of feature keys to draw (default: CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,repeat_region)',
        type=str,
        default="CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,repeat_region")
    parser.add_argument(
        '--feature_shape',
        help='Feature shape override (repeatable): TYPE=SHAPE where SHAPE is arrow or rectangle.',
        type=_parse_feature_shape_assignment_arg,
        action='append',
        default=[],
        metavar='TYPE=SHAPE',
    )
    parser.add_argument(
        '--block_stroke_color',
        help='Block stroke color (str; default: "gray")',
        type=str)
    parser.add_argument(
        '--block_stroke_width',
        help='Block stroke width (optional; float; default: 2 pt for genomes <= 50 kb, 0 pt for genomes >= 50 kb)',
        type=float)
    parser.add_argument(
        '--axis_stroke_color',
        help='Axis stroke color (str; default: auto: "lightgray", or "dimgray" with --ruler_on_axis)',
        type=str,
        default=None)
    parser.add_argument(
        '--axis_stroke_width',
        help='Axis stroke width (optional; float; default: 5 pt for genomes <= 50 kb, 2 pt for genomes >= 50 kb)',
        type=float)
    parser.add_argument(
        '--line_stroke_color',
        help='Line stroke color (optional; str; default: "lightgray")',
        type=str)
    parser.add_argument(
        '--line_stroke_width',
        help='Line stroke width (optional; float; default: 5 pt for genomes <= 50 kb, 1 pt for genomes >= 50 kb)',
        type=float)
    parser.add_argument(
        '--definition_font_size',
        help='Definition font size (optional; float; default: 24 pt for genomes <= 50 kb, 10 pt for genomes >= 50 kb)',
        type=float)
    parser.add_argument(
        '--plot_title',
        help='Shared plot title text (optional).',
        type=str,
        default="")
    parser.add_argument(
        '--plot_title_position',
        help='Shared plot title position ("center", "top", "bottom"; default: "bottom").',
        type=str,
        choices=["center", "top", "bottom"],
        default="bottom")
    parser.add_argument(
        '--plot_title_font_size',
        help='Shared plot title font size (optional; float; default: 32).',
        type=float)
    parser.add_argument(
        '--record_label',
        help='Optional top definition line (for example organism/strain; repeatable; order matches input records)',
        type=str,
        action='append',
        default=[])
    parser.add_argument(
        '--show_replicon',
        help='Show inferred replicon labels in linear record definitions (default: False).',
        action='store_true')
    parser.add_argument(
        '--hide_accession',
        help='Hide accession labels in linear record definitions (default: False).',
        action='store_true')
    parser.add_argument(
        '--hide_length',
        help='Hide length/coordinate labels in linear record definitions (default: False).',
        action='store_true')
    parser.add_argument(
        '--label_font_size',
        help='Label font size (optional; default: 24 pt for genomes <= 50 kb, 5 pt for genomes >= 50 kb)',
        type=float)
    parser.add_argument(
        '--label_placement',
        help='Linear label placement mode ("auto" or "above_feature"; default: "auto"). "above_feature" draws labels above features (or below negative-strand features when --separate_strands is used).',
        type=_parse_linear_label_placement,
        metavar="{auto,above_feature}",
    )
    parser.add_argument(
        '--label_rendering',
        help='Label rendering policy: "auto" embeds fitting labels and routes others externally; "embedded_only" drops external labels; "external_only" forces labels outside feature bodies. Cannot be combined with --label_placement above_feature except as "auto". Default: "auto".',
        choices=['auto', 'embedded_only', 'external_only'],
        default='auto',
        type=str)
    parser.add_argument(
        '--label_rotation',
        help='Linear label rotation in degrees (optional; float; default: 0). In above_feature mode, rotated labels start from the feature midpoint.',
        type=float,
    )
    parser.add_argument(
        '--linear_label_spacing',
        help='Linear label-to-label vertical spacing in px (optional; float; must be > 0).',
        type=float,
    )
    parser.add_argument(
        '--track_layout',
        help=(
            'Linear track layout mode ("above", "middle", or "below"; default: "middle"). '
            'Aliases: "spreadout" -> "above", "tuckin" -> "below".'
        ),
        type=_parse_linear_track_layout,
        metavar="{above,middle,below}",
        default="middle",
    )
    parser.add_argument(
        '--track_axis_gap',
        help=(
            "Gap between axis and nearest feature edge in pixels for above/below layouts. "
            "Use 'auto' to derive it from feature height."
        ),
        type=_parse_linear_track_axis_gap,
        default=None,
        metavar="AUTO|PX",
    )
    parser.add_argument(
        '--linear_track_order',
        help='Linear custom track shortcut order, for example features,depth,gc_content,gc_skew.',
        type=str,
        default="",
    )
    parser.add_argument(
        '--linear_track_slot',
        help='Linear custom track slot: <slot_id>:<renderer>@key=value,key=value. Repeat to add slots.',
        action='append',
        default=[],
        metavar='SLOT',
    )
    parser.add_argument(
        '--linear_track_axis_index',
        help='Axis boundary index for linear custom track slots.',
        type=int,
        default=None,
    )
    parser.add_argument(
        '--ruler_on_axis',
        help=(
            'Use each record axis as the ruler in linear mode. '
            'Effective only with --scale_style ruler and --track_layout above|below.'
        ),
        action='store_true',
    )
    parser.add_argument(
        '-f',
        '--format',
        help='Comma-separated list of output file formats (svg, png, pdf, eps, ps; default: svg; non-SVG requires CairoSVG).',
        type=str,
        default="svg")
    parser.add_argument(
        '-l',
        '--legend',
        help='Legend position (default: "right"; "right", "left", "top", "bottom", "none")',
        type=str,
        default="right")
    parser.add_argument(
            "--show_labels",
            help="Show labels: no argument or 'all' (all records), 'first' (first record only), 'none' (no labels). Default: 'none'",
            nargs='?',
            const="all",
            default="none",
            choices=["all", "first", "none"],
            type=str
        )
    parser.add_argument(
        '--resolve_overlaps',
        help='Resolve overlaps (experimental; default: False). ',
        action='store_true')
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
    parser.add_argument(
        '--label_table',
        help='Path to a TSV file defining post-filter label text overrides (optional)',
        type=str,
        default="")
    parser.add_argument(
        '--feature_table',
        help='Path to a TSV file defining per-feature visibility overrides (optional)',
        type=str,
        default="")
    parser.add_argument(
        '--feature_height',
        help='Feature vertical width (optional; float; default: 80 (pixels, 96 dpi) for genomes <= 50 kb, 20 for genomes >= 50 kb)',
        type=float),
    parser.add_argument(
        '--gc_height',
        help='GC content/skew vertical width (optional; float; default: 20 (pixels, 96 dpi))',
        type=float),
    parser.add_argument(
        '--comparison_height',
        help='Comparison block height (optional; float; optional; default: 60 (pixels, 96 dpi))',
        type=float)
    parser.add_argument(
        '--scale_style',
        help='Style for the length scale (default: "bar"; "bar", "ruler")',
        type=str,
        choices=["bar", "ruler"],
        default="bar")
    parser.add_argument(
        '--scale_stroke_color',
        help='Scale bar/ruler stroke color (optional; str; default: "black")',
        type=str)
    parser.add_argument(
        '--scale_stroke_width',
        help='Scale bar/ruler stroke width (optional; float; default: 3 (pt))',
        type=float)
    parser.add_argument(
        '--scale_font_size',
        help='Scale bar/ruler font size (optional; float; default: 24 (pt) for genomes <= 50 kb, 16 for genomes >= 50 kb).',
        type=float)
    parser.add_argument(
        '--ruler_label_font_size',
        help='Ruler label font size (optional; float). Overrides --scale_font_size when both are set.',
        type=float)
    parser.add_argument(
        '--ruler_label_color',
        help='Ruler label color (optional; str; default follows axis color when --ruler_on_axis is active, otherwise black).',
        type=str)
    parser.add_argument(
        '--scale_interval',
        help='Manual tick interval for "ruler" scale style (in bp). Overrides automatic calculation; optional',
        type=int)
    parser.add_argument(
        '--legend_box_size',
        help='Legend box size (optional; float; default: 24 (pixels, 96 dpi) for genomes <= 50 kb, 20 for genomes >= 50 kb).',
        type=float)
    parser.add_argument(
        '--legend_font_size',
        help='Legend font size (optional; float; default: 20 (pt) for genomes <= 50 kb, 16 for genomes >= 50 kb).',
        type=float)
    parser.add_argument(
        '--normalize_length',
        help='Normalize record length (experimental; default: False). ',
        action='store_true')
    parser.add_argument(
        '--region',
        help=(
            'Crop a region (repeatable). Format: record_id:start-end[:rc], #index:start-end[:rc], '
            'or file:record_selector:start-end[:rc]. '
            'Coordinates are 1-based inclusive. For multiple records without selectors, provide one spec per record in input order (file order, then record order within each file).'
        ),
        type=str,
        action='append',
        default=[])
    parser.add_argument(
        '--record_id',
        help=(
            'Select a record by ID or #index per input file (repeatable; order matches input files). '
            'Use an empty value to skip selection for a file.'
        ),
        type=str,
        action='append',
        default=[])
    parser.add_argument(
        '--reverse_complement',
        help=(
            'Reverse complement record per input file (repeatable; order matches input files). '
            'Accepted values: 1/0, true/false, yes/no.'
        ),
        type=str,
        action='append',
        default=[])
    args = parser.parse_args(args)
    validate_input_args(parser, args)
    validate_label_args(parser, args)
    if args.protein_blastp_mode != "none" and args.blast:
        parser.error("--protein_blastp_mode cannot be used with -b/--blast")
    if args.collinear_blocks and args.blast:
        parser.error("--collinear_blocks cannot be used with -b/--blast")
    if args.collinear_blocks and args.protein_blastp_mode != "none":
        parser.error("--collinear_blocks imports native blocks and cannot be used with --protein_blastp_mode")
    if args.save_collinear_blocks and not args.collinear_blocks and args.protein_blastp_mode != "collinear":
        parser.error("--save_collinear_blocks requires --protein_blastp_mode collinear or --collinear_blocks")
    if args.protein_blastp_max_hits <= 0:
        parser.error("--protein_blastp_max_hits must be > 0")
    if args.orthogroup_member_max_hits <= 0:
        parser.error("--orthogroup_member_max_hits must be > 0")
    if args.losatp_threads is not None and args.losatp_threads <= 0:
        parser.error("--losatp_threads must be > 0")
    if args.align_orthogroup_feature and args.protein_blastp_mode != "orthogroup" and not args.blast:
        parser.error("--align_orthogroup_feature requires --protein_blastp_mode orthogroup")
    if args.depth and args.depth_track:
        parser.error("--depth cannot be combined with --depth_track")
    if args.show_depth and not (args.depth or args.depth_track):
        parser.error("--show_depth requires --depth or --depth_track")
    if args.depth_height is not None and args.depth_height <= 0:
        parser.error("--depth_height must be > 0")
    if args.depth_window is not None and args.depth_window <= 0:
        parser.error("--depth_window must be > 0")
    if args.depth_step is not None and args.depth_step <= 0:
        parser.error("--depth_step must be > 0")
    if args.depth_min is not None and args.depth_min < 0:
        parser.error("--depth_min must be >= 0")
    if args.depth_max is not None and args.depth_max < 0:
        parser.error("--depth_max must be >= 0")
    if args.depth_min is not None and args.depth_max is not None and args.depth_min > args.depth_max:
        parser.error("--depth_min must be <= --depth_max")
    if args.depth_tick_interval is not None and args.depth_tick_interval <= 0:
        parser.error("--depth_tick_interval must be > 0")
    if args.depth_large_tick_interval is not None and args.depth_large_tick_interval <= 0:
        parser.error("--depth_large_tick_interval must be > 0")
    if args.depth_small_tick_interval is not None and args.depth_small_tick_interval <= 0:
        parser.error("--depth_small_tick_interval must be > 0")
    if args.depth_tick_font_size is not None and args.depth_tick_font_size <= 0:
        parser.error("--depth_tick_font_size must be > 0")
    for option_name in (
        "depth_track_height",
        "depth_track_large_tick_interval",
        "depth_track_small_tick_interval",
        "depth_track_tick_font_size",
    ):
        for option_value in getattr(args, option_name) or []:
            option_text = str(option_value).strip().lower()
            if option_text in {"", "auto", "none", "null", "-"}:
                continue
            try:
                numeric_option_value = float(option_value)
            except (TypeError, ValueError):
                parser.error(f"--{option_name} values must be numbers or auto")
            if numeric_option_value <= 0:
                parser.error(f"--{option_name} values must be > 0")
    for option_name in ("gc_content_min_percent", "gc_content_max_percent"):
        option_value = getattr(args, option_name)
        if option_value is not None and not math.isfinite(float(option_value)):
            parser.error(f"--{option_name} must be a finite number")
    if (
        args.gc_content_min_percent is not None
        and args.gc_content_max_percent is not None
        and args.gc_content_min_percent > args.gc_content_max_percent
    ):
        parser.error("--gc_content_min_percent must be <= --gc_content_max_percent")
    if args.gc_content_tick_interval is not None and args.gc_content_tick_interval <= 0:
        parser.error("--gc_content_tick_interval must be > 0")
    if args.gc_content_large_tick_interval is not None and args.gc_content_large_tick_interval <= 0:
        parser.error("--gc_content_large_tick_interval must be > 0")
    if args.gc_content_small_tick_interval is not None and args.gc_content_small_tick_interval <= 0:
        parser.error("--gc_content_small_tick_interval must be > 0")
    if args.gc_content_tick_font_size is not None and args.gc_content_tick_font_size <= 0:
        parser.error("--gc_content_tick_font_size must be > 0")
    if args.linear_track_order and args.linear_track_slot:
        parser.error("--linear_track_order cannot be combined with --linear_track_slot")
    if args.linear_track_axis_index is not None and not (args.linear_track_order or args.linear_track_slot):
        parser.error("--linear_track_axis_index requires --linear_track_order or --linear_track_slot")
    try:
        linear_track_slot_specs = None
        if args.linear_track_order:
            linear_track_slot_specs = linear_track_slots_from_order(
                args.linear_track_order,
                show_depth=bool(args.show_depth or args.depth or args.depth_track),
                depth_track_count=max(1, len(args.depth_track or [])),
                show_gc=bool(args.show_gc),
                show_skew=bool(args.show_skew),
                dinucleotide=str(args.nt or "GC").upper(),
                track_layout=args.track_layout,
            )
        elif args.linear_track_slot:
            linear_track_slot_specs = parse_linear_track_slots(args.linear_track_slot)
        if linear_track_slot_specs is not None and args.linear_track_axis_index is not None:
            normalize_linear_track_slots_with_axis(
                linear_track_slot_specs,
                args.linear_track_axis_index,
            )
        args.linear_track_slot_specs = linear_track_slot_specs
    except Exception as exc:
        parser.error(str(exc))
    if args.label_placement == "above_feature" and args.label_rendering != "auto":
        parser.error("--label_rendering embedded_only|external_only cannot be used with --label_placement above_feature")
    if args.collinear_min_anchors <= 0:
        parser.error("--collinear_min_anchors must be > 0")
    if args.collinear_max_unit_gap < 0:
        parser.error("--collinear_max_unit_gap must be >= 0")
    if args.collinear_block_merge_gap < 0:
        parser.error("--collinear_block_merge_gap must be >= 0")
    if args.collinear_singleton_merge_gap < 0:
        parser.error("--collinear_singleton_merge_gap must be >= 0")
    if args.collinear_max_diagonal_drift < 0:
        parser.error("--collinear_max_diagonal_drift must be >= 0")
    if args.collinear_max_conflicts_in_merge_gap < 0:
        parser.error("--collinear_max_conflicts_in_merge_gap must be >= 0")
    if args.collinear_max_paralog_links_per_orthogroup <= 0:
        parser.error("--collinear_max_paralog_links_per_orthogroup must be > 0")
    if args.collinear_gap_penalty < 0:
        parser.error("--collinear_gap_penalty must be >= 0")
    if args.collinear_nearby_duplicate_window < 0:
        parser.error("--collinear_nearby_duplicate_window must be >= 0")
    if args.collinear_constant_anchor_score <= 0:
        parser.error("--collinear_constant_anchor_score must be > 0")
    return args


def _record_major_depth_track_files_from_cli(
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


def linear_main(cmd_args) -> None:
    """
    Main function for generating linear genome diagrams.

    This function orchestrates the creation of linear genome diagrams by processing
    input GenBank files and optional BLAST comparison files. It leverages various
    configurations and Configurator provided via command-line arguments to customize the
    visualization, such as color schemes, genome feature selections, and layout options.

    Args:
        cmd_args (argparse.Namespace): Command-line arguments parsed by argparse, providing
                                       specifications for input files, output formats, and
                                       visualization Configurator.

    The function performs the following key steps:
    - Loading and validating input files and Configurator.
    - Configuring the visualization canvas and feature Configurator.
    - Executing the plotting process to generate the linear diagrams.
    - Handling output in specified file formats.

    The final output includes linear genome diagrams in user-specified file formats,
    illustrating genomic features and optional BLAST comparison results.
    """
    args: argparse.Namespace = _get_args(cmd_args)
    if '-i' in cmd_args or '--input' in cmd_args:
        logger.warning(
            "WARNING: The -i/--input option is deprecated and will be removed in a future version. Please use --gbk instead.")
    out_file_prefix: str = args.output
    blast_files: str = args.blast
    protein_blastp_mode: str = str(args.protein_blastp_mode or "none")
    losatp_bin: str = args.losatp_bin
    losatp_threads: int | None = args.losatp_threads
    protein_blastp_max_hits: int = args.protein_blastp_max_hits
    protein_blastp_candidate_limit: int | None = args.protein_blastp_candidate_limit
    orthogroup_membership_mode: str = str(args.orthogroup_membership_mode or ORTHOGROUP_INFERENCE_VERSION)
    orthogroup_member_max_hits: int = args.orthogroup_member_max_hits
    align_orthogroup_feature: str = str(args.align_orthogroup_feature or "").strip()
    collinear_unit_mode: str = str(args.collinear_unit_mode or "auto")
    collinear_anchor_mode: str = str(args.collinear_anchor_mode or "rbh")
    collinear_search_scope: str = str(args.collinear_search_scope or "adjacent")
    collinear_color_mode: str = str(args.collinear_color_mode or "orientation")
    collinear_blocks_path: str = str(args.collinear_blocks or "").strip()
    save_collinear_blocks_path: str = str(args.save_collinear_blocks or "").strip()
    collinearity_params = LosslessCollinearityParameters(
        min_anchors=args.collinear_min_anchors,
        max_unit_gap=args.collinear_max_unit_gap,
        max_diagonal_drift=args.collinear_max_diagonal_drift,
    )
    color_table_path: str = args.table
    strandedness: bool = args.separate_strands
    resolve_overlaps: bool = args.resolve_overlaps
    dinucleotide: str = args.nt.upper()
    show_gc: bool = args.show_gc
    gc_content_mode: str | None = args.gc_content_mode
    gc_content_min_percent: Optional[float] = args.gc_content_min_percent
    gc_content_max_percent: Optional[float] = args.gc_content_max_percent
    gc_content_show_axis: bool | None = args.gc_content_show_axis
    gc_content_show_ticks: bool | None = args.gc_content_show_ticks
    gc_content_tick_interval: Optional[float] = args.gc_content_tick_interval
    gc_content_large_tick_interval: Optional[float] = args.gc_content_large_tick_interval
    gc_content_small_tick_interval: Optional[float] = args.gc_content_small_tick_interval
    gc_content_tick_font_size: Optional[float] = args.gc_content_tick_font_size
    manual_window: int = args.window
    manual_step: int = args.step
    depth_files: list[str] | None = args.depth
    depth_track_groups: list[list[str]] | None = args.depth_track
    depth_track_labels: list[str] | None = list(args.depth_track_label or []) or None
    depth_track_colors: list[str] | None = list(args.depth_track_color or []) or None
    depth_track_heights: list[str] | None = list(args.depth_track_height or []) or None
    depth_track_large_tick_intervals: list[str] | None = list(args.depth_track_large_tick_interval or []) or None
    depth_track_small_tick_intervals: list[str] | None = list(args.depth_track_small_tick_interval or []) or None
    depth_track_tick_font_sizes: list[str] | None = list(args.depth_track_tick_font_size or []) or None
    show_depth: bool = bool(args.show_depth or depth_files or depth_track_groups)
    depth_color: str | None = args.depth_color
    depth_height: Optional[float] = args.depth_height
    depth_window: Optional[int] = args.depth_window
    depth_step: Optional[int] = args.depth_step
    depth_share_axis: bool = bool(args.share_depth_axis)
    depth_min: Optional[float] = args.depth_min
    depth_max: Optional[float] = args.depth_max
    depth_normalize: bool | None = args.depth_normalize
    depth_show_axis: bool | None = args.depth_show_axis
    depth_show_ticks: bool | None = args.depth_show_ticks
    depth_tick_interval: Optional[float] = args.depth_tick_interval
    depth_large_tick_interval: Optional[float] = args.depth_large_tick_interval
    depth_small_tick_interval: Optional[float] = args.depth_small_tick_interval
    depth_tick_font_size: Optional[float] = args.depth_tick_font_size
    align_center: bool = args.align_center
    keep_definition_left_aligned: bool = args.keep_definition_left_aligned
    evalue: float = args.evalue
    legend: str = args.legend
    gc_height: Optional[float] = args.gc_height
    show_skew: bool = args.show_skew
    bitscore: float = args.bitscore
    identity: float = args.identity
    alignment_length: int = args.alignment_length
    pairwise_match_style: str = args.pairwise_match_style
    show_labels: str = args.show_labels
    label_whitelist: str = args.label_whitelist
    label_blacklist: str = args.label_blacklist
    qualifier_priority_path: str = args.qualifier_priority
    label_table_path: str = args.label_table
    feature_table_path: str = args.feature_table
    selected_features_set: str = args.features.split(',')
    feature_shapes = parse_feature_shape_overrides(args.feature_shape)
    feature_height: Optional[float] = args.feature_height
    comparison_height: Optional[float] = args.comparison_height

    out_formats: list[str] = parse_formats(args.format)
    out_formats = handle_output_formats(out_formats)
    user_defined_default_colors: str = args.default_colors
    scale_style: str = args.scale_style
    scale_stroke_color: Optional[str] = args.scale_stroke_color
    scale_stroke_width: Optional[float] = args.scale_stroke_width
    scale_font_size: Optional[float] = args.scale_font_size
    ruler_label_font_size: Optional[float] = args.ruler_label_font_size
    effective_ruler_label_font_size: Optional[float] = (
        ruler_label_font_size if ruler_label_font_size is not None else scale_font_size
    )
    scale_label_color: Optional[str] = args.ruler_label_color
    scale_interval: Optional[int] = args.scale_interval
    legend_box_size: Optional[float] = args.legend_box_size
    legend_font_size: Optional[float] = args.legend_font_size
    normalize_length = args.normalize_length
    if alignment_length < 0:
        raise ValidationError("alignment_length must be >= 0")
    if blast_files or protein_blastp_mode != "none" or collinear_blocks_path:
        load_comparison = True
    else:
        load_comparison = False
    palette: str = args.palette
    default_colors: Optional[DataFrame] = load_default_colors(
        user_defined_default_colors, palette, load_comparison)
    color_table: Optional[DataFrame] = read_color_table(color_table_path)
    feature_table: Optional[DataFrame] = read_feature_visibility_file(feature_table_path)
    config_dict: dict = load_config_toml('gbdraw.data', 'config.toml')

    filtering_cfg = config_dict.setdefault("labels", {}).setdefault("filtering", {})
    if qualifier_priority_path:
        filtering_cfg["qualifier_priority_df"] = read_qualifier_priority_file(qualifier_priority_path)
    else:
        filtering_cfg["qualifier_priority_df"] = None
    if label_whitelist:
        filtering_cfg["whitelist_df"] = read_filter_list_file(label_whitelist)
    else:
        filtering_cfg["whitelist_df"] = None
    if label_table_path:
        filtering_cfg["label_override_df"] = read_label_override_file(label_table_path)
    else:
        filtering_cfg["label_override_df"] = None

    block_stroke_color: Optional[str] = args.block_stroke_color
    block_stroke_width: Optional[float] = args.block_stroke_width
    definition_font_size: Optional[float] = args.definition_font_size
    definition_show_replicon: bool = bool(args.show_replicon)
    definition_show_accession: bool = not bool(args.hide_accession)
    definition_show_length: bool = not bool(args.hide_length)
    plot_title: str = str(args.plot_title or "").strip()
    plot_title_position: str = str(args.plot_title_position or "bottom").strip().lower()
    plot_title_font_size: Optional[float] = args.plot_title_font_size
    label_font_size: Optional[float] = args.label_font_size
    label_rendering: str = args.label_rendering
    label_placement: Optional[str] = args.label_placement
    label_rotation: Optional[float] = args.label_rotation
    track_layout: str = args.track_layout
    track_axis_gap: Optional[float] = args.track_axis_gap
    linear_track_slot_specs = args.linear_track_slot_specs
    linear_track_axis_index: int | None = args.linear_track_axis_index
    ruler_on_axis: bool = bool(args.ruler_on_axis)
    if ruler_on_axis and not (scale_style == "ruler" and track_layout in {"above", "below"}):
        logger.warning(
            "WARNING: --ruler_on_axis is ignored unless --scale_style ruler and --track_layout above|below are set."
        )
        ruler_on_axis = False
    axis_stroke_color: str = (
        args.axis_stroke_color
        if args.axis_stroke_color is not None
        else ("dimgray" if ruler_on_axis else "lightgray")
    )
    if ruler_on_axis and scale_stroke_color is None:
        scale_stroke_color = axis_stroke_color
    if ruler_on_axis and scale_label_color is None:
        scale_label_color = axis_stroke_color
    axis_stroke_width: Optional[float] = args.axis_stroke_width
    line_stroke_color: Optional[str] = args.line_stroke_color
    line_stroke_width: Optional[float] = args.line_stroke_width
    if plot_title_font_size is not None and float(plot_title_font_size) <= 0:
        raise ValidationError("plot_title_font_size must be > 0")
    if args.linear_label_spacing is not None and float(args.linear_label_spacing) <= 0:
        raise ValidationError("linear_label_spacing must be > 0")
    config_dict = modify_config_dict(
        config_dict,
        block_stroke_color=block_stroke_color,
        block_stroke_width=block_stroke_width,
        linear_axis_stroke_color=axis_stroke_color,
        linear_axis_stroke_width=axis_stroke_width,
        linear_definition_font_size=definition_font_size,
        linear_definition_show_replicon=definition_show_replicon,
        linear_definition_show_accession=definition_show_accession,
        linear_definition_show_length=definition_show_length,
        label_font_size=label_font_size,
        linear_label_spacing=args.linear_label_spacing,
        label_rendering=label_rendering,
        label_placement=label_placement,
        label_rotation=label_rotation,
        line_stroke_color=line_stroke_color,
        line_stroke_width=line_stroke_width,
        show_gc=show_gc,
        gc_content_mode=gc_content_mode,
        gc_content_min_percent=gc_content_min_percent,
        gc_content_max_percent=gc_content_max_percent,
        gc_content_show_axis=gc_content_show_axis,
        gc_content_show_ticks=gc_content_show_ticks,
        gc_content_tick_interval=gc_content_tick_interval,
        gc_content_large_tick_interval=gc_content_large_tick_interval,
        gc_content_small_tick_interval=gc_content_small_tick_interval,
        gc_content_tick_font_size=gc_content_tick_font_size,
        show_skew=show_skew,
        show_depth=show_depth,
        depth_color=depth_color,
        depth_height=depth_height,
        depth_min=depth_min,
        depth_max=depth_max,
        depth_normalize=depth_normalize,
        depth_show_axis=depth_show_axis,
        depth_show_ticks=depth_show_ticks,
        depth_tick_interval=depth_tick_interval,
        depth_large_tick_interval=depth_large_tick_interval,
        depth_small_tick_interval=depth_small_tick_interval,
        depth_tick_font_size=depth_tick_font_size,
        depth_share_axis=depth_share_axis,
        show_labels=show_labels,
        resolve_overlaps=resolve_overlaps,
        linear_track_layout=track_layout,
        linear_track_axis_gap=track_axis_gap,
        linear_ruler_on_axis=ruler_on_axis,
        align_center=align_center,
        keep_definition_left_aligned=keep_definition_left_aligned,
        strandedness=strandedness,
        label_blacklist=label_blacklist,
        label_whitelist=label_whitelist,
        label_table=label_table_path,
        default_cds_height=feature_height,
        comparison_height=comparison_height,
        gc_height=gc_height,
        scale_style=scale_style,
        scale_stroke_color=scale_stroke_color,
        scale_label_color=scale_label_color,
        scale_stroke_width=scale_stroke_width,
        scale_font_size=scale_font_size,
        ruler_label_font_size=effective_ruler_label_font_size,
        scale_interval=scale_interval,
        legend_box_size=legend_box_size,
        legend_font_size=legend_font_size,
        pairwise_match_style=pairwise_match_style,
        normalize_length=normalize_length
        )

    def _normalize_list(values, target_len, fill_value=""):
        items = list(values or [])
        if len(items) > target_len:
            logger.error(
                "ERROR: Too many --record_id/--reverse_complement values (expected at most %s).", target_len
            )
            raise ValidationError(
                f"Too many --record_id/--reverse_complement values (expected at most {target_len})."
            )
        while len(items) < target_len:
            items.append(fill_value)
        return items

    def _parse_bool(value: str | None) -> bool:
        if value is None:
            return False
        text = str(value).strip().lower()
        if text in {"1", "true", "yes", "y", "on"}:
            return True
        if text in {"0", "false", "no", "n", "off", "", "none", "null", "-"}:
            return False
        logger.error("ERROR: Invalid reverse_complement value: %s", value)
        raise ValidationError(f"Invalid reverse_complement value: {value}")

    if args.gbk:
        file_count = len(args.gbk)
        record_selectors = _normalize_list(args.record_id, file_count, "")
        reverse_flags_raw = _normalize_list(args.reverse_complement, file_count, "")
        reverse_flags = [_parse_bool(v) for v in reverse_flags_raw]
        records = load_gbks(
            args.gbk,
            "linear",
            load_comparison,
            record_selectors=record_selectors,
            reverse_flags=reverse_flags,
        )
    elif args.gff and args.fasta:
        file_count = len(args.gff)
        record_selectors = _normalize_list(args.record_id, file_count, "")
        reverse_flags_raw = _normalize_list(args.reverse_complement, file_count, "")
        reverse_flags = [_parse_bool(v) for v in reverse_flags_raw]
        records = load_gff_fasta(
            args.gff,
            args.fasta,
            "linear",
            selected_features_set,
            keep_all_features=bool(feature_table_path),
            load_comparison=load_comparison,
            record_selectors=record_selectors,
            reverse_flags=reverse_flags,
        )
    else:
        logger.error("A critical error occurred with input file arguments.")
        raise ValidationError("Invalid input file arguments.")
    record_labels = args.record_label or []
    if record_labels:
        if len(record_labels) > len(records):
            logger.warning(
                "WARNING: More --record_label values were provided than records loaded; extra labels will be ignored."
            )
        for idx, label in enumerate(record_labels[: len(records)]):
            if label is None:
                continue
            label = str(label).strip()
            if not label:
                continue
            if getattr(records[idx], "annotations", None) is None:
                records[idx].annotations = {}
            records[idx].annotations["gbdraw_record_label"] = label
    region_specs = parse_region_specs(args.region)
    if region_specs:
        try:
            records = apply_region_specs(records, region_specs, log=logger)
        except ValueError as exc:
            logger.error(f"ERROR: {exc}")
            raise ValidationError(str(exc)) from exc
        if blast_files:
            logger.warning(
                "WARNING: Region cropping is enabled; ensure BLAST coordinates match the cropped regions (and reverse complements if specified)."
            )
    if protein_blastp_mode != "none" and len(records) < 2:
        raise ValidationError("--protein_blastp_mode requires at least two linear records.")
    if collinear_blocks_path and len(records) < 2:
        raise ValidationError("--collinear_blocks requires at least two linear records.")
    depth_track_files = _record_major_depth_track_files_from_cli(
        depth_track_groups,
        record_count=len(records),
    )
    collinearity_comparisons: list[DataFrame] | None = None
    if collinear_blocks_path:
        collinearity_result = parse_native_collinearity_tsv(
            collinear_blocks_path,
            records,
            params=collinearity_params,
            unit_mode=collinear_unit_mode,
        )
        collinearity_comparisons = convert_collinearity_blocks_to_comparisons(
            collinearity_result,
            records=records,
            color_mode=collinear_color_mode,
        )
    elif protein_blastp_mode == "collinear":
        collinearity_result = build_orthogroup_collinearity_blocks(
            records,
            losatp_bin=losatp_bin,
            losatp_threads=losatp_threads,
            candidate_limit=protein_blastp_candidate_limit,
            orthogroup_membership_mode=orthogroup_membership_mode,
            orthogroup_member_max_hits=orthogroup_member_max_hits,
            max_paralog_links_per_orthogroup=args.collinear_max_paralog_links_per_orthogroup,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
            params=collinearity_params,
            unit_mode=collinear_unit_mode,
            edge_mode=collinear_anchor_mode,
            search_scope=collinear_search_scope,
        )
        collinearity_comparisons = convert_collinearity_blocks_to_comparisons(
            collinearity_result,
            records=records,
            color_mode=collinear_color_mode,
        )
    if save_collinear_blocks_path:
        Path(save_collinear_blocks_path).write_text(
            write_native_collinearity_tsv(collinearity_result),
            encoding="utf-8",
        )
    # Use raw records to avoid collapsing lengths when IDs are duplicated.
    longest_genome: int = max(len(record.seq) for record in records)
    cfg = GbdrawConfig.from_dict(config_dict)
    window, step = calculate_window_step(longest_genome, cfg, manual_window, manual_step)

    canvas = assemble_linear_diagram_from_records(
        records=records,
        blast_files=blast_files,
        config_dict=config_dict,
        color_table=color_table,
        default_colors=default_colors,
        selected_features_set=selected_features_set,
        feature_table=feature_table,
        feature_shapes=feature_shapes or None,
        output_prefix=out_file_prefix,
        legend=legend,
        dinucleotide=dinucleotide,
        window=window,
        step=step,
        depth_window=depth_window,
        depth_step=depth_step,
        depth_files=depth_files,
        depth_track_files=depth_track_files,
        depth_track_labels=depth_track_labels,
        depth_track_colors=depth_track_colors,
        depth_track_heights=depth_track_heights,
        depth_track_large_tick_intervals=depth_track_large_tick_intervals,
        depth_track_small_tick_intervals=depth_track_small_tick_intervals,
        depth_track_tick_font_sizes=depth_track_tick_font_sizes,
        linear_track_slots=linear_track_slot_specs,
        linear_track_axis_index=linear_track_axis_index,
        plot_title=plot_title,
        plot_title_position=plot_title_position,
        plot_title_font_size=plot_title_font_size,
        protein_comparisons=collinearity_comparisons,
        protein_blastp_mode="none" if collinearity_comparisons is not None else protein_blastp_mode,
        losatp_bin=losatp_bin,
        losatp_threads=losatp_threads,
        protein_blastp_max_hits=protein_blastp_max_hits,
        protein_blastp_candidate_limit=protein_blastp_candidate_limit,
        orthogroup_membership_mode=orthogroup_membership_mode,
        orthogroup_member_max_hits=orthogroup_member_max_hits,
        collinear_max_paralog_links_per_orthogroup=args.collinear_max_paralog_links_per_orthogroup,
        align_orthogroup_feature=align_orthogroup_feature or None,
        pairwise_match_style=pairwise_match_style,
        evalue=evalue,
        bitscore=bitscore,
        identity=identity,
        alignment_length=alignment_length,
        cfg=cfg,
    )
    save_figure(canvas, out_formats)


if __name__ == "__main__":
    # This gets all arguments passed to the script, excluding the script name
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(handler)
    main_args = sys.argv[1:]
    if not main_args:
        main_args.append('--help')
    linear_main(main_args)
