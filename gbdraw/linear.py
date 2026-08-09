#!/usr/bin/env python
# coding: utf-8


import argparse
import copy
import logging
import math
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Mapping, Optional, Sequence
from .config.toml import load_config_toml
from .render.export import parse_formats
from .api.request_render import CurrentRequestArtifacts, render_request
from .api.session_compat import render_session_compatible_request
from .api.record_planning import (
    depth_track_inputs_from_cli,
    record_input_manifest_from_paths,
    record_input_manifest_from_table,
)
from .api.options import (
    AnnotationOptions,
    ColorOptions,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    LinearOutputOptions,
    LinearRequestTrackOptions,
)
from .layout.record_placement import (
    parse_record_row_position,
)
from .api.requests import (
    LinearDiagramRequest,
    RecordCardinality,
    RenderOutputRequest,
)
from .definition_line_styles import (
    parse_definition_line_style_assignment,
    parse_definition_line_style_overrides,
)
from .analysis.collinearity import (
    LosslessCollinearityParameters,
    normalize_collinearity_color_mode,
    normalize_collinearity_search_scope,
)
from .analysis.protein_colinearity import (
    ORTHOGROUP_INFERENCE_VERSION,
    PROTEIN_BLASTP_MODES,
    hydrate_protein_losat_tsv,
    is_protein_losat_cache_entry,
)
from .config.modify import modify_config_dict  # type: ignore[reportMissingImports]
from .config.models.objects import normalize_pairwise_match_style
from .features.shapes import parse_feature_shape_overrides
from .exceptions import ValidationError
from .mode_profiles import ComparisonThresholds, LINEAR_MODE_PROFILE
from .tracks import (
    linear_track_slots_from_order,
    normalize_linear_track_slots_with_axis,
    parse_linear_track_slots,
)


from .cli_utils.common import (
    _add_arrow_geometry_args,
    _add_block_stroke_args,
    _add_depth_axis_args,
    _add_depth_track_arg,
    _add_depth_track_label_color_args,
    _add_depth_track_tick_args,
    _add_feature_shape_arg,
    _add_format_arg,
    _add_overwrite_arg,
    _add_gc_skew_toggle_args,
    _add_gc_content_axis_args,
    _add_legend_size_args,
    _add_window_step_args,
    add_feature_args,
    add_input_args,
    add_label_args,
    setup_logging,
    validate_input_args,
    validate_label_args,
    handle_output_formats,
)
from .cli_utils.session import (
    DiagramRunResult,
    add_session_args,
    diagram_request_output_paths,
    make_rendered_svg,
    parse_session_pre_args,
    preflight_session_sidecar_if_requested,
    render_canonical_session_if_present,
    save_session_sidecar_if_requested,
)
from .render.track_slot_metadata import (
    build_track_slot_geometry_run_metadata,
    collect_track_slot_geometry_records,
)
from .render.output_paths import commit_staged_output_file, preflight_output_paths
from .session_io import load_session, session_to_cli_args


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


# Setup for the logging system
logger = logging.getLogger()
setup_logging()


def _linear_cli_record_cardinality(
    *,
    is_gff_source: bool,
    source_count: int,
    load_comparison: bool,
) -> RecordCardinality:
    """State the legacy CLI multi-entry selection rule as typed policy."""

    if load_comparison and (is_gff_source or source_count > 1):
        return RecordCardinality.FIRST
    return RecordCardinality.ALL


def _parse_linear_label_placement(value: str) -> str:
    normalized = str(value).strip().lower()
    if normalized in {"auto", "above_feature"}:
        return normalized
    raise argparse.ArgumentTypeError(
        "label placement must be 'auto' or 'above_feature'"
    )


def _parse_linear_track_layout(value: str) -> str:
    normalized = str(value).strip().lower()
    if normalized in {"above", "middle", "below"}:
        return normalized
    raise argparse.ArgumentTypeError(
        "track layout must be 'above', 'middle', or 'below'"
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
    try:
        return normalize_pairwise_match_style(value)
    except ValidationError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def _parse_collinear_color_mode(value: str) -> str:
    try:
        return normalize_collinearity_color_mode(value)
    except ValidationError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def _parse_collinear_search_scope(value: str) -> str:
    try:
        return normalize_collinearity_search_scope(value)
    except ValidationError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc




def _parse_definition_line_style_arg(value: str) -> str:
    try:
        parse_definition_line_style_assignment(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return value


def _parse_linear_multi_record_position(value: str) -> str:
    try:
        parse_record_row_position(value)
    except ValidationError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return str(value).strip()


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
    add_input_args(parser)
    parser.add_argument(
        "--records_table",
        metavar="TSV",
        help="TSV manifest for row-based input records and per-record options.",
        type=str)
    parser.add_argument(
        "--multi_record_position",
        metavar="SELECTOR@ROW",
        action="append",
        default=[],
        type=_parse_linear_multi_record_position,
        help="Place a record in a Linear row; repeat once for every loaded record.",
    )
    parser.add_argument(
        "--linear_record_gap",
        metavar="PX",
        type=float,
        default=24.0,
        help="Gap in pixels between records in the same Linear row (default: 24).",
    )
    parser.add_argument(
        "--comparisons_table",
        metavar="TSV",
        type=str,
        help="TSV manifest with blast, query, and subject columns.",
    )
    parser.add_argument(
        '-b',
        '--blast',
        help="input BLAST result file in tab-separated format (-outfmt 6 or 7) (optional)",
        type=str,
        nargs='*')
    parser.add_argument(
        '--losatp_bin',
        dest='losatp_bin',
        help='Native LOSAT executable for --protein_blastp_mode pairwise/orthogroup/collinear (default: losat).',
        type=str,
        default='losat')
    parser.add_argument(
        '--ncbi_blastp_bin',
        dest='ncbi_blastp_bin',
        help='NCBI BLAST+ blastp executable for --protein_blastp_mode pairwise/orthogroup/collinear (default: use automatic runtime resolution).',
        type=str,
        default=None)
    parser.add_argument(
        '--losatp_threads',
        dest='losatp_threads',
        help='Threads passed to the selected protein blastp runtime for --protein_blastp_mode pairwise/orthogroup/collinear (default: runtime default).',
        type=int,
        default=None)
    parser.add_argument(
        '--protein_blastp_mode',
        dest='protein_blastp_mode',
        help='Protein blastp comparison mode: none, pairwise adjacent ribbons, all-record similarity groups (orthogroup), or collinear blocks (default: none).',
        choices=PROTEIN_BLASTP_MODES,
        default='none')
    parser.add_argument(
        '--protein_blastp_max_hits',
        dest='protein_blastp_max_hits',
        help='Maximum distinct subject protein hits per query protein for pairwise protein blastp display links (default: 5).',
        type=int,
        default=5)
    parser.add_argument(
        '--protein_blastp_candidate_limit',
        dest='protein_blastp_candidate_limit',
        help="Optional protein blastp candidate cap per query; use 'none' for no cap (default: none).",
        type=_parse_optional_positive_int,
        default=None)
    parser.add_argument(
        '--protein_blastp_output',
        metavar='TSV',
        help=(
            'Write the raw protein-search evidence to one deterministic TSV. '
            'Runtime handles are replaced with user-visible protein IDs; requires '
            '--protein_blastp_mode.'
        ),
        type=str,
        default=None)
    parser.add_argument(
        '--align_orthogroup_feature',
        dest='align_orthogroup_feature',
        help='Align linear records by the gbdraw similarity group containing this feature SVG hash or protein ID.',
        type=str,
        default="")
    parser.add_argument(
        '--collinear_unit_mode',
        dest='collinear_unit_mode',
        help=argparse.SUPPRESS,
        choices=["auto", "cds", "locus"],
        default='auto')
    parser.add_argument(
        '--collinear_search_scope',
        dest='collinear_search_scope',
        help=(
            'Collinear protein blastp scope: adjacent displayed records/rows or all record pairs. '
            'With multi-record rows, adjacent searches every cross-record pair between '
            'neighboring rows; all uses every pair as grouping evidence but renders only '
            'pairs across adjacent rows '
            '(default: adjacent).'
        ),
        type=_parse_collinear_search_scope,
        choices=["adjacent", "all"],
        default='adjacent')
    parser.add_argument(
        '--collinear_min_anchors',
        dest='collinear_min_anchors',
        help='Minimum anchors/genes required for a rendered Collinear block; 1 allows singleton links (default: 1).',
        type=int,
        default=1)
    parser.add_argument(
        '--collinear_max_unit_gap',
        dest='collinear_max_unit_gap',
        help='Maximum unit gap between neighboring collinear anchors (default: 0).',
        type=int,
        default=0)
    parser.add_argument(
        '--collinear_max_diagonal_drift',
        dest='collinear_max_diagonal_drift',
        help='Maximum order-space diagonal drift allowed within or between collinear runs (default: 0).',
        type=int,
        default=0)
    parser.add_argument(
        '--collinear_max_conflicts_in_merge_gap',
        dest='collinear_max_conflicts_in_merge_gap',
        help=argparse.SUPPRESS,
        type=int,
        default=1)
    parser.add_argument(
        '--collinear_max_paralog_links_per_orthogroup',
        dest='collinear_max_paralog_links_per_orthogroup',
        help=argparse.SUPPRESS,
        type=int,
        default=2)
    parser.add_argument(
        '--collinear_color_mode',
        dest='collinear_color_mode',
        help='Collinear ribbon color mode: average_identity, orientation, or orientation_identity (default: orientation).',
        type=_parse_collinear_color_mode,
        choices=["average_identity", "orientation", "orientation_identity"],
        default='orientation')
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
    _add_overwrite_arg(parser)
    parser.add_argument(
        '-n',
        '--nt',
        help='dinucleotide skew (default: GC). ',
        type=str,
        default="GC")
    _add_window_step_args(parser)
    parser.add_argument(
        '--separate_strands',
        help='separate forward and reverse strands (default: False). Features of undefined strands are shown on the forward strand.',
        action='store_true')
    _add_gc_skew_toggle_args(
        parser,
        show_gc_default=LINEAR_MODE_PROFILE.show_gc,
        show_skew_default=LINEAR_MODE_PROFILE.show_skew,
    )
    _add_gc_content_axis_args(parser)
    _add_depth_track_arg(parser, mode="linear")
    _add_depth_track_label_color_args(parser)
    parser.add_argument(
        '--depth_track_height',
        metavar='PX',
        help='Linear depth track height(s) in px. Provide one value or one per --depth_track.',
        type=str,
        nargs='+')
    _add_depth_track_tick_args(parser)
    parser.add_argument(
        '--depth_color',
        help='Depth track fill color (optional; default: #4A90E2).',
        type=str)
    parser.add_argument(
        '--depth_height',
        help='Depth track height for linear mode (in px; must be > 0).',
        type=float)
    _add_depth_axis_args(parser)
    parser.add_argument(
        '--align_center',
        help='Align genomes to the center (default: False). ',
        action='store_true')
    parser.add_argument(
        '--keep_definition_left_aligned',
        dest='keep_definition_left_aligned',
        help=(
            'Keep linear record definitions in the left column. With multi-record rows, '
            'the leading record label becomes the row definition while remaining record '
            'text stays above its record (default: False).'
        ),
        action='store_true')
    parser.add_argument(
        '--evalue',
        help='evalue threshold (default=1e-2)',
        type=float,
        default=LINEAR_MODE_PROFILE.comparison.evalue)
    parser.add_argument(
        '--bitscore',
        help='bitscore threshold (default=50)',
        type=float,
        default=LINEAR_MODE_PROFILE.comparison.bitscore)
    parser.add_argument(
        '--identity',
        help='identity threshold (default=0)',
        type=float,
        default=LINEAR_MODE_PROFILE.comparison.identity)
    parser.add_argument(
        '--alignment_length',
        help='minimum BLAST alignment length threshold (default=0)',
        type=int,
        default=LINEAR_MODE_PROFILE.comparison.alignment_length)
    parser.add_argument(
        '--pairwise_match_style',
        dest='pairwise_match_style',
        help=(
            'Pairwise comparison link style: ribbon keeps straight filled ribbons; '
            'curve draws curved filled ribbons that preserve alignment spans.'
        ),
        type=_parse_pairwise_match_style,
        choices=["ribbon", "curve"],
        default="ribbon")
    add_feature_args(parser)
    _add_feature_shape_arg(parser)
    _add_arrow_geometry_args(parser)
    _add_block_stroke_args(parser)
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
        '--definition_line_style',
        dest='definition_line_style',
        help=(
            'Definition line style override (repeatable): '
            'LINE:weight=bold,color=#000000,size=12. '
            'Lines: name/species/record_label, subtitle/record_subtitle, replicon, accession, length/coordinates.'
        ),
        type=_parse_definition_line_style_arg,
        action='append',
        default=[],
        metavar='LINE:KEY=VALUE',
    )
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
        help='Optional top record-label line (for example organism/strain; repeatable; order matches input records)',
        type=str,
        action='append',
        default=[])
    parser.add_argument(
        '--record_subtitle',
        dest='record_subtitle',
        help='Optional second record-label line (repeatable; order matches input records)',
        type=str,
        action='append',
        default=[])
    parser.add_argument(
        '--show_replicon',
        help='Show inferred replicon labels in linear record-label blocks (default: False).',
        action='store_true')
    parser.add_argument(
        '--hide_accession',
        help='Hide accession labels in linear record-label blocks (default: False).',
        action='store_true')
    parser.add_argument(
        '--hide_length',
        help='Hide length/coordinate labels in linear record-label blocks (default: False).',
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
        help='Linear track layout mode ("above", "middle", or "below"; default: "middle").',
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
        help='Linear custom track slot: <slot_id>:<renderer>@key=value,key=value. Repeat to add slots. The annotations renderer requires set_id from --annotation_table; overlay annotations also require anchor_slot and layer.',
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
            'Effective only with a visible scale, --scale_style ruler, '
            'and --track_layout above|below.'
        ),
        action='store_true',
    )
    _add_format_arg(parser)
    parser.add_argument(
        '-l',
        '--legend',
        help='Legend position (default: "right"; "right", "left", "top", "bottom", "none")',
        type=str,
        default="right")
    parser.add_argument(
            "--show_labels",
            help="Show labels: no argument or 'all' (all records), 'first' (first record only), 'orthogroup_top' (topmost record containing each gbdraw similarity group), 'none' (no labels). Default: 'none'",
            nargs='?',
            const="all",
            default="none",
            choices=["all", "first", "orthogroup_top", "none"],
            type=str
        )
    parser.add_argument(
        '--resolve_overlaps',
        help='Resolve overlaps (experimental; default: False). ',
        action='store_true')
    add_label_args(parser)
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
        '--hide_scale',
        help='Hide the coordinate scale while retaining each record axis (default: False).',
        action='store_true')
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
    _add_legend_size_args(parser)
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
    add_session_args(parser)

    args = parser.parse_args(args)
    validate_input_args(parser, args)
    validate_label_args(parser, args)
    if args.records_table:
        for option_name in (
            "record_label",
            "record_subtitle",
            "record_id",
            "reverse_complement",
            "region",
        ):
            if getattr(args, option_name):
                parser.error(
                    f"--records_table cannot be combined with --{option_name}; "
                    f"use the records table {option_name} column instead."
                )
    if args.records_table and args.multi_record_position:
        parser.error(
            "--records_table cannot be combined with --multi_record_position; "
            "use row and column table columns instead."
        )
    if args.comparisons_table and args.blast:
        parser.error("--comparisons_table cannot be combined with -b/--blast")
    if not math.isfinite(args.linear_record_gap) or args.linear_record_gap < 0:
        parser.error("--linear_record_gap must be a finite non-negative number")
    if args.comparison_height is not None and (
        not math.isfinite(args.comparison_height) or args.comparison_height <= 0
    ):
        parser.error("--comparison_height must be a positive finite number")
    if args.protein_blastp_mode != "none" and args.blast:
        parser.error("--protein_blastp_mode cannot be used with -b/--blast")
    if args.protein_blastp_max_hits <= 0:
        parser.error("--protein_blastp_max_hits must be > 0")
    if args.protein_blastp_output:
        if args.protein_blastp_mode == "none":
            parser.error(
                "--protein_blastp_output requires --protein_blastp_mode"
            )
        if Path(args.protein_blastp_output).suffix.lower() != ".tsv":
            parser.error("--protein_blastp_output must use a .tsv suffix")
    if args.losatp_threads is not None and args.losatp_threads <= 0:
        parser.error("--losatp_threads must be > 0")
    if args.align_orthogroup_feature and args.protein_blastp_mode != "orthogroup" and not args.blast:
        parser.error("--align_orthogroup_feature requires --protein_blastp_mode orthogroup")
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
                show_depth=bool(args.depth_track),
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
    if args.collinear_max_diagonal_drift < 0:
        parser.error("--collinear_max_diagonal_drift must be >= 0")
    if args.collinear_max_conflicts_in_merge_gap < 0:
        parser.error("--collinear_max_conflicts_in_merge_gap must be >= 0")
    if args.collinear_max_paralog_links_per_orthogroup <= 0:
        parser.error("--collinear_max_paralog_links_per_orthogroup must be > 0")
    try:
        ComparisonThresholds(
            evalue=args.evalue,
            bitscore=args.bitscore,
            identity=args.identity,
            alignment_length=args.alignment_length,
        )
    except ValidationError as exc:
        parser.error(str(exc))
    return args




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
    session_request = parse_session_pre_args(cmd_args, mode="linear")
    if session_request is not None:
        with TemporaryDirectory(prefix="gbdraw-session-") as temp_dir:
            session = load_session(session_request.session_path)
            if render_canonical_session_if_present(
                session,
                mode="linear",
                output_override=session_request.output,
                format_override=session_request.format,
                overwrite=session_request.overwrite,
                save_session=session_request.save_session,
                session_output=session_request.session_output,
            ):
                return
            run_spec = session_to_cli_args(
                session,
                mode="linear",
                temp_dir=Path(temp_dir),
                output_override=session_request.output,
                format_override=session_request.format,
            )
            args = _get_args(list(run_spec.args))
            args.overwrite = session_request.overwrite
            args.save_session = session_request.save_session
            args.session_output = session_request.session_output
            args._gbdraw_source_session = session
            run_result = run_linear_from_namespace(args)
            save_session_sidecar_if_requested(
                save_session=session_request.save_session,
                session_output=session_request.session_output,
                output_prefix=args.output,
                run_result=run_result,
                source_session=session,
                cli_invocation_args=run_spec.cli_invocation_args,
                file_bindings=run_spec.file_bindings,
                overwrite=session_request.overwrite,
            )
        return

    args: argparse.Namespace = _get_args(cmd_args)
    run_result = run_linear_from_namespace(args)
    _write_protein_blastp_output(
        args.protein_blastp_output,
        run_result=run_result,
        overwrite=bool(args.overwrite),
    )
    save_session_sidecar_if_requested(
        save_session=bool(args.save_session or args.session_output),
        session_output=args.session_output,
        output_prefix=args.output,
        run_result=run_result,
        cmd_args=cmd_args,
        overwrite=bool(args.overwrite),
    )


def _protein_blastp_output_text(run_result: DiagramRunResult) -> str:
    entries = tuple(
        entry
        for entry in (run_result.losat_cache_entries or ())
        if is_protein_losat_cache_entry(entry)
    )
    manifest = run_result.protein_identity_manifest
    if not entries or manifest is None:
        raise ValidationError(
            "Protein-search output requested, but no raw protein evidence was produced."
        )

    sections = [
        "# gbdraw raw protein-search evidence",
        "# Non-comment lines use the 12 BLAST outfmt 6 columns.",
    ]
    for index, entry in enumerate(entries, start=1):
        filename = str(entry.get("filename") or f"losat_pair_{index}.tsv")
        if "\n" in filename or "\r" in filename:
            raise ValidationError("Protein-search evidence filename contains a newline.")
        hydrated = hydrate_protein_losat_tsv(entry, manifest).rstrip("\r\n")
        sections.extend((f"# entry {index}: {filename}", hydrated))
    return "\n".join(sections) + "\n"


def _write_protein_blastp_output(
    output: str | None,
    *,
    run_result: DiagramRunResult,
    overwrite: bool,
) -> Path | None:
    if not output:
        return None
    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(
        prefix=f".{output_path.name}.",
        dir=output_path.parent,
    ) as temp_name:
        staged_path = Path(temp_name) / output_path.name
        with staged_path.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(_protein_blastp_output_text(run_result))
        commit_staged_output_file(
            staged_path,
            output_path,
            overwrite=overwrite,
        )
    return output_path


def _preflight_protein_blastp_output(
    output: str | None,
    *,
    diagram_output_paths: Sequence[Path],
    session_output_path: Path | None,
    overwrite: bool,
) -> None:
    if not output:
        return
    output_path = Path(output)
    preflight_output_paths((output_path,), overwrite=overwrite)
    try:
        output_identity = output_path.resolve(strict=False)
        reserved_paths = tuple(diagram_output_paths) + (
            ((session_output_path,) if session_output_path is not None else ())
        )
        collides = any(
            path.resolve(strict=False) == output_identity for path in reserved_paths
        )
    except (OSError, ValueError) as exc:
        raise ValidationError(
            f"Could not resolve protein-search output path: {output_path}."
        ) from exc
    if collides:
        raise ValidationError(
            "Protein-search output path collides with a diagram or session output: "
            f"{output_path}."
        )


def run_linear_from_namespace(args: argparse.Namespace) -> DiagramRunResult:
    """Run linear rendering from an already parsed argparse namespace."""

    source_session = getattr(args, "_gbdraw_source_session", None)
    if not isinstance(source_session, Mapping):
        source_session = None
    out_file_prefix: str = args.output
    blast_files: list[str] | None = args.blast
    protein_blastp_mode: str = str(args.protein_blastp_mode or "none")
    losatp_bin: str = args.losatp_bin
    ncbi_blastp_bin: str | None = getattr(args, "ncbi_blastp_bin", None)
    losatp_threads: int | None = args.losatp_threads
    protein_blastp_max_hits: int = args.protein_blastp_max_hits
    protein_blastp_candidate_limit: int | None = args.protein_blastp_candidate_limit
    orthogroup_membership_mode: str = ORTHOGROUP_INFERENCE_VERSION
    orthogroup_member_max_hits: int = 5
    align_orthogroup_feature: str = str(args.align_orthogroup_feature or "").strip()
    collinear_unit_mode: str = str(args.collinear_unit_mode or "auto")
    collinear_anchor_mode: str = "rbh"
    collinear_search_scope: str = str(args.collinear_search_scope or "adjacent")
    collinear_color_mode: str = str(args.collinear_color_mode or "orientation")
    collinearity_params = LosslessCollinearityParameters(
        min_anchors=args.collinear_min_anchors,
        max_unit_gap=args.collinear_max_unit_gap,
        max_diagonal_drift=args.collinear_max_diagonal_drift,
        max_conflicts=args.collinear_max_conflicts_in_merge_gap,
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
    depth_track_groups: list[list[str]] | None = args.depth_track
    depth_track_labels: list[str] | None = list(args.depth_track_label or []) or None
    depth_track_colors: list[str] | None = list(args.depth_track_color or []) or None
    depth_track_heights: list[str] | None = list(args.depth_track_height or []) or None
    depth_track_large_tick_intervals: list[str] | None = list(args.depth_track_large_tick_interval or []) or None
    depth_track_small_tick_intervals: list[str] | None = list(args.depth_track_small_tick_interval or []) or None
    depth_track_tick_font_sizes: list[str] | None = list(args.depth_track_tick_font_size or []) or None
    show_depth: bool = bool(depth_track_groups)
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
    arrow_head_length_ratio: str | float | None = args.arrow_head_length_ratio
    arrow_shaft_width_ratio: float | None = args.arrow_shaft_width_ratio
    feature_height: Optional[float] = args.feature_height
    comparison_height: Optional[float] = args.comparison_height

    out_formats: list[str] = parse_formats(args.format)
    out_formats = handle_output_formats(out_formats)
    user_defined_default_colors: str = args.default_colors
    show_scale: bool = not bool(args.hide_scale)
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
    if blast_files or args.comparisons_table or protein_blastp_mode != "none":
        load_comparison = True
    else:
        load_comparison = False
    palette: str = args.palette
    annotation_options = (
        AnnotationOptions(table_file=args.annotation_table)
        if args.annotation_table
        else None
    )
    config_dict: dict = load_config_toml('gbdraw.data', 'config.toml')

    filtering_cfg = config_dict.setdefault("labels", {}).setdefault("filtering", {})

    block_stroke_color: Optional[str] = args.block_stroke_color
    block_stroke_width: Optional[float] = args.block_stroke_width
    definition_font_size: Optional[float] = args.definition_font_size
    definition_line_styles = parse_definition_line_style_overrides(args.definition_line_style)
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
    if ruler_on_axis and not show_scale:
        logger.warning(
            "WARNING: --ruler_on_axis is ignored when --hide_scale is set."
        )
        ruler_on_axis = False
    elif ruler_on_axis and not (scale_style == "ruler" and track_layout in {"above", "below"}):
        logger.warning(
            "WARNING: --ruler_on_axis is ignored unless --scale_style ruler and --track_layout above|below are set."
        )
        ruler_on_axis = False
    axis_stroke_color: str = (
        args.axis_stroke_color
        if args.axis_stroke_color is not None
        else (
            LINEAR_MODE_PROFILE.linear_ruler_axis_color or "dimgray"
            if ruler_on_axis
            else LINEAR_MODE_PROFILE.linear_axis_color or "lightgray"
        )
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
    filtering_override = dict(filtering_cfg)
    if label_blacklist is not None:
        filtering_override["blacklist_keywords"] = (
            [keyword.strip() for keyword in label_blacklist.split(",") if keyword.strip()]
            if isinstance(label_blacklist, str)
            else label_blacklist
        )
    gc_large_tick_interval = (
        gc_content_large_tick_interval
        if gc_content_large_tick_interval is not None
        else gc_content_tick_interval
    )
    override_candidates: dict[str, object | None] = {
        "objects.features.block_stroke_color": block_stroke_color,
        "objects.features.arrow_geometry.head_length_ratio": arrow_head_length_ratio,
        "objects.features.arrow_geometry.shaft_width_ratio": arrow_shaft_width_ratio,
        "objects.axis.linear.stroke_color": axis_stroke_color,
        "objects.definition.linear.show_replicon": definition_show_replicon,
        "objects.definition.linear.show_accession": definition_show_accession,
        "objects.definition.linear.show_length": definition_show_length,
        "labels.spacing.linear": args.linear_label_spacing,
        "labels.rendering": label_rendering,
        "labels.linear.placement": label_placement,
        "labels.linear.rotation": label_rotation,
        "objects.features.line_stroke_color": line_stroke_color,
        "canvas.show_gc": show_gc,
        "objects.gc_content.mode": gc_content_mode,
        "objects.gc_content.min_percent": gc_content_min_percent,
        "objects.gc_content.max_percent": gc_content_max_percent,
        "objects.gc_content.show_axis": gc_content_show_axis,
        "objects.gc_content.show_ticks": gc_content_show_ticks,
        "objects.gc_content.large_tick_interval": gc_large_tick_interval,
        "objects.gc_content.small_tick_interval": gc_content_small_tick_interval,
        "objects.gc_content.tick_font_size": gc_content_tick_font_size,
        "canvas.show_skew": show_skew,
        "canvas.show_depth": show_depth,
        "objects.depth.fill_color": depth_color,
        "canvas.linear.depth_height": depth_height,
        "objects.depth.min_depth": depth_min,
        "objects.depth.max_depth": depth_max,
        "objects.depth.normalize": depth_normalize,
        "objects.depth.show_axis": depth_show_axis,
        "objects.depth.show_ticks": depth_show_ticks,
        "objects.depth.large_tick_interval": depth_large_tick_interval,
        "objects.depth.small_tick_interval": depth_small_tick_interval,
        "objects.depth.tick_font_size": depth_tick_font_size,
        "objects.depth.share_axis": depth_share_axis,
        "labels.linear.scope": show_labels,
        "canvas.resolve_overlaps": resolve_overlaps,
        "canvas.linear.track_layout": track_layout,
        "canvas.linear.track_axis_gap": track_axis_gap,
        "canvas.linear.ruler_on_axis": ruler_on_axis,
        "canvas.linear.align_center": align_center,
        "canvas.linear.keep_definition_left_aligned": keep_definition_left_aligned,
        "canvas.strandedness": strandedness,
        "labels.filtering.raw": filtering_override,
        "canvas.linear.comparison_height": comparison_height,
        "canvas.linear.default_gc_height": gc_height,
        "objects.scale.show": show_scale,
        "objects.scale.style": scale_style,
        "objects.scale.stroke_color": scale_stroke_color,
        "objects.scale.label_color": scale_label_color,
        "objects.scale.stroke_width": scale_stroke_width,
        "objects.scale.interval": scale_interval,
        "objects.blast_match.style": pairwise_match_style,
        "canvas.linear.normalize_length": normalize_length,
    }
    for line_kind, properties in (definition_line_styles or {}).items():
        for property_name, value in properties.items():
            override_candidates[
                "objects.definition.linear.line_styles."
                f"{line_kind}.{property_name}"
            ] = value
    for width_path in (
        "objects.features.block_stroke_width.short",
        "objects.features.block_stroke_width.long",
    ):
        override_candidates[width_path] = block_stroke_width
    for width_path in (
        "objects.axis.linear.stroke_width.short",
        "objects.axis.linear.stroke_width.long",
    ):
        override_candidates[width_path] = axis_stroke_width
    for width_path in (
        "objects.features.line_stroke_width.short",
        "objects.features.line_stroke_width.long",
    ):
        override_candidates[width_path] = line_stroke_width
    for font_path in (
        "objects.definition.linear.font_size.short",
        "objects.definition.linear.font_size.long",
    ):
        override_candidates[font_path] = definition_font_size
    for font_path in (
        "labels.font_size.linear.short",
        "labels.font_size.linear.long",
    ):
        override_candidates[font_path] = label_font_size
    for height_path in (
        "canvas.linear.default_cds_height.short",
        "canvas.linear.default_cds_height.long",
    ):
        override_candidates[height_path] = feature_height
    for font_path in (
        "objects.scale.font_size.short",
        "objects.scale.font_size.long",
    ):
        override_candidates[font_path] = scale_font_size
    for font_path in (
        "objects.scale.ruler_label_font_size.short",
        "objects.scale.ruler_label_font_size.long",
    ):
        override_candidates[font_path] = effective_ruler_label_font_size
    for size_path in (
        "objects.legends.color_rect_size.short",
        "objects.legends.color_rect_size.long",
    ):
        override_candidates[size_path] = legend_box_size
    for size_path in (
        "objects.legends.font_size.short",
        "objects.legends.font_size.long",
    ):
        override_candidates[size_path] = legend_font_size
    config_dict = modify_config_dict(
        config_dict,
        {
            path: value
            for path, value in override_candidates.items()
            if value is not None
        },
    )

    if args.records_table:
        record_manifest = record_input_manifest_from_table(args.records_table)
        linear_positions: list[str] = []
    else:
        source_count = len(args.gbk or args.gff or ())
        cardinality = _linear_cli_record_cardinality(
            is_gff_source=not bool(args.gbk),
            source_count=source_count,
            load_comparison=load_comparison,
        )
        record_manifest = record_input_manifest_from_paths(
            gbk_paths=args.gbk,
            gff_paths=args.gff,
            fasta_paths=args.fasta,
            cardinalities=cardinality,
            selectors=args.record_id,
            reverse_flags=args.reverse_complement,
            labels=args.record_label,
            subtitles=args.record_subtitle,
            regions=args.region,
        )
        linear_positions = list(args.multi_record_position or [])
    if record_manifest.record_options.regions and blast_files:
        logger.warning(
            "WARNING: Region cropping is enabled; ensure BLAST coordinates "
            "match the cropped regions (and reverse complements if specified)."
        )
    has_table_placement = any(
        record.presentation.grid_row is not None
        for record in record_manifest.records
    )
    linear_layout = (
        LinearMultiRecordOptions(
            record_gap_px=float(args.linear_record_gap),
            multi_record_positions=(
                tuple(linear_positions) if linear_positions else None
            ),
        )
        if linear_positions or has_table_placement
        else None
    )
    depth_tracks = depth_track_inputs_from_cli(
        depth_track_groups,
        labels=depth_track_labels,
        colors=depth_track_colors,
        heights=depth_track_heights,
        large_tick_intervals=depth_track_large_tick_intervals,
        small_tick_intervals=depth_track_small_tick_intervals,
        tick_font_sizes=depth_track_tick_font_sizes,
    )

    canonical_config = copy.deepcopy(config_dict)
    request_path = Path(str(out_file_prefix))
    canonical_request = LinearDiagramRequest(
        records=record_manifest.records,
        options=LinearDiagramOptions(
            config=canonical_config,
            colors=ColorOptions(
                color_table_file=color_table_path or None,
                default_colors_palette=palette,
                default_colors_file=user_defined_default_colors or None,
            ),
            tracks=LinearRequestTrackOptions(
                linear_track_slots=linear_track_slot_specs,
                linear_track_axis_index=linear_track_axis_index,
            ),
            annotations=annotation_options,
            output=LinearOutputOptions(
                legend=legend,
                plot_title_position=plot_title_position,
            ),
            selected_features_set=tuple(selected_features_set),
            feature_visibility_table_file=feature_table_path or None,
            label_whitelist_file=label_whitelist or None,
            qualifier_priority_file=qualifier_priority_path or None,
            label_override_file=label_table_path or None,
            feature_shapes=feature_shapes or None,
            dinucleotide=dinucleotide,
            window=manual_window,
            step=manual_step,
            depth_window=depth_window,
            depth_step=depth_step,
            depth_tracks=depth_tracks,
            plot_title=plot_title or None,
            plot_title_font_size=plot_title_font_size,
            blast_files=tuple(blast_files) if blast_files else None,
            comparison_table_file=args.comparisons_table,
            protein_blastp_mode=protein_blastp_mode,
            pairwise_match_style=pairwise_match_style,
            collinearity_params=collinearity_params,
            collinearity_unit_mode=collinear_unit_mode,
            collinearity_anchor_mode=collinear_anchor_mode,
            collinearity_search_scope=collinear_search_scope,
            collinearity_color_mode=collinear_color_mode,
            losatp_bin=losatp_bin,
            ncbi_blastp_bin=ncbi_blastp_bin,
            losatp_threads=losatp_threads,
            protein_blastp_max_hits=protein_blastp_max_hits,
            protein_blastp_candidate_limit=protein_blastp_candidate_limit,
            orthogroup_membership_mode=orthogroup_membership_mode,
            orthogroup_member_max_hits=orthogroup_member_max_hits,
            collinear_max_paralog_links_per_orthogroup=args.collinear_max_paralog_links_per_orthogroup,
            align_orthogroup_feature=align_orthogroup_feature or None,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
        ),
        layout=linear_layout,
        record_options=record_manifest.record_options,
        output=RenderOutputRequest(
            output_prefix=request_path.name,
            output_directory=(
                request_path.parent if request_path.parent != Path(".") else None
            ),
            formats=tuple(out_formats),
            overwrite=args.overwrite,
        ),
    )
    diagram_output_paths = diagram_request_output_paths(canonical_request)
    session_output_path = preflight_session_sidecar_if_requested(
        save_session=bool(args.save_session or args.session_output),
        session_output=args.session_output,
        output_prefix=args.output,
        diagram_output_paths=diagram_output_paths,
        overwrite=bool(args.overwrite),
    )
    _preflight_protein_blastp_output(
        args.protein_blastp_output,
        diagram_output_paths=diagram_output_paths,
        session_output_path=session_output_path,
        overwrite=bool(args.overwrite),
    )
    legacy_protein_raw_candidates = None
    legacy_protein_derived_evidence = None
    include_feature_catalog = bool(args.save_session or args.session_output)
    if source_session is not None:
        render_result = render_session_compatible_request(
            canonical_request,
            source_session,
            include_feature_catalog=include_feature_catalog,
        )
        legacy_protein_raw_candidates = (
            render_result.legacy_protein_raw_candidates
        )
        legacy_protein_derived_evidence = (
            render_result.legacy_protein_derived_evidence
        )
    else:
        render_result = render_request(
            canonical_request,
            artifacts=CurrentRequestArtifacts(),
            include_feature_catalog=include_feature_catalog,
        )
    canvas = render_result.drawing
    interactive_context = render_result.interactive_context
    rendered_svg = make_rendered_svg(out_file_prefix, request_path.name)
    track_slot_geometry_records = collect_track_slot_geometry_records(
        canvas,
        result_index=0,
        result_name=rendered_svg.svg_path.name,
    )

    return DiagramRunResult(
        mode="linear",
        render_formats=tuple(out_formats),
        outputs=(rendered_svg,),
        feature_metadata=tuple(interactive_context.features) if interactive_context else (),
        biological_feature_metadata=(
            tuple(interactive_context.biological_features)
            if interactive_context
            else ()
        ),
        interactive_contexts=(interactive_context,),
        orthogroup_metadata=(
            tuple(interactive_context.orthogroups)
            if interactive_context is not None
            else None
        ),
        losat_cache_entries=render_result.losat_cache_entries,
        losat_derived_cache_entries=render_result.losat_derived_cache_entries,
        protein_identity_manifest=render_result.protein_identity_manifest,
        legacy_protein_raw_candidates=legacy_protein_raw_candidates,
        legacy_protein_derived_evidence=legacy_protein_derived_evidence,
        run_metadata=build_track_slot_geometry_run_metadata(
            mode="linear",
            records=track_slot_geometry_records,
        ),
        canonical_request=render_result.request,
    )


if __name__ == "__main__":
    # This gets all arguments passed to the script, excluding the script name
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(handler)
    main_args = sys.argv[1:]
    if not main_args:
        main_args.append('--help')
    linear_main(main_args)
