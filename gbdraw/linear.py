#!/usr/bin/env python
# coding: utf-8
"""The gbdraw linear command: its arguments, Session replay and the Linear request."""


import argparse
import copy
import logging
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal, Mapping, NamedTuple, Sequence, TypedDict, cast
from .config.toml import load_config_toml
from .api.request_render import (
    CurrentRequestArtifacts,
    RequestRenderResult,
    diagram_request_output_paths,
    render_request,
    stored_protein_result_request,
)
from .api.session_compat import render_session_compatible_request
from .api.record_planning import (
    RecordInputManifest,
    parse_cli_bool,
    record_input_manifest_from_paths,
    record_input_manifest_from_table,
)
from .api.options import (
    DepthTrackInput,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    LinearOutputOptions,
    LinearRequestTrackOptions,
    LinearWrapDefaults,
    LinearWrapTrackOrder,
    LINEAR_DEFAULT_PLOT_TITLE_POSITION,
    LosatRuntimeOptions,
    LosatSearchOptions,
)
from .layout.similarity_alignment import (
    SimilarityAlignmentReference,
    SimilarityAlignmentReferenceError,
)
from .layout.record_placement import (
    parse_record_row_position,
)
from .api.requests import (
    LinearDiagramRequest,
    RecordCardinality,
)
from .definition_line_styles import (
    parse_definition_line_style_assignment,
    parse_definition_line_style_overrides,
)
from .analysis.collinearity import (
    LosslessCollinearityParameters,
    COLLINEARITY_ANCHOR_MODES,
    COLLINEARITY_MERGE_ORIENTATIONS,
    normalize_collinearity_anchor_mode,
    normalize_collinearity_color_mode,
    normalize_collinearity_search_scope,
)
from .analysis.protein_colinearity import (
    hydrate_protein_losat_tsv,
    is_protein_losat_cache_entry,
)
from .config.modify import modify_config_dict
from .config.models.objects import normalize_pairwise_match_style
from .config.models.canvas import LINEAR_WRAP_LENGTH_PATH, linear_axis_ruler_applies
from .features.shapes import FeatureShape, parse_feature_shape_overrides
from .exceptions import ValidationError
from .mode_profiles import (
    LINEAR_MODE_PROFILE,
    LINEAR_WRAP_DEFAULT_TRACK_LAYOUT,
    effective_linear_track_leaves,
    linear_wrap_default_overrides,
    resolve_mode_profile_overrides,
)
from .tracks import (
    LinearTrackSlot,
    linear_track_slots_from_order,
    normalize_linear_track_slots_with_axis,
    parse_linear_track_slots,
)


from .cli_utils.common import (
    add_feature_density_args,
    depth_track_inputs_from_args,
    _add_depth_track_source_args,
    add_identity_scale_arg,
    _add_arrow_geometry_args,
    _add_block_stroke_args,
    _add_color_table_arg,
    _add_default_colors_arg,
    _add_depth_axis_args,
    _add_depth_color_arg,
    _add_depth_track_tick_args,
    _add_feature_shape_arg,
    _add_figure_size_args,
    _add_format_arg,
    _add_overwrite_arg,
    _add_gc_skew_toggle_args,
    _add_retired_option_args,
    _add_gc_content_axis_args,
    _add_font_args,
    _add_label_rendering_args,
    _add_legend_size_args,
    _add_palette_arg,
    add_analysis_args,
    add_feature_args,
    add_input_args,
    apply_record_display_cli_options,
    add_label_args,
    require_positive_args,
    setup_logging,
    validate_comparison_threshold_args,
    validate_depth_axis_args,
    validate_gc_content_axis_args,
    validate_input_args,
    validate_label_args,
    validate_output_choice_args,
    validate_scale_interval_arg,
    figure_size_overrides,
    legend_row_overrides,
    png_dpi,
    register_font_files,
)
from .cli_utils.diagram_run import (
    cli_output_formats,
    depth_overrides,
    gc_content_overrides,
    label_filtering_override,
    label_style_overrides,
    legend_size_overrides,
    log_render_notices,
    mode_option_fields,
    paired_overrides,
    render_output_request,
    stroke_width_overrides,
    style_overrides,
)
from .cli_utils.session import (
    DiagramRunResult,
    SessionCliRequest,
    add_session_args,
    legacy_replay_png_dpi,
    make_rendered_svg,
    parse_session_pre_args,
    preflight_session_sidecar_if_requested,
    render_canonical_session_if_present,
    save_session_sidecar_if_requested,
    read_stored_font_size_options,
)
from .api.run_metadata import RenderedResult, build_run_metadata
from .render.output_paths import preflight_output_paths
from .session import (
    SessionDocument,
    load_session_document,
    noting_session_png_dpi,
    require_search_frame_protein_results,
    session_drawing_artifacts,
)
from .session_drawings import SessionDrawingArtifacts
from .session_io import session_to_cli_args
from .cli_utils.losat_output import (
    add_losat_runtime_args,
    add_losatn_task_arg,
    parse_positive_int as _parse_positive_int,
    write_losat_output_files,
)
from .comparisons.identity_scale import LINEAR_IDENTITY_SCALE_PATH
from .comparisons.losat_jobs import losat_entry_filename, losat_pair_filename, unique_losat_filenames
from .comparisons.losat_runtime import COMPARISON_LOSAT_PROGRAMS, LOSAT_GENETIC_CODE_TEXT, NUCLEOTIDE_LOSAT_PROGRAMS
from .core.numeric import (
    argparse_float,
    argparse_integer,
    argparse_length_bp,
    argparse_scale_interval,
    parse_ascii_finite_float,
    parse_ascii_integer,
)


# ``--losat`` choices (design 3.2).
_LOSATP_MODES = ("similarity_groups", "collinear", "pairwise")
# ``--losat_output_dir`` file for LOSATP raw evidence (design 3.7, D13).
LOSATP_RAW_OUTPUT_NAME = "losatp.raw.tsv"
# ``--losat_output_dir`` manifest of LOSATN / TLOSATX edges; it is a valid
# ``--comparisons_table`` (design 3.7).
LOSAT_COMPARISONS_OUTPUT_NAME = "comparisons.tsv"
# CLI spelling of the typed LOSAT fields named in option/program diagnostics.
_CLI_LOSAT_FIELD_FLAGS = {
    "losatp_mode": "--losatp_mode",
    "losatn_task": "--losatn_task",
    "record_gencodes": "--losat_gencode (or the records table losat_gencode column)",
    "losatp_max_hits": "--losatp_max_hits",
    "losatp_max_target_seqs": "--losatp_max_target_seqs",
    "losatp_member_max_hits": "--losatp_member_max_hits",
}


def _parse_optional_positive_int(value: str) -> int | None:
    normalized = str(value or "").strip().lower()
    if normalized in {"", "none", "auto", "null"}:
        return None
    parsed = parse_ascii_integer(normalized)
    if parsed is None or parsed <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer or 'none'")
    return parsed


# Setup for the logging system
logger = logging.getLogger()

_LINEAR_OPTION_DEFAULTS = LinearDiagramOptions()
_COLLINEARITY_PARAMETER_DEFAULTS = LosslessCollinearityParameters()


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


def _parse_wrap_length(value: str) -> int | str:
    return "auto" if value.strip().lower() == "auto" else argparse_length_bp(value)


def _parse_linear_track_axis_gap(value: str) -> float | None:
    normalized = str(value).strip().lower()
    if normalized in {"", "auto"}:
        return None
    axis_gap = parse_ascii_finite_float(normalized)
    if axis_gap is None:
        raise argparse.ArgumentTypeError(
            "track axis gap must be a non-negative number or 'auto'"
        )
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


def _parse_collinear_anchor_mode(value: str) -> str:
    try:
        return normalize_collinearity_anchor_mode(value)
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


def _add_linear_input_args(parser: argparse.ArgumentParser) -> None:
    """The record inputs, the linear records table, and the rows of the records."""

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
        type=argparse_float,
        default=24.0,
        help="Gap in pixels between records in the same Linear row (default: 24).",
    )


def _add_comparison_input_args(parser: argparse.ArgumentParser) -> None:
    """The comparisons of the records: a comparisons table or BLAST files, and the orientation they set."""

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
        '--auto_orient',
        action='store_true',
        help=(
            'Reverse-complement a record when its reverse-strand nucleotide hits against '
            'the row above, as that row is drawn, weigh more than half of its hits '
            '(sum of alignment lengths, after the match thresholds). Reads -b/--blast, '
            '--comparisons_table, and losatn/tlosatx results, not protein comparisons; '
            '--reverse_complement 1 keeps its record reversed (default: off).'
        ),
    )


def _add_linear_losat_args(parser: argparse.ArgumentParser) -> None:
    """The LOSAT comparison between the records: the program and LOSATP display, the task, the genetic codes, the runtime, the protein hit caps, the output folder, and the Similarity group alignment."""

    parser.add_argument(
        '--losat',
        choices=COMPARISON_LOSAT_PROGRAMS,
        default=None,
        help=(
            'Run a LOSAT comparison between the records: losatn compares the '
            'nucleotide sequences (LOSAT blastn), tlosatx compares their '
            'translations (LOSAT tblastx), and losatp compares the CDS proteins '
            '(LOSAT blastp). Records in adjacent rows are compared unless '
            '--comparisons_table lists source=losat rows. Cannot be combined with '
            '-b/--blast (default: no LOSAT comparison).'
        ),
    )
    parser.add_argument(
        '--losatp_mode',
        choices=_LOSATP_MODES,
        default=None,
        help=(
            'LOSATP display with --losat losatp: similarity_groups (Similarity groups across all records), collinear '
            '(Collinear blocks), or pairwise (links between records in adjacent rows, or the pairs that the source=losat '
            'rows of --comparisons_table list) (default: pairwise with --comparisons_table, else similarity_groups).'
        ),
    )
    add_losatn_task_arg(parser)
    parser.add_argument(
        '--losat_gencode',
        metavar='CODE',
        type=_parse_positive_int,
        nargs='+',
        default=None,
        help=(
            'TLOSATX genetic code with --losat tlosatx: one for every record or one '
            'per record input; the records table losat_gencode column sets it per '
            f'row. LOSAT reads the NCBI tables {LOSAT_GENETIC_CODE_TEXT} '
            '(default: the runtime default, 1).'
        ),
    )
    add_losat_runtime_args(parser, blast_programs='blastn, tblastx, or blastp')
    parser.add_argument(
        '--losatp_max_hits',
        metavar='N',
        type=_parse_positive_int,
        default=None,
        help='Maximum distinct subject proteins per query protein in --losatp_mode pairwise links (default: 5).',
    )
    parser.add_argument(
        '--losatp_max_target_seqs',
        metavar='N',
        type=_parse_optional_positive_int,
        default=None,
        help=(
            "LOSATP -max_target_seqs per query; 'none' for no cap (default: none). "
            "In a within-record search the protein's hit to itself takes one of these slots."
        ),
    )
    parser.add_argument(
        '--losatp_member_max_hits',
        metavar='N',
        type=_parse_optional_positive_int,
        default=None,
        help=(
            "Maximum distinct subject proteins, not counting the protein itself, kept per "
            "query protein in each directed record-pair search for Similarity group and "
            "Collinear membership; 'none' for no cap (default: none)."
        ),
    )
    parser.add_argument(
        '--losat_output_dir',
        metavar='DIR',
        default=None,
        help=(
            'Write the raw LOSAT evidence to DIR; requires --losat. LOSATN and '
            'TLOSATX: one TSV per compared record pair and comparisons.tsv, which '
            '--comparisons_table accepts unchanged. LOSATP: losatp.raw.tsv with '
            'user-visible protein IDs.'
        ),
    )
    parser.add_argument(
        '--similarity_alignment_feature',
        metavar='ID',
        default="",
        help=(
            'Align linear records using this exact feature SVG hash or protein ID; '
            'Similarity Group IDs are not accepted. Requires --losat losatp with '
            '--losatp_mode similarity_groups.'
        ),
    )


def _add_collinear_args(parser: argparse.ArgumentParser) -> None:
    """The Collinear blocks of LOSATP: units, search scope, Similarity groups, anchors, merge rules, and colors."""

    parser.add_argument(
        '--collinear_unit_mode',
        dest='collinear_unit_mode',
        help=(
            'Units that order Collinear anchors: auto (the CDSs that share a locus ID form '
            'one unit; any other CDS is a unit of its own), cds (one unit per CDS), or locus '
            '(one unit per locus ID; every CDS needs one) (default: auto).'
        ),
        choices=["auto", "cds", "locus"],
        default=_LINEAR_OPTION_DEFAULTS.collinearity_unit_mode)
    parser.add_argument(
        '--collinear_search_scope',
        dest='collinear_search_scope',
        help=(
            'Record pairs whose protein searches supply Collinear similarity-group '
            'evidence: adjacent (records in neighboring rows) or all (every record pair). '
            'It changes the group evidence only, never the blocks, which always connect '
            'adjacent rows, so it has an effect only with --collinear_infer_orthogroups on '
            '(default: adjacent).'
        ),
        type=_parse_collinear_search_scope,
        choices=["adjacent", "all"],
        default=_LINEAR_OPTION_DEFAULTS.collinearity_search_scope)
    parser.add_argument(
        '--collinear_infer_orthogroups',
        choices=["on", "off"],
        default="on" if _LINEAR_OPTION_DEFAULTS.collinear_infer_orthogroups else "off",
        help=(
            'Add within-record protein searches and group proteins, paralogs included, '
            'into Similarity groups before Collinear blocks are built '
            f'(default: {"on" if _LINEAR_OPTION_DEFAULTS.collinear_infer_orthogroups else "off"}).'
        ),
    )
    parser.add_argument(
        '--collinear_anchor_mode',
        dest='collinear_anchor_mode',
        help=(
            'Between-record hits that become Collinear anchors: rbh (reciprocal best hits; '
            'with --collinear_infer_orthogroups on, the links of the Similarity groups) or '
            f'all (every hit that passes the filters) (default: {_LINEAR_OPTION_DEFAULTS.collinearity_anchor_mode}).'
        ),
        type=_parse_collinear_anchor_mode,
        choices=list(COLLINEARITY_ANCHOR_MODES),
        default=_LINEAR_OPTION_DEFAULTS.collinearity_anchor_mode)
    parser.add_argument(
        '--collinear_merge_orientation',
        dest='collinear_merge_orientation',
        help=(
            'How neighboring anchors of a Collinear block must agree: strand (the anchor '
            'strands follow the block orientation and the gene order) or order (the gene '
            f'order only) (default: {_COLLINEARITY_PARAMETER_DEFAULTS.merge_orientation}).'
        ),
        choices=list(COLLINEARITY_MERGE_ORIENTATIONS),
        default=_COLLINEARITY_PARAMETER_DEFAULTS.merge_orientation)
    parser.add_argument(
        '--collinear_min_anchors',
        dest='collinear_min_anchors',
        help='Minimum anchors/genes required for a rendered Collinear block; 1 allows singleton links (default: 1).',
        type=argparse_integer,
        default=_COLLINEARITY_PARAMETER_DEFAULTS.min_anchors)
    parser.add_argument(
        '--collinear_max_unit_gap',
        dest='collinear_max_unit_gap',
        help=(
            'Maximum number of units (see --collinear_unit_mode) between two neighboring '
            'anchors of a Collinear block, on both records; 0 joins adjacent units only. '
            'With --collinear_max_diagonal_drift 0, a gap is bridged only when both records '
            'skip the same number of units (default: 0).'
        ),
        type=argparse_integer,
        default=_COLLINEARITY_PARAMETER_DEFAULTS.max_unit_gap)
    parser.add_argument(
        '--collinear_max_diagonal_drift',
        dest='collinear_max_diagonal_drift',
        help=(
            'Maximum difference between the unit gaps on the two records for two '
            'neighboring anchors of a Collinear block. It needs --collinear_max_unit_gap '
            '1 or more; a value above that gap acts as the gap (default: 0).'
        ),
        type=argparse_integer,
        default=_COLLINEARITY_PARAMETER_DEFAULTS.max_diagonal_drift)
    parser.add_argument(
        '--collinear_max_conflicts_in_merge_gap',
        dest='collinear_max_conflicts_in_merge_gap',
        help=(
            'Maximum other anchors that may lie strictly inside the gap between two '
            'neighboring anchors of a Collinear block, on both records; those anchors '
            'stay in blocks of their own. It needs --collinear_max_unit_gap 1 or more '
            '(default: 1).'
        ),
        type=argparse_integer,
        default=_COLLINEARITY_PARAMETER_DEFAULTS.max_conflicts)
    parser.add_argument(
        '--collinear_max_paralog_links_per_orthogroup',
        dest='collinear_max_paralog_links_per_orthogroup',
        help=(
            'Maximum paralog links kept for each Similarity group, in --losatp_mode '
            'similarity_groups and with --collinear_infer_orthogroups on; it does not change '
            'Collinear blocks (default: 2).'
        ),
        type=argparse_integer,
        default=_LINEAR_OPTION_DEFAULTS.collinear_max_paralog_links_per_orthogroup)
    parser.add_argument(
        '--collinear_color_mode',
        dest='collinear_color_mode',
        help='Collinear ribbon color mode: average_identity, orientation, or orientation_identity (default: orientation).',
        type=_parse_collinear_color_mode,
        choices=["average_identity", "orientation", "orientation_identity"],
        default=_LINEAR_OPTION_DEFAULTS.collinearity_color_mode)


def _add_linear_depth_args(parser: argparse.ArgumentParser) -> None:
    """The depth tracks: sources, heights, ticks, feature density tracks, color, and the depth axis."""

    _add_depth_track_source_args(parser, mode="linear")
    parser.add_argument(
        '--depth_track_height',
        metavar='PX',
        help='Linear depth track height(s) in px. Provide one value or one per --depth_track.',
        type=str,
        nargs='+')
    _add_depth_track_tick_args(parser)
    add_feature_density_args(parser)
    _add_depth_color_arg(parser)
    parser.add_argument(
        '--depth_height',
        help='Depth track height for linear mode (in px; must be > 0).',
        type=argparse_float)
    _add_depth_axis_args(parser)


def _add_record_alignment_args(parser: argparse.ArgumentParser) -> None:
    """How the records line up: centered, at a feature, and where their definitions sit."""

    parser.add_argument(
        '--align_center',
        help='Align genomes to the center (default: False). ',
        action='store_true')
    parser.add_argument(
        '--anchor_feature',
        metavar='KEY=VALUE',
        default=None,
        help=(
            'Align the records at the start as drawn (its 5\' base when it reads rightward, '
            'else its 3\' base) of the feature this exact selector names in each record, for example gene=dnaA or locus_tag=b0001 '
            '(the feature_selector form of the placement and override tables). A complete circular '
            'record is rotated so that the feature starts at its left end, replacing its display start; '
            'a linear or cropped record moves sideways to that place (to the first record\'s feature '
            'when none is rotated). A record without the feature keeps its place and start and is '
            'reported. Wins over --align_center for the records it moves (default: off).'
        ),
    )
    parser.add_argument(
        '--keep_definition_left_aligned',
        dest='keep_definition_left_aligned',
        help=(
            'Keep linear record definitions in the left column. In a multi-record row, a '
            'label or subtitle that no record of the row contradicts becomes a row '
            'definition drawn once beside the row, and text that varies within the row '
            'stays above its own record; this flag aligns those row definitions in a '
            'shared left column. By default, definitions are centered within a shared '
            'column width and follow their row offset '
            '(default: False).'
        ),
        action='store_true')


def _add_pairwise_match_args(parser: argparse.ArgumentParser) -> None:
    """The pairwise matches: the thresholds, the identity scale, and the link style."""

    parser.add_argument(
        '--evalue',
        help='evalue threshold (default=1e-2)',
        type=argparse_float,
        default=LINEAR_MODE_PROFILE.comparison.evalue)
    parser.add_argument(
        '--bitscore',
        help='bitscore threshold (default=50)',
        type=argparse_float,
        default=LINEAR_MODE_PROFILE.comparison.bitscore)
    parser.add_argument(
        '--identity',
        help='identity threshold (default=0)',
        type=argparse_float,
        default=LINEAR_MODE_PROFILE.comparison.identity)
    parser.add_argument(
        '--alignment_length',
        help='minimum BLAST alignment length threshold (default=0)',
        type=argparse_integer,
        default=LINEAR_MODE_PROFILE.comparison.alignment_length)
    add_identity_scale_arg(parser, "pairwise matches")
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


def _add_linear_stroke_args(parser: argparse.ArgumentParser) -> None:
    """The axis and line strokes."""

    parser.add_argument(
        '--axis_stroke_color',
        help='Axis stroke color (str; default: auto: "lightgray", or "dimgray" with --ruler_on_axis)',
        type=str,
        default=None)
    parser.add_argument(
        '--axis_stroke_width',
        help='Axis stroke width (optional; float; default: 5 pt for genomes < 100 kb, 2 pt for genomes >= 100 kb)',
        type=argparse_float)
    parser.add_argument(
        '--line_stroke_color',
        help='Line stroke color (optional; str; default: "lightgray")',
        type=str)
    parser.add_argument(
        '--line_stroke_width',
        help='Line stroke width (optional; float; default: 5 pt for genomes < 100 kb, 1 pt for genomes >= 100 kb)',
        type=argparse_float)


def _add_linear_definition_args(parser: argparse.ArgumentParser) -> None:
    """The record definitions and the shared plot title."""

    parser.add_argument(
        '--definition_font_size',
        help='Definition font size (optional; float; default: 24 pt for genomes < 100 kb, 10 pt for genomes >= 100 kb)',
        type=argparse_float)
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
        type=argparse_float)
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
        help=(
            'Show one automatic name per record, using chromosome, plasmid, then '
            'organelle source qualifiers in that order. Explicit subtitles are '
            'independent (default: False).'
        ),
        action='store_true')
    parser.add_argument(
        '--hide_accession',
        help='Hide accession labels in linear record-label blocks (default: False).',
        action='store_true')
    parser.add_argument(
        '--hide_length',
        help='Hide length/coordinate labels in linear record-label blocks (default: False).',
        action='store_true')


def _add_linear_label_style_args(parser: argparse.ArgumentParser) -> None:
    """The feature label size, placement, rendering policy and style, rotation, and spacing."""

    parser.add_argument(
        '--label_font_size',
        help='Label font size (optional; default: 24 pt for genomes < 100 kb, 5 pt for genomes >= 100 kb)',
        type=argparse_float)
    parser.add_argument(
        '--label_placement',
        help='Linear label placement mode ("auto" or "above_feature"; default: "auto"). "above_feature" draws each label over its own feature (below negative-strand features when --separate_strands is used): a horizontal label starts one row off its feature with a leader, as in "auto"; a turned label (--label_rotation) sits on its feature.',
        type=_parse_linear_label_placement,
        metavar="{auto,above_feature}",
    )
    _add_label_rendering_args(parser, mode="linear")
    parser.add_argument(
        '--label_rotation',
        help='Linear label rotation in degrees (optional; float; default: 0). In above_feature mode, rotated labels start from the feature midpoint.',
        type=argparse_float,
    )
    parser.add_argument(
        '--linear_label_spacing',
        help='Linear label-to-label vertical spacing in px (optional; float; must be > 0).',
        type=argparse_float,
    )


def _add_linear_track_layout_args(parser: argparse.ArgumentParser) -> None:
    """The track layout: above, middle or below the axis, row wrapping, custom track slots, and the ruler on the axis."""

    parser.add_argument(
        '--track_layout',
        help='Linear track layout mode ("above", "middle", or "below"; default: "middle", or "above" when --wrap_length cuts a record).',
        type=_parse_linear_track_layout,
        metavar="{above,middle,below}",
    )
    parser.add_argument(
        '--wrap_length',
        help=(
            'Draw each record over several rows of this length: a length in bp with an optional unit '
            '(for example 50000, 50kb, 0.5Mb), or "auto" (at most 50 kb per row, by the longest record). '
            'With no record longer than one row, the option has no effect. Default: no wrapping.'
        ),
        type=_parse_wrap_length,
        metavar="AUTO|LENGTH",
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
        type=argparse_integer,
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


def _add_linear_legend_label_args(parser: argparse.ArgumentParser) -> None:
    """The legend position, which records show labels, overlap resolution, and the label lists and tables."""

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


def _add_linear_track_height_args(parser: argparse.ArgumentParser) -> None:
    """The feature, GC, and comparison heights."""

    parser.add_argument(
        '--feature_height',
        help='Feature vertical width (optional; float; default: 80 (pixels, 96 dpi) for genomes < 100 kb, 20 for genomes >= 100 kb)',
        type=argparse_float)
    parser.add_argument(
        '--gc_height',
        help='GC content/skew vertical width (optional; float; default: 20 (pixels, 96 dpi))',
        type=argparse_float)
    parser.add_argument(
        '--comparison_height',
        help='Comparison block height (optional; float; optional; default: 60 (pixels, 96 dpi))',
        type=argparse_float)


def _add_linear_scale_args(parser: argparse.ArgumentParser) -> None:
    """The length scale: bar or ruler, its strokes, labels, and interval."""

    parser.add_argument(
        '--hide_scale',
        help='Hide the coordinate scale while retaining each record axis (default: False).',
        action='store_true')
    parser.add_argument(
        '--scale_style',
        help='Style for the length scale ("bar", "ruler"; default: "bar", or "ruler" when --wrap_length cuts a record)',
        type=str,
        choices=["bar", "ruler"])
    parser.add_argument(
        '--scale_stroke_color',
        help='Scale bar/ruler stroke color (optional; str; default: "black")',
        type=str)
    parser.add_argument(
        '--scale_stroke_width',
        help='Scale bar/ruler stroke width (optional; float; default: 3 (pt))',
        type=argparse_float)
    parser.add_argument(
        '--scale_font_size',
        help='Scale bar/ruler font size (optional; float; default: 24 (pt) for genomes < 100 kb, 16 for genomes >= 100 kb).',
        type=argparse_float)
    parser.add_argument(
        '--ruler_label_font_size',
        help='Ruler label font size (optional; float). Overrides --scale_font_size when both are set.',
        type=argparse_float)
    parser.add_argument(
        '--ruler_label_color',
        help='Ruler label color (optional; str; default follows axis color when --ruler_on_axis is active, otherwise black).',
        type=str)
    parser.add_argument(
        '--scale_interval',
        help='Manual tick interval for "ruler" scale style (a length in bp with an optional unit, for example 5000 or 5kb; must be > 0). Overrides automatic calculation; optional',
        type=argparse_scale_interval)


def _add_record_selection_args(parser: argparse.ArgumentParser) -> None:
    """Length normalization, and the regions, records, and orientation drawn of each input."""

    parser.add_argument(
        '--normalize_length',
        help='Normalize record length (experimental; default: False). ',
        action='store_true')
    parser.add_argument(
        '--region',
        help=(
            'Crop a region (repeatable). Format: record_id:start-end[:origin][:rc], #index:start-end[:origin][:rc], '
            'or file:record_selector:start-end[:origin][:rc]. Coordinates are 1-based inclusive. A start greater than the end '
            'without :origin is the reversed span; with :origin it is a window across the origin of a circular record, from start through the origin to end. For multiple records without selectors, provide one spec per record in input order (file order, then record order within each file). Repeat --region for one record to draw several windows of it side by side in one row; windows may overlap.'
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


def _linear_parser() -> argparse.ArgumentParser:
    """The ``gbdraw linear`` parser, with its options in the order ``--help`` lists them."""

    parser = argparse.ArgumentParser(
        description='Generate  plot in PNG/PDF/SVG/PS/EPS.')
    _add_linear_input_args(parser)
    _add_comparison_input_args(parser)
    _add_linear_losat_args(parser)
    _add_retired_option_args(parser, mode="linear")
    _add_collinear_args(parser)
    _add_color_table_arg(parser)
    _add_palette_arg(parser)
    _add_default_colors_arg(parser)
    parser.add_argument(
        '-o',
        '--output',
        help='output file prefix (default: out)',
        type=str,
        default="out")
    _add_overwrite_arg(parser)
    add_analysis_args(parser, nt_subject="dinucleotide skew")
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
    _add_linear_depth_args(parser)
    _add_record_alignment_args(parser)
    _add_pairwise_match_args(parser)
    add_feature_args(parser)
    _add_feature_shape_arg(parser)
    _add_arrow_geometry_args(parser)
    _add_block_stroke_args(parser)
    _add_linear_stroke_args(parser)
    _add_linear_definition_args(parser)
    _add_linear_label_style_args(parser)
    _add_linear_track_layout_args(parser)
    _add_format_arg(parser)
    _add_figure_size_args(parser)
    _add_linear_legend_label_args(parser)
    _add_linear_track_height_args(parser)
    _add_linear_scale_args(parser)
    _add_legend_size_args(parser)
    _add_font_args(parser)
    _add_record_selection_args(parser)
    add_session_args(parser)
    return parser


def _validate_linear_args(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
    argv_tokens: Sequence[str],
    *,
    legacy_session: bool,
) -> None:
    """Refuse option values and combinations the linear command cannot draw, in this order, and set the values
    that follow from them (the LOSATP display, the effective track layout and scale style, the track slots).
    ``argv_tokens``: the arguments as typed; ``legacy_session``: the replay of a released legacy Session, which
    reads its stored font sizes and keeps the fallbacks it was drawn with."""

    if legacy_session:
        read_stored_font_size_options(parser, args)
    validate_input_args(parser, args)
    validate_label_args(parser, args)
    if not legacy_session:
        validate_output_choice_args(args)
    _validate_linear_inputs(parser, args)
    _validate_linear_losat(parser, args, argv_tokens)
    require_positive_args(parser, args, "depth_height")
    validate_depth_axis_args(parser, args, track_heights=True)
    validate_gc_content_axis_args(parser, args)
    validate_scale_interval_arg(parser, args, legacy_session=legacy_session)
    _resolve_linear_track_slots(parser, args)
    if args.label_placement == "above_feature" and args.label_rendering != "auto":
        parser.error("--label_rendering embedded_only|external_only cannot be used with --label_placement above_feature")
    _validate_collinear_args(parser, args)
    validate_comparison_threshold_args(parser, args)


def _validate_linear_inputs(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Refuse per-record options beside a records table, and comparison inputs, record gaps, and comparison
    heights that cannot be drawn."""

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
    if args.linear_record_gap < 0:
        parser.error("--linear_record_gap must be a finite non-negative number")
    if args.comparison_height is not None and args.comparison_height <= 0:
        parser.error("--comparison_height must be a positive finite number")


def _validate_linear_losat(
    parser: argparse.ArgumentParser, args: argparse.Namespace, argv_tokens: Sequence[str]
) -> None:
    """Refuse LOSAT options that do not fit the program or each other, and set the default LOSATP display."""

    if args.losat and args.blast:
        parser.error("--losat cannot be used with -b/--blast")
    if args.losatp_mode is not None and not args.losat:
        parser.error("--losatp_mode requires --losat losatp")
    if args.losat in NUCLEOTIDE_LOSAT_PROGRAMS:
        collinear_flags = sorted(
            {
                str(token).split("=", 1)[0]
                for token in argv_tokens
                if str(token).startswith("--collinear_")
            }
        )
        if collinear_flags:
            raise ValidationError(
                f"{collinear_flags[0]} applies to --losat losatp --losatp_mode "
                f"collinear, not --losat {args.losat}.",
                diagnostic={
                    "code": "COMPARISON_INPUT",
                    "reason": "LOSAT_OPTION_PROGRAM",
                    "field": collinear_flags[0].lstrip("-"),
                    "program": args.losat,
                },
            )
    if args.losat == "losatp" and args.losatp_mode is None:
        args.losatp_mode = "pairwise" if args.comparisons_table else "similarity_groups"
    if args.losat_output_dir and not args.losat:
        parser.error("--losat_output_dir requires --losat")
    if args.similarity_alignment_feature and args.anchor_feature:
        parser.error("--anchor_feature cannot be combined with --similarity_alignment_feature")
    if args.similarity_alignment_feature and (
        args.losat != "losatp" or args.losatp_mode != "similarity_groups"
    ):
        parser.error(
            "--similarity_alignment_feature requires --losat losatp with "
            "--losatp_mode similarity_groups"
        )


def _linear_track_slots_from_order(args: argparse.Namespace, track_layout: str) -> list[LinearTrackSlot]:
    """The track slots of ``--linear_track_order`` in one track layout."""

    return linear_track_slots_from_order(
        args.linear_track_order,
        show_depth=bool(args.depth_track),
        depth_track_count=max(1, len(args.depth_track or [])),
        show_gc=bool(args.show_gc),
        show_skew=bool(args.show_skew),
        dinucleotide=str(args.nt or "GC").upper(),
        track_layout=track_layout,
    )


def _resolve_linear_track_slots(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Set the effective track layout and scale style, and the track slots of ``--linear_track_order`` or
    ``--linear_track_slot`` (or, when the record lengths decide the layout, both orders for the planner)."""

    args.effective_track_layout, args.effective_scale_style = effective_linear_track_leaves(
        {
            path: value
            for path, value in (
                ("canvas.linear.track_layout", args.track_layout),
                ("objects.scale.style", args.scale_style),
            )
            if value is not None
        }
    )
    if args.linear_track_order and args.linear_track_slot:
        parser.error("--linear_track_order cannot be combined with --linear_track_slot")
    if args.linear_track_axis_index is not None and not (args.linear_track_order or args.linear_track_slot):
        parser.error("--linear_track_axis_index requires --linear_track_order or --linear_track_slot")
    try:
        linear_track_slot_specs = None
        args.wrap_track_order = None
        if args.linear_track_order and args.wrap_length is not None and args.track_layout is None:
            # The order follows the layout, which the record lengths decide: the planner chooses and checks the slots.
            args.wrap_track_order = LinearWrapTrackOrder(
                _linear_track_slots_from_order(args, LINEAR_WRAP_DEFAULT_TRACK_LAYOUT),
                _linear_track_slots_from_order(args, args.effective_track_layout),
                args.linear_track_axis_index,
            )
        elif args.linear_track_order:
            linear_track_slot_specs = _linear_track_slots_from_order(args, args.effective_track_layout)
        elif args.linear_track_slot:
            linear_track_slot_specs = parse_linear_track_slots(args.linear_track_slot)
        if linear_track_slot_specs is not None and args.linear_track_axis_index is not None:
            normalize_linear_track_slots_with_axis(linear_track_slot_specs, args.linear_track_axis_index)
        args.linear_track_slot_specs = linear_track_slot_specs
    except Exception as exc:
        parser.error(str(exc))


def _validate_collinear_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Refuse Collinear block counts out of range."""

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


def _get_args(args: Sequence[str] | None, *, _legacy_session: bool = False) -> argparse.Namespace:
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
    parser = _linear_parser()
    args_list = [str(token) for token in (args or [])]
    namespace = parser.parse_args(args)
    _validate_linear_args(parser, namespace, args_list, legacy_session=_legacy_session)
    return namespace


def replay_legacy_session(
    document: SessionDocument,
    session_request: SessionCliRequest,
) -> None:
    """Replay a Session 27-30, which has no canonical request, through CLI arguments."""

    with TemporaryDirectory(prefix="gbdraw-session-") as temp_dir:
        # Each caller loaded this document from a file and is its only holder;
        # the replay reads the parsed payload in place.
        session = document._data
        run_spec = session_to_cli_args(
            session,
            mode="linear",
            temp_dir=Path(temp_dir),
            output_override=session_request.output,
            format_override=session_request.format,
        )
        args = _get_args(list(run_spec.args), _legacy_session=True)
        if (
            args.region
            and any(parse_cli_bool(flag) for flag in args.reverse_complement)
            and (args.blast or args.comparisons_table)
        ):
            logger.warning(
                "This Session crops and reverse-complements a record that has comparisons. "
                "Since 0.14.0 the region is cropped before it is reverse-complemented, so the drawn "
                "window can differ from the one the stored comparison was searched on; re-run the "
                "comparison if its matches look misplaced."
            )
        args.overwrite = session_request.overwrite
        args.embed_fonts = session_request.embed_fonts
        png = legacy_replay_png_dpi(document, session_request, args.format)
        args.dpi = png.dpi
        args.save_session = session_request.save_session
        args.session_output = session_request.session_output
        args._gbdraw_source_session = (
            session_drawing_artifacts(document) if document.drawings else None
        )
        with noting_session_png_dpi(png):
            run_result = run_linear_from_namespace(args)
        save_session_sidecar_if_requested(
            save_session=session_request.save_session,
            session_output=session_request.session_output,
            output_prefix=args.output,
            run_result=run_result,
            source_session=session,
            source_resources=document.resources.sources,
            cli_invocation_args=run_spec.cli_invocation_args,
            file_bindings=run_spec.file_bindings,
            overwrite=session_request.overwrite,
        )


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
    setup_logging()
    session_request = parse_session_pre_args(cmd_args, mode="linear")
    if session_request is not None:
        # Before the Session is read, so its font names resolve to these files.
        register_font_files(session_request.font_files)
        document = load_session_document(session_request.session_path)
        if not render_canonical_session_if_present(
            document,
            mode="linear",
            drawing=session_request.drawings[0] if session_request.drawings else None,
            output_override=session_request.output,
            format_override=session_request.format,
            overwrite=session_request.overwrite,
            save_session=session_request.save_session,
            session_output=session_request.session_output,
            embed_fonts=session_request.embed_fonts,
            dpi=session_request.dpi,
        ):
            if session_request.drawings:
                document.drawing(session_request.drawings[0], mode="linear")
            replay_legacy_session(document, session_request)
        return

    args = _get_args(cmd_args)
    run_result = run_linear_from_namespace(args)
    _write_losat_output_dir(
        args.losat_output_dir,
        program=args.losat,
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


def _losatp_raw_output_text(run_result: DiagramRunResult) -> str:
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
        filename = losat_pair_filename(entry, index)
        if "\n" in filename or "\r" in filename:
            raise ValidationError("Protein-search evidence filename contains a newline.")
        hydrated = hydrate_protein_losat_tsv(entry, manifest).rstrip("\r\n")
        sections.extend((f"# entry {index}: {filename}", hydrated))
    return "\n".join(sections) + "\n"


def _nucleotide_output_files(run_result: DiagramRunResult) -> list[tuple[str, str]]:
    """Per-edge raw TSVs and a ``--comparisons_table`` manifest (design 3.7)."""

    # A Linear run result always carries its Linear request.
    request = cast(LinearDiagramRequest, run_result.canonical_request)
    comparisons = [
        comparison
        for comparison in getattr(request.options, "linear_comparisons", None) or ()
        if comparison.search_frame_text is not None
    ]
    names = [
        losat_entry_filename(entry)
        for entry in run_result.losat_cache_entries or ()
        if entry.get("identityKind") == "nucleotide" and entry.get("display") is not False
    ]
    if not comparisons or len(names) != len(comparisons):
        raise ValidationError(
            "LOSAT output requested, but no raw nucleotide evidence was produced.",
            diagnostic={"code": "LOSAT_RUNTIME", "reason": "OUTPUT"},
        )
    files: list[tuple[str, str]] = []
    rows = ["blast\tquery\tsubject"]
    unique_names = unique_losat_filenames(names, reserved=(LOSAT_COMPARISONS_OUTPUT_NAME,))
    for comparison, candidate in zip(comparisons, unique_names):
        files.append((candidate, str(comparison.search_frame_text)))
        rows.append(
            f"{candidate}\t#{comparison.query_record_index + 1}"
            f"\t#{comparison.subject_record_index + 1}"
        )
    files.append((LOSAT_COMPARISONS_OUTPUT_NAME, "\n".join(rows) + "\n"))
    return files


def _losat_output_path(output_dir: str | None, program: str | None) -> Path | None:
    """The fixed file of ``--losat_output_dir`` (edge TSV names depend on records)."""

    if not output_dir:
        return None
    name = LOSATP_RAW_OUTPUT_NAME if program == "losatp" else LOSAT_COMPARISONS_OUTPUT_NAME
    return Path(output_dir) / name


def _write_losat_output_dir(
    output_dir: str | None,
    *,
    program: str | None,
    run_result: DiagramRunResult,
    overwrite: bool,
) -> Path | None:
    output_path = _losat_output_path(output_dir, program)
    if output_path is None:
        return None
    files = (
        [(LOSATP_RAW_OUTPUT_NAME, _losatp_raw_output_text(run_result))]
        if program == "losatp"
        else _nucleotide_output_files(run_result)
    )
    write_losat_output_files(output_path.parent, files, overwrite=overwrite)
    return output_path


def _preflight_losat_output_dir(
    output_dir: str | None,
    *,
    program: str | None,
    diagram_output_paths: Sequence[Path],
    session_output_path: Path | None,
    overwrite: bool,
) -> None:
    output_path = _losat_output_path(output_dir, program)
    if output_path is None:
        return
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


def _cli_losat_search(
    args: argparse.Namespace,
    record_gencodes: Sequence[int | None] = (),
) -> LosatSearchOptions | None:
    """Map the CLI LOSAT options; the typed options validate them."""

    if not args.losat:
        return None
    try:
        return LosatSearchOptions(
            program=args.losat,
            losatp_mode=args.losatp_mode,
            losatn_task=args.losatn_task,
            record_gencodes=tuple(args.losat_gencode or record_gencodes),
            losatp_max_hits=5 if args.losatp_max_hits is None else args.losatp_max_hits,
            losatp_max_target_seqs=args.losatp_max_target_seqs,
            losatp_member_max_hits=args.losatp_member_max_hits,
            runtime=LosatRuntimeOptions(
                losat_executable=args.losat_bin,
                ncbi_blast_executable=args.ncbi_blast_bin,
                threads=args.losat_threads,
            ),
        )
    except ValidationError as exc:
        diagnostic = exc.diagnostic or {}
        field = str(diagnostic.get("field"))
        flag = _CLI_LOSAT_FIELD_FLAGS.get(field)
        if flag is None or diagnostic.get("reason") not in {"LOSAT_OPTION_PROGRAM", "GENETIC_CODE"}:
            raise
        message = (
            f"{flag} does not apply to --losat {args.losat}."
            if diagnostic.get("reason") == "LOSAT_OPTION_PROGRAM"
            else str(exc).replace(field, flag, 1)
        )
        raise ValidationError(message, diagnostic=diagnostic) from exc


class _CollinearityFields(TypedDict):
    """The collinearity fields of the Linear options: each mode as given, else the Linear default."""

    collinearity_params: LosslessCollinearityParameters
    collinearity_unit_mode: str
    collinearity_anchor_mode: str
    collinearity_search_scope: str
    collinearity_color_mode: str
    orthogroup_membership_mode: str
    collinear_infer_orthogroups: bool
    collinear_max_paralog_links_per_orthogroup: int


def _linear_collinearity_fields(args: argparse.Namespace) -> _CollinearityFields:
    """The collinearity fields; ``LosslessCollinearityParameters`` validates the numbers."""

    defaults = _LINEAR_OPTION_DEFAULTS
    unit_mode = str(args.collinear_unit_mode or defaults.collinearity_unit_mode)
    anchor_mode = str(args.collinear_anchor_mode or defaults.collinearity_anchor_mode)
    search_scope = str(args.collinear_search_scope or defaults.collinearity_search_scope)
    color_mode = str(args.collinear_color_mode or defaults.collinearity_color_mode)
    return {
        "collinearity_params": LosslessCollinearityParameters(
            min_anchors=args.collinear_min_anchors,
            max_unit_gap=args.collinear_max_unit_gap,
            max_diagonal_drift=args.collinear_max_diagonal_drift,
            max_conflicts=args.collinear_max_conflicts_in_merge_gap,
            merge_orientation=args.collinear_merge_orientation,
        ),
        "collinearity_unit_mode": unit_mode,
        "collinearity_anchor_mode": anchor_mode,
        "collinearity_search_scope": search_scope,
        "collinearity_color_mode": color_mode,
        "orthogroup_membership_mode": str(defaults.orthogroup_membership_mode),
        "collinear_infer_orthogroups": args.collinear_infer_orthogroups == "on",
        "collinear_max_paralog_links_per_orthogroup": args.collinear_max_paralog_links_per_orthogroup,
    }


def _linear_ruler_on_axis(args: argparse.Namespace) -> bool:
    """``--ruler_on_axis`` when the scale draws a ruler beside the axis; a warning when it is ignored."""

    if not args.ruler_on_axis:
        return False
    show_scale = not bool(args.hide_scale)
    if not show_scale:
        logger.warning(
            "--ruler_on_axis is ignored when --hide_scale is set."
        )
        return False
    if not linear_axis_ruler_applies(
        show_scale=show_scale, scale_style=args.effective_scale_style, track_layout=args.effective_track_layout
    ):
        # With a wrap length the record lengths decide whether it draws: rows always have the axis ruler.
        if args.wrap_length is None:
            logger.warning(
                "--ruler_on_axis is ignored unless --scale_style ruler and --track_layout above|below are set."
            )
        return False
    return True


def _linear_override_candidates(
    args: argparse.Namespace,
    filtering: Mapping[str, object],
    definition_line_styles: dict[str, dict[str, object]],
    ruler_on_axis: bool,
) -> dict[str, object | None]:
    """Every config leaf a Linear option sets, in the order the request records them (None: not given)."""

    ruler_label_font_size = args.ruler_label_font_size
    return {
        **style_overrides(args, mode="linear"),
        "objects.definition.linear.show_replicon": bool(args.show_replicon),
        "objects.definition.linear.show_accession": not bool(args.hide_accession),
        "objects.definition.linear.show_length": not bool(args.hide_length),
        "labels.spacing.linear": args.linear_label_spacing,
        **label_style_overrides(args, mode="linear"),
        LINEAR_IDENTITY_SCALE_PATH: list(args.identity_scale) if args.identity_scale else None,
        "labels.linear.placement": args.label_placement,
        "labels.linear.rotation": args.label_rotation,
        "objects.features.line_stroke_color": args.line_stroke_color,
        "canvas.show_gc": args.show_gc,
        **gc_content_overrides(args, gc_content_mode=args.gc_content_mode),
        "canvas.show_skew": args.show_skew,
        "canvas.show_depth": bool(args.depth_track),
        "objects.depth.fill_color": args.depth_color,
        "canvas.linear.depth_height": args.depth_height,
        **depth_overrides(args),
        "labels.linear.scope": args.show_labels,
        "canvas.resolve_overlaps": args.resolve_overlaps,
        "canvas.feature_overlap_tolerance_bp": args.feature_overlap_tolerance_bp,
        "canvas.linear.track_layout": args.track_layout,
        LINEAR_WRAP_LENGTH_PATH: args.wrap_length,
        **figure_size_overrides(args),
        **legend_row_overrides(args),
        "canvas.linear.track_axis_gap": args.track_axis_gap,
        "canvas.linear.ruler_on_axis": ruler_on_axis,
        "canvas.linear.align_center": args.align_center,
        "canvas.linear.keep_definition_left_aligned": args.keep_definition_left_aligned,
        "canvas.strandedness": args.separate_strands,
        "labels.filtering.raw": label_filtering_override(filtering, args.label_blacklist),
        "canvas.linear.comparison_height": args.comparison_height,
        "canvas.linear.default_gc_height": args.gc_height,
        "objects.scale.show": not bool(args.hide_scale),
        "objects.scale.style": args.scale_style,
        "objects.scale.stroke_color": args.scale_stroke_color,
        "objects.scale.label_color": args.ruler_label_color,
        "objects.scale.stroke_width": args.scale_stroke_width,
        "objects.scale.interval": args.scale_interval,
        "objects.blast_match.style": args.pairwise_match_style,
        "canvas.linear.normalize_length": args.normalize_length,
        **{
            f"objects.definition.linear.line_styles.{line_kind}.{property_name}": value
            for line_kind, properties in (definition_line_styles or {}).items()
            for property_name, value in properties.items()
        },
        **stroke_width_overrides(args, mode="linear"),
        **paired_overrides(
            {
                "objects.definition.linear.font_size": args.definition_font_size,
                "labels.font_size.linear": args.label_font_size,
                "canvas.linear.default_cds_height": args.feature_height,
                "objects.scale.font_size": args.scale_font_size,
                # --ruler_label_font_size falls back on --scale_font_size.
                "objects.scale.ruler_label_font_size": (
                    ruler_label_font_size if ruler_label_font_size is not None else args.scale_font_size
                ),
            }
        ),
        **legend_size_overrides(args),
    }


class _LinearConfig(NamedTuple):
    """The config the request draws with, and the overrides a wrap applies when it cuts a record."""

    config: dict
    wrap_default_overrides: dict[str, object] | None


def _linear_config(args: argparse.Namespace) -> _LinearConfig:
    """The packaged config with every config option given, completed by the Linear profile."""

    config_dict: dict = load_config_toml('gbdraw.data', 'config.toml')
    filtering = config_dict.setdefault("labels", {}).setdefault("filtering", {})
    definition_line_styles = parse_definition_line_style_overrides(args.definition_line_style)
    ruler_on_axis = _linear_ruler_on_axis(args)
    if args.linear_label_spacing is not None and float(args.linear_label_spacing) <= 0:
        raise ValidationError("linear_label_spacing must be > 0")
    override_candidates = _linear_override_candidates(args, filtering, definition_line_styles, ruler_on_axis)
    explicit_overrides = {path: value for path, value in override_candidates.items() if value is not None}
    # The Linear profile supplies the auto axis color, and the override owner
    # derives unset ruler colors from it, as for API and Web requests.
    config_dict = modify_config_dict(
        config_dict,
        resolve_mode_profile_overrides("linear", explicit_overrides),
    )
    return _LinearConfig(copy.deepcopy(config_dict), linear_wrap_default_overrides(explicit_overrides))


class _LinearRecords(NamedTuple):
    """The records, the LOSAT search with their genetic codes, and their placement."""

    manifest: RecordInputManifest
    losat_search: LosatSearchOptions | None
    layout: LinearMultiRecordOptions | None


def _linear_records(args: argparse.Namespace, losat_search: LosatSearchOptions | None) -> _LinearRecords:
    """Read the records (a records table, else the input files) and take the table's genetic codes for LOSAT."""

    if args.records_table:
        record_manifest = record_input_manifest_from_table(args.records_table)
        linear_positions: list[str] = []
    else:
        cardinality = _linear_cli_record_cardinality(
            is_gff_source=bool(args.gff),
            source_count=len(args.gbk or args.fasta or ()),
            load_comparison=bool(args.blast or args.comparisons_table or losat_search is not None),
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
    record_manifest = apply_record_display_cli_options(record_manifest, args)
    if record_manifest.losat_gencodes and args.losat:
        if args.losat_gencode:
            raise ValidationError(
                "Pass --losat_gencode or the records table losat_gencode column, not both.",
                diagnostic={"code": "COMPARISON_INPUT", "field": "record_gencodes"},
            )
        losat_search = _cli_losat_search(args, record_manifest.losat_gencodes)
    if record_manifest.record_options.regions and args.blast:
        logger.warning(
            "Region cropping is enabled; BLAST coordinates must refer "
            "to the cropped regions, counted on the source strand."
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
    return _LinearRecords(record_manifest, losat_search, linear_layout)


class _LinearInputs(NamedTuple):
    """What the command reads and parses from its options before it builds the request."""

    collinearity: _CollinearityFields
    feature_shapes: dict[str, FeatureShape]
    out_formats: list[str]
    config: _LinearConfig
    records: _LinearRecords
    depth_tracks: tuple[DepthTrackInput, ...] | None


def _read_linear_inputs(args: argparse.Namespace) -> _LinearInputs:
    """The LOSAT search, collinearity, feature shapes, formats, config, records and depth tracks, in this order."""

    losat_search = _cli_losat_search(args)
    collinearity = _linear_collinearity_fields(args)
    feature_shapes = parse_feature_shape_overrides(args.feature_shape)
    out_formats = cli_output_formats(args.format)
    if args.alignment_length < 0:
        raise ValidationError("alignment_length must be >= 0")
    config = _linear_config(args)
    records = _linear_records(args, losat_search)
    depth_tracks = depth_track_inputs_from_args(args, heights=args.depth_track_height)
    return _LinearInputs(collinearity, feature_shapes, out_formats, config, records, depth_tracks)


def _linear_request(
    args: argparse.Namespace,
    inputs: _LinearInputs,
    source_session: SessionDrawingArtifacts | None,
) -> LinearDiagramRequest:
    """The typed Linear request; the typed layer validates it."""

    similarity_alignment_feature = str(args.similarity_alignment_feature or "").strip()
    wrap_default_overrides = inputs.config.wrap_default_overrides
    return LinearDiagramRequest(
        records=inputs.records.manifest.records,
        options=LinearDiagramOptions(
            config=inputs.config.config,
            tracks=(
                None
                if args.wrap_track_order is not None
                else LinearRequestTrackOptions(
                    linear_track_slots=args.linear_track_slot_specs,
                    linear_track_axis_index=args.linear_track_axis_index,
                )
            ),
            output=LinearOutputOptions(
                legend=args.legend,
                # argparse restricts --plot_title_position to the Linear choices.
                plot_title_position=cast(
                    'Literal["center", "top", "bottom"]',
                    str(args.plot_title_position or LINEAR_DEFAULT_PLOT_TITLE_POSITION).strip().lower(),
                ),
            ),
            **mode_option_fields(args, feature_shapes=inputs.feature_shapes, depth_tracks=inputs.depth_tracks),
            blast_files=tuple(args.blast) if args.blast else None,
            comparison_table_file=args.comparisons_table,
            auto_orient=args.auto_orient,
            losat_search=inputs.records.losat_search,
            pairwise_match_style=args.pairwise_match_style,
            **inputs.collinearity,
            wrap_defaults=(
                None
                if wrap_default_overrides is None
                else LinearWrapDefaults(wrap_default_overrides, args.wrap_track_order)
            ),
            anchor_feature=args.anchor_feature,
        ),
        layout=inputs.records.layout,
        similarity_alignment=(
            SimilarityAlignmentReference(feature_id=similarity_alignment_feature)
            if similarity_alignment_feature and source_session is None
            else None
        ),
        record_options=inputs.records.manifest.record_options,
        output=render_output_request(
            args.output,
            inputs.out_formats,
            overwrite=args.overwrite,
            embed_fonts=args.embed_fonts,
            dpi=png_dpi(args.dpi, inputs.out_formats),
        ),
    )


def _render_linear(
    args: argparse.Namespace,
    request: LinearDiagramRequest,
    source_session: SessionDrawingArtifacts | None,
    *,
    writes_session: bool,
) -> RequestRenderResult:
    """Refuse an existing Session sidecar or LOSAT output, then render (a Session's drawing with its artifacts)."""

    diagram_output_paths = diagram_request_output_paths(request)
    session_output_path = preflight_session_sidecar_if_requested(
        save_session=writes_session,
        session_output=args.session_output,
        output_prefix=args.output,
        diagram_output_paths=diagram_output_paths,
        overwrite=bool(args.overwrite),
    )
    _preflight_losat_output_dir(
        args.losat_output_dir,
        program=args.losat,
        diagram_output_paths=diagram_output_paths,
        session_output_path=session_output_path,
        overwrite=bool(args.overwrite),
    )
    if source_session is None:
        try:
            # A Linear request renders to one result, never a batch.
            return cast(
                RequestRenderResult,
                render_request(
                    request,
                    artifacts=CurrentRequestArtifacts(),
                    include_feature_catalog=writes_session,
                ),
            )
        except SimilarityAlignmentReferenceError as exc:
            raise ValidationError(f"--similarity_alignment_feature {exc.detail}") from exc
    # A Linear request renders to one session-compatible result, never a batch.
    render_result = cast(
        RequestRenderResult,
        render_session_compatible_request(
            request,
            source_session,
            include_feature_catalog=writes_session,
        ),
    )
    if writes_session:
        require_search_frame_protein_results(render_result.migration_report)
    return render_result


def _linear_run_result(
    args: argparse.Namespace,
    render_result: RequestRenderResult,
    out_formats: list[str],
    *,
    from_session: bool,
    writes_session: bool,
) -> DiagramRunResult:
    """The run's SVG, feature metadata, LOSAT results and the request the Session saves, after the render's log."""

    log_render_notices((render_result,))
    interactive_context = render_result.interactive_context
    rendered_svg = make_rendered_svg(args.output, Path(str(args.output)).name)
    # A Session replay keeps the protein evidence the Session saved.
    migration_report = render_result.migration_report if from_session else None

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
        legacy_protein_raw_candidates=(
            migration_report.protein_raw_candidates if migration_report is not None else None
        ),
        legacy_protein_derived_evidence=(
            migration_report.protein_derived_evidence if migration_report is not None else None
        ),
        run_metadata=build_run_metadata("linear", [RenderedResult(render_result, rendered_svg.result_name)]),
        # A Session stores the protein results this render computed (OV-1068).
        canonical_request=(
            stored_protein_result_request(render_result.request, render_result.linear_metadata)
            if writes_session
            else render_result.request
        ),
    )


def run_linear_from_namespace(args: argparse.Namespace) -> DiagramRunResult:
    """Run linear rendering from an already parsed argparse namespace."""

    source_session = getattr(args, "_gbdraw_source_session", None)
    if not isinstance(source_session, SessionDrawingArtifacts):
        source_session = None
    inputs = _read_linear_inputs(args)
    request = _linear_request(args, inputs, source_session)
    writes_session = bool(args.save_session or args.session_output)
    render_result = _render_linear(args, request, source_session, writes_session=writes_session)
    return _linear_run_result(
        args,
        render_result,
        inputs.out_formats,
        from_session=source_session is not None,
        writes_session=writes_session,
    )


if __name__ == "__main__":
    # ``python -m gbdraw.linear``: linear_main sets up the log as ``gbdraw linear`` does.
    main_args = sys.argv[1:]
    if not main_args:
        main_args.append('--help')
    linear_main(main_args)
