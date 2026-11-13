#!/usr/bin/env python
# coding: utf-8
"""The gbdraw circular command: its arguments, Session replay and the Circular request."""

import sys
import argparse
import copy
import logging
import os
from pathlib import Path
from dataclasses import replace
from tempfile import TemporaryDirectory
from typing import Any, Literal, Mapping, NamedTuple, Sequence, TypedDict, cast
from .io.cli_tables import (
    CircularTrackTable,
    ConservationTableEntry,
    conservation_table_text,
    read_circular_track_table,
    require_listable_comparison_sequences,
)
from .config.models.objects import CONSERVATION_COLOR_SCHEME_NAMES
from .config.toml import load_config_toml
from .api.options import (
    CircularDiagramOptions,
    LosatRuntimeOptions,
    LosatSearchOptions,
    CircularMultiRecordOptions,
    CircularOutputOptions,
    ONE_CIRCLE_GAP_DEGREES_RANGE,
    CircularRequestTrackOptions,
    DepthTrackInput,
    RingGenomeFiles,
    _ring_genome_files,
)
from .api.request_render import (
    CircularBatchRenderResult,
    CircularBatchRequestPlan,
    RequestRenderResult,
    build_request_plan_diagram,
    diagram_request_output_paths,
    plan_request,
    render_prepared_request,
    render_request,
)
from .api.record_planning import (
    RecordInputManifest,
    record_input_manifest_from_paths,
    record_input_manifest_from_table,
)
from .api.requests import (
    CircularBatchOutputPolicy,
    CircularBatchRequest,
    CircularDiagramRequest,
    DiagramRequest,
    RecordCardinality,
)
from .config.modify import modify_config_dict
from .layout.record_placement import parse_record_row_position
from .features.shapes import FeatureShape, parse_feature_shape_overrides
from .exceptions import ValidationError
from .mode_profiles import CIRCULAR_MODE_PROFILE, resolve_mode_profile_overrides
from .tracks import (
    CircularTrackSlot,
    ScalarSpec,
    circular_track_slots_from_order,
    normalize_listed_circular_track_slots,
    parse_circular_track_slots,
)

from .cli_utils.losat_output import (
    add_losat_runtime_args,
    add_losatn_task_arg,
    parse_positive_int,
    write_losat_output_files,
)
from .comparisons.identity_scale import CIRCULAR_IDENTITY_SCALE_PATH
from .comparisons.losat_jobs import unique_losat_filenames
from .comparisons.losat_runtime import COMPARISON_LOSAT_PROGRAMS, LOSAT_GENETIC_CODE_TEXT
from .core.numeric import argparse_float, argparse_integer, argparse_scale_interval
from .render.output_paths import preflight_output_paths
from .cli_utils.common import (
    add_feature_density_args,
    depth_track_inputs_from_args,
    _add_depth_track_source_args,
    _add_arrow_geometry_args,
    _add_block_stroke_args,
    _add_comparison_filter_args,
    _add_depth_axis_args,
    _add_depth_color_arg,
    _add_depth_track_tick_args,
    _add_feature_shape_arg,
    _add_figure_size_args,
    _add_format_arg,
    _add_overwrite_arg,
    _add_retired_option_args,
    _add_gc_skew_toggle_args,
    _add_gc_content_axis_args,
    _add_font_args,
    _add_label_rendering_args,
    _add_legend_size_args,
    add_analysis_args,
    add_color_args,
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
    RenderedSvg,
    add_session_args,
    legacy_replay_png_dpi,
    diagram_request_rendered_svgs,
    parse_session_pre_args,
    preflight_session_sidecar_if_requested,
    render_canonical_session_if_present,
    resolve_session_sidecar_path,
    save_session_sidecar_if_requested,
    read_stored_font_size_options,
)
from .api.run_metadata import RenderedResult, build_run_metadata
from .session import SessionDocument, load_session_document, noting_session_png_dpi
from .session_io import session_to_cli_args

# Setup for the logging system
logger = logging.getLogger()

# ``--losat`` for similarity rings (design 3.3). losatp is accepted so that it
# reports COMPARISON_INPUT / RING_LOSAT_PROGRAM rather than a usage error.
# ``--losat_output_dir`` manifest of the rings; it is a valid --conservation_table.
LOSAT_CONSERVATION_OUTPUT_NAME = "conservation.tsv"
# CLI spelling of typed ring fields named in ring diagnostics (mapping only).
_CLI_RING_FIELD_FLAGS = {
    "conservation_blast_files": "--conservation_blast",
    "conservation_dataframes": "--conservation_blast",
    "conservation_reference": "--conservation_reference",
    "conservation_sequence_files": "--conservation_sequence",
    "conservation_labels": "--conservation_labels",
    "conservation_colors": "--conservation_colors",
    "conservation_losat_gencodes": "--conservation_losat_gencode",
    "record_gencodes": "--losat_gencode",
    "losatn_task": "--losatn_task",
}


def _circular_cli_record_cardinality() -> RecordCardinality:
    """Preserve the CLI contract that every direct input entry is rendered."""

    return RecordCardinality.ALL




def _parse_multi_record_position_arg(value: str) -> str:
    try:
        selector_text, row = parse_record_row_position(value)
    except ValidationError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return f"{selector_text}@{row}"


def conservation_colors_from_cli(values: Sequence[str] | None) -> list[str] | None:
    """``--conservation_colors``: ``auto`` is a ring without a color of its own, which gbdraw reads as ``''``."""

    return ['' if color.strip().lower() == 'auto' else color for color in values or ()] or None


def _add_conservation_input_args(parser: argparse.ArgumentParser) -> None:
    """The similarity ring inputs: BLAST tables, the ring manifest and the comparison genomes."""

    parser.add_argument(
        '--conservation_blast',
        metavar='BLAST',
        help='Precomputed BLAST outfmt 6/7 file(s) for circular similarity rings.',
        type=str,
        nargs='+')
    parser.add_argument(
        '--conservation_table',
        metavar='TSV',
        help=(
            'TSV manifest of the similarity rings: BLAST files, comparison genomes (the files of '
            'one genome separated by ";"), labels, and colors.'
        ),
        type=str)
    parser.add_argument(
        '--conservation_sequence',
        metavar='FILE',
        help=(
            'Comparison genome sequence file(s) (FASTA, GenBank, or DDBJ; one file is one '
            'genome): the ring queries with --losat, otherwise one per --conservation_blast '
            'for interactive span export. For a genome in several files use --conservation_genome.'
        ),
        type=str,
        nargs='+')
    parser.add_argument(
        '--conservation_genome',
        metavar='FILE',
        help=(
            'The files of one comparison genome (FASTA, GenBank, or DDBJ), such as a chromosome '
            'and its plasmids; repeat the option once per ring. The ring reads every file whole, '
            'in the given order, as one genome. Labels, colors, genetic codes, and '
            '--conservation_blast files align by occurrence. Cannot be combined with '
            '--conservation_sequence or --conservation_table.'
        ),
        type=str,
        action='append',
        nargs='+')


def _validate_conservation_input_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Refuse similarity ring inputs that cannot be given together, and empty ``--conservation_genome`` files."""

    if args.conservation_genome and args.conservation_sequence:
        parser.error("--conservation_genome cannot be combined with --conservation_sequence")
    if args.conservation_genome and args.conservation_table:
        parser.error("--conservation_genome cannot be combined with --conservation_table")
    if any("" in genome and genome != [""] for genome in args.conservation_genome or ()):
        parser.error("--conservation_genome files must not be empty; give '' alone for a ring without a genome")
    if args.conservation_table and args.conservation_blast:
        parser.error("--conservation_table cannot be combined with --conservation_blast")
    if args.conservation_table and args.conservation_sequence:
        parser.error("--conservation_table cannot be combined with --conservation_sequence")
    if args.conservation_table and args.conservation_labels:
        parser.error("--conservation_table cannot be combined with --conservation_labels")
    if args.conservation_table and args.conservation_colors:
        parser.error("--conservation_table cannot be combined with --conservation_colors")


def _cli_conservation_genomes(args: argparse.Namespace) -> list[str | list[str] | None] | None:
    """Each ring's comparison genome: one ``--conservation_genome`` occurrence ('' alone: none) or one file."""

    if args.conservation_genome:
        return [genome if genome != [""] else None for genome in args.conservation_genome]
    return list(args.conservation_sequence or []) or None


def _add_record_layout_args(parser: argparse.ArgumentParser) -> None:
    """The options that lay several records out: a grid canvas, or one after another on one circle
    (``CircularMultiRecordOptions``)."""

    parser.add_argument(
        '--multi_record_canvas',
        help='Place multiple records on one shared canvas using automatic grid layout (default: False).',
        action='store_true')
    parser.add_argument(
        '--one_circle',
        help=(
            'Draw all records one after another on one circle, each on its own arc with its own ticks and GC '
            'content (default: False). GC skew is off and GC content is drawn in percent unless set.'
        ),
        action='store_true')
    # The defaults are the API's (``CircularMultiRecordOptions``), which also checks the ranges.
    defaults = CircularMultiRecordOptions()
    low, high = ONE_CIRCLE_GAP_DEGREES_RANGE
    parser.add_argument(
        '--one_circle_gap_degrees',
        help=(
            f'Gap between records on one circle, in degrees ({low:g} to {high:g}; '
            f'default: {defaults.one_circle_gap_degrees}).'
        ),
        type=argparse_float,
        default=defaults.one_circle_gap_degrees)
    parser.add_argument(
        '--one_circle_min_length',
        metavar='BP',
        help=(
            'Records shorter than this are left off the circle and counted in the centre text, also with '
            f'--hide_definition (default: {defaults.one_circle_min_length}).'
        ),
        type=argparse_integer,
        default=defaults.one_circle_min_length)
    parser.add_argument(
        '--one_circle_no_record_band',
        help='Leave out the band that names each record on one circle.',
        action='store_true')


class _OneCircleLayoutFields(TypedDict):
    """The ``CircularMultiRecordOptions`` fields the one-circle options set."""

    record_layout: Literal["grid", "one_circle"]
    one_circle_gap_degrees: float
    one_circle_min_length: int
    one_circle_record_band: bool


def _one_circle_layout_fields(args: argparse.Namespace) -> _OneCircleLayoutFields:
    """The ``CircularMultiRecordOptions`` fields the one-circle options set."""

    if args.one_circle and args.multi_record_canvas:
        raise ValidationError(
            "--one_circle and --multi_record_canvas are two layouts of several records; choose one.",
            diagnostic={"code": "INPUT_INVALID", "reason": "FIELDS"},
        )
    return {
        "record_layout": "one_circle" if args.one_circle else "grid",
        "one_circle_gap_degrees": args.one_circle_gap_degrees,
        "one_circle_min_length": args.one_circle_min_length,
        "one_circle_record_band": not args.one_circle_no_record_band,
    }


def _profile_values(args: argparse.Namespace) -> dict[str, object]:
    """The CLI's values of the mode profile's config paths: as typed, else the profile's, one circle's on one
    circle (``resolve_mode_profile_overrides``, the API's rule). The CLI writes its whole config, so the API's
    profile step never reaches it."""

    typed: dict[str, object] = {}
    if args.show_skew_typed:
        typed["canvas.show_skew"] = args.show_skew
    if args.gc_content_mode is not None:
        typed["objects.gc_content.mode"] = args.gc_content_mode
    return resolve_mode_profile_overrides("circular", typed, one_circle=args.one_circle)


def _shown_skew(args: argparse.Namespace) -> bool:
    """Whether the figure draws GC skew: as typed, else the profile's default, which one circle turns off."""

    return bool(_profile_values(args)["canvas.show_skew"])


def _add_circular_input_output_args(parser: argparse.ArgumentParser) -> None:
    """The record inputs, the circular records table, and the output prefix."""

    add_input_args(parser)
    parser.add_argument(
        "--records_table",
        metavar="TSV",
        help="TSV manifest for row-based input records and circular placement metadata.",
        type=str)
    parser.add_argument(
        '-o',
        '--output',
        help='output file prefix (default: accession number of the sequence)',
        type=str)
    _add_overwrite_arg(parser)


def _add_species_strain_args(parser: argparse.ArgumentParser) -> None:
    """The species and strain of the center definition."""

    parser.add_argument(
        '--species',
        help='Species name (optional; e.g. "<i>Escherichia coli</i>", "<i>Ca.</i> Hepatoplasma crinochetorum")',
        type=str)
    parser.add_argument(
        '--strain',
        help='Strain/isolate name (optional; e.g. "K-12", "Av")',
        type=str)


def _add_circular_stroke_args(parser: argparse.ArgumentParser) -> None:
    """The axis and line strokes."""

    parser.add_argument(
        '--axis_stroke_color',
        help='Axis stroke color (str; default: "gray")',
        type=str)
    parser.add_argument(
        '--axis_stroke_width',
        help='Axis stroke width (optional; float; default: 3 pt for genomes < 75 kb, 1 pt for genomes >= 75 kb)',
        type=argparse_float)
    parser.add_argument(
        '--line_stroke_color',
        help='Line stroke color (str; default: "gray")',
        type=str)
    parser.add_argument(
        '--line_stroke_width',
        help='Line stroke width (optional; float; default: 5 pt for genomes < 75 kb, 1 pt for genomes >= 75 kb)',
        type=argparse_float)


def _add_circular_title_args(parser: argparse.ArgumentParser) -> None:
    """The center definition and the plot title."""

    parser.add_argument(
        '--definition_font_size',
        help='Definition font size (optional; default: 18)',
        type=argparse_float)
    parser.add_argument(
        '--hide_definition',
        help=(
            'Draw no center definition for any record (name, strain, replicon, accession, length, GC); '
            'the inside tracks keep no room for it. A plot title (--plot_title_position) is still drawn.'
        ),
        action='store_true')
    parser.add_argument(
        '--plot_title',
        help='Circular plot title shown when plot_title_position is top or bottom (optional; defaults to species + strain).',
        type=str)
    parser.add_argument(
        '--plot_title_font_size',
        help='Plot title font size for circular top/bottom title layout (optional; default: 32).',
        type=argparse_float)
    parser.add_argument(
        '--keep_full_definition_with_plot_title',
        help='Keep the full species/strain center label when a circular plot title is shown (default: False).',
        action='store_true')


def _add_circular_center_args(parser: argparse.ArgumentParser) -> None:
    """The display start (``--origin_feature``) and the room the center definition keeps inside the tracks."""

    parser.add_argument(
        '--origin_feature',
        metavar='KEY=VALUE',
        default=None,
        help=(
            'Rotate each record so that the feature this exact selector names starts at 12 o\'clock, '
            'for example gene=dnaA or locus_tag=b0001 (the feature_selector form of the placement and '
            'override tables): its start as drawn (its 5\' base when it reads clockwise, else its 3\' '
            'base) becomes the display start. A record without '
            'the feature keeps its start and is reported. Cannot be combined with '
            '--display_start_coordinate or a crop (default: off).'
        ),
    )
    parser.add_argument(
        '--center_reserved_radius',
        help='Override the center-label reserved radius for circular track packing (in px; must be >= 0).',
        type=argparse_float)
    parser.add_argument(
        '--no_track_narrowing',
        dest='track_narrowing',
        help=(
            'Never narrow the inside tracks for the center definition: when moving them outward does not clear it, '
            'keep their widths and draw the definition over them (default: narrow Auto widths down to '
            '--track_narrowing_min_px).'
        ),
        action='store_false',
        default=None)
    parser.add_argument(
        '--track_narrowing_min_px',
        help=(
            'The narrowest width (px; must be >= 0) an inside track with an Auto width takes when it is narrowed for '
            'the center definition (default: 30). A narrower track, or a typed width, is not narrowed.'
        ),
        type=argparse_float)


def _add_conservation_ring_args(parser: argparse.ArgumentParser) -> None:
    """The similarity rings: the reference side of the BLAST tables, the ring labels, colors, width and gap."""

    parser.add_argument(
        '--conservation_reference',
        help='BLAST side containing displayed circular reference coordinates.',
        choices=['query', 'subject', 'auto'],
        default='auto',
        type=str)
    parser.add_argument(
        '--conservation_labels',
        metavar='LABEL',
        help='Labels for similarity rings, aligned by logical source index.',
        type=str,
        nargs='+')
    parser.add_argument(
        '--conservation_colors',
        metavar='COLOR',
        help=(
            'Colors for similarity rings, aligned by logical source index. Accepts SVG color names, #RRGGBB, '
            'or auto (the ring takes the --conservation_color_scheme color of its position).'
        ),
        type=str,
        nargs='+')
    parser.add_argument(
        '--conservation_color_scheme',
        help=(
            'Colors of the similarity rings without a color of their own, by drawn position from the outermost '
            'ring: same (every ring the default color), ten (ten colors, repeating from ring 11), or per_five '
            '(one color per five rings, repeating from ring 31) (default: same).'
        ),
        choices=CONSERVATION_COLOR_SCHEME_NAMES,
        default=None)
    parser.add_argument(
        '--conservation_ring_width',
        help='Similarity ring width for circular mode (in px; must be > 0).',
        type=argparse_float)
    parser.add_argument(
        '--conservation_ring_gap',
        help='Similarity ring gap for circular mode (in px; must be > 0).',
        type=argparse_float)


def _add_ring_losat_args(parser: argparse.ArgumentParser) -> None:
    """The LOSAT search that builds the similarity rings: the program, the task, the genetic codes, the runtime and the output folder."""

    parser.add_argument(
        '--losat',
        choices=COMPARISON_LOSAT_PROGRAMS,
        default=None,
        help=(
            'Build the similarity rings with LOSAT: losatn compares nucleotide '
            'sequences (LOSAT blastn) and tlosatx their translations (LOSAT tblastx). '
            'Each --conservation_sequence or --conservation_genome genome is the query and '
            'every record of the reference file (uncropped, also records not drawn) is the '
            'subject. '
            'Cannot be combined with --conservation_blast (default: no LOSAT search).'
        ),
    )
    add_losatn_task_arg(parser)
    parser.add_argument(
        '--losat_gencode',
        metavar='CODE',
        type=parse_positive_int,
        default=None,
        help=(
            'TLOSATX genetic code of the reference file with --losat tlosatx; LOSAT '
            f'reads the NCBI tables {LOSAT_GENETIC_CODE_TEXT} (default: 1).'
        ),
    )
    parser.add_argument(
        '--conservation_losat_gencode',
        metavar='CODE',
        type=parse_positive_int,
        nargs='+',
        default=None,
        help=(
            'TLOSATX genetic code of each comparison genome with --losat tlosatx: one '
            'for every ring or one per comparison genome, in the same order; the '
            '--conservation_table losat_gencode column sets it per row. LOSAT reads '
            f'the NCBI tables {LOSAT_GENETIC_CODE_TEXT} (default: 1).'
        ),
    )
    add_losat_runtime_args(parser, blast_programs='blastn or tblastx')
    parser.add_argument(
        '--losat_output_dir',
        metavar='DIR',
        default=None,
        help=(
            'Write the ring tables of the LOSAT search to DIR; requires --losat. One TSV per ring, '
            'in the search frame of the drawn records, and conservation.tsv, which '
            '--conservation_table accepts unchanged.'
        ),
    )


def _add_multi_record_canvas_args(parser: argparse.ArgumentParser) -> None:
    """The grid canvas of several records: sizes, gaps, and placement."""

    parser.add_argument(
        '--multi_record_size_mode',
        help='Size mode for multi-record circular canvas ("auto", "linear", "equal"; default: "auto").',
        type=str,
        choices=['auto', 'linear', 'equal'],
        default='auto')
    parser.add_argument(
        '--multi_record_min_radius_ratio',
        help='Minimum radius ratio for multi-record scaling (0 < ratio <= 1; default: 0.55).',
        type=argparse_float,
        default=0.55)
    parser.add_argument(
        '--multi_record_column_gap_ratio',
        help='Additional horizontal gap ratio between visible content bounds in each multi-record row (>= 0; default: 0.10).',
        type=argparse_float,
        default=0.10)
    parser.add_argument(
        '--multi_record_row_gap_ratio',
        help='Additional gap ratio between multi-record row content bounds (>= 0; default: 0.05).',
        type=argparse_float,
        default=0.05)
    parser.add_argument(
        '--multi_record_position',
        help="Record placement for multi-record canvas (repeatable): <selector>@<row> where selector is #index or record_id and row starts at 1.",
        type=_parse_multi_record_position_arg,
        action='append',
        default=[])


def _add_circular_track_args(parser: argparse.ArgumentParser) -> None:
    """The plot title position, the strands, the track preset, and the overlap rule."""

    parser.add_argument(
        '--plot_title_position',
        help='Plot title position in circular mode ("none", "top", "bottom"; default: "none").',
        type=str,
        choices=['none', 'top', 'bottom'],
        default='none')
    parser.add_argument(
        '--separate_strands',
        help='Separate strands (default: False).',
        action='store_true')
    parser.add_argument(
        '--track_type',
        help='Circular track preset. Custom track slots inherit omitted geometry from this preset.',
        type=str,
        choices=['tuckin', 'middle', 'spreadout'],
        default="tuckin")
    parser.add_argument(
        '--resolve_overlaps',
        help='Resolve overlapping features by placing them on separate tracks (default: False). Useful for plasmid visualization.',
        action='store_true')


def _add_circular_label_mode_args(parser: argparse.ArgumentParser) -> None:
    """Which labels are drawn and how: ``--labels``, the rendering policy and label style, and the placement."""

    parser.add_argument(
        '--labels',
        help='Label placement mode: no argument or "out" (outside), "both" (outside+inside), or "none" (hidden). Default: "none".',
        nargs='?',
        const='out',
        default='none',
        choices=['none', 'out', 'both'],
        type=str)
    _add_label_rendering_args(parser, mode="circular")
    parser.add_argument(
        '--label_placement',
        help='External circular label placement: "horizontal" or "radial". Default: "horizontal".',
        choices=['horizontal', 'radial'],
        default='horizontal',
        type=str)


def _add_label_ellipse_args(parser: argparse.ArgumentParser) -> None:
    """The offsets of the outer and inner label ellipses."""

    parser.add_argument(
        '--outer_label_x_radius_offset',
        help='Outer label x-radius offset factor (float; default from config)',
        type=argparse_float)
    parser.add_argument(
        '--outer_label_y_radius_offset',
        help='Outer label y-radius offset factor (float; default from config)',
        type=argparse_float)
    parser.add_argument(
        '--inner_label_x_radius_offset',
        help='Inner label x-radius offset factor (float; default from config)',
        type=argparse_float)
    parser.add_argument(
        '--inner_label_y_radius_offset',
        help='Inner label y-radius offset factor (float; default from config)',
        type=argparse_float)


def _add_circular_scale_args(parser: argparse.ArgumentParser) -> None:
    """The scale ticks and their labels, the label spacing, and the feature track width."""

    parser.add_argument(
        '--scale_interval',
        help='Manual scale interval for circular mode (a length in bp with an optional unit, for example 5000 or 5kb; must be > 0). Overrides automatic calculation.',
        type=argparse_scale_interval)
    parser.add_argument(
        '--hide_scale',
        help='Hide the primary genome-coordinate scale while retaining the circular axis.',
        action='store_true')
    parser.add_argument(
        '--tick_label_font_size',
        help='Tick label font size for circular mode (optional; float; default: 14 (pt)).',
        type=argparse_float)
    parser.add_argument(
        '--circular_label_spacing',
        help='Circular label-to-label vertical spacing in px (optional; float; must be > 0).',
        type=argparse_float)
    parser.add_argument(
        '--feature_width',
        help='Feature track width for circular mode (in px; must be > 0).',
        type=argparse_float)


def _add_circular_track_slot_args(parser: argparse.ArgumentParser) -> None:
    """The custom track slots: an order, slot specs, or a table, and the axis index."""

    parser.add_argument(
        '--circular_track_order',
        help='Comma-separated circular slot order. Omitted built-in slot geometry inherits --track_type.',
        type=str)
    parser.add_argument(
        '--circular_track_slot',
        help='Circular track slot spec: <slot_id>:<renderer>@key=value,key=value. Can be repeated. Omitted built-in slot geometry inherits --track_type; use r, w, inner_gap_px, outer_gap_px, side, and z for explicit overrides. The annotations renderer requires set_id from --annotation_table; overlay annotations also require anchor_slot and layer. Implicit inside numeric/depth slots auto-compress and never move outside automatically.',
        action='append',
        default=[])
    parser.add_argument(
        '--circular_track_table',
        metavar='TSV',
        help='TSV manifest for circular track slots and axis placement.',
        type=str)
    parser.add_argument(
        '--circular_track_axis_index',
        type=argparse_integer,
        default=None,
        help='Axis boundary index for --circular_track_slot order. Slots before this index render outside; slots at or after it render inside.')


def _add_gc_track_size_args(parser: argparse.ArgumentParser) -> None:
    """The GC content and GC skew track sizes, and the GC content axis."""

    parser.add_argument(
        '--gc_content_width',
        help='GC content track width for circular mode (in px; must be > 0).',
        type=argparse_float)
    parser.add_argument(
        '--gc_content_radius',
        help='GC content track center radius for circular mode (as a ratio of base radius; must be > 0).',
        type=argparse_float)
    _add_gc_content_axis_args(parser)
    parser.add_argument(
        '--gc_skew_width',
        help='GC skew track width for circular mode (in px; must be > 0).',
        type=argparse_float)
    parser.add_argument(
        '--gc_skew_radius',
        help='GC skew track center radius for circular mode (as a ratio of base radius; must be > 0).',
        type=argparse_float)


def _circular_parser() -> argparse.ArgumentParser:
    """The ``gbdraw circular`` parser, with its options in the order ``--help`` lists them."""

    parser = argparse.ArgumentParser(
        description=(
            "Generate genome diagrams in PNG/PDF/SVG/PS/EPS. "
            "By default, diagrams for multiple entries are saved separately. "
            "Use --multi_record_canvas to place multiple records on one grid canvas."
        )
    )
    _add_circular_input_output_args(parser)
    add_color_args(parser)
    add_analysis_args(parser)
    _add_species_strain_args(parser)
    add_feature_args(parser)
    _add_feature_shape_arg(parser)
    _add_arrow_geometry_args(parser)
    _add_block_stroke_args(parser)
    _add_circular_stroke_args(parser)
    _add_circular_title_args(parser)
    _add_circular_center_args(parser)
    parser.add_argument(
        '--label_font_size',
        help='Label font size (optional; default: 14 (pt) for genomes < 75 kb, 8 for genomes >= 75 kb)',
        type=argparse_float)
    _add_format_arg(parser)
    _add_figure_size_args(parser)
    _add_gc_skew_toggle_args(
        parser,
        show_gc_default=CIRCULAR_MODE_PROFILE.show_gc,
        show_skew_default=CIRCULAR_MODE_PROFILE.show_skew,
        typed_dests=frozenset({"show_skew"}),  # One circle draws no skew unless asked (``_shown_skew``).
    )
    _add_depth_track_source_args(parser, mode="circular")
    _add_depth_track_tick_args(parser)
    add_feature_density_args(parser)
    _add_conservation_input_args(parser)
    _add_retired_option_args(parser, mode="circular")
    _add_conservation_ring_args(parser)
    _add_ring_losat_args(parser)
    _add_comparison_filter_args(parser)
    _add_depth_color_arg(parser)
    parser.add_argument(
        '--depth_width',
        help='Depth track width for circular mode (in px; must be > 0).',
        type=argparse_float)
    _add_depth_axis_args(parser)
    parser.add_argument(
        '-l',
        '--legend',
        help='Legend position (default: "right"; "left", "right", "top", "bottom", "upper_left", "upper_right", "lower_left", "lower_right", "none")',
        type=str,
        default="right")
    _add_record_layout_args(parser)
    _add_multi_record_canvas_args(parser)
    _add_circular_track_args(parser)
    _add_circular_label_mode_args(parser)
    add_label_args(parser)
    _add_label_ellipse_args(parser)
    _add_circular_scale_args(parser)
    _add_circular_track_slot_args(parser)
    _add_gc_track_size_args(parser)
    _add_legend_size_args(parser)
    _add_font_args(parser)
    add_session_args(parser)
    return parser


def _validate_circular_args(
    parser: argparse.ArgumentParser, args: argparse.Namespace, *, legacy_session: bool
) -> None:
    """Refuse option values and combinations the circular command cannot draw, in this order, and normalize
    ``--scale_interval`` and ``--circular_track_slot`` in place. ``legacy_session``: the replay of a released
    legacy Session, which reads its stored font sizes and keeps the fallbacks it was drawn with."""

    if legacy_session:
        read_stored_font_size_options(parser, args)
    validate_input_args(parser, args)
    validate_label_args(parser, args)
    if not legacy_session:
        validate_output_choice_args(args)
    if args.records_table and args.multi_record_position:
        parser.error("--records_table cannot be combined with --multi_record_position; use row and column table columns instead.")
    _validate_conservation_input_args(parser, args)
    if args.circular_track_table and args.circular_track_order:
        parser.error("--circular_track_table cannot be combined with --circular_track_order")
    if args.circular_track_table and args.circular_track_slot:
        parser.error("--circular_track_table cannot be combined with --circular_track_slot")
    if args.circular_track_table and args.circular_track_axis_index is not None:
        parser.error("--circular_track_table cannot be combined with --circular_track_axis_index")
    require_positive_args(parser, args, "feature_width", "gc_content_width", "gc_content_radius")
    validate_gc_content_axis_args(parser, args)
    validate_scale_interval_arg(parser, args, legacy_session=legacy_session)
    require_positive_args(
        parser,
        args,
        "gc_skew_width",
        "gc_skew_radius",
        "depth_width",
        "conservation_ring_width",
        "conservation_ring_gap",
    )
    if args.center_reserved_radius is not None and args.center_reserved_radius < 0:
        parser.error("--center_reserved_radius must be a finite number >= 0")
    if args.track_narrowing_min_px is not None and args.track_narrowing_min_px < 0:
        parser.error("--track_narrowing_min_px must be a finite number >= 0")
    validate_comparison_threshold_args(parser, args)
    _validate_ring_sources(parser, args)
    validate_depth_axis_args(parser, args)
    require_positive_args(parser, args, "circular_label_spacing")
    _normalize_circular_track_slots(parser, args, legacy_session=legacy_session)
    if args.multi_record_min_radius_ratio <= 0 or args.multi_record_min_radius_ratio > 1:
        parser.error("--multi_record_min_radius_ratio must be > 0 and <= 1")
    if args.multi_record_column_gap_ratio < 0:
        parser.error("--multi_record_column_gap_ratio must be a finite number >= 0")
    if args.multi_record_row_gap_ratio < 0:
        parser.error("--multi_record_row_gap_ratio must be a finite number >= 0")


def _validate_ring_sources(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Refuse a ring LOSAT program or output folder that cannot run, and ring labels, colors or genomes without
    the ring sources they align with."""

    if args.losat == "losatp":
        raise ValidationError(
            "Circular similarity rings run --losat losatn or --losat tlosatx; "
            "losatp is not available for rings.",
            diagnostic={
                "code": "COMPARISON_INPUT",
                "reason": "RING_LOSAT_PROGRAM",
                "field": "losat",
                "program": "losatp",
            },
        )
    if args.losat_output_dir and not args.losat:
        parser.error("--losat_output_dir requires --losat")
    ring_genomes = args.conservation_genome or args.conservation_sequence
    genome_flag = "--conservation_genome" if args.conservation_genome else "--conservation_sequence"
    ring_sources = args.conservation_blast if not args.losat else ring_genomes
    ring_source_flag = "--conservation_blast" if not args.losat else genome_flag
    if args.conservation_labels and not ring_sources and not args.conservation_table:
        parser.error(f"--conservation_labels requires {ring_source_flag}")
    if args.conservation_colors and not ring_sources and not args.conservation_table:
        parser.error(f"--conservation_colors requires {ring_source_flag}")
    if not args.losat:
        if ring_genomes and not args.conservation_blast:
            parser.error(f"{genome_flag} requires --conservation_blast or --losat")
        if ring_genomes and len(ring_genomes) != len(args.conservation_blast):
            parser.error(
                "--conservation_genome must give one genome per --conservation_blast"
                if args.conservation_genome
                else "--conservation_sequence must provide one source per --conservation_blast"
            )


def _normalize_circular_track_slots(
    parser: argparse.ArgumentParser, args: argparse.Namespace, *, legacy_session: bool
) -> None:
    """Refuse a track order with slot specs, or an axis index without either; parse ``--circular_track_slot``
    into slots in place, and check the axis index against the slots."""

    if args.circular_track_order and args.circular_track_slot:
        parser.error("--circular_track_order cannot be combined with --circular_track_slot")
    if args.circular_track_axis_index is not None and not (args.circular_track_order or args.circular_track_slot):
        parser.error("--circular_track_axis_index requires --circular_track_slot or --circular_track_order")
    circular_track_slots_for_validation: list[CircularTrackSlot] | None = None
    if args.circular_track_order:
        try:
            circular_track_slots_for_validation = circular_track_slots_from_order(args.circular_track_order)
        except Exception as exc:
            parser.error(str(exc))
    if args.circular_track_slot:
        try:
            circular_track_slots_for_validation = parse_circular_track_slots(
                args.circular_track_slot,
                _allow_legacy_transport=legacy_session,
            )
        except Exception as exc:
            parser.error(str(exc))
        args.circular_track_slot = circular_track_slots_for_validation
    if args.circular_track_axis_index is not None:
        try:
            normalize_listed_circular_track_slots(
                circular_track_slots_for_validation or [],
                args.circular_track_axis_index,
            )
        except Exception as exc:
            parser.error(str(exc))


def _get_args(
    args: Sequence[str] | None,
    *,
    _allow_legacy_track_transport: bool = False,
) -> argparse.Namespace:
    """
    Parses command-line arguments for generating circular genome diagrams.

    Args:
        args (list of str): Command-line arguments.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments.

    This function defines and parses the command-line arguments required for
    generating genome diagrams in a circular layout. It includes options for input
    GenBank files, color customization, output formats, and various Configurator for
    visualizing GC content, GC skew, and specific genomic features.
    """
    parser = _circular_parser()
    namespace = parser.parse_args(args)
    _validate_circular_args(
        parser, namespace, legacy_session=_allow_legacy_track_transport
    )
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
            mode="circular",
            temp_dir=Path(temp_dir),
            output_override=session_request.output,
            format_override=session_request.format,
        )
        args = _get_args(
            list(run_spec.args),
            _allow_legacy_track_transport=True,
        )
        if args.conservation_blast and not args.losat and not args.conservation_labels:
            # D-03, OV-202: an unlabelled ring keeps the file name its original run
            # drew, which the Session binds; the replay reads a temporary copy.
            names = {run_spec.args[binding.argIndex]: binding.name for binding in run_spec.file_bindings}
            args.conservation_labels = [names.get(path, Path(path).name) for path in args.conservation_blast]
        args.overwrite = session_request.overwrite
        args.embed_fonts = session_request.embed_fonts
        png = legacy_replay_png_dpi(document, session_request, args.format)
        args.dpi = png.dpi
        args.save_session = session_request.save_session
        args.session_output = session_request.session_output
        args._allow_legacy_track_transport = True
        args._require_canonical_session = bool(
            session_request.save_session
            or session_request.session_output
        )
        with noting_session_png_dpi(png):
            run_result = run_circular_from_namespace(args)
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


def circular_main(cmd_args) -> None:
    """
    Main function for generating circular genome diagrams.

    Args:
        cmd_args (argparse.Namespace): Command-line arguments parsed by _get_args function.

    This function orchestrates the creation of circular genome diagrams. It processes
    input GenBank files and applies user-specified visualization Configurator. The diagrams
    can include GC content, GC skew, and various genomic features. Outputs are saved
    in specified file formats, named using the accession number of the sequence.

    Key steps include:
    - Loading GenBank records and extracting relevant data.
    - Configuring the visualization Configurator and canvas.
    - Plotting the circular diagrams with genomic features and GC-related tracks.
    - Generating output files in specified formats.
    """
    setup_logging()
    session_request = parse_session_pre_args(cmd_args, mode="circular")
    if session_request is not None:
        # Before the Session is read, so its font names resolve to these files.
        register_font_files(session_request.font_files)
        document = load_session_document(session_request.session_path)
        if not render_canonical_session_if_present(
            document,
            mode="circular",
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
                document.drawing(session_request.drawings[0], mode="circular")
            replay_legacy_session(document, session_request)
        return

    args = _get_args(cmd_args)
    if args.losat_output_dir:
        preflight_output_paths(
            (Path(args.losat_output_dir) / LOSAT_CONSERVATION_OUTPUT_NAME,),
            overwrite=bool(args.overwrite),
        )
        # Before the search: a ring whose files the manifest cannot list stops the run here (D-123).
        _conservation_table_genomes(
            _ring_genome_files(_cli_conservation_genomes(args)) or (), Path(args.losat_output_dir)
        )
    run_result = run_circular_from_namespace(args)
    if args.losat_output_dir:
        write_losat_output_files(
            Path(args.losat_output_dir),
            _ring_output_files(run_result, Path(args.losat_output_dir)),
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


def _cli_ring_losat_search(args: argparse.Namespace) -> LosatSearchOptions | None:
    if not args.losat:
        return None
    try:
        return LosatSearchOptions(
            program=args.losat,
            losatn_task=args.losatn_task,
            record_gencodes=(args.losat_gencode,) if args.losat_gencode else (),
            runtime=LosatRuntimeOptions(
                losat_executable=args.losat_bin,
                ncbi_blast_executable=args.ncbi_blast_bin,
                threads=args.losat_threads,
            ),
        )
    except ValidationError as exc:
        raise _cli_ring_error(exc, args) from exc


def _cli_ring_error(exc: ValidationError, args: argparse.Namespace) -> ValidationError:
    """Name the CLI flag in a typed ring diagnostic (the typed layer validates)."""

    diagnostic = dict(exc.diagnostic or {})
    flags = {
        **_CLI_RING_FIELD_FLAGS,
        **({"conservation_sequence_files": "--conservation_genome"} if args.conservation_genome else {}),
    }
    flag = flags.get(str(diagnostic.get("field")))
    reason = diagnostic.get("reason")
    if flag is None or reason not in {"LOSAT_OPTION_PROGRAM", "RING_LOSAT_INPUT", "GENETIC_CODE"}:
        return exc
    if reason == "LOSAT_OPTION_PROGRAM":
        program = f"--losat {args.losat}" if args.losat else "no --losat"
        message = f"{flag} does not apply to {program}."
    elif flag == "--conservation_blast":
        message = (
            "--conservation_blast cannot be combined with --losat; the search "
            "builds the rings from --conservation_sequence."
        )
    elif flag == "--conservation_reference":
        message = (
            "--losat uses the reference file as the subject; "
            "--conservation_reference must be auto or subject."
        )
    else:
        message = str(exc)
        for typed_name, cli_flag in flags.items():
            message = message.replace(typed_name, cli_flag)
    return ValidationError(message, diagnostic=diagnostic)


def _cli_ring_options(*, _args: argparse.Namespace, **values) -> CircularDiagramOptions:
    try:
        return CircularDiagramOptions(**values)
    except ValidationError as exc:
        rewritten = _cli_ring_error(exc, _args)
        if rewritten is exc:
            raise
        raise rewritten from exc


def _conservation_table_genomes(
    genomes: Sequence[RingGenomeFiles], output_dir: Path
) -> tuple[RingGenomeFiles, ...]:
    """Each ring's files as the ``--losat_output_dir`` manifest names them, refused when it could not list them.

    A path is relative when next to the output directory, so both can move together.
    """

    resolved_dir = output_dir.resolve()

    def table_path(path: str) -> str:
        absolute = str(Path(path).resolve())
        try:
            relative = os.path.relpath(absolute, resolved_dir)
        except ValueError:
            relative = absolute
        return relative if relative.split(os.sep).count("..") <= 1 else absolute

    written = tuple(tuple(table_path(path) for path in genome) for genome in genomes)
    for genome in written:
        require_listable_comparison_sequences(genome)
    return written


def _ring_output_files(run_result: DiagramRunResult, output_dir: Path) -> list[tuple[str, str]]:
    """Per-ring tables in the search frame and a ``--conservation_table`` manifest (design 3.7)."""

    # A Circular run result always carries its Circular request.
    request = cast(
        CircularDiagramRequest | CircularBatchRequest, run_result.canonical_request
    )
    options = request.options
    rings = tuple(getattr(options, "conservation_search_results", None) or ())
    if not rings:
        raise ValidationError(
            "LOSAT output requested, but no ring table was produced.",
            diagnostic={"code": "LOSAT_RUNTIME", "reason": "OUTPUT"},
        )
    genomes = _conservation_table_genomes(options.conservation_genome_files, output_dir)
    labels = tuple(options.conservation_labels or ())
    colors = tuple(options.conservation_colors or ())
    files: list[tuple[str, str]] = []
    names = unique_losat_filenames(
        [ring.name for ring in rings], reserved=(LOSAT_CONSERVATION_OUTPUT_NAME,)
    )
    entries: list[ConservationTableEntry] = []
    for index, (ring, candidate) in enumerate(zip(rings, names)):
        files.append((candidate, ring.text))
        entries.append(
            ConservationTableEntry(
                blast=candidate,
                comparison_sequences=genomes[index] if index < len(genomes) else (),
                label=labels[index] if index < len(labels) else "",
                color=colors[index] if index < len(colors) else "",
            )
        )
    files.append((LOSAT_CONSERVATION_OUTPUT_NAME, conservation_table_text(entries)))
    return files


def _drawn_request(request: DiagramRequest) -> DiagramRequest:
    """The request a run's Session saves: one record drawn as separate diagrams is one single drawing.

    The command has no choice between one drawing and separate drawings; it draws each record of its input on its
    own canvas unless ``--multi_record_canvas`` asks for a grid. When the records come to one, the Session saves
    that drawing as the Web saves its one-record drawing, so the Web opens it on that record with its crop in the
    single-record fields (OV-3297, PD-OI-044). A grid of one record stays a grid. The figure and its file do not
    change: a batch already draws each of its records as a single-record diagram.
    """

    if not isinstance(request, CircularBatchRequest) or len(request.records) != 1:
        return request
    return CircularDiagramRequest(
        records=request.records,
        options=request.options,
        output=request.outputs[0],
        grouping="single",
        record_options=request.record_options,
    )


class _CircularTrackInput(NamedTuple):
    """The track list the options name: a ``--circular_track_table``, else ``--circular_track_slot`` and
    ``--circular_track_order``."""

    table: CircularTrackTable | None
    order: str | None
    slot_specs: Sequence[str | CircularTrackSlot]
    axis_index: int | None
    allow_legacy_transport: bool


def _circular_track_input(args: argparse.Namespace) -> _CircularTrackInput:
    """Read ``--circular_track_table`` (it wins over the slot and order options) or take the listed slots."""

    allow_legacy_transport = bool(getattr(args, "_allow_legacy_track_transport", False))
    table = (
        read_circular_track_table(args.circular_track_table, _allow_legacy_transport=allow_legacy_transport)
        if args.circular_track_table
        else None
    )
    if table:
        slot_specs = parse_circular_track_slots(table.slot_specs, _allow_legacy_transport=allow_legacy_transport)
        return _CircularTrackInput(table, None, slot_specs, table.axis_index, allow_legacy_transport)
    return _CircularTrackInput(
        None,
        args.circular_track_order,
        list(args.circular_track_slot or []),
        args.circular_track_axis_index,
        allow_legacy_transport,
    )


class _CircularInputs(NamedTuple):
    """What the command reads and parses from its options before it builds the config."""

    feature_shapes: dict[str, FeatureShape]
    record_manifest: RecordInputManifest
    multi_record_positions: list[str]
    show_skew: bool
    gc_content_mode: str | None
    tracks: _CircularTrackInput
    depth_tracks: tuple[DepthTrackInput, ...] | None


def _read_circular_inputs(args: argparse.Namespace) -> _CircularInputs:
    """The feature shapes, the records, the track list and the depth tracks, read in this order."""

    feature_shapes = parse_feature_shape_overrides(args.feature_shape)
    record_manifest = apply_record_display_cli_options(
        record_input_manifest_from_table(args.records_table)
        if args.records_table
        else record_input_manifest_from_paths(
            gbk_paths=args.gbk,
            gff_paths=args.gff,
            fasta_paths=args.fasta,
            cardinalities=_circular_cli_record_cardinality(),
        ),
        args,
    )
    multi_record_positions = (
        []
        if args.records_table
        else [str(position) for position in (args.multi_record_position or [])]
    )
    show_skew = _shown_skew(args)
    tracks = _circular_track_input(args)
    gc_content_mode = cast("str | None", _profile_values(args).get("objects.gc_content.mode"))
    depth_tracks = depth_track_inputs_from_args(args)
    if not args.multi_record_canvas and multi_record_positions:
        logger.info(
            "Ignoring --multi_record_position because --multi_record_canvas is disabled."
        )
    return _CircularInputs(
        feature_shapes,
        record_manifest,
        multi_record_positions,
        show_skew,
        gc_content_mode,
        tracks,
        depth_tracks,
    )


def _circular_config(args: argparse.Namespace, inputs: _CircularInputs) -> dict:
    """The packaged config with every config option the command was given (the CLI writes its whole config)."""

    config_dict: dict = load_config_toml('gbdraw.data', 'config.toml')
    filtering = config_dict.setdefault("labels", {}).setdefault("filtering", {})
    override_candidates: dict[str, object | None] = {
        **figure_size_overrides(args),
        **legend_row_overrides(args),
        **style_overrides(args, mode="circular"),
        "objects.features.line_stroke_color": args.line_stroke_color,
        "labels.circular.scope": "outer" if args.labels == "out" else args.labels,
        "canvas.circular.track_type": args.track_type,
        "canvas.circular.track_narrowing": args.track_narrowing,
        "canvas.circular.track_narrowing_min_px": args.track_narrowing_min_px,
        "canvas.strandedness": args.separate_strands,
        "canvas.resolve_overlaps": args.resolve_overlaps,
        "canvas.feature_overlap_tolerance_bp": args.feature_overlap_tolerance_bp,
        "canvas.show_gc": args.show_gc,
        **gc_content_overrides(args, gc_content_mode=inputs.gc_content_mode),
        "canvas.show_skew": inputs.show_skew,
        "canvas.show_depth": bool(args.depth_track),
        "objects.depth.fill_color": args.depth_color,
        **depth_overrides(args),
        "objects.definition.circular.font_size": args.definition_font_size,
        "objects.definition.circular.plot_title_font_size": args.plot_title_font_size,
        # Only a hidden definition is an override, so a default run's request is unchanged.
        "objects.definition.circular.show": False if args.hide_definition else None,
        "objects.conservation.color_scheme": args.conservation_color_scheme,
        CIRCULAR_IDENTITY_SCALE_PATH: list(args.identity_scale) if args.identity_scale else None,
        **label_style_overrides(args, mode="circular"),
        "labels.circular.placement": args.label_placement,
        "labels.filtering.raw": label_filtering_override(filtering, args.label_blacklist),
        "labels.unified_adjustment.outer_labels.x_radius_offset": args.outer_label_x_radius_offset,
        "labels.unified_adjustment.outer_labels.y_radius_offset": args.outer_label_y_radius_offset,
        "labels.unified_adjustment.inner_labels.x_radius_offset": args.inner_label_x_radius_offset,
        "labels.unified_adjustment.inner_labels.y_radius_offset": args.inner_label_y_radius_offset,
        "objects.scale.interval": args.scale_interval,
        "objects.scale.show": not bool(args.hide_scale),
        "objects.ticks.tick_labels.font_size": args.tick_label_font_size,
        "labels.spacing.circular": args.circular_label_spacing,
        **stroke_width_overrides(args, mode="circular"),
        **paired_overrides({"labels.font_size": args.label_font_size}),
        **legend_size_overrides(args),
    }
    config_dict = modify_config_dict(
        config_dict,
        {path: value for path, value in override_candidates.items() if value is not None},
    )
    return copy.deepcopy(config_dict)


def _with_shortcut_geometry(slot: CircularTrackSlot, args: argparse.Namespace) -> CircularTrackSlot:
    """The slot with the radius and width the geometry shortcut options give its track."""

    if slot.id == "features" and args.feature_width is not None:
        return replace(slot, width=ScalarSpec(float(args.feature_width), "px"))
    if str(slot.renderer) == "depth" and args.depth_width is not None:
        return replace(slot, width=ScalarSpec(float(args.depth_width), "px"))
    for slot_id, radius, width in (
        ("gc_content", args.gc_content_radius, args.gc_content_width),
        ("gc_skew", args.gc_skew_radius, args.gc_skew_width),
    ):
        if slot.id == slot_id:
            return replace(
                slot,
                radius=ScalarSpec(float(radius), "factor") if radius is not None else slot.radius,
                width=ScalarSpec(float(width), "px") if width is not None else slot.width,
            )
    return slot


def _circular_track_slots(
    args: argparse.Namespace, tracks: _CircularTrackInput, *, show_skew: bool, depth_track_count: int
) -> Sequence[str | CircularTrackSlot] | None:
    """The slots the request draws: the listed ones; else those of ``--circular_track_order``, or of the default
    order when a geometry shortcut is given, with the shortcut sizes; else None (the mode's default slots)."""

    hide_scale = bool(args.hide_scale)
    shortcut_geometry = (
        args.feature_width,
        args.depth_width,
        args.gc_content_width,
        args.gc_content_radius,
        args.gc_skew_width,
        args.gc_skew_radius,
    )
    shortcut_geometry_requested = any(value is not None for value in shortcut_geometry)
    if tracks.slot_specs and shortcut_geometry_requested:
        slot_source_option = "--circular_track_table" if tracks.table else "--circular_track_slot"
        raise ValidationError(
            f"Circular geometry shortcut options cannot be combined with {slot_source_option}; "
            "put r= and/or w= on the matching circular track slot."
        )

    slots: Sequence[str | CircularTrackSlot] | None = None
    if tracks.slot_specs:
        slots = tracks.slot_specs
    elif tracks.order or shortcut_geometry_requested:
        # A shortcut without a track list draws the default order, whose ticks --hide_scale leaves out.
        slots = circular_track_slots_from_order(
            tracks.order or "features,ticks,depth,gc_content,gc_skew",
            show_ticks=bool(tracks.order) or not hide_scale,
            show_depth=bool(args.depth_track),
            depth_track_count=max(1, depth_track_count),
            show_gc=args.show_gc,
            show_skew=show_skew,
            dinucleotide=args.nt.upper(),
        )

    if (
        hide_scale
        and (tracks.table or tracks.order or tracks.slot_specs)
        and any(
            slot.enabled and str(slot.renderer) == "ticks"
            for slot in parse_circular_track_slots(
                slots or (),
                _allow_legacy_transport=tracks.allow_legacy_transport,
            )
        )
    ):
        logger.warning(
            "--hide_scale does not suppress an enabled ticks slot in a "
            "user-explicit Circular track list; the explicit slot wins."
        )

    if slots is None or tracks.slot_specs:
        return slots
    return [_with_shortcut_geometry(slot, args) for slot in parse_circular_track_slots(slots)]


def _circular_request_options(
    args: argparse.Namespace,
    inputs: _CircularInputs,
    config: dict,
    track_slots: Sequence[str | CircularTrackSlot] | None,
) -> CircularDiagramOptions:
    """The typed Circular options; the typed layer validates them, and ring diagnostics name the CLI flags."""

    return _cli_ring_options(
        config=config,
        tracks=CircularRequestTrackOptions(
            circular_track_slots=track_slots,
            circular_track_axis_index=inputs.tracks.axis_index,
            center_reserved_radius=args.center_reserved_radius,
        ),
        output=CircularOutputOptions(
            legend=args.legend,
            # argparse restricts --plot_title_position to the Circular choices.
            plot_title_position=cast(
                'Literal["none", "top", "bottom"]', args.plot_title_position
            ),
        ),
        **mode_option_fields(args, feature_shapes=inputs.feature_shapes, depth_tracks=inputs.depth_tracks),
        conservation_blast_files=list(args.conservation_blast or []) or None,
        conservation_sequence_files=_cli_conservation_genomes(args),
        conservation_table_file=args.conservation_table,
        conservation_reference=args.conservation_reference,
        conservation_labels=list(args.conservation_labels or []) or None,
        conservation_colors=conservation_colors_from_cli(args.conservation_colors),
        conservation_ring_width=args.conservation_ring_width,
        conservation_ring_gap=args.conservation_ring_gap,
        losat_search=_cli_ring_losat_search(args),
        conservation_losat_gencodes=(
            tuple(args.conservation_losat_gencode)
            if args.conservation_losat_gencode
            else None
        ),
        _args=args,
        keep_full_definition_with_plot_title=args.keep_full_definition_with_plot_title,
        species=args.species or None,
        strain=args.strain or None,
        origin_feature=args.origin_feature,
    )


def _circular_request(
    args: argparse.Namespace,
    inputs: _CircularInputs,
    options: CircularDiagramOptions,
    out_formats: list[str],
    dpi: int | None,
) -> CircularDiagramRequest | CircularBatchRequest:
    """One diagram of every record (a grid canvas or one circle), else a batch of one diagram per record."""

    record_manifest = inputs.record_manifest
    if not (args.multi_record_canvas or args.one_circle):
        return CircularBatchRequest(
            records=record_manifest.records,
            options=options,
            output_policy=CircularBatchOutputPolicy(
                output_prefix=args.output,
                formats=tuple(out_formats),
                overwrite=args.overwrite,
                embed_fonts=args.embed_fonts,
                dpi=dpi,
            ),
            record_options=record_manifest.record_options,
        )
    layout = CircularMultiRecordOptions(
        # argparse restricts --multi_record_size_mode to these choices.
        multi_record_size_mode=cast('Literal["linear", "auto", "equal"]', args.multi_record_size_mode),
        multi_record_min_radius_ratio=args.multi_record_min_radius_ratio,
        multi_record_column_gap_ratio=args.multi_record_column_gap_ratio,
        multi_record_row_gap_ratio=args.multi_record_row_gap_ratio,
        multi_record_positions=inputs.multi_record_positions or None,
        **_one_circle_layout_fields(args),
    )
    output = render_output_request(
        args.output or "out",
        out_formats,
        overwrite=args.overwrite,
        embed_fonts=args.embed_fonts,
        dpi=dpi,
    )
    if args.output is None:
        # Without -o the file is named after the first record.
        output = replace(output, resolve_prefix_from_first_record=True)
    return CircularDiagramRequest(
        records=record_manifest.records,
        options=options,
        layout=layout,
        output=output,
        record_options=record_manifest.record_options,
        grouping="grid",
    )


def _render_circular(
    args: argparse.Namespace, request: CircularDiagramRequest | CircularBatchRequest
) -> RequestRenderResult | CircularBatchRenderResult:
    """Render the request; with a Session to save, refuse an existing sidecar or output before drawing."""

    if not (args.save_session or args.session_output):
        return render_request(request)
    request_plan = plan_request(request)
    materialized_request = request_plan.request
    sidecar_path = resolve_session_sidecar_path(
        explicit_path=args.session_output,
        output_prefix=args.output,
        outputs=diagram_request_rendered_svgs(materialized_request),
    )
    preflight_session_sidecar_if_requested(
        save_session=True,
        session_output=str(sidecar_path),
        output_prefix=None,
        diagram_output_paths=diagram_request_output_paths(
            materialized_request
        ),
        overwrite=bool(args.overwrite),
    )
    batch_outputs_preflighted = isinstance(
        request_plan,
        CircularBatchRequestPlan,
    )
    request_plan.preflight_outputs()
    prepared = build_request_plan_diagram(request_plan)
    return render_prepared_request(
        prepared,
        batch_outputs_preflighted=batch_outputs_preflighted,
        include_feature_catalog=True,
    )


def _circular_run_result(
    rendered: RequestRenderResult | CircularBatchRenderResult, out_formats: list[str]
) -> DiagramRunResult:
    """The run's SVGs, the feature metadata and interactive contexts of each, and the request the Session saves."""

    if isinstance(rendered, CircularBatchRenderResult):
        rendered_items = rendered.items
    elif isinstance(rendered, RequestRenderResult):
        rendered_items = (rendered,)
    else:  # pragma: no cover - render_request exhausts the Circular result union.
        raise ValidationError("Circular request renderer returned an unsupported result.")

    materialized_request = _drawn_request(rendered.request)

    outputs: list[RenderedSvg] = []
    run_results: list[RenderedResult] = []
    session_feature_metadata: list[Mapping[str, Any]] = []
    session_biological_feature_metadata: list[Mapping[str, Any]] = []
    session_interactive_contexts = []
    log_render_notices(rendered_items, warned=set())
    for rendered_item in rendered_items:
        if not rendered_item.output_paths:
            raise ValidationError("Circular request renderer did not produce an SVG output.")
        svg_path = Path(rendered_item.output_paths[0])
        # Every rendered item (also of a batch) comes from a single-diagram request.
        item_request = cast(CircularDiagramRequest, rendered_item.request)
        rendered_svg = RenderedSvg(
            output_prefix=str(svg_path.with_suffix("")),
            svg_path=svg_path,
            result_name=item_request.output.output_prefix,
        )
        outputs.append(rendered_svg)
        run_results.append(RenderedResult(rendered_item, rendered_svg.result_name))
        interactive_context = rendered_item.interactive_context
        session_interactive_contexts.append(interactive_context)
        if interactive_context is not None:
            session_feature_metadata.extend(interactive_context.features)
            session_biological_feature_metadata.extend(
                interactive_context.biological_features
            )

    return DiagramRunResult(
        mode="circular",
        render_formats=tuple(out_formats),
        outputs=tuple(outputs),
        feature_metadata=tuple(session_feature_metadata),
        biological_feature_metadata=tuple(session_biological_feature_metadata),
        interactive_contexts=tuple(session_interactive_contexts),
        run_metadata=build_run_metadata("circular", run_results),
        canonical_request=materialized_request,
        losat_cache_entries=rendered.losat_cache_entries,
    )


def run_circular_from_namespace(args: argparse.Namespace) -> DiagramRunResult:
    """Run circular rendering from an already parsed argparse namespace."""

    inputs = _read_circular_inputs(args)
    config = _circular_config(args, inputs)
    out_formats = cli_output_formats(args.format)
    dpi = png_dpi(args.dpi, out_formats)
    track_slots = _circular_track_slots(
        args,
        inputs.tracks,
        show_skew=inputs.show_skew,
        depth_track_count=len(inputs.depth_tracks or ()),
    )
    options = _circular_request_options(args, inputs, config, track_slots)
    request = _circular_request(args, inputs, options, out_formats, dpi)
    return _circular_run_result(_render_circular(args, request), out_formats)


if __name__ == "__main__":
    # ``python -m gbdraw.circular``: circular_main sets up the log as ``gbdraw circular`` does.
    main_args = sys.argv[1:]
    if not main_args:
        main_args.append('--help')
    circular_main(main_args)
