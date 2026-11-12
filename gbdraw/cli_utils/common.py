#!/usr/bin/env python
# coding: utf-8

"""Shared CLI argument definitions and utilities.

This module provides common argument parsers and validation logic
used by both circular and linear diagram commands.
"""

import argparse
import csv
import logging
import sys
from dataclasses import dataclass, replace
from typing import Sequence, get_args

from gbdraw.analysis.feature_density import (
    FeatureDensityTrack,
    parse_feature_density_words,
    read_feature_density_table,
    unique_feature_density_tracks,
)
from gbdraw.cli_utils.losat_output import parse_positive_int
from gbdraw.comparisons.identity_scale import IDENTITY_SCALE_RULE, identity_scale_from_config
from gbdraw.config.models.canvas import (
    CANVAS_PADDING_PATH,
    FIGURE_HEIGHT_PATH,
    FIGURE_WIDTH_PATH,
    canvas_padding_from_config,
)
from gbdraw.config.models.objects import (
    SCALAR_AXIS_LABEL_STYLES,
    LEGEND_HIDDEN_ROWS_PATH,
    LEGEND_ROW_NAMES_PATH,
    LEGEND_ROW_ORDER_PATH,
    legend_row_keys_from_config,
    legend_row_names_from_config,
)
from gbdraw.core.figure_length import argparse_figure_length
from gbdraw.core.numeric import JS_WHITESPACE, argparse_float, argparse_integer, ascii_int, parse_ascii_finite_float
from gbdraw.core.fonts import RegisteredFont, complete_font_family, font_family_list_for, register_font_file_faces
from gbdraw.exceptions import ValidationError
from gbdraw.api.options import DepthTrackInput
from gbdraw.api.record_planning import RecordInputManifest, depth_track_inputs_from_cli
from gbdraw.api.requests import RecordCardinality, RecordDisplayOptions
from gbdraw.features.shapes import parse_feature_shape_assignment
from gbdraw.labels.label_style import (
    DEFAULT_LABEL_LINE,
    DEFAULT_LABEL_PROTECTION,
    LABEL_LINES,
    LABEL_PROTECTIONS,
    LabelStyleMode,
)
from gbdraw.labels.policy import LabelRenderingPolicy
from gbdraw.mode_profiles import CIRCULAR_MODE_PROFILE, DEFAULT_FEATURE_TYPES, ComparisonThresholds
from gbdraw.render.export import require_cairosvg
from gbdraw.io.colors import BUILTIN_COLOR_TABLE_PREFIX, BUILTIN_COLOR_TABLES, require_palette_name
from gbdraw.render.formats import (
    SVG_FORMAT,
    is_cairosvg_format,
    normalize_format_token,
    parse_output_formats,
)

logger = logging.getLogger(__name__)

_OUTPUT_FORMAT_HELP = (
    "Comma-separated list of output file formats "
    "(svg, interactive_svg, png, pdf, eps, ps; default: svg; "
    "png/pdf/eps/ps require CairoSVG)."
)
_TRACK_VISIBILITY_FLAG_SPECS = (
    ("gc", "show_gc", "GC content"),
    ("skew", "show_skew", "GC skew"),
)
_DEPTH_TRACK_HELP_BY_MODE = {
    "circular": (
        "Repeatable logical depth track. Provide one file for all records "
        "or one file per record."
    ),
    "linear": (
        "Repeatable logical depth track. Provide one file for all records "
        "or one file per input record."
    ),
}
# The rule ``--label_rendering`` meets in one mode, said in its help before the default.
_LABEL_RENDERING_RULE_BY_MODE: dict[LabelStyleMode, str] = {
    "circular": "",
    "linear": 'Cannot be combined with --label_placement above_feature except as "auto". ',
}
# The per-track depth options whose values are numbers, ``auto``, or empty.
_DEPTH_TRACK_NUMBER_OPTIONS = (
    "depth_track_large_tick_interval",
    "depth_track_small_tick_interval",
    "depth_track_tick_font_size",
)


CLI_LOG_FORMAT = "%(levelname)s: %(message)s"


def cli_label(levelname: str, message: str) -> str:
    """Return the label of one CLI line, shared by the log formatter and the CLI's own prints."""
    return CLI_LOG_FORMAT % {"levelname": levelname, "message": message}


def setup_logging() -> None:
    """Log INFO and above to standard error, each line labelled with its level.

    The command's entry points call it; importing a module never does. Standard output stays for what a
    command exists to print (help, a version, a path).
    """
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    if not root_logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter(CLI_LOG_FORMAT))
        root_logger.addHandler(handler)


def add_input_args(parser: argparse.ArgumentParser) -> None:
    """Add input file arguments (--gbk, --gff, --fasta)."""
    parser.add_argument(
        "--feature_placement_table", metavar="TSV", default=None,
        help="Feature placements: record, feature_selector, placement, and optional level columns.",
    )
    parser.add_argument(
        "--feature_override_table", metavar="TSV", default=None,
        help=(
            "Per-feature edits: record, feature_selector, and feature_visibility, "
            "label_visibility, or label_text columns."
        ),
    )
    parser.add_argument(
        "--feature_overlap_tolerance_bp", metavar="BP", type=argparse_integer, default=0,
        help="Non-negative permitted feature overlap in base pairs (default: 0).",
    )
    parser.add_argument(
        "--record_topology", choices=("auto", "linear", "circular"), default=None,
        help="Topology override for exactly one direct-input record (default: auto).",
    )
    parser.add_argument(
        "--display_start_coordinate", type=argparse_integer, default=None,
        help="1-based source coordinate at the display start; requires one complete circular record.",
    )
    parser.add_argument(
        "--gbk",
        metavar="GBK_FILE",
        help='GenBank, DDBJ, or EMBL flat file, plain or gzip-compressed',
        type=str,
        nargs='*')
    parser.add_argument(
        "--gff",
        metavar="GFF3_FILE",
        help="GFF3 file, plain or gzip-compressed (instead of --gbk; requires --fasta)",
        type=str,
        nargs='*')
    parser.add_argument(
        "--fasta",
        metavar="FASTA_FILE",
        help=(
            "FASTA file, plain or gzip-compressed: with --gff, its sequences; alone (instead of --gbk), "
            "each nucleotide entry is a record without features"
        ),
        type=str,
        nargs='*')


def _add_window_step_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '-w',
        '--window',
        help=(
            'window size (optional; default: 1kb for genomes < 1Mb, 10kb for genomes <10Mb, 100kb for genomes <100Mb, 1Mb for genomes >=100Mb, '
            'and at most a tenth of the drawn length)'
        ),
        type=argparse_integer)
    parser.add_argument(
        '-s',
        '--step',
        help=(
            'step size (optional; default: 100 bp for genomes < 1Mb, 1kb for genomes <10Mb, 10kb for genomes <100Mb, 100kb for genomes >=100Mb, '
            'and at most a tenth of the Auto window)'
        ),
        type=argparse_integer)


def _add_feature_shape_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--feature_shape',
        help=(
            'Feature rendering override (repeatable): TYPE=SHAPE where SHAPE is '
            'arrow, rectangle, or underlay.'
        ),
        type=parse_feature_shape_assignment_arg,
        action='append',
        default=[],
        metavar='TYPE=SHAPE',
    )


def _add_arrow_geometry_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--arrow_head_length_ratio',
        help=(
            'Arrow head length divided by full feature thickness: auto or a '
            'positive finite number (default: auto).'
        ),
        type=parse_arrow_head_length_ratio_arg,
    )
    parser.add_argument(
        '--arrow_shaft_width_ratio',
        help=(
            'Arrow shaft thickness divided by full head thickness: a '
            'finite number in (0, 1] (default: 1.0).'
        ),
        type=parse_arrow_shaft_width_ratio_arg,
    )


def _add_block_stroke_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--block_stroke_color',
        help='Block stroke color (str; default: "gray")',
        type=str)
    parser.add_argument(
        '--block_stroke_width',
        help='Block stroke width (optional; float; default: 2 pt for genomes shorter than the short/long boundary (75 kb circular, 100 kb linear), 0 pt otherwise)',
        type=argparse_float)


def _add_format_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '-f',
        '--format',
        help=_OUTPUT_FORMAT_HELP,
        type=str,
        default="svg")


def add_dpi_arg(parser: argparse.ArgumentParser) -> None:
    """``--dpi``, the one definition for the diagram commands, ``gbdraw render`` and the ``--session`` replays."""

    parser.add_argument(
        '--dpi',
        help=(
            'PNG resolution in pixels per inch (positive integer; default: 300, or with --session the PNG DPI the '
            'Session saved; 96 draws one pixel per SVG pixel). PNG only. The PNG also states it in its pHYs chunk.'
        ),
        type=parse_positive_int,
    )


def _add_figure_size_args(parser: argparse.ArgumentParser) -> None:
    """``--figure_width``, ``--figure_height`` and ``--dpi``, the same for both diagram commands."""

    for side in ("width", "height"):
        parser.add_argument(
            f'--figure_{side}',
            help=(
                f'{side.capitalize()} of the whole figure, legend and title included, with a unit: mm, cm, in, or px '
                '(for example 180mm). The figure is laid out again to fit (the Linear sequence width or the Circular '
                'radius changes); fonts and line widths keep their size. Default: auto.'
            ),
            type=argparse_figure_length,
            metavar='LENGTH|auto',
        )
    add_dpi_arg(parser)
    parser.add_argument(
        '--canvas_padding',
        help=(
            'Room around the figure, outside its edge padding, in px: TOP,RIGHT,BOTTOM,LEFT (for example 0,40,0,0). '
            'With --figure_width or --figure_height the figure fits inside the requested size with the padding. '
            'Default: 0,0,0,0.'
        ),
        type=argparse_canvas_padding,
        metavar='TOP,RIGHT,BOTTOM,LEFT',
    )


def argparse_canvas_padding(text: str) -> list[float]:
    """argparse ``type=`` of ``--canvas_padding TOP,RIGHT,BOTTOM,LEFT`` (px), read by the ``canvas.padding`` leaf's
    validator."""

    parts = [part.strip() for part in str(text).split(",")]
    try:
        sides = [parse_ascii_finite_float(part) for part in parts]
        return list(canvas_padding_from_config(sides))
    except (TypeError, ValueError, ValidationError) as exc:
        raise argparse.ArgumentTypeError(
            f"expected four px numbers TOP,RIGHT,BOTTOM,LEFT, each 0 or more (for example 0,40,0,0); got {text!r}"
        ) from exc


def figure_size_overrides(args: argparse.Namespace) -> dict[str, object]:
    """The config overrides of ``--figure_width``, ``--figure_height`` and ``--canvas_padding`` (``None``: not
    given; a padding of zeros is no override, so the request is the one without it)."""

    padding = args.canvas_padding if args.canvas_padding is not None and any(args.canvas_padding) else None
    return {FIGURE_WIDTH_PATH: args.figure_width, FIGURE_HEIGHT_PATH: args.figure_height, CANVAS_PADDING_PATH: padding}


# The columns of a ``--legend_rows`` table, in order.
LEGEND_ROWS_COLUMNS = ("key", "name", "hidden", "position")


@dataclass(frozen=True)
class LegendRowsTable:
    """A ``--legend_rows`` table read as the three Legend row edit leaves (``objects.legends.*``)."""

    row_names: tuple[tuple[str, str], ...]
    hidden_rows: tuple[str, ...]
    row_order: tuple[str, ...]


def read_legend_rows_table(path: str) -> LegendRowsTable:
    """argparse ``type=`` of ``--legend_rows FILE``: a TSV with the header ``key name hidden position``.

    ``key`` is the Legend row's key, trimmed as the Web's executor trims a key (JavaScript's ``trim()``); ``name`` the name it is given
    (blank: Python's); ``hidden`` ``true`` or ``false`` (blank: false); ``position`` its 1-based place in the Legend
    in ASCII digits, on every shown row or on none (none: Python's order; the rows the table does not list follow the
    listed ones); a hidden row may leave it blank, as the Web's edited order lists the shown rows only. A cell is read as written: a quote is part of it. A UTF-8 byte order mark is skipped.
    """

    def refuse(message: str) -> argparse.ArgumentTypeError:
        return argparse.ArgumentTypeError(f"{path}: {message}")

    try:
        with open(path, encoding="utf-8-sig", newline="") as handle:
            lines = [
                line
                for line in csv.reader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
                if any(cell.strip() for cell in line)
            ]
    except OSError as exc:
        raise refuse(f"cannot read the file ({exc.strerror})") from exc
    if not lines or tuple(cell.strip().lower() for cell in lines[0]) != LEGEND_ROWS_COLUMNS:
        raise refuse("the first line must be the header " + "\t".join(LEGEND_ROWS_COLUMNS))
    names: list[tuple[str, str]] = []
    hidden: list[str] = []
    positions: list[tuple[int, str]] = []
    keys: list[str] = []
    unplaced_shown = 0
    for number, cells in enumerate(lines[1:], start=2):
        if len(cells) > len(LEGEND_ROWS_COLUMNS) or not cells[0].strip():
            raise refuse(f"line {number} must have a key and at most four columns")
        raw_key, name, hidden_text, position = (*cells, "", "", "")[:4]
        key = raw_key.strip(JS_WHITESPACE)
        keys.append(key)
        if name.strip():
            names.append((key, name))
        flag = hidden_text.strip().lower()
        if flag not in ("", "true", "false"):
            raise refuse(f"line {number}: hidden must be true or false; got {hidden_text!r}")
        if flag == "true":
            hidden.append(key)
        if position.strip():
            try:
                place = ascii_int(position.strip())
            except ValueError:
                place = 0
            if place < 1:
                raise refuse(f"line {number}: position must be a whole number from 1; got {position!r}")
            positions.append((place, key))
        elif flag != "true":
            unplaced_shown += 1
    rows = len(lines) - 1
    if len(set(keys)) != rows:
        raise refuse("each key must be on one row")
    if positions and unplaced_shown:
        raise refuse("give a position on every shown row or on none")
    if len({place for place, _ in positions}) != len(positions):
        raise refuse("each position must be on one row")
    try:
        return LegendRowsTable(
            row_names=legend_row_names_from_config(names),
            hidden_rows=legend_row_keys_from_config(hidden, path=LEGEND_HIDDEN_ROWS_PATH),
            row_order=legend_row_keys_from_config([key for _, key in sorted(positions)], path=LEGEND_ROW_ORDER_PATH),
        )
    except ValidationError as exc:
        raise refuse(str(exc)) from exc


def legend_row_overrides(args: argparse.Namespace) -> dict[str, object]:
    """The config overrides of ``--legend_rows`` (``None``: not given or empty), each leaf only when it is set."""

    table: LegendRowsTable | None = args.legend_rows
    if table is None:
        return {}
    return {
        LEGEND_ROW_NAMES_PATH: [list(pair) for pair in table.row_names] or None,
        LEGEND_HIDDEN_ROWS_PATH: list(table.hidden_rows) or None,
        LEGEND_ROW_ORDER_PATH: list(table.row_order) or None,
    }


def png_dpi(dpi: int | None, formats: Sequence[str]) -> int | None:
    """``--dpi``, with a warning when none of ``formats`` is a PNG."""

    if dpi is not None and "png" not in formats:
        logger.warning("--dpi applies to PNG output only; no PNG is written (-f).")
    return dpi


def _add_overwrite_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--overwrite',
        help='Replace existing output files (default: refuse to overwrite).',
        action='store_true')


class _RetiredOptionAction(argparse.Action):
    """Reject a retired flag and name its replacement (exit status 2)."""

    def __init__(self, option_strings, dest, *, retired, **kwargs) -> None:
        self.retired = retired
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None) -> None:
        parser.error(self.retired.message(values))


def _add_retired_option_args(
    parser: argparse.ArgumentParser,
    *,
    mode: str,
) -> None:
    """Register the retired flags of one mode from the shared retired table."""

    from gbdraw.session_io import RETIRED_CLI_OPTIONS

    for retired in RETIRED_CLI_OPTIONS.values():
        if mode not in retired.modes:
            continue
        parser.add_argument(
            retired.option,
            action=_RetiredOptionAction,
            retired=retired,
            nargs=retired.nargs,
            dest=f"_retired_{retired.option.lstrip('-')}",
            default=argparse.SUPPRESS,
            help=argparse.SUPPRESS,
        )


class _TypedToggleAction(argparse.Action):
    """A track toggle that also notes it was typed (``<dest>_typed``), so a layout with defaults of its own (one
    circle) still yields to what the user typed."""

    def __init__(self, option_strings, dest, *, const: bool, help: str | None = None) -> None:
        super().__init__(option_strings, dest, nargs=0, const=const, default=not const, help=help)

    def __call__(self, parser, namespace, values, option_string=None) -> None:
        setattr(namespace, self.dest, self.const)
        setattr(namespace, f"{self.dest}_typed", True)


def _add_gc_skew_toggle_args(
    parser: argparse.ArgumentParser,
    *,
    show_gc_default: bool,
    show_skew_default: bool,
    typed_dests: frozenset[str] = frozenset(),
) -> None:
    """``--gc``/``--no-gc`` and ``--skew``/``--no-skew``; a toggle whose dest is in ``typed_dests`` also notes
    whether it was typed (``_TypedToggleAction``)."""

    for option_name, destination, label in _TRACK_VISIBILITY_FLAG_SPECS:
        group = parser.add_mutually_exclusive_group()
        typed = destination in typed_dests
        for option, shown, verb in ((f'--{option_name}', True, 'Show'), (f'--no-{option_name}', False, 'Hide')):
            group.add_argument(
                option,
                dest=destination,
                help=f'{verb} the {label} track.',
                **({'action': _TypedToggleAction, 'const': shown} if typed else {'action': 'store_true' if shown else 'store_false'}),
            )
    parser.set_defaults(
        show_gc=bool(show_gc_default),
        show_skew=bool(show_skew_default),
        **{f"{dest}_typed": False for dest in sorted(typed_dests)},
    )


def _add_depth_track_arg(
    parser: argparse.ArgumentParser,
    *,
    mode: str,
) -> None:
    parser.add_argument(
        '--depth_track',
        metavar='DEPTH',
        help=_DEPTH_TRACK_HELP_BY_MODE[mode],
        type=str,
        nargs='+',
        action='append')


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


def _add_depth_axis_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--depth_window',
        help='Depth aggregation window size. Defaults to one tenth of the GC/skew window, with a 100 bp minimum.',
        type=argparse_integer)
    parser.add_argument(
        '--depth_step',
        help='Depth aggregation step size. Defaults to one tenth of the GC/skew step.',
        type=argparse_integer)
    parser.add_argument(
        '--share_depth_axis',
        help='Use one depth y-axis scale across records.',
        action='store_true')
    parser.add_argument(
        '--depth_min',
        help='Minimum depth for clipping/normalization (optional; must be >= 0).',
        type=argparse_float)
    parser.add_argument(
        '--depth_max',
        help='Maximum depth for clipping/normalization (optional; must be >= 0).',
        type=argparse_float)
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
        '--depth_axis_labels',
        help=(
            'Depth coverage axis labels: range (one [min-max] label for the axis; default) or numbered '
            '(a number on each large tick).'
        ),
        choices=SCALAR_AXIS_LABEL_STYLES)
    parser.add_argument(
        '--depth_large_tick_interval',
        help='Depth coverage large tick interval in x coverage units (optional; must be > 0).',
        type=argparse_float)
    parser.add_argument(
        '--depth_small_tick_interval',
        help='Depth coverage small tick interval in x coverage units (optional; must be > 0; hidden by default).',
        type=argparse_float)
    parser.add_argument(
        '--depth_tick_font_size',
        help='Depth coverage tick label font size (optional; must be > 0).',
        type=argparse_float)


def _add_depth_color_arg(parser: argparse.ArgumentParser) -> None:
    """``--depth_color`` of either mode."""
    parser.add_argument(
        '--depth_color',
        help='Depth track fill color (optional; default: #4A90E2).',
        type=str)


def require_positive_args(parser: argparse.ArgumentParser, args: argparse.Namespace, *dests: str) -> None:
    """Refuse a value <= 0 of each option ``--<dest>`` that was given, in this order."""

    for dest in dests:
        value = getattr(args, dest)
        if value is not None and value <= 0:
            parser.error(f"--{dest} must be > 0")


def validate_depth_axis_args(
    parser: argparse.ArgumentParser, args: argparse.Namespace, *, track_heights: bool = False
) -> None:
    """Refuse a depth range or tick interval out of bounds, and a per-track depth number that is not a number
    > 0 or ``auto`` (``track_heights``: also Linear's ``--depth_track_height``)."""

    if args.depth_min is not None and args.depth_min < 0:
        parser.error("--depth_min must be >= 0")
    if args.depth_max is not None and args.depth_max < 0:
        parser.error("--depth_max must be >= 0")
    if args.depth_min is not None and args.depth_max is not None and args.depth_min > args.depth_max:
        parser.error("--depth_min must be <= --depth_max")
    require_positive_args(parser, args, "depth_large_tick_interval", "depth_small_tick_interval")
    per_track = (("depth_track_height",) if track_heights else ()) + _DEPTH_TRACK_NUMBER_OPTIONS
    for option_name in per_track:
        for option_value in getattr(args, option_name) or []:
            option_text = str(option_value).strip().lower()
            if option_text in {"", "auto", "none", "null", "-"}:
                continue
            numeric_option_value = parse_ascii_finite_float(str(option_value))
            if numeric_option_value is None:
                parser.error(f"--{option_name} values must be numbers or auto")
            if numeric_option_value <= 0:
                parser.error(f"--{option_name} values must be > 0")


def _add_comparison_filter_args(parser: argparse.ArgumentParser) -> None:
    thresholds = CIRCULAR_MODE_PROFILE.comparison
    parser.add_argument(
        '--evalue',
        help='Maximum BLAST e-value retained for similarity rings (default: 1e-5).',
        type=argparse_float,
        default=thresholds.evalue)
    parser.add_argument(
        '--bitscore',
        help='Minimum BLAST bitscore retained for similarity rings (default: 50).',
        type=argparse_float,
        default=thresholds.bitscore)
    parser.add_argument(
        '--identity',
        help='Minimum BLAST identity percentage retained for similarity rings (default: 70).',
        type=argparse_float,
        default=thresholds.identity)
    parser.add_argument(
        '--alignment_length',
        help='Minimum BLAST alignment length retained for similarity rings (default: 0).',
        type=argparse_integer,
        default=thresholds.alignment_length)
    add_identity_scale_arg(parser, "similarity rings")


def validate_comparison_threshold_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Refuse match thresholds that ``ComparisonThresholds`` refuses, in either mode."""

    try:
        ComparisonThresholds(
            evalue=args.evalue,
            bitscore=args.bitscore,
            identity=args.identity,
            alignment_length=args.alignment_length,
        )
    except ValidationError as exc:
        parser.error(str(exc))


class _IdentityScaleAction(argparse.Action):
    """Read ``--identity_scale MIN MAX`` with the scale's one validator; a refusal names the option."""

    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: object,
        option_string: str | None = None,
    ) -> None:
        try:
            scale = identity_scale_from_config(values, path=self.dest)
        except ValidationError:
            parser.error(f"argument {option_string}: must be {IDENTITY_SCALE_RULE}")
        setattr(namespace, self.dest, scale)


def add_identity_scale_arg(parser: argparse.ArgumentParser, drawn: str) -> None:
    """``--identity_scale MIN MAX`` of either mode; ``drawn`` names what the scale colours."""

    parser.add_argument(
        '--identity_scale',
        help=f'Fix the identity colour scale of the {drawn} and their Legend to MIN..MAX percent '
        '(0 <= MIN < MAX <= 100), so figures with different --identity thresholds share one scale '
        '(default: from --identity to 100).',
        nargs=2,
        metavar=('MIN', 'MAX'),
        type=argparse_float,
        action=_IdentityScaleAction)


def _add_gc_content_axis_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--gc_content_mode',
        help='GC content display mode: deviation (current centered behavior) or percent (absolute 0-100%% area track).',
        choices=['deviation', 'percent'],
        default=None)
    parser.add_argument(
        '--gc_content_min_percent',
        help='Minimum GC percent for percent-mode clipping/axis (optional; finite number).',
        type=argparse_float)
    parser.add_argument(
        '--gc_content_max_percent',
        help='Maximum GC percent for percent-mode clipping/axis (optional; finite number).',
        type=argparse_float)
    parser.add_argument(
        '--gc_content_tick_interval',
        help='GC content percent-mode large tick interval (optional; must be > 0; alias for --gc_content_large_tick_interval).',
        type=argparse_float)
    parser.add_argument(
        '--gc_content_large_tick_interval',
        help='GC content percent-mode large tick interval (optional; must be > 0).',
        type=argparse_float)
    parser.add_argument(
        '--gc_content_small_tick_interval',
        help='GC content percent-mode small tick interval (optional; must be > 0; hidden by default).',
        type=argparse_float)
    parser.add_argument(
        '--gc_content_tick_font_size',
        help='GC content percent-mode tick label font size (optional; must be > 0).',
        type=argparse_float)
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
        '--gc_content_axis_labels',
        help=(
            'GC content percent-mode axis labels: range (one [min-max] label for the axis; default) or numbered '
            '(a number on each large tick).'
        ),
        choices=SCALAR_AXIS_LABEL_STYLES)


def validate_gc_content_axis_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Refuse a GC content percent range whose minimum is above its maximum, and a tick interval <= 0."""

    if (
        args.gc_content_min_percent is not None
        and args.gc_content_max_percent is not None
        and args.gc_content_min_percent > args.gc_content_max_percent
    ):
        parser.error("--gc_content_min_percent must be <= --gc_content_max_percent")
    require_positive_args(
        parser,
        args,
        "gc_content_tick_interval",
        "gc_content_large_tick_interval",
        "gc_content_small_tick_interval",
    )


def validate_scale_interval_arg(
    parser: argparse.ArgumentParser, args: argparse.Namespace, *, legacy_session: bool
) -> None:
    """Refuse ``--scale_interval`` <= 0; the replay of a released legacy Session reads it as automatic."""

    if args.scale_interval is not None and args.scale_interval <= 0:
        if not legacy_session:
            parser.error("--scale_interval must be > 0")
        args.scale_interval = None  # D-04: a released legacy Session drew <= 0 as automatic.


def _add_legend_size_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--legend_rows',
        help=(
            'Legend row edits: a TSV with the header "key\tname\thidden\tposition", one line per Legend row by '
            'its key; name renames the row (blank: keep its name), hidden true hides it, position (1, 2, ...) on '
            'every shown row or on none orders the rows. A key the Legend does not draw is ignored.'
        ),
        type=read_legend_rows_table,
        metavar='FILE',
    )
    parser.add_argument(
        '--legend_box_size',
        help='Legend box size (optional; float; default: 24 (pixels, 96 dpi) for genomes shorter than the short/long boundary (75 kb circular, 100 kb linear), 20 otherwise).',
        type=argparse_float)
    parser.add_argument(
        '--legend_font_size',
        help='Legend font size (optional; float; default: 20 (pt) for genomes shorter than the short/long boundary (75 kb circular, 100 kb linear), 16 otherwise).',
        type=argparse_float)


def _add_font_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--font_family',
        help=(
            'Font for all diagram text: one font family name (for example "DejaVu Sans"), '
            'which gets its generic family appended, or a CSS font-family list. Text is '
            'measured with the font file of the first family that is a --font_file font, '
            'a bundled Liberation font or alias, that fontconfig (fc-match) finds, or, without '
            'fontconfig, that is installed in the system font folders. A family found nowhere '
            'is replaced with the nearest Liberation font, with a warning. A Session saves the '
            "font name, not the font file (optional; default: 'Liberation Sans', 'Arial', "
            "'Helvetica', 'Nimbus Sans L', sans-serif; with --font_file, the family of the "
            "first file)."
        ),
        type=str)
    add_font_file_args(parser)


def add_font_file_args(parser: argparse.ArgumentParser) -> None:
    """``--font_file`` and ``--no_font_embedding``, for the mode commands and
    ``gbdraw render`` (also with ``--session``)."""
    parser.add_argument(
        '--font_file',
        metavar='PATH',
        action='append',
        default=[],
        help=(
            'A TrueType or OpenType font file (.ttf, .otf, .ttc, .otc) whose family '
            'the diagram text can use (--font_family, or the font a Session names); '
            'every face of a collection is used. Repeatable. Works with --session. '
            'A Session does not store the file (optional).'
        ),
        type=str)
    parser.add_argument(
        '--no_font_embedding',
        dest='embed_fonts',
        action='store_false',
        help=(
            'Do not embed fonts in SVG and interactive SVG output. By default, each font '
            'the SVG uses that gbdraw does not bundle is embedded as a subset of the '
            'characters it draws (the whole font when its license forbids subsetting), '
            'unless its license (OS/2 fsType) forbids embedding.'
        ))


def register_font_files(paths: Sequence[str]) -> list[RegisteredFont]:
    """Register every face of each ``--font_file``, in order."""
    return [face for path in paths for face in register_font_file_faces(path)]


def cli_font_family(args: argparse.Namespace) -> str | None:
    """Register each ``--font_file`` and return the ``font_family`` override:
    ``--font_family`` completed with its generic family, else the family of
    the first file's first face, else None (the configured default)."""
    registered = register_font_files(args.font_file)
    if args.font_family:
        return complete_font_family(args.font_family)
    if registered:
        return font_family_list_for(registered[0].family, registered[0].generic)
    return None


# -t/--table of both modes: a color table file or a built-in table by name.
COLOR_TABLE_HELP = (
    "color table: a TSV file, or a built-in table "
    f"({', '.join(BUILTIN_COLOR_TABLE_PREFIX + name for name in BUILTIN_COLOR_TABLES)}) (optional)"
)


def add_color_args(parser: argparse.ArgumentParser) -> None:
    """Add color-related arguments (-p, -t, -d)."""
    _add_palette_arg(parser)
    _add_color_table_arg(parser)
    _add_default_colors_arg(parser)


def _add_palette_arg(parser: argparse.ArgumentParser) -> None:
    """``-p/--palette`` of either mode."""
    parser.add_argument(
        "-p", "--palette",
        metavar="PALETTE",
        default="default",
        help="Palette name (default: default)",
        type=str)


def _add_color_table_arg(parser: argparse.ArgumentParser) -> None:
    """``-t/--table`` of either mode."""
    parser.add_argument(
        '-t',
        '--table',
        help=COLOR_TABLE_HELP,
        type=str,
        default="")


def _add_default_colors_arg(parser: argparse.ArgumentParser) -> None:
    """``-d/--default_colors`` of either mode."""
    parser.add_argument(
        '-d',
        '--default_colors',
        help='TSV file that overrides the color palette (optional)',
        type=str,
        default="")


def add_analysis_args(parser: argparse.ArgumentParser, *, nt_subject: str = "dinucleotide") -> None:
    """Add GC/skew analysis arguments (-n, -w, -s); ``nt_subject`` begins the help of ``-n``."""
    parser.add_argument(
        '-n',
        '--nt',
        help=f'{nt_subject}: two letters from A, C, G, T, and U; case-insensitive, U counts as T (default: GC).',
        type=str,
        default="GC")
    _add_window_step_args(parser)


def add_feature_args(parser: argparse.ArgumentParser) -> None:
    """Add feature selection arguments (-k)."""
    parser.add_argument(
        '-k',
        '--features',
        help='Comma-separated list of feature keys to draw (default: CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,repeat_region)',
        type=str,
        default=",".join(DEFAULT_FEATURE_TYPES))


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
        '--italic_gene_labels',
        help='Draw each feature label that the qualifier priority names from the gene qualifier in italics, '
        'also when --label_table replaces its text; other labels stay upright (default: off).',
        action='store_true')
    parser.add_argument(
        '--feature_visibility_table',
        dest='feature_table',
        help='Path to a TSV file defining per-feature visibility overrides (optional)',
        type=str,
        default="")
    parser.add_argument(
        '--annotation_table',
        '--annotation-table',
        dest='annotation_table',
        help='Path to a TSV file defining region annotation sets and marks (optional)',
        type=str,
        default='')


def _add_label_rendering_args(parser: argparse.ArgumentParser, *, mode: LabelStyleMode) -> None:
    """``--label_rendering``, ``--label_line`` and ``--label_protection``; the help names the mode's rule and
    defaults."""
    parser.add_argument(
        '--label_rendering',
        help='Label rendering policy: "auto" embeds fitting labels and routes others externally; "embedded_only" '
        'drops external labels; "external_only" forces labels outside feature bodies. '
        f'{_LABEL_RENDERING_RULE_BY_MODE[mode]}Default: "auto".',
        choices=list(get_args(LabelRenderingPolicy)),
        default='auto',
        type=str)
    parser.add_argument(
        '--label_line',
        help='Line under each horizontal feature label with a leader: "under" draws one that continues the leader '
        '(a turned label gets none); "none" draws none. '
        f'Default: "{DEFAULT_LABEL_LINE[mode]}".',
        choices=LABEL_LINES,
        type=str)
    parser.add_argument(
        '--label_protection',
        help='What keeps a leader passing behind a feature label off its letters: "halo" (a white outline of the letters), '
        f'"box" (a white box), or "none". Default: "{DEFAULT_LABEL_PROTECTION[mode]}".',
        choices=LABEL_PROTECTIONS,
        type=str)


def _cli_list(values: Sequence[str] | None) -> list[str] | None:
    return list(values or []) or None


def depth_track_inputs_from_args(
    args: argparse.Namespace, *, heights: Sequence[str] | None = None
) -> tuple[DepthTrackInput, ...] | None:
    """The depth tracks of ``--depth_track`` and its per-track options, in either mode (``heights``: Linear's)."""

    return depth_track_inputs_from_cli(
        args.depth_track,
        labels=_cli_list(args.depth_track_label),
        colors=_cli_list(args.depth_track_color),
        heights=_cli_list(heights),
        large_tick_intervals=_cli_list(args.depth_track_large_tick_interval),
        small_tick_intervals=_cli_list(args.depth_track_small_tick_interval),
        tick_font_sizes=_cli_list(args.depth_track_tick_font_size),
    )


def _add_depth_track_source_args(parser: argparse.ArgumentParser, *, mode: str) -> None:
    _add_depth_track_arg(parser, mode=mode)
    _add_depth_track_label_color_args(parser)


def add_feature_density_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        '--feature_density_track',
        metavar='KEY=VALUE',
        help='A feature density track, as key=value words: id=<name> type=<feature type, or * for every type> '
        '[qualifier=<name> [value=<regex>]] [min_length=<bp>] [max_length=<bp>] [metric=count|coverage] '
        '[style=heatmap|histogram] [window=<bp>] [step=<bp>] [color=<color>] [caption=<text>]. '
        'count: the features in each window; coverage: the percentage of each window the features cover. '
        'Window and step default to --window and --step, else to the automatic window and step of the longest '
        'record in the figure. Repeat the option for more tracks.',
        type=str,
        nargs='+',
        action='append')
    parser.add_argument(
        '--feature_density_table',
        help='Path to a TSV file of feature density tracks: a header line naming the columns (id, type, qualifier, '
        'value, min_length, max_length, metric, style, window, step, color, caption) and one track per line '
        '(optional)',
        type=str,
        default='')


def feature_density_tracks_from_args(args: argparse.Namespace) -> tuple[FeatureDensityTrack, ...] | None:
    """The tracks of ``--feature_density_table`` and then of each ``--feature_density_track``, or None."""

    table: str = args.feature_density_table
    tracks = list(read_feature_density_table(table)) if table else []
    for index, words in enumerate(args.feature_density_track or (), start=1):
        tracks.append(parse_feature_density_words(words, where=f"--feature_density_track #{index}"))
    return unique_feature_density_tracks(tracks) or None


def validate_input_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Validate input file argument combinations."""
    records_table = getattr(args, "records_table", None)
    if records_table and (
        getattr(args, "record_topology", None) is not None
        or getattr(args, "display_start_coordinate", None) is not None
    ):
        parser.error("--records_table cannot be combined with --record_topology or --display_start_coordinate.")
    if records_table and (args.gbk or args.gff or args.fasta):
        parser.error("--records_table cannot be used with --gbk, --gff, or --fasta.")
    if args.gbk and (args.gff or args.fasta):
        parser.error("--gbk cannot be used with --gff or --fasta.")
    if args.gff and not args.fasta:
        parser.error("--gff requires --fasta.")
    if not records_table and not args.gbk and not args.fasta:
        parser.error("Either --records_table, --gbk, or --fasta (with or without --gff) must be provided.")


def apply_record_display_cli_options(
    manifest: RecordInputManifest, args: argparse.Namespace,
) -> RecordInputManifest:
    """Adapt direct flags once, leaving biological validation to the planner."""
    topology = getattr(args, "record_topology", None)
    start = getattr(args, "display_start_coordinate", None)
    if topology is None and start is None:
        return manifest
    if args.records_table or len(manifest.records) != 1:
        raise ValidationError(
            "Display flags require exactly one direct-input record; use --records_table."
        )
    display = RecordDisplayOptions(
        is_circular=None if topology in {None, "auto"} else topology == "circular",
        start_coordinate=start,
    )
    return replace(manifest, records=(replace(
        manifest.records[0], display=display, cardinality=RecordCardinality.EXACTLY_ONE,
    ),))


def validate_output_choice_args(args: argparse.Namespace) -> None:
    """Refuse an unknown ``-f/--format`` token or ``-p/--palette`` name (OV-1244).

    The replay of a released legacy Session skips this check: its saved
    arguments keep the fallbacks they were drawn with.
    """

    parse_output_formats(args.format)
    require_palette_name(args.palette)


def validate_label_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    """Validate label argument combinations."""
    if args.label_whitelist and args.label_blacklist:
        parser.error("--label_whitelist and --label_blacklist are mutually exclusive.")


def parse_feature_shape_assignment_arg(value: str) -> str:
    try:
        parse_feature_shape_assignment(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return value


def parse_arrow_head_length_ratio_arg(value: str) -> str | float:
    text = str(value).strip().lower()
    if text == "auto":
        return "auto"
    parsed = parse_ascii_finite_float(text)
    if parsed is None or parsed <= 0:
        raise argparse.ArgumentTypeError(
            "arrow head length ratio must be 'auto' or a positive finite number"
        )
    return parsed


def parse_arrow_shaft_width_ratio_arg(value: str) -> float:
    parsed = parse_ascii_finite_float(str(value))
    if parsed is None or parsed <= 0 or parsed > 1:
        raise argparse.ArgumentTypeError(
            "arrow shaft width ratio must be a finite number in (0, 1]"
        )
    return parsed


def handle_output_formats(out_formats: list[str]) -> list[str]:
    """The output formats this run writes, checked before anything is drawn.

    Under WebAssembly the browser converts PNG, PDF, EPS and PS, so they are
    left out. Elsewhere a format CairoSVG cannot make raises
    ``ValidationError`` (``require_cairosvg``), as the Python API does.
    """
    out_formats = [normalize_format_token(fmt) for fmt in out_formats]
    if "pyodide" in sys.modules:
        if any(is_cairosvg_format(f) for f in out_formats):
            logger.info("Running in WebAssembly mode: Image conversion is handled by the browser.")
            return [f for f in out_formats if not is_cairosvg_format(f)] or [SVG_FORMAT]
        return out_formats
    require_cairosvg(out_formats)
    return out_formats


__all__ = [
    "add_analysis_args",
    "add_feature_density_args",
    "depth_track_inputs_from_args",
    "feature_density_tracks_from_args",
    "COLOR_TABLE_HELP",
    "add_color_args",
    "add_feature_args",
    "add_input_args",
    "add_label_args",
    "cli_label",
    "handle_output_formats",
    "setup_logging",
    "validate_input_args",
    "validate_label_args",
]
