#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
import logging
import math
from typing import Optional
from pandas import DataFrame  # type: ignore[reportMissingImports]
from .io.genome import load_gbks, load_gff_fasta
from .io.colors import load_default_colors, read_color_table
from .config.toml import load_config_toml
from .render.export import parse_formats, save_figure
from .api.diagram import (  # type: ignore[reportMissingImports]
    assemble_circular_diagram_from_record,
    assemble_circular_diagram_from_records,
)
from .config.modify import suppress_gc_content_and_skew, modify_config_dict  # type: ignore[reportMissingImports]
from .config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from .core.sequence import determine_output_file_prefix  # type: ignore[reportMissingImports]
from .labels.filtering import (
    read_filter_list_file,
    read_label_override_file,
    read_qualifier_priority_file,
)  # type: ignore[reportMissingImports]
from .io.record_select import parse_record_selector
from .features.shapes import parse_feature_shape_assignment, parse_feature_shape_overrides
from .features.visibility import read_feature_visibility_file
from .exceptions import ValidationError

from .cli_utils.common import (
    CAIROSVG_AVAILABLE,
    setup_logging,
    validate_input_args,
    validate_label_args,
    handle_output_formats,
    calculate_window_step,
)

# Setup for the logging system
logger = logging.getLogger()
setup_logging()


def _parse_feature_shape_assignment_arg(value: str) -> str:
    try:
        parse_feature_shape_assignment(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return value


def _parse_multi_record_position_arg(value: str) -> str:
    raw = str(value or "").strip()
    if not raw:
        raise argparse.ArgumentTypeError(
            "multi_record_position does not allow empty entries."
        )
    if "@" not in raw:
        raise argparse.ArgumentTypeError(
            f"multi_record_position entry '{raw}' must be in '<selector>@<row>' format."
        )
    selector_text, row_text = raw.rsplit("@", 1)
    selector_text = selector_text.strip()
    row_text = row_text.strip()
    if not selector_text:
        raise argparse.ArgumentTypeError(
            f"multi_record_position entry '{raw}' must include a selector before '@'."
        )
    if not row_text.isdigit() or int(row_text) <= 0:
        raise argparse.ArgumentTypeError(
            f"multi_record_position entry '{raw}' must use a positive integer row."
        )
    try:
        selector = parse_record_selector(selector_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    if selector is None:
        raise argparse.ArgumentTypeError(
            f"multi_record_position selector '{selector_text}' is invalid."
        )
    return f"{selector_text}@{int(row_text)}"


def _get_args(args) -> argparse.Namespace:
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
    parser = argparse.ArgumentParser(
        description=(
            "Generate genome diagrams in PNG/PDF/SVG/PS/EPS. "
            "By default, diagrams for multiple entries are saved separately. "
            "Use --multi_record_canvas to place multiple records on one grid canvas."
        )
    )
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
        '-o',
        '--output',
        help='output file prefix (default: accession number of the sequence)',
        type=str)
    parser.add_argument(
        "-p", "--palette",
        metavar="PALETTE",
        default="default",
        help="Palette name (default: default)",
        type=str
    )
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
    parser.add_argument(
        '--species',
        help='Species name (optional; e.g. "<i>Escherichia coli</i>", "<i>Ca.</i> Hepatoplasma crinochetorum")',
        type=str)
    parser.add_argument(
        '--strain',
        help='Strain/isolate name (optional; e.g. "K-12", "Av")',
        type=str)
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
        help='Axis stroke color (str; default: "gray")',
        type=str)
    parser.add_argument(
        '--axis_stroke_width',
        help='Axis stroke width (optional; float; default: 3 pt for genomes <= 50 kb, 1 pt for genomes >= 50 kb)',
        type=float)
    parser.add_argument(
        '--line_stroke_color',
        help='Line stroke color (str; default: "gray")',
        type=str)
    parser.add_argument(
        '--line_stroke_width',
        help='Line stroke width (optional; float; default: 5 pt for genomes <= 50 kb, 1 pt for genomes >= 50 kb)',
        type=float)
    parser.add_argument(
        '--definition_font_size',
        help='Definition font size (optional; default: 18)',
        type=float)
    parser.add_argument(
        '--plot_title',
        help='Circular plot title shown when plot_title_position is top or bottom (optional; defaults to species + strain).',
        type=str)
    parser.add_argument(
        '--plot_title_font_size',
        help='Plot title font size for circular top/bottom title layout (optional; default: 32).',
        type=float)
    parser.add_argument(
        '--label_font_size',
        help='Label font size (optional; default: 14 (pt) for genomes <= 50 kb, 8 for genomes >= 50 kb)',
        type=float)
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
    parser.add_argument(
        '--suppress_gc',
        help='Suppress GC content track (default: False).',
        action='store_true')
    parser.add_argument(
        '--suppress_skew',
        help='Suppress GC skew track (default: False).',
        action='store_true')
    parser.add_argument(
        '-l',
        '--legend',
        help='Legend position (default: "right"; "left", "right", "top", "bottom", "upper_left", "upper_right", "lower_left", "lower_right", "none")',
        type=str,
        default="right")
    parser.add_argument(
        '--multi_record_canvas',
        help='Place multiple records on one shared canvas using automatic grid layout (default: False).',
        action='store_true')
    parser.add_argument(
        '--multi_record_size_mode',
        help='Size mode for multi-record circular canvas ("auto", "linear", "equal"; "sqrt" is accepted as an alias of "auto"; default: "auto").',
        type=str,
        choices=['auto', 'linear', 'equal', 'sqrt'],
        default='auto')
    parser.add_argument(
        '--multi_record_min_radius_ratio',
        help='Minimum radius ratio for multi-record scaling (0 < ratio <= 1; default: 0.55).',
        type=float,
        default=0.55)
    parser.add_argument(
        '--multi_record_column_gap_ratio',
        help='Additional horizontal gap ratio between visible content bounds in each multi-record row (>= 0; default: 0.10).',
        type=float,
        default=0.10)
    parser.add_argument(
        '--multi_record_row_gap_ratio',
        help='Additional gap ratio between multi-record row content bounds (>= 0; default: 0.05).',
        type=float,
        default=0.05)
    parser.add_argument(
        '--multi_record_position',
        help="Record placement for multi-record canvas (repeatable): <selector>@<row> where selector is #index or record_id and row starts at 1.",
        type=_parse_multi_record_position_arg,
        action='append',
        default=[])
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
        help='Track type (default: "tuckin"; "tuckin", "middle", "spreadout")',
        type=str,
        default="tuckin")
    parser.add_argument(
        '--resolve_overlaps',
        help='Resolve overlapping features by placing them on separate tracks (default: False). Useful for plasmid visualization.',
        action='store_true')
    parser.add_argument(
        '--labels',
        help='Label placement mode: no argument or "out" (outside), "both" (outside+inside), or "none" (hidden). Default: "none".',
        nargs='?',
        const='out',
        default='none',
        choices=['none', 'out', 'both'],
        type=str)
    
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
        '--outer_label_x_radius_offset',
        help='Outer label x-radius offset factor (float; default from config)',
        type=float)
    parser.add_argument(
        '--outer_label_y_radius_offset',
        help='Outer label y-radius offset factor (float; default from config)',
        type=float)
    parser.add_argument(
        '--inner_label_x_radius_offset',
        help='Inner label x-radius offset factor (float; default from config)',
        type=float)
    parser.add_argument(
        '--inner_label_y_radius_offset',
        help='Inner label y-radius offset factor (float; default from config)',
        type=float)
    parser.add_argument(
        '--scale_interval',
        help='Manual scale interval for circular mode (in bp). Overrides automatic calculation.',
        type=int)
    parser.add_argument(
        '--feature_width',
        help='Feature track width for circular mode (in px; must be > 0).',
        type=float)
    parser.add_argument(
        '--gc_content_width',
        help='GC content track width for circular mode (in px; must be > 0).',
        type=float)
    parser.add_argument(
        '--gc_content_radius',
        help='GC content track center radius for circular mode (as a ratio of base radius; must be > 0).',
        type=float)
    parser.add_argument(
        '--gc_skew_width',
        help='GC skew track width for circular mode (in px; must be > 0).',
        type=float)
    parser.add_argument(
        '--gc_skew_radius',
        help='GC skew track center radius for circular mode (as a ratio of base radius; must be > 0).',
        type=float)
    parser.add_argument(
        '--legend_box_size',
        help='Legend box size (optional; float; default: 24 (pixels, 96 dpi) for genomes <= 50 kb, 20 for genomes >= 50 kb).',
        type=float)
    parser.add_argument(
        '--legend_font_size',
        help='Legend font size (optional; float; default: 20 (pt) for genomes <= 50 kb, 16 for genomes >= 50 kb).',
        type=float)
    
    args = parser.parse_args(args)
    validate_input_args(parser, args)
    validate_label_args(parser, args)
    if args.feature_width is not None and args.feature_width <= 0:
        parser.error("--feature_width must be > 0")
    if args.gc_content_width is not None and args.gc_content_width <= 0:
        parser.error("--gc_content_width must be > 0")
    if args.gc_content_radius is not None and args.gc_content_radius <= 0:
        parser.error("--gc_content_radius must be > 0")
    if args.gc_skew_width is not None and args.gc_skew_width <= 0:
        parser.error("--gc_skew_width must be > 0")
    if args.gc_skew_radius is not None and args.gc_skew_radius <= 0:
        parser.error("--gc_skew_radius must be > 0")
    if args.multi_record_min_radius_ratio <= 0 or args.multi_record_min_radius_ratio > 1:
        parser.error("--multi_record_min_radius_ratio must be > 0 and <= 1")
    if not math.isfinite(args.multi_record_column_gap_ratio) or args.multi_record_column_gap_ratio < 0:
        parser.error("--multi_record_column_gap_ratio must be a finite number >= 0")
    if not math.isfinite(args.multi_record_row_gap_ratio) or args.multi_record_row_gap_ratio < 0:
        parser.error("--multi_record_row_gap_ratio must be a finite number >= 0")
    return args


    
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
    args: argparse.Namespace = _get_args(cmd_args)   
    output_prefix = args.output
    dinucleotide: str = args.nt.upper()
    manual_window: int = args.window
    manual_step: int = args.step
    color_table_path: str = args.table
    selected_features_set: str = args.features.split(',')
    feature_shapes = parse_feature_shape_overrides(args.feature_shape)
    species: str = args.species
    strain: str = args.strain
    legend: str = args.legend
    multi_record_canvas: bool = args.multi_record_canvas
    multi_record_size_mode: str = args.multi_record_size_mode
    multi_record_min_radius_ratio: float = args.multi_record_min_radius_ratio
    multi_record_column_gap_ratio: float = args.multi_record_column_gap_ratio
    multi_record_row_gap_ratio: float = args.multi_record_row_gap_ratio
    multi_record_positions: list[str] = [str(position) for position in (args.multi_record_position or [])]
    plot_title: str = str(args.plot_title or "").strip()
    plot_title_position: str = args.plot_title_position
    definition_font_size: Optional[float] = args.definition_font_size
    plot_title_font_size: Optional[float] = args.plot_title_font_size
    label_font_size: Optional[float] = args.label_font_size
    suppress_gc: bool = args.suppress_gc
    suppress_skew: bool = args.suppress_skew
    labels_mode: str = args.labels
    show_labels: bool = labels_mode != "none"
    resolve_overlaps: bool = args.resolve_overlaps
    allow_inner_labels: bool = labels_mode == "both"
    label_whitelist: str = args.label_whitelist
    label_blacklist: str = args.label_blacklist
    qualifier_priority_path: str = args.qualifier_priority
    label_table_path: str = args.label_table
    feature_table_path: str = args.feature_table
    scale_interval: Optional[int] = args.scale_interval
    legend_box_size = args.legend_box_size
    legend_font_size = args.legend_font_size
    feature_width: Optional[float] = args.feature_width
    gc_content_width: Optional[float] = args.gc_content_width
    gc_content_radius: Optional[float] = args.gc_content_radius
    gc_skew_width: Optional[float] = args.gc_skew_width
    gc_skew_radius: Optional[float] = args.gc_skew_radius
    if plot_title_font_size is not None and float(plot_title_font_size) <= 0:
        raise ValidationError("plot_title_font_size must be > 0")
    if args.gbk:
        gb_records = load_gbks(args.gbk, "circular")
    elif args.gff and args.fasta:
        gb_records = load_gff_fasta(
            args.gff,
            args.fasta,
            "circular",
            selected_features_set,
            keep_all_features=bool(feature_table_path),
        )
    else:
        # This case should not be reached due to arg validation
        logger.error("Invalid input file configuration.")
        raise ValidationError("Invalid input file configuration.")

    outer_label_x_radius_offset: Optional[float] = args.outer_label_x_radius_offset
    outer_label_y_radius_offset: Optional[float] = args.outer_label_y_radius_offset
    inner_label_x_radius_offset: Optional[float] = args.inner_label_x_radius_offset
    inner_label_y_radius_offset: Optional[float] = args.inner_label_y_radius_offset
    if allow_inner_labels and not (suppress_gc and suppress_skew):

        suppress_gc = True 
        suppress_skew = True
        logger.warning(
            "WARNING: --labels both requires suppressing GC and skew tracks. Suppressing GC and skew tracks.")  # 

    user_defined_default_colors: str = args.default_colors
    block_stroke_color: Optional[str] = args.block_stroke_color
    block_stroke_width: Optional[float] = args.block_stroke_width
    axis_stroke_color: Optional[str] = args.axis_stroke_color
    axis_stroke_width: Optional[float] = args.axis_stroke_width
    line_stroke_color: Optional[str] = args.line_stroke_color
    line_stroke_width: Optional[float] = args.line_stroke_width   
    track_type: str = args.track_type
    strandedness = args.separate_strands
    scale_interval: Optional[int] = args.scale_interval
    if (
        not multi_record_canvas
        and bool(multi_record_positions)
    ):
        logger.info(
            "Ignoring --multi_record_position because --multi_record_canvas is disabled."
        )
    
    # Warn if resolve_overlaps is used with separate_strands
    if strandedness and resolve_overlaps:
        logger.warning(
            "WARNING: --resolve_overlaps is ignored when --separate_strands is enabled.")
        resolve_overlaps = False
    
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

    palette: str = args.palette
    default_colors: Optional[DataFrame] = load_default_colors(
        user_defined_default_colors, palette)
    
    color_table: Optional[DataFrame] = read_color_table(color_table_path)
    feature_table: Optional[DataFrame] = read_feature_visibility_file(feature_table_path)
    show_gc, show_skew = suppress_gc_content_and_skew(
        suppress_gc, suppress_skew)

    config_dict = modify_config_dict(
        config_dict, 
        block_stroke_color=block_stroke_color, 
        block_stroke_width=block_stroke_width,
        circular_axis_stroke_color=axis_stroke_color, 
        circular_axis_stroke_width=axis_stroke_width, 
        line_stroke_color=line_stroke_color, 
        line_stroke_width=line_stroke_width, 
        show_labels=show_labels, 
        track_type=track_type, 
        strandedness=strandedness, 
        resolve_overlaps=resolve_overlaps,
        show_gc=show_gc, 
        show_skew=show_skew, 
        allow_inner_labels=allow_inner_labels,
        circular_definition_font_size=definition_font_size,
        plot_title_font_size=plot_title_font_size,
        label_font_size=label_font_size,
        label_blacklist=label_blacklist,
        label_whitelist=label_whitelist,
        label_table=label_table_path,
        outer_label_x_radius_offset=outer_label_x_radius_offset,
        outer_label_y_radius_offset=outer_label_y_radius_offset,
        inner_label_x_radius_offset=inner_label_x_radius_offset,
        inner_label_y_radius_offset=inner_label_y_radius_offset,
        scale_interval=scale_interval,
        legend_box_size=legend_box_size,
        legend_font_size=legend_font_size
    )    

    out_formats: list[str] = parse_formats(args.format)
    out_formats = handle_output_formats(out_formats)
    record_count: int = 0



    cfg = GbdrawConfig.from_dict(config_dict)
    track_specs: list[str] = []
    if feature_width is not None:
        track_specs.append(f"features@w={float(feature_width):g}px")

    gc_content_spec_requested = (gc_content_width is not None) or (gc_content_radius is not None)
    if gc_content_spec_requested:
        if not show_gc:
            logger.warning(
                "WARNING: GC content track is suppressed. Ignoring --gc_content_width/--gc_content_radius."
            )
        else:
            gc_content_opts: list[str] = []
            if gc_content_radius is not None:
                gc_content_opts.append(f"r={float(gc_content_radius):g}")
            if gc_content_width is not None:
                gc_content_opts.append(f"w={float(gc_content_width):g}px")
            if gc_content_opts:
                track_specs.append(f"gc_content@{','.join(gc_content_opts)}")

    gc_skew_spec_requested = (gc_skew_width is not None) or (gc_skew_radius is not None)
    if gc_skew_spec_requested:
        if not show_skew:
            logger.warning(
                "WARNING: GC skew track is suppressed. Ignoring --gc_skew_width/--gc_skew_radius."
            )
        else:
            gc_skew_opts: list[str] = []
            if gc_skew_radius is not None:
                gc_skew_opts.append(f"r={float(gc_skew_radius):g}")
            if gc_skew_width is not None:
                gc_skew_opts.append(f"w={float(gc_skew_width):g}px")
            if gc_skew_opts:
                track_specs.append(f"gc_skew@{','.join(gc_skew_opts)}")

    track_specs_or_none = track_specs or None

    if multi_record_canvas and len(gb_records) > 1:
        first_accession = gb_records[0].id if gb_records else "out"
        outfile_prefix = output_prefix if output_prefix is not None else first_accession
        canvas = assemble_circular_diagram_from_records(
            gb_records,
            config_dict=config_dict,
            color_table=color_table,
            default_colors=default_colors,
            selected_features_set=selected_features_set,
            feature_table=feature_table,
            feature_shapes=feature_shapes or None,
            output_prefix=outfile_prefix,
            legend=legend,
            dinucleotide=dinucleotide,
            window=manual_window,
            step=manual_step,
            species=species,
            strain=strain,
            plot_title=plot_title,
            plot_title_position=plot_title_position,
            plot_title_font_size=plot_title_font_size,
            multi_record_size_mode=multi_record_size_mode,
            multi_record_min_radius_ratio=multi_record_min_radius_ratio,
            multi_record_column_gap_ratio=multi_record_column_gap_ratio,
            multi_record_row_gap_ratio=multi_record_row_gap_ratio,
            multi_record_positions=multi_record_positions or None,
            cfg=cfg,
            track_specs=track_specs_or_none,
        )
        save_figure(canvas, out_formats)
    else:
        for gb_record in gb_records:
            record_count += 1
            accession = gb_record.id
            seq_length = len(gb_record.seq)
            window, step = calculate_window_step(seq_length, cfg, manual_window, manual_step)

            outfile_prefix = determine_output_file_prefix(gb_records, output_prefix, record_count, accession)
            canvas = assemble_circular_diagram_from_record(
                gb_record,
                config_dict=config_dict,
                color_table=color_table,
                default_colors=default_colors,
                selected_features_set=selected_features_set,
                feature_table=feature_table,
                feature_shapes=feature_shapes or None,
                output_prefix=outfile_prefix,
                legend=legend,
                dinucleotide=dinucleotide,
                window=window,
                step=step,
                species=species,
                strain=strain,
                plot_title=plot_title,
                plot_title_position=plot_title_position,
                plot_title_font_size=plot_title_font_size,
                cfg=cfg,
                track_specs=track_specs_or_none,
            )
            save_figure(canvas, out_formats)

if __name__ == "__main__":
    # Entry point for the script when run as a standalone program.
    # Calls the main function `circular_main` with command-line arguments.
    # args are parsed command-line arguments
    # This gets all arguments passed to the script, excluding the script name
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(handler)
    main_args = sys.argv[1:]
    if not main_args:
        main_args.append('--help')
    circular_main(main_args)
