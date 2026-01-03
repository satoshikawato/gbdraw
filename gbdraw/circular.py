#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
import logging
from typing import Optional
from pandas import DataFrame  # type: ignore[reportMissingImports]
from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from .io.genome import load_gbks, load_gff_fasta
from .io.colors import load_default_colors, read_color_table
from .config.toml import load_config_toml
from .render.export import parse_formats, save_figure
from .api.diagram import assemble_circular_diagram_from_record  # type: ignore[reportMissingImports]
from .config.modify import suppress_gc_content_and_skew, modify_config_dict  # type: ignore[reportMissingImports]
from .config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from .core.sequence import determine_output_file_prefix  # type: ignore[reportMissingImports]
from .labels.filtering import read_qualifier_priority_file, read_filter_list_file  # type: ignore[reportMissingImports]

try:
    import cairosvg  # type: ignore[reportMissingImports]
    CAIROSVG_AVAILABLE = True
except (ImportError, OSError):
    CAIROSVG_AVAILABLE = False

# Setup for the logging system and sets the logging level to INFO
logger = logging.getLogger()
logger.setLevel(logging.INFO)
# Ensure handler is added if not already present
if not logger.handlers:
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(handler)


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
        description='Generate genome diagrams in PNG/PDF/SVG/PS/EPS. Diagrams for multiple entries are saved separately.')
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
        help='Legend position (default: "right"; "left", "right", "upper_left", "upper_right", "lower_left", "lower_right", "none")',
        type=str,
        default="right")
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
        '--show_labels',
        help='Show feature labels (default: False).',
        action='store_true')
    parser.add_argument(
        '--allow_inner_labels',
        help='Place labels inside the circle (default: False). If enabled, labels are placed both inside and outside the circle, and gc and skew tracks are not shown.',
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
        '--legend_box_size',
        help='Legend box size (optional; float; default: 24 (pixels, 96 dpi) for genomes <= 50 kb, 20 for genomes >= 50 kb).',
        type=float)
    parser.add_argument(
        '--legend_font_size',
        help='Legend font size (optional; float; default: 20 (pt) for genomes <= 50 kb, 16 for genomes >= 50 kb).',
        type=float)
    
    args = parser.parse_args(args)
    if args.gbk and (args.gff or args.fasta):
        parser.error("Error: --gbk cannot be used with --gff or --fasta.")
    
    # Check if both --gff and --fasta are provided together
    if args.gff and not args.fasta:
        parser.error("Error: --gff requires --fasta.")
    
    if args.fasta and not args.gff:
        parser.error("Error: --fasta requires --gff.")
        
    # Ensure that either --gbk or both --gff and --fasta are provided
    if not args.gbk and not (args.gff and args.fasta):
        parser.error("Error: Either --gbk or both --gff and --fasta must be provided.")
    if args.label_whitelist and args.label_blacklist:
        parser.error("Error: --label_whitelist and --label_blacklist are mutually exclusive.")
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
    species: str = args.species
    strain: str = args.strain
    legend: str = args.legend
    definition_font_size: Optional[float] = args.definition_font_size
    label_font_size: Optional[float] = args.label_font_size
    suppress_gc: bool = args.suppress_gc
    suppress_skew: bool = args.suppress_skew
    show_labels: bool = args.show_labels
    resolve_overlaps: bool = args.resolve_overlaps
    allow_inner_labels: bool = args.allow_inner_labels
    label_whitelist: str = args.label_whitelist
    label_blacklist: str = args.label_blacklist
    qualifier_priority_path: str = args.qualifier_priority
    scale_interval: Optional[int] = args.scale_interval
    legend_box_size = args.legend_box_size
    legend_font_size = args.legend_font_size
    if args.gbk:
        gb_records = load_gbks(args.gbk, "circular")
    elif args.gff and args.fasta:
        gb_records = load_gff_fasta(args.gff, args.fasta, "circular", selected_features_set)
    else:
        # This case should not be reached due to arg validation
        logger.error("Invalid input file configuration.")
        sys.exit(1)

    outer_label_x_radius_offset: Optional[float] = args.outer_label_x_radius_offset
    outer_label_y_radius_offset: Optional[float] = args.outer_label_y_radius_offset
    inner_label_x_radius_offset: Optional[float] = args.inner_label_x_radius_offset
    inner_label_y_radius_offset: Optional[float] = args.inner_label_y_radius_offset
    if allow_inner_labels and not show_labels:
        show_labels = True  # If inner labels are allowed, labels must be shown
        logger.warning(
            "WARNING: Inner labels are allowed, but labels are not shown. Enabling labels.")
    if allow_inner_labels and not (suppress_gc and suppress_skew):

        suppress_gc = True 
        suppress_skew = True
        logger.warning(
            "WARNING: Inner labels are allowed, but GC and skew tracks are not suppressed. Suppressing GC and skew tracks.")  # 

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

    palette: str = args.palette
    default_colors: Optional[DataFrame] = load_default_colors(
        user_defined_default_colors, palette)
    
    color_table: Optional[DataFrame] = read_color_table(color_table_path)
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
        label_font_size=label_font_size,
        label_blacklist=label_blacklist,
        label_whitelist=label_whitelist,
        outer_label_x_radius_offset=outer_label_x_radius_offset,
        outer_label_y_radius_offset=outer_label_y_radius_offset,
        inner_label_x_radius_offset=inner_label_x_radius_offset,
        inner_label_y_radius_offset=inner_label_y_radius_offset,
        scale_interval=scale_interval,
        legend_box_size=legend_box_size,
        legend_font_size=legend_font_size
    )    

    out_formats: list[str] = parse_formats(args.format)
    # Handle WebAssembly environment and CairoSVG availability
    if "pyodide" in sys.modules:
        if any(f != 'svg' for f in out_formats):
            logger.info("Running in WebAssembly mode: Output format constrained to SVG. (Image conversion is handled by the browser)")
            out_formats = ['svg']
    # Handle absence of CairoSVG
    elif not CAIROSVG_AVAILABLE:
        non_svg_formats = [f for f in out_formats if f != 'svg']
        if non_svg_formats:
            logger.warning(
                f"⚠️  CairoSVG is not installed. Cannot generate: {', '.join(non_svg_formats).upper()}\n"
                f"   Output restricted to SVG only.\n"
                f"   (To enable PNG/PDF, run: pip install gbdraw[export])"
            )
            out_formats = ['svg']
    record_count: int = 0



    cfg = GbdrawConfig.from_dict(config_dict)

    for gb_record in gb_records:
        record_count += 1 
        accession = gb_record.id
        seq_length = len(gb_record.seq)
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

        outfile_prefix = determine_output_file_prefix(gb_records, output_prefix, record_count, accession)
        canvas = assemble_circular_diagram_from_record(
            gb_record,
            config_dict=config_dict,
            color_table=color_table,
            default_colors=default_colors,
            selected_features_set=selected_features_set,
            output_prefix=outfile_prefix,
            legend=legend,
            dinucleotide=dinucleotide,
            window=window,
            step=step,
            species=species,
            strain=strain,
            cfg=cfg,
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