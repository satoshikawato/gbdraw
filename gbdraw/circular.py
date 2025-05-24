#!/usr/bin/env python
# coding: utf-8

import sys
import argparse
import logging
from typing import Optional
from pandas import DataFrame
from Bio.SeqRecord import SeqRecord
from .file_processing import load_gbks, load_default_colors, read_color_table, load_config_toml, parse_formats
from .circular_diagram_components import plot_circular_diagram
from .data_processing import skew_df
from .canvas_generator import CircularCanvasConfigurator
from .object_configurators import LegendDrawingConfigurator, GcSkewConfigurator,GcContentConfigurator, FeatureDrawingConfigurator
from .utility_functions import suppress_gc_content_and_skew, modify_config_dict, determine_output_file_prefix

# Setup the logging system. Configures a stream handler to output log messages to stdout.
# Default logging level is set to INFO.
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)


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
        description='Generate genome diagrams in PNG/PDF/SVG/PS/EPS. Diagrams for multiple entries are saved separately (hence the lack of output file name option).')
    parser.add_argument(
        '-i',
        '--input',
        help='Genbank/DDBJ flatfile (required)',
        type=str,
        required=True,
        nargs='*')
    parser.add_argument(
        '-o',
        '--output',
        help='output file prefix (default: accession number of the sequence)',
        type=str)
    parser.add_argument(
        "-p", "--palette",
        metavar="NAME",
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
        help='TSV file that specifies default color Configurator (optional; default: data/default_colors.tsv)',
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
        help='window size (default: 1000) ',
        type=int,
        default="1000")
    parser.add_argument(
        '-s',
        '--step',
        help='step size (default: 100) ',
        type=int,
        default="100")
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
        help='Comma-separated list of feature keys to draw (default: CDS,tRNA,rRNA,repeat_region)',
        type=str,
        default="CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,repeat_region,regulatory,rep_origin")
    parser.add_argument(
        '--block_stroke_color',
        help='Block stroke color (str; default: "gray")',
        type=str,
        default="gray")
    parser.add_argument(
        '--block_stroke_width',
        help='Block stroke width (float; default: 0)',
        type=float,
        default=0)
    parser.add_argument(
        '--line_stroke_color',
        help='Line stroke color (str; default: "gray")',
        type=str,
        default="gray")
    parser.add_argument(
        '--line_stroke_width',
        help='Line stroke width (float; default: 1.0)',
        type=float,
        default=1.0)
    parser.add_argument(
        '-f',
        '--format',
        help='Comma-separated list of output file formats (default: png)',
        type=str,
        default="png")
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
        help='Legend position (default: "right"; "left", "right", "upper_left", "upper_right", "lower_left", "lower_right")',
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
        '--show_labels',
        help='Show feature labels (default: False).',
        action='store_true')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args(args)
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
    input_file: str = args.input
    output_prefix = args.output
    dinucleotide: str = args.nt
    window: int = args.window
    step: int = args.step
    color_table_path: str = args.table
    selected_features_set: str = args.features.split(',')
    species: str = args.species
    strain: str = args.strain
    legend: str = args.legend
    suppress_gc: bool = args.suppress_gc
    suppress_skew: bool = args.suppress_skew
    show_labels: bool = args.show_labels
    user_defined_default_colors: str = args.default_colors
    block_stroke_color: str = args.block_stroke_color
    block_stroke_width: str = args.block_stroke_width
    line_stroke_color: str = args.line_stroke_color
    line_stroke_width: str = args.line_stroke_width   
    track_type: str = args.track_type
    strandedness = args.separate_strands
    config_dict: dict = load_config_toml('gbdraw.data', 'config.toml')
    palette: str = args.palette
    default_colors: Optional[DataFrame] = load_default_colors(
        user_defined_default_colors, palette)
    color_table: Optional[DataFrame] = read_color_table(color_table_path)
    show_gc, show_skew = suppress_gc_content_and_skew(
        suppress_gc, suppress_skew)
    config_dict = modify_config_dict(config_dict, block_stroke_color=block_stroke_color, block_stroke_width=block_stroke_width, line_stroke_color=line_stroke_color, line_stroke_width=line_stroke_width, show_labels=show_labels, track_type=track_type, strandedness=strandedness, show_gc=show_gc, show_skew=show_skew)
    out_formats: list[str] = parse_formats(args.format)
    record_count: int = 0
    gb_records: list[SeqRecord] = load_gbks(input_file, "circular")
    gc_config = GcContentConfigurator(
        window=window, step=step, dinucleotide=dinucleotide, config_dict=config_dict, default_colors_df=default_colors)
    skew_config = GcSkewConfigurator(
        window=window, step=step, dinucleotide=dinucleotide, config_dict=config_dict, default_colors_df=default_colors)
    feature_config = FeatureDrawingConfigurator(
        color_table=color_table, default_colors=default_colors, selected_features_set=selected_features_set, config_dict=config_dict)
    legend_config = LegendDrawingConfigurator(color_table=color_table, default_colors=default_colors, selected_features_set=selected_features_set, config_dict=config_dict, gc_config=gc_config, skew_config=skew_config, feature_config=feature_config)


    for gb_record in gb_records:
        record_count += 1 
        accession = str(gb_record.annotations["accessions"][0])  # type: ignore
        outfile_prefix = determine_output_file_prefix(gb_records, output_prefix, record_count, accession)
        gc_df: DataFrame = skew_df(gb_record, window, step, dinucleotide)
        canvas_config = CircularCanvasConfigurator(
            output_prefix=outfile_prefix, config_dict=config_dict, legend=legend, gb_record=gb_record
            )
        plot_circular_diagram(gb_record, canvas_config, gc_df, gc_config, skew_config,
                              feature_config, species, strain, config_dict, out_formats, legend_config)

if __name__ == "__main__":
    # Entry point for the script when run as a standalone program.
    # Calls the main function `circular_main` with command-line arguments.
    # args are parsed command-line arguments
    # This gets all arguments passed to the script, excluding the script name
    args = sys.argv[1:]
    circular_main(args)
