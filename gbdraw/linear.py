#!/usr/bin/env python
# coding: utf-8


import argparse
import logging
import sys
from Bio.SeqRecord import SeqRecord
from typing import Any, Optional
from pandas import DataFrame
from .file_processing import load_default_colors, load_gbks, read_color_table, load_config_toml, parse_formats, read_qualifier_priority_file
from .linear_diagram_components import plot_linear_diagram
from .utility_functions import create_dict_for_sequence_lengths, modify_config_dict
from .canvas_generator import LinearCanvasConfigurator
from .object_configurators import GcSkewConfigurator, LegendDrawingConfigurator, GcContentConfigurator, FeatureDrawingConfigurator, BlastMatchConfigurator

# Setup for the logging system. Configures a stream handler to direct log messages to stdout
# and sets the logging level to INFO for both the handler and the logger.
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.INFO)
logger.addHandler(handler)
logger.setLevel(logging.INFO)


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
        '-i',
        '--input',
        help='genbank (required)',
        type=str,
        required=True,
        nargs='*')
    parser.add_argument(
        '-b',
        '--blast',
        help="input BLAST result file in tab-separated format (-outfmt 6 or 7) (optional)",
        type=str,
        nargs='*')
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
        '--separate_strands',
        help='separate forward and reverse strands (default: False). Features of undefined strands are shown on the forward strand.',
        action='store_true')
    parser.add_argument(
        '--show_gc',
        help='plot GC content below genome (default: False). ',
        action='store_true')
    parser.add_argument(
        '--align_center',
        help='Align genomes to the center (default: False). ',
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
        '-k',
        '--features',
        help='Comma-separated list of feature keys to draw (default: CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,repeat_region)',
        type=str,
        default="CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,repeat_region")
    parser.add_argument(
        '--block_stroke_color',
        help='Block stroke color (str; default: "black")',
        type=str,
        default="black")
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
        '--label_font_size',
        help='Label font size (optional; default: 16 for short genomes, 5 for long genomes)',
        type=float)
    parser.add_argument(
        '-f',
        '--format',
        help='Comma-separated list of output file formats (default: png)',
        type=str,
        default="png")
    parser.add_argument(
        '-l',
        '--legend',
        help='Legend position (default: "right"; "right", "left", "none")',
        type=str,
        default="right")
    parser.add_argument("--show_labels", help="Show labels", action="store_true")
    parser.add_argument(
        '--resolve_overlaps',
        help='Resolve overlaps (experimental; default: False). ',
        action='store_true')
    # NEW: ラベルブラックリスト用のコマンドラインオプションを追加
    parser.add_argument(
        '--label_blacklist',
        help='Comma-separated keywords or path to a file for label blacklisting (optional)',
        type=str,
        default="")
        
    # NEW: ラベル優先順位ファイル用のコマンドラインオプションを追加
    parser.add_argument(
        '--qualifier_priority',
        help='Path to a TSV file defining qualifier priority for labels (optional)',
        type=str,
        default="")
    args = parser.parse_args(args)
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
    args: argparse.Namespace = _get_args(cmd_args)
    genbank_files: str = args.input
    out_file_prefix: str = args.output
    blast_files: str = args.blast
    color_table_path: str = args.table
    strandedness: bool = args.separate_strands
    resolve_overlaps: bool = args.resolve_overlaps
    dinucleotide: str = args.nt
    show_gc: bool = args.show_gc
    window: int = args.window
    step: int = args.step
    align_center: bool = args.align_center
    evalue: float = args.evalue
    legend: str = args.legend
    show_skew: bool = False
    bitscore: float = args.bitscore
    identity: float = args.identity
    show_labels: bool = args.show_labels
    # NEW: 新しいコマンドライン引数をローカル変数に格納
    label_blacklist: str = args.label_blacklist
    qualifier_priority_path: str = args.qualifier_priority
    selected_features_set: str = args.features.split(',')

    out_formats: list[str] = parse_formats(args.format)
    user_defined_default_colors: str = args.default_colors
    if blast_files:
        load_comparison = True
    else:
        load_comparison = False
    palette: str = args.palette
    default_colors: Optional[DataFrame] = load_default_colors(
        user_defined_default_colors, palette)
    color_table: Optional[DataFrame] = read_color_table(color_table_path)
    config_dict: dict = load_config_toml('gbdraw.data', 'config.toml')


    if qualifier_priority_path:
        qualifier_priority_df = read_qualifier_priority_file(qualifier_priority_path)
        config_dict['labels']['filtering']['qualifier_priority_df'] = qualifier_priority_df
    else:
        config_dict['labels']['filtering']['qualifier_priority_df'] = None


    block_stroke_color: str = args.block_stroke_color
    block_stroke_width: str = args.block_stroke_width
    line_stroke_color: str = args.line_stroke_color
    line_stroke_width: str = args.line_stroke_width       
    config_dict = modify_config_dict(
        config_dict, 
        block_stroke_color=block_stroke_color, 
        block_stroke_width=block_stroke_width, 
        line_stroke_color=line_stroke_color, 
        line_stroke_width=line_stroke_width, 
        show_gc=show_gc, 
        show_skew=show_skew, 
        align_center=align_center, 
        strandedness=strandedness,
        label_blacklist=label_blacklist
    )
    records: list[SeqRecord] = load_gbks(
        genbank_files, "linear", load_comparison)
    sequence_length_dict: dict[str,
                               int] = create_dict_for_sequence_lengths(records)
    longest_genome: int = max(sequence_length_dict.values())
    num_of_entries: int = len(sequence_length_dict)
    config_dict = modify_config_dict(config_dict, block_stroke_color=block_stroke_color, block_stroke_width=block_stroke_width, line_stroke_color=line_stroke_color, line_stroke_width=line_stroke_width, show_gc=show_gc, show_skew=show_skew, align_center=align_center, strandedness=strandedness, show_labels=show_labels, resolve_overlaps=resolve_overlaps)

    blast_config = BlastMatchConfigurator(
        evalue=evalue, bitscore=bitscore, identity=identity, sequence_length_dict=sequence_length_dict, config_dict=config_dict, default_colors_df=default_colors)
    canvas_config = LinearCanvasConfigurator(output_prefix=out_file_prefix,
                                            num_of_entries=num_of_entries, longest_genome=longest_genome, config_dict=config_dict, legend=legend)
    feature_config = FeatureDrawingConfigurator(
        color_table=color_table, default_colors=default_colors, selected_features_set=selected_features_set, config_dict=config_dict)
    gc_config = GcContentConfigurator(
        window=window, step=step, dinucleotide=dinucleotide, config_dict=config_dict, default_colors_df=default_colors)
    skew_config = GcSkewConfigurator(
        window=window, step=step, dinucleotide=dinucleotide, config_dict=config_dict, default_colors_df=default_colors)
    legend_config = LegendDrawingConfigurator(color_table=color_table, default_colors=default_colors, selected_features_set=selected_features_set, config_dict=config_dict, gc_config=gc_config, skew_config=skew_config, feature_config=feature_config)
    plot_linear_diagram(records, blast_files, canvas_config, blast_config,
                        feature_config, gc_config, config_dict, out_formats, legend_config, skew_config)


if __name__ == "__main__":
    # This gets all arguments passed to the script, excluding the script name
    main_args = sys.argv[1:]
    if not main_args:
        main_args.append('--help')
    linear_main(main_args)
