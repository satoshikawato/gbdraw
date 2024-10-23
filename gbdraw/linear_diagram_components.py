#!/usr/bin/env python
# coding: utf-8

import logging
import sys
from pandas import DataFrame
from typing import Tuple
from svgwrite import Drawing
from svgwrite.container import Group
from Bio.SeqRecord import SeqRecord
from .canvas_generator import LinearCanvasConfigurator
from .circular_object_groups import LegendGroup
from .linear_object_groups import SeqRecordGroup, DefinitionGroup, GcContentGroup, PairWiseMatchGroup, LengthBarGroup
from .file_processing import load_comparisons, save_figure
from .object_configurators import LegendDrawingConfigurator, GcContentConfigurator, FeatureDrawingConfigurator
from .data_processing import prepare_legend_table
from .utility_functions import check_feature_presence

# Logging setup
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)


def calculate_record_offsets(count: int, record: SeqRecord, canvas_config: LinearCanvasConfigurator) -> Tuple[float, float]:
    """
    Calculates the vertical and horizontal offsets for placing a genomic record in the canvas.

    Args:
        count (int): The sequential number of the record being processed.
        record (SeqRecord): The genomic record to be plotted.
        canvas_config (LinearCanvasConfigurator): Configuration object for the canvas.

    Returns:
        Tuple[float, float]: A tuple containing the vertical (offset) and horizontal (offset_x) offsets.

    The function computes the vertical offset based on the count of the record and various configuration Configurator.
    The horizontal offset is calculated to align the record in the center, considering the length of the record and the longest genome length.
    """
    offset: float = canvas_config.vertical_offset + (canvas_config.cds_height + canvas_config.vertical_padding +
                                                     canvas_config.gc_padding + canvas_config.comparison_height + canvas_config.vertical_padding) * (count - 1)
    offset_x: float = (canvas_config.alignment_width *
                       ((canvas_config.longest_genome - len(record.seq)) / canvas_config.longest_genome) / 2) if canvas_config.align_center else 0
    return offset, offset_x


def position_record_group(record_group: Group, offset: float, offset_x: float, canvas_config: LinearCanvasConfigurator) -> Group:
    """
    Translates the record group to its designated position on the canvas.

    Args:
        record_group (Group): The SVG group element containing the record visualization.
        offset (float): Vertical offset for the position.
        offset_x (float): Horizontal offset for the position.
        canvas_config (LinearCanvasConfigurator): Configuration object for the canvas.

    Returns:
        Group: The translated SVG group element.

    This function positions the record group on the canvas based on the calculated offsets.
    """
    record_group.translate(offset_x + canvas_config.horizontal_offset, offset)
    return record_group


def position_gc_content_group(gc_content_group: Group, offset: float, offset_x: float, canvas_config: LinearCanvasConfigurator) -> Group:
    """
    Positions the GC content group on the canvas.

    Args:
        gc_content_group (Group): The SVG group element containing the GC content visualization.
        offset (float): Vertical offset for the position.
        offset_x (float): Horizontal offset for the position.
        canvas_config (LinearCanvasConfigurator): Configuration object for the canvas.

    Returns:
        Group: The translated SVG group element.

    This function adjusts the position of the GC content visualization relative to the genomic record it corresponds to.
    """
    gc_content_group.translate(
        offset_x + canvas_config.horizontal_offset,
        offset + 2 * canvas_config.cds_padding + canvas_config.vertical_padding)
    return gc_content_group


def position_record_definition_group(record_definition_group: Group, offset: float, offset_x: float, canvas_config: LinearCanvasConfigurator) -> Group:
    """
    Places the record definition group in the correct position on the canvas.

    Args:
        record_definition_group (Group): The SVG group element containing the record's definition.
        offset (float): Vertical offset for the position.
        offset_x (float): Horizontal offset for the position.
        canvas_config (LinearCanvasConfigurator): Configuration object for the canvas.

    Returns:
        Group: The adjusted SVG group element.

    This function aligns the record's definition with the rest of the record's visualization on the canvas.
    """
    record_definition_group.translate(
        offset_x + canvas_config.horizontal_offset - 45, offset)
    return record_definition_group


def position_length_bar_group(total_height: float, vertical_offset: float) -> float:
    """
    Calculates the vertical position for the length bar on the canvas.

    Args:
        total_height (float): The total height of the canvas.
        vertical_offset (float): The vertical offset from the top of the canvas.

    Returns:
        float: The vertical position for placing the length bar.
    """
    return (int(float(total_height)) - 0.5 * vertical_offset)


def position_comparison_group(comparison_count: int, canvas_config: LinearCanvasConfigurator) -> float:
    """
    Calculates the vertical position for placing a comparison group on the canvas.

    Args:
        comparison_count (int): The sequential number of the comparison being processed.
        canvas_config (LinearCanvasConfigurator): Configuration object for the canvas.

    Returns:
        float: The vertical position for placing the comparison group.
    """
    return ((canvas_config.vertical_offset + 0.5 * canvas_config.cds_height +
             canvas_config.vertical_padding + 0.9 * canvas_config.gc_padding) +
            ((canvas_config.comparison_height + canvas_config.vertical_padding +
              canvas_config.cds_height + canvas_config.vertical_padding +
              canvas_config.gc_padding) * (comparison_count - 1)))


def add_record_group(canvas: Drawing, record: SeqRecord, offset: float, offset_x: float, canvas_config: LinearCanvasConfigurator, feature_config: FeatureDrawingConfigurator, config_dict: dict) -> Drawing:
    """
    Adds a record group to the linear canvas.

    This function creates a group for a SeqRecord visualization and positions it on the canvas based on the given offsets.

    Args:
        canvas (Drawing): The SVG drawing canvas.
        record (SeqRecord): The GenBank record to be visualized.
        offset (float): Vertical offset for the group's position.
        offset_x (float): Horizontal offset for the group's position.
        canvas_config (LinearCanvasConfigurator): Configuration for the linear canvas.
        feature_config (FeatureDrawingConfigurator): Configuration for feature drawing.
        config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
        Drawing: The updated SVG drawing with the record group added.
    """
    # Unpack parameters
    record_group: Group = SeqRecordGroup(
        gb_record=record, canvas_config=canvas_config, feature_config=feature_config, config_dict=config_dict).get_group()
    position_record_group(record_group, offset, offset_x, canvas_config)
    canvas.add(record_group)
    return canvas


def add_gc_content_group(canvas: Drawing, record: SeqRecord, offset: float, offset_x: float, canvas_config: LinearCanvasConfigurator, gc_config: GcContentConfigurator, config_dict: dict) -> Drawing:
    """
    Adds a GC content group to the linear canvas.

    This function creates a group for visualizing the GC content of a SeqRecord and positions it on the canvas based on the given offsets.

    Args:
        canvas (Drawing): The SVG drawing canvas.
        record (SeqRecord): The GenBank record with GC content data.
        offset (float): Vertical offset for the group's position.
        offset_x (float): Horizontal offset for the group's position.
        canvas_config (LinearCanvasConfigurator): Configuration for the linear canvas.
        gc_config (GcContentConfigurator): Configuration for the GC content representation.
        config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
        Drawing: The updated SVG drawing with the GC content group added.
    """
    gc_content_group: Group = GcContentGroup(gb_record=record, alignment_width=canvas_config.alignment_width, longest_record_len=canvas_config.longest_genome,
                                             track_height=canvas_config.cds_height, gc_config=gc_config, config_dict=config_dict).get_group()
    position_gc_content_group(
        gc_content_group, offset, offset_x, canvas_config)
    canvas.add(gc_content_group)
    return canvas


def add_record_definition_group(canvas: Drawing, record: SeqRecord, offset: float, offset_x: float, canvas_config: LinearCanvasConfigurator, config_dict: dict) -> Drawing:
    """
    Adds a record definition group to the linear canvas.

    This function creates a group that contains definition details of the SeqRecord, such as annotations, and positions it on the canvas.

    Args:
        canvas (Drawing): The SVG drawing canvas.
        record (SeqRecord): The GenBank record with definition details.
        offset (float): Vertical offset for the group's position.
        offset_x (float): Horizontal offset for the group's position.
        canvas_config (LinearCanvasConfigurator): Configuration for the linear canvas.
        config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
        Drawing: The updated SVG drawing with the record definition group added.
    """
    record_definition_group: Group = DefinitionGroup(
        record, config_dict).get_group()
    position_record_definition_group(
        record_definition_group, offset, offset_x, canvas_config)
    canvas.add(record_definition_group)
    return canvas


def add_comparison_on_linear_canvas(canvas: Drawing, comparisons, canvas_config: LinearCanvasConfigurator, blast_config, config_dict: dict) -> Drawing:
    """
    Adds comparison groups, such as pairwise matches, to the linear canvas.

    This function iterates over comparisons and creates groups for each pairwise match.
    These groups are positioned on the canvas based on the configuration settings.

    Args:
        canvas (Drawing): The SVG drawing canvas for adding the comparisons.
        comparisons: A collection of comparison data (e.g., pairwise matches).
        canvas_config (LinearCanvasConfigurator): Configuration settings for the linear canvas.
        blast_config: Configuration settings for BLAST comparisons.
        config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
        Drawing: The updated SVG drawing with the comparison groups added.
    """
    for comparison_count, comparison in enumerate(comparisons, start=1):
        match_group: Group = PairWiseMatchGroup(canvas_config, blast_config.sequence_length_dict,
                                                comparison, canvas_config.comparison_height, comparison_count, blast_config).get_group()
        offset: float = position_comparison_group(
            comparison_count, canvas_config)
        match_group.translate(canvas_config.horizontal_offset, offset)
        canvas.add(match_group)
    return canvas


def add_length_bar_on_linear_canvas(canvas: Drawing, canvas_config: LinearCanvasConfigurator, config_dict: dict) -> Drawing:
    """
    Adds a length bar to the linear canvas.

    This function creates and positions a length bar on the canvas, which serves as a scale
    reference for the genomic data. The length bar is adjusted according to the canvas configuration.

    Args:
        canvas (Drawing): The SVG drawing canvas for adding the length bar.
        canvas_config (LinearCanvasConfigurator): Configuration settings for the linear canvas.
        config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
        Drawing: The updated SVG drawing with the length bar added.
    """
    length_bar_group: Group = LengthBarGroup(
        canvas_config.fig_width, canvas_config.longest_genome, config_dict).get_group()
    offset_for_length_bar: float = position_length_bar_group(
        canvas_config.total_height, canvas_config.vertical_offset)
    length_bar_group.translate(
        canvas_config.horizontal_offset, offset_for_length_bar)
    canvas.add(length_bar_group)
    return canvas

def add_legends_on_linear_canvas(canvas: Drawing, canvas_config: LinearCanvasConfigurator, legend_config, legend_table):
    legend_group: Group = LegendGroup(canvas_config, legend_config, legend_table).get_group()
    offset_x = canvas_config.legend_offset_x
    offset_y = canvas_config.legend_offset_y
    legend_group.translate(offset_x, offset_y) 
    canvas.add(legend_group)
    return canvas

def add_records_on_linear_canvas(canvas: Drawing, records: list[SeqRecord], feature_config: FeatureDrawingConfigurator, gc_config: GcContentConfigurator, canvas_config: LinearCanvasConfigurator, config_dict: dict) -> Drawing:
    """
    Adds multiple SeqRecord groups to the linear canvas.

    This function iterates over a list of SeqRecords and adds each as a group to the canvas.
    It includes options to add associated features like GC content and record definitions.

    Args:
        canvas (Drawing): The SVG drawing canvas for adding the SeqRecords.
        records (list[SeqRecord]): A list of SeqRecords to be visualized.
        feature_config (FeatureDrawingConfigurator): Configuration for feature drawing.
        gc_config (GcContentConfigurator): Configuration for GC content representation.
        canvas_config (LinearCanvasConfigurator): Configuration settings for the linear canvas.
        config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
        Drawing: The updated SVG drawing with the SeqRecord groups added.
    """
    for count, record in enumerate(records, start=1):
        offset: float
        offset_x: float
        offset, offset_x = calculate_record_offsets(
            count, record, canvas_config)
        # Add record group
        add_record_group(canvas, record, offset, offset_x,
                         canvas_config, feature_config, config_dict)
        # Add record definition group
        add_record_definition_group(
            canvas, record, offset, offset_x, canvas_config, config_dict)
        # Add GC content group if configured to show
        if canvas_config.show_gc:
            add_gc_content_group(canvas, record, offset,
                                 offset_x, canvas_config, gc_config, config_dict)

    return canvas


def plot_linear_diagram(records: list[SeqRecord], blast_files, canvas_config: LinearCanvasConfigurator, blast_config, feature_config: FeatureDrawingConfigurator, gc_config: GcContentConfigurator, config_dict: dict, out_formats, legend_config, skew_config) -> None:
    """
    Plots a linear diagram of genomic records with optional BLAST comparison data.

    Args:
        records (list[SeqRecord]): List of genomic records to be plotted.
        blast_files: List of paths to BLAST output files for comparison.
        canvas_config (LinearCanvasConfigurator): Configuration object for the canvas.
        blast_config: Configuration for BLAST comparison visualization.
        feature_config (FeatureDrawingConfigurator): Configurator for feature drawing.
        gc_config (GcContentConfigurator): Configurator for GC content visualization.
        config_dict (dict): General configuration dictionary.
        out_formats: List of output formats for saving the diagram.

    This function orchestrates the creation of a linear diagram, including genomic records, GC content, BLAST comparisons, and a length bar. It handles the layout and positioning of these elements and saves the final diagram in specified formats.
    """
    # Initialize and create canvas



    features_present = check_feature_presence(records, feature_config.selected_features_set)
    legend_table = prepare_legend_table(gc_config, skew_config, feature_config, features_present)
    legend_config = legend_config.recalculate_legend_dimensions(legend_table)
    canvas_config.recalculate_canvas_dimensions(legend_config)
    canvas: Drawing = canvas_config.create_svg_canvas()


    # Add length bar
    canvas = add_legends_on_linear_canvas(canvas, canvas_config, legend_config, legend_table)
    canvas = add_length_bar_on_linear_canvas(
        canvas, canvas_config, config_dict)
    # Add BLAST pairwise matches (if specified)
    if blast_files is not None:
        comparisons: list[DataFrame] = load_comparisons(
            blast_files, blast_config)  # file_processing
        canvas = add_comparison_on_linear_canvas(
            canvas, comparisons, canvas_config, blast_config, config_dict)
    # Add records
    canvas = add_records_on_linear_canvas(
        canvas, records, feature_config, gc_config, canvas_config, config_dict)
    # Create SVG file
    save_figure(canvas, out_formats)
