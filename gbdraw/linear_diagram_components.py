#!/usr/bin/env python
# coding: utf-8

import logging
import sys
from pandas import DataFrame
from typing import Tuple, Optional
from svgwrite import Drawing
from svgwrite.container import Group
from Bio.SeqRecord import SeqRecord
from .canvas_generator import LinearCanvasConfigurator
from .linear_object_groups import SeqRecordGroup, DefinitionGroup, GcContentGroup, GcSkewGroup, PairWiseMatchGroup, LengthBarGroup, LegendGroup
from .file_processing import load_comparisons, save_figure
from .object_configurators import LegendDrawingConfigurator, GcContentConfigurator, GcSkewConfigurator, FeatureDrawingConfigurator
from .data_processing import prepare_legend_table
from .utility_functions import check_feature_presence
from .create_feature_objects import create_feature_dict, preprocess_color_tables
from .data_processing import prepare_label_list_linear
from .utility_functions import preprocess_label_filtering

# Logging setup
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)


def _precalculate_definition_widths(records: list[SeqRecord], config_dict: dict) -> float:
    """
    Pre-calculates the maximum definition width among all records.
    """
    max_definition_width = 0
    if not records:
        return 0

    for record in records:
        def_group = DefinitionGroup(record, config_dict)
        if def_group.definition_bounding_box_width > max_definition_width:
            max_definition_width = def_group.definition_bounding_box_width
    return max_definition_width

def _precalculate_label_dimensions(records: list[SeqRecord], feature_config: FeatureDrawingConfigurator, canvas_config: LinearCanvasConfigurator, config_dict: dict) -> tuple[float, dict, dict]:
    """Pre-calculates label placements for all records to determine the required canvas height."""
    if not canvas_config.show_labels:
        return 0, {}, {}

    max_required_height = 0
    all_labels_by_record = {}
    record_label_heights = {} # Store height required for labels per record

    for record in records:
        color_table, default_colors = preprocess_color_tables(
            feature_config.color_table, feature_config.default_colors
        )
        label_filtering = preprocess_label_filtering(config_dict['labels']['filtering'])
        feature_dict = create_feature_dict(
            record, color_table, feature_config.selected_features_set, default_colors,
            canvas_config.strandedness, canvas_config.resolve_overlaps, label_filtering
        )

        record_length = len(record.seq)
        genome_size_normalization_factor = record_length / canvas_config.longest_genome

        label_list = prepare_label_list_linear(
            feature_dict, record_length, canvas_config.alignment_width,
            genome_size_normalization_factor, canvas_config.cds_height,
            canvas_config.strandedness, config_dict
        )
        all_labels_by_record[record.id] = label_list

        min_y_coord_record = 0
        for label in label_list:
            if not label['is_embedded']:
                label_bottom_y = label['middle_y'] - (label['height_px'] / 2)
                if label_bottom_y < min_y_coord_record:
                    min_y_coord_record = label_bottom_y
        
        record_height = abs(min_y_coord_record)
        record_label_heights[record.id] = record_height

        if record_height > max_required_height:
            max_required_height = record_height

    return max_required_height, all_labels_by_record, record_label_heights

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
                                                     canvas_config.gc_padding + canvas_config.skew_padding + canvas_config.comparison_height + canvas_config.vertical_padding) * (count - 1)
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


def position_gc_content_group(gc_content_group: Group, offset_y: float, offset_x: float, canvas_config: LinearCanvasConfigurator) -> Group:
    """
    Positions the GC content group on the canvas.

    Args:
        gc_content_group (Group): The SVG group element containing the GC content visualization.
        offset_y (float): Vertical offset for the position.
        offset_x (float): Horizontal offset for the position.
        canvas_config (LinearCanvasConfigurator): Configuration object for the canvas.

    Returns:
        Group: The translated SVG group element.

    This function adjusts the position of the GC content visualization relative to the genomic record it corresponds to.
    """
    gc_content_group.translate(
        offset_x + canvas_config.horizontal_offset,
        offset_y + canvas_config.cds_padding + canvas_config.vertical_padding)
    return gc_content_group

def position_gc_skew_group(gc_skew_group: Group, offset_y: float, offset_x: float, canvas_config: LinearCanvasConfigurator) -> Group:
    """
    Positions the GC skew group on the canvas.
    """
    y_offset = offset_y + canvas_config.cds_padding + canvas_config.vertical_padding
    if canvas_config.show_gc:
        y_offset += canvas_config.gc_height
        
    gc_skew_group.translate(
        offset_x + canvas_config.horizontal_offset,
        y_offset
    )
    return gc_skew_group


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
        canvas_config.horizontal_offset - offset_x, offset)
    return record_definition_group


def position_length_bar_group(total_height: float, vertical_offset: float, vertical_padding: float) -> float:
    """
    Calculates the vertical position for the length bar on the canvas.
    """
    return total_height - (vertical_offset -  vertical_padding)

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
             canvas_config.vertical_padding + 0.9 * canvas_config.gc_padding + 0.9 * canvas_config.skew_padding) +
            ((canvas_config.comparison_height + canvas_config.vertical_padding +
              canvas_config.cds_height + canvas_config.vertical_padding +
              canvas_config.gc_padding + canvas_config.skew_padding) * (comparison_count - 1)))


def add_record_group(canvas: Drawing, record: SeqRecord, offset_y: float, offset_x: float, canvas_config: LinearCanvasConfigurator, feature_config: FeatureDrawingConfigurator, config_dict: dict, precalculated_labels: Optional[list]) -> Drawing:
    """Adds a record group to the linear canvas."""
    record_group: Group = SeqRecordGroup(
        gb_record=record, canvas_config=canvas_config, feature_config=feature_config, config_dict=config_dict, precalculated_labels=precalculated_labels
    ).get_group()
    position_record_group(record_group, offset_y, offset_x, canvas_config)
    canvas.add(record_group)
    return canvas


def add_gc_content_group(canvas: Drawing, record: SeqRecord, offset_y: float, offset_x: float, canvas_config: LinearCanvasConfigurator, gc_config: GcContentConfigurator, config_dict: dict) -> Drawing:
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
                                             track_height=canvas_config.gc_height, gc_config=gc_config, config_dict=config_dict).get_group()
    position_gc_content_group(
        gc_content_group, offset_y, offset_x, canvas_config)
    canvas.add(gc_content_group)
    return canvas

def add_gc_skew_group(canvas: Drawing, record: SeqRecord, offset: float, offset_x: float, canvas_config: LinearCanvasConfigurator, skew_config: GcSkewConfigurator, config_dict: dict) -> Drawing:
    """
    Adds a GC skew group to the linear canvas.

    This function creates a group for visualizing the GC skew of a SeqRecord and positions it on the canvas based on the given offsets.

    Args:
        canvas (Drawing): The SVG drawing canvas.
        record (SeqRecord): The GenBank record with GC skew data.
        offset (float): Vertical offset for the group's position.
        offset_x (float): Horizontal offset for the group's position.
        canvas_config (LinearCanvasConfigurator): Configuration for the linear canvas.
        skew_config (GcContentConfigurator): Configuration for the GC skew representation.
        config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
        Drawing: The updated SVG drawing with the GC skew group added.
    """
    gc_skew_group: Group = GcSkewGroup(gb_record=record, alignment_width=canvas_config.alignment_width, longest_record_len=canvas_config.longest_genome,
                                             track_height=canvas_config.skew_height, skew_config=skew_config, config_dict=config_dict).get_group()
    position_gc_skew_group(
        gc_skew_group, offset, offset_x, canvas_config)
    canvas.add(gc_skew_group)
    return canvas

def add_record_definition_group(canvas: Drawing, record: SeqRecord, record_offset_y: float, record_offset_x: float, canvas_config: LinearCanvasConfigurator, config_dict: dict) -> Drawing:
    """
    Adds a record definition group to the linear canvas.
    """
    definition_group_obj = DefinitionGroup(record, config_dict)

    definition_offset_x = definition_group_obj.definition_bounding_box_width / 2

    record_definition_group: Group = definition_group_obj.get_group()
    
    position_record_definition_group(
        record_definition_group, record_offset_y, (definition_offset_x - record_offset_x), canvas_config)
    
    canvas.add(record_definition_group)
    return canvas


def add_comparison_on_linear_canvas(canvas: Drawing, comparisons, canvas_config: LinearCanvasConfigurator, blast_config, config_dict: dict, records: list, comparison_offsets: list, actual_comparison_heights: list) -> Drawing:
    """
    Adds comparison groups at specified y-offsets with dynamic height.
    """
    for comparison_count, comparison in enumerate(comparisons, start=1):
        if comparison_count > len(comparison_offsets):
            break
        
        height = actual_comparison_heights[comparison_count - 1]
        offset = comparison_offsets[comparison_count - 1]

        match_group: Group = PairWiseMatchGroup(
            canvas_config, blast_config.sequence_length_dict,
            comparison, height, comparison_count, blast_config, records
        ).get_group()
        
        match_group.translate(canvas_config.horizontal_offset, offset)
        canvas.add(match_group)
    return canvas

def add_length_bar_on_linear_canvas(canvas: Drawing, canvas_config: LinearCanvasConfigurator, config_dict: dict) -> Drawing:
    """
    Adds a length bar to the linear canvas.
    (snip)
    """
    length_bar_group: Group = LengthBarGroup(
        canvas_config.fig_width, canvas_config.alignment_width, canvas_config.longest_genome, config_dict).get_group()
    
    offset_for_length_bar: float = position_length_bar_group(
        canvas_config.total_height, canvas_config.original_vertical_offset, canvas_config.vertical_padding)
    length_bar_group.translate(
        canvas_config.horizontal_offset, offset_for_length_bar)
    canvas.add(length_bar_group)
    return canvas

def add_legends_on_linear_canvas(canvas: Drawing, config_dict, canvas_config: LinearCanvasConfigurator, legend_config, legend_table):
    legend_group: Group = LegendGroup(config_dict, canvas_config, legend_config, legend_table).get_group()
    offset_x = canvas_config.legend_offset_x
    offset_y = canvas_config.legend_offset_y
    legend_group.translate(offset_x, offset_y) 
    canvas.add(legend_group)
    return canvas


def plot_linear_diagram(records: list[SeqRecord], blast_files, canvas_config: LinearCanvasConfigurator, blast_config, feature_config: FeatureDrawingConfigurator, gc_config: GcContentConfigurator, config_dict: dict, out_formats, legend_config, skew_config) -> None:
    """
    Plots a linear diagram of genomic records with optional BLAST comparison data.
    """
    max_def_width = _precalculate_definition_widths(records, config_dict)
    required_label_height, all_labels, record_label_heights = _precalculate_label_dimensions(
        records, feature_config, canvas_config, config_dict
    )
    
    if required_label_height > 0:
        if canvas_config.vertical_offset < required_label_height:
            canvas_config.vertical_offset = required_label_height
    else:
        canvas_config.vertical_offset = canvas_config.original_vertical_offset + canvas_config.cds_padding
    has_blast = bool(blast_files)
    record_ids = [r.id for r in records]
    record_offsets = []
    
    current_y = canvas_config.vertical_offset 
    for i, record_id in enumerate(record_ids):
        record_offsets.append(current_y)
        
        if i < len(record_ids) - 1:
            height_below_axis = canvas_config.cds_padding + canvas_config.gc_padding + canvas_config.skew_padding
            next_record_id = record_ids[i+1]
            next_label_height = record_label_heights.get(next_record_id, 0)
            if next_label_height > canvas_config.comparison_height:
                inter_record_space = height_below_axis + next_label_height + canvas_config.cds_padding
            else:
                inter_record_space = height_below_axis + canvas_config.comparison_height + canvas_config.cds_padding
            current_y += inter_record_space
        

    final_height = current_y + canvas_config.cds_padding + canvas_config.gc_padding + canvas_config.skew_padding + canvas_config.original_vertical_offset + canvas_config.vertical_padding
    canvas_config.total_height = int(final_height)

    features_present = check_feature_presence(records, feature_config.selected_features_set)
    legend_table = prepare_legend_table(gc_config, skew_config, feature_config, features_present, blast_config, has_blast)
    legend_config = legend_config.recalculate_legend_dimensions(legend_table)
    padding = canvas_config.canvas_padding * 2  
    required_legend_height = legend_config.legend_height + padding

    if required_legend_height > canvas_config.total_height:
        height_difference = required_legend_height - canvas_config.total_height
        canvas_config.total_height = int(required_legend_height)
        vertical_shift = height_difference / 2
        record_offsets = [offset + vertical_shift for offset in record_offsets]

    canvas_config.recalculate_canvas_dimensions(legend_config, max_def_width)
    canvas: Drawing = canvas_config.create_svg_canvas()

    if canvas_config.legend_position != 'none':
        canvas = add_legends_on_linear_canvas(canvas, config_dict, canvas_config, legend_config, legend_table)
    canvas = add_length_bar_on_linear_canvas(canvas, canvas_config, config_dict)
    
    if blast_files:
        comparisons = load_comparisons(blast_files, blast_config)
        comparison_offsets = []
        actual_comparison_heights = []
        for i in range(len(records) - 1):
            height_below_axis = canvas_config.cds_padding + canvas_config.gc_padding + canvas_config.skew_padding
            ribbon_start_y =  record_offsets[i] + height_below_axis
            comparison_offsets.append(ribbon_start_y)
            next_record_id = records[i+1].id
            next_label_height = record_label_heights.get(next_record_id, 0)
            ribbon_end_y = record_offsets[i+1] - canvas_config.cds_padding
            height = ribbon_end_y - ribbon_start_y
            actual_comparison_heights.append(height)

        canvas = add_comparison_on_linear_canvas(
            canvas, comparisons, canvas_config, blast_config, config_dict, records, 
            comparison_offsets, actual_comparison_heights
        )

    
    for count, record in enumerate(records, start=1):
        offset_y = record_offsets[count-1]
        offset_x = (canvas_config.alignment_width *
                    ((canvas_config.longest_genome - len(record.seq)) / canvas_config.longest_genome) / 2) if canvas_config.align_center else 0
        
        labels_for_record = all_labels.get(record.id)
        add_record_group(canvas, record, offset_y, offset_x, canvas_config, feature_config, config_dict, precalculated_labels=labels_for_record)
        add_record_definition_group(canvas, record, offset_y, offset_x, canvas_config, config_dict)
        if canvas_config.show_gc:
            add_gc_content_group(canvas, record, offset_y, offset_x, canvas_config, gc_config, config_dict)
        if canvas_config.show_skew:
            add_gc_skew_group(canvas, record, offset_y, offset_x, canvas_config, skew_config, config_dict)

    save_figure(canvas, out_formats)
