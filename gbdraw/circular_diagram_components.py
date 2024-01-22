#!/usr/bin/env python
# coding: utf-8

import logging
import sys
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from .canvas_generator import CircularCanvasConfigurator
from .circular_object_groups import GcContentGroup, GcSkewGroup, SeqRecordGroup, DefinitionGroup, TickGroup, AxisGroup
from .object_configurators import GcContentConfigurator, FeatureDrawingConfigurator


from .file_processing import save_figure
from svgwrite import Drawing
from svgwrite.container import Group
# Logging setup
logger = logging.getLogger()
handler = logging.StreamHandler(sys.stdout)


def center_group_on_canvas(group: Group, canvas_config: CircularCanvasConfigurator) -> Group:
    """
    Centers a given SVG group on the canvas based on the canvas configuration.

    Parameters:
    group (Group): The SVG group to be centered.
    canvas_config (CircularCanvasConfigurator): The configuration of the circular canvas.

    Returns:
    Group: The centered SVG group.
    """
    group.translate(canvas_config.offset_x, canvas_config.offset_y)
    return group


def add_gc_skew_group_on_canvas(canvas: Drawing, gb_record: SeqRecord, gc_df: DataFrame, canvas_config: CircularCanvasConfigurator, config_dict: dict) -> Drawing:
    """
    Adds the GC skew group to the canvas.

    Parameters:
    canvas (Drawing): The SVG drawing canvas.
    gb_record (SeqRecord): The GenBank record.
    gc_df (DataFrame): DataFrame containing GC content and skew information.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
    Drawing: The updated SVG drawing with the GC skew group added.
    """
    gc_skew_group: Group = GcSkewGroup(gb_record, gc_df, canvas_config.radius, canvas_config.track_width,
                                       config_dict, canvas_config.track_ids['skew_track']).get_group()
    gc_skew_group = center_group_on_canvas(gc_skew_group, canvas_config)
    canvas.add(gc_skew_group)
    return canvas


def add_gc_content_group_on_canvas(canvas: Drawing, gb_record: SeqRecord, gc_df: DataFrame, canvas_config: CircularCanvasConfigurator, gc_config: GcContentConfigurator, config_dict: dict) -> Drawing:
    """
    Adds the GC content group to the canvas.

    Parameters:
    canvas (Drawing): The SVG drawing canvas.
    gb_record (SeqRecord): The GenBank record.
    gc_df (DataFrame): DataFrame containing GC content and skew information.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    gc_config (GcContentConfigurator): Configuration for the GC content representation.
    config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
    Drawing: The updated SVG drawing with the GC content group added.
    """
    gc_content_group: Group = GcContentGroup(gb_record, gc_df, canvas_config.radius, canvas_config.track_width,
                                             gc_config, config_dict, canvas_config.track_ids['gc_track']).get_group()
    gc_content_group = center_group_on_canvas(gc_content_group, canvas_config)
    canvas.add(gc_content_group)
    return canvas


def add_record_definition_group_on_canvas(canvas: Drawing, gb_record: SeqRecord, canvas_config: CircularCanvasConfigurator, species: str, strain: str, config_dict: dict) -> Drawing:
    """
    Adds the record definition group to the canvas.

    Parameters:
    canvas (Drawing): The SVG drawing canvas.
    gb_record (SeqRecord): The GenBank record.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    species (str): Species name.
    strain (str): Strain name.
    config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
    Drawing: The updated SVG drawing with the record definition group added.
    """
    definition_group: Group = DefinitionGroup(
        gb_record, canvas_config, species=species, strain=strain, config_dict=config_dict).get_group()
    definition_group = center_group_on_canvas(definition_group, canvas_config)
    canvas.add(definition_group)
    return canvas


def add_record_group_on_canvas(canvas: Drawing, record: SeqRecord, canvas_config: CircularCanvasConfigurator, feature_config: FeatureDrawingConfigurator, config_dict: dict) -> Drawing:
    """
    Adds the record group to the canvas.

    Parameters:
    canvas (Drawing): The SVG drawing canvas.
    record (SeqRecord): The GenBank record.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    feature_config (FeatureDrawingConfigurator): Configuration for feature drawing.
    config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
    Drawing: The updated SVG drawing with the record group added.
    """
    # Unpack parameters
    record_group: Group = SeqRecordGroup(
        gb_record=record, canvas_config=canvas_config, feature_config=feature_config, config_dict=config_dict).get_group()
    record_group = center_group_on_canvas(record_group, canvas_config)
    canvas.add(record_group)
    return canvas


def add_axis_group_on_canvas(canvas: Drawing, canvas_config: CircularCanvasConfigurator, config_dict: dict) -> Drawing:
    """
    Adds the axis group to the canvas.

    Parameters:
    canvas (Drawing): The SVG drawing canvas.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
    Drawing: The updated SVG drawing with the axis group added.
    """
    axis_group: Group = AxisGroup(
        canvas_config.radius, config_dict).get_group()
    axis_group = center_group_on_canvas(axis_group, canvas_config)
    canvas.add(axis_group)
    return canvas


def add_tick_group_on_canvas(canvas: Drawing, gb_record: SeqRecord, canvas_config: CircularCanvasConfigurator, config_dict: dict) -> Drawing:
    """
    Adds the tick group to the canvas.

    Parameters:
    canvas (Drawing): The SVG drawing canvas.
    gb_record (SeqRecord): The GenBank record.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
    Drawing: The updated SVG drawing with the tick group added.
    """
    tick_group: Group = TickGroup(
        gb_record, canvas_config, config_dict).get_group()
    tick_group = center_group_on_canvas(tick_group, canvas_config)
    canvas.add(tick_group)
    return canvas


def add_record_on_circular_canvas(canvas: Drawing, gb_record: SeqRecord, canvas_config: CircularCanvasConfigurator, feature_config: FeatureDrawingConfigurator, gc_config: GcContentConfigurator, gc_df: DataFrame, species: str, strain: str, config_dict: dict) -> Drawing:
    """
    Adds various record-related groups to a circular canvas.

    Parameters:
    canvas (Drawing): The SVG drawing canvas.
    gb_record (SeqRecord): The GenBank record.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    feature_config (FeatureDrawingConfigurator): Configuration for feature drawing.
    gc_config (GcContentConfigurator): Configuration for GC content representation.
    gc_df (DataFrame): DataFrame containing GC content and skew information.
    species (str): Species name.
    strain (str): Strain name.
    config_dict (dict): Configuration dictionary for drawing parameters.

    Returns:
    Drawing: The updated SVG drawing with all record-related groups added.
    """
    canvas = add_axis_group_on_canvas(canvas, canvas_config, config_dict)
    # Add record group
    canvas = add_record_group_on_canvas(
        canvas, gb_record, canvas_config, feature_config, config_dict)
    # Add record definition group
    canvas = add_record_definition_group_on_canvas(
        canvas, gb_record, canvas_config, species, strain, config_dict)
    # Add record tick group
    canvas = add_tick_group_on_canvas(
        canvas, gb_record, canvas_config, config_dict)

    # Add GC content group if configured to show
    if canvas_config.show_gc:
        canvas = add_gc_content_group_on_canvas(
            canvas, gb_record, gc_df, canvas_config, gc_config, config_dict)
    if canvas_config.show_skew:
        canvas = add_gc_skew_group_on_canvas(
            canvas, gb_record, gc_df, canvas_config, config_dict)
    return canvas


def plot_circular_diagram(gb_record: SeqRecord, canvas_config: CircularCanvasConfigurator, gc_df: DataFrame, gc_config: GcContentConfigurator, feature_config: FeatureDrawingConfigurator, species: str, strain: str, config_dict: dict, out_formats: list) -> None:
    """
    Plots a circular diagram for a GenBank record.

    Parameters:
    gb_record (SeqRecord): The GenBank record to plot.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    gc_df (DataFrame): DataFrame containing GC content and skew information.
    gc_config (GcContentConfigurator): Configuration for GC content representation.
    feature_config (FeatureDrawingConfigurator): Configuration for feature drawing.
    species (str): Species name.
    strain (str): Strain name.
    config_dict (dict): Configuration dictionary for drawing parameters.
    out_formats (list): List of formats to save the output (e.g., ['png', 'svg']).

    Returns:
    None: The function saves the plotted diagram to specified output formats.
    """
    # Configure and create canvas
    canvas: Drawing = canvas_config.create_svg_canvas()

    canvas = add_record_on_circular_canvas(
        canvas, gb_record, canvas_config, feature_config, gc_config, gc_df, species, strain, config_dict)
    save_figure(canvas, out_formats)
