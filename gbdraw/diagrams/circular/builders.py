#!/usr/bin/env python
# coding: utf-8

"""Group-builder helpers for circular diagram assembly."""

from __future__ import annotations

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ...canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...configurators import (  # type: ignore[reportMissingImports]
    FeatureDrawingConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
)
from ...render.groups.circular import (  # type: ignore[reportMissingImports]
    AxisGroup,
    DefinitionGroup,
    GcContentGroup,
    GcSkewGroup,
    LabelsGroup,
    LegendGroup,
    SeqRecordGroup,
    TickGroup,
)
from .positioning import center_group_on_canvas, place_legend_on_canvas


def add_gc_skew_group_on_canvas(
    canvas: Drawing,
    gb_record: SeqRecord,
    gc_df: DataFrame,
    canvas_config: CircularCanvasConfigurator,
    skew_config: GcSkewConfigurator,
    config_dict: dict,
    *,
    track_width_override: float | None = None,
    norm_factor_override: float | None = None,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
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
    cfg = cfg or canvas_config._cfg
    skew_track_width = (
        float(track_width_override)
        if track_width_override is not None
        else canvas_config.radius * canvas_config.track_ratio * canvas_config.track_ratio_factors[2]
    )
    gc_skew_group: Group = GcSkewGroup(
        gb_record,
        gc_df,
        canvas_config.radius,
        skew_track_width,
        skew_config,
        config_dict,
        canvas_config.track_ids["skew_track"],
        norm_factor_override=norm_factor_override,
        cfg=cfg,
    ).get_group()
    gc_skew_group = center_group_on_canvas(gc_skew_group, canvas_config)
    canvas.add(gc_skew_group)
    return canvas


def add_gc_content_group_on_canvas(
    canvas: Drawing,
    gb_record: SeqRecord,
    gc_df: DataFrame,
    canvas_config: CircularCanvasConfigurator,
    gc_config: GcContentConfigurator,
    config_dict: dict,
    *,
    track_width_override: float | None = None,
    norm_factor_override: float | None = None,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
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
    cfg = cfg or canvas_config._cfg
    gc_content_track_width = (
        float(track_width_override)
        if track_width_override is not None
        else canvas_config.radius * canvas_config.track_ratio * canvas_config.track_ratio_factors[1]
    )
    gc_content_group: Group = GcContentGroup(
        gb_record,
        gc_df,
        canvas_config.radius,
        gc_content_track_width,
        gc_config,
        config_dict,
        canvas_config.track_ids["gc_track"],
        norm_factor_override=norm_factor_override,
        cfg=cfg,
    ).get_group()
    gc_content_group = center_group_on_canvas(gc_content_group, canvas_config)
    canvas.add(gc_content_group)
    return canvas


def add_record_definition_group_on_canvas(
    canvas: Drawing,
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    species: str,
    strain: str,
    config_dict: dict,
    *,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
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
        gb_record, canvas_config, species=species, strain=strain, config_dict=config_dict, cfg=cfg or canvas_config._cfg).get_group()
    definition_group = center_group_on_canvas(definition_group, canvas_config)
    canvas.add(definition_group)
    return canvas


def add_record_group_on_canvas(
    canvas: Drawing,
    record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    config_dict: dict,
    *,
    cfg: GbdrawConfig | None = None,
    precomputed_feature_dict: dict | None = None,
    precalculated_labels: list[dict] | None = None,
) -> Drawing:
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
        gb_record=record,
        canvas_config=canvas_config,
        feature_config=feature_config,
        config_dict=config_dict,
        cfg=cfg or canvas_config._cfg,
        precomputed_feature_dict=precomputed_feature_dict,
        precalculated_labels=precalculated_labels,
    ).get_group()
    # Calculate start and end points for the 60-degree arc

    record_group = center_group_on_canvas(record_group, canvas_config)
    canvas.add(record_group)

    return canvas


def add_axis_group_on_canvas(
    canvas: Drawing,
    canvas_config: CircularCanvasConfigurator,
    config_dict: dict,
    *,
    radius_override: float | None = None,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
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
        float(radius_override) if radius_override is not None else canvas_config.radius,
        config_dict,
        canvas_config,
        cfg=cfg or canvas_config._cfg,
    ).get_group()
    axis_group = center_group_on_canvas(axis_group, canvas_config)
    canvas.add(axis_group)
    return canvas


def add_tick_group_on_canvas(
    canvas: Drawing,
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    config_dict: dict,
    *,
    radius_override: float | None = None,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
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
        gb_record,
        canvas_config,
        config_dict,
        radius=radius_override,
        cfg=cfg or canvas_config._cfg,
    ).get_group()
    tick_group = center_group_on_canvas(tick_group, canvas_config)
    canvas.add(tick_group)
    return canvas


def add_labels_group_on_canvas(
    canvas: Drawing,
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    config_dict: dict,
    *,
    outer_arena: tuple[float, float] | None = None,
    cfg: GbdrawConfig | None = None,
    precomputed_feature_dict: dict | None = None,
    precalculated_labels: list[dict] | None = None,
) -> Drawing:
    """
    Adds the labels group to the canvas.

    Parameters:
    canvas (Drawing): The SVG drawing canvas.
    gb_record (SeqRecord): The GenBank record.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    feature_config (FeatureDrawingConfigurator): Configuration for feature drawing.
    config_dict (dict): Configuration dictionary for drawing parameters.
    outer_arena (tuple[float, float] | None): Optional outer arena bounds.

    Returns:
    Drawing: The updated SVG drawing with the labels group added.
    """
    labels_group: Group = LabelsGroup(
        gb_record=gb_record,
        canvas_config=canvas_config,
        feature_config=feature_config,
        config_dict=config_dict,
        outer_arena=outer_arena,
        cfg=cfg or canvas_config._cfg,
        precomputed_feature_dict=precomputed_feature_dict,
        precalculated_labels=precalculated_labels,
    ).get_group()
    labels_group = center_group_on_canvas(labels_group, canvas_config)
    canvas.add(labels_group)
    return canvas


def add_legend_group_on_canvas(canvas: Drawing, canvas_config: CircularCanvasConfigurator, legend_config, legend_table) -> Drawing:
    """
    Adds the legend group to the canvas.

    Parameters:
    canvas (Drawing): The SVG drawing canvas.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    legend_config: Configuration for the legend.
    legend_table: The legend table data.

    Returns:
    Drawing: The updated SVG drawing with the legend group added.
    """
    legend_group = LegendGroup(canvas_config, legend_config, legend_table).get_group()
    legend_group = place_legend_on_canvas(legend_group, canvas_config)
    canvas.add(legend_group)

    return canvas


__all__ = [
    "add_gc_skew_group_on_canvas",
    "add_gc_content_group_on_canvas",
    "add_record_definition_group_on_canvas",
    "add_record_group_on_canvas",
    "add_axis_group_on_canvas",
    "add_tick_group_on_canvas",
    "add_labels_group_on_canvas",
    "add_legend_group_on_canvas",
]
