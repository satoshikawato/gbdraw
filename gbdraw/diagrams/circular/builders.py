#!/usr/bin/env python
# coding: utf-8

"""Group-builder helpers for circular diagram assembly."""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable, Literal, ParamSpec, Sequence

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from svgwrite import Drawing
from svgwrite.container import Group

from ...auto_sizes import circular_track_width_px
from ...canvas import CircularCanvasConfigurator
from ...analysis.conservation import ConservationTrack
from ...features.factory import FeatureBuildResult
from ...labels.circular_types import CircularLabel
from ...layout.circular import CircularRecordRenderContext
from ...legend.table import LegendTable
from ...configurators import (
    FeatureDrawingConfigurator,
    DepthConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
    LegendMeasurement,
)
from ...render.groups.circular import (
    AxisGroup,
    CircularDefinitionStyle,
    ConservationGroup,
    DefinitionGroup,
    DepthGroup,
    GcContentGroup,
    GcSkewGroup,
    LabelsGroup,
    SeqRecordGroup,
    TickGroup,
)
from ...render.group_content import (
    CANVAS_SCOPE,
    BuildsGroup,
    GroupContent,
    GroupKey,
    GroupKind,
    build_group_content,
    place_content,
)
from ...render_stages import render_stage
from ...svg.ids import stable_svg_id
from ...render.groups.circular.legend import circular_legend_content
from .positioning import (
    center_group_on_canvas,
    place_definition_group_on_canvas,
)


from gbdraw.layout.record_coordinates import RecordDisplayTransform

if TYPE_CHECKING:
    from gbdraw.api.requests import RecordPresentation

P = ParamSpec("P")


def _namespace_slot_group(group: Group, slot_group_id: str) -> None:
    """Give a slot group its slot DOM id and namespace its renderer-owned child IDs."""

    subtree: list[Group] = []

    def collect(node: Group) -> None:
        subtree.append(node)
        for child in getattr(node, "elements", []) or []:
            collect(child)

    collect(group)
    id_map: dict[str, str] = {}
    for child_index, child in enumerate(subtree):
        child_attribs = getattr(child, "attribs", None)
        if not isinstance(child_attribs, dict):
            continue
        old_id = str(child_attribs.get("id", ""))
        if not old_id:
            continue
        new_id = (
            slot_group_id
            if child is group
            else stable_svg_id("track_slot_child", slot_group_id, old_id, child_index)
        )
        child_attribs["id"] = new_id
        id_map[old_id] = new_id

    for child in subtree:
        child_attribs = getattr(child, "attribs", None)
        if not isinstance(child_attribs, dict):
            continue
        for name, value in tuple(child_attribs.items()):
            if not isinstance(value, str):
                continue
            updated = value
            for old_id, new_id in id_map.items():
                updated = updated.replace(f"url(#{old_id})", f"url(#{new_id})")
                if updated == f"#{old_id}":
                    updated = f"#{new_id}"
            if updated != value:
                child_attribs[name] = updated


def slot_group_content(
    slot_group_id: str | None,
    group_class: Callable[P, BuildsGroup],
    /,
    *args: P.args,
    **kwargs: P.kwargs,
) -> GroupContent:
    """A render-stage function: build one track group and give it its slot IDs.

    With ``slot_group_id`` the group and its child IDs are namespaced inside
    the stage, before the content is kept, so a stored content is never
    changed afterwards.
    """

    group = group_class(*args, **kwargs).get_group()
    if slot_group_id is not None:
        _namespace_slot_group(group, slot_group_id)
    return GroupContent.of(group)


def _add_track_content(
    canvas: Drawing,
    content: GroupContent,
    canvas_config: CircularCanvasConfigurator,
    group_key: GroupKey,
    *,
    slot_id: str | None,
    slot_renderer: str | None,
) -> Drawing:
    group = center_group_on_canvas(place_content(canvas, content, group_key), canvas_config)
    if slot_id is not None and slot_renderer is not None:
        group.attribs["data-gbdraw-slot-id"] = str(slot_id)
        group.attribs["data-gbdraw-slot-renderer"] = str(slot_renderer)
    canvas.add(group)
    return canvas

def add_depth_group_on_canvas(
    canvas: Drawing,
    gb_record: SeqRecord,
    depth_df: DataFrame,
    canvas_config: CircularCanvasConfigurator,
    depth_config: DepthConfigurator,
    *,
    group_key: GroupKey,
    track_width_override: float | None = None,
    norm_factor_override: float,
    group_id: str | None = None,
    axis_group_id: str | None = None,
    slot_group_id: str | None = None,
    slot_id: str | None = None,
    slot_renderer: str | None = None,
    record_transform: RecordDisplayTransform | None = None,
) -> Drawing:
    """Adds the depth coverage group to the canvas."""
    depth_track_width = (
        float(track_width_override)
        if track_width_override is not None
        else circular_track_width_px(
            "depth",
            radius=canvas_config.radius,
            track_ratio=canvas_config.track_ratio,
            factors=canvas_config.track_ratio_factors,
        )
    )
    content = render_stage(
        "depth_content",
        slot_group_content,
        slot_group_id,
        DepthGroup,
        gb_record,
        depth_df,
        canvas_config.radius,
        depth_track_width,
        depth_config,
        norm_factor_override=norm_factor_override,
        group_id=group_id,
        axis_group_id=axis_group_id,
        record_transform=record_transform,
    )
    return _add_track_content(
        canvas, content, canvas_config, group_key, slot_id=slot_id, slot_renderer=slot_renderer
    )


def add_conservation_group_on_canvas(
    canvas: Drawing,
    gb_record: SeqRecord,
    conservation_track: ConservationTrack,
    canvas_config: CircularCanvasConfigurator,
    *,
    group_key: GroupKey,
    inner_radius_px: float,
    outer_radius_px: float,
    min_identity: float,
    group_id: str | None = None,
    slot_id: str | None = None,
    record_transform: RecordDisplayTransform | None = None,
    record_index: int = 0,
) -> Drawing:
    """Adds one circular conservation ring group to the canvas."""
    content = render_stage(
        "conservation_content",
        build_group_content,
        ConservationGroup,
        track=conservation_track,
        total_length=len(gb_record.seq),
        inner_radius_px=inner_radius_px,
        outer_radius_px=outer_radius_px,
        min_identity=min_identity,
        cfg=canvas_config.profile.config,
        group_id=group_id,
        slot_id=slot_id,
        record_transform=record_transform,
        record_index=record_index,
    )
    canvas.add(center_group_on_canvas(place_content(canvas, content, group_key), canvas_config))
    return canvas


def add_gc_skew_group_on_canvas(
    canvas: Drawing,
    gb_record: SeqRecord,
    gc_df: DataFrame,
    canvas_config: CircularCanvasConfigurator,
    skew_config: GcSkewConfigurator,
    *,
    group_key: GroupKey,
    track_width_override: float | None = None,
    norm_factor_override: float,
    group_id: str | None = None,
    slot_group_id: str | None = None,
    slot_id: str | None = None,
    slot_renderer: str | None = None,
    record_transform: RecordDisplayTransform | None = None,
) -> Drawing:
    """Adds the GC skew group to the canvas.

    With ``slot_group_id`` the group takes its slot DOM ids; ``slot_id`` and
    ``slot_renderer`` tag the placed group.
    """
    skew_track_width = (
        float(track_width_override)
        if track_width_override is not None
        else circular_track_width_px(
            "dinucleotide_skew",
            radius=canvas_config.radius,
            track_ratio=canvas_config.track_ratio,
            factors=canvas_config.track_ratio_factors,
        )
    )
    content = render_stage(
        "gc_skew_content",
        slot_group_content,
        slot_group_id,
        GcSkewGroup,
        gb_record,
        gc_df,
        canvas_config.radius,
        skew_track_width,
        skew_config,
        norm_factor_override=norm_factor_override,
        group_id=group_id,
        record_transform=record_transform,
    )
    return _add_track_content(
        canvas, content, canvas_config, group_key, slot_id=slot_id, slot_renderer=slot_renderer
    )


def add_gc_content_group_on_canvas(
    canvas: Drawing,
    gb_record: SeqRecord,
    gc_df: DataFrame,
    canvas_config: CircularCanvasConfigurator,
    gc_config: GcContentConfigurator,
    *,
    group_key: GroupKey,
    track_width_override: float | None = None,
    norm_factor_override: float,
    group_id: str | None = None,
    slot_group_id: str | None = None,
    slot_id: str | None = None,
    slot_renderer: str | None = None,
    record_transform: RecordDisplayTransform | None = None,
) -> Drawing:
    """Adds the GC content group to the canvas.

    With ``slot_group_id`` the group takes its slot DOM ids; ``slot_id`` and
    ``slot_renderer`` tag the placed group.
    """
    gc_content_track_width = (
        float(track_width_override)
        if track_width_override is not None
        else circular_track_width_px(
            "dinucleotide_content",
            radius=canvas_config.radius,
            track_ratio=canvas_config.track_ratio,
            factors=canvas_config.track_ratio_factors,
        )
    )
    content = render_stage(
        "gc_content_content",
        slot_group_content,
        slot_group_id,
        GcContentGroup,
        gb_record,
        gc_df,
        canvas_config.radius,
        gc_content_track_width,
        gc_config,
        norm_factor_override=norm_factor_override,
        group_id=group_id,
        record_transform=record_transform,
    )
    return _add_track_content(
        canvas, content, canvas_config, group_key, slot_id=slot_id, slot_renderer=slot_renderer
    )


def add_record_definition_group_on_canvas(
    canvas: Drawing,
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    species: str | None,
    strain: str | None,
    *,
    group_key: GroupKey,
    plot_title: str | None = None,
    definition_profile: str = "full",
    definition_position: str = "center",
    definition_group_id: str | None = None,
    record_index: int = 0,
    record_count: int = 1,
    record_transform: RecordDisplayTransform | None = None,
    species_line_count: int = 1,
    record_presentation: RecordPresentation | None = None,
) -> Drawing:
    """
    Adds the record definition group to the canvas.

    Parameters:
    canvas (Drawing): The SVG drawing canvas.
    gb_record (SeqRecord): The GenBank record.
    canvas_config (CircularCanvasConfigurator): Configuration for the circular canvas.
    species (str): Species name.
    strain (str): Strain name.
    Returns:
    Drawing: The updated SVG drawing with the record definition group added.
    """
    content = render_stage(
        "circular_definition_content",
        build_group_content,
        DefinitionGroup,
        gb_record,
        canvas_config.geometry(),
        style=CircularDefinitionStyle.of(canvas_config.profile.config, definition_profile),
        species=species,
        strain=strain,
        plot_title=plot_title,
        definition_profile=definition_profile,
        definition_group_id=definition_group_id,
        record_index=record_index,
        record_count=record_count,
        record_transform=record_transform,
        species_line_count=species_line_count,
        record_presentation=record_presentation,
    )
    definition_group = place_content(canvas, content, group_key)
    setattr(definition_group, "_gbdraw_local_bounds", content.local_bounds)
    definition_group = place_definition_group_on_canvas(
        definition_group,
        canvas_config,
        position=definition_position,
    )
    canvas.add(definition_group)
    return canvas


def add_record_group_on_canvas(
    canvas: Drawing,
    record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    *,
    group_key: GroupKey,
    feature_layers: FeatureBuildResult,
    render_context: CircularRecordRenderContext,
    precalculated_labels: tuple[CircularLabel, ...] = (),
    feature_track_ratio_factor_override: float | None = None,
    feature_anchor_radius_px: float | None = None,
    record_index: int = 0,
    group_id: str | None = None,
    slot_id: str | None = None,
    feature_dom_namespace: str | None = None,
) -> Drawing:
    """Adds the record group (features and embedded labels) to the canvas."""
    content = render_stage(
        "circular_record_content",
        build_group_content,
        SeqRecordGroup,
        gb_record=record,
        geometry=canvas_config.geometry(),
        feature_config=feature_config,
        feature_layers=feature_layers,
        render_context=render_context,
        precalculated_labels=precalculated_labels,
        feature_track_ratio_factor_override=feature_track_ratio_factor_override,
        feature_anchor_radius_px=feature_anchor_radius_px,
        record_index=record_index,
        group_id=group_id,
        slot_id=slot_id,
        feature_dom_namespace=feature_dom_namespace,
    )
    return _add_track_content(canvas, content, canvas_config, group_key, slot_id=None, slot_renderer=None)


def add_axis_group_on_canvas(
    canvas: Drawing,
    canvas_config: CircularCanvasConfigurator,
    *,
    group_key: GroupKey,
) -> Drawing:
    """Adds the axis circle at the canvas radius."""
    content = render_stage(
        "axis_content",
        build_group_content,
        AxisGroup,
        canvas_config.radius,
        canvas_config.geometry(),
        cfg=canvas_config.profile.config,
    )
    return _add_track_content(canvas, content, canvas_config, group_key, slot_id=None, slot_renderer=None)


def add_tick_group_on_canvas(
    canvas: Drawing,
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    *,
    group_key: GroupKey,
    radius_override: float | None = None,
    label_side: str = "legacy",
    tick_side: str = "legacy",
    tick_length_px: float | None = None,
    track_preset: str | None = None,
    group_id: str | None = None,
    slot_id: str | None = None,
    record_transform: RecordDisplayTransform | None = None,
) -> Drawing:
    """Adds the tick group to the canvas."""
    content = render_stage(
        "ticks_content",
        build_group_content,
        TickGroup,
        gb_record,
        canvas_config.geometry(),
        profile=canvas_config.profile,
        radius=radius_override,
        label_side=label_side,
        tick_side=tick_side,
        tick_length_px=tick_length_px,
        track_preset=track_preset,
        group_id=group_id,
        slot_id=slot_id,
        record_transform=record_transform,
    )
    return _add_track_content(canvas, content, canvas_config, group_key, slot_id=None, slot_renderer=None)


def add_labels_group_on_canvas(
    canvas: Drawing,
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    *,
    group_key: GroupKey,
    render_context: CircularRecordRenderContext,
    precalculated_labels: Sequence[CircularLabel],
    phase: Literal["leaders", "text"],
) -> Drawing:
    """Add the leader lines or the text of the placed external labels.

    Two stages, so a leader stroke change keeps the text and a font change
    keeps the leaders. The labels are this request's copies after the legend
    shift, keyed by value; the record assembler pre-keys them once for both
    phases, and a tuple is passed on as it is so that key is kept.
    """
    content = render_stage(
        "label_leaders_content" if phase == "leaders" else "label_text_content",
        build_group_content,
        LabelsGroup,
        len(gb_record.seq),
        canvas_config.geometry(),
        render_context,
        precalculated_labels=(
            precalculated_labels if isinstance(precalculated_labels, tuple) else tuple(precalculated_labels)
        ),
        phase=phase,
    )
    canvas.add(center_group_on_canvas(place_content(canvas, content, group_key), canvas_config))
    return canvas


def add_legend_group_on_canvas(
    canvas: Drawing,
    legend_measurement: LegendMeasurement,
    legend_table: LegendTable,
) -> Drawing:
    """
    Adds the legend group to the canvas.

    Parameters:
    canvas (Drawing): The SVG drawing canvas.
    legend_measurement: Immutable measurement for the legend.
    legend_table: The legend table data.

    Returns:
    Drawing: The updated SVG drawing with the legend group added.
    """
    content = render_stage("legend_content", circular_legend_content, legend_measurement, legend_table)
    canvas.add(place_content(canvas, content, GroupKey(CANVAS_SCOPE, GroupKind.LEGEND)))
    return canvas


__all__ = [
    "add_depth_group_on_canvas",
    "add_gc_skew_group_on_canvas",
    "add_gc_content_group_on_canvas",
    "add_record_definition_group_on_canvas",
    "add_record_group_on_canvas",
    "add_axis_group_on_canvas",
    "add_tick_group_on_canvas",
    "add_labels_group_on_canvas",
    "add_legend_group_on_canvas",
]
