#!/usr/bin/env python
# coding: utf-8

"""Circular diagram assembly (implementation).

This module was extracted from `gbdraw.circular_diagram_components` to improve cohesion.
"""

from __future__ import annotations

from gbdraw.features.placement import FeaturePlacementSlot, ResolvedRecordFeatureInputs

import logging
import math
import copy
from collections.abc import Hashable, Iterator
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field, replace
from functools import partial
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Callable, Collection, Generic, Literal, Mapping, NamedTuple, Sequence, TypeAlias, TypedDict, TypeVar, cast

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from svgwrite import Drawing
from svgwrite.container import Group

from ...auto_sizes import circular_track_width_px
from ...canvas import CircularCanvasConfigurator, CircularCanvasGeometry
from ...analysis.conservation import (
    ConservationLegendSource,
    ConservationTrack,
    conservation_legend_sources,
    named_conservation_tracks,
    conservation_slot_track,
    conservation_track_gradient_colors,
)
from ...analysis.feature_density import FeatureDensityLayer
from ..feature_density import (
    feature_density_color,
    feature_density_histogram_source,
    feature_density_legend,
    slot_feature_density_layer,
)
from ...analysis.depth_tracks import (
    DepthTrackData,
    depth_track_data_count,
    depth_track_legend_name,
    index_depth_track_row,
)
from ...config.models import (
    CircularRenderProfile,
    GbdrawConfig,
)
from ...configurators import (
    FeatureDrawingConfigurator,
    DepthConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
    LegendMeasurement,
)
from ...configurators.bundle import DiagramConfigurators
from ...configurators.gc import _slot_skew_config
from ...configurators.legend import measure_legend
from ...render_stages import _pre_keyed, render_stage
from ...config.modify import with_config_overrides, with_fit_row_overrides
from ...exceptions import ValidationError
from ...features.factory import FeatureBuildResult, create_feature_layers
from ...features.objects import FeatureObject
from ...labels.label_ink import circular_label_ink
from ...labels.label_style import FeatureLabelStyle, label_leader_width_px
from ...render.drawers.label_backing import feature_label_drawn_parts
from ...labels.circular import (
    CircularLabelProbeStopped,
    CircularLabelRadii,
    assign_leader_start_points,
    circular_label_fits,
    circular_label_font_size,
    circular_label_shape,
    minimum_bbox_gap_px,
    circular_label_radii,
    prepare_label_list,
    x_overlap,
    y_overlap,
)
from ...labels.circular_repair import make_way_for_region_labels
from ...labels.circular_capacity import CircularLabelCapacity, CircularLabelCapacityExceeded
from ...labels.circular_fit import CircularDefinitionFitProbe, CircularLabelFitProbe, CircularLabelLayout
from ...labels.circular_radial import RADIAL_GROWTH_TOLERANCE_PX
from ...labels.circular_types import CircularLabel, CircularLabelCandidate, CircularLabelCandidateCache
from ...labels.fit import (
    DEFINITION_FONT_PATH,
    FITS,
    LabelFit,
    LabelFitValue,
    LabelKind,
    LabelShape,
    Polygon,
    RegionLabelWarnings,
    region_rooms_left_out,
)
from ...labels.filtering import label_filter_maps
from ...layout.record_coordinates import RecordDisplayTransform
from ...layout.record_transform import bind_record_transforms, circular_display_length, unset_record_transform
from ...render.track_slot_metadata import (
    RadialExtent,
    CircularTrackSlotEntry,
    CircularTrackSlotRecord,
    RecordAxisInput,
    attach_feature_placement_targets,
    attach_track_slot_geometry,
    circular_axis_geometry,
)
from ...layout.circular import (
    CircularAxisLayout,
    CircularFeatureLaneDirection,
    CircularFeatureLayout,
    CircularRadialLayout,
    CircularRecordRenderContext,
    CircularResolvedSlot,
    CircularTickLayout,
    RadialBand,
    feature_radius_intervals,
    slots_feature_layout,
)
from ...layout.circular_depth_axis import CircularScalarAxes, DepthSlotWidth
from ...layout.circular_ink import CircularInk, CircularInkRing, box_meets_ink, stroke_pad_px
from ...layout.composition import (
    CompositionInputs,
    CompositionItem,
    CompositionPlan,
    CompositionRequest,
    DEFAULT_COMPOSITION_SPACING,
    LegendPlacement,
    plan_composition,
)
from ...layout.spatial import Aabb, BoxedLetters, bounding_aabb, union_aabbs
from ...comparisons.identity_scale import IdentityColorDomain, identity_color_domain
from ...legend.row_facts import LegendRowFacts, attach_legend_row_facts
from ...legend.sources import (
    DepthLegendSource,
    LegendSlotSource,
    LegendSources,
    NumericLegendSource,
    depth_legend_sources,
    feature_legend_presence,
    legend_slot_sources,
)
from ...legend.table import (
    LegendRowEdits,
    LegendTable,
    apply_legend_row_edits,
    prepare_legend_table,
    sync_feature_density_legend_entries,
)
from ...tracks import (
    CircularTrackSlot,
    ScalarSpec,
    default_circular_track_slots,
    normalize_circular_track_slots,
    normalize_circular_track_slots_with_axis,
)
from ...annotations import (
    AnnotationTrackParams,
    LabelSizing,
    LaneRoom,
    ResolvedAnnotationBundle,
    ResolvedAnnotationTrack,
    ResolvedRegionAnnotation,
    RoomBounds,
    annotation_track_params_from_mapping,
    feature_underlay_anchor_slot_id,
    feature_underlay_slot_id,
    layout_annotation_track,
    merge_feature_underlays,
    sync_annotation_legend_entries,
)
from ...annotations.planning import prepare_annotation_track_slots
from ...annotations.models import slots_in_paint_order
from ...render.drawers.circular.depth import DepthDrawer
from ...render.drawers.circular.feature_density import (
    FeatureDensityHeatmapGroup,
    feature_density_heatmap_ink,
    feature_density_paints_band_edges,
)
from ...render.drawers.circular.record_band import (
    RECORD_BAND_SLOT_ID,
    RECORD_BAND_WIDTH_EM,
    RecordBandArc,
    record_band_group,
)
from ...render.drawers.circular.gc_content import GcContentDrawer
from ...render.groups.circular.conservation import conservation_ring_ink, conservation_ring_paints_band_edges
from ...render.groups.circular.ticks import TickGroup
from ...render.drawers.circular.gc_skew import SkewDrawer
from ...render.drawers.circular.annotations import (
    AnnotationRing,
    CircularAnnotationLabel,
    CircularAnnotationLabels,
    CircularLabelObstacles,
    CircularLabelRoom,
    CircularLaneLabels,
    CircularObstacleLabel,
    LeaderSegment,
    RegionLabelPlace,
    circular_annotation_label_extent,
    circular_annotation_label_room,
    circular_annotation_lane_labels,
    circular_annotation_lanes_meeting,
    circular_annotation_mark_past_lane,
    circular_annotation_labels,
    circular_annotation_labels_across_centre,
    circular_annotation_track_content,
    circular_annotation_track_ink,
    circular_region_label_front_content,
)
from ...render.group_content import GroupContent, GroupKey, GroupKind, place_content
from ...layout.linear_multi_record import record_key_at
from ...render.composition import (
    CompositionPrimary,
    PlacedLegend,
    compose_drawing,
    place_legend,
    place_title,
)
from ...render.groups.circular.legend import circular_legend_content
from ...svg.ids import track_slot_svg_id
from .positioning import center_group_on_canvas
from ...tracks.circular import tick_sides_for_tick_label_layout
from ...tracks.parsing import depth_slot_track_index, slot_dinucleotide, slot_legend_label, slots_draw_features

from .builders import (
    add_axis_group_on_canvas,
    add_conservation_group_on_canvas,
    add_depth_group_on_canvas,
    add_feature_density_heatmap_on_canvas,
    add_gc_content_group_on_canvas,
    add_gc_skew_group_on_canvas,
    add_labels_group_on_canvas,
    add_record_definition_group_on_canvas,
    add_record_group_on_canvas,
    add_tick_group_on_canvas,
    tick_group_content,
)
from ...render.groups.circular.definition import (
    CircularDefinitionProfile,
    CircularDefinitionStyle,
    DefinitionGroup,
    DefinitionNameLines,
    DefinitionRecordFacts,
    definition_record_facts,
)
from ...svg.circular_tracks import circle_unit_arc
from .one_circle import one_circle_definition_facts, record_band_names
from .parts import (
    CircularOneCircle,
    CircularPlannedTracks,
    CircularPresentationInputs,
    CircularRecordData,
    CircularTrackInputs,
)
from .radial_layout import (
    AnnotationLabelRoom,
    cannot_fit_diagnostic,
    resolve_circular_radial_layout,
    track_layout_reason,
)
from .presets import (
    CircularPresetContext,
    circular_feature_lane_direction_for_preset,
    circular_radial_plan_for_preset,
    circular_track_slots_from_preset_order,
    normalize_circular_track_preset,
)

if TYPE_CHECKING:
    from gbdraw.api.requests import RecordPresentation


LEGEND_LABEL_MARGIN_PX = 2.0 * DEFAULT_COMPOSITION_SPACING.overlay_clearance_px
LABEL_NUDGE_STEP_PX = 6.0
MAX_LABEL_NUDGE_PX = 180.0
MIN_LABEL_ORDER_GAP_RAD = 1e-4
FEATURE_BAND_EPSILON = 1e-6
logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class _CircularPlotAssembly:
    """Uncomposed Circular plot geometry and its existing SVG targets."""

    drawing: Drawing
    canvas_config: CircularCanvasConfigurator
    legend_measurement: LegendMeasurement
    legend_table: LegendTable
    source_content_bounds: Aabb
    source_overlay_obstacles: tuple[Aabb, ...]
    primary_targets: tuple[Group, ...]
    legend_placement: LegendPlacement
    legend: PlacedLegend | None
    track_slot_geometry: CircularTrackSlotGeometry
    label_layout: CircularLabelLayout | None


@dataclass(frozen=True, slots=True)
class CircularAssemblyResult:
    """Final single-record Circular assembly with authoritative geometry."""

    drawing: Drawing
    legend_measurement: LegendMeasurement
    legend_table: LegendTable
    source_content_bounds: Aabb
    content_bounds: Aabb
    overlay_obstacles: tuple[Aabb, ...]
    primary_targets: tuple[Group, ...]
    legend: PlacedLegend | None
    title_target: Group | None
    track_slot_geometry: CircularTrackSlotGeometry
    composition_plan: CompositionPlan
    # The drawn external labels: their fit by kind and a probe for the fit search (None without them).
    label_layout: CircularLabelLayout | None = None


def plan_circular_annotation_slots(
    tracks: CircularTrackInputs,
    records: Sequence[SeqRecord],
    *,
    show_ticks: bool,
    show_depth: bool,
    show_gc: bool,
    show_skew: bool,
    depth_track_count: int,
    record_transforms: Sequence[RecordDisplayTransform] | None = None,
    record_indices: Sequence[int] | None = None,
) -> tuple[list[CircularTrackSlot] | None, ResolvedAnnotationBundle]:
    """Bind the track inputs' annotation sets to the Circular slot list that the records use. A list the
    engine built (``implicit_track_stack``: rings, depth tracks) is the default list the sets' slots join,
    not a reader's stack that leaves them out (OV-2820)."""

    slots = tracks.circular_track_slots
    engine_slots = slots if tracks.implicit_track_stack else None
    planned, bundle, _auto_slot_ids = prepare_annotation_track_slots(
        tracks.annotations,
        records,
        None if tracks.implicit_track_stack else slots,
        mode="circular",
        record_transforms=record_transforms,
        record_indices=record_indices,
        default_slots=lambda: list(engine_slots) if engine_slots is not None else default_circular_track_slots(
            show_features=True,
            show_ticks=show_ticks,
            show_depth=show_depth,
            depth_track_count=depth_track_count,
            show_gc=show_gc,
            show_skew=show_skew,
        ),
        slot_factory=CircularTrackSlot,
    )
    return (planned if planned is not None else slots), bundle


def _prepare_circular_annotation_tracks(
    gb_record: SeqRecord,
    tracks: CircularTrackInputs,
    *,
    canvas_config: CircularCanvasConfigurator,
    show_ticks: bool,
    show_depth: bool,
    show_gc: bool,
    show_skew: bool,
    depth_track_count: int,
    record_transform: RecordDisplayTransform | None = None,
) -> CircularPlannedTracks:
    record_index = tracks.annotation_record_index
    slots, bundle = plan_circular_annotation_slots(
        tracks,
        [gb_record],
        show_ticks=show_ticks,
        show_depth=show_depth,
        show_gc=show_gc,
        show_skew=show_skew,
        depth_track_count=depth_track_count,
        record_transforms=(record_transform,) if record_transform is not None else None,
        record_indices=(record_index,),
    )
    if not bundle.set_ids and not bundle.annotations:
        return CircularPlannedTracks(slots, bundle, {}, tracks.circular_track_axis_index, record_index)
    # The planner generates the slot list whenever annotations are present.
    assert slots is not None
    # The drawing's complete list: an overlay reference the given list left
    # open (a comparison ring) is checked here (OV-2810).
    try:
        normalized = normalize_circular_track_slots_with_axis(slots, tracks.circular_track_axis_index)
    except ValueError as exc:
        raise ValidationError(
            str(exc), diagnostic={"code": "TRACK_INVALID"}
        ) from exc
    layouts = _layout_circular_annotation_tracks(
        _CircularAnnotationPacking.of(slots, bundle, canvas_config, (record_index,), len(gb_record.seq)),
        sides={slot.id: slot.side for slot in normalized},
    )
    return CircularPlannedTracks(slots, bundle, layouts, tracks.circular_track_axis_index, record_index)


_Kept = TypeVar("_Kept")


@dataclass(eq=False)
class _ResultsByInputs(Generic[_Kept]):
    """Results of one pure measure by its inputs: a hashable part picks the entries, the rest is compared by value
    (``==``), so inputs that hold mappings are kept too."""

    _entries: dict[Hashable, list[tuple[object, _Kept]]] = field(default_factory=dict)

    def get(self, bucket: Hashable, inputs: object, measure: Callable[[], _Kept]) -> _Kept:
        entries = self._entries.setdefault(bucket, [])
        for seen, result in entries:
            if seen == inputs:
                return result
        result = measure()
        entries.append((inputs, result))
        return result


@dataclass(eq=False)
class CircularAnnotationMeasures:
    """The annotation measures one Circular figure's assemblies share (OV-2664): the slots' packings, their label
    rooms, the drawn labels' places, and the slots' drawn content and paint (OV-3111). The label fit check assembles
    the figure again for each row it verifies (a smaller definition font, for one), and those assemblies pack, place
    and draw the same slots with the same inputs, and count the definition lines against the same paint; a measure
    whose inputs equal an earlier one's takes its result. Kept while ``kept_annotation_measures`` runs, in a context
    variable, so no stage key takes it: a stage's inputs still decide its result."""

    tracks: _ResultsByInputs[ResolvedAnnotationTrack] = field(default_factory=_ResultsByInputs)
    rooms: _ResultsByInputs[_AnnotationLabelRoom] = field(default_factory=_ResultsByInputs)
    labels: _ResultsByInputs[CircularAnnotationLabels] = field(default_factory=_ResultsByInputs)
    contents: _ResultsByInputs[GroupContent] = field(default_factory=_ResultsByInputs)
    inks: _ResultsByInputs[tuple[CircularInk, ...]] = field(default_factory=_ResultsByInputs)


_KEPT_MEASURES: ContextVar[CircularAnnotationMeasures | None] = ContextVar("gbdraw_circular_annotation_measures", default=None)


@contextmanager
def kept_annotation_measures(measures: CircularAnnotationMeasures) -> Iterator[None]:
    """Let the Circular assemblies run inside share ``measures`` (``CircularAnnotationMeasures``)."""

    token = _KEPT_MEASURES.set(measures)
    try:
        yield
    finally:
        _KEPT_MEASURES.reset(token)


def _kept_measure(
    table: Callable[[CircularAnnotationMeasures], _ResultsByInputs[_Kept]],
    bucket: Hashable,
    inputs: object,
    measure: Callable[[], _Kept],
) -> _Kept:
    """``measure()``, or the result kept for equal ``bucket`` and ``inputs`` while measures are kept."""

    measures = _KEPT_MEASURES.get()
    return measure() if measures is None else table(measures).get(bucket, inputs, measure)


@dataclass(frozen=True)
class _CircularAnnotationPacking:
    """What packs one circle's annotation slots on a Circular canvas: its records' annotations, its slots as
    requested (an annotation slot without a width takes its required extent), and the axis radius and label font
    the packer measures with: the canvas's as planned, which radial labels that grow it later do not change; and how
    far the labels keep from each other where they are drawn (``CircularLabelObstacles.clearance``)."""

    annotations: tuple[ResolvedRegionAnnotation, ...]
    requested: tuple[CircularTrackSlot, ...]
    record_indices: tuple[int, ...]
    record_length: int
    axis_radius_px: float
    font_size: float
    font_family: str
    label_clearance_px: float

    @property
    def record_index(self) -> int:
        """The record whose labels the slots measure: the first record on the circle."""

        return self.record_indices[0]

    @classmethod
    def of(
        cls,
        slots: Sequence[CircularTrackSlot],
        bundle: ResolvedAnnotationBundle,
        canvas_config: CircularCanvasConfigurator,
        record_indices: Sequence[int],
        length: int,
    ) -> _CircularAnnotationPacking:
        """``record_indices``: the records drawn on the circle, the first one's labels measured; ``length``: the
        circle's length in bases (``circular_display_length``)."""

        config = canvas_config.profile.config
        return cls(
            tuple(item for item in bundle.annotations if item.record_index in record_indices),
            tuple(slots),
            tuple(int(index) for index in record_indices),
            length,
            float(canvas_config.radius),
            circular_label_font_size(config, length),
            str(config.objects.text.font_family),
            float(config.labels.spacing.circular),
        )

    def track(
        self,
        slot: CircularTrackSlot,
        *,
        inside: bool,
        lanes_from_px: float,
        rooms_at_least: Sequence[LaneRoom] = (),
        drawn: bool = True,
        checks_set_size: bool = True,
    ) -> ResolvedAnnotationTrack:
        """Pack one annotation slot, its labels measured beside their marks on the slot's side
        (``circular_annotation_label_extent``), the lanes stacking from ``lanes_from_px``: the slot's edge nearest
        the axis side it starts from; each gap keeps at least its room in ``rooms_at_least``. An inside slot adds no
        room that would carry its lanes past the centre (LBL-ROOM-STOP-AT-CENTRE-OMIT): ``lanes_from_px`` is all the
        room it has, whatever its set width (``_lane_band`` takes the nearer of the two for a clipped slot). Where the
        layout ``drawn`` the slot's lanes from ``lanes_from_px``, each lane keeps its labels' letters where the drawer
        draws it (``circular_annotation_lanes_meeting``, OV-2849); the planner's packing at the axis circle, made
        before the radial layout places the slot, is not drawn there. A slot whose lanes and rooms do not fit its set
        size raises unless not ``checks_set_size`` (``layout_annotation_track``). A packing with the same inputs is
        kept (``CircularAnnotationMeasures``)."""

        room_bounds = RoomBounds(tuple(rooms_at_least), float(lanes_from_px) if inside else None)
        return _kept_measure(
            lambda measures: measures.tracks,
            (slot.id, inside, lanes_from_px, room_bounds, drawn, checks_set_size),
            (self, slot),
            partial(
                self._track, slot, inside=inside, lanes_from_px=lanes_from_px, room_bounds=room_bounds, drawn=drawn,
                checks_set_size=checks_set_size,
            ),
        )

    def _track(
        self,
        slot: CircularTrackSlot,
        *,
        inside: bool,
        lanes_from_px: float,
        room_bounds: RoomBounds,
        drawn: bool,
        checks_set_size: bool,
    ) -> ResolvedAnnotationTrack:
        params = annotation_track_params_from_mapping(slot.params)
        # The two checks of where the lanes are drawn read one placement of the labels for the same track and extent.
        placements: _ResultsByInputs[CircularLaneLabels] = _ResultsByInputs()

        def placed(track: ResolvedAnnotationTrack, track_params: AnnotationTrackParams, extent_px: float) -> CircularLaneLabels:
            return placements.get(
                extent_px,
                (track, track_params),
                partial(
                    circular_annotation_lane_labels,
                    track,
                    track_params,
                    extent_px,
                    font_family=self.font_family,
                    record_length=self.record_length,
                    lanes_from_px=lanes_from_px,
                    outward=not inside,
                    clearance_px=self.label_clearance_px,
                ),
            )

        return layout_annotation_track(
            slot.id,
            params.set_id,
            self.annotations,
            record_lengths=dict.fromkeys(self.record_indices, self.record_length),
            params=params,
            available_extent_px=slot.width.resolve(self.axis_radius_px) if slot.width is not None else None,
            labels=LabelSizing(
                font_size=self.font_size,
                extent=partial(
                    circular_annotation_label_extent,
                    font_family=self.font_family,
                    record_length=self.record_length,
                    lanes_from_px=lanes_from_px,
                    outward=not inside,
                ),
                reach=circular_annotation_mark_past_lane,
                bp_per_px=dict.fromkeys(
                    self.record_indices, self.record_length / max(1.0, 2.0 * math.pi * self.axis_radius_px)
                ),
                meet=partial(circular_annotation_lanes_meeting, placed=placed) if drawn else None,
                across_centre=partial(circular_annotation_labels_across_centre, placed=placed) if drawn and inside else None,
            ),
            join_origin=True,
            room_bounds=room_bounds,
            checks_set_size=checks_set_size,
        )


@dataclass(frozen=True)
class _AnnotationTracks:
    """One record's annotation slots: how they are packed, each slot's lanes (``layouts``), and the edge each was
    last packed from where a radial layout draws it (``packed_at``; a slot not in it was packed at the axis circle).
    The ``circular_labels`` stage takes them and packs the slots again where its radial growth draws them."""

    packing: _CircularAnnotationPacking
    layouts: Mapping[str, ResolvedAnnotationTrack]
    packed_at: Mapping[str, float]

    def label_rooms(self) -> dict[str, _AnnotationLabelRoom]:
        """The room each annotation slot's labels need beyond its band: one room, with what it measured, for each
        slot lanes equal to those of an earlier packing while annotation measures are kept (``CircularAnnotationMeasures``)."""

        packing = self.packing
        return {
            slot_id: _kept_measure(
                lambda measures: measures.rooms,
                (packing.record_index, packing.record_length, packing.font_family),
                track,
                partial(_AnnotationLabelRoom, track, packing.record_index, packing.record_length, packing.font_family),
            )
            for slot_id, track in self.layouts.items()
        }


def _layout_circular_annotation_tracks(
    packing: _CircularAnnotationPacking, *, sides: Mapping[str, str]
) -> dict[str, ResolvedAnnotationTrack]:
    """Resolve the lanes of each annotation slot, its labels measured on the slot's side (``sides``, by slot id)
    beside marks at the axis circle: no slot has its radii before the radial layout, which then packs each slot
    beside the axis (inside or outside) again where it draws it (``_annotation_tracks_packed_where_drawn``). Only
    that packing checks such a slot's set size (OV-3024); a slot not packed again is checked here."""

    return {
        slot.id: packing.track(
            slot, inside=sides.get(slot.id) == "inside", lanes_from_px=packing.axis_radius_px, drawn=False,
            checks_set_size=sides.get(slot.id) not in {"inside", "outside"},
        )
        for slot in packing.requested
        if str(slot.renderer).strip().lower() == "annotations"
    }


def _sized_annotation_slots(slots: Sequence[CircularTrackSlot], tracks: _AnnotationTracks) -> list[CircularTrackSlot]:
    """``slots`` with each annotation slot requested without a width taking its required extent."""

    requested = {slot.id: slot for slot in tracks.packing.requested}
    return [
        replace(slot, width=ScalarSpec(tracks.layouts[slot.id].required_extent_px, "px"))
        if slot.id in tracks.layouts and requested[slot.id].width is None
        else slot
        for slot in slots
    ]


# How often the annotation slots are packed again where the radial layout draws them while their lanes' edges
# move, and the move under which an edge has settled.
_LANE_EDGE_PASSES = 4
_LANE_EDGE_SETTLED_PX = 0.01


class _AnnotationTracksWhereDrawn(NamedTuple):
    """The tracks laid out with each annotation slot packed where the layout draws it."""

    fit: _CentreDefinitionFit
    slots: list[CircularTrackSlot]
    tracks: _AnnotationTracks
    label_rooms: dict[str, _AnnotationLabelRoom]


def _place_annotation_tracks_where_drawn(
    place: Callable[[Sequence[CircularTrackSlot], Mapping[str, AnnotationLabelRoom]], _CentreDefinitionFit],
    slots: Sequence[CircularTrackSlot],
    tracks: _AnnotationTracks,
) -> _AnnotationTracksWhereDrawn:
    """Lay the tracks out (``place``) with each annotation slot packed where the layout draws it, so that its
    labels are measured there (OV-2662, OV-2663).

    ``slots`` are the slots to lay out: an annotation slot requested without a width takes its required extent. A
    slot whose packing changes where it is drawn lays the tracks out again, until the slots' edges settle. A gap
    that needed room where this loop packed the slot before keeps it in later packings: a slot centred on its radius
    widens toward the axis too, and its lanes' edge moving back as the room dropped would bring the crossing back.
    A packing made for another layout (the planner's at the axis circle, or one at a circle that radial labels
    grew since) keeps no room."""

    packed_at = dict(tracks.packed_at)
    packed_here: set[str] = set()
    passes = 0
    while True:
        sized = _sized_annotation_slots(slots, tracks)
        rooms = tracks.label_rooms()
        fit = place(sized, rooms)
        repacked = (
            _annotation_tracks_packed_where_drawn(tracks, fit.layout, packed_at, packed_here)
            if passes < _LANE_EDGE_PASSES
            else {}
        )
        if not repacked:
            return _AnnotationTracksWhereDrawn(fit, sized, replace(tracks, packed_at=packed_at), rooms)
        tracks = replace(tracks, layouts={**tracks.layouts, **repacked})
        passes += 1


def _annotation_tracks_packed_where_drawn(
    tracks: _AnnotationTracks,
    radial_layout: CircularRadialLayout,
    packed_at: dict[str, float],
    packed_here: set[str],
) -> dict[str, ResolvedAnnotationTrack]:
    """The annotation slots beside the axis packed again at the edge their lanes stack from in ``radial_layout``
    (the inner edge outside the axis, the outer edge inside it), where that edge moved from where they were packed
    (``packed_at``, which this updates), or where none was packed where it is drawn yet (the planner packed it at
    the axis circle, ``_CircularAnnotationPacking.track``): only the slots whose packing changed. A slot this loop
    packed before (``packed_here``, which this updates) keeps at least the rooms of that packing."""

    requested = {slot.id: slot for slot in tracks.packing.requested}
    repacked: dict[str, ResolvedAnnotationTrack] = {}
    for resolved in radial_layout.slots:
        if resolved.id not in tracks.layouts or resolved.side not in {"inside", "outside"}:
            continue
        inside = resolved.side == "inside"
        edge = resolved.draw_outer_radius_px if inside else resolved.draw_inner_radius_px
        packed = packed_at.get(resolved.id)
        if packed is not None and abs(edge - packed) <= _LANE_EDGE_SETTLED_PX:
            continue
        rooms_at_least = tracks.layouts[resolved.id].lane_rooms if resolved.id in packed_here else ()
        packed_at[resolved.id] = edge
        packed_here.add(resolved.id)
        layout = tracks.packing.track(
            requested[resolved.id], inside=inside, lanes_from_px=edge, rooms_at_least=rooms_at_least
        )
        if layout != tracks.layouts[resolved.id]:
            repacked[resolved.id] = layout
    return repacked


def _add_circular_feature_underlays(
    gb_records: Sequence[SeqRecord],
    underlay_features: Sequence[Sequence[FeatureObject]],
    slots: list[CircularTrackSlot],
    bundle: ResolvedAnnotationBundle,
    layouts: dict[str, ResolvedAnnotationTrack],
    *,
    canvas_config: CircularCanvasConfigurator,
    record_indices: Sequence[int],
    length: int,
) -> tuple[
    list[CircularTrackSlot],
    ResolvedAnnotationBundle,
    dict[str, ResolvedAnnotationTrack],
]:
    """The underlays of the records on the circle (one, or each record of one circle), as one private set that
    covers their feature track; ``length``: the circle's length in bases."""

    merged, set_id = merge_feature_underlays(
        bundle,
        underlay_features,
        gb_records,
        mode="circular",
        record_indices=record_indices,
    )
    anchor_id = (
        feature_underlay_anchor_slot_id(slots, (feature for features in underlay_features for feature in features))
        if set_id is not None
        else None
    )
    if anchor_id is None:
        return slots, bundle, layouts
    anchor = next(slot for slot in slots if str(slot.id) == anchor_id)
    slot_id = feature_underlay_slot_id(slots)
    underlay_slot = CircularTrackSlot(
        id=slot_id,
        renderer="annotations",
        side="overlay",
        z=int(anchor.z) - 1,
        params={
            "set_id": set_id,
            "marks": ("highlight",),
            "anchor_slot": anchor_id,
            "layer": "underlay",
            "cover_anchor": True,
            "padding_px": 0.0,
            "show_labels": False,
        },
    )
    # The underlay set is new (collision-safe id), so the slots laid out before keep their layouts.
    # The underlay draws no labels.
    underlay_layouts = _layout_circular_annotation_tracks(
        _CircularAnnotationPacking.of([underlay_slot], merged, canvas_config, record_indices, length),
        sides={},
    )
    return [underlay_slot, *slots], merged, {**underlay_layouts, **layouts}


def _circular_preset_for_layout(
    cfg: GbdrawConfig,
) -> str:
    return normalize_circular_track_preset(str(cfg.canvas.circular.track_type))


def _lane_direction_for_feature_slot(
    feature_slot: CircularTrackSlot | None,
    *,
    track_preset: str,
) -> CircularFeatureLaneDirection:
    if feature_slot is not None:
        raw = feature_slot.params.get("lane_direction", feature_slot.params.get("lanes"))
        if raw is not None:
            direction = str(raw).strip().lower()
            if direction in {"inside", "outside", "split"}:
                return cast(CircularFeatureLaneDirection, direction)
    return circular_feature_lane_direction_for_preset(
        track_preset
    )


class CircularTrackSlotGeometry(TypedDict):
    """The Circular track-slot geometry of a drawing (schema 1): one entry per drawn record."""

    schema: Literal[1]
    mode: Literal["circular"]
    source: Literal["resolved"]
    records: list[CircularTrackSlotRecord]


def circular_track_slot_geometry(records: list[CircularTrackSlotRecord]) -> CircularTrackSlotGeometry:
    """The Circular track-slot geometry of a drawing with these record entries."""

    return {
        "schema": 1,
        "mode": "circular",
        "source": "resolved",
        "records": records,
    }


def _serialize_circular_track_slot_geometry(
    *,
    drawn: Sequence[tuple[int, RecordAxisInput]],
    radial_layout: CircularRadialLayout,
    layout_slots: Sequence[CircularTrackSlot],
    base_radius_px: float,
) -> CircularTrackSlotGeometry:
    """The track-slot geometry of the drawn records (index and axis input): one record, or each record of a
    circle several records share."""
    source_by_index = {
        int(getattr(slot, "slot_index", index)): slot
        for index, slot in enumerate(layout_slots)
    }
    safe_base_radius = max(FEATURE_BAND_EPSILON, float(base_radius_px))
    slots_payload: list[CircularTrackSlotEntry] = []
    for resolved in radial_layout.slots:
        source_slot = source_by_index.get(int(resolved.slot_index))
        width_spec = getattr(source_slot, "width", None)
        width_factor = None
        if width_spec is not None and getattr(width_spec, "unit", None) == "factor":
            width_factor = float(width_spec.value)
        resolved_width = (
            float(resolved.resolved_width_px)
            if resolved.resolved_width_px is not None
            else float(resolved.draw_width_px)
        )
        slots_payload.append(
            {
                "slotIndex": int(resolved.slot_index),
                "slotId": str(resolved.id),
                "renderer": str(resolved.renderer),
                "side": str(resolved.side),
                "widthPx": resolved_width,
                "widthFactor": width_factor,
                "radiusFactor": float(resolved.center_radius_px) / safe_base_radius,
                "innerGapPx": float(resolved.inner_gap_px),
                "outerGapPx": float(resolved.outer_gap_px),
                "compressed": bool(resolved.compressed),
                "source": "resolved",
            }
        )
    feature_bands: list[RadialExtent] = [
        {"innerPx": slot.draw_inner_radius_px, "outerPx": slot.draw_outer_radius_px}
        for slot in radial_layout.slots
        if slot.renderer == "features" and slot.draw_band_px is not None
    ]
    return circular_track_slot_geometry(
        [
            {
                "recordIndex": int(record_index),
                "recordId": str(getattr(axis_input.record, "id", "") or ""),
                "recordLabel": str(getattr(axis_input.record, "id", "") or ""),
                "recordKey": axis_input.record_key,
                "axisRadiusPx": float(base_radius_px),
                "axis": circular_axis_geometry(axis_input, radius_px=float(base_radius_px), feature_bands=feature_bands),
                "slots": slots_payload,
            }
            for record_index, axis_input in drawn
        ]
    )


def _legend_bbox(legend_bounds: Aabb) -> tuple[float, float, float, float]:
    """Return legend bounds as a collision tuple."""
    return (
        legend_bounds.min_x,
        legend_bounds.min_y,
        legend_bounds.max_x,
        legend_bounds.max_y,
    )


def _bbox_overlaps(a: tuple[float, float, float, float], b: tuple[float, float, float, float]) -> bool:
    """Return True when two axis-aligned bboxes overlap."""
    return a[0] < b[2] and a[2] > b[0] and a[1] < b[3] and a[3] > b[1]


def _unwrap_angle_near_reference(angle_rad: float, reference_rad: float) -> float:
    """Project angle to the nearest 2*pi branch around reference."""
    return angle_rad + (2.0 * math.pi) * round((reference_rad - angle_rad) / (2.0 * math.pi))


def _label_target_unwrapped_angle(label: CircularLabel, total_length: int) -> float:
    """Target label angle derived from genomic midpoint, in unwrapped radians."""
    return (2.0 * math.pi * (float(label["middle"]) / float(total_length))) - (0.5 * math.pi)


def _label_unwrapped_angle(label: CircularLabel, total_length: int) -> float:
    """Current label angle projected near its target branch."""
    current_angle = math.atan2(float(label["start_y"]), float(label["start_x"]))
    target_angle = _label_target_unwrapped_angle(label, total_length)
    return _unwrap_angle_near_reference(current_angle, target_angle)


def _label_bbox_local(label: CircularLabel, total_length: int, margin_px: float) -> tuple[float, float, float, float]:
    """Approximate label bbox in local (centered) coordinates."""
    authoritative_aabb = label.get("aabb_local")
    if authoritative_aabb is not None and len(authoritative_aabb) == 4:
        half_margin = 0.5 * float(margin_px)
        return (
            float(authoritative_aabb[0]) - half_margin,
            float(authoritative_aabb[1]) - half_margin,
            float(authoritative_aabb[2]) + half_margin,
            float(authoritative_aabb[3]) + half_margin,
        )
    half_margin = 0.5 * margin_px
    start_x = float(label["start_x"])
    start_y = float(label["start_y"])
    width_px = float(label["width_px"])
    height_px = float(label["height_px"])
    is_inner = bool(label.get("is_inner", False))

    if not is_inner:
        if start_x > 0:
            min_x = start_x - half_margin
            max_x = start_x + width_px + half_margin
        else:
            min_x = start_x - width_px - half_margin
            max_x = start_x + half_margin
    else:
        if start_x > 0:
            min_x = start_x - width_px - half_margin
            max_x = start_x + half_margin
        else:
            min_x = start_x - half_margin
            max_x = start_x + width_px + half_margin

    angle_deg = (360.0 * (float(label["middle"]) / float(total_length))) % 360.0
    if not is_inner:
        if 0.0 <= angle_deg < 10.0 or angle_deg >= 350.0:
            max_y = start_y
            min_y = start_y - height_px - half_margin
        elif 170.0 <= angle_deg < 190.0:
            max_y = start_y + height_px + half_margin
            min_y = start_y - half_margin
        else:
            max_y = start_y + 0.5 * height_px + half_margin
            min_y = start_y - 0.5 * height_px - half_margin
    else:
        if 0.0 <= angle_deg < 10.0:
            max_y = start_y + height_px + half_margin
            min_y = start_y
        elif angle_deg >= 350.0:
            max_y = start_y
            min_y = start_y - height_px - half_margin
        else:
            max_y = start_y + 0.5 * height_px + half_margin
            min_y = start_y - 0.5 * height_px - half_margin
    return min_x, min_y, max_x, max_y


def _label_bbox_on_canvas(
    label: CircularLabel,
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    margin_px: float,
) -> tuple[float, float, float, float]:
    """Return label bbox on canvas as (min_x, min_y, max_x, max_y)."""
    min_x, min_y, max_x, max_y = _label_bbox_local(label, total_length, margin_px)
    return (
        min_x + float(canvas_config.offset_x),
        min_y + float(canvas_config.offset_y),
        max_x + float(canvas_config.offset_x),
        max_y + float(canvas_config.offset_y),
    )


def _move_label_text_anchor(
    label: CircularLabel,
    *,
    start_x: float,
    start_y: float,
) -> None:
    """Move label text geometry while keeping authoritative bounds coherent."""
    old_x = float(label.get("start_x", 0.0))
    old_y = float(label.get("start_y", 0.0))
    dx = float(start_x) - old_x
    dy = float(start_y) - old_y
    label["start_x"] = float(start_x)
    label["start_y"] = float(start_y)

    if "text_x" in label:
        label["text_x"] = float(label["text_x"]) + dx
    if "text_y" in label:
        label["text_y"] = float(label["text_y"]) + dy
    if "leader_start_x" in label:
        label["leader_start_x"] = float(label["leader_start_x"]) + dx
    if "leader_start_y" in label:
        label["leader_start_y"] = float(label["leader_start_y"]) + dy

    authoritative_aabb = label.get("aabb_local")
    if authoritative_aabb is not None and len(authoritative_aabb) == 4:
        label["aabb_local"] = (
            float(authoritative_aabb[0]) + dx,
            float(authoritative_aabb[1]) + dy,
            float(authoritative_aabb[2]) + dx,
            float(authoritative_aabb[3]) + dy,
        )
    oriented_corners = label.get("oriented_corners")
    if oriented_corners is not None:
        label["oriented_corners"] = tuple(
            (float(point[0]) + dx, float(point[1]) + dy)
            for point in oriented_corners
        )


def _external_label_bounds_on_canvas(
    labels: list[CircularLabel],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    *,
    margin_px: float = LEGEND_LABEL_MARGIN_PX,
) -> tuple[float, float, float, float] | None:
    """Return bounds of all non-embedded labels in canvas coordinates."""
    min_x: float | None = None
    min_y: float | None = None
    max_x: float | None = None
    max_y: float | None = None

    for label in labels:
        if label.get("is_embedded"):
            continue
        box_min_x, box_min_y, box_max_x, box_max_y = _label_bbox_on_canvas(
            label,
            total_length,
            canvas_config,
            margin_px=margin_px,
        )
        min_x = box_min_x if min_x is None else min(min_x, box_min_x)
        min_y = box_min_y if min_y is None else min(min_y, box_min_y)
        max_x = box_max_x if max_x is None else max(max_x, box_max_x)
        max_y = box_max_y if max_y is None else max(max_y, box_max_y)

    if min_x is None or min_y is None or max_x is None or max_y is None:
        return None
    return min_x, min_y, max_x, max_y


def _external_label_paint_bounds_local(label: CircularLabel, total_length: int, style: FeatureLabelStyle) -> Aabb:
    """Return what an outer or inner label paints in ``style`` (its turned letters, its underline or white box,
    and its leader at its label's own width) in circle-local coordinates, from the label's ink owner."""
    ink = circular_label_ink(label, total_length, style)
    box = ink.ink
    letters = tuple(
        ink.from_label_frame(corner)
        for corner in ((box.min_x, box.min_y), (box.max_x, box.min_y), (box.max_x, box.max_y), (box.min_x, box.max_y))
    )
    drawn = tuple(point for part in feature_label_drawn_parts(ink, style) for point in part)
    text_bounds = bounding_aabb((*letters, *drawn))
    leader_bounds = bounding_aabb(ink.leader, 0.5 * label_leader_width_px(label))
    if (
        text_bounds is None or leader_bounds is None or (combined := union_aabbs((text_bounds, leader_bounds))) is None
    ):  # pragma: no cover - placed labels have finite coordinates
        raise RuntimeError("external label paint bounds unexpectedly empty")
    return combined


def _external_label_obstacles_on_canvas(
    labels: Sequence[CircularLabel],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
) -> tuple[Aabb, ...]:
    """Return authoritative per-label overlay obstacles in plot coordinates."""
    style = canvas_config.profile.label_style
    return tuple(
        _external_label_paint_bounds_local(label, total_length, style).translated(
            float(canvas_config.offset_x), float(canvas_config.offset_y)
        )
        for label in labels
        if not label.get("is_embedded")
    )


# The drawn labels of each annotation slot of a record, by slot id.
_SlotLabels: TypeAlias = dict[str, tuple[CircularAnnotationLabel, ...]]


# A drawn region label of a record: its slot id and its annotation's id.
_RegionLabelKey: TypeAlias = tuple[str, str]


class _RegionLabelsNotClear(NamedTuple):
    """What a record's annotation slots could not draw clear (``CircularAnnotationLabels``): the labels left out
    because no letter fits their circle, and what each drawn label still overlaps (OV-3015); with the letters each
    drawn label keeps clear and their height, which the feature labels clear when they make way (CIR-RGN-YIELD).
    Both are keyed by slot id and annotation id, in slot order, then in each slot's order."""

    without_letters: tuple[CircularAnnotationLabel, ...] = ()
    overlapped: Mapping[_RegionLabelKey, frozenset[int]] = MappingProxyType({})
    letters: Mapping[_RegionLabelKey, BoxedLetters] = MappingProxyType({})
    # The places of the drawn labels whose letters meet a feature label's letters, by slot id (OV-3017).
    over_feature_letters: Mapping[str, tuple[RegionLabelPlace, ...]] = MappingProxyType({})

    @property
    def overlaps(self) -> int:
        """How many things the drawn labels still overlap, as the region label warning counts them."""

        return sum(len(met) for met in self.overlapped.values())


# The places along their regions that the feature labels made way for, by slot id (CIR-RGN-YIELD).
_RegionLabelsMadeWay: TypeAlias = Mapping[str, tuple[RegionLabelPlace, ...]]


def _circular_annotation_labels_by_slot(
    radial_layout: CircularRadialLayout,
    annotation_track_layouts: Mapping[str, ResolvedAnnotationTrack] | None,
    *,
    record_index: int,
    record_length: int,
    font_family: str,
    obstacles: CircularLabelObstacles,
    made_way: _RegionLabelsMadeWay,
) -> tuple[_SlotLabels, _RegionLabelsNotClear]:
    """The drawn labels of each annotation slot of a record, measured once for the drawer and the canvas bounds,
    and, for the label fit warning, the labels the slots leave out because no letter fits their circle and what the
    drawn labels still overlap (``CircularAnnotationLabels``).

    The letters clear the axis circle's paint too where they would meet it, and
    move along their regions off ``obstacles`` and the labels of earlier slots,
    except those drawn where the feature labels made way for them
    (``made_way``, by slot id).
    """

    axis_paint = radial_layout.axis.paint_band_px
    labels: dict[str, tuple[CircularAnnotationLabel, ...]] = {}
    without_letters: list[CircularAnnotationLabel] = []
    over_feature_letters: dict[str, tuple[RegionLabelPlace, ...]] = {}
    overlapped: dict[_RegionLabelKey, frozenset[int]] = {}
    letters: dict[_RegionLabelKey, BoxedLetters] = {}
    for slot in radial_layout.slots:
        layout = (annotation_track_layouts or {}).get(slot.id)
        if str(slot.renderer) != "annotations" or layout is None:
            continue
        inner, outer, side, params = (
            float(slot.draw_inner_radius_px),
            float(slot.draw_outer_radius_px),
            str(slot.side),
            annotation_track_params_from_mapping(slot.params),
        )
        obstacles = replace(obstacles, made_way=made_way.get(slot.id, ()))
        rings = _label_rings(axis_paint, slot.label_clear_band_px)
        measured = _kept_measure(
            lambda measures: measures.labels,
            (record_index, record_length, inner, outer, side, font_family, slot.spans_slots),
            (layout, params, rings, obstacles),
            partial(
                circular_annotation_labels,
                layout,
                record_index=record_index,
                record_length=record_length,
                inner_radius_px=inner,
                outer_radius_px=outer,
                side=side,
                font_family=font_family,
                params=params,
                rings=rings,
                obstacles=obstacles,
                covers_slots=slot.spans_slots,
            ),
        )
        labels[slot.id] = measured.drawn
        without_letters.extend(measured.without_letters)
        if measured.over_feature_letters:
            over_feature_letters[slot.id] = measured.over_feature_letters
        overlapped.update(((slot.id, annotation_id), met) for annotation_id, met in measured.overlapped.items())
        letters.update(((slot.id, annotation_id), drawn) for annotation_id, drawn in measured.letters.items())
        obstacles = obstacles.with_labels(measured.drawn)
    return labels, _RegionLabelsNotClear(tuple(without_letters), overlapped, letters, over_feature_letters)


def _feature_labels_making_way(
    labels: Sequence[CircularLabel],
    not_clear: _RegionLabelsNotClear,
    total_length: int,
    definition_text: Sequence[Aabb],
    style: FeatureLabelStyle,
    placed: Callable[[Sequence[CircularLabel], _RegionLabelsMadeWay], _RegionLabelsNotClear],
) -> tuple[list[CircularLabel], _RegionLabelsMadeWay] | None:
    """The feature labels moved off the letters of each region label drawn over them, and the places they made way
    for, by slot id (CIR-RGN-YIELD, OV-3017); None when they clear none. ``not_clear`` is what ``placed`` drew with
    the feature labels at ``labels``. The feature label repair moves them (``make_way_for_region_labels``, the one
    owner of their positions) off every drawn region label's letters, and keeps no move that raises what a region
    label overlaps where ``placed`` then draws it, as its warning counts it; a region label it cannot clear is drawn
    and warns (CIR-REGION-WARN). Every region label is found by its slot id and annotation id."""

    keys = list(not_clear.letters)
    places = {
        (slot_id, place.annotation_id): place for slot_id, over in not_clear.over_feature_letters.items() for place in over
    }
    blocked = [at for at, key in enumerate(keys) if key in places]
    if not labels or not blocked:
        return None

    def made_way_at(cleared: Sequence[int]) -> _RegionLabelsMadeWay:
        made_way: dict[str, tuple[RegionLabelPlace, ...]] = {}
        for slot_id, annotation_id in (keys[at] for at in cleared):
            made_way[slot_id] = (*made_way.get(slot_id, ()), places[slot_id, annotation_id])
        return made_way

    def region_overlaps(moved: Sequence[CircularLabel], cleared: Sequence[int]) -> tuple[frozenset[int], ...]:
        overlapped = placed(moved, made_way_at(cleared)).overlapped
        return tuple(overlapped[key] for key in keys)

    made_way = make_way_for_region_labels(
        labels,
        total_length,
        [not_clear.letters[key] for key in keys],
        blocked,
        overlapped=[not_clear.overlapped[key] for key in keys],
        region_overlaps=region_overlaps,
        definition_text=definition_text,
        drawn_parts=lambda label: _feature_label_drawn_parts(label, total_length, style),
    )
    if made_way is None:
        return None
    moved, cleared = made_way
    return moved, made_way_at(cleared)


@dataclass(eq=False)
class _RecordRegionLabels:
    """The annotation labels of one record, placed off its feature labels and leaders, and the places along their
    regions that the feature labels made way for (CIR-RGN-YIELD): the one owner of those places across the record's
    label passes, so every later placement draws those labels where the feature labels made way."""

    radial_layout: CircularRadialLayout
    annotation_track_layouts: Mapping[str, ResolvedAnnotationTrack] | None
    record_index: int
    record_length: int
    font_family: str
    canvas_config: CircularCanvasConfigurator
    made_way: _RegionLabelsMadeWay = field(default_factory=dict)

    def _placed(
        self, feature_labels: Sequence[CircularLabel], made_way: _RegionLabelsMadeWay
    ) -> tuple[_SlotLabels, _RegionLabelsNotClear]:
        """The drawn labels of each slot with the feature labels at ``feature_labels`` and those of ``made_way`` drawn
        where the feature labels made way for them (``_circular_annotation_labels_by_slot``)."""

        return _circular_annotation_labels_by_slot(
            self.radial_layout,
            self.annotation_track_layouts,
            record_index=self.record_index,
            record_length=self.record_length,
            font_family=self.font_family,
            obstacles=_feature_label_obstacles(feature_labels, self.record_length, self.canvas_config),
            made_way=made_way,
        )

    def place(self, feature_labels: Sequence[CircularLabel]) -> tuple[_SlotLabels, tuple[Aabb, ...], _RegionLabelsNotClear]:
        """The drawn labels of each slot, their ink bounds on the canvas (a legend or title avoids them), and what
        they could not draw clear (``_circular_annotation_labels_by_slot``)."""

        labels_by_slot, not_clear = self._placed(feature_labels, self.made_way)
        offset_x, offset_y = float(self.canvas_config.offset_x), float(self.canvas_config.offset_y)
        return labels_by_slot, tuple(
            label.ink_bounds.translated(offset_x, offset_y) for labels in labels_by_slot.values() for label in labels
        ), not_clear

    def make_way(
        self, feature_labels: Sequence[CircularLabel], not_clear: _RegionLabelsNotClear, definition_text: Sequence[Aabb]
    ) -> list[CircularLabel] | None:
        """The feature labels moved off the region labels that ``not_clear`` found still over their letters, keeping
        the places they made way for; None when they clear none (``_feature_labels_making_way``, which judges each
        move by where the next ``place`` draws the region labels)."""

        made_way = _feature_labels_making_way(
            feature_labels,
            not_clear,
            self.record_length,
            definition_text,
            self.canvas_config.profile.label_style,
            lambda moved, made_way: self._placed(moved, made_way)[1],
        )
        if made_way is None:
            return None
        moved, self.made_way = made_way
        return moved


def _region_label_warnings(
    radial_layout: CircularRadialLayout,
    tracks: Mapping[str, ResolvedAnnotationTrack],
    record_index: int,
    labels: Mapping[str, Sequence[CircularAnnotationLabel]],
    not_clear: _RegionLabelsNotClear,
) -> RegionLabelWarnings:
    """The record's region label counts that only warn, with the font sizes of the labels left out: the labels its
    tracks leave out for want of room (``region_rooms_left_out``), what the drawn ``labels`` still overlap (OV-3015,
    REGION-WARN as in Linear), the drawn labels that lie on another track (LP-D86), and the labels left out because
    no letter fits their circle (RGN-CUT-EMPTY-OMIT)."""

    return region_rooms_left_out(tracks.values(), record_index).adding(
        LabelFit(
            0,
            0,
            region_overlaps=not_clear.overlaps,
            region_track_overlaps=_region_labels_on_other_tracks(radial_layout, labels),
            region_cuts_left_out=len(not_clear.without_letters),
        ),
        tuple(label.font_size for label in not_clear.without_letters),
    )


# Letters lie on a track when they reach more than this into its band, as a label covers the definition text.
REGION_TRACK_MIN_OVERLAP_PX = 0.5


def _region_labels_on_other_tracks(
    radial_layout: CircularRadialLayout, labels: Mapping[str, Sequence[CircularAnnotationLabel]]
) -> int:
    """How many annotation labels' letters lie on the band another track takes.

    The layout keeps a room for these letters beyond an inside or outside
    annotation slot and beyond a span overlay (OV-1192), so they lie on a
    track only where set radii leave no room (OV-1207); the figure is then
    drawn and warns (LP-D86). A track takes the band the layout keeps for it
    (``own_reserved_band_px``): its drawing with its tick or axis labels and
    its own labels' room. Annotation overlays draw over other tracks and take
    none. An overlay the layout keeps no room for (``span_keeps_label_room``:
    a span whose outer edge is a plot overlay) has its labels drawn over the
    tracks around it whatever the radii, so they are not counted.
    """

    tracks = _track_bands(radial_layout)
    count = 0
    for slot in radial_layout.slots:
        own = slot.draw_band_px
        if own is None or (slot.side == "overlay" and not slot.span_keeps_label_room):
            continue
        others = [
            band for track, band in tracks if track is not slot and not _reaches_band(own.inner_px, own.outer_px, band)
        ]
        count += sum(
            1 for label in labels.get(slot.id, ()) if any(_reaches_band(*label.ink_radii, band) for band in others)
        )
    return count


def _track_bands(radial_layout: CircularRadialLayout) -> list[tuple[CircularResolvedSlot, RadialBand]]:
    """The band each drawn track takes (``own_reserved_band_px``: its drawing with its tick or axis labels and its
    own labels' room); annotation overlays draw over other tracks and take none."""

    return [
        (slot, band)
        for slot in radial_layout.slots
        if slot.draw_band_px is not None
        and not (slot.renderer == "annotations" and slot.side == "overlay")
        and (band := slot.own_reserved_band_px) is not None
    ]


def _reaches_band(inner: float, outer: float, band: RadialBand) -> bool:
    """Whether ink between radii ``inner`` and ``outer`` reaches more than the threshold into ``band``."""

    return min(outer - band.inner_px, band.outer_px - inner) > REGION_TRACK_MIN_OVERLAP_PX


# What one drawn GC content, GC skew or depth slot paints with: its data and its settings.
SlotPaintData: TypeAlias = (
    tuple[DataFrame, GcContentConfigurator] | tuple[DataFrame, GcSkewConfigurator] | tuple[DataFrame, DepthConfigurator]
)


@dataclass(frozen=True, eq=False)
class CircularTrackInkSources:
    """What a record's tracks paint from, besides the radial layout: the sources the drawing reads for each slot
    (``_draw_resolved_circular_slot``). The ``definition_track_overlaps`` stage is keyed by these fields."""

    gb_record: SeqRecord
    geometry: CircularCanvasGeometry
    profile: CircularRenderProfile
    record_transform: RecordDisplayTransform | None
    # Each drawn GC content, GC skew and depth slot's data and settings, by slot id (``_slot_paint_data``): only
    # what a drawn slot paints with, not a depth track's label.
    slot_data: Mapping[str, SlotPaintData]
    # The record's feature density layers by track id: a heatmap slot paints its runs.
    feature_density: Mapping[str, FeatureDensityLayer] | None
    conservation_tracks: Sequence[ConservationTrack] | None
    conservation_min_identity: float
    # The feature layers the radial layout measured (``_radial_layout_feature_dict``).
    feature_layers: FeatureBuildResult
    # Half the feature blocks' stroke (``stroke_pad_px``): not their arrow shape.
    feature_stroke_pad_px: float
    annotation_track_layouts: Mapping[str, ResolvedAnnotationTrack] | None
    annotation_record_index: int
    use_slot_tick_options: bool


def _record_track_ink(
    axis: CircularAxisLayout,
    slots: Sequence[CircularResolvedSlot],
    sources: CircularTrackInkSources,
    annotation_labels: Mapping[str, Sequence[CircularAnnotationLabel]],
) -> tuple[CircularInk, ...]:
    """What the record's tracks paint, from the owners that draw them (D-18), read with the same slot sources
    the drawing reads (``_draw_resolved_circular_slot``).

    It reads the radial layout's axis and resolved slots only, not the centre band or the outer content radius
    (OV-2111).

    The GC content, GC skew and depth areas and their axes come from their drawers, the feature blocks from the
    intervals the labels clear (``feature_radius_intervals``), the conservation hits and background from the
    ring's owner, the annotation marks and labels from the annotation track's owner, the tick marks and labels
    from the tick group, and the axis circle from its paint band. A spacer paints nothing.
    """

    gb_record = sources.gb_record
    canvas_radius = float(sources.geometry.radius)
    record_length = circular_display_length(gb_record, sources.record_transform)
    cfg = sources.profile.config
    record_transform = (
        sources.record_transform if sources.record_transform is not None else unset_record_transform(gb_record)
    )
    inks: list[CircularInk] = [CircularInkRing(axis.paint_band_px)]
    for slot in slots:
        renderer = str(slot.renderer)
        if slot.draw_band_px is None:
            continue
        norm_factor = _slot_norm_factor(slot, canvas_radius)
        slot_data = sources.slot_data.get(str(slot.id))
        if slot_data is not None:
            data_df, data_config = slot_data
            if renderer == "dinucleotide_content" and isinstance(data_config, GcContentConfigurator):
                inks.extend(
                    GcContentDrawer(data_config).ink(
                        canvas_radius, data_df, record_length, slot.draw_width_px, norm_factor,
                        data_config.dinucleotide, record_transform,
                    )
                )
            elif renderer == "dinucleotide_skew" and isinstance(data_config, GcSkewConfigurator):
                inks.extend(
                    SkewDrawer(data_config).ink(
                        canvas_radius, data_df, record_length, slot.draw_width_px, norm_factor,
                        data_config.dinucleotide, record_transform,
                    )
                )
            elif renderer in {"depth", "feature_density"} and isinstance(data_config, DepthConfigurator):
                inks.extend(
                    DepthDrawer(data_config).ink(
                        canvas_radius, data_df, record_length, float(slot.draw_width_px), norm_factor,
                        record_transform=record_transform,
                    )
                )
        elif renderer == "feature_density":
            layer = slot_feature_density_layer(str(slot.id), slot, sources.feature_density)
            if layer is not None:
                inks.extend(feature_density_heatmap_ink(layer, record_length, slot.draw_band_px, record_transform))
        elif renderer == "sequence_conservation":
            track = _slot_conservation_track(slot, sources.conservation_tracks)
            if track is not None:
                inks.extend(
                    conservation_ring_ink(
                        track,
                        total_length=record_length,
                        inner_radius_px=float(slot.draw_inner_radius_px),
                        outer_radius_px=float(slot.draw_outer_radius_px),
                        min_identity=float(sources.conservation_min_identity),
                        cfg=cfg,
                        record_transform=record_transform,
                    )
                )
        elif renderer == "annotations":
            annotation_layout = (sources.annotation_track_layouts or {}).get(slot.id)
            if annotation_layout is not None:
                ring = _annotation_ring(slot)
                params = annotation_track_params_from_mapping(slot.params)
                labels = tuple(annotation_labels.get(slot.id, ()))
                painted: tuple[CircularInk, ...] = _kept_measure(
                    lambda measures: measures.inks,
                    (slot.id, sources.annotation_record_index, record_length),
                    (annotation_layout, ring, params, labels),
                    partial(
                        circular_annotation_track_ink,
                        annotation_layout,
                        record_index=sources.annotation_record_index,
                        record_length=record_length,
                        ring=ring,
                        params=params,
                        labels=labels,
                    ),
                )
                inks.extend(painted)
        elif renderer == "ticks":
            tick_options = _slot_tick_options(slot, sources.use_slot_tick_options)
            if tick_options is not None:
                inks.extend(
                    TickGroup(
                        gb_record,
                        sources.geometry,
                        profile=sources.profile,
                        radius=float(slot.anchor_radius_px or 0.0),
                        record_transform=record_transform,
                        slot_id=str(slot.id),
                        **tick_options,
                    ).ink()
                )
        elif renderer == "features" and isinstance(slot.payload, CircularFeatureLayout):
            pad = sources.feature_stroke_pad_px
            inks.extend(
                CircularInkRing(RadialBand(inner - pad, outer + pad), start / record_length, end / record_length)
                for start, end, inner, outer in feature_radius_intervals(
                    _radial_layout_feature_dict(sources.feature_layers), record_length, slot.payload
                )
            )
    return tuple(inks)


def definition_lines_over_tracks(
    definition_lines: tuple[Aabb, ...],
    axis: CircularAxisLayout,
    slots: tuple[CircularResolvedSlot, ...],
    sources: CircularTrackInkSources,
    annotation_labels: Mapping[str, Sequence[CircularAnnotationLabel]],
) -> tuple[int, ...]:
    """Which centre definition lines (their places in ``definition_lines``) meet a track's paint of the record
    ``sources`` paints (stage ``definition_track_overlaps``, OV-2088); on a circle several records share, the lines
    that meet any record's paint are those of the records together.

    Every call builds every track's paint and counts against all of it; the stage key (the lines, the radial
    layout's axis and slots, the paint sources and the annotation labels) lets a Generate that changes none of
    them reuse the count. The centre band a species change moves beside a plot title (OV-1710) is not an input:
    the paint does not read it (OV-2111).
    """

    return _definition_lines_meeting(_record_track_ink(axis, slots, sources, annotation_labels), definition_lines)


def _definition_lines_meeting(inks: Sequence[CircularInk], definition_text: Sequence[Aabb]) -> tuple[int, ...]:
    """The places of the lines ``_definition_lines_on_tracks`` counts."""

    return tuple(
        index for index, box in enumerate(definition_text)
        if box_meets_ink(box, inks, threshold_px=REGION_TRACK_MIN_OVERLAP_PX)
    )


def _definition_lines_on_tracks(inks: Sequence[CircularInk], definition_text: Sequence[Aabb]) -> int:
    """How many lines of the centre definition text meet what a track paints (D-18).

    ``definition_text`` holds the lines' ink boxes centred on the record
    (``_centre_definition_ink``, the boxes the centre reserve holds); a line
    counts when its box reaches more than
    the threshold that region labels use (LP-D86) into a track's paint
    (``_record_track_ink``). The figure is drawn as before and warns.
    """

    return len(_definition_lines_meeting(inks, definition_text))


def _feature_label_drawn_parts(label: CircularLabel, total_length: int, style: FeatureLabelStyle) -> tuple[Polygon, ...]:
    """What ``style`` draws with a feature label besides its letters, its underline and white box
    (``feature_label_drawn_parts``); nothing for an embedded label, or when the style draws neither (its letters'
    ink is not measured then)."""

    if label.get("is_embedded") or not style.draws_from_ink:
        return ()
    return feature_label_drawn_parts(circular_label_ink(label, total_length, style), style)


def _leader_segments(label: CircularLabel, shape: LabelShape) -> tuple[LeaderSegment, ...]:
    """The straight parts of a feature label's drawn leader, each with half its label's leader width."""

    return tuple(
        LeaderSegment(start, end, 0.5 * label_leader_width_px(label)) for start, end in zip(shape.leader, shape.leader[1:])
    )


def _feature_label_leaders(labels: Sequence[CircularLabel], total_length: int) -> tuple[LeaderSegment, ...]:
    """The record's feature label leaders, drawn before every slot: the scale backs the numbers they run behind
    (OV-3301)."""

    return tuple(
        segment
        for label in labels
        if not label.get("is_embedded")
        for segment in _leader_segments(label, circular_label_shape(label, total_length))
    )


def _feature_label_obstacles(
    labels: Sequence[CircularLabel], total_length: int, canvas_config: CircularCanvasConfigurator
) -> CircularLabelObstacles:
    """The record's drawn feature labels and their leaders, as the fit count sees them, for annotation labels to
    keep the label spacing from (OV-1173). A label with a leader also covers what its style draws with its letters,
    its underline and its white box (OV-1358, the Linear rule of OV-1357)."""

    style = canvas_config.profile.label_style
    shapes = [
        (label, shape, _feature_label_drawn_parts(label, total_length, style))
        for label in labels
        if (shape := circular_label_shape(label, total_length)).parts
    ]
    return CircularLabelObstacles(
        labels=tuple(
            CircularObstacleLabel(
                (*shape.parts, *drawn),
                _leader_segments(label, shape),
                feature=True,
            )
            for label, shape, drawn in shapes
        ),
        clearance=float(canvas_config.profile.config.labels.spacing.circular),
    )


class _AnnotationLabelRoom:
    """The room one annotation slot's labels need beyond a band (``AnnotationLabelRoom``), and whether the band
    leaves out a label the slot's lanes draw (``drops_labels``), measured by the labels'
    own owner once for each band, axis ring, side, parameters and ``covers_slots`` flag it is asked for: the layouts
    of one drawing, the definition font scan's among them, ask for the same ones again and again (OV-2664). The
    memo is reused for equal lanes while annotation measures are kept (``CircularAnnotationMeasures``): within one
    assembly and across the assemblies of one figure. It is no stage input: a stage takes the tracks (``_AnnotationTracks``) and
    makes its own rooms."""

    __slots__ = ("track", "_record_index", "_record_length", "_font_family", "_measured")

    def __init__(self, track: ResolvedAnnotationTrack, record_index: int, record_length: int, font_family: str) -> None:
        self.track = track
        self._record_index = record_index
        self._record_length = record_length
        self._font_family = font_family
        self._measured: dict[
            tuple[RadialBand, tuple[RadialBand, ...], str, AnnotationTrackParams, bool], CircularLabelRoom
        ] = {}

    def __call__(
        self,
        band: RadialBand,
        axis: CircularAxisLayout,
        side: str,
        params: AnnotationTrackParams,
        covers_slots: bool,
        clears: RadialBand | None = None,
    ) -> float:
        return self._room(band, _label_rings(axis.paint_band_px, clears), side, params, covers_slots).px

    def drops_labels(
        self, band: RadialBand, axis: CircularAxisLayout, side: str, params: AnnotationTrackParams
    ) -> bool:
        return self._room(band, _label_rings(axis.paint_band_px, None), side, params, False).drops_labels

    def _room(
        self,
        band: RadialBand,
        rings: tuple[RadialBand, ...],
        side: str,
        params: AnnotationTrackParams,
        covers_slots: bool,
    ) -> CircularLabelRoom:
        key = (band, rings, side, params, covers_slots)
        room = self._measured.get(key)
        if room is None:
            room = self._measured[key] = circular_annotation_label_room(
                self.track,
                record_index=self._record_index,
                record_length=self._record_length,
                band=band,
                side=side,
                font_family=self._font_family,
                params=params,
                rings=rings,
                covers_slots=covers_slots,
            )
        return room


def _label_rings(axis_paint: RadialBand, clears: RadialBand | None) -> tuple[RadialBand, ...]:
    """The rings an annotation slot's letters keep ``label_offset`` beyond where they would meet them: the axis
    circle's paint, and the lanes of the features slot an overlay is anchored to when its labels keep off them
    (``CircularResolvedSlot.label_clear_band_px``, OV-3220)."""

    return (axis_paint,) if clears is None else (axis_paint, clears)


def _primary_radial_half_stroke_px(
    *,
    canvas_config: CircularCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    depth_config: DepthConfigurator | None,
) -> float:
    """Return conservative stroke safety for radial record/track geometry."""
    cfg = canvas_config.profile.config
    widths = [
        float(
            cfg.objects.axis.circular.stroke_width.for_length_param(
                canvas_config.length_param
            )
        ),
        float(feature_config.block_stroke_width),
        float(feature_config.line_stroke_width),
        float(gc_config.stroke_width),
        float(skew_config.stroke_width),
    ]
    if depth_config is not None:
        widths.append(float(depth_config.stroke_width))
    conservation = getattr(cfg.objects, "conservation", None)
    if conservation is not None:
        widths.append(float(getattr(conservation, "stroke_width", 0.0)))
    return 0.5 * max(widths, default=0.0)


def _collect_circular_primary_bounds(
    *,
    radial_layout: CircularRadialLayout,
    canvas_config: CircularCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    depth_config: DepthConfigurator | None,
    overlay_obstacles: Sequence[Aabb],
    definition_target: Group | None = None,
) -> Aabb:
    """Collect authoritative plot-space bounds before outer composition: the rings, and past them only the
    ``overlay_obstacles`` (the labels as drawn), not the room the slots keep for region labels (OV-2809)."""
    stroke_safety = _primary_radial_half_stroke_px(
        canvas_config=canvas_config,
        feature_config=feature_config,
        gc_config=gc_config,
        skew_config=skew_config,
        depth_config=depth_config,
    )
    radius = float(radial_layout.ring_outer_radius_px) + stroke_safety
    center_x = float(canvas_config.offset_x)
    center_y = float(canvas_config.offset_y)
    bounds: list[Aabb] = [
        Aabb(
            center_x - radius,
            center_y - radius,
            center_x + radius,
            center_y + radius,
        )
    ]
    bounds.extend(overlay_obstacles)
    if definition_target is not None:
        local_bounds = getattr(definition_target, "_gbdraw_local_bounds", None)
        translation = getattr(definition_target, "_gbdraw_plot_translation", None)
        if isinstance(local_bounds, Aabb) and translation is not None:
            bounds.append(
                local_bounds.translated(float(translation[0]), float(translation[1]))
            )
    combined = union_aabbs(bounds)
    if combined is None:  # pragma: no cover - radial bounds are always present
        raise RuntimeError("circular primary bounds unexpectedly empty")
    return combined


def _overlay_legend_bounds_for_collision(
    *,
    legend_measurement: LegendMeasurement,
    primary_bounds: Aabb,
    placement: LegendPlacement,
    record_obstacles: Sequence[Aabb],
) -> Aabb:
    """Return where the composer puts a corner legend before label adjustment.

    ``record_obstacles`` are the obstacles of the final plan that labels cannot
    clear (the record and its tracks), so the predicted corner is one the final
    plan can take; labels are left out because they move off the prediction
    (OV-1609).
    """
    if not placement.is_overlay:
        raise ValueError("collision bounds require an overlay legend placement")
    plan = render_stage(
        "composition_plan",
        plan_composition,
        CompositionRequest(
            primary=CompositionItem("primary", primary_bounds),
            legend=CompositionItem("legend", legend_measurement.local_bounds),
            legend_placement=placement,
            overlay_obstacles=tuple(record_obstacles),
        ),
    )
    primary = plan.placement_for("primary")
    legend = plan.placement_for("legend")
    if primary is None or legend is None:  # pragma: no cover - request has both items
        raise RuntimeError("overlay collision plan is incomplete")
    return legend.final_bounds.translated(-primary.dx, -primary.dy)


def _legend_collision_indices(
    labels: list[CircularLabel],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    legend_bounds: Aabb,
    margin_px: float = LEGEND_LABEL_MARGIN_PX,
) -> list[int]:
    """Return indices of labels that currently collide with the legend."""
    legend_box = _legend_bbox(legend_bounds)
    collided: list[int] = []
    for idx, label in enumerate(labels):
        if label.get("is_embedded"):
            continue
        label_box = _label_bbox_on_canvas(label, total_length, canvas_config, margin_px=margin_px)
        if _bbox_overlaps(label_box, legend_box):
            collided.append(idx)
    return collided


def _labels_collide_with_legend(
    labels: list[CircularLabel],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    legend_bounds: Aabb,
) -> bool:
    """Whether any external label overlaps with the legend bbox."""
    return bool(
        _legend_collision_indices(
            labels,
            total_length,
            canvas_config,
            legend_bounds,
        )
    )




def _legend_center_local(
    canvas_config: CircularCanvasConfigurator,
    legend_bounds: Aabb,
) -> tuple[float, float]:
    """Return legend center in local (circle-centered) coordinates."""
    legend_min_x, legend_min_y, legend_max_x, legend_max_y = _legend_bbox(legend_bounds)
    legend_center_x = 0.5 * (legend_min_x + legend_max_x)
    legend_center_y = 0.5 * (legend_min_y + legend_max_y)
    return legend_center_x - float(canvas_config.offset_x), legend_center_y - float(canvas_config.offset_y)


def _preferred_angular_shift_sign(
    label: CircularLabel,
    canvas_config: CircularCanvasConfigurator,
    legend_bounds: Aabb,
) -> int:
    """Return preferred angular direction (+1 ccw / -1 cw) to move away from legend."""
    start_x = float(label["start_x"])
    start_y = float(label["start_y"])
    radius = math.hypot(start_x, start_y)
    if radius <= 1e-6:
        return 1
    legend_local_x, legend_local_y = _legend_center_local(
        canvas_config,
        legend_bounds,
    )
    away_x = start_x - legend_local_x
    away_y = start_y - legend_local_y
    # Unit tangents at current point (CCW and CW).
    ccw_tx = -start_y / radius
    ccw_ty = start_x / radius
    cw_tx = -ccw_tx
    cw_ty = -ccw_ty
    ccw_score = ccw_tx * away_x + ccw_ty * away_y
    cw_score = cw_tx * away_x + cw_ty * away_y
    return 1 if ccw_score >= cw_score else -1


def _expand_shift_block_to_preserve_order(
    labels: list[CircularLabel],
    total_length: int,
    side_indices: list[int],
    center_pos: int,
    delta_angle: float,
) -> tuple[int, int]:
    """Expand move block so angular order cannot be inverted by this shift."""
    left = center_pos
    right = center_pos
    while True:
        expanded = False
        moved_left_idx = side_indices[left]
        moved_right_idx = side_indices[right]
        moved_left_angle = _label_unwrapped_angle(labels[moved_left_idx], total_length) + delta_angle
        moved_right_angle = _label_unwrapped_angle(labels[moved_right_idx], total_length) + delta_angle

        if left > 0:
            prev_idx = side_indices[left - 1]
            prev_angle = _label_unwrapped_angle(labels[prev_idx], total_length)
            if moved_left_angle <= prev_angle + MIN_LABEL_ORDER_GAP_RAD:
                left -= 1
                expanded = True

        if right < len(side_indices) - 1:
            next_idx = side_indices[right + 1]
            next_angle = _label_unwrapped_angle(labels[next_idx], total_length)
            if moved_right_angle >= next_angle - MIN_LABEL_ORDER_GAP_RAD:
                right += 1
                expanded = True

        if not expanded:
            break
    return left, right


def _try_shift_labels_away_from_legend(
    labels: list[CircularLabel],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    legend_bounds: Aabb,
) -> bool:
    """Resolve collisions by moving labels first (preferred strategy)."""
    if not labels:
        return True

    max_passes = 4
    for _ in range(max_passes):
        collided_indices = _legend_collision_indices(
            labels,
            total_length,
            canvas_config,
            legend_bounds,
        )
        if not collided_indices:
            return True

        changed = False
        legend_box = _legend_bbox(legend_bounds)
        for idx in collided_indices:
            label = labels[idx]
            if label.get("is_embedded"):
                continue
            start_x = float(label["start_x"])
            start_y = float(label["start_y"])
            radius = math.hypot(start_x, start_y)
            if radius <= 1e-6:
                continue
            moving_side = 1 if start_x >= 0 else -1
            preferred_sign = _preferred_angular_shift_sign(
                label,
                canvas_config,
                legend_bounds,
            )
            side_indices = [
                side_idx
                for side_idx, side_label in enumerate(labels)
                if not side_label.get("is_embedded") and ((float(side_label["start_x"]) >= 0) == (moving_side > 0))
            ]
            side_indices.sort(key=lambda side_idx: float(labels[side_idx]["middle"]))
            if idx not in side_indices:
                continue
            center_pos = side_indices.index(idx)

            placed = False
            shift_px = LABEL_NUDGE_STEP_PX
            while shift_px <= MAX_LABEL_NUDGE_PX and not placed:
                angle_delta = shift_px / radius
                for direction_sign in (preferred_sign, -preferred_sign):
                    signed_delta = direction_sign * angle_delta
                    block_left, block_right = _expand_shift_block_to_preserve_order(
                        labels, total_length, side_indices, center_pos, signed_delta
                    )
                    block_indices = side_indices[block_left : block_right + 1]

                    candidate_positions: dict[int, tuple[float, float]] = {}
                    candidate_labels: dict[int, CircularLabel] = {}
                    candidate_valid = True

                    for block_idx in block_indices:
                        block_label = labels[block_idx]
                        block_radius = math.hypot(float(block_label["start_x"]), float(block_label["start_y"]))
                        if block_radius <= 1e-6:
                            candidate_valid = False
                            break
                        block_angle = math.atan2(float(block_label["start_y"]), float(block_label["start_x"]))
                        candidate_angle = block_angle + signed_delta
                        candidate_x = block_radius * math.cos(candidate_angle)
                        candidate_y = block_radius * math.sin(candidate_angle)

                        block_sign = 1 if float(block_label["start_x"]) >= 0 else -1
                        if block_sign > 0 and candidate_x <= 0:
                            candidate_valid = False
                            break
                        if block_sign < 0 and candidate_x >= 0:
                            candidate_valid = False
                            break

                        candidate = block_label.copy()
                        _move_label_text_anchor(
                            candidate,
                            start_x=float(candidate_x),
                            start_y=float(candidate_y),
                        )
                        candidate_positions[block_idx] = (float(candidate_x), float(candidate_y))
                        candidate_labels[block_idx] = candidate

                    if not candidate_valid:
                        continue

                    for candidate in candidate_labels.values():
                        candidate_box = _label_bbox_on_canvas(
                            candidate, total_length, canvas_config, margin_px=LEGEND_LABEL_MARGIN_PX
                        )
                        if _bbox_overlaps(candidate_box, legend_box):
                            candidate_valid = False
                            break
                    if not candidate_valid:
                        continue

                    for block_idx, candidate in candidate_labels.items():
                        for peer_idx, peer_label in enumerate(labels):
                            if peer_idx == block_idx or peer_label.get("is_embedded"):
                                continue
                            peer = candidate_labels.get(peer_idx, peer_label)
                            min_gap_px = minimum_bbox_gap_px(candidate, peer)
                            if y_overlap(candidate, peer, total_length, min_gap_px) and x_overlap(
                                candidate, peer, minimum_margin=min_gap_px
                            ):
                                candidate_valid = False
                                break
                        if not candidate_valid:
                            break
                    if not candidate_valid:
                        continue

                    def _unwrapped_for_order(check_idx: int) -> float:
                        candidate = candidate_labels.get(check_idx)
                        if candidate is not None:
                            return _label_unwrapped_angle(candidate, total_length)
                        return _label_unwrapped_angle(labels[check_idx], total_length)

                    for order_pos in range(1, len(side_indices)):
                        prev_idx = side_indices[order_pos - 1]
                        curr_idx = side_indices[order_pos]
                        if _unwrapped_for_order(curr_idx) <= _unwrapped_for_order(prev_idx) + MIN_LABEL_ORDER_GAP_RAD:
                            candidate_valid = False
                            break
                    if not candidate_valid:
                        continue

                    for block_idx, (candidate_x, candidate_y) in candidate_positions.items():
                        _move_label_text_anchor(
                            labels[block_idx],
                            start_x=candidate_x,
                            start_y=candidate_y,
                        )
                    changed = True
                    placed = True
                    break
                shift_px += LABEL_NUDGE_STEP_PX

        if not changed:
            break

    return not _labels_collide_with_legend(
        labels,
        total_length,
        canvas_config,
        legend_bounds,
    )


def _resolve_label_legend_collisions(
    labels: list[CircularLabel],
    total_length: int,
    canvas_config: CircularCanvasConfigurator,
    legend_bounds: Aabb,
) -> None:
    """Shift horizontal labels only for the four corner-overlay placements."""
    try:
        placement = LegendPlacement(str(canvas_config.legend_position))
    except ValueError:
        return
    if not placement.is_overlay:
        return
    if (
        legend_bounds.width <= 0
        or legend_bounds.height <= 0
    ):
        return

    external_labels = [label for label in labels if not label.get("is_embedded")]
    if not external_labels:
        return
    if not _labels_collide_with_legend(
        external_labels,
        total_length,
        canvas_config,
        legend_bounds,
    ):
        return

    if any(label.get("placement") == "radial" for label in external_labels):
        return
    _try_shift_labels_away_from_legend(
        external_labels,
        total_length,
        canvas_config,
        legend_bounds,
    )


def _slot_width_ratio_factor(
    slot: CircularTrackSlot | None,
    *,
    base_radius_px: float,
    base_track_ratio: float,
) -> float | None:
    if slot is None or slot.width is None:
        return None
    denominator = float(base_radius_px) * float(base_track_ratio)
    if denominator <= 0:
        return None
    return float(slot.width.resolve(float(base_radius_px))) / denominator




def _definition_reserved_radius_px(
    facts: DefinitionRecordFacts,
    geometry: CircularCanvasGeometry,
    species: str | None,
    strain: str | None,
    *,
    cfg: GbdrawConfig,
    plot_title: str | None,
    definition_profile: CircularDefinitionProfile,
    record_presentation: RecordPresentation | None = None,
) -> float:
    style = CircularDefinitionStyle.of(cfg, definition_profile)
    definition = DefinitionGroup(
        facts,
        geometry,
        style=style,
        names=DefinitionNameLines.of(
            facts.metadata,
            style.profile,
            species=species,
            strain=strain,
            plot_title=plot_title,
            record_presentation=record_presentation,
        ),
    )

    # Track-slot reservation is a record-local compatibility contract. Keep it on
    # the letters each line paints where it is drawn (its ink box in the faces it is
    # drawn in, centred as drawn), independent from the exact painted bounds used by
    # outer composition, so composition metadata cannot reflow tracks. The reserve
    # reaches the ink box's farthest corner from the centre (OV-1822), which no
    # painted point passes; the track check reads the same boxes.
    max_extent = max(
        (
            math.hypot(max(abs(box.min_x), abs(box.max_x)), max(abs(box.min_y), abs(box.max_y)))
            for box in _centre_definition_ink(definition)
        ),
        default=0.0,
    )

    if max_extent <= FEATURE_BAND_EPSILON:
        return 0.0
    padding_px = max(8.0, 0.02 * float(geometry.radius))
    return float(max_extent + padding_px)


def _definition_clearance_px(layout: CircularRadialLayout) -> float:
    """How far the centre definition's band may reach before it enters a drawn track's band (``_track_bands``).

    This is the narrow row's reading of "the definition over a track": the
    definition's band (``_definition_reserved_radius_px``: the farthest corner
    of its lines' ink boxes, plus padding) against each track's band (its
    drawing, its tick or axis labels, and its labels' room). The notice reads
    the same ink boxes (``_centre_definition_ink``) against what the tracks
    paint (``_definition_lines_on_tracks`` over ``_record_track_ink``), so a
    band that clears the tracks' bands leaves no letters on a track (OV-1854).
    """

    return min((float(band.inner_px) for _slot, band in _track_bands(layout)), default=float("inf"))


@dataclass(frozen=True)
class _CentreDefinitionFit:
    """The tracks and the centre definition drawn among them (PD-OI-078).

    ``reserved_radius_px`` is the definition's band (None without a centre definition); radial labels keep out of
    it. ``narrow_radius_px`` is that band when it enters a track's band (``_definition_clearance_px``), the
    ``center_reserved_radius`` that narrows the inside tracks around it (the "narrow the inside tracks" row); 0 when
    the band clears the tracks' bands or the radius is explicit already.
    """

    layout: CircularRadialLayout
    reserved_radius_px: float | None
    narrow_radius_px: float = 0.0


def _resolve_radial_layout_fitting_definition(
    resolve_layout: Callable[[float | None], CircularRadialLayout],
    *,
    facts: DefinitionRecordFacts,
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    species: str | None,
    strain: str | None,
    plot_title: str | None,
    definition_profile: CircularDefinitionProfile,
    layout_definition_profile: CircularDefinitionProfile,
    definition_position: str,
    definition_shown: bool,
    center_reserved_radius: float | None,
    record_presentation: RecordPresentation | None,
) -> _CentreDefinitionFit:
    """Place the tracks around the centre definition, then draw it among them (PD-OI-078, rule L).

    The definition is never wrapped (DF-D25). The radial layout places the
    tracks around its band (``resolve_circular_radial_layout``): it moves the
    inside tracks outward, then narrows their Auto widths down to
    ``canvas.circular.track_narrowing_min_px``, and when that does not clear
    the band it keeps the widths and the definition is drawn over the tracks;
    the label fit check warns when its painted lines meet a track's ink
    (``_definition_lines_on_tracks``, PD-OI-097). An explicit
    ``center_reserved_radius`` keeps the tracks out of the centre (DF-D05).

    The tracks are placed around the band of ``layout_definition_profile``,
    the definition the record draws without a plot title or, for rows that
    pack on a Multi-Record Canvas, alone (OV-1710, OV-1711), so a title or the
    canvas changes only the centre text. The band and the narrow row read the
    drawn definition (``definition_profile``). A hidden definition
    (``definition_shown`` false) has no band.
    """

    if center_reserved_radius is not None:
        reserved = max(0.0, float(center_reserved_radius))
        return _CentreDefinitionFit(resolve_layout(reserved), reserved)
    if not definition_shown or str(definition_position).strip().lower() != "center":
        return _CentreDefinitionFit(resolve_layout(None), None)

    def band_px(profile: CircularDefinitionProfile) -> float:
        return _definition_reserved_radius_px(
            facts,
            geometry,
            species,
            strain,
            cfg=cfg,
            plot_title=plot_title,
            definition_profile=profile,
            record_presentation=record_presentation,
        )

    reserved = band_px(definition_profile)
    layout_band = reserved if layout_definition_profile == definition_profile else band_px(layout_definition_profile)
    if max(reserved, layout_band) <= FEATURE_BAND_EPSILON:
        return _CentreDefinitionFit(resolve_layout(None), None)
    layout = resolve_layout(float(layout_band) if layout_band > FEATURE_BAND_EPSILON else None)
    narrow_radius = 0.0 if reserved <= _definition_clearance_px(layout) + FEATURE_BAND_EPSILON else reserved
    return _CentreDefinitionFit(layout, reserved, narrow_radius)


def _slot_config_with_dinucleotide(config: Any, nt: str) -> Any:
    if str(getattr(config, "dinucleotide", "")).upper() == str(nt).upper():
        return config
    cloned = copy.copy(config)
    cloned.dinucleotide = str(nt).upper()
    return cloned






def _slot_gc_config_with_axis_font_size(
    config: GcContentConfigurator,
    nt: str,
    tick_font_size: float | None,
) -> GcContentConfigurator:
    slot_config = _slot_config_with_dinucleotide(config, nt)
    if tick_font_size is None or getattr(slot_config, "tick_font_size", None) is not None:
        return slot_config
    cloned = copy.copy(slot_config)
    cloned.tick_font_size = float(tick_font_size)
    return cloned


def _edge_painting_slot_ids(
    slots: Sequence[CircularTrackSlot], configurators: DiagramConfigurators, record: CircularRecordData, cfg: GbdrawConfig
) -> Iterator[str]:
    """The slots whose tracks paint their band's edges, each by its drawer's rule: a feature track (its blocks fill
    the height of their lanes, OV-3192), the GC content band in percent mode with its background or border, a
    conservation ring, a feature density heatmap."""

    gc_paints_edges = GcContentDrawer(configurators.gc).paints_percent_band_edges()
    for slot in slots:
        renderer = str(slot.renderer)
        if renderer == "features":
            paints = True
        elif renderer == "dinucleotide_content":
            paints = gc_paints_edges
        elif renderer == "sequence_conservation":
            track = conservation_slot_track(slot.params, record.conservation_tracks)
            paints = track is not None and conservation_ring_paints_band_edges(track, cfg)
        elif renderer == "feature_density":
            layer = slot_feature_density_layer(str(slot.id), slot, record.feature_density)
            paints = layer is not None and feature_density_paints_band_edges(layer)
        else:
            paints = False
        if paints:
            yield str(slot.id)


def _circular_scalar_axes(
    configurators: DiagramConfigurators,
    record: CircularRecordData,
    slots: Sequence[CircularTrackSlot],
    cfg: GbdrawConfig,
    *,
    show_depth_track: bool,
    depth_tracks_by_index: Mapping[int, DepthTrackData] | None,
) -> CircularScalarAxes:
    """The axes the radial layout keeps room for: those of the depth tracks when they are drawn, each in its own
    settings as its slot draws it (``_slot_depth_source``), and the GC content track's in percent mode; and which of
    ``record``'s ``slots`` paint their band's edges (``_edge_painting_slot_ids``)."""

    gc = configurators.gc
    return CircularScalarAxes.of(
        configurators.depth if show_depth_track else None,
        gc if str(gc.mode).strip().lower() == "percent" else None,
        None if depth_tracks_by_index is None else {index: track.config for index, track in depth_tracks_by_index.items()},
        edge_painting_slot_ids=_edge_painting_slot_ids(slots, configurators, record, cfg),
    )


def _gc_content_axis_font_size(
    configurators: DiagramConfigurators,
    scalar_axes: CircularScalarAxes,
    resolved_track_slots: Sequence[CircularResolvedSlot],
) -> float | None:
    """The tick font the Percent GC axis takes from the first drawn depth axis, as the radial layout kept its room
    for it."""

    gc = configurators.gc
    if str(gc.mode).strip().lower() != "percent":
        return None
    return scalar_axes.gc_tick_font_size(
        gc,
        (
            DepthSlotWidth(depth_slot_track_index(str(slot.id), slot.params), float(slot.requested_width_px))
            for slot in resolved_track_slots
            if slot.renderer == "depth" and slot.requested_width_px is not None
        ),
    )


def _slot_dataframe_for_nt(
    *,
    nt: str,
    default_df: DataFrame,
    default_nt: str,
    dinucleotide_dataframes: dict[str, DataFrame] | None,
) -> DataFrame:
    normalized_nt = str(nt).upper()
    if dinucleotide_dataframes and normalized_nt in dinucleotide_dataframes:
        return dinucleotide_dataframes[normalized_nt]
    if normalized_nt == str(default_nt).upper():
        return default_df
    return DataFrame()


def _slot_content_source(
    resolved_slot: CircularResolvedSlot,
    *,
    gc_df: DataFrame,
    gc_config: GcContentConfigurator,
    dataframes: dict[str, DataFrame] | None,
    tick_font_size_override: float | None,
) -> tuple[str, DataFrame | None, GcContentConfigurator]:
    """The dinucleotide, data (``None`` when there are none) and settings a content slot draws with."""

    default_nt = str(getattr(gc_config, "dinucleotide", "GC")).upper()
    nt = slot_dinucleotide(resolved_slot.params, default_nt)
    slot_df = _slot_dataframe_for_nt(nt=nt, default_df=gc_df, default_nt=default_nt, dinucleotide_dataframes=dataframes)
    config = _slot_gc_config_with_axis_font_size(gc_config, nt, tick_font_size_override)
    return nt, slot_df if f"{nt} content" in slot_df.columns else None, config


def _slot_skew_source(
    resolved_slot: CircularResolvedSlot,
    *,
    gc_df: DataFrame,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    dataframes: dict[str, DataFrame] | None,
) -> tuple[str, DataFrame | None, GcSkewConfigurator]:
    """The dinucleotide, data (``None`` when there are none) and settings a skew slot draws with."""

    default_nt = str(getattr(gc_config, "dinucleotide", "GC")).upper()
    nt = slot_dinucleotide(resolved_slot.params, default_nt)
    slot_df = _slot_dataframe_for_nt(nt=nt, default_df=gc_df, default_nt=default_nt, dinucleotide_dataframes=dataframes)
    config = _slot_skew_config(skew_config, resolved_slot, nt)
    return nt, slot_df if f"{nt} skew" in slot_df.columns else None, config


def _slot_norm_factor(resolved_slot: CircularResolvedSlot, canvas_radius: float) -> float:
    """The slot's anchor radius as a share of the canvas radius, as its drawer takes it."""

    # Radial layout always resolves a radius for a drawn slot.
    anchor_radius_px = resolved_slot.anchor_radius_px
    assert anchor_radius_px is not None
    return float(anchor_radius_px) / float(canvas_radius)


def _slot_depth_source(
    resolved_slot: CircularResolvedSlot,
    depth_df: DataFrame | None,
    depth_config: DepthConfigurator | None,
    depth_tracks_by_index: Mapping[int, DepthTrackData] | None,
) -> tuple[DataFrame | None, DepthConfigurator | None, DepthTrackData | None] | None:
    """The data and settings a depth slot draws with (``None`` where they are missing) and its track when there
    are several; ``None`` when the slot names a track there is not."""

    if depth_tracks_by_index is None:
        return depth_df, depth_config, None
    depth_track = depth_tracks_by_index.get(depth_slot_track_index(str(resolved_slot.id), resolved_slot.params))
    if depth_track is None:
        return None
    return depth_track.df, depth_track.config, depth_track


def _slot_paint_data(
    radial_layout: CircularRadialLayout,
    *,
    gc_df: DataFrame,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    dinucleotide_dataframes: dict[str, DataFrame] | None,
    gc_content_tick_font_size_override: float | None,
    depth_df: DataFrame | None,
    depth_config: DepthConfigurator | None,
    depth_tracks_by_index: Mapping[int, DepthTrackData] | None,
    feature_density: Mapping[str, FeatureDensityLayer] | None,
    profile: CircularRenderProfile,
) -> dict[str, SlotPaintData]:
    """The data and settings each drawn GC content, GC skew and depth slot of ``radial_layout`` paints with, by
    slot id, from the sources the drawing reads (``_slot_content_source``, ``_slot_skew_source``,
    ``_slot_depth_source``). A slot with no data is left out."""

    data: dict[str, SlotPaintData] = {}
    for slot in radial_layout.slots:
        renderer = str(slot.renderer)
        if slot.draw_band_px is None:
            continue
        if renderer == "dinucleotide_content":
            _, content_df, content_config = _slot_content_source(
                slot,
                gc_df=gc_df,
                gc_config=gc_config,
                dataframes=dinucleotide_dataframes,
                tick_font_size_override=gc_content_tick_font_size_override,
            )
            if content_df is not None:
                data[str(slot.id)] = (content_df, content_config)
        elif renderer == "dinucleotide_skew":
            _, skew_df, slot_skew_config = _slot_skew_source(
                slot,
                gc_df=gc_df,
                gc_config=gc_config,
                skew_config=skew_config,
                dataframes=dinucleotide_dataframes,
            )
            if skew_df is not None:
                data[str(slot.id)] = (skew_df, slot_skew_config)
        elif renderer == "depth":
            source = _slot_depth_source(slot, depth_df, depth_config, depth_tracks_by_index)
            if source is not None and source[0] is not None and source[1] is not None:
                data[str(slot.id)] = (source[0], source[1])
        elif renderer == "feature_density":
            # A histogram paints as a depth area; a heatmap band's ink comes from its runs (``_record_track_ink``).
            layer = slot_feature_density_layer(str(slot.id), slot, feature_density)
            if layer is not None and layer.track.style == "histogram":
                data[str(slot.id)] = feature_density_histogram_source(layer, profile)
    return data


def _annotation_ring(resolved_slot: CircularResolvedSlot) -> AnnotationRing:
    """Where an annotation slot's track is drawn, for its content and for the ink that content paints."""

    return AnnotationRing(
        inner_radius_px=float(resolved_slot.draw_inner_radius_px),
        outer_radius_px=float(resolved_slot.draw_outer_radius_px),
        side=str(resolved_slot.side),
        strand_bands_px=resolved_slot.strand_bands_px,
    )


def _slot_conservation_track(
    resolved_slot: CircularResolvedSlot, conservation_tracks: Sequence[ConservationTrack] | None
) -> ConservationTrack | None:
    """The conservation track a conservation slot draws."""

    return conservation_slot_track(resolved_slot.params, conservation_tracks)


class _TickOptions(TypedDict):
    label_side: str
    tick_side: str
    tick_length_px: float | None
    track_preset: str | None


def _slot_tick_options(resolved_slot: CircularResolvedSlot, use_slot_tick_options: bool) -> _TickOptions | None:
    """The sides, tick length and preset a tick slot draws with; ``None`` when it draws neither ticks nor
    labels."""

    tick_layout = resolved_slot.payload if isinstance(resolved_slot.payload, CircularTickLayout) else None
    if tick_layout is not None:
        label_side, tick_side = str(tick_layout.label_side), str(tick_layout.tick_side)
    else:
        label_side, tick_side = tick_sides_for_tick_label_layout(
            resolved_slot.params.get("tick_label_layout"),
            side=resolved_slot.side,
        )[:2]
    label_side, tick_side = label_side.strip().lower(), tick_side.strip().lower()
    if label_side == "none" and tick_side == "none":
        return None
    if not (use_slot_tick_options or tick_layout is not None):
        return {"label_side": "legacy", "tick_side": "legacy", "tick_length_px": None, "track_preset": None}
    return {
        "label_side": label_side,
        "tick_side": tick_side,
        "tick_length_px": (
            tick_layout.tick_length_px
            if tick_layout is not None
            else (
                float(resolved_slot.draw_width_px)
                if resolved_slot.explicit_width and resolved_slot.draw_width_px > 0
                else None
            )
        ),
        "track_preset": (
            tick_layout.track_preset
            if tick_layout is not None
            else normalize_circular_track_preset(
                str(resolved_slot.params.get("preset", resolved_slot.params.get("track_preset", "tuckin")))
            )
        ),
    }




class _TickGroupOptions(_TickOptions):
    radius_override: float
    slot_id: str


def _tick_group_options(
    resolved_slot: CircularResolvedSlot, *, use_slot_tick_options: bool
) -> _TickGroupOptions | None:
    """What the tick group of a tick slot is drawn with; ``None`` when it draws neither ticks nor labels."""

    tick_options = _slot_tick_options(resolved_slot, use_slot_tick_options)
    if tick_options is None:
        return None
    anchor_radius_px = resolved_slot.anchor_radius_px
    assert anchor_radius_px is not None
    return {
        **tick_options,
        "radius_override": float(anchor_radius_px),
        "slot_id": str(resolved_slot.id),
    }


def _with_slot_row(error: ValidationError, slot: CircularResolvedSlot) -> ValidationError:
    """A tick slot's scale labels that do not fit (``CANNOT_FIT``) as the layout's refusal of that slot's row."""

    if track_layout_reason(error) != "CANNOT_FIT":
        return error
    return ValidationError(str(error), diagnostic=cannot_fit_diagnostic(slot.slot_index, renderer=str(slot.renderer)))


def _refuse_unfit_scale_labels(
    radial_layout: CircularRadialLayout,
    units: Sequence[_DrawnUnit],
    canvas_config: CircularCanvasConfigurator,
    user_slot_mode: bool,
) -> None:
    """Refuse the figure, before its feature labels are placed, when a tick slot of the layout it draws has scale
    labels longer than their circle (OV-2559): the slot's row does not fit, and a radius too small for its scale
    labels fails at once, so a requested figure size tries it cheaply. It builds the tick groups the figure draws (a
    render stage: within a Worker transaction, drawing them reuses the build); drawing a slot refuses the same way.
    """

    for slot in radial_layout.slots:
        tick_group = (
            _tick_group_options(slot, use_slot_tick_options=user_slot_mode)
            if slot.renderer == "ticks"
            else None
        )
        if tick_group is None:
            continue
        try:
            for unit in units:
                tick_group_content(unit.record.gb_record, canvas_config, **tick_group, record_transform=unit.transform)
        except ValidationError as error:
            raise _with_slot_row(error, slot) from None


def _feature_band_overrides(
    layout: CircularRadialLayout, radius_px: float, track_ratio: float
) -> tuple[float, float] | None:
    """The feature row's anchor radius and its width as a factor of the configured feature width (``radius_px`` times
    ``track_ratio``); None without a feature row."""

    if layout.features is None:
        return None
    width_factor = float(layout.features.width_px) / max(FEATURE_BAND_EPSILON, radius_px * track_ratio)
    return float(layout.features.anchor_radius_px), width_factor


def _sync_legend_table_for_circular_slots(
    out: LegendTable,
    *,
    slots: Sequence[LegendSlotSource] | None,
    gc: NumericLegendSource,
    skew: NumericLegendSource,
    depth: DepthLegendSource | None,
    depth_tracks: Sequence[DepthLegendSource] | None = None,
    conservation_legend: Sequence[ConservationLegendSource] = (),
    conservation_colors: tuple[str, str] = ("", ""),
    conservation_domain: IdentityColorDomain | None = None,
) -> LegendTable:
    """Replace singleton numeric legend entries with slot-aware entries.

    ``depth`` is the lone depth track's swatch when its frame was built;
    ``conservation_colors`` the default gradient ends (min, max).
    """
    if slots is None:
        return out

    default_nt = str(gc.dinucleotide).upper()
    removable = {
        "Depth",
        f"{default_nt} content",
        f"{default_nt} content (+)",
        f"{default_nt} content (-)",
        f"{default_nt} skew",
        f"{default_nt} skew (+)",
        f"{default_nt} skew (-)",
    }
    for key in removable:
        out.discard_default_track(key)

    depth_by_index = index_depth_track_row(depth_tracks or ())
    depth_count = max(depth_by_index, default=-1) + 1
    for slot in slots:
        renderer = slot.renderer
        if renderer == "depth":
            if depth_by_index:
                track = depth_by_index.get(depth_slot_track_index(slot.id, slot.params))
                if track is None:
                    continue
                swatch = track
                label = slot_legend_label(slot, depth_track_legend_name(track, depth_count))
            elif depth is None:
                continue
            else:
                swatch = depth
                label = slot_legend_label(slot, "Depth")
            out.add_source(label, {
                "type": "solid",
                "fill": swatch.fill_color,
                "stroke": swatch.stroke_color,
                "width": swatch.stroke_width,
            })
        elif renderer == "dinucleotide_content":
            nt = slot_dinucleotide(slot.params, default_nt)
            out.add_numeric_source(slot_legend_label(slot, f"{nt} content"), gc)
        elif renderer == "dinucleotide_skew":
            nt = slot_dinucleotide(slot.params, default_nt)
            out.add_numeric_source(slot_legend_label(slot, f"{nt} skew"), skew.for_skew_slot(slot, nt))
    if conservation_legend:
        _add_conservation_ring_rows(
            out, conservation_legend, conservation_colors, conservation_domain or identity_color_domain(0.0, None)
        )
    return out


def _add_conservation_ring_rows(
    out: LegendTable,
    rings: Sequence[ConservationLegendSource],
    default_colors: tuple[str, str],
    domain: IdentityColorDomain,
) -> None:
    """One row per drawn ring, outermost first (PD-OI-110, scenario revision 2).

    When the drawn rings differ in color (a ring without a color counts as the configured max), each row carries
    its ring's color (``fill``) and its ring's own identity gradient, which its ramp draws; otherwise no row has a
    ``fill`` and every row carries the first drawn ring's gradient, which the one band draws.
    """

    def gradient(color: str | None) -> tuple[str, str]:
        return conservation_track_gradient_colors(
            color, default_min_color=default_colors[0], default_max_color=default_colors[1]
        )

    gradients = [gradient(ring.track_color) for ring in rings]
    one_color = len({max_color for _, max_color in gradients}) == 1
    for ring, (min_color, max_color) in zip(rings, gradients):
        ramp_min, ramp_max = gradients[0] if one_color else (min_color, max_color)
        out.add_source(ring.track_label, {
            "type": "ring",
            "source_index": ring.source_index,
            "track_index": ring.track_index,
            "fill": None if one_color else max_color,
            "stroke": "none",
            "width": 0,
            "min_color": ramp_min,
            "max_color": ramp_max,
            "min_value": domain.minimum,
            "max_value": domain.maximum,
        })


def _draw_circular_feature_density_slot(
    canvas: Drawing,
    resolved_slot: CircularResolvedSlot,
    feature_density: Mapping[str, FeatureDensityLayer] | None,
    *,
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    render_context: CircularRecordRenderContext,
    norm_factor_override: float,
    record_key: str,
    record_transform: RecordDisplayTransform | None,
) -> Drawing:
    """Draw one record's feature density track in its slot: a histogram as a depth area, a heatmap as a band."""

    layer = slot_feature_density_layer(str(resolved_slot.id), resolved_slot, feature_density)
    if layer is None:
        return canvas
    renderer = "feature_density"
    slot_group_id = track_slot_svg_id(resolved_slot.id, renderer=renderer)
    slot_key = GroupKey(record_key, GroupKind.DEPTH, str(resolved_slot.id))
    slot_id = str(resolved_slot.id)
    if layer.track.style == "histogram":
        frame, histogram_config = feature_density_histogram_source(layer, render_context.profile)
        return add_depth_group_on_canvas(
            canvas,
            gb_record,
            frame,
            canvas_config,
            histogram_config,
            group_key=slot_key,
            track_width_override=float(resolved_slot.draw_width_px),
            norm_factor_override=norm_factor_override,
            group_id=str(resolved_slot.id),
            slot_group_id=slot_group_id,
            slot_id=slot_id,
            slot_renderer=renderer,
            record_transform=record_transform,
        )
    heatmap = FeatureDensityHeatmapGroup(
        layer,
        circular_display_length(gb_record, record_transform),
        float(resolved_slot.draw_inner_radius_px),
        float(resolved_slot.draw_outer_radius_px),
        feature_density_color(layer, render_context.profile),
        group_id=slot_group_id,
        record_transform=record_transform,
    )
    return add_feature_density_heatmap_on_canvas(
        canvas, heatmap, canvas_config, group_key=slot_key, slot_id=slot_id, slot_renderer=renderer
    )


def _annotation_track_content(
    track: ResolvedAnnotationTrack,
    *,
    record_id: str,
    record_index: int,
    record_length: int,
    ring: AnnotationRing,
    font_family: str,
    params: AnnotationTrackParams,
    labels: Sequence[CircularAnnotationLabel],
    dom_group_id: str,
    semantic_slot_id: str,
) -> GroupContent:
    """``circular_annotation_track_content``, or the content drawn before with equal inputs while the figure's
    assemblies share their measures (``CircularAnnotationMeasures``): each label fit row's assembly draws the slot
    the figure drew."""

    shown = tuple(labels)
    return _kept_measure(
        lambda measures: measures.contents,
        (semantic_slot_id, dom_group_id, record_id, record_index, record_length, font_family),
        (track, ring, params, shown),
        partial(
            circular_annotation_track_content,
            track,
            record_id=record_id,
            record_index=record_index,
            record_length=record_length,
            ring=ring,
            font_family=font_family,
            params=params,
            labels=shown,
            dom_group_id=dom_group_id,
            semantic_slot_id=semantic_slot_id,
        ),
    )


def _draw_circular_tick_slot(
    canvas: Drawing,
    resolved_slot: CircularResolvedSlot,
    *,
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    precalculated_labels: list[CircularLabel] | None,
    use_slot_tick_options: bool,
    record_transform: RecordDisplayTransform | None,
    record_key: str,
) -> Drawing:
    """Draw one record's ticks slot, clear of the feature labels' leader lines; a row error names the slot."""

    tick_group = _tick_group_options(resolved_slot, use_slot_tick_options=use_slot_tick_options)
    if tick_group is None:
        return canvas
    leaders = _feature_label_leaders(precalculated_labels or (), len(gb_record.seq))
    key = GroupKey(record_key, GroupKind.TICKS, str(resolved_slot.id))
    try:
        return add_tick_group_on_canvas(
            canvas, gb_record, canvas_config, group_key=key, **tick_group, record_transform=record_transform,
            feature_leaders=leaders,
        )
    except ValidationError as error:
        raise _with_slot_row(error, resolved_slot) from None


def _draw_circular_annotation_slot(
    canvas: Drawing,
    resolved_slot: CircularResolvedSlot,
    *,
    gb_record: SeqRecord,
    record_length: int,
    canvas_config: CircularCanvasConfigurator,
    annotation_layout: ResolvedAnnotationTrack | None,
    annotation_labels: Sequence[CircularAnnotationLabel],
    annotation_record_index: int,
    record_key: str,
) -> Drawing:
    """Draw one record's annotation slot, then its backed labels in front of the record's leader lines, which are
    drawn before every slot (CIR-RGN-FRONT, OV-3156)."""

    if annotation_layout is None:
        return canvas
    dom_group_id = track_slot_svg_id(resolved_slot.id, renderer="annotations")
    font_family = str(canvas_config.profile.config.objects.text.font_family)
    content = render_stage(
        "annotation_content",
        _annotation_track_content,
        annotation_layout,
        record_id=str(gb_record.id),
        record_index=annotation_record_index,
        record_length=record_length,
        ring=_annotation_ring(resolved_slot),
        font_family=font_family,
        params=annotation_track_params_from_mapping(resolved_slot.params),
        labels=annotation_labels,
        dom_group_id=dom_group_id,
        semantic_slot_id=str(resolved_slot.id),
    )
    group = place_content(canvas, content, GroupKey(record_key, GroupKind.ANNOTATIONS, str(resolved_slot.id)))
    canvas.add(center_group_on_canvas(group, canvas_config))
    front = circular_region_label_front_content(
        annotation_layout,
        record_index=annotation_record_index,
        font_family=font_family,
        labels=annotation_labels,
        dom_group_id=dom_group_id,
    )
    if front is not None:
        front_group = place_content(canvas, front, GroupKey(record_key, GroupKind.LABEL_TEXT, str(resolved_slot.id)))
        canvas.add(center_group_on_canvas(front_group, canvas_config))
    return canvas


def _draw_resolved_circular_slot(
    canvas: Drawing,
    resolved_slot: CircularResolvedSlot,
    *,
    gb_record: SeqRecord,
    canvas_config: CircularCanvasConfigurator,
    feature_config: FeatureDrawingConfigurator,
    gc_df: DataFrame,
    gc_config: GcContentConfigurator,
    skew_config: GcSkewConfigurator,
    depth_df: DataFrame | None,
    depth_config: DepthConfigurator | None,
    depth_tracks_by_index: Mapping[int, DepthTrackData] | None,
    feature_density: Mapping[str, FeatureDensityLayer] | None,
    conservation_tracks: Sequence[ConservationTrack] | None,
    conservation_min_identity: float,
    dinucleotide_dataframes: dict[str, DataFrame] | None,
    gc_content_tick_font_size_override: float | None = None,
    feature_layers: FeatureBuildResult,
    render_context: CircularRecordRenderContext,
    precalculated_labels: list[CircularLabel] | None = None,
    use_slot_tick_options: bool = True,
    use_feature_anchor_override: bool = True,
    several_feature_slots: bool = False,
    annotation_track_layouts: dict[str, ResolvedAnnotationTrack] | None = None,
    annotation_labels: Mapping[str, Sequence[CircularAnnotationLabel]],
    annotation_record_index: int = 0,
    record_transform: RecordDisplayTransform | None = None,
    record_key: str,
) -> Drawing:
    """Draw one resolved circular slot of one record of the circle (the record alone, or one of several)."""
    cfg = render_context.profile.config
    renderer = str(resolved_slot.renderer)
    # Radial layout always resolves a radius; float(None) already failed here before.
    anchor_radius_px = resolved_slot.anchor_radius_px
    assert anchor_radius_px is not None
    norm_factor_override = _slot_norm_factor(resolved_slot, float(canvas_config.radius))
    if renderer == "spacer":
        return canvas

    if renderer == "annotations":
        return _draw_circular_annotation_slot(
            canvas,
            resolved_slot,
            gb_record=gb_record,
            record_length=circular_display_length(gb_record, record_transform),
            canvas_config=canvas_config,
            annotation_layout=(annotation_track_layouts or {}).get(resolved_slot.id),
            annotation_labels=annotation_labels.get(resolved_slot.id, ()),
            annotation_record_index=annotation_record_index,
            record_key=record_key,
        )

    if renderer == "features":
        base_width = circular_track_width_px(
            "features",
            radius=canvas_config.radius,
            track_ratio=canvas_config.track_ratio,
            factors=cfg.canvas.circular.track_ratio_factors[str(canvas_config.length_param)],
        )
        resolved_feature_width = float(resolved_slot.draw_width_px)
        if render_context.feature_layout is not None:
            resolved_feature_width = float(render_context.feature_layout.width_px)
        ratio_override = None
        if base_width > 0 and resolved_feature_width > 0:
            ratio_override = resolved_feature_width / max(
                FEATURE_BAND_EPSILON,
                float(canvas_config.radius) * float(canvas_config.track_ratio),
            )
        anchor_moved = use_feature_anchor_override or not math.isclose(
            float(anchor_radius_px),
            float(canvas_config.radius),
            rel_tol=1e-9,
            abs_tol=1e-9,
        )
        return add_record_group_on_canvas(
            canvas,
            gb_record,
            canvas_config,
            feature_config,
            group_key=GroupKey(record_key, GroupKind.FEATURES, str(resolved_slot.id)),
            feature_layers=feature_layers,
            render_context=render_context,
            # The record group draws the embedded labels; the legend shift moves only external ones.
            precalculated_labels=tuple(
                label for label in precalculated_labels or () if label.get("is_embedded")
            ),
            feature_track_ratio_factor_override=ratio_override,
            record_index=annotation_record_index,
            group_id=track_slot_svg_id(resolved_slot.id, renderer=renderer),
            slot_id=str(resolved_slot.id),
            # Two feature slots may draw one feature: its id takes the slot's.
            feature_dom_namespace=str(resolved_slot.id) if several_feature_slots else None,
            feature_anchor_radius_px=float(anchor_radius_px) if anchor_moved else None,
        )

    if renderer == "ticks":
        return _draw_circular_tick_slot(
            canvas,
            resolved_slot,
            gb_record=gb_record,
            canvas_config=canvas_config,
            precalculated_labels=precalculated_labels,
            use_slot_tick_options=use_slot_tick_options,
            record_transform=record_transform,
            record_key=record_key,
        )

    if renderer == "depth":
        depth_source = _slot_depth_source(resolved_slot, depth_df, depth_config, depth_tracks_by_index)
        if depth_source is None:
            return canvas
        selected_depth_df, selected_depth_config, _depth_track = depth_source
        if selected_depth_config is None or selected_depth_df is None:
            logger.warning("Skipping circular depth slot '%s' because depth data are unavailable.", resolved_slot.id)
            return canvas
        return add_depth_group_on_canvas(
            canvas,
            gb_record,
            selected_depth_df,
            canvas_config,
            selected_depth_config,
            group_key=GroupKey(record_key, GroupKind.DEPTH, str(resolved_slot.id)),
            track_width_override=float(resolved_slot.draw_width_px),
            norm_factor_override=norm_factor_override,
            group_id=str(resolved_slot.id),
            slot_group_id=track_slot_svg_id(resolved_slot.id, renderer=renderer),
            slot_id=str(resolved_slot.id),
            slot_renderer=renderer,
            record_transform=record_transform,
        )

    if renderer == "feature_density":
        return _draw_circular_feature_density_slot(
            canvas, resolved_slot, feature_density, gb_record=gb_record, canvas_config=canvas_config,
            render_context=render_context, norm_factor_override=norm_factor_override,
            record_key=record_key, record_transform=record_transform,
        )

    if renderer == "sequence_conservation":
        track = _slot_conservation_track(resolved_slot, conservation_tracks)
        if track is None:
            logger.warning(
                "Skipping circular conservation slot '%s' because conservation data are unavailable.",
                resolved_slot.id,
            )
            return canvas
        return add_conservation_group_on_canvas(
            canvas,
            gb_record,
            track,
            canvas_config,
            inner_radius_px=float(resolved_slot.draw_inner_radius_px),
            outer_radius_px=float(resolved_slot.draw_outer_radius_px),
            min_identity=float(conservation_min_identity),
            record_index=annotation_record_index,
            record_transform=record_transform,
            group_key=GroupKey(record_key, GroupKind.CONSERVATION, str(resolved_slot.id)),
            group_id=track_slot_svg_id(resolved_slot.id, renderer=renderer),
            slot_id=str(resolved_slot.id),
        )

    if renderer == "dinucleotide_content":
        nt, slot_df, slot_gc_config = _slot_content_source(
            resolved_slot,
            gc_df=gc_df,
            gc_config=gc_config,
            dataframes=dinucleotide_dataframes,
            tick_font_size_override=gc_content_tick_font_size_override,
        )
        if slot_df is None:
            logger.warning(
                "Skipping circular content slot '%s' because %s data are unavailable.",
                resolved_slot.id,
                nt,
            )
            return canvas
        return add_gc_content_group_on_canvas(
            canvas,
            gb_record,
            slot_df,
            canvas_config,
            slot_gc_config,
            group_key=GroupKey(record_key, GroupKind.GC_CONTENT, str(resolved_slot.id)),
            track_width_override=float(resolved_slot.draw_width_px),
            norm_factor_override=norm_factor_override,
            slot_id=str(resolved_slot.id),
            slot_renderer=renderer,
            group_id=str(resolved_slot.id),
            slot_group_id=track_slot_svg_id(resolved_slot.id, renderer=renderer),
            record_transform=record_transform,
        )

    if renderer == "dinucleotide_skew":
        nt, slot_df, slot_skew_config = _slot_skew_source(
            resolved_slot,
            gc_df=gc_df,
            gc_config=gc_config,
            skew_config=skew_config,
            dataframes=dinucleotide_dataframes,
        )
        if slot_df is None:
            logger.warning(
                "Skipping circular skew slot '%s' because %s data are unavailable.",
                resolved_slot.id,
                nt,
            )
            return canvas
        return add_gc_skew_group_on_canvas(
            canvas,
            gb_record,
            slot_df,
            canvas_config,
            slot_skew_config,
            group_key=GroupKey(record_key, GroupKind.GC_SKEW, str(resolved_slot.id)),
            track_width_override=float(resolved_slot.draw_width_px),
            norm_factor_override=norm_factor_override,
            slot_id=str(resolved_slot.id),
            slot_renderer=renderer,
            group_id=str(resolved_slot.id),
            slot_group_id=track_slot_svg_id(resolved_slot.id, renderer=renderer),
            record_transform=record_transform,
        )

    logger.warning("Skipping unsupported circular track slot renderer '%s'.", renderer)
    return canvas


def _default_outer_label_arena(
    *,
    geometry: CircularCanvasGeometry,
    cfg: GbdrawConfig,
    strandedness: bool,
) -> tuple[float, float]:
    """Return default (anchor_radius_px, arc_outer_radius_px) for outer labels.

    ``cfg`` gives the labels settings (a fit probe passes its own).
    """
    length_param = str(geometry.length_param)
    track_type = _circular_preset_for_layout(cfg)
    strands = "separate" if strandedness else "single"
    base_radius = float(geometry.radius)

    anchor_radius = float(cfg.labels.radius_factor[track_type][strands][length_param]) * base_radius
    offset_cfg = cfg.labels.unified_adjustment.outer_labels
    arc_x = (
        base_radius
        * float(cfg.labels.arc_x_radius_factor[track_type][strands][length_param])
        * float(offset_cfg.x_radius_offset)
    )
    arc_y = (
        base_radius
        * float(cfg.labels.arc_y_radius_factor[track_type][strands][length_param])
        * float(offset_cfg.y_radius_offset)
    )
    arc_outer = max(abs(arc_x), abs(arc_y))

    if arc_outer < anchor_radius:
        arc_outer = anchor_radius
    return anchor_radius, arc_outer




class _DefinitionTrackOverlap(NamedTuple):
    """The centre definition lines that reach into a track (``lines``, how many), and what the fit report offers to clear
    them: a probe of the tracks and text at another font size, and the radius that narrows the inside tracks."""

    lines: int = 0
    probe: CircularDefinitionFitProbe | None = None
    narrow_radius_px: float = 0.0


def _record_label_layout(
    gb_record: SeqRecord,
    profile: CircularRenderProfile,
    drawn_labels: Sequence[CircularLabel],
    place_labels: Callable[..., list[CircularLabel]] | None,
    *,
    record_index: int,
    definition_text: Sequence[Aabb],
    radial_refused: bool,
    capacity: Sequence[CircularLabelCapacity] = (),
    region_labels: Sequence[CircularAnnotationLabel] = (),
    region_warnings: RegionLabelWarnings = RegionLabelWarnings(FITS),
    definition_overlap: _DefinitionTrackOverlap = _DefinitionTrackOverlap(),
    total_length: int | None = None,
) -> CircularLabelLayout | None:
    """The drawn labels' fit by kind, and a probe that places horizontal outer labels with other settings.

    Region labels cut to fit their circle (LP-D87), and those the layout
    counts in ``region_warnings`` (lying on another track, LP-D86; left out
    for want of lane room, LBL-ROOM-STOP-AT-CENTRE-OMIT and LBL-ROOM-CLIP-FIT-OMIT,
    or because the cut would keep no letter, RGN-CUT-EMPTY-OMIT, with
    their font sizes), are a
    ``circular_region`` kind that only warns, and so are the
    lines of the centre definition text that reach into a track (D-12,
    ``circular_definition``). ``capacity`` names the sides whose labels need
    much more room than the map has and were not placed (CAPACITY-STOP). ``total_length`` is the length the
    labels are placed round: the circle's when several records share it.
    """

    definition_track_overlaps = definition_overlap.lines
    total_length = len(gb_record.seq) if total_length is None else total_length
    label_counts: dict[LabelKind, int] = {}
    for label in drawn_labels:
        kind: LabelKind = (
            "circular_embedded"
            if label.get("is_embedded")
            else "circular_inner"
            if label.get("is_inner")
            else "circular_outer"
        )
        label_counts[kind] = label_counts.get(kind, 0) + 1
    placed = bool(label_counts) and place_labels is not None
    cuts = sum(label.cut for label in region_labels)
    region_counts = replace(region_warnings.fit, region_cuts=cuts)
    if not placed and region_counts.fits and not definition_track_overlaps and not capacity:
        return None
    cfg = profile.config
    horizontal = profile.label_placement == "horizontal"

    probe: CircularLabelFitProbe | None = None
    if place_labels is not None and horizontal and label_counts.get("circular_outer"):
        placer = place_labels

        def probe_outer(overrides: Mapping[str, LabelFitValue]) -> LabelFit:
            probe_profile = replace(profile, config=with_config_overrides(cfg, overrides))
            try:
                labels = placer(probe_profile, probe=True)
            except CircularLabelProbeStopped as stopped:
                # The dense candidate already overlaps; crossings are not counted.
                return LabelFit(overlap_pairs=stopped.overlap_pairs, crossing_pairs=0)
            return circular_label_fits(labels, total_length).get("circular_outer", FITS)

        probe = probe_outer

    counts: dict[LabelKind, int] = dict(label_counts) if placed else {}
    counts.update((side.label_kind, side.record.label_count) for side in capacity)
    kind_fits = circular_label_fits(drawn_labels, total_length, definition_text=definition_text) if placed else {}
    if not region_counts.fits:
        counts["circular_region"] = len(region_labels) + region_counts.region_left_out
        kind_fits["circular_region"] = region_counts
    if definition_track_overlaps:
        counts["circular_definition"] = definition_track_overlaps
        kind_fits["circular_definition"] = LabelFit(0, 0, definition_track_overlaps=definition_track_overlaps)

    outer_offsets = cfg.labels.unified_adjustment.outer_labels
    font_size = circular_label_font_size(cfg, total_length)
    return CircularLabelLayout(
        record_id=str(gb_record.id),
        record_index=record_index,
        placement="horizontal" if horizontal else "radial",
        font_size=font_size,
        label_counts=counts,
        kind_fits=kind_fits,
        x_radius_offset=float(outer_offsets.x_radius_offset),
        y_radius_offset=float(outer_offsets.y_radius_offset),
        label_spacing=float(cfg.labels.spacing.circular),
        refused=radial_refused,
        capacity={side.label_kind: side.record for side in capacity},
        features_separated=bool(profile.resolve_overlaps),
        probe=probe,
        # The region labels' own sizes, drawn or left out; a label without a size of its own takes the record's.
        region_font_size=max(
            (*(label.font_size for label in region_labels), *region_warnings.left_out_font_sizes), default=font_size
        ),
        definition_font_size=float(cfg.objects.definition.circular.font_size) if definition_track_overlaps else 0.0,
        definition_probe=definition_overlap.probe if definition_track_overlaps else None,
        definition_narrow_radius_px=definition_overlap.narrow_radius_px if definition_track_overlaps else 0.0,
    )


def _centred_on_record(definition: DefinitionGroup, boxes: Sequence[Aabb]) -> tuple[Aabb, ...]:
    """Boxes of a centred definition's lines in the coordinates of the record's labels.

    ``place_definition_group_on_canvas`` centres the line boxes' union on the
    record's centre, the origin of the label coordinates.
    """

    bounds = definition.local_bounds
    dx = -0.5 * (bounds.min_x + bounds.max_x)
    dy = -0.5 * (bounds.min_y + bounds.max_y)
    return tuple(box.translated(dx, dy) for box in boxes)


def _centre_definition_text(definition: DefinitionGroup) -> tuple[Aabb, ...]:
    """The line boxes of a centred definition, which the labels keep clear of."""

    return _centred_on_record(definition, definition.line_bounds)


def _centre_definition_ink(definition: DefinitionGroup) -> tuple[Aabb, ...]:
    """The ink boxes of a centred definition's lines: what the centre reserve holds and the track check reads."""

    return _centred_on_record(definition, definition.line_inks)


def _radial_layout_feature_dict(feature_layers: FeatureBuildResult) -> dict[str, FeatureObject]:
    """The features the radial layout measures: foreground, then underlays."""

    layout_feature_dict = {
        f"foreground:{feature_id}": feature
        for feature_id, feature in feature_layers.foreground_features.items()
    }
    layout_feature_dict.update(
        {
            f"underlay:{feature_index}": feature
            for feature_index, feature in enumerate(feature_layers.underlay_features)
        }
    )
    return layout_feature_dict


@dataclass(frozen=True)
class _CircularLabelState:
    """What a placement of one record's labels reads besides the label settings.

    The geometry is the one the radial growth left; the fit probes of the
    label-fit check place the labels again with it.
    """

    feature_dict: Mapping[str, FeatureObject]
    total_length: int
    geometry: CircularCanvasGeometry
    strandedness: bool
    track_preset: str
    resolved_feature_anchor_radius_px: float | None
    outer_arena_enabled: bool
    feature_track_ratio_factor_override: float | None
    # What the labels read of the radial layout: the feature lanes and the room around them
    # (``circular_label_radii``: radial labels' content radius and inner reserve, OV-2103; the outer floor, OV-3301).
    feature_layout: CircularFeatureLayout | None
    label_radii: CircularLabelRadii
    feature_lane_direction: CircularFeatureLaneDirection
    record_transform: RecordDisplayTransform | None
    definition_text: tuple[Aabb, ...]
    label_filtering: Mapping[str, Any]
    # The record's ID: a drawn placement reads the labels' room first and names the record when it warns or stops.
    record_id: str


def _place_circular_labels(
    state: _CircularLabelState,
    label_profile: CircularRenderProfile,
    *,
    probe: bool = False,
    candidate_caches: dict[float, CircularLabelCandidateCache],
) -> list[CircularLabel]:
    """Place one record's labels with ``label_profile`` and the geometry of ``state``."""

    font_size = circular_label_font_size(label_profile.config, state.total_length)
    return prepare_label_list(
        state.feature_dict,
        state.total_length,
        (
            float(state.resolved_feature_anchor_radius_px)
            if state.resolved_feature_anchor_radius_px is not None
            else state.geometry.radius
        ),
        state.geometry.track_ratio,
        label_profile,
        # The arena follows the outer label offsets, which a probe may change.
        outer_arena=(
            _default_outer_label_arena(
                geometry=state.geometry,
                cfg=label_profile.config,
                strandedness=state.strandedness,
            )
            if state.outer_arena_enabled
            else None
        ),
        feature_track_ratio_factor_override=state.feature_track_ratio_factor_override,
        feature_layout=state.feature_layout,
        room=state.label_radii,
        track_preset=state.track_preset,
        feature_lane_direction=state.feature_lane_direction,
        _candidate_cache=candidate_caches.setdefault(font_size, {}),
        record_transform=state.record_transform,
        probe=probe,
        definition_text=state.definition_text,
        label_filtering=state.label_filtering,
        record_id=state.record_id,
    )


@dataclass(frozen=True)
class CircularLabelPlacement:
    """The drawn labels of one record and the geometry their radial growth left.

    Result of the ``circular_labels`` stage. The caller writes ``state.geometry.radius``
    into its canvas configurator and draws copies of ``labels``.
    """

    labels: tuple[CircularLabel, ...]
    state: _CircularLabelState
    # The layout the radial growth replanned, or ``None`` when the labels kept the one they were given.
    radial_layout: CircularRadialLayout | None
    feature_layers: FeatureBuildResult
    layout_slots: tuple[CircularTrackSlot, ...]
    feature_slot: CircularTrackSlot | None
    # The annotation slots as packed where the layout draws them.
    annotation_tracks: _AnnotationTracks
    radial_refused: bool
    # The sides whose labels need much more room than the map has: no label of the record was placed (CAPACITY-STOP).
    capacity: tuple[CircularLabelCapacity, ...]
    # Text measurements per label font size, for the fit probes.
    label_candidates: tuple[tuple[float, tuple[CircularLabelCandidate, ...]], ...]


def circular_labels(
    feature_layers: FeatureBuildResult,
    profile: CircularRenderProfile,
    *,
    total_length: int,
    geometry: CircularCanvasGeometry,
    resolved_slots: tuple[CircularResolvedSlot, ...],
    label_radii: CircularLabelRadii,
    layout_slots: tuple[CircularTrackSlot, ...],
    feature_slot: CircularTrackSlot | None,
    feature_lane_direction: CircularFeatureLaneDirection,
    circular_preset: str,
    track_preset: str,
    resolved_feature_anchor_radius_px: float | None,
    feature_track_ratio_factor_override: float | None,
    feature_width_override_requested: bool,
    record_features: ResolvedRecordFeatureInputs | None,
    definition_reserved_radius_px: float | None,
    preferred_anchor_slot_ids: Collection[str],
    scalar_axes: CircularScalarAxes,
    center_reserved_radius_explicit: bool,
    record_transform: RecordDisplayTransform | None,
    definition_text: tuple[Aabb, ...],
    label_filtering: Mapping[str, Any],
    annotation_tracks: _AnnotationTracks,
    record_id: str,
) -> CircularLabelPlacement:
    """Place one record's external and embedded labels; radial labels grow the radius (stage ``circular_labels``).

    The labels read the radial layout through its resolved slots and
    ``label_radii`` (``circular_label_radii``), not the centre band (OV-2103). The growth loop replans the radial layout
    with a larger radius until the radial labels fit or growing stops making
    room. It works on copies: the feature layers it is given and the canvas
    configurator are not changed.
    """

    cfg = profile.config
    candidate_caches: dict[float, CircularLabelCandidateCache] = {}
    replanned_layout: CircularRadialLayout | None = None

    def state() -> _CircularLabelState:
        return _CircularLabelState(
            feature_dict=feature_layers.foreground_features,
            total_length=total_length,
            geometry=geometry,
            strandedness=profile.strandedness,
            track_preset=track_preset,
            resolved_feature_anchor_radius_px=resolved_feature_anchor_radius_px,
            outer_arena_enabled=feature_width_override_requested,
            feature_track_ratio_factor_override=feature_track_ratio_factor_override,
            feature_layout=slots_feature_layout(resolved_slots),
            label_radii=label_radii,
            feature_lane_direction=feature_lane_direction,
            record_transform=record_transform,
            definition_text=definition_text,
            label_filtering=label_filtering,
            record_id=record_id,
        )

    slots = list(layout_slots)
    capacity: tuple[CircularLabelCapacity, ...] = ()
    try:
        labels = _place_circular_labels(state(), profile, candidate_caches=candidate_caches)
    except CircularLabelCapacityExceeded as exceeded:
        # CAPACITY-STOP: the figure stops at the label fit check; nothing of this record is placed.
        labels, capacity = [], exceeded.capacities
    # Radial labels that collide once the canvas stops growing are refused (OV-694, RADIAL-JOIN).
    radial_refused = False
    if profile.label_placement == "radial":
        previous_growth: float | None = None
        stagnant_passes = 0
        preflight_tracks_frozen = False
        while True:
            required_growth = max(
                (
                    float(label.get("required_radius_growth_px", 0.0))
                    for label in labels
                    if not label.get("is_embedded")
                ),
                default=0.0,
            )
            if required_growth <= RADIAL_GROWTH_TOLERANCE_PX:
                break
            if (
                previous_growth is not None
                and required_growth >= previous_growth - RADIAL_GROWTH_TOLERANCE_PX
            ):
                stagnant_passes += 1
            else:
                stagnant_passes = 0
            if stagnant_passes >= 2:
                # Growing the canvas no longer makes room (LP-D49 RADIAL-JOIN).
                radial_refused = True
                break
            previous_growth = required_growth

            if not preflight_tracks_frozen:
                frozen_slots: list[CircularTrackSlot] = []
                for slot_index, slot in enumerate(slots):
                    resolved_slot = next(
                        (
                            candidate
                            for candidate in resolved_slots
                            if int(candidate.slot_index) == int(slot_index)
                        ),
                        None,
                    )
                    # Outside rows are not frozen: they stay beyond the enlarged
                    # axis and the feature slot moved outside, in stack order.
                    if (
                        str(slot.renderer) != "features"
                        and resolved_slot is not None
                        and resolved_slot.anchor_radius_px is not None
                        and resolved_slot.side != "outside"
                    ):
                        frozen_slots.append(
                            replace(
                                slot,
                                radius=ScalarSpec(
                                    float(resolved_slot.anchor_radius_px),
                                    unit="px",
                                ),
                            )
                        )
                    elif str(slot.renderer) == "features" and slot.radius is None and slot.side == "inside":
                        # Radial inner labels need a real arena between features and
                        # the frozen inner tracks. Move the automatic tuck-in feature
                        # slot to the outside of the enlarged axis while preserving
                        # the remaining slot parameters.
                        feature_params = dict(slot.params)
                        feature_params["lane_direction"] = "outside"
                        frozen_slots.append(
                            replace(slot, side="outside", params=feature_params)
                        )
                    else:
                        frozen_slots.append(slot)
                slots = frozen_slots
                feature_slot = next(
                    (
                        slot
                        for slot in slots
                        if slot.enabled and str(slot.renderer) == "features"
                    ),
                    None,
                )
                if feature_slot is not None:
                    feature_lane_direction = _lane_direction_for_feature_slot(
                        feature_slot,
                        track_preset=circular_preset,
                    )
                # Reflow changed the physical slot. Rebuild occupancy through
                # the same planner before measuring any final geometry.
                from ...features.placement import plan_feature_placements

                if feature_slot is not None:
                    # The layers are a stage result: replan on copies of the features.
                    replanned = {
                        feature_id: copy.copy(feature)
                        for feature_id, feature in feature_layers.foreground_features.items()
                    }
                    feature_layers = replace(feature_layers, foreground_features=replanned)
                    plan_feature_placements(
                        replanned,
                        slot=FeaturePlacementSlot(
                            "circular", feature_lane_direction, profile.strandedness,
                        ),
                        record_features=record_features,
                        resolve_overlaps=profile.resolve_overlaps,
                        tolerance_bp=profile.feature_overlap_tolerance_bp,
                        genome_length=total_length,
                    )
                preflight_tracks_frozen = True

            geometry = replace(geometry, radius=float(geometry.radius) + required_growth)
            feature_track_ratio_factor_override = _slot_width_ratio_factor(
                feature_slot,
                base_radius_px=float(geometry.radius),
                base_track_ratio=float(geometry.track_ratio),
            )
            # The annotation slots are packed again where the grown circle draws them (OV-2663).
            where_drawn = _place_annotation_tracks_where_drawn(
                lambda sized, rooms: _CentreDefinitionFit(resolve_circular_radial_layout(
                    total_length=total_length, geometry=geometry, cfg=cfg, slots=sized,
                    feature_dict=_radial_layout_feature_dict(feature_layers),
                    definition_reserved_radius_px=definition_reserved_radius_px,
                    preferred_anchor_slot_ids=preferred_anchor_slot_ids, scalar_axes=scalar_axes,
                    center_reserved_radius_explicit=center_reserved_radius_explicit, annotation_label_rooms=rooms,
                ), definition_reserved_radius_px),
                slots, annotation_tracks,
            )
            replanned_layout, slots, annotation_tracks = where_drawn.fit.layout, where_drawn.slots, where_drawn.tracks
            resolved_slots = replanned_layout.slots
            label_radii = circular_label_radii(replanned_layout)
            band = _feature_band_overrides(replanned_layout, float(geometry.radius), float(geometry.track_ratio))
            if band is not None:
                resolved_feature_anchor_radius_px, feature_track_ratio_factor_override = band
            labels = _place_circular_labels(state(), profile, candidate_caches=candidate_caches)
        radial_refused = radial_refused or any(
            any(label.get("layout_collision_counts", ()))
            for label in labels
            if not label.get("is_embedded")
        )

    return CircularLabelPlacement(
        labels=tuple(labels),
        state=state(),
        radial_layout=replanned_layout,
        feature_layers=feature_layers,
        layout_slots=tuple(slots),
        feature_slot=feature_slot,
        annotation_tracks=annotation_tracks,
        radial_refused=radial_refused,
        capacity=capacity,
        label_candidates=tuple(
            (font_size, cache["candidates"])
            for font_size, cache in candidate_caches.items()
            if "candidates" in cache
        ),
    )


class _DrawnUnit(NamedTuple):
    """One record a Circular figure draws on its circle: the record alone, or one of several on one circle."""

    record: CircularRecordData
    transform: RecordDisplayTransform | None
    record_index: int
    record_key: str


def _drawn_units(
    record: CircularRecordData,
    planned: CircularPlannedTracks,
    presentation: CircularPresentationInputs,
    one_circle: CircularOneCircle | None,
) -> tuple[_DrawnUnit, ...]:
    """The records the circle draws, in drawing order: ``record`` alone, or the records of ``one_circle``."""

    if one_circle is None:
        key = str(record_key_at(record.gb_record, planned.annotation_record_index))
        return (_DrawnUnit(record, presentation.record_transform, planned.annotation_record_index, key),)
    return tuple(
        _DrawnUnit(unit.record, unit.transform, unit.record_index, str(record_key_at(unit.record.gb_record, unit.record_index)))
        for unit in one_circle.units
    )


class _CircleDrawing(NamedTuple):
    """What one Circular canvas draws round its circle: its records in drawing order, the length the drawing reads
    as the record's (the circle's display length when records share it), and what the definition says.

    The definition's facts are built here, once per drawing (OV-3043): the record's, or the facts of the circle's
    records. Every reader takes them: the band the tracks keep clear, the fit probes and the drawn group. A hidden
    definition reads no GC percent."""

    units: tuple[_DrawnUnit, ...]
    display_length: int
    definition_facts: DefinitionRecordFacts
    # Whether the centre definition is drawn: as configured, or, hidden, its count of the records a circle leaves out.
    definition_shown: bool

    @classmethod
    def of(
        cls,
        record: CircularRecordData,
        planned: CircularPlannedTracks,
        presentation: CircularPresentationInputs,
        one_circle: CircularOneCircle | None,
        definition_shown: bool,
    ) -> _CircleDrawing:
        units = _drawn_units(record, planned, presentation, one_circle)
        if one_circle is None:
            facts = definition_record_facts(record.gb_record, gc_percent=definition_shown)
            return cls(units, len(record.gb_record.seq), facts, definition_shown)
        _refuse_one_circle_inputs(planned.annotation_track_layouts)
        facts = one_circle_definition_facts(one_circle, definition_shown=definition_shown)
        return cls(units, one_circle.display_length, facts, definition_shown or facts.summary_only)


def _record_feature_layers(
    position: int,
    unit: _DrawnUnit,
    *,
    feature_config: FeatureDrawingConfigurator,
    profile: CircularRenderProfile,
    label_filtering: Mapping[str, object],
    split_overlaps_by_strand: bool,
    placement_slot: FeaturePlacementSlot | None,
    compute_label_text: bool,
) -> FeatureBuildResult:
    """The feature layers of the record a circle draws at ``position`` (stage ``feature_layers``); the record's
    feature input is the ``position``-th, one per record in drawing order (a single record reads the first)."""

    record_feature_inputs = feature_config.record_features or ()
    return render_stage(
        "feature_layers",
        create_feature_layers,
        unit.record.gb_record,
        feature_config.specific_color_rules,
        feature_config.selected_features_set,
        feature_config.default_color_map,
        profile.strandedness,
        profile.resolve_overlaps,
        label_filtering,
        split_overlaps_by_strand=split_overlaps_by_strand,
        record_features=record_feature_inputs[position] if position < len(record_feature_inputs) else None,
        placement_slot=placement_slot,
        feature_overlap_tolerance_bp=profile.feature_overlap_tolerance_bp,
        feature_shapes=feature_config.feature_shapes,
        feature_visibility_rules=feature_config.feature_visibility_rules,
        compute_label_text=compute_label_text,
        record_transform=unit.transform,
    )


class _CircleFeatureLayers(NamedTuple):
    """The feature layers of the records a circle draws: ``layers`` as the circle measures and labels them (one
    record's, or several records' merged), and each record's own (``unit_keys`` picks them from ``layers``;
    ``None`` when one record has the circle to itself)."""

    layers: FeatureBuildResult
    unit_layers: tuple[FeatureBuildResult, ...]
    unit_keys: tuple[tuple[str, ...], ...] | None

    @classmethod
    def build(
        cls, units: Sequence[_DrawnUnit], build_layers: Callable[[int, _DrawnUnit], FeatureBuildResult], *, shared: bool
    ) -> _CircleFeatureLayers:
        """Each unit's layers from ``build_layers`` (its place in drawing order, the unit), merged when the units
        share one circle (``shared``)."""

        unit_layers = tuple(build_layers(position, unit) for position, unit in enumerate(units))
        if not shared:
            return cls(unit_layers[0], unit_layers, None)
        merged, keys = _merged_feature_layers(unit_layers)
        return cls(merged, unit_layers, keys)


def _unit_track_ink_sources(
    unit: _DrawnUnit,
    layout: CircularRadialLayout,
    canvas_config: CircularCanvasConfigurator,
    configurators: DiagramConfigurators,
    *,
    feature_layers: FeatureBuildResult,
    gc_content_tick_font_size_override: float | None,
    depth_tracks_by_index: Mapping[int, DepthTrackData] | None,
    annotation_track_layouts: Mapping[str, ResolvedAnnotationTrack],
    use_slot_tick_options: bool,
) -> CircularTrackInkSources:
    """What the tracks of one record of the circle paint, as ``layout`` places them (``_record_track_ink``)."""

    record = unit.record
    return CircularTrackInkSources(
        gb_record=record.gb_record,
        geometry=canvas_config.geometry(),
        profile=canvas_config.profile,
        record_transform=unit.transform,
        slot_data=_slot_paint_data(
            layout,
            gc_df=record.gc_df,
            gc_config=configurators.gc,
            skew_config=configurators.skew,
            dinucleotide_dataframes=record.dinucleotide_dataframes,
            gc_content_tick_font_size_override=gc_content_tick_font_size_override,
            depth_df=record.depth_df,
            depth_config=configurators.depth,
            depth_tracks_by_index=depth_tracks_by_index,
            feature_density=record.feature_density,
            profile=canvas_config.profile,
        ),
        feature_density=record.feature_density,
        conservation_tracks=record.conservation_tracks,
        conservation_min_identity=record.conservation_min_identity,
        feature_layers=feature_layers,
        feature_stroke_pad_px=stroke_pad_px(
            configurators.feature.block_stroke_color, configurators.feature.block_stroke_width
        ),
        annotation_track_layouts=annotation_track_layouts,
        annotation_record_index=unit.record_index,
        use_slot_tick_options=use_slot_tick_options,
    )


def _unit_axis_inputs(units: Sequence[_DrawnUnit], *, shared: bool) -> list[tuple[int, RecordAxisInput]]:
    """The track-slot metadata's record index and axis input of each record the circle draws (index 0 for a record
    alone on its canvas, its index among the inputs on a shared circle)."""

    return [
        (
            unit.record_index if shared else 0,
            RecordAxisInput(
                unit.record.gb_record,
                unit.record_key,
                bind_record_transforms((unit.record.gb_record,), None if unit.transform is None else (unit.transform,))[0],
            ),
        )
        for unit in units
    ]


def _record_band_font_size(cfg: GbdrawConfig) -> float:
    """The band names the records in the tick labels' font size."""

    return float(cfg.objects.ticks.tick_labels.font_size)


def _record_band_slot(cfg: GbdrawConfig) -> CircularTrackSlot:
    """The room the record band keeps outside the other slots: a spacer the band is then drawn in."""

    return CircularTrackSlot(
        id=RECORD_BAND_SLOT_ID, renderer="spacer", side="outside",
        width=ScalarSpec(RECORD_BAND_WIDTH_EM * _record_band_font_size(cfg), "px"),
    )


def _add_record_band(
    canvas: Drawing,
    canvas_config: CircularCanvasConfigurator,
    units: Sequence[_DrawnUnit],
    radial_layout: CircularRadialLayout,
) -> None:
    """Draw the record band in the room its spacer kept (none when the band is off): each record's arc and name."""

    band_slot = next((slot for slot in radial_layout.slots if slot.id == RECORD_BAND_SLOT_ID), None)
    band_room = None if band_slot is None else band_slot.reserved_band_px
    if band_room is None:
        return
    band_arcs = []
    for unit in units:
        arc = None if unit.transform is None else circle_unit_arc(unit.transform)
        if arc is None:
            raise ValueError("A record on the record band must lie on a shared circle.")
        band_arcs.append(RecordBandArc(
            unit.record_key, str(unit.record.gb_record.id),
            record_band_names(unit.record.gb_record, unit.record_index + 1), *arc,
        ))
    cfg = canvas_config.profile.config
    canvas.add(center_group_on_canvas(record_band_group(
        tuple(band_arcs),
        inner_px=float(band_room.inner_px),
        outer_px=float(band_room.outer_px),
        font_family=str(cfg.objects.text.font_family),
        font_size=_record_band_font_size(cfg),
    ), canvas_config))


def _refuse_one_circle_inputs(annotation_track_layouts: Mapping[str, object]) -> None:
    """Refuse the region tracks several records on one circle do not draw yet (``_refuse_one_circle_tables`` in
    ``gbdraw/api/diagram.py`` refuses depth and conservation before they are read)."""

    if annotation_track_layouts:
        raise ValidationError(
            "Region tracks are not drawn when several records share one circle.",
            diagnostic={"code": "TRACK_INVALID", "reason": "ONE_CIRCLE_REGIONS"},
        )


def _merged_feature_layers(layers: Sequence[FeatureBuildResult]) -> tuple[FeatureBuildResult, tuple[tuple[str, ...], ...]]:
    """The features of several records on one circle as one set, each keyed by its record's place on the circle
    and its own key, and the keys of each record's features."""

    merged: dict[str, FeatureObject] = {}
    keys: list[tuple[str, ...]] = []
    for position, layer in enumerate(layers):
        unit_keys = tuple(f"{position}:{key}" for key in layer.foreground_features)
        merged.update(zip(unit_keys, layer.foreground_features.values()))
        keys.append(unit_keys)
    rules = frozenset(rule for layer in layers for rule in layer.used_color_rules)
    underlays = tuple(feature for layer in layers for feature in layer.underlay_features)
    return FeatureBuildResult(merged, underlays, rules), tuple(keys)


def _unit_feature_layers(
    layers: FeatureBuildResult, unit_keys: tuple[tuple[str, ...], ...] | None, position: int
) -> FeatureBuildResult:
    """The features one record of the circle draws: all of them when the record has the circle to itself."""

    if unit_keys is None:
        return layers
    features = layers.foreground_features
    return FeatureBuildResult(
        {key: features[key] for key in unit_keys[position] if key in features}, (), layers.used_color_rules
    )


def _walk_ids(element: object) -> Iterator[dict[str, object]]:
    attribs = getattr(element, "attribs", None)
    if isinstance(attribs, dict):
        yield attribs
    for child in getattr(element, "elements", []) or []:
        yield from _walk_ids(child)


class _UnitIds:
    """Keeps the ids of the records a shared circle draws apart: what a later record draws under an id an earlier
    one uses takes the record's suffix (``_rename_unit_ids``). A record alone keeps its ids."""

    def __init__(self, shared: bool) -> None:
        self.shared = shared
        self.used: set[str] = set()

    @contextmanager
    def drawing(self, canvas: Drawing, position: int, unit: _DrawnUnit) -> Iterator[None]:
        """Around drawing ``unit`` (at ``position`` in drawing order) on ``canvas``."""

        start = len(getattr(canvas, "elements", []))
        yield
        if self.shared:
            _rename_unit_ids(getattr(canvas, "elements", [])[start:], self.used, unit.record_index if position else 0)


def _rename_unit_ids(elements: Sequence[object], used_ids: set[str], record_index: int) -> None:
    """Give the ids a later record of the circle draws that an earlier one uses the record's suffix (as a
    Multi-Record Canvas names them), and point the record's references at the new names."""

    renamed: dict[str, str] = {}
    for element in elements:
        for attribs in _walk_ids(element):
            value = attribs.get("id")
            if not isinstance(value, str):
                continue
            if value in used_ids:
                candidate = f"{value}__record_{record_index + 1}"
                ordinal = 2
                while candidate in used_ids:
                    candidate = f"{value}__record_{record_index + 1}_{ordinal}"
                    ordinal += 1
                renamed.setdefault(value, candidate)
                attribs["id"] = candidate
                value = candidate
            used_ids.add(value)
    if not renamed:
        return
    for element in elements:
        for attribs in _walk_ids(element):
            for name, value in tuple(attribs.items()):
                if not isinstance(value, str) or "#" not in value:
                    continue
                for old, new in renamed.items():
                    value = value.replace(f"url(#{old})", f"url(#{new})")
                    if value == f"#{old}":
                        value = f"#{new}"
                attribs[name] = value


class _CircleSlotPlan(NamedTuple):
    """The slots a Circular canvas lays out round its circle, the slots its feature anchor prefers, and whether
    they show a depth track."""

    slots: list[CircularTrackSlot]
    preferred_anchor_slot_ids: frozenset[str]
    show_depth_track: bool


def _circle_slot_plan(
    canvas_config: CircularCanvasConfigurator,
    configurators: DiagramConfigurators,
    planned: CircularPlannedTracks,
    *,
    total_length: int,
    circular_preset: str,
    show_features: bool,
    record_band: bool,
) -> _CircleSlotPlan:
    """The slots the circle lays out over ``total_length``: the planned slots in the preset's order, or the
    preset's own when none were requested, with the tracks the slots or the configuration show; and the record
    band's room outside them when several records share the circle and the band is on (``record_band``)."""

    profile = canvas_config.profile
    cfg = profile.config
    requested = planned.slots
    depth_enabled = bool(profile.show_depth and configurators.depth is not None)
    if requested is not None:
        renderers = {str(slot.renderer) for slot in requested if slot.enabled}
        show_depth_track = bool("depth" in renderers and depth_enabled)
        show_gc_track = bool("dinucleotide_content" in renderers)
        show_skew_track = bool("dinucleotide_skew" in renderers)
        show_ticks_track = bool("ticks" in renderers)
    else:
        show_depth_track = bool(depth_enabled)
        show_gc_track = profile.show_gc
        show_skew_track = profile.show_skew
        show_ticks_track = bool(cfg.objects.scale.show)
    context = CircularPresetContext(
        cfg=cfg,
        geometry=canvas_config.geometry(),
        total_length=total_length,
        strandedness=profile.strandedness,
        show_features=show_features,
        show_ticks=show_ticks_track,
        show_depth=show_depth_track,
        show_gc=show_gc_track,
        show_skew=show_skew_track,
        track_ids=tuple(canvas_config.track_ids.items()),
        dinucleotide=str(getattr(configurators.gc, "dinucleotide", "GC")),
    )
    plan = (
        circular_radial_plan_for_preset(circular_preset, context)
        if requested is None
        else circular_track_slots_from_preset_order(
            list(requested), circular_preset, context, axis_index=planned.circular_track_axis_index
        )
    )
    band = [_record_band_slot(cfg)] if record_band else []
    return _CircleSlotPlan([*plan.slots, *band], plan.preferred_anchor_slot_ids, show_depth_track)


def add_record_on_circular_canvas(
    canvas: Drawing,
    canvas_config: CircularCanvasConfigurator,
    configurators: DiagramConfigurators,
    record: CircularRecordData,
    planned: CircularPlannedTracks,
    presentation: CircularPresentationInputs,
    legend_table: LegendTable,
    one_circle: CircularOneCircle | None = None,
) -> _CircularPlotAssembly:
    """Draw one record's tracks, labels, definition and legend on a Circular canvas, uncomposed.

    ``planned`` holds what annotation planning made of the record's track
    inputs: the slots, the annotations, the axis index and the record's index
    on its canvas. The requested slots are not passed.

    With ``one_circle`` several records share the circle (``record`` is the first): each record's features, ticks and tracks are drawn on its own
    arc, labels are placed once round the circle, and the centre definition speaks for the circle. One record is the one-record case of these steps.
    """
    annotation_track_layouts = planned.annotation_track_layouts
    primary_target_start = len(getattr(canvas, "elements", []))
    record_key = str(record_key_at(record.gb_record, planned.annotation_record_index))
    profile = canvas_config.profile
    cfg = profile.config
    drawing = _CircleDrawing.of(record, planned, presentation, one_circle, bool(cfg.objects.definition.circular.show))
    units, display_length = drawing.units, drawing.display_length
    depth_tracks_by_index = (
        index_depth_track_row(record.depth_tracks)
        if record.depth_tracks is not None
        else None
    )

    show_labels_base = profile.labels_enabled
    # A hidden definition is no definition (no band, no text to fit, no group) unless it draws a circle's count of the
    # records left out (``_CircleDrawing.definition_shown``). The plot title is drawn by the caller either way.
    definition_shown = drawing.definition_shown
    # The caller planned and laid out the annotation slots (_prepare_circular_annotation_tracks).
    effective_circular_track_slots = planned.slots
    resolved_annotations = planned.annotations
    user_slot_mode = effective_circular_track_slots is not None
    circular_preset = normalize_circular_track_preset(cfg.canvas.circular.track_type)
    show_features = slots_draw_features(effective_circular_track_slots)
    layout_slots, preferred_anchor_slot_ids, show_depth_track = _circle_slot_plan(
        canvas_config, configurators, planned, total_length=display_length, circular_preset=circular_preset,
        show_features=show_features, record_band=one_circle is not None and one_circle.record_band,
    )
    feature_slot = next(
        (
            slot
            for slot in layout_slots
            if slot.enabled and str(slot.renderer) == "features"
        ),
        None,
    )
    feature_lane_direction = _lane_direction_for_feature_slot(
        feature_slot,
        track_preset=circular_preset,
    )

    show_external_labels = show_labels_base and show_features
    split_overlaps_by_strand = (
        profile.resolve_overlaps
        and (not profile.strandedness)
        and feature_lane_direction == "split"
    )

    feature_track_ratio_factor_override: float | None = None
    if show_features:
        feature_track_ratio_factor_override = _slot_width_ratio_factor(
            feature_slot,
            base_radius_px=float(canvas_config.radius),
            base_track_ratio=float(canvas_config.track_ratio),
        )
    feature_width_override_requested = feature_track_ratio_factor_override is not None

    compute_label_text = show_external_labels
    label_filtering = (
        label_filter_maps(cfg.labels.filtering.as_dict())
        if compute_label_text
        else {}
    )
    record_feature_layers = partial(
        _record_feature_layers,
        feature_config=configurators.feature,
        profile=profile,
        label_filtering=label_filtering,
        split_overlaps_by_strand=split_overlaps_by_strand,
        placement_slot=(FeaturePlacementSlot(
            "circular", feature_lane_direction, profile.strandedness,
        ) if feature_slot is not None else None),
        compute_label_text=compute_label_text,
    )
    circle_layers = _CircleFeatureLayers.build(units, record_feature_layers, shared=one_circle is not None)
    feature_layers, unit_layers, unit_feature_keys = circle_layers
    # The layers the radial layout measures; the definition count reads their blocks (OV-2088).
    layout_feature_layers = feature_layers
    layout_feature_dict = _radial_layout_feature_dict(layout_feature_layers)
    if feature_layers.underlay_features:
        layout_slots, resolved_annotations, annotation_track_layouts = _add_circular_feature_underlays(
            [unit.record.gb_record for unit in units], [layer.underlay_features for layer in unit_layers],
            layout_slots, resolved_annotations, annotation_track_layouts, canvas_config=canvas_config,
            record_indices=[unit.record_index for unit in units], length=display_length,
        )
    resolved_track_slots: list[CircularResolvedSlot] = []
    resolved_feature_anchor_radius_px: float | None = None
    # Rows that pack from the feature track are placed around the full definition the record draws alone, so a
    # record summary on a Multi-Record Canvas changes only the center text (OV-1711). The preset lanes keep the
    # drawn definition there (D-08 covers the plot title).
    # The scalar axes the radial layout keeps room for; the layout sets the GC axis font they are drawn in.
    scalar_axes = _circular_scalar_axes(
        configurators, record, layout_slots, cfg, show_depth_track=show_depth_track, depth_tracks_by_index=depth_tracks_by_index
    )
    reserved_definition_profile: CircularDefinitionProfile = (
        "full"
        if presentation.definition_profile == "record_summary" and not preferred_anchor_slot_ids
        else presentation.definition_profile if presentation.layout_definition_profile is None else presentation.layout_definition_profile
    )
    # The tracks and the definition text drawn among them, for a config: the figure's own, or a smaller definition
    # font that the definition row probes (D-12); the slots and the room their annotation labels take.
    def place_around_definition(
        layout_cfg: GbdrawConfig, slots: Sequence[CircularTrackSlot], rooms: Mapping[str, AnnotationLabelRoom]
    ) -> _CentreDefinitionFit:
        return _resolve_radial_layout_fitting_definition(
            lambda reserved_radius_px: resolve_circular_radial_layout(
                total_length=display_length,
                geometry=canvas_config.geometry(),
                cfg=layout_cfg,
                slots=slots,
                feature_dict=layout_feature_dict,
                definition_reserved_radius_px=reserved_radius_px,
                preferred_anchor_slot_ids=preferred_anchor_slot_ids,
                scalar_axes=scalar_axes,
                center_reserved_radius_explicit=presentation.center_reserved_radius is not None,
                annotation_label_rooms=rooms,
            ),
            facts=drawing.definition_facts,
            geometry=canvas_config.geometry(),
            cfg=layout_cfg,
            species=presentation.species,
            strain=presentation.strain,
            plot_title=presentation.plot_title,
            definition_profile=presentation.definition_profile,
            layout_definition_profile=reserved_definition_profile,
            definition_position=presentation.definition_position,
            definition_shown=definition_shown,
            center_reserved_radius=presentation.center_reserved_radius,
            record_presentation=presentation.record_presentation,
        )

    def centre_definition(text_cfg: GbdrawConfig) -> DefinitionGroup | None:
        if not definition_shown or str(presentation.definition_position) != "center" or str(presentation.definition_profile) == "shared_common":
            return None
        facts = drawing.definition_facts
        style = CircularDefinitionStyle.of(text_cfg, presentation.definition_profile)
        return DefinitionGroup(
            facts,
            canvas_config.geometry(),
            style=style,
            names=DefinitionNameLines.of(
                facts.metadata,
                style.profile,
                species=presentation.species,
                strain=presentation.strain,
                plot_title=presentation.plot_title,
                record_presentation=presentation.record_presentation,
            ),
            record_index=int(planned.annotation_record_index),
            record_count=int(presentation.definition_record_count),
        )

    def definition_lines_on_tracks(
        layout: CircularRadialLayout,
        definition_lines: Sequence[Aabb],
        labels: Mapping[str, Sequence[CircularAnnotationLabel]],
    ) -> int:
        return len(frozenset().union(*(render_stage(
            "definition_track_overlaps",
            definition_lines_over_tracks,
            tuple(definition_lines),
            layout.axis,
            layout.slots,
            _unit_track_ink_sources(
                unit,
                layout,
                canvas_config,
                configurators,
                feature_layers=_unit_feature_layers(layout_feature_layers, unit_feature_keys, position),
                gc_content_tick_font_size_override=gc_content_tick_font_size_override,
                depth_tracks_by_index=depth_tracks_by_index,
                annotation_track_layouts=annotation_track_layouts,
                use_slot_tick_options=user_slot_mode,
            ),
            labels,
        ) for position, unit in enumerate(units))))

    def definition_clears_tracks(font_size: float) -> bool:
        """A hint for the definition row: the tracks and the text at ``font_size`` before the labels are placed."""

        probe_cfg = with_fit_row_overrides(cfg, {DEFINITION_FONT_PATH: font_size})
        try:
            probe = place_around_definition(probe_cfg, layout_slots, annotation_label_rooms)
        except ValidationError:
            return False
        probe_definition = centre_definition(probe_cfg)
        if probe_definition is None:
            return True
        return definition_lines_on_tracks(probe.layout, _centre_definition_ink(probe_definition), {}) == 0

    centre_fit, layout_slots, annotation_tracks, annotation_label_rooms = _place_annotation_tracks_where_drawn(
        partial(place_around_definition, cfg), layout_slots, _AnnotationTracks(_CircularAnnotationPacking.of(
            layout_slots, resolved_annotations, canvas_config, [unit.record_index for unit in units], display_length
        ), annotation_track_layouts, {}),
    )
    annotation_track_layouts = dict(annotation_tracks.layouts)
    radial_layout = centre_fit.layout
    if not (show_external_labels and profile.label_placement == "radial"):  # Radial labels can move the tracks.
        _refuse_unfit_scale_labels(radial_layout, units, canvas_config, user_slot_mode)
    definition_reserved_radius_px = centre_fit.reserved_radius_px
    band = _feature_band_overrides(radial_layout, float(canvas_config.radius), float(canvas_config.track_ratio))
    if band is not None:
        resolved_feature_anchor_radius_px, feature_track_ratio_factor_override = band
    # External labels: separate group (label arena). Embedded labels remain in the record group.
    # Add labels BEFORE features so leader lines appear behind features.
    precalculated_labels: list[CircularLabel] | None = None
    # LP-D49 DEF-WARN and the LP-D57 repair: labels over the centre definition text (not a plot title); D-12: the
    # text over a track.
    centre = centre_definition(cfg)
    definition_text = _centre_definition_text(centre) if centre is not None else ()
    definition_ink = _centre_definition_ink(centre) if centre is not None else ()
    radial_refused = False
    label_capacity: tuple[CircularLabelCapacity, ...] = ()
    label_state: _CircularLabelState | None = None
    # Text measurements per label font size, shared by the drawn placement and fit probes.
    label_candidate_caches: dict[float, CircularLabelCandidateCache] = {}
    if show_external_labels:
        # The stage key leaves out whether the capacity check runs (off only while a fit row is verified). A result
        # reused across that flag is still right: a verified row's settings differ from the drawn ones, so only a
        # later Generate with a found row's settings can reuse it, and those labels fit (no figure with a found row
        # reaches the stop level in the calibration); reusing it can skip a warning, nothing more.
        placement = render_stage(
            "circular_radial_labels" if profile.label_placement == "radial" else "circular_labels",
            circular_labels,
            feature_layers,
            profile,
            total_length=display_length,
            geometry=canvas_config.geometry(),
            # The labels read the centre band only as their inner reserve (OV-2103).
            resolved_slots=radial_layout.slots,
            label_radii=circular_label_radii(radial_layout),
            layout_slots=tuple(layout_slots),
            feature_slot=feature_slot,
            feature_lane_direction=feature_lane_direction,
            circular_preset=circular_preset,
            track_preset=_circular_preset_for_layout(cfg),
            resolved_feature_anchor_radius_px=resolved_feature_anchor_radius_px,
            feature_track_ratio_factor_override=feature_track_ratio_factor_override,
            feature_width_override_requested=feature_width_override_requested,
            record_features=None if one_circle is not None else next(iter(configurators.feature.record_features or ()), None),
            definition_reserved_radius_px=definition_reserved_radius_px,
            preferred_anchor_slot_ids=preferred_anchor_slot_ids,
            # The radial replan reads the scalar axes and ticks only.
            scalar_axes=scalar_axes,
            center_reserved_radius_explicit=presentation.center_reserved_radius is not None,
            record_transform=presentation.record_transform if one_circle is None else None,
            definition_text=definition_text,
            label_filtering=label_filtering,
            annotation_tracks=annotation_tracks,
            record_id=str(record.gb_record.id),
        )
        # The one write of the radius the labels grew (a stage result is never changed).
        canvas_config.radius = placement.state.geometry.radius
        label_state = placement.state
        if placement.radial_layout is not None:
            radial_layout = placement.radial_layout
        feature_layers = placement.feature_layers
        layout_slots = list(placement.layout_slots)
        feature_slot = placement.feature_slot
        annotation_tracks = placement.annotation_tracks
        annotation_track_layouts = dict(annotation_tracks.layouts)
        annotation_label_rooms = annotation_tracks.label_rooms()
        feature_lane_direction = label_state.feature_lane_direction
        resolved_feature_anchor_radius_px = label_state.resolved_feature_anchor_radius_px
        feature_track_ratio_factor_override = label_state.feature_track_ratio_factor_override
        radial_refused = placement.radial_refused
        label_capacity = placement.capacity
        label_candidate_caches.update(
            (font_size, {"candidates": candidates}) for font_size, candidates in placement.label_candidates
        )
        # Copies: the legend collision shift and the leader pass move them.
        precalculated_labels = [label.copy() for label in placement.labels]

    resolved_track_slots = list(radial_layout.slots)
    track_slot_geometry = _serialize_circular_track_slot_geometry(
        drawn=_unit_axis_inputs(units, shared=one_circle is not None),
        radial_layout=radial_layout,
        layout_slots=layout_slots,
        base_radius_px=float(canvas_config.radius),
    )
    attach_track_slot_geometry(canvas, track_slot_geometry)
    gc_content_tick_font_size_override = _gc_content_axis_font_size(configurators, scalar_axes, resolved_track_slots)

    # Annotation labels lie outside the radial bands too: the canvas holds them and a legend or title avoids them.
    # They move off the feature labels and leaders, which make way for one still over their letters (CIR-RGN-YIELD).
    region_labels = _RecordRegionLabels(
        radial_layout,
        annotation_track_layouts,
        planned.annotation_record_index,
        display_length,
        str(cfg.objects.text.font_family),
        canvas_config,
    )

    # A corner legend must also clear the record and its tracks, even when
    # every label has moved outside and the interior has no label obstacles.
    # The record is an obstacle whatever this record's own legend position:
    # a grid places its legend over records drawn without one, and the Web
    # moves a drawn legend to a corner with these obstacles. The label pass
    # predicts the corner with it (region labels too) before the centre definition is drawn.
    def record_obstacle(definition: Group | None = None) -> Aabb:
        return _collect_circular_primary_bounds(
            radial_layout=radial_layout,
            canvas_config=canvas_config,
            feature_config=configurators.feature,
            gc_config=configurators.gc,
            skew_config=configurators.skew,
            depth_config=configurators.depth,
            overlay_obstacles=annotation_label_obstacles,
            definition_target=definition,
        )

    annotation_labels, annotation_label_obstacles, region_labels_not_clear = region_labels.place(precalculated_labels or ())
    if (moved := region_labels.make_way(precalculated_labels or (), region_labels_not_clear, definition_text)) is not None:
        precalculated_labels = moved
        annotation_labels, annotation_label_obstacles, region_labels_not_clear = region_labels.place(moved)
    overlay_obstacles = (
        *_external_label_obstacles_on_canvas(precalculated_labels or (), display_length, canvas_config),
        *annotation_label_obstacles,
    )
    source_content_bounds = _collect_circular_primary_bounds(
        radial_layout=radial_layout,
        canvas_config=canvas_config,
        feature_config=configurators.feature,
        gc_config=configurators.gc,
        skew_config=configurators.skew,
        depth_config=configurators.depth,
        overlay_obstacles=overlay_obstacles,
    )
    legend_measurement = measure_legend(
        configurators.legend,
        legend_table,
        placement=canvas_config.legend_position,
        wrap_width=source_content_bounds.width,
    )
    legend_placement = LegendPlacement(str(canvas_config.legend_position))
    if (
        legend_placement.is_overlay
        and legend_measurement.local_bounds.width > 0
        and legend_measurement.local_bounds.height > 0
    ):
        legend_collision_bounds = _overlay_legend_bounds_for_collision(
            legend_measurement=legend_measurement,
            primary_bounds=source_content_bounds,
            placement=legend_placement,
            record_obstacles=(record_obstacle(),),
        )
        if precalculated_labels:
            _resolve_label_legend_collisions(
                precalculated_labels,
                display_length,
                canvas_config,
                legend_collision_bounds,
            )
            if profile.label_placement == "horizontal":
                assign_leader_start_points(
                    [
                        label
                        for label in precalculated_labels
                        if not label.get("is_embedded")
                    ],
                    display_length,
                )
            annotation_labels, annotation_label_obstacles, region_labels_not_clear = region_labels.place(
                precalculated_labels
            )
            overlay_obstacles = (
                *_external_label_obstacles_on_canvas(precalculated_labels, display_length, canvas_config),
                *annotation_label_obstacles,
            )
            source_content_bounds = _collect_circular_primary_bounds(
                radial_layout=radial_layout,
                canvas_config=canvas_config,
                feature_config=configurators.feature,
                gc_config=configurators.gc,
                skew_config=configurators.skew,
                depth_config=configurators.depth,
                overlay_obstacles=overlay_obstacles,
            )

    # The labels do not change after this: both label phases take one keyed copy.
    placed_labels: Sequence[CircularLabel] = (
        _pre_keyed(label.copy() for label in precalculated_labels or ()) if show_external_labels else ()
    )
    render_context = CircularRecordRenderContext(
        profile=profile,
        track_preset=circular_preset,
        feature_lane_direction=feature_lane_direction,
        radial_layout=radial_layout,
    )

    unit_ids = _UnitIds(one_circle is not None)
    for position, unit in enumerate(units):
        with unit_ids.drawing(canvas, position, unit):
            canvas = add_axis_group_on_canvas(
                canvas,
                canvas_config,
                group_key=GroupKey(unit.record_key, GroupKind.AXIS_CIRCLE),
                unit_arc=None if unit.transform is None else circle_unit_arc(unit.transform),
            )

    if show_external_labels:
        canvas = add_labels_group_on_canvas(
            canvas,
            record.gb_record,
            canvas_config,
            group_key=GroupKey(record_key, GroupKind.LABEL_LEADERS),
            profile=render_context.profile,
            precalculated_labels=placed_labels,
            phase="leaders",
            display_length=display_length,
        )

    drawable_slots = [
        slot for slot in resolved_track_slots
        if str(slot.renderer) != "feature_labels"
    ]
    several_feature_slots = sum(str(slot.renderer) == "features" for slot in drawable_slots) > 1
    for resolved_slot in slots_in_paint_order(drawable_slots):
        for position, unit in enumerate(units):
            with unit_ids.drawing(canvas, position, unit):
                canvas = _draw_resolved_circular_slot(
                    canvas,
                    resolved_slot,
                    gb_record=unit.record.gb_record,
                    canvas_config=canvas_config,
                    feature_config=configurators.feature,
                    gc_df=unit.record.gc_df,
                    gc_config=configurators.gc,
                    skew_config=configurators.skew,
                    depth_df=unit.record.depth_df,
                    depth_config=configurators.depth,
                    depth_tracks_by_index=depth_tracks_by_index,
                    feature_density=unit.record.feature_density,
                    conservation_tracks=unit.record.conservation_tracks,
                    conservation_min_identity=unit.record.conservation_min_identity,
                    dinucleotide_dataframes=unit.record.dinucleotide_dataframes,
                    gc_content_tick_font_size_override=gc_content_tick_font_size_override,
                    feature_layers=_unit_feature_layers(feature_layers, unit_feature_keys, position),
                    render_context=render_context.on_circle(None if unit.transform is None else unit.transform.circle),
                    precalculated_labels=precalculated_labels if position == 0 else None,
                    use_slot_tick_options=user_slot_mode,
                    use_feature_anchor_override=user_slot_mode,
                    several_feature_slots=several_feature_slots,
                    annotation_track_layouts=annotation_track_layouts,
                    annotation_labels=annotation_labels,
                    annotation_record_index=unit.record_index,
                    record_transform=unit.transform,
                    record_key=unit.record_key,
                )

    if one_circle is not None:
        _add_record_band(canvas, canvas_config, units, radial_layout)

    if show_external_labels:
        canvas = add_labels_group_on_canvas(
            canvas,
            record.gb_record,
            canvas_config,
            group_key=GroupKey(record_key, GroupKind.LABEL_TEXT),
            profile=render_context.profile,
            precalculated_labels=placed_labels,
            phase="text",
            display_length=display_length,
        )

    definition_target: Group | None = None
    if definition_shown:
        definition_kwargs: dict[str, Any] = {}
        definition_kwargs["record_index"] = int(planned.annotation_record_index)
        definition_kwargs["record_count"] = int(presentation.definition_record_count)
        if presentation.plot_title is not None:
            definition_kwargs["plot_title"] = presentation.plot_title
        if str(presentation.definition_profile) != "full":
            definition_kwargs["definition_profile"] = presentation.definition_profile
        if str(presentation.definition_position) != "center":
            definition_kwargs["definition_position"] = presentation.definition_position
        if presentation.definition_group_id is not None:
            definition_kwargs["definition_group_id"] = presentation.definition_group_id
        definition_target_start = len(getattr(canvas, "elements", []))
        canvas = add_record_definition_group_on_canvas(
            canvas,
            drawing.definition_facts,
            canvas_config,
            presentation.species,
            presentation.strain,
            group_key=GroupKey(record_key, GroupKind.DEFINITION),
            **definition_kwargs,
            record_presentation=presentation.record_presentation,
        )
        definition_target = next(
            (
                element
                for element in getattr(canvas, "elements", [])[definition_target_start:]
                if isinstance(element, Group)
            ),
            None,
        )
    # The check counts the labels as drawn: after the legend shift and the last leader pass.
    label_layout = _record_label_layout(
        record.gb_record,
        profile,
        precalculated_labels or [],
        (
            partial(_place_circular_labels, label_state, candidate_caches=label_candidate_caches)
            if label_state is not None
            else None
        ),
        # The index the record's definition group carries (``definition_kwargs``).
        record_index=int(planned.annotation_record_index),
        definition_text=definition_text,
        radial_refused=radial_refused,
        capacity=label_capacity,
        region_labels=tuple(label for labels in annotation_labels.values() for label in labels),
        region_warnings=_region_label_warnings(
            radial_layout, annotation_track_layouts, planned.annotation_record_index, annotation_labels, region_labels_not_clear
        ),
        definition_overlap=_DefinitionTrackOverlap(
            definition_lines_on_tracks(radial_layout, definition_ink, annotation_labels) if definition_ink else 0,
            definition_clears_tracks,
            centre_fit.narrow_radius_px,
        ),
        total_length=display_length,
    )
    source_content_bounds = _collect_circular_primary_bounds(
        radial_layout=radial_layout,
        canvas_config=canvas_config,
        feature_config=configurators.feature,
        gc_config=configurators.gc,
        skew_config=configurators.skew,
        depth_config=configurators.depth,
        overlay_obstacles=overlay_obstacles,
        definition_target=definition_target,
    )
    overlay_obstacles = (*overlay_obstacles, record_obstacle(definition_target))
    primary_targets = tuple(
        element
        for element in getattr(canvas, "elements", [])[primary_target_start:]
        if isinstance(element, Group)
    )
    legend = place_legend(canvas, circular_legend_content, legend_placement, legend_measurement, legend_table)
    attach_feature_placement_targets(
        track_slot_geometry["records"],
        FeaturePlacementSlot("circular", feature_lane_direction, profile.strandedness) if feature_slot is not None else None,
    )

    return _CircularPlotAssembly(
        drawing=canvas,
        canvas_config=canvas_config,
        legend_measurement=legend_measurement,
        legend_table=legend_table,
        source_content_bounds=source_content_bounds,
        source_overlay_obstacles=overlay_obstacles,
        primary_targets=primary_targets,
        legend_placement=legend_placement,
        legend=legend,
        track_slot_geometry=track_slot_geometry,
        label_layout=label_layout,
    )


def _compose_circular_plot(
    plot: _CircularPlotAssembly,
    composition: CompositionInputs,
) -> CircularAssemblyResult:
    """Compose one fully measured Circular plot exactly once, in the requested canvas size when one is set, with
    the canvas padding outside its edge padding."""
    title_target = place_title(plot.drawing, composition)
    plan = compose_drawing(
        plot.drawing,
        composition,
        CompositionPrimary(plot.source_content_bounds, plot.primary_targets, plot.source_overlay_obstacles),
        plot.legend_placement,
        plot.legend,
        title_target,
    )

    primary_placement = plan.placement_for("primary")
    if primary_placement is None:  # pragma: no cover - request always has primary
        raise RuntimeError("circular composition has no primary placement")
    final_obstacles = tuple(
        obstacle.translated(primary_placement.dx, primary_placement.dy)
        for obstacle in plot.source_overlay_obstacles
    )
    plot.canvas_config.total_width = plan.width
    plot.canvas_config.total_height = plan.height
    plot.canvas_config.offset_x = float(plot.canvas_config.offset_x) + primary_placement.dx
    plot.canvas_config.offset_y = float(plot.canvas_config.offset_y) + primary_placement.dy
    legend_plan_placement = plan.placement_for("legend")
    if legend_plan_placement is not None:
        plot.canvas_config.legend_offset_x = legend_plan_placement.dx
        plot.canvas_config.legend_offset_y = legend_plan_placement.dy

    result = CircularAssemblyResult(
        drawing=plot.drawing,
        legend_measurement=plot.legend_measurement,
        legend_table=plot.legend_table,
        source_content_bounds=plot.source_content_bounds,
        content_bounds=plan.primary_bounds,
        overlay_obstacles=final_obstacles,
        primary_targets=plot.primary_targets,
        legend=plot.legend,
        title_target=title_target,
        track_slot_geometry=plot.track_slot_geometry,
        composition_plan=plan,
        label_layout=plot.label_layout,
    )
    setattr(plot.drawing, "_gbdraw_circular_assembly_result", result)
    return result


def _circular_legend_slots(
    slots: tuple[LegendSlotSource, ...] | None,
    depth_tracks: tuple[DepthLegendSource, ...] | None,
    *,
    show_gc: bool,
    show_skew: bool,
    gc: NumericLegendSource,
) -> tuple[LegendSlotSource, ...] | None:
    """The slots whose order the legend's track rows follow: the explicit ones, else, with Depth series, the
    default stack the rings are drawn in (``default_circular_track_slots``), so the Depth rows stand between
    the feature rows and the GC rows with or without a stack (OV-2823)."""

    if slots is not None or not depth_tracks:
        return slots
    return legend_slot_sources(
        normalize_circular_track_slots(
            default_circular_track_slots(
                show_ticks=False,
                show_depth=True,
                depth_track_count=max(int(track.track_index) for track in depth_tracks) + 1,
                show_gc=show_gc,
                show_skew=show_skew,
                dinucleotide=str(gc.dinucleotide),
            )
        )
    )


def build_circular_legend_table(
    records: Sequence[SeqRecord],
    *,
    sources: LegendSources,
    show_gc: bool,
    show_skew: bool,
    show_depth: bool,
    depth: DepthLegendSource | None = None,
    has_depth_frame: bool = False,
    depth_tracks: tuple[DepthLegendSource, ...] | None = None,
    depth_track_count_value: int = 0,
    slots: tuple[LegendSlotSource, ...] | None = None,
    annotation_slots: tuple[LegendSlotSource, ...] | None = None,
    annotations: ResolvedAnnotationBundle | None = None,
    conservation_legend: Sequence[ConservationLegendSource] = (),
    conservation_colors: tuple[str, str] = ("", ""),
    conservation_domain: IdentityColorDomain | None = None,
    edits: LegendRowEdits = LegendRowEdits(),
) -> LegendTable:
    """Build the one Circular legend table for the records drawn on a canvas.

    Single-record and Multi-Record Canvas output both call this owner, so slot
    renderer parameters, region annotation and depth ``legend_label`` values,
    and conservation rows reach every Circular legend in the same order.
    ``circular_legend_table`` translates the drawing settings into these
    legend-only arguments. The request's Legend row edits are applied last.
    """

    features_present, used_color_rules, default_used_features = feature_legend_presence(
        records, sources.features, slots
    )
    singleton_depth = bool(
        show_depth
        and depth is not None
        and (has_depth_frame or int(depth_track_count_value) == 1)
    )
    legend_table = prepare_legend_table(
        sources.gc,
        sources.skew,
        sources.features,
        features_present,
        used_color_rules=used_color_rules,
        default_used_features=default_used_features,
        depth=depth if singleton_depth else None,
        show_gc=show_gc,
        show_skew=show_skew,
        show_depth=singleton_depth,
    )
    legend_table = _sync_legend_table_for_circular_slots(
        legend_table,
        slots=_circular_legend_slots(slots, depth_tracks, show_gc=show_gc, show_skew=show_skew, gc=sources.gc),
        gc=sources.gc,
        skew=sources.skew,
        depth=depth if has_depth_frame else None,
        depth_tracks=depth_tracks,
        conservation_legend=conservation_legend,
        conservation_colors=conservation_colors,
        conservation_domain=conservation_domain,
    )
    legend_table = sync_feature_density_legend_entries(legend_table, sources.feature_density, slots=slots)
    if annotations is not None:
        legend_table = sync_annotation_legend_entries(legend_table, annotations, annotation_slots)
    return apply_legend_row_edits(legend_table, edits)


def circular_legend_slot_sources(
    slots: Sequence[CircularTrackSlot] | None, axis_index: int | None
) -> tuple[LegendSlotSource, ...] | None:
    """The slot list after annotation planning, in the Legend's order: normalized with its Axis index, so each slot
    has the side it is drawn on (``legend_slot_sources``)."""

    return legend_slot_sources(normalize_circular_track_slots_with_axis(slots, axis_index)) if slots is not None else None


def circular_legend_table(
    records: Sequence[SeqRecord],
    sources: LegendSources,
    *,
    profile: CircularRenderProfile,
    depth_config: DepthConfigurator | None,
    depth_df: DataFrame | None,
    depth_tracks: Sequence[DepthTrackData] | None,
    depth_track_count_value: int,
    slots: tuple[LegendSlotSource, ...] | None,
    annotations: ResolvedAnnotationBundle | None,
    conservation_legend: Sequence[ConservationLegendSource],
    conservation_min_identity: float,
) -> LegendTable:
    """Call the ``legend_table`` stage with the legend-only parts of the drawing settings.

    ``slots`` are the drawing's slots in the Legend's order (``circular_legend_slot_sources``).
    """

    conservation = profile.config.objects.conservation
    return render_stage(
        "legend_table",
        build_circular_legend_table,
        list(records),
        sources=sources,
        show_gc=profile.show_gc,
        show_skew=profile.show_skew,
        show_depth=profile.show_depth,
        depth=DepthLegendSource.of_config(depth_config) if depth_config is not None else None,
        has_depth_frame=depth_df is not None,
        depth_tracks=depth_legend_sources(depth_tracks) if depth_tracks is not None else None,
        depth_track_count_value=depth_track_count_value,
        slots=slots,
        annotation_slots=slots if annotations is not None else None,
        annotations=annotations,
        conservation_legend=tuple(conservation_legend),
        conservation_colors=(conservation.min_color, conservation.max_color),
        conservation_domain=identity_color_domain(conservation_min_identity, conservation.identity_scale),
        edits=LegendRowEdits.of_config(profile.config.objects.legends),
    )


def _assemble_circular_diagram_result(
    canvas_config: CircularCanvasConfigurator,
    configurators: DiagramConfigurators,
    record: CircularRecordData,
    tracks: CircularTrackInputs,
    presentation: CircularPresentationInputs,
    composition: CompositionInputs,
    one_circle: CircularOneCircle | None = None,
) -> CircularAssemblyResult:
    """Assemble and compose one Circular record: plan its annotation tracks, build its legend, draw it, compose it.

    Returns the composed canvas and its authoritative geometry. The assembly shares the annotation measures its
    figure keeps (``kept_annotation_measures``), or keeps its own while it runs. With ``one_circle`` the circle
    draws its records one after another (``record`` is the first).
    """

    if _KEPT_MEASURES.get() is None:
        with kept_annotation_measures(CircularAnnotationMeasures()):
            return _assemble_circular_record(
                canvas_config, configurators, record, tracks, presentation, composition, one_circle
            )
    return _assemble_circular_record(canvas_config, configurators, record, tracks, presentation, composition, one_circle)


def _assemble_circular_record(
    canvas_config: CircularCanvasConfigurator,
    configurators: DiagramConfigurators,
    record: CircularRecordData,
    tracks: CircularTrackInputs,
    presentation: CircularPresentationInputs,
    composition: CompositionInputs,
    one_circle: CircularOneCircle | None,
) -> CircularAssemblyResult:

    profile = canvas_config.profile
    legend_records = [record.gb_record] if one_circle is None else [unit.record.gb_record for unit in one_circle.units]
    cfg = profile.config
    resolved_depth_track_count = (
        int(record.depth_track_count_value)
        if record.depth_track_count_value is not None
        else depth_track_data_count([record.depth_tracks or ()])
    )
    planned = _prepare_circular_annotation_tracks(
        record.gb_record,
        tracks,
        canvas_config=canvas_config,
        show_ticks=bool(cfg.objects.scale.show),
        show_depth=bool(profile.show_depth and configurators.depth is not None),
        show_gc=profile.show_gc,
        show_skew=profile.show_skew,
        depth_track_count=max(1, resolved_depth_track_count),
        record_transform=presentation.record_transform,
    )

    legend_sources = LegendSources.of(
        configurators.feature, configurators.gc, configurators.skew, feature_density_legend(record.feature_density, profile)
    )
    # The table the Legend draws; with no Legend, the one that names the rings (a Multi-Record Canvas names them
    # from its shared Legend before it draws each record).
    legend_table = LegendTable()
    unnamed = any(track.legend_name is None for track in record.conservation_tracks or ())
    if canvas_config.legend_position != "none" or unnamed:
        legend_table = circular_legend_table(
            legend_records,
            legend_sources,
            profile=profile,
            depth_config=configurators.depth,
            depth_df=record.depth_df,
            depth_tracks=record.depth_tracks,
            depth_track_count_value=resolved_depth_track_count,
            slots=circular_legend_slot_sources(planned.slots, planned.circular_track_axis_index),
            annotations=planned.annotations,
            conservation_legend=conservation_legend_sources(record.conservation_tracks),
            conservation_min_identity=record.conservation_min_identity,
        )
    if unnamed:
        record = replace(
            record, conservation_tracks=named_conservation_tracks(record.conservation_tracks or (), legend_table.ring_names())
        )
    if canvas_config.legend_position == "none":
        legend_table = LegendTable()
    canvas: Drawing = canvas_config.create_svg_canvas()
    attach_legend_row_facts(
        canvas,
        LegendRowFacts(
            records=legend_records,
            sources=legend_sources,
            legend_table=legend_table,
        ),
    )
    plot = add_record_on_circular_canvas(
        canvas,
        canvas_config,
        configurators,
        record,
        planned,
        presentation,
        legend_table,
        one_circle,
    )
    return _compose_circular_plot(plot, composition)


__all__ = [
    "CircularAnnotationMeasures",
    "CircularAssemblyResult",
    "add_record_on_circular_canvas",
    "build_circular_legend_table",
    "circular_legend_slot_sources",
    "kept_annotation_measures",
    "plan_circular_annotation_slots",
]

