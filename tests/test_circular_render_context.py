from dataclasses import FrozenInstanceError
from inspect import Parameter, signature

import pytest

from gbdraw.config.models import CircularRenderProfile, GbdrawConfig
from gbdraw.config.toml import load_config_toml
from gbdraw.layout.circular import (
    CircularAxisLayout,
    CircularFeatureLane,
    CircularFeatureLayout,
    CircularRadialLayout,
    CircularRecordRenderContext,
    CircularResolvedSlot,
    RadialBand,
)
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.canvas import CircularCanvasConfigurator, CircularCanvasGeometry
from gbdraw.render.groups.circular.axis import AxisGroup
from gbdraw.render.groups.circular.definition import DefinitionGroup
from gbdraw.render.groups.circular.labels import LabelsGroup
from gbdraw.render.groups.circular.seq_record import SeqRecordGroup
from gbdraw.render.groups.circular.ticks import TickGroup


def _profile() -> CircularRenderProfile:
    return CircularRenderProfile(
        GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    )


def _feature_radial_layout() -> CircularRadialLayout:
    feature_layout = CircularFeatureLayout(
        anchor_radius_px=95.0,
        width_px=10.0,
        lanes_by_track_id={
            0: CircularFeatureLane(
                track_id=0,
                strand_group="combined",
                inner_px=90.0,
                center_px=95.0,
                outer_px=100.0,
            )
        },
        primary_band_px=RadialBand(90.0, 100.0),
        all_band_px=RadialBand(90.0, 100.0),
    )
    feature_slot = CircularResolvedSlot(
        slot_index=0,
        id="features",
        renderer="features",
        side="inside",
        z=10,
        anchor_radius_px=95.0,
        anchor_offset_px=-5.0,
        requested_width_px=10.0,
        resolved_width_px=10.0,
        packing_band_px=RadialBand(90.0, 100.0),
        draw_band_px=RadialBand(90.0, 100.0),
        reserved_band_px=RadialBand(90.0, 100.0),
        inner_gap_px=0.0,
        outer_gap_px=0.0,
        params={"nested": {"values": [1, 2]}},
        payload=feature_layout,
    )
    return CircularRadialLayout(
        axis=CircularAxisLayout(radius_px=100.0, stroke_width_px=1.0),
        slots=(feature_slot,),
        definition_reserved_band_px=None,
        outer_content_radius_px=100.0,
    )


def test_circular_record_render_context_is_frozen_and_owns_feature_layout() -> None:
    radial_layout = _feature_radial_layout()
    context = CircularRecordRenderContext(
        profile=_profile(),
        track_preset="tuckin",
        feature_lane_direction="inside",
        radial_layout=radial_layout,
    )

    assert context.feature_layout is radial_layout.features
    with pytest.raises(FrozenInstanceError):
        context.track_preset = "spreadout"  # type: ignore[misc]


def test_circular_record_render_context_deep_freezes_layout_collections() -> None:
    context = CircularRecordRenderContext(
        profile=_profile(),
        track_preset="tuckin",
        feature_lane_direction="inside",
        radial_layout=_feature_radial_layout(),
    )
    feature_layout = context.feature_layout

    assert feature_layout is not None
    with pytest.raises(TypeError):
        feature_layout.lanes_by_track_id[1] = feature_layout.lane_for_track_id(0)  # type: ignore[index]
    with pytest.raises(TypeError):
        context.radial_layout.slots[0].params["nested"]["extra"] = True  # type: ignore[index]
    assert context.radial_layout.slots[0].params["nested"]["values"] == (1, 2)


@pytest.mark.parametrize("group_type", [SeqRecordGroup, LabelsGroup])
def test_circular_record_render_groups_require_explicit_context(group_type: type) -> None:
    render_context = signature(group_type).parameters["render_context"]

    assert render_context.default is Parameter.empty


@pytest.mark.parametrize("group_type", [AxisGroup, TickGroup, DefinitionGroup, LabelsGroup, SeqRecordGroup])
def test_circular_group_classes_read_the_canvas_geometry_not_the_configurator(group_type: type) -> None:
    parameters = signature(group_type).parameters

    assert "canvas_config" not in parameters
    assert parameters["geometry"].annotation in {CircularCanvasGeometry, "CircularCanvasGeometry"}


def test_circular_canvas_geometry_is_a_frozen_value_of_the_configurator() -> None:
    record = SeqRecord(Seq("ATGC" * 250), id="NC_TEST.1")
    canvas = CircularCanvasConfigurator("test", _profile(), "none", record)
    geometry = canvas.geometry()

    assert (geometry.radius, geometry.track_ratio, geometry.length_param) == (
        canvas.radius, canvas.track_ratio, canvas.length_param,
    )
    assert geometry.track_ratio_factors == tuple(canvas.track_ratio_factors)
    assert hash(geometry) == hash(canvas.geometry())
    with pytest.raises(FrozenInstanceError):
        geometry.radius = 1.0  # type: ignore[misc]
    canvas.radius += 10.0
    assert canvas.geometry().radius == geometry.radius + 10.0


def test_circular_track_groups_are_placed_contents_with_their_slot_ids_inside() -> None:
    """Axis, features, ticks, GC, skew and annotation groups are contents of render stages (PR 4).

    A slot group's DOM ids are given inside its stage, so the stored content
    already holds them and the placed group only adds its position and slot tags.
    """

    from pathlib import Path

    from Bio import SeqIO

    from gbdraw.api.diagram import assemble_circular_diagram_from_record
    from gbdraw.render.group_content import ContentGroup, GroupKind, group_manifest

    record = SeqIO.read(Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk", "genbank")
    drawing = assemble_circular_diagram_from_record(
        record,
        cfg=GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml")),
        circular_track_slots=["features:features", "ticks:ticks", "gc_content:dinucleotide_content", "gc_skew:dinucleotide_skew"],
    )
    kinds = {entry.key.kind for entry in group_manifest(drawing)}
    assert {GroupKind.AXIS_CIRCLE, GroupKind.FEATURES, GroupKind.TICKS, GroupKind.GC_CONTENT, GroupKind.GC_SKEW} <= kinds
    for element in drawing.elements:
        if isinstance(element, ContentGroup) and element.key.kind in {GroupKind.GC_CONTENT, GroupKind.GC_SKEW}:
            # The slot DOM id (not the slot id) is already the stored content's id.
            assert element.content.attributes["id"] == element.attribs["id"] != element.key.slot
            assert element.attribs["data-gbdraw-slot-id"] == element.key.slot
            assert "data-gbdraw-slot-id" not in element.content.attributes


def test_circular_conservation_ring_is_a_placed_content_of_its_track() -> None:
    """A conservation ring is the content of a render stage over its track (PR 4)."""

    from pathlib import Path

    from Bio import SeqIO

    from gbdraw.api.diagram import assemble_circular_diagram_from_record
    from gbdraw.render.group_content import ContentGroup, GroupKind

    examples = Path(__file__).parent.parent / "examples"
    record = SeqIO.read(examples / "MjeNMV.gb", "genbank")
    drawing = assemble_circular_diagram_from_record(
        record,
        cfg=GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml")),
        conservation_blast_files=[str(examples / "MjeNMV.MelaMJNV.tblastx.out")],
        conservation_labels=["MelaMJNV"],
    )
    rings = [
        element
        for element in drawing.elements
        if isinstance(element, ContentGroup) and element.key.kind is GroupKind.CONSERVATION
    ]
    assert len(rings) == 1
    assert rings[0].attribs["data-track-label"] == "MelaMJNV"
    assert rings[0].holds_content_children()
