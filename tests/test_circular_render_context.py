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
from gbdraw.render.groups.circular.labels import LabelsGroup
from gbdraw.render.groups.circular.seq_record import SeqRecordGroup


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
