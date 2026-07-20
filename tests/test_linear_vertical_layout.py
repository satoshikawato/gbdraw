from __future__ import annotations

from types import SimpleNamespace

import pytest

from gbdraw.layout.linear import (
    VerticalBand,
    measure_linear_feature_lanes,
    union_vertical_bands,
)
from gbdraw.canvas import LinearCanvasConfigurator
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.toml import load_config_toml
from gbdraw.diagrams.linear.track_slots import (
    LinearSlotFootprint,
    resolve_linear_record_vertical_plan,
    resolve_linear_track_layout,
)
from gbdraw.tracks import normalize_linear_track_slots_with_axis, parse_linear_track_slots


def _feature(track_id: int, strand: str) -> SimpleNamespace:
    return SimpleNamespace(feature_track_id=track_id, strand=strand)


def _base_layout(slot_specs: list[str]):
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_config = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=1200,
        config_dict=config_dict,
        legend="none",
        cfg=cfg,
    )
    slots = normalize_linear_track_slots_with_axis(parse_linear_track_slots(slot_specs), None)
    layout = resolve_linear_track_layout(slots, canvas_config=canvas_config, cfg=cfg)
    return layout, canvas_config


def test_vertical_band_geometry_operations() -> None:
    band = VerticalBand(-4.0, 6.0)

    assert band.height == pytest.approx(10.0)
    assert band.translate(3.0) == VerticalBand(-1.0, 9.0)
    assert band.expand(2.0, 4.0) == VerticalBand(-6.0, 10.0)
    assert band.union(VerticalBand(-8.0, 2.0)) == VerticalBand(-8.0, 6.0)
    assert band.intersects(VerticalBand(5.0, 8.0))
    assert not band.intersects(VerticalBand(6.0, 8.0))


@pytest.mark.parametrize(
    ("top_y", "bottom_y"),
    [
        (1.0, -1.0),
        (float("nan"), 1.0),
        (0.0, float("inf")),
    ],
)
def test_vertical_band_rejects_invalid_coordinates(top_y: float, bottom_y: float) -> None:
    with pytest.raises(ValueError):
        VerticalBand(top_y, bottom_y)


def test_union_vertical_bands_uses_explicit_empty_default() -> None:
    default = VerticalBand(0.0, 0.0)

    assert union_vertical_bands([], default=default) == default
    with pytest.raises(ValueError, match="at least one"):
        union_vertical_bands([])


def test_measure_middle_separated_feature_lanes_uses_renderer_factors() -> None:
    geometry = measure_linear_feature_lanes(
        {
            "positive": _feature(0, "positive"),
            "negative-1": _feature(-1, "negative"),
            "negative-2": _feature(-2, "negative"),
            "negative-5": _feature(-5, "negative"),
        },
        cds_height=20.0,
        separate_strands=True,
        track_layout="middle",
        stroke_width=2.0,
    )

    lanes = {(lane.strand_pool, lane.track_id): lane.band for lane in geometry.lanes}
    assert (lanes[("positive", 0)].top_y, lanes[("positive", 0)].bottom_y) == pytest.approx((-13.0, -1.0))
    assert (lanes[("negative", -1)].top_y, lanes[("negative", -1)].bottom_y) == pytest.approx((1.0, 13.0))
    assert (lanes[("negative", -2)].top_y, lanes[("negative", -2)].bottom_y) == pytest.approx((14.0, 26.0))
    assert (lanes[("negative", -5)].top_y, lanes[("negative", -5)].bottom_y) == pytest.approx((53.0, 65.0))
    assert (geometry.occupied_band.top_y, geometry.occupied_band.bottom_y) == pytest.approx((-13.0, 65.0))
    assert geometry.lane_for(
        strand="positive",
        track_id=0,
        separate_strands=True,
    ).positions == pytest.approx((-12.0, -7.0, -2.0))


@pytest.mark.parametrize(("layout", "sign"), [("above", -1.0), ("below", 1.0)])
def test_measure_feature_lanes_tracks_selected_side(layout: str, sign: float) -> None:
    geometry = measure_linear_feature_lanes(
        {"feature": _feature(0, "positive")},
        cds_height=20.0,
        separate_strands=False,
        track_layout=layout,
        axis_gap=6.0,
    )

    band = geometry.occupied_band
    assert band.top_y * sign >= 0.0 if sign > 0 else band.bottom_y <= 0.0
    if layout == "above":
        assert band.bottom_y == pytest.approx(-6.0)
    else:
        assert band.top_y == pytest.approx(6.0)


def test_record_planner_reserves_middle_feature_band_before_depth() -> None:
    layout, canvas_config = _base_layout(
        [
            "features:features@side=overlay",
            "depth:depth@side=below,h=20px",
        ]
    )
    feature_band = VerticalBand(-12.0, 30.0)
    plan = resolve_linear_record_vertical_plan(
        layout,
        axis_band=VerticalBand(-0.5, 0.5),
        footprints={
            "features": LinearSlotFootprint(feature_band, feature_band),
            "depth": LinearSlotFootprint(VerticalBand(0.0, 20.0), VerticalBand(0.0, 20.0)),
        },
    )

    features = plan.slot_by_id("features")
    depth = plan.slot_by_id("depth")
    assert not features.reserve_band.intersects(depth.reserve_band)
    assert depth.reserve_band.top_y - features.reserve_band.bottom_y == pytest.approx(
        canvas_config.vertical_padding
    )
    assert plan.record_body_band == VerticalBand(-12.0, depth.reserve_band.bottom_y)


def test_record_planner_missing_depth_has_no_paint_or_reserve() -> None:
    layout, _canvas_config = _base_layout(
        [
            "features:features@side=overlay",
            "depth:depth@side=below,h=20px",
        ]
    )
    plan = resolve_linear_record_vertical_plan(
        layout,
        axis_band=VerticalBand(-0.5, 0.5),
        footprints={
            "features": LinearSlotFootprint(VerticalBand(-5.0, 5.0), VerticalBand(-5.0, 5.0)),
            "depth": LinearSlotFootprint(
                paint_band=None,
                reserve_band=VerticalBand(0.0, 0.0),
                data_available=False,
            ),
        },
    )

    depth = plan.slot_by_id("depth")
    assert not depth.data_available
    assert depth.paint_band is None
    assert depth.reserve_band.height == pytest.approx(0.0)
    assert plan.record_body_band == VerticalBand(-5.0, 5.0)


def test_record_planner_comparison_band_uses_paint_not_empty_reserve() -> None:
    layout, _canvas_config = _base_layout(["features:features@side=overlay,h=40px"])
    plan = resolve_linear_record_vertical_plan(
        layout,
        axis_band=VerticalBand(-0.5, 0.5),
        footprints={
            "features": LinearSlotFootprint(
                paint_band=VerticalBand(-5.0, 5.0),
                reserve_band=VerticalBand(-20.0, 20.0),
            ),
        },
    )

    assert plan.record_body_band == VerticalBand(-20.0, 20.0)
    assert plan.comparison_exclusion_band == VerticalBand(-5.0, 5.0)


def test_record_planner_packs_same_side_feature_and_depth_in_slot_order() -> None:
    layout, canvas_config = _base_layout(
        [
            "features:features@side=below",
            "depth:depth@side=below,h=20px",
        ]
    )
    plan = resolve_linear_record_vertical_plan(
        layout,
        axis_band=VerticalBand(-0.5, 0.5),
        footprints={
            "features": LinearSlotFootprint(VerticalBand(5.0, 25.0), VerticalBand(5.0, 25.0)),
            "depth": LinearSlotFootprint(VerticalBand(0.0, 20.0), VerticalBand(0.0, 20.0)),
        },
    )

    features = plan.slot_by_id("features")
    depth = plan.slot_by_id("depth")
    assert depth.reserve_band.top_y - features.reserve_band.bottom_y >= canvas_config.vertical_padding


def test_record_planner_z_only_changes_paint_order_metadata() -> None:
    first, _canvas_config = _base_layout(
        [
            "content:gc_content@side=below,h=20px,z=0",
            "depth:depth@side=below,h=20px,z=5",
        ]
    )
    second, _canvas_config = _base_layout(
        [
            "content:gc_content@side=below,h=20px,z=9",
            "depth:depth@side=below,h=20px,z=-2",
        ]
    )
    footprints = {
        "content": LinearSlotFootprint(VerticalBand(-10.0, 10.0), VerticalBand(-10.0, 10.0)),
        "depth": LinearSlotFootprint(VerticalBand(0.0, 20.0), VerticalBand(0.0, 20.0)),
    }

    first_plan = resolve_linear_record_vertical_plan(
        first,
        axis_band=VerticalBand(-0.5, 0.5),
        footprints=footprints,
    )
    second_plan = resolve_linear_record_vertical_plan(
        second,
        axis_band=VerticalBand(-0.5, 0.5),
        footprints=footprints,
    )

    assert [slot.origin_y for slot in first_plan.slots] == pytest.approx(
        [slot.origin_y for slot in second_plan.slots]
    )
    assert [slot.z for slot in first_plan.slots] != [slot.z for slot in second_plan.slots]


@pytest.mark.parametrize("side", ["above", "below"])
def test_record_planner_non_overlay_reservations_do_not_intersect(side: str) -> None:
    layout, _canvas_config = _base_layout(
        [
            f"depth:depth@side={side},h=20px,spacing=7px",
            f"content:gc_content@side={side},h=18px,spacing=3px",
            f"skew:gc_skew@side={side},h=16px",
        ]
    )
    footprints = {
        "depth": LinearSlotFootprint(VerticalBand(0.0, 20.0), VerticalBand(-2.0, 22.0)),
        "content": LinearSlotFootprint(VerticalBand(-9.0, 9.0), VerticalBand(-10.0, 10.0)),
        "skew": LinearSlotFootprint(VerticalBand(-8.0, 8.0), VerticalBand(-9.0, 9.0)),
    }
    plan = resolve_linear_record_vertical_plan(
        layout,
        axis_band=VerticalBand(-0.5, 0.5),
        footprints=footprints,
    )
    ordered = sorted(plan.slots, key=lambda slot: slot.reserve_band.top_y)

    for inner, outer in zip(ordered, ordered[1:]):
        assert not inner.reserve_band.intersects(outer.reserve_band)


def test_record_planner_spacer_reserves_distance_without_paint() -> None:
    layout, _canvas_config = _base_layout(
        [
            "features:features@side=overlay",
            "gap:spacer@side=below,h=12px",
            "depth:depth@side=below,h=20px",
        ]
    )
    plan = resolve_linear_record_vertical_plan(
        layout,
        axis_band=VerticalBand(-0.5, 0.5),
        footprints={
            "features": LinearSlotFootprint(VerticalBand(-5.0, 5.0), VerticalBand(-5.0, 5.0)),
            "gap": LinearSlotFootprint(None, VerticalBand(0.0, 12.0), data_available=False),
            "depth": LinearSlotFootprint(VerticalBand(0.0, 20.0), VerticalBand(0.0, 20.0)),
        },
    )
    spacer = plan.slot_by_id("gap")
    depth = plan.slot_by_id("depth")

    assert spacer.paint_band is None
    assert spacer.reserve_band.height == pytest.approx(12.0)
    assert depth.reserve_band.top_y >= spacer.reserve_band.bottom_y


def test_record_planner_overlay_expansion_pushes_next_track_outward() -> None:
    layout, _canvas_config = _base_layout(
        [
            "features:features@side=overlay",
            "depth:depth@side=below,h=20px",
            "note:annotations@side=overlay,h=30px,anchor_slot=depth,set_id=note,z=1",
            "content:gc_content@side=below,h=20px",
        ]
    )
    plan = resolve_linear_record_vertical_plan(
        layout,
        axis_band=VerticalBand(-0.5, 0.5),
        footprints={
            "features": LinearSlotFootprint(VerticalBand(-5.0, 5.0), VerticalBand(-5.0, 5.0)),
            "depth": LinearSlotFootprint(VerticalBand(0.0, 20.0), VerticalBand(0.0, 20.0)),
            "note": LinearSlotFootprint(VerticalBand(-10.0, 40.0), VerticalBand(-10.0, 40.0)),
            "content": LinearSlotFootprint(VerticalBand(-10.0, 10.0), VerticalBand(-10.0, 10.0)),
        },
    )
    depth = plan.slot_by_id("depth")
    note = plan.slot_by_id("note")
    content = plan.slot_by_id("content")
    composite_bottom = max(depth.reserve_band.bottom_y, note.reserve_band.bottom_y)

    assert depth.reserve_band.intersects(note.reserve_band)
    assert content.reserve_band.top_y >= composite_bottom
