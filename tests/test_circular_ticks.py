from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gbdraw
from gbdraw.api.config import apply_config_overrides
from gbdraw.diagrams.circular.radial_layout import track_layout_reason
from gbdraw.auto_sizes import tick_and_window_tiers
from gbdraw.exceptions import ValidationError
from gbdraw.labels.label_style import ROW_LABEL_GAP_EM
from gbdraw.render.groups.circular import ticks as ticks_module
from gbdraw.render.groups.circular.ticks import TickGroup
from gbdraw.render.groups.linear.length_bar import format_linear_tick_label
from gbdraw.svg.circular_ticks import (
    CircularTickLabelGeometry,
    circular_tick_labels,
    generate_circular_tick_labels,
    get_circular_tick_path_radius_bounds,
    resolve_circular_tick_label_geometry,
    set_tick_label_anchor_value,
)

# The packaged labels.length_threshold.circular.
LENGTH_THRESHOLD = int(apply_config_overrides(None, None).labels.length_threshold.circular)


def _tick_label_kwargs() -> dict:
    return {
        "center_radius_px": 100.0,
        "total_len": 4_641_652,
        "size": "large",
        "font_size": 14.0,
        "font_family": "Arial",
        "track_type": "tuckin",
        "strandedness": True,
        "label_side": "outside",
        "tick_side": "inside",
        "tick_length_px": 10.0,
        "tick_width": 2.0,
        "length_threshold": LENGTH_THRESHOLD,
    }


def test_circular_tick_label_anchor_stays_middle_and_baseline_uses_edge_values_by_angle() -> None:
    assert set_tick_label_anchor_value(4_641_652, 1_000_000) == ("middle", "text-after-edge")
    assert set_tick_label_anchor_value(4_641_652, 2_500_000) == ("middle", "text-before-edge")
    assert set_tick_label_anchor_value(4_641_652, 3_000_000) == ("middle", "text-before-edge")
    assert set_tick_label_anchor_value(4_641_652, 4_500_000) == ("middle", "text-after-edge")


def test_circular_tick_label_geometry_uses_middle_anchor() -> None:
    base_kwargs = _tick_label_kwargs()

    right_side = resolve_circular_tick_label_geometry(
        **base_kwargs,
        tick=1_000_000,
        label_text="1.0 Mbp",
    )
    left_side = resolve_circular_tick_label_geometry(
        **base_kwargs,
        tick=3_000_000,
        label_text="3.0 Mbp",
    )
    top_side = resolve_circular_tick_label_geometry(
        **base_kwargs,
        tick=4_500_000,
        label_text="4.5 Mbp",
    )

    assert right_side.text_anchor == "middle"
    assert left_side.text_anchor == "middle"
    assert top_side.text_anchor == "middle"


def test_circular_tick_label_outside_uses_margin_from_tick_band() -> None:
    base_kwargs = _tick_label_kwargs()
    expected_gap = max(2.0, base_kwargs["font_size"] * 0.15) + (base_kwargs["tick_width"] / 2.0)
    _tick_inner, tick_outer = get_circular_tick_path_radius_bounds(
        center_radius_px=base_kwargs["center_radius_px"],
        total_len=base_kwargs["total_len"],
        size=base_kwargs["size"],
        track_type=base_kwargs["track_type"],
        strandedness=base_kwargs["strandedness"],
        tick_side=base_kwargs["tick_side"],
        tick_length_px=base_kwargs["tick_length_px"],
        length_threshold=LENGTH_THRESHOLD,
    )

    for tick, label_text in ((1_000_000, "1.0 Mbp"), (3_000_000, "3.0 Mbp")):
        geometry = resolve_circular_tick_label_geometry(
            **base_kwargs,
            tick=tick,
            label_text=label_text,
        )

        assert geometry.radial_inner_px - tick_outer == pytest.approx(expected_gap)


def test_circular_tick_label_inside_uses_legacy_center_offset_clearance() -> None:
    base_kwargs = {
        **_tick_label_kwargs(),
        "label_side": "inside",
        "tick_side": "outside",
    }
    expected_margin = max(2.0, base_kwargs["font_size"] * 0.15) + (base_kwargs["tick_width"] / 2.0)
    tick_inner, _tick_outer = get_circular_tick_path_radius_bounds(
        center_radius_px=base_kwargs["center_radius_px"],
        total_len=base_kwargs["total_len"],
        size=base_kwargs["size"],
        track_type=base_kwargs["track_type"],
        strandedness=base_kwargs["strandedness"],
        tick_side=base_kwargs["tick_side"],
        tick_length_px=base_kwargs["tick_length_px"],
        length_threshold=LENGTH_THRESHOLD,
    )

    path_radii = []
    for tick, label_text in ((1_000_000, "1.0 Mbp"), (3_000_000, "3.0 Mbp")):
        geometry = resolve_circular_tick_label_geometry(
            **base_kwargs,
            tick=tick,
            label_text=label_text,
        )

        path_radii.append(geometry.path_radius_px)
        expected_gap = expected_margin + (geometry.bbox_height_px / 4.0)
        assert tick_inner - geometry.radial_outer_px == pytest.approx(expected_gap)
    assert path_radii[0] == pytest.approx(path_radii[1])


def test_circular_tick_label_textpath_uses_middle_anchor() -> None:
    elements = generate_circular_tick_labels(
        100.0,
        4_641_652,
        "large",
        [1_000_000, 3_000_000, 4_500_000],
        "none",
        "black",
        14.0,
        "normal",
        "Arial",
        "tuckin",
        True,
        label_side="outside",
        tick_side="inside",
        tick_length_px=10.0,
        tick_width=2.0,
        length_threshold=LENGTH_THRESHOLD,
    )
    text_elements = [element.tostring() for element in elements if element.elementname == "text"]

    assert 'text-anchor="middle"' in text_elements[0]
    assert 'startOffset="50%"' in text_elements[0]
    assert 'dominant-baseline="text-after-edge"' in text_elements[0]
    assert 'text-anchor="middle"' in text_elements[1]
    assert 'startOffset="50%"' in text_elements[1]
    assert 'dominant-baseline="text-before-edge"' in text_elements[1]
    assert 'text-anchor="middle"' in text_elements[2]
    assert 'startOffset="50%"' in text_elements[2]
    assert 'dominant-baseline="text-after-edge"' in text_elements[2]


def test_circular_tick_labels_resolve_an_interval_below_one_kbp() -> None:
    # OV-201: a 500 bp interval read "0 kbp, 1 kbp, 1 kbp"; the labels now keep
    # the decimals the interval needs, as Linear tick labels do.
    elements = generate_circular_tick_labels(
        100.0,
        16_569,
        "large",
        [500, 1_000, 1_500, 2_250],
        "none",
        "black",
        14.0,
        "normal",
        "Arial",
        "tuckin",
        True,
        label_side="outside",
        tick_side="inside",
        tick_length_px=10.0,
        tick_width=2.0,
        tick_interval=250,
        length_threshold=LENGTH_THRESHOLD,
    )
    texts = [element.tostring() for element in elements if element.elementname == "text"]

    assert [">0.5 kbp<" in texts[0], ">1 kbp<" in texts[1], ">1.5 kbp<" in texts[2], ">2.25 kbp<" in texts[3]] == [
        True,
        True,
        True,
        True,
    ]


@pytest.mark.parametrize(
    ("position", "context_length", "interval", "expected"),
    [
        (2_250, 16_569, 250, "2.25 kbp"),
        (4_500, 16_569, 1_500, "4.5 kbp"),
        (1_200, 1_500, 100, "1.2 kbp"),
        (1_500_000, 4_641_652, 500_000, "1.5 Mbp"),
        (40_000, 48_502, 5_000, "40 kbp"),
    ],
)
def test_linear_tick_labels_resolve_their_interval(
    position: int, context_length: int, interval: int, expected: str
) -> None:
    # OV-201: Linear rounded 2250 to "2.2 kbp" (interval 250) and 4500 to
    # "4 kbp" (interval 1500); automatic intervals keep their labels.
    assert (
        format_linear_tick_label(position, context_length=context_length, tick_interval=interval)
        == expected
    )


@pytest.mark.parametrize(("radius", "label_side"), [(5.0, "outside"), (20.0, "inside")])
def test_a_tick_label_longer_than_its_circle_is_a_track_that_does_not_fit(radius, label_side) -> None:
    """OV-2559: a scale label whose path would wrap its circle more than once (outside labels at radius 5; inside
    labels at radius 20, whose circle is about 1 px) is refused as a track that does not fit, which a requested
    figure size reads as a radius too small."""

    with pytest.raises(ValidationError) as caught:
        generate_circular_tick_labels(
            radius,
            16_569,
            "large",
            [5_000],
            "none",
            "black",
            14.0,
            "normal",
            "Arial",
            "tuckin",
            True,
            label_side=label_side,
            tick_side="inside",
            tick_length_px=10.0,
            tick_width=2.0,
            tick_interval=5_000,
            length_threshold=LENGTH_THRESHOLD,
        )
    assert track_layout_reason(caught.value) == "CANNOT_FIT"
    assert "5 kbp" in str(caught.value)


def test_a_map_too_small_for_its_tick_labels_is_refused_as_a_track_that_does_not_fit() -> None:
    """OV-2559: the mitochondrial map with labels on both sides at radius 25 raised a bare ValueError from the arc
    of a tick label; it is now the track layout refusal, with the row of the tick slot."""

    records = gbdraw.read_genbank(Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk")
    with pytest.raises(ValidationError) as caught:
        gbdraw.draw_circular(
            records,
            options=gbdraw.CircularOptions(
                config_overrides={"canvas.circular.radius": 25, "labels.circular.scope": "both"}
            ),
        )
    assert track_layout_reason(caught.value) == "CANNOT_FIT"
    assert caught.value.diagnostic["slotIndex"] >= 0
    assert caught.value.diagnostic["slotTrack"] == "ticks"


def test_radial_labels_refuse_scale_labels_longer_than_their_circle_when_the_slot_is_drawn() -> None:
    """OV-2559 with radial labels, which can move the tracks: the early check is skipped, and drawing the tick slot
    refuses with its row instead of the bare ValueError of the arc."""

    records = gbdraw.read_genbank(Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk")
    slots = ["features:features", "ticks:ticks@tick_label_layout=label_only", "gc_content:dinucleotide_content"]
    with pytest.raises(ValidationError) as caught:
        gbdraw.draw_circular(
            records,
            options=gbdraw.CircularOptions(
                config_overrides={
                    "canvas.circular.radius": 25,
                    "labels.circular.scope": "both",
                    "labels.circular.placement": "radial",
                    "objects.definition.circular.show": False,
                },
                tracks=gbdraw.CircularTrackOptions(slots=slots),
            ),
        )
    assert track_layout_reason(caught.value) == "CANNOT_FIT"
    assert caught.value.diagnostic["slotIndex"] == 1
    assert caught.value.diagnostic["slotTrack"] == "ticks"


# OV-2747: a 112,934,049 bp record took 1 Mbp ticks, and its 112 labels merged into a black band.
_LONG_RECORD_BP = 112_934_049


@pytest.fixture(scope="module")
def long_record_tick_group() -> tuple[TickGroup, tuple, dict]:
    """The tick group of the default Circular figure of a 112,934,049 bp record, with the arguments it was made
    with: the default figure gives every long record the same tick slot."""

    made: list[tuple[TickGroup, tuple, dict]] = []
    record = SeqRecord(Seq("A" * _LONG_RECORD_BP), id="long", annotations={"molecule_type": "DNA"})
    original_init = TickGroup.__init__

    def recording_init(self: TickGroup, *args: object, **kwargs: object) -> None:
        original_init(self, *args, **kwargs)
        made.append((self, args, kwargs))

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(TickGroup, "__init__", recording_init)
        gbdraw.draw_circular(
            [record],
            options=gbdraw.CircularOptions(tracks=gbdraw.CircularTrackOptions(slots=["features:features", "ticks:ticks"])),
        )
    return made[-1]


def _label_gap_shortfall_px(group: TickGroup) -> float:
    """How far the closest two neighbouring large tick labels stay from a quarter em apart (OV-2748): their gap
    along the circle at the smaller of their path radii, minus ``ROW_LABEL_GAP_EM`` of the label font. Negative:
    they are closer, so they read as one. The labels are those ``TickGroup.ink()`` draws."""

    drawn: list[tuple[str, float, CircularTickLabelGeometry]] = []

    def recording_labels(*args: object, **kwargs: object) -> list[tuple[str, float, CircularTickLabelGeometry]]:
        labels = circular_tick_labels(*args, **kwargs)
        drawn.extend(labels)
        return labels

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(ticks_module, "circular_tick_labels", recording_labels)
        group.ink()
    arcs = sorted(
        (
            float(tick) / group.total_len - half_turn,
            float(tick) / group.total_len + half_turn,
            radius,
        )
        for _text, tick, geometry in drawn
        for radius in [max(1.0, geometry.path_radius_px)]
        for half_turn in [0.5 * geometry.bbox_width_px / radius / (2.0 * math.pi)]
    )
    assert len(arcs) > 1
    first_again = (arcs[0][0] + 1.0, arcs[0][1] + 1.0, arcs[0][2])
    gap_px = min(
        (right[0] - left[1]) * 2.0 * math.pi * min(left[2], right[2])
        for left, right in zip(arcs, [*arcs[1:], first_again])
    )
    return gap_px - ROW_LABEL_GAP_EM * float(group.font_size)


def _tick_group_of_length(captured: tuple[TickGroup, tuple, dict], length: int, interval: int | None = None) -> TickGroup:
    """The captured default tick group rebuilt for a record of ``length`` bp (the group reads only the length),
    with Auto ticks or a forced large ``interval``."""

    _group, args, kwargs = captured
    kwargs = {**kwargs, "record_transform": None}
    if interval is not None:
        profile = kwargs["profile"]
        cfg = profile.config
        scale = replace(cfg.objects.scale, interval=interval)
        kwargs["profile"] = replace(profile, config=replace(cfg, objects=replace(cfg.objects, scale=scale)))
    return TickGroup(SeqRecord(Seq(None, length=length), id="long"), *args[1:], **kwargs)


def test_the_tick_labels_of_a_113_mbp_record_keep_a_quarter_em_apart(long_record_tick_group) -> None:
    group = long_record_tick_group[0]

    assert _LONG_RECORD_BP // group.tick_large < 40
    assert _label_gap_shortfall_px(group) >= 0.0


_TIERS = tick_and_window_tiers(apply_config_overrides(None, None))["circularTicks"]


@pytest.mark.parametrize(("at_most", "large", "_small"), _TIERS["tiers"])
def test_the_tick_labels_of_every_tier_keep_a_quarter_em_apart_at_its_longest_record(
    long_record_tick_group, at_most, large, _small
) -> None:
    group = _tick_group_of_length(long_record_tick_group, at_most)

    assert group.tick_large == large
    assert _label_gap_shortfall_px(group) >= 0.0


@pytest.mark.parametrize(
    ("at_most", "large"), [(at_most, large) for at_most, large, _small in _TIERS["tiers"] if at_most > 10_000_000]
)
def test_a_tier_above_10_mbp_ends_where_its_tick_labels_come_closer_than_a_quarter_em(
    long_record_tick_group, at_most, large
) -> None:
    """The rows above 10 Mbp are measured, not chosen: each ends at the longest record whose labels keep a quarter
    em apart, so only figures whose labels would come closer take a coarser interval."""

    assert _tick_group_of_length(long_record_tick_group, at_most, large).label_stride == 1
    assert _tick_group_of_length(long_record_tick_group, at_most + 1, large).label_stride > 1


def test_the_ticks_above_the_last_tier_fit_records_of_10_gbp(long_record_tick_group) -> None:
    group = _tick_group_of_length(long_record_tick_group, 10_000_000_000)

    assert group.tick_large == _TIERS["above"][0]
    assert _label_gap_shortfall_px(group) >= 0.0


# OV-3020: the 5,498,578 bp E. coli O157:H7 Sakai chromosome with ten comparison rings drew its tick circle at radius
# 117 px, and its 0.5 Mbp labels printed into each other ("4.5 Mbp5.0 Mbp"). A small figure puts the ticks there too.
_SAKAI_BP = 5_498_578


def _drawn_tick_group(radius: float, slot: str = "ticks:ticks", interval: int | None = None) -> TickGroup:
    """The tick group that ``draw_circular`` draws for a featureless record of ``_SAKAI_BP`` at ``radius``."""

    made: list[TickGroup] = []
    original_init = TickGroup.__init__

    def recording_init(self: TickGroup, *args: object, **kwargs: object) -> None:
        original_init(self, *args, **kwargs)
        made.append(self)

    record = SeqRecord(Seq("ACGT" * (_SAKAI_BP // 4) + "AC"), id="sakai", annotations={"molecule_type": "DNA"})
    overrides: dict[str, object] = {"canvas.circular.radius": radius}
    if interval is not None:
        overrides["objects.scale.interval"] = interval
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(TickGroup, "__init__", recording_init)
        gbdraw.draw_circular(
            [record],
            options=gbdraw.CircularOptions(
                config_overrides=overrides, tracks=gbdraw.CircularTrackOptions(slots=["features:features", slot])
            ),
        )
    return made[-1]


def _drawn_label_texts(group: TickGroup) -> list[str]:
    return [
        text_path.text
        for text in group.get_group().elements
        if text.elementname == "text"
        for text_path in text.elements
    ]


@pytest.mark.parametrize(
    ("radius", "slot"),
    [
        (120.0, "ticks:ticks"),
        (90.0, "ticks:ticks"),
        (90.0, "ticks:ticks@tick_label_layout=label_in_tick_out"),
        (90.0, "ticks:ticks@tick_label_layout=label_out_tick_in"),
        (90.0, "ticks:ticks@tick_label_layout=label_only"),
    ],
)
def test_tick_labels_on_a_small_circle_leave_labels_out_instead_of_printing_into_each_other(radius, slot) -> None:
    """OV-3020: where the labels of every tick would come closer than a quarter em, every second (then fifth, tenth,
    ...) tick keeps its label; every tick keeps its mark at the interval."""

    group = _drawn_tick_group(radius, slot)

    assert group.tick_large == 500_000
    assert group.label_stride > 1
    assert _label_gap_shortfall_px(group) >= 0.0
    marks = [element for element in group.get_group().elements if element.elementname == "path" and element["stroke"] == "gray"]
    assert len(marks) == (0 if slot.endswith("label_only") else 11)
    texts = _drawn_label_texts(group)
    assert texts and all(float(text.split()[0]) * 1_000_000 % (500_000 * group.label_stride) == 0 for text in texts)


def test_a_typed_interval_keeps_its_ticks_and_leaves_out_labels_that_would_print_into_each_other() -> None:
    group = _drawn_tick_group(120.0, interval=250_000)

    assert group.tick_large == 250_000
    assert group.label_stride > 1
    assert _label_gap_shortfall_px(group) >= 0.0


def test_tick_labels_that_keep_apart_are_all_drawn() -> None:
    group = _drawn_tick_group(390.0)

    assert group.label_stride == 1
    assert _drawn_label_texts(group) == [f"{tick / 1_000_000} Mbp" for tick in range(500_000, _SAKAI_BP, 500_000)]
