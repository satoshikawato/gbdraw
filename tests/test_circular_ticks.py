from __future__ import annotations

import math
import re
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
    circular_tick_label_margin_px,
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


@pytest.mark.parametrize(
    ("label_side", "tick_side"),
    [("outside", "inside"), ("inside", "outside")],
    ids=["label_out_tick_in", "label_in_tick_out"],
)
def test_circular_tick_labels_stand_one_margin_from_their_ticks_on_either_side(label_side, tick_side) -> None:
    """Both label layouts keep one gap between a label's text box and the tick marks: the label margin
    (OV-3300; labels inside stood a quarter of their height farther)."""

    base_kwargs = {**_tick_label_kwargs(), "label_side": label_side, "tick_side": tick_side}
    margin = circular_tick_label_margin_px(base_kwargs["font_size"], base_kwargs["tick_width"])
    tick_inner, tick_outer = get_circular_tick_path_radius_bounds(
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
        geometry = resolve_circular_tick_label_geometry(**base_kwargs, tick=tick, label_text=label_text)

        gap = tick_inner - geometry.radial_outer_px if label_side == "inside" else geometry.radial_inner_px - tick_outer
        assert gap == pytest.approx(margin)


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


def _slot_arc_radii(svg: str, slot_id: str, path_id_prefix: str = "") -> list[float]:
    """The arc radii of the paths in slot ``slot_id``'s group, of those whose id starts with ``path_id_prefix``."""

    import xml.etree.ElementTree as ET

    (group,) = [
        element for element in ET.fromstring(svg).iter() if element.get("data-gbdraw-slot-id") == slot_id
    ]
    return [
        float(radius)
        for element in group.iter()
        if element.tag.rsplit("}", 1)[-1] == "path" and str(element.get("id") or "").startswith(path_id_prefix)
        for radius in re.findall(r"A ([0-9.]+),[0-9.]+ ", str(element.get("d") or ""))
    ]


def _slot_outer_edge_px(svg: str, slot_id: str) -> float:
    """Where the paint of slot ``slot_id`` ends outward: the largest arc radius in its group, or the farthest
    point of its lines."""

    import xml.etree.ElementTree as ET

    (group,) = [
        element for element in ET.fromstring(svg).iter() if element.get("data-gbdraw-slot-id") == slot_id
    ]
    points = [
        math.hypot(float(x), float(y))
        for element in group.iter()
        if element.tag.rsplit("}", 1)[-1] == "path"
        for x, y in re.findall(r"[ML] ?([-0-9.]+)[ ,]([-0-9.]+)", str(element.get("d") or ""))
    ]
    return max([*_slot_arc_radii(svg, slot_id), *points])


def _tick_label_inner_edge_px(svg: str) -> float:
    """The circle the scale's tick labels stand on inside the ticks: their paths (above the path at the top, below
    it at the bottom)."""

    return min(_slot_arc_radii(svg, "ticks", "circular_tick_label_path_"))


def _ring_below_the_ticks(ring: str, background: bool, slots: tuple[str, ...] | None = None) -> str:
    """A 600 kbp record drawn with ``ring`` (a conservation ring with or without its background, or a feature
    density heatmap or histogram) in the row just inside the scale, whose labels face it; or in the ``slots`` given,
    where ``RING`` stands for the ring's slot (followed by more of its parameters)."""

    from gbdraw.analysis.feature_density import FeatureDensityTrack
    from gbdraw.api.assembly_inputs import CircularAssemblyInputs
    from gbdraw.api.diagram import assemble_circular_diagram_from_record
    from gbdraw.features.query import FeatureQuery
    from gbdraw.io.comparisons import COMPARISON_COLUMNS
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    import pandas as pd

    length = 600_000
    record = SeqRecord(Seq("ACGT" * (length // 4)), id="rec1", name="rec1", description="test record")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "circular"
    record.features = [
        SeqFeature(SimpleLocation(start, start + 900, strand=1), type="CDS") for start in range(0, length, 2_000)
    ]
    hits = pd.DataFrame(
        [("q", "rec1", 90.0, 1000, 0, 0, 1, 1000, start + 1, start + 1000, 1e-20, 100.0) for start in range(0, length, 4_000)],
        columns=COMPARISON_COLUMNS,
    )
    slot = {
        "conservation": "ring:sequence_conservation@source_index=0",
        "heatmap": "ring:feature_density@track=cds",
        "histogram": "ring:feature_density@track=cds_bars",
    }[ring]
    return assemble_circular_diagram_from_record(
        record,
        CircularAssemblyInputs(
            cfg=apply_config_overrides(None, {"objects.conservation.show_background": background}),
            legend="none",
            conservation_dataframes=[hits] if ring == "conservation" else None,
            conservation_reference="subject" if ring == "conservation" else None,
            feature_density_tracks=(
                FeatureDensityTrack("cds", FeatureQuery("CDS")),
                FeatureDensityTrack("cds_bars", FeatureQuery("CDS"), style="histogram"),
            ),
            circular_track_slots=[
                slot + entry[len("RING"):] if entry.startswith("RING") else entry
                for entry in (slots or ("features:features", "ticks:ticks", "RING", "gc_content:dinucleotide_content"))
            ],
        ),
    ).tostring()


@pytest.mark.circular
@pytest.mark.parametrize(
    ("ring", "background"),
    [("conservation", True), ("conservation", False), ("heatmap", True)],
    ids=["conservation-background", "conservation-hits", "density-heatmap"],
)
def test_the_tick_labels_keep_their_margin_from_a_ring_that_paints_its_edge(ring, background) -> None:
    """OV-3149: the tick labels inside the ticks stood on the edge of a conservation ring (its background ring, or
    hits that fill its band's height) or of a feature density heatmap (cells that fill its band's height), as on a
    Percent GC band (OV-3116)."""

    cfg = apply_config_overrides(None, None)
    margin = circular_tick_label_margin_px(float(cfg.objects.ticks.tick_labels.font_size), float(cfg.objects.ticks.tick_width))
    svg = _ring_below_the_ticks(ring, background)

    assert _tick_label_inner_edge_px(svg) - _slot_outer_edge_px(svg, "ring") == pytest.approx(margin, abs=0.01)


@pytest.mark.circular
def test_the_tick_labels_keep_the_band_edge_of_a_density_histogram() -> None:
    """A histogram's bars rise to its band's outer edge only at the highest windows: the tick labels keep the
    packing's band edge beside it, as beside a GC content line."""

    svg = _ring_below_the_ticks("histogram", True)

    assert _tick_label_inner_edge_px(svg) - _slot_outer_edge_px(svg, "ring") == pytest.approx(0.0, abs=0.05)


@pytest.mark.circular
def test_the_tick_labels_keep_their_margin_from_a_feature_track_below_them() -> None:
    """OV-3192: with the scale stacked above the feature track, its labels inside the ticks stood on the feature
    track's outer edge, where the feature blocks end (the "p" of "kbp" touched the arrows). A feature track paints the
    edges of its band, as a Percent GC band does (OV-3116)."""

    cfg = apply_config_overrides(None, None)
    margin = circular_tick_label_margin_px(float(cfg.objects.ticks.tick_labels.font_size), float(cfg.objects.ticks.tick_width))
    svg = _ring_below_the_ticks(
        "heatmap",
        True,
        ("ticks:ticks@tick_label_layout=label_in_tick_out", "features:features", "gc_content:dinucleotide_content"),
    )

    assert _tick_label_inner_edge_px(svg) - _slot_outer_edge_px(svg, "features") == pytest.approx(margin, abs=0.01)


@pytest.mark.circular
def test_outward_tick_labels_keep_their_margin_from_a_feature_track_above_them() -> None:
    """OV-3192 on the other side of the axis: labels outside their ticks keep the larger of the slot gap and their
    margin from a feature track stacked outside them (OV-3150)."""

    cfg = apply_config_overrides(None, None)
    margin = circular_tick_label_margin_px(float(cfg.objects.ticks.tick_labels.font_size), float(cfg.objects.ticks.tick_width))
    svg = _ring_below_the_ticks(
        "heatmap",
        True,
        (
            "features:features@lane_direction=outside,side=outside",
            "ticks:ticks@side=outside,tick_label_layout=label_out_tick_in",
            "gc_content:dinucleotide_content",
        ),
    )

    # The packing keeps the slot gap (3.9 px at this radius) between the labels and a row packed outside them.
    assert min(_slot_arc_radii(svg, "features")) - _tick_label_outer_edge_px(svg) == pytest.approx(max(3.9, margin), abs=0.01)


def _tick_label_outer_edge_px(svg: str) -> float:
    """Where the scale's tick labels end outward when they stand outside the ticks: their paths plus a line box."""

    from gbdraw.core.text import drawn_text_box

    cfg = apply_config_overrides(None, None)
    _width, height = drawn_text_box("0 kbp", str(cfg.objects.text.font_family), float(cfg.objects.ticks.tick_labels.font_size))
    return max(_slot_arc_radii(svg, "ticks", "circular_tick_label_path_")) + float(height)


@pytest.mark.circular
@pytest.mark.parametrize("gap", [None, 1.0, 6.0], ids=["default-gap", "gap-below-margin", "gap-above-margin"])
def test_outward_tick_labels_keep_the_larger_of_the_slot_gap_and_their_margin_from_a_ring(gap) -> None:
    """OV-3150: below a ring that paints its edge, labels standing outside their ticks kept the slot gap plus their
    margin (7.0 px by default); the packing keeps the gap there, so the margin counts only where the gap is less."""

    cfg = apply_config_overrides(None, None)
    margin = circular_tick_label_margin_px(float(cfg.objects.ticks.tick_labels.font_size), float(cfg.objects.ticks.tick_width))
    ring_gap = "" if gap is None else f",inner_gap_px={gap}"
    tick_gap = "" if gap is None else f",outer_gap_px={gap}"
    svg = _ring_below_the_ticks(
        "conservation",
        True,
        (
            "features:features",
            f"RING{ring_gap}",
            f"ticks:ticks@tick_label_layout=label_out_tick_in{tick_gap}",
            "gc_content:dinucleotide_content",
        ),
    )
    clearance = min(_slot_arc_radii(svg, "ring")) - _tick_label_outer_edge_px(svg)
    # The default slot gap at this radius: what dev draws there.
    default_gap = 3.9
    expected = max(default_gap if gap is None else gap, margin)
    assert clearance == pytest.approx(expected, abs=0.01)
