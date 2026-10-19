"""The shared label backing (LP-D82 region labels, LP-D85 study prototype for feature labels)."""

from __future__ import annotations

import pytest

from gbdraw.layout.spatial import Aabb
from gbdraw.svg.elements import Text
from gbdraw.core.text import get_font_vertical_metrics
from gbdraw.labels.label_ink import (
    GAP_STOP_EM,
    UNDERLINE_GAP_EM,
    FeatureLabelInk,
    OverlineMode,
    label_underline_counted,
    label_underline_in_box,
    label_underline_on,
    line_half_width_em,
    rule_leader_width,
    study_switch,
    underline_reach_em,
)
from gbdraw.render.drawers.label_backing import (
    BACKING_PADDING_EM,
    backed_feature_labels,
    backing_box,
    feature_label_backing,
    feature_label_backing_mode,
    feature_label_backing_opacity,
    feature_label_halo_em,
    label_halos,
    leaders_through_text,
)


def test_backing_box_pads_and_stops_at_an_edge_on_its_side() -> None:
    ink = Aabb(0.0, 0.0, 50.0, 10.0)
    pad = BACKING_PADDING_EM * 10.0
    assert backing_box(ink, 10.0) == Aabb(-pad, -pad, 50.0 + pad, 10.0 + pad)
    assert backing_box(ink, 10.0, edge_y=11.0).max_y == 11.0
    assert backing_box(ink, 10.0, edge_y=-1.0).min_y == -1.0
    assert backing_box(ink, 10.0, edge_x=51.0) == Aabb(-pad, -pad, 51.0, 10.0 + pad)


def _label(ink: Aabb, leader: tuple[tuple[float, float], ...]) -> FeatureLabelInk:
    return FeatureLabelInk(ink=ink, font_size=10.0, rotation_deg=0.0, center=(ink.min_x, ink.min_y), leader=leader)


def test_through_backs_only_the_label_another_leader_passes_through() -> None:
    low = _label(Aabb(0.0, 90.0, 40.0, 100.0), ((20.0, 120.0), (20.0, 100.0)))
    high = _label(Aabb(0.0, 60.0, 40.0, 70.0), ((30.0, 120.0), (30.0, 70.0)))  # its leader crosses ``low``
    labels = [low, high]
    assert leaders_through_text(labels) == [(1, 0)]
    assert backed_feature_labels(labels, "through") == frozenset({0})
    assert backed_feature_labels(labels, "all") == frozenset({0, 1})
    assert backed_feature_labels(labels, "off") == frozenset()


def test_a_backing_stops_at_its_own_leader_end() -> None:
    label = _label(Aabb(0.0, 90.0, 40.0, 100.0), ((20.0, 120.0), (20.0, 100.5)))
    assert feature_label_backing(label).max_y == pytest.approx(100.5)


@pytest.mark.parametrize(("value", "mode"), [("", "off"), ("through", "through"), ("ALL", "all"), ("x", "off")])
def test_prototype_switch(monkeypatch: pytest.MonkeyPatch, value: str, mode: str) -> None:
    monkeypatch.setenv("GBDRAW_LABEL_BACKING", value)
    assert feature_label_backing_mode("linear") == mode


def test_underline_continues_the_leader_and_is_part_of_it() -> None:
    # A leader from below (Linear, label above its feature) meets the line under the letters in a T.
    label = FeatureLabelInk(
        ink=Aabb(0.0, 90.0, 40.0, 100.0), font_size=10.0, rotation_deg=0.0, center=(0.0, 90.0),
        leader=((20.0, 130.0), (20.0, 100.5)), under_y=103.0, over_y=88.0,
    )
    assert label.underline == ((20.0, 103.0), (0.0, 103.0), (40.0, 103.0))
    assert label.drawn == ((20.0, 130.0), (20.0, 103.0), (0.0, 103.0), (40.0, 103.0))
    # The backing stops one gap short of the line, so it never paints over the line's stroke.
    assert feature_label_backing(label).max_y == pytest.approx(100.0 + BACKING_PADDING_EM * 10.0)
    tight = FeatureLabelInk(
        ink=Aabb(0.0, 90.0, 40.0, 100.0), font_size=10.0, rotation_deg=0.0, center=(0.0, 90.0),
        leader=((20.0, 130.0), (20.0, 100.5)), under_y=101.0, over_y=88.0,
    )
    assert feature_label_backing(tight).max_y == pytest.approx(101.0 - UNDERLINE_GAP_EM * 10.0)
    # A leader from the side (Circular) meets the line at its near end; from above, the line is over the letters.
    side = FeatureLabelInk(
        ink=Aabb(0.0, 90.0, 40.0, 100.0), font_size=10.0, rotation_deg=0.0, center=(0.0, 95.0),
        leader=((-30.0, 60.0), (-10.0, 95.0), (-1.0, 95.0)), under_y=103.0, over_y=88.0,
    )
    assert side.drawn[-3:] == ((-10.0, 95.0), (0.0, 103.0), (40.0, 103.0))
    above = FeatureLabelInk(
        ink=Aabb(0.0, 90.0, 40.0, 100.0), font_size=10.0, rotation_deg=0.0, center=(0.0, 95.0),
        leader=((20.0, 50.0), (20.0, 89.5)), under_y=103.0, over_y=88.0,
    )
    assert above.underline is not None and above.underline[0] == (20.0, 88.0)
    # A leader from above but beside the letters (a Circular inner label near the top) still meets the line under them.
    beside = FeatureLabelInk(
        ink=Aabb(0.0, 90.0, 40.0, 100.0), font_size=10.0, rotation_deg=0.0, center=(0.0, 95.0),
        leader=((60.0, 50.0), (41.0, 95.0)), under_y=103.0, over_y=88.0,
    )
    assert beside.underline == ((40.0, 103.0), (40.0, 103.0), (0.0, 103.0))
    # A Circular leader keeps its bends: its last segment is re-aimed at the line's near end, never a T (LP-D90).
    circular = FeatureLabelInk(
        ink=Aabb(0.0, 90.0, 40.0, 100.0), font_size=10.0, rotation_deg=0.0, center=(0.0, 95.0),
        leader=((15.0, 140.0), (10.0, 120.0), (5.0, 100.5)), under_y=103.0, over_y=88.0, join_at_end=True,
    )
    assert circular.drawn == ((15.0, 140.0), (10.0, 120.0), (0.0, 103.0), (40.0, 103.0))
    assert len(circular.drawn) - 1 == len(circular.leader)
    # It meets the end on the side where the leader ends today, even when its bend lies nearer the other end.
    anchored = FeatureLabelInk(
        ink=Aabb(0.0, 90.0, 40.0, 100.0), font_size=10.0, rotation_deg=0.0, center=(0.0, 95.0),
        leader=((60.0, 140.0), (45.0, 110.0), (-0.5, 95.0)), under_y=103.0, over_y=88.0, join_at_end=True,
    )
    assert anchored.drawn == ((60.0, 140.0), (45.0, 110.0), (0.0, 103.0), (40.0, 103.0))
    # From above, re-aimed under the letters the leader would run through them: the line goes over them instead.
    over = FeatureLabelInk(
        ink=Aabb(0.0, 90.0, 40.0, 100.0), font_size=10.0, rotation_deg=0.0, center=(0.0, 95.0),
        leader=((60.0, 60.0), (30.0, 80.0), (40.5, 95.0)), under_y=103.0, over_y=88.0, join_at_end=True,
    )
    assert over.drawn[-2:] == ((40.0, 88.0), (0.0, 88.0))


def test_backing_opacity_switch(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("GBDRAW_LABEL_BACKING_OPACITY", raising=False)
    assert feature_label_backing_opacity() == 0.7
    monkeypatch.setenv("GBDRAW_LABEL_BACKING_OPACITY", "0.8")
    assert feature_label_backing_opacity() == 0.8


def test_a_study_switch_reads_the_mode_override_first(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GBDRAW_LABEL_UNDERLINE", "draw")
    monkeypatch.setenv("GBDRAW_LABEL_UNDERLINE_LINEAR", " Box ")
    assert study_switch("GBDRAW_LABEL_UNDERLINE", "linear") == "box"
    assert study_switch("GBDRAW_LABEL_UNDERLINE", "circular") == "draw"
    assert label_underline_in_box("linear") and not label_underline_counted("linear")
    assert label_underline_on("circular") and not label_underline_in_box("circular")
    monkeypatch.setenv("GBDRAW_LABEL_UNDERLINE_CIRCULAR", "box")
    assert not label_underline_in_box("circular")  # Circular draws a box-mode underline only


def test_the_underline_reach_with_a_central_baseline_is_half_the_text_plus_the_gap() -> None:
    metrics = get_font_vertical_metrics("Arial")
    reach = underline_reach_em("Arial", "central")
    assert metrics is not None and reach is not None
    assert reach == pytest.approx((metrics.ascent + metrics.descent) / 2.0 + UNDERLINE_GAP_EM)


def test_a_leader_from_above_gets_an_overline_a_gap_stop_or_its_old_end() -> None:
    def above(overline: OverlineMode, over_y: float) -> FeatureLabelInk:
        return FeatureLabelInk(
            ink=Aabb(0.0, 90.0, 40.0, 100.0), font_size=10.0, rotation_deg=0.0, center=(0.0, 95.0),
            leader=((20.0, 50.0), (20.0, 89.5)), under_y=103.0, over_y=over_y, overline=overline,
        )

    assert above("line", 88.0).underline == ((20.0, 88.0), (0.0, 88.0), (40.0, 88.0))
    # No line: the leader stops in the gap above the letters (over_y is set GAP_STOP_EM above them).
    gap = above("gap", 90.0 - GAP_STOP_EM * 10.0)
    assert gap.underline == ((20.0, 88.0), (20.0, 88.0), (20.0, 88.0))
    assert gap.drawn == ((20.0, 50.0), (20.0, 88.0))
    # No line and today's end.
    assert above("none", 88.0).underline is None
    assert above("none", 88.0).drawn == ((20.0, 50.0), (20.0, 89.5))


def test_a_halo_is_a_white_use_of_the_label_text(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GBDRAW_LABEL_HALO_LINEAR", "0.12")
    monkeypatch.setenv("GBDRAW_LABEL_HALO", "x")
    assert feature_label_halo_em("linear") == pytest.approx(0.12)
    assert feature_label_halo_em("circular") is None
    label = Text("tail protein", x=[0], y=[0], font_size=24.0)
    label.attribs["data-label-feature-id"] = "f1"
    other = Text("no font size")
    other.attribs["data-label-feature-id"] = "f2"
    twin = Text("tail protein", x=[0], y=[30], font_size=24.0)
    twin.attribs["data-label-feature-id"] = "f1"
    halos = label_halos([label, other, twin], 0.12)
    assert [label.attribs["id"], twin.attribs["id"]] == ["label_text_f1", "label_text_f1_1"]
    assert "id" not in other.attribs
    assert [h.attribs["xlink:href"] for h in halos] == ["#label_text_f1", "#label_text_f1_1"]
    assert halos[0].attribs["stroke"] == "#ffffff" and float(halos[0].attribs["stroke-width"]) == pytest.approx(2.88)
    assert "data-label-feature-id" not in halos[0].attribs  # the Web label editor binds labels by this attribute


def test_the_lp_d98_width_rule_and_the_gap_to_the_line_edge(monkeypatch: pytest.MonkeyPatch) -> None:
    assert rule_leader_width(10.0, "linear") is None and line_half_width_em(10.0, "linear") == 0.0
    monkeypatch.setenv("GBDRAW_LABEL_LEADER_WIDTH", "rule")
    assert rule_leader_width(10.0, "circular") == pytest.approx(1.0)  # the 1.0 px floor
    assert rule_leader_width(24.0, "linear") == pytest.approx(1.68)
    assert line_half_width_em(24.0, "linear") == 0.0  # gap still measured to the line's centre
    monkeypatch.setenv("GBDRAW_LABEL_LINE_GAP_LINEAR", "edge")
    assert line_half_width_em(24.0, "linear") == pytest.approx(0.035)
    assert line_half_width_em(24.0, "circular") == 0.0
