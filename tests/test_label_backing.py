"""The shared label backing (LP-D82 region labels, LP-D85 study prototype for feature labels)."""

from __future__ import annotations

import pytest

from gbdraw.layout.spatial import Aabb
from gbdraw.render.drawers.label_backing import (
    BACKING_PADDING_EM,
    FeatureLabelInk,
    backed_feature_labels,
    backing_box,
    feature_label_backing_mode,
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
    assert label.backing.max_y == pytest.approx(100.5)


@pytest.mark.parametrize(("value", "mode"), [("", "off"), ("through", "through"), ("ALL", "all"), ("x", "off")])
def test_prototype_switch(monkeypatch: pytest.MonkeyPatch, value: str, mode: str) -> None:
    monkeypatch.setenv("GBDRAW_LABEL_BACKING", value)
    assert feature_label_backing_mode() == mode
