from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
import pytest
from svgwrite.container import Group
from svgwrite.path import Path

from gbdraw.render.drawers.circular.gc_content import (
    GcContentDrawer as CircularGcContentDrawer,
)
from gbdraw.render.drawers.circular.gc_skew import SkewDrawer as CircularSkewDrawer
from gbdraw.render.drawers.linear.gc_content import (
    GcContentDrawer as LinearGcContentDrawer,
)
from gbdraw.render.drawers.linear.gc_skew import SkewDrawer as LinearSkewDrawer


def _gc_content_config(*, stroke_color: str, stroke_width: float) -> SimpleNamespace:
    return SimpleNamespace(
        high_fill_color="#10b981",
        low_fill_color="#10b981",
        fill_color="#10b981",
        stroke_color=stroke_color,
        stroke_width=stroke_width,
        fill_opacity=0.75,
        mode="deviation",
    )


def _gc_skew_config(*, stroke_color: str, stroke_width: float) -> SimpleNamespace:
    return SimpleNamespace(
        high_fill_color="#ef4444",
        low_fill_color="#3b82f6",
        stroke_color=stroke_color,
        stroke_width=stroke_width,
        fill_opacity=0.75,
    )


def _gc_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "GC content": [0.35, 0.65, 0.45],
            "GC skew": [0.25, -0.25, 0.1],
        },
        index=[0, 40, 80],
    )


def _main_paths(group: Group, *, stroke_color: str) -> list[Path]:
    return [
        element
        for element in group.elements
        if isinstance(element, Path) and element.attribs.get("stroke") == stroke_color
    ]


@pytest.mark.parametrize("stroke_width", [0.0, 2.5])
def test_gc_content_outline_contract_matches_in_both_modes(
    stroke_width: float,
) -> None:
    stroke_color = "#123456"
    frame = _gc_frame()
    circular_group = CircularGcContentDrawer(
        _gc_content_config(stroke_color=stroke_color, stroke_width=stroke_width)
    ).draw(100.0, Group(), frame, 120, 20.0, 0.8, "GC")
    linear_group = LinearGcContentDrawer(
        _gc_content_config(stroke_color=stroke_color, stroke_width=stroke_width)
    ).draw(Group(), frame, 120, 600.0, 1.0, 24.0, 0.0, 0.0, "GC")

    for group in (circular_group, linear_group):
        paths = _main_paths(group, stroke_color=stroke_color)
        assert len(paths) == 1
        assert paths[0].attribs["stroke-width"] == stroke_width


@pytest.mark.parametrize("stroke_width", [0.0, 2.5])
def test_gc_skew_outline_contract_matches_in_both_modes(
    stroke_width: float,
) -> None:
    stroke_color = "#654321"
    frame = _gc_frame()
    circular_group = CircularSkewDrawer(
        _gc_skew_config(stroke_color=stroke_color, stroke_width=stroke_width)
    ).draw(100.0, Group(), frame, 120, 20.0, 0.8, "GC")
    linear_group = LinearSkewDrawer(
        _gc_skew_config(stroke_color=stroke_color, stroke_width=stroke_width)
    ).draw(
        Group(),
        frame,
        120,
        600.0,
        1.0,
        24.0,
        0.0,
        0.0,
        "GC",
    )

    for group in (circular_group, linear_group):
        paths = _main_paths(group, stroke_color=stroke_color)
        assert len(paths) == 2
        assert {path.attribs["stroke-width"] for path in paths} == {stroke_width}


@pytest.mark.parametrize(
    ("drawer", "args"),
    [
        (CircularGcContentDrawer, (100.0, Group(), _gc_frame(), 120, 20.0, 0.8, "GC")),
        (LinearGcContentDrawer, (Group(), _gc_frame(), 120, 600.0, 1.0, 24.0, 0.0, 0.0, "GC")),
    ],
)
def test_none_stroke_remains_invisible_with_positive_width(drawer, args) -> None:
    group = drawer(
        _gc_content_config(stroke_color="none", stroke_width=3.0)
    ).draw(*args)

    paths = _main_paths(group, stroke_color="none")
    assert paths
    assert {path.attribs["stroke-width"] for path in paths} == {3.0}
