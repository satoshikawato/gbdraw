from __future__ import annotations

import pandas as pd
import pytest
from svgwrite.container import Group
from svgwrite.path import Path

from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.render.drawers.circular.numeric_series import circular_series_ink, draw_circular_series
from gbdraw.render.drawers.linear.numeric_series import draw_linear_series
from gbdraw.render.numeric_series import NumericSeriesPlot, gc_content_series, gc_skew_series
from gbdraw.svg.circular_tracks import CircularTrackPlace
from gbdraw.svg.linear_tracks import LinearTrackFrame
from tests.utils.series_configs import gc_content_config, gc_skew_config

ALONE = RecordDisplayTransform(120, 1, 1, None, True)


def _draw_circular(plot: NumericSeriesPlot, group_id: str = "gc_track") -> Group:
    return draw_circular_series(Group(), plot, CircularTrackPlace(100.0, 120, 20.0, 0.8, ALONE), group_id, None)


def _draw_linear(plot: NumericSeriesPlot, group_id: str = "gc_track") -> Group:
    return draw_linear_series(Group(), plot, LinearTrackFrame(0, 0, 120, 600.0, 24.0), group_id)


def _gc_content_config(*, stroke_color: str, stroke_width: float):
    return gc_content_config(
        fill_color="#10b981", stroke_color=stroke_color, stroke_width=stroke_width, fill_opacity=0.75, mode="deviation"
    )


def _gc_skew_config(*, stroke_color: str, stroke_width: float):
    return gc_skew_config(
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
    plot = gc_content_series(_gc_content_config(stroke_color=stroke_color, stroke_width=stroke_width), frame)
    circular_group, linear_group = _draw_circular(plot), _draw_linear(plot)

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
    plot = gc_skew_series(_gc_skew_config(stroke_color=stroke_color, stroke_width=stroke_width), frame)
    circular_group, linear_group = _draw_circular(plot), _draw_linear(plot)

    for group in (circular_group, linear_group):
        paths = _main_paths(group, stroke_color=stroke_color)
        assert len(paths) == 2
        assert {path.attribs["stroke-width"] for path in paths} == {stroke_width}


@pytest.mark.parametrize("draw", [_draw_circular, _draw_linear])
def test_none_stroke_remains_invisible_with_positive_width(draw) -> None:
    group = draw(gc_content_series(_gc_content_config(stroke_color="none", stroke_width=3.0), _gc_frame()))

    paths = _main_paths(group, stroke_color="none")
    assert paths
    assert {path.attribs["stroke-width"] for path in paths} == {3.0}


def test_a_record_on_a_shared_circle_inks_its_deviation_baseline_along_its_own_arc() -> None:
    # The ink model of GC content (deviation) and GC skew matches the drawn baseline, which closes along the
    # record's arc on a circle several records share, not round the whole circle.
    from gbdraw.layout.circular_ink import CircularInkRing
    from gbdraw.layout.record_coordinates import CirclePlace

    frame = _gc_frame()
    shared = RecordDisplayTransform(120, 1, 1, None, False).with_circle(CirclePlace(offset=30, length=300, scale_length=120), 120)
    alone = RecordDisplayTransform(120, 1, 1, None, True)
    for plot in (
        gc_content_series(_gc_content_config(stroke_color="black", stroke_width=1.0), frame),
        gc_skew_series(_gc_skew_config(stroke_color="black", stroke_width=1.0), frame),
    ):
        for transform, turns in ((shared, (0.1, 0.5)), (alone, (0.0, 1.0))):
            place = CircularTrackPlace(100.0, transform.length if transform.circle is None else 300, 20.0, 0.8, transform)
            rings = [ink for ink in circular_series_ink(plot, place) if isinstance(ink, CircularInkRing)]
            assert len(rings) == 1
            assert (rings[0].start_turn, rings[0].end_turn) == pytest.approx(turns)
