"""One series value per numeric track: the per-kind differences the drawers read are data on the value."""

from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
import pytest

from gbdraw.render.numeric_series import depth_series, gc_content_series, gc_skew_series

GC_FRAME = pd.DataFrame({"GC content": [0.4, 0.6, 0.5], "GC skew": [0.1, -0.1, 0.0]}, index=[0, 40, 80])
DEPTH_FRAME = pd.DataFrame({"position": [0, 40, 80], "depth": [5.0, 10.0, 20.0]})


def _paint(**extra: object) -> dict[str, object]:
    return {"stroke_color": "none", "stroke_width": 0.0, "fill_opacity": 1.0, "dinucleotide": "GC", **extra}


def _plots():
    content = SimpleNamespace(**_paint(fill_color="#000", mode="deviation"))
    skew = SimpleNamespace(**_paint(high_fill_color="#f00", low_fill_color="#00f"))
    depth = SimpleNamespace(**_paint(fill_color="#000", min_depth=None, max_depth=None, normalize=False, show_axis=False))
    return {
        "gc_content": gc_content_series(content, GC_FRAME),
        "gc_skew": gc_skew_series(skew, GC_FRAME),
        "depth": depth_series(depth, DEPTH_FRAME),
    }


@pytest.mark.parametrize(
    ("kind", "circular", "linear"),
    [("gc_content", False, False), ("gc_skew", False, True), ("depth", True, True)],
)
def test_a_stroke_that_paints_nothing_writes_its_width_as_each_kind_did(kind: str, circular: bool, linear: bool) -> None:
    paint = _plots()[kind].paint
    assert (paint.writes_stroke_width("circular"), paint.writes_stroke_width("linear")) == (circular, linear)


def test_the_kinds_keep_their_line_ends_fills_and_origins() -> None:
    plots = _plots()
    closes = {kind: plot.samples.closes_at_record_end for kind, plot in plots.items() if plot.samples is not None}
    assert closes == {"gc_content": True, "gc_skew": False, "depth": False}
    # Only skew fills below its baseline, with the clip named by its dinucleotide.
    assert {kind: plot.paint.below_fill for kind, plot in plots.items()} == {
        "gc_content": None, "gc_skew": "#00f", "depth": None,
    }
    assert plots["gc_skew"].clip_key == "GC"
    # svgwrite writes numbers as given: GC tracks place their Linear row at 0, depth at 0.0.
    origins = {kind: repr(plot.linear_origin) for kind, plot in plots.items()}
    assert origins == {"gc_content": "0", "gc_skew": "0", "depth": "0.0"}
