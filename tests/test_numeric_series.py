"""One series value per numeric track: the per-kind differences the drawers read are data on the value."""

from __future__ import annotations

import pandas as pd
import pytest

from gbdraw.render.numeric_series import depth_series, gc_content_series, gc_skew_series
from tests.utils.series_configs import depth_config, gc_content_config, gc_skew_config

GC_FRAME = pd.DataFrame({"GC content": [0.4, 0.6, 0.5], "GC skew": [0.1, -0.1, 0.0]}, index=[0, 40, 80])
DEPTH_FRAME = pd.DataFrame({"position": [0, 40, 80], "depth": [5.0, 10.0, 20.0]})


def _paint(**extra: object) -> dict[str, object]:
    return {"stroke_color": "none", "stroke_width": 0.0, "fill_opacity": 1.0, **extra}


def _plots():
    content = gc_content_config(**_paint(fill_color="#000", mode="deviation"))
    skew = gc_skew_config(**_paint(high_fill_color="#f00", low_fill_color="#00f"))
    depth = depth_config(**_paint(fill_color="#000", min_depth=None, max_depth=None, normalize=False, show_axis=False))
    return {
        "gc_content": gc_content_series(content, GC_FRAME),
        "gc_skew": gc_skew_series(skew, GC_FRAME),
        "depth": depth_series(depth, DEPTH_FRAME),
    }


@pytest.mark.parametrize("kind", ["gc_content", "gc_skew", "depth"])
def test_a_stroke_that_paints_nothing_writes_no_width_in_any_kind(kind: str) -> None:
    assert not _plots()[kind].paint.writes_stroke_width


def test_the_kinds_keep_their_line_ends_and_fills() -> None:
    plots = _plots()
    closes = {kind: plot.samples.closes_at_record_end for kind, plot in plots.items() if plot.samples is not None}
    assert closes == {"gc_content": True, "gc_skew": False, "depth": False}
    # Only skew fills below its baseline, with the clip named by its dinucleotide.
    assert {kind: plot.paint.below_fill for kind, plot in plots.items()} == {
        "gc_content": None, "gc_skew": "#00f", "depth": None,
    }
    assert plots["gc_skew"].clip_key == "GC"
