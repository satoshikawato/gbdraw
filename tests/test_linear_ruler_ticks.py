"""Linear ruler labels at the Auto tick interval do not overlap at the default width (OV-2747)."""

from __future__ import annotations

from dataclasses import replace

import pytest

from gbdraw.api.config import apply_config_overrides
from gbdraw.auto_sizes import tick_and_window_tiers
from gbdraw.canvas import LinearCanvasConfigurator
from gbdraw.config.models import LinearRenderProfile
from gbdraw.core.text import drawn_text_box
from gbdraw.render.groups.linear.length_bar import LengthBarGroup

_CFG = apply_config_overrides(None, {"objects.scale.style": "ruler"})
_TIERS = tick_and_window_tiers(_CFG)["linearTicks"]


def _ruler(length: int, interval: int | None = None) -> LengthBarGroup:
    """The ruler scale of a default Linear figure whose longest record is ``length`` bp, at the Auto interval or a
    forced one."""

    canvas = LinearCanvasConfigurator(
        num_of_entries=1, longest_genome=length, profile=LinearRenderProfile(_CFG), legend="none"
    )
    return LengthBarGroup(
        float(canvas.alignment_width),
        length,
        canvas.length_param,
        scale=replace(_CFG.objects.scale, interval=interval),
        font_family=_CFG.objects.text.font_family,
    )


def _smallest_label_gap_px(ruler: LengthBarGroup) -> float:
    """The smallest gap between the drawn boxes of two neighbouring ruler labels (negative: they overlap)."""

    boxes = []
    for element in ruler.scale_group.elements:
        if element.elementname == "text":
            width, _ = drawn_text_box(element.text, ruler.font_family, ruler.ruler_label_font_size, ruler.font_weight)
            x = float(element.attribs["x"])
            boxes.append((x - width / 2.0, x + width / 2.0))
    boxes.sort()
    assert len(boxes) > 1
    return min(right[0] - left[1] for left, right in zip(boxes, boxes[1:]))


def test_the_ruler_labels_of_a_113_mbp_record_do_not_overlap() -> None:
    assert _smallest_label_gap_px(_ruler(112_934_049)) >= 0.0


@pytest.mark.parametrize(("below", "interval"), _TIERS["tiers"])
def test_the_ruler_labels_of_every_tier_fit_at_its_longest_record(below: int, interval: int) -> None:
    ruler = _ruler(below - 1)

    assert _smallest_label_gap_px(ruler) >= 0.0


@pytest.mark.parametrize(
    ("below", "interval"), [(below, interval) for below, interval in _TIERS["tiers"] if below > 5_000_000]
)
def test_a_tier_above_5_mbp_ends_where_its_ruler_labels_start_to_overlap(below: int, interval: int) -> None:
    """The rows above 5 Mbp are measured, not chosen: each ends where its ruler labels stop fitting, so a figure
    whose labels fit keeps its interval. The axis ruler spaces its labels over one base less and fits a few bases
    longer, so the ruler scale sets the bound."""

    assert _smallest_label_gap_px(_ruler(below, interval)) < 0.0


def test_the_ruler_above_the_last_tier_fits_records_of_10_gbp() -> None:
    assert _smallest_label_gap_px(_ruler(10_000_000_000)) >= 0.0
