"""Linear ruler labels at the Auto tick interval do not overlap at the default width (OV-2747)."""

from __future__ import annotations

from dataclasses import replace

import pytest

from gbdraw.api.config import apply_config_overrides
from gbdraw.auto_sizes import tick_and_window_tiers
from gbdraw.canvas import LinearCanvasConfigurator
from gbdraw.config.models import LinearRenderProfile
from gbdraw.core.text import drawn_text_box
from gbdraw.labels.label_style import ROW_LABEL_GAP_EM
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


def _label_gap_shortfall_px(ruler: LengthBarGroup) -> float:
    """How far the closest two neighbouring ruler labels, the end label included, stay from a quarter em apart
    (OV-2748): the gap between their drawn boxes minus ``ROW_LABEL_GAP_EM`` of the label font. Negative: they
    are closer, so they read as one."""

    boxes = []
    for element in ruler.scale_group.elements:
        if element.elementname == "text":
            width, _ = drawn_text_box(element.text, ruler.font_family, ruler.ruler_label_font_size, ruler.font_weight)
            x = float(element.attribs["x"])
            boxes.append((x - width / 2.0, x + width / 2.0))
    boxes.sort()
    assert len(boxes) > 1
    gap_px = min(right[0] - left[1] for left, right in zip(boxes, boxes[1:]))
    return gap_px - ROW_LABEL_GAP_EM * float(ruler.ruler_label_font_size)


def test_the_ruler_labels_of_a_113_mbp_record_keep_a_quarter_em_apart() -> None:
    assert _label_gap_shortfall_px(_ruler(112_934_049)) >= 0.0


@pytest.mark.parametrize(("below", "interval"), _TIERS["tiers"])
def test_the_ruler_labels_of_every_tier_keep_a_quarter_em_apart_at_its_longest_record(
    below: int, interval: int
) -> None:
    assert _label_gap_shortfall_px(_ruler(below - 1)) >= 0.0


@pytest.mark.parametrize(
    ("below", "interval"), [(below, interval) for below, interval in _TIERS["tiers"] if below > 5_000_000]
)
def test_a_tier_above_5_mbp_ends_where_its_ruler_labels_come_closer_than_a_quarter_em(
    below: int, interval: int
) -> None:
    """The rows above 5 Mbp are measured, not chosen: each ends where its ruler labels would come closer than a
    quarter em, so only figures whose labels would read as one take a coarser interval. The axis ruler spaces its
    labels over one base less and keeps the gap a few bases longer, so the ruler scale sets the bound."""

    assert _ruler(below - 1, interval).label_stride == 1
    assert _ruler(below, interval).label_stride > 1


def test_the_ruler_above_the_last_tier_fits_records_of_10_gbp() -> None:
    assert _label_gap_shortfall_px(_ruler(10_000_000_000)) >= 0.0


def test_the_ruler_drops_an_end_label_closer_than_a_quarter_em() -> None:
    """At 4,606,000 bp the end label "4.6 Mbp" stood 2.7 px after "4.5 Mbp" and read as one with it; the end
    label now needs the quarter em the other labels keep (OV-2748)."""

    ruler = _ruler(4_606_000)

    assert [element.text for element in ruler.scale_group.elements if element.elementname == "text"][-1] == "4.5 Mbp"
    assert _label_gap_shortfall_px(ruler) >= 0.0


def _texts(ruler: LengthBarGroup) -> list[str]:
    return [element.text for element in ruler.scale_group.elements if element.elementname == "text"]


def _tick_lines(ruler: LengthBarGroup) -> int:
    return sum(1 for element in ruler.scale_group.elements if element.elementname == "line") - 1  # less the axis


def test_a_typed_interval_whose_labels_would_print_into_each_other_labels_every_second_tick() -> None:
    """OV-3051: lambda (48,502 bp) at a typed 1 kbp interval printed "1 kbp2 kbp3 kbp" along the ruler; every tick
    keeps its mark, and the labels step 0, 2, 4, ... kbp so they keep a quarter em apart (as Circular, OV-3020)."""

    ruler = _ruler(48_502, 1_000)

    assert ruler.label_stride == 2
    assert _tick_lines(ruler) == 49  # 0 to 48 kbp; the 48,502 bp end tick is too close to 48 kbp to draw
    assert _texts(ruler)[:4] == ["0", "2 kbp", "4 kbp", "6 kbp"]
    assert _label_gap_shortfall_px(ruler) >= 0.0


def test_ruler_labels_that_keep_apart_are_all_drawn() -> None:
    ruler = _ruler(48_502)

    assert ruler.label_stride == 1
    assert _texts(ruler)[:10] == ["0", *(f"{kbp} kbp" for kbp in range(5, 50, 5))]
