"""The drawn axis of a Circular scalar track (depth, GC content in percent): one owner for both drawers.

The axis stands up at the top of the circle from the band's baseline to its outer edge, with its ticks to the right.
A numbered tick label sits on the band, so the track's own area (the depth fill, the GC content area and band) runs
under it wherever the data reach the tick; each one gets a halo (``label_halo``: a white outline of its letters
drawn just before them), so the number stays readable on the data and the data stay visible around it. The range
label sits inside the band near its outer edge and is drawn as before.
"""

from __future__ import annotations

from collections.abc import Sequence

from gbdraw.svg.elements import Circle, Group, Line, Path, Text

from ....layout.circular_depth_axis import SCALAR_AXIS_LABEL_GAP_PX, CircularRangeLabel, ScalarAxisTick
from ....svg.circular_ticks import text_along_circle
from ....svg.ids import stable_svg_id
from ..label_backing import label_halo


def circular_scalar_axis_group(
    axis_group_id: str,
    baseline_radius: float,
    outer_radius: float,
    ticks: Sequence[ScalarAxisTick],
    range_label: CircularRangeLabel | None,
    *,
    stroke: str,
    stroke_width: float,
    font_size: float,
    font_family: str,
    baseline_circle: bool = True,
) -> Group:
    """The axis group ``axis_group_id``: the baseline circle (left out with ``baseline_circle`` false, where each
    record of a shared circle draws its own baseline arc), the axis line, the ``ticks`` with their numbers (each on
    a halo) and the ``range_label``."""

    # The numbered labels take their colour from the group: a halo (a ``<use>`` with no fill) draws no letters of
    # its own only when the text it repeats has no fill of its own, so a CairoSVG PDF carries each number once.
    numbered = any(tick.label is not None for tick in ticks)
    axis_group = Group(id=axis_group_id, fill=stroke) if numbered else Group(id=axis_group_id)
    if baseline_circle:
        axis_group.add(Circle(center=(0, 0), r=baseline_radius, fill="none", stroke=stroke, stroke_width=stroke_width))
    axis_group.add(Line(start=(0, -baseline_radius), end=(0, -outer_radius), stroke=stroke, stroke_width=stroke_width))
    for tick in ticks:
        axis_group.add(
            Line(start=(0, -tick.radius_px), end=(tick.length_px, -tick.radius_px), stroke=stroke, stroke_width=stroke_width)
        )
        if tick.label is not None:
            text = Text(
                tick.label,
                insert=(tick.length_px + SCALAR_AXIS_LABEL_GAP_PX, -tick.radius_px),
                font_size=font_size,
                font_family=font_family,
                text_anchor="start",
                dominant_baseline="central",
            )
            axis_group.add(label_halo(text, stable_svg_id("scalar_axis_label", axis_group_id, tick.label, tick.radius_px)))
            axis_group.add(text)
    if range_label is not None:
        for element in circular_range_label_elements(
            range_label, axis_group_id, stroke=stroke, font_size=font_size, font_family=font_family
        ):
            axis_group.add(element)
    return axis_group


def circular_range_label_elements(
    range_label: CircularRangeLabel, axis_group_id: str, *, stroke: str, font_size: float, font_family: str
) -> tuple[Path, Text]:
    """The ``range_label`` of the axis ``axis_group_id``: its path and its text along the circle, upright
    (``CircularRangeLabel.text_path``). The axis group draws it, or, on a circle several records share, the record
    whose arc holds it."""

    path = range_label.text_path(font_family)
    return text_along_circle(
        range_label.text,
        path.radius_px,
        path.start_deg,
        path.span_deg,
        id_parts=(axis_group_id,),
        fill=stroke,
        font_size=font_size,
        font_family=font_family,
    )


__all__ = ["circular_range_label_elements", "circular_scalar_axis_group"]
