"""The record band of a circle that several records share: one arc per record, neighbours told apart by their
fills, each named on its arc by the first of its names that fits."""

from __future__ import annotations

import math
from typing import NamedTuple

from gbdraw.svg.elements import Group, Path, Text

from ....core.text import drawn_text_box
from ....layout.text_geometry import reads_upside_down
from ....svg.circular_arcs import circle_arc, circle_arc_to, circle_point
from ....svg.path_data import path_number, path_point

# The spacer slot that keeps the band's room outside the tracks, and the band's width in name font sizes.
RECORD_BAND_SLOT_ID = "record_band"
RECORD_BAND_WIDTH_EM = 1.4
RECORD_BAND_FILLS = ("#9ca3af", "#d1d5db")
# The last arc of an odd count, which meets the first (the first tone) across 12 o'clock and the one before it (the
# second tone).
RECORD_BAND_LAST_FILL = "#b7bcc5"
RECORD_BAND_TEXT_COLOR = "#1f2937"
# The least room between the names of two neighbouring arcs, in pixels: a gap narrower than this keeps the rest
# inside the arcs (``record_band_name_end_pad_px``).
RECORD_BAND_NAME_CLEARANCE_PX = 2.0


class RecordBandArc(NamedTuple):
    """One record's arc on the band: its identity, the names it can show, longest first (``record_band_names``),
    and where it lies (start and span in degrees, the SVG frame)."""

    record_key: str
    record_id: str
    names: tuple[str, ...]
    start_deg: float
    span_deg: float


def _sector(inner_px: float, outer_px: float, start_deg: float, span_deg: float) -> str:
    return (
        f"{circle_arc(outer_px, start_deg, span_deg)} "
        f"L {path_point(*circle_point(inner_px, start_deg + span_deg))} "
        f"{circle_arc_to(inner_px, start_deg + span_deg, -span_deg)} Z"
    )


def _arc_px(span_deg: float, radius_px: float) -> float:
    return 2.0 * math.pi * radius_px * span_deg / 360.0


def _gap_deg(earlier: RecordBandArc, later: RecordBandArc) -> float:
    """The gap from the end of ``earlier`` to the start of ``later``, clockwise, in degrees."""

    gap = (later.start_deg - (earlier.start_deg + earlier.span_deg)) % 360.0
    return gap if gap <= 360.0 - (earlier.span_deg + later.span_deg) else 0.0


def record_band_name_end_pad_px(gap_deg: float, radius_px: float) -> float:
    """The room a name keeps from an end of its arc where the gap there is ``gap_deg`` wide: what the gap lacks of
    ``RECORD_BAND_NAME_CLEARANCE_PX`` between two neighbours' names, half on each side; none at a gap that wide."""

    return max(0.0, 0.5 * (RECORD_BAND_NAME_CLEARANCE_PX - _arc_px(gap_deg, radius_px)))


def record_band_name_fits(
    name: str, span_deg: float, *, radius_px: float, font_family: str, font_size: float, end_pad_px: float
) -> bool:
    """Whether ``name`` fits on an arc of ``span_deg`` at ``radius_px``: the box it is drawn in
    (``drawn_text_box``), centred on the arc, lies within it and ``end_pad_px`` from each end."""

    width, _height = drawn_text_box(name, font_family, font_size)
    return float(width) + 2.0 * end_pad_px <= _arc_px(span_deg, radius_px)


def record_band_fill(index: int, count: int) -> str:
    """The fill of arc ``index`` of the band's ``count`` arcs: the two tones in turn, and a third for the last arc
    of an odd count, so that any two neighbours differ, the last and the first across 12 o'clock included."""

    if count > 1 and count % 2 == 1 and index == count - 1:
        return RECORD_BAND_LAST_FILL
    return RECORD_BAND_FILLS[index % len(RECORD_BAND_FILLS)]


def record_band_group(
    arcs: tuple[RecordBandArc, ...],
    *,
    inner_px: float,
    outer_px: float,
    font_family: str,
    font_size: float,
) -> Group:
    """The band group: each record's sector, filled by ``record_band_fill``, and on it the first of its names that
    fits, upright and along the circle, apart from its neighbours' names across a narrow gap; no name only when
    none fits."""

    group = Group(id="record_band_arcs")
    middle = 0.5 * (inner_px + outer_px)
    for index, arc in enumerate(arcs):
        # The gaps to the neighbours on either side. Two arcs leave at most a turn less both their spans between
        # them, so a reading past that is arcs that meet to within rounding: a gap of 0. An arc alone has no
        # neighbour's name to keep apart from.
        pairs = ((arcs[index - 1], arc), (arc, arcs[(index + 1) % len(arcs)])) if len(arcs) > 1 else ()
        end_pad = max((
            record_band_name_end_pad_px(_gap_deg(earlier, later), middle) for earlier, later in pairs
        ), default=0.0)
        item = Group(id=f"record_band_{index + 1}")
        item.attribs["data-gbdraw-record-id"] = arc.record_id
        item.attribs["data-gbdraw-record-key"] = arc.record_key
        item.add(Path(d=_sector(inner_px, outer_px, arc.start_deg, arc.span_deg),
                      fill=record_band_fill(index, len(arcs)), stroke="none"))
        name = next((
            name for name in arc.names
            if record_band_name_fits(
                name, arc.span_deg, radius_px=middle, font_family=font_family, font_size=font_size, end_pad_px=end_pad
            )
        ), None)
        if name is not None:
            angle = arc.start_deg + 0.5 * arc.span_deg
            x, y = circle_point(middle, angle)
            # Tangent to the circle, turned to read upright on the lower half.
            rotation = angle + 90.0
            if reads_upside_down(angle):
                rotation -= 180.0
            item.add(Text(
                name,
                insert=(round(x, 3), round(y, 3)),
                transform=f"rotate({path_number(rotation)},{path_number(x)},{path_number(y)})",
                font_family=font_family,
                font_size=font_size,
                fill=RECORD_BAND_TEXT_COLOR,
                text_anchor="middle",
                dominant_baseline="central",
            ))
        group.add(item)
    return group


__all__ = [
    "RECORD_BAND_FILLS",
    "RECORD_BAND_LAST_FILL",
    "RECORD_BAND_NAME_CLEARANCE_PX",
    "RECORD_BAND_SLOT_ID",
    "RECORD_BAND_WIDTH_EM",
    "RecordBandArc",
    "record_band_fill",
    "record_band_group",
    "record_band_name_end_pad_px",
    "record_band_name_fits",
]
