"""The sub-structure of a wrapped record: rows are id-less groups inside its one record and track groups."""

from __future__ import annotations

from collections.abc import Callable, Sequence

from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.scalar_axis import ScalarRow, ScalarRows
from gbdraw.svg.elements import ClipPath, Group
from svgwrite.base import BaseElement

WRAP_ROW_ATTRIBUTE = "data-gbdraw-wrap-row"


def wrap_row_group(row: int, shift_y: float) -> Group:
    """One row's group, ``shift_y`` below the frame of the group that holds it."""

    group = Group()
    group.attribs[WRAP_ROW_ATTRIBUTE] = str(row)
    if shift_y:
        group.translate(0.0, shift_y)
    return group


def _without_ids(element: BaseElement) -> None:
    """Drop the ids a track's drawer gave a later row: the first row keeps the track's ids. A clip path keeps
    its own, which its row refers to."""

    for child in getattr(element, "elements", ()):
        if isinstance(child, ClipPath):
            continue
        child.attribs.pop("id", None)
        _without_ids(child)


def draw_track_rows(
    group: Group,
    transform: RecordDisplayTransform | None,
    record_length: int,
    width: float,
    row_shifts: Sequence[float],
    draw: Callable[[Group, int, float, ScalarRow | None], Group],
) -> Group:
    """Draw a scalar track: once, or, on a wrapped record, once per row at the one scale of its rows.

    ``draw`` takes the group, the bases and width of what it draws, and the row (``None`` unwrapped).
    """

    if transform is None or transform.wrap_length is None:
        return draw(group, record_length, width, None)
    samples = ScalarRows(transform)
    for row, shift_y in enumerate(row_shifts):
        row_length = transform.row_length(row)
        target = wrap_row_group(row, shift_y)
        draw(target, row_length, width * row_length / transform.axis_length, ScalarRow(samples, row))
        if row:
            _without_ids(target)
        if target.elements:
            group.add(target)
    return group


__all__ = ["WRAP_ROW_ATTRIBUTE", "draw_track_rows", "wrap_row_group"]
