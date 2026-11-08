"""Record length text and source-coordinate tick positions for record displays."""
from typing import NamedTuple

from Bio.SeqRecord import SeqRecord

from gbdraw.core.record_metadata import _read_coord_map, _source_period, _wrap_source_base

from .record_coordinates import RecordDisplayTransform


class RulerTick(NamedTuple):
    """A source base on the ruler: its row and its base index in that row."""

    base: int
    row: int
    offset: int


def record_length_label(length: int) -> str:
    """Return the definition length line, the same for every display start and mode."""
    return f"{length:,} bp"


def source_ruler_ticks(transform: RecordDisplayTransform, interval: int) -> list[RulerTick]:
    """The multiples of ``interval`` on the display axis, each on its row.

    A rotated axis shows every base and keeps every multiple. An affine axis (a crop shows only
    part of the source) keeps those strictly inside it, as its coordinate-map ruler does.
    """
    ticks = []
    for base in range(interval, transform.length, interval):
        index = transform.source_base_to_display_index(base)
        if transform.start_coordinate is None and not 0 < index < transform.drawn_length - 1:
            continue
        ticks.append(RulerTick(base, *transform.display_index_row(index)))
    return ticks


class MapTick(NamedTuple):
    """A ruler tick of a record drawn by its coordinate map: where it lies along the axis (0 at the axis start,
    1 at its end) and the base its label names."""

    position: float
    base: int


def record_map_span(record: SeqRecord, axis_length: int) -> tuple[int, int]:
    """The coordinates at the two ends of a record's axis, read through its coordinate map.

    Past the source end of a crop across the origin they run on (``source_period`` wraps them for labels).
    """

    base, step = _read_coord_map(record)
    return base, base + step * max(0, axis_length - 1)


def coordinate_map_ticks(
    record: SeqRecord, axis_length: int, interval: int, *, record_local: bool = False
) -> tuple[MapTick, ...]:
    """The ruler ticks of a record drawn by its coordinate map: the multiples of ``interval`` strictly inside it.

    A record-local ruler counts from 0 along the axis. The axis of a crop across
    the origin (D-28) counts on past the source end; its multiples restart at
    the origin, and each label names the source base. Ticks run along the axis.
    """

    start, end = (0, axis_length) if record_local else record_map_span(record, axis_length)
    period = None if record_local else _source_period(record)
    span = abs(end - start)
    if span <= 0 or interval <= 0:
        return ()
    low, high = min(start, end), max(start, end)
    # Each turn of the circle starts its multiples again; without a period the axis is one turn from 0.
    offsets = [0] if period is None else range(((low - 1) // period) * period, high, period)
    coordinates: list[int] = []
    for offset in offsets:
        multiple = max(1, (low - offset) // interval) * interval
        while offset + multiple < high and (period is None or multiple <= period):
            if offset + multiple > low:
                coordinates.append(offset + multiple)
            multiple += interval
    coordinates.sort(reverse=start > end)
    return tuple(MapTick(abs(value - start) / span, _wrap_source_base(value, period)) for value in coordinates)
