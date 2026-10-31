"""Record length text and source-coordinate tick positions for record displays."""
from typing import NamedTuple

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
