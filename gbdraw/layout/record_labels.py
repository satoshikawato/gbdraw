"""Source-coordinate tick positions for explicit record displays."""
from .record_coordinates import RecordDisplayTransform


def source_ruler_ticks(transform: RecordDisplayTransform, interval: int) -> list[tuple[int, int]]:
    return [(base, transform.source_base_to_display_index(base))
            for base in range(interval, transform.length, interval)]
