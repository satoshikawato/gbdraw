"""Record length text and source-coordinate tick positions for record displays."""
from .record_coordinates import RecordDisplayTransform


def record_length_label(length: int) -> str:
    """Return the definition length line, the same for every display start and mode."""
    return f"{length:,} bp"


def source_ruler_ticks(transform: RecordDisplayTransform, interval: int) -> list[tuple[int, int]]:
    return [(base, transform.source_base_to_display_index(base))
            for base in range(interval, transform.length, interval)]
