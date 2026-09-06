"""Source-coordinate text and tick positions for explicit record displays."""
from .record_coordinates import RecordDisplayTransform


def display_coordinate_label(transform: RecordDisplayTransform) -> str:
    start = transform.display_index_to_source_base(0)
    length = transform.length
    if transform.source_step == 1:
        spans = [(start, length)]
        if start > 1:
            spans.append((1, start - 1))
    else:
        spans = [(start, 1)]
        if start < length:
            spans.append((length, start + 1))
    return ", ".join(f"[{first:,}..{last:,}]" for first, last in spans) + " bp"


def source_ruler_ticks(transform: RecordDisplayTransform, interval: int) -> list[tuple[int, int]]:
    return [(base, transform.source_base_to_display_index(base))
            for base in range(interval, transform.length, interval)]
