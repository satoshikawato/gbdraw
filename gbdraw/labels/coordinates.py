"""Label anchors selected from the feature owner's projected coverage."""
from __future__ import annotations

from dataclasses import dataclass

from ..features.objects import FeatureObject
from ..layout.record_coordinates import RecordDisplayTransform


@dataclass(frozen=True)
class DisplayLabelSegment:
    start: float
    end: float
    middle: float
    span: float
    strand: str
    source_span: tuple[int, int]
    part_index: int
    source_middle: float | None


def display_label_segment(
    feature: FeatureObject, length: int, *, circular: bool = False,
    record_transform: RecordDisplayTransform | None = None,
) -> DisplayLabelSegment | None:
    """Choose one anchor, retaining biological parts for circular placement.

    Linear placement uses the longest continuous fragment, breaking ties by
    display start and original part order. Circular placement retains the
    source part across its artificial cut and the existing origin-part join.
    No source or local feature coordinates are rewritten here.
    """
    blocks = [part for part in feature.display_parts or () if part.kind == "block"]
    if not blocks:
        return None
    groups = [[part] for part in blocks]
    if circular:
        by_part = {}
        for part in blocks:
            by_part.setdefault(part.fragment.part_index, []).append(part)
        groups = list(by_part.values())
        if len(groups) == 2 and len({part.strand for part in blocks}) == 1:
            spans = [(min(p.fragment.local_start for p in group),
                      max(p.fragment.local_end for p in group)) for group in groups]
            if min(start for start, _ in spans) == 0 and max(end for _, end in spans) == length:
                groups = [groups[0] + groups[1]]

    group = min(groups, key=lambda items: (
        -sum(p.fragment.display_end - p.fragment.display_start for p in items),
        (0 if circular else min(p.fragment.display_start for p in items)),
        items[0].fragment.part_index,
    ))
    span = sum(p.fragment.display_end - p.fragment.display_start for p in group)
    remaining = span / 2.0
    anchor = group[-1].fragment
    for part in group:
        anchor = part.fragment
        size = anchor.display_end - anchor.display_start
        if remaining <= size:
            break
        remaining -= size
    middle = (anchor.display_start + remaining if anchor.orientation == 1
              else anchor.display_end - remaining)
    local_forward = anchor.orientation == 1
    ordered = group if local_forward else list(reversed(group))
    source_middle = None
    if record_transform is not None:
        distance = middle - anchor.display_start
        source_middle = (anchor.source_start + distance if record_transform.source_step == 1
                         else anchor.source_end - distance)
        if circular and middle == length:
            middle = record_transform.source_boundary_to_display_offset(int(source_middle))
    return DisplayLabelSegment(
        ordered[0].fragment.display_start, ordered[-1].fragment.display_end,
        middle, span, group[0].strand,
        (anchor.source_start, anchor.source_end), anchor.part_index, source_middle,
    )
