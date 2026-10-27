"""The order-preserving packer of Circular horizontal labels (CP-D04, CP-D05) and the counts it reads."""

from __future__ import annotations

from typing import cast

from gbdraw.labels.circular import _count_label_overlaps
from gbdraw.labels.circular_types import CircularLabel


def _box_label(start_x: float, start_y: float, width: float, **extra: float | bool) -> CircularLabel:
    """An outer label at 3 o'clock (middle baseline) whose text box runs from ``start_x`` to ``start_x + width``."""

    label = {
        "middle": 250.0,
        "start_x": start_x,
        "start_y": start_y,
        "width_px": width,
        "height_px": 10.0,
        "is_inner": False,
        **extra,
    }
    return cast(CircularLabel, label)


# ---- OV-730: the min-gap overlap count pads its broad phase by the largest pair gap ----------------------


def test_min_gap_overlap_count_finds_a_pair_wider_than_any_single_label_gap() -> None:
    """Labels 1 and 3 (not list neighbours) are 8 px apart; their pair gap is 9.2 px (5 px spacing + 4.2 px).

    Each label alone needs at most 5 px, so a broad phase padded by the single-label gap puts the two boxes in
    different grid cells and never tests the pair (OV-730). The pair-correct pad finds it.
    """

    spaced = {"_enforce_spacing_px": True, "_spacing_px": 5.0}
    labels = [
        _box_label(1000.0, 1000.0, 10.0, **spaced),
        _box_label(57.0, 0.0, 50.0, min_outer_start_radius_px=10.0),
        _box_label(2000.0, 1000.0, 10.0, **spaced),
        _box_label(115.0, 0.0, 50.0, **spaced),
    ]
    assert _count_label_overlaps(labels, 1000, use_min_gap=True) == 1
