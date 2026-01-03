#!/usr/bin/env python
# coding: utf-8

"""
Legacy label placement utilities (compatibility façade).

Historically this module contained both circular and linear label placement logic.
To improve cohesion, implementations were moved into:

- `gbdraw.labels.placement_legacy_circular`
- `gbdraw.labels.placement_legacy_linear`

The public API remains available from this module via re-exports.
"""

from .placement_legacy_circular import (  # type: ignore[reportMissingImports]
    improved_label_placement_fc,
    place_labels_on_arc_fc,
    prepare_label_list,
    rearrange_labels_fc,
    sort_labels,
    x_overlap,
    y_overlap,
)
from .placement_legacy_linear import (  # type: ignore[reportMissingImports]
    check_label_overlap,
    edit_available_tracks,
    find_lowest_available_track,
    prepare_label_list_linear,
)

__all__ = [
    "check_label_overlap",
    "edit_available_tracks",
    "find_lowest_available_track",
    "improved_label_placement_fc",
    "place_labels_on_arc_fc",
    "prepare_label_list",
    "prepare_label_list_linear",
    "rearrange_labels_fc",
    "sort_labels",
    "x_overlap",
    "y_overlap",
]


