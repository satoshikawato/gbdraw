#!/usr/bin/env python
# coding: utf-8

"""
Backward-compatible facade for legacy imports.

Historically `gbdraw.data_processing` contained several responsibilities:
  - sequence statistics (GC%, dinucleotide skew)
  - sliding window helpers
  - legend table preparation
  - label placement (circular/linear)

These have been split into:
  - `gbdraw.analysis.*`
  - `gbdraw.legend.*`
  - `gbdraw.labels.placement`

New code should import from those modules directly.
"""

from .analysis.gc import calculate_gc_percent  # noqa: F401
from .analysis.skew import (  # noqa: F401
    calculate_dinucleotide_skew,
    skew_df,
    sliding_window,
)
from .legend.table import prepare_legend_table  # noqa: F401
from .labels.placement import (  # noqa: F401
    calculate_angle_degrees,
    calculate_angle_for_y,
    calculate_coordinates,
    check_label_overlap,
    euclidean_distance,
    find_lowest_available_track,
    improved_label_placement_fc,
    place_labels_on_arc_fc,
    prepare_label_list,
    prepare_label_list_linear,
    rearrange_labels_fc,
    sort_labels,
    x_overlap,
    y_overlap,
)

__all__ = [
    # analysis/gc
    "calculate_gc_percent",
    # analysis/skew
    "calculate_dinucleotide_skew",
    "skew_df",
    "sliding_window",
    # legend
    "prepare_legend_table",
    # labels/placement
    "calculate_angle_degrees",
    "calculate_angle_for_y",
    "calculate_coordinates",
    "check_label_overlap",
    "euclidean_distance",
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


