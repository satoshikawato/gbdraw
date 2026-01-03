#!/usr/bin/env python
# coding: utf-8

"""
Backward-compatible facade for legacy imports.

Historically `gbdraw.utility_functions` accumulated a wide range of unrelated helpers.
These have been split into smaller, responsibility-focused modules:

  - Text / bbox helpers: `gbdraw.text`
  - Label filtering: `gbdraw.labels.filtering`
  - Label placement (legacy): `gbdraw.labels.placement_legacy`
  - Geometry helpers: `gbdraw.geometry`
  - Layout helpers: `gbdraw.layout.*`
  - Config mutation helpers: `gbdraw.config.modify`
  - Core helpers: `gbdraw.core.*`

New code should prefer importing from those modules directly.
"""

from .text import (  # noqa: F401
    calculate_bbox_dimensions,
    create_text_element,
    get_text_bbox_size_pixels,
    parse_mixed_content_text,
)
from .labels.filtering import (  # noqa: F401
    get_label_text,
    preprocess_label_filtering,
    read_filter_list_file,
    read_qualifier_priority_file,
)
from .labels.placement_legacy import (  # noqa: F401
    check_label_overlap,
    edit_available_tracks,
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
from .core.color import interpolate_color  # noqa: F401
from .layout.linear_coords import (  # noqa: F401
    normalize_position_linear,
    normalize_position_to_linear_track,
)
from .core.sequence import (  # noqa: F401
    check_feature_presence,
    create_dict_for_sequence_lengths,
    determine_length_parameter,
    determine_output_file_prefix,
    get_coordinates_of_longest_segment,
)
from .geometry import (  # noqa: F401
    calculate_angle_degrees,
    calculate_coordinates,
    deg_to_rad,
    euclidean_distance,
)
from .layout.common import calculate_cds_ratio  # noqa: F401
from .config.modify import (  # noqa: F401
    modify_config_dict,
    suppress_gc_content_and_skew,
    update_config_value,
)

__all__ = [
    # text
    "calculate_bbox_dimensions",
    "create_text_element",
    "get_text_bbox_size_pixels",
    "parse_mixed_content_text",
    # labels (filtering)
    "get_label_text",
    "preprocess_label_filtering",
    "read_filter_list_file",
    "read_qualifier_priority_file",
    # labels (placement legacy)
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
    # color
    "interpolate_color",
    # layout coords
    "normalize_position_linear",
    "normalize_position_to_linear_track",
    # core/sequence
    "check_feature_presence",
    "create_dict_for_sequence_lengths",
    "determine_length_parameter",
    "determine_output_file_prefix",
    "get_coordinates_of_longest_segment",
    # geometry
    "calculate_angle_degrees",
    "calculate_coordinates",
    "deg_to_rad",
    "euclidean_distance",
    # layout common
    "calculate_cds_ratio",
    # config modify
    "modify_config_dict",
    "suppress_gc_content_and_skew",
    "update_config_value",
]


