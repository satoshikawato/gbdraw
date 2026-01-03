#!/usr/bin/env python
# coding: utf-8

"""
Backward-compatible facade for legacy imports.

Historically `gbdraw.linear_path_drawer` included:
  - layout factors for linear feature tracks
  - SVG path generators for linear tracks/features
  - misc helpers for label/track placement

These have been split into:
  - `gbdraw.layout.linear`
  - `gbdraw.svg.linear_features`
  - `gbdraw.svg.linear_tracks`
  - (coordinate normalization is in `gbdraw.layout.linear_coords`)
"""

from .layout.linear import calculate_feature_position_factors_linear  # noqa: F401
from .svg.linear_features import (  # noqa: F401
    calculate_normalized_arrow_length,
    construct_arrowhead_path,
    create_arrowhead_path_linear,
    create_intron_path_linear,
    create_rectangle_path_linear,
    get_arrow_strand_positions,
    normalize_feature_positions,
)
from .svg.linear_tracks import (  # noqa: F401
    calculate_corrdinate,
    calculate_gc_content_path_desc,
    calculate_gc_skew_path_desc,
)

__all__ = [
    # layout
    "calculate_feature_position_factors_linear",
    # linear features
    "calculate_normalized_arrow_length",
    "construct_arrowhead_path",
    "create_arrowhead_path_linear",
    "create_intron_path_linear",
    "create_rectangle_path_linear",
    "get_arrow_strand_positions",
    "normalize_feature_positions",
    # linear tracks
    "calculate_corrdinate",
    "calculate_gc_content_path_desc",
    "calculate_gc_skew_path_desc",
]


