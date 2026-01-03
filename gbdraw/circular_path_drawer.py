#!/usr/bin/env python
# coding: utf-8

"""
Backward-compatible facade for legacy imports.

Historically `gbdraw.circular_path_drawer` included:
  - layout/radius factors for circular features
  - SVG path generation for circular tracks/features/ticks
  - some text helpers
  - re-exports of feature coordinate helpers

These have been split into:
  - `gbdraw.layout.circular`
  - `gbdraw.svg.circular_tracks`
  - `gbdraw.svg.circular_features`
  - `gbdraw.svg.circular_ticks`
  - `gbdraw.svg.text_path`
  - (feature coordinate helpers live in `gbdraw.features.coordinates`)

New code should import from the above modules directly.
"""

from .layout.circular import calculate_feature_position_factors_circular  # noqa: F401
from .svg.circular_tracks import (  # noqa: F401
    draw_circle_path,
    generate_circle_path_desc,
    generate_circular_gc_content_path_desc,
    generate_circular_gc_skew_path_desc,
)
from .svg.circular_features import (  # noqa: F401
    generate_circular_arrowhead_path,
    generate_circular_intron_path,
    generate_circular_rectangle_path,
)
from .svg.circular_ticks import (  # noqa: F401
    generate_circular_tick_labels,
    generate_circular_tick_paths,
    set_tick_label_anchor_value,
)
from .svg.text_path import generate_name_path, generate_text_path  # noqa: F401

# legacy re-exports (used by older code)
from .features.coordinates import get_exon_and_intron_coordinates  # noqa: F401

__all__ = [
    # layout
    "calculate_feature_position_factors_circular",
    # tracks
    "draw_circle_path",
    "generate_circle_path_desc",
    "generate_circular_gc_content_path_desc",
    "generate_circular_gc_skew_path_desc",
    # features
    "generate_circular_arrowhead_path",
    "generate_circular_intron_path",
    "generate_circular_rectangle_path",
    # ticks
    "generate_circular_tick_labels",
    "generate_circular_tick_paths",
    "set_tick_label_anchor_value",
    # text helpers
    "generate_name_path",
    "generate_text_path",
    # legacy re-exports
    "get_exon_and_intron_coordinates",
]


