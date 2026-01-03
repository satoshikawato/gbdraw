#!/usr/bin/env python
# coding: utf-8

"""
Backward-compatible facade for legacy imports.

Historically `gbdraw.create_feature_objects` mixed multiple responsibilities:
  - feature color table preprocessing
  - FeatureObject/GeneObject/RepeatObject construction from Biopython records
  - track assignment / overlap handling
  - exon/intron coordinate normalization helpers
  - small drawing helpers (arrow shoulder)

These have been split into:
  - `gbdraw.features.colors`
  - `gbdraw.features.factory`
  - `gbdraw.features.tracks`
  - `gbdraw.features.coordinates`
  - `gbdraw.svg.arrows`

New code should import from those modules directly.
"""

# classes (legacy import path support)
from .feature_objects import FeatureObject, GeneObject, RepeatObject  # noqa: F401

# colors
from .features.colors import (  # type: ignore[reportMissingImports]  # noqa: F401
    get_color,
    preprocess_color_tables,
)

# factory
from .features.factory import (  # type: ignore[reportMissingImports]  # noqa: F401
    create_feature_dict,
    create_feature_object,
    create_gene_object,
    create_repeat_object,
)

# tracks
from .features.tracks import (  # type: ignore[reportMissingImports]  # noqa: F401
    arrange_feature_tracks,
    calculate_feature_metrics,
    check_feature_overlap,
    find_best_track,
    get_feature_ends,
)

# coordinates
from .features.coordinates import (  # type: ignore[reportMissingImports]  # noqa: F401
    get_coordinate,
    get_exon_and_intron_coordinates,
    get_exon_coordinate,
    get_intron_coordinate,
    get_strand,
)

# drawing helpers
from .svg.arrows import set_arrow_shoulder  # type: ignore[reportMissingImports]  # noqa: F401

__all__ = [
    # classes
    "FeatureObject",
    "GeneObject",
    "RepeatObject",
    # colors
    "get_color",
    "preprocess_color_tables",
    # factory
    "create_feature_dict",
    "create_feature_object",
    "create_gene_object",
    "create_repeat_object",
    # tracks
    "arrange_feature_tracks",
    "calculate_feature_metrics",
    "check_feature_overlap",
    "find_best_track",
    "get_feature_ends",
    # coordinates
    "get_coordinate",
    "get_exon_and_intron_coordinates",
    "get_exon_coordinate",
    "get_intron_coordinate",
    "get_strand",
    # drawing helpers
    "set_arrow_shoulder",
]


