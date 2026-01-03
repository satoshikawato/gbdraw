#!/usr/bin/env python
# coding: utf-8

"""
Backward-compatible facade for legacy imports.

Historically `gbdraw.circular_feature_drawer` bundled multiple responsibilities:
  - GC content / GC skew track drawing
  - feature path generation and drawing
  - definition text drawing
  - label drawing

These have been split into:
  - `gbdraw.drawers.circular.*`
  - `gbdraw.svg.*` for pure SVG path factories
  - `gbdraw.layout.*` for positioning math

New code should import from `gbdraw.drawers.circular`.
"""

from .drawers.circular import (  # type: ignore[reportMissingImports]  # noqa: F401
    DefinitionDrawer,
    FeatureDrawer,
    FeaturePathGenerator,
    GcContentDrawer,
    LabelDrawer,
    SkewDrawer,
)

__all__ = [
    "DefinitionDrawer",
    "FeatureDrawer",
    "FeaturePathGenerator",
    "GcContentDrawer",
    "LabelDrawer",
    "SkewDrawer",
]


