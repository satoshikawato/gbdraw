#!/usr/bin/env python
# coding: utf-8

"""
Backward-compatible facade for legacy imports.

Historically `gbdraw.linear_feature_drawer` bundled multiple responsibilities:
  - feature drawing on linear canvas
  - label drawing on linear canvas
  - GC content / GC skew drawing on linear canvas

These have been split into:
  - `gbdraw.drawers.linear.*`
  - `gbdraw.svg.*` for pure SVG path factories

New code should import from `gbdraw.drawers.linear`.
"""

from .drawers.linear import (  # type: ignore[reportMissingImports]  # noqa: F401
    FeatureDrawer,
    FeaturePathGenerator,
    GcContentDrawer,
    LabelDrawer,
    SkewDrawer,
)

__all__ = [
    "FeatureDrawer",
    "FeaturePathGenerator",
    "GcContentDrawer",
    "LabelDrawer",
    "SkewDrawer",
]


