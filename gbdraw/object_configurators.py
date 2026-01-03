#!/usr/bin/env python
# coding: utf-8

"""
Object configurators (compatibility façade).

Historically `gbdraw.object_configurators` contained multiple configurator classes.
To improve cohesion and make library usage clearer, implementations were moved into:

- `gbdraw.configurators.*`

The original import paths remain available via re-exports.
"""

from .configurators.blast import BlastMatchConfigurator  # type: ignore[reportMissingImports]
from .configurators.features import FeatureDrawingConfigurator  # type: ignore[reportMissingImports]
from .configurators.gc import GcContentConfigurator, GcSkewConfigurator  # type: ignore[reportMissingImports]
from .configurators.legend import LegendDrawingConfigurator  # type: ignore[reportMissingImports]

__all__ = [
    "BlastMatchConfigurator",
    "FeatureDrawingConfigurator",
    "GcContentConfigurator",
    "GcSkewConfigurator",
    "LegendDrawingConfigurator",
]


