#!/usr/bin/env python
# coding: utf-8

"""
Canvas configurators (compatibility façade).

Historically this module contained both `CircularCanvasConfigurator` and `LinearCanvasConfigurator`.
To improve cohesion and make library usage clearer, implementations were moved into:

- `gbdraw.canvas.circular`
- `gbdraw.canvas.linear`

The original import paths remain available via re-exports.
"""

from .canvas.circular import CircularCanvasConfigurator  # type: ignore[reportMissingImports]
from .canvas.linear import LinearCanvasConfigurator  # type: ignore[reportMissingImports]

__all__ = [
    "CircularCanvasConfigurator",
    "LinearCanvasConfigurator",
]


