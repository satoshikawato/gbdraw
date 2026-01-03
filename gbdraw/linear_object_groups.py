#!/usr/bin/env python
# coding: utf-8

"""
Backward-compatible facade for legacy imports.

Historically `gbdraw.linear_object_groups` contained multiple responsibilities (and many classes)
for building SVG groups in the linear renderer.

These have been split into:
  - `gbdraw.groups.linear.*`

New code should import from `gbdraw.groups.linear`.
"""

from .groups.linear import (  # noqa: F401
    DefinitionGroup,
    GcContentGroup,
    GcSkewGroup,
    LegendGroup,
    LengthBarGroup,
    PairWiseMatchGroup,
    SeqRecordGroup,
)

__all__ = [
    "DefinitionGroup",
    "GcContentGroup",
    "GcSkewGroup",
    "LegendGroup",
    "LengthBarGroup",
    "PairWiseMatchGroup",
    "SeqRecordGroup",
]


