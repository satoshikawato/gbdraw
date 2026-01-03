#!/usr/bin/env python
# coding: utf-8

"""
Backward-compatible facade for legacy imports.

Historically `gbdraw.circular_object_groups` contained multiple responsibilities (and many classes)
for building SVG groups in the circular renderer.

These have been split into:
  - `gbdraw.groups.circular.*`

New code should import from `gbdraw.groups.circular`.
"""

from .groups.circular import (  # type: ignore[reportMissingImports]  # noqa: F401
    AxisGroup,
    DefinitionGroup,
    GcContentGroup,
    GcSkewGroup,
    LegendGroup,
    SeqRecordGroup,
    TickGroup,
)

__all__ = [
    "AxisGroup",
    "DefinitionGroup",
    "GcContentGroup",
    "GcSkewGroup",
    "LegendGroup",
    "SeqRecordGroup",
    "TickGroup",
]


