"""Diagram assembly internals.

This package holds the implementation of diagram composition (placing groups on canvases).
It is **not** the stable public API; see `gbdraw.api` for the public-ish layer.

Subpackages:
- circular: Circular diagram assembly
- linear: Linear diagram assembly
"""

from . import circular
from . import linear

__all__ = ["circular", "linear"]


