"""Canvas configuration helpers (internal).

These classes encapsulate canvas sizing/offset logic for circular and linear diagrams.
"""

from .circular import CircularCanvasConfigurator, CircularCanvasGeometry
from .linear import LinearCanvasConfigurator, LinearCanvasGeometry

__all__ = [
    "CircularCanvasConfigurator",
    "CircularCanvasGeometry",
    "LinearCanvasConfigurator",
    "LinearCanvasGeometry",
]


