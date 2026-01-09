"""SVG object group builders (internal).

Subpackages:
- circular: Group builders for circular diagrams
- linear: Group builders for linear diagrams
"""

from . import circular
from . import linear

__all__ = ["circular", "linear"]
