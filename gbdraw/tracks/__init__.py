"""Multi-track specification (foundation, internal).

This package introduces *data models* for describing arbitrary track layouts for both
linear and circular diagrams. The current renderer does not yet consume these models;
they are a foundation for the upcoming multi-track refactor.
"""

from .spec import (  # type: ignore[reportMissingImports]
    CircularTrackPlacement,
    LinearTrackPlacement,
    ScalarSpec,
    TrackSpec,
)
from .parser import (  # type: ignore[reportMissingImports]
    TrackSpecParseError,
    parse_track_spec,
    parse_track_specs,
)

__all__ = [
    "CircularTrackPlacement",
    "LinearTrackPlacement",
    "ScalarSpec",
    "TrackSpec",
    "TrackSpecParseError",
    "parse_track_spec",
    "parse_track_specs",
]


