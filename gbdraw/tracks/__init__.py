"""Multi-track specification (foundation, internal).

These data models describe arbitrary track layouts for both linear and circular
diagrams. They are still experimental. Currently only the circular assembler
consumes a subset of these specs; linear support will follow in a later release.
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


