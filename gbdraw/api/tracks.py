"""Track specification helpers (public API layer).

Track specs are experimental. At the moment they are only consumed by the
circular diagram assembler; linear support will follow in a later release.
"""

from gbdraw.tracks import (  # type: ignore[reportMissingImports]
    CircularTrackPlacement,
    LinearTrackPlacement,
    ScalarSpec,
    TrackSpec,
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


