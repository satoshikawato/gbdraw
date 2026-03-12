"""Track specification helpers (public API layer).

Track specs are experimental. At the moment they are only consumed by the
circular diagram assembler; linear support will follow in a later release.
"""

from gbdraw.tracks import (  # type: ignore[reportMissingImports]
    CircularCustomRule,
    CircularTrackPlacement,
    LinearTrackPlacement,
    ScalarSpec,
    TrackSpec,
    TrackSpecParseError,
    get_circular_feature_type_union,
    iter_compiled_custom_rules,
    load_circular_track_specs,
    normalize_circular_track_specs,
    parse_track_spec,
    parse_track_specs,
)

__all__ = [
    "CircularCustomRule",
    "CircularTrackPlacement",
    "LinearTrackPlacement",
    "ScalarSpec",
    "TrackSpec",
    "TrackSpecParseError",
    "get_circular_feature_type_union",
    "iter_compiled_custom_rules",
    "load_circular_track_specs",
    "normalize_circular_track_specs",
    "parse_track_spec",
    "parse_track_specs",
]


