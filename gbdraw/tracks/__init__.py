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
from .circular import (  # type: ignore[reportMissingImports]
    CircularCustomRule,
    get_circular_feature_type_union,
    iter_compiled_custom_rules,
    load_circular_track_specs,
    normalize_circular_track_specs,
)

__all__ = [
    "CircularTrackPlacement",
    "LinearTrackPlacement",
    "ScalarSpec",
    "TrackSpec",
    "TrackSpecParseError",
    "parse_track_spec",
    "parse_track_specs",
    "CircularCustomRule",
    "get_circular_feature_type_union",
    "iter_compiled_custom_rules",
    "load_circular_track_specs",
    "normalize_circular_track_specs",
]


