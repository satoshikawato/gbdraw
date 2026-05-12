"""Track specification helpers (public API layer).

Track specs are experimental. At the moment they are only consumed by the
circular diagram assembler; linear support will follow in a later release.
"""

from gbdraw.tracks import (  # type: ignore[reportMissingImports]
    CircularTrackLayoutContext,
    CircularTrackPlacement,
    CircularTrackRendererName,
    CircularTrackSlot,
    LinearTrackPlacement,
    ResolvedCircularTrackSlot,
    ScalarSpec,
    SUPPORTED_CIRCULAR_TRACK_RENDERERS,
    TrackSpec,
    TrackSpecParseError,
    circular_track_slots_from_order,
    default_circular_track_slots,
    parse_circular_track_slot,
    parse_circular_track_slots,
    parse_track_spec,
    parse_track_specs,
    resolve_circular_track_slots,
)

__all__ = [
    "CircularTrackLayoutContext",
    "CircularTrackPlacement",
    "CircularTrackRendererName",
    "CircularTrackSlot",
    "LinearTrackPlacement",
    "ResolvedCircularTrackSlot",
    "ScalarSpec",
    "SUPPORTED_CIRCULAR_TRACK_RENDERERS",
    "TrackSpec",
    "TrackSpecParseError",
    "circular_track_slots_from_order",
    "default_circular_track_slots",
    "parse_circular_track_slot",
    "parse_circular_track_slots",
    "parse_track_spec",
    "parse_track_specs",
    "resolve_circular_track_slots",
]


