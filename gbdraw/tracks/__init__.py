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
from .circular import (  # type: ignore[reportMissingImports]
    CircularTrackRendererName,
    CircularTrackSlot,
    SUPPORTED_CIRCULAR_TRACK_RENDERERS,
    circular_track_slots_from_order,
    default_circular_track_slots,
    parse_circular_track_slot,
    parse_circular_track_slots,
)
from ..diagrams.circular.slot_layout import (  # type: ignore[reportMissingImports]
    CircularSlotFootprint,
    CircularTrackLayoutContext,
    ResolvedCircularTrackSlot,
    resolve_circular_track_slots,
)
from .parser import (  # type: ignore[reportMissingImports]
    TrackSpecParseError,
    parse_track_spec,
    parse_track_specs,
)

__all__ = [
    "CircularTrackPlacement",
    "CircularTrackLayoutContext",
    "CircularTrackRendererName",
    "CircularSlotFootprint",
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


