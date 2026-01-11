"""Public, stable-ish library API for gbdraw.

This package is intended for pipeline/library usage (as opposed to CLI entry points).
We keep it thin and mostly re-export internal building blocks.
"""

from .canvas import CircularCanvasConfigurator, LinearCanvasConfigurator
from .configurators import (
    BlastMatchConfigurator,
    FeatureDrawingConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
    LegendDrawingConfigurator,
)
from .diagram import (  # type: ignore[reportMissingImports]
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
)
from .tracks import (  # type: ignore[reportMissingImports]
    CircularTrackPlacement,
    LinearTrackPlacement,
    ScalarSpec,
    TrackSpec,
    TrackSpecParseError,
    parse_track_spec,
    parse_track_specs,
)

__all__ = [
    # canvas
    "CircularCanvasConfigurator",
    "LinearCanvasConfigurator",
    # configurators
    "BlastMatchConfigurator",
    "FeatureDrawingConfigurator",
    "GcContentConfigurator",
    "GcSkewConfigurator",
    "LegendDrawingConfigurator",
    # diagrams
    "assemble_circular_diagram_from_record",
    "assemble_linear_diagram_from_records",
    # tracks (foundation)
    "CircularTrackPlacement",
    "LinearTrackPlacement",
    "ScalarSpec",
    "TrackSpec",
    "TrackSpecParseError",
    "parse_track_spec",
    "parse_track_specs",
]


