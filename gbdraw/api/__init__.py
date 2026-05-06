"""Public, stable-ish library API for gbdraw.

This package is the official entry point for pipeline/library usage (as opposed to CLI
entry points). It stays intentionally thin and mostly re-exports internal building
blocks through a stable namespace.
"""

from .canvas import CircularCanvasConfigurator, LinearCanvasConfigurator
from .config import GbdrawConfig, apply_config_overrides, load_default_config
from .configurators import (
    BlastMatchConfigurator,
    DepthConfigurator,
    FeatureDrawingConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
    LegendDrawingConfigurator,
)
from .diagram import (  # type: ignore[reportMissingImports]
    DEFAULT_SELECTED_FEATURES,
    assemble_circular_diagram_from_records,
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
    assemble_linear_diagram_from_rows,
    build_circular_diagram,
    build_linear_diagram,
)
from .io import (
    LinearInputRow,
    RegionSpec,
    RecordSelector,
    apply_region_specs,
    load_gbk_rows,
    load_gbks,
    load_gff_fasta_rows,
    load_gff_fasta,
    parse_record_selector,
    parse_record_selectors,
    parse_region_spec,
    parse_region_specs,
)
from gbdraw.analysis.collinearity import (  # type: ignore[reportMissingImports]
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityParameters,
    CollinearityResult,
    CollinearitySearchScope,
    LosslessCollinearityParameters,
    build_orthogroup_collinearity_blocks,
    iter_collinearity_search_pairs,
    normalize_collinearity_search_scope,
)
from .options import ColorOptions, DiagramOptions, OutputOptions, TrackOptions
from .render import parse_formats, render_to_bytes, save_figure, save_figure_to
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
    # config
    "GbdrawConfig",
    "apply_config_overrides",
    "load_default_config",
    # configurators
    "BlastMatchConfigurator",
    "DepthConfigurator",
    "FeatureDrawingConfigurator",
    "GcContentConfigurator",
    "GcSkewConfigurator",
    "LegendDrawingConfigurator",
    # diagrams
    "DEFAULT_SELECTED_FEATURES",
    "assemble_circular_diagram_from_records",
    "assemble_circular_diagram_from_record",
    "assemble_linear_diagram_from_records",
    "assemble_linear_diagram_from_rows",
    "build_circular_diagram",
    "build_linear_diagram",
    # io
    "LinearInputRow",
    "RegionSpec",
    "RecordSelector",
    "apply_region_specs",
    "load_gbk_rows",
    "load_gbks",
    "load_gff_fasta_rows",
    "load_gff_fasta",
    "parse_record_selector",
    "parse_record_selectors",
    "parse_region_spec",
    "parse_region_specs",
    # collinearity
    "CollinearityAnchor",
    "CollinearityBlock",
    "CollinearityParameters",
    "CollinearityResult",
    "CollinearitySearchScope",
    "LosslessCollinearityParameters",
    "build_orthogroup_collinearity_blocks",
    "iter_collinearity_search_pairs",
    "normalize_collinearity_search_scope",
    # options
    "ColorOptions",
    "DiagramOptions",
    "OutputOptions",
    "TrackOptions",
    # render
    "parse_formats",
    "render_to_bytes",
    "save_figure",
    "save_figure_to",
    # tracks (foundation)
    "CircularTrackPlacement",
    "LinearTrackPlacement",
    "ScalarSpec",
    "TrackSpec",
    "TrackSpecParseError",
    "parse_track_spec",
    "parse_track_specs",
]


