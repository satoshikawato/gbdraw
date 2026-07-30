"""Feature parsing / object construction helpers (internal)."""

from .colors import preprocess_color_tables, get_color
from .factory import (
    FeatureBuildResult,
    create_feature_dict,
    create_feature_layers,
    create_feature_object,
    create_gene_object,
    create_repeat_object,
)
from .tracks import (
    arrange_feature_tracks,
    calculate_feature_metrics,
    find_best_track,
    get_feature_ends,
)
from .coordinates import (
    get_coordinate,
    get_exon_and_intron_coordinates,
    get_exon_coordinate,
    get_intron_coordinate,
    get_strand,
)
from .shapes import (
    DEFAULT_DIRECTIONAL_FEATURE_TYPES,
    DEFAULT_FEATURE_RENDERINGS,
    FEATURE_RENDERING_VALUES,
    FeatureRendering,
    default_feature_rendering,
    normalize_feature_shape,
    normalize_feature_shape_overrides,
    parse_feature_shape_assignment,
    parse_feature_shape_overrides,
    resolve_directional_feature_types,
    resolve_feature_rendering,
    resolve_underlay_feature_types,
)
from .visibility import (
    compile_feature_visibility_rules,
    read_feature_visibility_file,
    resolve_candidate_feature_types,
    should_include_feature_in_analysis,
    should_render_feature,
)

__all__ = [
    # colors
    "get_color",
    "preprocess_color_tables",
    # factory
    "FeatureBuildResult",
    "create_feature_dict",
    "create_feature_layers",
    "create_feature_object",
    "create_gene_object",
    "create_repeat_object",
    # tracks
    "arrange_feature_tracks",
    "calculate_feature_metrics",
    "find_best_track",
    "get_feature_ends",
    # coordinates
    "get_coordinate",
    "get_exon_and_intron_coordinates",
    "get_exon_coordinate",
    "get_intron_coordinate",
    "get_strand",
    # shapes
    "DEFAULT_DIRECTIONAL_FEATURE_TYPES",
    "DEFAULT_FEATURE_RENDERINGS",
    "FEATURE_RENDERING_VALUES",
    "FeatureRendering",
    "default_feature_rendering",
    "normalize_feature_shape",
    "normalize_feature_shape_overrides",
    "parse_feature_shape_assignment",
    "parse_feature_shape_overrides",
    "resolve_directional_feature_types",
    "resolve_feature_rendering",
    "resolve_underlay_feature_types",
    # visibility
    "compile_feature_visibility_rules",
    "read_feature_visibility_file",
    "resolve_candidate_feature_types",
    "should_include_feature_in_analysis",
    "should_render_feature",
]
