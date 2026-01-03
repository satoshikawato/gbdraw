"""Feature parsing / object construction helpers (internal)."""

from .colors import preprocess_color_tables, get_color
from .factory import (
    create_feature_dict,
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

__all__ = [
    # colors
    "get_color",
    "preprocess_color_tables",
    # factory
    "create_feature_dict",
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
]


