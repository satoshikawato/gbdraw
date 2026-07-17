"""Region annotations for circular and linear diagrams."""

from .io import annotation_sets_from_dataframe, read_annotation_table
from .layout import (
    PlacedRegionAnnotation,
    ResolvedAnnotationTrack,
    assign_annotation_lanes,
    layout_annotation_track,
)
from .legend import sync_annotation_legend_entries
from .models import (
    AnnotationOptions,
    AnnotationSet,
    AnnotationTrackParams,
    CoordinateSpan,
    FeatureSelector,
    FeatureSpan,
    HatchStyle,
    RegionAnnotation,
    RegionAnnotationStyle,
    ResolvedAnnotationBundle,
    ResolvedRegionAnnotation,
    ResolutionWarning,
    annotation_track_params_from_mapping,
    effective_annotation_style,
    parse_feature_selector,
)
from .resolve import resolve_annotation_set, resolve_annotations

__all__ = [
    "AnnotationOptions",
    "AnnotationSet",
    "AnnotationTrackParams",
    "CoordinateSpan",
    "FeatureSelector",
    "FeatureSpan",
    "HatchStyle",
    "PlacedRegionAnnotation",
    "RegionAnnotation",
    "RegionAnnotationStyle",
    "ResolvedAnnotationBundle",
    "ResolvedAnnotationTrack",
    "ResolvedRegionAnnotation",
    "ResolutionWarning",
    "annotation_sets_from_dataframe",
    "annotation_track_params_from_mapping",
    "effective_annotation_style",
    "assign_annotation_lanes",
    "layout_annotation_track",
    "parse_feature_selector",
    "read_annotation_table",
    "resolve_annotation_set",
    "resolve_annotations",
    "sync_annotation_legend_entries",
]
