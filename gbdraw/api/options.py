"""Option bundles for the public API.

These dataclasses provide a lighter-weight entry point than the long list of
keyword arguments in assemble_* helpers. They are optional and additive.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Mapping, Sequence

from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.analysis.collinearity import (  # type: ignore[reportMissingImports]
    CollinearityBlock,
    CollinearityAnchorMode,
    CollinearityColorMode,
    CollinearityParameters,
    CollinearityResult,
    CollinearitySearchScope,
    LosslessCollinearityParameters,
)
from gbdraw.analysis.collinearity_units import CollinearityUnitMode  # type: ignore[reportMissingImports]
from gbdraw.analysis.protein_colinearity import OrthogroupResult  # type: ignore[reportMissingImports]
from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]
from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot  # type: ignore[reportMissingImports]
from gbdraw.annotations import AnnotationOptions


@dataclass(frozen=True)
class ColorOptions:
    """Color table and palette inputs."""

    color_table: DataFrame | None = None
    color_table_file: str | None = None
    default_colors: DataFrame | None = None
    default_colors_palette: str = "default"
    default_colors_file: str | None = None


@dataclass(frozen=True)
class TrackOptions:
    """Track layout options."""

    circular_track_slots: Sequence[str | CircularTrackSlot] | None = None
    circular_track_axis_index: int | None = None
    linear_track_slots: Sequence[str | LinearTrackSlot] | None = None
    linear_track_axis_index: int | None = None
    center_reserved_radius: float | None = None


@dataclass(frozen=True)
class OutputOptions:
    """High-level output/layout options used by the assemblers."""

    output_prefix: str = "out"
    legend: str = "right"
    plot_title_position: Literal["none", "center", "top", "bottom"] | None = None


@dataclass(frozen=True)
class CircularMultiRecordOptions:
    """Layout values used only by circular multi-record canvases."""

    multi_record_size_mode: Literal["linear", "auto", "equal", "sqrt"] = "auto"
    multi_record_min_radius_ratio: float = 0.55
    multi_record_column_gap_ratio: float = 0.10
    multi_record_row_gap_ratio: float = 0.05
    multi_record_positions: Sequence[str] | None = None


@dataclass(frozen=True)
class DiagramOptions:
    """Bundled options for diagram assembly helpers."""

    config: GbdrawConfig | dict | None = None
    config_overrides: Mapping[str, object] | None = None
    colors: ColorOptions | None = None
    tracks: TrackOptions | None = None
    annotations: AnnotationOptions | None = None
    output: OutputOptions | None = None
    selected_features_set: Sequence[str] | None = None
    feature_table: DataFrame | None = None
    feature_table_file: str | None = None
    feature_visibility_table: DataFrame | None = None
    feature_visibility_table_file: str | None = None
    label_whitelist_table: DataFrame | None = None
    label_whitelist_file: str | None = None
    qualifier_priority_table: DataFrame | None = None
    qualifier_priority_file: str | None = None
    label_override_table: DataFrame | None = None
    label_override_file: str | None = None
    feature_shapes: Mapping[str, str] | None = None
    dinucleotide: str = "GC"
    window: int | None = None
    step: int | None = None
    depth_window: int | None = None
    depth_step: int | None = None
    depth_table: DataFrame | None = None
    depth_file: str | None = None
    depth_tables: Sequence[DataFrame] | None = None
    depth_files: Sequence[str] | None = None
    depth_track_tables: Sequence[Sequence[DataFrame | None]] | None = None
    depth_track_files: Sequence[Sequence[str | None]] | None = None
    depth_track_labels: Sequence[str] | None = None
    depth_track_colors: Sequence[str] | None = None
    depth_track_heights: Sequence[float | str | None] | None = None
    depth_track_large_tick_intervals: Sequence[float | str | None] | None = None
    depth_track_small_tick_intervals: Sequence[float | str | None] | None = None
    depth_track_tick_font_sizes: Sequence[float | str | None] | None = None
    conservation_blast_files: Sequence[str] | None = None
    conservation_dataframes: Sequence[DataFrame] | None = None
    conservation_reference: Literal["query", "subject", "auto"] = "auto"
    conservation_labels: Sequence[str] | None = None
    conservation_colors: Sequence[str] | None = None
    conservation_ring_width: float | None = None
    conservation_ring_gap: float | None = None
    plot_title: str | None = None
    plot_title_font_size: float | None = None
    keep_full_definition_with_plot_title: bool = False
    species: str | None = None
    strain: str | None = None
    blast_files: Sequence[str] | None = None
    protein_comparisons: Sequence[DataFrame] | None = None
    orthogroups: OrthogroupResult | None = None
    protein_blastp_mode: Literal["none", "pairwise", "orthogroup", "collinear"] = "none"
    pairwise_match_style: Literal["ribbon", "curve"] = "ribbon"
    collinearity_blocks: CollinearityResult | Sequence[CollinearityBlock] | None = None
    collinearity_params: CollinearityParameters | LosslessCollinearityParameters | None = None
    collinearity_unit_mode: CollinearityUnitMode | str = "auto"
    collinearity_anchor_mode: CollinearityAnchorMode | str = "rbh"
    collinearity_search_scope: CollinearitySearchScope | str = "adjacent"
    collinearity_color_mode: CollinearityColorMode | str = "orientation"
    losatp_bin: str = "losat"
    ncbi_blastp_bin: str | None = None
    losatp_threads: int | None = None
    protein_blastp_max_hits: int = 5
    protein_blastp_candidate_limit: int | None = None
    orthogroup_membership_mode: Literal["anchor_core_v1"] | str = "anchor_core_v1"
    orthogroup_member_max_hits: int = 5
    collinear_max_paralog_links_per_orthogroup: int = 2
    align_orthogroup_feature: str | None = None
    evalue: float = 1e-5
    bitscore: float = 50.0
    identity: float = 70.0
    alignment_length: int = 0


_DEFAULT_DIAGRAM_OPTIONS = DiagramOptions()
_CIRCULAR_ONLY_DIAGRAM_OPTION_NAMES = (
    "conservation_blast_files",
    "conservation_dataframes",
    "conservation_reference",
    "conservation_labels",
    "conservation_colors",
    "conservation_ring_width",
    "conservation_ring_gap",
    "keep_full_definition_with_plot_title",
    "species",
    "strain",
)
_LINEAR_ONLY_DIAGRAM_OPTION_NAMES = (
    "depth_track_heights",
    "blast_files",
    "protein_comparisons",
    "orthogroups",
    "protein_blastp_mode",
    "pairwise_match_style",
    "collinearity_blocks",
    "collinearity_params",
    "collinearity_unit_mode",
    "collinearity_anchor_mode",
    "collinearity_search_scope",
    "collinearity_color_mode",
    "losatp_bin",
    "ncbi_blastp_bin",
    "losatp_threads",
    "protein_blastp_max_hits",
    "protein_blastp_candidate_limit",
    "orthogroup_membership_mode",
    "orthogroup_member_max_hits",
    "collinear_max_paralog_links_per_orthogroup",
    "align_orthogroup_feature",
)


def _validate_diagram_options_mode(
    options: DiagramOptions,
    *,
    mode: Literal["circular", "circular_multi", "linear"],
) -> None:
    incompatible_names = (
        _CIRCULAR_ONLY_DIAGRAM_OPTION_NAMES
        if mode == "linear"
        else _LINEAR_ONLY_DIAGRAM_OPTION_NAMES
    )
    non_default_names = [
        name
        for name in incompatible_names
        if getattr(options, name) != getattr(_DEFAULT_DIAGRAM_OPTIONS, name)
    ]
    if non_default_names:
        raise ValidationError(
            f"{mode} builder does not support non-default DiagramOptions fields: "
            f"{', '.join(non_default_names)}"
        )


__all__ = [
    "CircularMultiRecordOptions",
    "AnnotationOptions",
    "ColorOptions",
    "DiagramOptions",
    "OutputOptions",
    "TrackOptions",
]
