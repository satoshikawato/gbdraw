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
from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot  # type: ignore[reportMissingImports]


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
class DiagramOptions:
    """Bundled options for diagram assembly helpers."""

    config: GbdrawConfig | dict | None = None
    config_overrides: Mapping[str, object] | None = None
    colors: ColorOptions | None = None
    tracks: TrackOptions | None = None
    output: OutputOptions | None = None
    selected_features_set: Sequence[str] | None = None
    feature_table: DataFrame | None = None
    feature_table_file: str | None = None
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
    losatp_threads: int | None = None
    protein_blastp_max_hits: int = 5
    protein_blastp_candidate_limit: int | None = None
    align_orthogroup_feature: str | None = None
    evalue: float = 1e-5
    bitscore: float = 50.0
    identity: float = 70.0
    alignment_length: int = 0


__all__ = [
    "ColorOptions",
    "DiagramOptions",
    "OutputOptions",
    "TrackOptions",
]
