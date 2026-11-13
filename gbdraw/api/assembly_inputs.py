"""Typed inputs of the diagram assemblers in gbdraw.api.diagram.

``DiagramAssemblyInputs`` holds what both modes take, declared once;
``CircularAssemblyInputs`` and ``LinearAssemblyInputs`` add each mode's inputs
and defaults. ``CircularRecordContext`` holds the private per-record values
that only the single-record Circular assembler reads. Every field keeps the
name, type and default of the assembler keyword it replaced.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Mapping, Optional, Sequence

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw.analysis.collinearity import (
    CollinearityAnchorMode,
    CollinearityBlock,
    CollinearityColorMode,
    CollinearityResult,
    CollinearitySearchScope,
    LosslessCollinearityParameters,
)
from gbdraw.analysis.collinearity_units import CollinearityUnitMode
from gbdraw.analysis.conservation import ConservationInput, ConservationTrack
from gbdraw.analysis.depth_tracks import DepthTrackData, DepthTrackSpec, ResolvedDepthTracks
from gbdraw.analysis.feature_density import FeatureDensityLayer, FeatureDensityTrack
from gbdraw.analysis.protein_colinearity import (
    LosatpCacheManager,
    OrthogroupGraphResult,
    OrthogroupMembershipMode,
    OrthogroupResult,
    ProteinExtractionResult,
)
from gbdraw.annotations import ResolvedAnnotationBundle
from gbdraw.api.options import (
    CIRCULAR_DEFAULT_PLOT_TITLE_POSITION,
    LINEAR_DEFAULT_PLOT_TITLE_POSITION,
    AnnotationOptions,
    DepthTrackInput,
    LinearDiagramOptions,
    LosatSearchOptions,
)
from gbdraw.api.prepared import ResolvedFeatureInputs
from gbdraw.api.requests import RecordPresentation
from gbdraw.config.models import GbdrawConfig
from gbdraw.features.placement import ResolvedRecordFeatureInputs
from gbdraw.labels.linear import FeatureLabelRows
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.record_units import CircleUnits
from gbdraw.layout.record_placement import FeatureAnchorPlacement
from gbdraw.layout.similarity_alignment import SimilarityAlignmentPlan
from gbdraw.linear_comparison import LinearComparison
from gbdraw.mode_profiles import CIRCULAR_MODE_PROFILE, LINEAR_MODE_PROFILE
from gbdraw.render.groups.circular.definition import CircularDefinitionProfile
from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot


@dataclass(frozen=True, kw_only=True)
class DiagramAssemblyInputs:
    """The assembly inputs both modes take.

    If ``default_colors`` is None, the built-in default palette is loaded. If
    ``color_table`` is None and ``color_table_file`` is set, the file is read.
    If ``selected_features_set`` is None, the CLI default feature list is used.
    ``depth_tables`` and ``depth_files`` hold one entry per record of a
    multi-record figure. The ``_resolved_*`` values are inputs a builder read
    once (OV-720, OV-732); ``_resolved_depth_tracks`` replaces every depth
    input when set.
    """

    cfg: GbdrawConfig
    color_table: Optional[DataFrame] = None
    color_table_file: str | None = None
    default_colors: DataFrame | None = None
    default_colors_palette: str = "default"
    default_colors_file: str | None = None
    selected_features_set: Sequence[str] | None = None
    feature_visibility_table: DataFrame | None = None
    feature_visibility_table_file: str | None = None
    feature_shapes: Mapping[str, str] | None = None
    output_prefix: str = "out"
    legend: str = "right"
    dinucleotide: str = "GC"
    window: Optional[int] = None
    step: Optional[int] = None
    depth_window: Optional[int] = None
    depth_step: Optional[int] = None
    depth_tracks: Sequence[DepthTrackInput] | None = None
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
    # The feature density tracks a ``feature_density`` slot names by id; a stack the reader did not write draws
    # each of them after its features slot.
    feature_density_tracks: Sequence[FeatureDensityTrack] | None = None
    annotation_options: AnnotationOptions | None = None
    plot_title: str | None = None
    plot_title_font_size: float | None = None
    _resolved_feature_inputs: ResolvedFeatureInputs | None = None
    _resolved_record_features: tuple[ResolvedRecordFeatureInputs, ...] = ()
    _resolved_annotations: ResolvedAnnotationBundle | None = None
    _resolved_depth_tracks: ResolvedDepthTracks | None = None


@dataclass(frozen=True, kw_only=True)
class CircularAssemblyInputs(DiagramAssemblyInputs):
    """The Circular assembly inputs, with the Circular comparison thresholds."""

    conservation_blast_files: Sequence[ConservationInput] | None = None
    conservation_dataframes: Sequence[DataFrame] | None = None
    conservation_reference: Literal["query", "subject", "auto"] | str = "auto"
    conservation_labels: Sequence[str] | None = None
    conservation_colors: Sequence[str] | None = None
    conservation_ring_width: float | None = None
    conservation_ring_gap: float | None = None
    species: Optional[str] = None
    strain: Optional[str] = None
    plot_title_position: Literal["none", "top", "bottom"] = CIRCULAR_DEFAULT_PLOT_TITLE_POSITION
    keep_full_definition_with_plot_title: bool = False
    center_reserved_radius: float | None = None
    circular_track_slots: Sequence[str | CircularTrackSlot] | None = None
    # True when ``circular_track_slots`` is the stack a grid built for its records because the reader gave none.
    _implicit_track_stack: bool = False
    circular_track_axis_index: int | None = None
    evalue: float = CIRCULAR_MODE_PROFILE.comparison.evalue
    bitscore: float = CIRCULAR_MODE_PROFILE.comparison.bitscore
    identity: float = CIRCULAR_MODE_PROFILE.comparison.identity
    alignment_length: int = CIRCULAR_MODE_PROFILE.comparison.alignment_length


@dataclass(frozen=True, kw_only=True)
class LinearAssemblyInputs(DiagramAssemblyInputs):
    """The Linear assembly inputs, with the Linear comparison thresholds."""

    blast_files: Optional[Sequence[str]] = None
    linear_comparisons: Sequence[LinearComparison] | None = None
    protein_comparisons: Sequence[DataFrame] | None = None
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None = None
    losat_search: LosatSearchOptions | None = None
    pairwise_match_style: Literal["ribbon", "curve"] | str = LinearDiagramOptions.pairwise_match_style
    collinearity_blocks: CollinearityResult | Sequence[CollinearityBlock] | None = None
    collinearity_params: LosslessCollinearityParameters | None = None
    collinearity_unit_mode: CollinearityUnitMode | str = LinearDiagramOptions.collinearity_unit_mode
    collinearity_anchor_mode: CollinearityAnchorMode | str = LinearDiagramOptions.collinearity_anchor_mode
    collinearity_search_scope: CollinearitySearchScope | str = LinearDiagramOptions.collinearity_search_scope
    collinearity_color_mode: CollinearityColorMode | str = LinearDiagramOptions.collinearity_color_mode
    losatp_cache: LosatpCacheManager | None = None
    protein_extraction: ProteinExtractionResult | None = None
    orthogroup_membership_mode: OrthogroupMembershipMode | str = LinearDiagramOptions.orthogroup_membership_mode
    collinear_infer_orthogroups: bool = LinearDiagramOptions.collinear_infer_orthogroups
    collinear_max_paralog_links_per_orthogroup: int = LinearDiagramOptions.collinear_max_paralog_links_per_orthogroup
    similarity_alignment: SimilarityAlignmentPlan | None = None
    depth_track_heights: Sequence[float | str | None] | None = None
    linear_track_slots: Sequence[str | LinearTrackSlot] | None = None
    linear_track_axis_index: int | None = None
    plot_title_position: Literal["center", "top", "bottom"] = LINEAR_DEFAULT_PLOT_TITLE_POSITION
    evalue: float = LINEAR_MODE_PROFILE.comparison.evalue
    bitscore: float = LINEAR_MODE_PROFILE.comparison.bitscore
    identity: float = LINEAR_MODE_PROFILE.comparison.identity
    alignment_length: int = LINEAR_MODE_PROFILE.comparison.alignment_length
    _alignment_anchor_centers: Sequence[float | None] | None = None
    _feature_anchor: FeatureAnchorPlacement | None = None
    _feature_label_rows: FeatureLabelRows | None = None


@dataclass(frozen=True)
class CircularOneCircleInputs:
    """The records of a figure drawn on one circle: every record and the transform it is drawn by, by index,
    and where each drawn record lies on the circle (``circle_units``)."""

    records: tuple[SeqRecord, ...]
    transforms: tuple[RecordDisplayTransform, ...]
    layout: CircleUnits
    # Records shorter than this are left out of the circle and counted (``layout.hidden``).
    min_length_bp: int
    # Whether a band outside the tracks tells the records apart and names them.
    record_band: bool


@dataclass(frozen=True, kw_only=True)
class CircularRecordContext:
    """Per-record values a multi-record or builder caller hands one Circular record.

    ``precomputed_depth_tracks`` is None when the canvas has no depth input;
    ``[]`` means depth input without a cell for this record, which keeps the
    record's depth slot.
    """

    definition_profile: CircularDefinitionProfile = "full"
    layout_definition_profile: CircularDefinitionProfile | None = None
    precomputed_depth_df: DataFrame | None = None
    precomputed_depth_track_specs: Sequence[DepthTrackSpec] | None = None
    precomputed_depth_tracks: Sequence[DepthTrackData] | None = None
    precomputed_depth_track_count: int | None = None
    shared_depth_max: float | None = None
    # The record's feature density layers, on the scale of the canvas's shared legend; None: the record's own.
    feature_density: Mapping[str, FeatureDensityLayer] | None = None
    precomputed_conservation_tracks: Sequence[ConservationTrack] | None = None
    # Each ring's Legend name by ``track_index`` (``LegendTable.ring_names`` of a Multi-Record Canvas's shared
    # Legend); None: the record names its rings from its own Legend.
    conservation_ring_names: Mapping[int, str] | None = None
    annotation_record_index: int = 0
    definition_group_id: str | None = None
    record_transform: RecordDisplayTransform | None = None
    record_presentation: RecordPresentation | None = None
    requested_canvas: bool = True
    # Several records drawn one after another on one circle (this record is the first drawn).
    one_circle: CircularOneCircleInputs | None = None
