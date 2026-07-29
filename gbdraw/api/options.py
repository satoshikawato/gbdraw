"""Option bundles for the public API.

These dataclasses provide a lighter-weight entry point than the long list of
keyword arguments in assemble_* helpers. They are optional and additive.
"""

from __future__ import annotations

from dataclasses import dataclass, fields, replace
import math
from numbers import Integral, Real
from typing import Literal, Mapping, Sequence

from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.analysis.collinearity import (  # type: ignore[reportMissingImports]
    CollinearityBlock,
    CollinearityAnchorMode,
    CollinearityColorMode,
    CollinearityResult,
    CollinearitySearchScope,
    LosslessCollinearityParameters,
    normalize_collinearity_anchor_mode,
    normalize_collinearity_color_mode,
    normalize_collinearity_search_scope,
)
from gbdraw.analysis.collinearity_units import (  # type: ignore[reportMissingImports]
    CollinearityUnitMode,
    normalize_collinearity_unit_mode,
)
from gbdraw.analysis.conservation import (  # type: ignore[reportMissingImports]
    normalize_conservation_reference,
)
from gbdraw.analysis.protein_colinearity import (  # type: ignore[reportMissingImports]
    OrthogroupResult,
    normalize_orthogroup_membership_mode,
    normalize_protein_blastp_mode,
)
from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]
from gbdraw.linear_comparison import LinearComparison
from gbdraw.mode_profiles import (
    ComparisonThresholds,
    DiagramMode,
    get_mode_profile,
    resolve_mode_profile_overrides,
)
from gbdraw.tracks import (  # type: ignore[reportMissingImports]
    CircularTrackSlot,
    LinearTrackSlot,
    normalize_circular_track_slots_with_axis,
    normalize_linear_track_slots_with_axis,
    parse_circular_track_slots,
    parse_linear_track_slots,
)
from gbdraw.annotations import AnnotationOptions


@dataclass(frozen=True)
class ColorOptions:
    """Color table and palette inputs."""

    color_table: DataFrame | None = None
    color_table_file: str | None = None
    default_colors: DataFrame | None = None
    default_colors_palette: str = "default"
    default_colors_file: str | None = None


def _validate_track_configuration(
    values: Sequence[object] | None,
    *,
    axis_index: object,
    mode: Literal["circular", "linear"],
    expected_type: type[CircularTrackSlot] | type[LinearTrackSlot],
    slots_field_name: str,
    axis_field_name: str,
) -> int | None:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        if values is not None:
            raise ValidationError(f"{slots_field_name} must be a sequence.")
        parsed_slots = None
    else:
        if not all(isinstance(value, (str, expected_type)) for value in values):
            raise ValidationError(
                f"{slots_field_name} must contain strings or "
                f"{expected_type.__name__} values."
            )
        parser = (
            parse_circular_track_slots
            if mode == "circular"
            else parse_linear_track_slots
        )
        try:
            parsed_slots = parser(values)
        except (TypeError, ValueError) as exc:
            raise ValidationError(f"{slots_field_name}: {exc}") from exc

    if axis_index is None:
        return None
    if isinstance(axis_index, bool) or not isinstance(axis_index, Integral):
        raise ValidationError(f"{axis_field_name} must be an integer.")
    if parsed_slots is None:
        raise ValidationError(f"{axis_field_name} requires {slots_field_name}.")
    normalized_axis_index = int(axis_index)
    normalizer = (
        normalize_circular_track_slots_with_axis
        if mode == "circular"
        else normalize_linear_track_slots_with_axis
    )
    try:
        normalizer(parsed_slots, normalized_axis_index)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{axis_field_name}: {exc}") from exc
    return normalized_axis_index


def _validate_center_reserved_radius(value: object, *, field_name: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValidationError(f"{field_name} must be a finite number >= 0.")
    try:
        radius = float(value)
    except (OverflowError, TypeError, ValueError) as exc:
        raise ValidationError(f"{field_name} must be a finite number >= 0.") from exc
    if not math.isfinite(radius) or radius < 0:
        raise ValidationError(f"{field_name} must be a finite number >= 0.")
    return radius


def _validate_positive_real(value: object, *, field_name: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValidationError(f"{field_name} must be a finite number > 0 or None.")
    try:
        normalized = float(value)
    except (OverflowError, TypeError, ValueError) as exc:
        raise ValidationError(
            f"{field_name} must be a finite number > 0 or None."
        ) from exc
    if not math.isfinite(normalized) or normalized <= 0:
        raise ValidationError(f"{field_name} must be a finite number > 0 or None.")
    return normalized


def _validate_positive_int(
    value: object,
    *,
    field_name: str,
    allow_none: bool = False,
) -> int | None:
    if value is None and allow_none:
        return None
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) <= 0:
        suffix = " or None" if allow_none else ""
        raise ValidationError(f"{field_name} must be a positive integer{suffix}.")
    return int(value)


def _validate_sequence_elements(
    values: object,
    *,
    field_name: str,
    element_type: type,
    allow_none: bool = False,
) -> None:
    if values is None:
        return
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ValidationError(f"{field_name} must be a sequence.")
    for index, value in enumerate(values):
        if allow_none and value is None:
            continue
        if not isinstance(value, element_type):
            raise ValidationError(
                f"{field_name}[{index}] must be "
                f"{element_type.__name__}{' or None' if allow_none else ''}."
            )


def _validate_nested_sequence_elements(
    rows: object,
    *,
    field_name: str,
    element_type: type,
) -> None:
    if rows is None:
        return
    if isinstance(rows, (str, bytes)) or not isinstance(rows, Sequence):
        raise ValidationError(f"{field_name} must be a nested sequence.")
    for row_index, row in enumerate(rows):
        _validate_sequence_elements(
            row,
            field_name=f"{field_name}[{row_index}]",
            element_type=element_type,
            allow_none=True,
        )


@dataclass(frozen=True)
class CircularTrackOptions:
    """Circular track layout options."""

    circular_track_slots: Sequence[str | CircularTrackSlot] | None = None
    circular_track_axis_index: int | None = None
    center_reserved_radius: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "circular_track_axis_index",
            _validate_track_configuration(
                self.circular_track_slots,
                axis_index=self.circular_track_axis_index,
                mode="circular",
                expected_type=CircularTrackSlot,
                slots_field_name="circular_track_slots",
                axis_field_name="circular_track_axis_index",
            ),
        )
        object.__setattr__(
            self,
            "center_reserved_radius",
            _validate_center_reserved_radius(
                self.center_reserved_radius,
                field_name="center_reserved_radius",
            ),
        )


@dataclass(frozen=True)
class LinearTrackOptions:
    """Linear track layout options."""

    linear_track_slots: Sequence[str | LinearTrackSlot] | None = None
    linear_track_axis_index: int | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "linear_track_axis_index",
            _validate_track_configuration(
                self.linear_track_slots,
                axis_index=self.linear_track_axis_index,
                mode="linear",
                expected_type=LinearTrackSlot,
                slots_field_name="linear_track_slots",
                axis_field_name="linear_track_axis_index",
            ),
        )


@dataclass(frozen=True)
class TrackOptions:
    """Internal compatibility options for low-level assembly helpers."""

    circular_track_slots: Sequence[str | CircularTrackSlot] | None = None
    circular_track_axis_index: int | None = None
    linear_track_slots: Sequence[str | LinearTrackSlot] | None = None
    linear_track_axis_index: int | None = None
    center_reserved_radius: float | None = None

    def __post_init__(self) -> None:
        circular_axis_index = _validate_track_configuration(
            self.circular_track_slots,
            axis_index=self.circular_track_axis_index,
            mode="circular",
            expected_type=CircularTrackSlot,
            slots_field_name="circular_track_slots",
            axis_field_name="circular_track_axis_index",
        )
        linear_axis_index = _validate_track_configuration(
            self.linear_track_slots,
            axis_index=self.linear_track_axis_index,
            mode="linear",
            expected_type=LinearTrackSlot,
            slots_field_name="linear_track_slots",
            axis_field_name="linear_track_axis_index",
        )
        object.__setattr__(
            self,
            "circular_track_axis_index",
            circular_axis_index,
        )
        object.__setattr__(
            self,
            "linear_track_axis_index",
            linear_axis_index,
        )
        object.__setattr__(
            self,
            "center_reserved_radius",
            _validate_center_reserved_radius(
                self.center_reserved_radius,
                field_name="center_reserved_radius",
            ),
        )


@dataclass(frozen=True)
class CircularOutputOptions:
    """Circular legend and title placement options."""

    legend: str = "right"
    plot_title_position: Literal["none", "top", "bottom"] | None = None

    def __post_init__(self) -> None:
        if self.plot_title_position not in {None, "none", "top", "bottom"}:
            raise ValidationError(
                "Circular plot_title_position must be one of: none, top, bottom."
            )


@dataclass(frozen=True)
class LinearOutputOptions:
    """Linear legend and title placement options."""

    legend: str = "right"
    plot_title_position: Literal["center", "top", "bottom"] | None = None

    def __post_init__(self) -> None:
        if self.plot_title_position not in {None, "center", "top", "bottom"}:
            raise ValidationError(
                "Linear plot_title_position must be one of: center, top, bottom."
            )


@dataclass(frozen=True)
class OutputOptions:
    """Internal compatibility output options used by low-level assemblers."""

    legend: str = "right"
    plot_title_position: Literal["none", "center", "top", "bottom"] | None = None

    def __post_init__(self) -> None:
        if self.plot_title_position not in {
            None,
            "none",
            "center",
            "top",
            "bottom",
        }:
            raise ValidationError(
                "plot_title_position must be one of: none, center, top, bottom."
            )


@dataclass(frozen=True)
class CircularMultiRecordOptions:
    """Layout values used only by circular multi-record canvases."""

    multi_record_size_mode: Literal["linear", "auto", "equal"] = "auto"
    multi_record_min_radius_ratio: float = 0.55
    multi_record_column_gap_ratio: float = 0.10
    multi_record_row_gap_ratio: float = 0.05
    multi_record_positions: Sequence[str] | None = None

    def __post_init__(self) -> None:
        if self.multi_record_size_mode not in {"linear", "auto", "equal"}:
            raise ValidationError(
                "multi_record_size_mode must be one of: auto, linear, equal."
            )
        ratio_values = (
            self.multi_record_min_radius_ratio,
            self.multi_record_column_gap_ratio,
            self.multi_record_row_gap_ratio,
        )
        if any(
            isinstance(value, bool) or not isinstance(value, Real)
            for value in ratio_values
        ):
            raise ValidationError(
                "Circular multi-record ratios must be finite numbers."
            )
        try:
            min_radius_ratio = float(self.multi_record_min_radius_ratio)
            column_gap_ratio = float(self.multi_record_column_gap_ratio)
            row_gap_ratio = float(self.multi_record_row_gap_ratio)
        except (OverflowError, TypeError, ValueError) as exc:
            raise ValidationError(
                "Circular multi-record ratios must be finite numbers."
            ) from exc
        if (
            not math.isfinite(min_radius_ratio)
            or min_radius_ratio <= 0
            or min_radius_ratio > 1
        ):
            raise ValidationError(
                "multi_record_min_radius_ratio must be a finite number in (0, 1]."
            )
        for field_name, value in (
            ("multi_record_column_gap_ratio", column_gap_ratio),
            ("multi_record_row_gap_ratio", row_gap_ratio),
        ):
            if not math.isfinite(value) or value < 0:
                raise ValidationError(
                    f"{field_name} must be a finite number >= 0."
                )
        object.__setattr__(self, "multi_record_min_radius_ratio", min_radius_ratio)
        object.__setattr__(self, "multi_record_column_gap_ratio", column_gap_ratio)
        object.__setattr__(self, "multi_record_row_gap_ratio", row_gap_ratio)
        if self.multi_record_positions is not None:
            if isinstance(self.multi_record_positions, (str, bytes)) or not isinstance(
                self.multi_record_positions,
                Sequence,
            ):
                raise ValidationError(
                    "multi_record_positions must be a sequence of strings or None."
                )
            positions = tuple(self.multi_record_positions)
            if not all(isinstance(item, str) and item.strip() for item in positions):
                raise ValidationError(
                    "multi_record_positions must contain non-empty strings."
                )
            object.__setattr__(self, "multi_record_positions", positions)


@dataclass(frozen=True)
class LinearMultiRecordOptions:
    """Layout values used only by Linear multi-record rows."""

    record_gap_px: float = 24.0
    multi_record_positions: Sequence[str] | None = None

    def __post_init__(self) -> None:
        if isinstance(self.record_gap_px, bool) or not isinstance(
            self.record_gap_px,
            Real,
        ):
            raise ValidationError("record_gap_px must be a finite non-negative number.")
        try:
            value = float(self.record_gap_px)
        except (OverflowError, TypeError, ValueError) as exc:
            raise ValidationError(
                "record_gap_px must be a finite non-negative number."
            ) from exc
        if not math.isfinite(value) or value < 0:
            raise ValidationError("record_gap_px must be a finite non-negative number.")
        object.__setattr__(self, "record_gap_px", value)
        if self.multi_record_positions is not None:
            if isinstance(self.multi_record_positions, (str, bytes)) or not isinstance(
                self.multi_record_positions,
                Sequence,
            ):
                raise ValidationError("multi_record_positions must be a sequence of strings or None.")
            positions = tuple(self.multi_record_positions)
            if not all(isinstance(item, str) and item.strip() for item in positions):
                raise ValidationError("multi_record_positions must contain non-empty strings.")
            object.__setattr__(self, "multi_record_positions", positions)


@dataclass(frozen=True)
class DiagramOptions:
    """Internal compatibility options for low-level assembly helpers."""

    config: GbdrawConfig | dict | None = None
    config_overrides: Mapping[str, object] | None = None
    colors: ColorOptions | None = None
    tracks: TrackOptions | None = None
    annotations: AnnotationOptions | None = None
    output: OutputOptions | None = None
    selected_features_set: Sequence[str] | None = None
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
    linear_comparisons: Sequence[LinearComparison] | None = None
    protein_comparisons: Sequence[DataFrame] | None = None
    orthogroups: OrthogroupResult | None = None
    protein_blastp_mode: Literal["none", "pairwise", "orthogroup", "collinear"] = "none"
    protein_comparison_pairs: Sequence[tuple[int, int]] | None = None
    pairwise_match_style: Literal["ribbon", "curve"] = "ribbon"
    collinearity_blocks: CollinearityResult | Sequence[CollinearityBlock] | None = None
    collinearity_params: LosslessCollinearityParameters | None = None
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
    evalue: float | None = None
    bitscore: float | None = None
    identity: float | None = None
    alignment_length: int | None = None

    def __post_init__(self) -> None:
        nested_types = (
            ("colors", self.colors, ColorOptions),
            ("tracks", self.tracks, TrackOptions),
            ("annotations", self.annotations, AnnotationOptions),
            ("output", self.output, OutputOptions),
        )
        for field_name, value, expected_type in nested_types:
            if value is not None and not isinstance(value, expected_type):
                raise ValidationError(
                    f"{field_name} must be {expected_type.__name__} or None."
                )
        if self.config is not None and not isinstance(
            self.config,
            (GbdrawConfig, dict),
        ):
            raise ValidationError("config must be GbdrawConfig, dict, or None.")
        for field_name, value in (
            ("config_overrides", self.config_overrides),
            ("feature_shapes", self.feature_shapes),
        ):
            if value is not None and not isinstance(value, Mapping):
                raise ValidationError(f"{field_name} must be a mapping or None.")
        if self.selected_features_set is not None:
            if isinstance(self.selected_features_set, (str, bytes)) or not isinstance(
                self.selected_features_set,
                Sequence,
            ):
                raise ValidationError(
                    "selected_features_set must be a sequence of feature names or None."
                )
            if not all(
                isinstance(feature_name, str) and feature_name.strip()
                for feature_name in self.selected_features_set
            ):
                raise ValidationError(
                    "selected_features_set must contain non-empty strings."
                )
        if self.depth_table is not None and not isinstance(
            self.depth_table,
            DataFrame,
        ):
            raise ValidationError("depth_table must be DataFrame or None.")
        if self.depth_file is not None and not isinstance(self.depth_file, str):
            raise ValidationError("depth_file must be a string or None.")
        _validate_sequence_elements(
            self.depth_tables,
            field_name="depth_tables",
            element_type=DataFrame,
        )
        _validate_sequence_elements(
            self.depth_files,
            field_name="depth_files",
            element_type=str,
        )
        _validate_nested_sequence_elements(
            self.depth_track_tables,
            field_name="depth_track_tables",
            element_type=DataFrame,
        )
        _validate_nested_sequence_elements(
            self.depth_track_files,
            field_name="depth_track_files",
            element_type=str,
        )
        thresholds = ComparisonThresholds(
            evalue=0.0 if self.evalue is None else self.evalue,
            bitscore=0.0 if self.bitscore is None else self.bitscore,
            identity=0.0 if self.identity is None else self.identity,
            alignment_length=(
                0 if self.alignment_length is None else self.alignment_length
            ),
        )
        if self.evalue is not None:
            object.__setattr__(self, "evalue", thresholds.evalue)
        if self.bitscore is not None:
            object.__setattr__(self, "bitscore", thresholds.bitscore)
        if self.identity is not None:
            object.__setattr__(self, "identity", thresholds.identity)
        if self.alignment_length is not None:
            object.__setattr__(
                self,
                "alignment_length",
                thresholds.alignment_length,
            )


@dataclass(frozen=True)
class _ModeDiagramOptions:
    """Fields shared by the mode-specific typed request options."""

    config: GbdrawConfig | dict | None = None
    config_overrides: Mapping[str, object] | None = None
    colors: ColorOptions | None = None
    annotations: AnnotationOptions | None = None
    selected_features_set: Sequence[str] | None = None
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
    depth_track_large_tick_intervals: Sequence[float | str | None] | None = None
    depth_track_small_tick_intervals: Sequence[float | str | None] | None = None
    depth_track_tick_font_sizes: Sequence[float | str | None] | None = None
    plot_title: str | None = None
    plot_title_font_size: float | None = None
    evalue: float | None = None
    bitscore: float | None = None
    identity: float | None = None
    alignment_length: int | None = None

    def __post_init__(self) -> None:
        # Keep this conversion private while low-level builders still consume the
        # old shared object. It also centralizes validation during the migration.
        common_values = {
            item.name: getattr(self, item.name)
            for item in fields(_ModeDiagramOptions)
        }
        validated = DiagramOptions(**common_values)
        for item in fields(_ModeDiagramOptions):
            object.__setattr__(self, item.name, getattr(validated, item.name))


@dataclass(frozen=True)
class CircularDiagramOptions(_ModeDiagramOptions):
    """Options accepted by a Circular typed request."""

    tracks: CircularTrackOptions | None = None
    output: CircularOutputOptions | None = None
    conservation_blast_files: Sequence[str] | None = None
    conservation_fasta_files: Sequence[str | None] | None = None
    conservation_dataframes: Sequence[DataFrame] | None = None
    conservation_reference: Literal["query", "subject", "auto"] = "auto"
    conservation_labels: Sequence[str] | None = None
    conservation_colors: Sequence[str] | None = None
    conservation_ring_width: float | None = None
    conservation_ring_gap: float | None = None
    keep_full_definition_with_plot_title: bool = False
    species: str | None = None
    strain: str | None = None

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.tracks is not None and not isinstance(
            self.tracks,
            CircularTrackOptions,
        ):
            raise ValidationError(
                "tracks must be CircularTrackOptions or None."
            )
        if self.output is not None and not isinstance(
            self.output,
            CircularOutputOptions,
        ):
            raise ValidationError(
                "output must be CircularOutputOptions or None."
            )
        _validate_sequence_elements(
            self.conservation_dataframes,
            field_name="conservation_dataframes",
            element_type=DataFrame,
        )
        _validate_sequence_elements(
            self.conservation_blast_files,
            field_name="conservation_blast_files",
            element_type=str,
        )
        _validate_sequence_elements(
            self.conservation_fasta_files,
            field_name="conservation_fasta_files",
            element_type=str,
            allow_none=True,
        )
        _validate_sequence_elements(
            self.conservation_labels,
            field_name="conservation_labels",
            element_type=str,
        )
        _validate_sequence_elements(
            self.conservation_colors,
            field_name="conservation_colors",
            element_type=str,
        )
        object.__setattr__(
            self,
            "conservation_reference",
            normalize_conservation_reference(self.conservation_reference),
        )
        for field_name in ("conservation_ring_width", "conservation_ring_gap"):
            object.__setattr__(
                self,
                field_name,
                _validate_positive_real(
                    getattr(self, field_name),
                    field_name=field_name,
                ),
            )


@dataclass(frozen=True)
class LinearDiagramOptions(_ModeDiagramOptions):
    """Options accepted by a Linear typed request."""

    tracks: LinearTrackOptions | None = None
    output: LinearOutputOptions | None = None
    depth_track_heights: Sequence[float | str | None] | None = None
    blast_files: Sequence[str] | None = None
    linear_comparisons: Sequence[LinearComparison] | None = None
    protein_comparisons: Sequence[DataFrame] | None = None
    orthogroups: OrthogroupResult | None = None
    protein_blastp_mode: Literal[
        "none",
        "pairwise",
        "orthogroup",
        "collinear",
    ] = "none"
    protein_comparison_pairs: Sequence[tuple[int, int]] | None = None
    pairwise_match_style: Literal["ribbon", "curve"] = "ribbon"
    collinearity_blocks: (
        CollinearityResult | Sequence[CollinearityBlock] | None
    ) = None
    collinearity_params: LosslessCollinearityParameters | None = None
    collinearity_unit_mode: CollinearityUnitMode | str = "auto"
    collinearity_anchor_mode: CollinearityAnchorMode | str = "rbh"
    collinearity_search_scope: CollinearitySearchScope | str = "adjacent"
    collinearity_color_mode: CollinearityColorMode | str = "orientation"
    losatp_bin: str = "losat"
    ncbi_blastp_bin: str | None = None
    losatp_threads: int | None = None
    protein_blastp_max_hits: int = 5
    protein_blastp_candidate_limit: int | None = None
    orthogroup_membership_mode: Literal[
        "anchor_core_v1"
    ] | str = "anchor_core_v1"
    orthogroup_member_max_hits: int = 5
    collinear_max_paralog_links_per_orthogroup: int = 2
    align_orthogroup_feature: str | None = None

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.tracks is not None and not isinstance(
            self.tracks,
            LinearTrackOptions,
        ):
            raise ValidationError(
                "tracks must be LinearTrackOptions or None."
            )
        if self.output is not None and not isinstance(
            self.output,
            LinearOutputOptions,
        ):
            raise ValidationError(
                "output must be LinearOutputOptions or None."
            )
        _validate_sequence_elements(
            self.linear_comparisons,
            field_name="linear_comparisons",
            element_type=LinearComparison,
        )
        _validate_sequence_elements(
            self.protein_comparisons,
            field_name="protein_comparisons",
            element_type=DataFrame,
        )
        if self.collinearity_params is not None and not isinstance(
            self.collinearity_params,
            LosslessCollinearityParameters,
        ):
            raise ValidationError(
                "collinearity_params has an unsupported type."
            )
        if self.collinearity_params is not None:
            self.collinearity_params.validate()
        pairwise_match_style = str(self.pairwise_match_style).strip().lower()
        if pairwise_match_style not in {"ribbon", "curve"}:
            raise ValidationError(
                "pairwise_match_style must be one of: ribbon, curve."
            )
        object.__setattr__(self, "pairwise_match_style", pairwise_match_style)
        object.__setattr__(
            self,
            "protein_blastp_mode",
            normalize_protein_blastp_mode(self.protein_blastp_mode),
        )
        object.__setattr__(
            self,
            "collinearity_unit_mode",
            normalize_collinearity_unit_mode(str(self.collinearity_unit_mode)),
        )
        object.__setattr__(
            self,
            "collinearity_anchor_mode",
            normalize_collinearity_anchor_mode(str(self.collinearity_anchor_mode)),
        )
        object.__setattr__(
            self,
            "collinearity_search_scope",
            normalize_collinearity_search_scope(
                str(self.collinearity_search_scope)
            ),
        )
        object.__setattr__(
            self,
            "collinearity_color_mode",
            normalize_collinearity_color_mode(str(self.collinearity_color_mode)),
        )
        object.__setattr__(
            self,
            "orthogroup_membership_mode",
            normalize_orthogroup_membership_mode(
                str(self.orthogroup_membership_mode)
            ),
        )
        for field_name in (
            "protein_blastp_max_hits",
            "orthogroup_member_max_hits",
            "collinear_max_paralog_links_per_orthogroup",
        ):
            object.__setattr__(
                self,
                field_name,
                _validate_positive_int(
                    getattr(self, field_name),
                    field_name=field_name,
                ),
            )
        for field_name in ("losatp_threads", "protein_blastp_candidate_limit"):
            object.__setattr__(
                self,
                field_name,
                _validate_positive_int(
                    getattr(self, field_name),
                    field_name=field_name,
                    allow_none=True,
                ),
            )


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
    "linear_comparisons",
    "protein_comparisons",
    "orthogroups",
    "protein_blastp_mode",
    "protein_comparison_pairs",
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
    tracks = options.tracks
    if tracks is not None:
        nested_incompatible = (
            (
                "circular_track_slots",
                "circular_track_axis_index",
                "center_reserved_radius",
            )
            if mode == "linear"
            else ("linear_track_slots", "linear_track_axis_index")
        )
        nested_non_default = [
            name for name in nested_incompatible if getattr(tracks, name) is not None
        ]
        if nested_non_default:
            raise ValidationError(
                f"{mode} builder does not support TrackOptions fields: "
                f"{', '.join(nested_non_default)}"
            )
    if options.output is not None:
        title_position = options.output.plot_title_position
        if mode == "linear" and title_position == "none":
            raise ValidationError(
                "Linear plot_title_position must be one of: center, top, bottom."
            )
        if mode != "linear" and title_position == "center":
            raise ValidationError(
                "Circular plot_title_position must be one of: none, top, bottom."
            )


def _resolve_options_for_mode(
    options: DiagramOptions | CircularDiagramOptions | LinearDiagramOptions,
    *,
    mode: DiagramMode,
) -> DiagramOptions | CircularDiagramOptions | LinearDiagramOptions:
    profile = get_mode_profile(mode)
    thresholds = ComparisonThresholds(
        evalue=profile.comparison.evalue if options.evalue is None else options.evalue,
        bitscore=(
            profile.comparison.bitscore
            if options.bitscore is None
            else options.bitscore
        ),
        identity=(
            profile.comparison.identity
            if options.identity is None
            else options.identity
        ),
        alignment_length=(
            profile.comparison.alignment_length
            if options.alignment_length is None
            else options.alignment_length
        ),
    )

    config_overrides = options.config_overrides
    if options.config is None:
        config_overrides = resolve_mode_profile_overrides(mode, config_overrides)

    return replace(
        options,
        config_overrides=config_overrides,
        selected_features_set=(
            profile.feature_types
            if options.selected_features_set is None
            else options.selected_features_set
        ),
        evalue=thresholds.evalue,
        bitscore=thresholds.bitscore,
        identity=thresholds.identity,
        alignment_length=thresholds.alignment_length,
    )


def resolve_diagram_options_for_mode(
    options: DiagramOptions,
    *,
    mode: DiagramMode,
) -> DiagramOptions:
    """Resolve an internal shared option object through one mode profile."""

    if not isinstance(options, DiagramOptions):
        raise ValidationError("options must be DiagramOptions.")
    return _resolve_options_for_mode(options, mode=mode)


def resolve_circular_diagram_options(
    options: CircularDiagramOptions,
) -> CircularDiagramOptions:
    """Resolve omitted Circular request values through the Circular profile."""

    if not isinstance(options, CircularDiagramOptions):
        raise ValidationError("options must be CircularDiagramOptions.")
    return _resolve_options_for_mode(options, mode="circular")


def resolve_linear_diagram_options(
    options: LinearDiagramOptions,
) -> LinearDiagramOptions:
    """Resolve omitted Linear request values through the Linear profile."""

    if not isinstance(options, LinearDiagramOptions):
        raise ValidationError("options must be LinearDiagramOptions.")
    return _resolve_options_for_mode(options, mode="linear")


def _legacy_diagram_options(
    options: CircularDiagramOptions | LinearDiagramOptions,
) -> DiagramOptions:
    """Convert a strict request option object for an internal legacy builder."""

    values = {item.name: getattr(options, item.name) for item in fields(options)}
    values.pop("conservation_fasta_files", None)
    tracks = values.get("tracks")
    if isinstance(tracks, CircularTrackOptions):
        values["tracks"] = TrackOptions(
            circular_track_slots=tracks.circular_track_slots,
            circular_track_axis_index=tracks.circular_track_axis_index,
            center_reserved_radius=tracks.center_reserved_radius,
        )
    elif isinstance(tracks, LinearTrackOptions):
        values["tracks"] = TrackOptions(
            linear_track_slots=tracks.linear_track_slots,
            linear_track_axis_index=tracks.linear_track_axis_index,
        )
    output = values.get("output")
    if isinstance(output, (CircularOutputOptions, LinearOutputOptions)):
        values["output"] = OutputOptions(
            legend=output.legend,
            plot_title_position=output.plot_title_position,
        )
    return DiagramOptions(**values)


__all__ = [
    "CircularDiagramOptions",
    "CircularMultiRecordOptions",
    "CircularOutputOptions",
    "CircularTrackOptions",
    "LinearMultiRecordOptions",
    "LinearDiagramOptions",
    "LinearOutputOptions",
    "LinearTrackOptions",
    "AnnotationOptions",
    "ColorOptions",
    "resolve_circular_diagram_options",
    "resolve_linear_diagram_options",
]
