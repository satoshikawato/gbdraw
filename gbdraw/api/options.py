"""Mode-specific typed option bundles used by public requests and builders."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import math
from numbers import Integral, Real
from pathlib import Path
from types import MappingProxyType
from typing import Literal, Mapping, Sequence, cast

from pandas import DataFrame

from gbdraw.features.overrides import FeatureOverride, normalize_feature_override_inputs
from gbdraw.features.placement import FeaturePlacementOverride, normalize_feature_placements

from gbdraw.analysis.collinearity import (
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
from gbdraw.analysis.collinearity_units import (
    CollinearityUnitMode,
    normalize_collinearity_unit_mode,
)
from gbdraw.analysis.conservation import (
    ConservationSearchResult,
    normalize_conservation_reference,
)
from gbdraw.analysis.protein_colinearity import (
    OrthogroupResult,
    OrthogroupGraphResult,
    normalize_orthogroup_membership_mode,
)
from gbdraw.comparisons.losat_runtime import AUTOMATIC_LOSAT_BIN
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.models.objects import (
    normalize_pairwise_match_style,
)
from gbdraw.config.modify import validate_config_overrides
from gbdraw.exceptions import ValidationError
from gbdraw.features.shapes import (
    normalize_feature_shape_overrides,
)
from gbdraw.linear_comparison import LinearComparison
from gbdraw.mode_profiles import (
    ComparisonThresholds,
    DiagramMode,
    get_mode_profile,
    resolve_mode_profile_overrides,
    validate_dinucleotide,
)
from gbdraw.tracks import (
    CircularTrackSlot,
    LinearTrackSlot,
    normalize_circular_track_slots_with_axis,
    normalize_linear_track_slots_with_axis,
    parse_circular_track_slots,
    parse_linear_track_slots,
)
from gbdraw.annotations import AnnotationOptions


# Paths not listed here are shared by both renderers.
_MODE_CONFIG_OVERRIDE_PREFIXES: dict[DiagramMode, tuple[str, ...]] = {
    "circular": (
        "canvas.circular",
        "labels.circular",
        "labels.length_threshold.circular",
        "labels.font_size.short",
        "labels.font_size.long",
        "labels.radius_factor",
        "labels.inner_radius_factor",
        "labels.arc_x_radius_factor",
        "labels.arc_y_radius_factor",
        "labels.arc_center_x",
        "labels.arc_angle",
        "labels.inner_arc_x_radius_factor",
        "labels.inner_arc_y_radius_factor",
        "labels.inner_arc_center_x",
        "labels.inner_arc_angle",
        "labels.unified_adjustment",
        "labels.spacing.circular",
        "objects.axis.circular",
        "objects.conservation",
        "objects.definition.circular",
        "objects.ticks",
    ),
    "linear": (
        "canvas.linear",
        "labels.linear",
        "labels.length_threshold.linear",
        "labels.font_size.linear",
        "labels.spacing.linear",
        "objects.axis.linear",
        "objects.blast_match",
        "objects.definition.linear",
    ),
}


def _validate_mode_config_overrides(
    overrides: Mapping[str, object] | None,
    *,
    mode: DiagramMode,
) -> None:
    other_mode: DiagramMode = "linear" if mode == "circular" else "circular"
    other_prefixes = _MODE_CONFIG_OVERRIDE_PREFIXES[other_mode]
    wrong_paths = sorted(
        path
        for path in (overrides or {})
        if any(
            path == prefix or path.startswith(f"{prefix}.")
            for prefix in other_prefixes
        )
    )
    if not wrong_paths:
        return
    mode_name = mode.title()
    other_name = other_mode.title()
    target = (
        f"{other_name} label settings"
        if all(path.startswith(f"labels.{other_mode}.") for path in wrong_paths)
        else f"{other_name} settings"
    )
    raise ValidationError(
        f"{mode_name} config overrides cannot target {target}: "
        + ", ".join(wrong_paths)
        + "."
    )


@dataclass(frozen=True)
class ColorOptions:
    """Color table and palette inputs."""

    color_table: DataFrame | None = None
    color_table_file: str | None = None
    default_colors: DataFrame | None = None
    default_colors_palette: str = "default"
    default_colors_file: str | None = None


_DepthTrackSource = str | Path | DataFrame


@dataclass(frozen=True)
class DepthTrackInput:
    """Data and styling for one logical depth track in a typed request."""

    source: _DepthTrackSource | Sequence[_DepthTrackSource | None]
    label: str | None = None
    color: str | None = None
    height: float | None = None
    large_tick_interval: float | None = None
    small_tick_interval: float | None = None
    tick_font_size: float | None = None

    def __post_init__(self) -> None:
        source = self.source
        if isinstance(source, DataFrame):
            pass
        elif isinstance(source, (str, Path)):
            if not str(source).strip():
                raise ValidationError(
                    "DepthTrackInput.source path must not be empty."
                )
        elif isinstance(source, Sequence) and not isinstance(source, (str, bytes)):
            sources = tuple(source)
            if not sources:
                raise ValidationError(
                    "DepthTrackInput.source must include at least one source."
                )
            for index, item in enumerate(sources):
                if item is None or isinstance(item, DataFrame):
                    continue
                if isinstance(item, (str, Path)) and str(item).strip():
                    continue
                raise ValidationError(
                    "DepthTrackInput.source"
                    f"[{index}] must be a path, DataFrame, or None."
                )
            if not any(item is not None for item in sources):
                raise ValidationError(
                    "DepthTrackInput.source must include at least one non-None source."
                )
            object.__setattr__(self, "source", sources)
        else:
            raise ValidationError(
                "DepthTrackInput.source must be a path, DataFrame, or a sequence "
                "of per-record sources."
            )

        for field_name in ("label", "color"):
            value = getattr(self, field_name)
            if value is None:
                continue
            if not isinstance(value, str) or not value.strip():
                raise ValidationError(
                    f"DepthTrackInput.{field_name} must be a non-empty string or None."
                )
            object.__setattr__(self, field_name, value.strip())

        for field_name in (
            "height",
            "large_tick_interval",
            "small_tick_interval",
            "tick_font_size",
        ):
            object.__setattr__(
                self,
                field_name,
                _validate_positive_real(
                    getattr(self, field_name),
                    field_name=f"DepthTrackInput.{field_name}",
                ),
            )


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


def _invalid_input(field_name: str, reason: str) -> dict[str, str]:
    return {
        "code": "INPUT_INVALID",
        "field": field_name.rsplit(".", 1)[-1],
        "reason": reason,
    }


def _validate_positive_real(value: object, *, field_name: str) -> float | None:
    if value is None:
        return None
    normalized = math.nan
    if not isinstance(value, bool) and isinstance(value, Real):
        try:
            normalized = float(value)
        except (OverflowError, TypeError, ValueError):
            normalized = math.nan
    if not math.isfinite(normalized) or normalized <= 0:
        raise ValidationError(
            f"{field_name} must be a finite number > 0 or None.",
            diagnostic=_invalid_input(field_name, "POSITIVE_OR_AUTO"),
        )
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
        raise ValidationError(
            f"{field_name} must be a positive integer{suffix}.",
            diagnostic=_invalid_input(
                field_name,
                "POSITIVE_INTEGER_OR_AUTO" if allow_none else "POSITIVE_INTEGER",
            ),
        )
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
class CircularRequestTrackOptions:
    """Circular track layout options for typed requests."""

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
class LinearRequestTrackOptions:
    """Linear track layout options for typed requests."""

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


# Compatibility aliases for the original typed API names. The package-root
# classes with these names are separate beginner-facing option bundles.
CircularTrackOptions = CircularRequestTrackOptions
LinearTrackOptions = LinearRequestTrackOptions


@dataclass(frozen=True)
class CircularOutputOptions:
    """Circular legend and title placement options."""

    legend: str = "right"
    plot_title_position: Literal["none", "top", "bottom"] | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.legend, str) or self.legend not in {
            "left",
            "right",
            "top",
            "bottom",
            "upper_left",
            "upper_right",
            "lower_left",
            "lower_right",
            "none",
        }:
            raise ValidationError(
                "Circular legend must be one of: left, right, top, bottom, "
                "upper_left, upper_right, lower_left, lower_right, none."
            )
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
        if not isinstance(self.legend, str) or self.legend not in {
            "left",
            "right",
            "top",
            "bottom",
            "none",
        }:
            raise ValidationError(
                "Linear legend must be one of: left, right, top, bottom, none."
            )
        if self.plot_title_position not in {None, "center", "top", "bottom"}:
            raise ValidationError(
                "Linear plot_title_position must be one of: center, top, bottom."
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
                "multi_record_min_radius_ratio must be a finite number in (0, 1].",
                diagnostic=_invalid_input("multi_record_min_radius_ratio", "POSITIVE_UNIT_INTERVAL"),
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
class LinearRecordTranslation:
    """Persistent base translation for one displayed Linear record."""

    record_key: str
    x: float = 0.0
    y: float = 0.0

    def __post_init__(self) -> None:
        if (
            not isinstance(self.record_key, str)
            or not self.record_key.strip()
            or "\0" in self.record_key
        ):
            raise ValidationError(
                "record_key must be a non-empty string without NUL."
            )
        object.__setattr__(self, "record_key", self.record_key.strip())
        for field_name in ("x", "y"):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, Real):
                raise ValidationError(
                    f"{field_name} must be a finite number."
                )
            normalized = float(value)
            if not math.isfinite(normalized):
                raise ValidationError(
                    f"{field_name} must be a finite number."
                )
            object.__setattr__(self, field_name, normalized)


@dataclass(frozen=True)
class LinearMultiRecordOptions:
    """Layout values used only by Linear multi-record rows."""

    record_gap_px: float = 24.0
    multi_record_positions: Sequence[str] | None = None
    record_translations: Sequence[LinearRecordTranslation] = ()

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
        if isinstance(self.record_translations, (str, bytes)) or not isinstance(
            self.record_translations,
            Sequence,
        ):
            raise ValidationError(
                "record_translations must be a sequence of LinearRecordTranslation values."
            )
        translations = tuple(self.record_translations)
        if not all(
            isinstance(item, LinearRecordTranslation) for item in translations
        ):
            raise ValidationError(
                "record_translations must contain LinearRecordTranslation values."
            )
        keys = [item.record_key for item in translations]
        if len(set(keys)) != len(keys):
            raise ValidationError(
                "record_translations must not contain duplicate record keys."
            )
        object.__setattr__(self, "record_translations", translations)


LosatProgram = Literal["losatn", "tlosatx", "losatp"]
LosatpMode = Literal["similarity_groups", "collinear", "pairwise"]
LosatnTask = Literal["megablast", "blastn", "dc-megablast"]
_LOSAT_PROGRAMS: tuple[str, ...] = ("losatn", "tlosatx", "losatp")
LOSATN_TASKS: tuple[str, ...] = ("megablast", "blastn", "dc-megablast")
DEFAULT_LOSATN_TASK = "megablast"

# The one table between a typed LOSATP display mode and its persisted
# ``generatedProteinComparison.mode`` spelling, which is also the protein
# analysis spelling (design D3, D6). ``"none"`` keeps the search settings of
# protein evidence that a request already carries (a resolved request); no
# search runs.
LOSATP_MODE_WIRE: Mapping[str, str] = MappingProxyType(
    {
        "similarity_groups": "orthogroup",
        "collinear": "collinear",
        "pairwise": "pairwise",
        "none": "none",
    }
)


@dataclass(frozen=True)
class LosatRuntimeOptions:
    """Executable choice and thread count for LOSAT searches.

    ``None`` executables select the runtime automatically. The NCBI BLAST+
    executable is the one for the selected program.
    """

    losat_executable: str | None = None
    ncbi_blast_executable: str | None = None
    threads: int | None = None

    def __post_init__(self) -> None:
        for name in ("losat_executable", "ncbi_blast_executable"):
            value = getattr(self, name)
            if value is None:
                continue
            if not isinstance(value, str) or not value.strip() or "\0" in value:
                raise ValidationError(
                    f"{name} must be a non-empty string or None.",
                    diagnostic={"code": "COMPARISON_INPUT"},
                )
            normalized: str | None = value.strip()
            if name == "losat_executable" and normalized == AUTOMATIC_LOSAT_BIN:
                normalized = None
            object.__setattr__(self, name, normalized)
        if (
            self.losat_executable is not None
            and self.ncbi_blast_executable is not None
        ):
            raise ValidationError(
                "Pass either losat_executable or ncbi_blast_executable, not both.",
                diagnostic={"code": "COMPARISON_INPUT"},
            )
        object.__setattr__(
            self,
            "threads",
            _validate_positive_int(
                self.threads,
                field_name="threads",
                allow_none=True,
            ),
        )


@dataclass(frozen=True)
class LosatSearchOptions:
    """One LOSAT comparison search for a diagram.

    ``losatp_mode`` is required for ``losatp``. ``losatn_task`` applies to
    ``losatn`` (default ``megablast``). ``record_gencodes`` applies to
    ``tlosatx``: empty uses the runtime default table (1) for every record, one
    value applies to every record input, otherwise one value (or ``None``) per
    record input. ``pairs`` lists explicit ``(query, subject)`` record indexes;
    ``None`` compares adjacent rows.
    """

    program: LosatProgram
    pairs: Sequence[tuple[int, int]] | None = None
    losatp_mode: LosatpMode | Literal["none"] | None = None
    losatn_task: LosatnTask | None = None
    record_gencodes: Sequence[int | None] = ()
    losatp_max_hits: int = 5
    losatp_max_target_seqs: int | None = None
    losatp_member_max_hits: int | None = None
    runtime: LosatRuntimeOptions = field(default_factory=LosatRuntimeOptions)

    def _option_program_error(self, field_name: str) -> ValidationError:
        return ValidationError(
            f"{field_name} does not apply to LOSAT program {self.program!r}.",
            diagnostic={
                "code": "COMPARISON_INPUT",
                "reason": "LOSAT_OPTION_PROGRAM",
                "field": field_name,
                "program": self.program,
            },
        )

    def __post_init__(self) -> None:
        if self.program not in _LOSAT_PROGRAMS:
            raise ValidationError(
                "program must be one of: " + ", ".join(_LOSAT_PROGRAMS) + ".",
                diagnostic={"code": "COMPARISON_INPUT"},
            )
        if self.program == "losatp":
            self._validate_losatp()
        else:
            self._validate_nucleotide()
        if self.pairs is not None:
            object.__setattr__(self, "pairs", _normalize_record_pairs(self.pairs))
            if self.program == "losatp" and self.losatp_mode != "pairwise":
                raise ValidationError(
                    "pairs requires losatp_mode 'pairwise'.",
                    diagnostic={"code": "COMPARISON_INPUT", "reason": "LOSAT_PLAN"},
                )
        if not isinstance(self.runtime, LosatRuntimeOptions):
            raise ValidationError(
                "runtime must be LosatRuntimeOptions.",
                diagnostic={"code": "COMPARISON_INPUT"},
            )

    def _validate_losatp(self) -> None:
        for name in ("losatn_task", "record_gencodes"):
            if getattr(self, name):
                raise self._option_program_error(name)
        object.__setattr__(self, "record_gencodes", ())
        if self.losatp_mode is None:
            raise ValidationError(
                "losatp_mode is required when program is 'losatp'.",
                diagnostic={"code": "COMPARISON_INPUT"},
            )
        if self.losatp_mode not in LOSATP_MODE_WIRE:
            raise ValidationError(
                "losatp_mode must be one of: similarity_groups, collinear, pairwise.",
                diagnostic={"code": "COMPARISON_INPUT"},
            )
        object.__setattr__(
            self,
            "losatp_max_hits",
            _validate_positive_int(self.losatp_max_hits, field_name="losatp_max_hits"),
        )
        for name in ("losatp_max_target_seqs", "losatp_member_max_hits"):
            object.__setattr__(
                self,
                name,
                _validate_positive_int(
                    getattr(self, name),
                    field_name=name,
                    allow_none=True,
                ),
            )

    def _validate_nucleotide(self) -> None:
        if self.losatp_mode is not None:
            raise self._option_program_error("losatp_mode")
        for name, default in (
            ("losatp_max_hits", 5),
            ("losatp_max_target_seqs", None),
            ("losatp_member_max_hits", None),
        ):
            if getattr(self, name) != default:
                raise self._option_program_error(name)
        if self.program == "losatn":
            if self.record_gencodes:
                raise self._option_program_error("record_gencodes")
            task = DEFAULT_LOSATN_TASK if self.losatn_task is None else self.losatn_task
            if task not in LOSATN_TASKS:
                raise ValidationError(
                    "losatn_task must be one of: " + ", ".join(LOSATN_TASKS) + ".",
                    diagnostic={"code": "COMPARISON_INPUT", "field": "losatn_task"},
                )
            object.__setattr__(self, "losatn_task", task)
            object.__setattr__(self, "record_gencodes", ())
            return
        if self.losatn_task is not None:
            raise self._option_program_error("losatn_task")
        gencodes = self.record_gencodes
        if isinstance(gencodes, (str, bytes)) or not isinstance(gencodes, Sequence):
            raise ValidationError(
                "record_gencodes must be a sequence of positive integers or None.",
                diagnostic={"code": "COMPARISON_INPUT", "field": "record_gencodes"},
            )
        object.__setattr__(
            self,
            "record_gencodes",
            tuple(
                None
                if value is None
                else _validate_positive_int(value, field_name="record_gencodes")
                for value in gencodes
            ),
        )


def _normalize_record_pairs(
    pairs: object,
) -> tuple[tuple[int, int], ...]:
    if isinstance(pairs, (str, bytes)) or not isinstance(pairs, Sequence):
        raise ValidationError("pairs must contain integer index pairs", diagnostic={"code": "COMPARISON_INPUT"})
    normalized: list[tuple[int, int]] = []
    for pair in pairs:
        if (
            isinstance(pair, (str, bytes))
            or not isinstance(pair, Sequence)
            or len(pair) != 2
            or any(isinstance(item, bool) or not isinstance(item, Integral) for item in pair)
            or any(int(item) < 0 for item in pair)
        ):
            raise ValidationError("pairs must contain integer index pairs", diagnostic={"code": "COMPARISON_INPUT"})
        normalized.append((int(pair[0]), int(pair[1])))
    return tuple(normalized)


def losatp_analysis_mode(search: LosatSearchOptions | None) -> str:
    """Return the protein-analysis mode a LOSAT search requests.

    ``"none"`` means that no LOSATP search runs.
    """

    if search is None or search.program != "losatp":
        return "none"
    return LOSATP_MODE_WIRE[str(search.losatp_mode)]


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
    plot_title: str | None = None
    plot_title_font_size: float | None = None
    evalue: float | None = None
    bitscore: float | None = None
    identity: float | None = None
    alignment_length: int | None = None
    feature_placements: tuple[FeaturePlacementOverride, ...] = field(default=(), kw_only=True)
    feature_placement_table: DataFrame | None = field(default=None, kw_only=True)
    feature_placement_table_file: str | Path | None = field(default=None, kw_only=True)
    feature_overrides: tuple[FeatureOverride, ...] = field(default=(), kw_only=True)
    feature_override_table: DataFrame | None = field(default=None, kw_only=True)
    feature_override_table_file: str | Path | None = field(default=None, kw_only=True)

    def __post_init__(self) -> None:
        placements = normalize_feature_placements(self.feature_placements)
        object.__setattr__(self, "feature_placements", placements)
        object.__setattr__(self, "feature_overrides", normalize_feature_override_inputs(
            self.feature_overrides,
            table=self.feature_override_table,
            table_file=self.feature_override_table_file,
        ))
        if sum((bool(placements), self.feature_placement_table is not None,
                self.feature_placement_table_file is not None)) > 1:
            raise ValidationError("Feature placement exact/table/file inputs are mutually exclusive.")
        if self.feature_placement_table is not None and not isinstance(self.feature_placement_table, DataFrame):
            raise ValidationError("feature_placement_table must be a DataFrame or None.")
        if self.feature_placement_table_file is not None and (
            not isinstance(self.feature_placement_table_file, (str, Path))
            or not str(self.feature_placement_table_file).strip()
        ):
            raise ValidationError("feature_placement_table_file must identify a file.")
        nested_types = (
            ("colors", self.colors, ColorOptions),
            ("annotations", self.annotations, AnnotationOptions),
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
        if isinstance(self.config, dict):
            try:
                object.__setattr__(
                    self,
                    "config",
                    GbdrawConfig.from_dict(self.config),
                )
            except ValidationError:
                raise
            except (KeyError, TypeError, ValueError) as exc:
                raise ValidationError(f"Invalid configuration: {exc}") from exc
        for field_name, value in (
            ("config_overrides", self.config_overrides),
            ("feature_shapes", self.feature_shapes),
        ):
            if value is not None and not isinstance(value, Mapping):
                raise ValidationError(f"{field_name} must be a mapping or None.")
        validate_config_overrides(self.config_overrides)
        if self.feature_shapes is not None:
            try:
                normalized_shapes = normalize_feature_shape_overrides(
                    self.feature_shapes
                )
            except ValueError as exc:
                raise ValidationError(str(exc)) from exc
            object.__setattr__(self, "feature_shapes", normalized_shapes)
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
        if self.depth_tracks is not None:
            if (
                isinstance(self.depth_tracks, (str, bytes))
                or not isinstance(self.depth_tracks, Sequence)
            ):
                raise ValidationError(
                    "depth_tracks must be a sequence of DepthTrackInput values."
                )
            depth_tracks = tuple(self.depth_tracks)
            if not depth_tracks:
                raise ValidationError(
                    "depth_tracks must include at least one DepthTrackInput."
                )
            if not all(isinstance(track, DepthTrackInput) for track in depth_tracks):
                raise ValidationError(
                    "depth_tracks must contain DepthTrackInput values."
                )
            compatibility_fields = (
                "depth_table",
                "depth_file",
                "depth_tables",
                "depth_files",
                "depth_track_tables",
                "depth_track_files",
                "depth_track_labels",
                "depth_track_colors",
                "depth_track_heights",
                "depth_track_large_tick_intervals",
                "depth_track_small_tick_intervals",
                "depth_track_tick_font_sizes",
            )
            mixed = [
                field_name
                for field_name in compatibility_fields
                if getattr(self, field_name, None) is not None
            ]
            if mixed:
                raise ValidationError(
                    "depth_tracks cannot be combined with compatibility depth "
                    "inputs: "
                    + ", ".join(mixed)
                    + "."
                )
            object.__setattr__(self, "depth_tracks", depth_tracks)
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
        object.__setattr__(
            self,
            "dinucleotide",
            validate_dinucleotide(self.dinucleotide),
        )
        for field_name in ("window", "step", "depth_window", "depth_step"):
            object.__setattr__(
                self,
                field_name,
                _validate_positive_int(
                    getattr(self, field_name),
                    field_name=field_name,
                    allow_none=True,
                ),
            )
        object.__setattr__(
            self,
            "plot_title_font_size",
            _validate_positive_real(
                self.plot_title_font_size,
                field_name="plot_title_font_size",
            ),
        )
        thresholds = ComparisonThresholds(
            evalue=0.0 if self.evalue is None else self.evalue,
            bitscore=0.0 if self.bitscore is None else self.bitscore,
            identity=0.0 if self.identity is None else self.identity,
            alignment_length=(
                0 if self.alignment_length is None else self.alignment_length
            ),
        )
        for field_name in ("evalue", "bitscore", "identity", "alignment_length"):
            if getattr(self, field_name) is not None:
                object.__setattr__(
                    self,
                    field_name,
                    getattr(thresholds, field_name),
                )


@dataclass(frozen=True)
class CircularDiagramOptions(_ModeDiagramOptions):
    """Options accepted by a Circular typed request."""

    tracks: CircularRequestTrackOptions | None = None
    output: CircularOutputOptions | None = None
    conservation_blast_files: Sequence[str] | None = None
    conservation_sequence_files: Sequence[str | None] | None = None
    conservation_dataframes: Sequence[DataFrame] | None = None
    conservation_reference: Literal["query", "subject", "auto"] = "auto"
    conservation_labels: Sequence[str] | None = None
    conservation_colors: Sequence[str] | None = None
    conservation_ring_width: float | None = None
    conservation_ring_gap: float | None = None
    keep_full_definition_with_plot_title: bool = False
    species: str | None = None
    strain: str | None = None
    conservation_table_file: str | None = None
    losat_search: LosatSearchOptions | None = None
    conservation_losat_gencodes: Sequence[int] | None = None
    conservation_search_results: Sequence[ConservationSearchResult] | None = None

    def __post_init__(self) -> None:
        super().__post_init__()
        for placement in self.feature_placements:
            placement.target.validate_mode("circular")
        if self.depth_tracks is not None and any(
            track.height is not None for track in self.depth_tracks
        ):
            raise ValidationError(
                "DepthTrackInput.height is available only for Linear diagrams."
            )
        _validate_mode_config_overrides(
            self.config_overrides,
            mode="circular",
        )
        if self.tracks is not None and not isinstance(
            self.tracks,
            CircularRequestTrackOptions,
        ):
            raise ValidationError(
                "tracks must be CircularRequestTrackOptions or None."
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
            self.conservation_sequence_files,
            field_name="conservation_sequence_files",
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
        if self.conservation_table_file is not None:
            if (
                not isinstance(self.conservation_table_file, str)
                or not self.conservation_table_file.strip()
            ):
                raise ValidationError(
                    "conservation_table_file must be a non-empty string or None."
                )
            if any(
                value is not None
                for value in (
                    self.conservation_blast_files,
                    self.conservation_sequence_files,
                    self.conservation_dataframes,
                    self.conservation_labels,
                    self.conservation_colors,
                    self.conservation_losat_gencodes,
                    self.conservation_search_results,
                )
            ):
                raise ValidationError(
                    "conservation_table_file cannot be combined with direct "
                    "Circular comparison inputs."
                )
            object.__setattr__(
                self,
                "conservation_table_file",
                self.conservation_table_file.strip(),
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

        _validate_sequence_elements(
            self.conservation_search_results,
            field_name="conservation_search_results",
            element_type=ConservationSearchResult,
        )
        if self.conservation_search_results is not None and any(
            value is not None
            for value in (self.conservation_blast_files, self.conservation_dataframes)
        ):
            raise ValidationError(
                "conservation_search_results cannot be combined with "
                "conservation_blast_files or conservation_dataframes.",
                diagnostic={"code": "COMPARISON_INPUT", "field": "conservation_search_results"},
            )
        self._validate_ring_losat()

    def _validate_ring_losat(self) -> None:
        """Circular LOSATN / TLOSATX ring intent (design 3.3, 3.4)."""

        search = self.losat_search
        gencodes = self.conservation_losat_gencodes
        if gencodes is not None:
            if isinstance(gencodes, (str, bytes)) or not isinstance(gencodes, Sequence):
                raise ValidationError(
                    "conservation_losat_gencodes must be a sequence of positive integers.",
                    diagnostic={"code": "COMPARISON_INPUT", "field": "conservation_losat_gencodes"},
                )
            object.__setattr__(
                self,
                "conservation_losat_gencodes",
                tuple(
                    _validate_positive_int(value, field_name="conservation_losat_gencodes")
                    for value in gencodes
                ),
            )
            gencodes = self.conservation_losat_gencodes
        if search is None:
            if gencodes is not None:
                raise ValidationError(
                    "conservation_losat_gencodes requires losat_search with program 'tlosatx'.",
                    diagnostic={
                        "code": "COMPARISON_INPUT",
                        "reason": "LOSAT_OPTION_PROGRAM",
                        "field": "conservation_losat_gencodes",
                    },
                )
            return
        if not isinstance(search, LosatSearchOptions):
            raise ValidationError(
                "losat_search must be LosatSearchOptions or None.",
                diagnostic={"code": "COMPARISON_INPUT", "field": "losat_search"},
            )
        if search.program not in {"losatn", "tlosatx"}:
            raise ValidationError(
                "Circular similarity rings run LOSATN or TLOSATX; "
                f"{search.program} is not available for rings.",
                diagnostic={
                    "code": "COMPARISON_INPUT",
                    "reason": "RING_LOSAT_PROGRAM",
                    "field": "program",
                    "program": search.program,
                },
            )
        if search.pairs is not None:
            raise ValidationError(
                "Circular rings compare each comparison genome with the displayed "
                "records; losat_search.pairs applies to Linear diagrams only.",
                diagnostic={"code": "COMPARISON_INPUT", "reason": "RING_LOSAT_INPUT", "field": "pairs"},
            )
        if len(tuple(search.record_gencodes)) > 1:
            raise ValidationError(
                "A Circular ring search takes one reference translation table "
                f"(losat_search.record_gencodes); got {len(tuple(search.record_gencodes))}.",
                diagnostic={"code": "COMPARISON_INPUT", "field": "record_gencodes"},
            )
        if gencodes is not None and search.program != "tlosatx":
            raise ValidationError(
                "conservation_losat_gencodes applies to TLOSATX rings only.",
                diagnostic={
                    "code": "COMPARISON_INPUT",
                    "reason": "LOSAT_OPTION_PROGRAM",
                    "field": "conservation_losat_gencodes",
                    "program": search.program,
                },
            )

        def ring_input_error(message: str, field_name: str) -> ValidationError:
            return ValidationError(
                message,
                diagnostic={
                    "code": "COMPARISON_INPUT",
                    "reason": "RING_LOSAT_INPUT",
                    "field": field_name,
                },
            )

        for field_name in (
            "conservation_blast_files",
            "conservation_dataframes",
            "conservation_search_results",
        ):
            if getattr(self, field_name) is not None:
                raise ring_input_error(
                    f"{field_name} cannot be combined with a ring LOSAT search; the "
                    "search produces the ring rows from conservation_sequence_files.",
                    field_name,
                )
        if self.conservation_reference == "query":
            raise ring_input_error(
                "A ring LOSAT search uses the displayed records as the subject; "
                "conservation_reference must be 'auto' or 'subject'.",
                "conservation_reference",
            )
        if self.conservation_table_file is not None:
            return
        sequences = tuple(self.conservation_sequence_files or ())
        if not sequences or any(not value for value in sequences):
            raise ring_input_error(
                "A ring LOSAT search needs one comparison sequence file per ring "
                "(conservation_sequence_files).",
                "conservation_sequence_files",
            )
        for field_name in ("conservation_labels", "conservation_colors"):
            values = getattr(self, field_name)
            if values is not None and len(tuple(values)) != len(sequences):
                raise ring_input_error(
                    f"{field_name} must give one value per comparison sequence file "
                    f"({len(sequences)}); got {len(tuple(values))}.",
                    field_name,
                )
        if gencodes is not None and len(gencodes) not in {1, len(sequences)}:
            raise ring_input_error(
                "conservation_losat_gencodes must give one translation table for all "
                f"rings or one per comparison sequence file ({len(sequences)}); got "
                f"{len(gencodes)}.",
                "conservation_losat_gencodes",
            )


@dataclass(frozen=True)
class LinearDiagramOptions(_ModeDiagramOptions):
    """Options accepted by a Linear typed request."""

    tracks: LinearRequestTrackOptions | None = None
    output: LinearOutputOptions | None = None
    depth_track_heights: Sequence[float | str | None] | None = None
    blast_files: Sequence[str] | None = None
    linear_comparisons: Sequence[LinearComparison] | None = None
    protein_comparisons: Sequence[DataFrame] | None = None
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None = None
    losat_search: LosatSearchOptions | None = None
    pairwise_match_style: Literal["ribbon", "curve"] = "ribbon"
    collinearity_blocks: (
        CollinearityResult | Sequence[CollinearityBlock] | None
    ) = None
    collinearity_params: LosslessCollinearityParameters | None = None
    collinearity_unit_mode: CollinearityUnitMode | str = "auto"
    collinearity_anchor_mode: CollinearityAnchorMode | str = "rbh"
    collinearity_search_scope: CollinearitySearchScope | str = "adjacent"
    collinearity_color_mode: CollinearityColorMode | str = "orientation"
    orthogroup_membership_mode: Literal[
        "anchor_core_v1"
    ] | str = "anchor_core_v1"
    collinear_infer_orthogroups: bool = True
    collinear_max_paralog_links_per_orthogroup: int = 2
    comparison_table_file: str | None = None

    def __post_init__(self) -> None:
        super().__post_init__()
        for placement in self.feature_placements:
            placement.target.validate_mode("linear")
        _validate_mode_config_overrides(
            self.config_overrides,
            mode="linear",
        )
        if self.tracks is not None and not isinstance(
            self.tracks,
            LinearRequestTrackOptions,
        ):
            raise ValidationError(
                "tracks must be LinearRequestTrackOptions or None."
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
        if self.comparison_table_file is not None:
            if (
                not isinstance(self.comparison_table_file, str)
                or not self.comparison_table_file.strip()
            ):
                raise ValidationError(
                    "comparison_table_file must be a non-empty string or None."
                )
            if self.blast_files is not None or self.linear_comparisons is not None:
                raise ValidationError(
                    "comparison_table_file cannot be combined with blast_files or "
                    "linear_comparisons."
                )
            object.__setattr__(
                self,
                "comparison_table_file",
                self.comparison_table_file.strip(),
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
        object.__setattr__(
            self,
            "pairwise_match_style",
            normalize_pairwise_match_style(self.pairwise_match_style),
        )
        if self.losat_search is not None and not isinstance(
            self.losat_search,
            LosatSearchOptions,
        ):
            raise ValidationError(
                "losat_search must be LosatSearchOptions or None.",
                diagnostic={"code": "COMPARISON_INPUT"},
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
        if not isinstance(self.collinear_infer_orthogroups, bool):
            raise ValidationError("collinear_infer_orthogroups must be a boolean")
        object.__setattr__(
            self,
            "collinear_max_paralog_links_per_orthogroup",
            _validate_positive_int(
                self.collinear_max_paralog_links_per_orthogroup,
                field_name="collinear_max_paralog_links_per_orthogroup",
            ),
        )


def _resolve_options_for_mode(
    options: CircularDiagramOptions | LinearDiagramOptions,
    *,
    mode: DiagramMode,
) -> CircularDiagramOptions | LinearDiagramOptions:
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


def resolve_circular_diagram_options(
    options: CircularDiagramOptions,
) -> CircularDiagramOptions:
    """Resolve omitted Circular request values through the Circular profile."""

    if not isinstance(options, CircularDiagramOptions):
        raise ValidationError("options must be CircularDiagramOptions.")
    return cast(
        CircularDiagramOptions,
        _resolve_options_for_mode(options, mode="circular"),
    )


def resolve_linear_diagram_options(
    options: LinearDiagramOptions,
) -> LinearDiagramOptions:
    """Resolve omitted Linear request values through the Linear profile."""

    if not isinstance(options, LinearDiagramOptions):
        raise ValidationError("options must be LinearDiagramOptions.")
    return cast(
        LinearDiagramOptions,
        _resolve_options_for_mode(options, mode="linear"),
    )


__all__ = [
    "CircularDiagramOptions",
    "CircularMultiRecordOptions",
    "CircularOutputOptions",
    "CircularRequestTrackOptions",
    "CircularTrackOptions",
    "LinearMultiRecordOptions",
    "LinearRecordTranslation",
    "LinearDiagramOptions",
    "LinearOutputOptions",
    "LinearRequestTrackOptions",
    "LinearTrackOptions",
    "LosatProgram",
    "LosatRuntimeOptions",
    "LosatSearchOptions",
    "LosatpMode",
    "AnnotationOptions",
    "ColorOptions",
    "DepthTrackInput",
    "resolve_circular_diagram_options",
    "resolve_linear_diagram_options",
]
