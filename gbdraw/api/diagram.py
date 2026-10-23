"""Typed diagram builders and their internal assembly helpers.

The public builders create diagrams from already-parsed in-memory records without
going through CLI argument parsing. Lower-level assembler functions remain private
engine boundaries used by those builders.

The functions here return an `svgwrite.Drawing` (SVG canvas). Saving/conversion is
handled separately by `gbdraw.api.save_figure_to` or `gbdraw.api.render_to_bytes`.
"""

from __future__ import annotations

import copy
import logging
import math
import re
import time
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any, Literal, Mapping, Optional, Sequence, cast

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from svgwrite import Drawing
from svgwrite.base import BaseElement

from gbdraw.analysis.collinearity import (
    CollinearityAnchorMode,
    CollinearityBlock,
    CollinearityColorMode,
    CollinearityResult,
    CollinearitySearchScope,
    LosslessCollinearityParameters,
    build_orthogroup_collinearity_blocks,
    convert_collinearity_blocks_to_comparisons,
    convert_collinearity_blocks_to_pair_comparisons,
    normalize_collinearity_anchor_mode,
    normalize_collinearity_color_mode,
    normalize_collinearity_search_scope,
)
from gbdraw.analysis.collinearity_units import CollinearityUnitMode
from gbdraw.analysis.conservation import (
    conservation_legend_sources,
    ConservationLoadResult,
    ConservationReferenceMode,
    ConservationTrack,
    load_conservation_sources,
    normalize_conservation_reference,
    read_conservation_inputs,
    normalize_conservation_tracks_for_record,
)
from gbdraw.analysis.depth import depth_df as build_depth_df
from gbdraw.analysis.depth_tracks import (
    DepthTrackData,
    DepthTrackSpec,
    build_depth_track_dataframes,
    depth_track_count,
    depth_track_data_count,
    depth_track_indices,
    depth_track_inputs,
    index_depth_track_row,
    normalize_depth_tracks,
    representative_depth_tracks,
)
from gbdraw.analysis.depth_tracks import (
    depth_track_heights as _depth_track_heights_from_specs,
)
from gbdraw.analysis.gc import circular_dinucleotide_content_df
from gbdraw.analysis.protein_colinearity import (
    LosatpCacheManager,
    LosatpCacheState,
    OrthogroupGraphResult,
    OrthogroupMembershipMode,
    OrthogroupResult,
    ProteinBlastpMode,
    ProteinBlastpResult,
    ProteinExtractionResult,
    build_pairwise_protein_blastp_comparisons,
    build_rbh_orthogroup_protein_blastp_comparisons,
    normalize_orthogroup_membership_mode,
)
from gbdraw.analysis.skew import skew_df
from gbdraw.annotations import (
    ResolvedAnnotationBundle,
    read_annotation_options,
    resolve_annotations,
)
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.record_planning import ComparisonRecordBinding
from gbdraw.api.options import (
    AnnotationOptions,
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    ColorOptions,
    DepthTrackInput,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    LosatRuntimeOptions,
    LosatSearchOptions,
    losatp_analysis_mode,
    resolve_circular_diagram_options,
    resolve_linear_diagram_options,
)
from gbdraw.auto_sizes import (
    depth_sliding_window,
    determine_length_parameter,
    sliding_window,
)
from gbdraw.canvas import CircularCanvasConfigurator, LinearCanvasConfigurator
from gbdraw.comparisons.losat_runtime import AUTOMATIC_LOSAT_BIN, LosatCacheWrite
from gbdraw.config.models import (
    CircularRenderProfile,
    GbdrawConfig,
    LinearRenderProfile,
)
from gbdraw.config.models.labels import LabelsFilteringConfig
from gbdraw.config.models.objects import (
    normalize_pairwise_match_style,
)
from gbdraw.config.modify import with_config_overrides, with_fit_row_overrides
from gbdraw.configurators import (
    BlastMatchConfigurator,
    DepthConfigurator,
    FeatureDrawingConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
    LegendDrawingConfigurator,
)
from gbdraw.configurators.legend import measure_legend
from gbdraw.core.color import normalize_hex_color
from gbdraw.core.fonts import resolve_font_family
from gbdraw.core.numeric import compensated_sum
from gbdraw.core.sequence import create_dict_for_sequence_lengths
from gbdraw.diagrams.circular.assemble import (
    CircularAssemblyResult,
    _assemble_circular_diagram_result,
    circular_legend_table,
    plan_circular_annotation_slots,
)
from gbdraw.diagrams.circular.positioning import _parse_svg_number
from gbdraw.diagrams.circular.radial_layout import track_layout_reason
from gbdraw.diagrams.linear import assemble_linear_diagram
from gbdraw.exceptions import LabelLayoutError, LabelsDoNotFitError, ValidationError
from gbdraw.features.colors import build_default_color_map, build_specific_color_rules
from gbdraw.features.placement import ResolvedRecordFeatureInputs
from gbdraw.features.visibility import (
    compile_feature_visibility_rules,
    read_feature_visibility_file,
)
from gbdraw.io.colors import (
    check_user_color,
    merge_default_colors,
    read_default_color_overrides,
    read_color_table,
    resolve_color_to_hex,
)
from gbdraw.labels.circular_fit import (
    CircularLabelFitSearch,
    CircularLabelLayout,
    check_circular_label_fit,
    circular_label_outcome,
)
from gbdraw.labels.filtering import (
    DERIVED_LABEL_FILTERING_KEYS,
    read_filter_list_file,
    read_label_override_file,
    read_qualifier_priority_file,
)
from gbdraw.labels.fit import (
    FITS,
    WARN_SEARCH_SECONDS,
    LabelFitReport,
    LabelFitRow,
    LabelFitRowKind,
    LabelFitValue,
    LabelFitVerification,
)
from gbdraw.labels.linear import (
    LINEAR_LABEL_FIT_ATTR,
    FeatureLabelRows,
    LinearLabelFit,
    check_linear_label_fit,
    keeps_feature_label_rows,
)
from gbdraw.layout.composition import (
    CompositionItem,
    CompositionRequest,
    LegendPlacement,
    TitlePlacement,
    plan_composition,
)
from gbdraw.layout.linear_multi_record import record_key_at, record_pairs_between_adjacent_rows
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.record_placement import resolve_record_row_positions
from gbdraw.layout.similarity_alignment import SimilarityAlignmentPlan
from gbdraw.layout.spatial import Aabb
from gbdraw.legend.row_facts import LegendRowFacts, attach_legend_row_facts
from gbdraw.legend.sources import LegendSources
from gbdraw.legend.table import (
    LegendTable,
    pairwise_identity_legend_from_comparisons,
)
from gbdraw.linear_comparison import (
    ComparisonRecordIdWarning,
    LinearComparison,
    ProteinComparisonFrames,
)
from gbdraw.mode_profiles import (
    CIRCULAR_MODE_PROFILE,
    DEFAULT_FEATURE_TYPES,
    LINEAR_MODE_PROFILE,
    ComparisonThresholds,
)
from gbdraw.render.composition import (
    COMPOSITION_ROLE_ATTRIBUTE,
    apply_composition_plan,
)
from gbdraw.render.group_content import (
    CANVAS_SCOPE,
    ContentGroup,
    GroupContent,
    GroupKey,
    GroupKind,
    build_group_content,
    place_content,
)
from gbdraw.render.groups.circular import CircularDefinitionStyle, DefinitionGroup
from gbdraw.render.groups.circular.definition import CircularDefinitionProfile, definition_record_facts
from gbdraw.render.groups.circular.legend import circular_legend_content
from gbdraw.render_stages import render_stage
from gbdraw.svg.elements import Group
from gbdraw.svg.ids import definition_group_svg_id, instance_svg_id
from gbdraw.tracks import (
    CircularTrackSlot,
    LinearTrackSlot,
    ScalarSpec,
    circular_track_slots_from_order,
    default_circular_track_slots,
    normalize_circular_track_slots_with_axis,
    normalize_linear_track_slots_with_axis,
    parse_circular_track_slots,
    parse_linear_track_slots,
    parse_nonnegative_integer,
)
from gbdraw.tracks.parsing import slot_dinucleotide

from .prepared import ResolvedFeatureInputs

if TYPE_CHECKING:
    from gbdraw.analysis.conservation import ConservationInput
    from gbdraw.api.requests import RecordPresentation


DEFAULT_SELECTED_FEATURES = DEFAULT_FEATURE_TYPES

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class LinearDiagramMetadata:
    """Computed Linear analysis values consumed after diagram assembly."""

    protein_comparisons: tuple[DataFrame, ...] | None = None
    linear_comparisons: tuple[LinearComparison, ...] = ()
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None = None
    collinearity_result: CollinearityResult | None = None
    comparison_warnings: tuple[ComparisonRecordIdWarning, ...] = ()
    # Each record's search-frame binding, set by the plan's build when the
    # values above hold source-bound results: a Session stores them in F (R5).
    search_frame_bindings: tuple[ComparisonRecordBinding, ...] = ()


@dataclass(frozen=True)
class LinearDiagramBuildResult:
    """A Linear drawing paired with the analysis values that produced it."""

    drawing: Drawing
    metadata: LinearDiagramMetadata


@dataclass(frozen=True, slots=True)
class _CircularGridLayout:
    """Visible-bounds packing result for Circular record wrappers."""

    bounds: Aabb
    translations: tuple[tuple[float, float], ...]


_MULTI_RECORD_SUFFIXED_TOP_LEVEL_IDS = {
    "Axis",
    "tick",
    "labels",
    "gc_content",
    "skew",
    "gc_skew",
}
_MULTI_RECORD_SIZE_MODES = {"linear", "auto", "equal"}
_DEFINITION_POSITIONS = {"center", "top", "bottom"}
_LINEAR_PLOT_TITLE_POSITIONS = {"center", "top", "bottom"}
_CIRCULAR_PLOT_TITLE_POSITIONS = {"none", "top", "bottom"}
# Legacy fallback used only when a valid radius cannot be derived.
_MULTI_RECORD_GRID_GAP_PX = 16.0
_MULTI_RECORD_COLUMN_GAP_RATIO = 0.10
_MULTI_RECORD_ROW_GAP_RATIO = 0.05
_PLOT_TITLE_DEFAULT_FONT_SIZE = 32.0

def _resolve_single_circular_depth_options(
    options: CircularDiagramOptions,
) -> tuple[DataFrame | None, str | None]:
    singular_present = options.depth_table is not None or options.depth_file is not None
    plural_present = options.depth_tables is not None or options.depth_files is not None
    if singular_present and plural_present:
        raise ValidationError(
            "Single circular depth inputs cannot combine singular and plural forms."
        )
    if options.depth_table is not None and options.depth_file is not None:
        raise ValidationError("Pass either the depth table form or depth file form, not both.")
    if options.depth_tables is not None and options.depth_files is not None:
        raise ValidationError("Pass either the depth table form or depth file form, not both.")
    if options.depth_tables is not None:
        if len(options.depth_tables) != 1:
            raise ValidationError("Expected one depth table for one circular record.")
        return options.depth_tables[0], None
    if options.depth_files is not None:
        if len(options.depth_files) != 1:
            raise ValidationError("Expected one depth file for one circular record.")
        return None, options.depth_files[0]
    return options.depth_table, options.depth_file


def _resolve_optional_table(
    *,
    name: str,
    table: DataFrame | None,
    file_path: str | None,
    reader,
) -> tuple[bool, DataFrame | None]:
    if table is not None and file_path is not None:
        raise ValidationError(f"Pass either {name}_table or {name}_file, not both.")
    if table is not None:
        return True, table
    if file_path is not None:
        return True, reader(file_path)
    return False, None


# Requested font-family lists already warned about in this process.
_warned_font_families: set[str] = set()


def _with_available_font_family(cfg: GbdrawConfig) -> GbdrawConfig:
    """Draw and measure text only with font families that resolve.

    A family that is not available is replaced with the nearest one
    (``resolve_font_family``), with one warning per requested list and
    process naming both when the requested first family is not drawn.
    """

    resolution = resolve_font_family(cfg.objects.text.font_family)
    if not resolution.substituted:
        return cfg
    # The warning only when the drawn family is not the requested first one (OV-743).
    if resolution.message is not None and resolution.requested not in _warned_font_families:
        _warned_font_families.add(resolution.requested)
        logger.warning(resolution.message)
    text = replace(cfg.objects.text, font_family=resolution.effective)
    return replace(cfg, objects=replace(cfg.objects, text=text))


def _resolve_diagram_options_config(
    options: CircularDiagramOptions | LinearDiagramOptions,
) -> GbdrawConfig:
    """Resolve one typed config with available fonts, and losslessly attach
    explicit label tables."""

    cfg = _with_available_font_family(
        apply_config_overrides(options.config, options.config_overrides)
    )

    whitelist_given, whitelist = _resolve_optional_table(
        name="label_whitelist",
        table=options.label_whitelist_table,
        file_path=options.label_whitelist_file,
        reader=read_filter_list_file,
    )
    priority_given, priority = _resolve_optional_table(
        name="qualifier_priority",
        table=options.qualifier_priority_table,
        file_path=options.qualifier_priority_file,
        reader=read_qualifier_priority_file,
    )
    override_given, label_override = _resolve_optional_table(
        name="label_override",
        table=options.label_override_table,
        file_path=options.label_override_file,
        reader=read_label_override_file,
    )
    if not (whitelist_given or priority_given or override_given):
        return cfg

    filtering = copy.deepcopy(cfg.labels.filtering.as_dict())
    # Maps compiled by a previous render (for example, persisted in a CLI
    # Session config) would otherwise shadow the attached tables.
    for derived_key in DERIVED_LABEL_FILTERING_KEYS:
        filtering.pop(derived_key, None)
    if whitelist_given:
        filtering["whitelist_df"] = whitelist
    if priority_given:
        filtering["qualifier_priority_df"] = priority
    if override_given:
        filtering["label_override_df"] = label_override
    cfg = replace(
        cfg,
        labels=replace(
            cfg.labels,
            filtering=LabelsFilteringConfig.from_dict(filtering),
        ),
    )
    return cfg


def _parse_circular_track_slot_inputs(
    circular_track_slots: Sequence[str | CircularTrackSlot] | None,
) -> list[CircularTrackSlot] | None:
    if circular_track_slots is None:
        return None
    return parse_circular_track_slots(list(circular_track_slots))


def _validate_circular_track_axis_index(
    circular_track_axis_index: int | None,
    parsed_circular_track_slots: Sequence[CircularTrackSlot] | None,
) -> int | None:
    if circular_track_axis_index is None:
        return None
    if parsed_circular_track_slots is None:
        raise ValidationError("circular_track_axis_index requires circular_track_slots.")
    if not isinstance(circular_track_axis_index, int):
        raise ValidationError("circular_track_axis_index must be an integer.")
    if circular_track_axis_index < 0 or circular_track_axis_index > len(parsed_circular_track_slots):
        raise ValidationError(
            f"circular_track_axis_index must be between 0 and the number of circular track slots ({len(parsed_circular_track_slots)})."
        )
    try:
        normalize_circular_track_slots_with_axis(
            parsed_circular_track_slots,
            circular_track_axis_index,
        )
    except Exception as exc:
        raise ValidationError(str(exc)) from exc
    return circular_track_axis_index


def _slots_have_renderer(
    slots: Sequence[CircularTrackSlot | LinearTrackSlot] | None,
    renderer: str,
) -> bool:
    return any(slot.enabled and str(slot.renderer) == renderer for slot in (slots or []))


def _resolve_circular_track_visibility(
    cfg: GbdrawConfig,
    slots: Sequence[CircularTrackSlot] | None,
) -> tuple[bool, bool, bool]:
    """Resolve canvas track flags, including the horizontal inner-label policy."""

    if slots is not None:
        return (
            _slots_have_renderer(slots, "depth"),
            _slots_have_renderer(slots, "dinucleotide_content"),
            _slots_have_renderer(slots, "dinucleotide_skew"),
        )

    show_gc = bool(cfg.canvas.show_gc)
    show_skew = bool(cfg.canvas.show_skew)
    circular_labels = cfg.labels.circular
    if circular_labels.scope == "both" and circular_labels.placement == "horizontal":
        show_gc = False
        show_skew = False
    return bool(cfg.canvas.show_depth), show_gc, show_skew


def _circular_slots_define_renderer(
    slots: Sequence[CircularTrackSlot] | None,
    renderer: str,
) -> bool:
    return any(str(slot.renderer) == renderer for slot in (slots or []))


def _dinucleotides_from_circular_slots(
    slots: Sequence[CircularTrackSlot] | None,
    *,
    default_nt: str,
) -> set[str]:
    nts: set[str] = set()
    for slot in slots or []:
        if not slot.enabled or str(slot.renderer) not in {"dinucleotide_content", "dinucleotide_skew"}:
            continue
        nts.add(slot_dinucleotide(slot.params, default_nt))
    return nts


def _parse_linear_track_slot_inputs(
    linear_track_slots: Sequence[str | LinearTrackSlot] | None,
) -> list[LinearTrackSlot] | None:
    if linear_track_slots is None:
        return None
    return parse_linear_track_slots(list(linear_track_slots))


def _validate_linear_track_axis_index(
    linear_track_axis_index: int | None,
    parsed_linear_track_slots: Sequence[LinearTrackSlot] | None,
) -> int | None:
    if linear_track_axis_index is None:
        return None
    if parsed_linear_track_slots is None:
        raise ValidationError("linear_track_axis_index requires linear_track_slots.")
    if not isinstance(linear_track_axis_index, int):
        raise ValidationError("linear_track_axis_index must be an integer.")
    if linear_track_axis_index < 0 or linear_track_axis_index > len(parsed_linear_track_slots):
        raise ValidationError(
            f"linear_track_axis_index must be between 0 and the number of linear track slots ({len(parsed_linear_track_slots)})."
        )
    try:
        normalize_linear_track_slots_with_axis(
            parsed_linear_track_slots,
            linear_track_axis_index,
        )
    except Exception as exc:
        raise ValidationError(str(exc)) from exc
    return linear_track_axis_index



def _build_circular_dinucleotide_content_dataframes(
    record: SeqRecord,
    *,
    window: int,
    step: int,
    nts: set[str],
) -> dict[str, DataFrame]:
    return {
        nt: render_stage("gc_skew_data", circular_dinucleotide_content_df, record, window, step, nt)
        for nt in sorted(nts)
    }


def _build_circular_dinucleotide_skew_dataframes(
    record: SeqRecord,
    *,
    window: int,
    step: int,
    nts: set[str],
) -> dict[str, DataFrame]:
    return {nt: render_stage("gc_skew_data", skew_df, record, window, step, nt) for nt in sorted(nts)}


def _has_conservation_inputs(
    conservation_blast_files: Sequence[ConservationInput] | None,
    conservation_dataframes: Sequence[DataFrame] | None,
) -> bool:
    return bool(conservation_blast_files or conservation_dataframes)


def _load_conservation_result(
    *,
    conservation_blast_files: Sequence[ConservationInput] | None,
    conservation_dataframes: Sequence[DataFrame] | None,
    conservation_labels: Sequence[str] | None,
    conservation_colors: Sequence[str] | None,
    thresholds: ComparisonThresholds,
) -> ConservationLoadResult | None:
    """The filtered conservation sources, from the ``conservation_sources`` stage.

    Files are read outside the stage, so a staged file is keyed by its
    resource identity, never by its path.
    """

    if not _has_conservation_inputs(conservation_blast_files, conservation_dataframes):
        return None
    return render_stage(
        "conservation_sources",
        load_conservation_sources,
        blast_config=thresholds,
        conservation_files=read_conservation_inputs(conservation_blast_files),
        conservation_dataframes=conservation_dataframes,
        labels=conservation_labels,
        colors=conservation_colors,
    )


def _conservation_tracks_for_record(
    load_result: ConservationLoadResult | None,
    records: Sequence[SeqRecord],
    record: SeqRecord,
    conservation_reference: ConservationReferenceMode,
) -> tuple[ConservationTrack, ...]:
    """The rings of one displayed record, from the ``conservation_tracks`` stage."""

    if load_result is None:
        return ()
    return render_stage(
        "conservation_tracks",
        normalize_conservation_tracks_for_record,
        load_result,
        displayed_records=tuple(records),
        record=record,
        conservation_reference=conservation_reference,
    )


def _conservation_ring_width(
    cfg: GbdrawConfig,
    explicit_width: float | None,
) -> float | None:
    return float(explicit_width) if explicit_width is not None else cfg.objects.conservation.ring_width


def _conservation_ring_gap(
    cfg: GbdrawConfig,
    explicit_gap: float | None,
) -> float | None:
    return float(explicit_gap) if explicit_gap is not None else cfg.objects.conservation.ring_gap


def _conservation_slots_for_tracks(
    tracks: Sequence[ConservationTrack],
    *,
    cfg: GbdrawConfig,
    conservation_ring_width: float | None,
    conservation_ring_gap: float | None,
) -> list[CircularTrackSlot]:
    width = _conservation_ring_width(cfg, conservation_ring_width)
    gap = _conservation_ring_gap(cfg, conservation_ring_gap)
    return [
        CircularTrackSlot(
            id=f"conservation_{int(track.track_index)}",
            renderer="sequence_conservation",
            side="inside",
            width=ScalarSpec(float(width), "px") if width is not None else None,
            inner_gap_px=float(gap) if gap is not None else None,
            outer_gap_px=float(gap) if gap is not None else None,
            params={
                "track_index": int(track.track_index),
                "source_index": int(track.source_index),
                "label": track.track_label,
                "color": track.track_color or "",
            },
        )
        for track in tracks
    ]


def _insert_conservation_slots(
    base_slots: Sequence[CircularTrackSlot],
    conservation_slots: Sequence[CircularTrackSlot],
    *,
    axis_index: int | None,
) -> list[CircularTrackSlot]:
    if not conservation_slots:
        return list(base_slots)
    out: list[CircularTrackSlot]
    if axis_index is not None:
        out = list(base_slots)
        out[axis_index:axis_index] = conservation_slots
        return out
    out = []
    inserted = False
    for slot in base_slots:
        out.append(slot)
        if not inserted and str(slot.renderer) == "features":
            out.extend(conservation_slots)
            inserted = True
    if not inserted:
        out = list(conservation_slots) + out
    return out


def _apply_conservation_track_params_to_slots(
    base_slots: Sequence[CircularTrackSlot],
    tracks: Sequence[ConservationTrack],
    *,
    cfg: GbdrawConfig,
    conservation_ring_width: float | None,
    conservation_ring_gap: float | None,
) -> list[CircularTrackSlot]:
    if not tracks:
        return list(base_slots)
    width = _conservation_ring_width(cfg, conservation_ring_width)
    gap = _conservation_ring_gap(cfg, conservation_ring_gap)
    out: list[CircularTrackSlot] = []
    track_order = list(tracks)
    next_track_index = 0
    used_source_indexes: set[int] = set()
    for slot in base_slots:
        if str(slot.renderer) != "sequence_conservation":
            out.append(slot)
            continue

        params = dict(slot.params or {})
        raw_track_index = params.get("track_index")
        track: ConservationTrack | None = None
        raw_source_index = params.get("source_index")
        if raw_source_index is not None:
            try:
                wanted_source_index = int(raw_source_index)
                track = next(
                    item
                    for item in track_order
                    if int(item.source_index) == wanted_source_index
                )
            except Exception:
                track = None
        if track is None and raw_track_index is not None:
            try:
                wanted_track_index = int(raw_track_index)
                track = next(
                    item
                    for item in track_order
                    if int(item.track_index) == wanted_track_index
                )
            except Exception:
                track = None
        if track is None and next_track_index < len(track_order):
            while (
                next_track_index < len(track_order)
                and int(track_order[next_track_index].source_index) in used_source_indexes
            ):
                next_track_index += 1
            if next_track_index < len(track_order):
                track = track_order[next_track_index]
                next_track_index += 1

        if track is not None:
            used_source_indexes.add(int(track.source_index))
            params.setdefault("track_index", int(track.track_index))
            params.setdefault("source_index", int(track.source_index))
            params.setdefault("label", track.track_label)
            params.setdefault("color", track.track_color or "")
        out.append(
            replace(
                slot,
                side=slot.side,
                width=slot.width or (ScalarSpec(float(width), "px") if width is not None else None),
                inner_gap_px=(
                    slot.inner_gap_px
                    if slot.inner_gap_px is not None
                    else (float(gap) if gap is not None else None)
                ),
                outer_gap_px=(
                    slot.outer_gap_px
                    if slot.outer_gap_px is not None
                    else (float(gap) if gap is not None else None)
                ),
                params=params,
            )
        )
    return out


def _resolve_circular_window_step(
    record: SeqRecord,
    cfg: GbdrawConfig,
    *,
    window: int | None,
    step: int | None,
) -> tuple[int, int]:
    """Resolve circular window/step with record-length defaults."""
    auto_window, auto_step = sliding_window(len(record.seq), cfg)
    return (
        int(window) if window is not None else auto_window,
        int(step) if step is not None else auto_step,
    )


def _validate_depth_bounds(min_depth: float | None, max_depth: float | None) -> None:
    if min_depth is not None and float(min_depth) < 0:
        raise ValidationError("depth_min must be >= 0")
    if max_depth is not None and float(max_depth) < 0:
        raise ValidationError("depth_max must be >= 0")
    if min_depth is not None and max_depth is not None and float(min_depth) > float(max_depth):
        raise ValidationError("depth_min must be <= depth_max")


def _validate_depth_config(depth_config) -> None:
    _validate_depth_bounds(depth_config.min_depth, depth_config.max_depth)
    if (
        depth_config.large_tick_interval is not None
        and float(depth_config.large_tick_interval) <= 0
    ):
        raise ValidationError("depth_large_tick_interval must be > 0")
    if depth_config.small_tick_interval is not None and float(depth_config.small_tick_interval) <= 0:
        raise ValidationError("depth_small_tick_interval must be > 0")


def _validate_gc_content_config(gc_content_config) -> None:
    if str(getattr(gc_content_config, "mode", "deviation")) not in {"deviation", "percent"}:
        raise ValidationError("gc_content_mode must be one of: deviation, percent")
    min_percent = getattr(gc_content_config, "min_percent", None)
    max_percent = getattr(gc_content_config, "max_percent", None)
    if min_percent is not None and not math.isfinite(float(min_percent)):
        raise ValidationError("gc_content_min_percent must be a finite number")
    if max_percent is not None and not math.isfinite(float(max_percent)):
        raise ValidationError("gc_content_max_percent must be a finite number")
    if min_percent is not None and max_percent is not None and float(min_percent) > float(max_percent):
        raise ValidationError("gc_content_min_percent must be <= gc_content_max_percent")
    for attr, label in (
        ("large_tick_interval", "gc_content_large_tick_interval"),
        ("small_tick_interval", "gc_content_small_tick_interval"),
    ):
        value = getattr(gc_content_config, attr, None)
        if value is not None and (not math.isfinite(float(value)) or float(value) <= 0):
            raise ValidationError(f"{label} must be > 0")
    percent_background_opacity = getattr(gc_content_config, "percent_background_opacity", 1.0)
    if (
        not math.isfinite(float(percent_background_opacity))
        or float(percent_background_opacity) < 0.0
        or float(percent_background_opacity) > 1.0
    ):
        raise ValidationError("gc_content_percent_background_opacity must be between 0 and 1")
    percent_border_width = getattr(gc_content_config, "percent_border_width", 0.8)
    if not math.isfinite(float(percent_border_width)) or float(percent_border_width) < 0.0:
        raise ValidationError("gc_content_percent_border_width must be >= 0")


def _validate_positive_optional(name: str, value: int | None) -> None:
    if value is not None and int(value) <= 0:
        raise ValidationError(f"{name} must be > 0")


def _validate_positive_float_optional(name: str, value: float | None) -> None:
    if value is not None and (not math.isfinite(float(value)) or float(value) <= 0):
        raise ValidationError(f"{name} must be > 0")


def _validate_nonnegative_float_optional(name: str, value: float | None) -> None:
    if value is not None and (not math.isfinite(float(value)) or float(value) < 0):
        raise ValidationError(f"{name} must be a finite number >= 0")


def _resolve_depth_window_step(
    *,
    window: int,
    step: int,
    depth_window: int | None,
    depth_step: int | None,
) -> tuple[int, int]:
    """Resolve depth window/step from GC/skew settings.

    Depth uses a denser sampling step than GC/skew, but keeps at least a 100 bp
    aggregation window so the default track is not overly noisy.
    """

    _validate_positive_optional("depth_window", depth_window)
    _validate_positive_optional("depth_step", depth_step)
    auto_depth_window, auto_depth_step = depth_sliding_window(window, step)
    return (
        int(depth_window) if depth_window is not None else auto_depth_window,
        int(depth_step) if depth_step is not None else auto_depth_step,
    )


def _cfg_with_depth_scale_max(cfg: GbdrawConfig, max_depth: float | None) -> GbdrawConfig:
    if max_depth is None:
        return cfg
    return replace(
        cfg,
        objects=replace(
            cfg.objects,
            depth=replace(cfg.objects.depth, max_depth=float(max_depth)),
        ),
    )










def _depth_track_count_for_render(
    record_depth_tracks: Sequence[Sequence[DepthTrackSpec]] | None,
    precomputed_depth_tracks: Sequence[DepthTrackData] | None,
    precomputed_depth_df: DataFrame | None = None,
) -> int:
    if precomputed_depth_tracks:
        return depth_track_data_count([precomputed_depth_tracks])
    count = depth_track_count(record_depth_tracks)
    if count > 0:
        return count
    return 1 if precomputed_depth_df is not None else 0


def _default_circular_depth_slots_if_needed(
    *,
    parsed_circular_track_slots: list[CircularTrackSlot] | None,
    show_depth: bool,
    depth_track_count_value: int,
    show_ticks: bool,
    show_gc: bool,
    show_skew: bool,
    dinucleotide: str,
) -> list[CircularTrackSlot] | None:
    if parsed_circular_track_slots is not None:
        return parsed_circular_track_slots
    if not show_depth or depth_track_count_value <= 1:
        return None
    return default_circular_track_slots(
        show_ticks=show_ticks,
        show_depth=True,
        depth_track_count=depth_track_count_value,
        show_gc=show_gc,
        show_skew=show_skew,
        dinucleotide=dinucleotide,
    )


def _validate_depth_track_indices(
    slots: Sequence[CircularTrackSlot | LinearTrackSlot] | None,
    *,
    depth_track_count_value: int,
    available_indices: Sequence[int] | None = None,
) -> None:
    if not slots:
        return
    available_index_set = set(available_indices) if available_indices is not None else None
    for slot in slots:
        if not slot.enabled or str(slot.renderer) != "depth":
            continue
        raw_index = (slot.params or {}).get("track_index", 0)
        try:
            track_index = parse_nonnegative_integer(
                raw_index,
                field_name=f"Depth slot '{slot.id}' track_index",
            )
        except (TypeError, ValueError) as exc:
            raise ValidationError(
                f"Depth slot '{slot.id}' has invalid track_index={raw_index!r}."
            ) from exc
        if (
            track_index < 0
            or track_index >= depth_track_count_value
            or (available_index_set is not None and track_index not in available_index_set)
        ):
            available_text = (
                f"; globally available logical indices: {', '.join(map(str, sorted(available_index_set))) or 'none'}"
                if available_index_set is not None
                else ""
            )
            raise ValidationError(
                f"Depth slot '{slot.id}' track_index={track_index} is outside the available "
                f"depth track range 0..{max(0, depth_track_count_value - 1)}{available_text}."
            )


def _estimate_square_grid(record_count: int) -> tuple[int, int]:
    """Return grid dimensions (cols, rows) close to square."""
    if record_count <= 0:
        return 1, 1
    cols = int(math.ceil(math.sqrt(record_count)))
    rows = int(math.ceil(float(record_count) / float(cols)))
    return cols, rows


def _resolve_multi_record_default_row_counts(record_count: int) -> list[int]:
    """Resolve default near-square row counts for multi-record canvas layout."""
    if record_count <= 0:
        return []
    cols, rows = _estimate_square_grid(record_count)
    counts: list[int] = []
    for row in range(rows):
        start = row * cols
        end = min(record_count, start + cols)
        if end > start:
            counts.append(end - start)
    return counts


def _pack_circular_record_bounds(
    record_bounds: Sequence[Aabb],
    row_counts: Sequence[int],
    *,
    column_gap_px: float,
    row_gap_px: float,
) -> _CircularGridLayout:
    """Pack record content by authoritative visible bounds."""
    if not record_bounds:
        raise ValueError("record_bounds is empty")
    if (
        not math.isfinite(float(column_gap_px))
        or not math.isfinite(float(row_gap_px))
        or column_gap_px < 0.0
        or row_gap_px < 0.0
    ):
        raise ValueError("record grid gaps must be finite and non-negative")

    rows: list[tuple[int, ...]] = []
    cursor = 0
    for raw_count in row_counts:
        count = max(0, int(raw_count))
        row = tuple(range(cursor, min(cursor + count, len(record_bounds))))
        if row:
            rows.append(row)
            cursor += len(row)
        if cursor >= len(record_bounds):
            break
    if cursor < len(record_bounds):
        rows.append(tuple(range(cursor, len(record_bounds))))

    row_widths = tuple(
        compensated_sum(record_bounds[index].width for index in row)
        + max(0, len(row) - 1) * float(column_gap_px)
        for row in rows
    )
    row_heights = tuple(
        max(record_bounds[index].height for index in row)
        for row in rows
    )
    grid_width = max(row_widths)
    grid_height = compensated_sum(row_heights) + max(0, len(rows) - 1) * float(row_gap_px)

    translations: list[tuple[float, float] | None] = [None] * len(record_bounds)
    row_top = 0.0
    for row, row_width, row_height in zip(rows, row_widths, row_heights):
        content_left = 0.5 * (grid_width - row_width)
        for position, record_index in enumerate(row):
            bounds = record_bounds[record_index]
            translations[record_index] = (
                content_left - bounds.min_x,
                row_top + (0.5 * (row_height - bounds.height)) - bounds.min_y,
            )
            content_left += bounds.width
            if position < len(row) - 1:
                content_left += float(column_gap_px)
        row_top += row_height + float(row_gap_px)

    if any(translation is None for translation in translations):
        raise RuntimeError("circular record grid left a record unplaced")
    resolved_translations = tuple(
        translation for translation in translations if translation is not None
    )
    return _CircularGridLayout(
        bounds=Aabb(0.0, 0.0, grid_width, grid_height),
        translations=resolved_translations,
    )


def _strip_nested_composition_roles(element: object) -> None:
    """Remove single-record role markers before nesting them in a grid."""
    for child in _walk_svgwrite_elements(element):
        attribs = getattr(child, "attribs", None)
        if isinstance(attribs, dict):
            attribs.pop(COMPOSITION_ROLE_ATTRIBUTE, None)


def _require_circular_assembly_result(drawing: Drawing) -> CircularAssemblyResult:
    """Return the authoritative record result or fail at the assembly boundary."""
    result = getattr(drawing, "_gbdraw_circular_assembly_result", None)
    if not isinstance(result, CircularAssemblyResult):
        raise RuntimeError(
            "Circular record assembly did not expose authoritative layout bounds."
        )
    return result


def _suffix_fixed_top_level_group_id(element: object, record_index: int) -> None:
    """Suffix selected top-level group IDs to avoid collisions on merged canvas."""
    attribs = getattr(element, "attribs", None)
    if not isinstance(attribs, dict):
        return
    group_id = attribs.get("id")
    if group_id in _MULTI_RECORD_SUFFIXED_TOP_LEVEL_IDS:
        attribs["id"] = f"{group_id}_{record_index}"
    elif isinstance(group_id, str) and group_id.startswith("depth_"):
        if group_id.endswith("_axis"):
            attribs["id"] = f"{group_id[:-5]}_record_{record_index + 1}_axis"
        else:
            attribs["id"] = f"{group_id}_record_{record_index + 1}"
    for child in getattr(element, "elements", []) or []:
        child_attribs = getattr(child, "attribs", None)
        if not isinstance(child_attribs, dict):
            continue
        child_id = child_attribs.get("id")
        if isinstance(child_id, str) and child_id.startswith("depth_"):
            if child_id.endswith("_axis"):
                child_attribs["id"] = f"{child_id[:-5]}_record_{record_index + 1}_axis"
            else:
                child_attribs["id"] = f"{child_id}_record_{record_index + 1}"


_SVG_URL_REFERENCE_RE = re.compile(r"url\(#([^)]+)\)")


def _walk_svgwrite_elements(element: object):
    yield element
    for child in getattr(element, "elements", []) or []:
        yield from _walk_svgwrite_elements(child)


def _uniquify_copied_subtrees_ids(
    elements: Sequence[object],
    *,
    record_index: int,
    used_ids: set[str],
    bind_record_identity: bool = False,
) -> None:
    """Repair copied multi-record subtrees with one shared reference namespace."""

    renamed: dict[str, str] = {}
    for element in elements:
        for child in _walk_svgwrite_elements(element):
            attribs = getattr(child, "attribs", None)
            if not isinstance(attribs, dict):
                continue
            if bind_record_identity and attribs.get("data-gbdraw-feature-id"):
                stable_id = attribs["data-gbdraw-feature-id"]
                rendered_id = attribs.get("data-gbdraw-rendered-feature-id") or stable_id
                attribs["data-gbdraw-rendered-feature-id"] = instance_svg_id(rendered_id, f"record_{record_index + 1}")
                attribs["data-gbdraw-stable-feature-id"] = stable_id
                attribs["data-gbdraw-record-index"] = str(record_index)
            if bind_record_identity and attribs.get("data-label-feature-id"):
                attribs["data-label-feature-id"] = instance_svg_id(
                    attribs["data-label-feature-id"], f"record_{record_index + 1}"
                )
            raw_id = attribs.get("id")
            if not isinstance(raw_id, str) or not raw_id:
                continue
            candidate = raw_id
            if candidate in used_ids:
                base = f"{raw_id}__record_{record_index + 1}"
                candidate = base
                duplicate_index = 2
                while candidate in used_ids:
                    candidate = f"{base}_{duplicate_index}"
                    duplicate_index += 1
                attribs["id"] = candidate
                renamed.setdefault(raw_id, candidate)
            used_ids.add(candidate)

    if not renamed:
        return

    def replace_url_reference(match: re.Match[str]) -> str:
        referenced_id = match.group(1)
        return f"url(#{renamed.get(referenced_id, referenced_id)})"

    for element in elements:
        for child in _walk_svgwrite_elements(element):
            attribs = getattr(child, "attribs", None)
            if not isinstance(attribs, dict):
                continue
            for name, value in tuple(attribs.items()):
                if not isinstance(value, str):
                    continue
                rewritten = _SVG_URL_REFERENCE_RE.sub(replace_url_reference, value)
                if (
                    name in {"href", "xlink:href"}
                    or str(name).rsplit("}", 1)[-1] == "href"
                ) and rewritten.startswith("#"):
                    rewritten = f"#{renamed.get(rewritten[1:], rewritten[1:])}"
                if rewritten != value:
                    attribs[name] = rewritten


# One placed top-level group of a record drawing: its content and its attributes in document order.
PlacedRecordContent = tuple[GroupContent, tuple[tuple[str, object], ...]]


@dataclass(frozen=True)
class MergedCircularRecord:
    """Result of the ``record_merge`` stage: one record's ``record_<i>`` group content.

    ``used_ids`` holds the ids of this record and of every earlier one; the next
    record's merge takes this result, so its ids depend on all earlier records.
    """

    content: GroupContent
    used_ids: frozenset[str]


def placed_record_contents(
    drawing: Drawing,
) -> tuple[tuple[PlacedRecordContent, ...], tuple[tuple[int, int], ...]]:
    """The placed contents of a single-record Circular drawing and where each of its definitions comes from.

    A definition's position is (index of the placed content, index in its definitions).
    """

    placed: list[PlacedRecordContent] = []
    definitions: list[object] = []
    for element in drawing.elements:
        if _is_defs_element(element):
            definitions.extend(getattr(element, "elements", ()) or ())
            continue
        if not isinstance(element, ContentGroup):
            raise RuntimeError(f"a Circular record drawing holds a group that is not a content: {element!r}")
        placed.append((element.content, tuple(element.attribs.items())))
    owners = {
        id(definition): (content_index, definition_index)
        for content_index, (content, _attributes) in enumerate(placed)
        for definition_index, definition in enumerate(content.definitions)
    }
    try:
        positions = tuple(owners[id(definition)] for definition in definitions)
    except KeyError as error:
        raise RuntimeError("a Circular record drawing holds a definition that no content holds") from error
    return tuple(placed), positions


def merge_circular_record(
    placed: tuple[PlacedRecordContent, ...],
    definition_positions: tuple[tuple[int, int], ...],
    *,
    record_index: int,
    record_id: str,
    bind_record_identity: bool,
    previous: MergedCircularRecord | None,
) -> MergedCircularRecord:
    """The ``record_merge`` stage: copy one record's groups into its ``record_<i>`` group.

    The copies drop the single-record composition roles, suffix the fixed
    top-level ids, and rename ids that an earlier record (``previous``) or this
    one already uses, with the references to them. The stored contents are not
    changed: their children are deep copies (perf-0140 D-07 keeps the copy).
    """

    used_ids = set(previous.used_ids) if previous is not None else set()
    used_ids.add(f"record_{record_index}")
    copied_definitions: list[BaseElement] = [
        copy.deepcopy(placed[content_index][0].definitions[definition_index])
        for content_index, definition_index in definition_positions
    ]
    copied_elements: list[Group] = []
    for content, attributes in placed:
        group = Group(debug=False)
        if content.parameter is not None:
            group._parameter = content.parameter
        group.attribs = dict(attributes)
        group.elements = list(content.children)
        copied = copy.deepcopy(group)
        _strip_nested_composition_roles(copied)
        _suffix_fixed_top_level_group_id(copied, record_index)
        copied_elements.append(copied)
    for definition in copied_definitions:
        _strip_nested_composition_roles(definition)
    _uniquify_copied_subtrees_ids(
        (*copied_definitions, *copied_elements),
        record_index=record_index,
        used_ids=used_ids,
        bind_record_identity=bind_record_identity,
    )
    record_group = Group(id=f"record_{record_index}")
    record_group.attribs["data-gbdraw-record-id"] = record_id
    record_group.attribs["data-gbdraw-record-index"] = str(record_index)
    for copied in copied_elements:
        record_group.add(copied)
    return MergedCircularRecord(
        content=GroupContent.of(record_group, tuple(copied_definitions)),
        used_ids=frozenset(used_ids),
    )


def _is_defs_element(element: object) -> bool:
    """Return True if the element is an SVG <defs> container."""
    class_name = element.__class__.__name__.lower()
    if class_name == "defs":
        return True
    element_name = getattr(element, "elementname", None)
    if isinstance(element_name, str) and element_name.lower() == "defs":
        return True
    return False


def _resolve_multi_record_size_mode(mode: str) -> Literal["linear", "auto", "equal"]:
    """Normalize and validate multi-record size mode."""
    normalized = str(mode).strip().lower()
    if normalized not in _MULTI_RECORD_SIZE_MODES:
        raise ValidationError(
            "multi_record_size_mode must be one of: auto, linear, equal"
        )
    return cast(Literal["linear", "auto", "equal"], normalized)




def _resolve_linear_plot_title_position(
    position: str,
) -> Literal["center", "top", "bottom"]:
    normalized = str(position).strip().lower()
    if normalized not in _LINEAR_PLOT_TITLE_POSITIONS:
        raise ValidationError(
            "plot_title_position must be one of: center, top, bottom"
        )
    return cast(Literal["center", "top", "bottom"], normalized)


def _resolve_circular_plot_title_position(
    position: str,
) -> Literal["none", "top", "bottom"]:
    normalized = str(position).strip().lower()
    if normalized not in _CIRCULAR_PLOT_TITLE_POSITIONS:
        raise ValidationError(
            "plot_title_position must be one of: none, top, bottom"
        )
    return cast(Literal["none", "top", "bottom"], normalized)


def _resolve_plot_title_font_size(value: float | None) -> float:
    if value is None:
        return float(_PLOT_TITLE_DEFAULT_FONT_SIZE)
    if not math.isfinite(float(value)) or float(value) <= 0.0:
        raise ValidationError("plot_title_font_size must be a finite number > 0")
    return float(value)


def _normalize_plot_title(plot_title: str | None) -> str:
    return str(plot_title or "").strip()


def _circular_record_definition_profile(
    *,
    record_count: int,
    show_plot_title: bool,
    keep_full_definition_with_plot_title: bool,
) -> Literal["full", "record_summary"]:
    """The center definition a record draws: the full one, or the record summary.

    A shown plot title replaces the species lines unless the full definition is
    kept; without a title only a single record draws the full definition.
    """

    if show_plot_title:
        return "full" if keep_full_definition_with_plot_title else "record_summary"
    return "full" if record_count == 1 else "record_summary"


def _apply_circular_plot_title_font_size_override(
    *,
    cfg: GbdrawConfig,
    plot_title_font_size: float | None,
) -> GbdrawConfig:
    if plot_title_font_size is None:
        return cfg
    resolved_font_size = _resolve_plot_title_font_size(plot_title_font_size)
    return replace(
        cfg,
        objects=replace(
            cfg.objects,
            definition=replace(
                cfg.objects.definition,
                circular=replace(
                    cfg.objects.definition.circular,
                    plot_title_font_size=resolved_font_size,
                ),
            ),
        ),
    )


def _circular_plot_title_content(
    *,
    gb_record: SeqRecord,
    output_prefix: str,
    cfg: GbdrawConfig,
    profile: CircularRenderProfile,
    species: str | None,
    strain: str | None,
    plot_title: str | None,
    record_presentation: RecordPresentation | None,
) -> GroupContent:
    """The plot title content, keyed by the record's facts (not the record) and the title geometry."""

    plot_title_canvas_config = CircularCanvasConfigurator(
        output_prefix=output_prefix,
        profile=profile,
        legend="none",
        gb_record=gb_record,
    )
    return render_stage(
        "circular_plot_title_content",
        build_group_content,
        DefinitionGroup,
        definition_record_facts(gb_record, gc_percent=False),
        plot_title_canvas_config.geometry(),
        style=CircularDefinitionStyle.of(cfg, "shared_common"),
        species=species,
        strain=strain,
        plot_title=plot_title,
        definition_group_id="plot_title",
        record_presentation=record_presentation,
    )


def _validate_multi_record_min_radius_ratio(value: float) -> float:
    """Validate minimum radius ratio for multi-record scaling."""
    ratio = float(value)
    if not (0.0 < ratio <= 1.0):
        raise ValidationError("multi_record_min_radius_ratio must be > 0 and <= 1")
    return ratio


def _validate_multi_record_row_gap_ratio(value: float) -> float:
    """Validate row gap ratio for multi-record grid spacing."""
    ratio = float(value)
    if not math.isfinite(ratio) or ratio < 0.0:
        raise ValidationError("multi_record_row_gap_ratio must be a finite number >= 0")
    return ratio


def _validate_multi_record_column_gap_ratio(value: float) -> float:
    """Validate column gap ratio for multi-record grid spacing."""
    ratio = float(value)
    if not math.isfinite(ratio) or ratio < 0.0:
        raise ValidationError("multi_record_column_gap_ratio must be a finite number >= 0")
    return ratio


def _resolve_multi_record_scales(
    record_lengths: Sequence[int],
    *,
    mode: Literal["linear", "auto", "equal"],
    min_radius_ratio: float,
) -> list[float]:
    """Return per-record circular scales for multi-record canvas rendering."""
    if not record_lengths:
        return []
    if mode == "equal":
        return [1.0] * len(record_lengths)

    max_record_length = max(int(length) for length in record_lengths)
    if max_record_length <= 0:
        ratios = [1.0] * len(record_lengths)
    else:
        ratios = [
            max(0.0, float(length) / float(max_record_length))
            for length in record_lengths
        ]

    if mode == "linear":
        raw_scales = [float(ratio) for ratio in ratios]
    else:
        base_scales = [math.sqrt(float(ratio)) for ratio in ratios]
        below_min = [float(scale) for scale in base_scales if float(scale) < float(min_radius_ratio)]
        if len(below_min) >= 2:
            base_min = min(below_min)
            denominator = 1.0 - float(base_min)
            if denominator > 0.0:
                raw_scales = [
                    float(min_radius_ratio)
                    + (
                        (float(scale) - float(base_min))
                        * (1.0 - float(min_radius_ratio))
                        / denominator
                    )
                    for scale in base_scales
                ]
            else:
                raw_scales = [float(min_radius_ratio)] * len(base_scales)
        else:
            raw_scales = [float(scale) for scale in base_scales]

    return [
        max(float(min_radius_ratio), min(float(scale), 1.0))
        for scale in raw_scales
    ]


def _resolve_multi_record_grid_gap_px(max_record_radius_px: float, *, gap_ratio: float) -> float:
    """Return inter-record gap as a ratio of the largest record radius."""
    radius_px = float(max_record_radius_px)
    if radius_px > 0:
        return radius_px * float(gap_ratio)
    return float(_MULTI_RECORD_GRID_GAP_PX)


def _harmonize_multi_record_circular_style_cfg(
    cfg: GbdrawConfig,
    *,
    record_lengths: Sequence[int],
) -> GbdrawConfig:
    """Give a canvas with short and long records one size class: long.

    Every size-class-dependent value of a record is read through
    ``labels.length_threshold.circular``, so the canvas config classes every
    record as long instead of copying long values over short ones.
    """
    threshold = int(cfg.labels.length_threshold.circular)
    size_classes = {determine_length_parameter(int(length), threshold) for length in record_lengths}
    if len(size_classes) < 2:
        return cfg
    one_class = replace(cfg.labels.length_threshold, circular=0)
    return replace(cfg, labels=replace(cfg.labels, length_threshold=one_class))


def _scale_circular_cfg(cfg: GbdrawConfig, *, scale: float) -> GbdrawConfig:
    """Scale circular base geometry while preserving other config sections."""
    circular_cfg = cfg.canvas.circular
    scaled_width = replace(
        circular_cfg.width,
        with_labels=max(1, int(round(float(circular_cfg.width.with_labels) * float(scale)))),
        without_labels=max(1, int(round(float(circular_cfg.width.without_labels) * float(scale)))),
    )
    scaled_circular_cfg = replace(
        circular_cfg,
        height=max(1, int(round(float(circular_cfg.height) * float(scale)))),
        radius=max(1.0, float(circular_cfg.radius) * float(scale)),
        width=scaled_width,
    )
    scaled_canvas = replace(cfg.canvas, circular=scaled_circular_cfg)
    return replace(cfg, canvas=scaled_canvas)


def _check_feature_colors(table: DataFrame | None, label: str, *, by_row: bool) -> None:
    """Check every color of a -t/-d table once, whether it came from a file or a DataFrame.

    Rows of the specific color table carry their position; default colors are
    named by feature type (their rows merge the palette with the user's file).
    """

    if table is None or table.empty:
        return
    for row, (feature_type, color) in enumerate(zip(table["feature_type"], table["color"]), start=1):
        check_user_color(
            color,
            where=f"for feature type {feature_type!r} in the {label}" + (f" (row {row})" if by_row else ""),
            diagnostic={"code": "TABLE_INVALID", "field": "color", "reason": "COLOR", **({"row": row} if by_row else {})},
        )


def _check_caption_colors(color_table: DataFrame | None) -> None:
    """Check the colors of a color table with captions.

    A caption groups legend rows by color (``gbdraw.legend.table``), so every
    color of such a table is none, an SVG color name, #RGB, or #RRGGBB.
    """

    if color_table is None or color_table.empty or "caption" not in color_table:
        return
    for row, color in enumerate((str(value).strip() for value in color_table["color"]), start=1):
        if color.lower() == "none":
            continue
        try:
            normalize_hex_color(resolve_color_to_hex(color))
        except (ValueError, ValidationError) as exc:
            raise ValidationError(
                f"Invalid color {color!r} in the color table (row {row}). "
                "A color table with captions uses none, an SVG color name, #RGB, or #RRGGBB.",
                diagnostic={"code": "TABLE_INVALID", "field": "color", "reason": "COLOR", "row": row},
            ) from exc


def resolve_default_colors(
    default_colors: DataFrame | None,
    overrides: DataFrame | None,
    palette: str,
    load_comparison: bool,
) -> tuple[DataFrame, dict[str, str]]:
    """Return the default-color table and its map.

    A given table is used as it is; otherwise the packaged palette with the
    override rows applied. The ``default_colors`` render stage owns this, so
    the table's key names the overrides, the palette and ``load_comparison``.
    """

    if default_colors is None:
        default_colors = merge_default_colors(
            overrides,
            palette=palette or "default",
            load_comparison=load_comparison,
        )
    return default_colors, build_default_color_map(default_colors)


def resolve_feature_inputs(
    *,
    color_table: DataFrame | None,
    default_colors: DataFrame | None,
    feature_visibility_table: DataFrame | None,
    default_color_overrides: DataFrame | None = None,
    default_colors_palette: str = "default",
    load_comparison_colors: bool = False,
) -> ResolvedFeatureInputs:
    """Compile already-loaded feature inputs into one reusable value.

    Three render stages (design 2.6, M13): default colors (the given table, or
    the palette with the override rows), color rules from the color table,
    visibility rules from the visibility table. Errors come in that order.
    """

    default_colors, default_color_map = render_stage(
        "default_colors",
        resolve_default_colors,
        default_colors,
        default_color_overrides,
        default_colors_palette,
        load_comparison_colors,
    )
    _check_feature_colors(color_table, "color table", by_row=True)
    _check_feature_colors(default_colors, "default colors", by_row=False)
    _check_caption_colors(color_table)
    specific_color_rules = render_stage("color_rules", build_specific_color_rules, color_table)
    return ResolvedFeatureInputs(
        color_table=color_table,
        default_colors=default_colors,
        feature_visibility_table=feature_visibility_table,
        feature_visibility_rules=render_stage(
            "visibility_rules",
            compile_feature_visibility_rules,
            feature_visibility_table,
        ),
        specific_color_rules=specific_color_rules,
        default_color_map=default_color_map,
    )


def _resolve_feature_render_inputs(
    *,
    color_table: DataFrame | None,
    color_table_file: str | None,
    default_colors: DataFrame | None,
    default_colors_palette: str,
    default_colors_file: str | None,
    feature_visibility_table: DataFrame | None,
    feature_visibility_table_file: str | None,
    load_comparison_colors: bool,
    resolved_inputs: ResolvedFeatureInputs | None,
) -> ResolvedFeatureInputs:
    """Resolve direct-builder feature inputs or reuse a request planner result."""

    if resolved_inputs is not None:
        if not isinstance(resolved_inputs, ResolvedFeatureInputs):
            raise ValidationError(
                "_resolved_feature_inputs must be ResolvedFeatureInputs or None"
            )
        return resolved_inputs
    if color_table is None and color_table_file is not None:
        color_table = read_color_table(color_table_file)
    if (
        feature_visibility_table is None
        and feature_visibility_table_file is not None
    ):
        feature_visibility_table = read_feature_visibility_file(
            feature_visibility_table_file
        )
    return resolve_feature_inputs(
        color_table=color_table,
        default_colors=default_colors,
        feature_visibility_table=feature_visibility_table,
        default_color_overrides=(
            read_default_color_overrides(default_colors_file)
            if default_colors is None and default_colors_file
            else None
        ),
        default_colors_palette=default_colors_palette,
        load_comparison_colors=load_comparison_colors,
    )


def _feature_inputs_read_once(
    options: CircularDiagramOptions | LinearDiagramOptions,
    *,
    load_comparison_colors: bool,
    resolved_inputs: ResolvedFeatureInputs | None,
) -> ResolvedFeatureInputs:
    """Read a builder's color and visibility files once, for every assembly."""

    colors = options.colors or ColorOptions()
    return _resolve_feature_render_inputs(
        color_table=colors.color_table,
        color_table_file=colors.color_table_file,
        default_colors=colors.default_colors,
        default_colors_palette=colors.default_colors_palette,
        default_colors_file=colors.default_colors_file,
        feature_visibility_table=options.feature_visibility_table,
        feature_visibility_table_file=options.feature_visibility_table_file,
        load_comparison_colors=load_comparison_colors,
        resolved_inputs=resolved_inputs,
    )


def _linear_has_comparisons(
    *,
    blast_files: Sequence[str] | None,
    linear_comparisons: Sequence[LinearComparison] | None,
    protein_comparisons: Sequence[DataFrame] | None,
    collinearity_blocks: CollinearityResult | Sequence[CollinearityBlock] | None,
    losat_search: LosatSearchOptions | None,
) -> bool:
    return bool(
        blast_files
        or linear_comparisons
        or protein_comparisons
        or collinearity_blocks
        or losatp_analysis_mode(losat_search) != "none"
    )


def _web_safe_cache_filename(name: object, fallback: str = "losat") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name or "")).strip("_")
    return cleaned or fallback


def _web_normalize_cache_label(label: object, fallback: str) -> str:
    base = str(label or "").strip() or str(fallback or "")
    dotted = re.sub(r"\.+", ".", re.sub(r"[\s/]+", ".", base)).strip(".")
    return _web_safe_cache_filename(
        dotted,
        _web_safe_cache_filename(fallback, "losat"),
    )


def _linear_losat_cache_filenames(
    records: Sequence[SeqRecord],
    presentations: Sequence[RecordPresentation] | None,
    pairs: Sequence[tuple[int, int]] | None = None,
) -> tuple[str, ...]:
    """Raw TSV names of the displayed pairs (consecutive records by default)."""

    def label(index: int) -> str:
        record = records[index]
        presentation = presentations[index] if presentations is not None else None
        fallback = f"seq_{index + 1}"
        return _web_normalize_cache_label(
            str((presentation.label if presentation is not None else None) or "").strip()
            or str(record.id or "").strip()
            or fallback,
            fallback,
        )

    displayed = (
        pairs
        if pairs is not None
        else tuple((index, index + 1) for index in range(max(0, len(records) - 1)))
    )
    return tuple(
        f"{label(query)}.{label(subject)}.losatp.tsv" for query, subject in displayed
    )


def project_collinearity_frames(
    result: CollinearityResult,
    records: Sequence[SeqRecord],
    *,
    color_mode: CollinearityColorMode,
    search_scope: CollinearitySearchScope,
    comparison_pairs: Sequence[tuple[int, int]] | None,
) -> ProteinComparisonFrames:
    """The comparisons a Linear diagram draws for collinearity blocks.

    Without ``comparison_pairs`` the blocks give one ordered frame per adjacent
    record pair; with them, one comparison per listed pair.
    """

    if comparison_pairs is None:
        return ProteinComparisonFrames(
            ordered=tuple(
                convert_collinearity_blocks_to_comparisons(
                    result,
                    records=records,
                    color_mode=color_mode,
                    search_scope=search_scope,
                )
            )
        )
    pair_comparisons = convert_collinearity_blocks_to_pair_comparisons(
        result,
        records=records,
        color_mode=color_mode,
        search_scope=search_scope,
    )
    allowed_pairs = set(comparison_pairs)
    return ProteinComparisonFrames(
        pairs=tuple(
            LinearComparison(query_index, subject_index, matches)
            for (query_index, subject_index), matches in pair_comparisons.items()
            if (query_index, subject_index) in allowed_pairs
        )
    )


@dataclass(frozen=True)
class ProteinComparisonAnalysis:
    """A protein analysis, the comparisons it draws, and its raw-cache writes.

    ``comparisons`` is a direct field, so the stages that draw it can key it
    by this stage. The caller replays ``cache_writes`` (writes into a fork of
    the request's raw cache) on its cache on a build and on a hit, since
    ``session_entries()`` reads the display marks after the drawing.
    """

    result: ProteinBlastpResult | CollinearityResult
    comparisons: ProteinComparisonFrames
    cache_writes: tuple[LosatCacheWrite, ...]


def _invoke_protein_analysis_helper(
    mode: ProteinBlastpMode,
    records: Sequence[SeqRecord],
    *,
    losatp_bin: str,
    ncbi_blastp_bin: str | None,
    losatp_threads: int | None,
    pairwise_max_hits: int,
    candidate_limit: int | None,
    orthogroup_membership_mode: OrthogroupMembershipMode,
    orthogroup_member_max_hits: int | None,
    max_paralog_links_per_orthogroup: int,
    evalue: float,
    bitscore: float,
    identity: float,
    alignment_length: int,
    collinearity_params: LosslessCollinearityParameters | None,
    collinearity_unit_mode: CollinearityUnitMode | str,
    collinearity_anchor_mode: CollinearityAnchorMode,
    collinearity_search_scope: CollinearitySearchScope,
    collinearity_color_mode: CollinearityColorMode,
    collinear_infer_orthogroups: bool,
    collinearity_comparison_pairs: Sequence[tuple[int, int]] | None,
    cache_state: LosatpCacheState | None,
    protein_extraction: ProteinExtractionResult | None,
    feature_visibility_rules: list[dict[str, object]] | None,
    cache_filenames: Sequence[str] | None,
    pairwise_pairs: Sequence[tuple[int, int]] | None = None,
) -> ProteinComparisonAnalysis:
    """The ``protein_analysis`` stage: one analysis-helper call on a fork of the raw cache."""

    losatp_cache = (
        LosatpCacheManager.from_state(cache_state) if cache_state is not None else None
    )
    result: ProteinBlastpResult | CollinearityResult
    comparisons: ProteinComparisonFrames
    if mode == "pairwise":
        result = build_pairwise_protein_blastp_comparisons(
            records,
            pairs=pairwise_pairs,
            losatp_bin=losatp_bin,
            ncbi_blastp_bin=ncbi_blastp_bin,
            losatp_threads=losatp_threads,
            max_hits=pairwise_max_hits,
            candidate_limit=candidate_limit,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
            losatp_cache=losatp_cache,
            protein_extraction=protein_extraction,
            feature_visibility_rules=feature_visibility_rules,
            cache_filenames=cache_filenames,
        )
        comparisons = (
            ProteinComparisonFrames(
                pairs=tuple(
                    LinearComparison(query_index, subject_index, matches)
                    for (query_index, subject_index), matches in zip(
                        pairwise_pairs, result.comparisons, strict=True
                    )
                )
            )
            if pairwise_pairs is not None
            else ProteinComparisonFrames(ordered=tuple(result.comparisons))
        )
    elif mode == "orthogroup":
        result = build_rbh_orthogroup_protein_blastp_comparisons(
            records,
            losatp_bin=losatp_bin,
            ncbi_blastp_bin=ncbi_blastp_bin,
            losatp_threads=losatp_threads,
            candidate_limit=candidate_limit,
            orthogroup_membership_mode=orthogroup_membership_mode,
            orthogroup_member_max_hits=orthogroup_member_max_hits,
            max_related_edges_per_orthogroup=max_paralog_links_per_orthogroup,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
            losatp_cache=losatp_cache,
            protein_extraction=protein_extraction,
            feature_visibility_rules=feature_visibility_rules,
            cache_filenames=cache_filenames,
        )
        comparisons = ProteinComparisonFrames(ordered=tuple(result.comparisons))
    elif mode == "collinear":
        result = build_orthogroup_collinearity_blocks(
            records,
            infer_orthogroups=collinear_infer_orthogroups,
            losatp_bin=losatp_bin,
            ncbi_blastp_bin=ncbi_blastp_bin,
            losatp_threads=losatp_threads,
            candidate_limit=candidate_limit,
            orthogroup_membership_mode=orthogroup_membership_mode,
            orthogroup_member_max_hits=orthogroup_member_max_hits,
            max_paralog_links_per_orthogroup=max_paralog_links_per_orthogroup,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
            params=collinearity_params,
            unit_mode=collinearity_unit_mode,
            edge_mode=collinearity_anchor_mode,
            search_scope=collinearity_search_scope,
            comparison_pairs=collinearity_comparison_pairs,
            losatp_cache=losatp_cache,
            protein_extraction=protein_extraction,
            feature_visibility_rules=feature_visibility_rules,
            cache_filenames=cache_filenames,
        )
        comparisons = project_collinearity_frames(
            result,
            records,
            color_mode=collinearity_color_mode,
            search_scope=collinearity_search_scope,
            comparison_pairs=collinearity_comparison_pairs,
        )
    else:
        raise ValidationError(f"Unsupported protein analysis mode: {mode!r}.")
    return ProteinComparisonAnalysis(
        result, comparisons, losatp_cache.writes if losatp_cache is not None else ()
    )


def _shows_identity_labels(record_features: Sequence[ResolvedRecordFeatureInputs]) -> bool:
    """Whether an identity-addressed edit shows a label outside the display scope."""
    return any(
        override.label_visibility == "on"
        for inputs in record_features
        for override in inputs.overrides.values()
    )


def assemble_linear_diagram_from_records(
    records: Sequence[SeqRecord],
    *,
    cfg: GbdrawConfig,
    blast_files: Optional[Sequence[str]] = None,
    linear_comparisons: Sequence[LinearComparison] | None = None,
    layout: LinearMultiRecordOptions | None = None,
    protein_comparisons: Sequence[DataFrame] | None = None,
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None = None,
    losat_search: LosatSearchOptions | None = None,
    pairwise_match_style: Literal["ribbon", "curve"] | str = "ribbon",
    collinearity_blocks: CollinearityResult | Sequence[CollinearityBlock] | None = None,
    collinearity_params: LosslessCollinearityParameters | None = None,
    collinearity_unit_mode: CollinearityUnitMode | str = "auto",
    collinearity_anchor_mode: CollinearityAnchorMode | str = "rbh",
    collinearity_search_scope: CollinearitySearchScope | str = "adjacent",
    collinearity_color_mode: CollinearityColorMode | str = "orientation",
    losatp_cache: LosatpCacheManager | None = None,
    protein_extraction: ProteinExtractionResult | None = None,
    orthogroup_membership_mode: OrthogroupMembershipMode | str = "anchor_core_v1",
    collinear_infer_orthogroups: bool = True,
    collinear_max_paralog_links_per_orthogroup: int = 2,
    similarity_alignment: SimilarityAlignmentPlan | None = None,
    color_table: Optional[DataFrame] = None,
    color_table_file: str | None = None,
    default_colors: DataFrame | None = None,
    default_colors_palette: str = "default",
    default_colors_file: str | None = None,
    selected_features_set: Sequence[str] | None = None,
    feature_visibility_table: DataFrame | None = None,
    feature_visibility_table_file: str | None = None,
    feature_shapes: Mapping[str, str] | None = None,
    output_prefix: str = "out",
    legend: str = "right",
    dinucleotide: str = "GC",
    window: Optional[int] = None,
    step: Optional[int] = None,
    depth_window: Optional[int] = None,
    depth_step: Optional[int] = None,
    depth_tracks: Sequence[DepthTrackInput] | None = None,
    depth_table: DataFrame | None = None,
    depth_file: str | None = None,
    depth_tables: Sequence[DataFrame] | None = None,
    depth_files: Sequence[str] | None = None,
    depth_track_tables: Sequence[Sequence[DataFrame | None]] | None = None,
    depth_track_files: Sequence[Sequence[str | None]] | None = None,
    depth_track_labels: Sequence[str] | None = None,
    depth_track_colors: Sequence[str] | None = None,
    depth_track_heights: Sequence[float | str | None] | None = None,
    depth_track_large_tick_intervals: Sequence[float | str | None] | None = None,
    depth_track_small_tick_intervals: Sequence[float | str | None] | None = None,
    depth_track_tick_font_sizes: Sequence[float | str | None] | None = None,
    linear_track_slots: Sequence[str | LinearTrackSlot] | None = None,
    linear_track_axis_index: int | None = None,
    annotation_options: AnnotationOptions | None = None,
    plot_title: str | None = None,
    plot_title_position: Literal["center", "top", "bottom"] = "bottom",
    plot_title_font_size: float | None = None,
    evalue: float = LINEAR_MODE_PROFILE.comparison.evalue,
    bitscore: float = LINEAR_MODE_PROFILE.comparison.bitscore,
    identity: float = LINEAR_MODE_PROFILE.comparison.identity,
    alignment_length: int = LINEAR_MODE_PROFILE.comparison.alignment_length,
    _resolved_feature_inputs: ResolvedFeatureInputs | None = None,
    _resolved_record_features: tuple[ResolvedRecordFeatureInputs, ...] = (),
    _return_build_result: bool = False,
    _record_transforms: Sequence[RecordDisplayTransform] | None = None,
    _record_presentations: Sequence[RecordPresentation] | None = None,
    _alignment_anchor_centers: Sequence[float | None] | None = None,
    _resolved_annotations: ResolvedAnnotationBundle | None = None,
    _feature_label_rows: FeatureLabelRows | None = None,
) -> Drawing | LinearDiagramBuildResult:
    """Builds and assembles a linear diagram for the given records.

    This is a convenience wrapper that builds internal configurators/canvas objects and
    returns the assembled SVG canvas.
    If default_colors is None, it loads the built-in default palette.
    If color_table is None and color_table_file is provided, it is loaded.
    If selected_features_set is None, it uses the CLI default feature list.
    """
    if not isinstance(cfg, GbdrawConfig):
        raise ValidationError("cfg must be GbdrawConfig")
    if not records:
        raise ValidationError("records is empty")
    thresholds = ComparisonThresholds(
        evalue=evalue,
        bitscore=bitscore,
        identity=identity,
        alignment_length=alignment_length,
    )
    evalue = thresholds.evalue
    bitscore = thresholds.bitscore
    identity = thresholds.identity
    alignment_length = thresholds.alignment_length
    if layout is not None and not isinstance(layout, LinearMultiRecordOptions):
        raise ValidationError("layout must be LinearMultiRecordOptions or None")
    if linear_comparisons is not None and not all(
        isinstance(item, LinearComparison) for item in linear_comparisons
    ):
        raise ValidationError("linear_comparisons must contain LinearComparison values")
    if losat_search is not None and not isinstance(losat_search, LosatSearchOptions):
        raise ValidationError("losat_search must be LosatSearchOptions or None")
    if losat_search is not None and losat_search.program != "losatp":
        raise ValidationError(
            f"A {losat_search.program} search is resolved by the request planner; "
            "render it with gbdraw.api.render_request().",
            diagnostic={"code": "COMPARISON_INPUT", "reason": "LOSAT_PLAN"},
        )
    normalized_protein_blastp_mode = losatp_analysis_mode(losat_search)
    runtime = losat_search.runtime if losat_search is not None else LosatRuntimeOptions()
    normalized_protein_pairs: tuple[tuple[int, int], ...] | None = None
    if losat_search is not None and losat_search.pairs is not None:
        normalized_pairs: list[tuple[int, int]] = []
        for query_index, subject_index in losat_search.pairs:
            if query_index >= len(records) or subject_index >= len(records):
                raise ValidationError("pairs contains an out-of-range record index")
            if query_index == subject_index:
                raise ValidationError("pairs cannot compare a record to itself")
            normalized_pairs.append((query_index, subject_index))
        if len(set(normalized_pairs)) != len(normalized_pairs):
            raise ValidationError("pairs must not contain duplicates")
        normalized_protein_pairs = tuple(normalized_pairs)
        _ordered_indices, pair_rows = resolve_record_row_positions(
            records,
            layout.multi_record_positions if layout is not None else None,
        )
        for query_index, subject_index in normalized_protein_pairs:
            if abs(pair_rows[query_index] - pair_rows[subject_index]) != 1:
                raise ValidationError(
                    "pairs must connect records in adjacent rows: "
                    f"query=#{query_index + 1} row={pair_rows[query_index] + 1}, "
                    f"subject=#{subject_index + 1} row={pair_rows[subject_index] + 1}."
                )
    if normalized_protein_pairs is not None and linear_comparisons:
        raise ValidationError("Pass either losat_search pairs or linear_comparisons, not both")
    normalized_pairwise_match_style = normalize_pairwise_match_style(pairwise_match_style)
    normalized_collinearity_anchor_mode = normalize_collinearity_anchor_mode(
        str(collinearity_anchor_mode)
    )
    normalized_collinearity_search_scope = normalize_collinearity_search_scope(str(collinearity_search_scope))
    normalized_collinearity_color_mode = normalize_collinearity_color_mode(str(collinearity_color_mode))
    normalized_orthogroup_membership_mode = normalize_orthogroup_membership_mode(str(orthogroup_membership_mode))
    collinearity_comparison_pairs: tuple[tuple[int, int], ...] | None = None
    if (
        (normalized_protein_blastp_mode == "collinear" or collinearity_blocks is not None)
        and layout is not None
    ):
        _ordered_indices, collinearity_rows = resolve_record_row_positions(
            records,
            layout.multi_record_positions,
        )
        if len(set(collinearity_rows)) < len(records):
            collinearity_comparison_pairs = record_pairs_between_adjacent_rows(collinearity_rows)
    if int(collinear_max_paralog_links_per_orthogroup) <= 0:
        raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")
    if normalized_protein_blastp_mode != "none" and protein_comparisons is not None:
        raise ValidationError("Pass either a LOSATP losat_search or protein_comparisons, not both.")
    if collinearity_blocks is not None and (
        normalized_protein_blastp_mode != "none" or protein_comparisons is not None or blast_files
    ):
        raise ValidationError(
            "Pass collinearity_blocks without a LOSATP losat_search, protein_comparisons, or blast_files."
        )
    if normalized_protein_blastp_mode != "none" and blast_files:
        raise ValidationError("A LOSATP losat_search cannot be used with blast_files.")
    if normalized_protein_blastp_mode != "none" and len(records) < 2:
        raise ValidationError("A LOSATP losat_search requires at least two records")
    if similarity_alignment is not None and not isinstance(
        similarity_alignment, SimilarityAlignmentPlan
    ):
        raise ValidationError(
            "similarity_alignment must be SimilarityAlignmentPlan or None."
        )
    _validate_positive_optional("depth_window", depth_window)
    _validate_positive_optional("depth_step", depth_step)
    has_comparisons = _linear_has_comparisons(
        blast_files=blast_files,
        linear_comparisons=linear_comparisons,
        protein_comparisons=protein_comparisons,
        collinearity_blocks=collinearity_blocks,
        losat_search=losat_search,
    )
    resolved_feature_inputs = _resolve_feature_render_inputs(
        color_table=color_table,
        color_table_file=color_table_file,
        default_colors=default_colors,
        default_colors_palette=default_colors_palette,
        default_colors_file=default_colors_file,
        feature_visibility_table=feature_visibility_table,
        feature_visibility_table_file=feature_visibility_table_file,
        load_comparison_colors=has_comparisons,
        resolved_inputs=_resolved_feature_inputs,
    )
    color_table = resolved_feature_inputs.color_table
    default_colors = resolved_feature_inputs.default_colors
    feature_table = resolved_feature_inputs.feature_visibility_table
    feature_visibility_rules = resolved_feature_inputs.feature_visibility_rules

    if normalized_pairwise_match_style == "ribbon":
        normalized_pairwise_match_style = normalize_pairwise_match_style(
            str(cfg.objects.blast_match.style)
        )
    cfg = replace(
        cfg,
        objects=replace(
            cfg.objects,
            blast_match=replace(
                cfg.objects.blast_match,
                style=normalized_pairwise_match_style,
            ),
        ),
    )
    record_depth_tracks = normalize_depth_tracks(
        records,
        depth_tracks=depth_tracks,
        depth_table=depth_table,
        depth_file=depth_file,
        depth_tables=depth_tables,
        depth_files=depth_files,
        depth_track_tables=depth_track_tables,
        depth_track_files=depth_track_files,
        depth_track_labels=depth_track_labels,
        depth_track_colors=depth_track_colors,
        depth_track_heights=depth_track_heights,
        depth_track_large_tick_intervals=depth_track_large_tick_intervals,
        depth_track_small_tick_intervals=depth_track_small_tick_intervals,
        depth_track_tick_font_sizes=depth_track_tick_font_sizes,
    )
    parsed_linear_track_slots = _parse_linear_track_slot_inputs(linear_track_slots)
    resolved_linear_track_axis_index = _validate_linear_track_axis_index(
        linear_track_axis_index,
        parsed_linear_track_slots,
    )
    available_depth_track_count = depth_track_count(record_depth_tracks)
    if parsed_linear_track_slots is not None:
        show_depth = _slots_have_renderer(parsed_linear_track_slots, "depth")
        show_gc = _slots_have_renderer(parsed_linear_track_slots, "dinucleotide_content")
        show_skew = _slots_have_renderer(parsed_linear_track_slots, "dinucleotide_skew")
        if show_depth and record_depth_tracks is None:
            raise ValidationError("A linear depth track slot requires a depth_table, depth_file, or depth_track input.")
        if show_depth:
            _validate_depth_track_indices(
                parsed_linear_track_slots,
                depth_track_count_value=max(1, available_depth_track_count),
                available_indices=depth_track_indices(record_depth_tracks),
            )
        cfg = replace(
            cfg,
            canvas=replace(
                cfg.canvas,
                show_depth=show_depth,
                show_gc=show_gc,
                show_skew=show_skew,
            ),
        )
    else:
        if cfg.canvas.show_depth and record_depth_tracks is None:
            raise ValidationError("show_depth requires a depth_table, depth_file, or depth_track input.")
        show_depth = record_depth_tracks is not None
        if show_depth != bool(cfg.canvas.show_depth):
            cfg = replace(cfg, canvas=replace(cfg.canvas, show_depth=show_depth))
    _validate_gc_content_config(cfg.objects.gc_content)
    _validate_depth_config(cfg.objects.depth)
    profile = LinearRenderProfile(
        cfg, identity_label_overrides=_shows_identity_labels(_resolved_record_features)
    )

    if selected_features_set is None:
        selected_features_set = DEFAULT_SELECTED_FEATURES
    resolved_protein_frames: ProteinComparisonFrames | None = None
    resolved_orthogroups: OrthogroupResult | OrthogroupGraphResult | None = orthogroups
    resolved_collinearity_result: CollinearityResult | None = None
    losat_cache_filenames = _linear_losat_cache_filenames(records, _record_presentations)

    def invoke_protein_analysis(
        mode: ProteinBlastpMode,
        analysis_records: Sequence[SeqRecord],
        *,
        analysis_cache_filenames: Sequence[str] | None = losat_cache_filenames,
        pairwise_pairs: Sequence[tuple[int, int]] | None = None,
    ) -> ProteinComparisonAnalysis:
        analysis = render_stage(
            "protein_analysis",
            _invoke_protein_analysis_helper,
            mode,
            analysis_records,
            pairwise_pairs=pairwise_pairs,
            losatp_bin=runtime.losat_executable or AUTOMATIC_LOSAT_BIN,
            ncbi_blastp_bin=runtime.ncbi_blast_executable,
            losatp_threads=runtime.threads,
            pairwise_max_hits=(
                losat_search.losatp_max_hits if losat_search is not None else 5
            ),
            candidate_limit=(
                losat_search.losatp_max_target_seqs if losat_search is not None else None
            ),
            orthogroup_membership_mode=normalized_orthogroup_membership_mode,
            orthogroup_member_max_hits=(
                losat_search.losatp_member_max_hits if losat_search is not None else None
            ),
            max_paralog_links_per_orthogroup=int(
                collinear_max_paralog_links_per_orthogroup
            ),
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
            collinearity_params=collinearity_params,
            collinearity_unit_mode=collinearity_unit_mode,
            collinearity_anchor_mode=normalized_collinearity_anchor_mode,
            collinearity_search_scope=normalized_collinearity_search_scope,
            collinearity_color_mode=normalized_collinearity_color_mode,
            collinear_infer_orthogroups=collinear_infer_orthogroups,
            collinearity_comparison_pairs=collinearity_comparison_pairs,
            cache_state=losatp_cache.state() if losatp_cache is not None else None,
            protein_extraction=protein_extraction,
            feature_visibility_rules=feature_visibility_rules,
            cache_filenames=analysis_cache_filenames,
        )
        if losatp_cache is not None:
            losatp_cache.apply_writes(analysis.cache_writes)
        return analysis

    if protein_comparisons is not None:
        resolved_protein_frames = ProteinComparisonFrames(ordered=tuple(protein_comparisons))
    elif collinearity_blocks is not None:
        if isinstance(collinearity_blocks, CollinearityResult):
            collinearity_result = collinearity_blocks
        else:
            collinearity_result = CollinearityResult(blocks=tuple(collinearity_blocks))
        resolved_collinearity_result = collinearity_result
        if collinearity_result.orthogroups is not None:
            resolved_orthogroups = collinearity_result.orthogroups
        resolved_protein_frames = render_stage(
            "collinearity_projection",
            project_collinearity_frames,
            collinearity_result,
            tuple(records),
            color_mode=normalized_collinearity_color_mode,
            search_scope=normalized_collinearity_search_scope,
            comparison_pairs=collinearity_comparison_pairs,
        )
    elif normalized_protein_blastp_mode == "pairwise" and normalized_protein_pairs is not None:
        # One plan for every pair, so records of one file share a database.
        resolved_protein_frames = invoke_protein_analysis(
            "pairwise",
            records,
            analysis_cache_filenames=_linear_losat_cache_filenames(
                records, _record_presentations, normalized_protein_pairs
            ),
            pairwise_pairs=normalized_protein_pairs,
        ).comparisons
    elif normalized_protein_blastp_mode == "pairwise":
        resolved_protein_frames = invoke_protein_analysis("pairwise", records).comparisons
    elif normalized_protein_blastp_mode == "orthogroup":
        analysis = invoke_protein_analysis("orthogroup", records)
        resolved_protein_frames = analysis.comparisons
        resolved_orthogroups = cast(ProteinBlastpResult, analysis.result).orthogroups
    elif normalized_protein_blastp_mode == "collinear":
        analysis = invoke_protein_analysis("collinear", records)
        collinearity_result = cast(CollinearityResult, analysis.result)
        resolved_collinearity_result = collinearity_result
        resolved_orthogroups = collinearity_result.orthogroups
        resolved_protein_frames = analysis.comparisons
    # Ordered frames and explicit pairs, as the reply metadata and the canvas read them.
    resolved_protein_comparisons = (
        resolved_protein_frames.ordered if resolved_protein_frames is not None else None
    )
    resolved_linear_comparisons = (
        *(linear_comparisons or ()),
        *(resolved_protein_frames.pairs if resolved_protein_frames is not None else ()),
    )
    normalized_plot_title = str(plot_title or "").strip()
    normalized_plot_title_position = _resolve_linear_plot_title_position(
        str(plot_title_position)
    )
    resolved_plot_title_font_size = _resolve_plot_title_font_size(plot_title_font_size)

    seq_len_dict = create_dict_for_sequence_lengths(records)
    # Use raw records to avoid collapsing lengths when IDs are duplicated.
    longest_genome = max(len(record.seq) for record in records)

    # Match legacy CLI behavior: linear window/step are based on the longest genome.
    auto_window, auto_step = sliding_window(longest_genome, cfg)
    window = auto_window if window is None else window
    step = auto_step if step is None else step
    resolved_depth_window, resolved_depth_step = _resolve_depth_window_step(
        window=int(window),
        step=int(step),
        depth_window=depth_window,
        depth_step=depth_step,
    )

    blast_config = BlastMatchConfigurator(
        evalue=evalue,
        bitscore=bitscore,
        identity=identity,
        alignment_length=alignment_length,
        sequence_length_dict=seq_len_dict,
        profile=profile,
        default_colors_df=default_colors,
    )
    additional_color_modes = (
        {normalized_collinearity_color_mode}
        if collinearity_blocks is not None or normalized_protein_blastp_mode == "collinear"
        else None
    )
    protein_identity_legend = pairwise_identity_legend_from_comparisons(
        blast_config,
        resolved_protein_comparisons,
        additional_color_modes=additional_color_modes,
    )

    canvas_config = LinearCanvasConfigurator(
        num_of_entries=len(records),
        longest_genome=longest_genome,
        profile=profile,
        legend=legend,
        output_prefix=output_prefix,
        has_comparisons=bool(blast_files or resolved_linear_comparisons or resolved_protein_comparisons),
        depth_track_count=max(1, depth_track_count(record_depth_tracks)),
        depth_track_heights=_depth_track_heights_from_specs(record_depth_tracks),
    )
    feature_config = FeatureDrawingConfigurator(
        color_table=color_table,
        default_colors=default_colors,
        selected_features_set=list(selected_features_set),
        profile=profile,
        feature_table=feature_table,
        feature_shapes=feature_shapes,
        feature_visibility_rules=feature_visibility_rules,
        specific_color_rules=resolved_feature_inputs.specific_color_rules,
        default_color_map=resolved_feature_inputs.default_color_map,
        canvas_config=canvas_config,
        record_features=_resolved_record_features,
    )
    gc_config = GcContentConfigurator(
        window=window,
        step=step,
        dinucleotide=dinucleotide,
        profile=profile,
        default_colors_df=default_colors,
    )
    skew_config = GcSkewConfigurator(
        window=window,
        step=step,
        dinucleotide=dinucleotide,
        profile=profile,
        default_colors_df=default_colors,
    )
    depth_config = DepthConfigurator(
        window=resolved_depth_window,
        step=resolved_depth_step,
        profile=profile,
    ) if show_depth else None
    legend_config = LegendDrawingConfigurator(
        profile=profile,
        canvas_config=canvas_config,
    )

    canvas = assemble_linear_diagram(
        record_transforms=_record_transforms,
        record_presentations=_record_presentations,
        records=list(records),
        blast_files=list(blast_files) if blast_files else None,
        canvas_config=canvas_config,
        blast_config=blast_config,
        feature_config=feature_config,
        gc_config=gc_config,
        legend_config=legend_config,
        skew_config=skew_config,
        depth_config=depth_config,
        record_depth_tracks=record_depth_tracks,
        linear_track_slots=parsed_linear_track_slots,
        linear_track_axis_index=resolved_linear_track_axis_index,
        annotations=_resolved_annotations if _resolved_annotations is not None else annotation_options,
        plot_title=normalized_plot_title or None,
        plot_title_position=normalized_plot_title_position,
        plot_title_font_size=resolved_plot_title_font_size,
        protein_comparisons=resolved_protein_frames,
        pairwise_identity_legend=protein_identity_legend,
        linear_comparisons=list(linear_comparisons) if linear_comparisons else None,
        linear_layout=layout,
        orthogroups=resolved_orthogroups,
        similarity_alignment=similarity_alignment,
        alignment_anchor_centers=_alignment_anchor_centers,
        feature_label_rows=_feature_label_rows,
    )
    build_result = LinearDiagramBuildResult(
        drawing=canvas,
        metadata=LinearDiagramMetadata(
            protein_comparisons=resolved_protein_comparisons,
            linear_comparisons=resolved_linear_comparisons,
            orthogroups=resolved_orthogroups,
            collinearity_result=resolved_collinearity_result,
            comparison_warnings=getattr(canvas, "_gbdraw_comparison_record_id_warnings", ()),
        ),
    )
    return build_result if _return_build_result else canvas


def assemble_circular_diagram_from_record(
    gb_record: SeqRecord,
    *,
    cfg: GbdrawConfig,
    conservation_blast_files: Sequence[ConservationInput] | None = None,
    conservation_dataframes: Sequence[DataFrame] | None = None,
    conservation_reference: Literal["query", "subject", "auto"] | str = "auto",
    conservation_labels: Sequence[str] | None = None,
    conservation_colors: Sequence[str] | None = None,
    conservation_ring_width: float | None = None,
    conservation_ring_gap: float | None = None,
    color_table: Optional[DataFrame] = None,
    color_table_file: str | None = None,
    default_colors: DataFrame | None = None,
    default_colors_palette: str = "default",
    default_colors_file: str | None = None,
    selected_features_set: Sequence[str] | None = None,
    feature_visibility_table: DataFrame | None = None,
    feature_visibility_table_file: str | None = None,
    feature_shapes: Mapping[str, str] | None = None,
    output_prefix: str = "out",
    legend: str = "right",
    dinucleotide: str = "GC",
    window: Optional[int] = None,
    step: Optional[int] = None,
    depth_window: Optional[int] = None,
    depth_step: Optional[int] = None,
    depth_tracks: Sequence[DepthTrackInput] | None = None,
    depth_table: DataFrame | None = None,
    depth_file: str | None = None,
    depth_track_tables: Sequence[Sequence[DataFrame | None]] | None = None,
    depth_track_files: Sequence[Sequence[str | None]] | None = None,
    depth_track_labels: Sequence[str] | None = None,
    depth_track_colors: Sequence[str] | None = None,
    depth_track_large_tick_intervals: Sequence[float | str | None] | None = None,
    depth_track_small_tick_intervals: Sequence[float | str | None] | None = None,
    depth_track_tick_font_sizes: Sequence[float | str | None] | None = None,
    species: Optional[str] = None,
    strain: Optional[str] = None,
    plot_title: str | None = None,
    plot_title_position: Literal["none", "top", "bottom"] = "none",
    plot_title_font_size: float | None = None,
    keep_full_definition_with_plot_title: bool = False,
    center_reserved_radius: float | None = None,
    circular_track_slots: Sequence[str | CircularTrackSlot] | None = None,
    circular_track_axis_index: int | None = None,
    annotation_options: AnnotationOptions | None = None,
    evalue: float = CIRCULAR_MODE_PROFILE.comparison.evalue,
    bitscore: float = CIRCULAR_MODE_PROFILE.comparison.bitscore,
    identity: float = CIRCULAR_MODE_PROFILE.comparison.identity,
    alignment_length: int = CIRCULAR_MODE_PROFILE.comparison.alignment_length,
    _definition_profile: Literal["full", "record_summary", "shared_common"] = "full",
    _layout_definition_profile: CircularDefinitionProfile | None = None,
    _precomputed_depth_df: DataFrame | None = None,
    _precomputed_depth_track_specs: Sequence[DepthTrackSpec] | None = None,
    _precomputed_depth_tracks: Sequence[DepthTrackData] | None = None,
    _precomputed_depth_track_count: int | None = None,
    _shared_depth_max: float | None = None,
    _precomputed_conservation_tracks: Sequence[ConservationTrack] | None = None,
    _resolved_annotations: ResolvedAnnotationBundle | None = None,
    _annotation_record_index: int = 0,
    _definition_group_id: str | None = None,
    _resolved_feature_inputs: ResolvedFeatureInputs | None = None,
    _resolved_record_features: tuple[ResolvedRecordFeatureInputs, ...] = (),
    _record_transform: RecordDisplayTransform | None = None,
    _record_presentation: RecordPresentation | None = None,
) -> Drawing:
    """Builds and assembles a circular diagram for a single record.

    If default_colors is None, it loads the built-in default palette.
    If color_table is None and color_table_file is provided, it is loaded.
    If selected_features_set is None, it uses the CLI default feature list.
    The private ``_precomputed_depth_tracks`` is None when the canvas has no
    depth input; ``[]`` means depth input without a cell for this record, which
    keeps the record's depth slot.
    """
    if not isinstance(cfg, GbdrawConfig):
        raise ValidationError("cfg must be GbdrawConfig")
    thresholds = ComparisonThresholds(
        evalue=evalue,
        bitscore=bitscore,
        identity=identity,
        alignment_length=alignment_length,
    )
    evalue = thresholds.evalue
    bitscore = thresholds.bitscore
    identity = thresholds.identity
    alignment_length = thresholds.alignment_length
    _validate_positive_optional("depth_window", depth_window)
    _validate_positive_optional("depth_step", depth_step)
    _validate_positive_float_optional("conservation_ring_width", conservation_ring_width)
    _validate_positive_float_optional("conservation_ring_gap", conservation_ring_gap)
    _validate_nonnegative_float_optional("center_reserved_radius", center_reserved_radius)
    resolved_feature_inputs = _resolve_feature_render_inputs(
        color_table=color_table,
        color_table_file=color_table_file,
        default_colors=default_colors,
        default_colors_palette=default_colors_palette,
        default_colors_file=default_colors_file,
        feature_visibility_table=feature_visibility_table,
        feature_visibility_table_file=feature_visibility_table_file,
        load_comparison_colors=False,
        resolved_inputs=_resolved_feature_inputs,
    )
    color_table = resolved_feature_inputs.color_table
    default_colors = resolved_feature_inputs.default_colors
    feature_table = resolved_feature_inputs.feature_visibility_table

    cfg = _apply_circular_plot_title_font_size_override(
        cfg=cfg,
        plot_title_font_size=plot_title_font_size,
    )
    record_depth_tracks = (
        [list(_precomputed_depth_track_specs)]
        if _precomputed_depth_track_specs is not None
        else normalize_depth_tracks(
            [gb_record],
            depth_tracks=depth_tracks,
            depth_table=depth_table,
            depth_file=depth_file,
            depth_track_tables=depth_track_tables,
            depth_track_files=depth_track_files,
            depth_track_labels=depth_track_labels,
            depth_track_colors=depth_track_colors,
            depth_track_large_tick_intervals=depth_track_large_tick_intervals,
            depth_track_small_tick_intervals=depth_track_small_tick_intervals,
            depth_track_tick_font_sizes=depth_track_tick_font_sizes,
        )
    )
    precomputed_depth_tracks_provided = _precomputed_depth_tracks is not None
    precomputed_depth_track_list = list(_precomputed_depth_tracks or [])
    if (
        cfg.canvas.show_depth
        and record_depth_tracks is None
        and not precomputed_depth_tracks_provided
        and _precomputed_depth_df is None
    ):
        raise ValidationError("show_depth requires a depth_table, depth_file, or depth_track input.")
    show_depth_from_input = (
        record_depth_tracks is not None
        or precomputed_depth_tracks_provided
        or _precomputed_depth_df is not None
    )
    if show_depth_from_input != bool(cfg.canvas.show_depth):
        cfg = replace(cfg, canvas=replace(cfg.canvas, show_depth=show_depth_from_input))
    if cfg.objects.depth.share_axis and cfg.objects.depth.max_depth is None and _shared_depth_max is not None:
        cfg = _cfg_with_depth_scale_max(cfg, _shared_depth_max)
    _validate_gc_content_config(cfg.objects.gc_content)
    _validate_depth_config(cfg.objects.depth)

    if selected_features_set is None:
        selected_features_set = DEFAULT_SELECTED_FEATURES
    normalized_plot_title_position = _resolve_circular_plot_title_position(
        str(plot_title_position)
    )
    normalized_plot_title = _normalize_plot_title(plot_title)
    show_plot_title = False
    effective_definition_profile: Literal["full", "record_summary", "shared_common"] = (
        _definition_profile
    )
    if _definition_profile == "full":
        show_plot_title = normalized_plot_title_position != "none"
        effective_definition_profile = _circular_record_definition_profile(
            record_count=1,
            show_plot_title=show_plot_title,
            keep_full_definition_with_plot_title=keep_full_definition_with_plot_title,
        )

    parsed_circular_track_slots = _parse_circular_track_slot_inputs(circular_track_slots)
    circular_track_axis_index = _validate_circular_track_axis_index(
        circular_track_axis_index,
        parsed_circular_track_slots,
    )

    legend_effective = legend

    # Explicit slots override high-level show flags used by canvas sizing.
    show_depth, show_gc, show_skew = _resolve_circular_track_visibility(
        cfg,
        parsed_circular_track_slots,
    )
    if (
        parsed_circular_track_slots is not None
        and show_depth
        and not show_depth_from_input
    ):
        raise ValidationError("A circular depth track slot requires a depth_table, depth_file, or depth_track input.")
    available_depth_track_count = _depth_track_count_for_render(
        record_depth_tracks,
        precomputed_depth_track_list,
        _precomputed_depth_df,
    )
    if _precomputed_depth_track_count is not None:
        available_depth_track_count = int(_precomputed_depth_track_count)
    parsed_circular_track_slots = _default_circular_depth_slots_if_needed(
        parsed_circular_track_slots=parsed_circular_track_slots,
        show_depth=bool(show_depth and show_depth_from_input),
        depth_track_count_value=available_depth_track_count,
        show_ticks=bool(cfg.objects.scale.show),
        show_gc=bool(show_gc),
        show_skew=bool(show_skew),
        dinucleotide=dinucleotide,
    )
    canvas_cfg = cfg.canvas
    canvas_cfg = replace(
        canvas_cfg,
        show_depth=bool(show_depth and show_depth_from_input),
        show_gc=bool(show_gc),
        show_skew=bool(show_skew),
    )
    cfg = replace(cfg, canvas=canvas_cfg)
    profile = CircularRenderProfile(
        cfg, identity_label_overrides=_shows_identity_labels(_resolved_record_features)
    )

    conservation_mode = normalize_conservation_reference(conservation_reference)
    if _precomputed_conservation_tracks is not None:
        conservation_tracks = tuple(_precomputed_conservation_tracks)
    else:
        conservation_load_result = _load_conservation_result(
            conservation_blast_files=conservation_blast_files,
            conservation_dataframes=conservation_dataframes,
            conservation_labels=conservation_labels,
            conservation_colors=conservation_colors,
            thresholds=thresholds,
        )
        conservation_tracks = _conservation_tracks_for_record(
            conservation_load_result, [gb_record], gb_record, conservation_mode
        )
    if conservation_tracks:
        if _circular_slots_define_renderer(parsed_circular_track_slots, "sequence_conservation"):
            parsed_circular_track_slots = _apply_conservation_track_params_to_slots(
                parsed_circular_track_slots or (),
                conservation_tracks,
                cfg=cfg,
                conservation_ring_width=conservation_ring_width,
                conservation_ring_gap=conservation_ring_gap,
            )
        else:
            conservation_slots = _conservation_slots_for_tracks(
                conservation_tracks,
                cfg=cfg,
                conservation_ring_width=conservation_ring_width,
                conservation_ring_gap=conservation_ring_gap,
            )
            if parsed_circular_track_slots is None:
                parsed_circular_track_slots = circular_track_slots_from_order(
                    "features,ticks,depth,gc_content,gc_skew",
                    show_ticks=bool(cfg.objects.scale.show),
                    show_depth=show_depth,
                    depth_track_count=available_depth_track_count,
                    show_gc=show_gc,
                    show_skew=show_skew,
                    dinucleotide=dinucleotide,
                )
            parsed_circular_track_slots = _insert_conservation_slots(
                parsed_circular_track_slots,
                conservation_slots,
                axis_index=circular_track_axis_index,
            )

    # Match legacy CLI behavior: circular window/step are based on the record length.
    auto_window, auto_step = sliding_window(len(gb_record.seq), cfg)
    window = auto_window if window is None else window
    step = auto_step if step is None else step
    resolved_depth_window, resolved_depth_step = _resolve_depth_window_step(
        window=int(window),
        step=int(step),
        depth_window=depth_window,
        depth_step=depth_step,
    )

    gc_config = GcContentConfigurator(
        window=window,
        step=step,
        dinucleotide=dinucleotide,
        profile=profile,
        default_colors_df=default_colors,
    )
    skew_config = GcSkewConfigurator(
        window=window,
        step=step,
        dinucleotide=dinucleotide,
        profile=profile,
        default_colors_df=default_colors,
    )
    depth_config = (
        DepthConfigurator(
            window=resolved_depth_window,
            step=resolved_depth_step,
            profile=profile,
        )
        if cfg.canvas.show_depth
        else None
    )

    # Circular drawing expects precomputed dinucleotide dataframes, but only when needed.
    if parsed_circular_track_slots is not None:
        requested_nts = _dinucleotides_from_circular_slots(
            parsed_circular_track_slots,
            default_nt=dinucleotide,
        )
    else:
        requested_nts = {str(dinucleotide).upper()} if (cfg.canvas.show_gc or cfg.canvas.show_skew) else set()
    dinucleotide_content_dataframes = _build_circular_dinucleotide_content_dataframes(
        gb_record,
        window=int(window),
        step=int(step),
        nts=requested_nts,
    )
    dinucleotide_skew_dataframes = _build_circular_dinucleotide_skew_dataframes(
        gb_record,
        window=int(window),
        step=int(step),
        nts=requested_nts,
    )
    gc_df = dinucleotide_content_dataframes.get(str(dinucleotide).upper(), DataFrame())
    if cfg.canvas.show_depth and depth_config is not None:
        if precomputed_depth_track_list:
            resolved_depth_tracks = precomputed_depth_track_list
        elif _precomputed_depth_df is not None:
            resolved_depth_tracks = [
                DepthTrackData(
                    id="depth",
                    label="Depth",
                    df=_precomputed_depth_df,
                    config=depth_config,
                )
            ]
        elif record_depth_tracks is not None:
            resolved_depth_tracks = build_depth_track_dataframes(
                [gb_record],
                record_depth_tracks,
                base_config=depth_config,
                depth_df_builder=build_depth_df,
                window_steps=[(resolved_depth_window, resolved_depth_step)],
            )[0]
        else:
            resolved_depth_tracks = []
    else:
        resolved_depth_tracks = []
    resolved_depth_by_index = index_depth_track_row(resolved_depth_tracks)
    resolved_depth_track_zero = resolved_depth_by_index.get(0)
    resolved_depth_df = resolved_depth_track_zero.df if resolved_depth_track_zero is not None else None
    _validate_depth_track_indices(
        parsed_circular_track_slots,
        depth_track_count_value=max(1, available_depth_track_count),
        available_indices=(
            None
            if _precomputed_depth_track_count is not None
            else depth_track_indices([resolved_depth_tracks])
        ),
    )

    canvas_config = CircularCanvasConfigurator(
        output_prefix=output_prefix,
        profile=profile,
        legend=legend_effective,
        gb_record=gb_record,
    )
    feature_config = FeatureDrawingConfigurator(
        color_table=color_table,
        default_colors=default_colors,
        selected_features_set=list(selected_features_set),
        profile=profile,
        feature_table=feature_table,
        feature_shapes=feature_shapes,
        feature_visibility_rules=resolved_feature_inputs.feature_visibility_rules,
        specific_color_rules=resolved_feature_inputs.specific_color_rules,
        default_color_map=resolved_feature_inputs.default_color_map,
        canvas_config=canvas_config,
        record_features=_resolved_record_features,
    )
    legend_config = LegendDrawingConfigurator(
        profile=profile,
        canvas_config=canvas_config,
    )

    title_content: GroupContent | None = None
    title_placement = TitlePlacement.NONE
    if show_plot_title:
        title_content = _circular_plot_title_content(
            gb_record=gb_record,
            output_prefix=output_prefix,
            cfg=cfg,
            profile=profile,
            species=species,
            strain=strain,
            plot_title=normalized_plot_title or None,
            record_presentation=_record_presentation,
        )
        title_placement = TitlePlacement(normalized_plot_title_position)

    result = _assemble_circular_diagram_result(
        record_transform=_record_transform,
        record_presentation=_record_presentation,
        gb_record=gb_record,
        canvas_config=canvas_config,
        gc_df=gc_df,
        gc_config=gc_config,
        skew_config=skew_config,
        feature_config=feature_config,
        species=species,
        strain=strain,
        plot_title=(
            normalized_plot_title or None
            if effective_definition_profile == "shared_common"
            else None
        ),
        legend_config=legend_config,
        depth_df=resolved_depth_df,
        depth_config=depth_config,
        depth_tracks=resolved_depth_tracks,
        depth_track_count_value=available_depth_track_count,
        conservation_tracks=conservation_tracks,
        conservation_min_identity=float(identity),
        circular_track_slots=parsed_circular_track_slots,
        circular_track_axis_index=circular_track_axis_index,
        annotations=_resolved_annotations if _resolved_annotations is not None else annotation_options,
        annotation_record_index=_annotation_record_index,
        dinucleotide_content_dataframes=dinucleotide_content_dataframes,
        dinucleotide_skew_dataframes=dinucleotide_skew_dataframes,
        definition_position="center",
        definition_profile=effective_definition_profile,
        # The radial layout reserves the definition drawn without the plot title (OV-1710).
        layout_definition_profile=(
            _definition_profile if _layout_definition_profile is None else _layout_definition_profile
        ),
        definition_group_id=_definition_group_id,
        center_reserved_radius=center_reserved_radius,
        title_content=title_content,
        title_placement=title_placement,
    )
    _attach_circular_label_layouts(result.drawing, [result.label_layout])
    return result.drawing


def assemble_circular_diagram_from_records(
    records: Sequence[SeqRecord],
    *,
    cfg: GbdrawConfig,
    conservation_blast_files: Sequence[ConservationInput] | None = None,
    conservation_dataframes: Sequence[DataFrame] | None = None,
    conservation_reference: Literal["query", "subject", "auto"] | str = "auto",
    conservation_labels: Sequence[str] | None = None,
    conservation_colors: Sequence[str] | None = None,
    conservation_ring_width: float | None = None,
    conservation_ring_gap: float | None = None,
    color_table: Optional[DataFrame] = None,
    color_table_file: str | None = None,
    default_colors: DataFrame | None = None,
    default_colors_palette: str = "default",
    default_colors_file: str | None = None,
    selected_features_set: Sequence[str] | None = None,
    feature_visibility_table: DataFrame | None = None,
    feature_visibility_table_file: str | None = None,
    feature_shapes: Mapping[str, str] | None = None,
    output_prefix: str = "out",
    legend: str = "right",
    dinucleotide: str = "GC",
    window: Optional[int] = None,
    step: Optional[int] = None,
    depth_window: Optional[int] = None,
    depth_step: Optional[int] = None,
    depth_tracks: Sequence[DepthTrackInput] | None = None,
    depth_table: DataFrame | None = None,
    depth_file: str | None = None,
    depth_tables: Sequence[DataFrame] | None = None,
    depth_files: Sequence[str] | None = None,
    depth_track_tables: Sequence[Sequence[DataFrame | None]] | None = None,
    depth_track_files: Sequence[Sequence[str | None]] | None = None,
    depth_track_labels: Sequence[str] | None = None,
    depth_track_colors: Sequence[str] | None = None,
    depth_track_large_tick_intervals: Sequence[float | str | None] | None = None,
    depth_track_small_tick_intervals: Sequence[float | str | None] | None = None,
    depth_track_tick_font_sizes: Sequence[float | str | None] | None = None,
    species: Optional[str] = None,
    strain: Optional[str] = None,
    plot_title: str | None = None,
    plot_title_position: Literal["none", "top", "bottom"] = "none",
    plot_title_font_size: float | None = None,
    keep_full_definition_with_plot_title: bool = False,
    center_reserved_radius: float | None = None,
    multi_record_size_mode: Literal["linear", "auto", "equal"] = "auto",
    multi_record_min_radius_ratio: float = 0.55,
    multi_record_column_gap_ratio: float = _MULTI_RECORD_COLUMN_GAP_RATIO,
    multi_record_row_gap_ratio: float = _MULTI_RECORD_ROW_GAP_RATIO,
    multi_record_positions: Sequence[str] | None = None,
    circular_track_slots: Sequence[str | CircularTrackSlot] | None = None,
    circular_track_axis_index: int | None = None,
    annotation_options: AnnotationOptions | None = None,
    evalue: float = CIRCULAR_MODE_PROFILE.comparison.evalue,
    bitscore: float = CIRCULAR_MODE_PROFILE.comparison.bitscore,
    identity: float = CIRCULAR_MODE_PROFILE.comparison.identity,
    alignment_length: int = CIRCULAR_MODE_PROFILE.comparison.alignment_length,
    _resolved_feature_inputs: ResolvedFeatureInputs | None = None,
    _resolved_record_features: tuple[ResolvedRecordFeatureInputs, ...] = (),
    _record_transforms: Sequence[RecordDisplayTransform] | None = None,
    _record_presentations: Sequence[RecordPresentation] | None = None,
    _resolved_annotations: ResolvedAnnotationBundle | None = None,
) -> Drawing:
    """Build and assemble a circular diagram grid from multiple records."""
    if not isinstance(cfg, GbdrawConfig):
        raise ValidationError("cfg must be GbdrawConfig")
    if not records:
        raise ValidationError("records is empty")
    thresholds = ComparisonThresholds(
        evalue=evalue,
        bitscore=bitscore,
        identity=identity,
        alignment_length=alignment_length,
    )
    evalue = thresholds.evalue
    bitscore = thresholds.bitscore
    identity = thresholds.identity
    alignment_length = thresholds.alignment_length
    resolved_annotations = (
        _resolved_annotations if _resolved_annotations is not None else
        resolve_annotations(annotation_options, records, mode="circular", record_transforms=_record_transforms)
    )
    _validate_positive_optional("depth_window", depth_window)
    _validate_positive_optional("depth_step", depth_step)
    _validate_positive_float_optional("conservation_ring_width", conservation_ring_width)
    _validate_positive_float_optional("conservation_ring_gap", conservation_ring_gap)
    _validate_nonnegative_float_optional("center_reserved_radius", center_reserved_radius)

    normalized_multi_record_size_mode = _resolve_multi_record_size_mode(
        str(multi_record_size_mode)
    )
    normalized_multi_record_min_radius_ratio = _validate_multi_record_min_radius_ratio(
        float(multi_record_min_radius_ratio)
    )
    normalized_multi_record_column_gap_ratio = _validate_multi_record_column_gap_ratio(
        float(multi_record_column_gap_ratio)
    )
    normalized_multi_record_row_gap_ratio = _validate_multi_record_row_gap_ratio(
        float(multi_record_row_gap_ratio)
    )
    normalized_plot_title_position = _resolve_circular_plot_title_position(
        str(plot_title_position)
    )
    normalized_plot_title = _normalize_plot_title(plot_title)
    resolved_feature_inputs = _resolve_feature_render_inputs(
        color_table=color_table,
        color_table_file=color_table_file,
        default_colors=default_colors,
        default_colors_palette=default_colors_palette,
        default_colors_file=default_colors_file,
        feature_visibility_table=feature_visibility_table,
        feature_visibility_table_file=feature_visibility_table_file,
        load_comparison_colors=False,
        resolved_inputs=_resolved_feature_inputs,
    )
    color_table = resolved_feature_inputs.color_table
    default_colors = resolved_feature_inputs.default_colors
    feature_table = resolved_feature_inputs.feature_visibility_table

    cfg = _apply_circular_plot_title_font_size_override(
        cfg=cfg,
        plot_title_font_size=plot_title_font_size,
    )

    if selected_features_set is None:
        selected_features_set = DEFAULT_SELECTED_FEATURES

    records = list(records)
    record_depth_tracks = normalize_depth_tracks(
        records,
        depth_tracks=depth_tracks,
        depth_table=depth_table,
        depth_file=depth_file,
        depth_tables=depth_tables,
        depth_files=depth_files,
        depth_track_tables=depth_track_tables,
        depth_track_files=depth_track_files,
        depth_track_labels=depth_track_labels,
        depth_track_colors=depth_track_colors,
        depth_track_large_tick_intervals=depth_track_large_tick_intervals,
        depth_track_small_tick_intervals=depth_track_small_tick_intervals,
        depth_track_tick_font_sizes=depth_track_tick_font_sizes,
    )
    if cfg.canvas.show_depth and record_depth_tracks is None:
        raise ValidationError("show_depth requires depth_tables, depth_files, or depth_track input.")
    show_depth_from_input = record_depth_tracks is not None
    if show_depth_from_input != bool(cfg.canvas.show_depth):
        cfg = replace(cfg, canvas=replace(cfg.canvas, show_depth=show_depth_from_input))
    _validate_gc_content_config(cfg.objects.gc_content)
    _validate_depth_config(cfg.objects.depth)

    if not records:
        ordered_indices, row_counts = [], []
    elif multi_record_positions:
        resolved_order, rows_by_record = resolve_record_row_positions(
            records,
            multi_record_positions,
            _compatibility="circular",
        )
        ordered_indices = list(resolved_order)
        row_counts = [
            sum(rows_by_record[index] == row for index in ordered_indices)
            for row in sorted(set(rows_by_record))
        ]
    else:
        ordered_indices = list(range(len(records)))
        row_counts = _resolve_multi_record_default_row_counts(len(records))
    records = [records[idx] for idx in ordered_indices]
    if record_depth_tracks is not None:
        record_depth_tracks = [record_depth_tracks[idx] for idx in ordered_indices]
    if ordered_indices != list(range(len(ordered_indices))):
        display_index_by_input = {
            input_index: display_index
            for display_index, input_index in enumerate(ordered_indices)
        }
        resolved_annotations = replace(
            resolved_annotations,
            annotations=tuple(
                replace(
                    annotation,
                    record_index=display_index_by_input[annotation.record_index],
                )
                for annotation in resolved_annotations.annotations
            ),
        )

    parsed_circular_track_slots = _parse_circular_track_slot_inputs(circular_track_slots)
    circular_track_axis_index = _validate_circular_track_axis_index(
        circular_track_axis_index,
        parsed_circular_track_slots,
    )

    legend_effective = legend

    show_depth, show_gc, show_skew = _resolve_circular_track_visibility(
        cfg,
        parsed_circular_track_slots,
    )
    if (
        parsed_circular_track_slots is not None
        and show_depth
        and record_depth_tracks is None
    ):
        raise ValidationError("A circular depth track slot requires depth_tables, depth_files, or depth_track input.")
    available_depth_track_count = depth_track_count(record_depth_tracks)
    parsed_circular_track_slots = _default_circular_depth_slots_if_needed(
        parsed_circular_track_slots=parsed_circular_track_slots,
        show_depth=bool(show_depth and record_depth_tracks is not None),
        depth_track_count_value=available_depth_track_count,
        show_ticks=bool(cfg.objects.scale.show),
        show_gc=bool(show_gc),
        show_skew=bool(show_skew),
        dinucleotide=dinucleotide,
    )
    canvas_cfg = replace(
        cfg.canvas,
        show_depth=bool(show_depth and record_depth_tracks is not None),
        show_gc=bool(show_gc),
        show_skew=bool(show_skew),
    )
    cfg = replace(cfg, canvas=canvas_cfg)
    show_plot_title = normalized_plot_title_position != "none"
    record_definition_profile = _circular_record_definition_profile(
        record_count=len(records),
        show_plot_title=show_plot_title,
        keep_full_definition_with_plot_title=keep_full_definition_with_plot_title,
    )
    layout_definition_profile = _circular_record_definition_profile(
        record_count=len(records),
        show_plot_title=False,
        keep_full_definition_with_plot_title=keep_full_definition_with_plot_title,
    )
    record_lengths = [len(record.seq) for record in records]
    cfg = _harmonize_multi_record_circular_style_cfg(
        cfg,
        record_lengths=record_lengths,
    )
    profile = CircularRenderProfile(
        cfg, identity_label_overrides=_shows_identity_labels(_resolved_record_features)
    )
    conservation_mode = normalize_conservation_reference(conservation_reference)
    conservation_load_result = _load_conservation_result(
        conservation_blast_files=conservation_blast_files,
        conservation_dataframes=conservation_dataframes,
        conservation_labels=conservation_labels,
        conservation_colors=conservation_colors,
        thresholds=thresholds,
    )
    first_record_conservation_tracks = _conservation_tracks_for_record(
        conservation_load_result, records, records[0], conservation_mode
    )
    if first_record_conservation_tracks:
        if _circular_slots_define_renderer(parsed_circular_track_slots, "sequence_conservation"):
            parsed_circular_track_slots = _apply_conservation_track_params_to_slots(
                parsed_circular_track_slots or (),
                first_record_conservation_tracks,
                cfg=cfg,
                conservation_ring_width=conservation_ring_width,
                conservation_ring_gap=conservation_ring_gap,
            )
        else:
            conservation_slots = _conservation_slots_for_tracks(
                first_record_conservation_tracks,
                cfg=cfg,
                conservation_ring_width=conservation_ring_width,
                conservation_ring_gap=conservation_ring_gap,
            )
            if parsed_circular_track_slots is None:
                parsed_circular_track_slots = circular_track_slots_from_order(
                    "features,ticks,depth,gc_content,gc_skew",
                    show_ticks=bool(cfg.objects.scale.show),
                    show_depth=show_depth,
                    depth_track_count=max(1, depth_track_count(record_depth_tracks)),
                    show_gc=show_gc,
                    show_skew=show_skew,
                    dinucleotide=dinucleotide,
                )
            parsed_circular_track_slots = _insert_conservation_slots(
                parsed_circular_track_slots,
                conservation_slots,
                axis_index=circular_track_axis_index,
            )
    record_scales = _resolve_multi_record_scales(
        record_lengths,
        mode=normalized_multi_record_size_mode,
        min_radius_ratio=normalized_multi_record_min_radius_ratio,
    )
    # None means no depth input at all; [] means depth input without a cell
    # for that record, which still reserves the record's depth slot.
    record_depth_track_data: list[list[DepthTrackData]] | None = (
        None if record_depth_tracks is None else [[] for _ in records]
    )
    if cfg.canvas.show_depth and record_depth_tracks is not None:
        record_depth_window_steps: list[tuple[int, int]] = []
        for record in records:
            record_window, record_step = _resolve_circular_window_step(
                record,
                cfg,
                window=window,
                step=step,
            )
            record_depth_window, record_depth_step = _resolve_depth_window_step(
                window=record_window,
                step=record_step,
                depth_window=depth_window,
                depth_step=depth_step,
            )
            record_depth_window_steps.append((record_depth_window, record_depth_step))
        first_depth_window, first_depth_step = record_depth_window_steps[0]
        base_depth_config = DepthConfigurator(
            window=first_depth_window,
            step=first_depth_step,
            profile=profile,
        )
        record_depth_track_data = build_depth_track_dataframes(
            records,
            record_depth_tracks,
            base_config=base_depth_config,
            depth_df_builder=build_depth_df,
            window_steps=record_depth_window_steps,
        )
        _validate_depth_track_indices(
            parsed_circular_track_slots,
            depth_track_count_value=max(1, depth_track_data_count(record_depth_track_data)),
            available_indices=depth_track_indices(record_depth_track_data),
        )

    record_results: list[CircularAssemblyResult] = []
    record_radii_px: list[float] = []
    track_slot_geometry_records: list[dict[str, Any]] = []
    record_id_counts = Counter(str(record.id) for record in records)
    total_record_count = len(records)
    for record_index, (record, record_scale) in enumerate(zip(records, record_scales)):
        scaled_cfg = _scale_circular_cfg(cfg, scale=record_scale)
        record_radii_px.append(float(scaled_cfg.canvas.circular.radius))
        raw_record_id = str(record.id)
        definition_group_id = (
            definition_group_svg_id(
                raw_record_id,
                mode="circular",
                record_index=record_index,
                record_count=total_record_count,
            )
            if record_id_counts[raw_record_id] > 1
            else None
        )
        record_conservation_tracks = _conservation_tracks_for_record(
            conservation_load_result, records, record, conservation_mode
        )
        sub_canvas = assemble_circular_diagram_from_record(
            record,
            cfg=scaled_cfg,
            conservation_reference=conservation_mode,
            color_table=color_table,
            default_colors=default_colors,
            selected_features_set=list(selected_features_set),
            feature_visibility_table=feature_table,
            feature_shapes=feature_shapes,
            output_prefix=output_prefix,
            legend="none",
            dinucleotide=dinucleotide,
            window=window,
            step=step,
            depth_window=depth_window,
            depth_step=depth_step,
            species=species,
            strain=strain,
            plot_title=None,
            plot_title_position="none",
            center_reserved_radius=center_reserved_radius,
            circular_track_slots=parsed_circular_track_slots,
            circular_track_axis_index=circular_track_axis_index,
            _resolved_annotations=resolved_annotations,
            _annotation_record_index=record_index,
            _definition_group_id=definition_group_id,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
            _definition_profile=record_definition_profile,
            _layout_definition_profile=layout_definition_profile,
            _precomputed_depth_tracks=(
                record_depth_track_data[record_index]
                if record_depth_track_data is not None
                else None
            ),
            _precomputed_depth_track_count=(
                available_depth_track_count
                if record_depth_track_data is not None
                else None
            ),
            _precomputed_conservation_tracks=record_conservation_tracks,
            _resolved_feature_inputs=resolved_feature_inputs,
            _resolved_record_features=(_resolved_record_features[record_index],) if _resolved_record_features else (),
            _record_transform=(_record_transforms[record_index] if _record_transforms is not None else None),
            _record_presentation=(
                _record_presentations[record_index] if _record_presentations is not None else None
            ),
        )
        result = _require_circular_assembly_result(sub_canvas)
        record_results.append(result)
        sub_geometry = result.track_slot_geometry
        if isinstance(sub_geometry, Mapping):
            for record_payload in sub_geometry.get("records", []) or []:
                if not isinstance(record_payload, Mapping):
                    continue
                updated_record = dict(record_payload)
                updated_record["recordIndex"] = int(record_index)
                track_slot_geometry_records.append(updated_record)

    max_record_radius_px = max(
        [float(radius) for radius in record_radii_px if float(radius) > 0.0],
        default=float(cfg.canvas.circular.radius),
    )
    column_gap_px = _resolve_multi_record_grid_gap_px(
        max_record_radius_px,
        gap_ratio=normalized_multi_record_column_gap_ratio,
    )
    row_gap_px = _resolve_multi_record_grid_gap_px(
        max_record_radius_px,
        gap_ratio=normalized_multi_record_row_gap_ratio,
    )

    if not row_counts:
        row_counts = _resolve_multi_record_default_row_counts(len(record_results))
    grid_layout = _pack_circular_record_bounds(
        [result.content_bounds for result in record_results],
        row_counts,
        column_gap_px=column_gap_px,
        row_gap_px=row_gap_px,
    )

    legend_table = LegendTable()
    legend_target: Group | None = None
    legend_content: GroupContent | None = None
    legend_bounds: Aabb | None = None
    plot_title_content: GroupContent | None = None

    if show_plot_title:
        plot_title_content = _circular_plot_title_content(
            gb_record=records[0],
            output_prefix=output_prefix,
            cfg=cfg,
            profile=profile,
            species=species,
            strain=strain,
            plot_title=normalized_plot_title or None,
            record_presentation=_record_presentations[0] if _record_presentations else None,
        )

    legend_canvas_config = CircularCanvasConfigurator(
        output_prefix=output_prefix,
        profile=profile,
        legend=legend_effective,
        gb_record=records[0],
    )
    legend_window, legend_step = _resolve_circular_window_step(
        records[0],
        cfg,
        window=window,
        step=step,
    )
    gc_config = GcContentConfigurator(
        window=legend_window,
        step=legend_step,
        dinucleotide=dinucleotide,
        profile=profile,
        default_colors_df=default_colors,
    )
    skew_config = GcSkewConfigurator(
        window=legend_window,
        step=legend_step,
        dinucleotide=dinucleotide,
        profile=profile,
        default_colors_df=default_colors,
    )
    depth_config = (
        DepthConfigurator(
            window=legend_window,
            step=legend_step,
            profile=profile,
        )
        if cfg.canvas.show_depth
        else None
    )
    feature_config = FeatureDrawingConfigurator(
        color_table=color_table,
        default_colors=default_colors,
        selected_features_set=list(selected_features_set),
        profile=profile,
        feature_table=feature_table,
        feature_shapes=feature_shapes,
        feature_visibility_rules=resolved_feature_inputs.feature_visibility_rules,
        specific_color_rules=resolved_feature_inputs.specific_color_rules,
        default_color_map=resolved_feature_inputs.default_color_map,
        canvas_config=legend_canvas_config,
        record_features=_resolved_record_features,
    )
    legend_sources = LegendSources.of(feature_config, gc_config, skew_config)
    if legend_effective != "none":
        legend_slots, legend_annotations = plan_circular_annotation_slots(
            resolved_annotations,
            records,
            parsed_circular_track_slots,
            show_ticks=bool(cfg.objects.scale.show),
            show_depth=bool(profile.show_depth and depth_config is not None),
            show_gc=profile.show_gc,
            show_skew=profile.show_skew,
            depth_track_count=max(1, available_depth_track_count),
        )
        legend_depth_tracks = representative_depth_tracks(record_depth_track_data)
        legend_depth_track_zero = index_depth_track_row(legend_depth_tracks).get(0)
        legend_table = circular_legend_table(
            records,
            legend_sources,
            profile=profile,
            depth_config=depth_config,
            depth_df=(
                legend_depth_track_zero.df
                if legend_depth_track_zero is not None
                else None
            ),
            depth_tracks=legend_depth_tracks,
            depth_track_count_value=depth_track_data_count(record_depth_track_data),
            circular_track_slots=legend_slots,
            annotations=legend_annotations,
            conservation_legend=conservation_legend_sources(first_record_conservation_tracks),
            conservation_min_identity=float(identity),
        )
        if legend_table:
            legend_config = LegendDrawingConfigurator(
                profile=profile,
                canvas_config=legend_canvas_config,
            )
            legend_measurement = measure_legend(
                legend_config,
                legend_table,
                placement=legend_canvas_config.legend_position,
                wrap_width=grid_layout.bounds.width,
            )
            legend_bounds = legend_measurement.local_bounds
            legend_content = render_stage(
                "legend_content", circular_legend_content, legend_measurement, legend_table
            )

    legend_row_facts = LegendRowFacts(records=records, sources=legend_sources, legend_table=legend_table)

    merged_canvas = Drawing(
        filename=f"{output_prefix}.svg",
        size=(f"{grid_layout.bounds.width}px", f"{grid_layout.bounds.height}px"),
        viewBox=(
            f"0 0 {grid_layout.bounds.width} {grid_layout.bounds.height}"
        ),
        debug=False,
    )
    attach_legend_row_facts(merged_canvas, legend_row_facts)
    record_targets: list[Group] = []
    grid_overlay_obstacles: list[Aabb] = []
    bind_record_identity = len(records) > 1 and (
        len({record.id for record in records}) != len(records)
        or bool(
            _record_transforms
            and any(transform.start_coordinate is not None for transform in _record_transforms)
        )
    )
    merged_record: MergedCircularRecord | None = None

    for record_index, result in enumerate(record_results):
        placed, definition_positions = placed_record_contents(result.drawing)
        merged_record = render_stage(
            "record_merge",
            merge_circular_record,
            placed,
            definition_positions,
            record_index=record_index,
            record_id=str(records[record_index].id),
            bind_record_identity=bind_record_identity,
            previous=merged_record,
        )
        record_group = place_content(
            merged_canvas,
            merged_record.content,
            GroupKey(record_key_at(records[record_index], record_index), GroupKind.RECORD),
        )
        record_dx, record_dy = grid_layout.translations[record_index]
        record_group.translate(record_dx, record_dy)
        merged_canvas.add(record_group)
        record_targets.append(record_group)
        grid_overlay_obstacles.extend(
            obstacle.translated(record_dx, record_dy)
            for obstacle in result.overlay_obstacles
        )

    plot_title_group: Group | None = None
    plot_title_bounds: Aabb | None = None
    if plot_title_content is not None:
        plot_title_group = place_content(
            merged_canvas, plot_title_content, GroupKey(CANVAS_SCOPE, GroupKind.PLOT_TITLE)
        )
        plot_title_bounds = plot_title_content.local_bounds
        merged_canvas.add(plot_title_group)
    if legend_content is not None:
        legend_target = place_content(
            merged_canvas, legend_content, GroupKey(CANVAS_SCOPE, GroupKind.LEGEND)
        )
        merged_canvas.add(legend_target)

    legend_placement = LegendPlacement(str(legend_effective))
    title_placement = (
        TitlePlacement(normalized_plot_title_position)
        if plot_title_group is not None
        else TitlePlacement.NONE
    )
    composition_plan = render_stage(
        "composition_plan",
        plan_composition,
        CompositionRequest(
            primary=CompositionItem("primary", grid_layout.bounds),
            legend=(
                CompositionItem("legend", legend_bounds)
                if legend_target is not None and legend_bounds is not None
                else None
            ),
            title=(
                CompositionItem("title", plot_title_bounds)
                if plot_title_group is not None and plot_title_bounds is not None
                else None
            ),
            legend_placement=legend_placement,
            title_placement=title_placement,
            overlay_obstacles=tuple(grid_overlay_obstacles),
        ),
    )
    apply_composition_plan(
        merged_canvas,
        composition_plan,
        primary_targets=record_targets,
        legend_target=legend_target,
        legend_side=legend_placement,
        legend_reflow_metrics=(
            legend_measurement.reflow_metrics
            if legend_target is not None
            else None
        ),
        title_target=plot_title_group,
        title_side=title_placement,
    )

    _attach_circular_label_layouts(
        merged_canvas, [result.label_layout for result in record_results]
    )
    if track_slot_geometry_records:
        setattr(
            merged_canvas,
            "_gbdraw_track_slot_geometry",
            {
                "schema": 1,
                "mode": "circular",
                "source": "resolved",
                "records": track_slot_geometry_records,
            },
        )

    return merged_canvas


_CIRCULAR_LABEL_LAYOUTS_ATTR = "_gbdraw_circular_label_layouts"
_LABEL_FIT_ATTR = "_gbdraw_label_fit"


def _attach_circular_label_layouts(
    drawing: Drawing,
    layouts: Sequence[CircularLabelLayout | None],
) -> None:
    setattr(
        drawing,
        _CIRCULAR_LABEL_LAYOUTS_ATTR,
        tuple(layout for layout in layouts if layout is not None),
    )


def _circular_label_layouts(drawing: Drawing) -> tuple[CircularLabelLayout, ...]:
    layouts = getattr(drawing, _CIRCULAR_LABEL_LAYOUTS_ATTR, ())
    return tuple(layout for layout in layouts if isinstance(layout, CircularLabelLayout))


def _drawing_size(drawing: Drawing) -> tuple[float, float]:
    return (
        _parse_svg_number(drawing.attribs.get("width")),
        _parse_svg_number(drawing.attribs.get("height")),
    )


# The reason of a row whose settings leave the tracks no room (OV-1266, D-07).
_TRACKS_DO_NOT_FIT_REASON = "The tracks do not fit with this setting"


def _unfit_row_verification(error: ValidationError) -> LabelFitVerification | None:
    """A row's assembly error that means the row does not fit, or None when the error must be raised.

    A label setting the engine cannot draw (``LabelLayoutError``) does not fit, and so do tracks that no longer fit
    (``TRACK_LAYOUT``): the figure assembled without the row's settings, so the row's settings caused it. Any other
    error, such as a missing file, is raised (OV-733).
    """

    if isinstance(error, LabelLayoutError):
        return LabelFitVerification(outcome="stop", canvas_size=None, reason=str(error))
    if track_layout_reason(error) is not None:
        return LabelFitVerification(outcome="stop", canvas_size=None, reason=_TRACKS_DO_NOT_FIT_REASON)
    return None


def _checked_circular_drawing(
    assemble: Callable[[Mapping[str, LabelFitValue]], Drawing],
    *,
    warn_search_seconds: float,
) -> Drawing:
    """Assemble a Circular figure and apply the label fit check to it.

    Overlapping labels (or radial labels without room) raise
    ``LabelsDoNotFitError``; crossing leaders or labels over the definition
    text alone keep the drawing, warn, and attach the report and its search
    (rows searched for ``warn_search_seconds``).
    """

    started = time.monotonic()
    drawing = assemble({})
    assembly_seconds = time.monotonic() - started

    def verify(overrides: Mapping[str, LabelFitValue]) -> LabelFitVerification:
        try:
            verified = assemble(overrides)
        except ValidationError as error:
            unfit = _unfit_row_verification(error)
            if unfit is None:
                raise
            return unfit
        layouts = _circular_label_layouts(verified)
        return LabelFitVerification(
            outcome=circular_label_outcome(layouts),
            canvas_size=_drawing_size(verified),
            fit=sum((layout.fit for layout in layouts), FITS),
        )

    checked = check_circular_label_fit(
        _circular_label_layouts(drawing),
        verify,
        assembly_seconds=assembly_seconds,
        warn_search_seconds=warn_search_seconds,
    )
    if checked is not None:
        setattr(drawing, _LABEL_FIT_ATTR, checked)
    return drawing


def _label_fit(drawing: Drawing) -> tuple[LabelFitReport, CircularLabelFitSearch | None] | None:
    """A drawn figure's report and the search of its rows (a Linear figure's report has none)."""

    report, search = getattr(drawing, _LABEL_FIT_ATTR, (None, None))
    if isinstance(report, LabelFitReport) and (search is None or isinstance(search, CircularLabelFitSearch)):
        return report, search
    return None


def label_fit_report(drawing: Drawing) -> LabelFitReport | None:
    """The report of a figure drawn with labels that do not fit (crossing leaders, labels over
    the definition text, Linear region labels over other labels or lines), if any."""

    checked = _label_fit(drawing)
    return checked[0] if checked is not None else None


def label_fit_row(
    source: Drawing | LabelsDoNotFitError, kind: LabelFitRowKind, *, deadline: float | None = None
) -> LabelFitRow | None:
    """Compute one row of a report, until a ``time.monotonic()`` ``deadline`` if given.

    ``source`` is a figure drawn with a report or the error of a figure that
    was not drawn. Rows that the bounded search did not reach are
    ``not_computed`` in the report; this returns the computed row (rows found
    earlier cost nothing), or ``not_computed`` again when ``deadline`` came
    first.
    """

    if isinstance(source, LabelsDoNotFitError):
        return source.row_search.row(kind, deadline=deadline) if source.row_search is not None else None
    checked = _label_fit(source)
    return checked[1].row(kind, deadline=deadline) if checked is not None and checked[1] is not None else None


def _circular_builder_options(
    options: CircularDiagramOptions | None,
) -> CircularDiagramOptions:
    if options is None:
        return resolve_circular_diagram_options(CircularDiagramOptions())
    if isinstance(options, CircularDiagramOptions):
        return resolve_circular_diagram_options(options)
    raise ValidationError("options must be CircularDiagramOptions.")


def _linear_builder_options(
    options: LinearDiagramOptions | None,
) -> LinearDiagramOptions:
    if options is None:
        return resolve_linear_diagram_options(LinearDiagramOptions())
    if isinstance(options, LinearDiagramOptions):
        return resolve_linear_diagram_options(options)
    raise ValidationError("options must be LinearDiagramOptions.")


def build_circular_diagram(
    gb_record: SeqRecord,
    *,
    options: CircularDiagramOptions | None = None,
    _precomputed_depth_track_specs: Sequence[DepthTrackSpec] | None = None,
    _precomputed_depth_track_count: int | None = None,
    _resolved_feature_inputs: ResolvedFeatureInputs | None = None,
    _resolved_record_features: tuple[ResolvedRecordFeatureInputs, ...] = (),
    _record_transform: RecordDisplayTransform | None = None,
    _record_presentation: RecordPresentation | None = None,
    _resolved_annotations: ResolvedAnnotationBundle | None = None,
    _label_fit_warn_search_seconds: float = WARN_SEARCH_SECONDS,
) -> Drawing:
    """Build a circular diagram using mode-specific typed options."""

    options = _circular_builder_options(options)
    depth_table, depth_file = _resolve_single_circular_depth_options(options)
    output = options.output
    tracks = options.tracks

    cfg = _resolve_diagram_options_config(options)
    # Read every file once: the label fit rows reassemble the figure, on the
    # Web after its staged files are deleted (OV-720, OV-732).
    conservation_inputs = read_conservation_inputs(
        options.conservation_search_results or options.conservation_blast_files
    )
    feature_inputs = _feature_inputs_read_once(
        options, load_comparison_colors=False, resolved_inputs=_resolved_feature_inputs
    )
    annotations = read_annotation_options(options.annotations, mode="circular")
    depth_tracks = (
        None
        if _precomputed_depth_track_specs is not None
        else depth_track_inputs(
            normalize_depth_tracks(
                [gb_record],
                depth_tracks=options.depth_tracks,
                depth_table=depth_table,
                depth_file=depth_file,
                depth_track_tables=options.depth_track_tables,
                depth_track_files=options.depth_track_files,
                depth_track_labels=options.depth_track_labels,
                depth_track_colors=options.depth_track_colors,
                depth_track_large_tick_intervals=options.depth_track_large_tick_intervals,
                depth_track_small_tick_intervals=options.depth_track_small_tick_intervals,
                depth_track_tick_font_sizes=options.depth_track_tick_font_sizes,
            )
        )
    )

    def assemble(overrides: Mapping[str, LabelFitValue]) -> Drawing:
        return assemble_circular_diagram_from_record(
            gb_record,
            cfg=with_fit_row_overrides(cfg, overrides),
            selected_features_set=options.selected_features_set,
            feature_shapes=options.feature_shapes,
            output_prefix="out",
            legend=output.legend if output else "right",
            dinucleotide=options.dinucleotide,
            window=options.window,
            step=options.step,
            depth_window=options.depth_window,
            depth_step=options.depth_step,
            depth_tracks=depth_tracks,
            conservation_blast_files=conservation_inputs,
            conservation_dataframes=options.conservation_dataframes,
            conservation_reference=options.conservation_reference,
            conservation_labels=options.conservation_labels,
            conservation_colors=options.conservation_colors,
            conservation_ring_width=options.conservation_ring_width,
            conservation_ring_gap=options.conservation_ring_gap,
            # The builder options resolved omitted thresholds from the mode profile.
            evalue=cast(float, options.evalue),
            bitscore=cast(float, options.bitscore),
            identity=cast(float, options.identity),
            alignment_length=cast(int, options.alignment_length),
            species=options.species,
            strain=options.strain,
            plot_title=options.plot_title,
            plot_title_position=(
                output.plot_title_position
                if output and output.plot_title_position is not None
                else "none"
            ),
            plot_title_font_size=options.plot_title_font_size,
            keep_full_definition_with_plot_title=options.keep_full_definition_with_plot_title,
            center_reserved_radius=tracks.center_reserved_radius if tracks else None,
            circular_track_slots=tracks.circular_track_slots if tracks else None,
            circular_track_axis_index=tracks.circular_track_axis_index if tracks else None,
            annotation_options=annotations,
            _resolved_annotations=_resolved_annotations,
            _precomputed_depth_track_specs=_precomputed_depth_track_specs,
            _precomputed_depth_track_count=_precomputed_depth_track_count,
            _resolved_feature_inputs=feature_inputs,
            _resolved_record_features=_resolved_record_features,
            _record_transform=_record_transform,
            _record_presentation=_record_presentation,
        )

    return _checked_circular_drawing(
        assemble, warn_search_seconds=_label_fit_warn_search_seconds
    )


def _build_linear_diagram(
    records: Sequence[SeqRecord],
    *,
    options: LinearDiagramOptions | None = None,
    layout: LinearMultiRecordOptions | None = None,
    similarity_alignment: SimilarityAlignmentPlan | None = None,
    losatp_cache: LosatpCacheManager | None = None,
    protein_extraction: ProteinExtractionResult | None = None,
    _resolved_feature_inputs: ResolvedFeatureInputs | None = None,
    _resolved_record_features: tuple[ResolvedRecordFeatureInputs, ...] = (),
    _return_build_result: bool = False,
    _record_transforms: Sequence[RecordDisplayTransform] | None = None,
    _record_presentations: Sequence[RecordPresentation] | None = None,
    _alignment_anchor_centers: Sequence[float | None] | None = None,
    _resolved_annotations: ResolvedAnnotationBundle | None = None,
) -> Drawing | LinearDiagramBuildResult:
    """Build a linear diagram using mode-specific typed options."""

    options = _linear_builder_options(options)
    output = options.output
    tracks = options.tracks
    cfg = _resolve_diagram_options_config(options)
    normalized_collinearity_anchor_mode = normalize_collinearity_anchor_mode(
        str(options.collinearity_anchor_mode)
    )
    # Read every file once: the label fit row reassembles the figure (OV-732).
    feature_inputs = _feature_inputs_read_once(
        options,
        load_comparison_colors=_linear_has_comparisons(
            blast_files=options.blast_files,
            linear_comparisons=options.linear_comparisons,
            protein_comparisons=options.protein_comparisons,
            collinearity_blocks=options.collinearity_blocks,
            losat_search=options.losat_search,
        ),
        resolved_inputs=_resolved_feature_inputs,
    )
    annotations = read_annotation_options(options.annotations, mode="linear")
    depth_tracks = depth_track_inputs(
        normalize_depth_tracks(
            records,
            depth_tracks=options.depth_tracks,
            depth_table=options.depth_table,
            depth_file=options.depth_file,
            depth_tables=options.depth_tables,
            depth_files=options.depth_files,
            depth_track_tables=options.depth_track_tables,
            depth_track_files=options.depth_track_files,
            depth_track_labels=options.depth_track_labels,
            depth_track_colors=options.depth_track_colors,
            depth_track_heights=options.depth_track_heights,
            depth_track_large_tick_intervals=options.depth_track_large_tick_intervals,
            depth_track_small_tick_intervals=options.depth_track_small_tick_intervals,
            depth_track_tick_font_sizes=options.depth_track_tick_font_sizes,
        )
    )

    def assemble(
        overrides: Mapping[str, LabelFitValue], feature_label_rows: FeatureLabelRows | None
    ) -> Drawing | LinearDiagramBuildResult:
        return assemble_linear_diagram_from_records(
            records,
            cfg=with_config_overrides(cfg, overrides),
            blast_files=options.blast_files,
            linear_comparisons=options.linear_comparisons,
            layout=layout,
            protein_comparisons=options.protein_comparisons,
            orthogroups=options.orthogroups,
            losat_search=options.losat_search,
            pairwise_match_style=options.pairwise_match_style,
            collinearity_blocks=options.collinearity_blocks,
            collinearity_params=options.collinearity_params,
            collinearity_unit_mode=options.collinearity_unit_mode,
            collinearity_anchor_mode=normalized_collinearity_anchor_mode,
            collinearity_search_scope=options.collinearity_search_scope,
            collinearity_color_mode=options.collinearity_color_mode,
            losatp_cache=losatp_cache,
            protein_extraction=protein_extraction,
            orthogroup_membership_mode=options.orthogroup_membership_mode,
            collinear_infer_orthogroups=options.collinear_infer_orthogroups,
            collinear_max_paralog_links_per_orthogroup=options.collinear_max_paralog_links_per_orthogroup,
            similarity_alignment=similarity_alignment,
            selected_features_set=options.selected_features_set,
            feature_shapes=options.feature_shapes,
            output_prefix="out",
            legend=output.legend if output else "right",
            dinucleotide=options.dinucleotide,
            window=options.window,
            step=options.step,
            depth_window=options.depth_window,
            depth_step=options.depth_step,
            depth_tracks=depth_tracks,
            linear_track_slots=tracks.linear_track_slots if tracks else None,
            linear_track_axis_index=tracks.linear_track_axis_index if tracks else None,
            annotation_options=annotations,
            _resolved_annotations=_resolved_annotations,
            plot_title=options.plot_title,
            plot_title_position=(
                output.plot_title_position
                if output and output.plot_title_position is not None
                else "bottom"
            ),
            plot_title_font_size=options.plot_title_font_size,
            # The builder options resolved omitted thresholds from the mode profile.
            evalue=cast(float, options.evalue),
            bitscore=cast(float, options.bitscore),
            identity=cast(float, options.identity),
            alignment_length=cast(int, options.alignment_length),
            _resolved_feature_inputs=feature_inputs,
            _resolved_record_features=_resolved_record_features,
            _record_transforms=_record_transforms,
            _record_presentations=_record_presentations,
            _alignment_anchor_centers=_alignment_anchor_centers,
            _return_build_result=_return_build_result,
            _feature_label_rows=feature_label_rows,
        )

    return _checked_linear_build(assemble, features_separated=bool(cfg.canvas.resolve_overlaps))


def _linear_drawing(result: Drawing | LinearDiagramBuildResult) -> Drawing:
    return result.drawing if isinstance(result, LinearDiagramBuildResult) else result


def _linear_label_fit(drawing: Drawing) -> LinearLabelFit:
    label_fit = getattr(drawing, LINEAR_LABEL_FIT_ATTR, None)
    return label_fit if isinstance(label_fit, LinearLabelFit) else LinearLabelFit()


def _checked_linear_build(
    assemble: Callable[[Mapping[str, LabelFitValue], FeatureLabelRows | None], Drawing | LinearDiagramBuildResult],
    *,
    features_separated: bool,
) -> Drawing | LinearDiagramBuildResult:
    """Assemble a Linear figure and apply the label fit check to it.

    A figure whose region labels need feature labels in a new row (LP-D82
    ``RGN-YIELD-ROW``) is assembled once more with those rows, and that
    figure is kept when ``keeps_feature_label_rows`` accepts it. A figure
    drawn with labels that do not fit keeps its report (the Web notice).
    """

    def build(overrides: Mapping[str, LabelFitValue]) -> Drawing | LinearDiagramBuildResult:
        first = assemble(overrides, None)
        label_fit = _linear_label_fit(_linear_drawing(first))
        if label_fit.new_rows is None:
            return first
        second = assemble(overrides, label_fit.new_rows)
        return second if keeps_feature_label_rows(label_fit, _linear_label_fit(_linear_drawing(second))) else first

    result = build({})

    def verify(overrides: Mapping[str, LabelFitValue]) -> LabelFitVerification:
        try:
            verified = _linear_drawing(build(overrides))
        except ValidationError as error:
            unfit = _unfit_row_verification(error)
            if unfit is None:
                raise
            return unfit
        label_fit = _linear_label_fit(verified)
        return LabelFitVerification(label_fit.outcome, _drawing_size(verified), label_fit.fit)

    drawing = _linear_drawing(result)
    report = check_linear_label_fit(_linear_label_fit(drawing), verify, features_separated=features_separated)
    if report is not None:
        setattr(drawing, _LABEL_FIT_ATTR, (report, None))
    return result


def build_linear_diagram(
    records: Sequence[SeqRecord],
    *,
    options: LinearDiagramOptions | None = None,
    layout: LinearMultiRecordOptions | None = None,
    similarity_alignment: SimilarityAlignmentPlan | None = None,
    losatp_cache: LosatpCacheManager | None = None,
    protein_extraction: ProteinExtractionResult | None = None,
    _resolved_feature_inputs: ResolvedFeatureInputs | None = None,
    _resolved_record_features: tuple[ResolvedRecordFeatureInputs, ...] = (),
    _record_transforms: Sequence[RecordDisplayTransform] | None = None,
    _record_presentations: Sequence[RecordPresentation] | None = None,
    _alignment_anchor_centers: Sequence[float | None] | None = None,
    _resolved_annotations: ResolvedAnnotationBundle | None = None,
) -> Drawing:
    """Build a linear diagram using mode-specific typed options."""

    result = _build_linear_diagram(
        records,
        options=options,
        layout=layout,
        similarity_alignment=similarity_alignment,
        losatp_cache=losatp_cache,
        protein_extraction=protein_extraction,
        _resolved_feature_inputs=_resolved_feature_inputs,
        _resolved_annotations=_resolved_annotations,
        _resolved_record_features=_resolved_record_features,
        _record_transforms=_record_transforms,
        _record_presentations=_record_presentations,
        _alignment_anchor_centers=_alignment_anchor_centers,
    )
    return cast(Drawing, result)


def build_linear_diagram_result(
    records: Sequence[SeqRecord],
    *,
    options: LinearDiagramOptions | None = None,
    layout: LinearMultiRecordOptions | None = None,
    similarity_alignment: SimilarityAlignmentPlan | None = None,
    losatp_cache: LosatpCacheManager | None = None,
    protein_extraction: ProteinExtractionResult | None = None,
    _resolved_feature_inputs: ResolvedFeatureInputs | None = None,
    _resolved_record_features: tuple[ResolvedRecordFeatureInputs, ...] = (),
    _record_transforms: Sequence[RecordDisplayTransform] | None = None,
    _record_presentations: Sequence[RecordPresentation] | None = None,
    _alignment_anchor_centers: Sequence[float | None] | None = None,
    _resolved_annotations: ResolvedAnnotationBundle | None = None,
) -> LinearDiagramBuildResult:
    """Build a Linear drawing with its computed analysis metadata."""

    result = _build_linear_diagram(
        records,
        options=options,
        layout=layout,
        similarity_alignment=similarity_alignment,
        losatp_cache=losatp_cache,
        protein_extraction=protein_extraction,
        _resolved_feature_inputs=_resolved_feature_inputs,
        _resolved_annotations=_resolved_annotations,
        _resolved_record_features=_resolved_record_features,
        _record_transforms=_record_transforms,
        _record_presentations=_record_presentations,
        _alignment_anchor_centers=_alignment_anchor_centers,
        _return_build_result=True,
    )
    if not isinstance(result, LinearDiagramBuildResult):  # pragma: no cover
        raise ValidationError("The Linear diagram builder did not return metadata.")
    return result


def build_circular_multi_diagram(
    records: Sequence[SeqRecord],
    *,
    options: CircularDiagramOptions | None = None,
    layout: CircularMultiRecordOptions | None = None,
    _resolved_feature_inputs: ResolvedFeatureInputs | None = None,
    _resolved_record_features: tuple[ResolvedRecordFeatureInputs, ...] = (),
    _record_transforms: Sequence[RecordDisplayTransform] | None = None,
    _record_presentations: Sequence[RecordPresentation] | None = None,
    _resolved_annotations: ResolvedAnnotationBundle | None = None,
    _label_fit_warn_search_seconds: float = WARN_SEARCH_SECONDS,
) -> Drawing:
    """Build a circular grid using mode-specific typed options."""

    options = _circular_builder_options(options)
    if layout is not None and not isinstance(layout, CircularMultiRecordOptions):
        raise ValidationError(
            "layout must be CircularMultiRecordOptions or None."
        )
    layout = layout or CircularMultiRecordOptions()
    output = options.output
    tracks = options.tracks
    cfg = _resolve_diagram_options_config(options)
    # Read every file once: the label fit rows reassemble the figure, on the
    # Web after its staged files are deleted (OV-720, OV-732).
    conservation_inputs = read_conservation_inputs(
        options.conservation_search_results or options.conservation_blast_files
    )
    feature_inputs = _feature_inputs_read_once(
        options, load_comparison_colors=False, resolved_inputs=_resolved_feature_inputs
    )
    annotations = read_annotation_options(options.annotations, mode="circular")
    depth_tracks = depth_track_inputs(
        normalize_depth_tracks(
            records,
            depth_tracks=options.depth_tracks,
            depth_table=options.depth_table,
            depth_file=options.depth_file,
            depth_tables=options.depth_tables,
            depth_files=options.depth_files,
            depth_track_tables=options.depth_track_tables,
            depth_track_files=options.depth_track_files,
            depth_track_labels=options.depth_track_labels,
            depth_track_colors=options.depth_track_colors,
            depth_track_large_tick_intervals=options.depth_track_large_tick_intervals,
            depth_track_small_tick_intervals=options.depth_track_small_tick_intervals,
            depth_track_tick_font_sizes=options.depth_track_tick_font_sizes,
        )
    )

    def assemble(overrides: Mapping[str, LabelFitValue]) -> Drawing:
        return assemble_circular_diagram_from_records(
            records,
            cfg=with_fit_row_overrides(cfg, overrides),
            conservation_blast_files=conservation_inputs,
            conservation_dataframes=options.conservation_dataframes,
            conservation_reference=options.conservation_reference,
            conservation_labels=options.conservation_labels,
            conservation_colors=options.conservation_colors,
            conservation_ring_width=options.conservation_ring_width,
            conservation_ring_gap=options.conservation_ring_gap,
            selected_features_set=options.selected_features_set,
            feature_shapes=options.feature_shapes,
            output_prefix="out",
            legend=output.legend if output else "right",
            dinucleotide=options.dinucleotide,
            window=options.window,
            step=options.step,
            depth_window=options.depth_window,
            depth_step=options.depth_step,
            depth_tracks=depth_tracks,
            species=options.species,
            strain=options.strain,
            plot_title=options.plot_title,
            plot_title_position=(
                output.plot_title_position
                if output and output.plot_title_position is not None
                else "none"
            ),
            plot_title_font_size=options.plot_title_font_size,
            keep_full_definition_with_plot_title=options.keep_full_definition_with_plot_title,
            center_reserved_radius=tracks.center_reserved_radius if tracks else None,
            multi_record_size_mode=layout.multi_record_size_mode,
            multi_record_min_radius_ratio=layout.multi_record_min_radius_ratio,
            multi_record_column_gap_ratio=layout.multi_record_column_gap_ratio,
            multi_record_row_gap_ratio=layout.multi_record_row_gap_ratio,
            multi_record_positions=layout.multi_record_positions,
            circular_track_slots=tracks.circular_track_slots if tracks else None,
            circular_track_axis_index=tracks.circular_track_axis_index if tracks else None,
            annotation_options=annotations,
            _resolved_annotations=_resolved_annotations,
            # The builder options resolved omitted thresholds from the mode profile.
            evalue=cast(float, options.evalue),
            bitscore=cast(float, options.bitscore),
            identity=cast(float, options.identity),
            alignment_length=cast(int, options.alignment_length),
            _resolved_feature_inputs=feature_inputs,
            _resolved_record_features=_resolved_record_features,
            _record_transforms=_record_transforms,
            _record_presentations=_record_presentations,
        )

    return _checked_circular_drawing(
        assemble, warn_search_seconds=_label_fit_warn_search_seconds
    )


__all__ = [
    "DEFAULT_SELECTED_FEATURES",
    "LinearDiagramMetadata",
    "build_circular_diagram",
    "build_circular_multi_diagram",
    "build_linear_diagram",
]
