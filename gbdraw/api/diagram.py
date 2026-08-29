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
from collections import Counter
from dataclasses import dataclass, replace
from typing import Any, Optional, Sequence, Mapping, Literal, cast

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from gbdraw.analysis.depth import depth_df as build_depth_df  # type: ignore[reportMissingImports]
from gbdraw.analysis.gc import circular_dinucleotide_content_df  # type: ignore[reportMissingImports]
from gbdraw.analysis.depth_tracks import (  # type: ignore[reportMissingImports]
    DepthTrackData,
    DepthTrackSpec,
    build_depth_track_dataframes,
    depth_track_count,
    depth_track_data_count,
    depth_track_heights as _depth_track_heights_from_specs,
    depth_track_indices,
    index_depth_track_row,
    normalize_depth_tracks,
    representative_depth_tracks,
    sync_depth_track_legend_entries,
)
from gbdraw.analysis.conservation import (  # type: ignore[reportMissingImports]
    ConservationLoadResult,
    ConservationTrack,
    conservation_track_gradient_colors,
    load_conservation_sources,
    normalize_conservation_reference,
    normalize_conservation_tracks_for_record,
)
from gbdraw.analysis.protein_colinearity import (  # type: ignore[reportMissingImports]
    LosatpCacheManager,
    OrthogroupMembershipMode,
    OrthogroupResult,
    ProteinBlastpResult,
    ProteinExtractionResult,
    ProteinBlastpMode,
    build_pairwise_protein_blastp_comparisons,
    build_rbh_orthogroup_protein_blastp_comparisons,
    normalize_orthogroup_membership_mode,
    normalize_protein_blastp_mode,
)
from gbdraw.analysis.collinearity import (  # type: ignore[reportMissingImports]
    CollinearityBlock,
    CollinearityAnchorMode,
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
from gbdraw.config.models.objects import (  # type: ignore[reportMissingImports]
    normalize_pairwise_match_style,
)
from gbdraw.analysis.collinearity_units import CollinearityUnitMode  # type: ignore[reportMissingImports]
from gbdraw.analysis.skew import skew_df  # type: ignore[reportMissingImports]
from gbdraw.api.config import apply_config_overrides  # type: ignore[reportMissingImports]
from gbdraw.api.options import (  # type: ignore[reportMissingImports]
    AnnotationOptions,
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    DepthTrackInput,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    resolve_circular_diagram_options,
    resolve_linear_diagram_options,
)
from gbdraw.linear_comparison import LinearComparison
from gbdraw.layout.linear_multi_record import record_pairs_between_adjacent_rows
from gbdraw.layout.record_placement import resolve_record_row_positions
from gbdraw.canvas import CircularCanvasConfigurator, LinearCanvasConfigurator  # type: ignore[reportMissingImports]
from gbdraw.config.models import (  # type: ignore[reportMissingImports]
    CircularRenderProfile,
    GbdrawConfig,
    LinearRenderProfile,
)
from gbdraw.config.models.labels import LabelsFilteringConfig  # type: ignore[reportMissingImports]
from gbdraw.io.colors import load_default_colors, read_color_table  # type: ignore[reportMissingImports]
from gbdraw.labels.filtering import (  # type: ignore[reportMissingImports]
    read_filter_list_file,
    read_label_override_file,
    read_qualifier_priority_file,
)
from gbdraw.features.visibility import read_feature_visibility_file  # type: ignore[reportMissingImports]
from gbdraw.configurators import (  # type: ignore[reportMissingImports]
    BlastMatchConfigurator,
    DepthConfigurator,
    FeatureDrawingConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
    LegendDrawingConfigurator,
)
from gbdraw.core.sequence import create_dict_for_sequence_lengths, check_feature_presence  # type: ignore[reportMissingImports]
from gbdraw.diagrams.circular.assemble import (  # type: ignore[reportMissingImports]
    CircularAssemblyResult,
    _assemble_circular_diagram_result,
)
from gbdraw.diagrams.linear import assemble_linear_diagram  # type: ignore[reportMissingImports]
from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]
from gbdraw.layout.composition import (
    CompositionItem,
    CompositionRequest,
    LegendPlacement,
    TitlePlacement,
    plan_composition,
)
from gbdraw.layout.spatial import Aabb
from gbdraw.mode_profiles import (
    CIRCULAR_MODE_PROFILE,
    ComparisonThresholds,
    DEFAULT_FEATURE_TYPES,
    LINEAR_MODE_PROFILE,
)
from gbdraw.annotations import ResolvedAnnotationBundle, resolve_annotations
from gbdraw.features.colors import precompute_used_color_rules  # type: ignore[reportMissingImports]
from gbdraw.legend.table import (  # type: ignore[reportMissingImports]
    configure_pairwise_identity_legend_from_comparisons,
    prepare_legend_table,
)
from gbdraw.render.groups.circular import DefinitionGroup, LegendGroup  # type: ignore[reportMissingImports]
from gbdraw.render.composition import (
    COMPOSITION_ROLE_ATTRIBUTE,
    apply_composition_plan,
)
from gbdraw.svg.ids import definition_group_svg_id
from gbdraw.tracks import (  # type: ignore[reportMissingImports]
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

from .prepared import ResolvedFeatureInputs, resolve_feature_inputs

DEFAULT_SELECTED_FEATURES = DEFAULT_FEATURE_TYPES

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class LinearDiagramMetadata:
    """Computed Linear analysis values consumed after diagram assembly."""

    protein_comparisons: tuple[DataFrame, ...] | None = None
    linear_comparisons: tuple[LinearComparison, ...] = ()
    orthogroups: OrthogroupResult | None = None
    collinearity_result: CollinearityResult | None = None


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


def _resolve_diagram_options_config(
    options: CircularDiagramOptions | LinearDiagramOptions,
) -> GbdrawConfig:
    """Resolve one typed config and losslessly attach explicit label tables."""

    cfg = apply_config_overrides(options.config, options.config_overrides)

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
        params = slot.params or {}
        nt = str(params.get("nt", params.get("dinucleotide", default_nt)) or default_nt).upper()
        if len(nt) >= 2:
            nts.add(nt)
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
    return {nt: circular_dinucleotide_content_df(record, window, step, nt) for nt in sorted(nts)}


def _build_circular_dinucleotide_skew_dataframes(
    record: SeqRecord,
    *,
    window: int,
    step: int,
    nts: set[str],
) -> dict[str, DataFrame]:
    return {nt: skew_df(record, window, step, nt) for nt in sorted(nts)}


def _has_conservation_inputs(
    conservation_blast_files: Sequence[str] | None,
    conservation_dataframes: Sequence[DataFrame] | None,
) -> bool:
    return bool(conservation_blast_files or conservation_dataframes)


def _build_conservation_blast_config(
    records: Sequence[SeqRecord],
    *,
    evalue: float,
    bitscore: float,
    identity: float,
    alignment_length: int,
    default_colors: DataFrame,
    profile: CircularRenderProfile,
) -> BlastMatchConfigurator:
    return BlastMatchConfigurator(
        evalue=float(evalue),
        bitscore=float(bitscore),
        identity=float(identity),
        alignment_length=int(alignment_length),
        sequence_length_dict=create_dict_for_sequence_lengths(records),
        profile=profile,
        default_colors_df=default_colors,
    )


def _load_conservation_result(
    records: Sequence[SeqRecord],
    *,
    conservation_blast_files: Sequence[str] | None,
    conservation_dataframes: Sequence[DataFrame] | None,
    conservation_labels: Sequence[str] | None,
    conservation_colors: Sequence[str] | None,
    evalue: float,
    bitscore: float,
    identity: float,
    alignment_length: int,
    default_colors: DataFrame,
    profile: CircularRenderProfile,
) -> ConservationLoadResult | None:
    if not _has_conservation_inputs(conservation_blast_files, conservation_dataframes):
        return None
    if alignment_length < 0:
        raise ValidationError("alignment_length must be >= 0")
    blast_config = _build_conservation_blast_config(
        records,
        evalue=evalue,
        bitscore=bitscore,
        identity=identity,
        alignment_length=alignment_length,
        default_colors=default_colors,
        profile=profile,
    )
    return load_conservation_sources(
        blast_config=blast_config,
        conservation_files=conservation_blast_files,
        conservation_dataframes=conservation_dataframes,
        labels=conservation_labels,
        colors=conservation_colors,
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
    if axis_index is not None:
        out = list(base_slots)
        out[axis_index:axis_index] = conservation_slots
        return out
    out: list[CircularTrackSlot] = []
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
    resolved_window = window
    resolved_step = step
    seq_length = len(record.seq)

    if resolved_window is None:
        if seq_length < 1_000_000:
            resolved_window = cfg.objects.sliding_window.default[0]
        elif seq_length < 10_000_000:
            resolved_window = cfg.objects.sliding_window.up1m[0]
        else:
            resolved_window = cfg.objects.sliding_window.up10m[0]

    if resolved_step is None:
        if seq_length < 1_000_000:
            resolved_step = cfg.objects.sliding_window.default[1]
        elif seq_length < 10_000_000:
            resolved_step = cfg.objects.sliding_window.up1m[1]
        else:
            resolved_step = cfg.objects.sliding_window.up10m[1]

    return int(resolved_window), int(resolved_step)


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
    if depth_config.tick_font_size is not None and float(depth_config.tick_font_size) <= 0:
        raise ValidationError("depth_tick_font_size must be > 0")


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
        ("tick_font_size", "gc_content_tick_font_size"),
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
    return (
        int(depth_window) if depth_window is not None else max(100, int(window) // 10),
        int(depth_step) if depth_step is not None else max(1, int(step) // 10),
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
        sum(record_bounds[index].width for index in row)
        + max(0, len(row) - 1) * float(column_gap_px)
        for row in rows
    )
    row_heights = tuple(
        max(record_bounds[index].height for index in row)
        for row in rows
    )
    grid_width = max(row_widths)
    grid_height = sum(row_heights) + max(0, len(rows) - 1) * float(row_gap_px)

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
) -> None:
    """Repair copied multi-record subtrees with one shared reference namespace."""

    renamed: dict[str, str] = {}
    for element in elements:
        for child in _walk_svgwrite_elements(element):
            attribs = getattr(child, "attribs", None)
            if not isinstance(attribs, dict):
                continue
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


def _build_circular_plot_title_group(
    *,
    gb_record: SeqRecord,
    output_prefix: str,
    cfg: GbdrawConfig,
    profile: CircularRenderProfile,
    species: str | None,
    strain: str | None,
    plot_title: str | None,
) -> tuple[Group, Aabb]:
    plot_title_canvas_config = CircularCanvasConfigurator(
        output_prefix=output_prefix,
        profile=profile,
        legend="none",
        gb_record=gb_record,
    )
    definition = DefinitionGroup(
        gb_record=gb_record,
        canvas_config=plot_title_canvas_config,
        species=species,
        strain=strain,
        plot_title=plot_title,
        definition_profile="shared_common",
        definition_group_id="plot_title",
        cfg=cfg,
    )
    return definition.get_group(), definition.local_bounds


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


def _has_mixed_short_and_long_records(
    record_lengths: Sequence[int],
    *,
    length_threshold: int,
) -> bool:
    """Return whether record lengths span both short and long buckets."""
    has_short = False
    has_long = False
    threshold = int(length_threshold)
    for record_length in record_lengths:
        if int(record_length) < threshold:
            has_short = True
        else:
            has_long = True
        if has_short and has_long:
            return True
    return False


def _harmonize_multi_record_circular_style_cfg(
    cfg: GbdrawConfig,
    *,
    record_lengths: Sequence[int],
) -> GbdrawConfig:
    """Harmonize short-record feature/axis style to long settings for mixed multi-record canvases."""
    if not _has_mixed_short_and_long_records(
        record_lengths,
        length_threshold=cfg.labels.length_threshold.circular,
    ):
        return cfg

    circular_cfg = cfg.canvas.circular
    harmonized_track_ratio_factors = {
        str(key): [float(value) for value in list(values)]
        for key, values in circular_cfg.track_ratio_factors.items()
    }
    short_factors = list(harmonized_track_ratio_factors.get("short", []))
    long_factors = list(harmonized_track_ratio_factors.get("long", []))
    if short_factors and long_factors:
        for factor_index in (0, 1, 2):
            if factor_index < len(short_factors) and factor_index < len(long_factors):
                short_factors[factor_index] = float(long_factors[factor_index])
        harmonized_track_ratio_factors["short"] = short_factors

    harmonized_track_dict: dict[str, dict[str, dict[str, float]]] = {
        str(length_param): {
            str(track_type): {
                str(track_id): float(track_value)
                for track_id, track_value in track_values.items()
            }
            for track_type, track_values in track_type_values.items()
        }
        for length_param, track_type_values in circular_cfg.track_dict.items()
    }
    short_track_dict = harmonized_track_dict.get("short")
    long_track_dict = harmonized_track_dict.get("long")
    if short_track_dict is not None and long_track_dict is not None:
        for track_type, short_track_values in short_track_dict.items():
            long_track_values = long_track_dict.get(str(track_type))
            if long_track_values is None:
                continue
            for track_id in ("2", "3"):
                if track_id in short_track_values and track_id in long_track_values:
                    short_track_values[track_id] = float(long_track_values[track_id])

    harmonized_features = replace(
        cfg.objects.features,
        block_stroke_width=replace(
            cfg.objects.features.block_stroke_width,
            short=float(cfg.objects.features.block_stroke_width.long),
        ),
        line_stroke_width=replace(
            cfg.objects.features.line_stroke_width,
            short=float(cfg.objects.features.line_stroke_width.long),
        ),
    )
    harmonized_axis_circular = replace(
        cfg.objects.axis.circular,
        stroke_width=replace(
            cfg.objects.axis.circular.stroke_width,
            short=float(cfg.objects.axis.circular.stroke_width.long),
        ),
    )
    harmonized_axis = replace(cfg.objects.axis, circular=harmonized_axis_circular)
    harmonized_objects = replace(
        cfg.objects,
        features=harmonized_features,
        axis=harmonized_axis,
    )
    harmonized_label_font_size = replace(
        cfg.labels.font_size,
        short=float(cfg.labels.font_size.long),
    )
    harmonized_labels = replace(
        cfg.labels,
        font_size=harmonized_label_font_size,
    )
    harmonized_circular = replace(
        circular_cfg,
        track_ratio_factors=harmonized_track_ratio_factors,
        track_dict=harmonized_track_dict,
    )
    harmonized_canvas = replace(cfg.canvas, circular=harmonized_circular)
    return replace(cfg, canvas=harmonized_canvas, objects=harmonized_objects, labels=harmonized_labels)


def _resolve_multi_record_tick_track_channel_override(
    record_lengths: Sequence[int],
    *,
    length_threshold: int,
) -> Literal["long"] | None:
    """Return tick ratio channel override for mixed short/long multi-record canvases."""
    if _has_mixed_short_and_long_records(
        record_lengths,
        length_threshold=length_threshold,
    ):
        return "long"
    return None


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
    if default_colors is None:
        default_colors = load_default_colors(
            user_defined_default_colors=default_colors_file or "",
            palette=default_colors_palette or "default",
            load_comparison=load_comparison_colors,
        )
    return resolve_feature_inputs(
        color_table=color_table,
        default_colors=default_colors,
        feature_visibility_table=feature_visibility_table,
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
) -> tuple[str, ...]:
    def label(record: SeqRecord, fallback: str) -> str:
        annotations = getattr(record, "annotations", {}) or {}
        return (
            str(annotations.get("gbdraw_record_label") or "").strip()
            or str(record.id or "").strip()
            or fallback
        )

    return tuple(
        (
            f"{_web_normalize_cache_label(label(records[index], f'seq_{index + 1}'), f'seq_{index + 1}')}"
            f".{_web_normalize_cache_label(label(records[index + 1], f'seq_{index + 2}'), f'seq_{index + 2}')}"
            ".losatp.tsv"
        )
        for index in range(max(0, len(records) - 1))
    )


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
    orthogroup_member_max_hits: int,
    max_paralog_links_per_orthogroup: int,
    evalue: float,
    bitscore: float,
    identity: float,
    alignment_length: int,
    collinearity_params: LosslessCollinearityParameters | None,
    collinearity_unit_mode: CollinearityUnitMode | str,
    collinearity_anchor_mode: CollinearityAnchorMode,
    collinearity_search_scope: CollinearitySearchScope,
    collinearity_comparison_pairs: Sequence[tuple[int, int]] | None,
    losatp_cache: LosatpCacheManager | None,
    protein_extraction: ProteinExtractionResult | None,
    feature_visibility_rules: list[dict[str, object]] | None,
    cache_filenames: Sequence[str] | None,
) -> ProteinBlastpResult | CollinearityResult:
    """Translate resolved typed values into one real analysis-helper call."""

    if mode == "pairwise":
        return build_pairwise_protein_blastp_comparisons(
            records,
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
    if mode == "orthogroup":
        return build_rbh_orthogroup_protein_blastp_comparisons(
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
    if mode == "collinear":
        return build_orthogroup_collinearity_blocks(
            records,
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
    raise ValidationError(f"Unsupported protein analysis mode: {mode!r}.")


def assemble_linear_diagram_from_records(
    records: Sequence[SeqRecord],
    *,
    cfg: GbdrawConfig,
    blast_files: Optional[Sequence[str]] = None,
    linear_comparisons: Sequence[LinearComparison] | None = None,
    layout: LinearMultiRecordOptions | None = None,
    protein_comparisons: Sequence[DataFrame] | None = None,
    orthogroups: OrthogroupResult | None = None,
    protein_blastp_mode: ProteinBlastpMode | str = "none",
    protein_comparison_pairs: Sequence[tuple[int, int]] | None = None,
    pairwise_match_style: Literal["ribbon", "curve"] | str = "ribbon",
    collinearity_blocks: CollinearityResult | Sequence[CollinearityBlock] | None = None,
    collinearity_params: LosslessCollinearityParameters | None = None,
    collinearity_unit_mode: CollinearityUnitMode | str = "auto",
    collinearity_anchor_mode: CollinearityAnchorMode | str = "rbh",
    collinearity_search_scope: CollinearitySearchScope | str = "adjacent",
    collinearity_color_mode: CollinearityColorMode | str = "orientation",
    losatp_bin: str = "losat",
    ncbi_blastp_bin: str | None = None,
    losatp_threads: int | None = None,
    protein_blastp_max_hits: int = 5,
    protein_blastp_candidate_limit: int | None = None,
    losatp_cache: LosatpCacheManager | None = None,
    protein_extraction: ProteinExtractionResult | None = None,
    orthogroup_membership_mode: OrthogroupMembershipMode | str = "anchor_core_v1",
    orthogroup_member_max_hits: int = 5,
    collinear_max_paralog_links_per_orthogroup: int = 2,
    align_orthogroup_feature: str | None = None,
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
    _return_build_result: bool = False,
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
    normalized_protein_pairs: tuple[tuple[int, int], ...] | None = None
    if protein_comparison_pairs is not None:
        normalized_pairs: list[tuple[int, int]] = []
        for pair in protein_comparison_pairs:
            if (
                not isinstance(pair, Sequence)
                or len(pair) != 2
                or not all(isinstance(index, int) and not isinstance(index, bool) for index in pair)
            ):
                raise ValidationError("protein_comparison_pairs must contain integer index pairs")
            query_index, subject_index = int(pair[0]), int(pair[1])
            if query_index < 0 or subject_index < 0 or query_index >= len(records) or subject_index >= len(records):
                raise ValidationError("protein_comparison_pairs contains an out-of-range record index")
            if query_index == subject_index:
                raise ValidationError("protein_comparison_pairs cannot compare a record to itself")
            normalized_pairs.append((query_index, subject_index))
        if len(set(normalized_pairs)) != len(normalized_pairs):
            raise ValidationError("protein_comparison_pairs must not contain duplicates")
        normalized_protein_pairs = tuple(normalized_pairs)
        _ordered_indices, pair_rows = resolve_record_row_positions(
            records,
            layout.multi_record_positions if layout is not None else None,
        )
        for query_index, subject_index in normalized_protein_pairs:
            if abs(pair_rows[query_index] - pair_rows[subject_index]) != 1:
                raise ValidationError(
                    "protein_comparison_pairs must connect records in adjacent rows: "
                    f"query=#{query_index + 1} row={pair_rows[query_index] + 1}, "
                    f"subject=#{subject_index + 1} row={pair_rows[subject_index] + 1}."
                )
    normalized_protein_blastp_mode = normalize_protein_blastp_mode(protein_blastp_mode)
    if normalized_protein_pairs is not None and normalized_protein_blastp_mode != "pairwise":
        raise ValidationError("protein_comparison_pairs requires protein_blastp_mode='pairwise'")
    if normalized_protein_pairs is not None and linear_comparisons:
        raise ValidationError("Pass either protein_comparison_pairs or linear_comparisons, not both")
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
    if int(protein_blastp_max_hits) <= 0:
        raise ValidationError("protein_blastp_max_hits must be > 0")
    if int(orthogroup_member_max_hits) <= 0:
        raise ValidationError("orthogroup_member_max_hits must be > 0")
    if int(collinear_max_paralog_links_per_orthogroup) <= 0:
        raise ValidationError("collinear_max_paralog_links_per_orthogroup must be > 0")
    if losatp_threads is not None and int(losatp_threads) <= 0:
        raise ValidationError("losatp_threads must be > 0 or None")
    if protein_blastp_candidate_limit is not None and int(protein_blastp_candidate_limit) <= 0:
        raise ValidationError("protein_blastp_candidate_limit must be > 0 or None")
    if normalized_protein_blastp_mode != "none" and protein_comparisons is not None:
        raise ValidationError("Pass either protein_blastp_mode or protein_comparisons, not both.")
    if collinearity_blocks is not None and (
        normalized_protein_blastp_mode != "none" or protein_comparisons is not None or blast_files
    ):
        raise ValidationError(
            "Pass collinearity_blocks without protein_blastp_mode, protein_comparisons, or blast_files."
        )
    if normalized_protein_blastp_mode != "none" and blast_files:
        raise ValidationError("protein_blastp_mode cannot be used with blast_files.")
    if normalized_protein_blastp_mode != "none" and len(records) < 2:
        raise ValidationError("protein_blastp_mode requires at least two records")
    has_precomputed_comparisons = bool(
        blast_files
        or linear_comparisons
        or protein_comparisons is not None
        or collinearity_blocks is not None
    )
    if (
        align_orthogroup_feature
        and normalized_protein_blastp_mode != "orthogroup"
        and not has_precomputed_comparisons
    ):
        raise ValidationError(
            "align_orthogroup_feature requires protein_blastp_mode='orthogroup'."
        )
    _validate_positive_optional("depth_window", depth_window)
    _validate_positive_optional("depth_step", depth_step)
    has_comparisons = bool(
        blast_files
        or linear_comparisons
        or protein_comparisons
        or collinearity_blocks
        or normalized_protein_blastp_mode != "none"
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
    profile = LinearRenderProfile(cfg)

    if selected_features_set is None:
        selected_features_set = DEFAULT_SELECTED_FEATURES
    resolved_protein_comparisons: list[DataFrame] | None = None
    resolved_linear_comparisons: list[LinearComparison] = list(linear_comparisons or ())
    resolved_orthogroups: OrthogroupResult | None = orthogroups
    resolved_collinearity_result: CollinearityResult | None = None
    losat_cache_filenames = _linear_losat_cache_filenames(records)

    def project_collinearity_comparisons(
        result: CollinearityResult,
    ) -> list[DataFrame] | None:
        if collinearity_comparison_pairs is None:
            return convert_collinearity_blocks_to_comparisons(
                result,
                records=records,
                color_mode=normalized_collinearity_color_mode,
                search_scope=normalized_collinearity_search_scope,
            )
        pair_comparisons = convert_collinearity_blocks_to_pair_comparisons(
            result,
            records=records,
            color_mode=normalized_collinearity_color_mode,
            search_scope=normalized_collinearity_search_scope,
        )
        allowed_pairs = set(collinearity_comparison_pairs)
        resolved_linear_comparisons.extend(
            LinearComparison(query_index, subject_index, matches)
            for (query_index, subject_index), matches in pair_comparisons.items()
            if (query_index, subject_index) in allowed_pairs
        )
        return None

    def invoke_protein_analysis(
        mode: ProteinBlastpMode,
        analysis_records: Sequence[SeqRecord],
        *,
        analysis_extraction: ProteinExtractionResult | None = protein_extraction,
        analysis_cache_filenames: Sequence[str] | None = losat_cache_filenames,
    ) -> ProteinBlastpResult | CollinearityResult:
        return _invoke_protein_analysis_helper(
            mode,
            analysis_records,
            losatp_bin=losatp_bin,
            ncbi_blastp_bin=ncbi_blastp_bin,
            losatp_threads=losatp_threads,
            pairwise_max_hits=int(protein_blastp_max_hits),
            candidate_limit=protein_blastp_candidate_limit,
            orthogroup_membership_mode=normalized_orthogroup_membership_mode,
            orthogroup_member_max_hits=int(orthogroup_member_max_hits),
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
            collinearity_comparison_pairs=collinearity_comparison_pairs,
            losatp_cache=losatp_cache,
            protein_extraction=analysis_extraction,
            feature_visibility_rules=feature_visibility_rules,
            cache_filenames=analysis_cache_filenames,
        )

    if protein_comparisons is not None:
        resolved_protein_comparisons = list(protein_comparisons)
    elif collinearity_blocks is not None:
        if isinstance(collinearity_blocks, CollinearityResult):
            collinearity_result = collinearity_blocks
        else:
            collinearity_result = CollinearityResult(blocks=tuple(collinearity_blocks))
        resolved_collinearity_result = collinearity_result
        if collinearity_result.orthogroups is not None:
            resolved_orthogroups = collinearity_result.orthogroups
        resolved_protein_comparisons = project_collinearity_comparisons(
            collinearity_result
        )
    elif normalized_protein_blastp_mode == "pairwise":
        pair_inputs = normalized_protein_pairs
        if pair_inputs is not None:
            for query_index, subject_index in pair_inputs:
                pair_extraction = (
                    replace(
                        protein_extraction,
                        proteins_by_record=[
                            protein_extraction.proteins_by_record[query_index],
                            protein_extraction.proteins_by_record[subject_index],
                        ],
                    )
                    if protein_extraction is not None
                    else None
                )
                protein_blastp_result = cast(
                    ProteinBlastpResult,
                    invoke_protein_analysis(
                        "pairwise",
                        (records[query_index], records[subject_index]),
                        analysis_extraction=pair_extraction,
                        analysis_cache_filenames=_linear_losat_cache_filenames(
                            (records[query_index], records[subject_index])
                        ),
                    ),
                )
                resolved_linear_comparisons.append(
                    LinearComparison(
                        query_index,
                        subject_index,
                        protein_blastp_result.comparisons[0],
                    )
                )
        else:
            protein_blastp_result = cast(
                ProteinBlastpResult,
                invoke_protein_analysis("pairwise", records),
            )
            resolved_protein_comparisons = protein_blastp_result.comparisons
    elif normalized_protein_blastp_mode == "orthogroup":
        protein_blastp_result = cast(
            ProteinBlastpResult,
            invoke_protein_analysis("orthogroup", records),
        )
        resolved_protein_comparisons = protein_blastp_result.comparisons
        resolved_orthogroups = protein_blastp_result.orthogroups
    elif normalized_protein_blastp_mode == "collinear":
        collinearity_result = cast(
            CollinearityResult,
            invoke_protein_analysis("collinear", records),
        )
        resolved_collinearity_result = collinearity_result
        resolved_orthogroups = collinearity_result.orthogroups
        resolved_protein_comparisons = project_collinearity_comparisons(
            collinearity_result
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
    if window is None:
        if longest_genome < 1_000_000:
            window = cfg.objects.sliding_window.default[0]
        elif longest_genome < 10_000_000:
            window = cfg.objects.sliding_window.up1m[0]
        else:
            window = cfg.objects.sliding_window.up10m[0]

    if step is None:
        if longest_genome < 1_000_000:
            step = cfg.objects.sliding_window.default[1]
        elif longest_genome < 10_000_000:
            step = cfg.objects.sliding_window.up1m[1]
        else:
            step = cfg.objects.sliding_window.up10m[1]
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
    blast_config.collinearity_color_mode = normalized_collinearity_color_mode
    additional_color_modes = (
        {normalized_collinearity_color_mode}
        if collinearity_blocks is not None or normalized_protein_blastp_mode == "collinear"
        else None
    )
    configure_pairwise_identity_legend_from_comparisons(
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
        color_table=color_table,
        default_colors=default_colors,
        selected_features_set=list(selected_features_set),
        profile=profile,
        gc_config=gc_config,
        skew_config=skew_config,
        feature_config=feature_config,
        blast_config=blast_config,
        canvas_config=canvas_config,
    )

    canvas = assemble_linear_diagram(
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
        annotations=annotation_options,
        plot_title=normalized_plot_title or None,
        plot_title_position=normalized_plot_title_position,
        plot_title_font_size=resolved_plot_title_font_size,
        comparison_dataframes=resolved_protein_comparisons,
        linear_comparisons=resolved_linear_comparisons or None,
        linear_layout=layout,
        orthogroups=resolved_orthogroups,
        align_orthogroup_feature=align_orthogroup_feature,
    )
    build_result = LinearDiagramBuildResult(
        drawing=canvas,
        metadata=LinearDiagramMetadata(
            protein_comparisons=(
                tuple(resolved_protein_comparisons)
                if resolved_protein_comparisons is not None
                else None
            ),
            linear_comparisons=tuple(resolved_linear_comparisons),
            orthogroups=resolved_orthogroups,
            collinearity_result=resolved_collinearity_result,
        ),
    )
    return build_result if _return_build_result else canvas


def assemble_circular_diagram_from_record(
    gb_record: SeqRecord,
    *,
    cfg: GbdrawConfig,
    conservation_blast_files: Sequence[str] | None = None,
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
    _tick_track_channel_override: Literal["short", "long"] | None = None,
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
) -> Drawing:
    """Builds and assembles a circular diagram for a single record.

    If default_colors is None, it loads the built-in default palette.
    If color_table is None and color_table_file is provided, it is loaded.
    If selected_features_set is None, it uses the CLI default feature list.
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
        if show_plot_title and keep_full_definition_with_plot_title:
            effective_definition_profile = "full"
        elif show_plot_title:
            effective_definition_profile = "record_summary"
        else:
            effective_definition_profile = "full"

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
    profile = CircularRenderProfile(cfg)

    conservation_mode = normalize_conservation_reference(conservation_reference)
    if _precomputed_conservation_tracks is not None:
        conservation_tracks = tuple(_precomputed_conservation_tracks)
    else:
        conservation_load_result = _load_conservation_result(
            [gb_record],
            conservation_blast_files=conservation_blast_files,
            conservation_dataframes=conservation_dataframes,
            conservation_labels=conservation_labels,
            conservation_colors=conservation_colors,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
            default_colors=default_colors,
            profile=profile,
        )
        conservation_tracks = (
            normalize_conservation_tracks_for_record(
                conservation_load_result,
                displayed_records=[gb_record],
                record=gb_record,
                conservation_reference=conservation_mode,
            )
            if conservation_load_result is not None
            else ()
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

    seq_length = len(gb_record.seq)

    # Match legacy CLI behavior: circular window/step are based on the record length.
    if window is None:
        if seq_length < 1_000_000:
            window = cfg.objects.sliding_window.default[0]
        elif seq_length < 10_000_000:
            window = cfg.objects.sliding_window.up1m[0]
        else:
            window = cfg.objects.sliding_window.up10m[0]

    if step is None:
        if seq_length < 1_000_000:
            step = cfg.objects.sliding_window.default[1]
        elif seq_length < 10_000_000:
            step = cfg.objects.sliding_window.up1m[1]
        else:
            step = cfg.objects.sliding_window.up10m[1]
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
    )
    legend_config = LegendDrawingConfigurator(
        color_table=color_table,
        default_colors=default_colors,
        selected_features_set=list(selected_features_set),
        profile=profile,
        gc_config=gc_config,
        skew_config=skew_config,
        feature_config=feature_config,
        canvas_config=canvas_config,
    )

    title_target: Group | None = None
    title_bounds: Aabb | None = None
    title_placement = TitlePlacement.NONE
    if show_plot_title:
        title_target, title_bounds = _build_circular_plot_title_group(
            gb_record=gb_record,
            output_prefix=output_prefix,
            cfg=cfg,
            profile=profile,
            species=species,
            strain=strain,
            plot_title=normalized_plot_title or None,
        )
        title_placement = TitlePlacement(normalized_plot_title_position)

    result = _assemble_circular_diagram_result(
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
        annotations=_resolved_annotations or annotation_options,
        annotation_record_index=_annotation_record_index,
        dinucleotide_content_dataframes=dinucleotide_content_dataframes,
        dinucleotide_skew_dataframes=dinucleotide_skew_dataframes,
        definition_position="center",
        definition_profile=effective_definition_profile,
        definition_group_id=_definition_group_id,
        center_reserved_radius=center_reserved_radius,
        _tick_track_channel_override=_tick_track_channel_override,
        title_target=title_target,
        title_bounds=title_bounds,
        title_placement=title_placement,
    )
    return result.drawing


def assemble_circular_diagram_from_records(
    records: Sequence[SeqRecord],
    *,
    cfg: GbdrawConfig,
    conservation_blast_files: Sequence[str] | None = None,
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
    resolved_annotations = resolve_annotations(annotation_options, records, mode="circular")
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
    record_definition_profile: Literal["full", "record_summary"] = "record_summary"
    if len(records) == 1 and not show_plot_title:
        record_definition_profile = "full"
    elif show_plot_title and keep_full_definition_with_plot_title:
        record_definition_profile = "full"
    record_lengths = [len(record.seq) for record in records]
    cfg = _harmonize_multi_record_circular_style_cfg(
        cfg,
        record_lengths=record_lengths,
    )
    profile = CircularRenderProfile(cfg)
    conservation_mode = normalize_conservation_reference(conservation_reference)
    conservation_load_result = _load_conservation_result(
        records,
        conservation_blast_files=conservation_blast_files,
        conservation_dataframes=conservation_dataframes,
        conservation_labels=conservation_labels,
        conservation_colors=conservation_colors,
        evalue=evalue,
        bitscore=bitscore,
        identity=identity,
        alignment_length=alignment_length,
        default_colors=default_colors,
        profile=profile,
    )
    first_record_conservation_tracks = (
        normalize_conservation_tracks_for_record(
            conservation_load_result,
            displayed_records=records,
            record=records[0],
            conservation_reference=conservation_mode,
        )
        if conservation_load_result is not None
        else ()
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
    tick_track_channel_override = _resolve_multi_record_tick_track_channel_override(
        record_lengths,
        length_threshold=cfg.labels.length_threshold.circular,
    )
    record_scales = _resolve_multi_record_scales(
        record_lengths,
        mode=normalized_multi_record_size_mode,
        min_radius_ratio=normalized_multi_record_min_radius_ratio,
    )
    record_depth_track_data: list[list[DepthTrackData]] = [[] for _ in records]
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
        record_conservation_tracks = (
            normalize_conservation_tracks_for_record(
                conservation_load_result,
                displayed_records=records,
                record=record,
                conservation_reference=conservation_mode,
            )
            if conservation_load_result is not None
            else ()
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
            _tick_track_channel_override=tick_track_channel_override,
            _precomputed_depth_tracks=record_depth_track_data[record_index],
            _precomputed_depth_track_count=available_depth_track_count,
            _precomputed_conservation_tracks=record_conservation_tracks,
            _resolved_feature_inputs=resolved_feature_inputs,
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

    legend_table: dict = {}
    legend_target: Group | None = None
    legend_bounds: Aabb | None = None
    plot_title_group: Group | None = None
    plot_title_bounds: Aabb | None = None

    if show_plot_title:
        plot_title_group, plot_title_bounds = _build_circular_plot_title_group(
            gb_record=records[0],
            output_prefix=output_prefix,
            cfg=cfg,
            profile=profile,
            species=species,
            strain=strain,
            plot_title=normalized_plot_title or None,
        )

    if legend_effective != "none":
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
        )
        color_map = feature_config.specific_color_rules
        default_color_map = feature_config.default_color_map
        features_present = check_feature_presence(
            list(records),
            list(selected_features_set),
            feature_visibility_rules=feature_config.feature_visibility_rules,
            specific_color_rules=color_map,
        )
        used_color_rules, default_used_features = precompute_used_color_rules(
            list(records),
            color_map,
            default_color_map,
            set(feature_config.selected_features_set),
            feature_visibility_rules=feature_config.feature_visibility_rules,
        )
        legend_table = prepare_legend_table(
            gc_config,
            skew_config,
            feature_config,
            features_present,
            used_color_rules=used_color_rules,
            default_used_features=default_used_features,
            depth_config=depth_config if depth_track_data_count(record_depth_track_data) == 1 else None,
            show_gc=profile.show_gc,
            show_skew=profile.show_skew,
            show_depth=bool(
                profile.show_depth
                and depth_track_data_count(record_depth_track_data) == 1
            ),
        )
        if profile.show_depth:
            legend_table = sync_depth_track_legend_entries(
                legend_table,
                representative_depth_tracks(record_depth_track_data),
            )
        if first_record_conservation_tracks:
            if any(track.track_color for track in first_record_conservation_tracks):
                for track in first_record_conservation_tracks:
                    min_color, max_color = conservation_track_gradient_colors(
                        track.track_color,
                        default_min_color=cfg.objects.conservation.min_color,
                        default_max_color=cfg.objects.conservation.max_color,
                    )
                    legend_table[track.track_label] = {
                        "type": "gradient",
                        "min_color": min_color,
                        "max_color": max_color,
                        "stroke": "none",
                        "width": 0,
                        "min_value": float(identity),
                    }
            else:
                legend_table["Conservation identity"] = {
                    "type": "gradient",
                    "min_color": cfg.objects.conservation.min_color,
                    "max_color": cfg.objects.conservation.max_color,
                    "stroke": "none",
                    "width": 0,
                    "min_value": float(identity),
                }
        if legend_table:
            legend_config = LegendDrawingConfigurator(
                color_table=color_table,
                default_colors=default_colors,
                selected_features_set=list(selected_features_set),
                profile=profile,
                gc_config=gc_config,
                skew_config=skew_config,
                feature_config=feature_config,
                canvas_config=legend_canvas_config,
            )
            legend_measurement = legend_config.measure_legend(
                legend_table,
                placement=legend_canvas_config.legend_position,
                wrap_width=grid_layout.bounds.width,
            )
            legend_bounds = legend_measurement.local_bounds
            legend_target = LegendGroup(
                legend_canvas_config,
                legend_measurement,
                legend_table,
            ).get_group()

    merged_canvas = Drawing(
        filename=f"{output_prefix}.svg",
        size=(f"{grid_layout.bounds.width}px", f"{grid_layout.bounds.height}px"),
        viewBox=(
            f"0 0 {grid_layout.bounds.width} {grid_layout.bounds.height}"
        ),
        debug=False,
    )
    used_ids: set[str] = set()
    record_targets: list[Group] = []
    grid_overlay_obstacles: list[Aabb] = []

    for record_index, result in enumerate(record_results):
        record_group = Group(id=f"record_{record_index}", debug=False)
        record_group.attribs["data-gbdraw-record-id"] = str(records[record_index].id)
        record_group.attribs["data-gbdraw-record-index"] = str(record_index)
        record_dx, record_dy = grid_layout.translations[record_index]
        record_group.translate(record_dx, record_dy)
        used_ids.add(f"record_{record_index}")

        copied_definitions: list[object] = []
        copied_elements: list[object] = []
        for element in result.drawing.elements:
            copied = copy.deepcopy(element)
            if _is_defs_element(element):
                copied_definitions.extend(
                    list(getattr(copied, "elements", ()) or ())
                )
                continue
            _strip_nested_composition_roles(copied)
            _suffix_fixed_top_level_group_id(copied, record_index)
            copied_elements.append(copied)

        for definition in copied_definitions:
            _strip_nested_composition_roles(definition)

        _uniquify_copied_subtrees_ids(
            (*copied_definitions, *copied_elements),
            record_index=record_index,
            used_ids=used_ids,
        )
        for definition in copied_definitions:
            merged_canvas.defs.add(definition)
        for copied in copied_elements:
            record_group.add(copied)
        merged_canvas.add(record_group)
        record_targets.append(record_group)
        grid_overlay_obstacles.extend(
            obstacle.translated(record_dx, record_dy)
            for obstacle in result.overlay_obstacles
        )

    if plot_title_group is not None:
        merged_canvas.add(plot_title_group)
    if legend_target is not None:
        merged_canvas.add(legend_target)

    legend_placement = LegendPlacement(str(legend_effective))
    title_placement = (
        TitlePlacement(normalized_plot_title_position)
        if plot_title_group is not None
        else TitlePlacement.NONE
    )
    composition_plan = plan_composition(
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
        )
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
) -> Drawing:
    """Build a circular diagram using mode-specific typed options."""

    options = _circular_builder_options(options)
    depth_table, depth_file = _resolve_single_circular_depth_options(options)
    colors = options.colors
    output = options.output
    tracks = options.tracks

    cfg = _resolve_diagram_options_config(options)

    return assemble_circular_diagram_from_record(
        gb_record,
        cfg=cfg,
        color_table=colors.color_table if colors else None,
        color_table_file=colors.color_table_file if colors else None,
        default_colors=colors.default_colors if colors else None,
        default_colors_palette=colors.default_colors_palette if colors else "default",
        default_colors_file=colors.default_colors_file if colors else None,
        selected_features_set=options.selected_features_set,
        feature_visibility_table=options.feature_visibility_table,
        feature_visibility_table_file=options.feature_visibility_table_file,
        feature_shapes=options.feature_shapes,
        output_prefix="out",
        legend=output.legend if output else "right",
        dinucleotide=options.dinucleotide,
        window=options.window,
        step=options.step,
        depth_window=options.depth_window,
        depth_step=options.depth_step,
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
        conservation_blast_files=options.conservation_blast_files,
        conservation_dataframes=options.conservation_dataframes,
        conservation_reference=options.conservation_reference,
        conservation_labels=options.conservation_labels,
        conservation_colors=options.conservation_colors,
        conservation_ring_width=options.conservation_ring_width,
        conservation_ring_gap=options.conservation_ring_gap,
        evalue=options.evalue,
        bitscore=options.bitscore,
        identity=options.identity,
        alignment_length=options.alignment_length,
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
        annotation_options=options.annotations,
        _precomputed_depth_track_specs=_precomputed_depth_track_specs,
        _precomputed_depth_track_count=_precomputed_depth_track_count,
        _resolved_feature_inputs=_resolved_feature_inputs,
    )


def _build_linear_diagram(
    records: Sequence[SeqRecord],
    *,
    options: LinearDiagramOptions | None = None,
    layout: LinearMultiRecordOptions | None = None,
    losatp_cache: LosatpCacheManager | None = None,
    protein_extraction: ProteinExtractionResult | None = None,
    _resolved_feature_inputs: ResolvedFeatureInputs | None = None,
    _return_build_result: bool = False,
) -> Drawing | LinearDiagramBuildResult:
    """Build a linear diagram using mode-specific typed options."""

    options = _linear_builder_options(options)
    colors = options.colors
    output = options.output
    tracks = options.tracks
    cfg = _resolve_diagram_options_config(options)
    normalized_collinearity_anchor_mode = normalize_collinearity_anchor_mode(
        str(options.collinearity_anchor_mode)
    )

    return assemble_linear_diagram_from_records(
        records,
        cfg=cfg,
        blast_files=options.blast_files,
        linear_comparisons=options.linear_comparisons,
        layout=layout,
        protein_comparisons=options.protein_comparisons,
        orthogroups=options.orthogroups,
        protein_blastp_mode=options.protein_blastp_mode,
        protein_comparison_pairs=options.protein_comparison_pairs,
        pairwise_match_style=options.pairwise_match_style,
        collinearity_blocks=options.collinearity_blocks,
        collinearity_params=options.collinearity_params,
        collinearity_unit_mode=options.collinearity_unit_mode,
        collinearity_anchor_mode=normalized_collinearity_anchor_mode,
        collinearity_search_scope=options.collinearity_search_scope,
        collinearity_color_mode=options.collinearity_color_mode,
        losatp_bin=options.losatp_bin,
        ncbi_blastp_bin=options.ncbi_blastp_bin,
        losatp_threads=options.losatp_threads,
        protein_blastp_max_hits=options.protein_blastp_max_hits,
        protein_blastp_candidate_limit=options.protein_blastp_candidate_limit,
        losatp_cache=losatp_cache,
        protein_extraction=protein_extraction,
        orthogroup_membership_mode=options.orthogroup_membership_mode,
        orthogroup_member_max_hits=options.orthogroup_member_max_hits,
        collinear_max_paralog_links_per_orthogroup=options.collinear_max_paralog_links_per_orthogroup,
        align_orthogroup_feature=options.align_orthogroup_feature,
        color_table=colors.color_table if colors else None,
        color_table_file=colors.color_table_file if colors else None,
        default_colors=colors.default_colors if colors else None,
        default_colors_palette=colors.default_colors_palette if colors else "default",
        default_colors_file=colors.default_colors_file if colors else None,
        selected_features_set=options.selected_features_set,
        feature_visibility_table=options.feature_visibility_table,
        feature_visibility_table_file=options.feature_visibility_table_file,
        feature_shapes=options.feature_shapes,
        output_prefix="out",
        legend=output.legend if output else "right",
        dinucleotide=options.dinucleotide,
        window=options.window,
        step=options.step,
        depth_window=options.depth_window,
        depth_step=options.depth_step,
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
        linear_track_slots=tracks.linear_track_slots if tracks else None,
        linear_track_axis_index=tracks.linear_track_axis_index if tracks else None,
        annotation_options=options.annotations,
        plot_title=options.plot_title,
        plot_title_position=(
            output.plot_title_position
            if output and output.plot_title_position is not None
            else "bottom"
        ),
        plot_title_font_size=options.plot_title_font_size,
        evalue=options.evalue,
        bitscore=options.bitscore,
        identity=options.identity,
        alignment_length=options.alignment_length,
        _resolved_feature_inputs=_resolved_feature_inputs,
        _return_build_result=_return_build_result,
    )


def build_linear_diagram(
    records: Sequence[SeqRecord],
    *,
    options: LinearDiagramOptions | None = None,
    layout: LinearMultiRecordOptions | None = None,
    losatp_cache: LosatpCacheManager | None = None,
    protein_extraction: ProteinExtractionResult | None = None,
    _resolved_feature_inputs: ResolvedFeatureInputs | None = None,
) -> Drawing:
    """Build a linear diagram using mode-specific typed options."""

    result = _build_linear_diagram(
        records,
        options=options,
        layout=layout,
        losatp_cache=losatp_cache,
        protein_extraction=protein_extraction,
        _resolved_feature_inputs=_resolved_feature_inputs,
    )
    return cast(Drawing, result)


def build_linear_diagram_result(
    records: Sequence[SeqRecord],
    *,
    options: LinearDiagramOptions | None = None,
    layout: LinearMultiRecordOptions | None = None,
    losatp_cache: LosatpCacheManager | None = None,
    protein_extraction: ProteinExtractionResult | None = None,
    _resolved_feature_inputs: ResolvedFeatureInputs | None = None,
) -> LinearDiagramBuildResult:
    """Build a Linear drawing with its computed analysis metadata."""

    result = _build_linear_diagram(
        records,
        options=options,
        layout=layout,
        losatp_cache=losatp_cache,
        protein_extraction=protein_extraction,
        _resolved_feature_inputs=_resolved_feature_inputs,
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
) -> Drawing:
    """Build a circular grid using mode-specific typed options."""

    options = _circular_builder_options(options)
    if layout is not None and not isinstance(layout, CircularMultiRecordOptions):
        raise ValidationError(
            "layout must be CircularMultiRecordOptions or None."
        )
    layout = layout or CircularMultiRecordOptions()
    colors = options.colors
    output = options.output
    tracks = options.tracks
    cfg = _resolve_diagram_options_config(options)

    return assemble_circular_diagram_from_records(
        records,
        cfg=cfg,
        conservation_blast_files=options.conservation_blast_files,
        conservation_dataframes=options.conservation_dataframes,
        conservation_reference=options.conservation_reference,
        conservation_labels=options.conservation_labels,
        conservation_colors=options.conservation_colors,
        conservation_ring_width=options.conservation_ring_width,
        conservation_ring_gap=options.conservation_ring_gap,
        color_table=colors.color_table if colors else None,
        color_table_file=colors.color_table_file if colors else None,
        default_colors=colors.default_colors if colors else None,
        default_colors_palette=colors.default_colors_palette if colors else "default",
        default_colors_file=colors.default_colors_file if colors else None,
        selected_features_set=options.selected_features_set,
        feature_visibility_table=options.feature_visibility_table,
        feature_visibility_table_file=options.feature_visibility_table_file,
        feature_shapes=options.feature_shapes,
        output_prefix="out",
        legend=output.legend if output else "right",
        dinucleotide=options.dinucleotide,
        window=options.window,
        step=options.step,
        depth_window=options.depth_window,
        depth_step=options.depth_step,
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
        annotation_options=options.annotations,
        evalue=options.evalue,
        bitscore=options.bitscore,
        identity=options.identity,
        alignment_length=options.alignment_length,
        _resolved_feature_inputs=_resolved_feature_inputs,
    )


__all__ = [
    "DEFAULT_SELECTED_FEATURES",
    "LinearDiagramMetadata",
    "build_circular_diagram",
    "build_circular_multi_diagram",
    "build_linear_diagram",
]
