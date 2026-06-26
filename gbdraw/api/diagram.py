"""Diagram assembly helpers (public API layer).

This module provides convenience functions to build diagrams from in-memory objects
(e.g., already-parsed SeqRecords) without going through the CLI argument parsing.

The functions here return an `svgwrite.Drawing` (SVG canvas). Saving/conversion is
handled separately (see `gbdraw.render.export.save_figure`).
"""

from __future__ import annotations

import copy
import logging
import math
import re
import xml.etree.ElementTree as ET
from dataclasses import replace
from typing import Optional, Sequence, Mapping, Literal, cast

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from gbdraw.analysis.depth import depth_df as build_depth_df, read_depth_tsv  # type: ignore[reportMissingImports]
from gbdraw.analysis.gc import circular_dinucleotide_content_df  # type: ignore[reportMissingImports]
from gbdraw.analysis.depth_tracks import (  # type: ignore[reportMissingImports]
    DepthTrackData,
    DepthTrackSpec,
    build_depth_track_dataframes,
    depth_track_count,
    depth_track_data_count,
    normalize_depth_tracks,
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
    OrthogroupMembershipMode,
    OrthogroupResult,
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
    CollinearityParameters,
    CollinearityResult,
    CollinearitySearchScope,
    LosslessCollinearityParameters,
    build_orthogroup_collinearity_blocks,
    convert_collinearity_blocks_to_comparisons,
    normalize_collinearity_anchor_mode,
    normalize_collinearity_color_mode,
    normalize_collinearity_search_scope,
)
from gbdraw.analysis.collinearity_units import CollinearityUnitMode  # type: ignore[reportMissingImports]
from gbdraw.analysis.skew import skew_df  # type: ignore[reportMissingImports]
from gbdraw.api.config import apply_config_overrides  # type: ignore[reportMissingImports]
from gbdraw.api.options import DiagramOptions  # type: ignore[reportMissingImports]
from gbdraw.canvas import CircularCanvasConfigurator, LinearCanvasConfigurator  # type: ignore[reportMissingImports]
from gbdraw.canvas.circular import resolve_circular_side_legend_geometry  # type: ignore[reportMissingImports]
from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from gbdraw.config.modify import modify_config_dict  # type: ignore[reportMissingImports]
from gbdraw.config.toml import load_config_toml  # type: ignore[reportMissingImports]
from gbdraw.io.colors import load_default_colors, read_color_table  # type: ignore[reportMissingImports]
from gbdraw.io.record_select import parse_record_selector  # type: ignore[reportMissingImports]
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
from gbdraw.diagrams.circular import assemble_circular_diagram  # type: ignore[reportMissingImports]
from gbdraw.diagrams.circular.positioning import place_definition_group_on_size  # type: ignore[reportMissingImports]
from gbdraw.diagrams.linear import assemble_linear_diagram  # type: ignore[reportMissingImports]
from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]
from gbdraw.features.colors import preprocess_color_tables, precompute_used_color_rules  # type: ignore[reportMissingImports]
from gbdraw.legend.table import (  # type: ignore[reportMissingImports]
    configure_pairwise_identity_legend_from_comparisons,
    prepare_legend_table,
)
from gbdraw.render.groups.circular import DefinitionGroup, LegendGroup  # type: ignore[reportMissingImports]
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
)

DEFAULT_SELECTED_FEATURES = (
    "CDS",
    "rRNA",
    "tRNA",
    "tmRNA",
    "ncRNA",
    "misc_RNA",
    "repeat_region",
)

logger = logging.getLogger(__name__)

_MULTI_RECORD_SUFFIXED_TOP_LEVEL_IDS = {
    "Axis",
    "tick",
    "labels",
    "gc_content",
    "skew",
    "gc_skew",
}
_MULTI_RECORD_SIZE_MODES = {"linear", "auto", "equal", "sqrt"}
_DEFINITION_POSITIONS = {"center", "top", "bottom"}
_LINEAR_PLOT_TITLE_POSITIONS = {"center", "top", "bottom"}
_CIRCULAR_PLOT_TITLE_POSITIONS = {"none", "top", "bottom"}
# Legacy fallback used only when a valid radius cannot be derived.
_MULTI_RECORD_GRID_GAP_PX = 16.0
_MULTI_RECORD_COLUMN_GAP_RATIO = 0.10
_MULTI_RECORD_ROW_GAP_RATIO = 0.05
_MULTI_RECORD_LEGEND_EDGE_PADDING_PX = 20.0
_MULTI_RECORD_LEGEND_GRID_GAP_PX = 20.0
_MULTI_RECORD_LEGEND_TOP_EDGE_PADDING_PX = 32.0
_MULTI_RECORD_LEGEND_PLOT_TITLE_GAP_PX = 20.0
_MULTI_RECORD_PLOT_TITLE_BOTTOM_MARGIN_PX = 24.0
_PLOT_TITLE_DEFAULT_FONT_SIZE = 32.0
_SVG_NUMBER_PATTERN = re.compile(r"[-+]?(?:\d*\.?\d+)(?:[eE][-+]?\d+)?")
_SVG_TRANSLATE_PATTERN = re.compile(
    r"translate\(\s*([-+0-9.eE]+)(?:[\s,]+([-+0-9.eE]+))?\s*\)"
)


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


def _circular_slots_have_renderer(
    slots: Sequence[CircularTrackSlot] | None,
    renderer: str,
) -> bool:
    return any(slot.enabled and str(slot.renderer) == renderer for slot in (slots or []))


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


def _depth_track_heights_from_specs(
    record_depth_tracks: Sequence[Sequence[DepthTrackSpec]] | None,
) -> list[float | None]:
    if not record_depth_tracks:
        return []
    track_count = depth_track_count(record_depth_tracks)
    heights: list[float | None] = [None for _ in range(track_count)]
    for row in record_depth_tracks:
        for track_index, spec in enumerate(row):
            if track_index >= len(heights) or heights[track_index] is not None:
                continue
            if spec.height is not None:
                heights[track_index] = float(spec.height)
    return heights


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


def _linear_slots_have_renderer(
    slots: Sequence[LinearTrackSlot] | None,
    renderer: str,
) -> bool:
    return any(slot.enabled and str(slot.renderer) == renderer for slot in (slots or []))


def _linear_slots_define_renderer(
    slots: Sequence[LinearTrackSlot] | None,
    renderer: str,
) -> bool:
    return any(str(slot.renderer) == renderer for slot in (slots or []))


def _dinucleotides_from_linear_slots(
    slots: Sequence[LinearTrackSlot] | None,
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
    config_dict: dict,
    default_colors: DataFrame,
    cfg: GbdrawConfig,
) -> BlastMatchConfigurator:
    return BlastMatchConfigurator(
        evalue=float(evalue),
        bitscore=float(bitscore),
        identity=float(identity),
        alignment_length=int(alignment_length),
        sequence_length_dict=create_dict_for_sequence_lengths(records),
        config_dict=config_dict,
        default_colors_df=default_colors,
        cfg=cfg,
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
    config_dict: dict,
    default_colors: DataFrame,
    cfg: GbdrawConfig,
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
        config_dict=config_dict,
        default_colors=default_colors,
        cfg=cfg,
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
            spacing=ScalarSpec(float(gap), "px") if gap is not None else None,
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
) -> list[CircularTrackSlot]:
    if not conservation_slots:
        return list(base_slots)
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
                spacing=slot.spacing or (ScalarSpec(float(gap), "px") if gap is not None else None),
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
    if depth_config.tick_interval is not None and float(depth_config.tick_interval) <= 0:
        raise ValidationError("depth_tick_interval/depth_large_tick_interval must be > 0")
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


def _max_depth_from_dataframes(depth_dfs: Sequence[DataFrame | None]) -> float | None:
    max_values = [
        float(depth_df["depth"].max())
        for depth_df in depth_dfs
        if depth_df is not None and not depth_df.empty and "depth" in depth_df.columns
    ]
    return max(max_values) if max_values else None


def _load_depth_table(depth_table: DataFrame | None, depth_file: str | None) -> DataFrame | None:
    if depth_table is not None and depth_file is not None:
        raise ValidationError("Pass either depth_table or depth_file, not both.")
    if depth_table is not None:
        return depth_table
    if depth_file:
        return read_depth_tsv(depth_file)
    return None


def _load_depth_tables(
    *,
    records: Sequence[SeqRecord],
    depth_tables: Sequence[DataFrame] | None,
    depth_files: Sequence[str] | None,
) -> list[DataFrame | None] | None:
    if depth_tables is not None and depth_files is not None:
        raise ValidationError("Pass either depth_tables or depth_files, not both.")
    if depth_tables is None and depth_files is None:
        return None

    record_count = len(records)
    if depth_tables is not None:
        tables = list(depth_tables)
    else:
        depth_paths = list(depth_files or [])
        tables = [read_depth_tsv(path) for path in depth_paths]

    if len(tables) == 1:
        return [tables[0] for _ in range(record_count)]
    if len(tables) != record_count:
        raise ValidationError(
            f"Expected one depth table/file or one per record ({record_count}); got {len(tables)}."
        )
    return list(tables)


def _depth_track_count_from_data(depth_tracks: Sequence[Sequence[DepthTrackData]] | None) -> int:
    return max((len(row) for row in (depth_tracks or ())), default=0)


def _depth_track_count_for_render(
    record_depth_tracks: Sequence[Sequence[DepthTrackSpec]] | None,
    precomputed_depth_tracks: Sequence[DepthTrackData] | None,
    precomputed_depth_df: DataFrame | None = None,
) -> int:
    if precomputed_depth_tracks:
        return len(precomputed_depth_tracks)
    count = depth_track_count(record_depth_tracks)
    if count > 0:
        return count
    return 1 if precomputed_depth_df is not None else 0


def _default_circular_depth_slots_if_needed(
    *,
    parsed_circular_track_slots: list[CircularTrackSlot] | None,
    show_depth: bool,
    depth_track_count_value: int,
    show_gc: bool,
    show_skew: bool,
    dinucleotide: str,
) -> list[CircularTrackSlot] | None:
    if parsed_circular_track_slots is not None:
        return parsed_circular_track_slots
    if not show_depth or depth_track_count_value <= 1:
        return None
    return default_circular_track_slots(
        show_depth=True,
        depth_track_count=depth_track_count_value,
        show_gc=show_gc,
        show_skew=show_skew,
        dinucleotide=dinucleotide,
    )


def _validate_circular_depth_track_indices(
    slots: Sequence[CircularTrackSlot] | None,
    *,
    depth_track_count_value: int,
) -> None:
    if not slots:
        return
    for slot in slots:
        if not slot.enabled or str(slot.renderer) != "depth":
            continue
        raw_index = (slot.params or {}).get("track_index", 0)
        try:
            track_index = int(raw_index or 0)
        except (TypeError, ValueError) as exc:
            raise ValidationError(
                f"Depth slot '{slot.id}' has invalid track_index={raw_index!r}."
            ) from exc
        if track_index < 0 or track_index >= depth_track_count_value:
            raise ValidationError(
                f"Depth slot '{slot.id}' track_index={track_index} is outside the available "
                f"depth track range 0..{max(0, depth_track_count_value - 1)}."
            )


def _validate_linear_depth_track_indices(
    slots: Sequence[LinearTrackSlot] | None,
    *,
    depth_track_count_value: int,
) -> None:
    if not slots:
        return
    for slot in slots:
        if not slot.enabled or str(slot.renderer) != "depth":
            continue
        raw_index = (slot.params or {}).get("track_index", 0)
        try:
            track_index = int(raw_index or 0)
        except (TypeError, ValueError) as exc:
            raise ValidationError(
                f"Depth slot '{slot.id}' has invalid track_index={raw_index!r}."
            ) from exc
        if track_index < 0 or track_index >= depth_track_count_value:
            raise ValidationError(
                f"Depth slot '{slot.id}' track_index={track_index} is outside the available "
                f"depth track range 0..{max(0, depth_track_count_value - 1)}."
            )


def _parse_svg_length_px(value: object, *, default: float = 0.0) -> float:
    """Parse SVG length values like '1000px' into float pixels."""
    if value is None:
        return float(default)
    raw = str(value).strip()
    if raw.endswith("px"):
        raw = raw[:-2]
    try:
        return float(raw)
    except (TypeError, ValueError):
        return float(default)


def _svg_local_name(tag: str) -> str:
    if "}" in tag:
        return tag.rsplit("}", maxsplit=1)[1]
    return tag


def _parse_translate_xy(value: object) -> tuple[float, float]:
    text = str(value or "")
    match = _SVG_TRANSLATE_PATTERN.search(text)
    if match is None:
        return 0.0, 0.0
    x = _parse_svg_length_px(match.group(1), default=0.0)
    y = _parse_svg_length_px(match.group(2), default=0.0) if match.group(2) is not None else 0.0
    return float(x), float(y)


def _parse_svg_numbers(value: object) -> list[float]:
    text = str(value or "")
    values: list[float] = []
    for raw_value in _SVG_NUMBER_PATTERN.findall(text):
        try:
            values.append(float(raw_value))
        except (TypeError, ValueError):
            continue
    return values


def _estimate_element_local_x_bounds(element: ET.Element) -> tuple[float, float] | None:
    """Estimate an element local x-range from geometry attributes."""
    tag = _svg_local_name(str(element.tag)).lower()
    attribs = element.attrib

    if tag == "circle":
        cx = _parse_svg_length_px(attribs.get("cx"), default=0.0)
        radius = abs(_parse_svg_length_px(attribs.get("r"), default=0.0))
        return float(cx - radius), float(cx + radius)

    if tag == "ellipse":
        cx = _parse_svg_length_px(attribs.get("cx"), default=0.0)
        radius_x = abs(_parse_svg_length_px(attribs.get("rx"), default=0.0))
        return float(cx - radius_x), float(cx + radius_x)

    if tag == "line":
        x1 = _parse_svg_length_px(attribs.get("x1"), default=0.0)
        x2 = _parse_svg_length_px(attribs.get("x2"), default=0.0)
        return float(min(x1, x2)), float(max(x1, x2))

    if tag in {"rect", "image"}:
        x = _parse_svg_length_px(attribs.get("x"), default=0.0)
        width = _parse_svg_length_px(attribs.get("width"), default=0.0)
        return float(min(x, x + width)), float(max(x, x + width))

    if tag in {"polygon", "polyline"}:
        points = _parse_svg_numbers(attribs.get("points"))
        if len(points) >= 2:
            xs = points[0::2]
            return float(min(xs)), float(max(xs))

    if tag == "text":
        text_x = _parse_svg_numbers(attribs.get("x"))
        if text_x:
            return float(min(text_x)), float(max(text_x))

    if tag == "path":
        path_numbers = _parse_svg_numbers(attribs.get("d"))
        if path_numbers:
            max_abs = max(abs(value) for value in path_numbers)
            return -float(max_abs), float(max_abs)

    generic_numbers: list[float] = []
    for key in ("x", "x1", "x2", "cx", "points", "d"):
        generic_numbers.extend(_parse_svg_numbers(attribs.get(key)))
    if generic_numbers:
        max_abs = max(abs(value) for value in generic_numbers)
        return -float(max_abs), float(max_abs)
    return None


def _estimate_element_local_y_bounds(element: ET.Element) -> tuple[float, float] | None:
    """Estimate an element local y-range from geometry attributes."""
    tag = _svg_local_name(str(element.tag)).lower()
    attribs = element.attrib

    if tag == "circle":
        cy = _parse_svg_length_px(attribs.get("cy"), default=0.0)
        radius = abs(_parse_svg_length_px(attribs.get("r"), default=0.0))
        return float(cy - radius), float(cy + radius)

    if tag == "ellipse":
        cy = _parse_svg_length_px(attribs.get("cy"), default=0.0)
        radius_y = abs(_parse_svg_length_px(attribs.get("ry"), default=0.0))
        return float(cy - radius_y), float(cy + radius_y)

    if tag == "line":
        y1 = _parse_svg_length_px(attribs.get("y1"), default=0.0)
        y2 = _parse_svg_length_px(attribs.get("y2"), default=0.0)
        return float(min(y1, y2)), float(max(y1, y2))

    if tag in {"rect", "image"}:
        y = _parse_svg_length_px(attribs.get("y"), default=0.0)
        height = _parse_svg_length_px(attribs.get("height"), default=0.0)
        return float(min(y, y + height)), float(max(y, y + height))

    if tag in {"polygon", "polyline"}:
        points = _parse_svg_numbers(attribs.get("points"))
        if len(points) >= 2:
            ys = points[1::2]
            return float(min(ys)), float(max(ys))

    if tag == "text":
        text_y = _parse_svg_numbers(attribs.get("y"))
        if text_y:
            font_size = _parse_svg_length_px(attribs.get("font-size"), default=0.0)
            half_height = 0.5 * font_size if font_size > 0.0 else 0.0
            return float(min(text_y) - half_height), float(max(text_y) + half_height)

    if tag == "path":
        path_numbers = _parse_svg_numbers(attribs.get("d"))
        if path_numbers:
            max_abs = max(abs(value) for value in path_numbers)
            return -float(max_abs), float(max_abs)

    generic_numbers: list[float] = []
    for key in ("y", "y1", "y2", "cy", "points", "d"):
        generic_numbers.extend(_parse_svg_numbers(attribs.get(key)))
    if generic_numbers:
        max_abs = max(abs(value) for value in generic_numbers)
        return -float(max_abs), float(max_abs)
    return None


def _estimate_subcanvas_content_bounds(
    sub_canvas: Drawing,
) -> tuple[tuple[float, float, float, float], tuple[float, float, float, float] | None] | None:
    """Estimate content bounds within a sub-canvas viewBox, excluding defs and legend."""
    try:
        root = ET.fromstring(sub_canvas.tostring())
    except ET.ParseError:
        return None

    view_box = str(root.attrib.get("viewBox", "")).strip()
    view_box_parts = [part for part in view_box.split() if part]
    if len(view_box_parts) == 4:
        vb_x = _parse_svg_length_px(view_box_parts[0], default=0.0)
        vb_y = _parse_svg_length_px(view_box_parts[1], default=0.0)
        vb_w = _parse_svg_length_px(view_box_parts[2], default=0.0)
        vb_h = _parse_svg_length_px(view_box_parts[3], default=0.0)
    else:
        vb_x = 0.0
        vb_y = 0.0
        vb_w = _parse_svg_length_px(root.attrib.get("width"), default=0.0)
        vb_h = _parse_svg_length_px(root.attrib.get("height"), default=0.0)

    if vb_w <= 0.0 and vb_h <= 0.0:
        return None

    min_x: float | None = None
    max_x: float | None = None
    min_y: float | None = None
    max_y: float | None = None

    def _walk(node: ET.Element, inherited_tx: float, inherited_ty: float) -> None:
        nonlocal min_x, max_x, min_y, max_y
        tag = _svg_local_name(str(node.tag)).lower()
        if tag == "defs":
            return
        if tag == "g" and str(node.attrib.get("id", "")) == "legend":
            return

        local_tx, local_ty = _parse_translate_xy(node.attrib.get("transform"))
        total_tx = inherited_tx + float(local_tx)
        total_ty = inherited_ty + float(local_ty)

        bounds_x = _estimate_element_local_x_bounds(node)
        if bounds_x is not None:
            candidate_min_x = total_tx + float(bounds_x[0])
            candidate_max_x = total_tx + float(bounds_x[1])
            min_x = candidate_min_x if min_x is None else min(min_x, candidate_min_x)
            max_x = candidate_max_x if max_x is None else max(max_x, candidate_max_x)

        bounds_y = _estimate_element_local_y_bounds(node)
        if bounds_y is not None:
            candidate_min_y = total_ty + float(bounds_y[0])
            candidate_max_y = total_ty + float(bounds_y[1])
            min_y = candidate_min_y if min_y is None else min(min_y, candidate_min_y)
            max_y = candidate_max_y if max_y is None else max(max_y, candidate_max_y)

        for child in list(node):
            _walk(child, total_tx, total_ty)

    _walk(root, 0.0, 0.0)

    content_bounds: tuple[float, float, float, float] | None = None
    if min_x is not None and max_x is not None and min_y is not None and max_y is not None:
        content_bounds = (float(min_x), float(max_x), float(min_y), float(max_y))

    return (
        (float(vb_x), float(vb_x) + float(vb_w), float(vb_y), float(vb_y) + float(vb_h)),
        content_bounds,
    )


def _estimate_subcanvas_horizontal_insets(sub_canvas: Drawing) -> tuple[float, float]:
    """Estimate left/right insets between viewBox edges and diagram geometry."""
    parsed = _estimate_subcanvas_content_bounds(sub_canvas)
    if parsed is None:
        return 0.0, 0.0
    (view_min_x, view_max_x, _view_min_y, _view_max_y), content_bounds = parsed
    if content_bounds is None:
        return 0.0, 0.0
    content_min_x, content_max_x, _content_min_y, _content_max_y = content_bounds
    left_inset = max(0.0, float(content_min_x) - float(view_min_x))
    right_inset = max(0.0, float(view_max_x) - float(content_max_x))
    return float(left_inset), float(right_inset)


def _estimate_subcanvas_vertical_insets(sub_canvas: Drawing) -> tuple[float, float]:
    """Estimate top/bottom insets between viewBox edges and diagram geometry."""
    parsed = _estimate_subcanvas_content_bounds(sub_canvas)
    if parsed is None:
        return 0.0, 0.0
    (_view_min_x, _view_max_x, view_min_y, view_max_y), content_bounds = parsed
    if content_bounds is None:
        return 0.0, 0.0
    _content_min_x, _content_max_x, content_min_y, content_max_y = content_bounds
    top_inset = max(0.0, float(content_min_y) - float(view_min_y))
    bottom_inset = max(0.0, float(view_max_y) - float(content_max_y))
    return float(top_inset), float(bottom_inset)


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


def _resolve_multi_record_selector_index(
    records: Sequence[SeqRecord],
    selector_text: str,
) -> int:
    """Resolve one record selector to an index in records."""
    record_count = len(records)
    try:
        selector = parse_record_selector(selector_text)
    except ValueError as exc:
        raise ValidationError(str(exc)) from exc
    if selector is None:
        raise ValidationError(f"multi_record_position selector '{selector_text}' is invalid.")

    if selector.record_index is not None:
        idx = int(selector.record_index)
        if idx < 0 or idx >= record_count:
            raise ValidationError(
                f"multi_record_position selector '{selector_text}' is out of range for "
                f"{record_count} record(s)."
            )
        return idx

    target_record_id = str(selector.record_id or "")
    matches = [idx for idx, record in enumerate(records) if str(record.id) == target_record_id]
    if not matches:
        raise ValidationError(
            f"multi_record_position selector '{selector_text}' did not match any record ID."
        )
    if len(matches) > 1:
        raise ValidationError(
            f"multi_record_position selector '{selector_text}' matched multiple records. "
            "Use #index to disambiguate."
        )
    return int(matches[0])


def _parse_multi_record_position(value: str) -> tuple[str, int]:
    """Parse one <selector>@<row> token."""
    raw = str(value or "").strip()
    if not raw:
        raise ValidationError("multi_record_position does not allow empty entries.")
    if "@" not in raw:
        raise ValidationError(
            f"multi_record_position entry '{raw}' must be in '<selector>@<row>' format."
        )
    selector_text, row_text = raw.rsplit("@", 1)
    selector_text = selector_text.strip()
    row_text = row_text.strip()
    if not selector_text:
        raise ValidationError(
            f"multi_record_position entry '{raw}' must include a selector before '@'."
        )
    if not row_text or not row_text.isdigit() or int(row_text) <= 0:
        raise ValidationError(
            f"multi_record_position entry '{raw}' must use a positive integer row."
        )
    return selector_text, int(row_text)


def _resolve_multi_record_positions(
    records: Sequence[SeqRecord],
    positions: Sequence[str] | None,
) -> tuple[list[int], list[int]]:
    """Resolve explicit multi-record positions to ordered indices and row counts."""
    record_count = len(records)
    if record_count <= 0:
        return [], []
    if not positions:
        return list(range(record_count)), _resolve_multi_record_default_row_counts(record_count)

    seen_indices: set[int] = set()
    row_entries: dict[int, list[int]] = {}
    provided_entries = 0
    for raw_position in positions:
        selector_text, row_value = _parse_multi_record_position(str(raw_position))
        resolved_index = _resolve_multi_record_selector_index(records, selector_text)
        if resolved_index in seen_indices:
            raise ValidationError(
                f"multi_record_position selector '{selector_text}' was specified more than once."
            )
        seen_indices.add(resolved_index)
        row_entries.setdefault(int(row_value), []).append(resolved_index)
        provided_entries += 1

    if len(seen_indices) != record_count:
        raise ValidationError(
            f"multi_record_position must include each loaded record exactly once "
            f"(expected {record_count}, got {len(seen_indices)} unique selector(s))."
        )
    if provided_entries != record_count:
        raise ValidationError(
            f"multi_record_position must provide exactly {record_count} entry(ies)."
        )

    ordered_indices: list[int] = []
    row_counts: list[int] = []
    for _row_value in sorted(row_entries):
        indices = row_entries[_row_value]
        if not indices:
            continue
        ordered_indices.extend(indices)
        row_counts.append(len(indices))

    if len(ordered_indices) != record_count:
        raise ValidationError(
            "multi_record_position internal error: failed to resolve all records."
        )
    return ordered_indices, row_counts


def _group_local_vertical_bounds(group: Group) -> tuple[float, float]:
    """Return local vertical bounds for text elements in a group."""
    min_y: float | None = None
    max_y: float | None = None

    for element in getattr(group, "elements", []):
        attribs = getattr(element, "attribs", None)
        if not isinstance(attribs, dict):
            continue
        if "y" not in attribs:
            continue
        y_value = _parse_svg_length_px(attribs.get("y"), default=0.0)
        font_size = _parse_svg_length_px(attribs.get("font-size"), default=0.0)
        half_height = 0.5 * font_size if font_size > 0 else 0.0
        top = y_value - half_height
        bottom = y_value + half_height
        min_y = top if min_y is None else min(min_y, top)
        max_y = bottom if max_y is None else max(max_y, bottom)

    if min_y is None or max_y is None:
        return 0.0, 0.0
    return float(min_y), float(max_y)


def _set_group_translate(group: Group, *, x: float, y: float) -> None:
    """Set a group transform to a simple translate(x, y)."""
    attribs = getattr(group, "attribs", None)
    if not isinstance(attribs, dict):
        return
    attribs["transform"] = f"translate({float(x)}, {float(y)})"


def _center_record_definition_group_on_record_axis(
    record_group: Group,
    *,
    record_index: int,
    record_id: str,
) -> None:
    """Center one per-record definition group vertically on its record axis."""
    axis_group_id = f"Axis_{record_index}"
    definition_group_id = f"{str(record_id).replace(' ', '_')}_definition"
    axis_group: Group | None = None
    definition_group: Group | None = None

    for child in getattr(record_group, "elements", []):
        attribs = getattr(child, "attribs", None)
        if not isinstance(attribs, dict):
            continue
        child_id = str(attribs.get("id", ""))
        if child_id == axis_group_id:
            axis_group = child
        elif child_id == definition_group_id:
            definition_group = child

    if axis_group is None or definition_group is None:
        return

    _axis_x, axis_center_y = _parse_translate_xy(axis_group.attribs.get("transform"))
    definition_x, definition_y = _parse_translate_xy(definition_group.attribs.get("transform"))
    definition_min_y, definition_max_y = _group_local_vertical_bounds(definition_group)
    definition_local_center_y = (float(definition_min_y) + float(definition_max_y)) * 0.5
    definition_text_center_y = float(definition_y) + float(definition_local_center_y)
    delta_y = float(axis_center_y) - float(definition_text_center_y)

    if math.isclose(delta_y, 0.0, rel_tol=1e-9, abs_tol=1e-9):
        return
    _set_group_translate(
        definition_group,
        x=float(definition_x),
        y=float(definition_y) + float(delta_y),
    )


def _center_single_record_definition_group_on_axis(
    canvas: Drawing,
    *,
    record_id: str,
) -> None:
    """Center a top-level single-record definition group vertically on its axis."""
    axis_group: Group | None = None
    definition_group: Group | None = None
    definition_group_id = f"{str(record_id).replace(' ', '_')}_definition"

    for child in getattr(canvas, "elements", []):
        attribs = getattr(child, "attribs", None)
        if not isinstance(attribs, dict):
            continue
        child_id = str(attribs.get("id", ""))
        if child_id == "Axis":
            axis_group = child
        elif child_id == definition_group_id:
            definition_group = child

    if axis_group is None or definition_group is None:
        return

    _axis_x, axis_center_y = _parse_translate_xy(axis_group.attribs.get("transform"))
    definition_x, definition_y = _parse_translate_xy(definition_group.attribs.get("transform"))
    definition_min_y, definition_max_y = _group_local_vertical_bounds(definition_group)
    definition_local_center_y = (float(definition_min_y) + float(definition_max_y)) * 0.5
    definition_text_center_y = float(definition_y) + float(definition_local_center_y)
    delta_y = float(axis_center_y) - float(definition_text_center_y)

    if math.isclose(delta_y, 0.0, rel_tol=1e-9, abs_tol=1e-9):
        return
    _set_group_translate(
        definition_group,
        x=float(definition_x),
        y=float(definition_y) + float(delta_y),
    )


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
            "multi_record_size_mode must be one of: auto, linear, equal, sqrt (alias of auto)"
        )
    if normalized == "sqrt":
        normalized = "auto"
    return cast(Literal["linear", "auto", "equal"], normalized)


def _resolve_definition_position(
    position: str,
    *,
    argument_name: str = "definition_position",
) -> Literal["center", "top", "bottom"]:
    normalized = str(position).strip().lower()
    if normalized not in _DEFINITION_POSITIONS:
        raise ValidationError(
            f"{argument_name} must be one of: center, top, bottom"
    )
    return cast(Literal["center", "top", "bottom"], normalized)


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


def _resolve_pairwise_match_style(style: str) -> Literal["ribbon", "curve"]:
    normalized = str(style).strip().lower()
    if normalized not in {"ribbon", "curve"}:
        raise ValidationError("pairwise_match_style must be one of: ribbon, curve")
    return cast(Literal["ribbon", "curve"], normalized)


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
    config_dict: dict,
    cfg: GbdrawConfig,
    plot_title_font_size: float | None,
) -> tuple[dict, GbdrawConfig]:
    if plot_title_font_size is None:
        return config_dict, cfg
    resolved_font_size = _resolve_plot_title_font_size(plot_title_font_size)
    updated_config_dict = modify_config_dict(
        config_dict,
        plot_title_font_size=resolved_font_size,
    )
    updated_cfg = apply_config_overrides(
        cfg,
        {"plot_title_font_size": resolved_font_size},
    )
    return updated_config_dict, updated_cfg


def _sync_drawing_canvas_size(
    canvas: Drawing,
    *,
    width: float,
    height: float,
) -> None:
    canvas.attribs["width"] = f"{float(width)}px"
    canvas.attribs["height"] = f"{float(height)}px"
    canvas.attribs["viewBox"] = f"0 0 {float(width)} {float(height)}"


def _translate_canvas_top_level_groups(canvas: Drawing, *, dy: float) -> None:
    if abs(float(dy)) <= 1e-6:
        return
    for element in getattr(canvas, "elements", []):
        if _is_defs_element(element):
            continue
        translate = getattr(element, "translate", None)
        if callable(translate):
            translate(0, float(dy))


def _build_circular_plot_title_group(
    *,
    gb_record: SeqRecord,
    output_prefix: str,
    config_dict: dict,
    cfg: GbdrawConfig,
    species: str | None,
    strain: str | None,
    plot_title: str | None,
) -> tuple[Group, tuple[float, float]]:
    plot_title_canvas_config = CircularCanvasConfigurator(
        output_prefix=output_prefix,
        config_dict=config_dict,
        legend="none",
        gb_record=gb_record,
        cfg=cfg,
    )
    plot_title_group = DefinitionGroup(
        gb_record=gb_record,
        canvas_config=plot_title_canvas_config,
        config_dict=config_dict,
        species=species,
        strain=strain,
        plot_title=plot_title,
        definition_profile="shared_common",
        definition_group_id="plot_title",
        cfg=cfg,
    ).get_group()
    return plot_title_group, _group_local_vertical_bounds(plot_title_group)


def _add_single_record_plot_title_group(
    *,
    canvas: Drawing,
    canvas_config: CircularCanvasConfigurator,
    legend_config: LegendDrawingConfigurator,
    gb_record: SeqRecord,
    output_prefix: str,
    config_dict: dict,
    cfg: GbdrawConfig,
    species: str | None,
    strain: str | None,
    plot_title: str | None,
    plot_title_position: Literal["none", "top", "bottom"],
) -> None:
    plot_title_group, plot_title_local_bounds = _build_circular_plot_title_group(
        gb_record=gb_record,
        output_prefix=output_prefix,
        config_dict=config_dict,
        cfg=cfg,
        species=species,
        strain=strain,
        plot_title=plot_title,
    )
    plot_title_min_y, plot_title_max_y = plot_title_local_bounds
    plot_title_height = max(
        0.0,
        float(plot_title_max_y) - float(plot_title_min_y),
    )

    if (
        plot_title_position == "top"
        and str(canvas_config.legend_position).strip().lower() == "top"
        and float(legend_config.legend_height) > 0.0
    ):
        legend_top = (
            float(canvas_config.legend_offset_y)
            - (0.5 * float(legend_config.color_rect_size))
        )
        required_legend_top = (
            _MULTI_RECORD_PLOT_TITLE_BOTTOM_MARGIN_PX
            + plot_title_height
            + _MULTI_RECORD_LEGEND_PLOT_TITLE_GAP_PX
        )
        if legend_top < required_legend_top:
            shift_y = required_legend_top - legend_top
            _translate_canvas_top_level_groups(canvas, dy=shift_y)
            canvas_config.offset_y = float(canvas_config.offset_y) + float(shift_y)
            canvas_config.total_height = float(canvas_config.total_height) + float(shift_y)
            _sync_drawing_canvas_size(
                canvas,
                width=float(canvas_config.total_width),
                height=float(canvas_config.total_height),
            )

    if plot_title_position == "bottom":
        anchor_bottom = 0.0
        if str(canvas_config.legend_position).strip().lower() == "bottom":
            anchor_bottom = (
                float(canvas_config.legend_offset_y)
                - (0.5 * float(legend_config.color_rect_size))
                + float(legend_config.legend_height)
            )
        required_height = (
            anchor_bottom
            + _MULTI_RECORD_LEGEND_PLOT_TITLE_GAP_PX
            + plot_title_height
            + _MULTI_RECORD_PLOT_TITLE_BOTTOM_MARGIN_PX
        )
        if required_height > float(canvas_config.total_height):
            canvas_config.total_height = float(required_height)
            _sync_drawing_canvas_size(
                canvas,
                width=float(canvas_config.total_width),
                height=float(canvas_config.total_height),
            )

    plot_title_group = place_definition_group_on_size(
        plot_title_group,
        canvas_width=float(canvas_config.total_width),
        canvas_height=float(canvas_config.total_height),
        position=plot_title_position,
    )
    canvas.add(plot_title_group)


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


def assemble_linear_diagram_from_records(
    records: Sequence[SeqRecord],
    *,
    blast_files: Optional[Sequence[str]] = None,
    protein_comparisons: Sequence[DataFrame] | None = None,
    orthogroups: OrthogroupResult | None = None,
    protein_blastp_mode: ProteinBlastpMode | str = "none",
    pairwise_match_style: Literal["ribbon", "curve"] | str = "ribbon",
    collinearity_blocks: CollinearityResult | Sequence[CollinearityBlock] | None = None,
    collinearity_params: CollinearityParameters | LosslessCollinearityParameters | None = None,
    collinearity_unit_mode: CollinearityUnitMode | str = "auto",
    collinearity_anchor_mode: CollinearityAnchorMode | str = "rbh",
    collinearity_search_scope: CollinearitySearchScope | str = "adjacent",
    collinearity_color_mode: CollinearityColorMode | str = "orientation",
    losatp_bin: str = "losat",
    ncbi_blastp_bin: str | None = None,
    losatp_threads: int | None = None,
    protein_blastp_max_hits: int = 5,
    protein_blastp_candidate_limit: int | None = None,
    orthogroup_membership_mode: OrthogroupMembershipMode | str = "anchor_core_v1",
    orthogroup_member_max_hits: int = 5,
    collinear_max_paralog_links_per_orthogroup: int = 2,
    align_orthogroup_feature: str | None = None,
    config_dict: dict | None = None,
    config_overrides: Mapping[str, object] | None = None,
    color_table: Optional[DataFrame] = None,
    color_table_file: str | None = None,
    default_colors: DataFrame | None = None,
    default_colors_palette: str = "default",
    default_colors_file: str | None = None,
    selected_features_set: Sequence[str] | None = None,
    feature_table: DataFrame | None = None,
    feature_table_file: str | None = None,
    feature_shapes: Mapping[str, str] | None = None,
    output_prefix: str = "out",
    legend: str = "right",
    dinucleotide: str = "GC",
    window: Optional[int] = None,
    step: Optional[int] = None,
    depth_window: Optional[int] = None,
    depth_step: Optional[int] = None,
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
    plot_title: str | None = None,
    plot_title_position: Literal["center", "top", "bottom"] = "bottom",
    plot_title_font_size: float | None = None,
    evalue: float = 1e-5,
    bitscore: float = 50.0,
    identity: float = 70.0,
    alignment_length: int = 0,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
    """Builds and assembles a linear diagram for the given records.

    This is a convenience wrapper that builds internal configurators/canvas objects and
    returns the assembled SVG canvas.
    If config_dict is None, it loads gbdraw.data/config.toml.
    If config_overrides is provided, modify_config_dict is applied.
    If default_colors is None, it loads the built-in default palette.
    If color_table is None and color_table_file is provided, it is loaded.
    If selected_features_set is None, it uses the CLI default feature list.
    """
    if not records:
        raise ValidationError("records is empty")
    if alignment_length < 0:
        raise ValidationError("alignment_length must be >= 0")
    normalized_protein_blastp_mode = normalize_protein_blastp_mode(protein_blastp_mode)
    normalized_pairwise_match_style = _resolve_pairwise_match_style(pairwise_match_style)
    normalize_collinearity_anchor_mode(str(collinearity_anchor_mode))
    normalized_collinearity_anchor_mode = "rbh"
    normalized_collinearity_search_scope = normalize_collinearity_search_scope(str(collinearity_search_scope))
    normalized_collinearity_color_mode = normalize_collinearity_color_mode(str(collinearity_color_mode))
    normalized_orthogroup_membership_mode = normalize_orthogroup_membership_mode(str(orthogroup_membership_mode))
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
    has_precomputed_comparisons = bool(blast_files or protein_comparisons is not None or collinearity_blocks is not None)
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
    if color_table is None and color_table_file is not None:
        color_table = read_color_table(color_table_file)
    if feature_table is None and feature_table_file is not None:
        feature_table = read_feature_visibility_file(feature_table_file)

    if default_colors is None:
        has_comparisons = bool(
            blast_files
            or protein_comparisons
            or collinearity_blocks
            or normalized_protein_blastp_mode != "none"
        )
        default_colors = load_default_colors(
            user_defined_default_colors=default_colors_file or "",
            palette=default_colors_palette or "default",
            load_comparison=has_comparisons,
        )

    if config_dict is None:
        config_dict = load_config_toml("gbdraw.data", "config.toml")
    if config_overrides:
        if cfg is not None:
            raise ValueError(
                "config_overrides cannot be used with cfg; pass cfg=None or apply overrides before."
            )
        config_dict = modify_config_dict(config_dict, **config_overrides)
    pairwise_style_from_overrides = (
        config_overrides is not None
        and normalized_pairwise_match_style == "ribbon"
        and "pairwise_match_style" in config_overrides
    )
    if pairwise_style_from_overrides:
        normalized_pairwise_match_style = _resolve_pairwise_match_style(
            str(config_overrides["pairwise_match_style"])
        )
    else:
        config_dict = modify_config_dict(config_dict, pairwise_match_style=normalized_pairwise_match_style)
    if cfg is not None:
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
    else:
        cfg = GbdrawConfig.from_dict(config_dict)
    record_depth_tracks = normalize_depth_tracks(
        records,
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
        show_depth = _linear_slots_have_renderer(parsed_linear_track_slots, "depth")
        show_gc = _linear_slots_have_renderer(parsed_linear_track_slots, "dinucleotide_content")
        show_skew = _linear_slots_have_renderer(parsed_linear_track_slots, "dinucleotide_skew")
        if show_depth and record_depth_tracks is None:
            raise ValidationError("A linear depth track slot requires a depth_table, depth_file, or depth_track input.")
        if show_depth:
            _validate_linear_depth_track_indices(
                parsed_linear_track_slots,
                depth_track_count_value=max(1, available_depth_track_count),
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

    if selected_features_set is None:
        selected_features_set = DEFAULT_SELECTED_FEATURES
    resolved_protein_comparisons: list[DataFrame] | None = None
    resolved_orthogroups: OrthogroupResult | None = orthogroups
    if protein_comparisons is not None:
        resolved_protein_comparisons = list(protein_comparisons)
    elif collinearity_blocks is not None:
        if isinstance(collinearity_blocks, CollinearityResult):
            collinearity_result = collinearity_blocks
        else:
            collinearity_result = CollinearityResult(blocks=tuple(collinearity_blocks))
        if collinearity_result.orthogroups is not None:
            resolved_orthogroups = collinearity_result.orthogroups
        resolved_protein_comparisons = convert_collinearity_blocks_to_comparisons(
            collinearity_result,
            records=records,
            color_mode=normalized_collinearity_color_mode,
        )
    elif normalized_protein_blastp_mode == "pairwise":
        protein_blastp_result = build_pairwise_protein_blastp_comparisons(
            records,
            losatp_bin=losatp_bin,
            ncbi_blastp_bin=ncbi_blastp_bin,
            losatp_threads=losatp_threads,
            max_hits=int(protein_blastp_max_hits),
            candidate_limit=protein_blastp_candidate_limit,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
        )
        resolved_protein_comparisons = protein_blastp_result.comparisons
    elif normalized_protein_blastp_mode == "orthogroup":
        protein_blastp_result = build_rbh_orthogroup_protein_blastp_comparisons(
            records,
            losatp_bin=losatp_bin,
            ncbi_blastp_bin=ncbi_blastp_bin,
            losatp_threads=losatp_threads,
            candidate_limit=protein_blastp_candidate_limit,
            orthogroup_membership_mode=normalized_orthogroup_membership_mode,
            orthogroup_member_max_hits=int(orthogroup_member_max_hits),
            max_related_edges_per_orthogroup=int(collinear_max_paralog_links_per_orthogroup),
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
        )
        resolved_protein_comparisons = protein_blastp_result.comparisons
        resolved_orthogroups = protein_blastp_result.orthogroups
    elif normalized_protein_blastp_mode == "collinear":
        collinearity_result = build_orthogroup_collinearity_blocks(
            records,
            losatp_bin=losatp_bin,
            ncbi_blastp_bin=ncbi_blastp_bin,
            losatp_threads=losatp_threads,
            candidate_limit=protein_blastp_candidate_limit,
            orthogroup_membership_mode=normalized_orthogroup_membership_mode,
            orthogroup_member_max_hits=int(orthogroup_member_max_hits),
            max_paralog_links_per_orthogroup=int(collinear_max_paralog_links_per_orthogroup),
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
            params=collinearity_params,
            unit_mode=collinearity_unit_mode,
            edge_mode=normalized_collinearity_anchor_mode,
            search_scope=normalized_collinearity_search_scope,
        )
        resolved_orthogroups = collinearity_result.orthogroups
        resolved_protein_comparisons = convert_collinearity_blocks_to_comparisons(
            collinearity_result,
            records=records,
            color_mode=normalized_collinearity_color_mode,
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
        config_dict=config_dict,
        default_colors_df=default_colors,
        cfg=cfg,
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
        config_dict=config_dict,
        legend=legend,
        output_prefix=output_prefix,
        cfg=cfg,
        has_comparisons=bool(blast_files or resolved_protein_comparisons),
        depth_track_count=max(1, depth_track_count(record_depth_tracks)),
        depth_track_heights=_depth_track_heights_from_specs(record_depth_tracks),
    )
    feature_config = FeatureDrawingConfigurator(
        color_table=color_table,
        default_colors=default_colors,
        selected_features_set=list(selected_features_set),
        feature_table=feature_table,
        feature_shapes=feature_shapes,
        config_dict=config_dict,
        canvas_config=canvas_config,
        cfg=cfg,
    )
    gc_config = GcContentConfigurator(
        window=window,
        step=step,
        dinucleotide=dinucleotide,
        config_dict=config_dict,
        default_colors_df=default_colors,
        cfg=cfg,
    )
    skew_config = GcSkewConfigurator(
        window=window,
        step=step,
        dinucleotide=dinucleotide,
        config_dict=config_dict,
        default_colors_df=default_colors,
        cfg=cfg,
    )
    depth_config = DepthConfigurator(
        window=resolved_depth_window,
        step=resolved_depth_step,
        config_dict=config_dict,
        cfg=cfg,
    ) if show_depth else None
    legend_config = LegendDrawingConfigurator(
        color_table=color_table,
        default_colors=default_colors,
        selected_features_set=list(selected_features_set),
        config_dict=config_dict,
        gc_config=gc_config,
        skew_config=skew_config,
        feature_config=feature_config,
        blast_config=blast_config,
        canvas_config=canvas_config,
        cfg=cfg,
    )

    canvas = assemble_linear_diagram(
        records=list(records),
        blast_files=list(blast_files) if blast_files else None,
        canvas_config=canvas_config,
        blast_config=blast_config,
        feature_config=feature_config,
        gc_config=gc_config,
        config_dict=config_dict,
        legend_config=legend_config,
        skew_config=skew_config,
        depth_config=depth_config,
        record_depth_tracks=record_depth_tracks,
        linear_track_slots=parsed_linear_track_slots,
        linear_track_axis_index=resolved_linear_track_axis_index,
        plot_title=normalized_plot_title or None,
        plot_title_position=normalized_plot_title_position,
        plot_title_font_size=resolved_plot_title_font_size,
        comparison_dataframes=resolved_protein_comparisons,
        orthogroups=resolved_orthogroups,
        align_orthogroup_feature=align_orthogroup_feature,
        cfg=cfg,
    )
    try:
        setattr(canvas, "_gbdraw_orthogroups", resolved_orthogroups)
    except Exception:
        pass
    return canvas


def assemble_circular_diagram_from_record(
    gb_record: SeqRecord,
    *,
    conservation_blast_files: Sequence[str] | None = None,
    conservation_dataframes: Sequence[DataFrame] | None = None,
    conservation_reference: Literal["query", "subject", "auto"] | str = "auto",
    conservation_labels: Sequence[str] | None = None,
    conservation_colors: Sequence[str] | None = None,
    conservation_ring_width: float | None = None,
    conservation_ring_gap: float | None = None,
    config_dict: dict | None = None,
    config_overrides: Mapping[str, object] | None = None,
    color_table: Optional[DataFrame] = None,
    color_table_file: str | None = None,
    default_colors: DataFrame | None = None,
    default_colors_palette: str = "default",
    default_colors_file: str | None = None,
    selected_features_set: Sequence[str] | None = None,
    feature_table: DataFrame | None = None,
    feature_table_file: str | None = None,
    feature_shapes: Mapping[str, str] | None = None,
    output_prefix: str = "out",
    legend: str = "right",
    dinucleotide: str = "GC",
    window: Optional[int] = None,
    step: Optional[int] = None,
    depth_window: Optional[int] = None,
    depth_step: Optional[int] = None,
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
    evalue: float = 1e-5,
    bitscore: float = 50.0,
    identity: float = 70.0,
    alignment_length: int = 0,
    _definition_profile: Literal["full", "record_summary", "shared_common"] = "full",
    _tick_track_channel_override: Literal["short", "long"] | None = None,
    _precomputed_depth_df: DataFrame | None = None,
    _precomputed_depth_tracks: Sequence[DepthTrackData] | None = None,
    _shared_depth_max: float | None = None,
    _precomputed_conservation_tracks: Sequence[ConservationTrack] | None = None,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
    """Builds and assembles a circular diagram for a single record.

    If config_dict is None, it loads gbdraw.data/config.toml.
    If config_overrides is provided, modify_config_dict is applied.
    If default_colors is None, it loads the built-in default palette.
    If color_table is None and color_table_file is provided, it is loaded.
    If selected_features_set is None, it uses the CLI default feature list.
    """
    _validate_positive_optional("depth_window", depth_window)
    _validate_positive_optional("depth_step", depth_step)
    _validate_positive_float_optional("conservation_ring_width", conservation_ring_width)
    _validate_positive_float_optional("conservation_ring_gap", conservation_ring_gap)
    _validate_nonnegative_float_optional("center_reserved_radius", center_reserved_radius)
    if color_table is None and color_table_file is not None:
        color_table = read_color_table(color_table_file)
    if feature_table is None and feature_table_file is not None:
        feature_table = read_feature_visibility_file(feature_table_file)

    if default_colors is None:
        default_colors = load_default_colors(
            user_defined_default_colors=default_colors_file or "",
            palette=default_colors_palette or "default",
            load_comparison=False,
        )

    if config_dict is None:
        config_dict = load_config_toml("gbdraw.data", "config.toml")
    if config_overrides:
        if cfg is not None:
            raise ValueError(
                "config_overrides cannot be used with cfg; pass cfg=None or apply overrides before."
            )
        config_dict = modify_config_dict(config_dict, **config_overrides)
    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    config_dict, cfg = _apply_circular_plot_title_font_size_override(
        config_dict=config_dict,
        cfg=cfg,
        plot_title_font_size=plot_title_font_size,
    )
    record_depth_tracks = normalize_depth_tracks(
        [gb_record],
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
    precomputed_depth_track_list = list(_precomputed_depth_tracks or [])
    if cfg.canvas.show_depth and record_depth_tracks is None and not precomputed_depth_track_list and _precomputed_depth_df is None:
        raise ValidationError("show_depth requires a depth_table, depth_file, or depth_track input.")
    show_depth_from_input = (
        record_depth_tracks is not None
        or bool(precomputed_depth_track_list)
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
            config_dict=config_dict,
            default_colors=default_colors,
            cfg=cfg,
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

    parsed_circular_track_slots = _parse_circular_track_slot_inputs(circular_track_slots)
    circular_track_axis_index = _validate_circular_track_axis_index(
        circular_track_axis_index,
        parsed_circular_track_slots,
    )

    legend_effective = legend

    # Explicit slots override high-level show flags used by canvas sizing.
    if parsed_circular_track_slots is not None:
        show_depth = _circular_slots_have_renderer(parsed_circular_track_slots, "depth")
        show_gc = _circular_slots_have_renderer(parsed_circular_track_slots, "dinucleotide_content")
        show_skew = _circular_slots_have_renderer(parsed_circular_track_slots, "dinucleotide_skew")
        if show_depth and not show_depth_from_input:
            raise ValidationError("A circular depth track slot requires a depth_table, depth_file, or depth_track input.")
    else:
        show_depth = cfg.canvas.show_depth
        show_gc = cfg.canvas.show_gc
        show_skew = cfg.canvas.show_skew
    available_depth_track_count = _depth_track_count_for_render(
        record_depth_tracks,
        precomputed_depth_track_list,
        _precomputed_depth_df,
    )
    parsed_circular_track_slots = _default_circular_depth_slots_if_needed(
        parsed_circular_track_slots=parsed_circular_track_slots,
        show_depth=bool(show_depth and show_depth_from_input),
        depth_track_count_value=available_depth_track_count,
        show_gc=bool(show_gc),
        show_skew=bool(show_skew),
        dinucleotide=dinucleotide,
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
                    show_depth=show_depth,
                    depth_track_count=available_depth_track_count,
                    show_gc=show_gc,
                    show_skew=show_skew,
                    dinucleotide=dinucleotide,
                )
            parsed_circular_track_slots = _insert_conservation_slots(
                parsed_circular_track_slots,
                conservation_slots,
            )
    canvas_cfg = cfg.canvas
    canvas_cfg = replace(
        canvas_cfg,
        show_depth=bool(show_depth and show_depth_from_input),
        show_gc=bool(show_gc),
        show_skew=bool(show_skew),
    )
    cfg = replace(cfg, canvas=canvas_cfg)

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
        config_dict=config_dict,
        default_colors_df=default_colors,
        cfg=cfg,
    )
    skew_config = GcSkewConfigurator(
        window=window,
        step=step,
        dinucleotide=dinucleotide,
        config_dict=config_dict,
        default_colors_df=default_colors,
        cfg=cfg,
    )
    depth_config = (
        DepthConfigurator(
            window=resolved_depth_window,
            step=resolved_depth_step,
            config_dict=config_dict,
            cfg=cfg,
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
            resolved_depth_df = resolved_depth_tracks[0].df if resolved_depth_tracks else None
        elif _precomputed_depth_df is not None:
            resolved_depth_df = _precomputed_depth_df
            resolved_depth_tracks = [
                DepthTrackData(
                    id="depth",
                    label="Depth",
                    df=resolved_depth_df,
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
            resolved_depth_df = resolved_depth_tracks[0].df if resolved_depth_tracks else None
        else:
            resolved_depth_df = None
            resolved_depth_tracks = []
    else:
        resolved_depth_df = None
        resolved_depth_tracks = []
    _validate_circular_depth_track_indices(
        parsed_circular_track_slots,
        depth_track_count_value=max(1, len(resolved_depth_tracks)),
    )

    canvas_config = CircularCanvasConfigurator(
        output_prefix=output_prefix, config_dict=config_dict, legend=legend_effective, gb_record=gb_record, cfg=cfg
    )
    feature_config = FeatureDrawingConfigurator(
        color_table=color_table,
        default_colors=default_colors,
        selected_features_set=list(selected_features_set),
        feature_table=feature_table,
        feature_shapes=feature_shapes,
        config_dict=config_dict,
        canvas_config=canvas_config,
        cfg=cfg,
    )
    legend_config = LegendDrawingConfigurator(
        color_table=color_table,
        default_colors=default_colors,
        selected_features_set=list(selected_features_set),
        config_dict=config_dict,
        gc_config=gc_config,
        skew_config=skew_config,
        feature_config=feature_config,
        canvas_config=canvas_config,
        cfg=cfg,
    )

    canvas = assemble_circular_diagram(
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
        config_dict=config_dict,
        legend_config=legend_config,
        depth_df=resolved_depth_df,
        depth_config=depth_config,
        depth_tracks=resolved_depth_tracks,
        conservation_tracks=conservation_tracks,
        conservation_min_identity=float(identity),
        cfg=cfg,
        circular_track_slots=parsed_circular_track_slots,
        circular_track_axis_index=circular_track_axis_index,
        dinucleotide_content_dataframes=dinucleotide_content_dataframes,
        dinucleotide_skew_dataframes=dinucleotide_skew_dataframes,
        definition_position="center",
        definition_profile=effective_definition_profile,
        center_reserved_radius=center_reserved_radius,
        _tick_track_channel_override=_tick_track_channel_override,
    )
    if show_plot_title:
        _center_single_record_definition_group_on_axis(
            canvas,
            record_id=str(gb_record.id),
        )
        _add_single_record_plot_title_group(
            canvas=canvas,
            canvas_config=canvas_config,
            legend_config=legend_config,
            gb_record=gb_record,
            output_prefix=output_prefix,
            config_dict=config_dict,
            cfg=cfg,
            species=species,
            strain=strain,
            plot_title=normalized_plot_title or None,
            plot_title_position=normalized_plot_title_position,
        )
    return canvas


def assemble_circular_diagram_from_records(
    records: Sequence[SeqRecord],
    *,
    conservation_blast_files: Sequence[str] | None = None,
    conservation_dataframes: Sequence[DataFrame] | None = None,
    conservation_reference: Literal["query", "subject", "auto"] | str = "auto",
    conservation_labels: Sequence[str] | None = None,
    conservation_colors: Sequence[str] | None = None,
    conservation_ring_width: float | None = None,
    conservation_ring_gap: float | None = None,
    config_dict: dict | None = None,
    config_overrides: Mapping[str, object] | None = None,
    color_table: Optional[DataFrame] = None,
    color_table_file: str | None = None,
    default_colors: DataFrame | None = None,
    default_colors_palette: str = "default",
    default_colors_file: str | None = None,
    selected_features_set: Sequence[str] | None = None,
    feature_table: DataFrame | None = None,
    feature_table_file: str | None = None,
    feature_shapes: Mapping[str, str] | None = None,
    output_prefix: str = "out",
    legend: str = "right",
    dinucleotide: str = "GC",
    window: Optional[int] = None,
    step: Optional[int] = None,
    depth_window: Optional[int] = None,
    depth_step: Optional[int] = None,
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
    multi_record_size_mode: Literal["linear", "auto", "equal", "sqrt"] = "auto",
    multi_record_min_radius_ratio: float = 0.55,
    multi_record_column_gap_ratio: float = _MULTI_RECORD_COLUMN_GAP_RATIO,
    multi_record_row_gap_ratio: float = _MULTI_RECORD_ROW_GAP_RATIO,
    multi_record_positions: Sequence[str] | None = None,
    circular_track_slots: Sequence[str | CircularTrackSlot] | None = None,
    circular_track_axis_index: int | None = None,
    evalue: float = 1e-5,
    bitscore: float = 50.0,
    identity: float = 70.0,
    alignment_length: int = 0,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
    """Build and assemble a circular diagram grid from multiple records."""
    if not records:
        raise ValidationError("records is empty")
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

    if len(records) == 1:
        if (depth_table is not None or depth_file is not None) and (
            depth_tables is not None or depth_files is not None
        ):
            raise ValidationError("Use depth_table/depth_file or depth_tables/depth_files, not both.")
        single_depth_table = depth_table
        single_depth_file = depth_file
        if depth_tables is not None:
            if len(depth_tables) != 1:
                raise ValidationError("Expected one depth table for one circular record.")
            single_depth_table = depth_tables[0]
        if depth_files is not None:
            if len(depth_files) != 1:
                raise ValidationError("Expected one depth file for one circular record.")
            single_depth_file = depth_files[0]
        return assemble_circular_diagram_from_record(
            records[0],
            conservation_blast_files=conservation_blast_files,
            conservation_dataframes=conservation_dataframes,
            conservation_reference=conservation_reference,
            conservation_labels=conservation_labels,
            conservation_colors=conservation_colors,
            conservation_ring_width=conservation_ring_width,
            conservation_ring_gap=conservation_ring_gap,
            config_dict=config_dict,
            config_overrides=config_overrides,
            color_table=color_table,
            color_table_file=color_table_file,
            default_colors=default_colors,
            default_colors_palette=default_colors_palette,
            default_colors_file=default_colors_file,
            selected_features_set=selected_features_set,
            feature_table=feature_table,
            feature_table_file=feature_table_file,
            feature_shapes=feature_shapes,
            output_prefix=output_prefix,
            legend=legend,
            dinucleotide=dinucleotide,
            window=window,
            step=step,
            depth_window=depth_window,
            depth_step=depth_step,
            depth_table=single_depth_table,
            depth_file=single_depth_file,
            depth_track_tables=depth_track_tables,
            depth_track_files=depth_track_files,
            depth_track_labels=depth_track_labels,
            depth_track_colors=depth_track_colors,
            depth_track_large_tick_intervals=depth_track_large_tick_intervals,
            depth_track_small_tick_intervals=depth_track_small_tick_intervals,
            depth_track_tick_font_sizes=depth_track_tick_font_sizes,
            species=species,
            strain=strain,
            plot_title=normalized_plot_title or None,
            plot_title_position=normalized_plot_title_position,
            plot_title_font_size=plot_title_font_size,
            keep_full_definition_with_plot_title=keep_full_definition_with_plot_title,
            center_reserved_radius=center_reserved_radius,
            circular_track_slots=circular_track_slots,
            circular_track_axis_index=circular_track_axis_index,
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
            cfg=cfg,
        )

    if color_table is None and color_table_file is not None:
        color_table = read_color_table(color_table_file)
    if feature_table is None and feature_table_file is not None:
        feature_table = read_feature_visibility_file(feature_table_file)

    if default_colors is None:
        default_colors = load_default_colors(
            user_defined_default_colors=default_colors_file or "",
            palette=default_colors_palette or "default",
            load_comparison=False,
        )

    if config_dict is None:
        config_dict = load_config_toml("gbdraw.data", "config.toml")
    if config_overrides:
        if cfg is not None:
            raise ValueError(
                "config_overrides cannot be used with cfg; pass cfg=None or apply overrides before."
            )
        config_dict = modify_config_dict(config_dict, **config_overrides)
    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    config_dict, cfg = _apply_circular_plot_title_font_size_override(
        config_dict=config_dict,
        cfg=cfg,
        plot_title_font_size=plot_title_font_size,
    )

    if selected_features_set is None:
        selected_features_set = DEFAULT_SELECTED_FEATURES

    records = list(records)
    record_depth_tracks = normalize_depth_tracks(
        records,
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

    ordered_indices, row_counts = _resolve_multi_record_positions(
        records, multi_record_positions
    )
    records = [records[idx] for idx in ordered_indices]
    if record_depth_tracks is not None:
        record_depth_tracks = [record_depth_tracks[idx] for idx in ordered_indices]

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
        config_dict=config_dict,
        default_colors=default_colors,
        cfg=cfg,
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

    parsed_circular_track_slots = _parse_circular_track_slot_inputs(circular_track_slots)
    circular_track_axis_index = _validate_circular_track_axis_index(
        circular_track_axis_index,
        parsed_circular_track_slots,
    )

    legend_effective = legend

    if parsed_circular_track_slots is not None:
        show_depth = _circular_slots_have_renderer(parsed_circular_track_slots, "depth")
        show_gc = _circular_slots_have_renderer(parsed_circular_track_slots, "dinucleotide_content")
        show_skew = _circular_slots_have_renderer(parsed_circular_track_slots, "dinucleotide_skew")
        if show_depth and record_depth_tracks is None:
            raise ValidationError("A circular depth track slot requires depth_tables, depth_files, or depth_track input.")
    else:
        show_depth = cfg.canvas.show_depth
        show_gc = cfg.canvas.show_gc
        show_skew = cfg.canvas.show_skew
    available_depth_track_count = depth_track_count(record_depth_tracks)
    parsed_circular_track_slots = _default_circular_depth_slots_if_needed(
        parsed_circular_track_slots=parsed_circular_track_slots,
        show_depth=bool(show_depth and record_depth_tracks is not None),
        depth_track_count_value=available_depth_track_count,
        show_gc=bool(show_gc),
        show_skew=bool(show_skew),
        dinucleotide=dinucleotide,
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
                    show_depth=show_depth,
                    depth_track_count=max(1, depth_track_count(record_depth_tracks)),
                    show_gc=show_gc,
                    show_skew=show_skew,
                    dinucleotide=dinucleotide,
                )
            parsed_circular_track_slots = _insert_conservation_slots(
                parsed_circular_track_slots,
                conservation_slots,
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
    if show_plot_title and keep_full_definition_with_plot_title:
        record_definition_profile = "full"
    record_lengths = [len(record.seq) for record in records]
    cfg = _harmonize_multi_record_circular_style_cfg(
        cfg,
        record_lengths=record_lengths,
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
            config_dict=config_dict,
            cfg=cfg,
        )
        record_depth_track_data = build_depth_track_dataframes(
            records,
            record_depth_tracks,
            base_config=base_depth_config,
            depth_df_builder=build_depth_df,
            window_steps=record_depth_window_steps,
        )
        _validate_circular_depth_track_indices(
            parsed_circular_track_slots,
            depth_track_count_value=max(1, depth_track_data_count(record_depth_track_data)),
        )

    canvases: list[Drawing] = []
    widths: list[float] = []
    heights: list[float] = []
    record_radii_px: list[float] = []
    for record_index, (record, record_scale) in enumerate(zip(records, record_scales)):
        scaled_cfg = _scale_circular_cfg(cfg, scale=record_scale)
        record_radii_px.append(float(scaled_cfg.canvas.circular.radius))
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
            conservation_reference=conservation_mode,
            config_dict=config_dict,
            color_table=color_table,
            default_colors=default_colors,
            selected_features_set=list(selected_features_set),
            feature_table=feature_table,
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
            evalue=evalue,
            bitscore=bitscore,
            identity=identity,
            alignment_length=alignment_length,
            _definition_profile=record_definition_profile,
            _tick_track_channel_override=tick_track_channel_override,
            _precomputed_depth_tracks=record_depth_track_data[record_index],
            _precomputed_conservation_tracks=record_conservation_tracks,
            cfg=scaled_cfg,
        )
        canvases.append(sub_canvas)
        widths.append(_parse_svg_length_px(sub_canvas.attribs.get("width"), default=0.0))
        heights.append(_parse_svg_length_px(sub_canvas.attribs.get("height"), default=0.0))

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
        row_counts = _resolve_multi_record_default_row_counts(len(canvases))
    row_record_indices: list[list[int]] = []
    record_row_by_index: dict[int, int] = {}
    cursor = 0
    for row_count in row_counts:
        if cursor >= len(canvases):
            break
        row_indices: list[int] = []
        for _ in range(max(0, int(row_count))):
            if cursor >= len(canvases):
                break
            row_indices.append(cursor)
            record_row_by_index[cursor] = len(row_record_indices)
            cursor += 1
        if row_indices:
            row_record_indices.append(row_indices)
    rows = len(row_record_indices)
    row_heights: list[float] = [0.0] * rows
    record_sizes: list[tuple[float, float]] = [
        (float(width_px), float(height_px))
        for width_px, height_px in zip(widths, heights)
    ]
    for row, row_indices in enumerate(row_record_indices):
        for record_index in row_indices:
            if record_index >= len(record_sizes):
                continue
            _record_width, record_height = record_sizes[record_index]
            row_heights[row] = max(row_heights[row], float(record_height))

    record_vertical_insets: list[tuple[float, float]] = [
        _estimate_subcanvas_vertical_insets(sub_canvas) for sub_canvas in canvases
    ]
    row_content_tops: list[float] = [0.0] * rows
    row_content_bottoms: list[float] = [0.0] * rows
    for row, row_indices in enumerate(row_record_indices):
        if not row_indices:
            continue
        cell_height = float(row_heights[row]) if row < len(row_heights) else 0.0
        row_content_top: float | None = None
        row_content_bottom: float | None = None
        for record_index in row_indices:
            sub_height = float(record_sizes[record_index][1])
            top_inset = max(0.0, float(record_vertical_insets[record_index][0]))
            bottom_inset = max(0.0, float(record_vertical_insets[record_index][1]))
            if sub_height > 0.0:
                if top_inset > sub_height:
                    top_inset = sub_height
                if top_inset + bottom_inset > sub_height:
                    bottom_inset = max(0.0, sub_height - top_inset)
            content_height = max(0.0, sub_height - top_inset - bottom_inset)
            cell_offset_y = (cell_height - sub_height) * 0.5
            content_top = cell_offset_y + top_inset
            content_bottom = content_top + content_height
            row_content_top = content_top if row_content_top is None else min(row_content_top, content_top)
            row_content_bottom = content_bottom if row_content_bottom is None else max(row_content_bottom, content_bottom)
        row_content_tops[row] = float(row_content_top) if row_content_top is not None else 0.0
        row_content_bottoms[row] = float(row_content_bottom) if row_content_bottom is not None else 0.0

    row_offsets: list[float] = [0.0] * rows
    for row in range(1, rows):
        previous_row = row - 1
        required_shift = (
            row_content_bottoms[previous_row]
            + float(row_gap_px)
            - row_content_tops[row]
        )
        row_offsets[row] = float(row_offsets[previous_row]) + max(0.0, float(required_shift))
    row_physical_widths: list[float] = []
    row_total_widths: list[float] = []
    row_outer_left_paddings: list[float] = []
    row_outer_right_paddings: list[float] = []
    record_horizontal_insets: list[tuple[float, float]] = [
        _estimate_subcanvas_horizontal_insets(sub_canvas) for sub_canvas in canvases
    ]
    record_content_widths: list[float] = []
    for record_index, (sub_width, _sub_height) in enumerate(record_sizes):
        left_inset = max(0.0, float(record_horizontal_insets[record_index][0]))
        right_inset = max(0.0, float(record_horizontal_insets[record_index][1]))
        if sub_width > 0.0:
            if left_inset > sub_width:
                left_inset = sub_width
            if left_inset + right_inset > sub_width:
                right_inset = max(0.0, sub_width - left_inset)
        record_horizontal_insets[record_index] = (float(left_inset), float(right_inset))
        record_content_widths.append(
            max(0.0, float(sub_width) - float(left_inset) - float(right_inset))
        )
    apply_row_margin_symmetry = (
        str(legend_effective).strip().lower() in {"none", "top", "bottom"}
    )
    for row_indices in row_record_indices:
        if not row_indices:
            row_physical_widths.append(0.0)
            row_total_widths.append(0.0)
            row_outer_left_paddings.append(0.0)
            row_outer_right_paddings.append(0.0)
            continue
        row_content_width = sum(
            float(record_content_widths[index]) for index in row_indices
        )
        row_content_width += max(0, len(row_indices) - 1) * column_gap_px

        first_index = row_indices[0]
        last_index = row_indices[-1]
        row_left = max(0.0, float(record_horizontal_insets[first_index][0]))
        row_right = max(0.0, float(record_horizontal_insets[last_index][1]))
        row_physical_width = float(row_left) + float(row_content_width) + float(row_right)
        row_physical_widths.append(float(row_physical_width))

        extra_left = 0.0
        extra_right = 0.0
        if apply_row_margin_symmetry:
            target_margin = max(row_left, row_right)
            extra_left = max(0.0, target_margin - row_left)
            extra_right = max(0.0, target_margin - row_right)

        row_outer_left_paddings.append(float(extra_left))
        row_outer_right_paddings.append(float(extra_right))
        row_total_widths.append(float(row_physical_width + extra_left + extra_right))

    grid_width = max(row_total_widths, default=0.0)
    grid_height = max(
        (
            float(row_offsets[row]) + float(row_heights[row])
            for row in range(rows)
        ),
        default=0.0,
    )
    record_offsets_x: dict[int, float] = {}
    for row, row_indices in enumerate(row_record_indices):
        if not row_indices:
            continue
        row_physical_width = (
            row_physical_widths[row] if row < len(row_physical_widths) else 0.0
        )
        row_total_width = (
            row_total_widths[row] if row < len(row_total_widths) else row_physical_width
        )
        row_outer_left_padding = (
            row_outer_left_paddings[row] if row < len(row_outer_left_paddings) else 0.0
        )
        first_index = row_indices[0]
        content_cursor_x = (
            max(0.0, (grid_width - row_total_width) * 0.5)
            + row_outer_left_padding
            + float(record_horizontal_insets[first_index][0])
        )
        for position, record_index in enumerate(row_indices):
            left_inset = float(record_horizontal_insets[record_index][0])
            content_width = float(record_content_widths[record_index])
            record_offsets_x[record_index] = content_cursor_x - left_inset
            content_cursor_x += content_width
            if position < len(row_indices) - 1:
                content_cursor_x += column_gap_px

    total_width = grid_width
    total_height = grid_height
    grid_origin_x = 0.0
    grid_origin_y = 0.0
    legend_offset_x = 0.0
    legend_offset_y = 0.0
    legend_config: LegendDrawingConfigurator | None = None
    legend_canvas_config: CircularCanvasConfigurator | None = None
    legend_table: dict = {}
    legend_local_top = 0.0
    legend_local_bottom = 0.0
    plot_title_group: Group | None = None
    plot_title_local_bounds = (0.0, 0.0)

    if show_plot_title:
        shared_canvas_config = CircularCanvasConfigurator(
            output_prefix=output_prefix,
            config_dict=config_dict,
            legend="none",
            gb_record=records[0],
            cfg=cfg,
        )
        plot_title_group = DefinitionGroup(
            gb_record=records[0],
            canvas_config=shared_canvas_config,
            config_dict=config_dict,
            species=species,
            strain=strain,
            plot_title=normalized_plot_title or None,
            definition_profile="shared_common",
            definition_group_id="plot_title",
            cfg=cfg,
        ).get_group()
        plot_title_local_bounds = _group_local_vertical_bounds(plot_title_group)

    if legend_effective != "none":
        legend_canvas_config = CircularCanvasConfigurator(
            output_prefix=output_prefix,
            config_dict=config_dict,
            legend=legend_effective,
            gb_record=records[0],
            cfg=cfg,
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
            config_dict=config_dict,
            default_colors_df=default_colors,
            cfg=cfg,
        )
        skew_config = GcSkewConfigurator(
            window=legend_window,
            step=legend_step,
            dinucleotide=dinucleotide,
            config_dict=config_dict,
            default_colors_df=default_colors,
            cfg=cfg,
        )
        depth_config = (
            DepthConfigurator(
                window=legend_window,
                step=legend_step,
                config_dict=config_dict,
                cfg=cfg,
            )
            if cfg.canvas.show_depth
            else None
        )
        feature_config = FeatureDrawingConfigurator(
            color_table=color_table,
            default_colors=default_colors,
            selected_features_set=list(selected_features_set),
            feature_table=feature_table,
            feature_shapes=feature_shapes,
            config_dict=config_dict,
            canvas_config=legend_canvas_config,
            cfg=cfg,
        )
        features_present = check_feature_presence(
            list(records),
            list(selected_features_set),
            feature_visibility_rules=feature_config.feature_visibility_rules,
        )
        color_map, default_color_map = preprocess_color_tables(
            feature_config.color_table,
            feature_config.default_colors,
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
        )
        if cfg.canvas.show_depth:
            first_depth_row = next((row for row in record_depth_track_data if row), [])
            legend_table = sync_depth_track_legend_entries(legend_table, first_depth_row)
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
                config_dict=config_dict,
                gc_config=gc_config,
                skew_config=skew_config,
                feature_config=feature_config,
                canvas_config=legend_canvas_config,
                cfg=cfg,
            )
            legend_config = legend_config.recalculate_legend_dimensions(
                legend_table, legend_canvas_config
            )
            legend_local_top = -0.5 * float(legend_config.color_rect_size)
            legend_local_bottom = legend_local_top + float(legend_config.legend_height)
            side_inner_gap = 0.0
            side_edge_margin = 0.0
            side_reserved_width = 0.0
            if legend_effective in {"left", "right"}:
                side_inner_gap, side_edge_margin, side_reserved_width = (
                    resolve_circular_side_legend_geometry(
                        canvas_height=float(total_height),
                        legend_width=float(legend_config.legend_width),
                        color_rect_size=float(legend_config.color_rect_size),
                    )
                )

            if legend_effective == "right":
                total_width = grid_width + side_reserved_width
                legend_offset_x = grid_width + side_inner_gap
                legend_offset_y = (total_height - legend_config.legend_height) / 2.0
            elif legend_effective == "left":
                total_width = grid_width + side_reserved_width
                grid_origin_x = side_reserved_width
                legend_offset_x = side_edge_margin
                legend_offset_y = (total_height - legend_config.legend_height) / 2.0
            elif legend_effective == "top":
                legend_offset_y = _MULTI_RECORD_LEGEND_TOP_EDGE_PADDING_PX - legend_local_top
                legend_bottom = legend_offset_y + legend_local_bottom
                grid_origin_y = legend_bottom + _MULTI_RECORD_LEGEND_GRID_GAP_PX
                total_height = grid_origin_y + grid_height
                legend_offset_x = (total_width - legend_config.legend_width) / 2.0
            elif legend_effective == "bottom":
                legend_offset_x = (total_width - legend_config.legend_width) / 2.0
                legend_offset_y = grid_height + _MULTI_RECORD_LEGEND_GRID_GAP_PX - legend_local_top
                legend_bottom = legend_offset_y + legend_local_bottom
                total_height = max(
                    total_height,
                    legend_bottom + _MULTI_RECORD_LEGEND_EDGE_PADDING_PX,
                )
            elif legend_effective == "upper_left":
                legend_offset_x = 0.025 * total_width
                legend_offset_y = 0.05 * total_height
            elif legend_effective == "upper_right":
                legend_offset_x = 0.85 * total_width
                legend_offset_y = 0.05 * total_height
            elif legend_effective == "lower_left":
                legend_offset_x = 0.025 * total_width
                legend_offset_y = 0.78 * total_height
            elif legend_effective == "lower_right":
                legend_offset_x = 0.875 * total_width
                legend_offset_y = 0.75 * total_height

    if (
        show_plot_title
        and plot_title_group is not None
        and normalized_plot_title_position == "bottom"
    ):
        plot_title_min_y, plot_title_max_y = plot_title_local_bounds
        plot_title_height = max(0.0, float(plot_title_max_y) - float(plot_title_min_y))

        records_bottom = float(grid_origin_y) + float(grid_height)
        anchor_bottom = records_bottom
        if legend_effective == "bottom" and legend_config is not None:
            legend_bottom = legend_offset_y + legend_local_bottom
            anchor_bottom = max(anchor_bottom, float(legend_bottom))

        required_height = (
            anchor_bottom
            + _MULTI_RECORD_LEGEND_PLOT_TITLE_GAP_PX
            + plot_title_height
            + _MULTI_RECORD_PLOT_TITLE_BOTTOM_MARGIN_PX
        )
        total_height = max(total_height, required_height)

    merged_canvas = Drawing(
        filename=f"{output_prefix}.svg",
        size=(f"{total_width}px", f"{total_height}px"),
        viewBox=f"0 0 {total_width} {total_height}",
        debug=False,
    )

    for record_index, sub_canvas in enumerate(canvases):
        row = int(record_row_by_index.get(record_index, 0))
        default_size = record_sizes[record_index] if record_index < len(record_sizes) else (0.0, 0.0)
        sub_width = _parse_svg_length_px(sub_canvas.attribs.get("width"), default=default_size[0])
        sub_height = _parse_svg_length_px(sub_canvas.attribs.get("height"), default=default_size[1])
        cell_height = row_heights[row] if row < len(row_heights) else sub_height
        cell_origin_x = grid_origin_x + record_offsets_x.get(record_index, 0.0)
        cell_origin_y = grid_origin_y + (row_offsets[row] if row < len(row_offsets) else 0.0)
        cell_offset_x = cell_origin_x
        cell_offset_y = cell_origin_y + (cell_height - sub_height) * 0.5
        record_group = Group(id=f"record_{record_index}")
        record_group.translate(cell_offset_x, cell_offset_y)

        for element in sub_canvas.elements:
            if _is_defs_element(element):
                continue
            copied = copy.deepcopy(element)
            _suffix_fixed_top_level_group_id(copied, record_index)
            record_group.add(copied)
        _center_record_definition_group_on_record_axis(
            record_group,
            record_index=record_index,
            record_id=str(records[record_index].id),
        )
        merged_canvas.add(record_group)

    if plot_title_group is not None:
        plot_title_group = place_definition_group_on_size(
            plot_title_group,
            canvas_width=float(total_width),
            canvas_height=float(total_height),
            position=normalized_plot_title_position,
        )
        merged_canvas.add(plot_title_group)

    if (
        legend_effective != "none"
        and legend_table
        and legend_config is not None
        and legend_canvas_config is not None
    ):
        legend_group = LegendGroup(
            legend_canvas_config,
            legend_config,
            legend_table,
        ).get_group()
        legend_group.translate(legend_offset_x, legend_offset_y)
        merged_canvas.add(legend_group)

    return merged_canvas


def build_circular_diagram(
    gb_record: SeqRecord,
    *,
    options: DiagramOptions | None = None,
) -> Drawing:
    """Build a circular diagram using bundled DiagramOptions."""

    options = options or DiagramOptions()
    colors = options.colors
    output = options.output
    tracks = options.tracks
    tracks = options.tracks

    config_dict: dict | None = None
    cfg: GbdrawConfig | None = None
    config_overrides = options.config_overrides
    if isinstance(options.config, GbdrawConfig):
        if config_overrides:
            cfg = apply_config_overrides(options.config, config_overrides)
            config_overrides = None
        else:
            cfg = options.config
    elif isinstance(options.config, dict):
        config_dict = options.config

    return assemble_circular_diagram_from_record(
        gb_record,
        config_dict=config_dict,
        config_overrides=config_overrides,
        color_table=colors.color_table if colors else None,
        color_table_file=colors.color_table_file if colors else None,
        default_colors=colors.default_colors if colors else None,
        default_colors_palette=colors.default_colors_palette if colors else "default",
        default_colors_file=colors.default_colors_file if colors else None,
        selected_features_set=options.selected_features_set,
        feature_table=options.feature_table,
        feature_table_file=options.feature_table_file,
        feature_shapes=options.feature_shapes,
        output_prefix=output.output_prefix if output else "out",
        legend=output.legend if output else "right",
        dinucleotide=options.dinucleotide,
        window=options.window,
        step=options.step,
        depth_window=options.depth_window,
        depth_step=options.depth_step,
        depth_table=(
            options.depth_table
            if options.depth_table is not None
            else (options.depth_tables[0] if options.depth_tables else None)
        ),
        depth_file=(
            options.depth_file
            if options.depth_file is not None
            else (options.depth_files[0] if options.depth_files else None)
        ),
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
        cfg=cfg,
    )


def build_linear_diagram(
    records: Sequence[SeqRecord],
    *,
    options: DiagramOptions | None = None,
) -> Drawing:
    """Build a linear diagram using bundled DiagramOptions."""

    options = options or DiagramOptions()
    colors = options.colors
    output = options.output
    tracks = options.tracks
    config_dict: dict | None = None
    cfg: GbdrawConfig | None = None
    config_overrides = options.config_overrides
    if isinstance(options.config, GbdrawConfig):
        if config_overrides:
            cfg = apply_config_overrides(options.config, config_overrides)
            config_overrides = None
        else:
            cfg = options.config
    elif isinstance(options.config, dict):
        config_dict = options.config
    normalize_collinearity_anchor_mode(str(options.collinearity_anchor_mode))

    return assemble_linear_diagram_from_records(
        records,
        blast_files=options.blast_files,
        protein_comparisons=options.protein_comparisons,
        orthogroups=options.orthogroups,
        protein_blastp_mode=options.protein_blastp_mode,
        pairwise_match_style=options.pairwise_match_style,
        collinearity_blocks=options.collinearity_blocks,
        collinearity_params=options.collinearity_params,
        collinearity_unit_mode=options.collinearity_unit_mode,
        collinearity_anchor_mode="rbh",
        collinearity_search_scope=options.collinearity_search_scope,
        collinearity_color_mode=options.collinearity_color_mode,
        losatp_bin=options.losatp_bin,
        ncbi_blastp_bin=options.ncbi_blastp_bin,
        losatp_threads=options.losatp_threads,
        protein_blastp_max_hits=options.protein_blastp_max_hits,
        protein_blastp_candidate_limit=options.protein_blastp_candidate_limit,
        orthogroup_membership_mode=options.orthogroup_membership_mode,
        orthogroup_member_max_hits=options.orthogroup_member_max_hits,
        collinear_max_paralog_links_per_orthogroup=options.collinear_max_paralog_links_per_orthogroup,
        align_orthogroup_feature=options.align_orthogroup_feature,
        config_dict=config_dict,
        config_overrides=config_overrides,
        color_table=colors.color_table if colors else None,
        color_table_file=colors.color_table_file if colors else None,
        default_colors=colors.default_colors if colors else None,
        default_colors_palette=colors.default_colors_palette if colors else "default",
        default_colors_file=colors.default_colors_file if colors else None,
        selected_features_set=options.selected_features_set,
        feature_table=options.feature_table,
        feature_table_file=options.feature_table_file,
        feature_shapes=options.feature_shapes,
        output_prefix=output.output_prefix if output else "out",
        legend=output.legend if output else "right",
        dinucleotide=options.dinucleotide,
        window=options.window,
        step=options.step,
        depth_window=options.depth_window,
        depth_step=options.depth_step,
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
        cfg=cfg,
    )


__all__ = [
    "DEFAULT_SELECTED_FEATURES",
    "assemble_circular_diagram_from_record",
    "assemble_circular_diagram_from_records",
    "assemble_linear_diagram_from_records",
    "build_circular_diagram",
    "build_linear_diagram",
]
