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
from dataclasses import replace
from typing import Optional, Sequence, Mapping, Literal, cast

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from gbdraw.analysis.skew import skew_df  # type: ignore[reportMissingImports]
from gbdraw.api.config import apply_config_overrides  # type: ignore[reportMissingImports]
from gbdraw.api.options import DiagramOptions  # type: ignore[reportMissingImports]
from gbdraw.canvas import CircularCanvasConfigurator, LinearCanvasConfigurator  # type: ignore[reportMissingImports]
from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from gbdraw.config.modify import modify_config_dict  # type: ignore[reportMissingImports]
from gbdraw.config.toml import load_config_toml  # type: ignore[reportMissingImports]
from gbdraw.io.colors import load_default_colors, read_color_table  # type: ignore[reportMissingImports]
from gbdraw.features.visibility import read_feature_visibility_file  # type: ignore[reportMissingImports]
from gbdraw.configurators import (  # type: ignore[reportMissingImports]
    BlastMatchConfigurator,
    FeatureDrawingConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
    LegendDrawingConfigurator,
)
from gbdraw.core.sequence import create_dict_for_sequence_lengths, check_feature_presence  # type: ignore[reportMissingImports]
from gbdraw.diagrams.circular import assemble_circular_diagram  # type: ignore[reportMissingImports]
from gbdraw.diagrams.linear import assemble_linear_diagram  # type: ignore[reportMissingImports]
from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]
from gbdraw.features.colors import preprocess_color_tables, precompute_used_color_rules  # type: ignore[reportMissingImports]
from gbdraw.legend.table import prepare_legend_table  # type: ignore[reportMissingImports]
from gbdraw.render.groups.circular import LegendGroup  # type: ignore[reportMissingImports]
from gbdraw.tracks import TrackSpec, parse_track_specs  # type: ignore[reportMissingImports]

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

_SUPPORTED_CIRCULAR_TRACK_KINDS = {
    "features",
    "gc_content",
    "gc_skew",
    "ticks",
    "axis",
    "legend",
    "labels",
}

_MULTI_RECORD_SUFFIXED_TOP_LEVEL_IDS = {
    "Axis",
    "tick",
    "labels",
    "gc_content",
    "skew",
    "gc_skew",
}
_MULTI_RECORD_SIZE_MODES = {"linear", "sqrt", "equal"}


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


def _estimate_square_grid(record_count: int) -> tuple[int, int]:
    """Return grid dimensions (cols, rows) close to square."""
    if record_count <= 0:
        return 1, 1
    cols = int(math.ceil(math.sqrt(record_count)))
    rows = int(math.ceil(float(record_count) / float(cols)))
    return cols, rows


def _suffix_fixed_top_level_group_id(element: object, record_index: int) -> None:
    """Suffix selected top-level group IDs to avoid collisions on merged canvas."""
    attribs = getattr(element, "attribs", None)
    if not isinstance(attribs, dict):
        return
    group_id = attribs.get("id")
    if group_id in _MULTI_RECORD_SUFFIXED_TOP_LEVEL_IDS:
        attribs["id"] = f"{group_id}_{record_index}"


def _is_defs_element(element: object) -> bool:
    """Return True if the element is an SVG <defs> container."""
    class_name = element.__class__.__name__.lower()
    if class_name == "defs":
        return True
    element_name = getattr(element, "elementname", None)
    if isinstance(element_name, str) and element_name.lower() == "defs":
        return True
    return False


def _resolve_multi_record_size_mode(mode: str) -> Literal["linear", "sqrt", "equal"]:
    """Normalize and validate multi-record size mode."""
    normalized = str(mode).strip().lower()
    if normalized not in _MULTI_RECORD_SIZE_MODES:
        raise ValidationError(
            "multi_record_size_mode must be one of: linear, sqrt, equal"
        )
    return cast(Literal["linear", "sqrt", "equal"], normalized)


def _validate_multi_record_min_radius_ratio(value: float) -> float:
    """Validate minimum radius ratio for multi-record scaling."""
    ratio = float(value)
    if not (0.0 < ratio <= 1.0):
        raise ValidationError("multi_record_min_radius_ratio must be > 0 and <= 1")
    return ratio


def _resolve_multi_record_scale(
    record_length: int,
    max_record_length: int,
    *,
    mode: Literal["linear", "sqrt", "equal"],
    min_radius_ratio: float,
) -> float:
    """Return per-record circular scale for multi-record canvas rendering."""
    if mode == "equal":
        return 1.0

    if max_record_length <= 0:
        ratio = 1.0
    else:
        ratio = max(0.0, float(record_length) / float(max_record_length))

    if mode == "linear":
        scale = ratio
    else:
        scale = math.sqrt(ratio)

    return max(float(min_radius_ratio), min(float(scale), 1.0))


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
    evalue: float = 1e-5,
    bitscore: float = 50.0,
    identity: float = 70.0,
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
    if color_table is None and color_table_file is not None:
        color_table = read_color_table(color_table_file)
    if feature_table is None and feature_table_file is not None:
        feature_table = read_feature_visibility_file(feature_table_file)

    if default_colors is None:
        default_colors = load_default_colors(
            user_defined_default_colors=default_colors_file or "",
            palette=default_colors_palette or "default",
            load_comparison=bool(blast_files),
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

    if selected_features_set is None:
        selected_features_set = DEFAULT_SELECTED_FEATURES

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

    blast_config = BlastMatchConfigurator(
        evalue=evalue,
        bitscore=bitscore,
        identity=identity,
        sequence_length_dict=seq_len_dict,
        config_dict=config_dict,
        default_colors_df=default_colors,
        cfg=cfg,
    )

    canvas_config = LinearCanvasConfigurator(
        num_of_entries=len(records),
        longest_genome=longest_genome,
        config_dict=config_dict,
        legend=legend,
        output_prefix=output_prefix,
        cfg=cfg,
        has_comparisons=bool(blast_files),
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

    return assemble_linear_diagram(
        records=list(records),
        blast_files=list(blast_files) if blast_files else None,
        canvas_config=canvas_config,
        blast_config=blast_config,
        feature_config=feature_config,
        gc_config=gc_config,
        config_dict=config_dict,
        legend_config=legend_config,
        skew_config=skew_config,
        cfg=cfg,
    )


def assemble_circular_diagram_from_record(
    gb_record: SeqRecord,
    *,
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
    species: Optional[str] = None,
    strain: Optional[str] = None,
    track_specs: Sequence[str | TrackSpec] | None = None,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
    """Builds and assembles a circular diagram for a single record.

    If config_dict is None, it loads gbdraw.data/config.toml.
    If config_overrides is provided, modify_config_dict is applied.
    If default_colors is None, it loads the built-in default palette.
    If color_table is None and color_table_file is provided, it is loaded.
    If selected_features_set is None, it uses the CLI default feature list.
    """
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

    if selected_features_set is None:
        selected_features_set = DEFAULT_SELECTED_FEATURES

    parsed_track_specs: list[TrackSpec] | None = None
    if track_specs is not None:
        parsed: list[TrackSpec] = []
        raw: list[str] = []
        for item in track_specs:
            if isinstance(item, TrackSpec):
                parsed.append(item)
            else:
                raw.append(str(item))
        if raw:
            parsed.extend(parse_track_specs(raw, mode="circular"))
        parsed_track_specs = parsed

    if parsed_track_specs:
        for ts in parsed_track_specs:
            if ts.mode != "circular":
                raise ValidationError(
                    f"TrackSpec mode '{ts.mode}' is not supported for circular diagrams."
                )
            if str(ts.kind) not in _SUPPORTED_CIRCULAR_TRACK_KINDS:
                logger.warning(
                    "TrackSpec kind '%s' is not supported for circular diagrams yet; it will be ignored.",
                    ts.kind,
                )

    ts_by_kind = {str(ts.kind): ts for ts in (parsed_track_specs or [])}
    legend_ts = ts_by_kind.get("legend")
    legend_effective = "none" if (legend_ts is not None and not legend_ts.show) else legend

    # Allow track specs to override high-level show flags (used by canvas sizing and track IDs).
    show_gc = ts_by_kind.get("gc_content").show if "gc_content" in ts_by_kind else cfg.canvas.show_gc
    show_skew = ts_by_kind.get("gc_skew").show if "gc_skew" in ts_by_kind else cfg.canvas.show_skew
    canvas_cfg = cfg.canvas
    canvas_cfg = replace(canvas_cfg, show_gc=bool(show_gc), show_skew=bool(show_skew))
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

    # Circular drawing expects the precomputed GC/skew dataframe, but only when needed.
    gc_df = skew_df(gb_record, window, step, dinucleotide) if (cfg.canvas.show_gc or cfg.canvas.show_skew) else DataFrame()

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

    return assemble_circular_diagram(
        gb_record=gb_record,
        canvas_config=canvas_config,
        gc_df=gc_df,
        gc_config=gc_config,
        skew_config=skew_config,
        feature_config=feature_config,
        species=species,
        strain=strain,
        config_dict=config_dict,
        legend_config=legend_config,
        cfg=cfg,
        track_specs=parsed_track_specs,
    )


def assemble_circular_diagram_from_records(
    records: Sequence[SeqRecord],
    *,
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
    species: Optional[str] = None,
    strain: Optional[str] = None,
    multi_record_size_mode: Literal["linear", "sqrt", "equal"] = "sqrt",
    multi_record_min_radius_ratio: float = 0.55,
    track_specs: Sequence[str | TrackSpec] | None = None,
    cfg: GbdrawConfig | None = None,
) -> Drawing:
    """Build and assemble a circular diagram grid from multiple records."""
    if not records:
        raise ValidationError("records is empty")

    normalized_multi_record_size_mode = _resolve_multi_record_size_mode(
        str(multi_record_size_mode)
    )
    normalized_multi_record_min_radius_ratio = _validate_multi_record_min_radius_ratio(
        float(multi_record_min_radius_ratio)
    )

    if len(records) == 1:
        return assemble_circular_diagram_from_record(
            records[0],
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
            species=species,
            strain=strain,
            track_specs=track_specs,
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

    if selected_features_set is None:
        selected_features_set = DEFAULT_SELECTED_FEATURES

    parsed_track_specs: list[TrackSpec] | None = None
    if track_specs is not None:
        parsed: list[TrackSpec] = []
        raw: list[str] = []
        for item in track_specs:
            if isinstance(item, TrackSpec):
                parsed.append(item)
            else:
                raw.append(str(item))
        if raw:
            parsed.extend(parse_track_specs(raw, mode="circular"))
        parsed_track_specs = parsed

    if parsed_track_specs:
        for ts in parsed_track_specs:
            if ts.mode != "circular":
                raise ValidationError(
                    f"TrackSpec mode '{ts.mode}' is not supported for circular diagrams."
                )
            if str(ts.kind) not in _SUPPORTED_CIRCULAR_TRACK_KINDS:
                logger.warning(
                    "TrackSpec kind '%s' is not supported for circular diagrams yet; it will be ignored.",
                    ts.kind,
                )

    ts_by_kind = {str(ts.kind): ts for ts in (parsed_track_specs or [])}
    legend_ts = ts_by_kind.get("legend")
    legend_effective = "none" if (legend_ts is not None and not legend_ts.show) else legend

    show_gc = ts_by_kind.get("gc_content").show if "gc_content" in ts_by_kind else cfg.canvas.show_gc
    show_skew = ts_by_kind.get("gc_skew").show if "gc_skew" in ts_by_kind else cfg.canvas.show_skew
    canvas_cfg = replace(cfg.canvas, show_gc=bool(show_gc), show_skew=bool(show_skew))
    cfg = replace(cfg, canvas=canvas_cfg)
    max_record_length = max(len(record.seq) for record in records)

    canvases: list[Drawing] = []
    widths: list[float] = []
    heights: list[float] = []
    for record in records:
        record_scale = _resolve_multi_record_scale(
            len(record.seq),
            max_record_length,
            mode=normalized_multi_record_size_mode,
            min_radius_ratio=normalized_multi_record_min_radius_ratio,
        )
        scaled_cfg = _scale_circular_cfg(cfg, scale=record_scale)
        sub_canvas = assemble_circular_diagram_from_record(
            record,
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
            species=species,
            strain=strain,
            track_specs=parsed_track_specs,
            cfg=scaled_cfg,
        )
        canvases.append(sub_canvas)
        widths.append(_parse_svg_length_px(sub_canvas.attribs.get("width"), default=0.0))
        heights.append(_parse_svg_length_px(sub_canvas.attribs.get("height"), default=0.0))

    cell_width = max(widths) if widths else 0.0
    cell_height = max(heights) if heights else 0.0
    cols, rows = _estimate_square_grid(len(canvases))
    grid_width = float(cols) * float(cell_width)
    grid_height = float(rows) * float(cell_height)

    total_width = grid_width
    total_height = grid_height
    grid_origin_x = 0.0
    legend_offset_x = 0.0
    legend_offset_y = 0.0
    legend_config: LegendDrawingConfigurator | None = None
    legend_canvas_config: CircularCanvasConfigurator | None = None
    legend_table: dict = {}

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
        )
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

            if legend_effective == "right":
                total_width = grid_width + (legend_config.legend_width * 1.1)
                legend_offset_x = grid_width + (legend_config.legend_width * 0.05)
                legend_offset_y = (total_height - legend_config.legend_height) / 2.0
            elif legend_effective == "left":
                total_width = grid_width + (legend_config.legend_width * 1.1)
                grid_origin_x = legend_config.legend_width * 1.1
                legend_offset_x = legend_config.legend_width * 0.05
                legend_offset_y = (total_height - legend_config.legend_height) / 2.0
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

    merged_canvas = Drawing(
        filename=f"{output_prefix}.svg",
        size=(f"{total_width}px", f"{total_height}px"),
        viewBox=f"0 0 {total_width} {total_height}",
        debug=False,
    )

    for record_index, sub_canvas in enumerate(canvases):
        row = record_index // cols
        col = record_index % cols
        sub_width = _parse_svg_length_px(sub_canvas.attribs.get("width"), default=cell_width)
        sub_height = _parse_svg_length_px(sub_canvas.attribs.get("height"), default=cell_height)
        cell_offset_x = grid_origin_x + float(col) * cell_width + (cell_width - sub_width) * 0.5
        cell_offset_y = float(row) * cell_height + (cell_height - sub_height) * 0.5
        record_group = Group(id=f"record_{record_index}")
        record_group.translate(cell_offset_x, cell_offset_y)

        for element in sub_canvas.elements:
            if _is_defs_element(element):
                continue
            copied = copy.deepcopy(element)
            _suffix_fixed_top_level_group_id(copied, record_index)
            record_group.add(copied)
        merged_canvas.add(record_group)

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
        species=options.species,
        strain=options.strain,
        track_specs=tracks.track_specs if tracks else None,
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

    if tracks and tracks.track_specs:
        logger.warning(
            "Track specs are not supported for linear diagrams yet; ignoring track_specs."
        )

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

    return assemble_linear_diagram_from_records(
        records,
        blast_files=options.blast_files,
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
        evalue=options.evalue,
        bitscore=options.bitscore,
        identity=options.identity,
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


