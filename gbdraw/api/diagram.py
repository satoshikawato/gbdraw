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
from gbdraw.diagrams.circular.positioning import place_definition_group_on_size  # type: ignore[reportMissingImports]
from gbdraw.diagrams.linear import assemble_linear_diagram  # type: ignore[reportMissingImports]
from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]
from gbdraw.features.colors import preprocess_color_tables, precompute_used_color_rules  # type: ignore[reportMissingImports]
from gbdraw.legend.table import prepare_legend_table  # type: ignore[reportMissingImports]
from gbdraw.render.groups.circular import DefinitionGroup, LegendGroup  # type: ignore[reportMissingImports]
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
    "definition",
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
_DEFINITION_POSITIONS = {"center", "top", "bottom"}
_MULTI_RECORD_DEFINITION_MODES = {"shared", "legacy"}
# Legacy fallback used only when a valid radius cannot be derived.
_MULTI_RECORD_GRID_GAP_PX = 16.0
_MULTI_RECORD_GRID_GAP_RATIO = 0.10
_MULTI_RECORD_LEGEND_EDGE_PADDING_PX = 20.0
_MULTI_RECORD_LEGEND_GRID_GAP_PX = 20.0
_MULTI_RECORD_LEGEND_TOP_EDGE_PADDING_PX = 32.0
_MULTI_RECORD_LEGEND_SHARED_GAP_PX = 20.0
_MULTI_RECORD_SHARED_BOTTOM_MARGIN_PX = 24.0
_SVG_NUMBER_PATTERN = re.compile(r"[-+]?(?:\d*\.?\d+)(?:[eE][-+]?\d+)?")
_SVG_TRANSLATE_PATTERN = re.compile(
    r"translate\(\s*([-+0-9.eE]+)(?:[\s,]+([-+0-9.eE]+))?\s*\)"
)


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


def _estimate_subcanvas_horizontal_insets(sub_canvas: Drawing) -> tuple[float, float]:
    """Estimate left/right insets between viewBox edges and diagram geometry."""
    try:
        root = ET.fromstring(sub_canvas.tostring())
    except ET.ParseError:
        return 0.0, 0.0

    view_box = str(root.attrib.get("viewBox", "")).strip()
    view_box_parts = [part for part in view_box.split() if part]
    if len(view_box_parts) == 4:
        vb_x = _parse_svg_length_px(view_box_parts[0], default=0.0)
        vb_w = _parse_svg_length_px(view_box_parts[2], default=0.0)
    else:
        vb_x = 0.0
        vb_w = _parse_svg_length_px(root.attrib.get("width"), default=0.0)

    if vb_w <= 0.0:
        return 0.0, 0.0

    min_x: float | None = None
    max_x: float | None = None

    def _walk(node: ET.Element, inherited_tx: float) -> None:
        nonlocal min_x, max_x
        tag = _svg_local_name(str(node.tag)).lower()
        if tag == "defs":
            return
        if tag == "g" and str(node.attrib.get("id", "")) == "legend":
            return

        local_tx, _local_ty = _parse_translate_xy(node.attrib.get("transform"))
        total_tx = inherited_tx + float(local_tx)

        bounds = _estimate_element_local_x_bounds(node)
        if bounds is not None:
            candidate_min = total_tx + float(bounds[0])
            candidate_max = total_tx + float(bounds[1])
            min_x = candidate_min if min_x is None else min(min_x, candidate_min)
            max_x = candidate_max if max_x is None else max(max_x, candidate_max)

        for child in list(node):
            _walk(child, total_tx)

    _walk(root, 0.0)

    if min_x is None or max_x is None:
        return 0.0, 0.0

    view_min_x = float(vb_x)
    view_max_x = float(vb_x) + float(vb_w)
    left_inset = max(0.0, float(min_x) - view_min_x)
    right_inset = max(0.0, view_max_x - float(max_x))
    return float(left_inset), float(right_inset)


def _estimate_square_grid(record_count: int) -> tuple[int, int]:
    """Return grid dimensions (cols, rows) close to square."""
    if record_count <= 0:
        return 1, 1
    cols = int(math.ceil(math.sqrt(record_count)))
    rows = int(math.ceil(float(record_count) / float(cols)))
    return cols, rows


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


def _resolve_multi_record_definition_mode(
    mode: str,
) -> Literal["shared", "legacy"]:
    normalized = str(mode).strip().lower()
    if normalized not in _MULTI_RECORD_DEFINITION_MODES:
        raise ValidationError(
            "multi_record_definition_mode must be one of: shared, legacy"
        )
    return cast(Literal["shared", "legacy"], normalized)


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


def _resolve_multi_record_grid_gap_px(max_record_radius_px: float) -> float:
    """Return inter-record gap as a ratio of the largest record radius."""
    radius_px = float(max_record_radius_px)
    if radius_px > 0:
        return radius_px * _MULTI_RECORD_GRID_GAP_RATIO
    return float(_MULTI_RECORD_GRID_GAP_PX)


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
    definition_position: Literal["center", "top", "bottom"] = "center",
    track_specs: Sequence[str | TrackSpec] | None = None,
    _definition_profile: Literal["full", "record_summary", "shared_common"] = "full",
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
    normalized_definition_position = _resolve_definition_position(
        str(definition_position),
        argument_name="definition_position",
    )

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
        definition_position=normalized_definition_position,
        definition_profile=_definition_profile,
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
    definition_position: Literal["center", "top", "bottom"] = "center",
    multi_record_definition_mode: Literal["shared", "legacy"] = "shared",
    shared_definition_position: Literal["center", "top", "bottom"] = "bottom",
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
    normalized_definition_position = _resolve_definition_position(
        str(definition_position),
        argument_name="definition_position",
    )
    normalized_multi_record_definition_mode = _resolve_multi_record_definition_mode(
        str(multi_record_definition_mode)
    )
    normalized_shared_definition_position = _resolve_definition_position(
        str(shared_definition_position),
        argument_name="shared_definition_position",
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
            definition_position=normalized_definition_position,
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
    use_shared_definition = normalized_multi_record_definition_mode == "shared"
    record_definition_profile: Literal["full", "record_summary"] = (
        "record_summary" if use_shared_definition else "full"
    )
    record_definition_position: Literal["center", "top", "bottom"] = (
        "center" if use_shared_definition else normalized_definition_position
    )
    max_record_length = max(len(record.seq) for record in records)

    canvases: list[Drawing] = []
    widths: list[float] = []
    heights: list[float] = []
    record_radii_px: list[float] = []
    for record in records:
        record_scale = _resolve_multi_record_scale(
            len(record.seq),
            max_record_length,
            mode=normalized_multi_record_size_mode,
            min_radius_ratio=normalized_multi_record_min_radius_ratio,
        )
        scaled_cfg = _scale_circular_cfg(cfg, scale=record_scale)
        record_radii_px.append(float(scaled_cfg.canvas.circular.radius))
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
            definition_position=record_definition_position,
            track_specs=parsed_track_specs,
            _definition_profile=record_definition_profile,
            cfg=scaled_cfg,
        )
        canvases.append(sub_canvas)
        widths.append(_parse_svg_length_px(sub_canvas.attribs.get("width"), default=0.0))
        heights.append(_parse_svg_length_px(sub_canvas.attribs.get("height"), default=0.0))

    max_record_radius_px = max(
        [float(radius) for radius in record_radii_px if float(radius) > 0.0],
        default=float(cfg.canvas.circular.radius),
    )
    grid_gap_px = _resolve_multi_record_grid_gap_px(max_record_radius_px)

    cols, rows = _estimate_square_grid(len(canvases))
    row_heights: list[float] = [0.0] * rows
    row_record_indices: list[list[int]] = [[] for _ in range(rows)]
    record_sizes: list[tuple[float, float]] = []
    for record_index, (width_px, height_px) in enumerate(zip(widths, heights)):
        row = record_index // cols
        width_value = float(width_px)
        height_value = float(height_px)
        row_heights[row] = max(row_heights[row], height_value)
        row_record_indices[row].append(record_index)
        record_sizes.append((width_value, height_value))

    def _grid_offsets(lengths: list[float], gap: float) -> list[float]:
        offsets: list[float] = []
        cursor = 0.0
        for index, length in enumerate(lengths):
            offsets.append(cursor)
            cursor += float(length)
            if index < len(lengths) - 1:
                cursor += float(gap)
        return offsets

    row_offsets = _grid_offsets(row_heights, grid_gap_px)
    row_widths: list[float] = []
    row_total_widths: list[float] = []
    row_outer_left_paddings: list[float] = []
    row_outer_right_paddings: list[float] = []
    record_horizontal_insets: list[tuple[float, float]] = [
        _estimate_subcanvas_horizontal_insets(sub_canvas) for sub_canvas in canvases
    ]
    apply_row_margin_symmetry = (
        str(legend_effective).strip().lower() in {"none", "top", "bottom"}
    )
    for row_indices in row_record_indices:
        if not row_indices:
            row_widths.append(0.0)
            row_total_widths.append(0.0)
            row_outer_left_paddings.append(0.0)
            row_outer_right_paddings.append(0.0)
            continue
        row_width = sum(record_sizes[index][0] for index in row_indices)
        row_width += max(0, len(row_indices) - 1) * grid_gap_px
        row_widths.append(float(row_width))

        extra_left = 0.0
        extra_right = 0.0
        if apply_row_margin_symmetry:
            first_index = row_indices[0]
            last_index = row_indices[-1]
            row_left = max(0.0, float(record_horizontal_insets[first_index][0]))
            row_right = max(0.0, float(record_horizontal_insets[last_index][1]))
            target_margin = max(row_left, row_right)
            extra_left = max(0.0, target_margin - row_left)
            extra_right = max(0.0, target_margin - row_right)

        row_outer_left_paddings.append(float(extra_left))
        row_outer_right_paddings.append(float(extra_right))
        row_total_widths.append(float(row_width + extra_left + extra_right))

    grid_width = max(row_total_widths, default=0.0)
    grid_height = (
        (row_offsets[-1] + row_heights[-1]) if row_offsets else 0.0
    )
    record_offsets_x: dict[int, float] = {}
    for row, row_indices in enumerate(row_record_indices):
        if not row_indices:
            continue
        row_width = row_widths[row] if row < len(row_widths) else 0.0
        row_total_width = (
            row_total_widths[row] if row < len(row_total_widths) else row_width
        )
        row_outer_left_padding = (
            row_outer_left_paddings[row] if row < len(row_outer_left_paddings) else 0.0
        )
        cursor_x = max(0.0, (grid_width - row_total_width) * 0.5) + row_outer_left_padding
        for position, record_index in enumerate(row_indices):
            record_offsets_x[record_index] = cursor_x
            cursor_x += record_sizes[record_index][0]
            if position < len(row_indices) - 1:
                cursor_x += grid_gap_px

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
    shared_definition_group: Group | None = None
    shared_definition_local_bounds = (0.0, 0.0)

    if use_shared_definition:
        shared_canvas_config = CircularCanvasConfigurator(
            output_prefix=output_prefix,
            config_dict=config_dict,
            legend="none",
            gb_record=records[0],
            cfg=cfg,
        )
        shared_definition_group = DefinitionGroup(
            gb_record=records[0],
            canvas_config=shared_canvas_config,
            config_dict=config_dict,
            species=species,
            strain=strain,
            definition_profile="shared_common",
            definition_group_id="shared_definition",
            cfg=cfg,
        ).get_group()
        shared_definition_local_bounds = _group_local_vertical_bounds(shared_definition_group)

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
            legend_local_top = -0.5 * float(legend_config.color_rect_size)
            legend_local_bottom = legend_local_top + float(legend_config.legend_height)

            if legend_effective == "right":
                total_width = grid_width + (legend_config.legend_width * 1.1)
                legend_offset_x = grid_width + (legend_config.legend_width * 0.05)
                legend_offset_y = (total_height - legend_config.legend_height) / 2.0
            elif legend_effective == "left":
                total_width = grid_width + (legend_config.legend_width * 1.1)
                grid_origin_x = legend_config.legend_width * 1.1
                legend_offset_x = legend_config.legend_width * 0.05
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
        use_shared_definition
        and shared_definition_group is not None
        and normalized_shared_definition_position == "bottom"
    ):
        shared_min_y, shared_max_y = shared_definition_local_bounds
        shared_height = max(0.0, float(shared_max_y) - float(shared_min_y))

        records_bottom = float(grid_origin_y) + float(grid_height)
        anchor_bottom = records_bottom
        if legend_effective == "bottom" and legend_config is not None:
            legend_bottom = legend_offset_y + legend_local_bottom
            anchor_bottom = max(anchor_bottom, float(legend_bottom))

        required_height = (
            anchor_bottom
            + _MULTI_RECORD_LEGEND_SHARED_GAP_PX
            + shared_height
            + _MULTI_RECORD_SHARED_BOTTOM_MARGIN_PX
        )
        total_height = max(total_height, required_height)

    merged_canvas = Drawing(
        filename=f"{output_prefix}.svg",
        size=(f"{total_width}px", f"{total_height}px"),
        viewBox=f"0 0 {total_width} {total_height}",
        debug=False,
    )

    for record_index, sub_canvas in enumerate(canvases):
        row = record_index // cols
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
        merged_canvas.add(record_group)

    if shared_definition_group is not None:
        shared_definition_group = place_definition_group_on_size(
            shared_definition_group,
            canvas_width=float(total_width),
            canvas_height=float(total_height),
            position=normalized_shared_definition_position,
        )
        merged_canvas.add(shared_definition_group)

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
        definition_position=(output.definition_position if output else "center"),
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


