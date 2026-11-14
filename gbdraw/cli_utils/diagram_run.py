"""What ``gbdraw circular`` and ``gbdraw linear`` derive alike from their parsed options.

Each mode's ``run_*_from_namespace`` reads the options it alone has; the config overrides, the request fields, the
output request and the log of a render that both modes take the same way from the same options come from here.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Mapping, Sequence, TypedDict

from gbdraw.analysis.feature_density import FeatureDensityTrack
from gbdraw.annotations.models import AnnotationOptions
from gbdraw.api.input_summary import InputWarningKey, log_input_summary
from gbdraw.api.options import ColorOptions, DepthTrackInput
from gbdraw.api.request_render import RequestRenderResult
from gbdraw.api.requests import RenderOutputRequest
from gbdraw.cli_utils.common import cli_font_family, feature_density_tracks_from_args, handle_output_formats
from gbdraw.features.overrides import log_feature_identity_notices
from gbdraw.features.shapes import FeatureShape
from gbdraw.labels.label_style import LabelStyleMode
from gbdraw.render.export import parse_formats

logger = logging.getLogger()


def cli_output_formats(format_option: str) -> list[str]:
    """The formats ``-f/--format`` names that this run writes (``handle_output_formats``)."""

    return handle_output_formats(parse_formats(format_option))


def render_output_request(
    output_prefix: str,
    formats: list[str],
    *,
    overwrite: bool,
    embed_fonts: bool,
    dpi: int | None,
) -> RenderOutputRequest:
    """Split a CLI path prefix into the typed output directory and basename."""

    prefix_path = Path(str(output_prefix))
    output_directory = (
        prefix_path.parent
        if prefix_path.parent != Path(".")
        else None
    )
    return RenderOutputRequest(
        output_prefix=prefix_path.name,
        output_directory=output_directory,
        formats=tuple(formats),
        overwrite=overwrite,
        embed_fonts=embed_fonts,
        dpi=dpi,
    )


def paired_overrides(values: Mapping[str, object]) -> dict[str, object]:
    """Each config leaf that has a value for short and for long records set to one option: ``<path>.short``, then
    ``<path>.long``, in the order of ``values``."""

    return {f"{path}.{length}": value for path, value in values.items() for length in ("short", "long")}


def style_overrides(args: argparse.Namespace, *, mode: LabelStyleMode) -> dict[str, object]:
    """The font family (registering each ``--font_file``), the feature stroke color, the arrow geometry, and the axis
    stroke color of ``mode``."""

    return {
        "objects.text.font_family": cli_font_family(args),
        "objects.features.block_stroke_color": args.block_stroke_color,
        "objects.features.arrow_geometry.head_length_ratio": args.arrow_head_length_ratio,
        "objects.features.arrow_geometry.shaft_width_ratio": args.arrow_shaft_width_ratio,
        f"objects.axis.{mode}.stroke_color": args.axis_stroke_color,
    }


def label_style_overrides(args: argparse.Namespace, *, mode: LabelStyleMode) -> dict[str, object]:
    """``--label_rendering``, ``--label_line``, ``--label_protection`` and ``--italic_gene_labels`` of ``mode``."""

    return {
        "labels.rendering": args.label_rendering,
        f"labels.{mode}.line": args.label_line,
        f"labels.{mode}.protection": args.label_protection,
        # Only the italic style is an override, so a default run's request is unchanged.
        f"labels.{mode}.italic_gene_labels": True if args.italic_gene_labels else None,
    }


def label_filtering_override(
    filtering: Mapping[str, object], label_blacklist: str | Sequence[str] | None
) -> dict[str, object]:
    """``labels.filtering`` of the config, with the keywords of ``--label_blacklist`` (comma-separated) when given."""

    override = dict(filtering)
    if label_blacklist is not None:
        override["blacklist_keywords"] = (
            [keyword.strip() for keyword in label_blacklist.split(",") if keyword.strip()]
            if isinstance(label_blacklist, str)
            else label_blacklist
        )
    return override


def gc_content_overrides(args: argparse.Namespace, *, gc_content_mode: str | None) -> dict[str, object]:
    """``objects.gc_content``: the mode as the command resolves it, then the axis, its ticks and its labels."""

    large_tick_interval = args.gc_content_large_tick_interval
    return {
        "objects.gc_content.mode": gc_content_mode,
        "objects.gc_content.min_percent": args.gc_content_min_percent,
        "objects.gc_content.max_percent": args.gc_content_max_percent,
        "objects.gc_content.show_axis": args.gc_content_show_axis,
        "objects.gc_content.show_ticks": args.gc_content_show_ticks,
        "objects.gc_content.axis_labels": args.gc_content_axis_labels,
        # --gc_content_tick_interval is the alias the large interval falls back on.
        "objects.gc_content.large_tick_interval": (
            large_tick_interval if large_tick_interval is not None else args.gc_content_tick_interval
        ),
        "objects.gc_content.small_tick_interval": args.gc_content_small_tick_interval,
        "objects.gc_content.tick_font_size": args.gc_content_tick_font_size,
    }


def depth_overrides(args: argparse.Namespace) -> dict[str, object]:
    """``objects.depth``: the axis range, normalization, axis, ticks, labels and shared axis (each mode sets the fill
    color and the track size before them)."""

    return {
        "objects.depth.min_depth": args.depth_min,
        "objects.depth.max_depth": args.depth_max,
        "objects.depth.normalize": args.depth_normalize,
        "objects.depth.show_axis": args.depth_show_axis,
        "objects.depth.show_ticks": args.depth_show_ticks,
        "objects.depth.axis_labels": args.depth_axis_labels,
        "objects.depth.large_tick_interval": args.depth_large_tick_interval,
        "objects.depth.small_tick_interval": args.depth_small_tick_interval,
        "objects.depth.tick_font_size": args.depth_tick_font_size,
        "objects.depth.share_axis": bool(args.share_depth_axis),
    }


def stroke_width_overrides(args: argparse.Namespace, *, mode: LabelStyleMode) -> dict[str, object]:
    """``--block_stroke_width``, ``--axis_stroke_width`` (of ``mode``) and ``--line_stroke_width``."""

    return paired_overrides(
        {
            "objects.features.block_stroke_width": args.block_stroke_width,
            f"objects.axis.{mode}.stroke_width": args.axis_stroke_width,
            "objects.features.line_stroke_width": args.line_stroke_width,
        }
    )


def legend_size_overrides(args: argparse.Namespace) -> dict[str, object]:
    """``--legend_box_size`` and ``--legend_font_size``."""

    return paired_overrides(
        {
            "objects.legends.color_rect_size": args.legend_box_size,
            "objects.legends.font_size": args.legend_font_size,
        }
    )


class ModeOptionFields(TypedDict):
    """The fields of either mode's diagram options that both commands set alike from their options."""

    colors: ColorOptions
    annotations: AnnotationOptions | None
    selected_features_set: tuple[str, ...]
    feature_visibility_table_file: str | None
    feature_placement_table_file: str | None
    feature_override_table_file: str | None
    label_whitelist_file: str | None
    qualifier_priority_file: str | None
    label_override_file: str | None
    feature_shapes: dict[str, FeatureShape] | None
    dinucleotide: str
    window: int | None
    step: int | None
    depth_window: int | None
    depth_step: int | None
    depth_tracks: tuple[DepthTrackInput, ...] | None
    feature_density_tracks: tuple[FeatureDensityTrack, ...] | None
    plot_title: str | None
    plot_title_font_size: float | None
    evalue: float
    bitscore: float
    identity: float
    alignment_length: int


def mode_option_fields(
    args: argparse.Namespace,
    *,
    feature_shapes: dict[str, FeatureShape],
    depth_tracks: tuple[DepthTrackInput, ...] | None,
) -> ModeOptionFields:
    """The shared option fields; ``feature_shapes`` and ``depth_tracks`` are read earlier, when each mode reads them.

    The feature density tracks are read here (``--feature_density_table`` is a file).
    """

    plot_title = str(args.plot_title or "").strip()
    return {
        "colors": ColorOptions(
            color_table_file=args.table or None,
            default_colors_palette=args.palette,
            default_colors_file=args.default_colors or None,
        ),
        "annotations": AnnotationOptions(table_file=args.annotation_table) if args.annotation_table else None,
        "selected_features_set": tuple(args.features.split(",")),
        "feature_visibility_table_file": args.feature_table or None,
        "feature_placement_table_file": args.feature_placement_table,
        "feature_override_table_file": args.feature_override_table,
        "label_whitelist_file": args.label_whitelist or None,
        "qualifier_priority_file": args.qualifier_priority or None,
        "label_override_file": args.label_table or None,
        "feature_shapes": feature_shapes or None,
        "dinucleotide": args.nt.upper(),
        "window": args.window,
        "step": args.step,
        "depth_window": args.depth_window,
        "depth_step": args.depth_step,
        "depth_tracks": depth_tracks,
        "feature_density_tracks": feature_density_tracks_from_args(args),
        "plot_title": plot_title or None,
        "plot_title_font_size": args.plot_title_font_size,
        "evalue": args.evalue,
        "bitscore": args.bitscore,
        "identity": args.identity,
        "alignment_length": args.alignment_length,
    }


def log_render_notices(
    rendered_items: Sequence[RequestRenderResult], *, warned: set[InputWarningKey] | None = None
) -> None:
    """Log each Result's annotation warnings, feature notices, input summary and, for a wrapped Linear figure, its wrap.

    ``warned``: the input warnings logged so far; with it, a warning of an input the Results of a batch share (one
    comparison file every Result reads) is logged once. Without it, every warning is logged.
    """

    for rendered_item in rendered_items:
        for warning in rendered_item.annotation_warnings:
            logger.warning("%s: %s/%s record #%s (%s): %s", warning.code,
                warning.set_id, warning.annotation_id, warning.record_index + 1,
                warning.record_id, warning.message)
        log_feature_identity_notices(rendered_item.feature_identity_notices)
        log_input_summary(rendered_item.input_summary, warned=warned)
        if rendered_item.linear_metadata is not None and rendered_item.linear_metadata.wrap is not None:
            logger.info("%s", rendered_item.linear_metadata.wrap.message)
