"""Diagram assembly helpers (public API layer).

This module provides convenience functions to build diagrams from in-memory objects
(e.g., already-parsed SeqRecords) without going through the CLI argument parsing.

The functions here return an `svgwrite.Drawing` (SVG canvas). Saving/conversion is
handled separately (see `gbdraw.render.export.save_figure`).
"""

from __future__ import annotations

import logging
from dataclasses import replace
from typing import Optional, Sequence, Mapping

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]

from gbdraw.analysis.skew import skew_df  # type: ignore[reportMissingImports]
from gbdraw.api.config import apply_config_overrides  # type: ignore[reportMissingImports]
from gbdraw.api.options import DiagramOptions  # type: ignore[reportMissingImports]
from gbdraw.canvas import CircularCanvasConfigurator, LinearCanvasConfigurator  # type: ignore[reportMissingImports]
from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from gbdraw.config.modify import modify_config_dict  # type: ignore[reportMissingImports]
from gbdraw.config.toml import load_config_toml  # type: ignore[reportMissingImports]
from gbdraw.io.colors import load_default_colors, read_color_table  # type: ignore[reportMissingImports]
from gbdraw.configurators import (  # type: ignore[reportMissingImports]
    BlastMatchConfigurator,
    FeatureDrawingConfigurator,
    GcContentConfigurator,
    GcSkewConfigurator,
    LegendDrawingConfigurator,
)
from gbdraw.core.sequence import create_dict_for_sequence_lengths  # type: ignore[reportMissingImports]
from gbdraw.diagrams.circular import assemble_circular_diagram  # type: ignore[reportMissingImports]
from gbdraw.diagrams.linear import assemble_linear_diagram  # type: ignore[reportMissingImports]
from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]
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
    "assemble_linear_diagram_from_records",
    "build_circular_diagram",
    "build_linear_diagram",
]


