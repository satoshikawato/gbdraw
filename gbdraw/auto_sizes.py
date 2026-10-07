"""Length- and mode-dependent Auto sizes.

A setting left at Auto has no explicit value: its value follows the drawn record
lengths and the diagram mode. This module owns that knowledge for every
renderer and surface: the size class (short or long), the sliding-window and
scale-tick tiers, the derived Depth window, the scalar-axis tick font, the
circular track widths, and the Auto value of every such setting of one drawing.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import reduce
from typing import Literal

from gbdraw.config.models import GbdrawConfig
from gbdraw.mode_profiles import DiagramMode


SizeClass = Literal["short", "long"]
AutoValue = float | int

# (at most bp, large tick, small tick); longer records use the last pair.
_CIRCULAR_TICK_TIERS: tuple[tuple[int, int, int], ...] = (
    (30_000, 1_000, 100),
    (50_000, 5_000, 1_000),
    (150_000, 10_000, 1_000),
    (1_000_000, 50_000, 10_000),
    (10_000_000, 500_000, 100_000),
)
_CIRCULAR_TICKS_ABOVE_TIERS = (1_000_000, 200_000)
# (below bp, tick interval); longer records use the last interval.
_LINEAR_TICK_TIERS: tuple[tuple[int, int], ...] = (
    (2_000, 100),
    (20_000, 1_000),
    (50_000, 5_000),
    (150_000, 10_000),
    (250_000, 50_000),
    (1_000_000, 100_000),
    (2_000_000, 200_000),
    (5_000_000, 500_000),
)
_LINEAR_TICK_ABOVE_TIERS = 1_000_000
_SCALAR_TICK_FONT_RANGE_PX = (5.0, 8.0)
_SCALAR_TICK_FONT_THICKNESS_FRACTION: Mapping[DiagramMode, float] = {
    "circular": 0.22,
    "linear": 0.7,
}

# Short/long config pairs, by the mode that draws them.
_SHARED_LENGTH_PAIRS = (
    "objects.features.block_stroke_width",
    "objects.features.line_stroke_width",
    "objects.legends.color_rect_size",
    "objects.legends.font_size",
)
_MODE_LENGTH_PAIRS: Mapping[DiagramMode, tuple[str, ...]] = {
    "circular": (
        "labels.font_size",
        "labels.stroke_width",
        "objects.axis.circular.stroke_width",
    ),
    "linear": (
        "canvas.linear.arrow_length_parameter",
        "canvas.linear.default_cds_height",
        "objects.axis.linear.stroke_width",
        "objects.definition.linear.font_size",
        "objects.scale.font_size",
        "objects.scale.ruler_label_font_size",
    ),
}
# The settings that give the track thickness the tick fonts follow.
_TRACK_THICKNESS_SETTINGS: Mapping[DiagramMode, Mapping[str, str]] = {
    "circular": {"gc_content": "circularSlot.gc_content.width", "depth": "circularSlot.depth.width"},
    "linear": {"gc_content": "canvas.linear.default_gc_height", "depth": "canvas.linear.depth_height"},
}
_CIRCULAR_SLOT_RENDERERS = {
    "features": "features",
    "gc_content": "dinucleotide_content",
    "gc_skew": "dinucleotide_skew",
    "depth": "depth",
}


def determine_length_parameter(record_length: int, length_threshold: int) -> SizeClass:
    """Return the size class of one record length."""
    if record_length < length_threshold:
        return "short"
    return "long"


def sliding_window(length: int, cfg: GbdrawConfig) -> tuple[int, int]:
    """Return the Auto GC content and skew (window, step) for a record length."""
    tiers = cfg.objects.sliding_window
    if length < 1_000_000:
        window, step = tiers.default
    elif length < 10_000_000:
        window, step = tiers.up1m
    else:
        window, step = tiers.up10m
    return int(window), int(step)


def depth_sliding_window(window: int, step: int) -> tuple[int, int]:
    """Return the Auto Depth (window, step) from the GC (window, step).

    Depth samples ten times denser than GC, with at least a 100 bp window so the
    default track is not overly noisy.
    """
    return max(100, int(window) // 10), max(1, int(step) // 10)


def circular_tick_intervals(length: int) -> tuple[int, int]:
    """Return the Auto (large, small) tick intervals of a circular record."""
    for at_most, large, small in _CIRCULAR_TICK_TIERS:
        if length <= at_most:
            return large, small
    return _CIRCULAR_TICKS_ABOVE_TIERS


def linear_tick_interval(length: int) -> int:
    """Return the Auto scale tick interval for a linear length."""
    for below, interval in _LINEAR_TICK_TIERS:
        if length < below:
            return interval
    return _LINEAR_TICK_ABOVE_TIERS


def linear_label_font_size(record_lengths: Sequence[int], cfg: GbdrawConfig) -> float:
    """Return the Auto Linear label font: the largest size over the labelled records."""
    threshold = int(cfg.labels.length_threshold.linear)
    return max(
        cfg.labels.font_size.linear.for_length_param(
            determine_length_parameter(record_length, threshold)
        )
        for record_length in record_lengths
    )


def circular_track_width_px(
    renderer: str,
    *,
    radius: float,
    track_ratio: float,
    factors: Sequence[float],
) -> float:
    """Return the Auto radial width of a circular track drawn by ``renderer``.

    ``factors`` are the ``canvas.circular.track_ratio_factors`` of the record's
    size class.
    """
    base = float(radius) * float(track_ratio)
    if renderer in {"features", "sequence_conservation"}:
        return base * float(factors[0])
    if renderer == "depth":
        return base * float(factors[1]) * 0.5
    if renderer == "dinucleotide_skew":
        return base * float(factors[2])
    return base * float(factors[1])


def scalar_axis_tick_font_size(mode: DiagramMode, track_thickness_px: float) -> float:
    """Return the Auto tick font of a GC content or Depth axis for its track thickness."""
    low, high = _SCALAR_TICK_FONT_RANGE_PX
    return max(
        low,
        min(high, float(track_thickness_px) * _SCALAR_TICK_FONT_THICKNESS_FRACTION[mode]),
    )


@dataclass(frozen=True)
class DrawnExtent:
    """The records one drawing draws: its mode and lengths after crops, in drawing order."""

    mode: DiagramMode
    record_lengths: tuple[int, ...]

    def __post_init__(self) -> None:
        if self.mode not in _MODE_LENGTH_PAIRS:
            raise ValueError(f"Unknown diagram mode: {self.mode!r}")
        if not self.record_lengths or min(self.record_lengths) < 1:
            raise ValueError("A drawn extent needs at least one record of length >= 1")

    def size_class(self, cfg: GbdrawConfig) -> SizeClass:
        """Return the drawing's size class: that of its longest record."""
        threshold = int(getattr(cfg.labels.length_threshold, self.mode))
        return determine_length_parameter(max(self.record_lengths), threshold)


def auto_setting_values(
    extent: DrawnExtent,
    cfg: GbdrawConfig,
    *,
    track_thickness_px: Mapping[str, float] | None = None,
) -> dict[str, AutoValue]:
    """Return the Auto value of every length- or mode-dependent setting of a drawing.

    Keys are config override paths without ``.short``/``.long``, the diagram
    options ``window``, ``step``, ``depth_window`` and ``depth_step``, and
    ``circularSlot.<slot>.width`` for the Auto circular track widths. The Linear
    GC content and Depth heights are listed too: they do not change with the
    length, but the tick fonts follow them. Any other setting has the same Auto
    value at every length in that mode.

    ``cfg`` must hold no explicit size values (normally the packaged defaults).
    The size class, windows and tick intervals follow the longest record; in a
    Circular canvas the records then share the long styles. The GC content and
    Depth tick fonts follow the track thickness: the Auto thickness unless
    ``track_thickness_px`` gives the drawing's own (keys ``gc_content`` and
    ``depth``).
    """
    mode = extent.mode
    longest = max(extent.record_lengths)
    size = extent.size_class(cfg)
    values: dict[str, AutoValue] = {
        path: float(getattr(reduce(getattr, path.split("."), cfg), size))
        for path in (*_SHARED_LENGTH_PAIRS, *_MODE_LENGTH_PAIRS[mode])
    }
    window, step = sliding_window(longest, cfg)
    depth_window, depth_step = depth_sliding_window(window, step)
    values.update(window=window, step=step, depth_window=depth_window, depth_step=depth_step)
    if mode == "circular":
        circular = cfg.canvas.circular
        for slot_id, renderer in _CIRCULAR_SLOT_RENDERERS.items():
            values[f"circularSlot.{slot_id}.width"] = circular_track_width_px(
                renderer,
                radius=circular.radius,
                track_ratio=circular.track_ratio,
                factors=circular.track_ratio_factors[size],
            )
        values["objects.scale.interval"] = circular_tick_intervals(longest)[0]
    else:
        values["labels.font_size.linear"] = linear_label_font_size(extent.record_lengths, cfg)
        values["objects.scale.interval"] = linear_tick_interval(longest)
        values["canvas.linear.default_gc_height"] = float(cfg.canvas.linear.default_gc_height)
        values["canvas.linear.depth_height"] = float(cfg.canvas.linear.depth_height)
    thickness = {
        track: float(values[setting]) for track, setting in _TRACK_THICKNESS_SETTINGS[mode].items()
    }
    thickness.update(track_thickness_px or {})
    values["objects.gc_content.tick_font_size"] = scalar_axis_tick_font_size(mode, thickness["gc_content"])
    values["objects.depth.tick_font_size"] = scalar_axis_tick_font_size(mode, thickness["depth"])
    return values


__all__ = [
    "AutoValue",
    "DrawnExtent",
    "SizeClass",
    "auto_setting_values",
    "circular_tick_intervals",
    "circular_track_width_px",
    "depth_sliding_window",
    "determine_length_parameter",
    "linear_label_font_size",
    "linear_tick_interval",
    "scalar_axis_tick_font_size",
    "sliding_window",
]
