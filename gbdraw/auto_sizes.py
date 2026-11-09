"""Length- and mode-dependent Auto sizes.

A setting left at Auto has no explicit value: its value follows the drawn record
lengths and the diagram mode. This module owns that knowledge for every
renderer and surface: the size class (short or long), the sliding-window and
scale-tick tiers, the derived Depth window, the scalar-axis tick font, the
circular track widths, and the Auto value of every such setting of one drawing.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from functools import reduce
from typing import Literal, Protocol

from gbdraw.config.models import GbdrawConfig
from gbdraw.config.models.labels import LINEAR_LABEL_SCOPE_RECORDS, LinearLabelScope
from gbdraw.mode_profiles import DiagramMode


SizeClass = Literal["short", "long"]
AutoValue = float | int

# (at most bp, large tick, small tick); longer records use the last pair. Above 10 Mbp the large ticks step
# 1-2-5 Mbp, and each row ends at the longest record whose neighbouring large tick labels, at their drawn
# widths, keep a quarter em (labels.label_style.ROW_LABEL_GAP_EM) between them in the default figure (OV-2747,
# OV-2748; tests/test_circular_ticks.py measures every bound again).
_CIRCULAR_TICK_TIERS: tuple[tuple[int, int, int], ...] = (
    (30_000, 1_000, 100),
    (50_000, 5_000, 1_000),
    (150_000, 10_000, 1_000),
    (1_000_000, 50_000, 10_000),
    (10_000_000, 500_000, 100_000),
    (33_839_331, 1_000_000, 200_000),
    (67_678_663, 2_000_000, 400_000),
    (150_284_490, 5_000_000, 1_000_000),
    (300_568_981, 10_000_000, 2_000_000),
    (601_137_962, 20_000_000, 4_000_000),
    (1_351_751_222, 50_000_000, 10_000_000),
    (2_703_502_444, 100_000_000, 20_000_000),
    (5_407_004_889, 200_000_000, 40_000_000),
)
# Its labels keep the gap on records up to 12.2 Gbp.
_CIRCULAR_TICKS_ABOVE_TIERS = (500_000_000, 100_000_000)
# (below bp, tick interval); longer records use the last interval. The scale bar takes its length from them
# too. Above 5 Mbp the interval steps 1-2-5 Mbp, and each row ends where neighbouring ruler labels, at their
# drawn widths and the default width, would come closer than a quarter em (labels.label_style.ROW_LABEL_GAP_EM;
# OV-2747, OV-2748; tests/test_linear_ruler_ticks.py measures every bound again).
_LINEAR_TICK_TIERS: tuple[tuple[int, int], ...] = (
    (2_000, 100),
    (20_000, 1_000),
    (50_000, 5_000),
    (150_000, 10_000),
    (250_000, 50_000),
    (1_000_000, 100_000),
    (2_000_000, 200_000),
    (5_000_000, 500_000),
    (46_484_180, 1_000_000),
    (92_968_360, 2_000_000),
    (201_210_407, 5_000_000),
    (402_420_813, 10_000_000),
    (804_841_626, 20_000_000),
    (1_773_897_378, 50_000_000),
    (3_547_794_755, 100_000_000),
    (7_095_589_510, 200_000_000),
)
# Its ruler labels keep the gap on records up to 15.8 Gbp.
_LINEAR_TICK_ABOVE_TIERS = 500_000_000
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
_WINDOW_SETTINGS = ("window", "step", "depth_window", "depth_step")
# The settings whose Auto value follows the drawn extent but is no length pair.
_MODE_EXTENT_SETTINGS: Mapping[DiagramMode, tuple[str, ...]] = {
    "circular": (
        *(f"circularSlot.{slot_id}.width" for slot_id in _CIRCULAR_SLOT_RENDERERS),
        "objects.scale.interval",
    ),
    "linear": (
        "labels.font_size.linear",
        "objects.scale.interval",
        "canvas.linear.default_gc_height",
        "canvas.linear.depth_height",
    ),
}
# The tick font of each scalar track, keyed as _TRACK_THICKNESS_SETTINGS.
_TICK_FONT_SETTINGS = {
    "gc_content": "objects.gc_content.tick_font_size",
    "depth": "objects.depth.tick_font_size",
}
ADAPTABLE_SETTINGS: Mapping[DiagramMode, tuple[str, ...]] = {
    mode: (
        *_SHARED_LENGTH_PAIRS,
        *pairs,
        *_WINDOW_SETTINGS,
        *_MODE_EXTENT_SETTINGS[mode],
        *_TICK_FONT_SETTINGS.values(),
    )
    for mode, pairs in _MODE_LENGTH_PAIRS.items()
}
"""The settings of each mode that can return to Auto, keyed as ``auto_setting_values`` keys them."""


def determine_length_parameter(record_length: int, length_threshold: int) -> SizeClass:
    """Return the size class of one record length."""
    if record_length < length_threshold:
        return "short"
    return "long"


# The Auto GC window is at most this fraction of the drawn length, and the Auto step at most this fraction of the
# window, so a record shorter than the tier's ten windows keeps about a hundred samples (OV-2864).
_SLIDING_WINDOW_DIVISOR = 10
# The Auto Depth window and step are the GC ones divided by this, the window at least _DEPTH_MIN_WINDOW bp.
_DEPTH_WINDOW_DIVISOR = 10
_DEPTH_MIN_WINDOW = 100


def sliding_window(length: int, cfg: GbdrawConfig) -> tuple[int, int]:
    """Return the Auto GC content and skew (window, step) for a drawn length: the tier's values, the window at most
    a tenth of the length and the step at most a tenth of the window, each at least 1 bp."""
    tier_window, tier_step = next(
        ((window, step) for below, window, step in sliding_window_tiers(cfg) if length < below),
        (int(cfg.objects.sliding_window.up10m[0]), int(cfg.objects.sliding_window.up10m[1])),
    )
    window = max(1, min(tier_window, int(length) // _SLIDING_WINDOW_DIVISOR))
    return window, max(1, min(tier_step, window // _SLIDING_WINDOW_DIVISOR))


def sliding_window_tiers(cfg: GbdrawConfig) -> tuple[tuple[int, int, int], ...]:
    """Return the Auto GC (below bp, window, step) tiers; longer records use ``up10m``."""
    tiers = cfg.objects.sliding_window
    return (
        (1_000_000, int(tiers.default[0]), int(tiers.default[1])),
        (10_000_000, int(tiers.up1m[0]), int(tiers.up1m[1])),
    )


def tick_and_window_tiers(cfg: GbdrawConfig) -> dict[str, object]:
    """Return the length tiers of the Auto ticks and sliding window, for the Web hints.

    Each ``tiers`` row is checked in order against the record length (circular
    ticks: at most; the others: below); a longer record takes ``above``.
    """
    up10m = cfg.objects.sliding_window.up10m
    return {
        "slidingWindow": {
            "tiers": [list(row) for row in sliding_window_tiers(cfg)],
            "above": [int(up10m[0]), int(up10m[1])],
            "divisor": _SLIDING_WINDOW_DIVISOR,
        },
        "depthWindow": {"divisor": _DEPTH_WINDOW_DIVISOR, "minWindow": _DEPTH_MIN_WINDOW},
        "circularTicks": {
            "tiers": [list(row) for row in _CIRCULAR_TICK_TIERS],
            "above": list(_CIRCULAR_TICKS_ABOVE_TIERS),
        },
        "linearTicks": {
            "tiers": [list(row) for row in _LINEAR_TICK_TIERS],
            "above": _LINEAR_TICK_ABOVE_TIERS,
        },
    }


def depth_sliding_window(window: int, step: int) -> tuple[int, int]:
    """Return the Auto Depth (window, step) from the GC (window, step).

    Depth samples ten times denser than GC, with at least a 100 bp window so the
    default track is not overly noisy.
    """
    return max(_DEPTH_MIN_WINDOW, int(window) // _DEPTH_WINDOW_DIVISOR), max(1, int(step) // _DEPTH_WINDOW_DIVISOR)


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


def linear_label_size_records(label_scope: LinearLabelScope) -> Literal["all", "first"]:
    """Return which records set the Auto Linear label font.

    They are the records in the label display scope. A scope that selects no
    record leaves the labels shown by per-feature overrides, which take the size
    of all records.
    """
    selected = LINEAR_LABEL_SCOPE_RECORDS[label_scope]
    return "all" if selected == "none" else selected


def linear_label_size_lengths(record_lengths: Sequence[int], label_scope: LinearLabelScope) -> Sequence[int]:
    """Return the lengths of the records that set the Auto Linear label font, in drawing order."""
    return record_lengths[:1] if linear_label_size_records(label_scope) == "first" else record_lengths


def linear_label_font_size(record_lengths: Sequence[int], cfg: GbdrawConfig) -> float:
    """Return the Auto Linear label font: that of the size class of the longest record in the label scope."""
    return cfg.labels.font_size.linear.for_length_param(
        determine_length_parameter(max(record_lengths), int(cfg.labels.length_threshold.linear))
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


def scalar_axis_tick_font_rule() -> dict[str, object]:
    """The rule of ``scalar_axis_tick_font_size`` as data, for the Web's Auto notes: the font is the track
    thickness times the mode's fraction, kept within the range."""
    return {
        "rangePx": list(_SCALAR_TICK_FONT_RANGE_PX),
        "thicknessFraction": dict(_SCALAR_TICK_FONT_THICKNESS_FRACTION),
    }


class TickFontSource(Protocol):
    """A scalar axis's tick font size: set, or none for the Auto size of its track."""

    @property
    def tick_font_size(self) -> float | None: ...


class TickedAxisSource(TickFontSource, Protocol):
    """A scalar axis that may draw ticks: whether it draws its axis and ticks, and their font size."""

    @property
    def show_axis(self) -> bool: ...

    @property
    def show_ticks(self) -> bool: ...


def depth_axis_tick_font_size_px(axis: TickFontSource, track_width_px: float) -> float:
    """The tick font of a Circular scalar axis on a track ``track_width_px`` wide: its own size, else the Auto
    size of that width."""

    if axis.tick_font_size is not None:
        return max(0.0, float(axis.tick_font_size))
    return scalar_axis_tick_font_size("circular", track_width_px)


def gc_axis_tick_font_size_of_depth(
    gc: TickFontSource, depth: TickedAxisSource | None, depth_track_width_px: float | None
) -> float | None:
    """The tick font a Percent GC axis without its own size takes from the depth axis beside it, so the two axes
    read alike: the depth axis's size on its track ``depth_track_width_px`` wide (as requested, before any
    narrowing). None when the GC axis sets its own size or no depth axis with ticks is drawn; the GC axis then
    keeps its own size or the Auto size of its own track. The radial layout keeps the GC axis's room in this font,
    the drawing draws it in it, and ``auto_setting_values`` gives it as the GC axis's Auto tick font."""

    if (
        gc.tick_font_size is not None
        or depth is None
        or not (depth.show_axis and depth.show_ticks)
        or depth_track_width_px is None
        or depth_track_width_px <= 0
    ):
        return None
    return depth_axis_tick_font_size_px(depth, depth_track_width_px)


def gc_axis_tick_font_size_of_depth_rows(
    gc: TickFontSource, depth_rows: Iterable[tuple[TickedAxisSource | None, float]]
) -> float | None:
    """The tick font a Percent GC axis ``gc`` takes from the depth rows beside it, each the axis of the depth track
    it draws (None for a track the record does not draw) and its width as requested: that of the first row, in stack
    order, that draws a track (``gc_axis_tick_font_size_of_depth``). The one reading of this rule: the drawing
    (``CircularScalarAxes.gc_tick_font_size``) and the Auto values (``auto_setting_values``) both call it."""

    for axis, width_px in depth_rows:
        if axis is not None and width_px > 0:
            return gc_axis_tick_font_size_of_depth(gc, axis, width_px)
    return None


@dataclass(frozen=True)
class DrawnDepthRow:
    """A depth row of a Circular drawing: the depth track it draws (``track_index``) and that track's axis, in its
    own settings (its tick font among them), and its width.

    ``width_px`` is a width the row states that no size setting holds (a share of the radius, or a width in px of
    a row other than the Depth row); None for the Auto width, or, when ``depth_width``, the Depth width. A
    ``depth_width`` row is the Depth row: its width is the setting ``circularSlot.depth.width``, which a drawing
    keeps by hand or returns to Auto.
    """

    axis: TickedAxisSource
    track_index: int
    width_px: float | None = None
    depth_width: bool = False


@dataclass(frozen=True)
class _DepthAxisInFont:
    """A depth row's axis drawn in the drawing's own Depth tick font."""

    show_axis: bool
    show_ticks: bool
    tick_font_size: float | None


@dataclass(frozen=True)
class DrawnExtent:
    """The records one drawing draws: its mode and lengths after crops, in drawing order.

    ``wrap_row_length`` is a wrapped drawing's row length: sizes and ticks follow what one row shows, the
    windows the whole record. ``depth_rows`` are the depth rows a Circular drawing draws, in stack order.
    """

    mode: DiagramMode
    record_lengths: tuple[int, ...]
    wrap_row_length: int | None = None
    depth_rows: tuple[DrawnDepthRow, ...] = ()

    def __post_init__(self) -> None:
        if self.mode not in _MODE_LENGTH_PAIRS:
            raise ValueError(f"Unknown diagram mode: {self.mode!r}")
        if not self.record_lengths or min(self.record_lengths) < 1:
            raise ValueError("A drawn extent needs at least one record of length >= 1")

    def size_class(self, cfg: GbdrawConfig) -> SizeClass:
        """Return the drawing's size class: that of its longest record."""
        threshold = int(getattr(cfg.labels.length_threshold, self.mode))
        return determine_length_parameter(self.drawn_axis_length, threshold)

    @property
    def drawn_axis_length(self) -> int:
        """The bases on the longest drawn axis."""
        longest = max(self.record_lengths)
        return longest if self.wrap_row_length is None else min(longest, self.wrap_row_length)


def auto_setting_values(
    extent: DrawnExtent,
    cfg: GbdrawConfig,
    *,
    track_thickness_px: Mapping[str, float] | None = None,
    depth_tick_font_size: float | None = None,
) -> dict[str, AutoValue]:
    """Return the Auto value of every length- or mode-dependent setting of a drawing.

    Keys are config override paths without ``.short``/``.long``, the diagram
    options ``window``, ``step``, ``depth_window`` and ``depth_step``, and
    ``circularSlot.<slot>.width`` for the Auto circular track widths. The Linear
    GC content and Depth heights are listed too: they do not change with the
    length, but the tick fonts follow them. Any other setting has the same Auto
    value at every length in that mode.

    ``cfg`` must hold no explicit size values (normally the packaged defaults).
    The size class and tick intervals follow the longest drawn axis and the windows the
    longest record; in a
    Circular canvas the records then share the long styles. The GC content and
    Depth tick fonts follow the track thickness: the Auto thickness unless
    ``track_thickness_px`` gives the drawing's own (keys ``gc_content`` and
    ``depth``); a Circular drawing's Depth font, that of the first depth
    row of series 0 (``extent.depth_rows``). The drawing gives a Percent GC
    axis the tick font of its first depth row
    (``gc_axis_tick_font_size_of_depth_rows``): the font of its track, else
    the drawing's own Depth tick font (``depth_tick_font_size``), else the
    Auto font of its width. ``cfg`` then holds the drawing's other settings
    (the GC content mode among them).
    """
    mode = extent.mode
    longest = max(extent.record_lengths)
    size = extent.size_class(cfg)
    values: dict[str, AutoValue] = {
        path: float(getattr(reduce(getattr, path.split("."), cfg), size))
        for path in (*_SHARED_LENGTH_PAIRS, *_MODE_LENGTH_PAIRS[mode])
    }
    window, step = sliding_window(longest, cfg)
    values.update(zip(_WINDOW_SETTINGS, (window, step, *depth_sliding_window(window, step))))
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
        values["labels.font_size.linear"] = linear_label_font_size((extent.drawn_axis_length,), cfg)
        values["objects.scale.interval"] = linear_tick_interval(extent.drawn_axis_length)
        values["canvas.linear.default_gc_height"] = float(cfg.canvas.linear.default_gc_height)
        values["canvas.linear.depth_height"] = float(cfg.canvas.linear.depth_height)
    thickness = {
        track: float(values[setting]) for track, setting in _TRACK_THICKNESS_SETTINGS[mode].items()
    }
    thickness.update(track_thickness_px or {})
    auto_depth_width = float(values[_TRACK_THICKNESS_SETTINGS[mode]["depth"]])
    rows = tuple(
        (
            row,
            row.width_px if row.width_px is not None else thickness["depth"] if row.depth_width else auto_depth_width,
        )
        for row in extent.depth_rows
    )
    # The Depth font's Auto value is that of the first row of the first series, the one the font sets.
    thickness["depth"] = next((width for row, width in rows if row.track_index == 0), thickness["depth"])
    for track, setting in _TICK_FONT_SETTINGS.items():
        values[setting] = scalar_axis_tick_font_size(mode, thickness[track])
    gc_content = cfg.objects.gc_content
    if mode == "circular" and gc_content.mode == "percent":

        def drawn_axis(axis: TickedAxisSource) -> TickedAxisSource:
            if axis.tick_font_size is not None or depth_tick_font_size is None:
                return axis
            return _DepthAxisInFont(axis.show_axis, axis.show_ticks, depth_tick_font_size)

        of_depth = gc_axis_tick_font_size_of_depth_rows(
            gc_content, ((drawn_axis(row.axis), width) for row, width in rows)
        )
        if of_depth is not None:
            values[_TICK_FONT_SETTINGS["gc_content"]] = of_depth
    return values


@dataclass(frozen=True)
class SettingReset:
    """A hand-set size that returns to Auto because its Auto value differs.

    An Auto value of ``None`` means that drawing's mode does not draw the setting.
    """

    setting: str
    explicit: AutoValue
    source_auto: AutoValue | None
    target_auto: AutoValue | None


@dataclass(frozen=True)
class SettingAdaptation:
    """The hand-set sizes a new drawing keeps, and those that return to Auto."""

    kept: Mapping[str, AutoValue]
    reset: tuple[SettingReset, ...]


def adapt_explicit_settings(
    explicit: Mapping[str, AutoValue],
    *,
    source: DrawnExtent,
    target: DrawnExtent,
    cfg: GbdrawConfig,
) -> SettingAdaptation:
    """Keep a hand-set size only where its Auto value is the same in ``target``.

    ``explicit`` is keyed as ``auto_setting_values`` is. A setting whose Auto
    value changes with the size class, a tier or the mode returns to Auto and
    is listed in ``reset``; the other sizes are listed in ``kept``. A value
    that neither drawing lists is not a size and is in neither. The tick fonts
    are compared at the track thickness each drawing ends with, and its depth
    rows at the Depth width it ends with.
    """

    def thickness(mode: DiagramMode, values: Mapping[str, AutoValue]) -> dict[str, float]:
        return {
            track: float(values[setting])
            for track, setting in _TRACK_THICKNESS_SETTINGS[mode].items()
            if setting in values
        }

    depth_font = _TICK_FONT_SETTINGS["depth"]
    source_auto = auto_setting_values(
        source,
        cfg,
        track_thickness_px=thickness(source.mode, explicit),
        depth_tick_font_size=explicit.get(depth_font),
    )
    target_auto = auto_setting_values(target, cfg)
    unchanged = {
        setting: value
        for setting, value in explicit.items()
        if source_auto.get(setting) == target_auto.get(setting)
    }
    # The thickness settings never depend on a thickness, so this pass only
    # moves the tick fonts to the thickness and the Depth font the target keeps.
    target_auto = auto_setting_values(
        target,
        cfg,
        track_thickness_px=thickness(target.mode, unchanged),
        depth_tick_font_size=unchanged.get(depth_font),
    )
    kept: dict[str, AutoValue] = {}
    reset: list[SettingReset] = []
    for setting, value in explicit.items():
        if setting not in source_auto and setting not in target_auto:
            continue
        if source_auto.get(setting) == target_auto.get(setting):
            kept[setting] = value
        else:
            reset.append(
                SettingReset(setting, value, source_auto.get(setting), target_auto.get(setting))
            )
    return SettingAdaptation(kept=kept, reset=tuple(reset))


__all__ = [
    "ADAPTABLE_SETTINGS",
    "AutoValue",
    "DrawnDepthRow",
    "DrawnExtent",
    "SettingAdaptation",
    "SettingReset",
    "SizeClass",
    "adapt_explicit_settings",
    "auto_setting_values",
    "circular_tick_intervals",
    "circular_track_width_px",
    "TickFontSource",
    "TickedAxisSource",
    "depth_axis_tick_font_size_px",
    "depth_sliding_window",
    "determine_length_parameter",
    "gc_axis_tick_font_size_of_depth",
    "gc_axis_tick_font_size_of_depth_rows",
    "linear_label_font_size",
    "linear_label_size_lengths",
    "linear_label_size_records",
    "linear_tick_interval",
    "scalar_axis_tick_font_rule",
    "scalar_axis_tick_font_size",
    "sliding_window",
]
