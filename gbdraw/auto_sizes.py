"""Length- and mode-dependent Auto sizes.

A setting left at Auto has no explicit value: its value follows the drawn record
lengths and the diagram mode. This module owns that knowledge for every
renderer and surface: the size class (short or long), the sliding-window and
scale-tick tiers, the derived Depth window, the scalar-axis tick font, the
circular track widths, and the Auto value of every such setting of one drawing.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from functools import reduce
from typing import Literal, Protocol, get_args

from gbdraw.config.models import GbdrawConfig
from gbdraw.config.models.labels import LINEAR_LABEL_SCOPE_RECORDS, LinearLabelScope
from gbdraw.mode_profiles import DiagramMode


SizeClass = Literal["short", "long"]
SIZE_CLASSES: tuple[SizeClass, ...] = get_args(SizeClass)
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
# The settings whose Auto value each record of a drawing takes for itself (``record_auto_values``).
_RECORD_SETTINGS = (*_WINDOW_SETTINGS, "objects.scale.interval")
# The Auto value a feature density track takes from the window and step settings, keyed in ``record_auto_values``.
_DENSITY_WINDOW_KEYS = {"window": "feature_density.window", "step": "feature_density.step"}
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
        (int(cfg.objects.sliding_window.up100m[0]), int(cfg.objects.sliding_window.up100m[1])),
    )
    window = max(1, min(tier_window, int(length) // _SLIDING_WINDOW_DIVISOR))
    return window, max(1, min(tier_step, window // _SLIDING_WINDOW_DIVISOR))


def sliding_window_tiers(cfg: GbdrawConfig) -> tuple[tuple[int, int, int], ...]:
    """Return the Auto GC (below bp, window, step) tiers; records of 100 Mbp or longer use ``up100m``."""
    tiers = cfg.objects.sliding_window
    return (
        (1_000_000, int(tiers.default[0]), int(tiers.default[1])),
        (10_000_000, int(tiers.up1m[0]), int(tiers.up1m[1])),
        (100_000_000, int(tiers.up10m[0]), int(tiers.up10m[1])),
    )


def tick_and_window_tiers(cfg: GbdrawConfig) -> dict[str, object]:
    """Return the length tiers of the Auto ticks and sliding window, for the Web hints.

    Each ``tiers`` row is checked in order against the record length (circular
    ticks: at most; the others: below); a longer record takes ``above``.
    """
    up100m = cfg.objects.sliding_window.up100m
    return {
        "slidingWindow": {
            "tiers": [list(row) for row in sliding_window_tiers(cfg)],
            "above": [int(up100m[0]), int(up100m[1])],
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


def record_sliding_windows(
    mode: DiagramMode, record_lengths: Sequence[int], cfg: GbdrawConfig
) -> tuple[tuple[int, int], ...]:
    """Return the Auto GC content and skew (window, step) each record of a drawing draws with, in drawing order.

    A Circular record draws at its own length (each record is a circle of its own); the records of a Linear
    drawing share one axis, so each takes the longest record's. The one reading of this rule: the renderer, the
    Auto values of a drawing (``record_auto_values``) and the region derivation (``adapt_explicit_settings``).
    """
    lengths = (max(record_lengths),) * len(record_lengths) if mode == "linear" else tuple(record_lengths)
    return tuple(sliding_window(length, cfg) for length in lengths)


def feature_density_windows(
    record_lengths: Sequence[int], cfg: GbdrawConfig, *, window: int | None = None, step: int | None = None
) -> tuple[tuple[int, int], ...]:
    """Return the (window, step) each record of one diagram draws its feature density tracks with when a track
    sets none, in drawing order.

    The records of one diagram, in either mode, share one window, so a track has one window and one Legend scale:
    the drawing's hand-set ``window`` and ``step``, else the longest record's Auto values. Each diagram of a
    Circular batch is a diagram of one record. The one reading of this rule: the Circular and Linear assemblers and
    the Auto values of a drawing (``record_auto_values``).
    """
    auto_window, auto_step = sliding_window(max(record_lengths), cfg)
    drawn = (auto_window if window is None else int(window), auto_step if step is None else int(step))
    return (drawn,) * len(record_lengths)


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
    if renderer in {"depth", "feature_density"}:
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
    ``density_window_settings`` names the settings (``window``, ``step``) that a feature density track the
    drawing draws takes, being one that sets none of its own. ``separate_records``: each record is a diagram of its
    own (a Circular batch), so it takes the sizes and the feature density window of its own length.
    """

    mode: DiagramMode
    record_lengths: tuple[int, ...]
    wrap_row_length: int | None = None
    depth_rows: tuple[DrawnDepthRow, ...] = ()
    density_window_settings: tuple[str, ...] = ()
    separate_records: bool = False

    def __post_init__(self) -> None:
        if self.mode not in _MODE_LENGTH_PAIRS:
            raise ValueError(f"Unknown diagram mode: {self.mode!r}")
        if not self.record_lengths or min(self.record_lengths) < 1:
            raise ValueError("A drawn extent needs at least one record of length >= 1")
        if not set(self.density_window_settings) <= set(_DENSITY_WINDOW_KEYS):
            raise ValueError(f"Unknown density window settings: {self.density_window_settings!r}")

    def size_class(self, cfg: GbdrawConfig) -> SizeClass:
        """Return the drawing's size class: that of its longest record."""
        threshold = int(getattr(cfg.labels.length_threshold, self.mode))
        return determine_length_parameter(self.drawn_axis_length, threshold)

    @property
    def drawn_axis_length(self) -> int:
        """The bases on the longest drawn axis."""
        longest = max(self.record_lengths)
        return longest if self.wrap_row_length is None else min(longest, self.wrap_row_length)


def record_auto_values(extent: DrawnExtent, cfg: GbdrawConfig) -> tuple[Mapping[str, int], ...]:
    """Return the Auto window, step, Depth window, Depth step and scale interval each record of ``extent`` draws
    with, in drawing order.

    A Circular record takes them at its own length (``record_sliding_windows``, ``circular_tick_intervals``); the
    records of a Linear drawing share the longest record's windows and the drawn axis's interval. Where a feature
    density track takes the drawing's window or step (``extent.density_window_settings``), the record also lists
    the Auto value that track draws with (``feature_density_windows``, per diagram: each record of a batch is one),
    keyed ``feature_density.window`` and ``feature_density.step``.
    """
    records: list[Mapping[str, int]] = []
    windows = record_sliding_windows(extent.mode, extent.record_lengths, cfg)
    density = (
        tuple(window for length in extent.record_lengths for window in feature_density_windows((length,), cfg))
        if extent.separate_records
        else feature_density_windows(extent.record_lengths, cfg)
    )
    for length, (window, step), density_values in zip(extent.record_lengths, windows, density, strict=True):
        interval = (
            circular_tick_intervals(length)[0]
            if extent.mode == "circular"
            else linear_tick_interval(extent.drawn_axis_length)
        )
        record = dict(zip(_RECORD_SETTINGS, (window, step, *depth_sliding_window(window, step), interval)))
        for setting, value in zip(("window", "step"), density_values):
            if setting in extent.density_window_settings:
                record[_DENSITY_WINDOW_KEYS[setting]] = value
        records.append(record)
    return tuple(records)


def _shared_record_value(records: Iterable[Mapping[str, int]], setting: str) -> int | None:
    """The value of ``setting`` that every record takes, or None when they differ."""
    values = {record[setting] for record in records}
    return values.pop() if len(values) == 1 else None


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
    The size class follows the longest drawn axis; in a Circular canvas the
    records then share the long styles. The windows and the scale interval are
    those every record takes (``record_auto_values``): a setting that the
    records of a Circular drawing take at different values is left out. The
    GC content and Depth tick fonts follow the track thickness: the Auto thickness unless
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
    size = extent.size_class(cfg)
    values: dict[str, AutoValue] = {
        path: float(getattr(reduce(getattr, path.split("."), cfg), size))
        for path in (*_SHARED_LENGTH_PAIRS, *_MODE_LENGTH_PAIRS[mode])
    }
    records = record_auto_values(extent, cfg)

    def put_shared(*settings: str) -> None:
        for setting in settings:
            shared = _shared_record_value(records, setting)
            if shared is not None:
                values[setting] = shared

    put_shared(*_WINDOW_SETTINGS)
    if mode == "circular":
        circular = cfg.canvas.circular
        for slot_id, renderer in _CIRCULAR_SLOT_RENDERERS.items():
            values[f"circularSlot.{slot_id}.width"] = circular_track_width_px(
                renderer,
                radius=circular.radius,
                track_ratio=circular.track_ratio,
                factors=circular.track_ratio_factors[size],
            )
        put_shared("objects.scale.interval")
    else:
        values["labels.font_size.linear"] = linear_label_font_size((extent.drawn_axis_length,), cfg)
        put_shared("objects.scale.interval")
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

    An Auto value of ``None`` means that drawing's mode does not draw the setting, or that the records compared
    take different values: for a setting each record takes for itself (``record_auto_values``), or in the
    diagrams of a batch source, each at its own size class. ``explicit`` is the value the source diagram drew: the
    hand-set one, or its Auto value where a short/long pair is set by hand only for the other size class.
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
    explicit: Mapping[str, Mapping[SizeClass, AutoValue]],
    *,
    source: DrawnExtent,
    target: DrawnExtent,
    cfg: GbdrawConfig,
    record_sources: Sequence[int | None] | None = None,
) -> SettingAdaptation:
    """Keep a hand-set size only where its Auto value is the same in ``target``.

    ``explicit`` holds, for each hand-set setting keyed as
    ``auto_setting_values`` is, the value a diagram of each size class draws:
    a short/long pair gives a class only the value of its own leaf, and a
    class it gives none draws its Auto value. A setting whose Auto
    value changes with the size class, a tier or the mode returns to Auto and
    is listed in ``reset``; the other sizes are listed in ``kept``. A value
    that neither drawing lists is not a size and is in neither. The tick fonts
    are compared at the track thickness each drawing ends with, and its depth
    rows at the Depth width it ends with.

    A setting each record takes for itself (the windows and the scale
    interval, ``record_auto_values``) is compared record by record:
    ``record_sources`` gives, for each record of ``target``, the index of the
    record of ``source`` it comes from (None for a record from elsewhere, and
    no list when none comes from ``source``); a record from elsewhere compares
    with every record of ``source``. It is kept only where every record keeps
    its Auto value. A batch source (``separate_records``) is a diagram per
    record, each at its own size class: the other sizes are compared with the
    diagram of each record the target's records come from. ``kept`` and
    ``reset`` name the value the first of those source diagrams, in the
    source's order, drew at its own size class (``DrawnExtent.size_class``):
    its hand-set value, else its Auto value; a size the source's mode does
    not draw names its hand-set value of the first class that sets it.
    """
    if record_sources is None:
        record_sources = (None,) * len(target.record_lengths)
    if len(record_sources) != len(target.record_lengths):
        raise ValueError("record_sources needs one entry per target record")
    source_records = record_auto_values(source, cfg)
    target_records = record_auto_values(target, cfg)
    compared = tuple(
        (target_record, source_records if index is None else (source_records[index],))
        for target_record, index in zip(target_records, record_sources, strict=True)
    )
    # The source records the target's records come from: the values a reset names.
    from_source = tuple(record for _, records in compared for record in records)
    # The source diagrams the target's records come from: a batch draws one per record.
    source_diagrams: tuple[DrawnExtent, ...] = (source,)
    if source.separate_records:
        used = {index for index in record_sources if index is not None}
        source_diagrams = tuple(
            replace(source, record_lengths=(length,), separate_records=False)
            for index, length in enumerate(source.record_lengths)
            if None in record_sources or index in used
        )

    def hand_set(size: SizeClass) -> dict[str, AutoValue]:
        """The hand-set values a diagram of ``size`` draws."""
        return {setting: values[size] for setting, values in explicit.items() if size in values}

    # The hand-set values each source diagram compared drew, at its own size class; the target draws its own
    # class's.
    source_explicit = tuple(hand_set(diagram.size_class(cfg)) for diagram in source_diagrams)

    def record_keys(setting: str) -> tuple[str, ...]:
        """The per-record Auto values a setting sets: its own and, for the window and step, a density track's."""
        return (setting, _DENSITY_WINDOW_KEYS[setting]) if setting in _DENSITY_WINDOW_KEYS else (setting,)

    def same_record_auto(key: str) -> bool:
        # A density key is compared where both records draw such a track.
        return all(
            target_record[key] == source_record[key]
            for target_record, records in compared
            for source_record in records
            if key in target_record and key in source_record
        )

    def changed_key(setting: str) -> str | None:
        return next((key for key in record_keys(setting) if not same_record_auto(key)), None)

    def same_auto(setting: str, target_values: Mapping[str, AutoValue]) -> bool:
        if setting in _RECORD_SETTINGS:
            return changed_key(setting) is None
        return all(values.get(setting) == target_values.get(setting) for values in source_autos)

    def autos_of(setting: str) -> tuple[AutoValue | None, AutoValue | None]:
        """The source and target Auto values a reset names: those of the per-record value that changed, and the
        value the source diagrams share (None where they differ)."""
        if setting not in _RECORD_SETTINGS:
            shared = {values.get(setting) for values in source_autos}
            return (shared.pop() if len(shared) == 1 else None), target_auto.get(setting)
        key = changed_key(setting) or setting
        return (
            _shared_record_value((record for record in from_source if key in record), key),
            _shared_record_value((record for record in target_records if key in record), key),
        )

    def thickness(mode: DiagramMode, values: Mapping[str, AutoValue]) -> dict[str, float]:
        return {
            track: float(values[setting])
            for track, setting in _TRACK_THICKNESS_SETTINGS[mode].items()
            if setting in values
        }

    depth_font = _TICK_FONT_SETTINGS["depth"]
    source_autos = tuple(
        auto_setting_values(
            diagram,
            cfg,
            track_thickness_px=thickness(source.mode, values),
            depth_tick_font_size=values.get(depth_font),
        )
        for diagram, values in zip(source_diagrams, source_explicit, strict=True)
    )
    target_auto = auto_setting_values(target, cfg)
    target_explicit = hand_set(target.size_class(cfg))
    unchanged = {setting: value for setting, value in target_explicit.items() if same_auto(setting, target_auto)}
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

    def drawn(setting: str) -> AutoValue:
        """The value the first source diagram compared drew (see the docstring)."""
        value = source_explicit[0].get(setting, source_autos[0].get(setting))
        if value is None:
            return next(explicit[setting][size] for size in SIZE_CLASSES if size in explicit[setting])
        return value

    for setting in explicit:
        value = drawn(setting)
        listed = setting in target_auto or any(setting in values for values in source_autos)
        if not listed and setting not in _RECORD_SETTINGS:
            continue
        if same_auto(setting, target_auto):
            kept[setting] = value
        else:
            reset.append(
                SettingReset(setting, value, *autos_of(setting))
            )
    return SettingAdaptation(kept=kept, reset=tuple(reset))


__all__ = [
    "ADAPTABLE_SETTINGS",
    "AutoValue",
    "DrawnDepthRow",
    "DrawnExtent",
    "SettingAdaptation",
    "SettingReset",
    "SIZE_CLASSES",
    "SizeClass",
    "adapt_explicit_settings",
    "auto_setting_values",
    "circular_tick_intervals",
    "circular_track_width_px",
    "TickFontSource",
    "TickedAxisSource",
    "depth_axis_tick_font_size_px",
    "depth_sliding_window",
    "feature_density_windows",
    "determine_length_parameter",
    "gc_axis_tick_font_size_of_depth",
    "gc_axis_tick_font_size_of_depth_rows",
    "linear_label_font_size",
    "linear_label_size_lengths",
    "linear_label_size_records",
    "linear_tick_interval",
    "record_auto_values",
    "record_sliding_windows",
    "scalar_axis_tick_font_rule",
    "scalar_axis_tick_font_size",
    "sliding_window",
]
