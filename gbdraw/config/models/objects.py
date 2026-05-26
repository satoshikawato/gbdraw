from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Literal, Mapping, Optional

from gbdraw.exceptions import ValidationError

from .common import ShortLongFloatConfig


PairwiseMatchStyle = Literal["ribbon", "curve"]
_PAIRWISE_MATCH_STYLES: tuple[PairwiseMatchStyle, ...] = ("ribbon", "curve")
GcContentMode = Literal["deviation", "percent"]
_GC_CONTENT_MODES: tuple[GcContentMode, ...] = ("deviation", "percent")


def _normalize_pairwise_match_style(value: Any) -> PairwiseMatchStyle:
    normalized = str(value if value is not None else "ribbon").strip().lower()
    if normalized not in _PAIRWISE_MATCH_STYLES:
        raise ValidationError("pairwise_match_style must be one of: ribbon, curve")
    return normalized  # type: ignore[return-value]


def _normalize_curve_tension(value: Any) -> float:
    tension = float(0.5 if value is None else value)
    if not math.isfinite(tension) or tension < 0.0 or tension > 1.0:
        raise ValidationError("curve_tension must be a finite number in [0.0, 1.0]")
    return tension


def _normalize_gc_content_mode(value: Any) -> GcContentMode:
    normalized = str(value if value is not None else "deviation").strip().lower()
    if normalized not in _GC_CONTENT_MODES:
        raise ValidationError("gc_content_mode must be one of: deviation, percent")
    return normalized  # type: ignore[return-value]


@dataclass(frozen=True)
class ObjectsTextConfig:
    font_family: str

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsTextConfig":
        return cls(font_family=str(d["font_family"]))


@dataclass(frozen=True)
class ObjectsGcContentConfig:
    stroke_color: str
    stroke_width: float
    fill_opacity: float
    mode: GcContentMode
    min_percent: float | None
    max_percent: float | None
    show_axis: bool
    show_ticks: bool
    large_tick_interval: float | None
    small_tick_interval: float | None
    tick_font_size: float | None

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsGcContentConfig":
        min_percent = _optional_gc_content_bound(d.get("min_percent", 0))
        max_percent = _optional_gc_content_bound(d.get("max_percent", 100))
        if min_percent is not None and max_percent is not None and min_percent > max_percent:
            raise ValidationError("gc_content_min_percent must be <= gc_content_max_percent")
        return cls(
            stroke_color=str(d["stroke_color"]),
            stroke_width=float(d["stroke_width"]),
            fill_opacity=float(d["fill_opacity"]),
            mode=_normalize_gc_content_mode(d.get("mode", "deviation")),
            min_percent=min_percent,
            max_percent=max_percent,
            show_axis=_bool_from_config(d.get("show_axis", True), default=True),
            show_ticks=_bool_from_config(d.get("show_ticks", True), default=True),
            large_tick_interval=_optional_gc_content_positive(d.get("large_tick_interval", 20)),
            small_tick_interval=_optional_gc_content_positive(d.get("small_tick_interval", None)),
            tick_font_size=_optional_gc_content_positive(d.get("tick_font_size", None)),
        )


@dataclass(frozen=True)
class ObjectsGcSkewConfig:
    stroke_color: str
    stroke_width: float
    fill_opacity: float

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsGcSkewConfig":
        return cls(
            stroke_color=str(d["stroke_color"]),
            stroke_width=float(d["stroke_width"]),
            fill_opacity=float(d["fill_opacity"]),
        )


def _optional_depth_bound(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"", "auto", "none", "null"}:
            return None
        return float(normalized)
    return float(value)


def _optional_depth_positive(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"", "auto", "none", "null"}:
            return None
        return float(normalized)
    return float(value)


def _optional_gc_content_bound(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"", "auto", "none", "null"}:
            return None
        value = normalized
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValidationError("gc_content_min_percent/gc_content_max_percent must be finite numbers")
    return parsed


def _optional_gc_content_positive(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"", "auto", "none", "null"}:
            return None
        value = normalized
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise ValidationError("gc_content tick intervals and tick_font_size must be positive finite numbers")
    return parsed


def _bool_from_config(value: Any, *, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"false", "0", "no", "off"}:
            return False
        if normalized in {"true", "1", "yes", "on"}:
            return True
    return bool(value)


@dataclass(frozen=True)
class ObjectsDepthConfig:
    fill_color: str
    stroke_color: str
    stroke_width: float
    fill_opacity: float
    normalize: bool
    min_depth: float | None
    max_depth: float | None
    show_axis: bool
    show_ticks: bool
    tick_interval: float | None
    large_tick_interval: float | None
    small_tick_interval: float | None
    tick_font_size: float | None
    share_axis: bool

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsDepthConfig":
        legacy_tick_interval = _optional_depth_positive(d.get("tick_interval", None))
        large_tick_interval = _optional_depth_positive(d.get("large_tick_interval", None))
        if large_tick_interval is None:
            large_tick_interval = legacy_tick_interval
        return cls(
            fill_color=str(d.get("fill_color", "#4A90E2")),
            stroke_color=str(d.get("stroke_color", "none")),
            stroke_width=float(d.get("stroke_width", 0)),
            fill_opacity=float(d.get("fill_opacity", 0.8)),
            normalize=_bool_from_config(d.get("normalize", False), default=False),
            min_depth=_optional_depth_bound(d.get("min_depth", None)),
            max_depth=_optional_depth_bound(d.get("max_depth", None)),
            show_axis=_bool_from_config(d.get("show_axis", True), default=True),
            show_ticks=_bool_from_config(d.get("show_ticks", True), default=True),
            tick_interval=large_tick_interval,
            large_tick_interval=large_tick_interval,
            small_tick_interval=_optional_depth_positive(d.get("small_tick_interval", None)),
            tick_font_size=_optional_depth_positive(d.get("tick_font_size", None)),
            share_axis=_bool_from_config(d.get("share_axis", False), default=False),
        )


@dataclass(frozen=True)
class ObjectsBlastMatchConfig:
    fill_opacity: float
    stroke_color: str
    stroke_width: float
    style: PairwiseMatchStyle | str
    curve_tension: float

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsBlastMatchConfig":
        return cls(
            fill_opacity=float(d["fill_opacity"]),
            stroke_color=str(d["stroke_color"]),
            stroke_width=float(d["stroke_width"]),
            style=_normalize_pairwise_match_style(d.get("style", "ribbon")),
            curve_tension=_normalize_curve_tension(d.get("curve_tension", 0.5)),
        )


def _optional_auto_positive_float(value: Any, *, field_name: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"", "auto", "none", "null"}:
            return None
        value = normalized
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise ValidationError(f"{field_name} must be a positive finite number or auto")
    return parsed


def _normalize_conservation_reference(value: Any) -> str:
    normalized = str(value if value is not None else "auto").strip().lower()
    if normalized not in {"query", "subject", "auto"}:
        raise ValidationError("conservation reference must be one of: query, subject, auto")
    return normalized


@dataclass(frozen=True)
class ObjectsConservationConfig:
    fill_opacity: float
    stroke_color: str
    stroke_width: float
    min_color: str
    max_color: str
    background_color: str
    background_opacity: float
    show_background: bool
    ring_width: float | None
    ring_gap: float | None
    reference: str

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsConservationConfig":
        return cls(
            fill_opacity=float(d.get("fill_opacity", 1.0)),
            stroke_color=str(d.get("stroke_color", "none")),
            stroke_width=float(d.get("stroke_width", 0)),
            min_color=str(d.get("min_color", "#f0f1f5")),
            max_color=str(d.get("max_color", "#8b9cc1")),
            background_color=str(d.get("background_color", "#f2f2f2")),
            background_opacity=float(d.get("background_opacity", 0.35)),
            show_background=_bool_from_config(d.get("show_background", True), default=True),
            ring_width=_optional_auto_positive_float(d.get("ring_width", None), field_name="conservation ring_width"),
            ring_gap=_optional_auto_positive_float(d.get("ring_gap", None), field_name="conservation ring_gap"),
            reference=_normalize_conservation_reference(d.get("reference", "auto")),
        )


@dataclass(frozen=True)
class ObjectsFeaturesConfig:
    features_drawn: list[str]
    block_stroke_color: str
    line_stroke_color: str
    font_weight: str
    block_stroke_width: ShortLongFloatConfig
    line_stroke_width: ShortLongFloatConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsFeaturesConfig":
        return cls(
            features_drawn=list(d.get("features_drawn", [])),
            block_stroke_color=str(d["block_stroke_color"]),
            line_stroke_color=str(d["line_stroke_color"]),
            font_weight=str(d.get("font_weight", "normal")),
            block_stroke_width=ShortLongFloatConfig.from_dict(d["block_stroke_width"]),
            line_stroke_width=ShortLongFloatConfig.from_dict(d["line_stroke_width"]),
        )


@dataclass(frozen=True)
class ObjectsLegendsConfig:
    font_size: ShortLongFloatConfig
    font_weight: str
    text_anchor: str
    color_rect_size: ShortLongFloatConfig
    dominant_baseline: str
    stroke: str
    fill: str

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsLegendsConfig":
        return cls(
            font_size=ShortLongFloatConfig.from_dict(d["font_size"]),
            font_weight=str(d["font_weight"]),
            text_anchor=str(d["text_anchor"]),
            color_rect_size=ShortLongFloatConfig.from_dict(d["color_rect_size"]),
            dominant_baseline=str(d["dominant_baseline"]),
            stroke=str(d.get("stroke", "none")),
            fill=str(d.get("fill", "black")),
        )


@dataclass(frozen=True)
class ObjectsSlidingWindowConfig:
    default: tuple[int, int]
    up1m: tuple[int, int]
    up10m: tuple[int, int]

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsSlidingWindowConfig":
        def _pair(v: Any) -> tuple[int, int]:
            a = list(v)
            return int(a[0]), int(a[1])

        return cls(default=_pair(d["default"]), up1m=_pair(d["up1m"]), up10m=_pair(d["up10m"]))


@dataclass(frozen=True)
class ObjectsTickLabelsConfig:
    stroke: str
    fill: str
    font_size: float
    font_weight: str

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsTickLabelsConfig":
        return cls(
            stroke=str(d["stroke"]),
            fill=str(d["fill"]),
            font_size=float(d["font_size"]),
            font_weight=str(d["font_weight"]),
        )


@dataclass(frozen=True)
class ObjectsTicksConfig:
    tick_width: float
    tick_labels: ObjectsTickLabelsConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsTicksConfig":
        return cls(tick_width=float(d["tick_width"]), tick_labels=ObjectsTickLabelsConfig.from_dict(d["tick_labels"]))


@dataclass(frozen=True)
class ObjectsAxisCircularConfig:
    stroke_color: str
    stroke_width: ShortLongFloatConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsAxisCircularConfig":
        return cls(stroke_color=str(d["stroke_color"]), stroke_width=ShortLongFloatConfig.from_dict(d["stroke_width"]))


@dataclass(frozen=True)
class ObjectsAxisLinearConfig:
    stroke_color: str
    stroke_width: ShortLongFloatConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsAxisLinearConfig":
        return cls(stroke_color=str(d["stroke_color"]), stroke_width=ShortLongFloatConfig.from_dict(d["stroke_width"]))


@dataclass(frozen=True)
class ObjectsAxisConfig:
    circular: ObjectsAxisCircularConfig
    linear: ObjectsAxisLinearConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsAxisConfig":
        return cls(
            circular=ObjectsAxisCircularConfig.from_dict(d["circular"]),
            linear=ObjectsAxisLinearConfig.from_dict(d["linear"]),
        )


@dataclass(frozen=True)
class ObjectsDefinitionLinearConfig:
    stroke: str
    fill: str
    font_weight: str
    text_anchor: str
    dominant_baseline: str
    interval: int
    show_replicon: bool
    show_accession: bool
    show_length: bool
    font_size: ShortLongFloatConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsDefinitionLinearConfig":
        return cls(
            stroke=str(d["stroke"]),
            fill=str(d["fill"]),
            font_weight=str(d["font_weight"]),
            text_anchor=str(d["text_anchor"]),
            dominant_baseline=str(d["dominant_baseline"]),
            interval=int(d["interval"]),
            show_replicon=bool(d.get("show_replicon", False)),
            show_accession=bool(d.get("show_accession", True)),
            show_length=bool(d.get("show_length", True)),
            font_size=ShortLongFloatConfig.from_dict(d["font_size"]),
        )


@dataclass(frozen=True)
class ObjectsDefinitionCircularConfig:
    interval: int
    font_size: float
    plot_title_font_size: float

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsDefinitionCircularConfig":
        return cls(
            interval=int(d["interval"]),
            font_size=float(d["font_size"]),
            plot_title_font_size=float(d.get("plot_title_font_size", d["font_size"])),
        )


@dataclass(frozen=True)
class ObjectsDefinitionConfig:
    linear: ObjectsDefinitionLinearConfig
    circular: ObjectsDefinitionCircularConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsDefinitionConfig":
        return cls(
            linear=ObjectsDefinitionLinearConfig.from_dict(d["linear"]),
            circular=ObjectsDefinitionCircularConfig.from_dict(d["circular"]),
        )


@dataclass(frozen=True)
class ObjectsScaleConfig:
    style: str
    stroke_color: str
    label_color: str
    stroke_width: float
    font_weight: str
    font_size: ShortLongFloatConfig
    ruler_label_font_size: ShortLongFloatConfig
    interval: Optional[int]

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsScaleConfig":
        interval = d.get("interval")
        font_size = ShortLongFloatConfig.from_dict(d["font_size"])
        return cls(
            style=str(d.get("style", "bar")),
            stroke_color=str(d["stroke_color"]),
            label_color=str(d.get("label_color", "black")),
            stroke_width=float(d["stroke_width"]),
            font_weight=str(d["font_weight"]),
            font_size=font_size,
            ruler_label_font_size=(
                ShortLongFloatConfig.from_dict(d["ruler_label_font_size"])
                if "ruler_label_font_size" in d
                else font_size
            ),
            interval=int(interval) if interval is not None else None,
        )


@dataclass(frozen=True)
class ObjectsConfig:
    scale: ObjectsScaleConfig
    sliding_window: ObjectsSlidingWindowConfig
    text: ObjectsTextConfig
    gc_content: ObjectsGcContentConfig
    gc_skew: ObjectsGcSkewConfig
    depth: ObjectsDepthConfig
    blast_match: ObjectsBlastMatchConfig
    conservation: ObjectsConservationConfig
    features: ObjectsFeaturesConfig
    legends: ObjectsLegendsConfig
    ticks: ObjectsTicksConfig
    axis: ObjectsAxisConfig
    definition: ObjectsDefinitionConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsConfig":
        return cls(
            scale=ObjectsScaleConfig.from_dict(d["scale"]),
            sliding_window=ObjectsSlidingWindowConfig.from_dict(d["sliding_window"]),
            text=ObjectsTextConfig.from_dict(d["text"]),
            gc_content=ObjectsGcContentConfig.from_dict(d["gc_content"]),
            gc_skew=ObjectsGcSkewConfig.from_dict(d["gc_skew"]),
            depth=ObjectsDepthConfig.from_dict(d.get("depth", {})),
            blast_match=ObjectsBlastMatchConfig.from_dict(d["blast_match"]),
            conservation=ObjectsConservationConfig.from_dict(d.get("conservation", {})),
            features=ObjectsFeaturesConfig.from_dict(d["features"]),
            legends=ObjectsLegendsConfig.from_dict(d["legends"]),
            ticks=ObjectsTicksConfig.from_dict(d["ticks"]),
            axis=ObjectsAxisConfig.from_dict(d["axis"]),
            definition=ObjectsDefinitionConfig.from_dict(d["definition"]),
        )


__all__ = ["ObjectsConfig"]
