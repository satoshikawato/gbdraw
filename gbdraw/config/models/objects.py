from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Optional

from .common import ShortLongFloatConfig


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

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsGcContentConfig":
        return cls(
            stroke_color=str(d["stroke_color"]),
            stroke_width=float(d["stroke_width"]),
            fill_opacity=float(d["fill_opacity"]),
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


@dataclass(frozen=True)
class ObjectsBlastMatchConfig:
    fill_opacity: float
    stroke_color: str
    stroke_width: float

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsBlastMatchConfig":
        return cls(
            fill_opacity=float(d["fill_opacity"]),
            stroke_color=str(d["stroke_color"]),
            stroke_width=float(d["stroke_width"]),
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
            font_size=ShortLongFloatConfig.from_dict(d["font_size"]),
        )


@dataclass(frozen=True)
class ObjectsDefinitionCircularConfig:
    interval: int
    font_size: float

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsDefinitionCircularConfig":
        return cls(interval=int(d["interval"]), font_size=float(d["font_size"]))


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
    stroke_width: float
    font_weight: str
    font_size: ShortLongFloatConfig
    interval: Optional[int]

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "ObjectsScaleConfig":
        interval = d.get("interval")
        return cls(
            style=str(d.get("style", "bar")),
            stroke_color=str(d["stroke_color"]),
            stroke_width=float(d["stroke_width"]),
            font_weight=str(d["font_weight"]),
            font_size=ShortLongFloatConfig.from_dict(d["font_size"]),
            interval=int(interval) if interval is not None else None,
        )


@dataclass(frozen=True)
class ObjectsConfig:
    scale: ObjectsScaleConfig
    sliding_window: ObjectsSlidingWindowConfig
    text: ObjectsTextConfig
    gc_content: ObjectsGcContentConfig
    gc_skew: ObjectsGcSkewConfig
    blast_match: ObjectsBlastMatchConfig
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
            blast_match=ObjectsBlastMatchConfig.from_dict(d["blast_match"]),
            features=ObjectsFeaturesConfig.from_dict(d["features"]),
            legends=ObjectsLegendsConfig.from_dict(d["legends"]),
            ticks=ObjectsTicksConfig.from_dict(d["ticks"]),
            axis=ObjectsAxisConfig.from_dict(d["axis"]),
            definition=ObjectsDefinitionConfig.from_dict(d["definition"]),
        )


__all__ = ["ObjectsConfig"]


