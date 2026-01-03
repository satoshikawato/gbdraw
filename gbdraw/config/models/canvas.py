from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal, Mapping


@dataclass(frozen=True)
class CircularCanvasWidthConfig:
    with_labels: int
    without_labels: int

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "CircularCanvasWidthConfig":
        return cls(with_labels=int(d["with_labels"]), without_labels=int(d["without_labels"]))


@dataclass(frozen=True)
class CircularCanvasConfig:
    height: int
    radius: float
    track_ratio: float
    track_type: str
    allow_inner_labels: bool
    width: CircularCanvasWidthConfig
    track_ratio_factors: dict[str, list[float]]
    track_dict: dict[str, dict[str, dict[str, float]]]

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "CircularCanvasConfig":
        track_dict_raw = d.get("track_dict", {})
        track_dict: dict[str, dict[str, dict[str, float]]] = {}
        for length_param, by_track_type in dict(track_dict_raw).items():
            track_dict[str(length_param)] = {}
            for track_type, by_track_id in dict(by_track_type).items():
                track_dict[str(length_param)][str(track_type)] = {
                    str(k): float(v) for k, v in dict(by_track_id).items()
                }

        return cls(
            height=int(d["height"]),
            radius=float(d["radius"]),
            track_ratio=float(d["track_ratio"]),
            track_type=str(d["track_type"]),
            allow_inner_labels=bool(d["allow_inner_labels"]),
            width=CircularCanvasWidthConfig.from_dict(d["width"]),
            track_ratio_factors={
                str(k): [float(x) for x in list(v)]
                for k, v in dict(d["track_ratio_factors"]).items()
            },
            track_dict=track_dict,
        )


@dataclass(frozen=True)
class LinearCanvasDefaultCdsHeightConfig:
    short: float
    long: float

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LinearCanvasDefaultCdsHeightConfig":
        return cls(short=float(d["short"]), long=float(d["long"]))


@dataclass(frozen=True)
class LinearCanvasArrowLengthParameterConfig:
    short: float
    long: float

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LinearCanvasArrowLengthParameterConfig":
        return cls(short=float(d["short"]), long=float(d["long"]))


@dataclass(frozen=True)
class LinearCanvasConfig:
    width: int
    vertical_offset: float
    horizontal_offset: float
    vertical_padding: float
    comparison_height: float
    canvas_padding: float
    default_gc_height: float
    align_center: bool
    normalize_length: bool
    default_cds_height: LinearCanvasDefaultCdsHeightConfig
    arrow_length_parameter: LinearCanvasArrowLengthParameterConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LinearCanvasConfig":
        return cls(
            width=int(d["width"]),
            vertical_offset=float(d["vertical_offset"]),
            horizontal_offset=float(d["horizontal_offset"]),
            vertical_padding=float(d["vertical_padding"]),
            comparison_height=float(d["comparison_height"]),
            canvas_padding=float(d["canvas_padding"]),
            default_gc_height=float(d["default_gc_height"]),
            align_center=bool(d["align_center"]),
            normalize_length=bool(d["normalize_length"]),
            default_cds_height=LinearCanvasDefaultCdsHeightConfig.from_dict(d["default_cds_height"]),
            arrow_length_parameter=LinearCanvasArrowLengthParameterConfig.from_dict(d["arrow_length_parameter"]),
        )


@dataclass(frozen=True)
class CanvasConfig:
    dpi: int
    show_gc: bool
    show_skew: bool
    # - bool: typical config.toml / circular CLI usage
    # - str: linear CLI supports mode strings ("all"/"first"/"none")
    show_labels: bool | Literal["all", "first", "none"]
    strandedness: bool
    resolve_overlaps: bool
    circular: CircularCanvasConfig
    linear: LinearCanvasConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "CanvasConfig":
        raw_show_labels = d.get("show_labels", False)
        if raw_show_labels is None:
            show_labels: bool | Literal["all", "first", "none"] = False
        elif isinstance(raw_show_labels, str):
            show_labels = raw_show_labels
        else:
            show_labels = bool(raw_show_labels)
        return cls(
            dpi=int(d["dpi"]),
            show_gc=bool(d["show_gc"]),
            show_skew=bool(d["show_skew"]),
            show_labels=show_labels,
            strandedness=bool(d["strandedness"]),
            resolve_overlaps=bool(d["resolve_overlaps"]),
            circular=CircularCanvasConfig.from_dict(d["circular"]),
            linear=LinearCanvasConfig.from_dict(d["linear"]),
        )


__all__ = [
    "CanvasConfig",
    "CircularCanvasConfig",
    "LinearCanvasConfig",
]


