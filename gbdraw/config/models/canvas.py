from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Literal, Mapping, TypeAlias

from gbdraw.exceptions import ValidationError

LinearTrackLayoutMode: TypeAlias = Literal["above", "middle", "below"]


def _positive_finite_float(value: Any, *, field_name: str) -> float:
    if isinstance(value, bool):
        raise ValidationError(f"{field_name} must be a positive finite number")
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{field_name} must be a positive finite number") from exc
    if not math.isfinite(parsed) or parsed <= 0:
        raise ValidationError(f"{field_name} must be a positive finite number")
    return parsed


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
    width: CircularCanvasWidthConfig
    track_ratio_factors: dict[str, list[float]]
    track_dict: dict[str, dict[str, dict[str, float]]]

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "CircularCanvasConfig":
        if "show_labels" in d:
            raise ValidationError(
                "canvas.circular.show_labels is no longer accepted; use "
                "labels.circular.scope"
            )
        if "allow_inner_labels" in d:
            raise ValidationError(
                "canvas.circular.allow_inner_labels is no longer accepted; use "
                "labels.circular.scope"
            )
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
    track_spacing: float
    comparison_height: float
    canvas_padding: float
    definition_gap: float
    default_gc_height: float
    depth_height: float
    depth_padding: float
    track_layout: LinearTrackLayoutMode
    track_axis_gap: float | None
    ruler_on_axis: bool
    align_center: bool
    keep_definition_left_aligned: bool
    normalize_length: bool
    default_cds_height: LinearCanvasDefaultCdsHeightConfig
    arrow_length_parameter: LinearCanvasArrowLengthParameterConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LinearCanvasConfig":
        if "show_labels" in d:
            raise ValidationError(
                "canvas.linear.show_labels is no longer accepted; use "
                "labels.linear.scope"
            )
        track_layout_value = d.get("track_layout", "middle")
        if not isinstance(track_layout_value, str):
            raise ValidationError(
                "canvas.linear.track_layout must be one of: above, middle, below"
            )
        track_layout_raw = track_layout_value.strip().lower()
        if track_layout_raw not in {"above", "middle", "below"}:
            raise ValidationError(
                "canvas.linear.track_layout must be one of: above, middle, below"
            )
        track_layout: Literal["above", "middle", "below"] = track_layout_raw  # type: ignore[assignment]
        track_axis_gap_raw = d.get("track_axis_gap", "auto")
        track_axis_gap: float | None
        if track_axis_gap_raw is None:
            track_axis_gap = None
        elif isinstance(track_axis_gap_raw, str):
            normalized_gap = track_axis_gap_raw.strip().lower()
            if normalized_gap in {"", "auto"}:
                track_axis_gap = None
            else:
                try:
                    parsed_gap = float(normalized_gap)
                except ValueError:
                    parsed_gap = -1.0
                track_axis_gap = parsed_gap if parsed_gap >= 0.0 else None
        else:
            try:
                parsed_gap = float(track_axis_gap_raw)
            except (TypeError, ValueError):
                parsed_gap = -1.0
            track_axis_gap = parsed_gap if parsed_gap >= 0.0 else None
        return cls(
            width=int(d["width"]),
            vertical_offset=float(d["vertical_offset"]),
            horizontal_offset=float(d["horizontal_offset"]),
            vertical_padding=float(d["vertical_padding"]),
            track_spacing=max(
                0.0,
                float(d.get("track_spacing", 0.0)),
            ),
            comparison_height=_positive_finite_float(
                d["comparison_height"], field_name="comparison_height"
            ),
            canvas_padding=float(d["canvas_padding"]),
            definition_gap=max(0.0, float(d.get("definition_gap", 20))),
            default_gc_height=float(d["default_gc_height"]),
            depth_height=float(d.get("depth_height", d.get("default_gc_height", 20))),
            depth_padding=float(d.get("depth_padding", d.get("vertical_padding", 8))),
            track_layout=track_layout,
            track_axis_gap=track_axis_gap,
            ruler_on_axis=bool(d.get("ruler_on_axis", False)),
            align_center=bool(d["align_center"]),
            keep_definition_left_aligned=bool(d.get("keep_definition_left_aligned", False)),
            normalize_length=bool(d["normalize_length"]),
            default_cds_height=LinearCanvasDefaultCdsHeightConfig.from_dict(d["default_cds_height"]),
            arrow_length_parameter=LinearCanvasArrowLengthParameterConfig.from_dict(d["arrow_length_parameter"]),
        )


@dataclass(frozen=True)
class CanvasConfig:
    dpi: int
    show_gc: bool
    show_skew: bool
    show_depth: bool
    strandedness: bool
    resolve_overlaps: bool
    circular: CircularCanvasConfig
    linear: LinearCanvasConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "CanvasConfig":
        if "show_labels" in d:
            raise ValidationError(
                "canvas.show_labels is no longer accepted; use "
                "labels.circular.scope and labels.linear.scope"
            )
        return cls(
            dpi=int(d["dpi"]),
            show_gc=bool(d["show_gc"]),
            show_skew=bool(d["show_skew"]),
            show_depth=bool(d.get("show_depth", False)),
            strandedness=bool(d["strandedness"]),
            resolve_overlaps=bool(d["resolve_overlaps"]),
            circular=CircularCanvasConfig.from_dict(d["circular"]),
            linear=LinearCanvasConfig.from_dict(d["linear"]),
        )


__all__ = [
    "CanvasConfig",
    "CircularCanvasConfig",
    "LinearCanvasConfig",
    "LinearTrackLayoutMode",
]
