from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from .common import ShortLongFloatConfig  # type: ignore[reportMissingImports]


@dataclass(frozen=True)
class LabelsFilteringConfig:
    raw: dict[str, Any]
    blacklist_keywords: list[str]
    qualifier_priority: dict[str, list[str]]

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LabelsFilteringConfig":
        raw = d if isinstance(d, dict) else dict(d)
        blacklist = list(raw.get("blacklist_keywords", []))
        qualifier_priority_raw = dict(raw.get("qualifier_priority", {}))
        qualifier_priority = {k: list(v) for k, v in qualifier_priority_raw.items()}
        return cls(raw=raw, blacklist_keywords=blacklist, qualifier_priority=qualifier_priority)

    def as_dict(self) -> dict[str, Any]:
        return self.raw


@dataclass(frozen=True)
class LabelsLengthThresholdConfig:
    circular: int
    linear: int

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LabelsLengthThresholdConfig":
        return cls(circular=int(d["circular"]), linear=int(d["linear"]))


@dataclass(frozen=True)
class LabelsFontSizeConfig:
    short: float
    long: float
    linear: ShortLongFloatConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LabelsFontSizeConfig":
        return cls(
            short=float(d["short"]),
            long=float(d["long"]),
            linear=ShortLongFloatConfig.from_dict(d["linear"]),
        )

    def for_length_param(self, length_param: str) -> float:
        return getattr(self, length_param)


@dataclass(frozen=True)
class LabelsStrokeColorConfig:
    label_stroke_color: str

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LabelsStrokeColorConfig":
        return cls(label_stroke_color=str(d["label_stroke_color"]))


@dataclass(frozen=True)
class LabelsArcOffsetConfig:
    x_radius_offset: float
    y_radius_offset: float

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LabelsArcOffsetConfig":
        return cls(
            x_radius_offset=float(d.get("x_radius_offset", 1.0)),
            y_radius_offset=float(d.get("y_radius_offset", 1.0)),
        )


@dataclass(frozen=True)
class LabelsUnifiedAdjustmentConfig:
    outer_labels: LabelsArcOffsetConfig
    inner_labels: LabelsArcOffsetConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LabelsUnifiedAdjustmentConfig":
        return cls(
            outer_labels=LabelsArcOffsetConfig.from_dict(d.get("outer_labels", {})),
            inner_labels=LabelsArcOffsetConfig.from_dict(d.get("inner_labels", {})),
        )


@dataclass(frozen=True)
class LabelsConfig:
    filtering: LabelsFilteringConfig
    length_threshold: LabelsLengthThresholdConfig
    font_size: LabelsFontSizeConfig
    stroke_width: ShortLongFloatConfig
    stroke_color: LabelsStrokeColorConfig
    radius_factor: dict[str, dict[str, dict[str, float]]]
    inner_radius_factor: dict[str, dict[str, dict[str, float]]]
    arc_x_radius_factor: dict[str, dict[str, dict[str, float]]]
    arc_y_radius_factor: dict[str, dict[str, dict[str, float]]]
    arc_center_x: dict[str, dict[str, float]]
    arc_angle: dict[str, dict[str, float]]
    inner_arc_x_radius_factor: dict[str, dict[str, dict[str, float]]]
    inner_arc_y_radius_factor: dict[str, dict[str, dict[str, float]]]
    inner_arc_center_x: dict[str, dict[str, float]]
    inner_arc_angle: dict[str, dict[str, float]]
    unified_adjustment: LabelsUnifiedAdjustmentConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LabelsConfig":
        def _factor_3level(table: Mapping[str, Any]) -> dict[str, dict[str, dict[str, float]]]:
            out: dict[str, dict[str, dict[str, float]]] = {}
            for k1, v1 in dict(table).items():
                out[str(k1)] = {}
                for k2, v2 in dict(v1).items():
                    out[str(k1)][str(k2)] = {str(k3): float(v3) for k3, v3 in dict(v2).items()}
            return out

        def _factor_2level(table: Mapping[str, Any]) -> dict[str, dict[str, float]]:
            return {
                str(k): {str(k2): float(v2) for k2, v2 in dict(v).items()}
                for k, v in dict(table).items()
            }

        return cls(
            filtering=LabelsFilteringConfig.from_dict(d.get("filtering", {})),
            length_threshold=LabelsLengthThresholdConfig.from_dict(d["length_threshold"]),
            font_size=LabelsFontSizeConfig.from_dict(d["font_size"]),
            stroke_width=ShortLongFloatConfig.from_dict(d["stroke_width"]),
            stroke_color=LabelsStrokeColorConfig.from_dict(d["stroke_color"]),
            radius_factor=_factor_3level(d["radius_factor"]),
            inner_radius_factor=_factor_3level(d["inner_radius_factor"]),
            arc_x_radius_factor=_factor_3level(d["arc_x_radius_factor"]),
            arc_y_radius_factor=_factor_3level(d["arc_y_radius_factor"]),
            arc_center_x=_factor_2level(d["arc_center_x"]),
            arc_angle=_factor_2level(d.get("arc_angle", {})),
            inner_arc_x_radius_factor=_factor_3level(d["inner_arc_x_radius_factor"]),
            inner_arc_y_radius_factor=_factor_3level(d["inner_arc_y_radius_factor"]),
            inner_arc_center_x=_factor_2level(d["inner_arc_center_x"]),
            inner_arc_angle=_factor_2level(d.get("inner_arc_angle", {})),
            unified_adjustment=LabelsUnifiedAdjustmentConfig.from_dict(d.get("unified_adjustment", {})),
        )


__all__ = [
    "LabelsConfig",
    "LabelsFilteringConfig",
    "LabelsLengthThresholdConfig",
    "LabelsFontSizeConfig",
    "LabelsStrokeColorConfig",
    "LabelsUnifiedAdjustmentConfig",
]


