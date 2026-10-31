from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal, Mapping, cast

from gbdraw.core.numeric import ascii_float, ascii_int

from gbdraw.exceptions import ValidationError
from gbdraw.labels.label_style import (
    DEFAULT_LABEL_LINE,
    DEFAULT_LABEL_PROTECTION,
    LabelLine,
    LabelProtection,
    normalize_label_line,
    normalize_label_protection,
)
from gbdraw.labels.policy import LabelRenderingPolicy, normalize_label_rendering
from gbdraw.labels.circular_types import (
    CircularLabelPlacement,
    normalize_circular_label_placement,
)

from .common import ShortLongFloatConfig

CircularLabelScope = Literal["none", "outer", "both"]
LinearLabelScope = Literal["none", "all", "first", "orthogroup_top"]
LinearLabelScopeRecords = Literal["none", "all", "first"]
def _config_bool(value: object, *, path: str) -> bool:
    """``value`` as a switch; anything but ``true`` or ``false`` is refused, naming ``path``."""

    if not isinstance(value, bool):
        raise ValidationError(
            f"{path} must be true or false", diagnostic={"code": "INPUT_INVALID", "configPath": path}
        )
    return value


# The records whose labels each Linear label display scope selects.
LINEAR_LABEL_SCOPE_RECORDS: Mapping[LinearLabelScope, LinearLabelScopeRecords] = {
    "none": "none",
    "all": "all",
    "first": "first",
    "orthogroup_top": "all",
}


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
        return cls(circular=ascii_int(d["circular"]), linear=ascii_int(d["linear"]))


@dataclass(frozen=True)
class LabelsFontSizeConfig:
    short: float
    long: float
    linear: ShortLongFloatConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LabelsFontSizeConfig":
        return cls(
            short=ascii_float(d["short"]),
            long=ascii_float(d["long"]),
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
            x_radius_offset=ascii_float(d.get("x_radius_offset", 1.0)),
            y_radius_offset=ascii_float(d.get("y_radius_offset", 1.0)),
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
class LabelsLinearConfig:
    scope: LinearLabelScope
    placement: Literal["auto", "above_feature"]
    rotation: float
    # The feature label style (gbdraw.labels.label_style); a stored config without them takes the defaults.
    line: LabelLine
    protection: LabelProtection
    italic_gene_labels: bool

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LabelsLinearConfig":
        raw = d if isinstance(d, dict) else dict(d)
        scope_raw = raw.get("scope", "none")
        if not isinstance(scope_raw, str):
            raise ValidationError(
                "labels.linear.scope must be one of: "
                "none, all, first, orthogroup_top"
            )
        scope = scope_raw.strip().lower()
        if scope not in {"none", "all", "first", "orthogroup_top"}:
            raise ValidationError(
                "labels.linear.scope must be one of: "
                "none, all, first, orthogroup_top"
            )
        placement_value = raw.get("placement", "auto")
        if not isinstance(placement_value, str):
            raise ValidationError(
                "labels.linear.placement must be one of: auto, above_feature"
            )
        placement_raw = placement_value.strip().lower()
        if placement_raw not in {"auto", "above_feature"}:
            raise ValidationError(
                "labels.linear.placement must be one of: auto, above_feature"
            )
        placement = cast(Literal["auto", "above_feature"], placement_raw)
        rotation = ascii_float(raw.get("rotation", 0.0))
        return cls(
            scope=cast(LinearLabelScope, scope),
            placement=placement,
            rotation=rotation,
            line=normalize_label_line(raw.get("line", DEFAULT_LABEL_LINE["linear"]), path="labels.linear.line"),
            protection=normalize_label_protection(
                raw.get("protection", DEFAULT_LABEL_PROTECTION["linear"]), path="labels.linear.protection"
            ),
            italic_gene_labels=_config_bool(
                raw.get("italic_gene_labels", False), path="labels.linear.italic_gene_labels"
            ),
        )


@dataclass(frozen=True)
class LabelsCircularConfig:
    scope: CircularLabelScope
    placement: CircularLabelPlacement
    # The feature label style (gbdraw.labels.label_style); a stored config without them takes the defaults.
    line: LabelLine
    protection: LabelProtection
    italic_gene_labels: bool

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LabelsCircularConfig":
        raw = d if isinstance(d, dict) else dict(d)
        scope_raw = raw.get("scope", "none")
        if not isinstance(scope_raw, str):
            raise ValidationError(
                "labels.circular.scope must be one of: none, outer, both"
            )
        scope = scope_raw.strip().lower()
        if scope not in {"none", "outer", "both"}:
            raise ValidationError(
                "labels.circular.scope must be one of: none, outer, both"
            )
        return cls(
            scope=cast(CircularLabelScope, scope),
            placement=normalize_circular_label_placement(raw.get("placement", "horizontal")),
            line=normalize_label_line(raw.get("line", DEFAULT_LABEL_LINE["circular"]), path="labels.circular.line"),
            protection=normalize_label_protection(
                raw.get("protection", DEFAULT_LABEL_PROTECTION["circular"]), path="labels.circular.protection"
            ),
            italic_gene_labels=_config_bool(
                raw.get("italic_gene_labels", False), path="labels.circular.italic_gene_labels"
            ),
        )


@dataclass(frozen=True)
class LabelsSpacingConfig:
    circular: float
    linear: float

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LabelsSpacingConfig":
        raw = d if isinstance(d, dict) else dict(d)
        return cls(
            circular=ascii_float(raw.get("circular", 3.0)),
            linear=ascii_float(raw.get("linear", 3.0)),
        )


@dataclass(frozen=True)
class LabelsConfig:
    rendering: LabelRenderingPolicy
    filtering: LabelsFilteringConfig
    length_threshold: LabelsLengthThresholdConfig
    font_size: LabelsFontSizeConfig
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
    circular: LabelsCircularConfig
    linear: LabelsLinearConfig
    spacing: LabelsSpacingConfig

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> "LabelsConfig":
        def _factor_3level(table: Mapping[str, Any]) -> dict[str, dict[str, dict[str, float]]]:
            out: dict[str, dict[str, dict[str, float]]] = {}
            for k1, v1 in dict(table).items():
                out[str(k1)] = {}
                for k2, v2 in dict(v1).items():
                    out[str(k1)][str(k2)] = {str(k3): ascii_float(v3) for k3, v3 in dict(v2).items()}
            return out

        def _factor_2level(table: Mapping[str, Any]) -> dict[str, dict[str, float]]:
            return {
                str(k): {str(k2): ascii_float(v2) for k2, v2 in dict(v).items()}
                for k, v in dict(table).items()
            }

        return cls(
            rendering=normalize_label_rendering(d.get("rendering", "auto")),
            filtering=LabelsFilteringConfig.from_dict(d.get("filtering", {})),
            length_threshold=LabelsLengthThresholdConfig.from_dict(d["length_threshold"]),
            font_size=LabelsFontSizeConfig.from_dict(d["font_size"]),
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
            circular=LabelsCircularConfig.from_dict(d.get("circular", {})),
            linear=LabelsLinearConfig.from_dict(d.get("linear", {})),
            spacing=LabelsSpacingConfig.from_dict(d.get("spacing", {})),
        )


__all__ = [
    "LabelsConfig",
    "CircularLabelScope",
    "LINEAR_LABEL_SCOPE_RECORDS",
    "LinearLabelScope",
    "LinearLabelScopeRecords",
    "LabelsFilteringConfig",
    "LabelsLengthThresholdConfig",
    "LabelsFontSizeConfig",
    "LabelsStrokeColorConfig",
    "LabelsUnifiedAdjustmentConfig",
    "LabelsCircularConfig",
    "LabelsLinearConfig",
    "LabelsSpacingConfig",
]

