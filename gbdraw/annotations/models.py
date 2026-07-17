"""Immutable public contracts for region annotation tracks."""

from __future__ import annotations

from dataclasses import dataclass, field
from math import isfinite
from types import MappingProxyType
from typing import Literal, Mapping, TypeAlias

from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.core.color import normalize_hex_color
from gbdraw.exceptions import ValidationError
from gbdraw.io.colors import resolve_color_to_hex
from gbdraw.io.record_select import RecordSelector


def _finite(value: object, name: str, *, positive: bool = False, nonnegative: bool = False) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{name} must be numeric.") from exc
    if not isfinite(number):
        raise ValidationError(f"{name} must be finite.")
    if positive and number <= 0:
        raise ValidationError(f"{name} must be positive.")
    if nonnegative and number < 0:
        raise ValidationError(f"{name} must be non-negative.")
    return number


def _color(value: object | None, name: str, *, optional: bool = False) -> str | None:
    if value is None and optional:
        return None
    text = str(value or "").strip()
    if not text:
        if optional:
            return None
        raise ValidationError(f"{name} cannot be empty.")
    try:
        return normalize_hex_color(resolve_color_to_hex(text))
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{name} is not a valid color: {text!r}.") from exc


@dataclass(frozen=True)
class FeatureSelector:
    """One stable feature selector used by :class:`FeatureSpan`.

    ``kind=value`` strings support ``hash``, ``location``, ``record_location``,
    ``type`` and arbitrary qualifier names. A bare value matches a stable hash
    or any qualifier value exactly.
    """

    value: str
    key: str | None = None

    def __post_init__(self) -> None:
        value = str(self.value).strip()
        key = str(self.key).strip() if self.key is not None else None
        if not value:
            raise ValidationError("Feature selector value cannot be empty.")
        if key == "":
            raise ValidationError("Feature selector key cannot be empty.")
        object.__setattr__(self, "value", value)
        object.__setattr__(self, "key", key)


def parse_feature_selector(value: str | FeatureSelector) -> FeatureSelector:
    if isinstance(value, FeatureSelector):
        return value
    text = str(value).strip()
    if not text:
        raise ValidationError("Feature selector cannot be empty.")
    if "=" in text:
        key, selector_value = text.split("=", 1)
        return FeatureSelector(value=selector_value, key=key)
    return FeatureSelector(value=text)


@dataclass(frozen=True)
class CoordinateSpan:
    record: RecordSelector | None
    start: int
    end: int
    coordinate_space: Literal["source", "local"] = "source"
    wraps_origin: bool = False
    out_of_bounds: Literal["clip", "skip", "error"] = "clip"

    def __post_init__(self) -> None:
        if self.record is not None and not isinstance(self.record, RecordSelector):
            raise ValidationError("CoordinateSpan record must be a RecordSelector.")
        if isinstance(self.start, bool) or isinstance(self.end, bool):
            raise ValidationError("Annotation coordinates must be integers.")
        try:
            start, end = int(self.start), int(self.end)
        except (TypeError, ValueError) as exc:
            raise ValidationError("Annotation coordinates must be integers.") from exc
        if start < 1 or end < 1:
            raise ValidationError("Annotation coordinates are 1-based and must be >= 1.")
        if self.coordinate_space not in {"source", "local"}:
            raise ValidationError("coordinate_space must be 'source' or 'local'.")
        if self.out_of_bounds not in {"clip", "skip", "error"}:
            raise ValidationError("out_of_bounds must be 'clip', 'skip', or 'error'.")
        if self.wraps_origin and start <= end:
            raise ValidationError("wraps_origin=True requires start > end.")
        if not self.wraps_origin and start > end:
            raise ValidationError("start must be <= end unless wraps_origin=True.")
        object.__setattr__(self, "start", start)
        object.__setattr__(self, "end", end)


@dataclass(frozen=True)
class FeatureSpan:
    record: RecordSelector | None
    selectors: tuple[FeatureSelector, ...]
    envelope: Literal["outer_bounds", "segments"] = "outer_bounds"
    circular_path: Literal["shortest", "forward", "reverse"] = "shortest"

    def __post_init__(self) -> None:
        if self.record is not None and not isinstance(self.record, RecordSelector):
            raise ValidationError("FeatureSpan record must be a RecordSelector.")
        selectors = tuple(parse_feature_selector(item) for item in self.selectors)
        if not selectors:
            raise ValidationError("FeatureSpan requires at least one selector.")
        if self.envelope not in {"outer_bounds", "segments"}:
            raise ValidationError("envelope must be 'outer_bounds' or 'segments'.")
        if self.circular_path not in {"shortest", "forward", "reverse"}:
            raise ValidationError("circular_path must be 'shortest', 'forward', or 'reverse'.")
        object.__setattr__(self, "selectors", selectors)


RegionTarget: TypeAlias = CoordinateSpan | FeatureSpan


@dataclass(frozen=True)
class HatchStyle:
    angle: float = 45.0
    spacing: float = 5.0
    color: str = "#666666"
    width: float = 1.0
    cross: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "angle", _finite(self.angle, "hatch angle"))
        object.__setattr__(self, "spacing", _finite(self.spacing, "hatch spacing", positive=True))
        object.__setattr__(self, "width", _finite(self.width, "hatch width", positive=True))
        object.__setattr__(self, "color", _color(self.color, "hatch color"))
        object.__setattr__(self, "cross", bool(self.cross))


@dataclass(frozen=True)
class RegionAnnotationStyle:
    stroke: str = "#404040"
    stroke_width: float = 1.5
    stroke_dasharray: tuple[float, ...] = ()
    line_cap: Literal["none", "tick", "arrow"] = "tick"
    fill: str | None = None
    fill_opacity: float = 0.2
    hatch: HatchStyle | None = None
    label_color: str = "#202020"
    label_font_size: float | None = None
    label_orientation: Literal["auto", "horizontal", "tangent", "radial", "arc"] = "auto"
    label_position: Literal["center", "start", "end"] = "center"
    label_offset: float = 4.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "stroke", _color(self.stroke, "annotation stroke"))
        object.__setattr__(self, "stroke_width", _finite(self.stroke_width, "stroke_width", positive=True))
        dasharray = tuple(_finite(item, "stroke_dasharray", positive=True) for item in self.stroke_dasharray)
        object.__setattr__(self, "stroke_dasharray", dasharray)
        if self.line_cap not in {"none", "tick", "arrow"}:
            raise ValidationError("line_cap must be 'none', 'tick', or 'arrow'.")
        object.__setattr__(self, "fill", _color(self.fill, "annotation fill", optional=True))
        opacity = _finite(self.fill_opacity, "fill_opacity")
        if opacity < 0 or opacity > 1:
            raise ValidationError("fill_opacity must be between 0 and 1.")
        object.__setattr__(self, "fill_opacity", opacity)
        if self.hatch is not None and not isinstance(self.hatch, HatchStyle):
            raise ValidationError("hatch must be HatchStyle or None.")
        object.__setattr__(self, "label_color", _color(self.label_color, "label color"))
        if self.label_font_size is not None:
            object.__setattr__(self, "label_font_size", _finite(self.label_font_size, "label_font_size", positive=True))
        if self.label_orientation not in {"auto", "horizontal", "tangent", "radial", "arc"}:
            raise ValidationError("Unsupported annotation label_orientation.")
        if self.label_position not in {"center", "start", "end"}:
            raise ValidationError("label_position must be 'center', 'start', or 'end'.")
        object.__setattr__(self, "label_offset", _finite(self.label_offset, "label_offset", nonnegative=True))


@dataclass(frozen=True)
class RegionAnnotation:
    id: str
    target: RegionTarget
    label: str = ""
    mark: Literal["line", "bracket", "band"] = "bracket"
    lane: int | None = None
    style: RegionAnnotationStyle | None = None
    legend_label: str | None = None
    metadata: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        annotation_id = str(self.id).strip()
        if not annotation_id:
            raise ValidationError("Annotation id cannot be empty.")
        if not isinstance(self.target, (CoordinateSpan, FeatureSpan)):
            raise ValidationError("Annotation target must be CoordinateSpan or FeatureSpan.")
        if self.mark not in {"line", "bracket", "band"}:
            raise ValidationError("Annotation mark must be 'line', 'bracket', or 'band'.")
        if self.lane is not None:
            if isinstance(self.lane, bool) or int(self.lane) < 0:
                raise ValidationError("Annotation lane must be a non-negative integer.")
            object.__setattr__(self, "lane", int(self.lane))
        if self.style is not None and not isinstance(self.style, RegionAnnotationStyle):
            raise ValidationError("Annotation style must be RegionAnnotationStyle or None.")
        metadata = {str(key): str(value) for key, value in dict(self.metadata or {}).items()}
        object.__setattr__(self, "id", annotation_id)
        object.__setattr__(self, "label", str(self.label or ""))
        object.__setattr__(self, "legend_label", None if self.legend_label is None else str(self.legend_label))
        object.__setattr__(self, "metadata", MappingProxyType(metadata))


@dataclass(frozen=True)
class AnnotationSet:
    id: str
    annotations: tuple[RegionAnnotation, ...]
    default_style: RegionAnnotationStyle = field(default_factory=RegionAnnotationStyle)
    legend_label: str | None = None

    def __post_init__(self) -> None:
        set_id = str(self.id).strip()
        if not set_id:
            raise ValidationError("Annotation set id cannot be empty.")
        annotations = tuple(self.annotations)
        if any(not isinstance(item, RegionAnnotation) for item in annotations):
            raise ValidationError("AnnotationSet contains an unsupported item.")
        ids = [item.id for item in annotations]
        if len(set(ids)) != len(ids):
            duplicate = next(item for item in ids if ids.count(item) > 1)
            raise ValidationError(f"Duplicate annotation id {duplicate!r} in set {set_id!r}.")
        if not isinstance(self.default_style, RegionAnnotationStyle):
            raise ValidationError("AnnotationSet default_style must be RegionAnnotationStyle.")
        object.__setattr__(self, "id", set_id)
        object.__setattr__(self, "annotations", annotations)
        object.__setattr__(self, "legend_label", None if self.legend_label is None else str(self.legend_label))


@dataclass(frozen=True)
class AnnotationOptions:
    sets: tuple[AnnotationSet, ...] = ()
    table: DataFrame | None = None
    table_file: str | None = None

    def __post_init__(self) -> None:
        sets = tuple(self.sets)
        if sum(bool(value) for value in (sets, self.table is not None, self.table_file)) > 1:
            raise ValidationError("Specify annotation sets, table, or table_file, not more than one.")
        ids = [item.id for item in sets]
        if len(set(ids)) != len(ids):
            raise ValidationError("Annotation set ids must be unique.")
        object.__setattr__(self, "sets", sets)


@dataclass(frozen=True)
class AnnotationTrackParams:
    set_id: str
    lane_gap_px: float = 3.0
    padding_px: float = 2.0
    overflow: Literal["error", "compress", "clip"] = "error"
    show_labels: bool = True
    style_override: RegionAnnotationStyle | None = None
    anchor_slot: str | None = None
    layer: Literal["underlay", "foreground"] = "foreground"

    def __post_init__(self) -> None:
        set_id = str(self.set_id).strip()
        if not set_id:
            raise ValidationError("Annotation track set_id is required.")
        object.__setattr__(self, "set_id", set_id)
        object.__setattr__(self, "lane_gap_px", _finite(self.lane_gap_px, "lane_gap_px", nonnegative=True))
        object.__setattr__(self, "padding_px", _finite(self.padding_px, "padding_px", nonnegative=True))
        if self.overflow not in {"error", "compress", "clip"}:
            raise ValidationError("overflow must be 'error', 'compress', or 'clip'.")
        if self.layer not in {"underlay", "foreground"}:
            raise ValidationError("layer must be 'underlay' or 'foreground'.")
        if self.style_override is not None and not isinstance(self.style_override, RegionAnnotationStyle):
            raise ValidationError("style_override must be RegionAnnotationStyle or None.")
        if self.anchor_slot is not None:
            anchor = str(self.anchor_slot).strip()
            if not anchor:
                raise ValidationError("anchor_slot cannot be empty.")
            object.__setattr__(self, "anchor_slot", anchor)


def annotation_track_params_from_mapping(params: Mapping[str, object]) -> AnnotationTrackParams:
    """Normalize raw TrackSlot params into the typed annotation contract."""

    raw = {str(key).strip().lower(): value for key, value in dict(params or {}).items()}
    allowed = {
        "set_id",
        "lane_gap_px",
        "padding_px",
        "overflow",
        "show_labels",
        "style_override",
        "anchor_slot",
        "layer",
    }
    unknown = set(raw) - allowed
    if unknown:
        raise ValidationError(
            f"Unknown annotation track parameter(s): {', '.join(sorted(unknown))}."
        )
    show_labels_raw = raw.get("show_labels", True)
    if isinstance(show_labels_raw, str):
        normalized = show_labels_raw.strip().lower()
        if normalized not in {"1", "0", "true", "false", "yes", "no", "on", "off"}:
            raise ValidationError("show_labels must be true or false.")
        show_labels = normalized in {"1", "true", "yes", "on"}
    else:
        show_labels = bool(show_labels_raw)
    style_override = raw.get("style_override")
    if isinstance(style_override, Mapping):
        try:
            style_override = RegionAnnotationStyle(**dict(style_override))
        except TypeError as exc:
            raise ValidationError(f"Invalid annotation style_override: {exc}") from exc
    return AnnotationTrackParams(
        set_id=str(raw.get("set_id", "")),
        lane_gap_px=raw.get("lane_gap_px", 3.0),  # type: ignore[arg-type]
        padding_px=raw.get("padding_px", 2.0),  # type: ignore[arg-type]
        overflow=str(raw.get("overflow", "error")).strip().lower(),  # type: ignore[arg-type]
        show_labels=show_labels,
        style_override=style_override,  # type: ignore[arg-type]
        anchor_slot=(None if raw.get("anchor_slot") is None else str(raw["anchor_slot"])),
        layer=str(raw.get("layer", "foreground")).strip().lower(),  # type: ignore[arg-type]
    )


@dataclass(frozen=True)
class ResolvedRegionAnnotation:
    id: str
    set_id: str
    record_index: int
    segments: tuple[tuple[int, int], ...]
    midpoint_bp: float
    span_bp: int
    label: str
    mark: str
    lane: int | None
    style: RegionAnnotationStyle
    style_is_annotation_override: bool = False
    legend_label: str | None = None
    metadata: Mapping[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class ResolutionWarning:
    code: str
    set_id: str
    annotation_id: str
    message: str


@dataclass(frozen=True)
class ResolvedAnnotationBundle:
    annotations: tuple[ResolvedRegionAnnotation, ...]
    warnings: tuple[ResolutionWarning, ...] = ()
    set_ids: tuple[str, ...] = ()


def effective_annotation_style(
    annotation: ResolvedRegionAnnotation,
    params: AnnotationTrackParams,
) -> RegionAnnotationStyle:
    """Apply set, track, and annotation style precedence."""

    if annotation.style_is_annotation_override:
        return annotation.style
    return params.style_override or annotation.style


__all__ = [
    "AnnotationOptions",
    "AnnotationSet",
    "AnnotationTrackParams",
    "CoordinateSpan",
    "FeatureSelector",
    "FeatureSpan",
    "HatchStyle",
    "RegionAnnotation",
    "RegionAnnotationStyle",
    "RegionTarget",
    "ResolvedAnnotationBundle",
    "ResolvedRegionAnnotation",
    "ResolutionWarning",
    "annotation_track_params_from_mapping",
    "effective_annotation_style",
    "parse_feature_selector",
]
