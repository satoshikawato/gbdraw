"""Immutable public contracts for region annotation tracks."""

from __future__ import annotations

from gbdraw.layout.record_coordinates import DisplayFragment

from dataclasses import dataclass, field, replace
from math import isfinite
from types import MappingProxyType
from typing import Container, Iterable, Literal, Mapping, NamedTuple, Protocol, Sequence, TypeAlias, TypeVar, cast

from pandas import DataFrame

from gbdraw.core.color import HEX_COLOR_DIGITS, normalize_hex_color
from gbdraw.core.numeric import ascii_float
from gbdraw.exceptions import ValidationError
from gbdraw.io.colors import COLOR_NAMES, resolve_color_to_hex
from gbdraw.io.record_select import RecordSelector


def _finite(value: object, name: str, *, positive: bool = False, nonnegative: bool = False) -> float:
    try:
        # TypeError/ValueError from non-numeric input is reported below.
        number = ascii_float(cast("float | str", value))
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{name} must be numeric.") from exc
    if not isfinite(number):
        raise ValidationError(f"{name} must be finite.")
    if positive and number <= 0:
        raise ValidationError(f"{name} must be positive.")
    if nonnegative and number < 0:
        raise ValidationError(f"{name} must be non-negative.")
    return number


# What `_color` reads, as data (the canonical request description): outer
# spaces, then a hex of 3 or 6 digits (normalize_hex_color) or one of
# COLOR_NAMES in any case (named_color_hex); an optional color also reads a
# blank text as unset. The name pattern spells each letter in both cases, so
# Python's re and the Web's RegExp read it alike without a flag.
HEX_COLOR_TEXT_PATTERN = rf"^\s*#(?:{HEX_COLOR_DIGITS})\s*$"
COLOR_NAME_TEXT_PATTERN = r"^\s*(?:{})\s*$".format(
    "|".join("".join(f"[{letter}{letter.upper()}]" for letter in name) for name in COLOR_NAMES)
)
BLANK_TEXT_PATTERN = r"^\s*$"

# Where a Linear region label stands along its region (``ResolvedRegionAnnotation.label_anchor``).
LabelPosition: TypeAlias = Literal["center", "start", "end"]


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
    except (TypeError, ValueError, ValidationError) as exc:
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
            raise ValidationError(
                "Feature selector value cannot be empty.",
                diagnostic={"code": "INPUT_INVALID", "reason": "REQUIRED"},
            )
        if key == "":
            raise ValidationError(
                "Feature selector key cannot be empty.",
                diagnostic={"code": "INPUT_INVALID", "reason": "REQUIRED"},
            )
        object.__setattr__(self, "value", value)
        object.__setattr__(self, "key", key)


def parse_feature_selector(value: str | FeatureSelector) -> FeatureSelector:
    if isinstance(value, FeatureSelector):
        return value
    text = str(value).strip()
    if not text:
        raise ValidationError(
            "Feature selector cannot be empty.",
            diagnostic={"code": "INPUT_INVALID", "reason": "REQUIRED"},
        )
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
        # An infinite coordinate (`1e400` in JSON) has no integer (OV-3390).
        except (TypeError, ValueError, OverflowError) as exc:
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
        _validate_feature_envelope(self)
        object.__setattr__(self, "selectors", selectors)


def _validate_feature_envelope(target: FeatureSpan | FeatureIdentitySpan) -> None:
    if target.envelope not in {"outer_bounds", "segments"}:
        raise ValidationError("envelope must be 'outer_bounds' or 'segments'.")
    if target.circular_path not in {"shortest", "forward", "reverse"}:
        raise ValidationError("circular_path must be 'shortest', 'forward', or 'reverse'.")


@dataclass(frozen=True)
class FeatureIdentitySpan:
    """One feature by original-source identity, as in feature overrides (design Q4).

    It keeps naming the same feature after crop, reverse complement, reordering
    and record duplication. The typed request planner resolves it.
    """

    record_key: str
    biological_feature_id: str
    envelope: Literal["outer_bounds", "segments"] = "outer_bounds"
    circular_path: Literal["shortest", "forward", "reverse"] = "shortest"

    def __post_init__(self) -> None:
        # Deferred: gbdraw.features.source imports this package.
        from gbdraw.features.source import FeatureIdentity

        identity = FeatureIdentity(self.record_key, self.biological_feature_id)
        object.__setattr__(self, "record_key", identity.record_key)
        object.__setattr__(self, "biological_feature_id", identity.biological_feature_id)
        _validate_feature_envelope(self)


RegionTarget: TypeAlias = CoordinateSpan | FeatureSpan | FeatureIdentitySpan


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
class AnnotationGradient:
    """A fill that runs along a region from ``start_color`` to ``end_color``.

    The colours are the two ends of a scale in the order of what the region
    stands for (a search hit: the query's start, then its end).
    ``query_start`` and ``query_end`` are where the region's two ends lie on
    that scale, as fractions from 0 (the query's start) to 1 (its end), in
    that order: a hit over the second half of its query runs from the
    middle colour to ``end_color``. ``direction`` says which way that order
    runs on the record's source strand: ``forward`` from the region's lower
    source coordinate, ``reverse`` from its higher one. A ``band`` or
    ``highlight`` mark draws it in place of its fill; other marks have no fill.
    """

    start_color: str
    end_color: str
    direction: Literal["forward", "reverse"] = "forward"
    query_start: float = 0.0
    query_end: float = 1.0

    def __post_init__(self) -> None:
        for name in ("start_color", "end_color"):
            text = str(getattr(self, name) or "").strip()
            try:
                color = normalize_hex_color(resolve_color_to_hex(text)) if text else None
            except (TypeError, ValueError, ValidationError):
                color = None
            if color is None:
                raise ValidationError(
                    f"Annotation gradient {name} is not a valid color: {text!r}.",
                    diagnostic={"code": "INPUT_INVALID", "reason": "COLOR"},
                )
            object.__setattr__(self, name, color)
        if self.direction not in {"forward", "reverse"}:
            raise ValidationError(
                "Annotation gradient direction must be 'forward' or 'reverse'.",
                diagnostic={"code": "INPUT_INVALID", "reason": "FIELDS"},
            )
        for name in ("query_start", "query_end"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not isfinite(value):
                raise ValidationError(
                    f"Annotation gradient {name} must be a finite number.",
                    diagnostic={"code": "INPUT_INVALID", "reason": "FINITE"},
                )
            if value < 0:
                raise ValidationError(
                    f"Annotation gradient {name} must be 0 or more.",
                    diagnostic={"code": "INPUT_INVALID", "reason": "NONNEGATIVE"},
                )
            object.__setattr__(self, name, float(value))
        if not self.query_start <= self.query_end <= 1:
            raise ValidationError(
                "Annotation gradient query positions run from query_start to query_end, at most 1.",
                diagnostic={"code": "INPUT_INVALID", "reason": "ORDER"},
            )


# The fill a ``highlight`` draws with when its style has none; the Web starts a new style on it too
# (``describe_draft_annotation_style`` in gbdraw/session_request_codec.py).
DEFAULT_HIGHLIGHT_FILL = "#94a3b8"


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
    # None: the font size of the feature labels of the drawing (annotation_label_font_size).
    label_font_size: float | None = None
    label_orientation: Literal["auto", "horizontal", "tangent", "radial", "arc"] = "auto"
    label_position: LabelPosition = "center"
    label_offset: float = 4.0
    # A fill along the region in place of ``fill`` (band and highlight marks).
    gradient: AnnotationGradient | None = None

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
        if self.gradient is not None and not isinstance(self.gradient, AnnotationGradient):
            raise ValidationError(
                "gradient must be AnnotationGradient or None.", diagnostic={"code": "INPUT_INVALID", "reason": "FIELDS"}
            )
        if self.gradient is not None and self.hatch is not None:
            raise ValidationError(
                "An annotation style has a hatch or a gradient, not both.",
                diagnostic={"code": "INPUT_INVALID", "reason": "FIELDS"},
            )
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
    mark: Literal["line", "bracket", "band", "highlight"] = "bracket"
    lane: int | None = None
    style: RegionAnnotationStyle | None = None
    legend_label: str | None = None
    metadata: Mapping[str, str] = field(default_factory=dict)
    # The strand the region lies on, on the record's source (a sequence search's hit, PD-OI-117): on a feature
    # track that separates the strands its highlight covers that strand's lane. None: no strand.
    strand: Literal["+", "-"] | None = None

    def __post_init__(self) -> None:
        annotation_id = str(self.id).strip()
        if not annotation_id:
            raise ValidationError("Annotation id cannot be empty.")
        if not isinstance(self.target, (CoordinateSpan, FeatureSpan, FeatureIdentitySpan)):
            raise ValidationError(
                "Annotation target must be CoordinateSpan, FeatureSpan or FeatureIdentitySpan."
            )
        if self.mark not in {"line", "bracket", "band", "highlight"}:
            raise ValidationError(
                "Annotation mark must be 'line', 'bracket', 'band', or 'highlight'."
            )
        if self.lane is not None:
            try:
                lane = None if isinstance(self.lane, bool) else int(self.lane)
            except OverflowError as exc:  # an infinite lane (`1e400` in JSON) has no integer (OV-3390)
                raise ValidationError("Annotation lane must be a non-negative integer.") from exc
            if lane is None or lane < 0:
                raise ValidationError("Annotation lane must be a non-negative integer.")
            object.__setattr__(self, "lane", lane)
        if self.style is not None and not isinstance(self.style, RegionAnnotationStyle):
            raise ValidationError("Annotation style must be RegionAnnotationStyle or None.")
        if self.strand not in {None, "+", "-"}:
            raise ValidationError(
                "Annotation strand must be '+', '-', or None.", diagnostic={"code": "INPUT_INVALID", "reason": "FIELDS"}
            )
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
    # A sequence search's set (PD-OI-117): its highlights keep the automatic overlay, in front of the features, also
    # in an explicit track stack that names no slot for the set (``prepare_annotation_track_slots``).
    keep_highlight_overlay: bool = False

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
        if not isinstance(self.keep_highlight_overlay, bool):
            raise ValidationError(
                "AnnotationSet keep_highlight_overlay must be a boolean.",
                diagnostic={"code": "INPUT_INVALID", "reason": "FIELDS"},
            )
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
    marks: tuple[str, ...] | None = None
    lane_gap_px: float = 3.0
    padding_px: float = 2.0
    overflow: Literal["error", "compress", "clip"] = "error"
    show_labels: bool = True
    style_override: RegionAnnotationStyle | None = None
    anchor_slot: str | None = None
    layer: Literal["underlay", "foreground"] = "foreground"
    cover_anchor: bool = False
    # With to_slot, an overlay spans from its anchor slot's band to this slot's band.
    to_slot: str | None = None

    def __post_init__(self) -> None:
        set_id = str(self.set_id).strip()
        if not set_id:
            raise ValidationError("Annotation track set_id is required.")
        object.__setattr__(self, "set_id", set_id)
        if self.marks is not None:
            raw_marks = (
                self.marks.replace(",", "|").split("|")
                if isinstance(self.marks, str)
                else self.marks
            )
            marks = tuple(dict.fromkeys(str(mark).strip().lower() for mark in raw_marks if str(mark).strip()))
            unknown_marks = set(marks) - {"line", "bracket", "band", "highlight"}
            if unknown_marks:
                raise ValidationError(
                    f"Unsupported annotation mark filter(s): {', '.join(sorted(unknown_marks))}."
                )
            if not marks:
                raise ValidationError("Annotation track marks filter cannot be empty.")
            object.__setattr__(self, "marks", marks)
        object.__setattr__(self, "lane_gap_px", _finite(self.lane_gap_px, "lane_gap_px", nonnegative=True))
        object.__setattr__(self, "padding_px", _finite(self.padding_px, "padding_px", nonnegative=True))
        if self.overflow not in {"error", "compress", "clip"}:
            raise ValidationError("overflow must be 'error', 'compress', or 'clip'.")
        if self.layer not in {"underlay", "foreground"}:
            raise ValidationError("layer must be 'underlay' or 'foreground'.")
        if self.style_override is not None and not isinstance(self.style_override, RegionAnnotationStyle):
            raise ValidationError("style_override must be RegionAnnotationStyle or None.")
        for name in ("anchor_slot", "to_slot"):
            reference = getattr(self, name)
            if reference is None:
                continue
            slot_id = str(reference).strip()
            if not slot_id:
                raise ValidationError(f"{name} cannot be empty.")
            object.__setattr__(self, name, slot_id)


def overlay_labels_clear_anchor(
    side: str, params: AnnotationTrackParams, feature_slot_ids: Container[str]
) -> str | None:
    """The features slot whose glyphs an annotation slot's labels keep off, in either mode: the anchor of an overlay
    that spans nothing (no ``to_slot``, no ``cover_anchor``) when ``feature_slot_ids`` holds it (OV-3078, OV-3219,
    OV-3220); None for every other slot. Each mode places the labels clear of that slot's lanes."""

    if side != "overlay" or params.to_slot is not None or params.cover_anchor:
        return None
    return params.anchor_slot if params.anchor_slot is not None and params.anchor_slot in feature_slot_ids else None


class LinearAnnotationLabelSide(NamedTuple):
    """Where a Linear annotation slot puts its labels: above its marks (else below them), and the features slot whose
    lanes they stand beyond, on that side (None: beside their marks)."""

    above: bool
    beyond_slot: str | None = None


def linear_annotation_label_side(
    side: str, params: AnnotationTrackParams, feature_slot_sides: Mapping[str, str]
) -> LinearAnnotationLabelSide:
    """Where a Linear annotation slot on ``side`` puts its labels; ``feature_slot_sides`` is the side of each drawn
    features slot, by slot id.

    An overlay anchored on a features slot drawn on the axis that spans no other slot (no ``to_slot``) puts its
    labels below the lanes on either layer, where feature labels do not stand (OV-3223): beyond them when it keeps
    off its anchor's glyphs (``overlay_labels_clear_anchor``), and below its covering band with ``cover_anchor``.
    An overlay that keeps off the glyphs of lanes drawn above or below the axis stands its labels on its marks'
    other side, away from the lanes. Every other slot puts its labels above its marks on the ``above`` side and on
    an underlay overlay, below them elsewhere.
    """

    anchor = overlay_labels_clear_anchor(side, params, feature_slot_sides)
    if (
        side == "overlay"
        and params.to_slot is None
        and params.anchor_slot is not None
        and feature_slot_sides.get(params.anchor_slot) == "overlay"
    ):
        return LinearAnnotationLabelSide(False, anchor)
    if anchor is not None:
        return LinearAnnotationLabelSide(feature_slot_sides[anchor] == "below")
    return LinearAnnotationLabelSide(side == "above" or (side == "overlay" and params.layer == "underlay"))


def annotation_track_params_from_mapping(params: Mapping[str, object]) -> AnnotationTrackParams:
    """Normalize raw TrackSlot params into the typed annotation contract."""

    raw = {str(key).strip().lower(): value for key, value in dict(params or {}).items()}
    allowed = {
        "set_id",
        "marks",
        "lane_gap_px",
        "padding_px",
        "overflow",
        "show_labels",
        "style_override",
        "anchor_slot",
        "layer",
        "cover_anchor",
        "to_slot",
    }
    unknown = set(raw) - allowed
    if unknown:
        raise ValidationError(
            f"Unknown annotation track parameter(s): {', '.join(sorted(unknown))}."
        )
    def parse_bool(raw_value: object, field_name: str) -> bool:
        if isinstance(raw_value, str):
            normalized = raw_value.strip().lower()
            if normalized not in {"1", "0", "true", "false", "yes", "no", "on", "off"}:
                raise ValidationError(f"{field_name} must be true or false.")
            return normalized in {"1", "true", "yes", "on"}
        return bool(raw_value)

    show_labels = parse_bool(raw.get("show_labels", True), "show_labels")
    cover_anchor = parse_bool(raw.get("cover_anchor", False), "cover_anchor")
    style_override = raw.get("style_override")
    if isinstance(style_override, Mapping):
        try:
            style_override = RegionAnnotationStyle(**dict(style_override))
        except TypeError as exc:
            raise ValidationError(f"Invalid annotation style_override: {exc}") from exc
    # The casts below only satisfy the checker: AnnotationTrackParams.__post_init__
    # validates every one of these values at run time.
    return AnnotationTrackParams(
        set_id=str(raw.get("set_id", "")),
        marks=cast("tuple[str, ...] | None", raw.get("marks")),
        lane_gap_px=cast(float, raw.get("lane_gap_px", 3.0)),
        padding_px=cast(float, raw.get("padding_px", 2.0)),
        overflow=cast(
            'Literal["error", "compress", "clip"]',
            str(raw.get("overflow", "error")).strip().lower(),
        ),
        show_labels=show_labels,
        style_override=cast("RegionAnnotationStyle | None", style_override),
        anchor_slot=(None if raw.get("anchor_slot") is None else str(raw["anchor_slot"])),
        layer=cast(
            'Literal["underlay", "foreground"]',
            str(raw.get("layer", "foreground")).strip().lower(),
        ),
        cover_anchor=cover_anchor,
        to_slot=(None if raw.get("to_slot") is None else str(raw["to_slot"])),
    )


class PaintedSlot(Protocol):
    """A resolved track slot as the paint order reads it, in either mode."""

    @property
    def slot_index(self) -> int: ...

    @property
    def renderer(self) -> str: ...

    @property
    def side(self) -> str: ...

    @property
    def z(self) -> int: ...

    @property
    def params(self) -> Mapping[str, object]: ...


PaintedSlotT = TypeVar("PaintedSlotT", bound=PaintedSlot)


def slots_in_paint_order(slots: Iterable[PaintedSlotT]) -> list[PaintedSlotT]:
    """The order a drawing paints its slots in, the same in Circular and Linear.

    Slots paint by ``z``. Within one ``z``, underlay overlays paint first, then
    the other slots, then foreground overlays, so an overlay lies on its layer's
    side of every slot of its ``z`` it overlaps, whatever order the slots were
    written in. Slots of one rank keep their written order (OV-3184, OV-3185).
    """

    def rank(slot: PaintedSlotT) -> int:
        if str(slot.renderer) != "annotations" or str(slot.side) != "overlay":
            return 1
        return 0 if annotation_track_params_from_mapping(slot.params).layer == "underlay" else 2

    return sorted(slots, key=lambda slot: (int(slot.z), rank(slot), int(slot.slot_index)))


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
    display_parts: tuple[DisplayFragment, ...] | None = None
    # Set with ``display_parts`` from the record's transform (``RecordDisplayTransform.circle``): the record is one of
    # several drawn on one circle, so its arc ends at its display cut, where the gap to the next record starts.
    shares_circle: bool = False
    # The span goes on across the origin of a circular record: from its last segment, which ends at the record's
    # end, into its first, which starts at 0 (``segments_cross_origin``).
    crosses_origin: bool = False

    # How the record's local axis runs on its source: -1 on a reverse-complemented record.
    source_step: Literal[-1, 1] = 1
    # ``RegionAnnotation.strand``, on the record's source.
    strand: Literal["+", "-"] | None = None

    @property
    def display_strand(self) -> Literal["+", "-"] | None:
        """The strand the region lies on as the record is drawn: a reverse-complemented record draws its
        source's plus strand as its minus strand, as it does its features."""

        if self.strand is None or self.source_step == 1:
            return self.strand
        return "-" if self.strand == "+" else "+"

    def __post_init__(self) -> None:
        if self.crosses_origin and (len(self.segments) < 2 or not any(start == 0 for start, _end in self.segments)):
            raise ValueError(
                f"Annotation {self.id!r} crosses the origin, so it needs at least two segments and one that "
                f"starts at 0; got {self.segments!r}."
            )

    @property
    def geometry_segments(self) -> tuple[tuple[int, int], ...]:
        if self.display_parts is None:
            return self.segments
        return tuple((part.display_start, part.display_end) for part in self.display_parts)

    def label_place(self, *, join_origin: bool) -> LabelPlace:
        """Where the region's centred label stands (display bp), and the drawn part it stands on.

        When its mark draws every part whole (``mark_parts``: each part capped at both ends, so no display cut,
        row boundary or unjoined origin cuts it), the label stands halfway along the region (``midpoint_bp``)
        as the display shows it. Else it stands on the centre of its longest part (PD-OI-119): where only the
        origin cuts it (no display seam), two equal parts go to the one where the region starts.
        """
        parts = mark_parts(self, join_origin=join_origin)
        if not all(part.capped_start and part.capped_end for part in parts):
            seamless = self.display_parts is None or not any(
                part.artificial_start or part.artificial_end for part in self.display_parts
            )
            longest = _longest_part(
                parts, prefer=self.label_anchor("start", join_origin=join_origin).part if seamless else None
            )
            return LabelPlace((longest.start + longest.end) / 2.0, longest)
        longest = _longest_part(parts)
        midpoint = self.midpoint_bp
        if self.display_parts is None:
            return LabelPlace(midpoint, next((part for part in parts if part.start <= midpoint <= part.end), longest))
        fragment = next(part for part in self.display_parts if part.local_start <= midpoint <= part.local_end)
        # A display fragment runs the same way as its record-local span.
        bp = fragment.display_start + (midpoint - fragment.local_start)
        return LabelPlace(
            bp, next((part for part in parts if part.row == fragment.row and part.start <= bp <= part.end), longest)
        )

    def label_midpoint_bp(self, *, join_origin: bool) -> float:
        """Where the region's centred label stands (``label_place``)."""
        return self.label_place(join_origin=join_origin).bp

    def label_anchor(self, position: LabelPosition, *, join_origin: bool) -> LabelPlace:
        """Where a Linear label at ``position`` stands and the drawn part it stands on: centred at
        ``label_place``; at the start, the region's first base, at the end its last, as the region runs
        (``region_run``) and as the display shows them (PD-OI-119): on the drawn part (``mark_parts``) that holds
        that base, also where a circle joins it to another part."""

        if position == "center":
            return self.label_place(join_origin=join_origin)
        run = region_run(self.segments, crosses_origin=self.crosses_origin)
        first = position == "start"
        local = run[0][0] if first else run[-1][1]
        # A part runs on from its record-local start; one joined over the origin runs on past the record's end.
        record_end = max(end for _start, end in self.segments) if self.crosses_origin else 0
        for part in mark_parts(self, join_origin=join_origin):
            offset = local - part.local_start
            offset += record_end if offset < 0 else 0
            if (0 <= offset < part.end - part.start) if first else (0 < offset <= part.end - part.start):
                return LabelPlace(float(part.start + offset), part)
        raise AssertionError(f"Annotation {self.id!r}: no drawn part holds its {position} base {local}.")

    def label_stretch(self, position: LabelPosition, *, join_origin: bool) -> tuple[LabelPlace, tuple[int, int]]:
        """Where a label at ``position`` stands (``label_anchor``) and the stretch of its row it may cover without
        leaving its mark: the region's drawn extent on that row where the mark draws every part whole, else the
        part it stands on."""
        place = self.label_anchor(position, join_origin=join_origin)
        parts = mark_parts(self, join_origin=join_origin)
        if all(part.capped_start and part.capped_end for part in parts):
            row = [part for part in parts if part.row == place.part.row]
            return place, (min(part.start for part in row), max(part.end for part in row))
        return place, (place.part.start, place.part.end)

    @property
    def geometry_rows(self) -> tuple[int, ...]:
        """The rows of a wrapped record the annotation draws on, in order."""
        if self.display_parts is None:
            return (0,)
        return tuple(sorted({part.row for part in self.display_parts}))

    def row_segments(self, row: int) -> tuple[tuple[int, int], ...]:
        """The drawn segments on one row, local to it."""
        if self.display_parts is None:
            return self.segments
        return tuple((part.display_start, part.display_end) for part in self.display_parts if part.row == row)


@dataclass(frozen=True)
class MarkPart:
    """One drawn part of a mark, in display bp on its ``row``, and whether it draws each of its ends.

    An end it does not draw is a cut, where the region goes on: a display cut that ends the record's drawing (a
    Linear bar, a record that shares one circle), a row boundary, or the origin of a span across it. A cut gets no
    edge stroke, cap or arrow head. A part joined across the origin or a circle's display cut ends past the
    record's length. ``local_start`` is where the part starts on its record (record-local bp); the part runs on
    from there as it runs on the display, past the record's end when it is joined across the origin.
    """

    start: int
    end: int
    capped_start: bool
    capped_end: bool
    local_start: int
    row: int = 0


class LabelPlace(NamedTuple):
    """Where a region's label stands (display bp) and the drawn part it stands on."""

    bp: float
    part: MarkPart


def _longest_part(parts: Sequence[MarkPart], *, prefer: MarkPart | None = None) -> MarkPart:
    """The longest drawn part; ties to ``prefer``, then to the earlier row, then the earlier start."""

    return min(parts, key=lambda part: (part.start - part.end, part != prefer, part.row, part.start))


def origin_sides(bounds: Sequence[tuple[int, int]]) -> tuple[int, int]:
    """The two sides at the origin of a span across it (``crosses_origin``), as indices into its record-local
    ``bounds``: the one that ends at the record's end, then the one that starts at 0."""

    # The resolver sets crosses_origin only when a segment starts at 0 and the last one ends at the record's end.
    record_end = max(end for _start, end in bounds)
    before = next(index for index, (_start, end) in enumerate(bounds) if end == record_end)
    after = next(index for index, (start, _end) in enumerate(bounds) if start == 0)
    return before, after


def region_run(segments: Sequence[tuple[int, int]], *, crosses_origin: bool) -> tuple[tuple[int, int], ...]:
    """A span's segments in the order the region runs (PD-OI-119): ascending, or across the origin from its side
    that ends at the record's end over the origin into its side that starts at 0 (``origin_sides``), then its
    other segments in ascending order. Its midpoint, its gradient and its start and end labels read this run."""

    ordered = tuple(sorted(segments))
    if not crosses_origin:
        return ordered
    before, after = origin_sides(ordered)
    return (ordered[before], ordered[after], *(item for index, item in enumerate(ordered) if index not in {before, after}))


def _joined(head: MarkPart, tail: MarkPart) -> MarkPart:
    """``head`` going on into ``tail`` where the drawing meets them: one part, ``tail``'s length past ``head``'s end."""

    return replace(head, end=head.end + tail.end - tail.start, capped_end=tail.capped_end)


def mark_parts(annotation: ResolvedRegionAnnotation, *, join_origin: bool) -> tuple[MarkPart, ...]:
    """The parts a mark draws, in the order of its geometry segments: the one owner of where a mark is cut, which
    its drawing and its label read.

    A record that closes its own circle (``join_origin``, a Circular drawing, and not ``shares_circle``) has no end
    at its display cut, the start of a rotated Circular drawing: a part that ends there goes on into the part that
    starts there with its next record-local base, as one part. A Linear record ends at its display cut, and so does
    a record that shares one circle with other records (the gap to the next record starts there), so their parts
    stay cut there.

    A span across the origin (``crosses_origin``) has two sides there: its part that ends at the record's end and
    its part that starts at 0. They become one part where they meet in the drawing: always on a record that closes
    its own circle, else where one ends on the row where the other starts (a record drawn from a display start).
    Where they do not meet, they leave those ends open.
    """

    closes_circle = join_origin and not annotation.shares_circle

    fragments = annotation.display_parts
    parts = [
        MarkPart(start, end, True, True, start)
        for start, end in annotation.segments
    ] if fragments is None else [
        MarkPart(
            part.display_start, part.display_end, not part.artificial_start, not part.artificial_end, part.local_start,
            part.row,
        )
        for part in fragments
    ]
    if closes_circle and fragments is not None:
        # On a circle every artificial seam is the display cut (no rows): join each part cut there to the part
        # that starts at its next base.
        tails = {
            head: next(
                index for index, other in enumerate(fragments)
                if other.artificial_start and other.local_start == fragment.local_end
            )
            for head, fragment in enumerate(fragments)
            if fragment.artificial_end
        }
        parts = [
            _joined(part, parts[tails[index]]) if index in tails else part
            for index, part in enumerate(parts)
            if index not in tails.values()
        ]
    if not annotation.crosses_origin:
        return tuple(parts)
    before, after = origin_sides([(part.local_start, part.local_start + part.end - part.start) for part in parts])
    if closes_circle or (parts[before].row == parts[after].row and parts[before].end == parts[after].start):
        parts[before] = _joined(parts[before], parts[after])
        del parts[after]
        return tuple(parts)
    parts[before] = replace(parts[before], capped_end=False)
    parts[after] = replace(parts[after], capped_start=False)
    return tuple(parts)


# Why a region annotation is not drawn whole; a Session stores it in ``runMetadata.annotationWarnings``.
AnnotationWarningCode = Literal[
    "feature_selector_unmatched", "empty_span", "out_of_bounds_skipped", "out_of_bounds_clipped"
]


@dataclass(frozen=True)
class ResolutionWarning:
    code: AnnotationWarningCode
    set_id: str
    annotation_id: str
    message: str
    record_id: str = ""
    record_index: int = 0
    missing_count: int = 0


@dataclass(frozen=True)
class ResolvedAnnotationBundle:
    annotations: tuple[ResolvedRegionAnnotation, ...]
    warnings: tuple[ResolutionWarning, ...] = ()
    set_ids: tuple[str, ...] = ()
    # The sets whose highlights keep their automatic overlay in an explicit track stack
    # (``AnnotationSet.keep_highlight_overlay``).
    highlight_overlay_set_ids: tuple[str, ...] = ()


def effective_annotation_style(
    annotation: ResolvedRegionAnnotation,
    params: AnnotationTrackParams,
) -> RegionAnnotationStyle:
    """Apply set, track, and annotation style precedence."""

    if annotation.style_is_annotation_override:
        return annotation.style
    return params.style_override or annotation.style


def annotation_solid_fill(mark: str, style: RegionAnnotationStyle) -> str:
    """The fill of a ``band`` or ``highlight`` mark without hatching: its style's own, else
    ``DEFAULT_HIGHLIGHT_FILL`` for a highlight and ``none`` for a band."""

    return style.fill or (DEFAULT_HIGHLIGHT_FILL if mark == "highlight" else "none")


def annotation_label_font_size(style: RegionAnnotationStyle, feature_label_font_size: float) -> float:
    """A region label's font size: its style's own, else ``feature_label_font_size``, the font size the
    feature labels of its drawing use."""

    return float(feature_label_font_size if style.label_font_size is None else style.label_font_size)


__all__ = [
    "AnnotationGradient",
    "AnnotationOptions",
    "AnnotationSet",
    "AnnotationTrackParams",
    "CoordinateSpan",
    "DEFAULT_HIGHLIGHT_FILL",
    "FeatureSelector",
    "FeatureSpan",
    "HatchStyle",
    "LabelPlace",
    "MarkPart",
    "RegionAnnotation",
    "RegionAnnotationStyle",
    "RegionTarget",
    "ResolvedAnnotationBundle",
    "ResolvedRegionAnnotation",
    "AnnotationWarningCode",
    "ResolutionWarning",
    "annotation_label_font_size",
    "annotation_solid_fill",
    "annotation_track_params_from_mapping",
    "effective_annotation_style",
    "LinearAnnotationLabelSide",
    "linear_annotation_label_side",
    "overlay_labels_clear_anchor",
    "mark_parts",
    "origin_sides",
    "parse_feature_selector",
    "region_run",
    "slots_in_paint_order",
]
