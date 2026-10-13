"""Whether drawn labels fit, the one decision point, and the report type.

Labels fit when no two labels overlap and no two leader polylines cross
(LP-D43 ``DEF-CROSS``); a label over the centre definition text counts like a
crossing (LP-D49 ``DEF-WARN``). The counts take plain shapes: each label kind
has one adapter next to its owner that turns the drawn labels into
:class:`LabelShape` values, and one record's labels of every kind are counted
together (LP-D48 ``SCOPE-ALL-LABELS``).
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Literal, Protocol, TypeAlias

from ..exceptions import LabelLeadersCrossWarning, LabelsDoNotFitError
from ..layout.spatial import Aabb, candidate_aabb_pairs
from ..layout.text_geometry import Point, convex_polygons_intersect, segments_intersect

LabelKind: TypeAlias = Literal["circular_outer", "circular_inner", "circular_embedded", "linear"]
LabelFitOutcome: TypeAlias = Literal["fits", "warn", "stop"]
LabelFitRowKind: TypeAlias = Literal[
    "outward_both",
    "outward_vertical",
    "outward_horizontal",
    "smaller_font",
    "tighter_rows",
    "combination",
    "labels_out",
    "radial",
    "separate_features",
]
# ``found``: ``values`` fit; ``not_found``: no value within the lever's limit
# (``reason`` says which); ``not_computed``: the bounded search stopped first.
LabelFitRowStatus: TypeAlias = Literal["found", "not_found", "not_computed"]
# A row sets numbers (offsets, font size, spacing), words (label scope,
# placement) and switches (separate overlapping features).
LabelFitValue: TypeAlias = bool | float | str
Polygon: TypeAlias = tuple[Point, ...]
Polyline: TypeAlias = Sequence[Point]

# A pair of labels of two kinds counts once, for the first kind in this order.
_KIND_ORDER: tuple[LabelKind, ...] = ("circular_inner", "circular_outer", "circular_embedded", "linear")
# A label covers the definition text when they overlap by more than this in x and y.
DEFINITION_TEXT_MIN_OVERLAP_PX = 0.5

# The settings the rows change, by config path (the identity the CLI, the Web
# request mapping and Run Info share).
X_OFFSET_PATH = "labels.unified_adjustment.outer_labels.x_radius_offset"
Y_OFFSET_PATH = "labels.unified_adjustment.outer_labels.y_radius_offset"
# ``--label_font_size`` sets both size classes; a row writes them the same way.
FONT_SIZE_PATHS = ("labels.font_size.short", "labels.font_size.long")
SPACING_PATH = "labels.spacing.circular"
SCOPE_PATH = "labels.circular.scope"
PLACEMENT_PATH = "labels.circular.placement"
SEPARATE_FEATURES_PATH = "canvas.resolve_overlaps"
# The switch rows' values: every label outside (``--labels out``), radial
# labels, and overlapping features on separate tracks (``--resolve_overlaps``).
ALL_OUTSIDE_SCOPE = "outer"
RADIAL_PLACEMENT = "radial"
SEPARATE_FEATURES_VALUES: Mapping[str, LabelFitValue] = {SEPARATE_FEATURES_PATH: True}

# How long the rows are searched: a figure drawn anyway (LP-D48
# SEARCH-BOUNDED-5S), and, when nothing was drawn, the switches and the font
# scan of inner and radial labels together (SWITCHES-FIRST, OD-F15).
WARN_SEARCH_SECONDS = 5.0
SWITCH_SEARCH_SECONDS = 30.0


@dataclass(frozen=True)
class LabelShape:
    """One drawn label as the fit count sees it.

    ``parts`` are convex polygons in the record's drawing coordinates whose
    union is the label (its margins included); ``leader`` is the drawn leader
    polyline, empty when the label has none.
    """

    kind: LabelKind
    parts: tuple[Polygon, ...]
    leader: tuple[Point, ...] = ()


@dataclass(frozen=True)
class LabelFit:
    """Overlapping label pairs, crossing leader pairs and labels over the definition text.

    The labels fit when all three are 0.
    """

    overlap_pairs: int
    crossing_pairs: int
    definition_overlaps: int = 0

    @property
    def fits(self) -> bool:
        return self.overlap_pairs == 0 and self.crossing_pairs == 0 and self.definition_overlaps == 0

    def __add__(self, other: LabelFit) -> LabelFit:
        return LabelFit(
            self.overlap_pairs + other.overlap_pairs,
            self.crossing_pairs + other.crossing_pairs,
            self.definition_overlaps + other.definition_overlaps,
        )


FITS = LabelFit(0, 0)


def _parts_aabb(parts: Sequence[Polygon]) -> Aabb:
    return Aabb(
        min(point[0] for part in parts for point in part),
        min(point[1] for part in parts for point in part),
        max(point[0] for part in parts for point in part),
        max(point[1] for part in parts for point in part),
    )


def _closed_polygons_overlap(first: Polygon, second: Polygon) -> bool:
    """Separating-axis test that counts touching polygons; exact on axis-aligned boxes."""

    for polygon in (first, second):
        for index, point in enumerate(polygon):
            following = polygon[(index + 1) % len(polygon)]
            axis = (point[1] - following[1], following[0] - point[0])
            if axis == (0.0, 0.0):
                continue
            if axis[0] == 0.0 or axis[1] == 0.0:
                # Unit axes keep box edges exact (the margins decide touching boxes).
                axis = (
                    math.copysign(1.0, axis[0]) if axis[0] else 0.0,
                    math.copysign(1.0, axis[1]) if axis[1] else 0.0,
                )
            first_projection = [x * axis[0] + y * axis[1] for x, y in first]
            second_projection = [x * axis[0] + y * axis[1] for x, y in second]
            if max(first_projection) < min(second_projection) or max(second_projection) < min(
                first_projection
            ):
                return False
    return True


def _open_polygons_overlap(first: Polygon, second: Polygon) -> bool:
    """The rotated-text test: touching does not count (the radial placement's own rule)."""

    first_box = _parts_aabb((first,))
    second_box = _parts_aabb((second,))
    epsilon = 1e-9
    if (
        first_box.max_x <= second_box.min_x + epsilon
        or second_box.max_x <= first_box.min_x + epsilon
        or first_box.max_y <= second_box.min_y + epsilon
        or second_box.max_y <= first_box.min_y + epsilon
    ):
        return False
    return convex_polygons_intersect(first, second, touching=False)


def overlapping_pairs(
    shapes: Sequence[Sequence[Polygon]],
    *,
    touching: bool,
) -> list[tuple[int, int]]:
    """Index pairs of shapes (each a union of convex parts) that overlap, in index order.

    ``touching`` counts shapes that only touch (horizontal labels, whose
    boxes carry their margins); without it they do not (rotated labels).
    """

    if not shapes:
        return []
    boxes = [_parts_aabb(parts) for parts in shapes]
    extent = max(max(box.width, box.max_y - box.min_y) for box in boxes)
    overlap = _closed_polygons_overlap if touching else _open_polygons_overlap
    return [
        (left, right)
        for left, right in sorted(
            candidate_aabb_pairs(boxes, bucket_size=max(16.0, min(256.0, extent)))
        )
        if any(overlap(first, second) for first in shapes[left] for second in shapes[right])
    ]


def _polyline_aabb(polyline: Polyline) -> Aabb:
    return Aabb(
        min(point[0] for point in polyline),
        min(point[1] for point in polyline),
        max(point[0] for point in polyline),
        max(point[1] for point in polyline),
    )


def crossing_pairs(polylines: Sequence[Polyline]) -> list[tuple[int, int]]:
    """Index pairs of polylines that cross (touching does not count); empty polylines never cross."""

    indexed = [index for index, polyline in enumerate(polylines) if len(polyline) >= 2]
    segments = [list(zip(polylines[index], polylines[index][1:])) for index in indexed]
    pairs: list[tuple[int, int]] = []
    for left, right in sorted(
        candidate_aabb_pairs([_polyline_aabb(polylines[index]) for index in indexed], bucket_size=64.0)
    ):
        if any(
            segments_intersect(first, second, touching=False)
            for first in segments[left]
            for second in segments[right]
        ):
            pairs.append((indexed[left], indexed[right]))
    return pairs


def polyline_crossing_count(polylines: Sequence[Polyline]) -> int:
    """Count polyline pairs that cross (touching does not count)."""

    return len(crossing_pairs(polylines))


def _covers_text(parts: Sequence[Polygon], text_box: Aabb) -> bool:
    label_box = _parts_aabb(parts)
    return (
        min(label_box.max_x, text_box.max_x) - max(label_box.min_x, text_box.min_x)
        > DEFINITION_TEXT_MIN_OVERLAP_PX
        and min(label_box.max_y, text_box.max_y) - max(label_box.min_y, text_box.min_y)
        > DEFINITION_TEXT_MIN_OVERLAP_PX
    )


def label_fits_by_kind(
    shapes: Sequence[LabelShape],
    *,
    touching: bool = True,
    definition_text: Sequence[Aabb] = (),
) -> dict[LabelKind, LabelFit]:
    """The fit of one record's labels of every kind together, split by kind.

    A pair of labels of two kinds counts once, for the kind that comes first
    in inner, outer, embedded, Linear order. Labels with a leader are tested
    against the ``definition_text`` boxes (embedded labels sit on their
    features). Every kind present has an entry.
    """

    counts: dict[LabelKind, list[int]] = {}
    for shape in shapes:
        counts.setdefault(shape.kind, [0, 0, 0])

    def owner(left: int, right: int) -> LabelKind:
        kinds = (shapes[left].kind, shapes[right].kind)
        return next(kind for kind in _KIND_ORDER if kind in kinds)

    for left, right in overlapping_pairs([shape.parts for shape in shapes], touching=touching):
        counts[owner(left, right)][0] += 1
    for left, right in crossing_pairs([shape.leader for shape in shapes]):
        counts[owner(left, right)][1] += 1
    if definition_text:
        for shape in shapes:
            if shape.leader and any(_covers_text(shape.parts, box) for box in definition_text):
                counts[shape.kind][2] += 1
    return {
        kind: LabelFit(overlaps, crossings, covered)
        for kind, (overlaps, crossings, covered) in counts.items()
    }


def label_fit_outcome(fits: Sequence[LabelFit], *, refused: bool = False) -> LabelFitOutcome:
    """The one decision point: overlaps stop; crossings and covered definition text alone warn.

    ``refused`` is a placement that could not place its labels at all (radial
    labels without room, LP-D49 ``RADIAL-JOIN``); it stops like an overlap.
    """

    if refused or any(fit.overlap_pairs > 0 for fit in fits):
        return "stop"
    if any(not fit.fits for fit in fits):
        return "warn"
    return "fits"


@dataclass(frozen=True)
class LabelFitRecord:
    """One record's labels of one kind that do not fit.

    ``refused`` marks radial labels that could not be placed in the room the
    tracks leave (their counts may then be 0).
    """

    record_id: str
    label_kind: LabelKind
    label_count: int
    font_size: float
    fit: LabelFit
    refused: bool = False


@dataclass(frozen=True)
class LabelFitRow:
    """One lever: its mildest fitting value, or why there is none.

    ``values`` maps config paths to the values to set and ``canvas_size`` is
    the (width, height) of the verified figure; both are set only when
    ``status`` is ``found``.
    """

    kind: LabelFitRowKind
    status: LabelFitRowStatus
    values: Mapping[str, LabelFitValue] | None = None
    canvas_size: tuple[float, float] | None = None
    reason: str | None = None


@dataclass(frozen=True)
class LabelFitVerification:
    """The full assembly with a row's settings.

    ``fits`` when every record fits (no overlaps, crossings or covered
    definition text); ``canvas_size`` is the figure's size, ``fit`` its
    labels' fit over every record, and ``reason`` why the assembly failed
    when it raised.
    """

    fits: bool
    canvas_size: tuple[float, float] | None
    fit: LabelFit | None = None
    reason: str | None = None


def switch_row(kind: LabelFitRowKind, values: Mapping[str, LabelFitValue], verification: LabelFitVerification) -> LabelFitRow:
    """A row with one setting to turn on, from its full assembly."""

    if verification.fits:
        return LabelFitRow(kind, "found", values, verification.canvas_size)
    if verification.reason is not None:
        reason = verification.reason
    elif verification.fit is not None and not verification.fit.fits:
        reason = f"Labels still do not fit ({fit_counts(verification.fit)})"
    else:
        reason = "Labels still do not fit"
    return LabelFitRow(kind, "not_found", reason=reason)


@dataclass(frozen=True)
class LabelFitReport:
    """Failing records and the rows that apply to them, in the dialog order.

    ``drawn`` is True when the figure was drawn anyway (crossing leaders or
    covered definition text only).
    """

    records: tuple[LabelFitRecord, ...]
    rows: tuple[LabelFitRow, ...]
    drawn: bool

    @property
    def overlap_pairs(self) -> int:
        return sum(record.fit.overlap_pairs for record in self.records)

    @property
    def crossing_pairs(self) -> int:
        return sum(record.fit.crossing_pairs for record in self.records)

    @property
    def definition_overlaps(self) -> int:
        return sum(record.fit.definition_overlaps for record in self.records)


class LabelFitRows(Protocol):
    """Computes one row of a report on request, without a time limit."""

    def row(self, kind: LabelFitRowKind) -> LabelFitRow: ...


# --- wording --------------------------------------------------------------------------------------

_CLI_OPTIONS: Mapping[str, str] = {
    X_OFFSET_PATH: "--outer_label_x_radius_offset",
    Y_OFFSET_PATH: "--outer_label_y_radius_offset",
    FONT_SIZE_PATHS[0]: "--label_font_size",
    SPACING_PATH: "--circular_label_spacing",
    SCOPE_PATH: "--labels",
    PLACEMENT_PATH: "--label_placement",
    SEPARATE_FEATURES_PATH: "--resolve_overlaps",
}
# The CLI spells the all-outside scope ``out``.
_CLI_WORDS: Mapping[str, str] = {ALL_OUTSIDE_SCOPE: "out"}
ROW_TITLES: Mapping[LabelFitRowKind, str] = {
    "outward_both": "outward, both directions",
    "outward_vertical": "outward, vertical only",
    "outward_horizontal": "outward, horizontal only",
    "smaller_font": "smaller font",
    "tighter_rows": "tighter rows",
    "combination": "combination",
    "labels_out": "all labels outside",
    "radial": "radial labels",
    "separate_features": "separate overlapping features",
}
_ROW_NAMES: Mapping[LabelFitRowKind, str] = {
    "outward_both": "outward in both directions",
    "outward_vertical": "outward vertical only",
    "outward_horizontal": "outward horizontal only",
    "smaller_font": "a smaller font",
    "tighter_rows": "tighter rows",
    "combination": "the combination",
    "labels_out": "all labels outside",
    "radial": "radial labels",
    "separate_features": "separating overlapping features",
}
_KIND_WORDS: Mapping[LabelKind, str] = {
    "circular_outer": "outer",
    "circular_inner": "inner",
    "circular_embedded": "embedded",
    "linear": "Linear",
}


def _counted(count: int, singular: str, plural: str) -> str:
    return f"{count} {singular if count == 1 else plural}"


def fit_counts(fit: LabelFit) -> str:
    counts = (
        f"{_counted(fit.overlap_pairs, 'label pair overlaps', 'label pairs overlap')}, "
        f"{_counted(fit.crossing_pairs, 'leader pair crosses', 'leader pairs cross')}"
    )
    if fit.definition_overlaps:
        counts += (
            f", {_counted(fit.definition_overlaps, 'label covers', 'labels cover')} the definition text"
        )
    return counts


def label_fit_record_counts(record: LabelFitRecord) -> str:
    """A record's counts as the report states them, or why radial labels were refused."""

    if record.refused:
        return "they need more room than the circular tracks leave"
    return fit_counts(record.fit)


def label_kind_words(record: LabelFitRecord) -> str:
    """The words naming a record's label kind (``outer``, ``radial inner``, ...)."""

    words = _KIND_WORDS[record.label_kind]
    return f"radial {words}" if record.refused else words


def _cli_value(value: LabelFitValue) -> str:
    return _CLI_WORDS.get(value, value) if isinstance(value, str) else f"{value:g}"


def _cli_argument(option: str, value: LabelFitValue) -> str:
    # A switch is a flag (``--resolve_overlaps``); rows only turn switches on.
    return option if isinstance(value, bool) else f"{option} {_cli_value(value)}"


def _cli_arguments(values: Mapping[str, LabelFitValue]) -> str:
    return " ".join(
        _cli_argument(option, values[path])
        for path, option in _CLI_OPTIONS.items()
        if path in values
    )


def _joined(names: Sequence[str]) -> str:
    if len(names) <= 1:
        return "".join(names)
    return ", ".join(names[:-1]) + " and " + names[-1]


def _first_line(report: LabelFitReport) -> str:
    outcome = "Drawn anyway." if report.drawn else "Nothing was drawn."
    records = report.records
    if len(records) == 1:
        (record,) = records
        return (
            f"{record.label_count} {label_kind_words(record)} labels do not fit at {record.font_size:g} pt "
            f"({label_fit_record_counts(record)}). {outcome}"
        )
    kinds = {label_kind_words(record) for record in records}
    if len(kinds) == 1 and len({record.record_id for record in records}) == len(records):
        listed = "; ".join(
            f"{record.record_id} ({record.label_count} labels at {record.font_size:g} pt: "
            f"{label_fit_record_counts(record)})"
            for record in records
        )
        words = next(iter(kinds))
        return f"{words[0].upper()}{words[1:]} labels do not fit in {len(records)} records: {listed}. {outcome}"
    listed = "; ".join(
        f"{record.record_id}: {record.label_count} {label_kind_words(record)} labels at "
        f"{record.font_size:g} pt ({label_fit_record_counts(record)})"
        for record in records
    )
    return f"Labels do not fit: {listed}. {outcome}"


def format_label_fit_report(report: LabelFitReport) -> str:
    """The CLI and Python text of the error and the warning (OD-F5, OD-F8)."""

    lines = [_first_line(report)]
    fitting = [(row, row.values) for row in report.rows if row.values is not None]
    if fitting:
        arguments = [_cli_arguments(values) for _, values in fitting]
        width = max(len(argument) for argument in arguments) + 3
        lines.append("Run again with one of:")
        lines.extend(
            f"  {argument:<{width}}({ROW_TITLES[row.kind]})"
            for argument, (row, _) in zip(arguments, fitting)
        )
    missing = [row.kind for row in report.rows if row.status == "not_found"]
    skipped = [row.kind for row in report.rows if row.status == "not_computed"]
    if missing and (fitting or skipped):
        names = _joined([_ROW_NAMES[kind] for kind in missing])
        verb = "do" if len(missing) > 1 or missing in (["tighter_rows"], ["radial"]) else "does"
        lines.append(f"{names[0].upper()}{names[1:]} {verb} not fit this map.")
    elif missing:
        lines.append("No single change of the label settings fits this map.")
    if skipped:
        seconds = WARN_SEARCH_SECONDS if report.drawn else SWITCH_SEARCH_SECONDS
        lines.append(
            f"Not searched within {seconds:g} s: {_joined([_ROW_NAMES[kind] for kind in skipped])}."
        )
    return "\n".join(lines)


def raise_or_warn(report: LabelFitReport, *, row_search: LabelFitRows | None = None) -> None:
    """Stop with ``LabelsDoNotFitError`` when nothing was drawn, else warn with the same words."""

    message = format_label_fit_report(report)
    if not report.drawn:
        raise LabelsDoNotFitError(message, report=report, row_search=row_search)
    warnings.warn(LabelLeadersCrossWarning(message, report=report), stacklevel=4)


__all__ = [
    "ALL_OUTSIDE_SCOPE",
    "DEFINITION_TEXT_MIN_OVERLAP_PX",
    "FITS",
    "FONT_SIZE_PATHS",
    "PLACEMENT_PATH",
    "RADIAL_PLACEMENT",
    "ROW_TITLES",
    "SCOPE_PATH",
    "SEPARATE_FEATURES_PATH",
    "SEPARATE_FEATURES_VALUES",
    "SWITCH_SEARCH_SECONDS",
    "SPACING_PATH",
    "WARN_SEARCH_SECONDS",
    "X_OFFSET_PATH",
    "Y_OFFSET_PATH",
    "LabelFit",
    "LabelFitOutcome",
    "LabelFitRecord",
    "LabelFitReport",
    "LabelFitRow",
    "LabelFitRowKind",
    "LabelFitRowStatus",
    "LabelFitVerification",
    "LabelFitRows",
    "LabelFitValue",
    "LabelKind",
    "LabelShape",
    "Polygon",
    "Polyline",
    "crossing_pairs",
    "fit_counts",
    "format_label_fit_report",
    "label_fit_outcome",
    "label_fit_record_counts",
    "label_fits_by_kind",
    "label_kind_words",
    "overlapping_pairs",
    "polyline_crossing_count",
    "raise_or_warn",
    "switch_row",
]
