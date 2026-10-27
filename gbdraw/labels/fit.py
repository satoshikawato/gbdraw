"""Whether drawn labels fit, the one decision point, and the report type.

Labels fit when no two labels overlap and no two leader polylines cross
(LP-D43 ``DEF-CROSS``); a label over the centre definition text counts like a
crossing (LP-D49 ``DEF-WARN``). A Linear region (annotation) label over other
labels or lines is counted by its placer and only warns (LP-D68), and so does a
Circular region label cut to fit its circle (LP-D87) or lying on another track
because set radii leave it no room (LP-D86). The centre definition text of a Circular record
whose lines meet what a track paints only warns too (D-12, ``circular_definition``). The counts take plain shapes: each label kind
has one adapter next to its owner that turns the drawn labels into
:class:`LabelShape` values, and one record's labels of every kind are counted
together (LP-D48 ``SCOPE-ALL-LABELS``).
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import Literal, Protocol, TypeAlias, get_args

from ..exceptions import LabelLeadersCrossWarning, LabelsDoNotFitError
from ..layout.spatial import (
    Aabb,
    GridCell,
    bounding_aabb,
    candidate_aabb_partners,
    cell_partners,
    grid_bucket_size,
    segment_cells,
)
from ..layout.text_geometry import Point, convex_polygons_intersect, segments_intersect

LabelKind: TypeAlias = Literal[
    "circular_outer",
    "circular_inner",
    "circular_embedded",
    "circular_region",
    "circular_definition",
    "linear",
    "linear_region",
]
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
    "smaller_definition_font",
    "narrow_tracks",
]
# The rows that answer for the Circular centre definition text over a track (D-12), not for the labels.
DEFINITION_ROW_KINDS: frozenset[LabelFitRowKind] = frozenset({"smaller_definition_font", "narrow_tracks"})
# ``found``: ``values`` fit; ``found_with_notice``: ``values`` leave no
# overlapping labels but leave crossing leaders or labels over the definition
# text, so the figure draws with the notice (LP-D57; offered only when nothing
# was drawn); ``not_found``: no value within the lever's limit (``reason`` says
# which); ``not_computed``: the bounded search stopped first.
LabelFitRowStatus: TypeAlias = Literal["found", "found_with_notice", "not_found", "not_computed"]
# A row sets numbers (offsets, font size, spacing), words (label scope,
# placement) and switches (separate overlapping features).
LabelFitValue: TypeAlias = bool | float | str
# What makes room for labels that need much more room than the map has (CAPACITY-STOP): label fewer features, shorten
# the label text, use a smaller font, or draw a Linear map.
LabelCapacityRemedy: TypeAlias = Literal["filter", "label_text", "font", "linear"]
Polygon: TypeAlias = tuple[Point, ...]
Polyline: TypeAlias = Sequence[Point]

# A pair of labels of two kinds counts once, for the first kind in this order.
_KIND_ORDER: tuple[LabelKind, ...] = (
    "circular_inner",
    "circular_outer",
    "circular_embedded",
    "circular_region",
    "circular_definition",
    "linear",
    "linear_region",
)
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
# ``--definition_font_size``: the Circular centre definition's font (its line interval follows it).
DEFINITION_FONT_PATH = "objects.definition.circular.font_size"
# ``--center_reserved_radius``: not a config path but the Circular track option (``CircularTrackOptions``) that keeps
# the inside tracks out of the centre; the "narrow the inside tracks" row writes it (PD-OI-097).
CENTER_RESERVED_RADIUS_PATH = "tracks.center_reserved_radius"


def split_row_option_values(values: Mapping[str, LabelFitValue]) -> tuple[dict[str, LabelFitValue], float | None]:
    """A row's settings as drawing again takes them: the config overrides, and the ``center_reserved_radius`` option
    when the row sets it (``CENTER_RESERVED_RADIUS_PATH``)."""

    config = {path: value for path, value in values.items() if path != CENTER_RESERVED_RADIUS_PATH}
    radius = values.get(CENTER_RESERVED_RADIUS_PATH)
    if radius is None:
        return config, None
    if isinstance(radius, (bool, str)):
        raise ValueError(f"{CENTER_RESERVED_RADIUS_PATH} takes a number, not {radius!r}")
    return config, float(radius)


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
    """Overlapping label pairs, crossing leader pairs, labels over the definition text,
    Linear region labels over other labels or lines (``region_overlaps``: pairs of a
    region label and what it covers, LP-D68), Circular region labels cut to fit
    their circle (``region_cuts``, LP-D87), Circular region labels whose letters
    lie on another track's band (``region_track_overlaps``, LP-D86), and lines of a
    Circular record's centre definition text that meet what a track paints, its
    ink (``definition_track_overlaps``, D-12).

    The labels fit when all seven are 0.
    """

    overlap_pairs: int
    crossing_pairs: int
    definition_overlaps: int = 0
    region_overlaps: int = 0
    region_cuts: int = 0
    region_track_overlaps: int = 0
    definition_track_overlaps: int = 0

    @property
    def fits(self) -> bool:
        return (
            self.overlap_pairs == 0
            and self.crossing_pairs == 0
            and self.definition_overlaps == 0
            and self.region_overlaps == 0
            and self.region_cuts == 0
            and self.region_track_overlaps == 0
            and self.definition_track_overlaps == 0
        )

    def __add__(self, other: LabelFit) -> LabelFit:
        return LabelFit(
            self.overlap_pairs + other.overlap_pairs,
            self.crossing_pairs + other.crossing_pairs,
            self.definition_overlaps + other.definition_overlaps,
            self.region_overlaps + other.region_overlaps,
            self.region_cuts + other.region_cuts,
            self.region_track_overlaps + other.region_track_overlaps,
            self.definition_track_overlaps + other.definition_track_overlaps,
        )


FITS = LabelFit(0, 0)


def _extent(points: Sequence[Point]) -> Aabb:
    """The box of a label's points (``bounding_aabb``); a label always has finite points."""

    box = bounding_aabb(points)
    if box is None:
        raise ValueError("a label shape needs finite points")
    return box


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

    first_box = _extent(first)
    second_box = _extent(second)
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
    boxes = [_extent([point for part in parts for point in part]) for parts in shapes]
    return [
        (left, right)
        for left, partners in enumerate(candidate_aabb_partners(boxes, bucket_size=grid_bucket_size(boxes)))
        for right in sorted(partners)
        if shapes_overlap(shapes[left], shapes[right], touching=touching)
    ]


def shapes_overlap(first: Sequence[Polygon], second: Sequence[Polygon], *, touching: bool) -> bool:
    """Whether two shapes (each a union of convex parts) overlap; ``touching`` as in :func:`overlapping_pairs`."""

    overlap = _closed_polygons_overlap if touching else _open_polygons_overlap
    return any(overlap(left, right) for left in first for right in second)


_CROSSING_BUCKET_PX = 64.0


def crossing_pairs(polylines: Sequence[Polyline]) -> list[tuple[int, int]]:
    """Index pairs of polylines that cross (touching does not count); empty polylines never cross."""

    indexed = [index for index, polyline in enumerate(polylines) if len(polyline) >= 2]

    def covers() -> Iterator[set[GridCell]]:
        # The buckets each polyline's segments pass through, not its box: a long radial
        # leader's box would make every leader in that part of the map a candidate (OV-1946).
        for index in indexed:
            polyline = polylines[index]
            yield set().union(
                *(segment_cells(start, end, _CROSSING_BUCKET_PX) for start, end in zip(polyline, polyline[1:]))
            )

    return [
        (indexed[left], indexed[right])
        for left, partners in enumerate(cell_partners(covers(), covers()))
        for right in sorted(partners)
        if polylines_cross(polylines[indexed[left]], polylines[indexed[right]])
    ]


def polylines_cross(first: Polyline, second: Polyline) -> bool:
    """Whether two polylines cross (touching does not count)."""

    # Loops, not a generator: the label repair tests about 1.3 M leader pairs on a record (N4 phage).
    for left in zip(first, first[1:]):
        for right in zip(second, second[1:]):
            if segments_intersect(left, right, touching=False):
                return True
    return False


def polyline_crossing_count(polylines: Sequence[Polyline]) -> int:
    """Count polyline pairs that cross (touching does not count)."""

    return len(crossing_pairs(polylines))


def covers_definition_text(shape: LabelShape, definition_text: Sequence[Aabb]) -> bool:
    """Whether a label with a leader covers a line of the definition text (embedded labels never do)."""

    if not shape.leader or not definition_text:
        return False
    label_box = _extent([point for part in shape.parts for point in part])
    return any(
        min(label_box.max_x, text_box.max_x) - max(label_box.min_x, text_box.min_x)
        > DEFINITION_TEXT_MIN_OVERLAP_PX
        and min(label_box.max_y, text_box.max_y) - max(label_box.min_y, text_box.min_y)
        > DEFINITION_TEXT_MIN_OVERLAP_PX
        for text_box in definition_text
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
    for shape in shapes:
        if covers_definition_text(shape, definition_text):
            counts[shape.kind][2] += 1
    return {
        kind: LabelFit(overlaps, crossings, covered)
        for kind, (overlaps, crossings, covered) in counts.items()
    }


def label_fit_outcome(fits: Sequence[LabelFit], *, refused: bool = False) -> LabelFitOutcome:
    """The one decision point: overlaps stop; crossings, covered definition text and region labels over others alone warn.

    ``refused`` is a placement that could not place its labels at all (radial
    labels without room, LP-D49 ``RADIAL-JOIN``); it stops like an overlap.
    """

    if refused or any(fit.overlap_pairs > 0 for fit in fits):
        return "stop"
    if any(not fit.fits for fit in fits):
        return "warn"
    return "fits"


@dataclass(frozen=True)
class LabelCapacityRecord:
    """The room one side of a record's labels needs, as a multiple of the room around the map (``ratio``).

    ``remedies`` are the ways to make room that apply, in the order they are offered.
    """

    label_count: int
    font_size: float
    ratio: float
    remedies: tuple[LabelCapacityRemedy, ...]


@dataclass(frozen=True)
class LabelFitRecord:
    """One record's labels of one kind that do not fit.

    ``refused`` marks radial labels that could not be placed in the room the
    tracks leave (their counts may then be 0). ``capacity`` marks labels that
    need much more room than the map has, so they were not placed (CAPACITY-STOP;
    their counts are 0). ``record_index`` is the index of the record's drawn
    groups (``data-gbdraw-record-index``), the record's identity where two
    records share an ID; Circular records carry it, Linear records do not yet
    (None).
    """

    record_id: str
    label_kind: LabelKind
    label_count: int
    font_size: float
    fit: LabelFit
    refused: bool = False
    record_index: int | None = None
    capacity: LabelCapacityRecord | None = None


@dataclass(frozen=True)
class LabelFitRow:
    """One lever: its mildest fitting value, or why there is none.

    ``values`` maps config paths to the values to set and ``canvas_size`` is
    the (width, height) of the verified figure; both are set only when
    ``status`` is ``found`` or ``found_with_notice``. ``remaining`` is what a
    ``found_with_notice`` row leaves: crossing leaders and labels over the
    definition text.
    """

    kind: LabelFitRowKind
    status: LabelFitRowStatus
    values: Mapping[str, LabelFitValue] | None = None
    canvas_size: tuple[float, float] | None = None
    reason: str | None = None
    remaining: LabelFit | None = None


@dataclass(frozen=True)
class LabelFitVerification:
    """The full assembly with a row's settings.

    ``outcome`` is :func:`label_fit_outcome` over every record (an assembly
    that raised stops); ``canvas_size`` is the figure's size, ``fit`` its
    labels' fit over every record, and ``reason`` why the assembly failed
    when it raised.
    """

    outcome: LabelFitOutcome
    canvas_size: tuple[float, float] | None
    fit: LabelFit | None = None
    reason: str | None = None

    @property
    def fits(self) -> bool:
        return self.outcome == "fits"


def verified_row(
    kind: LabelFitRowKind,
    values: Mapping[str, LabelFitValue],
    verification: LabelFitVerification,
    *,
    offer_notice: bool,
) -> LabelFitRow | None:
    """The row a verified value gives: ``found`` when it fits; with ``offer_notice``,
    ``found_with_notice`` when it only warns (LP-D57); otherwise None."""

    if verification.outcome == "fits":
        return LabelFitRow(kind, "found", values, verification.canvas_size)
    if offer_notice and verification.outcome == "warn":
        return LabelFitRow(
            kind, "found_with_notice", values, verification.canvas_size, remaining=verification.fit
        )
    return None


def switch_row(
    kind: LabelFitRowKind,
    values: Mapping[str, LabelFitValue],
    verification: LabelFitVerification,
    *,
    offer_notice: bool,
) -> LabelFitRow:
    """A row with one setting to turn on, from its full assembly."""

    row = verified_row(kind, values, verification, offer_notice=offer_notice)
    if row is not None:
        return row
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

    @property
    def region_overlaps(self) -> int:
        return sum(record.fit.region_overlaps for record in self.records)

    @property
    def region_cuts(self) -> int:
        return sum(record.fit.region_cuts for record in self.records)

    @property
    def region_track_overlaps(self) -> int:
        return sum(record.fit.region_track_overlaps for record in self.records)

    @property
    def definition_track_overlaps(self) -> int:
        return sum(record.fit.definition_track_overlaps for record in self.records)

    @property
    def capacity_remedies(self) -> tuple[LabelCapacityRemedy, ...]:
        """The ways to make room for the records that were not placed (CAPACITY-STOP), in the order offered."""

        offered = {remedy for record in self.records if record.capacity for remedy in record.capacity.remedies}
        return tuple(remedy for remedy in get_args(LabelCapacityRemedy) if remedy in offered)


class LabelFitRows(Protocol):
    """Computes one row of a report on request; a ``time.monotonic()`` deadline may stop it (``not_computed``)."""

    def row(self, kind: LabelFitRowKind, *, deadline: float | None = None) -> LabelFitRow: ...


# --- wording --------------------------------------------------------------------------------------

_CLI_OPTIONS: Mapping[str, str] = {
    X_OFFSET_PATH: "--outer_label_x_radius_offset",
    Y_OFFSET_PATH: "--outer_label_y_radius_offset",
    FONT_SIZE_PATHS[0]: "--label_font_size",
    SPACING_PATH: "--circular_label_spacing",
    SCOPE_PATH: "--labels",
    PLACEMENT_PATH: "--label_placement",
    SEPARATE_FEATURES_PATH: "--resolve_overlaps",
    DEFINITION_FONT_PATH: "--definition_font_size",
    CENTER_RESERVED_RADIUS_PATH: "--center_reserved_radius",
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
    "smaller_definition_font": "smaller definition font",
    "narrow_tracks": "narrow the inside tracks",
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
    "smaller_definition_font": "a smaller definition font",
    "narrow_tracks": "narrowing the inside tracks",
}
_KIND_WORDS: Mapping[LabelKind, str] = {
    "circular_outer": "outer",
    "circular_inner": "inner",
    "circular_embedded": "embedded",
    "circular_region": "region",
    "circular_definition": "definition",
    "linear": "Linear",
    "linear_region": "region",
}


def _counted(count: int, singular: str, plural: str) -> str:
    return f"{count} {singular if count == 1 else plural}"


def fit_counts(fit: LabelFit) -> str:
    counts: list[str] = []
    # Region labels (LP-D68, LP-D87, LP-D86) and the definition text over a track (D-12) are counted on their own
    # record, which has no other counts.
    if fit.overlap_pairs or fit.crossing_pairs or not (
        fit.region_overlaps or fit.region_cuts or fit.region_track_overlaps or fit.definition_track_overlaps
    ):
        counts.append(_counted(fit.overlap_pairs, "label pair overlaps", "label pairs overlap"))
        counts.append(_counted(fit.crossing_pairs, "leader pair crosses", "leader pairs cross"))
    if fit.definition_overlaps:
        counts.append(f"{_counted(fit.definition_overlaps, 'label covers', 'labels cover')} the definition text")
    if fit.region_overlaps:
        counts.append(_counted(fit.region_overlaps, "region label overlap", "region label overlaps"))
    if fit.region_cuts:
        counts.append(_counted(fit.region_cuts, "region label cut", "region labels cut"))
    if fit.region_track_overlaps:
        counts.append(_counted(fit.region_track_overlaps, "region label on another track", "region labels on other tracks"))
    if fit.definition_track_overlaps:
        counts.append(_counted(fit.definition_track_overlaps, "line overlaps a track", "lines overlap a track"))
    return ", ".join(counts)


def row_note(row: LabelFitRow) -> str | None:
    """What a ``found_with_notice`` row leaves of what it answers for, in the counts' words: crossing leaders and
    labels over the definition text for a label row (LP-D57), definition lines over a track for a definition row."""

    remaining = row.remaining
    if row.status != "found_with_notice" or remaining is None:
        return None
    if row.kind in DEFINITION_ROW_KINDS:
        if not remaining.definition_track_overlaps:
            return None
        return (
            _counted(remaining.definition_track_overlaps, "definition line still overlaps", "definition lines still overlap")
            + " a track"
        )
    notes: list[str] = []
    if remaining.crossing_pairs:
        notes.append(_counted(remaining.crossing_pairs, "leader pair still crosses", "leader pairs still cross"))
    if remaining.definition_overlaps:
        notes.append(
            _counted(remaining.definition_overlaps, "label still covers", "labels still cover")
            + " the definition text"
        )
    return " and ".join(notes) or None


def label_fit_notice(report: LabelFitReport) -> tuple[str, ...]:
    """The notice bar's sentences over a drawn figure (the Owner's Q8 wording)."""

    sentences: list[str] = []
    if report.crossing_pairs:
        sentences.append(f"Leader lines cross ({_counted(report.crossing_pairs, 'pair', 'pairs')}).")
    if report.definition_overlaps:
        sentences.append(f"Labels cover the definition text ({report.definition_overlaps}).")
    if report.region_overlaps:
        sentences.append(f"Region labels overlap other labels or lines ({report.region_overlaps}).")
    if report.region_cuts:
        sentences.append(f"Region labels are cut to fit their circle ({report.region_cuts}).")
    if report.region_track_overlaps:
        sentences.append(
            f"Region labels lie on other tracks: the set track radii leave them no room ({report.region_track_overlaps})."
        )
    if report.definition_track_overlaps:
        sentences.append(
            f"The definition text overlaps a track ({_counted(report.definition_track_overlaps, 'line', 'lines')})."
        )
    return tuple(sentences)


_RADIAL_REFUSAL = "need more room than the circular tracks leave"
# The CLI's ways to make room (CAPACITY-STOP); the Web names its own controls for the same remedies.
_CLI_CAPACITY_REMEDIES: Mapping[LabelCapacityRemedy, str] = {
    "filter": "label fewer features with --label_blacklist or --label_whitelist;",
    "label_text": "shorten the label text with --qualifier_priority (for example gene names only);",
    "font": "use a smaller --label_font_size;",
    "linear": "or draw a Linear map (gbdraw linear).",
}


def _room_needed(capacity: LabelCapacityRecord) -> str:
    return f"need about {capacity.ratio:.1f} times the room around the map"


def label_capacity_sentence(record_id: str, label_kind: LabelKind, capacity: LabelCapacityRecord, *, placed: bool) -> str:
    """The sentence naming a record's labels that need much more room than the map has: ``placed`` False when they
    were not placed (CAPACITY-STOP), True when they are placed after this warning (CAPACITY-WARN)."""

    words = _KIND_WORDS[label_kind]
    labels = _counted(capacity.label_count, f"{words} label", f"{words} labels")
    outcome = "may not fit. Placing them can take minutes." if placed else "were not placed."
    return f"{record_id}: {labels} at {capacity.font_size:g} pt {_room_needed(capacity)}, so they {outcome}"
# What else shortens a centre definition that overlaps a track: its text (an abbreviation typed as the species).
DEFINITION_TEXT_REMEDY = "Or shorten the definition text with --species or --strain (for example, an abbreviation)."


def label_fit_error_counts(report: LabelFitReport) -> str:
    """The counts sentence of the Web error over a figure that was not drawn.

    A refused record's counts are those of the layout it refused, so a report
    with one states the refusal instead (OV-1194).
    """

    if any(record.refused for record in report.records):
        return f"Radial labels {_RADIAL_REFUSAL}."
    largest: dict[LabelKind, LabelCapacityRecord] = {}
    for record in report.records:
        if record.capacity is not None and (
            record.label_kind not in largest or record.capacity.ratio > largest[record.label_kind].ratio
        ):
            largest[record.label_kind] = record.capacity
    if largest:
        return " ".join(
            f"{_KIND_WORDS[kind].capitalize()} labels {_room_needed(capacity)}." for kind, capacity in largest.items()
        )
    return f"Overlapping label pairs: {report.overlap_pairs}. Crossing leader pairs: {report.crossing_pairs}."


def label_fit_record_counts(record: LabelFitRecord) -> str:
    """A record's counts as the report states them, or why radial labels were refused."""

    if record.refused:
        return f"they {_RADIAL_REFUSAL}"
    if record.capacity is not None:
        return f"they {_room_needed(record.capacity)}"
    return fit_counts(record.fit)


def label_kind_words(record: LabelFitRecord) -> str:
    """The words naming a record's label kind (``outer``, ``radial inner``, ...)."""

    words = _KIND_WORDS[record.label_kind]
    return f"radial {words}" if record.refused else words


def label_kind_nouns(record: LabelFitRecord) -> tuple[str, str]:
    """The singular and plural noun of a record's counted items: ``outer label(s)``, ``definition line(s)``."""

    words = label_kind_words(record)
    noun = "line" if record.label_kind == "circular_definition" else "label"
    return f"{words} {noun}", f"{words} {noun}s"


def label_kind_title(record: LabelFitRecord) -> str:
    """A record's counted items as the Web dialog lists them (``outer labels``, ``definition line``)."""

    singular, plural = label_kind_nouns(record)
    return singular if record.label_count == 1 else plural


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
        # The record is named: a run that writes each record to its own file reports each on its own (OV-1836).
        (record,) = records
        singular, plural = label_kind_nouns(record)
        return (
            f"{record.record_id}: {_counted(record.label_count, f'{singular} does', f'{plural} do')} not fit at "
            f"{record.font_size:g} pt ({label_fit_record_counts(record)}). {outcome}"
        )
    kinds = {label_kind_nouns(record) for record in records}
    if len(kinds) == 1 and len({record.record_id for record in records}) == len(records):
        singular, plural = next(iter(kinds))
        noun = singular.rsplit(" ", 1)[-1]
        listed = "; ".join(
            f"{record.record_id} ({_counted(record.label_count, noun, f'{noun}s')} at {record.font_size:g} pt: "
            f"{label_fit_record_counts(record)})"
            for record in records
        )
        return f"{plural[0].upper()}{plural[1:]} do not fit in {len(records)} records: {listed}. {outcome}"
    listed = "; ".join(
        f"{record.record_id}: {_counted(record.label_count, *label_kind_nouns(record))} at "
        f"{record.font_size:g} pt ({label_fit_record_counts(record)})"
        for record in records
    )
    return f"Labels do not fit: {listed}. {outcome}"


def format_label_fit_report(report: LabelFitReport, *, search_seconds: float | None = None) -> str:
    """The CLI and Python text of the error and the warning (OD-F5, OD-F8).

    ``search_seconds`` is the bound the row search used (OV-708); by default
    the warning's or the error's default bound.
    """

    if any(record.capacity for record in report.records):
        return _capacity_report(report)
    lines = [_first_line(report)]
    fitting = [(row, row.values) for row in report.rows if row.values is not None]
    if fitting:
        arguments = [_cli_arguments(values) for _, values in fitting]
        width = max(len(argument) for argument in arguments) + 3
        lines.append("Run again with one of:")
        lines.extend(
            f"  {argument:<{width}}({'; '.join(filter(None, (ROW_TITLES[row.kind], row_note(row))))})"
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
        seconds = search_seconds
        if seconds is None:
            seconds = WARN_SEARCH_SECONDS if report.drawn else SWITCH_SEARCH_SECONDS
        lines.append(
            f"Not searched within {seconds:g} s: {_joined([_ROW_NAMES[kind] for kind in skipped])}."
        )
    if any(record.label_kind == "circular_definition" for record in report.records):
        # The Web dialog's "Edit Species / Strain…" button (PD-OI-097).
        lines.append(DEFINITION_TEXT_REMEDY)
    return "\n".join(lines)


def _capacity_report(report: LabelFitReport) -> str:
    """The error over records whose labels were not placed (CAPACITY-STOP): one sentence per such record, the other
    failing records as the error states them, and the ways to make room."""

    lines = [
        label_capacity_sentence(record.record_id, record.label_kind, record.capacity, placed=False)
        for record in report.records
        if record.capacity is not None
    ]
    others = tuple(record for record in report.records if record.capacity is None)
    if others:
        lines.append(_first_line(LabelFitReport(others, (), drawn=False)))
    else:
        lines[-1] += " Nothing was drawn."
    lines.append("To make room:")
    lines.extend(f"  - {_CLI_CAPACITY_REMEDIES[remedy]}" for remedy in report.capacity_remedies)
    return "\n".join(lines)


def raise_or_warn(
    report: LabelFitReport, *, row_search: LabelFitRows | None = None, search_seconds: float | None = None
) -> None:
    """Stop with ``LabelsDoNotFitError`` when nothing was drawn, else warn with the same words."""

    message = format_label_fit_report(report, search_seconds=search_seconds)
    if not report.drawn:
        raise LabelsDoNotFitError(message, report=report, row_search=row_search)
    warnings.warn(LabelLeadersCrossWarning(message, report=report), stacklevel=4)


__all__ = [
    "ALL_OUTSIDE_SCOPE",
    "CENTER_RESERVED_RADIUS_PATH",
    "DEFINITION_FONT_PATH",
    "DEFINITION_ROW_KINDS",
    "DEFINITION_TEXT_REMEDY",
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
    "LabelCapacityRecord",
    "LabelCapacityRemedy",
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
    "covers_definition_text",
    "crossing_pairs",
    "fit_counts",
    "format_label_fit_report",
    "label_fit_error_counts",
    "label_capacity_sentence",
    "label_fit_notice",
    "label_fit_outcome",
    "label_fit_record_counts",
    "label_fits_by_kind",
    "label_kind_nouns",
    "label_kind_title",
    "label_kind_words",
    "overlapping_pairs",
    "polyline_crossing_count",
    "polylines_cross",
    "raise_or_warn",
    "row_note",
    "shapes_overlap",
    "split_row_option_values",
    "switch_row",
    "verified_row",
]
