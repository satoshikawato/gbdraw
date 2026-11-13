"""Circular label fit check and the search for settings that fit.

The labels' fit and the one decision point live in :mod:`gbdraw.labels.fit`:
overlapping labels stop the run (LP-D42), crossing leaders or labels over the
centre definition text alone draw as before and warn (LP-D46 ``CROSS-WARN``,
LP-D49 ``DEF-WARN``), and radial labels without room stop like overlaps
(LP-D49 ``RADIAL-JOIN``). Every outcome carries one report: the failing
records by side and the rows of settings that apply to them, each the
mildest value of one lever that makes every record fit, verified with the
full assembly.

Outer labels have six rows, searched with fit probes (fit is monotone along
each lever). Inner labels and radial labels have the font row, scanned down
with the full assembly because fit is not monotone along the font there,
and two switches, all labels outside and radial labels (LP-D49 ``INNER-3``),
which are verified first (``SWITCHES-FIRST``). Embedded labels that overlap
also get the switch that puts overlapping features on separate tracks
(LP-D52 ``E-ROW``). When the figure is drawn anyway the search stops after
``WARN_SEARCH_SECONDS`` (LP-D48 ``SEARCH-BOUNDED-5S``); when nothing is
drawn, the switches and the font scan of inner and radial labels share
``SWITCH_SEARCH_SECONDS`` (OD-F15), and outer-only rows are all computed. A
row not reached is computed later with :meth:`CircularLabelFitSearch.row`.
When nothing is drawn and a lever has no value that fits, its mildest value
that leaves no overlapping labels is offered with what it leaves (LP-D57).

A record's centre definition text whose painted lines meet a track's ink is
drawn and warns too (D-12). A drawn figure with it gets two more rows: a
smaller definition font, scanned down with the full assembly like the inner
labels' font, and narrowing the inside tracks, offered when the definition's
band enters a track's band and no ``center_reserved_radius`` is set. That row
sets ``center_reserved_radius`` to the band the definition takes on one line
(PD-OI-097). The tracks narrow for the definition on their own, each Auto
width down to ``canvas.circular.track_narrowing_min_px`` (PD-OI-078, rule L);
the narrow row is offered when that did not clear the definition, and its
radius narrows them below that floor. "Over a track" is read two ways: the notice and the rows'
verification read the lines against the ink, the narrow row reads the band
against the tracks' bands. The band is a square extent, not the line boxes'
corners, so the narrowed figure may still have letters on a track and then
draws with the notice. Each row is checked for what it
fixes: the label rows for the labels, the definition rows for the definition
text (neither may stop the figure).
"""

from __future__ import annotations

import math
import time
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from typing import Literal, TypeAlias

from .fit import (
    ALL_OUTSIDE_SCOPE,
    CENTER_RESERVED_RADIUS_PATH,
    DEFINITION_FONT_PATH,
    DEFINITION_ROW_KINDS,
    FITS,
    OFFERED_ROW_STATUSES,
    REGION_WARNING_FIELDS,
    FONT_LIMIT_PT,
    FONT_SIZE_PATHS,
    CombinationLimitReason,
    LabelCapacityRecord,
    PLACEMENT_PATH,
    RADIAL_PLACEMENT,
    SCOPE_PATH,
    SEPARATE_FEATURES_VALUES,
    SPACING_PATH,
    SWITCH_SEARCH_SECONDS,
    WARN_SEARCH_SECONDS,
    X_OFFSET_PATH,
    Y_OFFSET_PATH,
    LabelFit,
    LabelFitRecord,
    LabelFitReport,
    LabelFitRow,
    LabelFitRowKind,
    LabelFitRowReason,
    LabelFitOutcome,
    LabelFitValue,
    LabelFitVerification,
    LabelKind,
    LimitReason,
    PlainReason,
    label_fit_outcome,
    raise_or_warn,
    shown_capacity_ratio,
    switch_row,
    verified_row,
)
from ..layout.figure_size import LargerFigure

OFFSET_STEP = 0.025
OFFSET_LIMIT = 2.0
FONT_STEP_PT = 0.5
TIGHTER_ROW_GRID_PX = (2.5, 2.0, 1.5, 1.0, 0.5, 0.25, 0.1)

# The dialog order of the rows; a report has the rows that apply to its records.
ROW_KINDS: tuple[LabelFitRowKind, ...] = (
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
)
_OUTER_ROWS: frozenset[LabelFitRowKind] = frozenset(ROW_KINDS[:6])
_INNER_ROWS: frozenset[LabelFitRowKind] = frozenset({"smaller_font", "labels_out", "radial"})
_RADIAL_ROWS: frozenset[LabelFitRowKind] = frozenset({"smaller_font", "labels_out"})
# The search order of outer rows: vertical bounds both directions and the
# combination, and the combination starts from the font row's value.
_SEARCH_ORDER: tuple[LabelFitRowKind, ...] = (
    "outward_vertical",
    "outward_both",
    "outward_horizontal",
    "smaller_font",
    "tighter_rows",
    "combination",
    "smaller_definition_font",
    "narrow_tracks",
)
# With inner, embedded or radial labels failing, the switches come first
# (SWITCHES-FIRST); separating features applies only to overlapping embedded labels.
_SWITCHES_FIRST_ORDER: tuple[LabelFitRowKind, ...] = (
    "separate_features",
    "labels_out",
    "radial",
    *_SEARCH_ORDER,
)
# The kinds no label row moves: region labels only warn (LP-D86, LP-D87); the definition text has its own rows.
_UNMOVED_LABEL_KINDS: frozenset[LabelKind] = frozenset({"circular_region", "circular_definition"})
_SWITCH_VALUES: Mapping[LabelFitRowKind, Mapping[str, LabelFitValue]] = {
    "labels_out": {SCOPE_PATH: ALL_OUTSIDE_SCOPE},
    "radial": {PLACEMENT_PATH: RADIAL_PLACEMENT},
    "separate_features": SEPARATE_FEATURES_VALUES,
}
_OUTER_LEVER_REASON = PlainReason("outer_labels_only")
_FONT_LIMIT_REASON = LimitReason("font_limit", FONT_LIMIT_PT)

CircularLabelPlacement: TypeAlias = Literal["horizontal", "radial"]
CircularLabelFitProbe: TypeAlias = Callable[[Mapping[str, LabelFitValue]], LabelFit]
# Whether a record's definition text clears every track at a definition font size (D-12), from the tracks laid out
# before the labels are placed: a hint that the full assembly verifies.
CircularDefinitionFitProbe: TypeAlias = Callable[[float], bool]


@dataclass(frozen=True)
class CircularLabelLayout:
    """One record's drawn labels: their fit by kind, and a probe for the outer rows.

    ``refused`` marks radial labels that could not be placed in the room the
    tracks leave; ``capacity`` the sides whose horizontal labels need much more
    room than the map has, so the record's labels were not placed
    (CAPACITY-STOP); ``features_separated`` that overlapping features already
    have their own tracks (``--resolve_overlaps``). ``probe`` (horizontal records with outer labels) places the
    record's labels again with config overrides and returns the outer labels'
    fit; it may stop early on a layout that cannot fit, so its answer is a
    hint that the full assembly verifies.
    """

    record_id: str
    # The index of the record's drawn groups (``data-gbdraw-record-index``): the record's identity on the figure.
    record_index: int
    placement: CircularLabelPlacement
    font_size: float
    label_counts: Mapping[LabelKind, int]
    kind_fits: Mapping[LabelKind, LabelFit]
    x_radius_offset: float
    y_radius_offset: float
    label_spacing: float
    refused: bool = False
    capacity: Mapping[LabelKind, LabelCapacityRecord] = field(default_factory=dict)
    features_separated: bool = False
    probe: CircularLabelFitProbe | None = field(default=None, compare=False, repr=False)
    # The largest font of the record's region labels, when some are cut or lie on another track (``circular_region``).
    region_font_size: float = 0.0
    # The centre definition's font, when its painted lines meet a track's ink (``circular_definition``).
    definition_font_size: float = 0.0
    definition_probe: CircularDefinitionFitProbe | None = field(default=None, compare=False, repr=False)
    # The ``center_reserved_radius`` that keeps the inside tracks out of the definition's one-line band, when that
    # band enters a track's band (not the ink the notice reads) and no radius is set (``circular_definition``); 0
    # otherwise.
    definition_narrow_radius_px: float = 0.0

    @property
    def fit(self) -> LabelFit:
        return sum(self.kind_fits.values(), FITS)

    @property
    def fits(self) -> bool:
        return not self.refused and not self.capacity and self.fit.fits

    @property
    def failing_kinds(self) -> tuple[LabelKind, ...]:
        """The sides that do not fit; a refused radial record names its inner or outer side."""

        failing = tuple(kind for kind, fit in self.kind_fits.items() if not fit.fits)
        failing += tuple(kind for kind in self.capacity if kind not in failing)
        if any(kind not in _UNMOVED_LABEL_KINDS for kind in failing) or not self.refused:
            return failing
        return (
            *failing,
            "circular_inner" if self.label_counts.get("circular_inner") else "circular_outer",
        )

    @property
    def failing_feature_kinds(self) -> tuple[LabelKind, ...]:
        """The failing sides the label rows can move: region labels that are cut (LP-D87) or lie on other tracks
        (LP-D86) only warn, with no row, and the definition text over a track has its own row (D-12)."""

        return tuple(kind for kind in self.failing_kinds if kind not in _UNMOVED_LABEL_KINDS)

    @property
    def definition_only(self) -> bool:
        """Only the definition text fails: the label rows do not search this record."""

        return self.failing_kinds == ("circular_definition",)

    @property
    def outer_only(self) -> bool:
        """Only horizontal outer labels fail: the six outer rows search with probes."""

        return (
            self.placement == "horizontal"
            and not self.refused
            and self.failing_feature_kinds == ("circular_outer",)
            and self.probe is not None
        )


CircularLabelFitVerifier: TypeAlias = Callable[[Mapping[str, LabelFitValue]], LabelFitVerification]


def labels_not_placed(layouts: Sequence[CircularLabelLayout]) -> bool:
    """Whether any record's labels need much more room than the map has, so they were not placed (CAPACITY-STOP)."""

    return any(layout.capacity for layout in layouts)


def circular_label_outcome(layouts: Sequence[CircularLabelLayout]) -> LabelFitOutcome:
    """The decision point over a figure's records; refused radial labels and labels over capacity stop."""

    return label_fit_outcome(
        [layout.fit for layout in layouts],
        refused=labels_not_placed(layouts) or any(layout.refused for layout in layouts),
    )


class _SearchStopped(Exception):
    """The bounded search reached its deadline before a new probe."""


class _RowDeferred(_SearchStopped):
    """In the first pass, a row whose value depends on a row still walking waits for the second pass (D-118)."""


def _offset_grid(current: float) -> list[float]:
    grid: list[float] = []
    step = 1
    while True:
        value = round(current + OFFSET_STEP * step, 3)
        if value >= OFFSET_LIMIT - 1e-9:
            break
        grid.append(value)
        step += 1
    if current < OFFSET_LIMIT - 1e-9:
        grid.append(OFFSET_LIMIT)
    return grid


def _font_grid(current: float) -> list[float]:
    grid: list[float] = []
    step = 1
    while True:
        value = round(current - FONT_STEP_PT * step, 3)
        if value <= FONT_LIMIT_PT + 1e-9:
            break
        grid.append(value)
        step += 1
    if current > FONT_LIMIT_PT + 1e-9:
        grid.append(FONT_LIMIT_PT)
    return grid


def _spacing_grid(current: float) -> list[float]:
    return [value for value in TIGHTER_ROW_GRID_PX if value < current - 1e-9]


def _font_values(font_size: float) -> dict[str, LabelFitValue]:
    return {path: font_size for path in FONT_SIZE_PATHS}


def _round_half_pt(value: float) -> float:
    return math.floor(value / FONT_STEP_PT + 0.5) * FONT_STEP_PT


def _bisect(fits: Callable[[float], bool], grid: Sequence[float], low: int, high: int) -> int:
    """Smallest fitting index in [low, high]; ``grid[high]`` fits."""

    while low < high:
        middle = (low + high) // 2
        if fits(grid[middle]):
            high = middle
        else:
            low = middle + 1
    return high


def _mildest_fitting_index(
    grid: Sequence[float],
    fits: Callable[[float], bool],
    bound: int | None = None,
    *,
    limit_first: bool = True,
    downward: bool = False,
) -> int | None:
    """Smallest index whose value fits, assuming fit is monotone along the grid.

    ``bound`` is an index another row suggests; when it fits, the answer lies at
    or below it. ``downward`` then gallops down from it (for a row whose value
    is usually the bound's), else an exponential search runs up from the
    mildest value. Without a fitting bound, ``limit_first`` probes the limit
    first, so a row that rarely has a value costs one probe; otherwise the
    limit is probed only when nothing milder fits. A bisection inside the
    bracket finds the answer. Probes that do not fit stop early on overlaps, so
    the search spends them near the mild end.
    """

    if not grid:
        return None
    last = len(grid) - 1
    high: int | None = None
    if bound is not None and bound < last and fits(grid[bound]):
        last = high = bound
        if downward:
            step = 1
            while high - step >= 0:
                if not fits(grid[high - step]):
                    return _bisect(fits, grid, high - step + 1, high)
                high -= step
                step *= 2
            return _bisect(fits, grid, 0, high)
    elif limit_first:
        if not fits(grid[last]):
            return None
        high = last
    low = 0
    span = 1
    while True:
        probe_index = low + span - 1
        if probe_index >= last:
            break
        if fits(grid[probe_index]):
            high = probe_index
            break
        low = probe_index + 1
        span *= 2
    if high is None:
        if not fits(grid[last]):
            return None
        high = last
    return _bisect(fits, grid, low, high)


def _radius_rounded_up(radius_px: float) -> float:
    """A radius the "narrow the inside tracks" row writes: up to the next 0.01 px, so the band stays clear."""

    return math.ceil(round(float(radius_px) * 100.0, 6)) / 100.0


def _settings_key(values: Mapping[str, LabelFitValue]) -> tuple[tuple[str, LabelFitValue], ...]:
    return tuple(sorted(values.items()))


@dataclass(frozen=True)
class _Lever:
    grid: list[float]
    values_for: Callable[[float], dict[str, LabelFitValue]]
    reason: LabelFitRowReason
    bound: float | None = None
    limit_first: bool = True
    downward: bool = False


def _applicable_rows(
    layouts: Sequence[CircularLabelLayout], *, drawn: bool
) -> frozenset[LabelFitRowKind]:
    """The rows of the failing sides: outer labels six, inner labels three, radial labels two.

    Overlapping embedded labels add the switch that separates overlapping
    features, unless they are separated already. A drawn figure whose
    definition lines meet a track's ink adds the smaller definition font, and
    narrowing the inside tracks when the definition's band also enters a
    track's band and no ``center_reserved_radius`` is set (when nothing was
    drawn, the labels come first).
    """

    kinds: set[LabelFitRowKind] = set()
    for layout in layouts:
        if drawn and "circular_definition" in layout.failing_kinds:
            kinds.add("smaller_definition_font")
            if layout.definition_narrow_radius_px > 0.0:
                kinds.add("narrow_tracks")
        embedded = layout.kind_fits.get("circular_embedded", FITS)
        if embedded.overlap_pairs and not layout.features_separated:
            kinds.add("separate_features")
        failing = layout.failing_feature_kinds
        if not failing:
            continue
        if layout.placement == "radial":
            kinds |= _RADIAL_ROWS
            continue
        if "circular_outer" in failing:
            kinds |= _OUTER_ROWS
        if any(kind != "circular_outer" for kind in failing):
            kinds |= _INNER_ROWS
    if not any(layout.label_counts.get("circular_inner") for layout in layouts):
        kinds.discard("labels_out")
    return frozenset(kinds)


class CircularLabelFitSearch:
    """The rows for one figure's failing records.

    When only horizontal outer labels fail, each row is the mildest grid value
    that every failing record's probe accepts (the strongest of their
    answers), verified with the full assembly and stepped one grid value
    stronger until verification fits. With ``offer_notice_rows`` (nothing
    was drawn), a lever without a value that fits searches again for the
    mildest value that leaves no overlap (LP-D57). Otherwise the switches are one full
    assembly each, the font is scanned down with the full assembly (the first
    value that fits), and the outer levers have no value. Probe and
    verification answers are kept, so rows share them and a row asked for
    again, or after its deadline, costs only what is new.
    """

    def __init__(
        self,
        layouts: Sequence[CircularLabelLayout],
        verify: CircularLabelFitVerifier,
        *,
        offer_notice_rows: bool,
    ) -> None:
        self.layouts = tuple(layouts)
        self._offer_notice_rows = offer_notice_rows
        self.row_kinds = _applicable_rows(self.layouts, drawn=not offer_notice_rows)
        # The label rows read only the records whose labels fail; a definition-only record has its own row.
        label_layouts = [layout for layout in self.layouts if not layout.definition_only]
        self._outer_only = all(layout.outer_only for layout in label_layouts)
        label_layouts = label_layouts or list(self.layouts)
        self._verify_assembly = verify
        # The region label warnings each record already has: a label row leaves out only these (D-112, D-115).
        self._record_fits = {layout.record_index: layout.fit for layout in self.layouts}
        # A verification starts only if one more, as long as the last, ends before the deadline; the first one
        # starts whenever the deadline has not passed (OV-2527: a full first assembly is no estimate of it).
        self._verification_seconds = 0.0
        # Multi-record rows start from the strongest current value among the failing records.
        self._x_offset = max(layout.x_radius_offset for layout in label_layouts)
        self._y_offset = max(layout.y_radius_offset for layout in label_layouts)
        self._font_size = min(layout.font_size for layout in label_layouts)
        self._definition_font_size = min(
            (layout.definition_font_size for layout in self.layouts if "circular_definition" in layout.failing_kinds),
            default=0.0,
        )
        self._spacing = min(layout.label_spacing for layout in label_layouts)
        self._probe_answers: list[dict[tuple[tuple[str, LabelFitValue], ...], LabelFit]] = [
            {} for _ in self.layouts
        ]
        self._verifications: dict[
            tuple[tuple[str, LabelFitValue], ...], LabelFitVerification
        ] = {}
        self._rows: dict[LabelFitRowKind, tuple[LabelFitRow, float | None]] = {}
        # A row found before the deadline that a later ask may make milder, whose walk past a region label warning
        # (OV-3110) is not finished, or that starts from such a row, is not kept in ``_rows``.
        self._unfinished = False
        # ``rows`` takes every row's first value before any row walks on past it (D-110).
        self._first_values_only = False
        self._deadline: float | None = None

    def row(self, kind: LabelFitRowKind, *, deadline: float | None = None) -> LabelFitRow:
        """One row; with a ``time.monotonic()`` deadline it may be ``not_computed``."""

        self._deadline = deadline
        self._unfinished = False
        try:
            return self._row(kind)[0]
        except _SearchStopped:
            return LabelFitRow(kind, "not_computed")
        finally:
            self._deadline = None

    @property
    def outer_only(self) -> bool:
        """Only horizontal outer labels fail: the rows search with probes."""

        return self._outer_only

    def rows(self, *, deadline: float | None = None) -> tuple[LabelFitRow, ...]:
        """The rows that apply, in the dialog order, searched until ``deadline``.

        Rows that leave crossings, labels over the definition text, or a
        region label warning come after the others (LP-D57, OV-3105). Every
        row's first value comes first; then the rows whose figure keeps a
        region label warning walk on to a value whose whole figure fits while
        the deadline allows (OV-3110, D-110). Only the combination, whose font
        comes from the smaller-font row, waits for that row's walk; a row the
        deadline leaves out is not searched (D-118).
        """

        order = _SEARCH_ORDER if self._outer_only else _SWITCHES_FIRST_ORDER
        kinds = [kind for kind in order if kind in self.row_kinds]
        self._first_values_only = True
        try:
            found = {kind: self.row(kind, deadline=deadline) for kind in kinds}
        finally:
            self._first_values_only = False
        for kind in kinds:
            if kind not in self._rows:
                row = self.row(kind, deadline=deadline)
                found[kind] = found[kind] if row.status == "not_computed" else row
        rows = [found[kind] for kind in ROW_KINDS if kind in found]
        return tuple(sorted(rows, key=lambda row: row.status in OFFERED_ROW_STATUSES[1:]))

    # --- probes ---

    def _check_deadline(self, expected_seconds: float = 0.0) -> None:
        if self._deadline is not None and time.monotonic() + expected_seconds >= self._deadline:
            raise _SearchStopped

    def _probe_fit(self, index: int, values: Mapping[str, LabelFitValue]) -> LabelFit:
        answers = self._probe_answers[index]
        key = _settings_key(values)
        if key not in answers:
            self._check_deadline()
            probe = self.layouts[index].probe
            if probe is None:  # pragma: no cover - outer rows search only records with a probe
                raise ValueError("an outer row searched a record without a probe")
            answers[key] = probe(values)
        return answers[key]

    def _verify(self, kind: LabelFitRowKind, values: Mapping[str, LabelFitValue]) -> LabelFitVerification:
        key = _settings_key(values)
        if key not in self._verifications:
            self._check_deadline(self._verification_seconds)
            started = time.monotonic()
            self._verifications[key] = self._verify_assembly(values)
            self._verification_seconds = time.monotonic() - started
        return _row_verification(kind, self._verifications[key], self._record_fits)

    # --- rows ---

    def _row(self, kind: LabelFitRowKind) -> tuple[LabelFitRow, float | None]:
        if kind in self._rows:
            return self._rows[kind]
        asking_row_unfinished = self._unfinished
        self._unfinished = False
        row, value = self._compute(kind)
        computed = self._with_region_notice(row), value
        finished = not self._unfinished
        if finished:
            self._rows[kind] = computed
        # A row that starts from an unfinished row is unfinished too: a later ask may bound it differently.
        self._unfinished = asking_row_unfinished or not finished
        return computed

    def _region_warned(self, row: LabelFitRow) -> LabelFit | None:
        """The figure's fit when a label row fits the labels it moves but its figure keeps region label warnings the
        figure already had, which no label row moves (OV-3105, D-112); otherwise None."""

        if row.status != "found" or row.kind in DEFINITION_ROW_KINDS or row.values is None:
            return None
        fit = self._verifications[_settings_key(row.values)].fit
        return fit if fit is not None and fit.region_warnings else None

    def _with_region_notice(self, row: LabelFitRow) -> LabelFitRow:
        """A row whose figure keeps a region label warning draws with the notice; it bounds no row that starts from
        it (D-118)."""

        fit = self._region_warned(row)
        return row if fit is None else replace(row, status="found_with_notice", remaining=fit)

    def _clear_of_regions(
        self,
        kind: LabelFitRowKind,
        found: tuple[LabelFitRow, float],
        later_values: Sequence[float],
        values_for: Callable[[float], dict[str, LabelFitValue]],
    ) -> tuple[LabelFitRow, float | None]:
        """A row whose first value keeps a region label warning walks on, as when it answered for the whole figure,
        to the first later value of its lever whose whole figure fits; without one up to the lever's limit, its first
        value stays and draws with the notice (OV-3110, D-117). Only a value whose whole figure fits bounds the rows
        that start from this one, as on dev; a row offered with the notice bounds none (D-118). At the deadline the
        first value is the row, kept out of the cache so a later ask walks on."""

        if self._region_warned(found[0]) is None:
            return found
        if self._first_values_only:
            self._unfinished = True
            return found[0], None
        try:
            for value in later_values:
                values = values_for(value)
                row = verified_row(kind, values, self._verify(kind, values), offer_notice=False)
                if row is not None and self._region_warned(row) is None:
                    return row, value
        except _SearchStopped:
            self._unfinished = True
        return found[0], None

    def _compute(self, kind: LabelFitRowKind) -> tuple[LabelFitRow, float | None]:
        if kind in _SWITCH_VALUES:
            return self._switch(kind), None
        if kind == "smaller_definition_font":
            return self._definition_font_scan(), None
        if kind == "narrow_tracks":
            return self._narrow_tracks(), None
        if self._outer_only:
            return self._search(kind, self._lever(kind))
        if kind == "smaller_font":
            return self._font_scan()
        return LabelFitRow(kind, "not_found", reason=_OUTER_LEVER_REASON), None

    def _switch(self, kind: LabelFitRowKind) -> LabelFitRow:
        values = _SWITCH_VALUES[kind]
        return switch_row(kind, values, self._verify(kind, values), offer_notice=self._offer_notice_rows)

    def _font_scan(self) -> tuple[LabelFitRow, float | None]:
        """The largest smaller font that fits, trying each step down (fit is not monotone).

        Without one, the largest that leaves no overlap, when such rows are offered.
        """

        notice: LabelFitRow | None = None
        grid = _font_grid(self._font_size)
        for index, value in enumerate(grid):
            values = _font_values(value)
            row = verified_row(
                "smaller_font", values, self._verify("smaller_font", values), offer_notice=self._offer_notice_rows
            )
            if row is not None and row.status == "found":
                return self._clear_of_regions("smaller_font", (row, value), grid[index + 1 :], _font_values)
            notice = notice or row
        if notice is not None:
            return notice, None
        return LabelFitRow("smaller_font", "not_found", reason=_FONT_LIMIT_REASON), None

    def _definition_font_scan(self) -> LabelFitRow:
        """The largest smaller definition font at which the definition text clears every track (D-12).

        The reserve and the tracks inside follow the font, so fit need not be monotone: the records' probes
        (the tracks before the labels are placed) take the first step down at which every record clears, and the
        full assembly verifies it, stepping down until it fits. Only a drawn figure offers this row.
        """

        grid = _font_grid(self._definition_font_size)
        probes = [
            layout.definition_probe for layout in self.layouts if "circular_definition" in layout.failing_kinds
        ]
        start = 0
        if all(probe is not None for probe in probes):

            def clears(value: float) -> bool:
                self._check_deadline()
                return all(probe(value) for probe in probes if probe is not None)

            start = next((index for index, value in enumerate(grid) if clears(value)), len(grid))

        def verified(index: int) -> LabelFitRow | None:
            values: dict[str, LabelFitValue] = {DEFINITION_FONT_PATH: grid[index]}
            return verified_row(
                "smaller_definition_font", values, self._verify("smaller_definition_font", values), offer_notice=False
            )

        for index in range(start, len(grid)):
            row = verified(index)
            if row is None:
                continue
            # The labels may widen the circle after the probe, so a milder step may fit too. At the deadline the
            # step already verified is the row, kept out of the cache so a later ask walks on.
            try:
                while index > 0 and (milder := verified(index - 1)) is not None:
                    index, row = index - 1, milder
            except _SearchStopped:
                self._unfinished = True
            return row
        return LabelFitRow("smaller_definition_font", "not_found", reason=_FONT_LIMIT_REASON)

    def _narrow_tracks(self) -> LabelFitRow:
        """The inside tracks kept out of the definition's one-line band (``center_reserved_radius``, PD-OI-097).

        One radius serves every record (a Multi-Record Canvas shares it), the largest band. The letters may still
        meet a track, as they did when the tracks narrowed on their own: the row then draws with the notice.
        """

        radius = max(
            (layout.definition_narrow_radius_px for layout in self.layouts if "circular_definition" in layout.failing_kinds),
            default=0.0,
        )
        values: dict[str, LabelFitValue] = {CENTER_RESERVED_RADIUS_PATH: _radius_rounded_up(radius)}
        return switch_row("narrow_tracks", values, self._verify("narrow_tracks", values), offer_notice=True)

    def _lever(self, kind: LabelFitRowKind) -> _Lever:
        y_offset = self._y_offset
        if kind == "outward_vertical":
            return _Lever(
                _offset_grid(y_offset),
                lambda value: {Y_OFFSET_PATH: value},
                LimitReason("y_offset_limit", OFFSET_LIMIT),
                limit_first=False,
            )
        if kind == "outward_both":
            return _Lever(
                _offset_grid(max(self._x_offset, y_offset)),
                lambda value: {X_OFFSET_PATH: value, Y_OFFSET_PATH: value},
                LimitReason("xy_offset_limit", OFFSET_LIMIT),
                bound=self._row("outward_vertical")[1],
                downward=True,
            )
        if kind == "outward_horizontal":
            return _Lever(
                _offset_grid(self._x_offset),
                lambda value: {X_OFFSET_PATH: value},
                LimitReason("x_offset_limit", OFFSET_LIMIT),
            )
        if kind == "smaller_font":
            return _Lever(
                _font_grid(self._font_size),
                _font_values,
                _FONT_LIMIT_REASON,
                limit_first=False,
            )
        if kind == "tighter_rows":
            return _Lever(
                _spacing_grid(self._spacing),
                lambda value: {SPACING_PATH: value},
                LimitReason("spacing_limit", TIGHTER_ROW_GRID_PX[-1]),
            )
        # Combination (OD-F4, OD-F7): half the font row's drop, at least one step,
        # then the smallest Y offset at that font; the current Y offset comes first.
        font_value = self._row("smaller_font")[1]
        if self._first_values_only and "smaller_font" not in self._rows:
            # The font row is still walking: its value decides this row's font, so this row waits (D-118).
            raise _RowDeferred
        drop = max(FONT_STEP_PT, _round_half_pt((self._font_size - (font_value or FONT_LIMIT_PT)) / 2.0))
        combined_font = max(FONT_LIMIT_PT, self._font_size - drop)

        def combined_values(value: float) -> dict[str, LabelFitValue]:
            values = _font_values(combined_font)
            if value != y_offset:
                values[Y_OFFSET_PATH] = value
            return values

        return _Lever(
            [y_offset, *_offset_grid(y_offset)],
            combined_values,
            CombinationLimitReason(combined_font, OFFSET_LIMIT),
            bound=self._row("outward_vertical")[1],
        )

    def _search(self, kind: LabelFitRowKind, lever: _Lever) -> tuple[LabelFitRow, float | None]:
        """The row and its grid value (None without a value that fits).

        The mildest value that fits; without one, when such rows are offered,
        the mildest value that leaves no overlap (LP-D57).
        """

        for offer_notice in (False, True) if self._offer_notice_rows else (False,):
            found = self._mildest_row(kind, lever, offer_notice=offer_notice)
            if found is not None:
                row, value = found
                # Only a value that fits bounds the rows that start from this one.
                return row, value if row.status == "found" else None
        return LabelFitRow(kind, "not_found", reason=lever.reason), None

    def _mildest_row(
        self, kind: LabelFitRowKind, lever: _Lever, *, offer_notice: bool
    ) -> tuple[LabelFitRow, float | None] | None:
        """The mildest value every record's probe accepts, stepped stronger until verified.

        Probes accept a fit, or with ``offer_notice`` a layout without overlaps.
        """

        def accepts(fit: LabelFit) -> bool:
            return label_fit_outcome([fit]) != "stop" if offer_notice else fit.fits

        grid = lever.grid
        bound_index = (
            next((index for index, value in enumerate(grid) if value >= lever.bound - 1e-9), None)
            if lever.bound is not None
            else None
        )
        strongest = 0
        for index in range(len(self.layouts)):

            def record_fits(value: float, index: int = index) -> bool:
                return accepts(self._probe_fit(index, lever.values_for(value)))

            mildest = _mildest_fitting_index(
                grid,
                record_fits,
                bound_index,
                limit_first=lever.limit_first,
                downward=lever.downward,
            )
            if mildest is None:
                return None
            strongest = max(strongest, mildest)
        for index in range(strongest, len(grid)):
            values = lever.values_for(grid[index])
            row = verified_row(kind, values, self._verify(kind, values), offer_notice=offer_notice)
            if row is not None:
                return self._clear_of_regions(kind, (row, grid[index]), grid[index + 1 :], lever.values_for)
        return None


def _kept_region_warnings(record_fits: Mapping[int, LabelFit], had: Mapping[int, LabelFit]) -> LabelFit:
    """The region label warnings each record keeps of those it had: per record and field, at most its count in
    ``had`` (D-112, D-115). A warning that moves to another record, or one more of a kind, is not kept."""

    kept = FITS
    for index, fit in record_fits.items():
        before = had.get(index, FITS)
        kept += LabelFit(0, 0, **{name: min(fit.count(name), before.count(name)) for name in REGION_WARNING_FIELDS})
    return kept


def _row_verification(
    kind: LabelFitRowKind, verification: LabelFitVerification, had: Mapping[int, LabelFit]
) -> LabelFitVerification:
    """A full assembly as one row reads it: the label rows answer for the labels they move, not for the region
    label warnings each record already had (``had``, by record index, up to its count of each, D-112, D-115;
    ``_UNMOVED_LABEL_KINDS``; ``_with_region_notice`` says which warning stays, OV-3105), so a value that adds one
    does not fit; the definition rows answer for the definition text over a track (D-12); an assembly that stops is
    read as it is."""

    fit = verification.fit
    if verification.outcome == "stop" or fit is None:
        return verification
    if kind in DEFINITION_ROW_KINDS:
        outcome: LabelFitOutcome = "fits" if fit.definition_track_overlaps == 0 else "warn"
    else:
        outcome = label_fit_outcome([replace(
            fit.without_region_warnings(_kept_region_warnings(verification.record_fits, had)),
            definition_track_overlaps=0,
        )])
    return replace(verification, outcome=outcome)


def _kind_font_size(layout: CircularLabelLayout, kind: LabelKind) -> float:
    if kind == "circular_region":
        return layout.region_font_size
    if kind == "circular_definition":
        return layout.definition_font_size
    return layout.font_size


def _fit_records(layouts: Sequence[CircularLabelLayout]) -> tuple[LabelFitRecord, ...]:
    return tuple(
        LabelFitRecord(
            record_id=layout.record_id,
            label_kind=kind,
            label_count=int(layout.label_counts.get(kind, 0)),
            font_size=_kind_font_size(layout, kind),
            fit=layout.kind_fits.get(kind, FITS),
            refused=layout.refused and kind not in _UNMOVED_LABEL_KINDS,
            record_index=layout.record_index,
            capacity=layout.capacity.get(kind),
        )
        for layout in layouts
        for kind in layout.failing_kinds
    )


def check_circular_label_fit(
    layouts: Sequence[CircularLabelLayout],
    verify: CircularLabelFitVerifier,
    *,
    warn_search_seconds: float = WARN_SEARCH_SECONDS,
    larger_figure: Callable[[float | None], LargerFigure] | None = None,
) -> tuple[LabelFitReport, CircularLabelFitSearch] | None:
    """Raise, warn, or pass for a drawn figure's records.

    Returns ``None`` when every record fits. Raises ``LabelsDoNotFitError``
    when labels overlap or radial labels have no room; its rows are computed
    (outer-only rows all; the switches and the font scan of inner, embedded and
    radial labels for ``SWITCH_SEARCH_SECONDS``) and ``row_search`` computes a
    row not reached. When only leaders cross, labels cover the
    definition text, or the definition text reaches into a track (D-12), warns with ``LabelLeadersCrossWarning`` and returns the
    report (rows searched for ``warn_search_seconds``; 0 leaves every row
    ``not_computed``) and the search that computes the other rows.
    A verification that would overrun a deadline, if it takes as long as the
    last one, is not started; the first one starts while any time is left. A record whose labels need
    much more room than the map has (CAPACITY-STOP) raises at once: the report
    lists every failing record, has no rows, and no row is searched. With a
    requested figure size, ``larger_figure`` gives the size those records need
    from the largest room they need (the ratio the message states); any other
    report names the requested sides without an estimate (None), offered when
    no row fits (NO-ROW-WAYS).
    """

    outcome = circular_label_outcome(layouts)
    if outcome == "fits":
        return None
    failing = tuple(layout for layout in layouts if not layout.fits)
    if labels_not_placed(failing):
        records = _fit_records(failing)
        ratio = max(shown_capacity_ratio(record.capacity) for record in records if record.capacity is not None)
        raise_or_warn(
            LabelFitReport(
                records=records,
                rows=(),
                drawn=False,
                larger_figure=None if larger_figure is None else larger_figure(ratio),
            )
        )
    drawn = outcome == "warn"
    search = CircularLabelFitSearch(failing, verify, offer_notice_rows=not drawn)
    search_seconds: float | None = (
        warn_search_seconds if drawn else None if search.outer_only else SWITCH_SEARCH_SECONDS
    )
    deadline = None if search_seconds is None else time.monotonic() + search_seconds
    report = LabelFitReport(
        records=_fit_records(failing),
        rows=search.rows(deadline=deadline),
        drawn=drawn,
        larger_figure=None if larger_figure is None else larger_figure(None),
    )
    raise_or_warn(report, row_search=search, search_seconds=search_seconds)
    return report, search


__all__ = [
    "ROW_KINDS",
    "CircularDefinitionFitProbe",
    "CircularLabelFitProbe",
    "CircularLabelFitSearch",
    "CircularLabelFitVerifier",
    "CircularLabelLayout",
    "CircularLabelPlacement",
    "check_circular_label_fit",
    "circular_label_outcome",
    "labels_not_placed",
]
