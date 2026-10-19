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
"""

from __future__ import annotations

import math
import time
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal, TypeAlias

from .fit import (
    ALL_OUTSIDE_SCOPE,
    FITS,
    FONT_SIZE_PATHS,
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
    LabelFitOutcome,
    LabelFitValue,
    LabelFitVerification,
    LabelKind,
    label_fit_outcome,
    raise_or_warn,
    switch_row,
    verified_row,
)

OFFSET_STEP = 0.025
OFFSET_LIMIT = 2.0
FONT_STEP_PT = 0.5
FONT_LIMIT_PT = 6.0
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
)
# With inner, embedded or radial labels failing, the switches come first
# (SWITCHES-FIRST); separating features applies only to overlapping embedded labels.
_SWITCHES_FIRST_ORDER: tuple[LabelFitRowKind, ...] = (
    "separate_features",
    "labels_out",
    "radial",
    *_SEARCH_ORDER,
)
_SWITCH_VALUES: Mapping[LabelFitRowKind, Mapping[str, LabelFitValue]] = {
    "labels_out": {SCOPE_PATH: ALL_OUTSIDE_SCOPE},
    "radial": {PLACEMENT_PATH: RADIAL_PLACEMENT},
    "separate_features": SEPARATE_FEATURES_VALUES,
}
_OUTER_LEVER_REASON = "Does not move inner, embedded, or radial labels"

CircularLabelPlacement: TypeAlias = Literal["horizontal", "radial"]
CircularLabelFitProbe: TypeAlias = Callable[[Mapping[str, LabelFitValue]], LabelFit]


@dataclass(frozen=True)
class CircularLabelLayout:
    """One record's drawn labels: their fit by kind, and a probe for the outer rows.

    ``refused`` marks radial labels that could not be placed in the room the
    tracks leave; ``features_separated`` that overlapping features already
    have their own tracks (``--resolve_overlaps``). ``probe`` (horizontal records with outer labels) places the
    record's labels again with config overrides and returns the outer labels'
    fit; it may stop early on a layout that cannot fit, so its answer is a
    hint that the full assembly verifies.
    """

    record_id: str
    placement: CircularLabelPlacement
    font_size: float
    label_counts: Mapping[LabelKind, int]
    kind_fits: Mapping[LabelKind, LabelFit]
    x_radius_offset: float
    y_radius_offset: float
    label_spacing: float
    refused: bool = False
    features_separated: bool = False
    probe: CircularLabelFitProbe | None = field(default=None, compare=False, repr=False)
    # The largest font of the record's region labels, when some are cut or lie on another track (``circular_region``).
    region_font_size: float = 0.0

    @property
    def fit(self) -> LabelFit:
        return sum(self.kind_fits.values(), FITS)

    @property
    def fits(self) -> bool:
        return not self.refused and self.fit.fits

    @property
    def failing_kinds(self) -> tuple[LabelKind, ...]:
        """The sides that do not fit; a refused radial record names its inner or outer side."""

        failing = tuple(kind for kind, fit in self.kind_fits.items() if not fit.fits)
        if any(kind != "circular_region" for kind in failing) or not self.refused:
            return failing
        return (
            *failing,
            "circular_inner" if self.label_counts.get("circular_inner") else "circular_outer",
        )

    @property
    def failing_feature_kinds(self) -> tuple[LabelKind, ...]:
        """The failing sides the rows can move: region labels that are cut (LP-D87) or lie on other tracks (LP-D86)
        only warn, with no row."""

        return tuple(kind for kind in self.failing_kinds if kind != "circular_region")

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


def circular_label_outcome(layouts: Sequence[CircularLabelLayout]) -> LabelFitOutcome:
    """The decision point over a figure's records; refused radial labels stop."""

    return label_fit_outcome(
        [layout.fit for layout in layouts], refused=any(layout.refused for layout in layouts)
    )


class _SearchStopped(Exception):
    """The bounded search reached its deadline before a new probe."""


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


def _settings_key(values: Mapping[str, LabelFitValue]) -> tuple[tuple[str, LabelFitValue], ...]:
    return tuple(sorted(values.items()))


@dataclass(frozen=True)
class _Lever:
    grid: list[float]
    values_for: Callable[[float], dict[str, LabelFitValue]]
    reason: str
    bound: float | None = None
    limit_first: bool = True
    downward: bool = False


def _applicable_rows(layouts: Sequence[CircularLabelLayout]) -> frozenset[LabelFitRowKind]:
    """The rows of the failing sides: outer labels six, inner labels three, radial labels two.

    Overlapping embedded labels add the switch that separates overlapping
    features, unless they are separated already.
    """

    kinds: set[LabelFitRowKind] = set()
    for layout in layouts:
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
        assembly_seconds: float = 0.0,
    ) -> None:
        self.layouts = tuple(layouts)
        self._offer_notice_rows = offer_notice_rows
        self.row_kinds = _applicable_rows(self.layouts)
        self._outer_only = all(layout.outer_only for layout in self.layouts)
        self._verify_assembly = verify
        # A verification starts only if one more assembly fits before the deadline.
        self._assembly_seconds = float(assembly_seconds)
        # Multi-record rows start from the strongest current value among the failing records.
        self._x_offset = max(layout.x_radius_offset for layout in self.layouts)
        self._y_offset = max(layout.y_radius_offset for layout in self.layouts)
        self._font_size = min(layout.font_size for layout in self.layouts)
        self._spacing = min(layout.label_spacing for layout in self.layouts)
        self._probe_answers: list[dict[tuple[tuple[str, LabelFitValue], ...], LabelFit]] = [
            {} for _ in self.layouts
        ]
        self._verifications: dict[
            tuple[tuple[str, LabelFitValue], ...], LabelFitVerification
        ] = {}
        self._rows: dict[LabelFitRowKind, tuple[LabelFitRow, float | None]] = {}
        self._deadline: float | None = None

    def row(self, kind: LabelFitRowKind, *, deadline: float | None = None) -> LabelFitRow:
        """One row; with a ``time.monotonic()`` deadline it may be ``not_computed``."""

        self._deadline = deadline
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

        Rows that leave crossings or labels over the definition text come
        after the others (LP-D57).
        """

        found: dict[LabelFitRowKind, LabelFitRow] = {}
        order = _SEARCH_ORDER if self._outer_only else _SWITCHES_FIRST_ORDER
        for kind in order:
            if kind in self.row_kinds:
                found[kind] = self.row(kind, deadline=deadline)
        rows = [found[kind] for kind in ROW_KINDS if kind in found]
        return tuple(sorted(rows, key=lambda row: row.status == "found_with_notice"))

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

    def _verify(self, values: Mapping[str, LabelFitValue]) -> LabelFitVerification:
        key = _settings_key(values)
        if key not in self._verifications:
            self._check_deadline(self._assembly_seconds)
            started = time.monotonic()
            self._verifications[key] = self._verify_assembly(values)
            self._assembly_seconds = time.monotonic() - started
        return self._verifications[key]

    # --- rows ---

    def _row(self, kind: LabelFitRowKind) -> tuple[LabelFitRow, float | None]:
        if kind not in self._rows:
            self._rows[kind] = self._compute(kind)
        return self._rows[kind]

    def _compute(self, kind: LabelFitRowKind) -> tuple[LabelFitRow, float | None]:
        if kind in _SWITCH_VALUES:
            return self._switch(kind), None
        if self._outer_only:
            return self._search(kind, self._lever(kind))
        if kind == "smaller_font":
            return self._font_scan()
        return LabelFitRow(kind, "not_found", reason=_OUTER_LEVER_REASON), None

    def _switch(self, kind: LabelFitRowKind) -> LabelFitRow:
        values = _SWITCH_VALUES[kind]
        return switch_row(kind, values, self._verify(values), offer_notice=self._offer_notice_rows)

    def _font_scan(self) -> tuple[LabelFitRow, float | None]:
        """The largest smaller font that fits, trying each step down (fit is not monotone).

        Without one, the largest that leaves no overlap, when such rows are offered.
        """

        notice: LabelFitRow | None = None
        for value in _font_grid(self._font_size):
            values = _font_values(value)
            row = verified_row("smaller_font", values, self._verify(values), offer_notice=self._offer_notice_rows)
            if row is not None and row.status == "found":
                return row, value
            notice = notice or row
        if notice is not None:
            return notice, None
        return LabelFitRow(
            "smaller_font", "not_found", reason=f"Does not fit down to {FONT_LIMIT_PT:g} pt"
        ), None

    def _lever(self, kind: LabelFitRowKind) -> _Lever:
        y_offset = self._y_offset
        if kind == "outward_vertical":
            return _Lever(
                _offset_grid(y_offset),
                lambda value: {Y_OFFSET_PATH: value},
                f"Does not fit up to Y offset {OFFSET_LIMIT:g}",
                limit_first=False,
            )
        if kind == "outward_both":
            return _Lever(
                _offset_grid(max(self._x_offset, y_offset)),
                lambda value: {X_OFFSET_PATH: value, Y_OFFSET_PATH: value},
                f"Does not fit up to X and Y offset {OFFSET_LIMIT:g}",
                bound=self._row("outward_vertical")[1],
                downward=True,
            )
        if kind == "outward_horizontal":
            return _Lever(
                _offset_grid(self._x_offset),
                lambda value: {X_OFFSET_PATH: value},
                f"Does not fit up to X offset {OFFSET_LIMIT:g}",
            )
        if kind == "smaller_font":
            return _Lever(
                _font_grid(self._font_size),
                _font_values,
                f"Does not fit down to {FONT_LIMIT_PT:g} pt",
                limit_first=False,
            )
        if kind == "tighter_rows":
            return _Lever(
                _spacing_grid(self._spacing),
                lambda value: {SPACING_PATH: value},
                f"Does not fit down to {TIGHTER_ROW_GRID_PX[-1]:g} px label spacing",
            )
        # Combination (OD-F4, OD-F7): half the font row's drop, at least one step,
        # then the smallest Y offset at that font; the current Y offset comes first.
        font_value = self._row("smaller_font")[1]
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
            f"Does not fit at {combined_font:g} pt up to Y offset {OFFSET_LIMIT:g}",
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
    ) -> tuple[LabelFitRow, float] | None:
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
        for value in grid[strongest:]:
            values = lever.values_for(value)
            row = verified_row(kind, values, self._verify(values), offer_notice=offer_notice)
            if row is not None:
                return row, value
        return None


def _fit_records(layouts: Sequence[CircularLabelLayout]) -> tuple[LabelFitRecord, ...]:
    return tuple(
        LabelFitRecord(
            record_id=layout.record_id,
            label_kind=kind,
            label_count=int(layout.label_counts.get(kind, 0)),
            font_size=layout.region_font_size if kind == "circular_region" else layout.font_size,
            fit=layout.kind_fits.get(kind, FITS),
            refused=layout.refused and kind != "circular_region",
        )
        for layout in layouts
        for kind in layout.failing_kinds
    )


def check_circular_label_fit(
    layouts: Sequence[CircularLabelLayout],
    verify: CircularLabelFitVerifier,
    *,
    warn_search_seconds: float = WARN_SEARCH_SECONDS,
    assembly_seconds: float = 0.0,
) -> tuple[LabelFitReport, CircularLabelFitSearch] | None:
    """Raise, warn, or pass for a drawn figure's records.

    Returns ``None`` when every record fits. Raises ``LabelsDoNotFitError``
    when labels overlap or radial labels have no room; its rows are computed
    (outer-only rows all; the switches and the font scan of inner, embedded and
    radial labels for ``SWITCH_SEARCH_SECONDS``) and ``row_search`` computes a
    row not reached. When only leaders cross or labels cover the
    definition text, warns with ``LabelLeadersCrossWarning`` and returns the
    report (rows searched for ``warn_search_seconds``; 0 leaves every row
    ``not_computed``) and the search that computes the other rows.
    ``assembly_seconds`` is the time of one full assembly, so a verification
    that would overrun a deadline is not started.
    """

    outcome = circular_label_outcome(layouts)
    if outcome == "fits":
        return None
    failing = tuple(layout for layout in layouts if not layout.fits)
    drawn = outcome == "warn"
    search = CircularLabelFitSearch(
        failing, verify, offer_notice_rows=not drawn, assembly_seconds=assembly_seconds
    )
    search_seconds: float | None = (
        warn_search_seconds if drawn else None if search.outer_only else SWITCH_SEARCH_SECONDS
    )
    deadline = None if search_seconds is None else time.monotonic() + search_seconds
    report = LabelFitReport(
        records=_fit_records(failing), rows=search.rows(deadline=deadline), drawn=drawn
    )
    raise_or_warn(report, row_search=search, search_seconds=search_seconds)
    return report, search


__all__ = [
    "ROW_KINDS",
    "CircularLabelFitProbe",
    "CircularLabelFitSearch",
    "CircularLabelFitVerifier",
    "CircularLabelLayout",
    "CircularLabelPlacement",
    "check_circular_label_fit",
    "circular_label_outcome",
]
