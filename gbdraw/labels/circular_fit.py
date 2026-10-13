"""Circular outer-label fit check and the search for settings that fit.

The labels' fit and the one decision point live in :mod:`gbdraw.labels.fit`:
overlapping labels stop the run (LP-D42), crossing leaders alone draw as
before and warn (LP-D46 ``CROSS-WARN``). Both outcomes carry one report: the
failing records and six rows of settings, each the mildest value of one lever
that fits, found with fit probes and verified with the full assembly. When the
figure is drawn anyway the search stops after ``WARN_SEARCH_SECONDS``
(LP-D48 ``SEARCH-BOUNDED-5S``); a row it did not reach can be computed later
with :meth:`CircularLabelFitSearch.row`.
"""

from __future__ import annotations

import math
import time
import warnings
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import TypeAlias

from ..exceptions import LabelLeadersCrossWarning, LabelsDoNotFitError
from .fit import (
    LabelFit,
    LabelFitRecord,
    LabelFitReport,
    LabelFitRow,
    LabelFitRowKind,
    label_fit_outcome,
)

X_OFFSET_PATH = "labels.unified_adjustment.outer_labels.x_radius_offset"
Y_OFFSET_PATH = "labels.unified_adjustment.outer_labels.y_radius_offset"
# ``--label_font_size`` sets both size classes; a row writes them the same way.
FONT_SIZE_PATHS = ("labels.font_size.short", "labels.font_size.long")
SPACING_PATH = "labels.spacing.circular"

OFFSET_STEP = 0.025
OFFSET_LIMIT = 2.0
FONT_STEP_PT = 0.5
FONT_LIMIT_PT = 6.0
TIGHTER_ROW_GRID_PX = (2.5, 2.0, 1.5, 1.0, 0.5, 0.25, 0.1)
WARN_SEARCH_SECONDS = 5.0

# The dialog order of the rows.
ROW_KINDS: tuple[LabelFitRowKind, ...] = (
    "outward_both",
    "outward_vertical",
    "outward_horizontal",
    "smaller_font",
    "tighter_rows",
    "combination",
)
# The search order: vertical bounds both directions and the combination, and
# the combination starts from the font row's value.
_SEARCH_ORDER: tuple[LabelFitRowKind, ...] = (
    "outward_vertical",
    "outward_both",
    "outward_horizontal",
    "smaller_font",
    "tighter_rows",
    "combination",
)

CircularLabelFitProbe: TypeAlias = Callable[[Mapping[str, float]], LabelFit]


@dataclass(frozen=True)
class CircularOuterLabelLayout:
    """One record's drawn horizontal outer labels and a probe for changed settings.

    ``probe`` places the record's labels again with config overrides (dotted
    config paths) and returns their fit. It may stop early on a layout that
    cannot fit; its answer is a hint that the full assembly verifies.
    """

    record_id: str
    outer_label_count: int
    font_size: float
    x_radius_offset: float
    y_radius_offset: float
    label_spacing: float
    fit: LabelFit
    probe: CircularLabelFitProbe = field(compare=False, repr=False)


@dataclass(frozen=True)
class CircularLabelVerification:
    """The full assembly with a row's settings: whether every record fits, and its canvas size."""

    fits: bool
    canvas_size: tuple[float, float]


CircularLabelFitVerifier: TypeAlias = Callable[[Mapping[str, float]], CircularLabelVerification]


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


def _font_values(font_size: float) -> dict[str, float]:
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


def _settings_key(values: Mapping[str, float]) -> tuple[tuple[str, float], ...]:
    return tuple(sorted(values.items()))


@dataclass(frozen=True)
class _Lever:
    grid: list[float]
    values_for: Callable[[float], dict[str, float]]
    reason: str
    bound: float | None = None
    limit_first: bool = True
    downward: bool = False


class CircularLabelFitSearch:
    """The rows for one figure's failing records.

    Each row is the mildest grid value that every failing record's probe
    accepts (the strongest of their answers), verified with the full assembly
    and stepped one grid value stronger until verification fits. Probe and
    verification answers are kept, so rows share them and a row asked for
    again costs nothing.
    """

    def __init__(
        self,
        layouts: Sequence[CircularOuterLabelLayout],
        verify: CircularLabelFitVerifier,
    ) -> None:
        self.layouts = tuple(layouts)
        self._verify_assembly = verify
        # Multi-record rows start from the strongest current value among the failing records.
        self._x_offset = max(layout.x_radius_offset for layout in self.layouts)
        self._y_offset = max(layout.y_radius_offset for layout in self.layouts)
        self._font_size = min(layout.font_size for layout in self.layouts)
        self._spacing = min(layout.label_spacing for layout in self.layouts)
        self._probe_answers: list[dict[tuple[tuple[str, float], ...], bool]] = [
            {} for _ in self.layouts
        ]
        self._verifications: dict[tuple[tuple[str, float], ...], CircularLabelVerification] = {}
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

    def rows(self, *, deadline: float | None = None) -> tuple[LabelFitRow, ...]:
        """All six rows in the dialog order."""

        found = {kind: self.row(kind, deadline=deadline) for kind in _SEARCH_ORDER}
        return tuple(found[kind] for kind in ROW_KINDS)

    # --- probes ---

    def _check_deadline(self) -> None:
        if self._deadline is not None and time.monotonic() >= self._deadline:
            raise _SearchStopped

    def _probe_fits(self, index: int, values: Mapping[str, float]) -> bool:
        answers = self._probe_answers[index]
        key = _settings_key(values)
        if key not in answers:
            self._check_deadline()
            answers[key] = self.layouts[index].probe(values).fits
        return answers[key]

    def _verify(self, values: Mapping[str, float]) -> CircularLabelVerification:
        key = _settings_key(values)
        if key not in self._verifications:
            self._check_deadline()
            self._verifications[key] = self._verify_assembly(values)
        return self._verifications[key]

    # --- rows ---

    def _row(self, kind: LabelFitRowKind) -> tuple[LabelFitRow, float | None]:
        if kind not in self._rows:
            self._rows[kind] = self._search(kind, self._lever(kind))
        return self._rows[kind]

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

        def combined_values(value: float) -> dict[str, float]:
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
        """The row and its grid value (None without a value)."""

        not_found = (LabelFitRow(kind, "not_found", reason=lever.reason), None)
        grid = lever.grid
        bound_index = (
            next((index for index, value in enumerate(grid) if value >= lever.bound - 1e-9), None)
            if lever.bound is not None
            else None
        )
        strongest = 0
        for index in range(len(self.layouts)):

            def record_fits(value: float, index: int = index) -> bool:
                return self._probe_fits(index, lever.values_for(value))

            mildest = _mildest_fitting_index(
                grid,
                record_fits,
                bound_index,
                limit_first=lever.limit_first,
                downward=lever.downward,
            )
            if mildest is None:
                return not_found
            strongest = max(strongest, mildest)
        for verified_index in range(strongest, len(grid)):
            values = lever.values_for(grid[verified_index])
            verification = self._verify(values)
            if verification.fits:
                row = LabelFitRow(kind, "found", values, verification.canvas_size)
                return row, grid[verified_index]
        return not_found


def check_circular_label_fit(
    layouts: Sequence[CircularOuterLabelLayout],
    verify: CircularLabelFitVerifier,
) -> tuple[LabelFitReport, CircularLabelFitSearch] | None:
    """Raise, warn, or pass for a drawn figure's records.

    Returns ``None`` when every record fits. Raises ``LabelsDoNotFitError``
    with every row when labels overlap. When only leaders cross, warns with
    ``LabelLeadersCrossWarning`` and returns the report (rows searched for
    ``WARN_SEARCH_SECONDS``) and the search that computes the other rows.
    """

    outcome = label_fit_outcome([layout.fit for layout in layouts])
    if outcome == "fits":
        return None
    failing = tuple(layout for layout in layouts if not layout.fit.fits)
    search = CircularLabelFitSearch(failing, verify)
    drawn = outcome == "warn"
    report = LabelFitReport(
        records=tuple(
            LabelFitRecord(
                record_id=layout.record_id,
                label_kind="circular_outer",
                label_count=layout.outer_label_count,
                font_size=layout.font_size,
                fit=layout.fit,
            )
            for layout in failing
        ),
        rows=search.rows(deadline=time.monotonic() + WARN_SEARCH_SECONDS if drawn else None),
        drawn=drawn,
    )
    message = format_circular_label_fit_report(report)
    if not drawn:
        raise LabelsDoNotFitError(message, report=report)
    warnings.warn(LabelLeadersCrossWarning(message, report=report), stacklevel=3)
    return report, search


# --- wording --------------------------------------------------------------------------------------

_CLI_OPTIONS: Mapping[str, str] = {
    X_OFFSET_PATH: "--outer_label_x_radius_offset",
    Y_OFFSET_PATH: "--outer_label_y_radius_offset",
    FONT_SIZE_PATHS[0]: "--label_font_size",
    SPACING_PATH: "--circular_label_spacing",
}
_ROW_TITLES: Mapping[LabelFitRowKind, str] = {
    "outward_both": "outward, both directions",
    "outward_vertical": "outward, vertical only",
    "outward_horizontal": "outward, horizontal only",
    "smaller_font": "smaller font",
    "tighter_rows": "tighter rows",
    "combination": "combination",
}
_ROW_NAMES: Mapping[LabelFitRowKind, str] = {
    "outward_both": "outward in both directions",
    "outward_vertical": "outward vertical only",
    "outward_horizontal": "outward horizontal only",
    "smaller_font": "a smaller font",
    "tighter_rows": "tighter rows",
    "combination": "the combination",
}


def _counted(count: int, singular: str, plural: str) -> str:
    return f"{count} {singular if count == 1 else plural}"


def _fit_counts(fit: LabelFit) -> str:
    return (
        f"{_counted(fit.overlap_pairs, 'label pair overlaps', 'label pairs overlap')}, "
        f"{_counted(fit.crossing_pairs, 'leader pair crosses', 'leader pairs cross')}"
    )


def _cli_arguments(values: Mapping[str, float]) -> str:
    return " ".join(
        f"{option} {values[path]:g}" for path, option in _CLI_OPTIONS.items() if path in values
    )


def _joined(names: Sequence[str]) -> str:
    if len(names) <= 1:
        return "".join(names)
    return ", ".join(names[:-1]) + " and " + names[-1]


def format_circular_label_fit_report(report: LabelFitReport) -> str:
    """The CLI and Python text of the error and the warning (OD-F5, OD-F8)."""

    outcome = "Drawn anyway." if report.drawn else "Nothing was drawn."
    if len(report.records) == 1:
        (record,) = report.records
        first = (
            f"{record.label_count} outer labels do not fit at {record.font_size:g} pt "
            f"({_fit_counts(record.fit)}). {outcome}"
        )
    else:
        records = "; ".join(
            f"{record.record_id} ({record.label_count} labels at {record.font_size:g} pt: "
            f"{_fit_counts(record.fit)})"
            for record in report.records
        )
        first = f"Outer labels do not fit in {len(report.records)} records: {records}. {outcome}"
    lines = [first]
    fitting = [(row, row.values) for row in report.rows if row.values is not None]
    if fitting:
        arguments = [_cli_arguments(values) for _, values in fitting]
        width = max(len(argument) for argument in arguments) + 3
        lines.append("Run again with one of:")
        lines.extend(
            f"  {argument:<{width}}({_ROW_TITLES[row.kind]})"
            for argument, (row, _) in zip(arguments, fitting)
        )
    missing = [row.kind for row in report.rows if row.status == "not_found"]
    skipped = [row.kind for row in report.rows if row.status == "not_computed"]
    if missing and (fitting or skipped):
        names = _joined([_ROW_NAMES[kind] for kind in missing])
        verb = "do" if len(missing) > 1 or missing == ["tighter_rows"] else "does"
        lines.append(f"{names[0].upper()}{names[1:]} {verb} not fit this map.")
    elif missing:
        lines.append("No single change of the label offsets, font size or label spacing fits this map.")
    if skipped:
        lines.append(
            f"Not searched within {WARN_SEARCH_SECONDS:g} s: "
            f"{_joined([_ROW_NAMES[kind] for kind in skipped])}."
        )
    return "\n".join(lines)


__all__ = [
    "FONT_SIZE_PATHS",
    "ROW_KINDS",
    "SPACING_PATH",
    "WARN_SEARCH_SECONDS",
    "X_OFFSET_PATH",
    "Y_OFFSET_PATH",
    "CircularLabelFitProbe",
    "CircularLabelFitSearch",
    "CircularLabelFitVerifier",
    "CircularLabelVerification",
    "CircularOuterLabelLayout",
    "check_circular_label_fit",
    "format_circular_label_fit_report",
]
