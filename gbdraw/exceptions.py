"""Custom exceptions for gbdraw."""

from __future__ import annotations

from typing import TYPE_CHECKING, Mapping

if TYPE_CHECKING:
    from gbdraw.labels.fit import LabelCapacityRecord, LabelFitReport, LabelFitRows, LabelKind


class GbdrawError(Exception):
    """Base class for gbdraw exceptions.

    ``diagnostic`` is the producer-owned failure meaning: a bounded mapping with
    ``code`` plus optional ``reason``, ``field``, and integer locators. It never
    carries user wording; the Web adapter validates it against its vocabulary.
    """

    diagnostic: Mapping[str, object] | None = None

    def __init__(
        self,
        *args: object,
        diagnostic: Mapping[str, object] | None = None,
    ) -> None:
        super().__init__(*args)
        if diagnostic is not None:
            self.diagnostic = dict(diagnostic)


class ConfigError(GbdrawError):
    """Raised when required configuration data is missing or invalid."""


class InputFileError(GbdrawError):
    """Raised when an input file is missing or unreadable."""


class ParseError(GbdrawError):
    """Raised when input data cannot be parsed."""


class ValidationError(GbdrawError):
    """Raised when user input or data fails validation."""


class LabelLayoutError(ValidationError):
    """A label setting the engine cannot draw (radial labels that find no room).

    The label fit counts it as a setting that does not fit, as it counts a
    ``TRACK_LAYOUT`` error of a row's assembly (the row's settings leave the
    tracks no room). Any other error of an assembly is raised.
    """


class ComparisonIdentityError(ValueError, GbdrawError):
    """Comparison endpoints conflict; remains catchable as ValueError.

    A conflict the user corrects in the input names it with ``diagnostic``.
    """

    def __init__(
        self,
        message: str,
        *,
        reason: str,
        diagnostic: Mapping[str, object] | None = None,
    ):
        super().__init__(message)
        self.reason = reason
        if diagnostic is not None:
            self.diagnostic = dict(diagnostic)


class ExportError(GbdrawError):
    """Raised when an explicitly requested library export cannot be generated."""


# The Web error's context key of the room each side's labels need when they were not placed (CAPACITY-STOP).
_CAPACITY_CONTEXT_KEYS: Mapping[LabelKind, str] = {
    "circular_outer": "labelCapacityOuter",
    "circular_inner": "labelCapacityInner",
}


def _label_fit_diagnostic(report: LabelFitReport) -> dict[str, object]:
    # The code alone names the failure; the shared CANNOT_FIT reason carries track wording.
    from gbdraw.labels.fit import label_fit_error_counts, largest_capacities, shown_capacity_ratio

    largest = largest_capacities(report)
    if largest:
        # Labels that need much more room than the map has were not placed: there is no fix to choose (CAPACITY-STOP).
        # The Web words the room each side needs from its value (OV-2179).
        return {
            "code": "LABEL_LAYOUT",
            "reason": "LABEL_CAPACITY",
            **{_CAPACITY_CONTEXT_KEYS[kind]: shown_capacity_ratio(capacity) for kind, capacity in largest.items()},
        }
    return {"code": "LABEL_LAYOUT", "labelFitCounts": label_fit_error_counts(report)}


class LabelsDoNotFitError(GbdrawError):
    """Labels would overlap, or radial labels have no room, so nothing was drawn.

    ``report`` names the failing records and the settings that fit;
    ``row_search``, when given, computes a row the report left
    ``not_computed``.
    """

    def __init__(
        self,
        message: str,
        *,
        report: LabelFitReport,
        row_search: LabelFitRows | None = None,
    ) -> None:
        super().__init__(message, diagnostic=_label_fit_diagnostic(report))
        self.report = report
        self.row_search = row_search


class LabelLeadersCrossWarning(UserWarning):
    """Label leaders cross, labels cover the definition text, or the Circular definition text overlaps a track;
    the figure was drawn anyway.

    ``report`` names the records and the settings that remove the crossings.
    """

    def __init__(self, message: str, *, report: LabelFitReport) -> None:
        super().__init__(message)
        self.report = report


class LabelsMayNotFitWarning(UserWarning):
    """A record's Circular labels need much of the room around the map, so they may not fit and placing them can take
    minutes (CAPACITY-WARN); warned before they are placed.

    ``capacity`` is the room they need; ``record_id`` and ``label_kind`` name them.
    """

    def __init__(self, message: str, *, record_id: str, label_kind: LabelKind, capacity: LabelCapacityRecord) -> None:
        super().__init__(message)
        self.record_id = record_id
        self.label_kind = label_kind
        self.capacity = capacity


__all__ = [
    "ComparisonIdentityError",
    "ConfigError",
    "ExportError",
    "GbdrawError",
    "InputFileError",
    "LabelLeadersCrossWarning",
    "LabelsDoNotFitError",
    "LabelsMayNotFitWarning",
    "ParseError",
    "ValidationError",
]
