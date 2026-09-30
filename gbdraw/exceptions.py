"""Custom exceptions for gbdraw."""

from __future__ import annotations

from typing import Mapping


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


class ComparisonIdentityError(ValueError, GbdrawError):
    """Comparison endpoints conflict; remains catchable as ValueError."""

    def __init__(self, message: str, *, reason: str):
        super().__init__(message)
        self.reason = reason


class ExportError(GbdrawError):
    """Raised when an explicitly requested library export cannot be generated."""


__all__ = [
    "ComparisonIdentityError",
    "ConfigError",
    "ExportError",
    "GbdrawError",
    "InputFileError",
    "ParseError",
    "ValidationError",
]
