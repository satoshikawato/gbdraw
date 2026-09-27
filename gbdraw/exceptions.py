"""Custom exceptions for gbdraw."""


class GbdrawError(Exception):
    """Base class for gbdraw exceptions."""


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
