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


__all__ = [
    "ConfigError",
    "GbdrawError",
    "InputFileError",
    "ParseError",
    "ValidationError",
]
