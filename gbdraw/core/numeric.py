"""Numeric helpers whose results must not depend on the Python version or on Unicode digits."""

from __future__ import annotations

import argparse
import math
import re
from collections.abc import Iterable
from fractions import Fraction
from typing import SupportsFloat, SupportsIndex, SupportsInt

# Python's ``\d`` and ``str.isdigit()`` accept other Unicode digits and ``int()`` reads the
# decimal ones (``int("١") == 1``). Every integer gbdraw reads from user or file input is
# ASCII; use these instead of ``\d``, ``isdigit()`` or a bare ``int()``. Floats read from
# user or file input go through :func:`ascii_float` instead of a bare ``float()``.
ASCII_NONNEGATIVE_INTEGER_PATTERN = "[0-9]+"
_ASCII_NONNEGATIVE_INTEGER_RE = re.compile(ASCII_NONNEGATIVE_INTEGER_PATTERN)
_ASCII_INTEGER_RE = re.compile(f"[+-]?{ASCII_NONNEGATIVE_INTEGER_PATTERN}")
# JavaScript's ``Number.MAX_SAFE_INTEGER``: the largest integer a JavaScript reader of the same value keeps exactly.
JS_MAX_SAFE_INTEGER = 2**53 - 1
# The white space and line terminators JavaScript's ``String.prototype.trim()``
# removes, for Python readers that read text as the Web does.
JS_WHITESPACE = (
    "\t\n\v\f\r \u00a0\u1680\u2000\u2001\u2002\u2003\u2004\u2005\u2006"
    "\u2007\u2008\u2009\u200a\u2028\u2029\u202f\u205f\u3000\ufeff"
)


def compensated_sum(values: Iterable[float]) -> float:
    """Sum floats as CPython 3.12+ ``sum()`` does, on every supported Python.

    Python 3.12 made the built-in ``sum()`` of floats compensated (Neumaier's
    variant of Kahan summation, ``Python/bltinmodule.c``); 3.10 and 3.11 add
    naively, so the same layout came out one ulp apart. The Web port of the
    Legend layout (``pythonFloatSum`` in ``gbdraw/web/js/services/legend-layout.js``)
    computes exactly this: the first value added to ``0``, each further value
    compensated, and the compensation added at the end when it is nonzero and
    finite. Every value is taken as a float, as the port does.
    """
    iterator = iter(values)
    try:
        first = next(iterator)
    except StopIteration:
        return 0.0
    total = 0 + float(first)
    compensation = 0.0
    for value in iterator:
        item = float(value)
        partial = total + item
        if abs(total) >= abs(item):
            compensation += (total - partial) + item
        else:
            compensation += (item - partial) + total
        total = partial
    if compensation and math.isfinite(compensation):
        total += compensation
    return total


def scaled_tick_text(position: int, divisor: int, interval: int | None) -> str:
    """Return ``position / divisor`` with the decimals that ``interval`` needs.

    The decimals are the fewest (at most 6) that write ``interval / divisor``
    exactly, so a 250 bp interval in kbp gives ``0.25``, ``0.5``, ``0.75``;
    trailing zeros are dropped. Without an interval the value is rounded.
    """
    decimals = 0
    step = abs(int(interval or 0))
    while step and decimals < 6 and (step * 10**decimals) % divisor:
        decimals += 1
    value = float(position) / float(divisor)
    if decimals == 0:
        return f"{value:.0f}"
    return f"{value:.{decimals}f}".rstrip("0").rstrip(".")


def _parse_ascii(pattern: re.Pattern[str], text: str) -> int | None:
    if pattern.fullmatch(text) is None:
        return None
    try:
        return int(text)
    except ValueError:  # more digits than int() converts (4300)
        return None


def parse_ascii_nonnegative_integer(text: str) -> int | None:
    """Return the integer ``text`` spells in ASCII digits, else ``None``.

    The whole text must be digits: no sign, space, underscore, or non-ASCII
    digit. A number too long for ``int()`` is ``None`` as well.
    """
    return _parse_ascii(_ASCII_NONNEGATIVE_INTEGER_RE, text)


def parse_ascii_integer(text: str) -> int | None:
    """Like :func:`parse_ascii_nonnegative_integer`, with an optional ``+`` or ``-`` sign."""
    return _parse_ascii(_ASCII_INTEGER_RE, text)


def argparse_integer(text: str) -> int:
    """argparse ``type=`` of an integer option: what ``int()`` read, in ASCII digits only.

    An optional sign and surrounding spaces are accepted, as ``int()`` did; other
    digits and underscores (``"1_0"``) are refused with ``int()``'s own message.
    """
    value = parse_ascii_integer(text.strip())
    if value is None:
        raise argparse.ArgumentTypeError(f"invalid int value: {text!r}")
    return value


def _has_other_digit_or_underscore(value: object) -> bool:
    return isinstance(value, str) and (
        "_" in value or (not value.isascii() and any(char.isdecimal() and not char.isascii() for char in value))
    )


def ascii_int(value: str | SupportsInt | SupportsIndex) -> int:
    """``int(value)``, except that text with another digit or an underscore is refused.

    The integer half of :func:`ascii_float`, for readers that take a number or
    its text (config leaves): such text raises ``int()``'s own ``ValueError``.
    """

    if _has_other_digit_or_underscore(value):
        raise ValueError(f"invalid literal for int() with base 10: {value!r}")
    return int(value)


def ascii_float(value: str | SupportsFloat | SupportsIndex) -> float:
    """``float(value)``, except that text with another digit or an underscore is refused.

    ``float()`` reads the decimal digits of every script (``float("١.٥") == 1.5``)
    and digit-grouping underscores (``float("1_0") == 10.0``); such text raises
    ``float()``'s own ``ValueError``. Everything else ``float()`` accepts is
    kept: a sign, surrounding whitespace, an exponent, ``.5`` and ``5.``, and
    ``inf``, ``Infinity`` and ``nan`` in any case. A value that is not text goes
    to ``float()`` unchanged.
    """

    if _has_other_digit_or_underscore(value):
        raise ValueError(f"could not convert string to float: {value!r}")
    return float(value)


def parse_ascii_float(text: str) -> float | None:
    """The float :func:`ascii_float` reads from ``text``, or ``None``."""

    try:
        return ascii_float(text)
    except ValueError:
        return None


def parse_ascii_finite_float(text: str) -> float | None:
    """The finite float :func:`ascii_float` reads from ``text``, or ``None``.

    ``nan``, ``inf`` and their other spellings are ``None``: no gbdraw option
    or table number takes a non-finite value.
    """

    value = parse_ascii_float(text)
    return value if value is not None and math.isfinite(value) else None


# A length in bp: a decimal count and an optional unit, in ASCII. ``float()`` is not used: the
# value must be an exact whole number of bases, and its spellings (``1e3``) are not lengths.
_LENGTH_RE = re.compile(r"\s*([0-9]+(?:\.[0-9]*)?|\.[0-9]+)\s*(bp|k|kb|kbp|m|mb|mbp)?\s*", re.IGNORECASE | re.ASCII)
_LENGTH_UNIT_BP = {"": 1, "bp": 1, "k": 10**3, "kb": 10**3, "kbp": 10**3, "m": 10**6, "mb": 10**6, "mbp": 10**6}


def parse_length_bp(text: str) -> int | None:
    """The positive whole number of bases ``text`` spells, with an optional unit, else ``None``.

    ``50000``, ``50kb``, ``50 kbp``, ``50k``, ``0.5Mb``; the unit's case does not matter. A value
    that is not a whole base (``0.0004kb``), is zero, or exceeds what a JavaScript reader keeps
    exactly is ``None``, as are a sign, an exponent, and other digits.
    """

    match = _LENGTH_RE.fullmatch(text) if isinstance(text, str) else None
    if match is None or len(match[1]) > 32:
        return None
    value = Fraction(match[1]) * _LENGTH_UNIT_BP[(match[2] or "").lower()]
    return int(value) if value.denominator == 1 and 0 < value <= JS_MAX_SAFE_INTEGER else None


def argparse_length_bp(text: str) -> int:
    """argparse ``type=`` of a length option (:func:`parse_length_bp`)."""

    value = parse_length_bp(text)
    if value is None:
        raise argparse.ArgumentTypeError(
            f"invalid length: {text!r} (a positive whole number of bases, with an optional bp, kb, or Mb unit)"
        )
    return value


def argparse_scale_interval(text: str) -> int:
    """argparse ``type=`` of ``--scale_interval``: a length, or any integer as the option always read.

    The caller refuses a value <= 0 or reads it as a released legacy Session drew it.
    """

    value = parse_ascii_integer(text.strip())
    return value if value is not None else argparse_length_bp(text)


def argparse_float(text: str) -> float:
    """argparse ``type=`` of every float option: a finite ``float()`` in ASCII only.

    Other digits and underscores (``"1_0"``) are refused with argparse's own
    message for a ``float`` option; ``nan`` and ``inf`` are refused as not finite.
    """

    value = parse_ascii_float(text)
    if value is None:
        raise argparse.ArgumentTypeError(f"invalid float value: {text!r}")
    if not math.isfinite(value):
        raise argparse.ArgumentTypeError(f"invalid float value: {text!r} (must be a finite number)")
    return value


__all__ = [
    "ASCII_NONNEGATIVE_INTEGER_PATTERN",
    "JS_MAX_SAFE_INTEGER",
    "JS_WHITESPACE",
    "argparse_float",
    "argparse_integer",
    "argparse_length_bp",
    "argparse_scale_interval",
    "ascii_float",
    "ascii_int",
    "compensated_sum",
    "parse_ascii_finite_float",
    "parse_ascii_float",
    "parse_ascii_integer",
    "parse_ascii_nonnegative_integer",
    "parse_length_bp",
    "scaled_tick_text",
]
