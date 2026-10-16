"""Numeric helpers whose results must not depend on the Python version or on Unicode digits."""

from __future__ import annotations

import argparse
import math
import re
from collections.abc import Iterable
from typing import SupportsFloat, SupportsIndex

# Python's ``\d`` and ``str.isdigit()`` accept other Unicode digits and ``int()`` reads the
# decimal ones (``int("١") == 1``). Every integer gbdraw reads from user or file input is
# ASCII; use these instead of ``\d``, ``isdigit()`` or a bare ``int()``. Floats read from
# user or file input go through :func:`ascii_float` instead of a bare ``float()``.
ASCII_NONNEGATIVE_INTEGER_PATTERN = "[0-9]+"
_ASCII_NONNEGATIVE_INTEGER_RE = re.compile(ASCII_NONNEGATIVE_INTEGER_PATTERN)
_ASCII_INTEGER_RE = re.compile(f"[+-]?{ASCII_NONNEGATIVE_INTEGER_PATTERN}")
# JavaScript's ``Number.MAX_SAFE_INTEGER``: the largest integer a JavaScript reader of the same value keeps exactly.
JS_MAX_SAFE_INTEGER = 2**53 - 1


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


def ascii_float(value: str | SupportsFloat | SupportsIndex) -> float:
    """``float(value)``, except that text with another digit or an underscore is refused.

    ``float()`` reads the decimal digits of every script (``float("١.٥") == 1.5``)
    and digit-grouping underscores (``float("1_0") == 10.0``); such text raises
    ``float()``'s own ``ValueError``. Everything else ``float()`` accepts is
    kept: a sign, surrounding whitespace, an exponent, ``.5`` and ``5.``, and
    ``inf``, ``Infinity`` and ``nan`` in any case. A value that is not text goes
    to ``float()`` unchanged.
    """

    if isinstance(value, str) and (
        "_" in value or (not value.isascii() and any(char.isdecimal() and not char.isascii() for char in value))
    ):
        raise ValueError(f"could not convert string to float: {value!r}")
    return float(value)


def parse_ascii_float(text: str) -> float | None:
    """The float :func:`ascii_float` reads from ``text``, or ``None``."""

    try:
        return ascii_float(text)
    except ValueError:
        return None


def argparse_float(text: str) -> float:
    """argparse ``type=`` of a float option: what ``float()`` read, in ASCII only.

    Other digits and underscores (``"1_0"``) are refused with argparse's own
    message for a ``float`` option.
    """

    value = parse_ascii_float(text)
    if value is None:
        raise argparse.ArgumentTypeError(f"invalid float value: {text!r}")
    return value


__all__ = [
    "ASCII_NONNEGATIVE_INTEGER_PATTERN",
    "JS_MAX_SAFE_INTEGER",
    "argparse_float",
    "argparse_integer",
    "ascii_float",
    "compensated_sum",
    "parse_ascii_float",
    "parse_ascii_integer",
    "parse_ascii_nonnegative_integer",
    "scaled_tick_text",
]
