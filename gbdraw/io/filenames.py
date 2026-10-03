"""File names of Session resources and written evidence files."""

from __future__ import annotations

import re
from typing import Callable, Collection, Sequence

_SAFE_FILENAME_RE = re.compile(r"[^A-Za-z0-9._-]+")


def safe_embedded_filename(name: object, *, fallback: str = "embedded-file") -> str:
    """Return a basename-only filename safe for materializing embedded content."""

    raw_name = str(name or "").replace("\\", "/").split("/")[-1].strip()
    cleaned = _SAFE_FILENAME_RE.sub("_", raw_name).strip("._")
    return cleaned or fallback


def unique_filename(
    name: str,
    used: Collection[str],
    *,
    key: Callable[[str], str] = str,
) -> str:
    """``name``, or ``<stem>.2.<ext>``, ``<stem>.3.<ext>``, ... when ``key(name)`` is used.

    The number goes before the last extension (``X.fna`` -> ``X.2.fna``); a
    name without one gets it at the end (``X`` -> ``X.2``).
    """

    stem, dot, extension = name.rpartition(".")
    candidate, ordinal = name, 1
    while key(candidate) in used:
        ordinal += 1
        candidate = f"{stem}.{ordinal}.{extension}" if dot and stem else f"{name}.{ordinal}"
    return candidate


def unique_filenames(names: Sequence[str], *, reserved: Sequence[str] = ()) -> tuple[str, ...]:
    """Apply :func:`unique_filename` to ``names`` in order, after ``reserved``."""

    used = set(reserved)
    result = []
    for name in names:
        candidate = unique_filename(name, used)
        used.add(candidate)
        result.append(candidate)
    return tuple(result)


__all__ = ["safe_embedded_filename", "unique_filename", "unique_filenames"]
