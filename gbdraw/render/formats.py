"""Shared output format parsing and filename helpers."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Iterable, Sequence

SVG_FORMAT = "svg"
INTERACTIVE_SVG_FORMAT = "interactive_svg"
PNG_FORMAT = "png"
PDF_FORMAT = "pdf"
EPS_FORMAT = "eps"
PS_FORMAT = "ps"

CAIROSVG_FORMATS = (PNG_FORMAT, PDF_FORMAT, EPS_FORMAT, PS_FORMAT)
BASE_SVG_FORMATS = (SVG_FORMAT, INTERACTIVE_SVG_FORMAT)
ACCEPTED_FORMATS = (*BASE_SVG_FORMATS, *CAIROSVG_FORMATS)

_FORMAT_ALIASES = {
    "interactive-svg": INTERACTIVE_SVG_FORMAT,
}


@dataclass(frozen=True)
class FormatClassification:
    """Classified format tokens used by export paths."""

    requested: tuple[str, ...]
    base_svg: tuple[str, ...]
    cairosvg: tuple[str, ...]
    interactive: bool


def normalize_format_token(value: object) -> str:
    """Normalize one user-facing format token to its canonical form."""

    token = str(value).strip().lower().lstrip(".")
    return _FORMAT_ALIASES.get(token, token)


def dedupe_formats(formats: Iterable[str]) -> list[str]:
    """Return formats in first-seen order with duplicates removed."""

    result: list[str] = []
    seen: set[str] = set()
    for fmt in formats:
        if fmt in seen:
            continue
        seen.add(fmt)
        result.append(fmt)
    return result


def parse_format_string(
    out_formats: str,
    *,
    logger: logging.Logger | None = None,
    fallback: Sequence[str] = (PNG_FORMAT,),
) -> list[str]:
    """Parse a comma-separated output format string."""

    parsed: list[str] = []
    for raw in str(out_formats).split(","):
        fmt = normalize_format_token(raw)
        if fmt in ACCEPTED_FORMATS:
            parsed.append(fmt)
        else:
            if logger is not None:
                logger.warning(
                    "WARNING: Unaccepted/unrecognized output file format: %s",
                    fmt,
                )

    accepted = dedupe_formats(parsed)
    if not accepted:
        if logger is not None:
            logger.warning(
                "WARNING: No valid output file format was specified; generate a PNG file."
            )
        accepted = list(fallback)
    return accepted


def normalize_format_sequence(
    formats: Sequence[str] | str | None,
    *,
    logger: logging.Logger | None = None,
    none_default: Sequence[str] = (SVG_FORMAT,),
    fallback: Sequence[str] = (PNG_FORMAT,),
) -> list[str]:
    """Normalize API format input while preserving existing fallback behavior."""

    if formats is None:
        return list(none_default)
    if isinstance(formats, str):
        return parse_format_string(formats, logger=logger, fallback=fallback)

    parsed: list[str] = []
    for raw in formats:
        if raw is None:
            continue
        fmt = normalize_format_token(raw)
        if not fmt:
            continue
        if fmt not in ACCEPTED_FORMATS:
            if logger is not None:
                logger.warning(
                    "WARNING: Unaccepted/unrecognized output file format: %s",
                    fmt,
                )
            continue
        parsed.append(fmt)

    accepted = dedupe_formats(parsed)
    if not accepted:
        if logger is not None:
            logger.warning(
                "WARNING: No valid output file format was specified; generate a PNG file."
            )
        accepted = list(fallback)
    return accepted


def classify_formats(formats: Sequence[str]) -> FormatClassification:
    """Classify normalized output formats."""

    requested = tuple(dedupe_formats(normalize_format_token(fmt) for fmt in formats))
    base_svg = tuple(fmt for fmt in requested if fmt in BASE_SVG_FORMATS)
    cairosvg = tuple(fmt for fmt in requested if fmt in CAIROSVG_FORMATS)
    return FormatClassification(
        requested=requested,
        base_svg=base_svg,
        cairosvg=cairosvg,
        interactive=INTERACTIVE_SVG_FORMAT in requested,
    )


def is_cairosvg_format(fmt: str) -> bool:
    return normalize_format_token(fmt) in CAIROSVG_FORMATS


def resolve_format_output_path(base_prefix: str, fmt: str) -> str:
    """Resolve an output path for one canonical format."""

    normalized = normalize_format_token(fmt)
    if normalized == INTERACTIVE_SVG_FORMAT:
        return f"{base_prefix}.interactive.svg"
    if normalized in ACCEPTED_FORMATS:
        return f"{base_prefix}.{normalized}"
    raise ValueError(f"Unsupported output format: {fmt}")


def resolve_output_paths(
    base_prefix: str,
    formats: Sequence[str],
    *,
    include_base_svg: bool = True,
) -> list[str]:
    """Resolve output paths for requested formats, including the base SVG once."""

    paths: list[str] = []
    if include_base_svg:
        paths.append(resolve_format_output_path(base_prefix, SVG_FORMAT))
    for fmt in dedupe_formats(normalize_format_token(value) for value in formats):
        if include_base_svg and fmt == SVG_FORMAT:
            continue
        paths.append(resolve_format_output_path(base_prefix, fmt))
    return dedupe_formats(paths)


__all__ = [
    "ACCEPTED_FORMATS",
    "BASE_SVG_FORMATS",
    "CAIROSVG_FORMATS",
    "EPS_FORMAT",
    "FormatClassification",
    "INTERACTIVE_SVG_FORMAT",
    "PDF_FORMAT",
    "PNG_FORMAT",
    "PS_FORMAT",
    "SVG_FORMAT",
    "classify_formats",
    "dedupe_formats",
    "is_cairosvg_format",
    "normalize_format_sequence",
    "normalize_format_token",
    "parse_format_string",
    "resolve_format_output_path",
    "resolve_output_paths",
]
