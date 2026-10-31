"""The owner of a Linear wrap length: its automatic rule and what it refuses.

Wrapping is a display format of one record: rows live in the record display
transform, never as records. Only a record longer than one row is cut, and every
record of a figure uses one row length, so all rows draw at one scale.
"""

from __future__ import annotations

from collections.abc import Sequence
from math import ceil

from gbdraw.config.models.canvas import LinearWrapLength
from gbdraw.exceptions import ValidationError


AUTO_WRAP_CAP_BP = 50_000
AUTO_WRAP_FLOOR_BP = 10_000
AUTO_WRAP_MAX_ROWS = 5
# The automatic rule's candidates, from the floor to the cap, on the 1-2-5 ladder.
_AUTO_WRAP_LADDER_BP = (10_000, 20_000, AUTO_WRAP_CAP_BP)


def auto_wrap_length(longest_bp: int) -> int | None:
    """The automatic row length for a figure whose longest record is ``longest_bp``.

    The smallest ladder value that keeps that record within ``AUTO_WRAP_MAX_ROWS``, and the cap
    when none does. A record no longer than the floor stays one row (``None``).
    """

    if longest_bp <= AUTO_WRAP_FLOOR_BP:
        return None
    for length in _AUTO_WRAP_LADDER_BP:
        if ceil(longest_bp / length) <= AUTO_WRAP_MAX_ROWS:
            return length
    return AUTO_WRAP_CAP_BP


def resolve_wrap_length(requested: LinearWrapLength, record_lengths: Sequence[int]) -> int | None:
    """The row length to cut at, or ``None`` when no record is longer than one row."""

    if requested == "none":
        return None
    longest = max(record_lengths)
    length = auto_wrap_length(longest) if requested == "auto" else requested
    return length if length is not None and longest > length else None


def wrap_blocking_inputs(
    *,
    blast_files: object = None,
    comparison_table_file: object = None,
    linear_comparisons: object = None,
    protein_comparisons: object = None,
    orthogroups: object = None,
    collinearity_blocks: object = None,
    losat_search: object = None,
    similarity_alignment: object = None,
) -> tuple[str, ...]:
    """The names of the given comparison inputs, none of which a wrapped figure can draw."""

    given = {
        "blast_files": blast_files,
        "comparison_table_file": comparison_table_file,
        "linear_comparisons": linear_comparisons,
        "protein_comparisons": protein_comparisons,
        "orthogroups": orthogroups,
        "collinearity_blocks": collinearity_blocks,
        "losat_search": losat_search,
        "similarity_alignment": similarity_alignment,
    }
    return tuple(name for name, value in given.items() if value is not None and value != () and value != [])


def reject_wrap_with_comparisons(requested: LinearWrapLength, blocking: Sequence[str]) -> None:
    """Refuse comparison inputs while wrapping, before any search runs."""

    if requested != "none" and blocking:
        raise ValidationError(
            "--wrap_length cannot be combined with comparison inputs "
            f"({', '.join(blocking)}): a wrapped figure draws no comparisons.",
            diagnostic={"code": "COMPARISON_INPUT", "field": "wrap_length", "reason": "WRAP"},
        )


def reject_wrap_in_one_grid_row(requested: LinearWrapLength, rows_by_record: Sequence[int]) -> None:
    """Refuse several records in one grid row while wrapping: a row holds one record's rows."""

    if requested != "none" and len(set(rows_by_record)) != len(rows_by_record):
        raise ValidationError(
            "--wrap_length cannot draw several records in one grid row; "
            "give each record its own row, or do not wrap.",
            diagnostic={"code": "INPUT_INVALID", "field": "wrap_length", "reason": "WRAP_GRID_ROW"},
        )


__all__ = [
    "AUTO_WRAP_CAP_BP",
    "AUTO_WRAP_FLOOR_BP",
    "AUTO_WRAP_MAX_ROWS",
    "auto_wrap_length",
    "reject_wrap_in_one_grid_row",
    "reject_wrap_with_comparisons",
    "resolve_wrap_length",
    "wrap_blocking_inputs",
]
