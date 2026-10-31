"""The owner of a Linear wrap length: its automatic rule, whether a record is cut, and what that refuses.

Wrapping is a display format of one record: rows live in the record display
transform, never as records. Only a record longer than one row is cut, and every
record of a figure uses one row length, so all rows draw at one scale.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from Bio.SeqRecord import SeqRecord

from gbdraw.config.models.canvas import LinearWrapLength
from gbdraw.config.models.root import GbdrawConfig
from gbdraw.core.record_metadata import _read_coord_map
from gbdraw.exceptions import ValidationError
from gbdraw.layout.linear import VerticalBand
from gbdraw.layout.record_coordinates import RecordDisplayTransform, wrapped_row_count


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
        if wrapped_row_count(longest_bp, length) <= AUTO_WRAP_MAX_ROWS:
            return length
    return AUTO_WRAP_CAP_BP


def resolve_wrap_length(requested: LinearWrapLength, record_lengths: Sequence[int]) -> int | None:
    """The row length to cut at, or ``None`` when no record is longer than one row.

    This is the one test of "a record is cut". Everything the option does beyond the option's own validation (rows,
    their defaults, its refusals) follows it, so with nothing to cut the option has no effect at all.
    """

    if requested == "none":
        return None
    longest = max(record_lengths)
    length = auto_wrap_length(longest) if requested == "auto" else requested
    return length if length is not None and longest > length else None


def _length_words(length: int) -> str:
    """A row length for a person: ``10 kb``, ``1.5 Mb``, ``800 bp``."""

    if length >= 1_000_000:
        return f"{length / 1_000_000:g} Mb"
    return f"{length / 1_000:g} kb" if length >= 1_000 else f"{length} bp"


def wrap_length_option(length: int) -> str:
    """A row length as the option a person types: ``--wrap_length 100kb``."""

    if length >= 1_000_000:
        return f"--wrap_length {length / 1_000_000:g}Mb"
    return f"--wrap_length {length / 1_000:g}kb" if length >= 1_000 else f"--wrap_length {length}"


def longer_row_lengths(row_length: int, longest_bp: int) -> tuple[int, ...]:
    """The 1-2-5 row lengths above ``row_length``, up to the first that holds the longest record on one row."""

    lengths: list[int] = []
    decade = 10 ** (len(str(row_length)) - 1)
    while True:
        for length in (decade, 2 * decade, 5 * decade):
            if length > row_length:
                lengths.append(length)
                if length >= longest_bp:
                    return tuple(lengths)
        decade *= 10


@dataclass(frozen=True)
class LinearWrapSummary:
    """What a figure drew: its row length, what was asked for, and each record's length; rows follow from them."""

    row_length: int
    requested: LinearWrapLength
    record_lengths: tuple[int, ...]

    def rows_at(self, row_length: int) -> int:
        return sum(wrapped_row_count(length, row_length) for length in self.record_lengths)

    @property
    def row_count(self) -> int:
        return self.rows_at(self.row_length)

    @property
    def record_rows(self) -> tuple[int, ...]:
        return tuple(wrapped_row_count(length, self.row_length) for length in self.record_lengths)

    @property
    def message(self) -> str:
        auto = " (auto)" if self.requested == "auto" else ""
        return f"Wrapping at {_length_words(self.row_length)} per row{auto}: {self.row_count} rows"


def drawn_axis_length(record: SeqRecord, transform: RecordDisplayTransform | None) -> int:
    """The bases one drawn row of a record shows: the only source of a record's axis length.

    A record's window, definition text and source coordinates read ``len(record.seq)`` instead.
    """

    return len(record.seq) if transform is None or transform.wrap_length is None else transform.axis_length


def _unset_transform(record: SeqRecord) -> RecordDisplayTransform:
    """The transform of a record no plan bound one to, from its coordinate map."""

    base, step = _read_coord_map(record)
    length = len(record.seq)
    if step == 1:
        return RecordDisplayTransform(base + length - 1, base, 1, None, False)
    return RecordDisplayTransform(max(base, length), base, -1, None, False)


def resolve_wrap(
    requested: LinearWrapLength,
    records: Sequence[SeqRecord],
    transforms: Sequence[RecordDisplayTransform] | None,
) -> tuple[Sequence[RecordDisplayTransform] | None, LinearWrapSummary | None]:
    """The one place a wrap length becomes transforms: only a record longer than the row is cut."""

    lengths = [len(record.seq) for record in records]
    row_length = resolve_wrap_length(requested, lengths)
    if row_length is None:
        return transforms, None
    bound = transforms if transforms is not None else [_unset_transform(record) for record in records]
    wrapped = tuple(
        transform.with_wrap(row_length, length) if length > row_length else transform
        for transform, length in zip(bound, lengths, strict=True)
    )
    return wrapped, LinearWrapSummary(row_length, requested, tuple(lengths))


def stack_row_offsets(row_bands: Sequence[VerticalBand], gap: float) -> tuple[float, ...]:
    """Each row's axis offset from the first row's axis, each row taking its own height (``row_bands``,
    local to its axis) with ``gap`` between neighbours."""

    offsets = [0.0]
    for above, below in zip(row_bands, row_bands[1:]):
        offsets.append(offsets[-1] + above.bottom_y - below.top_y + gap)
    return tuple(offsets)


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


def refuse_undrawable_wrap(
    config: GbdrawConfig, rows_by_record: Sequence[int], blocking: Sequence[str]
) -> None:
    """Refuse what rows cannot draw. Called only when a record is cut: with nothing to cut, the option has no effect."""

    if blocking:
        raise ValidationError(
            "wrap_length cannot cut a record of a figure with comparison inputs "
            f"({', '.join(blocking)}): a wrapped figure draws no comparisons.",
            diagnostic={"code": "COMPARISON_INPUT", "field": "wrap_length", "reason": "WRAP"},
        )
    linear = config.canvas.linear
    scale = config.objects.scale
    if linear.track_layout == "middle":
        raise _wrap_conflict(
            "WRAP_TRACK_LAYOUT",
            "track_layout middle cannot draw a record cut by wrap_length: "
            "rows draw their tracks above or below the axis.",
        )
    if scale.show and scale.style != "ruler":
        raise _wrap_conflict(
            "WRAP_SCALE_STYLE",
            f"scale_style {scale.style} cannot draw a record cut by wrap_length: each row draws a ruler on its axis.",
        )
    if linear.normalize_length:
        raise _wrap_conflict(
            "WRAP_NORMALIZE_LENGTH",
            "normalize_length cannot draw a record cut by wrap_length: every row draws at one scale.",
        )
    if len(set(rows_by_record)) != len(rows_by_record):
        raise _wrap_conflict(
            "WRAP_GRID_ROW",
            "wrap_length cannot cut a record in one grid row with others; give each record its own row, or do not wrap.",
        )


def refuse_feature_track_on_the_axis(feature_track_layout: str) -> None:
    """Refuse a wrapped figure whose features stand in the middle, whatever put them there: that is the middle layout
    under another name, and a ruler on the axis of the rows could not stand there."""

    if feature_track_layout == "middle":
        raise _wrap_conflict(
            "WRAP_TRACK_LAYOUT",
            "A feature track on the axis (side=overlay) cannot draw a record cut by wrap_length: "
            "rows draw their tracks above or below the axis.",
        )


def _wrap_conflict(reason: str, message: str) -> ValidationError:
    return ValidationError(message, diagnostic={"code": "INPUT_INVALID", "field": "wrap_length", "reason": reason})


__all__ = [
    "AUTO_WRAP_CAP_BP",
    "AUTO_WRAP_FLOOR_BP",
    "AUTO_WRAP_MAX_ROWS",
    "LinearWrapSummary",
    "auto_wrap_length",
    "drawn_axis_length",
    "longer_row_lengths",
    "refuse_feature_track_on_the_axis",
    "refuse_undrawable_wrap",
    "resolve_wrap",
    "resolve_wrap_length",
    "stack_row_offsets",
    "wrap_length_option",
    "wrap_blocking_inputs",
]
