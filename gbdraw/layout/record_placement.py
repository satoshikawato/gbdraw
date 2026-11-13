"""Shared record-row placement grammar and selector resolution."""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Literal, Sequence

from Bio.SeqRecord import SeqRecord

from gbdraw.core.numeric import parse_ascii_nonnegative_integer
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import RecordSelector, drawn_record_units, parse_record_selector


_PlacementCompatibility = Literal["linear", "circular"]

# The row of a ``<selector>@<row>`` Record Order token, the text after its last ``@``, as
# :func:`parse_record_row_position` reads it: spaces (``str.strip``'s, which are ``\s``'s) around it, an optional
# ``+``, and ASCII digits, leading zeros allowed; the row is that number when it is above zero. The Web reads every
# token's row with this pattern from the generated contracts (``RECORD_ROW_PATTERN``, OV-3313).
RECORD_ROW_PATTERN = r"\s*\+?([0-9]+)\s*"
_RECORD_ROW = re.compile(RECORD_ROW_PATTERN)

# What ``anchor_feature`` does to one record: ``rotated``, a complete circular record turned so that the feature
# starts at its left end; ``aligned``, any other record moved sideways so that its feature starts below the
# reference's; ``reference``, the record the others line up with when none is rotated, which keeps its place;
# ``unmatched``, a record without the feature, which keeps its place and its start.
FeatureAnchorStatus = Literal["reference", "rotated", "aligned", "unmatched"]


@dataclass(frozen=True)
class FeatureAnchorPlacement:
    """Where ``anchor_feature`` puts each record, decided once by the planner (``resolve_anchor_feature``).

    ``starts`` are where the feature starts in each drawn record, in display bp (0.0 for a rotated record; ``None``
    for an unmatched one). ``reference_index`` is the record the others line up with: the first rotated record,
    else the first record with the feature. The build reads them to give each record its x
    (``_final_record_translations``); a rotated reference is drawn at x 0.
    """

    starts: tuple[float | None, ...]
    statuses: tuple[FeatureAnchorStatus, ...]
    reference_index: int

    def __post_init__(self) -> None:
        # The planner's invariants, not a user's input: one start and status per record, an unmatched record
        # without a start, and a reference that is rotated or the reference.
        assert len(self.starts) == len(self.statuses)
        assert all((status == "unmatched") == (start is None) for start, status in zip(self.starts, self.statuses))
        assert self.statuses[self.reference_index] in ("reference", "rotated")

    @property
    def rotates(self) -> bool:
        """Whether the anchor rotated a record (then the request stores its new start)."""
        return "rotated" in self.statuses


def parse_record_row_position(value: str) -> tuple[str, int]:
    """Parse one ``<selector>@<row>`` placement token, in both modes.

    The row is the text after the last ``@``, a positive integer with an optional ``+``; the selector is the
    text before it, a ``#N`` that :func:`parse_record_selector` reads or a record ID (which may hold ``@``).
    """

    raw = str(value or "").strip()
    if not raw:
        raise ValidationError("multi_record_position does not allow empty entries.")
    if "@" not in raw:
        raise ValidationError(
            f"multi_record_position entry '{raw}' must be in '<selector>@<row>' format."
        )
    selector, row_text = (part.strip() for part in raw.rsplit("@", 1))
    if not selector:
        raise ValidationError(
            f"multi_record_position entry '{raw}' must include a selector before '@'."
        )
    match = _RECORD_ROW.fullmatch(row_text)
    row = parse_ascii_nonnegative_integer(match.group(1)) if match else None
    if not row:
        raise ValidationError(
            f"multi_record_position entry '{raw}' must use a positive integer row."
        )
    if selector.startswith("#"):
        _record_index_selector(selector)
    return selector, row


def record_order_selector(selector: str, records: Sequence[SeqRecord]) -> RecordSelector:
    """The record a Record Order selector names among the drawn ``records``, the one reader of it with or without
    windows, in both modes: a ``#N`` as :func:`parse_record_selector` reads it; any other text the record with that
    ID; a word :func:`parse_record_selector` reads as no selector (``none``, ``-``) that no record has as its ID, the
    one drawn record when only one is drawn (D-211). Any other text names no record, and resolving it refuses it."""

    if selector.startswith("#"):
        return _record_index_selector(selector)
    if (
        parse_record_selector(selector) is None
        and all(record.id != selector for record in records)
        and len(drawn_record_units(records)) == 1
    ):
        return RecordSelector(raw=selector, record_id=None, record_index=0)
    return RecordSelector(raw=selector, record_id=selector, record_index=None)


def _record_index_selector(selector: str) -> RecordSelector:
    """A ``#N`` Record Order selector as :func:`parse_record_selector` reads it, refused as a ``ValidationError``."""

    try:
        parsed = parse_record_selector(selector)
    except ValueError as exc:
        raise ValidationError(str(exc)) from exc
    # ``parse_record_selector`` refuses or reads every text that starts with ``#``.
    assert parsed is not None and parsed.record_index is not None
    return parsed


def _resolve_record_selector_index(records: Sequence[SeqRecord], selector: str) -> int:
    record_count = len(records)
    record_index = record_order_selector(selector, records).record_index
    if record_index is not None:
        if record_index < record_count:
            return record_index
        raise ValidationError(
            f"multi_record_position selector '{selector}' is out of range for "
            f"{record_count} loaded record(s)."
        )
    matches = [
        index
        for index, record in enumerate(records)
        if str(record.id) == selector
    ]
    if not matches:
        raise ValidationError(
            f"multi_record_position selector '{selector}' did not match any record ID."
        )
    if len(matches) > 1:
        raise ValidationError(
            f"multi_record_position selector '{selector}' matched multiple records. Use #index to disambiguate."
        )
    return matches[0]


def resolve_record_row_positions(
    records: Sequence[SeqRecord],
    positions: Sequence[str] | None,
    *,
    _compatibility: _PlacementCompatibility = "linear",
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Resolve selectors to input indices and contiguous zero-based rows."""

    record_count = len(records)
    if not positions:
        indices = tuple(range(record_count))
        return indices, indices
    if _compatibility == "linear" and len(positions) != record_count:
        raise ValidationError(
            f"multi_record_position must provide exactly {record_count} entry(ies)."
        )

    resolved: list[tuple[int, int, int]] = []
    seen: set[int] = set()
    for token_index, raw in enumerate(positions):
        selector, row = parse_record_row_position(str(raw))
        record_index = _resolve_record_selector_index(records, selector)
        if record_index in seen:
            raise ValidationError(
                f"multi_record_position selector '{selector}' was specified more than once."
            )
        seen.add(record_index)
        resolved.append((row, token_index, record_index))

    if len(seen) != record_count:
        if _compatibility == "circular":
            raise ValidationError(
                "multi_record_position must include each loaded record exactly once "
                f"(expected {record_count}, got {len(seen)} unique selector(s))."
            )
        raise ValidationError(
            "multi_record_position must include each loaded record exactly once."
        )
    if _compatibility == "circular" and len(positions) != record_count:
        raise ValidationError(
            f"multi_record_position must provide exactly {record_count} entry(ies)."
        )

    row_map = {
        value: index
        for index, value in enumerate(sorted({item[0] for item in resolved}))
    }
    ordered = sorted(resolved, key=lambda item: (row_map[item[0]], item[1]))
    rows_by_index = [0] * record_count
    for input_row, _token_index, record_index in resolved:
        rows_by_index[record_index] = row_map[input_row]
    return tuple(item[2] for item in ordered), tuple(rows_by_index)


__all__ = [
    "FeatureAnchorPlacement",
    "FeatureAnchorStatus",
    "parse_record_row_position",
    "resolve_record_row_positions",
]
