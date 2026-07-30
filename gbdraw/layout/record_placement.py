"""Shared record-row placement grammar and selector resolution."""

from __future__ import annotations

from typing import Literal, Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import parse_record_selector


_PlacementCompatibility = Literal["linear", "circular"]


def parse_record_row_position(
    value: str,
    *,
    _compatibility: _PlacementCompatibility = "linear",
) -> tuple[str, int]:
    """Parse one ``<selector>@<row>`` placement token."""

    raw = str(value or "").strip() if _compatibility == "circular" else str(value).strip()
    missing_separator = "@" not in raw
    if not raw:
        raise ValidationError("multi_record_position does not allow empty entries.")
    if missing_separator or (_compatibility == "linear" and raw.count("@") != 1):
        raise ValidationError(
            f"multi_record_position entry '{raw}' must be in '<selector>@<row>' format."
        )
    selector, row_text = (part.strip() for part in raw.rsplit("@", 1))
    if not selector:
        raise ValidationError(
            f"multi_record_position entry '{raw}' must include a selector before '@'."
        )
    try:
        row = int(row_text)
    except ValueError as exc:
        raise ValidationError(
            f"multi_record_position entry '{raw}' must use a positive integer row."
        ) from exc
    if (
        row <= 0
        or (_compatibility == "circular" and not row_text.isdigit())
    ):
        raise ValidationError(
            f"multi_record_position entry '{raw}' must use a positive integer row."
        )
    if _compatibility == "circular":
        try:
            parsed_selector = parse_record_selector(selector)
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc
        if parsed_selector is None:
            raise ValidationError(
                f"multi_record_position selector '{selector}' is invalid."
            )
    return selector, row


def _resolve_record_selector_index(
    records: Sequence[SeqRecord],
    selector: str,
    *,
    compatibility: _PlacementCompatibility,
) -> int:
    record_count = len(records)
    if compatibility == "linear":
        if selector.startswith("#"):
            try:
                record_index = int(selector[1:]) - 1
            except ValueError as exc:
                raise ValidationError(
                    f"multi_record_position selector '{selector}' is invalid."
                ) from exc
            if 0 <= record_index < record_count:
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
    else:
        try:
            parsed = parse_record_selector(selector)
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc
        if parsed is None:
            raise ValidationError(
                f"multi_record_position selector '{selector}' is invalid."
            )
        if parsed.record_index is not None:
            record_index = int(parsed.record_index)
            if 0 <= record_index < record_count:
                return record_index
            raise ValidationError(
                f"multi_record_position selector '{selector}' is out of range for "
                f"{record_count} record(s)."
            )
        matches = [
            index
            for index, record in enumerate(records)
            if str(record.id) == str(parsed.record_id)
        ]
    if not matches:
        raise ValidationError(
            f"multi_record_position selector '{selector}' did not match any record ID."
        )
    if len(matches) > 1:
        suffix = (
            "matched multiple records. Use #index to disambiguate."
            if compatibility == "circular"
            else "matched multiple records; use a #index selector."
        )
        raise ValidationError(
            f"multi_record_position selector '{selector}' {suffix}"
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
        selector, row = parse_record_row_position(
            str(raw),
            _compatibility=_compatibility,
        )
        record_index = _resolve_record_selector_index(
            records,
            selector,
            compatibility=_compatibility,
        )
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


__all__ = ["parse_record_row_position", "resolve_record_row_positions"]
