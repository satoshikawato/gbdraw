"""Shared CLI helpers for headered TSV table arguments."""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path, PurePath, PureWindowsPath
from typing import Mapping, Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError


_NONE_LIKE = {"", "auto", "none", "null", "-"}
_DISPLAY_SELECTOR_PREFIXES = {"input", "record", "name", "file", "accession"}
_SOURCE_SELECTOR_PREFIXES = {"record", "name", "accession"}


@dataclass(frozen=True)
class HeaderedTableRow:
    """One logical data row from a headered CLI TSV file."""

    row_number: int
    cells: Mapping[str, str]

    def cell(self, column: str) -> str:
        return str(self.cells.get(column, "")).strip()


def table_error(
    table_name: str,
    message: str,
    *,
    row_number: int | None = None,
    column: str | None = None,
) -> ValidationError:
    parts = [str(table_name)]
    if row_number is not None:
        parts.append(f"row {row_number}")
    if column:
        parts.append(f"column '{column}'")
    return ValidationError(f"{', '.join(parts)}: {message}")


def _is_comment_line(raw_line: str) -> bool:
    stripped = raw_line.lstrip()
    return (
        stripped == "#"
        or stripped.startswith("# ")
        or stripped.startswith("#\t")
        or stripped.startswith("##")
    )


def _parse_tsv_line(raw_line: str) -> list[str]:
    return next(csv.reader([raw_line.rstrip("\r\n")], delimiter="\t"))


def read_headered_tsv_table(
    path: str,
    required: Sequence[str],
    optional: Sequence[str],
    table_name: str,
) -> list[HeaderedTableRow]:
    """Read a headered TSV table and reject unknown columns."""

    table_path = Path(path)
    if not table_path.is_file():
        raise ValidationError(f"{table_name}: file does not exist or is not accessible: {path}")

    logical_lines: list[tuple[int, str]] = []
    with table_path.open("r", encoding="utf-8", newline="") as handle:
        for row_number, raw_line in enumerate(handle, start=1):
            if not raw_line.strip() or _is_comment_line(raw_line):
                continue
            logical_lines.append((row_number, raw_line))

    if not logical_lines:
        raise ValidationError(f"{table_name}: table is empty.")

    header_row_number, header_line = logical_lines[0]
    headers = [cell.strip() for cell in _parse_tsv_line(header_line)]
    if not headers or any(not header for header in headers):
        raise table_error(table_name, "header contains an empty column name", row_number=header_row_number)
    duplicated = sorted({header for header in headers if headers.count(header) > 1})
    if duplicated:
        raise table_error(
            table_name,
            f"duplicate column(s): {', '.join(duplicated)}",
            row_number=header_row_number,
        )

    required_set = set(required)
    optional_set = set(optional)
    allowed = required_set | optional_set
    missing = sorted(required_set.difference(headers))
    if missing:
        raise table_error(
            table_name,
            f"missing required column(s): {', '.join(missing)}",
            row_number=header_row_number,
        )
    unknown = sorted(set(headers).difference(allowed))
    if unknown:
        raise table_error(
            table_name,
            f"unknown column(s): {', '.join(unknown)}",
            row_number=header_row_number,
        )

    rows: list[HeaderedTableRow] = []
    for row_number, raw_line in logical_lines[1:]:
        values = _parse_tsv_line(raw_line)
        if len(values) > len(headers):
            raise table_error(
                table_name,
                f"expected at most {len(headers)} cell(s), got {len(values)}",
                row_number=row_number,
            )
        padded = values + [""] * (len(headers) - len(values))
        rows.append(
            HeaderedTableRow(
                row_number=row_number,
                cells={header: str(value).strip() for header, value in zip(headers, padded)},
            )
        )
    return rows


def resolve_table_path(table_path: str, value: str) -> str:
    """Resolve a table cell path relative to the table file directory."""

    text = str(value or "").strip()
    if not text:
        return ""
    path = Path(text)
    if path.is_absolute():
        return str(path)
    return str(Path(table_path).resolve().parent / path)


def parse_optional_float_cell(
    value: object,
    field_name: str,
    *,
    positive: bool = False,
    nonnegative: bool = False,
    table_name: str = "table",
    row_number: int | None = None,
    column: str | None = None,
) -> float | None:
    text = str(value or "").strip()
    if text.lower() in _NONE_LIKE:
        return None
    try:
        parsed = float(text)
    except (TypeError, ValueError) as exc:
        raise table_error(
            table_name,
            f"{field_name} must be a number",
            row_number=row_number,
            column=column or field_name,
        ) from exc
    if not math.isfinite(parsed):
        raise table_error(
            table_name,
            f"{field_name} must be finite",
            row_number=row_number,
            column=column or field_name,
        )
    if positive and parsed <= 0:
        raise table_error(
            table_name,
            f"{field_name} must be > 0",
            row_number=row_number,
            column=column or field_name,
        )
    if nonnegative and parsed < 0:
        raise table_error(
            table_name,
            f"{field_name} must be >= 0",
            row_number=row_number,
            column=column or field_name,
        )
    return parsed


def parse_optional_int_cell(
    value: object,
    field_name: str,
    *,
    positive: bool = False,
    nonnegative: bool = False,
    table_name: str = "table",
    row_number: int | None = None,
    column: str | None = None,
) -> int | None:
    text = str(value or "").strip()
    if text.lower() in _NONE_LIKE:
        return None
    try:
        parsed = int(text)
    except (TypeError, ValueError) as exc:
        raise table_error(
            table_name,
            f"{field_name} must be an integer",
            row_number=row_number,
            column=column or field_name,
        ) from exc
    if positive and parsed <= 0:
        raise table_error(
            table_name,
            f"{field_name} must be > 0",
            row_number=row_number,
            column=column or field_name,
        )
    if nonnegative and parsed < 0:
        raise table_error(
            table_name,
            f"{field_name} must be >= 0",
            row_number=row_number,
            column=column or field_name,
        )
    return parsed


def parse_optional_bool_cell(
    value: object,
    field_name: str,
    *,
    table_name: str = "table",
    row_number: int | None = None,
    column: str | None = None,
) -> bool | None:
    text = str(value or "").strip().lower()
    if text == "":
        return None
    if text in {"1", "true", "yes", "on"}:
        return True
    if text in {"0", "false", "no", "off"}:
        return False
    raise table_error(
        table_name,
        f"{field_name} must be one of 1, 0, true, false, yes, no, on, off",
        row_number=row_number,
        column=column or field_name,
    )


def _add_alias(
    index: dict[str, dict[str, set[int]]],
    namespace: str,
    alias: object,
    record_index: int,
) -> None:
    text = str(alias or "").strip()
    if not text:
        return
    index.setdefault(namespace, {}).setdefault(text, set()).add(record_index)


def _source_file_aliases(record: SeqRecord) -> set[str]:
    annotations = getattr(record, "annotations", None) or {}
    aliases: set[str] = set()
    for key in ("gbdraw_source_file", "gbdraw_source_basename"):
        raw = annotations.get(key)
        if raw:
            aliases.add(str(raw))
    expanded: set[str] = set()
    for alias in aliases:
        expanded.add(alias)
        try:
            expanded.add(PurePath(alias).name)
        except Exception:
            pass
        try:
            expanded.add(PureWindowsPath(alias).name)
        except Exception:
            pass
    return {alias for alias in expanded if alias}


def _accession_aliases(record: SeqRecord) -> set[str]:
    aliases = {str(getattr(record, "id", "") or ""), str(getattr(record, "name", "") or "")}
    annotations = getattr(record, "annotations", None) or {}
    for key in ("accession", "accessions", "sequence_version", "version"):
        raw = annotations.get(key)
        if isinstance(raw, (list, tuple, set)):
            aliases.update(str(item) for item in raw)
        elif raw:
            aliases.add(str(raw))
    return {alias.strip() for alias in aliases if str(alias).strip()}


def _explicit_selector_parts(selector: str, prefixes: set[str]) -> tuple[str | None, str]:
    text = str(selector or "").strip()
    if ":" not in text:
        return None, text
    prefix, value = text.split(":", 1)
    normalized_prefix = prefix.strip().lower()
    if normalized_prefix in prefixes:
        return normalized_prefix, value.strip()
    return None, text


def _format_selector_matches(matches: Mapping[int, set[str]]) -> str:
    parts = []
    for record_index, labels in sorted(matches.items()):
        parts.append(f"#{record_index + 1} ({', '.join(sorted(labels))})")
    return "; ".join(parts)


class DisplayRecordContext:
    """Resolve selectors against records that will be displayed in the diagram."""

    def __init__(self, records: Sequence[SeqRecord], *, table_name: str = "table") -> None:
        self.records = list(records)
        self.table_name = table_name
        self._index: dict[str, dict[str, set[int]]] = {}
        for idx, record in enumerate(self.records):
            annotations = getattr(record, "annotations", None) or {}
            _add_alias(self._index, "input", annotations.get("gbdraw_input_id"), idx)
            _add_alias(self._index, "record", getattr(record, "id", ""), idx)
            _add_alias(self._index, "name", getattr(record, "name", ""), idx)
            for alias in _source_file_aliases(record):
                _add_alias(self._index, "file", alias, idx)
            for alias in _accession_aliases(record):
                _add_alias(self._index, "accession", alias, idx)

    def resolve_display_record_selector(
        self,
        selector: str,
        *,
        row_number: int | None = None,
        column: str = "record_id",
    ) -> int:
        text = str(selector or "").strip()
        if not text:
            raise table_error(
                self.table_name,
                "record selector is empty",
                row_number=row_number,
                column=column,
            )
        if text.startswith("#"):
            idx_text = text[1:].strip()
            if not idx_text.isdigit() or int(idx_text) <= 0:
                raise table_error(
                    self.table_name,
                    f"invalid displayed-record index selector '{text}'",
                    row_number=row_number,
                    column=column,
                )
            idx = int(idx_text) - 1
            if idx >= len(self.records):
                raise table_error(
                    self.table_name,
                    f"displayed-record selector '{text}' is out of range; loaded {len(self.records)} record(s)",
                    row_number=row_number,
                    column=column,
                )
            return idx

        namespace, value = _explicit_selector_parts(text, _DISPLAY_SELECTOR_PREFIXES)
        matches: dict[int, set[str]] = {}
        if namespace:
            namespace_matches = self._index.get(namespace, {}).get(value, set())
            for idx in namespace_matches:
                matches.setdefault(idx, set()).add(f"{namespace}:{value}")
        else:
            for ns, alias_map in self._index.items():
                for idx in alias_map.get(value, set()):
                    matches.setdefault(idx, set()).add(f"{ns}:{value}")

        if not matches:
            raise table_error(
                self.table_name,
                f"record selector '{text}' did not match any displayed record",
                row_number=row_number,
                column=column,
            )
        if len(matches) > 1:
            detail = _format_selector_matches(matches)
            raise table_error(
                self.table_name,
                f"record selector '{text}' is ambiguous; matches {detail}. Use an explicit namespace such as input: or record:",
                row_number=row_number,
                column=column,
            )
        return next(iter(matches))


class SourceRecordContext:
    """Resolve selectors against records loaded from one input-table source row."""

    def __init__(self, records: Sequence[SeqRecord], *, table_name: str = "input_table") -> None:
        self.records = list(records)
        self.table_name = table_name
        self._index: dict[str, dict[str, set[int]]] = {}
        for idx, record in enumerate(self.records):
            _add_alias(self._index, "record", getattr(record, "id", ""), idx)
            _add_alias(self._index, "name", getattr(record, "name", ""), idx)
            for alias in _accession_aliases(record):
                _add_alias(self._index, "accession", alias, idx)

    def resolve_source_record_selector(
        self,
        selector: str,
        *,
        row_number: int | None = None,
        column: str = "record_id",
    ) -> int:
        text = str(selector or "").strip()
        if not text:
            raise table_error(
                self.table_name,
                "source-record selector is empty",
                row_number=row_number,
                column=column,
            )
        if text.startswith("#"):
            idx_text = text[1:].strip()
            if not idx_text.isdigit() or int(idx_text) <= 0:
                raise table_error(
                    self.table_name,
                    f"invalid row-local record index selector '{text}'",
                    row_number=row_number,
                    column=column,
                )
            idx = int(idx_text) - 1
            if idx >= len(self.records):
                raise table_error(
                    self.table_name,
                    f"row-local selector '{text}' is out of range; loaded {len(self.records)} record(s) from this source",
                    row_number=row_number,
                    column=column,
                )
            return idx

        namespace, value = _explicit_selector_parts(text, _SOURCE_SELECTOR_PREFIXES)
        if namespace == "input":
            raise table_error(
                self.table_name,
                "input: selectors are not valid before displayed records exist",
                row_number=row_number,
                column=column,
            )
        matches: dict[int, set[str]] = {}
        if namespace:
            namespace_matches = self._index.get(namespace, {}).get(value, set())
            for idx in namespace_matches:
                matches.setdefault(idx, set()).add(f"{namespace}:{value}")
        else:
            for ns, alias_map in self._index.items():
                for idx in alias_map.get(value, set()):
                    matches.setdefault(idx, set()).add(f"{ns}:{value}")

        if not matches:
            raise table_error(
                self.table_name,
                f"source-record selector '{text}' did not match any record in this row's source file",
                row_number=row_number,
                column=column,
            )
        if len(matches) > 1:
            detail = _format_selector_matches(matches)
            raise table_error(
                self.table_name,
                f"source-record selector '{text}' is ambiguous; matches {detail}. Use #index or an explicit namespace",
                row_number=row_number,
                column=column,
            )
        return next(iter(matches))


__all__ = [
    "DisplayRecordContext",
    "HeaderedTableRow",
    "SourceRecordContext",
    "parse_optional_bool_cell",
    "parse_optional_float_cell",
    "parse_optional_int_cell",
    "read_headered_tsv_table",
    "resolve_table_path",
    "table_error",
]
