"""CLI- and session-independent request value objects.

These models describe materialized record inputs and mode-specific render intent.
They intentionally do not parse session documents or expose CLI argument names.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path, PureWindowsPath
from typing import Literal, Sequence, TypeAlias

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import RecordSelector
from gbdraw.io.regions import RegionSpec
from gbdraw.render.formats import ACCEPTED_FORMATS, normalize_format_token

from .options import (
    CircularMultiRecordOptions,
    DiagramOptions,
    _validate_diagram_options_mode,
)


def _materialized_path(value: str | Path, *, field_name: str) -> Path:
    if not isinstance(value, (str, Path)):
        raise ValidationError(f"{field_name} must identify a materialized file.")
    raw = str(value).strip()
    if not raw or raw in {".", ".."}:
        raise ValidationError(f"{field_name} must identify a materialized file.")
    path = Path(raw)
    if not path.name:
        raise ValidationError(f"{field_name} must identify a materialized file.")
    return path


@dataclass(frozen=True)
class GenBankInputSource:
    """One materialized GenBank input file."""

    path: str | Path

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "path",
            _materialized_path(self.path, field_name="GenBank path"),
        )


@dataclass(frozen=True)
class GffFastaInputSource:
    """One materialized GFF3/FASTA input pair."""

    gff_path: str | Path
    fasta_path: str | Path

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "gff_path",
            _materialized_path(self.gff_path, field_name="GFF3 path"),
        )
        object.__setattr__(
            self,
            "fasta_path",
            _materialized_path(self.fasta_path, field_name="FASTA path"),
        )


@dataclass(frozen=True)
class InMemoryRecordSource:
    """One already-parsed sequence record."""

    record: SeqRecord

    def __post_init__(self) -> None:
        if not isinstance(self.record, SeqRecord):
            raise ValidationError("In-memory record source must contain a SeqRecord.")


RecordInputSource: TypeAlias = (
    GenBankInputSource | GffFastaInputSource | InMemoryRecordSource
)


@dataclass(frozen=True)
class RecordPresentation:
    """Per-record presentation and placement values."""

    label: str | None = None
    subtitle: str | None = None
    reverse_complement: bool = False
    grid_row: int | None = None
    grid_column: int | None = None

    def __post_init__(self) -> None:
        for field_name in ("label", "subtitle"):
            value = getattr(self, field_name)
            if value is not None and not isinstance(value, str):
                raise ValidationError(f"Record {field_name} must be a string or None.")
            if isinstance(value, str):
                object.__setattr__(self, field_name, value.strip() or None)
        if not isinstance(self.reverse_complement, bool):
            raise ValidationError("reverse_complement must be a boolean.")
        for field_name in ("grid_row", "grid_column"):
            value = getattr(self, field_name)
            if value is not None and (
                not isinstance(value, int) or isinstance(value, bool) or value < 1
            ):
                raise ValidationError(f"{field_name} must be a positive integer or None.")
        if self.grid_column is not None and self.grid_row is None:
            raise ValidationError("grid_column requires grid_row.")


@dataclass(frozen=True)
class RecordInput:
    """A record source plus selection, region, and presentation metadata."""

    source: RecordInputSource
    selector: RecordSelector | None = None
    region: RegionSpec | None = None
    presentation: RecordPresentation = field(default_factory=RecordPresentation)

    def __post_init__(self) -> None:
        if not isinstance(
            self.source,
            (GenBankInputSource, GffFastaInputSource, InMemoryRecordSource),
        ):
            raise ValidationError(
                "Record source must be GenBank, GFF3/FASTA, or an in-memory SeqRecord."
            )
        if not isinstance(self.presentation, RecordPresentation):
            raise ValidationError("Record presentation has an unsupported type.")
        if self.selector is not None:
            if not isinstance(self.selector, RecordSelector):
                raise ValidationError("Record selector has an unsupported type.")
            if (self.selector.record_id is None) == (self.selector.record_index is None):
                raise ValidationError("Record selector must identify one record ID or index.")
            if self.selector.record_index is not None and self.selector.record_index < 0:
                raise ValidationError("Record selector index must be non-negative.")
        if self.region is not None:
            if not isinstance(self.region, RegionSpec):
                raise ValidationError("Record region has an unsupported type.")
            if self.region.file_selector is not None:
                raise ValidationError(
                    "A materialized RecordInput region must not contain a file selector."
                )
            if self.region.start < 1 or self.region.end < self.region.start:
                raise ValidationError("Record region coordinates are invalid.")
            if self.region.record_index is not None and self.region.record_index < 0:
                raise ValidationError("Record region selector index must be non-negative.")
            region_selects_record = (
                self.region.record_id is not None or self.region.record_index is not None
            )
            if self.selector is not None and region_selects_record:
                raise ValidationError(
                    "Specify the record selector either on RecordInput or in its region, not both."
                )
            if self.region.reverse_complement and self.presentation.reverse_complement:
                raise ValidationError(
                    "Reverse complement is set in both the region and record presentation."
                )


@dataclass(frozen=True)
class RenderOutputRequest:
    """Output destination and format policy for a typed render request."""

    output_prefix: str = "out"
    output_directory: str | Path | None = None
    formats: Sequence[str] | str = ("svg",)
    overwrite: bool = False
    interactive_metadata_policy: Literal["auto", "required", "omit"] = "auto"

    def __post_init__(self) -> None:
        if not isinstance(self.output_prefix, str):
            raise ValidationError("output_prefix must be a string.")
        prefix = str(self.output_prefix).strip()
        if (
            not prefix
            or Path(prefix).name != prefix
            or PureWindowsPath(prefix).name != prefix
            or prefix in {".", ".."}
        ):
            raise ValidationError(
                "output_prefix must be a non-empty filename prefix without directories."
            )
        object.__setattr__(self, "output_prefix", prefix)

        if self.output_directory is not None:
            if not isinstance(self.output_directory, (str, Path)):
                raise ValidationError("output_directory must be a path or None.")
            raw_directory = str(self.output_directory).strip()
            if not raw_directory:
                raise ValidationError("output_directory must not be empty.")
            object.__setattr__(self, "output_directory", Path(raw_directory))

        if not isinstance(self.formats, (str, Sequence)):
            raise ValidationError("formats must be a string or sequence of strings.")
        raw_formats = self.formats.split(",") if isinstance(self.formats, str) else self.formats
        normalized_formats: list[str] = []
        for raw_format in raw_formats:
            normalized = normalize_format_token(raw_format)
            if normalized not in ACCEPTED_FORMATS:
                raise ValidationError(f"Unsupported output format: {raw_format}")
            if normalized not in normalized_formats:
                normalized_formats.append(normalized)
        if not normalized_formats:
            raise ValidationError("At least one output format is required.")
        object.__setattr__(self, "formats", tuple(normalized_formats))

        if not isinstance(self.overwrite, bool):
            raise ValidationError("overwrite must be a boolean.")
        if self.interactive_metadata_policy not in {"auto", "required", "omit"}:
            raise ValidationError(
                "interactive_metadata_policy must be 'auto', 'required', or 'omit'."
            )
        if (
            self.interactive_metadata_policy == "required"
            and "interactive_svg" not in normalized_formats
        ):
            raise ValidationError(
                "interactive metadata can be required only when interactive_svg is requested."
            )


def _request_records(records: Sequence[RecordInput]) -> tuple[RecordInput, ...]:
    try:
        normalized = tuple(records)
    except TypeError as exc:
        raise ValidationError("A diagram request requires a record input sequence.") from exc
    if not normalized:
        raise ValidationError("A diagram request requires at least one record input.")
    if not all(isinstance(record, RecordInput) for record in normalized):
        raise ValidationError("Diagram request records must be RecordInput values.")
    return normalized


def _validate_circular_placements(
    records: Sequence[RecordInput],
    *,
    layout: CircularMultiRecordOptions | None,
) -> None:
    placements = [record.presentation for record in records]
    has_row = any(item.grid_row is not None for item in placements)
    has_column = any(item.grid_column is not None for item in placements)
    if not has_row:
        return
    if layout is None:
        raise ValidationError("Grid placement requires a circular multi-record request.")
    if any(item.grid_row is None for item in placements):
        raise ValidationError("If one circular record has grid_row, every record must have it.")
    if has_column and any(item.grid_column is None for item in placements):
        raise ValidationError("If one circular record has grid_column, every record must have it.")
    occupied = [
        (item.grid_row, item.grid_column)
        for item in placements
        if item.grid_column is not None
    ]
    if len(set(occupied)) != len(occupied):
        raise ValidationError("Circular record grid placements must be unique.")
    if layout.multi_record_positions:
        raise ValidationError(
            "Specify circular placement in RecordPresentation or layout, not both."
        )


@dataclass(frozen=True)
class CircularDiagramRequest:
    """Materialized record inputs and options for a Circular render."""

    records: Sequence[RecordInput]
    options: DiagramOptions = field(default_factory=DiagramOptions)
    layout: CircularMultiRecordOptions | None = None
    output: RenderOutputRequest = field(default_factory=RenderOutputRequest)

    def __post_init__(self) -> None:
        records = _request_records(self.records)
        object.__setattr__(self, "records", records)
        if not isinstance(self.options, DiagramOptions):
            raise ValidationError("Circular request options must be DiagramOptions.")
        if self.layout is not None and not isinstance(
            self.layout, CircularMultiRecordOptions
        ):
            raise ValidationError("Circular request layout has an unsupported type.")
        if len(records) > 1 and self.layout is None:
            object.__setattr__(self, "layout", CircularMultiRecordOptions())
        if not isinstance(self.output, RenderOutputRequest):
            raise ValidationError("Circular request output has an unsupported type.")
        _validate_diagram_options_mode(
            self.options,
            mode="circular_multi" if self.layout is not None else "circular",
        )
        _validate_circular_placements(records, layout=self.layout)


@dataclass(frozen=True)
class LinearDiagramRequest:
    """Materialized record inputs and options for a Linear render."""

    records: Sequence[RecordInput]
    options: DiagramOptions = field(default_factory=DiagramOptions)
    output: RenderOutputRequest = field(default_factory=RenderOutputRequest)

    def __post_init__(self) -> None:
        records = _request_records(self.records)
        object.__setattr__(self, "records", records)
        if not isinstance(self.options, DiagramOptions):
            raise ValidationError("Linear request options must be DiagramOptions.")
        if not isinstance(self.output, RenderOutputRequest):
            raise ValidationError("Linear request output has an unsupported type.")
        if any(record.presentation.grid_row is not None for record in records):
            raise ValidationError("Grid placement is supported only by circular multi-record requests.")
        _validate_diagram_options_mode(self.options, mode="linear")


DiagramRequest: TypeAlias = CircularDiagramRequest | LinearDiagramRequest


__all__ = [
    "CircularDiagramRequest",
    "DiagramRequest",
    "GenBankInputSource",
    "GffFastaInputSource",
    "InMemoryRecordSource",
    "LinearDiagramRequest",
    "RecordInput",
    "RecordInputSource",
    "RecordPresentation",
    "RenderOutputRequest",
]
