"""CLI- and session-independent request value objects.

These models describe unresolved or materialized record inputs and render intent.
They intentionally do not parse session documents or expose CLI argument names.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path, PureWindowsPath
from typing import Literal, Sequence, TypeAlias

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import RecordSelector
from gbdraw.io.regions import RegionSpec
from gbdraw.render.formats import ACCEPTED_FORMATS, normalize_format_token
from gbdraw.render.output_paths import is_windows_reserved_filename_component

from .options import (
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    resolve_circular_diagram_options,
    resolve_linear_diagram_options,
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


class RecordCardinality(str, Enum):
    """How many displayed records one input source is allowed to resolve."""

    EXACTLY_ONE = "exactly_one"
    FIRST = "first"
    ALL = "all"


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
    record_key: str | None = None
    cardinality: RecordCardinality = RecordCardinality.EXACTLY_ONE

    def __post_init__(self) -> None:
        if not isinstance(
            self.source,
            (GenBankInputSource, GffFastaInputSource, InMemoryRecordSource),
        ):
            raise ValidationError(
                "Record source must be GenBank, GFF3/FASTA, or an in-memory SeqRecord."
            )
        if not isinstance(self.cardinality, RecordCardinality):
            raise ValidationError(
                "Record cardinality must be a RecordCardinality value."
            )
        if not isinstance(self.presentation, RecordPresentation):
            raise ValidationError("Record presentation has an unsupported type.")
        if self.record_key is not None:
            if not isinstance(self.record_key, str) or not self.record_key.strip():
                raise ValidationError("record_key must be a non-empty string or None.")
            object.__setattr__(self, "record_key", self.record_key.strip())
        if self.selector is not None:
            if not isinstance(self.selector, RecordSelector):
                raise ValidationError("Record selector has an unsupported type.")
            if (self.selector.record_id is None) == (self.selector.record_index is None):
                raise ValidationError("Record selector must identify one record ID or index.")
            if self.selector.record_index is not None and self.selector.record_index < 0:
                raise ValidationError("Record selector index must be non-negative.")
            if self.cardinality is not RecordCardinality.EXACTLY_ONE:
                raise ValidationError(
                    "A record selector identifies one record and therefore requires "
                    "RecordCardinality.EXACTLY_ONE."
                )
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
            if (
                region_selects_record
                and self.cardinality is not RecordCardinality.EXACTLY_ONE
            ):
                raise ValidationError(
                    "A selector-qualified region identifies one record and therefore "
                    "requires RecordCardinality.EXACTLY_ONE."
                )
            if (
                not region_selects_record
                and self.cardinality is RecordCardinality.ALL
            ):
                raise ValidationError(
                    "A selectorless per-input region cannot be combined with "
                    "RecordCardinality.ALL; use request record_options for "
                    "collection-level regions."
                )
            if self.region.reverse_complement and self.presentation.reverse_complement:
                raise ValidationError(
                    "Reverse complement is set in both the region and record presentation."
                )


@dataclass(frozen=True)
class RecordCollectionOptions:
    """Transforms applied after all unresolved inputs have expanded."""

    regions: Sequence[RegionSpec] = ()
    labels: Sequence[str] = ()
    subtitles: Sequence[str] = ()

    def __post_init__(self) -> None:
        for field_name in ("regions", "labels", "subtitles"):
            value = getattr(self, field_name)
            if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
                raise ValidationError(
                    f"Record collection {field_name} must be a sequence."
                )
        regions = tuple(self.regions)
        if not all(isinstance(region, RegionSpec) for region in regions):
            raise ValidationError(
                "Record collection regions must contain RegionSpec values."
            )
        object.__setattr__(self, "regions", regions)
        for field_name in ("labels", "subtitles"):
            raw_values = tuple(getattr(self, field_name))
            if not all(isinstance(value, str) for value in raw_values):
                raise ValidationError(
                    f"Record collection {field_name} must contain strings."
                )
            values = tuple(value.strip() for value in raw_values)
            object.__setattr__(self, field_name, values)


def _normalized_formats(formats: Sequence[str] | str) -> tuple[str, ...]:
    if not isinstance(formats, (str, Sequence)):
        raise ValidationError("formats must be a string or sequence of strings.")
    raw_formats = formats.split(",") if isinstance(formats, str) else formats
    normalized_formats: list[str] = []
    for raw_format in raw_formats:
        normalized = normalize_format_token(raw_format)
        if normalized not in ACCEPTED_FORMATS:
            raise ValidationError(f"Unsupported output format: {raw_format}")
        if normalized not in normalized_formats:
            normalized_formats.append(normalized)
    if not normalized_formats:
        raise ValidationError("At least one output format is required.")
    return tuple(normalized_formats)


def _validate_interactive_metadata_policy(
    policy: object,
    formats: Sequence[str],
) -> None:
    if policy not in {"auto", "required", "omit"}:
        raise ValidationError(
            "interactive_metadata_policy must be 'auto', 'required', or 'omit'."
        )
    if policy == "required" and "interactive_svg" not in formats:
        raise ValidationError(
            "interactive metadata can be required only when interactive_svg is requested."
        )


def _contains_ascii_control(value: str) -> bool:
    return any(ord(character) < 32 or ord(character) == 127 for character in value)


@dataclass(frozen=True)
class RenderOutputRequest:
    """Output destination and format policy for a typed render request."""

    output_prefix: str = "out"
    output_directory: str | Path | None = None
    formats: Sequence[str] | str = ("svg",)
    overwrite: bool = False
    interactive_metadata_policy: Literal["auto", "required", "omit"] = "auto"
    resolve_prefix_from_first_record: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.output_prefix, str):
            raise ValidationError("output_prefix must be a string.")
        raw_prefix = str(self.output_prefix)
        prefix = raw_prefix.strip()
        if (
            not prefix
            or Path(prefix).name != prefix
            or PureWindowsPath(prefix).name != prefix
            or prefix in {".", ".."}
            or _contains_ascii_control(raw_prefix)
            or is_windows_reserved_filename_component(prefix)
        ):
            raise ValidationError(
                "output_prefix must be a valid portable filename prefix "
                "without directories."
            )
        object.__setattr__(self, "output_prefix", prefix)

        if self.output_directory is not None:
            if not isinstance(self.output_directory, (str, Path)):
                raise ValidationError("output_directory must be a path or None.")
            raw_directory = str(self.output_directory)
            directory = raw_directory.strip()
            if not directory:
                raise ValidationError("output_directory must not be empty.")
            if _contains_ascii_control(raw_directory):
                raise ValidationError(
                    "output_directory must not contain control characters."
                )
            object.__setattr__(self, "output_directory", Path(directory))

        normalized_formats = _normalized_formats(self.formats)
        object.__setattr__(self, "formats", normalized_formats)

        if not isinstance(self.overwrite, bool):
            raise ValidationError("overwrite must be a boolean.")
        if not isinstance(self.resolve_prefix_from_first_record, bool):
            raise ValidationError(
                "resolve_prefix_from_first_record must be a boolean."
            )
        _validate_interactive_metadata_policy(
            self.interactive_metadata_policy,
            normalized_formats,
        )


@dataclass(frozen=True)
class CircularBatchOutputPolicy:
    """Resolve one collision-free output name per displayed Circular record."""

    output_prefix: str | Path | None = None
    formats: Sequence[str] | str = ("svg",)
    overwrite: bool = False
    interactive_metadata_policy: Literal["auto", "required", "omit"] = "auto"

    def __post_init__(self) -> None:
        if self.output_prefix is not None:
            if not isinstance(self.output_prefix, (str, Path)):
                raise ValidationError(
                    "Circular batch output_prefix must be a path or None."
                )
            raw_prefix = str(self.output_prefix)
            prefix = raw_prefix.strip()
            if (
                not prefix
                or prefix in {".", ".."}
                or _contains_ascii_control(raw_prefix)
            ):
                raise ValidationError(
                    "Circular batch output_prefix must identify a filename prefix."
                )
            path = Path(prefix)
            if not path.name:
                raise ValidationError(
                    "Circular batch output_prefix must identify a filename prefix."
                )
            object.__setattr__(self, "output_prefix", path)
        normalized_formats = _normalized_formats(self.formats)
        object.__setattr__(self, "formats", normalized_formats)
        if not isinstance(self.overwrite, bool):
            raise ValidationError("overwrite must be a boolean.")
        _validate_interactive_metadata_policy(
            self.interactive_metadata_policy,
            normalized_formats,
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
    keys = [record.record_key for record in normalized if record.record_key is not None]
    if len(set(keys)) != len(keys):
        raise ValidationError("Diagram request record keys must be unique.")
    return normalized


def _validate_circular_placements(
    records: Sequence[RecordInput],
    *,
    layout: CircularMultiRecordOptions | None,
) -> None:
    if any(
        record.cardinality is RecordCardinality.ALL
        for record in records
    ) and any(
        record.presentation.grid_row is not None for record in records
    ):
        raise ValidationError(
            "RecordCardinality.ALL cannot use per-input grid placement; "
            "placement is resolved only after record expansion."
        )
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


def _validate_linear_placements(
    records: Sequence[RecordInput],
    *,
    layout: LinearMultiRecordOptions | None,
) -> None:
    if any(
        record.cardinality is RecordCardinality.ALL
        for record in records
    ) and any(
        record.presentation.grid_row is not None for record in records
    ):
        raise ValidationError(
            "RecordCardinality.ALL cannot use per-input grid placement; "
            "placement is resolved only after record expansion."
        )
    placements = [record.presentation for record in records]
    has_row = any(item.grid_row is not None for item in placements)
    has_column = any(item.grid_column is not None for item in placements)
    if not has_row:
        return
    if layout is None:
        raise ValidationError("Grid placement requires a Linear layout.")
    if any(item.grid_row is None for item in placements):
        raise ValidationError("If one Linear record has grid_row, every record must have it.")
    if has_column and any(item.grid_column is None for item in placements):
        raise ValidationError("If one Linear record has grid_column, every record must have it.")
    occupied = [
        (item.grid_row, item.grid_column)
        for item in placements
        if item.grid_column is not None
    ]
    if len(set(occupied)) != len(occupied):
        raise ValidationError("Linear record grid placements must be unique.")
    if layout.multi_record_positions:
        raise ValidationError(
            "Specify Linear placement in RecordPresentation or layout, not both."
        )


def _resolve_circular_request_options(
    options: object,
) -> CircularDiagramOptions:
    if not isinstance(options, CircularDiagramOptions):
        raise ValidationError(
            "Circular request options must be CircularDiagramOptions."
        )
    return resolve_circular_diagram_options(options)


@dataclass(frozen=True)
class CircularDiagramRequest:
    """Possibly unresolved record inputs for one Circular diagram.

    ``grouping`` is normalized eagerly so the request states whether the one
    diagram is a regular single-record render or a grid.  A grid may contain
    one record; separate-diagram collections use :class:`CircularBatchRequest`.
    Request planning expands cardinality and returns an exact-one materialized
    projection suitable for canonical schema 5.
    """

    records: Sequence[RecordInput]
    options: CircularDiagramOptions = field(default_factory=CircularDiagramOptions)
    layout: CircularMultiRecordOptions | None = None
    output: RenderOutputRequest = field(default_factory=RenderOutputRequest)
    grouping: Literal["single", "grid"] | None = None
    record_options: RecordCollectionOptions = field(
        default_factory=RecordCollectionOptions
    )

    def __post_init__(self) -> None:
        records = _request_records(self.records)
        object.__setattr__(self, "records", records)
        object.__setattr__(
            self,
            "options",
            _resolve_circular_request_options(self.options),
        )
        if self.layout is not None and not isinstance(
            self.layout, CircularMultiRecordOptions
        ):
            raise ValidationError("Circular request layout has an unsupported type.")
        grouping = self.grouping
        if grouping is None:
            grouping = (
                "grid"
                if self.layout is not None
                or len(records) > 1
                or any(
                    record.cardinality is RecordCardinality.ALL
                    for record in records
                )
                else "single"
            )
        if grouping not in {"single", "grid"}:
            raise ValidationError("Circular request grouping must be 'single' or 'grid'.")
        if grouping == "single" and len(records) != 1:
            raise ValidationError(
                "A single Circular request requires exactly one record."
            )
        if (
            grouping == "single"
            and records[0].cardinality is RecordCardinality.ALL
        ):
            raise ValidationError(
                "A single Circular request cannot use RecordCardinality.ALL."
            )
        if grouping == "single" and self.layout is not None:
            raise ValidationError(
                "A single Circular request cannot define a grid layout."
            )
        if grouping == "grid" and self.layout is None:
            object.__setattr__(self, "layout", CircularMultiRecordOptions())
        object.__setattr__(self, "grouping", grouping)
        if not isinstance(self.output, RenderOutputRequest):
            raise ValidationError("Circular request output has an unsupported type.")
        if not isinstance(self.record_options, RecordCollectionOptions):
            raise ValidationError(
                "Circular request record_options has an unsupported type."
            )
        _validate_circular_placements(records, layout=self.layout)


@dataclass(frozen=True)
class CircularBatchRequest:
    """Possibly unresolved inputs for separate one-record Circular diagrams."""

    records: Sequence[RecordInput]
    options: CircularDiagramOptions = field(default_factory=CircularDiagramOptions)
    outputs: Sequence[RenderOutputRequest] = ()
    output_policy: CircularBatchOutputPolicy | None = None
    record_options: RecordCollectionOptions = field(
        default_factory=RecordCollectionOptions
    )

    def __post_init__(self) -> None:
        records = _request_records(self.records)
        object.__setattr__(self, "records", records)
        object.__setattr__(
            self,
            "options",
            _resolve_circular_request_options(self.options),
        )
        if isinstance(self.outputs, (str, bytes)) or not isinstance(
            self.outputs,
            Sequence,
        ):
            raise ValidationError(
                "Circular batch outputs must be a sequence of RenderOutputRequest values."
            )
        outputs = tuple(self.outputs)
        if bool(outputs) == bool(self.output_policy):
            raise ValidationError(
                "Circular batch requires either explicit outputs or one output_policy."
            )
        if outputs and len(outputs) != len(records):
            raise ValidationError(
                "Circular batch requires one resolved output request per "
                "single-record input."
            )
        if not all(isinstance(output, RenderOutputRequest) for output in outputs):
            raise ValidationError(
                "Circular batch outputs must be RenderOutputRequest values."
            )
        resolved_targets = {
            (
                Path(output.output_directory or "."),
                output.output_prefix,
            )
            for output in outputs
        }
        if len(resolved_targets) != len(outputs):
            raise ValidationError(
                "Circular batch output prefixes must be unique within each directory."
            )
        if outputs and any(
            record.cardinality is RecordCardinality.ALL for record in records
        ):
            raise ValidationError(
                "Circular batch inputs using RecordCardinality.ALL require an "
                "output_policy so outputs can be resolved after expansion."
            )
        if self.output_policy is not None and not isinstance(
            self.output_policy,
            CircularBatchOutputPolicy,
        ):
            raise ValidationError(
                "Circular batch output_policy has an unsupported type."
            )
        if not isinstance(self.record_options, RecordCollectionOptions):
            raise ValidationError(
                "Circular batch record_options has an unsupported type."
            )
        if any(record.presentation.grid_row is not None for record in records):
            raise ValidationError(
                "Circular batch records cannot define grid placement."
            )
        object.__setattr__(self, "outputs", outputs)

    @property
    def grouping(self) -> Literal["batch"]:
        return "batch"


@dataclass(frozen=True)
class LinearDiagramRequest:
    """Possibly unresolved record inputs and options for a Linear render."""

    records: Sequence[RecordInput]
    options: LinearDiagramOptions = field(default_factory=LinearDiagramOptions)
    layout: LinearMultiRecordOptions | None = None
    output: RenderOutputRequest = field(default_factory=RenderOutputRequest)
    record_options: RecordCollectionOptions = field(
        default_factory=RecordCollectionOptions
    )

    def __post_init__(self) -> None:
        records = _request_records(self.records)
        object.__setattr__(self, "records", records)
        if not isinstance(self.options, LinearDiagramOptions):
            raise ValidationError(
                "Linear request options must be LinearDiagramOptions."
            )
        object.__setattr__(
            self,
            "options",
            resolve_linear_diagram_options(self.options),
        )
        if self.layout is not None and not isinstance(self.layout, LinearMultiRecordOptions):
            raise ValidationError("Linear request layout has an unsupported type.")
        if not isinstance(self.output, RenderOutputRequest):
            raise ValidationError("Linear request output has an unsupported type.")
        if not isinstance(self.record_options, RecordCollectionOptions):
            raise ValidationError(
                "Linear request record_options has an unsupported type."
            )
        if any(record.presentation.grid_row is not None for record in records) and self.layout is None:
            raise ValidationError(
                "Grid placement requires LinearMultiRecordOptions; without it, grid "
                "placement is supported only by circular multi-record requests."
            )
        _validate_linear_placements(records, layout=self.layout)


DiagramRequest: TypeAlias = (
    CircularDiagramRequest | CircularBatchRequest | LinearDiagramRequest
)


__all__ = [
    "CircularBatchRequest",
    "CircularBatchOutputPolicy",
    "CircularDiagramRequest",
    "DiagramRequest",
    "GenBankInputSource",
    "GffFastaInputSource",
    "InMemoryRecordSource",
    "LinearDiagramRequest",
    "RecordCardinality",
    "RecordCollectionOptions",
    "RecordInput",
    "RecordInputSource",
    "RecordPresentation",
    "RenderOutputRequest",
]
