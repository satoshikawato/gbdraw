"""Resolve typed record inputs once while retaining their source provenance."""

from __future__ import annotations

import copy
import logging
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Callable, Literal, Sequence

import pandas as pd
from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError
from gbdraw.io.cli_tables import (
    read_comparisons_table,
    read_conservation_table,
    read_records_table,
)
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.io.genome import load_gbks, load_gff_fasta
from gbdraw.io.record_select import (
    RecordSelector,
    parse_record_selector,
    reverse_records,
)
from gbdraw.io.regions import RegionSpec, apply_region_specs, parse_region_spec
from gbdraw.layout.record_placement import resolve_record_row_positions
from gbdraw.linear_comparison import LinearComparison

from .options import (
    CircularDiagramOptions,
    DepthTrackInput,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
)
from .requests import (
    CircularBatchOutputPolicy,
    GenBankInputSource,
    GffFastaInputSource,
    InMemoryRecordSource,
    RecordCardinality,
    RecordCollectionOptions,
    RecordInput,
    RecordInputSource,
    RecordPresentation,
    RenderOutputRequest,
)

logger = logging.getLogger(__name__)

GenBankLoader = Callable[..., list[SeqRecord]]
GffFastaLoader = Callable[..., list[SeqRecord]]


@dataclass(frozen=True)
class ResolvedRecordProvenance:
    """Stable source coordinates for one displayed record."""

    resolved_index: int
    input_index: int
    source_record_index: int
    source_record_count: int
    source_record_id: str
    source_kind: Literal["genbank", "gff_fasta", "memory"]
    source_paths: tuple[str, ...]
    record_key: str
    cardinality: RecordCardinality
    selector: RecordSelector | None
    region: RegionSpec | None
    presentation: RecordPresentation


@dataclass(frozen=True)
class ResolvedRecordCollection:
    """Flattened records and one provenance entry per displayed record."""

    records: tuple[SeqRecord, ...]
    provenance: tuple[ResolvedRecordProvenance, ...]

    def __post_init__(self) -> None:
        if not self.records:
            raise ValidationError("A resolved record collection cannot be empty.")
        if len(self.records) != len(self.provenance):
            raise ValidationError(
                "Resolved record provenance must align with displayed records."
            )


@dataclass(frozen=True)
class RecordInputManifest:
    """Unresolved typed record inputs derived from CLI path syntax."""

    records: tuple[RecordInput, ...]
    record_options: RecordCollectionOptions
    source_paths: tuple[str, ...]
    multi_record_positions: tuple[str, ...] = ()


def _expanded_cli_track_values(
    values: Sequence[object] | None,
    *,
    track_count: int,
    field_name: str,
) -> list[object | None]:
    items = list(values or ())
    if not items:
        return [None] * track_count
    if len(items) == 1:
        return items * track_count
    if len(items) != track_count:
        raise ValidationError(
            f"{field_name} count must be one or equal to the number of "
            f"depth tracks ({track_count})."
        )
    return items


def _optional_positive_float(value: object | None, *, field_name: str) -> float | None:
    if value is None or str(value).strip().lower() in {
        "",
        "auto",
        "none",
        "null",
        "-",
    }:
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{field_name} values must be numeric or auto.") from exc
    if numeric <= 0:
        raise ValidationError(f"{field_name} values must be > 0.")
    return numeric


def depth_track_inputs_from_cli(
    groups: Sequence[Sequence[str]] | None,
    *,
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
    heights: Sequence[object] | None = None,
    large_tick_intervals: Sequence[object] | None = None,
    small_tick_intervals: Sequence[object] | None = None,
    tick_font_sizes: Sequence[object] | None = None,
) -> tuple[DepthTrackInput, ...] | None:
    """Translate logical CLI depth groups without knowing record cardinality."""

    if not groups:
        return None
    track_count = len(groups)
    label_values = _expanded_cli_track_values(
        labels,
        track_count=track_count,
        field_name="depth_track_labels",
    )
    color_values = _expanded_cli_track_values(
        colors,
        track_count=track_count,
        field_name="depth_track_colors",
    )
    height_values = _expanded_cli_track_values(
        heights,
        track_count=track_count,
        field_name="depth_track_heights",
    )
    large_values = _expanded_cli_track_values(
        large_tick_intervals,
        track_count=track_count,
        field_name="depth_track_large_tick_intervals",
    )
    small_values = _expanded_cli_track_values(
        small_tick_intervals,
        track_count=track_count,
        field_name="depth_track_small_tick_intervals",
    )
    font_values = _expanded_cli_track_values(
        tick_font_sizes,
        track_count=track_count,
        field_name="depth_track_tick_font_sizes",
    )
    inputs: list[DepthTrackInput] = []
    for index, group in enumerate(groups):
        sources = tuple(
            None
            if str(path).strip().lower() in {"", "-", "none", "null"}
            else str(path)
            for path in group
        )
        if not sources or not any(source is not None for source in sources):
            raise ValidationError(
                f"--depth_track #{index + 1} must include at least one file."
            )
        source: object = sources[0] if len(sources) == 1 else sources
        inputs.append(
            DepthTrackInput(
                source=source,
                label=(
                    str(label_values[index]).strip()
                    if label_values[index] is not None
                    else None
                ),
                color=(
                    str(color_values[index]).strip()
                    if color_values[index] is not None
                    else None
                ),
                height=_optional_positive_float(
                    height_values[index],
                    field_name="depth_track_heights",
                ),
                large_tick_interval=_optional_positive_float(
                    large_values[index],
                    field_name="depth_track_large_tick_intervals",
                ),
                small_tick_interval=_optional_positive_float(
                    small_values[index],
                    field_name="depth_track_small_tick_intervals",
                ),
                tick_font_size=_optional_positive_float(
                    font_values[index],
                    field_name="depth_track_tick_font_sizes",
                ),
            )
        )
    return tuple(inputs)


def _source_key(source: RecordInputSource) -> tuple[object, ...]:
    if isinstance(source, GenBankInputSource):
        return ("genbank", str(source.path))
    if isinstance(source, GffFastaInputSource):
        return ("gff_fasta", str(source.gff_path), str(source.fasta_path))
    if isinstance(source, InMemoryRecordSource):
        return ("memory", id(source.record))
    raise ValidationError("Unsupported record input source.")


def _source_details(
    source: RecordInputSource,
) -> tuple[Literal["genbank", "gff_fasta", "memory"], tuple[str, ...]]:
    if isinstance(source, GenBankInputSource):
        return "genbank", (str(source.path),)
    if isinstance(source, GffFastaInputSource):
        return "gff_fasta", (str(source.gff_path), str(source.fasta_path))
    if isinstance(source, InMemoryRecordSource):
        return "memory", ()
    raise ValidationError("Unsupported record input source.")


def _load_source_records(
    source: RecordInputSource,
    *,
    gff_candidate_features: Sequence[str] | None,
    gff_keep_all_features: bool,
    genbank_loader: GenBankLoader,
    gff_loader: GffFastaLoader,
) -> tuple[SeqRecord, ...]:
    if isinstance(source, GenBankInputSource):
        records = genbank_loader([str(source.path)])
    elif isinstance(source, GffFastaInputSource):
        records = gff_loader(
            [str(source.gff_path)],
            [str(source.fasta_path)],
            selected_features_set=gff_candidate_features,
            keep_all_features=gff_keep_all_features,
        )
    elif isinstance(source, InMemoryRecordSource):
        records = [source.record]
    else:  # pragma: no cover - RecordInput validates this union.
        raise ValidationError("Unsupported record input source.")
    if not records:
        raise ValidationError("A record input source resolved to no records.")
    return tuple(records)


def _selector_from_region(region: RegionSpec | None) -> RecordSelector | None:
    if region is None:
        return None
    if region.record_id is not None:
        return RecordSelector(
            raw=region.record_id,
            record_id=region.record_id,
            record_index=None,
        )
    if region.record_index is not None:
        return RecordSelector(
            raw=f"#{region.record_index + 1}",
            record_id=None,
            record_index=region.record_index,
        )
    return None


def _selected_source_indexes(
    records: Sequence[SeqRecord],
    selector: RecordSelector | None,
) -> list[int]:
    if selector is None:
        return list(range(len(records)))
    if selector.record_index is not None:
        index = selector.record_index
        if index < 0 or index >= len(records):
            raise ValidationError(
                f"Record selector {selector.label()} is out of range "
                f"(loaded {len(records)} record(s))."
            )
        return [index]
    record_id = selector.record_id or ""
    matches = [
        index for index, record in enumerate(records) if record.id == record_id
    ]
    if not matches:
        raise ValidationError(
            f"Record selector '{record_id}' did not match any record ID."
        )
    if len(matches) > 1:
        raise ValidationError(
            f"Record selector '{record_id}' matched multiple records. "
            "Use #index to disambiguate."
        )
    return matches


def _cardinality_indexes(
    indexes: Sequence[int],
    *,
    cardinality: RecordCardinality,
    input_index: int,
) -> tuple[int, ...]:
    selected = tuple(indexes)
    if cardinality is RecordCardinality.EXACTLY_ONE:
        if len(selected) != 1:
            raise ValidationError(
                f"RecordInput #{input_index + 1} requires exactly one record; "
                f"resolved {len(selected)}. Add a selector, use "
                "RecordCardinality.FIRST, or explicitly use RecordCardinality.ALL."
            )
        return selected
    if not selected:
        raise ValidationError(
            f"RecordInput #{input_index + 1} resolved no records."
        )
    if cardinality is RecordCardinality.FIRST:
        return selected[:1]
    if cardinality is RecordCardinality.ALL:
        return selected
    raise ValidationError("Unsupported record cardinality.")


def _unqualified_region(region: RegionSpec) -> RegionSpec:
    return replace(
        region,
        file_selector=None,
        record_id=None,
        record_index=None,
    )


def _apply_presentation(
    record: SeqRecord,
    presentation: RecordPresentation,
    *,
    record_key: str,
) -> None:
    if getattr(record, "annotations", None) is None:
        record.annotations = {}
    if presentation.label:
        record.annotations["gbdraw_record_label"] = presentation.label
    if presentation.subtitle:
        record.annotations["gbdraw_record_subtitle"] = presentation.subtitle
    record.annotations["gbdraw_record_key"] = record_key


def _apply_provenance_annotations(
    record: SeqRecord,
    provenance: ResolvedRecordProvenance,
) -> None:
    if getattr(record, "annotations", None) is None:
        record.annotations = {}
    record.annotations.update(
        {
            "gbdraw_input_index": provenance.input_index,
            "gbdraw_source_record_index": provenance.source_record_index,
            "gbdraw_source_record_count": provenance.source_record_count,
            "gbdraw_source_record_id": provenance.source_record_id,
            "gbdraw_record_key": provenance.record_key,
            "gbdraw_record_cardinality": provenance.cardinality.value,
            "gbdraw_source_kind": provenance.source_kind,
            "gbdraw_source_paths": provenance.source_paths,
        }
    )
    if provenance.source_paths:
        record.annotations["gbdraw_source_file"] = provenance.source_paths[0]
        record.annotations["gbdraw_source_basename"] = Path(
            provenance.source_paths[0]
        ).name


def _apply_collection_options(
    records: list[SeqRecord],
    provenance: Sequence[ResolvedRecordProvenance],
    options: RecordCollectionOptions,
) -> list[SeqRecord]:
    if options.regions:
        try:
            records = apply_region_specs(records, options.regions, log=logger)
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc
    for field_name, annotation_key in (
        ("labels", "gbdraw_record_label"),
        ("subtitles", "gbdraw_record_subtitle"),
    ):
        values = getattr(options, field_name)
        if len(values) > len(records):
            logger.warning(
                "WARNING: More record %s were provided than records resolved; "
                "extra values will be ignored.",
                field_name,
            )
        for index, value in enumerate(values[: len(records)]):
            if value:
                records[index].annotations[annotation_key] = value
    for record, item in zip(records, provenance, strict=True):
        _apply_provenance_annotations(record, item)
    return records


def resolve_record_inputs(
    record_inputs: Sequence[RecordInput],
    *,
    record_options: RecordCollectionOptions | None = None,
    gff_candidate_features: Sequence[str] | None,
    gff_keep_all_features: bool,
    genbank_loader: GenBankLoader = load_gbks,
    gff_loader: GffFastaLoader = load_gff_fasta,
) -> ResolvedRecordCollection:
    """Load each unique source once, then apply typed selection and transforms."""

    inputs = tuple(record_inputs)
    if not inputs:
        raise ValidationError("A request requires at least one RecordInput.")
    cache: dict[tuple[object, ...], tuple[SeqRecord, ...]] = {}
    records: list[SeqRecord] = []
    provenance: list[ResolvedRecordProvenance] = []
    for input_index, record_input in enumerate(inputs):
        key = _source_key(record_input.source)
        raw_records = cache.get(key)
        if raw_records is None:
            raw_records = _load_source_records(
                record_input.source,
                gff_candidate_features=gff_candidate_features,
                gff_keep_all_features=gff_keep_all_features,
                genbank_loader=genbank_loader,
                gff_loader=gff_loader,
            )
            cache[key] = raw_records
        selector = record_input.selector or _selector_from_region(record_input.region)
        source_indexes = _cardinality_indexes(
            _selected_source_indexes(raw_records, selector),
            cardinality=record_input.cardinality,
            input_index=input_index,
        )
        source_kind, source_paths = _source_details(record_input.source)
        expands = len(source_indexes) > 1
        for source_record_index in source_indexes:
            record = copy.deepcopy(raw_records[source_record_index])
            record = reverse_records(
                (record,),
                record_input.presentation.reverse_complement,
                log=logger,
            )[0]
            if record_input.region is not None:
                try:
                    record = apply_region_specs(
                        (record,),
                        (_unqualified_region(record_input.region),),
                        log=logger,
                    )[0]
                except ValueError as exc:
                    raise ValidationError(str(exc)) from exc
            base_key = record_input.record_key or f"record-{input_index + 1}"
            record_key = (
                f"{base_key}:{source_record_index + 1}"
                if expands
                else base_key
            )
            item = ResolvedRecordProvenance(
                resolved_index=len(records),
                input_index=input_index,
                source_record_index=source_record_index,
                source_record_count=len(raw_records),
                source_record_id=str(raw_records[source_record_index].id),
                source_kind=source_kind,
                source_paths=source_paths,
                record_key=record_key,
                cardinality=record_input.cardinality,
                selector=selector,
                region=record_input.region,
                presentation=record_input.presentation,
            )
            _apply_presentation(
                record,
                record_input.presentation,
                record_key=record_key,
            )
            _apply_provenance_annotations(record, item)
            records.append(record)
            provenance.append(item)
    records = _apply_collection_options(
        records,
        provenance,
        record_options or RecordCollectionOptions(),
    )
    return ResolvedRecordCollection(tuple(records), tuple(provenance))


def _normalized_cli_values(
    values: Sequence[object] | None,
    *,
    count: int,
    field_name: str,
    default: object,
) -> list[object]:
    normalized = list(values or ())
    if len(normalized) > count:
        raise ValidationError(
            f"Too many {field_name} values (expected at most {count})."
        )
    normalized.extend(default for _ in range(count - len(normalized)))
    return normalized


def _parse_cli_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value or "").strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "no", "n", "off", "", "none", "null", "-"}:
        return False
    raise ValidationError(f"Invalid reverse_complement value: {value}")


def record_input_manifest_from_paths(
    *,
    gbk_paths: Sequence[str] | None = None,
    gff_paths: Sequence[str] | None = None,
    fasta_paths: Sequence[str] | None = None,
    cardinalities: Sequence[RecordCardinality] | RecordCardinality,
    selectors: Sequence[str] | None = None,
    reverse_flags: Sequence[object] | None = None,
    labels: Sequence[str] | None = None,
    subtitles: Sequence[str] | None = None,
    regions: Sequence[str] | None = None,
) -> RecordInputManifest:
    """Translate CLI path lists into unresolved typed record inputs."""

    gbks = tuple(gbk_paths or ())
    gffs = tuple(gff_paths or ())
    fastas = tuple(fasta_paths or ())
    if bool(gbks) == bool(gffs):
        raise ValidationError(
            "Specify either GenBank paths or GFF3/FASTA paths."
        )
    if gffs and len(gffs) != len(fastas):
        raise ValidationError("GFF3 and FASTA path counts must match.")
    count = len(gbks) if gbks else len(gffs)
    if isinstance(cardinalities, RecordCardinality):
        cardinality_values = [cardinalities] * count
    else:
        cardinality_values = list(cardinalities)
        if len(cardinality_values) != count:
            raise ValidationError(
                "Record cardinality count must match the number of input sources."
            )
    if not all(
        isinstance(cardinality, RecordCardinality)
        for cardinality in cardinality_values
    ):
        raise ValidationError(
            "Record cardinalities must be RecordCardinality values."
        )
    selector_values = _normalized_cli_values(
        selectors,
        count=count,
        field_name="record_id",
        default="",
    )
    reverse_values = _normalized_cli_values(
        reverse_flags,
        count=count,
        field_name="reverse_complement",
        default=False,
    )
    records: list[RecordInput] = []
    source_paths: list[str] = []
    for index in range(count):
        source: RecordInputSource
        if gbks:
            source = GenBankInputSource(gbks[index])
            source_paths.append(gbks[index])
        else:
            source = GffFastaInputSource(gffs[index], fastas[index])
            source_paths.append(gffs[index])
        selector = parse_record_selector(str(selector_values[index] or ""))
        cardinality = cardinality_values[index]
        if selector is not None:
            cardinality = RecordCardinality.EXACTLY_ONE
        records.append(
            RecordInput(
                source=source,
                cardinality=cardinality,
                selector=selector,
                presentation=RecordPresentation(
                    reverse_complement=_parse_cli_bool(reverse_values[index]),
                ),
                record_key=f"record-{index + 1}",
            )
        )
    return RecordInputManifest(
        records=tuple(records),
        record_options=RecordCollectionOptions(
            regions=tuple(parse_region_spec(value) for value in (regions or ())),
            labels=tuple(labels or ()),
            subtitles=tuple(subtitles or ()),
        ),
        source_paths=tuple(source_paths),
    )


def record_input_manifest_from_table(path: str) -> RecordInputManifest:
    """Translate one records table into exact-one unresolved row inputs."""

    table = read_records_table(path)
    records: list[RecordInput] = []
    source_paths: list[str] = []
    for index, row in enumerate(table.rows):
        if table.input_kind == "gbk":
            source: RecordInputSource = GenBankInputSource(row.gbk)
            source_paths.append(row.gbk)
        else:
            source = GffFastaInputSource(row.gff, row.fasta)
            source_paths.append(row.gff)
        records.append(
            RecordInput(
                source=source,
                cardinality=RecordCardinality.EXACTLY_ONE,
                selector=parse_record_selector(row.record_id),
                region=parse_region_spec(row.region) if row.region else None,
                presentation=RecordPresentation(
                    label=row.record_label or None,
                    subtitle=row.record_subtitle or None,
                    reverse_complement=row.reverse_complement,
                    grid_row=row.row,
                    grid_column=row.column,
                ),
                record_key=f"record-{index + 1}",
            )
        )
    return RecordInputManifest(
        records=tuple(records),
        record_options=RecordCollectionOptions(),
        source_paths=tuple(source_paths),
        multi_record_positions=tuple(table.multi_record_positions()),
    )


def resolve_circular_options(
    options: CircularDiagramOptions,
) -> CircularDiagramOptions:
    """Resolve an optional Circular comparison manifest once."""

    if options.conservation_table_file is None:
        return options
    table = read_conservation_table(options.conservation_table_file)
    return replace(
        options,
        conservation_table_file=None,
        conservation_blast_files=tuple(table.conservation_blast_files),
        conservation_fasta_files=(
            tuple(table.comparison_fasta_files)
            if table.comparison_fasta_files is not None
            else None
        ),
        conservation_labels=(
            tuple(table.labels) if table.labels is not None else None
        ),
        conservation_colors=(
            tuple(table.colors) if table.colors is not None else None
        ),
    )


def _resolve_linear_comparison_selector(
    records: Sequence[SeqRecord],
    selector: str,
    *,
    table_path: str,
    row_number: int,
    column: str,
) -> int:
    if selector.startswith("#"):
        try:
            index = int(selector[1:]) - 1
        except ValueError as exc:
            raise ValidationError(
                f"{table_path}: row {row_number}, column {column!r}: "
                f"invalid record selector {selector!r}."
            ) from exc
        if 0 <= index < len(records):
            return index
    else:
        matches = [
            index
            for index, record in enumerate(records)
            if str(record.id) == selector
        ]
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ValidationError(
                f"{table_path}: row {row_number}, column {column!r}: "
                f"selector {selector!r} matched multiple record IDs; use #index."
            )
    raise ValidationError(
        f"{table_path}: row {row_number}, column {column!r}: unresolved "
        f"record selector {selector!r}."
    )


def resolve_linear_options(
    options: LinearDiagramOptions,
    *,
    records: Sequence[SeqRecord],
    layout: LinearMultiRecordOptions | None,
) -> LinearDiagramOptions:
    """Resolve an optional comparison manifest against displayed records."""

    positions = layout.multi_record_positions if layout is not None else None
    _ordered, rows_by_record = resolve_record_row_positions(records, positions)
    multi_record_rows = len(set(rows_by_record)) < len(records)
    if multi_record_rows and options.blast_files:
        raise ValidationError(
            "-b/--blast is ambiguous when a Linear row contains multiple records; "
            "use a comparison table with explicit query and subject selectors."
        )
    if options.comparison_table_file is None:
        return options
    table = read_comparisons_table(options.comparison_table_file)
    comparisons: list[LinearComparison] = []
    for row in table.rows:
        query_index = _resolve_linear_comparison_selector(
            records,
            row.query,
            table_path=table.table_path,
            row_number=row.row_number,
            column="query",
        )
        subject_index = _resolve_linear_comparison_selector(
            records,
            row.subject,
            table_path=table.table_path,
            row_number=row.row_number,
            column="subject",
        )
        query_row = int(rows_by_record[query_index])
        subject_row = int(rows_by_record[subject_index])
        if query_index == subject_index:
            raise ValidationError(
                f"{table.table_path}: row {row.row_number}, column 'subject': "
                "query and subject resolved to the same record."
            )
        if query_row == subject_row or abs(query_row - subject_row) != 1:
            topology = (
                "different rows"
                if query_row == subject_row
                else "adjacent rows"
            )
            raise ValidationError(
                f"{table.table_path}: row {row.row_number}, column 'subject': "
                f"query row {query_row + 1} and subject row "
                f"{subject_row + 1} must be in {topology}."
            )
        try:
            matches = pd.read_csv(
                row.blast,
                sep="\t",
                comment="#",
                names=COMPARISON_COLUMNS,
            )
        except (OSError, UnicodeError, pd.errors.ParserError) as exc:
            raise ValidationError(
                f"{table.table_path}: row {row.row_number}, column 'blast': "
                f"could not parse {row.blast}."
            ) from exc
        comparisons.append(
            LinearComparison(query_index, subject_index, matches)
        )
    return replace(
        options,
        comparison_table_file=None,
        linear_comparisons=tuple(comparisons),
    )


def resolve_circular_batch_outputs(
    policy: CircularBatchOutputPolicy,
    records: Sequence[SeqRecord],
) -> tuple[RenderOutputRequest, ...]:
    """Resolve one output request per displayed record after expansion."""

    prefix = str(policy.output_prefix) if policy.output_prefix is not None else None
    if prefix is None:
        raw_prefixes: list[str] = []
        used: set[str] = set()
        for record in records:
            base = resolve_implicit_record_output_prefix(record.id)
            candidate = base
            suffix = 2
            while candidate in used:
                candidate = f"{base}_{suffix}"
                suffix += 1
            used.add(candidate)
            raw_prefixes.append(candidate)
    elif len(records) == 1:
        raw_prefixes = [prefix]
    else:
        raw_prefixes = [
            f"{prefix}_{index}" for index in range(1, len(records) + 1)
        ]
    outputs: list[RenderOutputRequest] = []
    for raw_prefix in raw_prefixes:
        path = Path(raw_prefix)
        outputs.append(
            RenderOutputRequest(
                output_prefix=path.name,
                output_directory=(
                    path.parent if path.parent != Path(".") else None
                ),
                formats=policy.formats,
                overwrite=policy.overwrite,
                interactive_metadata_policy=policy.interactive_metadata_policy,
            )
        )
    return tuple(outputs)


def resolve_implicit_record_output_prefix(record_id: object) -> str:
    """Require a record-derived output name to be one filename component."""

    raw_record_id = str(record_id)
    try:
        output = RenderOutputRequest(output_prefix=raw_record_id)
    except ValidationError as exc:
        raise ValidationError(
            f"Record ID {raw_record_id!r} cannot be used as an implicit output "
            "filename prefix. Specify an explicit output prefix."
        ) from exc
    return output.output_prefix


__all__ = [
    "RecordInputManifest",
    "ResolvedRecordCollection",
    "ResolvedRecordProvenance",
    "depth_track_inputs_from_cli",
    "record_input_manifest_from_paths",
    "record_input_manifest_from_table",
    "resolve_circular_batch_outputs",
    "resolve_circular_options",
    "resolve_implicit_record_output_prefix",
    "resolve_linear_options",
    "resolve_record_inputs",
]
