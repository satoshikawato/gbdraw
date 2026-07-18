"""I/O helpers for the public API layer."""

from __future__ import annotations

from typing import Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]
from gbdraw.features.visibility import resolve_candidate_feature_types  # type: ignore[reportMissingImports]
from gbdraw.io.genome import (  # type: ignore[reportMissingImports]
    load_gbks as _load_gbks,
    load_gff_fasta as _load_gff_fasta,
)
from gbdraw.io.cli_tables import (  # type: ignore[reportMissingImports]
    CircularTrackTable,
    ConservationTable,
    ConservationTableRow,
    ComparisonTable,
    ComparisonTableRow,
    RecordsTable,
    RecordsTableRow,
    TablePathDependency,
    read_circular_track_table,
    read_conservation_table,
    read_comparisons_table,
    read_records_table,
)
from gbdraw.labels.filtering import (  # type: ignore[reportMissingImports]
    read_filter_list_file,
    read_label_override_file,
    read_qualifier_priority_file,
)
from gbdraw.io.record_select import (  # type: ignore[reportMissingImports]
    RecordSelector,
    parse_record_selector as _parse_record_selector,
    parse_record_selectors as _parse_record_selectors,
)
from gbdraw.io.regions import (  # type: ignore[reportMissingImports]
    RegionSpec,
    apply_region_specs as _apply_region_specs,
    parse_region_spec as _parse_region_spec,
    parse_region_specs as _parse_region_specs,
)


def load_gbks(
    gbk_list: Sequence[str],
    mode: str,
    load_comparison: bool = False,
    record_selectors: list[str] | None = None,
    reverse_flags: list[bool] | None = None,
) -> list[SeqRecord]:
    """Load GenBank files with consistent API exceptions."""

    try:
        return _load_gbks(
            list(gbk_list),
            mode=mode,
            load_comparison=load_comparison,
            record_selectors=record_selectors,
            reverse_flags=reverse_flags,
        )
    except ValueError as exc:
        raise ValidationError(str(exc)) from exc


def load_gff_fasta(
    gff_list: Sequence[str],
    fasta_list: Sequence[str],
    mode: str,
    selected_features_set=None,
    keep_all_features: bool = False,
    load_comparison: bool = False,
    record_selectors: list[str] | None = None,
    reverse_flags: list[bool] | None = None,
    color_table: DataFrame | None = None,
    feature_visibility_table: DataFrame | None = None,
) -> list[SeqRecord]:
    """Load paired GFF3 + FASTA files with consistent API exceptions."""

    try:
        candidate_features = selected_features_set
        resolved_keep_all_features = keep_all_features
        if not keep_all_features and (color_table is not None or feature_visibility_table is not None):
            candidate_types, table_requires_all = resolve_candidate_feature_types(
                selected_features_set,
                color_table=color_table,
                feature_visibility_table=feature_visibility_table,
            )
            candidate_features = candidate_types
            resolved_keep_all_features = table_requires_all
        return _load_gff_fasta(
            list(gff_list),
            list(fasta_list),
            mode=mode,
            selected_features_set=candidate_features,
            keep_all_features=resolved_keep_all_features,
            load_comparison=load_comparison,
            record_selectors=record_selectors,
            reverse_flags=reverse_flags,
        )
    except ValueError as exc:
        raise ValidationError(str(exc)) from exc


def parse_record_selector(text: str | None) -> RecordSelector | None:
    """Parse a record selector string, raising ValidationError on invalid input."""

    try:
        return _parse_record_selector(text)
    except ValueError as exc:
        raise ValidationError(str(exc)) from exc


def parse_record_selectors(specs: Sequence[str] | None) -> list[RecordSelector]:
    """Parse multiple record selectors, raising ValidationError on invalid input."""

    try:
        return _parse_record_selectors(specs)
    except ValueError as exc:
        raise ValidationError(str(exc)) from exc


def parse_region_spec(spec: str) -> RegionSpec:
    """Parse a region spec string, raising ValidationError on invalid input."""

    try:
        return _parse_region_spec(spec)
    except ValueError as exc:
        raise ValidationError(str(exc)) from exc


def parse_region_specs(specs: Sequence[str] | None) -> list[RegionSpec]:
    """Parse multiple region specs, raising ValidationError on invalid input."""

    try:
        return _parse_region_specs(specs)
    except ValueError as exc:
        raise ValidationError(str(exc)) from exc


def apply_region_specs(
    records: Sequence[SeqRecord],
    specs: Sequence[RegionSpec],
) -> list[SeqRecord]:
    """Apply region specs, raising ValidationError on invalid input."""

    try:
        return _apply_region_specs(records, specs)
    except ValueError as exc:
        raise ValidationError(str(exc)) from exc


def read_label_whitelist_table(path: str) -> DataFrame | None:
    """Read a label whitelist table using the CLI-compatible validator."""

    return read_filter_list_file(path)


def read_qualifier_priority_table(path: str) -> DataFrame | None:
    """Read a label qualifier-priority table."""

    return read_qualifier_priority_file(path)


def read_label_override_table(path: str) -> DataFrame | None:
    """Read a label override table."""

    return read_label_override_file(path)


__all__ = [
    "CircularTrackTable",
    "ConservationTable",
    "ConservationTableRow",
    "ComparisonTable",
    "ComparisonTableRow",
    "RecordsTable",
    "RecordsTableRow",
    "RegionSpec",
    "RecordSelector",
    "TablePathDependency",
    "apply_region_specs",
    "load_gbks",
    "load_gff_fasta",
    "parse_record_selector",
    "parse_record_selectors",
    "parse_region_spec",
    "parse_region_specs",
    "read_circular_track_table",
    "read_conservation_table",
    "read_comparisons_table",
    "read_label_override_table",
    "read_label_whitelist_table",
    "read_qualifier_priority_table",
    "read_records_table",
]
