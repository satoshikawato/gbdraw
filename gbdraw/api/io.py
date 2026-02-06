"""I/O helpers for the public API layer."""

from __future__ import annotations

from typing import Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]
from gbdraw.io.genome import (  # type: ignore[reportMissingImports]
    load_gbks as _load_gbks,
    load_gff_fasta as _load_gff_fasta,
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
    load_comparison: bool = False,
    record_selectors: list[str] | None = None,
    reverse_flags: list[bool] | None = None,
) -> list[SeqRecord]:
    """Load paired GFF3 + FASTA files with consistent API exceptions."""

    try:
        return _load_gff_fasta(
            list(gff_list),
            list(fasta_list),
            mode=mode,
            selected_features_set=selected_features_set,
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


__all__ = [
    "RegionSpec",
    "RecordSelector",
    "apply_region_specs",
    "load_gbks",
    "load_gff_fasta",
    "parse_record_selector",
    "parse_record_selectors",
    "parse_region_spec",
    "parse_region_specs",
]
