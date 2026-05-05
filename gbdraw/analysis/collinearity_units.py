#!/usr/bin/env python
# coding: utf-8

"""Collinearity unit resolution for native block calling."""

from __future__ import annotations

from dataclasses import dataclass
import logging
from typing import Literal, Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from gbdraw.analysis.protein_colinearity import CdsProtein, ProteinExtractionResult
from gbdraw.exceptions import ValidationError

logger = logging.getLogger(__name__)

CollinearityUnitMode = Literal["auto", "cds", "locus"]


@dataclass(frozen=True)
class CollinearityUnit:
    """A block-calling unit mapped to one representative CDS feature."""

    unit_id: str
    unit_kind: Literal["locus", "cds"]
    record_index: int
    record_id: str
    order: int
    representative_protein_id: str
    representative_feature_svg_id: str
    start: int
    end: int
    strand: int | None
    locus_id: str | None
    display_name: str
    cds_members: tuple[str, ...]
    aliases: tuple[str, ...]


@dataclass(frozen=True)
class CollinearityUnitIndex:
    """Lookup tables for collinearity units."""

    units_by_record: list[list[CollinearityUnit]]
    unit_by_id: dict[str, CollinearityUnit]
    unit_by_protein_id: dict[str, CollinearityUnit]
    aliases_by_record: list[dict[str, str]]
    ambiguous_aliases_by_record: list[set[str]]


def normalize_collinearity_unit_mode(mode: str | None) -> CollinearityUnitMode:
    normalized = str(mode or "auto").strip().lower()
    if normalized not in {"auto", "cds", "locus"}:
        raise ValidationError("collinear_unit_mode must be one of: auto, cds, locus")
    return normalized  # type: ignore[return-value]


def _geneid_db_xref(protein: CdsProtein) -> str | None:
    for value in getattr(protein, "db_xref", ()) or ():
        text = str(value).strip()
        if text.startswith("GeneID:") and text != "GeneID:":
            return text
    return None


def strong_locus_id(protein: CdsProtein) -> str | None:
    """Return the stable locus identifier used for locus-level collapse."""

    candidates = (
        getattr(protein, "gene_parent_id", None),
        getattr(protein, "locus_tag", None),
        getattr(protein, "gene_id", None),
        _geneid_db_xref(protein),
    )
    for candidate in candidates:
        text = str(candidate or "").strip()
        if text:
            return text
    return None


def _representative_sort_key(protein: CdsProtein) -> tuple[int, bool, int, int, str]:
    span = max(0, int(protein.end) - int(protein.start))
    return (
        -int(protein.protein_length),
        not bool(protein.source_protein_id),
        -span,
        int(protein.feature_index),
        str(protein.protein_id),
    )


def _unit_sort_key(proteins: Sequence[CdsProtein]) -> tuple[int, int, int, str]:
    starts = [int(protein.start) for protein in proteins]
    ends = [int(protein.end) for protein in proteins]
    feature_indices = [int(protein.feature_index) for protein in proteins]
    protein_ids = [str(protein.protein_id) for protein in proteins]
    return (min(starts), min(ends), min(feature_indices), min(protein_ids))


def _unit_strand(proteins: Sequence[CdsProtein]) -> int | None:
    strands = {protein.strand for protein in proteins if protein.strand in {-1, 1}}
    if len(strands) == 1:
        return strands.pop()
    return None


def _add_alias(aliases: list[str], value: object | None) -> None:
    text = str(value or "").strip()
    if text and text not in aliases:
        aliases.append(text)


def _unit_aliases(
    *,
    unit_id: str,
    proteins: Sequence[CdsProtein],
    representative: CdsProtein,
    locus_id: str | None,
    display_name: str,
) -> tuple[str, ...]:
    aliases: list[str] = []
    _add_alias(aliases, unit_id)
    _add_alias(aliases, representative.protein_id)
    _add_alias(aliases, representative.source_protein_id)
    _add_alias(aliases, representative.feature_svg_id)
    _add_alias(aliases, locus_id)
    _add_alias(aliases, display_name)
    for protein in proteins:
        _add_alias(aliases, protein.protein_id)
        _add_alias(aliases, protein.source_protein_id)
        _add_alias(aliases, protein.feature_svg_id)
        _add_alias(aliases, getattr(protein, "locus_tag", None))
        _add_alias(aliases, getattr(protein, "gene_id", None))
        _add_alias(aliases, getattr(protein, "old_locus_tag", None))
        _add_alias(aliases, getattr(protein, "gene", None))
    return tuple(aliases)


def _record_ids_from_records(
    records: Sequence[SeqRecord] | None,
    extraction: ProteinExtractionResult,
) -> list[str]:
    if records is not None:
        return [str(record.id) for record in records]
    record_ids: list[str] = []
    for record_proteins in extraction.proteins_by_record:
        if record_proteins:
            record_ids.append(str(record_proteins[0].record_id))
        else:
            record_ids.append("")
    return record_ids


def build_collinearity_unit_index(
    extraction: ProteinExtractionResult,
    *,
    records: Sequence[SeqRecord] | None = None,
    mode: CollinearityUnitMode | str = "auto",
) -> CollinearityUnitIndex:
    """Resolve CDS proteins into ordered collinearity units."""

    normalized_mode = normalize_collinearity_unit_mode(str(mode))
    record_ids = _record_ids_from_records(records, extraction)
    units_by_record: list[list[CollinearityUnit]] = []
    unit_by_id: dict[str, CollinearityUnit] = {}
    unit_by_protein_id: dict[str, CollinearityUnit] = {}
    global_unit_index = 0
    missing_locus_examples: list[str] = []
    mixed_fallback_counts: list[tuple[str, int]] = []
    collapsed_counts: list[tuple[str, int]] = []

    for record_index, proteins in enumerate(extraction.proteins_by_record):
        sorted_proteins = sorted(proteins, key=lambda protein: (int(protein.start), int(protein.end), int(protein.feature_index), str(protein.protein_id)))
        grouped: list[tuple[Literal["locus", "cds"], str | None, list[CdsProtein]]] = []

        if normalized_mode == "cds":
            grouped = [("cds", strong_locus_id(protein), [protein]) for protein in sorted_proteins]
        else:
            locus_groups: dict[str, list[CdsProtein]] = {}
            fallback_proteins: list[CdsProtein] = []
            for protein in sorted_proteins:
                locus_id = strong_locus_id(protein)
                if locus_id:
                    locus_groups.setdefault(locus_id, []).append(protein)
                elif normalized_mode == "locus":
                    missing_locus_examples.append(f"{protein.record_id}:{protein.protein_id}")
                else:
                    fallback_proteins.append(protein)
            if missing_locus_examples and normalized_mode == "locus":
                examples = ", ".join(missing_locus_examples[:5])
                raise ValidationError(
                    "collinear_unit_mode='locus' requires stable locus identifiers for all CDS proteins; "
                    f"missing examples: {examples}"
                )
            grouped.extend(
                ("locus", locus_id, members)
                for locus_id, members in locus_groups.items()
            )
            grouped.extend(("cds", strong_locus_id(protein), [protein]) for protein in fallback_proteins)

        grouped = sorted(grouped, key=lambda item: _unit_sort_key(item[2]))
        if normalized_mode == "auto":
            fallback_count = sum(1 for unit_kind, _locus_id, _members in grouped if unit_kind == "cds")
            collapsed_count = sum(len(members) - 1 for unit_kind, _locus_id, members in grouped if unit_kind == "locus" and len(members) > 1)
            if fallback_count:
                mixed_fallback_counts.append((record_ids[record_index] if record_index < len(record_ids) else str(record_index + 1), fallback_count))
            if collapsed_count:
                collapsed_counts.append((record_ids[record_index] if record_index < len(record_ids) else str(record_index + 1), collapsed_count))

        record_units: list[CollinearityUnit] = []
        for order, (unit_kind, locus_id, members) in enumerate(grouped):
            representative = min(members, key=_representative_sort_key)
            global_unit_index += 1
            unit_id = f"gbd_r{record_index + 1:04d}_unit{global_unit_index:06d}"
            display_name = str(locus_id or representative.label or representative.protein_id)
            start = min(int(protein.start) for protein in members)
            end = max(int(protein.end) for protein in members)
            unit = CollinearityUnit(
                unit_id=unit_id,
                unit_kind=unit_kind,
                record_index=record_index,
                record_id=record_ids[record_index] if record_index < len(record_ids) else representative.record_id,
                order=order,
                representative_protein_id=representative.protein_id,
                representative_feature_svg_id=str(representative.feature_svg_id or ""),
                start=start,
                end=end,
                strand=_unit_strand(members),
                locus_id=locus_id,
                display_name=display_name,
                cds_members=tuple(
                    protein.protein_id
                    for protein in sorted(
                        members,
                        key=lambda protein: (
                            int(protein.start),
                            int(protein.end),
                            int(protein.feature_index),
                            str(protein.protein_id),
                        ),
                    )
                ),
                aliases=(),
            )
            unit = CollinearityUnit(
                **{
                    **unit.__dict__,
                    "aliases": _unit_aliases(
                        unit_id=unit.unit_id,
                        proteins=members,
                        representative=representative,
                        locus_id=locus_id,
                        display_name=display_name,
                    ),
                }
            )
            record_units.append(unit)
            unit_by_id[unit.unit_id] = unit
            for protein in members:
                unit_by_protein_id[str(protein.protein_id)] = unit
        units_by_record.append(record_units)

    aliases_by_record: list[dict[str, str]] = []
    ambiguous_aliases_by_record: list[set[str]] = []
    for record_units in units_by_record:
        alias_targets: dict[str, set[str]] = {}
        for unit in record_units:
            for alias in unit.aliases:
                alias_targets.setdefault(alias, set()).add(unit.unit_id)
        aliases_by_record.append(
            {
                alias: next(iter(unit_ids))
                for alias, unit_ids in alias_targets.items()
                if len(unit_ids) == 1
            }
        )
        ambiguous_aliases_by_record.append(
            {alias for alias, unit_ids in alias_targets.items() if len(unit_ids) > 1}
        )

    if mixed_fallback_counts:
        preview = ", ".join(f"{record_id}: {count}" for record_id, count in mixed_fallback_counts[:5])
        logger.warning(
            "WARNING: collinear_unit_mode auto used CDS units for CDS proteins without stable locus IDs (%s).",
            preview,
        )
    if collapsed_counts:
        preview = ", ".join(f"{record_id}: {count}" for record_id, count in collapsed_counts[:5])
        logger.info(
            "INFO: Collapsed multiple CDS proteins into locus collinearity units (%s).",
            preview,
        )

    return CollinearityUnitIndex(
        units_by_record=units_by_record,
        unit_by_id=unit_by_id,
        unit_by_protein_id=unit_by_protein_id,
        aliases_by_record=aliases_by_record,
        ambiguous_aliases_by_record=ambiguous_aliases_by_record,
    )


__all__ = [
    "CollinearityUnit",
    "CollinearityUnitIndex",
    "CollinearityUnitMode",
    "build_collinearity_unit_index",
    "normalize_collinearity_unit_mode",
    "strong_locus_id",
]
