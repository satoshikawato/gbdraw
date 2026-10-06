"""Immutable original-source identities for record-aligned request planning."""

from __future__ import annotations

import csv
from collections import defaultdict
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Literal, overload

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame, isna

from gbdraw.core.record_metadata import (
    _feature_source_index_map,
    _iter_source_features,
    _mapped_feature_location_parts,
    _read_coord_map,
    _source_feature_index,
    _source_feature_location_parts,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import RecordSelector, parse_record_selector, select_record
from .ids import compute_feature_hash_from_location_parts, disambiguate_feature_ids
from .selector_values import (
    matches_feature_selector,
    normalize_qualifier_values,
    normalize_strand_token,
)


@dataclass(frozen=True)
class SourceFeatureIdentity:
    biological_feature_id: str
    source_feature_index: int
    feature_type: str
    location_parts: tuple[tuple[int, int, int | None], ...]
    qualifiers: tuple[tuple[str, tuple[str, ...]], ...]
    stable_feature_id: str

    def matches(self, *, key: str | None, value: str, record_id: str) -> bool:
        # hash= also names one of identical features by its `<hash>~<n>` ID.
        if key is not None and key.lower() == "hash" and value == self.biological_feature_id:
            return True
        start = min(part[0] for part in self.location_parts)
        end = max(part[1] for part in self.location_parts)
        strand = normalize_strand_token(self.location_parts[0][2])
        location = f"{start}..{end}"
        return matches_feature_selector(
            key=key,
            value=value,
            feature_type=self.feature_type,
            feature_hash=self.stable_feature_id,
            location=location,
            record_location=f"{record_id}:{location}:{strand}",
            qualifiers=dict(self.qualifiers),
        )


def build_source_feature_catalog(
    record: SeqRecord,
) -> tuple[SourceFeatureIdentity, ...]:
    """Capture before request crop/visibility; never mutate or reread the source."""
    base, step = _read_coord_map(record)
    entries = []
    for ordinal, feature in enumerate(_iter_source_features(record.features)):
        index = _source_feature_index(feature)
        parts = _source_feature_location_parts(
            feature
        ) or _mapped_feature_location_parts(
            feature,
            coord_base=base,
            coord_step=step,
        )
        if not parts:
            continue
        feature_type = str(feature.type)
        stable_id = compute_feature_hash_from_location_parts(
            feature_type, parts, record_id=record.id
        )
        entries.append(
            SourceFeatureIdentity(
                biological_feature_id=stable_id,
                source_feature_index=ordinal if index is None else index,
                feature_type=feature_type,
                location_parts=parts,
                qualifiers=tuple(
                    (str(key), tuple(normalize_qualifier_values(value)))
                    for key, value in (feature.qualifiers or {}).items()
                ),
                stable_feature_id=stable_id,
            )
        )
    if len({entry.source_feature_index for entry in entries}) != len(entries):
        raise ValidationError(
            "Source feature catalog contains duplicate source feature indexes."
        )
    ids = disambiguate_feature_ids(
        (str(record.id), entry.stable_feature_id, entry.source_feature_index)
        for entry in entries
    )
    return tuple(
        replace(entry, biological_feature_id=identity)
        for entry, identity in zip(entries, ids, strict=True)
    )


@dataclass(frozen=True)
class FeatureIdentity:
    """One original-source feature of one request record instance."""

    record_key: str
    biological_feature_id: str

    def __post_init__(self) -> None:
        for name in ("record_key", "biological_feature_id"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip() or "\0" in value:
                raise ValidationError(
                    f"{name} must be a non-empty identity without NUL.",
                    diagnostic={"code": "FEATURE_IDENTITY"},
                )
            object.__setattr__(self, name, value.strip())


@dataclass(frozen=True)
class IdentityBinding:
    """Where a request identity is in the drawn record of its record instance.

    ``present``: the drawn record has the feature. ``crop_excluded``: the drawn
    record is cropped and lacks it (outside the crop, or loading removed it).
    ``absent``: the uncropped drawn record lacks it; loading removed it (for
    example a GFF type filter). ``unresolved``: the original source has no such
    feature.
    """

    record_index: int
    source_feature_index: int | None
    status: Literal["present", "crop_excluded", "absent", "unresolved"]
    # The drawn feature when present, so consumers never resolve it again.
    feature: SeqFeature | None = field(default=None, compare=False, repr=False)


def _record_key_indexes(
    records: Sequence[SeqRecord],
    record_keys: Sequence[str],
    *aligned: Sequence[object],
) -> dict[str, int]:
    if any(len(values) != len(records) for values in (record_keys, *aligned)):
        raise ValidationError(
            "Feature identity context must align with records/provenance.",
            diagnostic={"code": "FEATURE_IDENTITY"},
        )
    indexes = {key: index for index, key in enumerate(record_keys)}
    if len(indexes) != len(record_keys):
        raise ValidationError(
            "Feature identity context contains duplicate record keys.",
            diagnostic={"code": "FEATURE_IDENTITY"},
        )
    return indexes


def resolve_feature_identities(
    *,
    records: Sequence[SeqRecord],
    record_keys: Sequence[str],
    source_catalogs: Sequence[tuple[SourceFeatureIdentity, ...]],
    identities: Iterable[FeatureIdentity],
) -> Mapping[FeatureIdentity, IdentityBinding]:
    """Bind identities to the drawn records of one request, aligned with provenance.

    Only the original-source catalog decides whether an identity exists; the
    runtime source ordinal decides whether the drawn record still has it. A record
    key outside the request is an error.
    """
    indexes = _record_key_indexes(records, record_keys, source_catalogs)
    wanted: dict[int, list[FeatureIdentity]] = defaultdict(list)
    for identity in identities:
        if identity.record_key not in indexes:
            raise ValidationError(
                f"Unknown feature identity record key {identity.record_key!r}.",
                diagnostic={"code": "FEATURE_IDENTITY"},
            )
        wanted[indexes[identity.record_key]].append(identity)
    bindings: dict[FeatureIdentity, IdentityBinding] = {}
    for record_index, items in wanted.items():
        record = records[record_index]
        known = {
            entry.biological_feature_id: entry.source_feature_index
            for entry in source_catalogs[record_index]
        }
        ordinals = _feature_source_index_map(record.features)
        drawn = {}
        for feature in _iter_source_features(record.features):
            index = _source_feature_index(feature)
            drawn[ordinals[id(feature)] if index is None else index] = feature
        missing: Literal["crop_excluded", "absent"] = (
            "crop_excluded" if record.annotations.get("gbdraw_region_applied") else "absent"
        )
        for identity in items:
            index = known.get(identity.biological_feature_id)
            feature = None if index is None else drawn.get(index)
            status: Literal["present", "crop_excluded", "absent", "unresolved"] = (
                "unresolved" if index is None else "present" if feature is not None else missing
            )
            bindings[identity] = IdentityBinding(record_index, index, status, feature)
    return bindings


def _identity_table_cell(value: object, *, verbatim: bool) -> str:
    if value is None or (not isinstance(value, str) and isna(value)):
        return ""
    if isinstance(value, float) and value.is_integer():
        value = int(value)  # pandas promotes an integer column with blank cells to float.
    return str(value) if verbatim else str(value).strip()


def read_identity_table(
    table: DataFrame | str | Path,
    *,
    table_name: str,
    columns: Sequence[str],
    required: Collection[str],
    verbatim: Collection[str] = (),
) -> list[dict[str, str]]:
    """Read an identity-addressed table: UTF-8 TSV (BOM allowed) or a DataFrame.

    Every allowed column is present in each row and a blank cell is ``""``. Cells
    are trimmed except ``verbatim`` columns. Row numbers count the header as 1.
    """
    def invalid(message: str, **context: object) -> ValidationError:
        return ValidationError(message, diagnostic={"code": "TABLE_INVALID", **context})

    if isinstance(table, DataFrame):
        header = list(table.columns)
        raw_rows = table.to_dict("records")
    else:
        try:
            with open(table, encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                header = list(reader.fieldnames or [])
                raw_rows = list(reader)
        except (OSError, UnicodeDecodeError, csv.Error) as exc:
            raise invalid(f"Cannot read {table_name}: {exc}") from exc
    if len(set(header)) != len(header) or set(header) - set(columns):
        raise invalid(f"{table_name} has duplicate or unknown columns; use {', '.join(columns)}.")
    missing = sorted(set(required) - set(header))
    if missing:
        raise invalid(
            f"{table_name} requires {' and '.join(missing)} column{'s' if len(missing) > 1 else ''}."
        )
    rows = []
    for row_number, raw in enumerate(raw_rows, start=2):
        if None in raw:
            raise invalid(f"{table_name} row {row_number} has more values than columns.", row=row_number)
        if any(isinstance(value, (Mapping, list, tuple, set)) for value in raw.values()):
            raise invalid(f"{table_name} row {row_number}: cells must be scalar values.", row=row_number)
        rows.append({
            column: _identity_table_cell(raw.get(column), verbatim=column in verbatim)
            for column in columns
        })
    return rows


def _selects_no_record(records: Sequence[SeqRecord], selector: RecordSelector | None) -> bool:
    if selector is None:
        return not records
    if selector.record_index is not None:
        return selector.record_index >= len(records)
    return all(record.id != selector.record_id for record in records)


@overload
def resolve_identity_table_rows(
    rows: Sequence[tuple[str, str]],
    *,
    table: str,
    records: Sequence[SeqRecord],
    record_keys: Sequence[str],
    source_record_ids: Sequence[str],
    source_catalogs: Sequence[tuple[SourceFeatureIdentity, ...]],
    unmatched: None = None,
) -> tuple[FeatureIdentity, ...]: ...


@overload
def resolve_identity_table_rows(
    rows: Sequence[tuple[str, str]],
    *,
    table: str,
    records: Sequence[SeqRecord],
    record_keys: Sequence[str],
    source_record_ids: Sequence[str],
    source_catalogs: Sequence[tuple[SourceFeatureIdentity, ...]],
    unmatched: list[int] | None = None,
) -> tuple[FeatureIdentity | None, ...]: ...


def resolve_identity_table_rows(
    rows: Sequence[tuple[str, str]],
    *,
    table: str,
    records: Sequence[SeqRecord],
    record_keys: Sequence[str],
    source_record_ids: Sequence[str],
    source_catalogs: Sequence[tuple[SourceFeatureIdentity, ...]],
    unmatched: list[int] | None = None,
) -> tuple[FeatureIdentity | None, ...]:
    """Resolve ``(record, feature_selector)`` cells to one identity per row.

    The record selector must select exactly one record and the feature selector
    exactly one original-source feature of it. Row numbers count the header as 1.
    With ``unmatched`` (the Web's Load Feature Edits TSV, Owner Q3 = A), a row
    that selects no record or no feature adds its row number there and resolves
    to ``None``; every other defect still fails.
    """
    # Deferred: the annotations package imports the feature factory, which uses this module.
    from gbdraw.annotations.models import parse_feature_selector

    _record_key_indexes(records, record_keys, source_record_ids, source_catalogs)
    identities: list[FeatureIdentity | None] = []
    seen: set[FeatureIdentity] = set()
    # hash= rows (Run Info writes one per edited feature) look up an index built
    # once per record instead of scanning the catalog per row; same matches.
    hash_indexes: dict[int, dict[str, list[SourceFeatureIdentity]]] = {}
    for row_number, (record_cell, selector_cell) in enumerate(rows, start=2):
        try:
            record_selector = parse_record_selector(record_cell)
            if unmatched is not None and _selects_no_record(records, record_selector):
                unmatched.append(row_number)
                identities.append(None)
                continue
            selected = select_record(records, record_selector)
            if len(selected) != 1:
                raise ValueError("Record selector must match exactly one record; use #index.")
            index = next(i for i, record in enumerate(records) if record is selected[0])
            selector = parse_feature_selector(selector_cell)
            if selector.key is not None and selector.key.lower() == "hash":
                if index not in hash_indexes:
                    by_hash: dict[str, list[SourceFeatureIdentity]] = defaultdict(list)
                    for entry in source_catalogs[index]:
                        for value in dict.fromkeys((entry.biological_feature_id, entry.stable_feature_id)):
                            by_hash[value].append(entry)
                    hash_indexes[index] = by_hash
                matched = hash_indexes[index].get(selector.value, [])
            else:
                matched = [
                    entry
                    for entry in source_catalogs[index]
                    if entry.matches(
                        key=selector.key, value=selector.value, record_id=source_record_ids[index]
                    )
                ]
            if not matched and unmatched is not None:
                unmatched.append(row_number)
                identities.append(None)
                continue
            if len(matched) != 1:
                raise ValueError(
                    "Feature selector must match exactly one source feature; "
                    f"matched {len(matched)}."
                )
            identity = FeatureIdentity(record_keys[index], matched[0].biological_feature_id)
            if identity in seen:
                raise ValueError("Duplicate resolved feature identity.")
        except (ValueError, ValidationError) as exc:
            raise ValidationError(
                f"{table} row {row_number}: {exc}",
                diagnostic={"code": "TABLE_INVALID", "row": row_number},
            ) from exc
        seen.add(identity)
        identities.append(identity)
    return tuple(identities)
