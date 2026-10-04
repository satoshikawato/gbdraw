#!/usr/bin/env python
# coding: utf-8

"""Orthogroup-aware label eligibility for linear diagrams."""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple, Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]

from ...analysis.protein_colinearity import OrthogroupResult, OrthogroupGraphResult  # type: ignore[reportMissingImports]
from ...core.record_metadata import (
    _mapped_feature_location_parts,
    _read_coord_map,
    _source_feature_index,
    _source_feature_location_parts,
)
from ...exceptions import ValidationError
from ...features.ids import compute_feature_hash_from_location_parts


@dataclass(frozen=True)
class OrthogroupLabelMember:
    orthogroup_id: str
    record_index: int
    protein_id: str
    feature_svg_id: str
    view_feature_svg_id: str
    feature_index: int


class OrthogroupLabelEligibility(NamedTuple):
    member_ids_by_record: dict[int, set[str | int]]
    top_member_ids_by_record: dict[int, set[str | int]]


def _row_value(row: object, column: str, default: object = "") -> object:
    return getattr(row, column, default)


def _row_str(row: object, column: str) -> str:
    value = _row_value(row, column, "")
    if value is None:
        return ""
    text = str(value).strip()
    return "" if text.lower() == "nan" else text


def _row_int(row: object, column: str, default: int = -1) -> int:
    value = _row_value(row, column, default)
    try:
        return int(value)
    except (TypeError, ValueError):
        return int(default)


def _member_from_row(row: object, role: str) -> OrthogroupLabelMember | None:
    orthogroup_id = _row_str(row, "orthogroup_id")
    if not orthogroup_id:
        return None
    return OrthogroupLabelMember(
        orthogroup_id=orthogroup_id,
        record_index=_row_int(row, f"{role}_record_index"),
        protein_id=_row_str(row, f"{role}_protein_id"),
        feature_svg_id=_row_str(row, f"{role}_feature_svg_id"),
        view_feature_svg_id=(
            _row_str(row, f"{role}_view_feature_svg_id")
            or _row_str(row, f"{role}_feature_svg_id")
        ),
        feature_index=_row_int(row, f"{role}_feature_index"),
    )


def _collect_label_members(
    comparisons: Sequence[DataFrame],
) -> dict[str, list[OrthogroupLabelMember]]:
    members_by_key: dict[tuple[str, int, int | str], OrthogroupLabelMember] = {}
    for comparison in comparisons:
        if comparison is None or comparison.empty or "orthogroup_id" not in comparison.columns:
            continue
        for row in comparison.itertuples(index=False):
            for role in ("query", "subject"):
                member = _member_from_row(row, role)
                if member is None or member.record_index < 0:
                    continue
                if member.feature_index >= 0:
                    key = ("source", member.record_index, member.feature_index)
                elif member.feature_svg_id:
                    key = ("stable", member.record_index, member.feature_svg_id)
                else:
                    continue
                existing = members_by_key.get(key)
                if (
                    existing is not None
                    and existing.orthogroup_id != member.orthogroup_id
                ):
                    raise ValidationError(
                        "Orthogroup alignment rows assign one feature to conflicting orthogroups."
                    )
                if (
                    existing is not None
                    and existing.feature_svg_id
                    and member.feature_svg_id
                    and existing.feature_svg_id != member.feature_svg_id
                ):
                    raise ValidationError(
                        "Orthogroup alignment rows contain conflicting stable feature identity."
                    )
                if existing is None:
                    members_by_key[key] = member

    members_by_orthogroup: dict[str, list[OrthogroupLabelMember]] = {}
    for member in members_by_key.values():
        members_by_orthogroup.setdefault(member.orthogroup_id, []).append(member)
    for orthogroup_id, members in members_by_orthogroup.items():
        members_by_orthogroup[orthogroup_id] = sorted(
            members,
            key=lambda item: (
                item.record_index,
                item.feature_index,
                item.protein_id,
                item.feature_svg_id,
            ),
        )
    return members_by_orthogroup


def _collect_label_members_from_orthogroups(
    orthogroups: OrthogroupResult | OrthogroupGraphResult,
) -> dict[str, list[OrthogroupLabelMember]]:
    members_by_orthogroup: dict[str, list[OrthogroupLabelMember]] = {}
    for orthogroup_id, members in orthogroups.orthogroups.items():
        group_members: list[OrthogroupLabelMember] = []
        for member in members:
            group_members.append(
                OrthogroupLabelMember(
                    orthogroup_id=orthogroup_id,
                    record_index=int(member.record_index),
                    protein_id=str(member.protein_id or ""),
                    feature_svg_id=str(member.feature_svg_id or ""),
                    view_feature_svg_id="",
                    feature_index=int(getattr(member, "feature_index", -1)),
                )
            )
        members_by_orthogroup[orthogroup_id] = sorted(
            group_members,
            key=lambda item: (
                item.record_index,
                item.feature_index,
                item.protein_id,
                item.feature_svg_id,
            ),
        )
    return members_by_orthogroup


def _features_by_source_index(record: SeqRecord) -> dict[int, object]:
    by_source_index: dict[int, object] = {}
    fallback_index = 0

    def walk(features: object) -> None:
        nonlocal fallback_index
        for feature in features or ():  # type: ignore[union-attr]
            source_index = _source_feature_index(feature)
            resolved_index = fallback_index if source_index is None else source_index
            fallback_index += 1
            existing = by_source_index.get(resolved_index)
            if existing is not None and existing is not feature:
                raise ValidationError(
                    "Linear record contains duplicate source feature indexes."
                )
            by_source_index[resolved_index] = feature
            walk(getattr(feature, "sub_features", None))

    walk(record.features)
    return by_source_index


def _biological_feature_svg_id(record: SeqRecord, feature: object) -> str:
    parts = _source_feature_location_parts(feature)
    if parts is None:
        coord_base, coord_step = _read_coord_map(record)
        parts = _mapped_feature_location_parts(
            feature,
            coord_base=coord_base,
            coord_step=coord_step,
        )
    if not parts:
        return ""
    return compute_feature_hash_from_location_parts(
        str(getattr(feature, "type", "") or ""),
        parts,
        record_id=record.id,
    )


def _member_label_identity(
    member: OrthogroupLabelMember,
    *,
    record_features: dict[int, dict[int, object]] | None,
    records: Sequence[SeqRecord] | None,
) -> str | int:
    if (
        record_features is not None
        and records is not None
        and member.feature_index >= 0
        and 0 <= member.record_index < len(records)
    ):
        feature = record_features.get(member.record_index, {}).get(
            member.feature_index
        )
        if feature is None:
            raise ValidationError(
                "Orthogroup member source feature index is absent from its record."
            )
        biological_id = _biological_feature_svg_id(
            records[member.record_index],
            feature,
        )
        if member.feature_svg_id and biological_id != member.feature_svg_id:
            raise ValidationError(
                "Orthogroup member source feature index conflicts with its stable feature ID."
            )
        return member.feature_index
    return member.view_feature_svg_id or member.feature_svg_id


def build_orthogroup_label_eligibility(
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None = None,
    comparisons: Sequence[DataFrame] | None = None,
    records: Sequence[SeqRecord] | None = None,
) -> OrthogroupLabelEligibility:
    """Return all orthogroup feature IDs and the top-record IDs eligible for labels."""

    members_by_orthogroup = (
        _collect_label_members_from_orthogroups(orthogroups)
        if orthogroups is not None
        else _collect_label_members(comparisons or [])
    )
    member_ids_by_record: dict[int, set[str | int]] = {}
    top_member_ids_by_record: dict[int, set[str | int]] = {}
    record_features = (
        {
            record_index: _features_by_source_index(record)
            for record_index, record in enumerate(records)
        }
        if records is not None
        else None
    )
    for members in members_by_orthogroup.values():
        members_with_ids = [
            (member, _member_label_identity(
                member,
                record_features=record_features,
                records=records,
            ))
            for member in members
        ]
        members_with_ids = [
            item
            for item in members_with_ids
            if item[1] is not None and item[1] != ""
        ]
        if not members_with_ids:
            continue
        for member, identity in members_with_ids:
            member_ids_by_record.setdefault(member.record_index, set()).add(identity)
        top_record_index = min(member.record_index for member, _identity in members_with_ids)
        for member, identity in members_with_ids:
            if member.record_index == top_record_index:
                top_member_ids_by_record.setdefault(member.record_index, set()).add(identity)
    return OrthogroupLabelEligibility(member_ids_by_record, top_member_ids_by_record)


def orthogroup_label_sets_for_record(
    eligibility: OrthogroupLabelEligibility | None,
    record_index: int,
) -> tuple[set[str | int] | None, set[str | int] | None]:
    if eligibility is None:
        return None, None
    return (
        eligibility.member_ids_by_record.get(record_index, set()),
        eligibility.top_member_ids_by_record.get(record_index, set()),
    )


__all__ = [
    "OrthogroupLabelEligibility",
    "build_orthogroup_label_eligibility",
    "orthogroup_label_sets_for_record",
]
