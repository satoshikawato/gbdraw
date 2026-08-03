#!/usr/bin/env python
# coding: utf-8

"""Orthogroup-aware horizontal alignment helpers for linear diagrams."""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple, Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]

from ...analysis.protein_colinearity import OrthogroupResult  # type: ignore[reportMissingImports]
from ...canvas import LinearCanvasConfigurator  # type: ignore[reportMissingImports]
from ...core.record_metadata import (
    _mapped_feature_location_parts,
    _read_coord_map,
    _source_feature_index,
    _source_feature_location_parts,
)
from ...exceptions import ValidationError
from ...features.ids import compute_feature_hash_from_location_parts


@dataclass(frozen=True)
class OrthogroupAlignmentMember:
    orthogroup_id: str
    record_index: int
    protein_id: str
    source_protein_id: str
    feature_svg_id: str
    view_feature_svg_id: str
    feature_index: int
    center: float
    bitscore: float
    evalue: float
    identity: float
    representative: bool = False


@dataclass(frozen=True)
class OrthogroupAlignmentCanvasExtents:
    horizontal_shift: float
    width_extension: float
    min_left: float
    max_right: float

    @property
    def ruler_offset_x(self) -> float:
        return self.min_left

    @property
    def ruler_width(self) -> float:
        return max(0.0, self.max_right - self.min_left)


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


def _row_float(row: object, column: str, default: float = 0.0) -> float:
    value = _row_value(row, column, default)
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def _row_bool(row: object, column: str) -> bool:
    value = _row_value(row, column, False)
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _center_from_row(row: object, role: str) -> float:
    if role == "query":
        start = _row_float(row, "qstart", 0.0)
        end = _row_float(row, "qend", 0.0)
    else:
        start = _row_float(row, "sstart", 0.0)
        end = _row_float(row, "send", 0.0)
    return (min(start, end) + max(start, end)) / 2.0


def _member_from_row(row: object, role: str) -> OrthogroupAlignmentMember | None:
    orthogroup_id = _row_str(row, "orthogroup_id")
    if not orthogroup_id:
        return None
    return OrthogroupAlignmentMember(
        orthogroup_id=orthogroup_id,
        record_index=_row_int(row, f"{role}_record_index"),
        protein_id=_row_str(row, f"{role}_protein_id"),
        source_protein_id=_row_str(row, f"{role}_source_protein_id"),
        feature_svg_id=_row_str(row, f"{role}_feature_svg_id"),
        view_feature_svg_id=(
            _row_str(row, f"{role}_view_feature_svg_id")
            or _row_str(row, f"{role}_feature_svg_id")
        ),
        feature_index=_row_int(row, f"{role}_feature_index"),
        center=_center_from_row(row, role),
        bitscore=_row_float(row, "bitscore", 0.0),
        evalue=_row_float(row, "evalue", float("inf")),
        identity=_row_float(row, "identity", 0.0),
        representative=_row_bool(row, f"{role}_orthogroup_representative"),
    )


def _is_better_member(
    candidate: OrthogroupAlignmentMember,
    current: OrthogroupAlignmentMember | None,
) -> bool:
    if current is None:
        return True
    return (
        candidate.bitscore > current.bitscore
        or (candidate.bitscore == current.bitscore and candidate.evalue < current.evalue)
        or (
            candidate.bitscore == current.bitscore
            and candidate.evalue == current.evalue
            and candidate.identity > current.identity
        )
        or (
            candidate.bitscore == current.bitscore
            and candidate.evalue == current.evalue
            and candidate.identity == current.identity
            and candidate.protein_id < current.protein_id
        )
    )


def _collect_alignment_members(
    comparisons: Sequence[DataFrame],
) -> dict[str, list[OrthogroupAlignmentMember]]:
    members_by_key: dict[tuple[str, int, int | str], OrthogroupAlignmentMember] = {}
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
                if member.representative and existing is not None and not existing.representative:
                    members_by_key[key] = member
                elif _is_better_member(member, existing):
                    members_by_key[key] = member

    members_by_orthogroup: dict[str, list[OrthogroupAlignmentMember]] = {}
    for member in members_by_key.values():
        members_by_orthogroup.setdefault(member.orthogroup_id, []).append(member)
    for orthogroup_id, members in members_by_orthogroup.items():
        members_by_orthogroup[orthogroup_id] = sorted(
            members,
            key=lambda item: (item.record_index, item.center, item.protein_id, item.feature_svg_id),
        )
    return members_by_orthogroup


def _collect_alignment_members_from_orthogroups(
    orthogroups: OrthogroupResult,
) -> dict[str, list[OrthogroupAlignmentMember]]:
    members_by_orthogroup: dict[str, list[OrthogroupAlignmentMember]] = {}
    for orthogroup_id, members in orthogroups.orthogroups.items():
        group_members: list[OrthogroupAlignmentMember] = []
        for member in members:
            center = (min(float(member.start + 1), float(member.end)) + max(float(member.start + 1), float(member.end))) / 2.0
            group_members.append(
                OrthogroupAlignmentMember(
                    orthogroup_id=orthogroup_id,
                    record_index=int(member.record_index),
                    protein_id=str(member.protein_id or ""),
                    source_protein_id=str(member.source_protein_id or ""),
                    feature_svg_id=str(member.feature_svg_id or ""),
                    view_feature_svg_id="",
                    feature_index=int(getattr(member, "feature_index", -1)),
                    center=center,
                    bitscore=0.0,
                    evalue=float("inf"),
                    identity=0.0,
                    representative=bool(member.representative),
                )
            )
        members_by_orthogroup[orthogroup_id] = sorted(
            group_members,
            key=lambda item: (item.record_index, item.center, item.protein_id, item.feature_svg_id),
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
    member: OrthogroupAlignmentMember,
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
    orthogroups: OrthogroupResult | None = None,
    comparisons: Sequence[DataFrame] | None = None,
    records: Sequence[SeqRecord] | None = None,
) -> OrthogroupLabelEligibility:
    """Return all orthogroup feature IDs and the top-record IDs eligible for labels."""

    members_by_orthogroup = (
        _collect_alignment_members_from_orthogroups(orthogroups)
        if orthogroups is not None
        else _collect_alignment_members(comparisons or [])
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


def _target_matches(member: OrthogroupAlignmentMember, target: str) -> bool:
    return target in {
        member.orthogroup_id,
        member.protein_id,
        member.source_protein_id,
        member.feature_svg_id,
        member.view_feature_svg_id,
    }


def _resolve_target_member(
    members_by_orthogroup: dict[str, list[OrthogroupAlignmentMember]],
    target: str,
) -> tuple[str, OrthogroupAlignmentMember]:
    if target in members_by_orthogroup and members_by_orthogroup[target]:
        members = members_by_orthogroup[target]
        representatives = [member for member in members if member.representative]
        return target, (representatives[0] if representatives else members[0])

    matches = [
        (orthogroup_id, member)
        for orthogroup_id, members in members_by_orthogroup.items()
        for member in members
        if _target_matches(member, target)
    ]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise ValidationError(
            "align_orthogroup_feature matched multiple orthogroup members; "
            "use an unambiguous orthogroup, runtime protein, or feature ID."
        )
    raise ValidationError(
        "align_orthogroup_feature did not match any LOSATP blastp orthogroup member."
    )


def _base_record_offset_x(
    record: SeqRecord,
    canvas_config: LinearCanvasConfigurator,
) -> float:
    if bool(canvas_config.normalize_length):
        return 0.0
    if not bool(canvas_config.align_center):
        return 0.0
    longest = max(1.0, float(canvas_config.longest_genome))
    return float(canvas_config.alignment_width) * ((longest - float(len(record.seq))) / longest) / 2.0


def _rendered_center_x(
    member: OrthogroupAlignmentMember,
    records: Sequence[SeqRecord],
    canvas_config: LinearCanvasConfigurator,
) -> float:
    if member.record_index < 0 or member.record_index >= len(records):
        raise ValidationError(
            "align_orthogroup_feature references an orthogroup member outside the rendered records."
        )
    record = records[member.record_index]
    if bool(canvas_config.normalize_length):
        record_length = max(1.0, float(len(record.seq)))
        return float(canvas_config.alignment_width) * (member.center / record_length)
    longest = max(1.0, float(canvas_config.longest_genome))
    return _base_record_offset_x(record, canvas_config) + (
        float(canvas_config.alignment_width) * (member.center / longest)
    )


def _rendered_record_width(
    record: SeqRecord,
    canvas_config: LinearCanvasConfigurator,
) -> float:
    if bool(canvas_config.normalize_length):
        return float(canvas_config.alignment_width)
    longest = max(1.0, float(canvas_config.longest_genome))
    return float(canvas_config.alignment_width) * (float(len(record.seq)) / longest)


def calculate_orthogroup_alignment_canvas_adjustment(
    records: Sequence[SeqRecord],
    canvas_config: LinearCanvasConfigurator,
    record_offsets_x: dict[int, float],
) -> tuple[float, float]:
    """Return (horizontal_shift, width_extension) needed to keep aligned records on canvas."""

    extents = calculate_orthogroup_alignment_canvas_extents(
        records,
        canvas_config,
        record_offsets_x,
    )
    return extents.horizontal_shift, extents.width_extension


def calculate_orthogroup_alignment_canvas_extents(
    records: Sequence[SeqRecord],
    canvas_config: LinearCanvasConfigurator,
    record_offsets_x: dict[int, float],
) -> OrthogroupAlignmentCanvasExtents:
    """Return the aligned record bounds and canvas growth needed to fit them."""

    alignment_width = float(canvas_config.alignment_width)

    if not record_offsets_x:
        return OrthogroupAlignmentCanvasExtents(
            horizontal_shift=0.0,
            width_extension=0.0,
            min_left=0.0,
            max_right=alignment_width,
        )

    bounds: list[tuple[float, float]] = []
    for record_index, record in enumerate(records):
        local_left = (
            _base_record_offset_x(record, canvas_config)
            + float(record_offsets_x.get(record_index, 0.0))
        )
        local_right = local_left + _rendered_record_width(record, canvas_config)
        bounds.append((local_left, local_right))

    if not bounds:
        return OrthogroupAlignmentCanvasExtents(
            horizontal_shift=0.0,
            width_extension=0.0,
            min_left=0.0,
            max_right=alignment_width,
        )

    min_left = min(left for left, _ in bounds)
    max_right = max(right for _, right in bounds)

    horizontal_shift = max(0.0, -min_left)
    width_extension = max(0.0, max_right + horizontal_shift - alignment_width)
    return OrthogroupAlignmentCanvasExtents(
        horizontal_shift=horizontal_shift,
        width_extension=width_extension,
        min_left=min_left,
        max_right=max_right,
    )


def calculate_orthogroup_alignment_offsets(
    records: Sequence[SeqRecord],
    comparisons: Sequence[DataFrame],
    canvas_config: LinearCanvasConfigurator,
    align_orthogroup_feature: str | None,
    *,
    orthogroups: OrthogroupResult | None = None,
) -> dict[int, float]:
    """Return per-record x offsets that align representatives to the selected member."""

    target = str(align_orthogroup_feature or "").strip()
    if not target:
        return {}

    members_by_orthogroup = (
        _collect_alignment_members_from_orthogroups(orthogroups)
        if orthogroups is not None
        else _collect_alignment_members(comparisons)
    )
    if not members_by_orthogroup:
        raise ValidationError(
            "align_orthogroup_feature requires LOSATP blastp orthogroup metadata."
        )

    orthogroup_id, anchor_member = _resolve_target_member(members_by_orthogroup, target)
    target_members = members_by_orthogroup.get(orthogroup_id, [])
    if not target_members:
        return {}

    representative_by_record: dict[int, OrthogroupAlignmentMember] = {}
    for member in target_members:
        current = representative_by_record.get(member.record_index)
        if member.representative and (current is None or not current.representative):
            representative_by_record[member.record_index] = member
        elif current is None or (
            not current.representative and _is_better_member(member, current)
        ):
            representative_by_record[member.record_index] = member

    # The selected member is the visual anchor, even when a paralog in the same
    # record has a stronger score.
    representative_by_record[anchor_member.record_index] = anchor_member

    anchor_center_x = _rendered_center_x(anchor_member, records, canvas_config)
    offsets: dict[int, float] = {}
    for record_index, representative in representative_by_record.items():
        if record_index < 0 or record_index >= len(records):
            continue
        representative_center_x = _rendered_center_x(representative, records, canvas_config)
        offsets[record_index] = anchor_center_x - representative_center_x
    return offsets


__all__ = [
    "OrthogroupAlignmentCanvasExtents",
    "OrthogroupAlignmentMember",
    "OrthogroupLabelEligibility",
    "build_orthogroup_label_eligibility",
    "calculate_orthogroup_alignment_canvas_adjustment",
    "calculate_orthogroup_alignment_canvas_extents",
    "calculate_orthogroup_alignment_offsets",
    "orthogroup_label_sets_for_record",
]
