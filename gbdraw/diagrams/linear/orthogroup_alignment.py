#!/usr/bin/env python
# coding: utf-8

"""Orthogroup-aware horizontal alignment helpers for linear diagrams."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]

from ...analysis.protein_colinearity import OrthogroupResult  # type: ignore[reportMissingImports]
from ...canvas import LinearCanvasConfigurator  # type: ignore[reportMissingImports]
from ...exceptions import ValidationError


@dataclass(frozen=True)
class OrthogroupAlignmentMember:
    orthogroup_id: str
    record_index: int
    protein_id: str
    source_protein_id: str
    feature_svg_id: str
    center: float
    bitscore: float
    evalue: float
    identity: float
    representative: bool = False


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
    members_by_key: dict[tuple[int, str, str, str], OrthogroupAlignmentMember] = {}
    for comparison in comparisons:
        if comparison is None or comparison.empty or "orthogroup_id" not in comparison.columns:
            continue
        for row in comparison.itertuples(index=False):
            for role in ("query", "subject"):
                member = _member_from_row(row, role)
                if member is None or member.record_index < 0:
                    continue
                key = (
                    member.record_index,
                    member.protein_id,
                    member.source_protein_id,
                    member.feature_svg_id,
                )
                existing = members_by_key.get(key)
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


def _target_matches(member: OrthogroupAlignmentMember, target: str) -> bool:
    return target in {
        member.orthogroup_id,
        member.protein_id,
        member.source_protein_id,
        member.feature_svg_id,
    }


def _resolve_target_member(
    members_by_orthogroup: dict[str, list[OrthogroupAlignmentMember]],
    target: str,
) -> tuple[str, OrthogroupAlignmentMember]:
    if target in members_by_orthogroup and members_by_orthogroup[target]:
        members = members_by_orthogroup[target]
        representatives = [member for member in members if member.representative]
        return target, (representatives[0] if representatives else members[0])

    for orthogroup_id, members in members_by_orthogroup.items():
        for member in members:
            if _target_matches(member, target):
                return orthogroup_id, member
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

    if not record_offsets_x:
        return 0.0, 0.0

    alignment_width = float(canvas_config.alignment_width)
    min_left = 0.0
    max_right = alignment_width
    for record_index, record in enumerate(records):
        local_left = (
            _base_record_offset_x(record, canvas_config)
            + float(record_offsets_x.get(record_index, 0.0))
        )
        local_right = local_left + _rendered_record_width(record, canvas_config)
        min_left = min(min_left, local_left)
        max_right = max(max_right, local_right)

    horizontal_shift = max(0.0, -min_left)
    width_extension = max(0.0, (max_right - min_left) - alignment_width)
    return horizontal_shift, width_extension


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
    "OrthogroupAlignmentMember",
    "calculate_orthogroup_alignment_canvas_adjustment",
    "calculate_orthogroup_alignment_offsets",
]
