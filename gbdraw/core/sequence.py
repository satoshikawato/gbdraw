#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, List, Union

from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition
from Bio.SeqRecord import SeqRecord

from ..features.overrides import feature_override_lookup
from ..features.visibility import should_render_feature
from .record_metadata import _source_feature_anchor_profile


def create_dict_for_sequence_lengths(records: Sequence[SeqRecord]) -> dict[str, int]:
    return {record.id: len(record.seq) for record in records}


def check_feature_presence(
    records: Union[List[SeqRecord], SeqRecord],
    features_list: List[str],
    feature_visibility_rules=None,
    specific_color_rules=None,
    record_features=(),
) -> list[str]:
    if isinstance(records, SeqRecord):
        records = [records]

    features_present: list[str] = []
    seen_feature_types: set[str] = set()

    for index, record in enumerate(records):
        override_of = feature_override_lookup(
            record, record_features[index].overrides if record_features else None
        )
        for feature in record.features:
            if not should_render_feature(
                feature,
                features_list,
                feature_visibility_rules=feature_visibility_rules,
                record_id=record.id,
                specific_color_rules=specific_color_rules,
                feature_override=override_of(feature),
            ):
                continue
            if feature.type in seen_feature_types:
                continue
            seen_feature_types.add(feature.type)
            features_present.append(feature.type)
    return features_present


def get_coordinates_of_longest_segment(feature_object):
    coords = feature_object.coordinates
    if not coords:
        return None, -1

    longest_segment_info = None
    max_length = -1

    for coord in coords:
        try:
            start, end = int(coord[2]), int(coord[3])
            length = abs(end - start)
            if length > max_length:
                max_length = length
                longest_segment_info = coord
        except (IndexError, TypeError, ValueError):
            continue

    return longest_segment_info, max_length


def _first_qualifier(qualifiers: dict, key: str) -> str | None:
    values = qualifiers.get(key)
    if isinstance(values, (list, tuple)):
        values = values[0] if values else None
    text = "" if values is None else str(values).strip()
    return text or None


def _cds_reading_frame(qualifiers: dict) -> int:
    """Return /codon_start, or the GFF3 phase of the 5'-most part plus one."""

    raw = _first_qualifier(qualifiers, "codon_start")
    name, shift = "codon_start", 0
    if raw is None:
        raw = _first_qualifier(qualifiers, "phase")
        name, shift = "phase", 1
    if raw is None:
        return 1
    try:
        frame = int(raw) + shift
    except ValueError:
        raise ValueError(f"{name} is invalid: {raw}") from None
    if frame not in (1, 2, 3):
        raise ValueError(f"{name} is outside {1 - shift}..{3 - shift}: {raw}")
    return frame


def _cds_five_prime_is_complete(feature: Any) -> bool:
    location = feature.location
    first_part = (getattr(location, "parts", None) or [location])[0]
    five_prime = first_part.end if first_part.strand == -1 else first_part.start
    if not isinstance(five_prime, ExactPosition):
        return False
    # GFF3 start_range/end_range name source columns, so use the source strand.
    source_strand = _source_feature_anchor_profile(feature).strand
    range_key = "end_range" if source_strand == "-" else "start_range"
    return range_key not in (feature.qualifiers or {})


def translate_cds(
    feature: Any,
    nucleotide_sequence: object,
    *,
    require_whole_codons: bool = False,
) -> str:
    """Translate a CDS without /translation from its extracted nucleotides.

    The reading frame is /codon_start or, for GFF3, the 5' part phase plus one;
    the table is /transl_table (default 1). As in INSDC /translation, which
    pseudo CDS do not carry, the first residue of a non-pseudo CDS is M when
    its 5' end is complete, the frame is 1, and the first codon is a start
    codon of the table. A trailing partial codon and one terminal stop are
    dropped. Raises ValueError for an invalid frame or table.
    """

    qualifiers = feature.qualifiers or {}
    raw_table = _first_qualifier(qualifiers, "transl_table") or "1"
    try:
        table_id = int(raw_table)
        start_codons = unambiguous_dna_by_id[table_id].start_codons
    except (KeyError, ValueError):
        raise ValueError(f"transl_table is invalid: {raw_table}") from None
    frame = _cds_reading_frame(qualifiers)
    coding = str(nucleotide_sequence).upper()[frame - 1 :]
    partial_codon = len(coding) % 3
    if partial_codon and require_whole_codons:
        raise ValueError("coding sequence length is not divisible by 3")
    coding = coding[: len(coding) - partial_codon]
    if not coding:
        raise ValueError("coding sequence has no complete codon")
    # A table id keeps Biopython's ambiguous table, so NNN becomes X.
    protein = str(Seq(coding).translate(table=table_id))
    if (
        frame == 1
        and coding[:3] in start_codons
        and "pseudo" not in qualifiers
        and "pseudogene" not in qualifiers
        and _cds_five_prime_is_complete(feature)
    ):
        protein = "M" + protein[1:]
    return protein[:-1] if protein.endswith("*") else protein


__all__ = [
    "check_feature_presence",
    "create_dict_for_sequence_lengths",
    "get_coordinates_of_longest_segment",
    "translate_cds",
]
