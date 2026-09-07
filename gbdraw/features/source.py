"""Immutable original-source identities for record-aligned request planning."""

from __future__ import annotations

from dataclasses import dataclass, replace

from Bio.SeqRecord import SeqRecord

from gbdraw.core.record_metadata import (
    _iter_source_features,
    _mapped_feature_location_parts,
    _read_coord_map,
    _source_feature_index,
    _source_feature_location_parts,
)
from gbdraw.exceptions import ValidationError
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
