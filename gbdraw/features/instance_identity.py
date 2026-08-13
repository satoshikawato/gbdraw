"""Exact, record-scoped identities for source feature instances."""

from __future__ import annotations

import base64
import hashlib
import re
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType

from Bio.SeqRecord import SeqRecord

from gbdraw.core.record_metadata import (
    _feature_source_index_map,
    _mapped_feature_location_parts,
    _read_coord_map,
    _source_feature_index,
    _source_feature_location_parts,
)
from gbdraw.exceptions import ValidationError

from .ids import compute_feature_hash_from_location_parts


FEATURE_INSTANCE_HASH_QUALIFIER = "__gbdraw_instance_hash__"
FEATURE_INSTANCE_HASH_PREFIX = "fi1_"
FEATURE_INSTANCE_HASH_PATTERN = re.compile(r"^fi1_[a-z2-7]{26}$")

_FEATURE_INSTANCE_HASH_DOMAIN = b"gbdraw-feature-instance-v1\0"


def _length_frame(value: str) -> bytes:
    encoded = str(value).encode("utf-8")
    if len(encoded) > 0xFFFFFFFF:
        raise ValidationError("Feature instance identity text is too large.")
    return len(encoded).to_bytes(4, byteorder="big") + encoded


def compute_feature_instance_hash(
    record_key: str,
    biological_feature_id: str,
) -> str:
    """Return the frozen v1 token for one canonical feature-instance pair.

    The SHA-256 input is the NUL-terminated ASCII domain followed by each UTF-8
    value framed with its unsigned 32-bit big-endian byte length.
    """

    normalized_record_key = str(record_key).strip()
    normalized_feature_id = str(biological_feature_id).strip()
    if not normalized_record_key or not normalized_feature_id:
        raise ValidationError(
            "Feature instance identity requires a non-empty record key and "
            "biological feature ID."
        )
    digest = hashlib.sha256(
        _FEATURE_INSTANCE_HASH_DOMAIN
        + _length_frame(normalized_record_key)
        + _length_frame(normalized_feature_id)
    ).digest()[:16]
    encoded = base64.b32encode(digest).decode("ascii").rstrip("=").lower()
    return f"{FEATURE_INSTANCE_HASH_PREFIX}{encoded}"


def is_feature_instance_hash(value: object) -> bool:
    """Return whether *value* is a canonical v1 feature-instance token."""

    return bool(FEATURE_INSTANCE_HASH_PATTERN.fullmatch(str(value or "")))


def validate_feature_instance_hash(value: object) -> str:
    """Return a valid token or raise a public input diagnostic."""

    token = str(value or "")
    if not is_feature_instance_hash(token):
        raise ValidationError(
            f"{FEATURE_INSTANCE_HASH_QUALIFIER} requires a case-sensitive "
            "fi1_ token with 26 lower-case base32 characters."
        )
    return token


def _iter_features(features: object):
    for feature in features or ():  # type: ignore[union-attr]
        yield feature
        yield from _iter_features(getattr(feature, "sub_features", None))


def _record_key(record: SeqRecord, index: int) -> str:
    annotations = getattr(record, "annotations", None) or {}
    value = annotations.get("gbdraw_record_key", f"record-{index + 1}")
    return "" if value is None else str(value).strip()


def _stable_feature_id(
    record: SeqRecord,
    feature: object,
    *,
    source_feature_index: int,
) -> str:
    parts = _source_feature_location_parts(feature)
    if parts is None:
        coord_base, coord_step = _read_coord_map(record)
        parts = _mapped_feature_location_parts(
            feature,
            coord_base=coord_base,
            coord_step=coord_step,
        )
    if not parts:
        return f"feature-{source_feature_index}"
    return compute_feature_hash_from_location_parts(
        str(getattr(feature, "type", "") or ""),
        parts,
        record_id=str(record.id) if record.id is not None else None,
    )


@dataclass(frozen=True)
class FeatureInstanceIdentity:
    """Canonical identity assigned to one source feature occurrence."""

    record_key: str
    record_id: str
    source_feature_index: int
    stable_feature_id: str
    biological_feature_id: str
    instance_hash: str


@dataclass(frozen=True)
class FeatureInstanceIdentityPlan:
    """Immutable identity lookup shared by rendering and metadata consumers."""

    record_keys: tuple[str, ...]
    identities: tuple[FeatureInstanceIdentity, ...]
    _by_feature_object_id: Mapping[int, FeatureInstanceIdentity] = field(
        repr=False,
        compare=False,
    )
    _record_key_by_object_id: Mapping[int, str] = field(
        repr=False,
        compare=False,
    )

    def identity_for_feature(self, feature: object) -> FeatureInstanceIdentity:
        """Resolve one feature object or fail closed on an unknown instance."""

        identity = self._by_feature_object_id.get(id(feature))
        if identity is None:
            raise ValidationError(
                "Feature instance identity context does not contain the source "
                "feature being resolved. Regenerate the diagram from one shared "
                "record plan."
            )
        return identity

    def record_key_for_record(self, record: SeqRecord) -> str:
        """Return the normalized key assigned to a materialized record object."""

        record_key = self._record_key_by_object_id.get(id(record))
        if record_key is None:
            raise ValidationError(
                "Feature instance identity context does not contain the record "
                "being resolved."
            )
        return record_key


def build_feature_instance_identity_plan(
    records: Sequence[SeqRecord] | SeqRecord,
) -> FeatureInstanceIdentityPlan:
    """Build exact identities for every top-level and nested source Feature."""

    record_list = [records] if isinstance(records, SeqRecord) else list(records)
    record_keys = tuple(
        _record_key(record, index) for index, record in enumerate(record_list)
    )
    if any(not record_key for record_key in record_keys):
        raise ValidationError("Feature instance record keys must be non-empty.")
    if len(set(record_keys)) != len(record_keys):
        raise ValidationError("Feature instance record keys must be unique.")

    candidates: list[tuple[SeqRecord, object, str, int, str]] = []
    seen_feature_objects: set[int] = set()
    for record, record_key in zip(record_list, record_keys, strict=True):
        source_indexes = _feature_source_index_map(record.features)
        used_source_indexes: set[int] = set()
        for feature in _iter_features(record.features):
            feature_object_id = id(feature)
            if feature_object_id in seen_feature_objects:
                raise ValidationError(
                    "Feature instance identity is ambiguous because one Feature "
                    "object occurs more than once."
                )
            seen_feature_objects.add(feature_object_id)
            source_feature_index = _source_feature_index(feature)
            if source_feature_index is None:
                source_feature_index = source_indexes[feature_object_id]
            if source_feature_index in used_source_indexes:
                raise ValidationError(
                    f"Feature instance source index {source_feature_index} is "
                    f"duplicated within record key {record_key!r}."
                )
            used_source_indexes.add(source_feature_index)
            stable_feature_id = _stable_feature_id(
                record,
                feature,
                source_feature_index=source_feature_index,
            )
            candidates.append(
                (
                    record,
                    feature,
                    record_key,
                    source_feature_index,
                    stable_feature_id,
                )
            )

    collision_counts = Counter(
        (record_key, stable_feature_id)
        for _record, _feature, record_key, _source_index, stable_feature_id in candidates
    )
    identities: list[FeatureInstanceIdentity] = []
    by_feature_object_id: dict[int, FeatureInstanceIdentity] = {}
    used: set[tuple[str, str]] = set()
    for record, feature, record_key, source_index, stable_feature_id in candidates:
        biological_feature_id = stable_feature_id
        if collision_counts[(record_key, stable_feature_id)] > 1:
            biological_feature_id = f"{stable_feature_id}~{source_index}"
        canonical = (record_key, biological_feature_id)
        if canonical in used:
            raise ValidationError(
                "Feature instance identity generated a duplicate canonical "
                f"feature pair {canonical!r}."
            )
        used.add(canonical)
        identity = FeatureInstanceIdentity(
            record_key=record_key,
            record_id=str(record.id or ""),
            source_feature_index=source_index,
            stable_feature_id=stable_feature_id,
            biological_feature_id=biological_feature_id,
            instance_hash=compute_feature_instance_hash(*canonical),
        )
        identities.append(identity)
        by_feature_object_id[id(feature)] = identity

    record_key_by_object_id = {
        id(record): record_key
        for record, record_key in zip(record_list, record_keys, strict=True)
    }
    return FeatureInstanceIdentityPlan(
        record_keys=record_keys,
        identities=tuple(identities),
        _by_feature_object_id=MappingProxyType(by_feature_object_id),
        _record_key_by_object_id=MappingProxyType(record_key_by_object_id),
    )


__all__ = [
    "FEATURE_INSTANCE_HASH_PATTERN",
    "FEATURE_INSTANCE_HASH_PREFIX",
    "FEATURE_INSTANCE_HASH_QUALIFIER",
    "FeatureInstanceIdentity",
    "FeatureInstanceIdentityPlan",
    "build_feature_instance_identity_plan",
    "compute_feature_instance_hash",
    "is_feature_instance_hash",
    "validate_feature_instance_hash",
]
