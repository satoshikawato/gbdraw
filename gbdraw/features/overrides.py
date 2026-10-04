"""Per-feature visibility and label edits addressed by original-source identity."""

from __future__ import annotations

import logging
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Literal

from gbdraw.core.record_metadata import _feature_source_index_map, _source_feature_index
from gbdraw.exceptions import ValidationError

from .source import FeatureIdentity, IdentityBinding, SourceFeatureIdentity

logger = logging.getLogger(__name__)

FeatureVisibilityMode = Literal["on", "off", "exclude_matching"]
LabelVisibilityMode = Literal["on", "off"]
_ROW_FIELDS = frozenset(
    {"recordKey", "biologicalFeatureId", "featureVisibility", "labelVisibility", "labelText"}
)


def _invalid(message: str) -> ValidationError:
    return ValidationError(message, diagnostic={"code": "FEATURE_IDENTITY"})


@dataclass(frozen=True)
class FeatureOverride:
    """One feature's edits; ``None`` keeps the rule-based result for that part.

    ``label_text`` alone changes the text of a label that is shown anyway; it
    never shows a label by itself.
    """

    record_key: str
    biological_feature_id: str
    feature_visibility: FeatureVisibilityMode | None = None
    label_visibility: LabelVisibilityMode | None = None
    label_text: str | None = None

    def __post_init__(self) -> None:
        identity = self.identity
        object.__setattr__(self, "record_key", identity.record_key)
        object.__setattr__(self, "biological_feature_id", identity.biological_feature_id)
        if self.feature_visibility not in (None, "on", "off", "exclude_matching"):
            raise _invalid("feature_visibility must be on, off, exclude_matching or None.")
        if self.label_visibility not in (None, "on", "off"):
            raise _invalid("label_visibility must be on, off or None.")
        if self.label_text is not None and (
            not isinstance(self.label_text, str)
            or not self.label_text.strip()
            or any(char in self.label_text for char in "\t\r\n\0")
        ):
            raise _invalid("label_text must be one non-blank line without tabs or NUL.")
        if self.edits == ():
            raise _invalid("A feature override must set at least one edit.")

    @property
    def identity(self) -> FeatureIdentity:
        return FeatureIdentity(self.record_key, self.biological_feature_id)

    @property
    def edits(self) -> tuple[str, ...]:
        """The edit kinds this row sets, in a fixed order."""
        return tuple(
            kind
            for kind, value in (
                ("feature_visibility", self.feature_visibility),
                ("label_visibility", self.label_visibility),
                ("label_text", self.label_text),
            )
            if value is not None
        )

    @classmethod
    def from_mapping(cls, value: Mapping[str, object]) -> FeatureOverride:
        if not isinstance(value, Mapping) or set(value) != _ROW_FIELDS:
            raise _invalid("Unknown or missing feature override fields.")
        return cls(
            value["recordKey"],
            value["biologicalFeatureId"],
            value["featureVisibility"],
            value["labelVisibility"],
            value["labelText"],
        )

    def to_mapping(self) -> dict[str, str | None]:
        return {
            "recordKey": self.record_key,
            "biologicalFeatureId": self.biological_feature_id,
            "featureVisibility": self.feature_visibility,
            "labelVisibility": self.label_visibility,
            "labelText": self.label_text,
        }


def normalize_feature_overrides(
    values: Sequence[FeatureOverride],
) -> tuple[FeatureOverride, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise _invalid("feature_overrides must be a sequence of FeatureOverride values.")
    if not all(isinstance(item, FeatureOverride) for item in values):
        raise _invalid("feature_overrides must contain FeatureOverride values.")
    if len({item.identity for item in values}) != len(values):
        raise _invalid("Duplicate feature override identity.")
    return tuple(sorted(values, key=lambda item: (item.record_key, item.biological_feature_id)))


@dataclass(frozen=True)
class ResolvedFeatureOverride:
    """A request row bound to its drawn feature."""

    feature_visibility: FeatureVisibilityMode | None
    label_visibility: LabelVisibilityMode | None
    label_text: str | None
    # The label of last resort: "<type> <start>..<end>" in source coordinates.
    default_label: str


@dataclass(frozen=True)
class FeatureIdentityNotice:
    """An identity-addressed edit that this diagram does not draw (design Q4, 3.4).

    ``crop_excluded`` and ``absent`` edits stay dormant and apply again when the
    feature is drawn; ``unresolved`` edits name no feature of the source.
    """

    record_key: str
    biological_feature_id: str
    status: Literal["crop_excluded", "absent", "unresolved"]
    kinds: tuple[str, ...]
    record_index: int

    @property
    def message(self) -> str:
        return {
            "crop_excluded": (
                "The cropped record does not have the feature: it is outside the crop, or loading"
                " removed it (for example a GFF3 type filter); the edit stays dormant."
            ),
            "absent": "Loading removed the feature (for example a GFF3 type filter); the edit stays dormant.",
            "unresolved": "The source record has no such feature; the edit was not applied.",
        }[self.status]


def log_feature_identity_notices(notices: Sequence[FeatureIdentityNotice]) -> None:
    """Write one CLI log line per notice, beside the annotation warnings."""
    for notice in notices:
        logger.warning(
            "feature_identity_%s: feature %s of record #%s (%s), %s: %s",
            notice.status, notice.biological_feature_id, notice.record_index + 1,
            notice.record_key, ", ".join(notice.kinds), notice.message,
        )


def bind_feature_overrides(
    overrides: Sequence[FeatureOverride],
    bindings: Mapping[FeatureIdentity, IdentityBinding],
    source_catalogs: Sequence[tuple[SourceFeatureIdentity, ...]],
) -> tuple[dict[int, ResolvedFeatureOverride], ...]:
    """Key each present row by record index and drawn source feature index."""
    resolved: tuple[dict[int, ResolvedFeatureOverride], ...] = tuple(
        {} for _ in source_catalogs
    )
    # One pass per record catalog, not one per row.
    entries: dict[int, dict[int, SourceFeatureIdentity]] = {}
    for row in overrides:
        binding = bindings[row.identity]
        if binding.status != "present":
            continue
        if binding.record_index not in entries:
            entries[binding.record_index] = {
                item.source_feature_index: item
                for item in source_catalogs[binding.record_index]
            }
        entry = entries[binding.record_index][int(binding.source_feature_index)]
        start = min(part[0] for part in entry.location_parts)
        end = max(part[1] for part in entry.location_parts)
        resolved[binding.record_index][int(binding.source_feature_index)] = (
            ResolvedFeatureOverride(
                row.feature_visibility,
                row.label_visibility,
                row.label_text,
                f"{entry.feature_type} {start}..{end}",
            )
        )
    return resolved


def feature_override_lookup(
    record: object,
    overrides: Mapping[int, ResolvedFeatureOverride] | None,
) -> Callable[[object], ResolvedFeatureOverride | None]:
    """Return the override of each top-level or nested feature of ``record``."""
    if not overrides:
        return lambda _feature: None
    ordinals = _feature_source_index_map(getattr(record, "features", ()))

    def lookup(feature: object) -> ResolvedFeatureOverride | None:
        index = _source_feature_index(feature)
        return overrides.get(ordinals.get(id(feature)) if index is None else index)

    return lookup


__all__ = [
    "FeatureIdentityNotice",
    "FeatureOverride",
    "ResolvedFeatureOverride",
    "bind_feature_overrides",
    "feature_override_lookup",
    "log_feature_identity_notices",
    "normalize_feature_overrides",
]
