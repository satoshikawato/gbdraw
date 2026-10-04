"""Shared requested placement, source binding and resolved physical lane planning."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from numbers import Integral
from pathlib import Path
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from .objects import FeatureObject

from pandas import DataFrame

from gbdraw.exceptions import ValidationError
from .overrides import (
    FeatureIdentityNotice,
    FeatureOverride,
    ResolvedFeatureOverride,
    bind_feature_overrides,
    normalize_feature_overrides,
)
from .shapes import resolve_feature_rendering
from .source import (
    FeatureIdentity,
    IdentityBinding,
    SourceFeatureIdentity,
    read_identity_table,
    resolve_feature_identities,
    resolve_identity_table_rows,
)
from .visibility import should_render_feature


@dataclass(frozen=True)
class FeaturePlacementTarget:
    kind: Literal["main", "lane"]
    side: Literal["outward", "inward", "above", "below"] | None = None
    level: int | None = None

    def __post_init__(self) -> None:
        if self.kind == "main":
            if self.side is not None or self.level is not None:
                raise ValidationError("Main placement must omit side and level.")
        elif self.kind == "lane":
            if self.side not in ("outward", "inward", "above", "below"):
                raise ValidationError(
                    "Lane placement requires a supported directional side."
                )
            if (
                isinstance(self.level, bool)
                or not isinstance(self.level, Integral)
                or self.level != 1
            ):
                raise ValidationError("Only feature lane level 1 is supported.")
            object.__setattr__(self, "level", int(self.level))
        else:
            raise ValidationError(
                "Placement kind must be main or lane; Auto is override absence."
            )

    def validate_mode(self, mode: str) -> None:
        if mode not in ("circular", "linear"):
            raise ValidationError("Placement mode must be circular or linear.")
        allowed = ("outward", "inward") if mode == "circular" else ("above", "below")
        if self.kind == "lane" and self.side not in allowed:
            raise ValidationError(
                f"Placement side {self.side!r} is unsupported in {mode} mode."
            )

    @classmethod
    def from_mapping(cls, value: Mapping[str, object]) -> FeaturePlacementTarget:
        """Validate the future nested shape without activating a canonical codec."""
        if not isinstance(value, Mapping):
            raise ValidationError("Placement must be an object.")
        expected = (
            {"kind"} if value.get("kind") == "main" else {"kind", "side", "level"}
        )
        if set(value) != expected:
            raise ValidationError(
                "Unknown or missing placement fields; main omits side/level and lane requires both."
            )
        return cls(**value)


@dataclass(frozen=True)
class FeaturePlacementOverride:
    record_key: str
    biological_feature_id: str
    target: FeaturePlacementTarget

    def __post_init__(self) -> None:
        identity = self.identity
        object.__setattr__(self, "record_key", identity.record_key)
        object.__setattr__(self, "biological_feature_id", identity.biological_feature_id)
        if not isinstance(self.target, FeaturePlacementTarget):
            raise ValidationError("target must be FeaturePlacementTarget.")

    @property
    def identity(self) -> FeatureIdentity:
        return FeatureIdentity(self.record_key, self.biological_feature_id)

    @classmethod
    def from_mapping(cls, value: Mapping[str, object]) -> FeaturePlacementOverride:
        if not isinstance(value, Mapping) or set(value) != {
            "recordKey",
            "biologicalFeatureId",
            "placement",
        }:
            raise ValidationError("Unknown or missing exact placement fields.")
        return cls(
            value["recordKey"],
            value["biologicalFeatureId"],
            FeaturePlacementTarget.from_mapping(value["placement"]),
        )


def normalize_feature_placements(
    values: Sequence[FeaturePlacementOverride],
) -> tuple[FeaturePlacementOverride, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise ValidationError(
            "feature_placements must be a sequence of FeaturePlacementOverride values."
        )
    if not all(isinstance(item, FeaturePlacementOverride) for item in values):
        raise ValidationError(
            "feature_placements must contain FeaturePlacementOverride values."
        )
    if len({item.identity for item in values}) != len(values):
        raise ValidationError("Duplicate feature placement identity.")
    return tuple(
        sorted(values, key=lambda item: (item.record_key, item.biological_feature_id))
    )


@dataclass(frozen=True)
class ResolvedFeaturePlacement:
    biological_feature_id: str
    source_feature_index: int
    target: FeaturePlacementTarget
    status: Literal["foreground", "hidden", "underlay", "crop_excluded"]
    # Row of the request's canonical feature_placements; a failure names it.
    placement_index: int


@dataclass(frozen=True)
class ResolvedRecordFeatureInputs:
    """Identity-addressed feature inputs of one record instance, aligned with records."""

    record_key: str
    placements: tuple[ResolvedFeaturePlacement, ...] = ()
    # Source feature index -> that feature's bound visibility and label edits.
    overrides: Mapping[int, ResolvedFeatureOverride] = field(default_factory=dict)

    @property
    def foreground(self) -> tuple[ResolvedFeaturePlacement, ...]:
        return tuple(item for item in self.placements if item.status == "foreground")


@dataclass(frozen=True)
class RecordFeatureResolution:
    """One request's identity-addressed feature inputs, resolved once."""

    records: tuple[ResolvedRecordFeatureInputs, ...]
    notices: tuple[FeatureIdentityNotice, ...]
    # Every requested identity, including annotation targets.
    bindings: Mapping[FeatureIdentity, IdentityBinding]


_PLACEMENT_TABLE = "Feature placement table"


def _table_target(row: Mapping[str, str], mode: str) -> FeaturePlacementTarget | None:
    token, level = row["placement"].lower(), row["level"]
    if token in ("auto", "main"):
        if level:
            raise ValidationError("Auto/Main placement must omit level.")
        return None if token == "auto" else FeaturePlacementTarget("main")
    if token not in ("outward", "inward", "above", "below"):
        raise ValidationError(f"Unsupported placement token {token!r}.")
    if level not in ("", "1"):
        raise ValidationError("Only feature lane level 1 is supported.")
    target = FeaturePlacementTarget("lane", token, 1)
    target.validate_mode(mode)
    return target


def read_feature_placement_table(
    table: DataFrame | str | Path,
    *,
    mode: str,
    records: Sequence,
    record_keys: Sequence[str],
    source_record_ids: Sequence[str],
    source_catalogs: Sequence[tuple[SourceFeatureIdentity, ...]],
) -> tuple[FeaturePlacementOverride, ...]:
    """Resolve a feature placement table to exact placements; ``auto`` rows add none."""
    rows = read_identity_table(
        table, table_name=_PLACEMENT_TABLE,
        columns=("record", "feature_selector", "placement", "level"),
        required=("feature_selector", "placement"),
    )
    targets = []
    for row_number, row in enumerate(rows, start=2):
        try:
            targets.append(_table_target(row, mode))
        except (ValueError, ValidationError) as exc:
            raise ValidationError(
                f"{_PLACEMENT_TABLE} row {row_number}: {exc}",
                diagnostic={"code": "TABLE_INVALID", "row": row_number},
            ) from exc
    identities = resolve_identity_table_rows(
        [(row["record"], row["feature_selector"]) for row in rows],
        table=_PLACEMENT_TABLE,
        records=records,
        record_keys=record_keys,
        source_record_ids=source_record_ids,
        source_catalogs=source_catalogs,
    )
    return normalize_feature_placements([
        FeaturePlacementOverride(identity.record_key, identity.biological_feature_id, target)
        for identity, target in zip(identities, targets, strict=True)
        if target is not None
    ])


def resolve_record_feature_inputs(
    *,
    records: Sequence,
    record_keys: Sequence[str],
    source_catalogs: Sequence[tuple[SourceFeatureIdentity, ...]],
    placements: Sequence[FeaturePlacementOverride],
    mode: str,
    feature_overrides: Sequence[FeatureOverride] = (),
    target_identities: Iterable[FeatureIdentity] = (),
    selected_features: Sequence[str],
    feature_visibility_rules: list | None,
    specific_color_rules: Mapping,
    feature_shapes: Mapping | None,
) -> RecordFeatureResolution:
    """Bind placements, feature overrides and annotation targets in one pass.

    Tables are already exact rows (``read_feature_placement_table``,
    ``read_feature_override_table``). Placement keeps its present/dormant
    classification; an edit whose identity is not drawn becomes a notice instead
    of failing the render (Owner Q3 = A).
    """
    exact = normalize_feature_placements(placements)
    for item in exact:
        item.target.validate_mode(mode)
    rows = normalize_feature_overrides(feature_overrides)
    bindings = resolve_feature_identities(
        records=records,
        record_keys=record_keys,
        source_catalogs=source_catalogs,
        identities={
            *(item.identity for item in exact),
            *(row.identity for row in rows),
            *target_identities,
        },
    )
    overrides = bind_feature_overrides(rows, bindings, source_catalogs)
    resolved: list[list[ResolvedFeaturePlacement]] = [[] for _ in records]
    # Nested source features are never foreground placement units.
    top_level = {
        index: {id(feature) for feature in records[index].features}
        for index in {bindings[item.identity].record_index for item in exact}
    }
    for placement_index, item in enumerate(exact):
        binding = bindings[item.identity]
        if binding.status == "unresolved":
            continue
        record = records[binding.record_index]
        feature = binding.feature
        if binding.status != "present":
            status = "crop_excluded" if binding.status == "crop_excluded" else "hidden"
        elif id(feature) not in top_level[binding.record_index] or not should_render_feature(
            feature,
            selected_features,
            feature_visibility_rules=feature_visibility_rules,
            record_id=record.id,
            specific_color_rules=specific_color_rules,
            feature_override=overrides[binding.record_index].get(binding.source_feature_index),
        ):
            status = "hidden"
        elif resolve_feature_rendering(feature.type, feature_shapes) == "underlay":
            status = "underlay"
        else:
            status = "foreground"
        resolved[binding.record_index].append(
            ResolvedFeaturePlacement(
                item.biological_feature_id,
                binding.source_feature_index,
                item.target,
                status,
                placement_index,
            )
        )
    kinds: dict[FeatureIdentity, list[str]] = {}
    for item in exact:
        kinds.setdefault(item.identity, []).append("placement")
    for row in rows:
        kinds.setdefault(row.identity, []).extend(row.edits)
    notices = tuple(
        FeatureIdentityNotice(
            identity.record_key, identity.biological_feature_id,
            bindings[identity].status, tuple(edit_kinds), bindings[identity].record_index,
        )
        for identity, edit_kinds in sorted(
            kinds.items(),
            key=lambda item: (
                bindings[item[0]].record_index, item[0].biological_feature_id
            ),
        )
        if bindings[identity].status != "present"
    )
    return RecordFeatureResolution(
        tuple(
            ResolvedRecordFeatureInputs(key, tuple(items), record_overrides)
            for key, items, record_overrides in zip(record_keys, resolved, overrides, strict=True)
        ),
        notices,
        bindings,
    )


@dataclass(frozen=True)
class FeaturePlacementSlot:
    """Resolved feature-slot direction, independent of presets and resolver state."""

    mode: Literal["circular", "linear"]
    direction: Literal["split", "inside", "outside", "overlay", "above", "below"]
    separate_strands: bool

    def __post_init__(self) -> None:
        allowed = {
            "circular": ("split", "inside", "outside"),
            "linear": ("overlay", "above", "below"),
        }
        if self.mode not in allowed or self.direction not in allowed[self.mode]:
            raise ValidationError("Invalid resolved feature slot mode/direction.")

    @property
    def bidirectional(self) -> bool:
        return self.direction == "split" or (
            self.direction == "overlay" and not self.separate_strands
        )

    def supported_targets(self) -> list[dict[str, object]]:
        """Expose requested targets using the same resolved-slot validity owner."""
        targets: list[dict[str, object]] = [{"kind": "main"}]
        if self.bidirectional:
            for side in (("outward", "inward") if self.mode == "circular" else ("above", "below")):
                targets.append({"kind": "lane", "side": side, "level": 1})
        return targets

    def validate_target(
        self, target: FeaturePlacementTarget, *, record_key: str = "",
        placement: ResolvedFeaturePlacement | None = None,
    ) -> None:
        target.validate_mode(self.mode)
        if target.kind == "lane" and not self.bidirectional:
            # R6: name the feature (its request row for the Web) and the lanes
            # that accept it; Auto and Main remain valid in every slot.
            feature = (
                f" for feature {placement.biological_feature_id!r} in record {record_key!r}"
                if placement is not None else ""
            )
            remedy = (
                "a split feature slot" if self.mode == "circular"
                else "an overlay feature slot without separate strands"
            )
            raise ValidationError(
                f"Directional placement {target.side!r}{feature} is unsupported in resolved "
                f"{self.mode} slot {self.direction!r} (separate strands={self.separate_strands}); "
                f"use auto or main placement, or {remedy}.",
                diagnostic={
                    "code": "FEATURE_PLACEMENT",
                    "reason": "SPLIT_LANES" if self.mode == "circular" else "OVERLAY_LANES",
                    **({} if placement is None else {"placementIndex": placement.placement_index}),
                },
            )


@dataclass(frozen=True)
class FeaturePlacementAssignment:
    """One physical lane for every block, connector and display fragment of a feature.

    Identity remains on the runtime feature's source_feature_index and record context.
    This value is derived, never a requested target or a persistence authority.
    """

    side: Literal["main", "outward", "inward", "above", "below"]
    strand_pool: Literal["combined", "positive", "negative"]
    level: int
    requested_target: FeaturePlacementTarget | None = None


def _assignment_for_track(
    slot: FeaturePlacementSlot, track_id: int,
    target: FeaturePlacementTarget | None,
) -> FeaturePlacementAssignment:
    negative_side = "inward" if slot.mode == "circular" else "below"
    positive_side = "outward" if slot.mode == "circular" else "above"
    pool = ("negative" if track_id < 0 else "positive") if slot.separate_strands else "combined"
    level = abs(track_id) - (1 if pool == "negative" else 0)
    if slot.direction in ("inside", "below"):
        side = negative_side
    elif slot.direction in ("outside", "above"):
        side = positive_side
    elif pool != "combined":
        side = negative_side if pool == "negative" else positive_side
    elif track_id == 0:
        side = "main"
    else:
        side = negative_side if track_id < 0 else positive_side
    return FeaturePlacementAssignment(side, pool, level, target)


def placement_track_id(
    assignment: FeaturePlacementAssignment, slot: FeaturePlacementSlot,
) -> int:
    """Derive the existing allocator/renderer carrier from the resolved assignment."""
    if assignment.strand_pool == "negative":
        return -(assignment.level + 1)
    if slot.bidirectional and assignment.side in ("inward", "below"):
        return -assignment.level
    return assignment.level


def plan_feature_placements(
    feature_dict: dict[str, FeatureObject],
    *,
    slot: FeaturePlacementSlot,
    record_features: ResolvedRecordFeatureInputs | None = None,
    resolve_overlaps: bool,
    tolerance_bp: int = 0,
    genome_length: int | None = None,
) -> dict[str, FeaturePlacementAssignment]:
    """Reserve exact foreground bindings first, then run the shared Auto allocator.

    Occupancy uses local outer envelopes. Circular display projection is an isometry,
    so display fragments must not replace these metrics or create placement units.
    """
    from gbdraw.config.models.canvas import validate_feature_overlap_tolerance
    from .tracks import arrange_feature_tracks, _strand_pool

    tolerance_bp = validate_feature_overlap_tolerance(tolerance_bp)
    overrides = record_features.placements if record_features is not None else ()
    for item in overrides:
        slot.validate_target(item.target, record_key=record_features.record_key, placement=item)
    by_source = {feature.source_feature_index: key for key, feature in feature_dict.items()}
    targets = {}
    fixed = {}
    for item in sorted(overrides, key=lambda item: item.biological_feature_id):
        if item.status != "foreground":
            continue
        key = by_source.get(item.source_feature_index)
        if key is None:
            raise ValidationError(
                f"Foreground placement binding {item.biological_feature_id!r} "
                f"has no runtime source ordinal {item.source_feature_index}."
            )
        target = item.target
        if target.kind == "main":
            track_id = (
                -1 if slot.separate_strands
                and _strand_pool(feature_dict[key].strand) == "negative" else 0
            )
        else:
            track_id = -(1 + int(slot.separate_strands)) if target.side in ("inward", "below") else 1
        targets[key] = target
        fixed[key] = track_id
    arrange_feature_tracks(
        feature_dict, slot.separate_strands, resolve_overlaps,
        # Linear Auto retains its existing all-above order in overlay geometry.
        split_overlaps_by_strand=slot.mode == "circular" and slot.bidirectional,
        genome_length=genome_length, fixed_tracks=fixed, tolerance_bp=tolerance_bp,
    )
    assignments = {}
    for key, feature in feature_dict.items():
        assignment = _assignment_for_track(slot, feature.feature_track_id, targets.get(key))
        assignments[key] = assignment
        feature.placement = assignment
        feature.feature_track_id = placement_track_id(assignment, slot)
    return assignments
