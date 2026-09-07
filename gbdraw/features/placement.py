"""Shared requested placement, source binding and resolved physical lane planning."""

from __future__ import annotations

import csv
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from numbers import Integral, Real
from pathlib import Path
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from .objects import FeatureObject

from pandas import DataFrame, isna

from gbdraw.annotations.models import parse_feature_selector
from gbdraw.core.record_metadata import (
    _feature_source_index_map,
    _iter_source_features,
    _source_feature_index,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import parse_record_selector, select_record
from .shapes import resolve_feature_rendering
from .source import SourceFeatureIdentity
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
        for name in ("record_key", "biological_feature_id"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip() or "\0" in value:
                raise ValidationError(
                    f"{name} must be a non-empty identity without NUL."
                )
            object.__setattr__(self, name, value.strip())
        if not isinstance(self.target, FeaturePlacementTarget):
            raise ValidationError("target must be FeaturePlacementTarget.")

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
    if len({(item.record_key, item.biological_feature_id) for item in values}) != len(
        values
    ):
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


@dataclass(frozen=True)
class ResolvedPlacementInputs:
    """Source bindings for one record instance, aligned with records/provenance."""

    record_key: str
    overrides: tuple[ResolvedFeaturePlacement, ...] = ()

    @property
    def foreground(self) -> tuple[ResolvedFeaturePlacement, ...]:
        return tuple(item for item in self.overrides if item.status == "foreground")


def _placement_table_rows(table: DataFrame | str | Path) -> list[dict]:
    if isinstance(table, DataFrame):
        columns = list(table.columns)
        rows = table.to_dict("records")
    else:
        try:
            with open(table, encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                columns = reader.fieldnames or []
                rows = list(reader)
        except (OSError, csv.Error) as exc:
            raise ValidationError(
                f"Cannot read feature placement table: {exc}"
            ) from exc
    if len(set(columns)) != len(columns) or set(columns) - {
        "record",
        "feature_selector",
        "placement",
        "level",
    }:
        raise ValidationError("Duplicate or unknown feature placement table columns.")
    if not {"feature_selector", "placement"} <= set(columns):
        raise ValidationError(
            "Feature placement table requires feature_selector and placement columns."
        )
    if any(None in row for row in rows):
        raise ValidationError(
            "Feature placement table row has more values than columns."
        )
    return rows


def _cell(value: object) -> str:
    if isinstance(value, (Mapping, list, tuple, set)):
        raise ValidationError("Feature placement table cells must be scalar values.")
    if value is None or (not isinstance(value, str) and isna(value)):
        return ""
    return str(value).strip()


def _table_target(row: dict, mode: str) -> FeaturePlacementTarget | None:
    token, level = _cell(row.get("placement")).lower(), _cell(row.get("level"))
    raw_level = row.get("level")
    if isinstance(raw_level, Real) and not isinstance(raw_level, bool) and raw_level == 1:
        level = "1"  # pandas promotes an integer column with blank cells to float.
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


def resolve_placement_inputs(
    *,
    records: Sequence,
    record_keys: Sequence[str],
    source_record_ids: Sequence[str],
    source_catalogs: Sequence[tuple[SourceFeatureIdentity, ...]],
    overrides: Sequence[FeaturePlacementOverride],
    mode: str,
    table: DataFrame | str | Path | None = None,
    selected_features: Sequence[str],
    feature_visibility_rules: list | None,
    specific_color_rules: Mapping,
    feature_shapes: Mapping | None,
) -> tuple[tuple[FeaturePlacementOverride, ...], tuple[ResolvedPlacementInputs, ...]]:
    """Materialize exactly-one source selectors and classify present/dormant intent."""
    if not (
        len(records)
        == len(record_keys)
        == len(source_record_ids)
        == len(source_catalogs)
    ):
        raise ValidationError(
            "Placement source context must align with records/provenance."
        )
    if len(set(record_keys)) != len(record_keys):
        raise ValidationError("Placement context contains duplicate record keys.")
    exact = normalize_feature_placements(overrides)
    if table is not None:
        if exact:
            raise ValidationError(
                "Exact placements and placement table inputs are mutually exclusive."
            )
        rows, seen = [], set()
        for row_number, row in enumerate(_placement_table_rows(table), start=2):
            try:
                target = _table_target(row, mode)
                selector = parse_record_selector(_cell(row.get("record")))
                selected = select_record(records, selector)
                if len(selected) != 1:
                    raise ValidationError(
                        "Placement record selector must match exactly one record; use #index."
                    )
                index = next(
                    index
                    for index, record in enumerate(records)
                    if record is selected[0]
                )
                feature_selector = parse_feature_selector(
                    _cell(row.get("feature_selector"))
                )
                matched = [
                    entry
                    for entry in source_catalogs[index]
                    if entry.matches(
                        key=feature_selector.key,
                        value=feature_selector.value,
                        record_id=source_record_ids[index],
                    )
                ]
                if len(matched) != 1:
                    raise ValidationError(
                        f"Feature placement selector must match exactly one source feature; matched {len(matched)}."
                    )
                identity = record_keys[index], matched[0].biological_feature_id
                if identity in seen:
                    raise ValidationError(
                        "Duplicate resolved feature placement identity."
                    )
                seen.add(identity)
                if target is not None:
                    rows.append(FeaturePlacementOverride(*identity, target))
            except (ValueError, ValidationError) as exc:
                raise ValidationError(
                    f"Feature placement table row {row_number}: {exc}"
                ) from exc
        exact = normalize_feature_placements(rows)
    by_key = {key: [] for key in record_keys}
    for item in exact:
        item.target.validate_mode(mode)
        if item.record_key not in by_key:
            raise ValidationError(f"Unknown placement record key {item.record_key!r}.")
        by_key[item.record_key].append(item)
    aligned = []
    for record, key, catalog in zip(records, record_keys, source_catalogs, strict=True):
        known = {entry.biological_feature_id: entry for entry in catalog}
        runtime = {}
        indexes = _feature_source_index_map(record.features)
        for feature in _iter_source_features(record.features):
            index = _source_feature_index(feature)
            runtime[indexes[id(feature)] if index is None else index] = feature
        foreground_candidates = {id(feature) for feature in record.features}
        resolved = []
        for item in by_key[key]:
            source = known.get(item.biological_feature_id)
            if source is None:
                raise ValidationError(
                    f"Unknown/stale placement identity ({key!r}, {item.biological_feature_id!r})."
                )
            feature = runtime.get(source.source_feature_index)
            if feature is None:
                status = (
                    "crop_excluded"
                    if record.annotations.get("gbdraw_region_applied")
                    else "hidden"
                )
            elif id(feature) not in foreground_candidates or not should_render_feature(
                feature,
                selected_features,
                feature_visibility_rules=feature_visibility_rules,
                record_id=record.id,
                specific_color_rules=specific_color_rules,
            ):
                status = "hidden"
            elif resolve_feature_rendering(feature.type, feature_shapes) == "underlay":
                status = "underlay"
            else:
                status = "foreground"
            resolved.append(
                ResolvedFeaturePlacement(
                    item.biological_feature_id,
                    source.source_feature_index,
                    item.target,
                    status,
                )
            )
        aligned.append(ResolvedPlacementInputs(key, tuple(resolved)))
    return exact, tuple(aligned)


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
        return not self.separate_strands and self.direction in ("split", "overlay")

    def supported_targets(self) -> list[dict[str, object]]:
        """Expose requested targets using the same resolved-slot validity owner."""
        targets: list[dict[str, object]] = [{"kind": "main"}]
        if self.bidirectional:
            for side in (("outward", "inward") if self.mode == "circular" else ("above", "below")):
                targets.append({"kind": "lane", "side": side, "level": 1})
        return targets

    def validate_target(self, target: FeaturePlacementTarget) -> None:
        target.validate_mode(self.mode)
        if target.kind == "lane" and not self.bidirectional:
            raise ValidationError(
                f"Directional placement is unsupported in resolved {self.mode} "
                f"slot {self.direction!r} (separate strands={self.separate_strands})."
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
    placement_inputs: ResolvedPlacementInputs | None = None,
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
    overrides = placement_inputs.overrides if placement_inputs is not None else ()
    for item in overrides:
        slot.validate_target(item.target)
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
            track_id = -1 if target.side in ("inward", "below") else 1
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
