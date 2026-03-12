from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Pattern, Sequence

from ..exceptions import ValidationError
from .spec import CircularTrackPlacement, TrackSpec


_TRACK_FILE_ALLOWED_KINDS = {"features", "gc_content", "gc_skew", "analysis", "custom"}
_CIRCULAR_SINGLETON_KINDS = {"features", "gc_content", "gc_skew"}
_CIRCULAR_CUSTOM_STRAND_MODES = {"all", "positive", "negative"}
_CIRCULAR_ANALYSIS_METRICS = {"content", "skew"}


@dataclass(frozen=True)
class CircularCustomRule:
    qualifier: str
    pattern: str
    regex: Pattern[str]


def _normalize_track_id(value: Any) -> str:
    track_id = str(value or "").strip()
    if not track_id:
        raise ValidationError("Circular track id must be a non-empty string.")
    return track_id


def _normalize_kind(value: Any) -> str:
    kind = str(value or "").strip().lower()
    if not kind:
        raise ValidationError("Circular track kind must be a non-empty string.")
    return kind


def _normalize_bool(value: Any, *, field_name: str) -> bool:
    if isinstance(value, bool):
        return value
    raise ValidationError(f"Circular track '{field_name}' must be a boolean.")


def _normalize_feature_types(value: Any, *, field_name: str, allow_missing: bool) -> list[str] | None:
    if value is None:
        return None if allow_missing else []
    if not isinstance(value, list):
        raise ValidationError(f"Circular track '{field_name}' must be an array of strings.")

    normalized: list[str] = []
    seen: set[str] = set()
    for item in value:
        feature_type = str(item or "").strip()
        if not feature_type:
            raise ValidationError(f"Circular track '{field_name}' entries must be non-empty strings.")
        if feature_type in seen:
            continue
        seen.add(feature_type)
        normalized.append(feature_type)
    return normalized


def _parse_placement_mapping(value: Any) -> CircularTrackPlacement | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise ValidationError("Circular track placement must be an object.")

    allowed_fields = {"radius", "inner_radius", "outer_radius", "width", "z"}
    unknown_fields = sorted(str(key) for key in value.keys() if str(key) not in allowed_fields)
    if unknown_fields:
        raise ValidationError(
            f"Unsupported circular track placement fields: {', '.join(unknown_fields)}."
        )

    placement_kwargs: dict[str, Any] = {}
    radius_value = value.get("radius")
    inner_radius_value = value.get("inner_radius")
    outer_radius_value = value.get("outer_radius")
    width_value = value.get("width")
    z_value = value.get("z", 0)

    try:
        from .spec import ScalarSpec

        if radius_value is not None:
            placement_kwargs["radius"] = ScalarSpec.parse(str(radius_value))
        if inner_radius_value is not None:
            placement_kwargs["inner_radius"] = ScalarSpec.parse(str(inner_radius_value))
        if outer_radius_value is not None:
            placement_kwargs["outer_radius"] = ScalarSpec.parse(str(outer_radius_value))
        if width_value is not None:
            placement_kwargs["width"] = ScalarSpec.parse(str(width_value))
        placement_kwargs["z"] = int(z_value)
    except ValueError as exc:
        raise ValidationError(f"Invalid circular track placement value: {exc}") from exc
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"Invalid circular track placement z value: {z_value!r}") from exc

    return CircularTrackPlacement(**placement_kwargs)


def _normalize_custom_rules(value: Any) -> list[dict[str, str]]:
    if value is None:
        return []
    if not isinstance(value, list):
        raise ValidationError("Custom circular track rules must be an array.")

    normalized_rules: list[dict[str, str]] = []
    for idx, rule in enumerate(value):
        if not isinstance(rule, Mapping):
            raise ValidationError(f"Custom circular track rule #{idx + 1} must be an object.")
        qualifier = str(rule.get("qualifier") or "").strip()
        pattern = str(rule.get("pattern") or "").strip()
        if not qualifier:
            raise ValidationError(f"Custom circular track rule #{idx + 1} qualifier must be non-empty.")
        try:
            re.compile(pattern, re.IGNORECASE)
        except re.error as exc:
            raise ValidationError(
                f"Invalid regex in custom circular track rule #{idx + 1}: {pattern!r} ({exc})"
            ) from exc
        normalized_rules.append({"qualifier": qualifier, "pattern": pattern})
    return normalized_rules


def _normalize_dinucleotide(value: Any, *, field_name: str) -> str:
    dinucleotide = str(value or "").strip().upper()
    if len(dinucleotide) != 2 or not dinucleotide.isalpha():
        raise ValidationError(f"Circular track {field_name} must be a two-letter nucleotide pair.")
    return dinucleotide


def _normalize_params_for_kind(kind: str, params: Any) -> dict[str, Any] | None:
    if params is None:
        return None
    if not isinstance(params, Mapping):
        raise ValidationError("Circular track params must be an object.")

    params_dict = dict(params)
    if kind == "custom":
        caption = str(params_dict.get("caption") or "").strip()
        if not caption:
            raise ValidationError("Custom circular track params.caption must be a non-empty string.")
        feature_types = _normalize_feature_types(
            params_dict.get("feature_types"),
            field_name="params.feature_types",
            allow_missing=False,
        )
        if not feature_types:
            raise ValidationError("Custom circular track params.feature_types must be non-empty.")
        strand_mode = str(params_dict.get("strand_mode") or "").strip().lower() or "all"
        if strand_mode not in _CIRCULAR_CUSTOM_STRAND_MODES:
            raise ValidationError(
                "Custom circular track params.strand_mode must be one of: all, positive, negative."
            )
        normalized_rules = _normalize_custom_rules(params_dict.get("rules"))
        return {
            "caption": caption,
            "feature_types": feature_types,
            "strand_mode": strand_mode,
            "rules": normalized_rules,
            "match_all": True,
        }

    if kind == "analysis":
        metric = str(params_dict.get("metric") or "").strip().lower()
        if metric not in _CIRCULAR_ANALYSIS_METRICS:
            raise ValidationError("Analysis circular track params.metric must be one of: content, skew.")
        dinucleotide = _normalize_dinucleotide(
            params_dict.get("dinucleotide"),
            field_name="params.dinucleotide",
        )
        caption = str(params_dict.get("caption") or "").strip()
        return {
            "caption": caption,
            "metric": metric,
            "dinucleotide": dinucleotide,
        }

    if kind == "features":
        feature_types = _normalize_feature_types(
            params_dict.get("feature_types"),
            field_name="params.feature_types",
            allow_missing=True,
        )
        normalized: dict[str, Any] = {}
        if feature_types is not None:
            if not feature_types:
                raise ValidationError("Circular features track params.feature_types must be non-empty when provided.")
            normalized["feature_types"] = feature_types
        for key, value in params_dict.items():
            if key == "feature_types":
                continue
            normalized[str(key)] = value
        return normalized or None

    return dict(params_dict) or None


def normalize_circular_track_specs(track_specs: Iterable[TrackSpec]) -> list[TrackSpec]:
    normalized_specs: list[TrackSpec] = []
    seen_ids: set[str] = set()
    seen_singletons: dict[str, str] = {}

    for track_spec in track_specs:
        if track_spec.mode != "circular":
            raise ValidationError(
                f"TrackSpec mode '{track_spec.mode}' is not supported for circular diagrams."
            )

        track_id = _normalize_track_id(track_spec.id)
        if track_id in seen_ids:
            raise ValidationError(f"Duplicate circular track id: {track_id!r}.")
        seen_ids.add(track_id)

        kind = _normalize_kind(track_spec.kind)
        placement = track_spec.placement
        params = _normalize_params_for_kind(kind, track_spec.params)

        if kind in _CIRCULAR_SINGLETON_KINDS:
            duplicate_id = seen_singletons.get(kind)
            if duplicate_id is not None:
                raise ValidationError(
                    f"Duplicate circular built-in track kind '{kind}' for ids {duplicate_id!r} and {track_id!r}."
                )
            seen_singletons[kind] = track_id

        normalized_specs.append(
            TrackSpec(
                id=track_id,
                kind=kind,  # type: ignore[arg-type]
                mode="circular",
                placement=placement,
                show=bool(track_spec.show),
                params=params,
            )
        )

    return normalized_specs


def load_circular_track_specs(track_file: str | Path) -> list[TrackSpec]:
    track_path = Path(track_file)
    try:
        payload = json.loads(track_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValidationError(f"Invalid JSON in circular track file '{track_path}': {exc}") from exc

    if not isinstance(payload, list):
        raise ValidationError("Circular track file payload must be an ordered JSON array.")

    parsed_specs: list[TrackSpec] = []
    for idx, entry in enumerate(payload):
        if not isinstance(entry, Mapping):
            raise ValidationError(f"Circular track entry #{idx + 1} must be an object.")

        track_id = _normalize_track_id(entry.get("id"))
        kind = _normalize_kind(entry.get("kind"))
        if kind not in _TRACK_FILE_ALLOWED_KINDS:
            raise ValidationError(
                f"Circular track file entry #{idx + 1} has unsupported kind {kind!r}."
            )
        show = _normalize_bool(entry.get("show"), field_name="show")
        placement = _parse_placement_mapping(entry.get("placement"))
        params = _normalize_params_for_kind(kind, entry.get("params"))

        parsed_specs.append(
            TrackSpec(
                id=track_id,
                kind=kind,  # type: ignore[arg-type]
                mode="circular",
                placement=placement,
                show=show,
                params=params,
            )
        )

    return normalize_circular_track_specs(parsed_specs)


def iter_compiled_custom_rules(track_spec: TrackSpec) -> tuple[CircularCustomRule, ...]:
    if str(track_spec.kind) != "custom":
        return tuple()
    params = dict(track_spec.params or {})
    compiled_rules: list[CircularCustomRule] = []
    for rule in list(params.get("rules") or []):
        qualifier = str(dict(rule).get("qualifier") or "").strip()
        pattern = str(dict(rule).get("pattern") or "").strip()
        if not qualifier:
            continue
        compiled_rules.append(
            CircularCustomRule(
                qualifier=qualifier,
                pattern=pattern,
                regex=re.compile(pattern, re.IGNORECASE),
            )
        )
    return tuple(compiled_rules)


def get_circular_feature_type_union(
    track_specs: Sequence[TrackSpec] | None,
    *,
    fallback_feature_types: Sequence[str],
) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    saw_features_row = False
    if track_specs:
        for track_spec in track_specs:
            kind = str(track_spec.kind)
            params = dict(track_spec.params or {})
            if kind == "features":
                saw_features_row = True
            if not track_spec.show:
                continue
            if kind == "custom":
                feature_types = list(params.get("feature_types") or [])
            elif kind == "features":
                feature_types = list(params.get("feature_types") or fallback_feature_types)
            else:
                continue
            for feature_type in feature_types:
                normalized = str(feature_type or "").strip()
                if not normalized or normalized in seen:
                    continue
                seen.add(normalized)
                ordered.append(normalized)
    if not saw_features_row:
        for feature_type in fallback_feature_types:
            normalized = str(feature_type or "").strip()
            if not normalized or normalized in seen:
                continue
            seen.add(normalized)
            ordered.append(normalized)
    return ordered


__all__ = [
    "CircularCustomRule",
    "get_circular_feature_type_union",
    "iter_compiled_custom_rules",
    "load_circular_track_specs",
    "normalize_circular_track_specs",
]
