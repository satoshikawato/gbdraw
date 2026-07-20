#!/usr/bin/env python
# coding: utf-8

"""Feature rendering helpers.

Public option names retain ``feature_shape`` for compatibility, while this module
models the value as a rendering policy because ``underlay`` is a layer choice.
"""

from __future__ import annotations

from typing import Iterable, Literal, Mapping, cast


FeatureRendering = Literal["arrow", "rectangle", "underlay"]
FeatureShape = FeatureRendering

FEATURE_RENDERING_VALUES: frozenset[str] = frozenset(
    {"arrow", "rectangle", "underlay"}
)

DEFAULT_FEATURE_RENDERINGS: dict[str, FeatureRendering] = {
    "CDS": "arrow",
    "rRNA": "arrow",
    "tRNA": "arrow",
    "tmRNA": "arrow",
    "ncRNA": "arrow",
    "misc_RNA": "arrow",
    "repeat_region": "underlay",
}

DEFAULT_DIRECTIONAL_FEATURE_TYPES: frozenset[str] = frozenset(
    {"CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA"}
)


def normalize_feature_shape(value: str) -> FeatureShape:
    normalized = str(value).strip().lower()
    if normalized not in FEATURE_RENDERING_VALUES:
        raise ValueError(
            f"invalid feature shape '{value}': expected 'arrow', 'rectangle', or 'underlay'"
        )
    return cast(FeatureShape, normalized)


def default_feature_rendering(feature_type: str) -> FeatureRendering:
    """Return the fresh-schema rendering default for one feature type."""

    return DEFAULT_FEATURE_RENDERINGS.get(str(feature_type), "rectangle")


def resolve_feature_rendering(
    feature_type: str,
    feature_shapes: Mapping[str, str] | None,
) -> FeatureRendering:
    """Resolve one feature type against strict explicit overrides and defaults."""

    normalized_shapes = normalize_feature_shape_overrides(feature_shapes)
    return normalized_shapes.get(
        str(feature_type),
        default_feature_rendering(str(feature_type)),
    )


def parse_feature_shape_assignment(raw: str) -> tuple[str, FeatureShape]:
    text = str(raw).strip()
    if "=" not in text:
        raise ValueError(
            "feature shape must be in TYPE=SHAPE format (e.g. CDS=arrow)"
        )
    feature_type_raw, shape_raw = text.split("=", 1)
    feature_type = feature_type_raw.strip()
    if not feature_type:
        raise ValueError("feature type in TYPE=SHAPE must not be empty")
    shape = normalize_feature_shape(shape_raw)
    return feature_type, shape


def parse_feature_shape_overrides(
    assignments: Iterable[str] | None,
) -> dict[str, FeatureShape]:
    overrides: dict[str, FeatureShape] = {}
    if assignments is None:
        return overrides
    for raw in assignments:
        feature_type, shape = parse_feature_shape_assignment(raw)
        overrides[feature_type] = shape
    return overrides


def normalize_feature_shape_overrides(
    feature_shapes: Mapping[str, str] | None,
) -> dict[str, FeatureShape]:
    overrides: dict[str, FeatureShape] = {}
    if feature_shapes is None:
        return overrides
    for feature_type_raw, shape_raw in feature_shapes.items():
        feature_type = str(feature_type_raw).strip()
        if not feature_type:
            raise ValueError("feature type in feature_shapes must not be empty")
        overrides[feature_type] = normalize_feature_shape(str(shape_raw))
    return overrides


def resolve_directional_feature_types(
    feature_shapes: Mapping[str, str] | None,
    *,
    base_directional_types: Iterable[str] = DEFAULT_DIRECTIONAL_FEATURE_TYPES,
) -> set[str]:
    directional_types = {str(feature_type) for feature_type in base_directional_types}
    normalized_shapes = normalize_feature_shape_overrides(feature_shapes)
    for feature_type, shape in normalized_shapes.items():
        if shape == "arrow":
            directional_types.add(feature_type)
        else:
            directional_types.discard(feature_type)
    return directional_types


def resolve_underlay_feature_types(
    feature_shapes: Mapping[str, str] | None,
) -> set[str]:
    """Return default and explicitly configured underlay feature types."""

    underlay_types = {
        feature_type
        for feature_type, rendering in DEFAULT_FEATURE_RENDERINGS.items()
        if rendering == "underlay"
    }
    for feature_type, rendering in normalize_feature_shape_overrides(feature_shapes).items():
        if rendering == "underlay":
            underlay_types.add(feature_type)
        else:
            underlay_types.discard(feature_type)
    return underlay_types


__all__ = [
    "DEFAULT_DIRECTIONAL_FEATURE_TYPES",
    "DEFAULT_FEATURE_RENDERINGS",
    "FEATURE_RENDERING_VALUES",
    "FeatureRendering",
    "FeatureShape",
    "default_feature_rendering",
    "normalize_feature_shape",
    "normalize_feature_shape_overrides",
    "parse_feature_shape_assignment",
    "parse_feature_shape_overrides",
    "resolve_directional_feature_types",
    "resolve_feature_rendering",
    "resolve_underlay_feature_types",
]
