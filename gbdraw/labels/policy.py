from __future__ import annotations

from typing import Literal

from gbdraw.exceptions import ValidationError


LabelRenderingPolicy = Literal["auto", "embedded_only", "external_only"]


def normalize_label_rendering(value: object) -> LabelRenderingPolicy:
    if not isinstance(value, str):
        raise ValidationError(
            "labels.rendering must be one of: auto, embedded_only, external_only"
        )
    normalized = value.strip().lower()
    if normalized not in {"auto", "embedded_only", "external_only"}:
        raise ValidationError(
            "labels.rendering must be one of: auto, embedded_only, external_only"
        )
    return normalized  # type: ignore[return-value]


__all__ = ["LabelRenderingPolicy", "normalize_label_rendering"]
