from __future__ import annotations

from typing import Literal


LabelRenderingPolicy = Literal["auto", "embedded_only", "external_only"]


def normalize_label_rendering(value: object) -> LabelRenderingPolicy:
    normalized = str(value or "auto").strip().lower()
    if normalized in {"embedded_only", "external_only"}:
        return normalized
    return "auto"


__all__ = ["LabelRenderingPolicy", "normalize_label_rendering"]
