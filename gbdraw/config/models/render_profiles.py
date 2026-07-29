from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

from gbdraw.exceptions import ValidationError
from gbdraw.labels.circular_types import CircularLabelPlacement
from gbdraw.labels.policy import LabelRenderingPolicy

from .canvas import LinearTrackLayoutMode
from .labels import CircularLabelScope, LinearLabelScope
from .root import GbdrawConfig


@dataclass(frozen=True)
class _RenderProfile:
    config: GbdrawConfig
    show_gc: bool = field(init=False)
    show_skew: bool = field(init=False)
    show_depth: bool = field(init=False)
    strandedness: bool = field(init=False)
    resolve_overlaps: bool = field(init=False)
    label_rendering: LabelRenderingPolicy = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.config, GbdrawConfig):
            raise ValidationError("render profile config must be GbdrawConfig")
        canvas = self.config.canvas
        object.__setattr__(self, "show_gc", bool(canvas.show_gc))
        object.__setattr__(self, "show_skew", bool(canvas.show_skew))
        object.__setattr__(self, "show_depth", bool(canvas.show_depth))
        object.__setattr__(self, "strandedness", bool(canvas.strandedness))
        object.__setattr__(self, "resolve_overlaps", bool(canvas.resolve_overlaps))
        object.__setattr__(self, "label_rendering", self.config.labels.rendering)


@dataclass(frozen=True)
class CircularRenderProfile(_RenderProfile):
    """Resolved circular render settings shared by downstream constructors."""

    label_scope: CircularLabelScope = field(init=False)
    label_placement: CircularLabelPlacement = field(init=False)
    labels_enabled: bool = field(init=False)
    inner_labels_enabled: bool = field(init=False)

    def __post_init__(self) -> None:
        super().__post_init__()
        scope = self.config.labels.circular.scope
        object.__setattr__(self, "label_scope", scope)
        object.__setattr__(self, "label_placement", self.config.labels.circular.placement)
        object.__setattr__(self, "labels_enabled", scope != "none")
        object.__setattr__(self, "inner_labels_enabled", scope == "both")


@dataclass(frozen=True)
class LinearRenderProfile(_RenderProfile):
    """Resolved linear render settings shared by downstream constructors."""

    label_scope: LinearLabelScope = field(init=False)
    label_placement: Literal["auto", "above_feature"] = field(init=False)
    label_rotation: float = field(init=False)
    track_layout: LinearTrackLayoutMode = field(init=False)
    ruler_on_axis: bool = field(init=False)

    def __post_init__(self) -> None:
        super().__post_init__()
        labels = self.config.labels.linear
        canvas = self.config.canvas.linear
        object.__setattr__(self, "label_scope", labels.scope)
        object.__setattr__(self, "label_placement", labels.placement)
        object.__setattr__(self, "label_rotation", labels.rotation)
        object.__setattr__(self, "track_layout", canvas.track_layout)
        object.__setattr__(self, "ruler_on_axis", bool(canvas.ruler_on_axis))


RenderProfile = CircularRenderProfile | LinearRenderProfile

__all__ = [
    "CircularRenderProfile",
    "LinearRenderProfile",
    "RenderProfile",
]
