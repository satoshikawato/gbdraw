from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

from gbdraw.exceptions import ValidationError
from gbdraw.labels.circular_types import CircularLabelPlacement
from gbdraw.labels.filtering import has_forced_label_overrides
from gbdraw.labels.label_style import FeatureLabelStyle
from gbdraw.labels.policy import LabelRenderingPolicy

from .canvas import LinearTrackLayoutMode
from .labels import LINEAR_LABEL_SCOPE_RECORDS, CircularLabelScope, LinearLabelScope
from .root import GbdrawConfig


@dataclass(frozen=True)
class _RenderProfile:
    config: GbdrawConfig
    # An identity-addressed feature override shows a label (label_visibility "on").
    identity_label_overrides: bool = False
    # A per-feature override shows its label outside the label display scope.
    forced_labels: bool = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.config, GbdrawConfig):
            raise ValidationError("render profile config must be GbdrawConfig")
        object.__setattr__(
            self,
            "forced_labels",
            self.identity_label_overrides
            or has_forced_label_overrides(self.config.labels.filtering.as_dict()),
        )

    # A render stage keys a profile field by value, whether it reads it or not,
    # and the config by the sections and leaves it declares (gbdraw.render_stages).
    # The settings below are read through the config, so changing one rebuilds
    # only the stages that declare it (``canvas.strandedness``, or all of
    # ``canvas``), not every content that takes the profile (OV-1309).
    # ``forced_labels`` stays a field: it is derived from all of
    # ``labels.filtering``, which a stage that reads it would have to declare.
    @property
    def strandedness(self) -> bool:
        return bool(self.config.canvas.strandedness)

    @property
    def resolve_overlaps(self) -> bool:
        return bool(self.config.canvas.resolve_overlaps)

    @property
    def feature_overlap_tolerance_bp(self) -> int:
        return self.config.canvas.feature_overlap_tolerance_bp

    @property
    def label_rendering(self) -> LabelRenderingPolicy:
        return self.config.labels.rendering

    @property
    def show_gc(self) -> bool:
        return bool(self.config.canvas.show_gc)

    @property
    def show_skew(self) -> bool:
        return bool(self.config.canvas.show_skew)

    @property
    def show_depth(self) -> bool:
        return bool(self.config.canvas.show_depth)


@dataclass(frozen=True)
class CircularRenderProfile(_RenderProfile):
    """Resolved circular render settings shared by downstream constructors."""

    label_scope: CircularLabelScope = field(init=False)
    label_placement: CircularLabelPlacement = field(init=False)
    labels_enabled: bool = field(init=False)
    # Only per-feature overrides label features when the scope is "none".
    labels_overrides_only: bool = field(init=False)
    inner_labels_enabled: bool = field(init=False)
    label_style: FeatureLabelStyle = field(init=False)

    def __post_init__(self) -> None:
        super().__post_init__()
        labels = self.config.labels.circular
        scope = labels.scope
        object.__setattr__(self, "label_scope", scope)
        object.__setattr__(self, "label_placement", self.config.labels.circular.placement)
        object.__setattr__(self, "labels_enabled", scope != "none" or self.forced_labels)
        object.__setattr__(self, "labels_overrides_only", scope == "none")
        object.__setattr__(self, "inner_labels_enabled", scope == "both")
        object.__setattr__(self, "label_style", FeatureLabelStyle(labels.line, labels.protection))


@dataclass(frozen=True)
class LinearRenderProfile(_RenderProfile):
    """Resolved linear render settings shared by downstream constructors."""

    label_scope: LinearLabelScope = field(init=False)
    label_placement: Literal["auto", "above_feature"] = field(init=False)
    label_rotation: float = field(init=False)
    track_layout: LinearTrackLayoutMode = field(init=False)
    ruler_on_axis: bool = field(init=False)
    label_style: FeatureLabelStyle = field(init=False)

    def __post_init__(self) -> None:
        super().__post_init__()
        labels = self.config.labels.linear
        canvas = self.config.canvas.linear
        object.__setattr__(self, "label_scope", labels.scope)
        object.__setattr__(self, "label_placement", labels.placement)
        object.__setattr__(self, "label_rotation", labels.rotation)
        object.__setattr__(self, "label_style", FeatureLabelStyle(labels.line, labels.protection))
        object.__setattr__(self, "track_layout", canvas.track_layout)
        object.__setattr__(self, "ruler_on_axis", bool(canvas.ruler_on_axis))

    def labels_in_scope(self, record_index: int) -> bool:
        """Return whether the label display scope selects this record's labels."""
        selected = LINEAR_LABEL_SCOPE_RECORDS[self.label_scope]
        return selected == "all" or (selected == "first" and record_index == 0)

    def record_has_labels(self, record_index: int) -> bool:
        """Return whether this record can draw labels, in scope or by override."""
        return self.labels_in_scope(record_index) or self.forced_labels


RenderProfile = CircularRenderProfile | LinearRenderProfile

__all__ = [
    "CircularRenderProfile",
    "LinearRenderProfile",
    "RenderProfile",
]
