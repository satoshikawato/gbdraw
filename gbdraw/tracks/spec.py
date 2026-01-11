from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Literal

ScalarUnit = Literal["px", "factor"]


@dataclass(frozen=True)
class ScalarSpec:
    """A scalar value with an explicit unit.

    - unit="px": absolute pixels
    - unit="factor": multiply by a caller-provided reference value (e.g., canvas radius)
    """

    value: float
    unit: ScalarUnit = "factor"

    @classmethod
    def parse(cls, raw: str) -> "ScalarSpec":
        s = str(raw).strip()
        if not s:
            raise ValueError("empty scalar")

        if s.endswith("px"):
            return cls(value=float(s[:-2]), unit="px")
        if s.endswith("%"):
            return cls(value=float(s[:-1]) / 100.0, unit="factor")
        return cls(value=float(s), unit="factor")

    def resolve(self, reference: float) -> float:
        return self.value if self.unit == "px" else (self.value * reference)


@dataclass(frozen=True)
class CircularTrackPlacement:
    """Placement parameters for circular tracks (unresolved)."""

    # Common options. Resolution is done later with a reference radius.
    radius: ScalarSpec | None = None
    inner_radius: ScalarSpec | None = None
    outer_radius: ScalarSpec | None = None
    width: ScalarSpec | None = None
    z: int = 0


@dataclass(frozen=True)
class LinearTrackPlacement:
    """Placement parameters for linear tracks (unresolved)."""

    y: ScalarSpec | None = None
    height: ScalarSpec | None = None
    z: int = 0


TrackKind = Literal[
    # built-ins
    "features",
    "gc_content",
    "gc_skew",
    "ticks",
    "axis",
    "legend",
    "labels",
    "pairwise_match",
    # future / user-defined
    "custom",
]

LayoutMode = Literal["linear", "circular"]


@dataclass(frozen=True)
class TrackSpec:
    """A single track specification.

    This object is intentionally general and keeps unknown options in `params`.
    """

    id: str
    kind: TrackKind
    mode: LayoutMode
    placement: CircularTrackPlacement | LinearTrackPlacement | None = None
    show: bool = True
    params: dict[str, Any] | None = None

    def with_param(self, key: str, value: Any) -> "TrackSpec":
        merged = dict(self.params or {})
        merged[str(key)] = value
        return TrackSpec(
            id=self.id,
            kind=self.kind,
            mode=self.mode,
            placement=self.placement,
            show=self.show,
            params=merged,
        )


__all__ = [
    "CircularTrackPlacement",
    "LayoutMode",
    "LinearTrackPlacement",
    "ScalarSpec",
    "ScalarUnit",
    "TrackKind",
    "TrackSpec",
]


