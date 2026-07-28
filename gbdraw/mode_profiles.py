"""Versioned defaults and shared validation for Circular and Linear modes."""

from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral, Real
from typing import Literal, Mapping, cast

from gbdraw.exceptions import ValidationError


DiagramMode = Literal["circular", "linear"]
MODE_PROFILE_VERSION = 1
DEFAULT_FEATURE_TYPES = (
    "CDS",
    "rRNA",
    "tRNA",
    "tmRNA",
    "ncRNA",
    "misc_RNA",
    "repeat_region",
)


def _finite_nonnegative(value: object, *, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValidationError(f"{field_name} must be a finite number >= 0.")
    try:
        normalized = float(value)
    except (OverflowError, TypeError, ValueError) as exc:
        raise ValidationError(
            f"{field_name} must be a finite number >= 0."
        ) from exc
    if not math.isfinite(normalized) or normalized < 0:
        raise ValidationError(f"{field_name} must be a finite number >= 0.")
    return normalized


def _identity_percent(value: object) -> float:
    normalized = _finite_nonnegative(value, field_name="identity")
    if normalized > 100:
        raise ValidationError("identity must be a finite number in [0, 100].")
    return normalized


def _alignment_length(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValidationError("alignment_length must be >= 0 and an integer.")
    normalized = int(value)
    if normalized < 0:
        raise ValidationError("alignment_length must be >= 0 and an integer.")
    return normalized


@dataclass(frozen=True)
class ComparisonThresholds:
    """Validated comparison and conservation filter thresholds."""

    evalue: float
    bitscore: float
    identity: float
    alignment_length: int = 0

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "evalue",
            _finite_nonnegative(self.evalue, field_name="evalue"),
        )
        object.__setattr__(
            self,
            "bitscore",
            _finite_nonnegative(self.bitscore, field_name="bitscore"),
        )
        object.__setattr__(self, "identity", _identity_percent(self.identity))
        object.__setattr__(
            self,
            "alignment_length",
            _alignment_length(self.alignment_length),
        )


@dataclass(frozen=True)
class ModeProfile:
    """Resolved defaults for one fresh drawing request."""

    mode: DiagramMode
    comparison: ComparisonThresholds
    show_gc: bool
    show_skew: bool
    feature_types: tuple[str, ...]
    linear_axis_color: str | None
    linear_ruler_axis_color: str | None

    @property
    def config_overrides(self) -> dict[str, object]:
        overrides: dict[str, object] = {
            "show_gc": self.show_gc,
            "show_skew": self.show_skew,
        }
        if self.linear_axis_color is not None:
            overrides["linear_axis_stroke_color"] = self.linear_axis_color
        return overrides


CIRCULAR_MODE_PROFILE = ModeProfile(
    mode="circular",
    comparison=ComparisonThresholds(
        evalue=1e-5,
        bitscore=50.0,
        identity=70.0,
        alignment_length=0,
    ),
    show_gc=True,
    show_skew=True,
    feature_types=DEFAULT_FEATURE_TYPES,
    linear_axis_color=None,
    linear_ruler_axis_color=None,
)

LINEAR_MODE_PROFILE = ModeProfile(
    mode="linear",
    comparison=ComparisonThresholds(
        evalue=1e-2,
        bitscore=50.0,
        identity=0.0,
        alignment_length=0,
    ),
    show_gc=False,
    show_skew=False,
    feature_types=DEFAULT_FEATURE_TYPES,
    linear_axis_color="lightgray",
    linear_ruler_axis_color="dimgray",
)

MODE_PROFILES: dict[DiagramMode, ModeProfile] = {
    "circular": CIRCULAR_MODE_PROFILE,
    "linear": LINEAR_MODE_PROFILE,
}


def get_mode_profile(mode: DiagramMode | str) -> ModeProfile:
    """Return the immutable defaults for a supported drawing mode."""

    normalized_mode = str(mode)
    if normalized_mode not in MODE_PROFILES:
        raise ValidationError(f"Unsupported drawing mode: {mode!r}.")
    try:
        return MODE_PROFILES[cast(DiagramMode, normalized_mode)]
    except KeyError as exc:  # pragma: no cover - guarded above
        raise ValidationError(f"Unsupported drawing mode: {mode!r}.") from exc


def resolve_mode_profile_overrides(
    mode: DiagramMode | str,
    explicit_overrides: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Merge fresh-request profile values with explicit legacy config overrides."""

    profile = get_mode_profile(mode)
    resolved = profile.config_overrides
    explicit = {
        key: value
        for key, value in dict(explicit_overrides or {}).items()
        if value is not None
    }
    if (
        profile.mode == "linear"
        and explicit.get("linear_ruler_on_axis") is True
        and "linear_axis_stroke_color" not in explicit
    ):
        resolved["linear_axis_stroke_color"] = profile.linear_ruler_axis_color
    resolved.update(explicit)
    return resolved


def mode_profiles_payload() -> dict[str, object]:
    """Return the JSON-compatible profile representation used by the web app."""

    return {
        "version": MODE_PROFILE_VERSION,
        "featureTypes": list(DEFAULT_FEATURE_TYPES),
        "modes": {
            mode: {
                "comparison": {
                    "evalue": profile.comparison.evalue,
                    "bitscore": profile.comparison.bitscore,
                    "identity": profile.comparison.identity,
                    "alignmentLength": profile.comparison.alignment_length,
                },
                "tracks": {
                    "gc": profile.show_gc,
                    "skew": profile.show_skew,
                },
                "linearAxisColor": profile.linear_axis_color,
                "linearRulerAxisColor": profile.linear_ruler_axis_color,
            }
            for mode, profile in MODE_PROFILES.items()
        },
    }


__all__ = [
    "CIRCULAR_MODE_PROFILE",
    "ComparisonThresholds",
    "DEFAULT_FEATURE_TYPES",
    "DiagramMode",
    "LINEAR_MODE_PROFILE",
    "MODE_PROFILE_VERSION",
    "MODE_PROFILES",
    "ModeProfile",
    "get_mode_profile",
    "mode_profiles_payload",
    "resolve_mode_profile_overrides",
]
