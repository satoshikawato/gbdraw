"""Versioned defaults and shared validation for Circular and Linear modes."""

from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral, Real
from types import MappingProxyType
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


@dataclass(frozen=True)
class ComparisonThresholdDomain:
    """Accepted values of one comparison threshold, published to the Web."""

    minimum: float
    maximum: float | None
    integer: bool
    reason: str
    requirement: str


# Owner of the comparison-threshold domains. The Web evaluates the generated
# copy only where it uses a threshold before Python (LOSAT post-processing).
COMPARISON_THRESHOLD_DOMAINS: Mapping[str, ComparisonThresholdDomain] = MappingProxyType({
    "evalue": ComparisonThresholdDomain(0, None, False, "NONNEGATIVE", "a finite number >= 0"),
    "bitscore": ComparisonThresholdDomain(0, None, False, "NONNEGATIVE", "a finite number >= 0"),
    "identity": ComparisonThresholdDomain(0, 100, False, "PERCENT", "a finite number between 0 and 100"),
    "alignment_length": ComparisonThresholdDomain(0, None, True, "NONNEGATIVE_INTEGER", "an integer >= 0"),
})


def validate_comparison_threshold(field_name: str, value: object) -> float | int:
    """Return one threshold normalized by its domain or raise a typed error."""

    domain = COMPARISON_THRESHOLD_DOMAINS[field_name]
    normalized: float | int | None = None
    if domain.integer:
        if not isinstance(value, bool) and isinstance(value, Integral):
            normalized = int(value)
    elif not isinstance(value, bool) and isinstance(value, Real):
        try:
            normalized = float(value)
        except (OverflowError, TypeError, ValueError):
            normalized = None
    if (
        normalized is None
        or not math.isfinite(normalized)
        or normalized < domain.minimum
        or (domain.maximum is not None and normalized > domain.maximum)
    ):
        raise ValidationError(
            f"{field_name} must be {domain.requirement}.",
            diagnostic={"code": "INPUT_INVALID", "field": field_name, "reason": domain.reason},
        )
    return normalized


_DINUCLEOTIDE_BASES = frozenset("ACGTU")


def validate_dinucleotide(value: object, *, field_name: str = "dinucleotide") -> str:
    """Return an upper-case pair from A, C, G, T, and U (U is counted as T)."""

    normalized = value.strip().upper() if isinstance(value, str) else ""
    if len(normalized) != 2 or not set(normalized) <= _DINUCLEOTIDE_BASES:
        raise ValidationError(
            f"{field_name} must be two letters from A, C, G, T, and U.",
            diagnostic={"code": "INPUT_INVALID", "field": "dinucleotide", "reason": "DINUCLEOTIDE"},
        )
    return normalized


@dataclass(frozen=True)
class ComparisonThresholds:
    """Validated comparison and conservation filter thresholds."""

    evalue: float
    bitscore: float
    identity: float
    alignment_length: int = 0

    def __post_init__(self) -> None:
        for field_name in COMPARISON_THRESHOLD_DOMAINS:
            object.__setattr__(
                self,
                field_name,
                validate_comparison_threshold(field_name, getattr(self, field_name)),
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
            "canvas.show_gc": self.show_gc,
            "canvas.show_skew": self.show_skew,
        }
        if self.linear_axis_color is not None:
            overrides["objects.axis.linear.stroke_color"] = self.linear_axis_color
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
    """Merge fresh-request profile values with explicit canonical config overrides."""

    profile = get_mode_profile(mode)
    resolved = profile.config_overrides
    explicit = dict(explicit_overrides or {})
    if (
        profile.mode == "linear"
        and explicit.get("canvas.linear.ruler_on_axis") is True
        and explicit.get("objects.scale.show", True) is True
        and "objects.axis.linear.stroke_color" not in explicit
    ):
        resolved["objects.axis.linear.stroke_color"] = profile.linear_ruler_axis_color
    resolved.update(explicit)
    return resolved


def _camel_threshold_name(name: str) -> str:
    head, *tail = name.split("_")
    return head + "".join(part.title() for part in tail)


def mode_profiles_payload() -> dict[str, object]:
    """Return the JSON-compatible profile representation used by the web app."""

    return {
        "version": MODE_PROFILE_VERSION,
        "comparisonDomains": {
            _camel_threshold_name(name): {
                "minimum": domain.minimum,
                "maximum": domain.maximum,
                "integer": domain.integer,
                "reason": domain.reason,
            }
            for name, domain in COMPARISON_THRESHOLD_DOMAINS.items()
        },
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
    "COMPARISON_THRESHOLD_DOMAINS",
    "ComparisonThresholdDomain",
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
    "validate_comparison_threshold",
    "validate_dinucleotide",
]
