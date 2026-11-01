"""Versioned defaults and shared validation for Circular and Linear modes."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import math
import re
from numbers import Integral, Real
from types import MappingProxyType
from typing import Literal, Mapping, cast, get_args

from gbdraw.config.models.canvas import (
    LINEAR_AXIS_RULER_COLOR_PATHS,
    LINEAR_WRAP_LENGTH_PATH,
    linear_axis_ruler_applies,
    linear_wrap_length_from_config,
)
from gbdraw.config.models.labels import LinearLabelScope
from gbdraw.config.toml import load_config_toml
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


_DINUCLEOTIDE_BASES = "".join(sorted("ACGTU"))
# Two bases in either case, with outer spaces (the request description names it).
DINUCLEOTIDE_PATTERN = rf"^\s*[{_DINUCLEOTIDE_BASES}{_DINUCLEOTIDE_BASES.lower()}]{{2}}\s*$"


def validate_dinucleotide(value: object, *, field_name: str = "dinucleotide") -> str:
    """Return an upper-case pair from A, C, G, T, and U (U is counted as T)."""

    if not isinstance(value, str) or re.search(DINUCLEOTIDE_PATTERN, value) is None:
        raise ValidationError(
            f"{field_name} must be two letters from A, C, G, T, and U.",
            diagnostic={"code": "INPUT_INVALID", "field": "dinucleotide", "reason": "DINUCLEOTIDE"},
        )
    return value.strip().upper()


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


# What a request that cuts a record and leaves these unset draws (the other values are refused then).
LINEAR_WRAP_DEFAULT_TRACK_LAYOUT = "above"
_LINEAR_WRAP_DEFAULTS: Mapping[str, object] = MappingProxyType({
    "canvas.linear.track_layout": LINEAR_WRAP_DEFAULT_TRACK_LAYOUT,
    "objects.scale.style": "ruler",
})


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
        and "objects.axis.linear.stroke_color" not in explicit
    ):
        # A fresh request starts from the packaged config, so its values
        # decide the leaves the request leaves unset.
        decisive = {**_packaged_axis_ruler_leaves(), **explicit}
        if linear_axis_ruler_applies(
            show_scale=decisive["objects.scale.show"],
            scale_style=decisive["objects.scale.style"],
            track_layout=decisive["canvas.linear.track_layout"],
        ):
            resolved["objects.axis.linear.stroke_color"] = profile.linear_ruler_axis_color
    resolved.update(explicit)
    return resolved


def linear_wrap_default_overrides(explicit_overrides: Mapping[str, object] | None = None) -> dict[str, object] | None:
    """What a fresh request that asks to wrap draws with if a record is cut; ``None`` when it does not ask.

    They are overrides for the planner to apply when it knows the record lengths. A figure with nothing to cut
    never reads them, so it draws what the same request without a wrap length draws. The request's own leaves stay
    out: those it names as something rows cannot draw are refused, not defaulted. The scale's colors it names are
    restated, because the axis ruler rule sets the ones its overrides leave out.
    """

    explicit = dict(explicit_overrides or {})
    if linear_wrap_length_from_config(explicit.get(LINEAR_WRAP_LENGTH_PATH)) == "none":
        return None
    defaults: dict[str, object] = {
        **{path: value for path, value in _LINEAR_WRAP_DEFAULTS.items() if path not in explicit},
        # Fixed while wrapping: the option's other value has nothing to draw (the CLI always sends one).
        "canvas.linear.ruler_on_axis": True,
        **{path: explicit[path] for path in LINEAR_AXIS_RULER_COLOR_PATHS if explicit.get(path) is not None},
    }
    if "objects.axis.linear.stroke_color" not in explicit:
        decisive = {**_packaged_axis_ruler_leaves(), **explicit, **defaults}
        if linear_axis_ruler_applies(
            show_scale=decisive["objects.scale.show"],
            scale_style=decisive["objects.scale.style"],
            track_layout=decisive["canvas.linear.track_layout"],
        ):
            defaults["objects.axis.linear.stroke_color"] = LINEAR_MODE_PROFILE.linear_ruler_axis_color
    return defaults


def effective_linear_track_leaves(
    explicit_overrides: Mapping[str, object] | None = None,
) -> tuple[str, str]:
    """The track layout and scale style of a fresh Linear request: explicit, then packaged."""

    decisive = {**_packaged_axis_ruler_leaves(), **dict(explicit_overrides or {})}
    return str(decisive["canvas.linear.track_layout"]), str(decisive["objects.scale.style"])


@lru_cache(maxsize=1)
def _packaged_axis_ruler_leaves() -> dict[str, object]:
    raw = load_config_toml("gbdraw.data", "config.toml")
    return {
        "objects.scale.show": raw["objects"]["scale"]["show"],
        "objects.scale.style": raw["objects"]["scale"]["style"],
        "canvas.linear.track_layout": raw["canvas"]["linear"]["track_layout"],
    }


def _camel_threshold_name(name: str) -> str:
    head, *tail = name.split("_")
    return head + "".join(part.title() for part in tail)


def mode_profiles_payload() -> dict[str, object]:
    """Return the JSON-compatible profile representation used by the web app."""

    # gbdraw.auto_sizes takes DiagramMode from this module.
    from gbdraw.auto_sizes import linear_label_size_records, tick_and_window_tiers
    from gbdraw.comparisons.identity_scale import IDENTITY_SCALE_PATHS
    from gbdraw.config.models import GbdrawConfig
    from gbdraw.diagrams.circular.presets import circular_auto_slot_sizes, circular_track_preset_table
    from gbdraw.labels.label_style import (
        DEFAULT_ITALIC_GENE_LABELS,
        DEFAULT_LABEL_LINE,
        DEFAULT_LABEL_PROTECTION,
        LABEL_LINES,
        LABEL_PROTECTIONS,
    )

    packaged = load_config_toml("gbdraw.data", "config.toml")
    packaged_cfg = GbdrawConfig.from_dict(packaged)
    length_threshold = packaged["labels"]["length_threshold"]
    return {
        "version": MODE_PROFILE_VERSION,
        # The packaged text font list: a resolved config that holds it names no font.
        "defaultFontFamily": str(packaged["objects"]["text"]["font_family"]),
        # The packaged short/long boundary of each mode, the records that set
        # the Auto Linear label font for each label display scope (rules:
        # gbdraw/auto_sizes.py), and the Auto sizes of a Circular row
        # (gbdraw/diagrams/circular/presets.py).
        "autoSizes": {
            "circularSlots": circular_auto_slot_sizes(packaged_cfg),
            "linearLabelSizeRecords": {
                scope: linear_label_size_records(scope) for scope in get_args(LinearLabelScope)
            },
            "sizeClassThresholdBp": {mode: int(length_threshold[mode]) for mode in MODE_PROFILES},
            **tick_and_window_tiers(packaged_cfg),
        },
        # The Circular track presets in menu order and their feature lane directions.
        "circularTrackPresets": circular_track_preset_table(),
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
        # The values of labels.<mode>.line and .protection.
        "labelStyleChoices": {"line": list(LABEL_LINES), "protection": list(LABEL_PROTECTIONS)},
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
                # The packaged feature label style of the mode (labels.<mode>.line, .protection,
                # .italic_gene_labels).
                "labelStyle": {
                    "line": DEFAULT_LABEL_LINE[mode],
                    "protection": DEFAULT_LABEL_PROTECTION[mode],
                    "italicGeneLabels": DEFAULT_ITALIC_GENE_LABELS[mode],
                },
                # The config leaf of the mode's fixed identity colour scale (an ``[MIN, MAX]`` override).
                "identityScalePath": IDENTITY_SCALE_PATHS[mode],
            }
            for mode, profile in MODE_PROFILES.items()
        },
    }


__all__ = [
    "LINEAR_WRAP_DEFAULT_TRACK_LAYOUT",
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
    "effective_linear_track_leaves",
    "get_mode_profile",
    "linear_wrap_default_overrides",
    "mode_profiles_payload",
    "resolve_mode_profile_overrides",
    "validate_comparison_threshold",
    "validate_dinucleotide",
]
