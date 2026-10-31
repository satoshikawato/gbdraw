"""The identity colour scale of comparisons: the one owner of its domain.

A comparison colours each pairwise match (Linear) or similarity-ring hit
(Circular) by where its percent identity falls in a domain, and its Legend
labels the two ends of that domain. By default the domain runs from the
mode's identity threshold to 100%, so it moves with the threshold. A fixed
scale (``objects.blast_match.identity_scale`` in Linear,
``objects.conservation.identity_scale`` in Circular; CLI ``--identity_scale
MIN MAX``) keeps the same colours and Legend labels across figures whatever
their thresholds.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from numbers import Real
from types import MappingProxyType
from typing import Mapping

from gbdraw.exceptions import ValidationError

# A fixed scale: (minimum, maximum) percent identity, 0 <= minimum < maximum <= 100.
IdentityScale = tuple[float, float]

LINEAR_IDENTITY_SCALE_PATH = "objects.blast_match.identity_scale"
CIRCULAR_IDENTITY_SCALE_PATH = "objects.conservation.identity_scale"
# The leaf each diagram mode reads; the Web reads this map from the generated mode profiles.
IDENTITY_SCALE_PATHS: Mapping[str, str] = MappingProxyType(
    {"circular": CIRCULAR_IDENTITY_SCALE_PATH, "linear": LINEAR_IDENTITY_SCALE_PATH}
)
# The text a config file or an override reads as "no fixed scale" (threshold to 100%).
_AUTO_SPELLINGS = frozenset({"", "auto", "none", "null"})


def identity_scale_from_config(value: object, *, path: str = LINEAR_IDENTITY_SCALE_PATH) -> IdentityScale | None:
    """The one validator of an identity scale leaf, for a config file, an override, and the CLI.

    ``None`` and ``"auto"`` are no fixed scale; two finite numbers with
    0 <= minimum < maximum <= 100 are one.
    """

    if value is None or (isinstance(value, str) and value.strip().lower() in _AUTO_SPELLINGS):
        return None
    if (
        isinstance(value, (list, tuple))
        and len(value) == 2
        and all(isinstance(item, Real) and not isinstance(item, bool) for item in value)
    ):
        minimum, maximum = float(value[0]), float(value[1])
        if math.isfinite(minimum) and math.isfinite(maximum) and 0.0 <= minimum < maximum <= 100.0:
            return (minimum, maximum)
    raise ValidationError(
        f'{path} must be "auto" or two percent identities [MIN, MAX] with 0 <= MIN < MAX <= 100; got {value!r}.',
        diagnostic={"code": "INPUT_INVALID", "configPath": path},
    )


@dataclass(frozen=True)
class IdentityColorDomain:
    """The percent identities at the two ends of a comparison's colour ramp."""

    minimum: float
    maximum: float

    def factor(self, identity: float) -> float:
        """Where ``identity`` falls on the ramp, clamped to 0..1.

        A domain without width (a threshold of 100% without a fixed scale)
        colours a match at its maximum or above with the end colour and
        every other match with the start colour.
        """

        if self.minimum >= self.maximum:
            return 1.0 if float(identity) >= self.maximum else 0.0
        factor = (float(identity) - self.minimum) / (self.maximum - self.minimum)
        return max(0.0, min(1.0, factor))


def identity_color_domain(threshold: float, scale: IdentityScale | None) -> IdentityColorDomain:
    """The colour domain of a comparison: its fixed scale, or its threshold to 100%."""

    if scale is not None:
        return IdentityColorDomain(float(scale[0]), float(scale[1]))
    return IdentityColorDomain(float(threshold), 100.0)


__all__ = [
    "CIRCULAR_IDENTITY_SCALE_PATH",
    "IdentityColorDomain",
    "IDENTITY_SCALE_PATHS",
    "IdentityScale",
    "LINEAR_IDENTITY_SCALE_PATH",
    "identity_color_domain",
    "identity_scale_from_config",
]
