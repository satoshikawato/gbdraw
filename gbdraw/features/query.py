"""Feature queries: which features of a record a rule or a track reads.

One owner of the feature condition the visibility table and the feature density tracks share: a feature type (exact,
or ``*`` for every type), an optional qualifier whose values a regular expression searches (case-insensitively),
and an optional range of feature length in bases. The qualifier may be one of the pseudo-qualifiers the visibility
table has always read (``hash``, ``location``, ``record_location``), which name values the feature carries rather
than a GenBank qualifier.

The color rules keep their own matcher in ``gbdraw.features.selector_values`` (``iter_specific_color_rules``): they
rank exact types before ``*`` and index long rule lists by their literal patterns.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cache

from ..exceptions import ParseError, ValidationError
from .selector_values import (
    _matches_constraint,
    get_feature_hash,
    get_feature_location_str,
    get_feature_qualifiers,
    get_feature_record_location_str,
    get_feature_type,
    get_qualifier_values,
)

# The qualifiers a query reads from the feature's drawn values instead of its GenBank qualifiers.
PSEUDO_QUALIFIERS = ("hash", "location", "record_location")

# A feature's drawn values for the pseudo-qualifiers, when the caller has them instead of the feature (the Web).
FeatureSelectorValues = Mapping[str, object]


# The diagnostic of a refused feature query: the Web words it as an invalid input.
_INVALID_QUERY = {"code": "INPUT_INVALID"}


@dataclass(frozen=True)
class FeatureQuery:
    """Which features a rule or track reads.

    ``feature_type`` is an exact feature type or ``*``. ``qualifier`` names a GenBank qualifier (any case) or a
    pseudo-qualifier; ``value`` is a regular expression searched in its values, case-insensitively. A qualifier
    without a value matches the features that carry that qualifier. ``min_length`` and ``max_length`` bound the
    bases the feature covers (the sum of its parts), inclusive.
    """

    feature_type: str = "*"
    qualifier: str | None = None
    value: str | None = None
    min_length: int | None = None
    max_length: int | None = None

    def __post_init__(self) -> None:
        feature_type = str(self.feature_type or "").strip()
        if not feature_type:
            raise ValidationError("A feature query needs a feature type (or * for every type).", diagnostic=_INVALID_QUERY)
        object.__setattr__(self, "feature_type", feature_type)
        qualifier = str(self.qualifier).strip() if self.qualifier is not None else ""
        object.__setattr__(self, "qualifier", qualifier or None)
        value = str(self.value).strip() if self.value is not None else ""
        if value and not qualifier:
            raise ValidationError(f"The feature query value '{value}' needs a qualifier to search.", diagnostic=_INVALID_QUERY)
        object.__setattr__(self, "value", value or None)
        if value:
            try:
                _compiled(value)
            except re.error as error:
                raise ParseError(f"Invalid regex '{value}' ({error})") from error
        for name in ("min_length", "max_length"):
            bound = getattr(self, name)
            if bound is not None and (isinstance(bound, bool) or int(bound) != bound or int(bound) < 0):
                raise ValidationError(f"The feature query {name} must be a whole number of bases >= 0; got {bound!r}.", diagnostic=_INVALID_QUERY)
        if self.min_length is not None and self.max_length is not None and self.min_length > self.max_length:
            raise ValidationError(
                f"The feature query min_length ({self.min_length}) is above its max_length ({self.max_length}).",
                diagnostic=_INVALID_QUERY,
            )

    @property
    def pattern(self) -> re.Pattern[str] | None:
        """The compiled value pattern, compiled once per distinct value."""

        return None if self.value is None else _compiled(self.value)

    @property
    def has_length_bounds(self) -> bool:
        return self.min_length is not None or self.max_length is not None

    def matches_type(self, feature_type: str | None) -> bool:
        return _matches_constraint(self.feature_type, feature_type)

    def _qualifier_values(
        self, feature: object, record_id: str | None, selector: FeatureSelectorValues | None
    ) -> list[str]:
        key = str(self.qualifier).lower()
        if key in PSEUDO_QUALIFIERS:
            if selector is not None:
                value = selector.get(key)
            elif key == "record_location":
                value = get_feature_record_location_str(feature, record_id)
            elif key == "hash":
                value = get_feature_hash(feature, record_id)
            else:
                value = get_feature_location_str(feature)
            return [str(value)] if value else []
        return get_qualifier_values(get_feature_qualifiers(feature), str(self.qualifier))

    def matches_qualifier(
        self, feature: object, *, record_id: str | None = None, selector: FeatureSelectorValues | None = None
    ) -> bool:
        """Whether the feature's qualifier values match (true for a query without a qualifier).

        ``selector`` gives the drawn values of the pseudo-qualifiers when the caller has them instead of the feature.
        """

        if self.qualifier is None:
            return True
        values = self._qualifier_values(feature, record_id, selector)
        pattern = self.pattern
        if pattern is None:
            return bool(values) or str(self.qualifier).lower() in {
                str(key).lower() for key in get_feature_qualifiers(feature)
            }
        return any(pattern.search(value) for value in values)

    def matches_length(self, feature: object) -> bool:
        """Whether the bases the feature covers (the sum of its parts) lie in the length bounds."""

        if not self.has_length_bounds:
            return True
        location = getattr(feature, "location", None)
        if location is None:
            return False
        length = len(location)
        if self.min_length is not None and length < self.min_length:
            return False
        return self.max_length is None or length <= self.max_length

    def matches(
        self, feature: object, *, record_id: str | None = None, selector: FeatureSelectorValues | None = None
    ) -> bool:
        """Whether ``feature`` is one the query reads."""

        return (
            self.matches_type(get_feature_type(feature))
            and self.matches_length(feature)
            and self.matches_qualifier(feature, record_id=record_id, selector=selector)
        )

    def describe(self) -> str:
        """The query in words, as a Legend caption shows it: ``repeat_region · rpt_type = long_terminal_repeat``."""

        parts = ["any feature" if self.feature_type == "*" else self.feature_type]
        if self.qualifier is not None:
            parts.append(self.qualifier if self.value is None else f"{self.qualifier} = {self.value}")
        if self.min_length is not None and self.max_length is not None:
            parts.append(f"length {format_bases(self.min_length)}–{format_bases(self.max_length)}")
        elif self.min_length is not None:
            parts.append(f"length ≥ {format_bases(self.min_length)}")
        elif self.max_length is not None:
            parts.append(f"length ≤ {format_bases(self.max_length)}")
        return " · ".join(parts)


@cache
def _compiled(value: str) -> re.Pattern[str]:
    return re.compile(value, re.IGNORECASE)


def format_bases(value: int | float) -> str:
    """A length in bases in words: ``500 bp``, ``5 kbp``, ``2.5 Mbp``."""

    number = float(value)
    for unit, size in (("Mbp", 1_000_000), ("kbp", 1_000)):
        if number >= size:
            scaled = number / size
            return f"{scaled:g} {unit}"
    return f"{number:g} bp"


__all__ = ["FeatureQuery", "PSEUDO_QUALIFIERS", "format_bases"]
