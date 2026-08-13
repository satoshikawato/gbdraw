"""Durable semantic selectors used by schema-6 Feature fill rules."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from urllib.parse import quote, unquote

from Bio.SeqRecord import SeqRecord

from gbdraw.core.record_metadata import _source_feature_index
from gbdraw.exceptions import ValidationError


FEATURE_SEMANTIC_SCOPE_QUALIFIER = "__gbdraw_semantic_scope__"
FEATURE_SEMANTIC_SCOPE_PREFIX = "fs1:"
FEATURE_SEMANTIC_SCOPE_KINDS = frozenset(
    {
        "feature-type",
        "base-legend-caption",
        "rendered-label",
        "source-annotation-label",
        "similarity-group",
    }
)

_SOURCE_LABEL_QUALIFIERS = ("product", "gene", "locus_tag", "note")
_URI_COMPONENT_SAFE = "-_.!~*'()"


def encode_feature_semantic_selector(kind: object, value: object) -> str:
    """Return one canonical opaque selector literal."""

    normalized_kind = str(kind or "").strip()
    normalized_value = str(value or "").strip()
    if normalized_kind not in FEATURE_SEMANTIC_SCOPE_KINDS:
        raise ValidationError(
            f"Unsupported Feature semantic selector kind: {normalized_kind!r}."
        )
    if not normalized_value:
        raise ValidationError("Feature semantic selector values must be non-empty.")
    encoded = quote(normalized_value, safe=_URI_COMPONENT_SAFE, encoding="utf-8")
    return f"{FEATURE_SEMANTIC_SCOPE_PREFIX}{normalized_kind}:{encoded}"


def parse_feature_semantic_selector(value: object) -> tuple[str, str]:
    """Decode and validate one canonical schema-6 selector literal."""

    literal = str(value or "")
    if not literal.startswith(FEATURE_SEMANTIC_SCOPE_PREFIX):
        raise ValidationError(
            f"{FEATURE_SEMANTIC_SCOPE_QUALIFIER} requires an fs1_compatible literal."
        )
    remainder = literal[len(FEATURE_SEMANTIC_SCOPE_PREFIX) :]
    kind, separator, encoded = remainder.partition(":")
    try:
        decoded = unquote(encoded, encoding="utf-8", errors="strict")
    except UnicodeDecodeError as exc:
        raise ValidationError(
            f"Invalid UTF-8 in {FEATURE_SEMANTIC_SCOPE_QUALIFIER}."
        ) from exc
    if not separator or encode_feature_semantic_selector(kind, decoded) != literal:
        raise ValidationError(
            f"{FEATURE_SEMANTIC_SCOPE_QUALIFIER} requires a canonical fs1 literal."
        )
    return kind, decoded


def is_feature_semantic_selector(value: object) -> bool:
    try:
        parse_feature_semantic_selector(value)
    except ValidationError:
        return False
    return True


def _normalized_values(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    values = value if isinstance(value, (list, tuple, set)) else (value,)
    return tuple(str(item).strip() for item in values if str(item).strip())


def source_annotation_label(feature: object) -> str:
    """Return the source-label grouping value used by the Web Feature editor."""

    qualifiers = getattr(feature, "qualifiers", None) or {}
    if not isinstance(qualifiers, Mapping):
        return ""
    for wanted in _SOURCE_LABEL_QUALIFIERS:
        key = next(
            (key for key in qualifiers if str(key).lower() == wanted),
            None,
        )
        if key is None:
            continue
        values = _normalized_values(qualifiers[key])
        if values:
            return values[0]
    return ""


def _iter_features(features: object):
    for feature in features or ():  # type: ignore[union-attr]
        yield feature
        yield from _iter_features(getattr(feature, "sub_features", None))


@dataclass(frozen=True)
class FeatureSemanticSelectorValues:
    rendered_label: str
    source_annotation_label: str
    similarity_group_ids: tuple[str, ...] = ()


@dataclass(frozen=True)
class FeatureSemanticSelectorContext:
    """Immutable semantic values keyed by source Feature object identity."""

    _by_feature_object_id: Mapping[int, FeatureSemanticSelectorValues] = field(
        repr=False,
        compare=False,
    )

    def values_for_feature(self, feature: object) -> FeatureSemanticSelectorValues:
        values = self._by_feature_object_id.get(id(feature))
        if values is None:
            raise ValidationError(
                "Feature semantic selector context does not contain the source "
                "Feature being resolved."
            )
        return values


def _orthogroups_by_feature(
    records: Sequence[SeqRecord],
    orthogroups: object | None,
) -> dict[int, set[str]]:
    groups_by_feature: dict[int, set[str]] = {}
    if orthogroups is None:
        return groups_by_feature

    features_by_source_index: dict[tuple[int, int], object] = {}
    for record_index, record in enumerate(records):
        for fallback_index, feature in enumerate(_iter_features(record.features)):
            source_index = _source_feature_index(feature)
            features_by_source_index[
                (record_index, fallback_index if source_index is None else source_index)
            ] = feature

    groups = getattr(orthogroups, "orthogroups", None) or {}
    for raw_group_id, members in groups.items():
        group_id = str(raw_group_id or "").strip()
        if not group_id:
            continue
        for member in members or ():
            try:
                key = (int(member.record_index), int(member.feature_index))
            except (AttributeError, TypeError, ValueError):
                continue
            feature = features_by_source_index.get(key)
            if feature is not None:
                groups_by_feature.setdefault(id(feature), set()).add(group_id)
    return groups_by_feature


def build_feature_semantic_selector_context(
    records: Sequence[SeqRecord] | SeqRecord,
    *,
    label_filtering: Mapping[str, object] | None = None,
    orthogroups: object | None = None,
) -> FeatureSemanticSelectorContext:
    """Build semantic selector values once for drawing, legend, and metadata."""

    # Imported at execution time because label filtering uses selector-value
    # helpers, while those helpers also consume this semantic context.
    from gbdraw.labels.filtering import get_label_text

    record_list = [records] if isinstance(records, SeqRecord) else list(records)
    filtering = dict(label_filtering or {})
    groups_by_feature = _orthogroups_by_feature(record_list, orthogroups)
    values_by_feature: dict[int, FeatureSemanticSelectorValues] = {}
    for record in record_list:
        for feature in _iter_features(record.features):
            values_by_feature[id(feature)] = FeatureSemanticSelectorValues(
                rendered_label=get_label_text(
                    feature,
                    filtering,
                    record_id=str(record.id or ""),
                ),
                source_annotation_label=source_annotation_label(feature),
                similarity_group_ids=tuple(
                    sorted(groups_by_feature.get(id(feature), ()))
                ),
            )
    return FeatureSemanticSelectorContext(
        _by_feature_object_id=MappingProxyType(values_by_feature)
    )


def feature_semantic_rules_require_context(color_map: Mapping[str, object]) -> bool:
    """Return whether rules contain rendered-label or Similarity selectors."""

    for rules_by_qualifier in color_map.values():
        if not isinstance(rules_by_qualifier, Mapping):
            continue
        for rule in rules_by_qualifier.get(FEATURE_SEMANTIC_SCOPE_QUALIFIER, ()):
            if not rule:
                continue
            selector = rule[0]
            if (
                isinstance(selector, tuple)
                and len(selector) == 2
                and selector[0] in {"rendered-label", "similarity-group"}
            ):
                return True
    return False


def feature_semantic_selector_matches(
    feature: object,
    selector: tuple[str, str],
    *,
    context: FeatureSemanticSelectorContext | None = None,
    base_legend_caption: str | None = None,
) -> bool:
    """Return whether one parsed selector matches a source Feature."""

    kind, value = selector
    if kind == "feature-type":
        return str(getattr(feature, "type", "") or "") == value
    if kind == "source-annotation-label":
        return source_annotation_label(feature) == value
    if kind == "base-legend-caption":
        return str(base_legend_caption or "").casefold() == value.casefold()
    if context is None:
        raise ValidationError(
            f"{FEATURE_SEMANTIC_SCOPE_QUALIFIER} selector {kind!r} requires "
            "shared render context. Regenerate from the canonical request."
        )
    values = context.values_for_feature(feature)
    if kind == "rendered-label":
        return values.rendered_label == value
    if kind == "similarity-group":
        return value in values.similarity_group_ids
    raise ValidationError(f"Unsupported Feature semantic selector kind: {kind!r}.")


__all__ = [
    "FEATURE_SEMANTIC_SCOPE_KINDS",
    "FEATURE_SEMANTIC_SCOPE_PREFIX",
    "FEATURE_SEMANTIC_SCOPE_QUALIFIER",
    "FeatureSemanticSelectorContext",
    "FeatureSemanticSelectorValues",
    "build_feature_semantic_selector_context",
    "encode_feature_semantic_selector",
    "feature_semantic_selector_matches",
    "feature_semantic_rules_require_context",
    "is_feature_semantic_selector",
    "parse_feature_semantic_selector",
    "source_annotation_label",
]
