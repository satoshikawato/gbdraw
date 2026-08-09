"""Compact feature-catalog normalization for browser render responses."""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Mapping, Sequence
import math
import re
import xml.etree.ElementTree as ET
from typing import Any

from gbdraw.exceptions import GbdrawError
from gbdraw.render.interactive_svg import (
    InteractiveSvgContext,
    _add_class_token,
    _collect_rendered_features,
    _compact_wire_value,
    _element_match_id_status,
    _feature_payloads,
    _feature_record_index_status,
    _feature_rendered_id_status,
    _feature_source_index_status,
    _feature_stable_id,
    _feature_stable_id_status,
    _first_text,
    _is_match_candidate,
    _normalize_string_array,
)

FEATURE_CATALOG_SCHEMA = 3

_BIOLOGICAL_ALIAS_KEYS = {
    "id",
    "svg_id",
    "svgId",
    "stable_svg_id",
    "stableSvgId",
    "stable_feature_id",
    "stableFeatureId",
    "stableFeatureSvgId",
    "stable_feature_svg_id",
    "featureSvgId",
    "feature_svg_id",
    "rendered_svg_id",
    "renderedSvgId",
    "rendered_feature_svg_id",
    "renderedFeatureSvgId",
}
_ORTHOGROUP_FEATURE_KEYS = {
    "orthogroup_id",
    "orthogroup_member_count",
    "orthogroup_record_coverage",
    "orthogroup_representative",
    "orthogroup_member",
    "orthogroupId",
    "orthogroupMemberCount",
    "orthogroupRecordCoverage",
    "orthogroupRepresentative",
    "orthogroupMember",
}
_MATCH_REDUNDANT_SUFFIXES = (
    "feature_svg_id",
    "stable_feature_svg_id",
    "feature_index",
    "protein_id",
    "locus_id",
    "display_name",
)
_ORTHOGROUP_COLLECTION_COUNTS = (
    ("orthologEdges", "orthologEdgeCount"),
    ("orthologPaths", "orthologPathCount"),
    ("relatedEdges", "relatedEdgeCount"),
)
_SEQUENCE_WHITESPACE_RE = re.compile(r"\s")
_PRESENTATION_KIND_BY_SCOPE = {
    "adjacent_local": "collinear_gene_group",
    "global_collinear": "orthogroup",
}
_PRESENTATION_BY_SEARCH_SCOPE = {
    "adjacent": ("adjacent_local", "collinear_gene_group"),
    "all": ("global_collinear", "orthogroup"),
}


def _sequence(value: object | None) -> list[object]:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return list(value)
    return []


def _text(value: object | None) -> str:
    return _first_text(value)


_MATCH_ATTRIBUTES = {
    "query_record_index": "data-query-record-index",
    "subject_record_index": "data-subject-record-index",
    "qstart": "data-qstart",
    "qend": "data-qend",
    "sstart": "data-sstart",
    "send": "data-send",
    "identity": "data-identity",
    "alignment_length": "data-alignment-length",
    "evalue": "data-evalue",
    "bitscore": "data-bitscore",
    "mismatches": "data-mismatches",
    "gap_opens": "data-gap-opens",
    "source_index": "data-source-index",
    "track_index": "data-track-index",
    "track_label": "data-track-label",
    "reference_side": "data-reference-side",
    "reference_record_id": "data-reference-record-id",
    "block_kind": "data-collinearity-block-kind",
    "group_scope": "data-group-scope",
    "collinear_group_scope": "data-collinear-group-scope",
    "group_kind": "data-group-kind",
    "block_color_mode": "data-collinearity-color-mode",
    "block_score": "data-collinearity-block-score",
    "block_evalue": "data-collinearity-block-evalue",
    "anchor_index": "data-collinearity-anchor-index",
    "anchor_count": "data-collinearity-anchor-count",
    "query_feature_svg_id": "data-query-feature-svg-id",
    "subject_feature_svg_id": "data-subject-feature-svg-id",
    "query_stable_feature_svg_id": "data-query-stable-feature-svg-id",
    "subject_stable_feature_svg_id": "data-subject-stable-feature-svg-id",
    "query_feature_index": "data-query-feature-index",
    "subject_feature_index": "data-subject-feature-index",
    "query_protein_id": "data-query-protein-id",
    "subject_protein_id": "data-subject-protein-id",
    "query_locus_id": "data-query-locus-id",
    "subject_locus_id": "data-subject-locus-id",
    "query_display_name": "data-query-display-name",
    "subject_display_name": "data-subject-display-name",
}


def _match_kind(element: ET.Element) -> str:
    value = _text(element.get("data-match-kind")).lower()
    if value in {"pairwise", "orthogroup", "collinear", "homology"}:
        return value
    if element.get("data-collinearity-block-id"):
        return "collinear"
    if element.get("data-orthogroup-id"):
        return "orthogroup"
    return "pairwise"


def _metadata_values(value: object | None) -> list[str]:
    return list(
        dict.fromkeys(
            part.strip() for part in _text(value).split(";") if part.strip()
        )
    )


def _match_payload(element: ET.Element, index: int) -> dict[str, object]:
    match_id, valid = _element_match_id_status(element)
    if not valid:
        raise GbdrawError("Rendered match contains conflicting match ID aliases.")
    match_id = match_id or f"match_{index + 1}"
    element.set("data-gbdraw-match-id", match_id)
    element.set("data-gbdraw-pairwise-match-id", match_id)
    payload: dict[str, object] = {
        "id": match_id,
        "match_kind": _match_kind(element),
        "orthogroup_ids": _metadata_values(element.get("data-orthogroup-id")),
        "collinearity_block_id": _text(element.get("data-collinearity-block-id")),
        "fill": _first_text(element.get("fill"), "#94a3b8"),
        "query_record_id": _first_text(
            element.get("data-query-record-id"), element.get("data-query")
        ),
        "subject_record_id": _first_text(
            element.get("data-subject-record-id"), element.get("data-subject")
        ),
        "orientation": _first_text(
            element.get("data-collinearity-orientation"),
            element.get("data-orientation"),
        ),
    }
    payload.update(
        {
            field: _text(element.get(attribute))
            for field, attribute in _MATCH_ATTRIBUTES.items()
        }
    )
    return dict(_compact_wire_value(payload) or {})


def _match_payloads(
    root: ET.Element,
    features: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    rendered_features = {
        _text(feature.get("svg_id")): feature
        for feature in features
        if _text(feature.get("svg_id"))
    }
    payloads: list[dict[str, object]] = []
    seen: set[str] = set()
    for element in root.iter():
        if not _is_match_candidate(element):
            continue
        payload = _match_payload(element, len(payloads))
        match_id = _text(payload.get("id"))
        if not match_id or match_id in seen:
            raise GbdrawError("Rendered SVG contains invalid or duplicate match IDs.")
        seen.add(match_id)
        for role in ("query", "subject"):
            rendered = rendered_features.get(
                _text(payload.get(f"{role}_feature_svg_id"))
            )
            if rendered is None:
                continue
            payload.setdefault(f"{role}_record_index", rendered.get("record_idx"))
            payload.setdefault(
                f"{role}_stable_feature_svg_id", _feature_stable_id(rendered)
            )
            payload.setdefault(f"{role}_feature_index", rendered.get("feature_index"))
        payloads.append(dict(_compact_wire_value(payload) or {}))
        element.set("data-gbdraw-interactive-match", "true")
        _add_class_token(element, "gbdraw-interactive-pairwise-match")
    return payloads


def _sequence_sources_for_matches(
    matches: Sequence[Mapping[str, object]],
    sources: Sequence[Mapping[str, object]],
) -> list[Mapping[str, object]]:
    if not matches:
        return []
    linear_indexes: set[int] = set()
    circular_record_ids: set[str] = set()
    comparison_source_indexes: set[int] = set()
    for match in matches:
        if _text(match.get("match_kind")) == "homology":
            reference_side = _text(match.get("reference_side"))
            circular_record_ids.add(_text(match.get(f"{reference_side}_record_id")))
            try:
                comparison_source_indexes.add(int(str(match.get("source_index"))))
            except (TypeError, ValueError):
                pass
            continue
        for role in ("query", "subject"):
            try:
                linear_indexes.add(int(str(match.get(f"{role}_record_index"))))
            except (TypeError, ValueError):
                pass

    selected: list[Mapping[str, object]] = []
    for source in sources:
        origin = _text(source.get("origin"))
        try:
            if (
                origin == "linear-record"
                and int(str(source.get("recordIndex"))) in linear_indexes
            ):
                selected.append(source)
            elif (
                origin == "homology-comparison"
                and int(str(source.get("sourceIndex")))
                in comparison_source_indexes
            ):
                selected.append(source)
        except (TypeError, ValueError):
            continue
        if origin == "circular-reference":
            aliases = {
                _text(source.get("recordId")),
                *_normalize_string_array(source.get("aliases")),
            }
            if aliases & circular_record_ids:
                selected.append(source)
    return selected


def _text_alias_status(
    payload: Mapping[str, object],
    keys: Sequence[str],
) -> tuple[str, bool]:
    values = {
        _text(payload.get(key))
        for key in keys
        if key in payload and _text(payload.get(key))
    }
    return (next(iter(values), ""), len(values) <= 1)


def _group_presentation(
    payload: Mapping[str, object],
) -> tuple[str, str] | None:
    presentation_scope_keys = ("presentationScope", "presentation_scope")
    collinear_scope_keys = (
        "collinearGroupScope",
        "collinear_group_scope",
    )
    kind_keys = ("groupKind", "group_kind")
    all_keys = presentation_scope_keys + collinear_scope_keys + kind_keys
    if not any(key in payload for key in all_keys):
        return None
    presentation_values = [
        _text(payload.get(key))
        for key in presentation_scope_keys
        if key in payload
    ]
    collinear_values = [
        _text(payload.get(key))
        for key in collinear_scope_keys
        if key in payload
    ]
    kind_values = [_text(payload.get(key)) for key in kind_keys if key in payload]
    presentation_scope = next(iter(presentation_values), "")
    collinear_scope = next(iter(collinear_values), "")
    group_kind = next(iter(kind_values), "")
    expected_kind = _PRESENTATION_KIND_BY_SCOPE.get(presentation_scope)
    if (
        not presentation_values
        or not collinear_values
        or not kind_values
        or any(not value for value in presentation_values + collinear_values + kind_values)
        or len(set(presentation_values)) != 1
        or len(set(collinear_values)) != 1
        or len(set(kind_values)) != 1
        or not presentation_scope
        or collinear_scope != presentation_scope
        or not expected_kind
        or group_kind != expected_kind
    ):
        raise GbdrawError(
            "Feature catalog orthogroup presentation metadata is incomplete or conflicting."
        )
    return presentation_scope, group_kind


def _presentation_from_search_scope(
    search_scope: object | None,
) -> tuple[str, str] | None:
    normalized = _text(search_scope).lower()
    if not normalized:
        return None
    presentation = _PRESENTATION_BY_SEARCH_SCOPE.get(normalized)
    if presentation is None:
        raise GbdrawError(
            "Feature catalog collinearity search scope must be 'adjacent' or 'all'."
        )
    return presentation


def _first_qualifier_value(value: object | None) -> str:
    values = _sequence(value)
    if not values:
        values = [value] if value is not None else []
    return next(
        (
            raw
            for entry in values
            if (raw := str(entry)) and raw.strip()
        ),
        "",
    )


def _record_index(feature: Mapping[str, object]) -> int | None:
    record_index, valid = _feature_record_index_status(feature)
    if not valid:
        raise GbdrawError(
            "Feature metadata contains conflicting or invalid record-index aliases."
        )
    return record_index


def _stable_feature_id(feature: Mapping[str, object]) -> str:
    stable_id, valid = _feature_stable_id_status(feature)
    if not valid:
        raise GbdrawError(
            "Feature metadata contains conflicting stable feature ID aliases."
        )
    return stable_id


def _feature_index(feature: Mapping[str, object]) -> int | None:
    feature_index, valid = _feature_source_index_status(feature)
    if not valid:
        raise GbdrawError(
            "Feature metadata contains conflicting or invalid source feature-index aliases."
        )
    return feature_index


def _record_key_alias(feature: Mapping[str, object]) -> tuple[str, bool]:
    values = {
        _text(feature.get(key))
        for key in ("record_key", "recordKey")
        if key in feature and _text(feature.get(key))
    }
    return (next(iter(values), ""), len(values) <= 1)


def _biological_feature_id_alias(
    feature: Mapping[str, object],
) -> tuple[str, bool]:
    values = {
        _text(feature.get(key))
        for key in ("biological_feature_id", "biologicalFeatureId")
        if key in feature and _text(feature.get(key))
    }
    return (next(iter(values), ""), len(values) <= 1)


def _rendered_svg_id(feature: Mapping[str, object]) -> str:
    rendered_id, rendered_valid = _feature_rendered_id_status(feature)
    svg_aliases = {
        _text(feature.get(key))
        for key in ("svg_id", "svgId")
        if key in feature and _text(feature.get(key))
    }
    values = set(svg_aliases)
    if rendered_id:
        values.add(rendered_id)
    if not rendered_valid or len(values) != 1:
        raise GbdrawError(
            "Rendered feature metadata contains missing or conflicting SVG ID aliases."
        )
    return next(iter(values))


def _validated_source_identity(
    feature: Mapping[str, object],
    record_keys: Sequence[str],
    *,
    description: str,
) -> tuple[int, str, int | None]:
    record_index = _record_index(feature)
    stable_id = _stable_feature_id(feature)
    feature_index = _feature_index(feature)
    record_key, record_key_valid = _record_key_alias(feature)
    if (
        record_index is None
        or not 0 <= record_index < len(record_keys)
        or not record_key_valid
        or (record_key and record_key != record_keys[record_index])
        or (not stable_id and feature_index is None)
    ):
        raise GbdrawError(
            f"{description} must identify one biological feature by a "
            "record-scoped stable feature ID or source feature index and a "
            "valid record index/key."
        )
    return record_index, stable_id, feature_index


def _record_keys(context: InteractiveSvgContext) -> list[str]:
    keys = [_text(value) for value in context.record_keys]
    if not keys:
        source_features = (
            context.biological_features
            if context.biological_features
            else context.features
        )
        record_indexes = []
        for feature in source_features:
            if not isinstance(feature, Mapping):
                continue
            record_index = _record_index(feature)
            if record_index is None:
                raise GbdrawError(
                    "Feature catalog source features require an explicit "
                    "nonnegative record index."
                )
            record_indexes.append(record_index)
        record_count = max(record_indexes, default=-1) + 1
        if set(record_indexes) != set(range(record_count)):
            raise GbdrawError(
                "Feature catalog source record indexes must be contiguous from zero."
            )
        keys = [
            f"record-{record_index + 1}"
            for record_index in range(record_count)
        ]
    if any(not key for key in keys) or len(set(keys)) != len(keys):
        raise GbdrawError("Feature catalog record keys must be non-empty and unique.")
    return keys


def _record_key(record_keys: Sequence[str], record_index: int | None) -> str:
    if record_index is None or not 0 <= record_index < len(record_keys):
        raise GbdrawError(
            "Feature catalog source identity has an invalid record index."
        )
    return str(record_keys[record_index])


def _identity_base(feature: Mapping[str, object]) -> str:
    stable_id = _stable_feature_id(feature)
    if stable_id:
        return stable_id
    feature_index = _feature_index(feature)
    if feature_index is not None:
        return f"feature-{feature_index}"
    raise GbdrawError(
        "Feature catalog source identity requires a stable feature ID or "
        "source feature index."
    )


def _normalized_biological_features(
    context: InteractiveSvgContext,
    record_keys: Sequence[str],
) -> tuple[
    list[dict[str, object]],
    list[tuple[Mapping[str, object], dict[str, object]]],
]:
    rendered_source_identities = {
        (_record_index(feature), _stable_feature_id(feature))
        for feature in context.features
        if isinstance(feature, Mapping)
        and _text(feature.get("type")).lower() == "source"
    }
    source_features = (
        context.biological_features
        if context.biological_features
        else context.features
    )
    candidates = [
        feature
        for feature in source_features
        if isinstance(feature, Mapping)
        and (
            _text(feature.get("type")).lower() != "source"
            or (
                (
                    _record_index(feature),
                    _stable_feature_id(feature),
                )
                in rendered_source_identities
            )
        )
    ]
    source_identities = [
        _validated_source_identity(
            feature,
            record_keys,
            description="Feature catalog source feature",
        )
        for feature in candidates
    ]
    identity_bases = [
        (_record_key(record_keys, record_index), _identity_base(feature))
        for feature, (record_index, _stable_id, _source_index) in zip(
            candidates,
            source_identities,
            strict=True,
        )
    ]
    collisions = Counter(identity_bases)
    collision_indexes: dict[tuple[str, str], list[int | None]] = defaultdict(list)
    for identity, (_record_index_value, _stable_id, feature_index) in zip(
        identity_bases,
        source_identities,
        strict=True,
    ):
        collision_indexes[identity].append(feature_index)
    for identity, indexes in collision_indexes.items():
        if len(indexes) > 1 and (
            any(index is None for index in indexes)
            or len(set(indexes)) != len(indexes)
        ):
            raise GbdrawError(
                "Feature catalog contains duplicate biological source identity "
                f"{identity!r}."
            )
    used: set[tuple[str, str]] = set()
    normalized: list[dict[str, object]] = []
    indexed: list[tuple[Mapping[str, object], dict[str, object]]] = []

    for feature, (record_key, identity_base), source_identity in zip(
        candidates,
        identity_bases,
        source_identities,
        strict=True,
    ):
        _record_index_value, stable_id, feature_index = source_identity
        biological_feature_id = identity_base
        if collisions[(record_key, identity_base)] > 1:
            biological_feature_id = f"{identity_base}~{feature_index}"
        if (record_key, biological_feature_id) in used:
            raise GbdrawError(
                "Feature catalog generated a duplicate canonical biological "
                "feature identity."
            )
        used.add((record_key, biological_feature_id))

        payload = dict(feature)
        for key in _BIOLOGICAL_ALIAS_KEYS | _ORTHOGROUP_FEATURE_KEYS:
            payload.pop(key, None)
        payload.pop("selector", None)
        for fasta_key in (
            "nucleotide_fasta",
            "nucleotideFasta",
            "amino_acid_fasta",
            "aminoAcidFasta",
        ):
            payload.pop(fasta_key, None)
        payload["biologicalFeatureId"] = biological_feature_id
        payload["recordKey"] = record_key
        for record_index_key in ("record_idx", "recordIndex", "record_index"):
            payload.pop(record_index_key, None)
        for feature_index_key in (
            "feature_index",
            "featureIndex",
            "source_feature_index",
            "sourceFeatureIndex",
        ):
            payload.pop(feature_index_key, None)
        if feature_index is not None:
            payload["sourceFeatureIndex"] = feature_index
        if stable_id and stable_id != biological_feature_id:
            payload["stableFeatureId"] = stable_id

        qualifiers = payload.get("qualifiers")
        if isinstance(qualifiers, Mapping):
            normalized_qualifiers = {
                str(key): list(value)
                if isinstance(value, Sequence) and not isinstance(value, (str, bytes))
                else value
                for key, value in qualifiers.items()
            }
            amino_acid_sequence = _text(
                payload.get("amino_acid_sequence")
                or payload.get("aminoAcidSequence")
            )
            translation = normalized_qualifiers.get("translation")
            if (
                amino_acid_sequence
                and isinstance(translation, list)
                and len(translation) == 1
                and translation[0] == amino_acid_sequence
            ):
                normalized_qualifiers.pop("translation", None)
                payload["translationFromAminoAcidSequence"] = True
            payload["qualifiers"] = normalized_qualifiers
            for field in (
                "protein_id",
                "locus_tag",
                "gene_id",
                "old_locus_tag",
                "gene",
                "product",
            ):
                if payload.get(field) == _first_qualifier_value(
                    normalized_qualifiers.get(field)
                ):
                    payload.pop(field, None)
            if payload.get("note") == _first_qualifier_value(
                normalized_qualifiers.get("note")
            )[:50]:
                payload.pop("note", None)

        protein_id = _text(
            payload.get("protein_id")
            or payload.get("proteinId")
            or (
                payload.get("qualifiers", {}).get("protein_id")
                if isinstance(payload.get("qualifiers"), Mapping)
                else None
            )
        )
        source_protein_id = _text(
            payload.get("source_protein_id")
            or payload.get("sourceProteinId")
        )
        if protein_id and source_protein_id == protein_id:
            payload.pop("source_protein_id", None)
            payload.pop("sourceProteinId", None)

        parts = _sequence(
            payload.get("location_parts") or payload.get("locationParts")
        )
        if len(parts) == 1 and isinstance(parts[0], Mapping):
            part = parts[0]
            start = _integer_or_none(payload.get("start"))
            end = _integer_or_none(payload.get("end"))
            if (
                start is not None
                and end is not None
                and _integer_or_none(part.get("start")) == start
                and _integer_or_none(part.get("end")) == end
                and _text(part.get("strand")) == _text(payload.get("strand"))
                and _text(part.get("display")) == f"{start + 1}..{end}"
            ):
                payload.pop("location_parts", None)
                payload.pop("locationParts", None)

        compact = dict(_compact_wire_value(payload) or {})
        normalized.append(compact)
        indexed.append((feature, compact))
    return normalized, indexed


class _BiologicalFeatureIndex:
    def __init__(
        self,
        indexed_features: Sequence[
            tuple[Mapping[str, object], Mapping[str, object]]
        ],
    ) -> None:
        self._by_record_feature_index: dict[
            tuple[int, int],
            list[Mapping[str, object]],
        ] = defaultdict(list)
        self._by_record_stable: dict[
            tuple[int, str],
            list[Mapping[str, object]],
        ] = defaultdict(list)
        self._by_canonical: dict[
            tuple[str, str],
            list[Mapping[str, object]],
        ] = defaultdict(list)
        for raw, normalized in indexed_features:
            record_index = _record_index(raw)
            stable_id = _stable_feature_id(raw)
            feature_index = _feature_index(raw)
            if record_index is not None and feature_index is not None:
                self._by_record_feature_index[
                    (record_index, feature_index)
                ].append(
                    normalized
                )
            if stable_id:
                if record_index is not None:
                    self._by_record_stable[(record_index, stable_id)].append(
                        normalized
                    )
            canonical = _feature_reference(normalized)
            if all(canonical):
                self._by_canonical[canonical].append(normalized)

    @staticmethod
    def _unique(
        candidates: Sequence[Mapping[str, object]],
    ) -> Mapping[str, object] | None:
        return candidates[0] if len(candidates) == 1 else None

    def resolve(
        self,
        *,
        record_index: int | None,
        stable_id: str = "",
        feature_index: int | None = None,
    ) -> Mapping[str, object] | None:
        if stable_id:
            if record_index is None:
                return None
            stable_candidates = self._by_record_stable.get(
                (record_index, stable_id),
                [],
            )
            if feature_index is None:
                return self._unique(stable_candidates)
            indexed = self._unique(
                self._by_record_feature_index.get(
                    (record_index, feature_index),
                    [],
                )
            )
            if indexed is None:
                return None
            return (
                indexed
                if any(candidate is indexed for candidate in stable_candidates)
                else None
            )
        if record_index is not None and feature_index is not None:
            return self._unique(
                self._by_record_feature_index.get(
                    (record_index, feature_index),
                    [],
                )
            )
        return None

    def resolve_canonical(
        self,
        *,
        record_key: str,
        biological_feature_id: str,
    ) -> Mapping[str, object] | None:
        if not record_key or not biological_feature_id:
            return None
        return self._unique(
            self._by_canonical.get(
                (record_key, biological_feature_id),
                [],
            )
        )


def _feature_reference(
    feature: Mapping[str, object] | None,
) -> tuple[str, str]:
    if feature is None:
        return "", ""
    return (
        _text(feature.get("recordKey")),
        _text(feature.get("biologicalFeatureId")),
    )


def _resolve_reference(
    feature: Mapping[str, object],
    record_keys: Sequence[str],
    biological_index: _BiologicalFeatureIndex,
    *,
    description: str,
) -> Mapping[str, object] | None:
    record_index = _record_index(feature)
    stable_id = _stable_feature_id(feature)
    feature_index = _feature_index(feature)
    record_key, record_key_valid = _record_key_alias(feature)
    biological_feature_id, biological_id_valid = (
        _biological_feature_id_alias(feature)
    )
    if not record_key_valid or not biological_id_valid:
        raise GbdrawError(f"{description} contains conflicting canonical aliases.")

    has_canonical_identity = bool(record_key or biological_feature_id)
    if has_canonical_identity and not (record_key and biological_feature_id):
        raise GbdrawError(
            f"{description} must supply record key and biological feature ID together."
        )
    has_source_identity = bool(
        record_index is not None or stable_id or feature_index is not None
    )
    references: list[Mapping[str, object]] = []
    if has_source_identity:
        record_index, stable_id, feature_index = _validated_source_identity(
            feature,
            record_keys,
            description=description,
        )
        source_reference = biological_index.resolve(
            record_index=record_index,
            stable_id=stable_id,
            feature_index=feature_index,
        )
        if source_reference is None:
            raise GbdrawError(
                f"{description} identity fields do not resolve to the same "
                "biological feature by a record-scoped stable feature ID or "
                "source feature index."
            )
        references.append(source_reference)
    if has_canonical_identity:
        canonical_reference = biological_index.resolve_canonical(
            record_key=record_key,
            biological_feature_id=biological_feature_id,
        )
        if canonical_reference is None:
            raise GbdrawError(
                f"{description} contains an unresolved canonical feature reference."
            )
        references.append(canonical_reference)
    if not references:
        return None
    reference = references[0]
    if any(candidate is not reference for candidate in references[1:]):
        raise GbdrawError(
            f"{description} identity fields do not resolve to the same biological feature."
        )
    return reference


def _normalized_rendered_features(
    rendered_features: Sequence[Mapping[str, object]],
    record_keys: Sequence[str],
    biological_index: _BiologicalFeatureIndex,
) -> tuple[list[dict[str, object]], dict[str, Mapping[str, object]]]:
    normalized: list[dict[str, object]] = []
    references_by_svg_id: dict[str, Mapping[str, object]] = {}
    for feature in rendered_features:
        svg_id = _rendered_svg_id(feature)
        reference = _resolve_reference(
            feature,
            record_keys,
            biological_index,
            description=f"Rendered feature {svg_id!r}",
        )
        record_key, biological_feature_id = _feature_reference(reference)
        if not svg_id or not record_key or not biological_feature_id:
            raise GbdrawError(
                "Rendered feature metadata does not resolve to one biological "
                f"feature: {svg_id or '<missing SVG ID>'}."
            )
        if svg_id in references_by_svg_id:
            raise GbdrawError(
                f"Rendered feature SVG ID {svg_id!r} is duplicated."
            )
        payload = {
            "svgId": svg_id,
            "recordKey": record_key,
            "biologicalFeatureId": biological_feature_id,
            "fillColor": _text(
                feature.get("fill_color") or feature.get("fillColor")
            ),
        }
        normalized.append(dict(_compact_wire_value(payload) or {}))
        references_by_svg_id[svg_id] = reference
    return normalized, references_by_svg_id


def _normalized_orthogroups(
    orthogroups: Sequence[Mapping[str, object]],
    record_keys: Sequence[str],
    biological_index: _BiologicalFeatureIndex,
    presentation_by_id: Mapping[str, tuple[str, str]] | None = None,
    default_presentation: tuple[str, str] | None = None,
) -> list[dict[str, object]]:
    normalized: list[dict[str, object]] = []
    known_group_ids: set[str] = set()
    member_specific_fields = (
        ("representative", ("representative",)),
        ("role", ("role",)),
        ("confidence", ("confidence",)),
        ("assignmentReason", ("assignmentReason", "assignment_reason")),
    )
    for group in orthogroups:
        if not isinstance(group, Mapping):
            continue
        group_id, group_id_valid = _text_alias_status(
            group,
            ("id", "orthogroupId", "orthogroup_id"),
        )
        if not group_id_valid or not group_id or group_id in known_group_ids:
            raise GbdrawError(
                "Feature catalog orthogroups require one consistent, unique ID."
            )
        known_group_ids.add(group_id)
        payload = dict(group)
        for key in ("id", "orthogroupId", "orthogroup_id"):
            payload.pop(key, None)
        payload["id"] = group_id
        explicit_presentation = _group_presentation(payload)
        for key in (
            "presentationScope",
            "presentation_scope",
            "collinearGroupScope",
            "collinear_group_scope",
            "groupKind",
            "group_kind",
        ):
            payload.pop(key, None)
        raw_members = _sequence(payload.pop("members", None))
        for collection_key, count_key in _ORTHOGROUP_COLLECTION_COUNTS:
            collection = _sequence(payload.pop(collection_key, None))
            if collection:
                payload[count_key] = len(collection)
        members: list[dict[str, object]] = []
        known_member_references: set[tuple[str, str]] = set()
        for member in raw_members:
            if not isinstance(member, Mapping):
                continue
            reference = _resolve_reference(
                member,
                record_keys,
                biological_index,
                description=(
                    f"Orthogroup {group_id!r} member"
                ),
            )
            record_key, biological_feature_id = _feature_reference(reference)
            if not record_key or not biological_feature_id:
                raise GbdrawError(
                    f"Orthogroup {group_id!r} member metadata must identify "
                    "exactly one biological feature by a record-scoped stable "
                    "feature ID or source feature index."
                )
            canonical_reference = (record_key, biological_feature_id)
            if canonical_reference in known_member_references:
                raise GbdrawError(
                    f"Orthogroup {group_id!r} contains duplicate biological members."
                )
            known_member_references.add(canonical_reference)
            compact_member = {
                "recordKey": record_key,
                "biologicalFeatureId": biological_feature_id,
                **{
                    output_key: next(
                        (
                            member.get(source_key)
                            for source_key in source_keys
                            if source_key in member
                        ),
                        None,
                    )
                    for output_key, source_keys in member_specific_fields
                },
            }
            members.append(dict(_compact_wire_value(compact_member) or {}))
        compact_group = dict(_compact_wire_value(payload) or {})
        presentations = {
            presentation
            for presentation in (
                explicit_presentation,
                default_presentation,
                (presentation_by_id or {}).get(group_id),
            )
            if presentation is not None
        }
        if len(presentations) > 1:
            raise GbdrawError(
                f"Orthogroup {group_id!r} has conflicting collinearity presentation metadata."
            )
        if presentations:
            presentation_scope, group_kind = next(iter(presentations))
            compact_group["presentationScope"] = presentation_scope
            compact_group["collinearGroupScope"] = presentation_scope
            compact_group["groupKind"] = group_kind
        compact_group["members"] = members
        normalized.append(compact_group)
    return normalized


def _orthogroup_presentation_by_id(
    matches: Sequence[Mapping[str, object]],
) -> dict[str, tuple[str, str]]:
    presentation_by_id: dict[str, tuple[str, str]] = {}
    for match in matches:
        scope_keys = (
            "collinear_group_scope",
            "collinearGroupScope",
            "group_scope",
            "groupScope",
        )
        kind_keys = ("group_kind", "groupKind")
        if not any(key in match for key in scope_keys + kind_keys):
            continue
        scope_values = [_text(match.get(key)) for key in scope_keys if key in match]
        kind_values = [_text(match.get(key)) for key in kind_keys if key in match]
        scope = next(iter(scope_values), "")
        kind = next(iter(kind_values), "")
        expected_kind = _PRESENTATION_KIND_BY_SCOPE.get(scope)
        if (
            not scope_values
            or not kind_values
            or any(not value for value in scope_values + kind_values)
            or len(set(scope_values)) != 1
            or len(set(kind_values)) != 1
            or expected_kind is None
            or kind != expected_kind
        ):
            raise GbdrawError(
                "Collinearity comparison metadata has an invalid presentation scope."
            )
        for raw_group_id in _sequence(match.get("orthogroup_ids")):
            group_id = _text(raw_group_id)
            if not group_id:
                continue
            presentation = (scope, kind)
            existing = presentation_by_id.get(group_id)
            if existing is not None and existing != presentation:
                raise GbdrawError(
                    f"Orthogroup {group_id!r} has conflicting collinearity presentation metadata."
                )
            presentation_by_id[group_id] = presentation
    return presentation_by_id


def _validate_match_orthogroup_references(
    matches: Sequence[Mapping[str, object]],
    orthogroups: Sequence[Mapping[str, object]],
) -> None:
    known_group_ids = {_text(group.get("id")) for group in orthogroups}
    for match in matches:
        match_id = _text(match.get("id")) or "<unnamed>"
        group_ids = [
            _text(group_id)
            for group_id in _sequence(match.get("orthogroup_ids"))
            if _text(group_id)
        ]
        missing = [group_id for group_id in group_ids if group_id not in known_group_ids]
        if missing:
            raise GbdrawError(
                f"Comparison match {match_id!r} references missing orthogroup "
                f"metadata: {', '.join(missing)}."
            )


def _match_endpoint_values(
    payload: Mapping[str, object],
    keys: Sequence[str],
    *,
    description: str,
    allow_duplicates: bool = False,
) -> list[str]:
    aliases: list[list[str]] = []
    for key in keys:
        if key not in payload:
            continue
        raw = _text(payload.get(key))
        values = [part.strip() for part in raw.split(";")]
        if (
            not raw
            or any(not value for value in values)
            or (
                not allow_duplicates
                and len(set(values)) != len(values)
            )
        ):
            raise GbdrawError(
                f"{description} contains invalid or duplicate endpoint IDs."
            )
        aliases.append(values)
    if not aliases:
        return []
    if any(values != aliases[0] for values in aliases[1:]):
        raise GbdrawError(f"{description} contains conflicting endpoint ID aliases.")
    return aliases[0]


def _match_endpoint_indexes(
    payload: Mapping[str, object],
    keys: Sequence[str],
    *,
    description: str,
) -> list[int]:
    values = _match_endpoint_values(
        payload,
        keys,
        description=description,
        allow_duplicates=True,
    )
    if any(not value.isascii() or not value.isdigit() for value in values):
        raise GbdrawError(
            f"{description} contains an invalid source feature index."
        )
    return [int(value) for value in values]


def _normalized_matches(
    matches: Sequence[Mapping[str, object]],
    rendered_references: Mapping[str, Mapping[str, object]],
    record_keys: Sequence[str],
    biological_index: _BiologicalFeatureIndex,
) -> list[dict[str, object]]:
    normalized: list[dict[str, object]] = []
    for match in matches:
        payload = dict(match)
        for role in ("query", "subject"):
            match_id = _text(payload.get("id")) or "<unnamed>"
            description = f"Comparison match {match_id!r} {role}"
            source_base = {
                "record_index": payload.get(f"{role}_record_index"),
                "recordIndex": payload.get(f"{role}RecordIndex"),
                "record_key": payload.get(f"{role}_record_key"),
                "recordKey": payload.get(f"{role}RecordKey"),
                "biological_feature_id": payload.get(
                    f"{role}_biological_feature_id"
                ),
                "biologicalFeatureId": payload.get(
                    f"{role}BiologicalFeatureId"
                ),
            }
            source_base = {
                key: value
                for key, value in source_base.items()
                if value is not None and str(value).strip()
            }
            stable_ids = _match_endpoint_values(
                payload,
                (
                    f"{role}_stable_feature_svg_id",
                    f"{role}StableFeatureSvgId",
                    f"{role}_stable_feature_id",
                    f"{role}StableFeatureId",
                ),
                description=description,
                allow_duplicates=True,
            )
            rendered_ids = _match_endpoint_values(
                payload,
                (
                    f"{role}_feature_svg_id",
                    f"{role}FeatureSvgId",
                    f"{role}_rendered_feature_svg_id",
                    f"{role}RenderedFeatureSvgId",
                ),
                description=description,
            )
            feature_indexes = _match_endpoint_indexes(
                payload,
                (
                    f"{role}_feature_index",
                    f"{role}FeatureIndex",
                    f"{role}_source_feature_index",
                    f"{role}SourceFeatureIndex",
                ),
                description=description,
            )
            canonical_keys = {
                "record_key",
                "recordKey",
                "biological_feature_id",
                "biologicalFeatureId",
            }
            endpoint_lengths = {
                len(values)
                for values in (stable_ids, rendered_ids, feature_indexes)
                if values
            }
            if len(endpoint_lengths) > 1:
                raise GbdrawError(
                    f"{description} endpoint identity arrays are not tuple-aligned."
                )
            endpoint_count = next(iter(endpoint_lengths), 0)
            if endpoint_count > 1 and any(
                key in source_base for key in canonical_keys
            ):
                raise GbdrawError(
                    f"{description} cannot combine plural endpoints with a "
                    "singular canonical identity."
                )
            source_references: list[Mapping[str, object]] = []
            if stable_ids or feature_indexes:
                for endpoint_index in range(endpoint_count):
                    endpoint_source = dict(source_base)
                    if stable_ids:
                        endpoint_source["stable_feature_id"] = stable_ids[
                            endpoint_index
                        ]
                    if feature_indexes:
                        endpoint_source["feature_index"] = feature_indexes[
                            endpoint_index
                        ]
                    reference = _resolve_reference(
                        endpoint_source,
                        record_keys,
                        biological_index,
                        description=description,
                    )
                    if reference is None:
                        raise GbdrawError(
                            f"{description} metadata does not resolve every stable endpoint."
                        )
                    source_references.append(reference)
            elif any(key in source_base for key in canonical_keys):
                reference = _resolve_reference(
                    source_base,
                    record_keys,
                    biological_index,
                    description=description,
                )
                if reference is None:
                    raise GbdrawError(
                        f"{description} metadata does not resolve to exactly one "
                        "biological feature."
                    )
                source_references.append(reference)
            rendered_endpoint_references: list[Mapping[str, object]] = []
            for rendered_id in rendered_ids:
                reference = rendered_references.get(rendered_id)
                if reference is None:
                    raise GbdrawError(
                        f"{description} rendered endpoint does not resolve to a "
                        "biological feature."
                    )
                rendered_endpoint_references.append(reference)
            if source_references and rendered_endpoint_references and [
                _feature_reference(reference)
                for reference in source_references
            ] != [
                _feature_reference(reference)
                for reference in rendered_endpoint_references
            ]:
                raise GbdrawError(
                    f"{description} identity fields do not resolve to the same "
                    "biological features: the rendered endpoints do not match "
                    "their canonical source identities."
                )
            references = source_references or rendered_endpoint_references
            canonical_references = [
                _feature_reference(reference)
                for reference in references
            ]
            if any(not all(reference) for reference in canonical_references):
                raise GbdrawError(
                    f"{description} contains an unresolved biological endpoint."
                )
            if len(set(canonical_references)) != len(canonical_references):
                raise GbdrawError(
                    f"{description} contains duplicate biological endpoints."
                )
            if canonical_references:
                payload[f"{role}FeatureReferences"] = [
                    {
                        "recordKey": record_key,
                        "biologicalFeatureId": biological_feature_id,
                    }
                    for record_key, biological_feature_id in canonical_references
                ]
            if len(canonical_references) == 1:
                record_key, biological_feature_id = canonical_references[0]
                payload[f"{role}RecordKey"] = record_key
                payload[f"{role}BiologicalFeatureId"] = biological_feature_id
            else:
                payload.pop(f"{role}RecordKey", None)
                payload.pop(f"{role}BiologicalFeatureId", None)
            for suffix in _MATCH_REDUNDANT_SUFFIXES:
                payload.pop(f"{role}_{suffix}", None)
            for suffix in (
                "FeatureIndex",
                "SourceFeatureIndex",
                "StableFeatureSvgId",
                "StableFeatureId",
                "FeatureSvgId",
                "RenderedFeatureSvgId",
            ):
                payload.pop(f"{role}{suffix}", None)
        normalized.append(dict(_compact_wire_value(payload) or {}))
    return normalized


def _normalized_annotations(
    root: ET.Element,
    context: InteractiveSvgContext,
) -> list[dict[str, object]]:
    context_by_key = {
        (
            _text(item.get("record_index")),
            _text(item.get("track_id")),
            _text(item.get("set_id")),
            _text(item.get("id")),
        ): dict(item)
        for item in context.annotations
        if isinstance(item, Mapping)
    }
    annotations: list[dict[str, object]] = []
    for element in root.iter():
        annotation_id = _text(element.get("data-gbdraw-annotation-id"))
        if not annotation_id:
            continue
        key = (
            _text(element.get("data-gbdraw-record-index")),
            _text(element.get("data-gbdraw-annotation-track-id")),
            _text(element.get("data-gbdraw-annotation-set-id")),
            annotation_id,
        )
        payload = context_by_key.get(
            key,
            context_by_key.get((key[0], "", key[2], key[3]), {}),
        ).copy()
        payload.update(
            {
                "dom_id": _text(element.get("id")),
                "id": annotation_id,
                "set_id": key[2],
                "track_id": key[1],
                "record_id": _text(element.get("data-gbdraw-record-id")),
                "record_index": (
                    int(key[0]) if key[0].isdigit() else key[0]
                ),
                "mark": _text(
                    element.get("data-gbdraw-annotation-mark")
                    or payload.get("mark")
                ),
                "label": _text(
                    element.get("data-gbdraw-annotation-label")
                    or payload.get("label")
                ),
            }
        )
        annotations.append(dict(_compact_wire_value(payload) or {}))
    return annotations


def _deduplicated_sequence_sources(
    matches: Sequence[Mapping[str, object]],
    context: InteractiveSvgContext,
) -> list[dict[str, object]]:
    selected = _sequence_sources_for_matches(matches, context.sequence_sources)
    normalized: list[dict[str, object]] = []
    by_key: dict[str, dict[str, object]] = {}
    for source in selected:
        if not isinstance(source, Mapping):
            continue
        payload = dict(_compact_wire_value(dict(source)) or {})
        key = _text(payload.get("key"))
        if not key:
            raise GbdrawError(
                "Feature catalog sequence sources require a non-empty key."
            )
        existing = by_key.get(key)
        if existing is not None:
            if existing != payload:
                raise GbdrawError(
                    f"Feature catalog sequence source key {key!r} is duplicated "
                    "with conflicting content."
                )
            continue
        by_key[key] = payload
        normalized.append(payload)
    return normalized


def _integer_or_none(value: object) -> int | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number) or not number.is_integer():
        return None
    return int(number)


def _reverse_complement_sequence(sequence: str) -> str:
    complements = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
        "U": "A",
        "R": "Y",
        "Y": "R",
        "S": "S",
        "W": "W",
        "K": "M",
        "M": "K",
        "B": "V",
        "D": "H",
        "H": "D",
        "V": "B",
        "N": "N",
        "-": "-",
    }
    return "".join(
        complements.get(base, "N")
        for base in reversed("".join(sequence.split()).upper())
    )


def _canonical_sequence_source(
    source: Mapping[str, object],
) -> str | None:
    sequence = source.get("sequence")
    if (
        not isinstance(sequence, str)
        or not sequence
        or _SEQUENCE_WHITESPACE_RE.search(sequence) is not None
    ):
        return None
    return sequence


def canonical_catalog_sequence_sources(
    sequence_sources: Sequence[object],
) -> tuple[str | None, ...]:
    """Validate whole-record sequence strings once per catalog item."""

    return tuple(
        _canonical_sequence_source(source)
        if isinstance(source, Mapping)
        else None
        for source in sequence_sources
    )


_UNVALIDATED_SOURCE_SEQUENCE = object()


def _reconstruct_nucleotide_sequence(
    feature: Mapping[str, object],
    source: Mapping[str, object],
    *,
    source_sequence: object = _UNVALIDATED_SOURCE_SEQUENCE,
) -> str | None:
    if source_sequence is _UNVALIDATED_SOURCE_SEQUENCE:
        source_sequence = _canonical_sequence_source(source)
    if not isinstance(source_sequence, str):
        return None
    raw_parts = _sequence(
        feature.get("location_parts") or feature.get("locationParts")
    )
    parts = raw_parts or [feature]
    reconstructed: list[str] = []
    for part in parts:
        if not isinstance(part, Mapping):
            return None
        start = _integer_or_none(part.get("start"))
        end = _integer_or_none(part.get("end"))
        if (
            start is None
            or end is None
            or start < 0
            or end < start
            or end > len(source_sequence)
        ):
            return None
        sequence = source_sequence[start:end]
        strand = _text(part.get("strand") or feature.get("strand"))
        try:
            reverse = strand == "-" or float(strand) < 0
        except ValueError:
            reverse = strand == "-"
        reconstructed.append(
            _reverse_complement_sequence(sequence) if reverse else sequence
        )
    return "".join(reconstructed)


def materialize_catalog_nucleotide_sequence(
    feature: Mapping[str, object],
    sequence_sources: Sequence[Mapping[str, object]],
    *,
    canonical_sources: Sequence[str | None] | None = None,
) -> str:
    """Return inline sequence data or reconstruct an explicit schema-3 reference."""

    if "sequenceSourceIndex" in feature:
        if (
            "nucleotide_sequence" in feature
            or "nucleotideSequence" in feature
        ):
            return ""
        source_index = feature.get("sequenceSourceIndex")
        if (
            type(source_index) is not int
            or not 0 <= source_index < len(sequence_sources)
        ):
            return ""
        if canonical_sources is None:
            canonical_sources = canonical_catalog_sequence_sources(
                sequence_sources
            )
        if source_index >= len(canonical_sources):
            return ""
        source = sequence_sources[source_index]
        if not isinstance(source, Mapping):
            return ""
        reconstructed = _reconstruct_nucleotide_sequence(
            feature,
            source,
            source_sequence=canonical_sources[source_index],
        )
        return reconstructed or ""

    inline = _text(
        feature.get("nucleotide_sequence")
        or feature.get("nucleotideSequence")
    )
    return inline


def _reference_reconstructible_nucleotide_sequences(
    features: Sequence[dict[str, object]],
    record_keys: Sequence[str],
    sequence_sources: Sequence[Mapping[str, object]],
) -> None:
    """Reference a record source only after exact sequence reconstruction."""

    record_index_by_key = {
        str(record_key): record_index
        for record_index, record_key in enumerate(record_keys)
    }
    canonical_sources = canonical_catalog_sequence_sources(sequence_sources)
    if any(source is None for source in canonical_sources):
        raise GbdrawError(
            "Feature catalog contains an invalid nucleotide sequence source."
        )
    for feature in features:
        inline = _text(
            feature.get("nucleotide_sequence")
            or feature.get("nucleotideSequence")
        )
        if not inline:
            continue
        record_index = record_index_by_key.get(_text(feature.get("recordKey")))
        if record_index is None:
            continue
        for source_index, source in enumerate(sequence_sources):
            if (
                _text(source.get("origin"))
                not in {"linear-record", "circular-reference"}
                or _integer_or_none(source.get("recordIndex")) != record_index
                or _reconstruct_nucleotide_sequence(
                    feature,
                    source,
                    source_sequence=canonical_sources[source_index],
                )
                != inline
            ):
                continue
            feature.pop("nucleotide_sequence", None)
            feature.pop("nucleotideSequence", None)
            feature["sequenceSourceIndex"] = source_index
            break


def build_feature_catalog_item(
    svg_source: str,
    context: InteractiveSvgContext,
    *,
    result_index: int,
    result_name: str,
) -> dict[str, Any]:
    """Normalize one prepared render context without reading source files."""

    if not isinstance(context, InteractiveSvgContext):
        raise GbdrawError("Feature catalog generation requires prepared context.")
    try:
        root = ET.fromstring(svg_source)
    except ET.ParseError as exc:
        raise GbdrawError(
            f"Malformed SVG source for feature catalog: {exc}"
        ) from exc
    if root.tag.rsplit("}", 1)[-1] != "svg":
        raise GbdrawError("Feature catalog generation expected an SVG root.")

    record_keys = _record_keys(context)
    biological_features, indexed_features = _normalized_biological_features(
        context,
        record_keys,
    )
    biological_index = _BiologicalFeatureIndex(indexed_features)
    rendered_payloads = _feature_payloads(root, context)
    rendered_entries = _collect_rendered_features(root)
    for rendered_payload in rendered_payloads:
        svg_id = _text(rendered_payload.get("svg_id"))
        if svg_id in rendered_entries:
            rendered_payload["fill_color"] = rendered_entries[svg_id].fill
    rendered_features, rendered_references = _normalized_rendered_features(
        rendered_payloads,
        record_keys,
        biological_index,
    )
    raw_matches = _match_payloads(root, rendered_payloads)
    orthogroups = list(context.orthogroups)
    known_group_ids: set[str] = set()
    for group in orthogroups:
        if isinstance(group, Mapping):
            group_id, _ = _text_alias_status(
                group,
                ("id", "orthogroupId", "orthogroup_id"),
            )
            known_group_ids.add(group_id)
    for match in raw_matches:
        if (
            _text(match.get("match_kind")) != "orthogroup"
            or "collinear_group_scope" in match
        ):
            continue
        group_scope = _text(match.pop("group_scope", None))
        group_kind = _text(match.pop("group_kind", None))
        if not group_scope and not group_kind:
            continue
        if group_scope not in {"cross_record", "record_local"} or group_kind != "orthogroup":
            raise GbdrawError(
                "Orthogroup comparison metadata has an invalid membership scope."
            )
        for group_id in _sequence(match.get("orthogroup_ids")):
            normalized_group_id = _text(group_id)
            if not normalized_group_id or normalized_group_id in known_group_ids:
                continue
            orthogroups.append(
                {"id": normalized_group_id, "scope": group_scope, "members": []}
            )
            known_group_ids.add(normalized_group_id)
    presentation_by_id = _orthogroup_presentation_by_id(raw_matches)
    sequence_sources = _deduplicated_sequence_sources(raw_matches, context)
    _reference_reconstructible_nucleotide_sequences(
        biological_features,
        record_keys,
        sequence_sources,
    )
    normalized_orthogroups = _normalized_orthogroups(
        orthogroups,
        record_keys,
        biological_index,
        presentation_by_id,
        _presentation_from_search_scope(
            context.collinearity_search_scope
        ),
    )
    normalized_matches = _normalized_matches(
        raw_matches,
        rendered_references,
        record_keys,
        biological_index,
    )
    _validate_match_orthogroup_references(
        normalized_matches,
        normalized_orthogroups,
    )
    item: dict[str, Any] = {
        "resultIndex": int(result_index),
        "resultName": str(result_name),
        "recordKeys": record_keys,
        "features": rendered_features,
        "biologicalFeatures": biological_features,
        "orthogroups": normalized_orthogroups,
        "annotations": _normalized_annotations(root, context),
        "comparisonMatches": normalized_matches,
    }
    if sequence_sources:
        item["sequenceSources"] = sequence_sources
    return item


def build_feature_catalog(
    items: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Wrap normalized logical-result items in the schema-3 envelope."""

    return {
        "schema": FEATURE_CATALOG_SCHEMA,
        "items": [dict(item) for item in items],
    }


def select_feature_catalog_item(
    catalog: Mapping[str, object],
    *,
    result_index: int,
    result_name: str,
) -> dict[str, object]:
    """Return one validated schema-3 item matched to a logical Result."""

    if (
        not isinstance(catalog, Mapping)
        or catalog.get("schema") != FEATURE_CATALOG_SCHEMA
    ):
        raise GbdrawError("Feature catalog must use schema 3.")
    items = catalog.get("items")
    if not isinstance(items, Sequence) or isinstance(items, (str, bytes)):
        raise GbdrawError("Feature catalog items must be an array.")
    matches = [
        item
        for item in items
        if isinstance(item, Mapping)
        and item.get("resultIndex") == result_index
        and _text(item.get("resultName")) == result_name
    ]
    if len(matches) != 1:
        raise GbdrawError(
            "Feature catalog must contain exactly one item for "
            f"Result {result_index} ({result_name})."
        )

    item = matches[0]
    required_arrays = (
        "recordKeys",
        "features",
        "biologicalFeatures",
        "orthogroups",
        "annotations",
        "comparisonMatches",
    )
    if any(
        not isinstance(item.get(key), Sequence)
        or isinstance(item.get(key), (str, bytes))
        for key in required_arrays
    ):
        raise GbdrawError("Feature catalog item is missing required arrays.")
    if (
        "sequenceSources" in item
        and (
            not isinstance(item.get("sequenceSources"), Sequence)
            or isinstance(item.get("sequenceSources"), (str, bytes))
        )
    ):
        raise GbdrawError("Feature catalog sequenceSources must be an array.")

    record_keys = [_text(value) for value in item["recordKeys"]]
    if any(not key for key in record_keys) or len(set(record_keys)) != len(
        record_keys
    ):
        raise GbdrawError(
            "Feature catalog record keys must be non-empty and unique."
        )
    known_features: set[tuple[str, str]] = set()
    biological_source_indexes: dict[
        tuple[str, str],
        list[int | None],
    ] = defaultdict(list)
    sequence_sources = item.get("sequenceSources")
    if not isinstance(sequence_sources, Sequence) or isinstance(
        sequence_sources,
        (str, bytes),
    ):
        sequence_sources = ()
    known_sequence_source_keys: set[str] = set()
    for source in sequence_sources:
        key = _text(source.get("key")) if isinstance(source, Mapping) else ""
        if not key or key in known_sequence_source_keys:
            raise GbdrawError(
                "Feature catalog sequence sources require unique non-empty keys."
            )
        known_sequence_source_keys.add(key)
    canonical_sources = canonical_catalog_sequence_sources(sequence_sources)
    if any(source is None for source in canonical_sources):
        raise GbdrawError(
            "Feature catalog contains an invalid nucleotide sequence "
            "source reference."
        )
    for feature in item["biologicalFeatures"]:
        if not isinstance(feature, Mapping):
            raise GbdrawError(
                "Feature catalog biologicalFeatures must contain objects."
            )
        reference = (
            _text(feature.get("recordKey")),
            _text(feature.get("biologicalFeatureId")),
        )
        if (
            not all(reference)
            or reference[0] not in record_keys
            or reference in known_features
        ):
            raise GbdrawError(
                "Feature catalog contains an invalid biological feature "
                "reference."
            )
        if "sourceFeatureIndex" in feature and (
            type(feature.get("sourceFeatureIndex")) is not int
            or feature["sourceFeatureIndex"] < 0
        ):
            raise GbdrawError(
                "Feature catalog contains an invalid source feature index."
            )
        stable_id = _text(feature.get("stableFeatureId")) or reference[1]
        source_feature_index = feature.get("sourceFeatureIndex")
        biological_source_indexes[(reference[0], stable_id)].append(
            source_feature_index
            if type(source_feature_index) is int
            else None
        )
        if "translationFromAminoAcidSequence" in feature:
            qualifiers = feature.get("qualifiers")
            amino_acid_sequence = (
                feature.get("amino_acid_sequence")
                if "amino_acid_sequence" in feature
                else feature.get("aminoAcidSequence")
            )
            if (
                feature.get("translationFromAminoAcidSequence") is not True
                or not isinstance(amino_acid_sequence, str)
                or not amino_acid_sequence.strip()
                or (
                    isinstance(qualifiers, Mapping)
                    and "translation" in qualifiers
                )
            ):
                raise GbdrawError(
                    "Feature catalog contains an invalid translation "
                    "source reference."
                )
        if "sequenceSourceIndex" in feature:
            source_index = feature.get("sequenceSourceIndex")
            source = (
                sequence_sources[source_index]
                if (
                    type(source_index) is int
                    and 0 <= source_index < len(sequence_sources)
                    and isinstance(sequence_sources[source_index], Mapping)
                )
                else None
            )
            if (
                "nucleotide_sequence" in feature
                or "nucleotideSequence" in feature
                or source is None
                or _text(source.get("origin"))
                not in {"linear-record", "circular-reference"}
                or type(source.get("recordIndex")) is not int
                or source.get("recordIndex") != record_keys.index(reference[0])
                or not _reconstruct_nucleotide_sequence(
                    feature,
                    source,
                    source_sequence=canonical_sources[source_index],
                )
            ):
                raise GbdrawError(
                    "Feature catalog contains an invalid nucleotide sequence "
                    "source reference."
                )
        known_features.add(reference)

    for source_identity, source_indexes in biological_source_indexes.items():
        if len(source_indexes) <= 1:
            continue
        if (
            any(index is None for index in source_indexes)
            or len(set(source_indexes)) != len(source_indexes)
        ):
            raise GbdrawError(
                "Feature catalog duplicate stable source identities require "
                "unique source feature indexes: "
                f"{source_identity!r}."
            )

    rendered_ids: set[str] = set()
    for feature in item["features"]:
        if not isinstance(feature, Mapping):
            raise GbdrawError(
                "Feature catalog features must contain objects."
            )
        svg_id = _text(feature.get("svgId"))
        reference = (
            _text(feature.get("recordKey")),
            _text(feature.get("biologicalFeatureId")),
        )
        if (
            not svg_id
            or svg_id in rendered_ids
            or reference not in known_features
        ):
            raise GbdrawError(
                "Feature catalog contains an invalid rendered feature "
                "reference."
            )
        rendered_ids.add(svg_id)

    known_group_ids: set[str] = set()
    presentation_by_group_id: dict[str, tuple[str, str] | None] = {}
    for group in item["orthogroups"]:
        if not isinstance(group, Mapping):
            raise GbdrawError("Feature catalog orthogroups must contain objects.")
        group_id, group_id_valid = _text_alias_status(
            group,
            ("id", "orthogroupId", "orthogroup_id"),
        )
        if not group_id_valid or not group_id or group_id in known_group_ids:
            raise GbdrawError(
                "Feature catalog orthogroups require one consistent, unique ID."
            )
        known_group_ids.add(group_id)
        presentation_by_group_id[group_id] = _group_presentation(group)
        members = group.get("members")
        if not isinstance(members, Sequence) or isinstance(
            members, (str, bytes)
        ):
            raise GbdrawError(
                "Feature catalog orthogroup members must be an array."
            )
        known_member_references: set[tuple[str, str]] = set()
        for member in members:
            if not isinstance(member, Mapping):
                raise GbdrawError(
                    "Feature catalog orthogroup members must contain objects."
                )
            reference = (
                _text(member.get("recordKey")),
                _text(member.get("biologicalFeatureId")),
            )
            if reference not in known_features:
                raise GbdrawError(
                    "Feature catalog contains an unresolved orthogroup member."
                )
            if reference in known_member_references:
                raise GbdrawError(
                    f"Orthogroup {group_id!r} contains duplicate biological members."
                )
            known_member_references.add(reference)

    known_match_ids: set[str] = set()
    for match in item["comparisonMatches"]:
        if not isinstance(match, Mapping):
            raise GbdrawError(
                "Feature catalog comparisonMatches must contain objects."
            )
        match_ids = {
            _text(match.get(key))
            for key in ("id", "matchId", "match_id")
            if key in match and _text(match.get(key))
        }
        if len(match_ids) != 1:
            raise GbdrawError(
                "Feature catalog comparison matches require one consistent match ID."
            )
        match_id = next(iter(match_ids))
        if match_id in known_match_ids:
            raise GbdrawError(
                "Feature catalog contains duplicate comparison match IDs."
            )
        known_match_ids.add(match_id)
        group_ids = [
            _text(group_id)
            for group_id in _sequence(match.get("orthogroup_ids"))
        ]
        if (
            any(not group_id for group_id in group_ids)
            or len(set(group_ids)) != len(group_ids)
            or any(group_id not in known_group_ids for group_id in group_ids)
        ):
            raise GbdrawError(
                "Feature catalog comparison match references invalid orthogroups."
            )
        match_presentations = _orthogroup_presentation_by_id((match,))
        if any(
            presentation_by_group_id.get(group_id) != presentation
            for group_id, presentation in match_presentations.items()
        ):
            raise GbdrawError(
                "Feature catalog comparison match has conflicting orthogroup presentation metadata."
            )
        for role in ("query", "subject"):
            singular_reference = (
                _text(match.get(f"{role}RecordKey")),
                _text(match.get(f"{role}BiologicalFeatureId")),
            )
            raw_references = match.get(f"{role}FeatureReferences")
            references: list[tuple[str, str]] = []
            if raw_references is not None:
                if not isinstance(raw_references, Sequence) or isinstance(
                    raw_references,
                    (str, bytes),
                ):
                    raise GbdrawError(
                        "Feature catalog comparison endpoint references must be an array."
                    )
                for raw_reference in raw_references:
                    if not isinstance(raw_reference, Mapping):
                        raise GbdrawError(
                            "Feature catalog comparison endpoint references must contain objects."
                        )
                    references.append(
                        (
                            _text(raw_reference.get("recordKey")),
                            _text(raw_reference.get("biologicalFeatureId")),
                        )
                    )
                if (
                    not references
                    or len(set(references)) != len(references)
                    or any(reference not in known_features for reference in references)
                ):
                    raise GbdrawError(
                        "Feature catalog contains invalid comparison endpoint references."
                    )
            if any(singular_reference) and singular_reference not in known_features:
                raise GbdrawError(
                    "Feature catalog contains an unresolved comparison match."
                )
            if (
                any(singular_reference)
                and (
                    not all(singular_reference)
                    or (
                        references
                        and references != [singular_reference]
                    )
                )
            ):
                raise GbdrawError(
                    "Feature catalog comparison endpoint aliases are conflicting."
                )
            if len(references) == 1 and singular_reference != references[0]:
                raise GbdrawError(
                    "Feature catalog singular comparison endpoint is missing or conflicting."
                )
            if len(references) > 1 and any(singular_reference):
                raise GbdrawError(
                    "Feature catalog plural comparison endpoints cannot use singular aliases."
                )
    return dict(item)


__all__ = [
    "FEATURE_CATALOG_SCHEMA",
    "build_feature_catalog",
    "build_feature_catalog_item",
    "canonical_catalog_sequence_sources",
    "materialize_catalog_nucleotide_sequence",
    "select_feature_catalog_item",
]
