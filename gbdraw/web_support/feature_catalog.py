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
    _collect_rendered_features,
    _compact_wire_value,
    _feature_record_index,
    _feature_payloads,
    _feature_stable_id,
    _first_text,
    _int_or_none,
    _match_payloads,
    _sequence_sources_for_matches,
)

FEATURE_CATALOG_SCHEMA = 3

_BIOLOGICAL_ALIAS_KEYS = {
    "id",
    "svg_id",
    "stable_svg_id",
    "stable_feature_id",
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


def _sequence(value: object | None) -> list[object]:
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        return list(value)
    return []


def _text(value: object | None) -> str:
    return _first_text(value)


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
    return _feature_record_index(feature)


def _stable_feature_id(feature: Mapping[str, object]) -> str:
    return _feature_stable_id(feature)


def _feature_index(feature: Mapping[str, object]) -> int | None:
    for key in ("feature_index", "featureIndex"):
        parsed = _int_or_none(feature.get(key))
        if parsed is not None:
            return parsed
    return None


def _record_keys(context: InteractiveSvgContext) -> list[str]:
    keys = [_text(value) for value in context.record_keys]
    if not keys:
        source_features = (
            context.biological_features
            if context.biological_features
            else context.features
        )
        record_indexes = [
            record_index
            for feature in source_features
            if isinstance(feature, Mapping)
            and (record_index := _record_index(feature)) is not None
        ]
        record_count = max(record_indexes, default=-1) + 1
        if source_features and record_count == 0:
            record_count = 1
        keys = [
            f"record-{record_index + 1}"
            for record_index in range(record_count)
        ]
    if any(not key for key in keys) or len(set(keys)) != len(keys):
        raise GbdrawError("Feature catalog record keys must be non-empty and unique.")
    return keys


def _record_key(record_keys: Sequence[str], record_index: int | None) -> str:
    if record_index is None:
        return str(record_keys[0]) if len(record_keys) == 1 else ""
    if 0 <= record_index < len(record_keys):
        return str(record_keys[record_index])
    return f"record-{record_index + 1}"


def _identity_base(
    feature: Mapping[str, object],
    *,
    ordinal: int,
) -> str:
    stable_id = _stable_feature_id(feature)
    if stable_id:
        return stable_id
    feature_index = _feature_index(feature)
    if feature_index is not None:
        return f"feature-{feature_index}"
    return f"feature-{ordinal}"


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
    identity_bases = [
        (
            _record_key(record_keys, _record_index(feature)),
            _identity_base(feature, ordinal=ordinal),
        )
        for ordinal, feature in enumerate(candidates)
    ]
    collisions = Counter(identity_bases)
    used: set[tuple[str, str]] = set()
    normalized: list[dict[str, object]] = []
    indexed: list[tuple[Mapping[str, object], dict[str, object]]] = []

    for ordinal, (feature, (record_key, identity_base)) in enumerate(
        zip(candidates, identity_bases, strict=True)
    ):
        feature_index = _feature_index(feature)
        biological_feature_id = identity_base
        if collisions[(record_key, identity_base)] > 1:
            discriminator = feature_index if feature_index is not None else ordinal
            biological_feature_id = f"{identity_base}~{discriminator}"
        unique_id = biological_feature_id
        duplicate_index = 2
        while (record_key, unique_id) in used:
            unique_id = f"{biological_feature_id}~{duplicate_index}"
            duplicate_index += 1
        used.add((record_key, unique_id))

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
        payload["biologicalFeatureId"] = unique_id
        payload["recordKey"] = record_key
        for record_index_key in ("record_idx", "recordIndex", "record_index"):
            payload.pop(record_index_key, None)
        payload.pop("feature_index", None)
        payload.pop("featureIndex", None)
        stable_id = _stable_feature_id(feature)
        if stable_id and stable_id != unique_id:
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
            Mapping[str, object],
        ] = {}
        self._by_record_stable: dict[
            tuple[int, str],
            list[Mapping[str, object]],
        ] = defaultdict(list)
        self._by_record_protein: dict[
            tuple[int, str],
            list[Mapping[str, object]],
        ] = defaultdict(list)
        self._by_stable: dict[str, list[Mapping[str, object]]] = defaultdict(list)
        for raw, normalized in indexed_features:
            record_index = _record_index(raw)
            stable_id = _stable_feature_id(raw)
            feature_index = _feature_index(raw)
            if record_index is not None and feature_index is not None:
                self._by_record_feature_index.setdefault(
                    (record_index, feature_index),
                    normalized,
                )
            if stable_id:
                self._by_stable[stable_id].append(normalized)
                if record_index is not None:
                    self._by_record_stable[(record_index, stable_id)].append(
                        normalized
                    )
            if record_index is None:
                continue
            for key in (
                "protein_id",
                "source_protein_id",
                "proteinId",
                "sourceProteinId",
            ):
                protein_id = _text(raw.get(key))
                if protein_id:
                    self._by_record_protein[
                        (record_index, protein_id)
                    ].append(normalized)

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
        protein_id: str = "",
    ) -> Mapping[str, object] | None:
        if record_index is not None and stable_id:
            resolved = self._unique(
                self._by_record_stable.get((record_index, stable_id), [])
            )
            if resolved is not None:
                return resolved
        if record_index is not None and feature_index is not None:
            resolved = self._by_record_feature_index.get(
                (record_index, feature_index)
            )
            if resolved is not None:
                return resolved
        if record_index is not None and protein_id:
            resolved = self._unique(
                self._by_record_protein.get((record_index, protein_id), [])
            )
            if resolved is not None:
                return resolved
        if stable_id:
            return self._unique(self._by_stable.get(stable_id, []))
        return None


def _feature_reference(
    feature: Mapping[str, object] | None,
) -> tuple[str, str]:
    if feature is None:
        return "", ""
    return (
        _text(feature.get("recordKey")),
        _text(feature.get("biologicalFeatureId")),
    )


def _normalized_rendered_features(
    rendered_features: Sequence[Mapping[str, object]],
    biological_index: _BiologicalFeatureIndex,
) -> tuple[list[dict[str, object]], dict[str, Mapping[str, object]]]:
    normalized: list[dict[str, object]] = []
    references_by_svg_id: dict[str, Mapping[str, object]] = {}
    for feature in rendered_features:
        svg_id = _text(
            feature.get("svg_id")
            or feature.get("rendered_feature_svg_id")
            or feature.get("renderedFeatureSvgId")
        )
        reference = biological_index.resolve(
            record_index=_record_index(feature),
            stable_id=_stable_feature_id(feature),
            feature_index=_feature_index(feature),
        )
        record_key, biological_feature_id = _feature_reference(reference)
        if not svg_id or not record_key or not biological_feature_id:
            raise GbdrawError(
                "Rendered feature metadata does not resolve to one biological "
                f"feature: {svg_id or '<missing SVG ID>'}."
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
    context: InteractiveSvgContext,
    biological_index: _BiologicalFeatureIndex,
) -> list[dict[str, object]]:
    normalized: list[dict[str, object]] = []
    member_specific_fields = (
        ("representative", ("representative",)),
        ("role", ("role",)),
        ("confidence", ("confidence",)),
        ("assignmentReason", ("assignmentReason", "assignment_reason")),
    )
    for group in context.orthogroups:
        if not isinstance(group, Mapping):
            continue
        payload = dict(group)
        raw_members = _sequence(payload.pop("members", None))
        for collection_key, count_key in _ORTHOGROUP_COLLECTION_COUNTS:
            collection = _sequence(payload.pop(collection_key, None))
            if collection:
                payload[count_key] = len(collection)
        members: list[dict[str, object]] = []
        for member in raw_members:
            if not isinstance(member, Mapping):
                continue
            reference = biological_index.resolve(
                record_index=_record_index(member),
                stable_id=_stable_feature_id(member),
                feature_index=_feature_index(member),
                protein_id=_text(
                    member.get("proteinId")
                    or member.get("protein_id")
                    or member.get("sourceProteinId")
                    or member.get("source_protein_id")
                ),
            )
            record_key, biological_feature_id = _feature_reference(reference)
            if not record_key or not biological_feature_id:
                raise GbdrawError(
                    "Orthogroup member metadata does not resolve to one "
                    "biological feature."
                )
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
        compact_group["members"] = members
        normalized.append(compact_group)
    return normalized


def _normalized_matches(
    matches: Sequence[Mapping[str, object]],
    rendered_references: Mapping[str, Mapping[str, object]],
    biological_index: _BiologicalFeatureIndex,
) -> list[dict[str, object]]:
    normalized: list[dict[str, object]] = []
    for match in matches:
        payload = dict(match)
        for role in ("query", "subject"):
            feature_svg_id = _text(payload.get(f"{role}_feature_svg_id"))
            record_index = _int_or_none(payload.get(f"{role}_record_index"))
            reference = rendered_references.get(feature_svg_id)
            if reference is None:
                reference = biological_index.resolve(
                    record_index=record_index,
                    stable_id=feature_svg_id,
                    protein_id=_text(payload.get(f"{role}_protein_id")),
                )
            record_key, biological_feature_id = _feature_reference(reference)
            if not record_key or not biological_feature_id:
                continue
            payload[f"{role}RecordKey"] = record_key
            payload[f"{role}BiologicalFeatureId"] = biological_feature_id
            for suffix in _MATCH_REDUNDANT_SUFFIXES:
                payload.pop(f"{role}_{suffix}", None)
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
    seen: set[str] = set()
    for source in selected:
        if not isinstance(source, Mapping):
            continue
        identity = _text(source.get("key")) or repr(sorted(source.items()))
        if identity in seen:
            continue
        seen.add(identity)
        normalized.append(dict(_compact_wire_value(dict(source)) or {}))
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
        biological_index,
    )
    raw_matches = _match_payloads(root, rendered_payloads, context.orthogroups)
    sequence_sources = _deduplicated_sequence_sources(raw_matches, context)
    _reference_reconstructible_nucleotide_sequences(
        biological_features,
        record_keys,
        sequence_sources,
    )
    item: dict[str, Any] = {
        "resultIndex": int(result_index),
        "resultName": str(result_name),
        "recordKeys": record_keys,
        "features": rendered_features,
        "biologicalFeatures": biological_features,
        "orthogroups": _normalized_orthogroups(context, biological_index),
        "annotations": _normalized_annotations(root, context),
        "comparisonMatches": _normalized_matches(
            raw_matches,
            rendered_references,
            biological_index,
        ),
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
    sequence_sources = item.get("sequenceSources")
    if not isinstance(sequence_sources, Sequence) or isinstance(
        sequence_sources,
        (str, bytes),
    ):
        sequence_sources = ()
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

    for group in item["orthogroups"]:
        if not isinstance(group, Mapping):
            raise GbdrawError("Feature catalog orthogroups must contain objects.")
        members = group.get("members")
        if not isinstance(members, Sequence) or isinstance(
            members, (str, bytes)
        ):
            raise GbdrawError(
                "Feature catalog orthogroup members must be an array."
            )
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

    for match in item["comparisonMatches"]:
        if not isinstance(match, Mapping):
            raise GbdrawError(
                "Feature catalog comparisonMatches must contain objects."
            )
        for role in ("query", "subject"):
            reference = (
                _text(match.get(f"{role}RecordKey")),
                _text(match.get(f"{role}BiologicalFeatureId")),
            )
            if any(reference) and reference not in known_features:
                raise GbdrawError(
                    "Feature catalog contains an unresolved comparison match."
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
