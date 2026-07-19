"""Standalone interactive SVG enrichment for CLI and API exports."""

from __future__ import annotations

import json
import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from importlib import resources
from typing import Literal, Mapping, Sequence

from gbdraw.exceptions import GbdrawError

SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
EV_NS = "http://www.w3.org/2001/xml-events"

INTERACTIVE_METADATA_ID = "gbdraw-interactive-feature-metadata"
INTERACTIVE_STYLE_ID = "gbdraw-interactive-feature-style"
INTERACTIVE_SCRIPT_ID = "gbdraw-interactive-feature-script"
INTERACTIVE_GLOW_FILTER_ID = "gbdraw-interactive-feature-glow"
INTERACTIVE_MATCH_GLOW_FILTER_ID = "gbdraw-interactive-feature-match-glow"
INTERACTIVE_SCHEMA = "gbdraw-interactive-feature-popup-v2"

_FEATURE_ELEMENT_SUFFIX_RE = re.compile(r"__(?:part|line)\d+$")
_FEATURE_CONNECTOR_SUFFIX_RE = re.compile(r"__line\d+$")
_FEATURE_RECORD_SUFFIX_RE = re.compile(r"_record_\d+$")
_ASSET_IDS = {
    INTERACTIVE_METADATA_ID,
    INTERACTIVE_STYLE_ID,
    INTERACTIVE_SCRIPT_ID,
    INTERACTIVE_GLOW_FILTER_ID,
    INTERACTIVE_MATCH_GLOW_FILTER_ID,
    "gbdraw-viewport-controls",
    "gbdraw-feature-search-controls",
    "gbdraw-feature-popup",
    "gbdraw-feature-hover-popup",
}

ET.register_namespace("", SVG_NS)
ET.register_namespace("xlink", XLINK_NS)
ET.register_namespace("ev", EV_NS)


@dataclass(frozen=True)
class InteractiveSvgContext:
    """Optional rich context for standalone interactive SVG export."""

    features: Sequence[Mapping[str, object]] = ()
    popup_mode: Literal["rich", "simple"] = "rich"
    orthogroups: Sequence[Mapping[str, object]] = ()
    legend_entries: Sequence[Mapping[str, object]] = ()
    current_colors: Mapping[str, str] = field(default_factory=dict)
    annotations: Sequence[Mapping[str, object]] = ()
    sequence_sources: Sequence[Mapping[str, object]] = ()


@dataclass
class _RenderedFeatureEntry:
    svg_id: str
    element: ET.Element
    fill: str


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _svg_tag(name: str) -> str:
    return f"{{{SVG_NS}}}{name}"


def _add_class_token(element: ET.Element, token: str) -> None:
    tokens = [item for item in str(element.get("class") or "").split() if item]
    if token not in tokens:
        tokens.append(token)
        element.set("class", " ".join(tokens))


def _remove_class_token(element: ET.Element, token: str) -> None:
    tokens = [item for item in str(element.get("class") or "").split() if item != token]
    if tokens:
        element.set("class", " ".join(tokens))
    else:
        element.attrib.pop("class", None)


def _compact_wire_value(value: object) -> object | None:
    if isinstance(value, Mapping):
        compact = {
            str(key): normalized
            for key, entry in value.items()
            if (normalized := _compact_wire_value(entry)) is not None
        }
        return compact or None
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        items = [normalized for entry in value if (normalized := _compact_wire_value(entry)) is not None]
        return items or None
    if value is None or value == "" or value is False:
        return None
    return value


def _normalize_feature_id(value: object | None) -> str:
    return _FEATURE_ELEMENT_SUFFIX_RE.sub("", str(value or "").strip())


def _element_feature_id(element: ET.Element) -> str:
    return _normalize_feature_id(
        element.get("data-gbdraw-feature-id") or element.get("id") or ""
    )


def _feature_svg_id_candidates(feature: Mapping[str, object]) -> list[str]:
    candidates: list[str] = []
    seen: set[str] = set()
    for value in (
        feature.get("svg_id"),
        feature.get("svgId"),
        feature.get("feature_svg_id"),
        feature.get("featureSvgId"),
        feature.get("stable_svg_id"),
        feature.get("stableSvgId"),
        feature.get("stable_feature_id"),
        feature.get("stableFeatureId"),
    ):
        text = _normalize_feature_id(value)
        if not text:
            continue
        for candidate in (text, _FEATURE_RECORD_SUFFIX_RE.sub("", text)):
            if candidate and candidate not in seen:
                seen.add(candidate)
                candidates.append(candidate)
    return candidates


def _is_feature_candidate(element: ET.Element) -> bool:
    if _local_name(element.tag) not in {"path", "polygon", "rect"}:
        return False
    element_id = str(element.get("id") or "")
    return bool(element.get("data-gbdraw-feature-id") or element_id.startswith("f"))


def _is_feature_fill_target(element: ET.Element) -> bool:
    explicit_part = str(element.get("data-gbdraw-feature-part") or "").strip()
    if explicit_part:
        return explicit_part == "block"
    return _FEATURE_CONNECTOR_SUFFIX_RE.search(str(element.get("id") or "")) is None


def _is_match_candidate(element: ET.Element) -> bool:
    if _local_name(element.tag) != "path":
        return False
    return any(
        element.get(attr)
        for attr in (
            "data-gbdraw-match-id",
            "data-gbdraw-pairwise-match-id",
            "data-match-kind",
            "data-pairwise-match-style",
        )
    )


def _collect_rendered_features(root: ET.Element) -> dict[str, _RenderedFeatureEntry]:
    entries: dict[str, _RenderedFeatureEntry] = {}
    for element in root.iter():
        if not _is_feature_candidate(element):
            continue
        svg_id = _element_feature_id(element)
        if not svg_id:
            continue
        existing = entries.get(svg_id)
        if existing is not None and (
            _is_feature_fill_target(existing.element) or not _is_feature_fill_target(element)
        ):
            continue
        entries[svg_id] = _RenderedFeatureEntry(
            svg_id=svg_id,
            element=element,
            fill=str(element.get("fill") or "#94a3b8"),
        )
    return entries


def _first_text(*values: object) -> str:
    for value in values:
        if isinstance(value, (list, tuple)):
            text = _first_text(*value)
            if text:
                return text
            continue
        text = str("" if value is None else value).strip()
        if text:
            return text
    return ""


def _strand_symbol(value: object) -> str:
    text = _first_text(value)
    if text == "1":
        return "+"
    if text == "-1":
        return "-"
    return text


def _normalize_string_array(value: object | None) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        return [str(item) for item in value if item is not None]
    if value is None or value == "":
        return []
    return [str(value)]


def _normalize_qualifier_map(qualifiers: object | None) -> dict[str, list[str]]:
    if not isinstance(qualifiers, Mapping):
        return {}
    normalized: dict[str, list[str]] = {}
    for key, value in qualifiers.items():
        normalized_key = str(key or "").strip()
        values = _normalize_string_array(value)
        if normalized_key and values:
            normalized[normalized_key] = values
    return normalized


def _first_qualifier_value(feature: Mapping[str, object] | None, key: str) -> str:
    if feature is None:
        return ""
    qualifiers = feature.get("qualifiers")
    if not isinstance(qualifiers, Mapping):
        return ""
    values = _normalize_string_array(qualifiers.get(key))
    return next((value for value in values if value.strip()), "")


def _number_or_none(value: object) -> float | int | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number):
        return None
    return int(number) if number.is_integer() else number


def _int_or_none(value: object) -> int | None:
    number = _number_or_none(value)
    return int(number) if isinstance(number, (int, float)) else None


def _build_feature_location(feature: Mapping[str, object]) -> str:
    start = _number_or_none(feature.get("start"))
    end = _number_or_none(feature.get("end"))
    start_text = str(int(start) + 1) if isinstance(start, (int, float)) else _first_text(feature.get("start"))
    end_text = str(int(end)) if isinstance(end, (int, float)) else _first_text(feature.get("end"))
    strand = _first_text(feature.get("strand"))
    range_text = f"{start_text}..{end_text}"
    return f"{range_text} ({strand})" if strand else range_text


def _get_feature_label(feature: Mapping[str, object]) -> str:
    for candidate in (
        feature.get("label"),
        feature.get("gene"),
        feature.get("locus_tag"),
        _first_qualifier_value(feature, "gene"),
        _first_qualifier_value(feature, "locus_tag"),
        _first_qualifier_value(feature, "product"),
        feature.get("product"),
        feature.get("type"),
        feature.get("svg_id"),
    ):
        text = _first_text(candidate)
        if text:
            return text
    return "Feature"


def _get_search_labels(
    feature: Mapping[str, object],
    fallback_label: str,
    display_label: str,
) -> list[str]:
    labels = [
        display_label,
        fallback_label,
        feature.get("label"),
        feature.get("gene"),
        feature.get("locus_tag"),
        feature.get("product"),
        _first_qualifier_value(feature, "gene"),
        _first_qualifier_value(feature, "locus_tag"),
        _first_qualifier_value(feature, "product"),
        feature.get("svg_id"),
    ]
    out: list[str] = []
    seen: set[str] = set()
    for value in labels:
        text = _first_text(value)
        key = text.lower()
        if not text or key in seen:
            continue
        seen.add(key)
        out.append(text)
    return out


def _normalize_location_parts(parts: object | None) -> list[dict[str, object]]:
    if not isinstance(parts, Sequence) or isinstance(parts, (str, bytes)):
        return []
    normalized: list[dict[str, object]] = []
    for part in parts:
        if not isinstance(part, Mapping):
            continue
        start = _number_or_none(part.get("start"))
        end = _number_or_none(part.get("end"))
        display = _first_text(part.get("display"))
        if not display:
            start_text = str(int(start) + 1) if isinstance(start, (int, float)) else ""
            end_text = str(int(end)) if isinstance(end, (int, float)) else ""
            display = f"{start_text}..{end_text}"
        if not display or display == "..":
            continue
        normalized.append(
            {
                "start": int(start) if isinstance(start, (int, float)) else None,
                "end": int(end) if isinstance(end, (int, float)) else None,
                "strand": _first_text(part.get("strand")),
                "display": display,
            }
        )
    return normalized


def _normalize_header_text(value: object) -> str:
    return re.sub(r"\s+", " ", _first_text(value)).strip()


def _normalize_fasta_id(value: object) -> str:
    return re.sub(r"\s+", "_", _normalize_header_text(value).lstrip(">"))


def _sequence_text(value: object) -> str:
    return re.sub(r"\s+", "", _first_text(value)).strip()


def _wrap_fasta_sequence(sequence: object, width: int = 60) -> str:
    text = _sequence_text(sequence)
    if not text:
        return ""
    width = max(1, int(width or 60))
    return "\n".join(text[index : index + width] for index in range(0, len(text), width))


def _format_fasta_entry(*, id_text: str, description: str, sequence: object) -> str:
    normalized_id = _normalize_fasta_id(id_text) or "sequence"
    wrapped = _wrap_fasta_sequence(sequence)
    if not wrapped:
        return ""
    description = _normalize_header_text(description)
    header = f">{normalized_id} {description}" if description else f">{normalized_id}"
    return f"{header}\n{wrapped}"


def _camelize_key(key: str) -> str:
    return re.sub(r"_([a-z])", lambda match: match.group(1).upper(), key)


def _get_feature_values(feature: Mapping[str, object], key: str) -> list[str]:
    for candidate_key in (key, _camelize_key(key)):
        if candidate_key in feature:
            values = _normalize_string_array(feature.get(candidate_key))
            if values:
                return values
    normalized_key = key.lower()
    qualifiers = feature.get("qualifiers")
    if isinstance(qualifiers, Mapping):
        for candidate_key, value in qualifiers.items():
            if str(candidate_key).lower() == normalized_key:
                return _normalize_string_array(value)
    return []


def _first_feature_text(feature: Mapping[str, object], keys: Sequence[str]) -> str:
    for key in keys:
        for value in _get_feature_values(feature, key):
            text = _normalize_header_text(value)
            if text:
                return text
    return ""


def _feature_description(feature: Mapping[str, object]) -> str:
    label = _first_feature_text(
        feature,
        ("product", "protein", "gene", "locus_tag", "name", "standard_name", "note", "type"),
    )
    organism = _first_feature_text(feature, ("organism", "source_organism", "record_organism"))
    if not organism:
        return label
    if not label:
        return f"[{organism}]"
    return label if re.search(r"\[[^\]]+\]\s*$", label) else f"{label} [{organism}]"


def _format_location_token(start_raw: object, end_raw: object, strand_raw: object) -> str:
    start = _int_or_none(start_raw)
    end = _int_or_none(end_raw)
    if start is None or end is None:
        return ""
    start_one_based = max(1, start + 1)
    end_one_based = max(start_one_based, end)
    return (
        f"c{end_one_based}-{start_one_based}"
        if _first_text(strand_raw) == "-"
        else f"{start_one_based}-{end_one_based}"
    )


def _feature_location_token(feature: Mapping[str, object]) -> str:
    parts = feature.get("location_parts")
    if isinstance(parts, Sequence) and not isinstance(parts, (str, bytes)):
        tokens = [
            _format_location_token(
                part.get("start"),
                part.get("end"),
                _first_text(part.get("strand"), feature.get("strand")),
            )
            for part in parts
            if isinstance(part, Mapping)
        ]
        tokens = [token for token in tokens if token]
        if len(tokens) == 1:
            return tokens[0]
        if len(tokens) > 1:
            return f"join({','.join(tokens)})"
    direct = _format_location_token(feature.get("start"), feature.get("end"), feature.get("strand"))
    if direct:
        return direct
    return _normalize_fasta_id(feature.get("location")) or "feature"


def _build_feature_sequence_fastas(feature: Mapping[str, object]) -> dict[str, str]:
    description = _feature_description(feature)
    record_id = _first_feature_text(feature, ("record_id", "recordId")) or "record"
    nucleotide_id = f"{record_id}:{_feature_location_token(feature)}"
    protein_id = _first_feature_text(
        feature,
        (
            "source_protein_id",
            "protein_id",
            "locus_tag",
            "gene",
            "gene_id",
            "old_locus_tag",
            "gff_id",
            "transcript_id",
            "name",
        ),
    ) or nucleotide_id
    return {
        "nucleotideFasta": _format_fasta_entry(
            id_text=nucleotide_id,
            description=description,
            sequence=_first_text(feature.get("nucleotide_sequence"), feature.get("nucleotideSequence")),
        ),
        "aminoAcidFasta": _format_fasta_entry(
            id_text=protein_id,
            description=description,
            sequence=_first_text(feature.get("amino_acid_sequence"), feature.get("aminoAcidSequence")),
        ),
    }


def _feature_orthogroup_entry(feature: Mapping[str, object]) -> dict[str, object] | None:
    orthogroup_id = _first_text(feature.get("orthogroupId"), feature.get("orthogroup_id"))
    if not orthogroup_id:
        return None
    return {
        "orthogroupId": orthogroup_id,
        "orthogroupMemberCount": _int_or_none(
            _first_text(feature.get("orthogroupMemberCount"), feature.get("orthogroup_member_count"))
        )
        or 0,
        "orthogroupRecordCoverage": _int_or_none(
            _first_text(feature.get("orthogroupRecordCoverage"), feature.get("orthogroup_record_coverage"))
        )
        or 0,
        "proteinId": _first_text(feature.get("proteinId"), feature.get("protein_id")),
        "sourceProteinId": _first_text(feature.get("sourceProteinId"), feature.get("source_protein_id")),
        "orthogroupRepresentative": bool(
            feature.get("orthogroupRepresentative") or feature.get("orthogroup_representative")
        ),
        "orthogroupMember": feature.get("orthogroupMember") or feature.get("orthogroup_member"),
    }


def _normalize_orthogroup_member(member: object | None) -> dict[str, object] | None:
    if not isinstance(member, Mapping):
        return None
    record_index = _int_or_none(_first_text(member.get("recordIndex"), member.get("record_index")))
    feature_index = _number_or_none(_first_text(member.get("featureIndex"), member.get("feature_index")))
    start = _number_or_none(member.get("start"))
    end = _number_or_none(member.get("end"))
    return {
        "orthogroup_id": _first_text(member.get("orthogroupId"), member.get("orthogroup_id")),
        "protein_id": _first_text(member.get("proteinId"), member.get("protein_id")),
        "source_protein_id": _first_text(member.get("sourceProteinId"), member.get("source_protein_id")),
        "record_index": record_index,
        "record_id": _first_text(member.get("recordId"), member.get("record_id")),
        "feature_index": feature_index,
        "label": _first_text(member.get("label")),
        "feature_svg_id": _first_text(member.get("featureSvgId"), member.get("feature_svg_id")),
        "start": start,
        "end": end,
        "strand": _strand_symbol(member.get("strand")),
        "representative": bool(member.get("representative")),
        "gene": _first_text(member.get("gene")),
        "locus_tag": _first_text(member.get("locusTag"), member.get("locus_tag")),
        "gene_id": _first_text(member.get("geneId"), member.get("gene_id")),
        "old_locus_tag": _first_text(member.get("oldLocusTag"), member.get("old_locus_tag")),
        "product": _first_text(member.get("product")),
        "note": _first_text(member.get("note")),
    }


def _normalize_orthogroup_candidate(candidate: object | None) -> dict[str, object] | None:
    if not isinstance(candidate, Mapping):
        return None
    return {
        "text": _first_text(candidate.get("text")),
        "source": _first_text(candidate.get("source")),
        "member_count": _int_or_none(
            _first_text(candidate.get("memberCount"), candidate.get("member_count"))
        )
        or 0,
        "record_coverage_count": _int_or_none(
            _first_text(candidate.get("recordCoverageCount"), candidate.get("record_coverage_count"))
        )
        or 0,
        "representative_count": _int_or_none(
            _first_text(candidate.get("representativeCount"), candidate.get("representative_count"))
        )
        or 0,
        "score": _number_or_none(candidate.get("score")) or 0,
    }


def _orthogroup_payloads(
    features: Sequence[Mapping[str, object]],
    context: InteractiveSvgContext,
) -> list[dict[str, object]]:
    needed_ids = {
        _first_text(feature.get("orthogroup_id"), feature.get("orthogroupId"))
        for feature in features
        if _first_text(feature.get("orthogroup_id"), feature.get("orthogroupId"))
    }
    if not needed_ids:
        return []
    payloads: list[dict[str, object]] = []
    for group in context.orthogroups or []:
        if not isinstance(group, Mapping):
            continue
        group_id = _first_text(group.get("id"))
        if group_id not in needed_ids:
            continue
        raw_members = group.get("members")
        raw_member_list = (
            raw_members
            if isinstance(raw_members, Sequence) and not isinstance(raw_members, (str, bytes))
            else []
        )
        members = [
            normalized
            for normalized in (
                _normalize_orthogroup_member(member)
                for member in raw_member_list
            )
            if normalized is not None
        ]
        member_count = _int_or_none(_first_text(group.get("member_count"), group.get("memberCount"))) or len(members)
        record_coverage_fallback = len(
            {
                member.get("record_index")
                for member in members
                if isinstance(member.get("record_index"), int)
            }
        )
        record_coverage_count = (
            _int_or_none(
                _first_text(group.get("record_coverage_count"), group.get("recordCoverageCount"))
            )
            or record_coverage_fallback
            or 0
        )
        raw_candidates = group.get("nameCandidates") or group.get("name_candidates")
        raw_candidate_list = (
            raw_candidates
            if isinstance(raw_candidates, Sequence) and not isinstance(raw_candidates, (str, bytes))
            else []
        )
        candidates = [
            normalized
            for normalized in (
                _normalize_orthogroup_candidate(candidate)
                for candidate in raw_candidate_list
            )
            if normalized is not None
        ]
        display_name = _first_text(
            group.get("display_name"),
            group.get("displayName"),
            group.get("name"),
            group_id,
        )
        payloads.append(
            {
                "id": group_id,
                "name": _first_text(group.get("name")),
                "display_name": display_name,
                "description": _first_text(group.get("description")),
                "name_confidence": _first_text(
                    group.get("nameConfidence"), group.get("name_confidence")
                ),
                "name_candidates": candidates,
                "member_count": member_count,
                "record_coverage_count": record_coverage_count,
                "members": members,
            }
        )
    return payloads


def _fallback_feature_payload(svg_id: str, entry: _RenderedFeatureEntry) -> dict[str, object]:
    fill_color = entry.fill
    stable_svg_id = _first_text(entry.element.get("data-gbdraw-stable-feature-id"), svg_id)
    label = _first_text(entry.element.get("data-label"), entry.element.get("id"), svg_id)
    search_labels = []
    for value in (label, svg_id):
        if value and value not in search_labels:
            search_labels.append(value)
    return {
        "svg_id": svg_id,
        "stable_svg_id": stable_svg_id,
        "label": label,
        "display_label": label,
        "search_labels": search_labels,
        "record_id": "",
        "record_idx": None,
        "type": "Feature",
        "start": None,
        "end": None,
        "strand": "",
        "location": "",
        "locus_tag": "",
        "gene_id": "",
        "old_locus_tag": "",
        "fill_color": fill_color,
        "orthogroup_id": "",
        "orthogroup_member_count": 0,
        "orthogroup_record_coverage": 0,
        "protein_id": "",
        "source_protein_id": "",
        "orthogroup_representative": False,
        "qualifiers": {},
        "location_parts": [],
        "nucleotide_sequence": "",
        "amino_acid_sequence": "",
        "sequence_warnings": [],
        "nucleotide_fasta": "",
        "amino_acid_fasta": "",
        "orthogroup_member": None,
    }


def _feature_payloads(
    root: ET.Element,
    context: InteractiveSvgContext,
) -> list[dict[str, object]]:
    rendered = _collect_rendered_features(root)
    if not rendered:
        return []

    payloads: list[dict[str, object]] = []
    seen: set[str] = set()
    popup_mode: Literal["rich", "simple"] = (
        "simple" if context.popup_mode == "simple" else "rich"
    )
    for feature in context.features:
        svg_id = next(
            (
                candidate
                for candidate in _feature_svg_id_candidates(feature)
                if candidate in rendered and candidate not in seen
            ),
            "",
        )
        if not svg_id or svg_id not in rendered or svg_id in seen:
            continue
        seen.add(svg_id)
        fallback_label = _get_feature_label(feature)
        display_label = fallback_label
        orthogroup_entry = _feature_orthogroup_entry(feature)
        orthogroup_member = _normalize_orthogroup_member(
            orthogroup_entry.get("orthogroupMember")
            if isinstance(orthogroup_entry, Mapping)
            else None
        )
        qualifiers = _normalize_qualifier_map(feature.get("qualifiers"))
        amino_acid_sequence = _first_text(
            feature.get("amino_acid_sequence"),
            feature.get("aminoAcidSequence"),
        )
        translation = next(
            (value for value in _normalize_string_array(qualifiers.get("translation")) if value.strip()),
            "",
        )
        payload = {
            "svg_id": svg_id,
            "stable_svg_id": _first_text(
                feature.get("stable_svg_id"),
                feature.get("stableFeatureSvgId"),
                feature.get("stable_feature_id"),
                rendered[svg_id].element.get("data-gbdraw-stable-feature-id"),
                svg_id,
            ),
            "label": fallback_label,
            "display_label": display_label,
            "search_labels": _get_search_labels(feature, fallback_label, display_label),
            "record_id": _first_text(feature.get("record_id")),
            "record_idx": _int_or_none(feature.get("record_idx")),
            "type": _first_text(feature.get("type")),
            "start": _int_or_none(feature.get("start")),
            "end": _int_or_none(feature.get("end")),
            "strand": _first_text(feature.get("strand")),
            "location": _build_feature_location(feature),
            "locus_tag": _first_text(feature.get("locus_tag"), feature.get("locusTag")),
            "gene_id": _first_text(feature.get("gene_id"), feature.get("geneId")),
            "old_locus_tag": _first_text(feature.get("old_locus_tag"), feature.get("oldLocusTag")),
            "orthogroup_id": _first_text(
                orthogroup_entry.get("orthogroupId")
                if isinstance(orthogroup_entry, Mapping)
                else ""
            ),
            "orthogroup_member_count": _int_or_none(
                orthogroup_entry.get("orthogroupMemberCount")
                if isinstance(orthogroup_entry, Mapping)
                else None
            )
            or 0,
            "orthogroup_record_coverage": _int_or_none(
                orthogroup_entry.get("orthogroupRecordCoverage")
                if isinstance(orthogroup_entry, Mapping)
                else None
            )
            or 0,
            "protein_id": _first_text(
                orthogroup_entry.get("proteinId")
                if isinstance(orthogroup_entry, Mapping)
                else "",
                feature.get("proteinId"),
                feature.get("protein_id"),
            ),
            "source_protein_id": _first_text(
                orthogroup_entry.get("sourceProteinId")
                if isinstance(orthogroup_entry, Mapping)
                else "",
                feature.get("sourceProteinId"),
                feature.get("source_protein_id"),
                _first_qualifier_value(feature, "protein_id"),
            ),
            "orthogroup_representative": bool(
                orthogroup_entry.get("orthogroupRepresentative")
                if isinstance(orthogroup_entry, Mapping)
                else False
            ),
        }
        if popup_mode == "rich":
            payload.update(
                {
                    "qualifiers": qualifiers,
                    "location_parts": _normalize_location_parts(feature.get("location_parts")),
                    "nucleotide_sequence": _first_text(
                        feature.get("nucleotide_sequence"),
                        feature.get("nucleotideSequence"),
                    ),
                    "amino_acid_sequence": (
                        amino_acid_sequence if amino_acid_sequence != translation else ""
                    ),
                    "sequence_warnings": _normalize_string_array(feature.get("sequence_warnings")),
                    "orthogroup_member": orthogroup_member,
                }
            )
        payloads.append(dict(_compact_wire_value(payload) or {}))

    for svg_id, entry in rendered.items():
        if svg_id in seen:
            continue
        seen.add(svg_id)
        payloads.append(dict(_compact_wire_value(_fallback_feature_payload(svg_id, entry)) or {}))

    feature_ids = {str(feature.get("svg_id") or "").strip() for feature in payloads}
    for element in root.iter():
        if not _is_feature_candidate(element):
            continue
        if _element_feature_id(element) not in feature_ids:
            continue
        element.set("data-gbdraw-interactive-feature", "true")
        _add_class_token(element, "gbdraw-interactive-feature")
    return payloads


def _match_kind(element: ET.Element) -> str:
    value = _first_text(element.get("data-match-kind")).lower()
    if value in {"pairwise", "orthogroup", "collinear", "homology"}:
        return value
    if element.get("data-collinearity-block-id"):
        return "collinear"
    if element.get("data-orthogroup-id"):
        return "orthogroup"
    return "pairwise"


_MATCH_TITLES = {
    "pairwise": "Pairwise match",
    "orthogroup": "Orthogroup match",
    "collinear": "Collinearity block",
    "homology": "Homology ring match",
}


def _add_match_row(rows: list[list[str]], label: str, value: object) -> None:
    text = _first_text(value)
    if text:
        rows.append([label, text])


def _interval_text(start: object, end: object) -> str:
    start_text = _first_text(start)
    end_text = _first_text(end)
    if start_text and end_text:
        return f"{start_text}..{end_text}"
    return start_text or end_text


def _split_metadata_values(value: object | None) -> list[str]:
    return [
        part.strip()
        for part in _first_text(value).split(";")
        if part.strip()
    ]


def _unique_metadata_values(value: object | None) -> list[str]:
    seen: set[str] = set()
    values: list[str] = []
    for item in _split_metadata_values(value):
        if item in seen:
            continue
        seen.add(item)
        values.append(item)
    return values


_GENERATED_PROTEIN_ID_RE = re.compile(
    r"^(?:gbd_r\d+_cds\d+|p_.+_\d+_\d+_-?\d+_[0-9a-f]{12}(?:_\d+)?)$",
    re.IGNORECASE,
)
_GENERATED_UNIT_ID_RE = re.compile(r"^gbd_r\d+_unit\d+$", re.IGNORECASE)


def _is_internal_display_id(value: object) -> bool:
    text = _first_text(value)
    return bool(text and (_GENERATED_PROTEIN_ID_RE.match(text) or _GENERATED_UNIT_ID_RE.match(text)))


def _add_unique_display_text(values: list[str], value: object) -> None:
    text = _first_text(value)
    if not text or _is_internal_display_id(text) or text in values:
        return
    values.append(text)


def _first_non_internal_display_text(*values: object) -> str:
    for value in values:
        text = _first_text(value)
        if text and not _is_internal_display_id(text):
            return text
    return ""


def _feature_lookup(features: Sequence[Mapping[str, object]]) -> dict[str, Mapping[str, object]]:
    lookup: dict[str, Mapping[str, object]] = {}
    for feature in features:
        svg_id = _first_text(feature.get("svg_id"))
        if svg_id and svg_id not in lookup:
            lookup[svg_id] = feature
    return lookup


def _orthogroup_lookup(
    orthogroups: Sequence[Mapping[str, object]],
) -> dict[str, Mapping[str, object]]:
    lookup: dict[str, Mapping[str, object]] = {}
    for group in orthogroups:
        group_id = _first_text(group.get("id"), group.get("orthogroupId"), group.get("orthogroup_id"))
        if group_id and group_id not in lookup:
            lookup[group_id] = group
    return lookup


def _member_feature_svg_id(member: Mapping[str, object]) -> str:
    return _first_text(member.get("featureSvgId"), member.get("feature_svg_id"))


def _feature_orthogroup_id(feature: Mapping[str, object] | None) -> str:
    if feature is None:
        return ""
    member = feature.get("orthogroup_member") or feature.get("orthogroupMember")
    member = member if isinstance(member, Mapping) else {}
    return _first_text(
        feature.get("orthogroup_id"),
        feature.get("orthogroupId"),
        member.get("orthogroupId"),
        member.get("orthogroup_id"),
    )


def _feature_product(feature: Mapping[str, object] | None) -> str:
    if feature is None:
        return ""
    member = feature.get("orthogroup_member") or feature.get("orthogroupMember")
    member = member if isinstance(member, Mapping) else {}
    return _first_text(
        feature.get("product"),
        member.get("product"),
        feature.get("note"),
        feature.get("label"),
        feature.get("display_label"),
        feature.get("displayLabel"),
    )


def _feature_to_member(feature: Mapping[str, object] | None) -> dict[str, object] | None:
    if feature is None:
        return None
    member = feature.get("orthogroup_member") or feature.get("orthogroupMember")
    member = member if isinstance(member, Mapping) else {}
    return {
        "recordId": _first_text(member.get("recordId"), member.get("record_id"), feature.get("record_id")),
        "start": member.get("start") if member.get("start") is not None else feature.get("start"),
        "end": member.get("end") if member.get("end") is not None else feature.get("end"),
        "strand": _first_text(member.get("strand"), feature.get("strand")),
        "proteinId": _first_text(
            member.get("proteinId"),
            member.get("protein_id"),
            feature.get("protein_id"),
            feature.get("proteinId"),
        ),
        "sourceProteinId": _first_text(
            member.get("sourceProteinId"),
            member.get("source_protein_id"),
            feature.get("source_protein_id"),
            feature.get("sourceProteinId"),
        ),
        "product": _first_text(member.get("product"), _feature_product(feature)),
        "note": _first_text(member.get("note"), feature.get("note")),
        "featureSvgId": _first_text(member.get("featureSvgId"), member.get("feature_svg_id"), feature.get("svg_id")),
    }


def _build_fallback_orthogroup(
    *,
    orthogroup_id: str,
    query_feature: Mapping[str, object] | None,
    subject_feature: Mapping[str, object] | None,
    features: Sequence[Mapping[str, object]],
) -> dict[str, object] | None:
    group_id = _first_text(orthogroup_id)
    if not group_id:
        return None
    matching_features = [
        feature
        for feature in features
        if _feature_orthogroup_id(feature) == group_id
    ]
    fallback_features = matching_features or [
        feature for feature in (query_feature, subject_feature) if feature is not None
    ]
    if not fallback_features:
        return None
    members = [
        member
        for member in (_feature_to_member(feature) for feature in fallback_features)
        if member is not None
    ]
    first_feature = fallback_features[0]
    record_coverage_fallback = len(
        {
            _first_text(
                feature.get("record_id"),
                feature.get("recordId"),
                feature.get("record_idx"),
                feature.get("recordIndex"),
            )
            for feature in fallback_features
            if _first_text(
                feature.get("record_id"),
                feature.get("recordId"),
                feature.get("record_idx"),
                feature.get("recordIndex"),
            )
        }
    )
    return {
        "id": group_id,
        "name": _feature_product(first_feature),
        "member_count": _first_text(
            first_feature.get("orthogroup_member_count"),
            first_feature.get("orthogroupMemberCount"),
            len(members),
        ),
        "record_coverage_count": _first_text(
            first_feature.get("orthogroup_record_coverage"),
            first_feature.get("orthogroupRecordCoverage"),
            record_coverage_fallback,
        ),
        "members": members,
    }


def _group_has_feature_svg_id(group: Mapping[str, object], feature_svg_id: str) -> bool:
    ids = set(_split_metadata_values(feature_svg_id))
    if not ids:
        return False
    members = group.get("members")
    if not isinstance(members, Sequence) or isinstance(members, (str, bytes)):
        return False
    return any(
        isinstance(member, Mapping) and _member_feature_svg_id(member) in ids
        for member in members
    )


def _find_orthogroup(
    orthogroups: Sequence[Mapping[str, object]],
    by_id: Mapping[str, Mapping[str, object]],
    orthogroup_id: str,
    query_feature_svg_id: str,
    subject_feature_svg_id: str,
) -> Mapping[str, object] | None:
    if orthogroup_id and orthogroup_id in by_id:
        return by_id[orthogroup_id]
    for group in orthogroups:
        if _group_has_feature_svg_id(group, query_feature_svg_id) or _group_has_feature_svg_id(
            group,
            subject_feature_svg_id,
        ):
            return group
    return None


def _member_location_text(member: Mapping[str, object]) -> str:
    try:
        start = int(member.get("start")) + 1
    except (TypeError, ValueError):
        start = _first_text(member.get("start"))
    try:
        end = int(member.get("end"))
    except (TypeError, ValueError):
        end = _first_text(member.get("end"))
    if start and end:
        text = f"{start}..{end}"
    else:
        text = _first_text(start, end)
    strand = _strand_symbol(member.get("strand"))
    return f"{text} ({strand})" if text and strand else text


def _display_protein_id(
    feature: Mapping[str, object] | None,
    member: Mapping[str, object] | None,
    fallback: object = "",
) -> str:
    qualifiers = feature.get("qualifiers") if feature else {}
    protein_qualifier = ""
    locus_tag = ""
    gene_id = ""
    old_locus_tag = ""
    id_qualifier = ""
    name_qualifier = ""
    parent_qualifier = ""
    gene_qualifier = ""
    if isinstance(qualifiers, Mapping):
        protein_qualifier = _first_text(qualifiers.get("protein_id"))
        locus_tag = _first_text(qualifiers.get("locus_tag"))
        gene_id = _first_text(qualifiers.get("gene_id"))
        old_locus_tag = _first_text(qualifiers.get("old_locus_tag"))
        id_qualifier = _first_text(qualifiers.get("ID"))
        name_qualifier = _first_text(qualifiers.get("Name"))
        parent_qualifier = _first_text(qualifiers.get("Parent"))
        gene_qualifier = _first_text(qualifiers.get("gene"))
    return _first_text(
        feature.get("source_protein_id") if feature else "",
        feature.get("sourceProteinId") if feature else "",
        member.get("sourceProteinId") if member else "",
        member.get("source_protein_id") if member else "",
        protein_qualifier,
        feature.get("locus_tag") if feature else "",
        feature.get("locusTag") if feature else "",
        locus_tag,
        member.get("locusTag") if member else "",
        member.get("locus_tag") if member else "",
        feature.get("gene_id") if feature else "",
        feature.get("geneId") if feature else "",
        gene_id,
        member.get("geneId") if member else "",
        member.get("gene_id") if member else "",
        feature.get("old_locus_tag") if feature else "",
        feature.get("oldLocusTag") if feature else "",
        old_locus_tag,
        member.get("oldLocusTag") if member else "",
        member.get("old_locus_tag") if member else "",
        feature.get("ID") if feature else "",
        id_qualifier,
        feature.get("Name") if feature else "",
        name_qualifier,
        feature.get("Parent") if feature else "",
        parent_qualifier,
        feature.get("gene") if feature else "",
        gene_qualifier,
        member.get("gene") if member else "",
        feature.get("protein_id") if feature else "",
        feature.get("proteinId") if feature else "",
        member.get("proteinId") if member else "",
        member.get("protein_id") if member else "",
        fallback,
    )


def _orthogroup_display_name(group: Mapping[str, object] | None) -> str:
    if group is None:
        return ""
    return _first_text(group.get("display_name"), group.get("displayName"), group.get("name"))


def _build_match_member_rows(
    group: Mapping[str, object] | None,
) -> list[dict[str, object]]:
    if group is None:
        return []
    members = group.get("members")
    if not isinstance(members, Sequence) or isinstance(members, (str, bytes)):
        return []
    group_id = _first_text(group.get("id"), group.get("orthogroupId"), group.get("orthogroup_id"))
    display_name = _orthogroup_display_name(group)
    rows: list[dict[str, object]] = []
    for member in members:
        if not isinstance(member, Mapping):
            continue
        feature_svg_id = _member_feature_svg_id(member)
        row = {
            "featureSvgId": feature_svg_id,
            "feature_svg_id": feature_svg_id,
            "orthogroupId": group_id,
            "orthogroup_id": group_id,
            "displayName": display_name,
            "display_name": display_name,
            "record": _first_text(member.get("recordId"), member.get("record_id")),
            "coordinates": _member_location_text(member),
            "proteinId": _display_protein_id(None, member),
            "productOrNote": _first_text(member.get("product"), member.get("note")),
        }
        if (
            row["record"]
            or row["coordinates"]
            or row["proteinId"]
            or row["productOrNote"]
            or row["featureSvgId"]
        ):
            rows.append(row)
    return rows


def _get_group_member_for_feature_svg_id(
    group: Mapping[str, object] | None,
    feature_svg_id: str,
) -> Mapping[str, object] | None:
    if group is None:
        return None
    members = group.get("members")
    if not isinstance(members, Sequence) or isinstance(members, (str, bytes)):
        return None
    for member in members:
        if isinstance(member, Mapping) and _member_feature_svg_id(member) == feature_svg_id:
            return member
    return None


def _resolve_block_member_labels(
    group: Mapping[str, object] | None,
    feature_svg_ids: str,
    features_by_id: Mapping[str, Mapping[str, object]],
) -> str:
    if group is None:
        return ""
    values: list[str] = []
    for feature_svg_id in _unique_metadata_values(feature_svg_ids):
        feature = features_by_id.get(feature_svg_id)
        member = _get_group_member_for_feature_svg_id(group, feature_svg_id)
        if member is None:
            continue
        text = _display_protein_id(feature, member)
        _add_unique_display_text(values, text)
    return "; ".join(values)


def _orthogroup_title(orthogroup_id: str, display_name: str) -> str:
    group_id = _first_text(orthogroup_id)
    name = _first_text(display_name)
    if group_id and name:
        return f"{group_id}:{name}"
    return group_id or name or _MATCH_TITLES["orthogroup"]


def _feature_location_text(feature: Mapping[str, object] | None) -> str:
    if feature is None:
        return ""
    direct = _first_text(feature.get("location"))
    if direct and direct != "..":
        return direct
    built = _build_feature_location(feature)
    return built if built and built != ".." else ""


def _feature_row_locus_id(feature: Mapping[str, object] | None, fallback_locus_id: object) -> str:
    if feature is None:
        return _first_text(fallback_locus_id)
    return _first_text(
        feature.get("locusTag"),
        feature.get("locus_tag"),
        _first_qualifier_value(feature, "locus_tag"),
        feature.get("geneId"),
        feature.get("gene_id"),
        _first_qualifier_value(feature, "gene_id"),
        fallback_locus_id,
    )


def _feature_row_display_name(
    feature: Mapping[str, object] | None,
    fallback_display_name: object,
) -> str:
    if feature is None:
        return _first_text(fallback_display_name)
    return _first_text(
        fallback_display_name,
        feature.get("displayLabel"),
        feature.get("display_label"),
        feature.get("label"),
        feature.get("gene"),
        _first_qualifier_value(feature, "gene"),
        feature.get("locus_tag"),
        feature.get("locusTag"),
        _first_qualifier_value(feature, "locus_tag"),
        _feature_product(feature),
    )


def _list_get(values: Sequence[str], index: int) -> str:
    return values[index] if index < len(values) else ""


def _resolve_feature_section_protein_ids(
    *,
    feature: Mapping[str, object] | None,
    feature_svg_ids: str,
    features_by_id: Mapping[str, Mapping[str, object]],
    group: Mapping[str, object] | None,
    fallback_protein_ids: str,
    locus_id: str,
    display_name: str,
) -> str:
    values: list[str] = []

    def add_feature_protein_id(
        candidate_feature: Mapping[str, object] | None,
        member: Mapping[str, object] | None = None,
    ) -> None:
        _add_unique_display_text(values, _display_protein_id(candidate_feature, member, ""))

    ids = _split_metadata_values(feature_svg_ids)
    for feature_svg_id in ids:
        candidate_feature = features_by_id.get(feature_svg_id)
        member = _get_group_member_for_feature_svg_id(group, feature_svg_id)
        add_feature_protein_id(candidate_feature, member)
    if feature is not None:
        add_feature_protein_id(
            feature,
            _get_group_member_for_feature_svg_id(group, ids[0])
            if len(ids) == 1
            else None,
        )
    if not values:
        for value in _split_metadata_values(locus_id):
            _add_unique_display_text(values, value)
    if not values:
        for value in _split_metadata_values(display_name):
            _add_unique_display_text(values, value)
    if not values:
        for value in _split_metadata_values(fallback_protein_ids):
            _add_unique_display_text(values, value)
    return "; ".join(values)


def _build_feature_list_rows(
    *,
    feature_svg_ids: str,
    features_by_id: Mapping[str, Mapping[str, object]],
    group: Mapping[str, object] | None,
    record_id: str,
    interval: str,
    protein_id: str,
    locus_id: str,
    display_name: str,
) -> list[dict[str, object]]:
    feature_ids = _unique_metadata_values(feature_svg_ids)
    protein_ids = _split_metadata_values(protein_id)
    locus_ids = _split_metadata_values(locus_id)
    display_names = _split_metadata_values(display_name)
    count = max(len(feature_ids), len(protein_ids), len(locus_ids), len(display_names))
    if count == 0:
        return []

    rows: list[dict[str, object]] = []
    for index in range(count):
        svg_id = _list_get(feature_ids, index)
        feature = features_by_id.get(svg_id) if svg_id else None
        member = _get_group_member_for_feature_svg_id(group, svg_id) if svg_id else None
        fallback_protein_id = _first_non_internal_display_text(
            _list_get(locus_ids, index),
            _list_get(display_names, index),
            _list_get(protein_ids, index),
        )
        resolved_protein_id = _display_protein_id(feature, member, "")
        display_protein_id = _first_text(
            "" if _is_internal_display_id(resolved_protein_id) else resolved_protein_id,
            fallback_protein_id,
            resolved_protein_id,
            _list_get(protein_ids, index),
        )
        row_record = _first_text(
            feature.get("record_id") if feature else "",
            feature.get("recordId") if feature else "",
            member.get("recordId") if member else "",
            member.get("record_id") if member else "",
            record_id,
        )
        row_location = _first_text(
            _feature_location_text(feature),
            _member_location_text(member) if member else "",
            interval if count == 1 else "",
        )
        row_locus_id = _feature_row_locus_id(feature, _list_get(locus_ids, index))
        row_display_name = _feature_row_display_name(feature, _list_get(display_names, index))
        product = _feature_product(feature)
        label = _first_text(
            display_protein_id,
            row_display_name,
            row_locus_id,
            product,
            svg_id,
            f"Feature {index + 1}",
        )
        copy_text = "\t".join(
            [
                row_record,
                row_location,
                display_protein_id,
                row_locus_id,
                row_display_name,
                product,
            ]
        )
        row = {
            "key": f"{svg_id or 'feature'}-{index}",
            "svg_id": svg_id,
            "svgId": svg_id,
            "can_open": bool(feature and feature.get("svg_id")),
            "canOpen": bool(feature and feature.get("svg_id")),
            "label": label,
            "record": row_record,
            "location": row_location,
            "protein_id": display_protein_id,
            "proteinId": display_protein_id,
            "locus_id": row_locus_id,
            "locusId": row_locus_id,
            "display_name": row_display_name,
            "displayName": row_display_name,
            "product": product,
            "type": _first_text(feature.get("type") if feature else ""),
            "copy_text": copy_text,
            "copyText": copy_text,
        }
        if (
            row["svg_id"]
            or row["record"]
            or row["location"]
            or row["protein_id"]
            or row["locus_id"]
            or row["display_name"]
            or row["product"]
        ):
            rows.append(row)
    return rows


def _build_match_feature_section(
    *,
    title: str,
    feature: Mapping[str, object] | None,
    record_id: str,
    interval: str,
    protein_id: str,
    locus_id: str,
    display_name: str,
    feature_svg_ids: str,
    features_by_id: Mapping[str, Mapping[str, object]],
    group: Mapping[str, object] | None,
) -> dict[str, object]:
    rows: list[list[str]] = []
    feature_rows = _build_feature_list_rows(
        feature_svg_ids=feature_svg_ids,
        features_by_id=features_by_id,
        group=group,
        record_id=record_id,
        interval=interval,
        protein_id=protein_id,
        locus_id=locus_id,
        display_name=display_name,
    )
    display_protein_ids = _resolve_feature_section_protein_ids(
        feature=feature,
        feature_svg_ids=feature_svg_ids,
        features_by_id=features_by_id,
        group=group,
        fallback_protein_ids=protein_id,
        locus_id=locus_id,
        display_name=display_name,
    )
    _add_match_row(rows, "Record", _first_text(feature.get("record_id") if feature else "", record_id))
    _add_match_row(rows, "Location", _first_text(feature.get("location") if feature else "", interval))
    _add_match_row(rows, "Protein IDs" if ";" in display_protein_ids else "Protein ID", display_protein_ids)
    _add_match_row(rows, "Type", feature.get("type") if feature else "")
    _add_match_row(rows, "Locus ID", locus_id)
    _add_match_row(rows, "Display name", display_name)
    return _match_section(
        title,
        rows,
        feature_rows=feature_rows,
        featureRows=feature_rows,
    )


def _build_orthogroup_detail_rows(
    *,
    orthogroup_id: str,
    display_name: str,
    description: str,
    member_count: str,
    record_coverage: str,
) -> list[list[str]]:
    rows: list[list[str]] = []
    _add_match_row(rows, "Orthogroup ID", orthogroup_id)
    _add_match_row(rows, "Display name", display_name)
    _add_match_row(rows, "Description", description)
    _add_match_row(rows, "Members", member_count)
    _add_match_row(rows, "Record coverage", record_coverage)
    return rows


def _build_block_orthogroups(
    *,
    orthogroup_ids: Sequence[str],
    orthogroups_by_id: Mapping[str, Mapping[str, object]],
    features_by_id: Mapping[str, Mapping[str, object]],
    query_feature_svg_id: str,
    subject_feature_svg_id: str,
) -> list[dict[str, object]]:
    block_orthogroups: list[dict[str, object]] = []
    for orthogroup_id in orthogroup_ids:
        group = orthogroups_by_id.get(orthogroup_id)
        display_name = _orthogroup_display_name(group)
        description = _first_text(group.get("description")) if group else ""
        member_count = (
            _first_text(
                group.get("member_count"),
                group.get("memberCount"),
                len(group.get("members")) if isinstance(group.get("members"), Sequence) else "",
            )
            if group
            else ""
        )
        record_coverage = (
            _first_text(group.get("record_coverage_count"), group.get("recordCoverage"))
            if group
            else ""
        )
        member_rows = _build_match_member_rows(group)
        block_orthogroups.append(
            {
                "id": orthogroup_id,
                "display_name": display_name,
                "displayName": display_name,
                "description": description,
                "member_count": member_count,
                "memberCount": member_count,
                "record_coverage": record_coverage,
                "recordCoverage": record_coverage,
                "query_member": _resolve_block_member_labels(
                    group,
                    query_feature_svg_id,
                    features_by_id,
                ),
                "subject_member": _resolve_block_member_labels(
                    group,
                    subject_feature_svg_id,
                    features_by_id,
                ),
                "detail_rows": _build_orthogroup_detail_rows(
                    orthogroup_id=orthogroup_id,
                    display_name=display_name,
                    description=description,
                    member_count=member_count,
                    record_coverage=record_coverage,
                ),
                "member_rows": member_rows,
                "member_copy_text": _member_copy_text(member_rows),
            }
        )
    return block_orthogroups


def _member_copy_text(member_rows: Sequence[Mapping[str, object]]) -> str:
    if not member_rows:
        return ""
    return "\n".join(
        [
            "Record\tCoordinates (+/-)\tProtein ID\tProduct / note",
            *[
                "\t".join(
                    [
                        _first_text(row.get("record")),
                        _first_text(row.get("coordinates")),
                        _first_text(row.get("proteinId")),
                        _first_text(row.get("productOrNote")),
                    ]
                )
                for row in member_rows
            ],
        ]
    )


def _match_section(title: str, rows: Sequence[Sequence[str]], **extras: object) -> dict[str, object]:
    section: dict[str, object] = {
        "title": title,
        "rows": [
            [str(row[0]), str(row[1])]
            for row in rows
            if len(row) >= 2 and _first_text(row[1])
        ],
    }
    section.update(extras)
    return section


def _match_payload_v1(
    element: ET.Element,
    index: int,
    *,
    features: Sequence[Mapping[str, object]],
    features_by_id: Mapping[str, Mapping[str, object]],
    orthogroups: Sequence[Mapping[str, object]],
    orthogroups_by_id: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
    match_id = _first_text(element.get("data-gbdraw-pairwise-match-id"))
    if not match_id:
        match_id = f"pairwise_match_{index + 1}"
        element.set("data-gbdraw-pairwise-match-id", match_id)

    kind = _match_kind(element)
    orthogroup_id = _first_text(element.get("data-orthogroup-id"))
    orthogroup_ids = _unique_metadata_values(orthogroup_id)
    block_id = _first_text(element.get("data-collinearity-block-id"))
    query_feature_svg_id = _first_text(element.get("data-query-feature-svg-id"))
    subject_feature_svg_id = _first_text(element.get("data-subject-feature-svg-id"))
    query_feature = features_by_id.get(query_feature_svg_id)
    subject_feature = features_by_id.get(subject_feature_svg_id)
    group = (
        None
        if kind == "collinear"
        else _find_orthogroup(
            orthogroups,
            orthogroups_by_id,
            orthogroup_id,
            query_feature_svg_id,
            subject_feature_svg_id,
        )
        or _build_fallback_orthogroup(
            orthogroup_id=orthogroup_id,
            query_feature=query_feature,
            subject_feature=subject_feature,
            features=features,
        )
    )
    block_orthogroups: list[dict[str, object]] = []
    if kind == "collinear":
        block_orthogroups = _build_block_orthogroups(
            orthogroup_ids=orthogroup_ids,
            orthogroups_by_id=orthogroups_by_id,
            features_by_id=features_by_id,
            query_feature_svg_id=query_feature_svg_id,
            subject_feature_svg_id=subject_feature_svg_id,
        )

    query_interval = _interval_text(element.get("data-qstart"), element.get("data-qend"))
    subject_interval = _interval_text(element.get("data-sstart"), element.get("data-send"))
    summary_rows: list[list[str]] = []
    _add_match_row(
        summary_rows,
        "Query record",
        _first_text(element.get("data-query-record-id"), element.get("data-query")),
    )
    _add_match_row(
        summary_rows,
        "Subject record",
        _first_text(element.get("data-subject-record-id"), element.get("data-subject")),
    )
    _add_match_row(summary_rows, "Query interval", query_interval)
    _add_match_row(summary_rows, "Subject interval", subject_interval)
    _add_match_row(
        summary_rows,
        "Orientation",
        _first_text(element.get("data-collinearity-orientation"), element.get("data-orientation")),
    )

    alignment_rows: list[list[str]] = []
    _add_match_row(alignment_rows, "Identity", element.get("data-identity"))
    _add_match_row(alignment_rows, "Alignment length", element.get("data-alignment-length"))
    _add_match_row(alignment_rows, "E-value", element.get("data-evalue"))
    _add_match_row(alignment_rows, "Bit score", element.get("data-bitscore"))
    _add_match_row(alignment_rows, "Mismatches", element.get("data-mismatches"))
    _add_match_row(alignment_rows, "Gap opens", element.get("data-gap-opens"))

    orthogroup_display_name = _first_text(
        group.get("display_name") if group else "",
        group.get("displayName") if group else "",
        group.get("name") if group else "",
        _feature_product(query_feature),
        _feature_product(subject_feature),
    )
    orthogroup_rows: list[list[str]] = []
    if kind != "collinear":
        _add_match_row(orthogroup_rows, "Orthogroup ID", orthogroup_id)
        _add_match_row(orthogroup_rows, "Display name", orthogroup_display_name)
        _add_match_row(orthogroup_rows, "Description", group.get("description") if group else "")
        _add_match_row(
            orthogroup_rows,
            "Members",
            _first_text(group.get("member_count"), group.get("memberCount")) if group else "",
        )
        _add_match_row(
            orthogroup_rows,
            "Record coverage",
            _first_text(group.get("record_coverage_count"), group.get("recordCoverage")) if group else "",
        )
    member_rows = _build_match_member_rows(group)

    block_rows: list[list[str]] = []
    _add_match_row(block_rows, "Block ID", block_id)
    _add_match_row(block_rows, "Kind", element.get("data-collinearity-block-kind"))
    _add_match_row(block_rows, "Orientation", element.get("data-collinearity-orientation"))
    _add_match_row(block_rows, "Color mode", element.get("data-collinearity-color-mode"))
    if kind == "collinear":
        _add_match_row(block_rows, "Average identity", element.get("data-identity"))
        _add_match_row(block_rows, "Aligned length", element.get("data-alignment-length"))
    _add_match_row(block_rows, "Block score", element.get("data-collinearity-block-score"))
    _add_match_row(block_rows, "Block e-value", element.get("data-collinearity-block-evalue"))
    _add_match_row(
        block_rows,
        "Anchor",
        " / ".join(
            value
            for value in (
                _first_text(element.get("data-collinearity-anchor-index")),
                _first_text(element.get("data-collinearity-anchor-count")),
            )
            if value
        ),
    )

    if kind == "orthogroup":
        summary_section = _match_section(
            "Summary",
            orthogroup_rows,
            member_rows=member_rows,
            member_copy_text=_member_copy_text(member_rows),
        )
        hover_rows: list[list[str]] = []
        _add_match_row(hover_rows, "Kind", kind)
        _add_match_row(hover_rows, "Orthogroup", orthogroup_id)
        _add_match_row(hover_rows, "Display name", orthogroup_display_name)
        _add_match_row(
            hover_rows,
            "Members",
            _first_text(group.get("member_count"), group.get("memberCount")) if group else "",
        )
        return {
            "id": match_id,
            "title": _orthogroup_title(orthogroup_id, orthogroup_display_name),
            "subtitle": "",
            "match_kind": kind,
            "orthogroup_id": orthogroup_id,
            "collinearity_block_id": block_id,
            "fill": _first_text(element.get("fill"), "#94a3b8"),
            "sections": (
                [summary_section]
                if summary_section["rows"] or summary_section.get("member_rows")
                else []
            ),
            "hover_rows": hover_rows,
        }

    sections: list[dict[str, object]] = [_match_section("Summary", summary_rows)]
    if kind != "collinear":
        sections.append(_match_section("Alignment", alignment_rows))
    if orthogroup_rows or kind == "orthogroup":
        sections.append(
            _match_section(
                "Orthogroup",
                orthogroup_rows,
                member_rows=member_rows,
                member_copy_text=_member_copy_text(member_rows),
            )
        )
    if kind == "collinear":
        block_orthogroup_rows: list[list[str]] = []
        _add_match_row(
            block_orthogroup_rows,
            "Number of orthogroups covered",
            str(len(orthogroup_ids)),
        )
        sections.append(
            _match_section(
                "Orthogroups covered",
                block_orthogroup_rows,
                block_orthogroups=block_orthogroups,
            )
        )
    if block_rows or kind == "collinear":
        sections.append(_match_section("Collinearity", block_rows))
    sections.append(
        _build_match_feature_section(
            title="Query",
            feature=query_feature,
            record_id=_first_text(element.get("data-query-record-id")),
            interval=query_interval,
            protein_id=_first_text(element.get("data-query-protein-id")),
            locus_id=_first_text(element.get("data-query-locus-id")),
            display_name=_first_text(element.get("data-query-display-name")),
            feature_svg_ids=query_feature_svg_id,
            features_by_id=features_by_id,
            group=group,
        )
    )
    sections.append(
        _build_match_feature_section(
            title="Subject",
            feature=subject_feature,
            record_id=_first_text(element.get("data-subject-record-id")),
            interval=subject_interval,
            protein_id=_first_text(element.get("data-subject-protein-id")),
            locus_id=_first_text(element.get("data-subject-locus-id")),
            display_name=_first_text(element.get("data-subject-display-name")),
            feature_svg_ids=subject_feature_svg_id,
            features_by_id=features_by_id,
            group=group,
        )
    )

    def find_row(section_title: str, row_label: str) -> str:
        for section in sections:
            if section.get("title") != section_title:
                continue
            rows = section.get("rows")
            if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)):
                continue
            for row in rows:
                if isinstance(row, Sequence) and not isinstance(row, (str, bytes)) and len(row) >= 2 and row[0] == row_label:
                    return _first_text(row[1])
        return ""

    hover_rows: list[list[str]] = []
    _add_match_row(hover_rows, "Kind", kind)
    _add_match_row(
        hover_rows,
        "Identity",
        _first_text(
            find_row("Alignment", "Identity"),
            find_row("Collinearity", "Average identity"),
        ),
    )
    _add_match_row(hover_rows, "Query", find_row("Summary", "Query interval"))
    _add_match_row(hover_rows, "Subject", find_row("Summary", "Subject interval"))
    if kind == "collinear":
        _add_match_row(hover_rows, "Orthogroups", str(len(block_orthogroups) or len(orthogroup_ids)))
    else:
        _add_match_row(hover_rows, "Orthogroup", orthogroup_id)
    _add_match_row(hover_rows, "Block", block_id)

    return {
        "id": match_id,
        "title": _MATCH_TITLES.get(kind, _MATCH_TITLES["pairwise"]),
        "subtitle": _first_text(
            block_id,
            orthogroup_id,
            match_id,
        ),
        "match_kind": kind,
        "orthogroup_id": orthogroup_id,
        "collinearity_block_id": block_id,
        "block_orthogroup_count": len(block_orthogroups)
        or (len(orthogroup_ids) if kind == "collinear" else 0),
        "block_orthogroups": block_orthogroups,
        "fill": _first_text(element.get("fill"), "#94a3b8"),
        "sections": [
            section
            for section in sections
            if section["rows"]
            or (
                isinstance(section.get("feature_rows"), Sequence)
                and not isinstance(section.get("feature_rows"), (str, bytes))
                and len(section.get("feature_rows")) > 0
            )
        ],
        "hover_rows": hover_rows,
    }


def _match_payload_v2(element: ET.Element, index: int) -> dict[str, object]:
    match_id = _first_text(
        element.get("data-gbdraw-match-id"),
        element.get("data-gbdraw-pairwise-match-id"),
    )
    if not match_id:
        match_id = f"match_{index + 1}"
    element.set("data-gbdraw-match-id", match_id)
    payload = {
        "id": match_id,
        "match_kind": _match_kind(element),
        "orthogroup_ids": _unique_metadata_values(element.get("data-orthogroup-id")),
        "collinearity_block_id": _first_text(element.get("data-collinearity-block-id")),
        "fill": _first_text(element.get("fill"), "#94a3b8"),
        "query_record_id": _first_text(
            element.get("data-query-record-id"), element.get("data-query")
        ),
        "query_record_index": _first_text(element.get("data-query-record-index")),
        "subject_record_id": _first_text(
            element.get("data-subject-record-id"), element.get("data-subject")
        ),
        "subject_record_index": _first_text(element.get("data-subject-record-index")),
        "qstart": _first_text(element.get("data-qstart")),
        "qend": _first_text(element.get("data-qend")),
        "sstart": _first_text(element.get("data-sstart")),
        "send": _first_text(element.get("data-send")),
        "orientation": _first_text(
            element.get("data-collinearity-orientation"), element.get("data-orientation")
        ),
        "identity": _first_text(element.get("data-identity")),
        "alignment_length": _first_text(element.get("data-alignment-length")),
        "evalue": _first_text(element.get("data-evalue")),
        "bitscore": _first_text(element.get("data-bitscore")),
        "mismatches": _first_text(element.get("data-mismatches")),
        "gap_opens": _first_text(element.get("data-gap-opens")),
        "source_index": _first_text(element.get("data-source-index")),
        "track_index": _first_text(element.get("data-track-index")),
        "track_label": _first_text(element.get("data-track-label")),
        "reference_side": _first_text(element.get("data-reference-side")),
        "reference_record_id": _first_text(element.get("data-reference-record-id")),
        "block_kind": _first_text(element.get("data-collinearity-block-kind")),
        "collinear_group_scope": _first_text(element.get("data-collinear-group-scope")),
        "group_kind": _first_text(element.get("data-group-kind")),
        "block_color_mode": _first_text(element.get("data-collinearity-color-mode")),
        "block_score": _first_text(element.get("data-collinearity-block-score")),
        "block_evalue": _first_text(element.get("data-collinearity-block-evalue")),
        "anchor_index": _first_text(element.get("data-collinearity-anchor-index")),
        "anchor_count": _first_text(element.get("data-collinearity-anchor-count")),
        "query_feature_svg_id": _first_text(element.get("data-query-feature-svg-id")),
        "subject_feature_svg_id": _first_text(element.get("data-subject-feature-svg-id")),
        "query_protein_id": _first_text(element.get("data-query-protein-id")),
        "subject_protein_id": _first_text(element.get("data-subject-protein-id")),
        "query_locus_id": _first_text(element.get("data-query-locus-id")),
        "subject_locus_id": _first_text(element.get("data-subject-locus-id")),
        "query_display_name": _first_text(element.get("data-query-display-name")),
        "subject_display_name": _first_text(element.get("data-subject-display-name")),
    }
    return dict(_compact_wire_value(payload) or {})


def _match_payloads(
    root: ET.Element,
    features: Sequence[Mapping[str, object]],
    orthogroups: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    payloads: list[dict[str, object]] = []
    for element in root.iter():
        if not _is_match_candidate(element):
            continue
        payloads.append(
            _match_payload_v2(element, len(payloads))
        )
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
        kind = _first_text(match.get("match_kind"))
        if kind == "homology":
            reference_side = _first_text(match.get("reference_side"))
            circular_record_ids.add(
                _first_text(match.get(f"{reference_side}_record_id"))
            )
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
        origin = _first_text(source.get("origin"))
        if origin == "linear-record":
            try:
                if int(str(source.get("recordIndex"))) in linear_indexes:
                    selected.append(source)
            except (TypeError, ValueError):
                continue
        elif origin == "circular-reference":
            aliases = {
                _first_text(source.get("recordId")),
                *(_normalize_string_array(source.get("aliases"))),
            }
            if aliases & circular_record_ids:
                selected.append(source)
        elif origin == "homology-comparison":
            try:
                if int(str(source.get("sourceIndex"))) in comparison_source_indexes:
                    selected.append(source)
            except (TypeError, ValueError):
                continue
    return selected


def _extract_template_literal(source: str, name: str) -> str:
    pattern = re.compile(rf"export\s+const\s+{re.escape(name)}\s*=\s*`")
    match = pattern.search(source)
    if match is None:
        raise GbdrawError(f"Missing standalone interactive SVG asset: {name}")
    i = match.end()
    chars: list[str] = []
    while i < len(source):
        char = source[i]
        if char == "`":
            return "".join(chars)
        if char == "\\":
            if i + 1 >= len(source):
                raise GbdrawError(f"Malformed standalone interactive SVG asset: {name}")
            next_char = source[i + 1]
            if next_char in {"`", "\\"}:
                chars.append(next_char)
            elif next_char == "n":
                chars.append("\n")
            elif next_char == "r":
                chars.append("\r")
            elif next_char == "t":
                chars.append("\t")
            elif next_char == "$" and i + 2 < len(source) and source[i + 2] == "{":
                chars.append("${")
                i += 1
            else:
                chars.append("\\")
                chars.append(next_char)
            i += 2
            continue
        chars.append(char)
        i += 1
    raise GbdrawError(f"Malformed standalone interactive SVG asset: {name}")


def _load_standalone_assets() -> tuple[str, str]:
    try:
        source = (
            resources.files("gbdraw.web")
            .joinpath("js/services/standalone-interactivity-assets.js")
            .read_text(encoding="utf-8")
        )
    except Exception as exc:
        raise GbdrawError("Missing standalone interactive SVG runtime/style package asset.") from exc
    return (
        _extract_template_literal(source, "STANDALONE_INTERACTIVE_STYLE"),
        _extract_template_literal(source, "STANDALONE_INTERACTIVE_SCRIPT"),
    )


def _remove_existing_assets(root: ET.Element) -> None:
    parent_by_child = {child: parent for parent in root.iter() for child in list(parent)}
    for element in list(root.iter()):
        if element.get("id") not in _ASSET_IDS:
            continue
        parent = parent_by_child.get(element)
        if parent is not None:
            parent.remove(element)


def _ensure_defs(root: ET.Element) -> ET.Element:
    for child in list(root):
        if _local_name(child.tag) == "defs":
            return child
    defs = ET.Element(_svg_tag("defs"))
    root.insert(0, defs)
    return defs


def _append_glow_filter(
    defs: ET.Element,
    *,
    filter_id: str,
    color: str,
    opacity: str,
    blur_std_deviation: str,
    slope: str,
    extent: str,
) -> None:
    filter_element = ET.SubElement(
        defs,
        _svg_tag("filter"),
        {
            "id": filter_id,
            "x": f"-{extent}%",
            "y": f"-{extent}%",
            "width": f"{100 + (float(extent) * 2):g}%",
            "height": f"{100 + (float(extent) * 2):g}%",
            "color-interpolation-filters": "sRGB",
        },
    )
    component_transfer = ET.SubElement(
        filter_element,
        _svg_tag("feComponentTransfer"),
        {"in": "SourceGraphic", "result": "gbdrawBrightenedFeature"},
    )
    for channel in ("R", "G", "B"):
        ET.SubElement(
            component_transfer,
            _svg_tag(f"feFunc{channel}"),
            {"type": "linear", "slope": slope},
        )
    ET.SubElement(
        filter_element,
        _svg_tag("feGaussianBlur"),
        {
            "in": "SourceAlpha",
            "stdDeviation": blur_std_deviation,
            "result": "gbdrawFeatureGlowBlur",
        },
    )
    ET.SubElement(
        filter_element,
        _svg_tag("feFlood"),
        {
            "flood-color": color,
            "flood-opacity": opacity,
            "result": "gbdrawFeatureGlowColor",
        },
    )
    ET.SubElement(
        filter_element,
        _svg_tag("feComposite"),
        {
            "in": "gbdrawFeatureGlowColor",
            "in2": "gbdrawFeatureGlowBlur",
            "operator": "in",
            "result": "gbdrawFeatureGlow",
        },
    )
    merge = ET.SubElement(filter_element, _svg_tag("feMerge"))
    ET.SubElement(merge, _svg_tag("feMergeNode"), {"in": "gbdrawFeatureGlow"})
    ET.SubElement(merge, _svg_tag("feMergeNode"), {"in": "gbdrawBrightenedFeature"})


def _ensure_glow_filters(root: ET.Element) -> None:
    defs = _ensure_defs(root)
    _append_glow_filter(
        defs,
        filter_id=INTERACTIVE_GLOW_FILTER_ID,
        color="#2563eb",
        opacity="0.85",
        blur_std_deviation="3",
        slope="1.2",
        extent="35",
    )
    _append_glow_filter(
        defs,
        filter_id=INTERACTIVE_MATCH_GLOW_FILTER_ID,
        color="#fbbf24",
        opacity="0.32",
        blur_std_deviation="1.5",
        slope="1.04",
        extent="25",
    )


def _parse_viewbox(value: str | None) -> str | None:
    if not value:
        return None
    parts = [part for part in re.split(r"[\s,]+", value.strip()) if part]
    if len(parts) < 4:
        return None
    try:
        numbers = [float(part) for part in parts[:4]]
    except ValueError:
        return None
    if numbers[2] <= 0 or numbers[3] <= 0:
        return None
    return " ".join(f"{value:g}" for value in numbers)


def _float_attr(value: str | None) -> float | None:
    if not value:
        return None
    match = re.match(r"\s*([0-9]+(?:\.[0-9]+)?)", value)
    if match is None:
        return None
    parsed = float(match.group(1))
    return parsed if parsed > 0 else None


def _original_viewbox(root: ET.Element) -> str:
    viewbox = _parse_viewbox(root.get("data-gbdraw-original-viewbox")) or _parse_viewbox(
        root.get("viewBox")
    )
    if viewbox:
        return viewbox
    width = _float_attr(root.get("width"))
    height = _float_attr(root.get("height"))
    if width is not None and height is not None:
        return f"0 0 {width:g} {height:g}"
    return "0 0 900 650"


def _set_style_properties(root: ET.Element) -> None:
    properties: dict[str, str] = {}
    for item in str(root.get("style") or "").split(";"):
        if ":" not in item:
            continue
        key, value = item.split(":", 1)
        properties[key.strip().lower()] = value.strip()
    properties.update(
        {
            "width": "100vw",
            "height": "100vh",
            "display": "block",
            "background": "#ffffff",
        }
    )
    root.set("style", "; ".join(f"{key}: {value}" for key, value in properties.items()))


def _apply_viewport_root(root: ET.Element) -> None:
    root.set("width", "100vw")
    root.set("height", "100vh")
    root.set("preserveAspectRatio", "xMidYMid meet")
    _set_style_properties(root)


def enrich_svg(
    svg_source: str,
    context: InteractiveSvgContext | None = None,
) -> str:
    """Return a standalone interactive SVG source string."""

    context = context or InteractiveSvgContext()
    try:
        root = ET.fromstring(svg_source)
    except ET.ParseError as exc:
        raise GbdrawError(f"Malformed SVG source for interactive export: {exc}") from exc
    if _local_name(root.tag) != "svg":
        raise GbdrawError("Interactive SVG export expected an SVG root element.")

    popup_mode: Literal["rich", "simple"] = (
        "simple" if context.popup_mode == "simple" else "rich"
    )
    original_viewbox = _original_viewbox(root)
    original_width = root.get("data-gbdraw-original-width") or root.get("width") or "900px"
    original_height = root.get("data-gbdraw-original-height") or root.get("height") or "650px"

    style_source, script_source = _load_standalone_assets()
    _remove_existing_assets(root)
    _remove_class_token(root, "gbdraw-feature-search-active")
    _remove_class_token(root, "gbdraw-feature-search-updating")
    for element in root.iter():
        _remove_class_token(element, "gbdraw-interactive-feature--dimmed")
        _remove_class_token(element, "gbdraw-interactive-feature--match")
        _remove_class_token(element, "gbdraw-interactive-feature--active-match")
    features = _feature_payloads(root, context)
    orthogroups = _orthogroup_payloads(features, context)
    matches = _match_payloads(root, features, orthogroups)
    sequence_sources = _sequence_sources_for_matches(matches, context.sequence_sources)
    annotation_context = {
        (
            str(item.get("record_index", "")),
            str(item.get("track_id", "")),
            str(item.get("set_id", "")),
            str(item.get("id", "")),
        ): dict(item)
        for item in context.annotations
    }
    annotations: list[dict[str, object]] = []
    for element in root.iter():
        annotation_id = str(element.get("data-gbdraw-annotation-id") or "")
        if not annotation_id:
            continue
        key = (
            str(element.get("data-gbdraw-record-index") or ""),
            str(element.get("data-gbdraw-annotation-track-id") or ""),
            str(element.get("data-gbdraw-annotation-set-id") or ""),
            annotation_id,
        )
        payload = annotation_context.get(key, annotation_context.get((key[0], "", key[2], key[3]), {})).copy()
        payload.update(
            {
                "dom_id": str(element.get("id") or ""),
                "id": annotation_id,
                "set_id": key[2],
                "track_id": key[1],
                "record_id": str(element.get("data-gbdraw-record-id") or ""),
                "record_index": int(key[0]) if key[0].isdigit() else key[0],
                "mark": str(element.get("data-gbdraw-annotation-mark") or payload.get("mark") or ""),
                "label": str(element.get("data-gbdraw-annotation-label") or payload.get("label") or ""),
            }
        )
        annotations.append(payload)

    try:
        metadata_payload = json.dumps(
            _compact_wire_value(
                {
                    "schema": INTERACTIVE_SCHEMA,
                    "popup_mode": popup_mode,
                    "features": features,
                    "orthogroups": orthogroups,
                    "matches": matches,
                    "annotations": annotations,
                    "sequence_sources": sequence_sources,
                }
            ),
            ensure_ascii=False,
            separators=(",", ":"),
        )
    except (TypeError, ValueError) as exc:
        raise GbdrawError("Interactive SVG metadata payload could not be serialized.") from exc

    _ensure_glow_filters(root)

    metadata = ET.SubElement(root, _svg_tag("metadata"))
    metadata.set("id", INTERACTIVE_METADATA_ID)
    metadata.set("data-schema", INTERACTIVE_SCHEMA)
    metadata.set("data-popup-mode", popup_mode)
    metadata.text = metadata_payload

    style = ET.SubElement(root, _svg_tag("style"))
    style.set("id", INTERACTIVE_STYLE_ID)
    style.set("type", "text/css")
    style.text = style_source

    script = ET.SubElement(root, _svg_tag("script"))
    script.set("id", INTERACTIVE_SCRIPT_ID)
    script.set("type", "application/ecmascript")
    script.text = script_source

    root.set("data-gbdraw-original-viewbox", original_viewbox)
    root.set("data-gbdraw-original-width", original_width)
    root.set("data-gbdraw-original-height", original_height)
    root.set("data-gbdraw-interactive-svg", "true")
    _apply_viewport_root(root)
    return ET.tostring(root, encoding="unicode")


__all__ = ["InteractiveSvgContext", "enrich_svg"]
