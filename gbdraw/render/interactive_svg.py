"""Standalone interactive SVG enrichment for CLI and API exports."""

from __future__ import annotations

import json
import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field, replace
from importlib import resources
from typing import Callable, Literal, Mapping, Sequence

from gbdraw.exceptions import GbdrawError

SVG_NS = "http://www.w3.org/2000/svg"
XLINK_NS = "http://www.w3.org/1999/xlink"
EV_NS = "http://www.w3.org/2001/xml-events"

INTERACTIVE_METADATA_ID = "gbdraw-interactive-feature-metadata"
INTERACTIVE_STYLE_ID = "gbdraw-interactive-feature-style"
INTERACTIVE_SCRIPT_ID = "gbdraw-interactive-feature-script"
INTERACTIVE_GLOW_FILTER_ID = "gbdraw-interactive-feature-glow"
INTERACTIVE_MATCH_GLOW_FILTER_ID = "gbdraw-interactive-feature-match-glow"
INTERACTIVE_SCHEMA = 3

_FEATURE_ELEMENT_SUFFIX_RE = re.compile(r"__(?:part|line)\d+$")
_FEATURE_CONNECTOR_SUFFIX_RE = re.compile(r"__line\d+$")
_FEATURE_RECORD_SUFFIX_RE = re.compile(r"_record_(\d+)$")
_FEATURE_INSTANCE_SUFFIX_RE = re.compile(r"__instance_.+_[0-9a-f]{16}$")
_NONNEGATIVE_INTEGER_RE = re.compile(r"\d+")
_FEATURE_RECORD_INDEX_KEYS = (
    "record_idx",
    "record_index",
    "recordIndex",
    "fileIdx",
    "file_idx",
)
_FEATURE_STABLE_ID_KEYS = (
    "stable_feature_id",
    "stable_svg_id",
    "stableFeatureSvgId",
    "stable_feature_svg_id",
    "stableFeatureId",
    "stableSvgId",
    "feature_svg_id",
    "featureSvgId",
)
_FEATURE_RENDERED_ID_KEYS = (
    "rendered_feature_svg_id",
    "renderedFeatureSvgId",
    "rendered_svg_id",
    "renderedSvgId",
)
_FEATURE_SOURCE_INDEX_KEYS = (
    "feature_index",
    "featureIndex",
    "source_feature_index",
    "sourceFeatureIndex",
)
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
    biological_features: Sequence[Mapping[str, object]] = ()
    record_keys: Sequence[str] = ()
    collinearity_search_scope: Literal["adjacent", "all"] | None = None


@dataclass
class _RenderedFeatureEntry:
    svg_id: str
    element: ET.Element
    fill: str
    stable_id: str
    stable_id_supplied: bool
    record_index: int | None
    record_index_supplied: bool
    source_index: int | None


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
        element.get("data-gbdraw-rendered-feature-id")
        or element.get("data-gbdraw-feature-id")
        or element.get("id")
        or ""
    )


def _feature_rendered_id_status(
    feature: Mapping[str, object],
) -> tuple[str, bool]:
    value, valid = _consistent_text_alias(
        feature,
        _FEATURE_RENDERED_ID_KEYS,
        normalize=_normalize_feature_id,
    )
    return value, valid


def _feature_stable_id(feature: Mapping[str, object]) -> str:
    stable_id, valid = _feature_stable_id_status(feature)
    return stable_id if valid else ""


def _consistent_text_alias(
    payload: Mapping[str, object],
    keys: Sequence[str],
    *,
    normalize: Callable[[object], str] = lambda value: str(value).strip(),
) -> tuple[str, bool]:
    values = {
        normalize(payload.get(key))
        for key in keys
        if key in payload
        and payload.get(key) is not None
        and str(payload.get(key)).strip()
    }
    values.discard("")
    if len(values) > 1:
        return "", False
    return (next(iter(values), ""), True)


def _feature_stable_id_status(
    feature: Mapping[str, object],
) -> tuple[str, bool]:
    stable_id, valid = _consistent_text_alias(
        feature,
        _FEATURE_STABLE_ID_KEYS,
        normalize=lambda value: _FEATURE_RECORD_SUFFIX_RE.sub(
            "",
            _normalize_feature_id(value),
        ),
    )
    if not valid or stable_id:
        return stable_id, valid
    rendered_id, rendered_valid = _feature_rendered_id_status(feature)
    if not rendered_valid:
        return "", False
    if rendered_id:
        return "", True
    return _consistent_text_alias(
        feature,
        ("svg_id", "svgId"),
        normalize=lambda value: _FEATURE_RECORD_SUFFIX_RE.sub(
            "",
            _normalize_feature_id(value),
        ),
    )


def _strict_nonnegative_integer(value: object) -> int | None:
    if isinstance(value, bool):
        return None
    text = str(value).strip()
    if not _NONNEGATIVE_INTEGER_RE.fullmatch(text):
        return None
    return int(text)


def _consistent_nonnegative_integer_alias(
    payload: Mapping[str, object],
    keys: Sequence[str],
) -> tuple[int | None, bool]:
    values: set[int] = set()
    for key in keys:
        if key not in payload:
            continue
        raw = payload.get(key)
        if raw is None or str(raw).strip() == "":
            continue
        parsed = _strict_nonnegative_integer(raw)
        if parsed is None:
            return None, False
        values.add(parsed)
    if len(values) > 1:
        return None, False
    return (next(iter(values), None), True)


def _feature_record_index_status(
    feature: Mapping[str, object],
) -> tuple[int | None, bool]:
    return _consistent_nonnegative_integer_alias(
        feature,
        _FEATURE_RECORD_INDEX_KEYS,
    )


def _feature_source_index_status(
    feature: Mapping[str, object],
) -> tuple[int | None, bool]:
    return _consistent_nonnegative_integer_alias(
        feature,
        _FEATURE_SOURCE_INDEX_KEYS,
    )


def _feature_record_index(feature: Mapping[str, object]) -> int | None:
    value, valid = _feature_record_index_status(feature)
    return value if valid else None


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


def _element_match_id_status(element: ET.Element) -> tuple[str, bool]:
    return _consistent_text_alias(
        {
            "match_id": element.get("data-gbdraw-match-id"),
            "pairwise_match_id": element.get("data-gbdraw-pairwise-match-id"),
        },
        ("match_id", "pairwise_match_id"),
    )


def _catalog_match_id_status(match: Mapping[str, object]) -> tuple[str, bool]:
    return _consistent_text_alias(match, ("id", "matchId", "match_id"))


def _rendered_element_identity(
    rendered_id: str,
    element: ET.Element,
) -> tuple[str, bool, int | None, bool, int | None]:
    raw_record_index = element.get("data-gbdraw-record-index")
    record_index_supplied = raw_record_index not in {None, ""}
    record_index = (
        _strict_nonnegative_integer(raw_record_index)
        if raw_record_index not in {None, ""}
        else None
    )
    if raw_record_index not in {None, ""} and record_index is None:
        raise GbdrawError(
            f"Rendered feature {rendered_id!r} has an invalid record index."
        )
    record_suffix = _FEATURE_RECORD_SUFFIX_RE.search(rendered_id)
    suffix_record_index = (
        int(record_suffix.group(1)) - 1
        if record_suffix is not None
        else None
    )
    if suffix_record_index is not None and suffix_record_index < 0:
        raise GbdrawError(
            f"Rendered feature {rendered_id!r} has an invalid record suffix."
        )
    if (
        record_index is not None
        and suffix_record_index is not None
        and record_index != suffix_record_index
    ):
        raise GbdrawError(
            f"Rendered feature {rendered_id!r} has conflicting record identity."
        )
    if record_index is None:
        record_index = suffix_record_index

    raw_source_indexes = [
        element.get("data-gbdraw-feature-index"),
        element.get("data-gbdraw-source-feature-index"),
    ]
    source_indexes = {
        parsed
        for raw in raw_source_indexes
        if raw not in {None, ""}
        for parsed in [_strict_nonnegative_integer(raw)]
        if parsed is not None
    }
    if any(
        raw not in {None, ""} and _strict_nonnegative_integer(raw) is None
        for raw in raw_source_indexes
    ) or len(source_indexes) > 1:
        raise GbdrawError(
            f"Rendered feature {rendered_id!r} has invalid source feature identity."
        )
    raw_stable_id = _normalize_feature_id(
        element.get("data-gbdraw-stable-feature-id")
    )
    stable_id = raw_stable_id or _FEATURE_RECORD_SUFFIX_RE.sub("", rendered_id)
    return (
        stable_id,
        bool(raw_stable_id),
        record_index,
        record_index_supplied,
        next(iter(source_indexes), None),
    )


def _collect_rendered_features(root: ET.Element) -> dict[str, _RenderedFeatureEntry]:
    entries: dict[str, _RenderedFeatureEntry] = {}
    for element in root.iter():
        if not _is_feature_candidate(element):
            continue
        svg_id = _element_feature_id(element)
        if not svg_id:
            continue
        (
            stable_id,
            stable_id_supplied,
            record_index,
            record_index_supplied,
            source_index,
        ) = _rendered_element_identity(svg_id, element)
        existing = entries.get(svg_id)
        if existing is not None:
            if (
                (
                    existing.stable_id_supplied
                    and stable_id_supplied
                    and existing.stable_id != stable_id
                )
                or (
                    existing.record_index_supplied
                    and record_index_supplied
                    and existing.record_index != record_index
                )
                or (
                    existing.source_index is not None
                    and source_index is not None
                    and existing.source_index != source_index
                )
            ):
                raise GbdrawError(
                    f"Rendered feature SVG ID {svg_id!r} has conflicting DOM identity."
                )
            if existing.stable_id_supplied:
                stable_id = existing.stable_id
            stable_id_supplied = existing.stable_id_supplied or stable_id_supplied
            record_index = existing.record_index if existing.record_index is not None else record_index
            record_index_supplied = (
                existing.record_index_supplied or record_index_supplied
            )
            source_index = existing.source_index if existing.source_index is not None else source_index
            if _is_feature_fill_target(existing.element) or not _is_feature_fill_target(element):
                element = existing.element
        entries[svg_id] = _RenderedFeatureEntry(
            svg_id=svg_id,
            element=element,
            fill=str(element.get("fill") or "#94a3b8"),
            stable_id=stable_id,
            stable_id_supplied=stable_id_supplied,
            record_index=record_index,
            record_index_supplied=record_index_supplied,
            source_index=source_index,
        )
    return entries


def _rendered_space_stable_id_candidates(rendered_id: str) -> set[str]:
    """Return stable-looking bases encoded by one rendered-space handle."""

    base = _FEATURE_RECORD_SUFFIX_RE.sub("", _normalize_feature_id(rendered_id))
    candidates = {base} if base else set()
    instance_base = _FEATURE_INSTANCE_SUFFIX_RE.sub("", base)
    if instance_base:
        candidates.add(instance_base)
    return candidates


def _rendered_entry_agrees_with_identity(
    entry: _RenderedFeatureEntry,
    *,
    rendered_id: str,
    record_index: int | None,
    stable_id: str,
    source_index: int | None = None,
    allow_rendered_space_stable_id: bool = False,
) -> bool:
    if (
        record_index is not None
        and entry.record_index is not None
        and record_index != entry.record_index
    ):
        return False
    if (
        stable_id
        and entry.stable_id_supplied
        and entry.stable_id != stable_id
        and not (
            allow_rendered_space_stable_id
            and entry.stable_id
            in _rendered_space_stable_id_candidates(rendered_id)
        )
    ):
        return False
    return not (
        source_index is not None
        and entry.source_index is not None
        and source_index != entry.source_index
    )


def _rendered_feature_id_index(
    root: ET.Element,
) -> dict[tuple[str, int | None, str], str]:
    """Index only feature IDs that are present in the final SVG DOM."""

    index: dict[tuple[str, int | None, str], str] = {}
    ambiguous: set[tuple[str, int | None, str]] = set()

    def add(key: tuple[str, int | None, str], rendered_id: str) -> None:
        if key in ambiguous:
            return
        existing = index.get(key)
        if existing is not None and existing != rendered_id:
            index.pop(key, None)
            ambiguous.add(key)
            return
        index[key] = rendered_id

    for rendered_id, entry in _collect_rendered_features(root).items():
        if entry.stable_id:
            add(("stable", entry.record_index, entry.stable_id), rendered_id)
        add(("rendered", entry.record_index, rendered_id), rendered_id)
    return index


def _resolve_rendered_feature_id(
    rendered_index: Mapping[tuple[str, int | None, str], str],
    record_index: int | None,
    stable_feature_id: str,
) -> str:
    if not stable_feature_id:
        return ""
    return _first_text(
        rendered_index.get(("stable", record_index, stable_feature_id))
    )


def _resolve_explicit_rendered_feature_id(
    rendered_index: Mapping[tuple[str, int | None, str], str],
    record_index: int | None,
    rendered_feature_id: str,
) -> str:
    if not rendered_feature_id:
        return ""
    return _first_text(
        rendered_index.get(("rendered", record_index, rendered_feature_id))
    )


def _resolve_explicit_recordless_rendered_id(
    rendered_index: Mapping[tuple[str, int | None, str], str],
    rendered_feature_id: str,
) -> str:
    """Bind an exact rendered handle when legacy DOM omits record metadata."""

    if not rendered_feature_id:
        return ""
    return _first_text(
        rendered_index.get(("rendered", None, rendered_feature_id))
    )


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
    return _strict_nonnegative_integer(value)


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


def _normalize_orthogroup_member(
    member: object | None,
    *,
    rendered_feature_svg_id: str = "",
) -> dict[str, object] | None:
    if not isinstance(member, Mapping):
        return None
    record_index = _int_or_none(_first_text(member.get("recordIndex"), member.get("record_index")))
    feature_index = _number_or_none(_first_text(member.get("featureIndex"), member.get("feature_index")))
    start = _number_or_none(member.get("start"))
    end = _number_or_none(member.get("end"))
    stable_feature_svg_id = _first_text(
        member.get("stableFeatureSvgId"),
        member.get("stable_feature_svg_id"),
        member.get("featureSvgId"),
        member.get("feature_svg_id"),
    )
    payload = {
        "orthogroup_id": _first_text(member.get("orthogroupId"), member.get("orthogroup_id")),
        "protein_id": _first_text(member.get("proteinId"), member.get("protein_id")),
        "source_protein_id": _first_text(member.get("sourceProteinId"), member.get("source_protein_id")),
        "record_index": record_index,
        "record_id": _first_text(member.get("recordId"), member.get("record_id")),
        "feature_index": feature_index,
        "label": _first_text(member.get("label")),
        "feature_svg_id": stable_feature_svg_id,
        "stable_feature_svg_id": stable_feature_svg_id,
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
    if rendered_feature_svg_id:
        payload["rendered_feature_svg_id"] = rendered_feature_svg_id
    return payload


def _fallback_feature_payload(
    svg_id: str,
    entry: _RenderedFeatureEntry,
    *,
    feature_index: int,
) -> dict[str, object]:
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
        "record_idx": _int_or_none(
            entry.element.get("data-gbdraw-record-index")
        ),
        "feature_index": feature_index,
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
    rendered_index = _rendered_feature_id_index(root)

    payloads: list[dict[str, object]] = []
    seen: set[str] = set()
    stable_ids_by_rendered_id: dict[str, str] = {}
    popup_mode: Literal["rich", "simple"] = (
        "simple" if context.popup_mode == "simple" else "rich"
    )
    for feature in context.features:
        record_index, record_index_valid = _feature_record_index_status(feature)
        stable_feature_id, stable_feature_id_valid = _feature_stable_id_status(
            feature
        )
        rendered_feature_id, rendered_feature_id_valid = (
            _feature_rendered_id_status(feature)
        )
        generic_svg_id, generic_svg_id_valid = _consistent_text_alias(
            feature,
            ("svg_id", "svgId"),
            normalize=_normalize_feature_id,
        )
        feature_index, feature_index_valid = _feature_source_index_status(feature)
        if not all(
            (
                record_index_valid,
                stable_feature_id_valid,
                rendered_feature_id_valid,
                generic_svg_id_valid,
                feature_index_valid,
            )
        ):
            raise GbdrawError(
                "Feature metadata contains conflicting or invalid identity aliases."
            )
        rendered_handle = rendered_feature_id or generic_svg_id
        mapped_svg_id = (
            _first_text(
                _resolve_explicit_rendered_feature_id(
                    rendered_index,
                    record_index,
                    rendered_handle,
                ),
                _resolve_explicit_recordless_rendered_id(
                    rendered_index,
                    rendered_handle,
                ),
            )
            if rendered_handle
            else _resolve_rendered_feature_id(
                rendered_index,
                record_index,
                stable_feature_id,
            )
        )
        mapped_entry = rendered.get(mapped_svg_id)
        if mapped_entry is not None and not _rendered_entry_agrees_with_identity(
            mapped_entry,
            rendered_id=mapped_svg_id,
            record_index=record_index,
            stable_id=stable_feature_id,
            source_index=feature_index,
            allow_rendered_space_stable_id=True,
        ):
            raise GbdrawError(
                f"Feature metadata identity does not agree with rendered SVG ID "
                f"{mapped_svg_id!r}."
            )
        if mapped_svg_id and mapped_svg_id in seen:
            raise GbdrawError(
                f"Multiple feature metadata entries resolve to rendered SVG ID "
                f"{mapped_svg_id!r}."
            )
        svg_id = (
            mapped_svg_id
            if mapped_svg_id in rendered and mapped_svg_id not in seen
            else ""
        )
        if not svg_id or svg_id not in rendered or svg_id in seen:
            continue
        seen.add(svg_id)
        if stable_feature_id:
            stable_ids_by_rendered_id[svg_id] = stable_feature_id
        fallback_label = _get_feature_label(feature)
        display_label = fallback_label
        orthogroup_entry = _feature_orthogroup_entry(feature)
        orthogroup_member = _normalize_orthogroup_member(
            orthogroup_entry.get("orthogroupMember")
            if isinstance(orthogroup_entry, Mapping)
            else None,
            rendered_feature_svg_id=svg_id,
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
            "rendered_feature_svg_id": svg_id,
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
            "record_idx": record_index,
            "feature_index": feature_index,
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

    for element in root.iter():
        stable_feature_id = stable_ids_by_rendered_id.get(_element_feature_id(element))
        if stable_feature_id:
            element.set("data-gbdraw-stable-feature-id", stable_feature_id)

    for fallback_index, (svg_id, entry) in enumerate(rendered.items()):
        if svg_id in seen:
            continue
        seen.add(svg_id)
        payloads.append(
            dict(
                _compact_wire_value(
                    _fallback_feature_payload(
                        svg_id,
                        entry,
                        feature_index=fallback_index,
                    )
                )
                or {}
            )
        )

    feature_ids = {str(feature.get("svg_id") or "").strip() for feature in payloads}
    for element in root.iter():
        if not _is_feature_candidate(element):
            continue
        if _element_feature_id(element) not in feature_ids:
            continue
        element.set("data-gbdraw-interactive-feature", "true")
        _add_class_token(element, "gbdraw-interactive-feature")
    return payloads


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


def _validate_catalog_feature_bindings(
    root: ET.Element,
    catalog_item: Mapping[str, object],
) -> dict[str, _RenderedFeatureEntry]:
    """Require one exact, record-scoped catalog binding per rendered feature."""

    record_keys_value = catalog_item.get("recordKeys")
    biological_value = catalog_item.get("biologicalFeatures")
    rendered_value = catalog_item.get("features")
    if any(
        not isinstance(value, Sequence) or isinstance(value, (str, bytes))
        for value in (record_keys_value, biological_value, rendered_value)
    ):
        raise GbdrawError("Interactive feature catalog is missing feature arrays.")

    record_keys = [str(value).strip() for value in record_keys_value]
    if any(not value for value in record_keys) or len(set(record_keys)) != len(
        record_keys
    ):
        raise GbdrawError(
            "Interactive feature catalog record keys must be non-empty and unique."
        )
    record_indexes = {record_key: index for index, record_key in enumerate(record_keys)}

    biological_identities: dict[tuple[str, str], tuple[str, int | None]] = {}
    for feature in biological_value:
        if not isinstance(feature, Mapping):
            raise GbdrawError(
                "Interactive feature catalog biological features must be objects."
            )
        record_key, record_key_valid = _consistent_text_alias(
            feature,
            ("recordKey", "record_key"),
        )
        biological_id, biological_id_valid = _consistent_text_alias(
            feature,
            ("biologicalFeatureId", "biological_feature_id"),
        )
        stable_id, stable_id_valid = _consistent_text_alias(
            feature,
            _FEATURE_STABLE_ID_KEYS,
        )
        source_index, source_index_valid = _feature_source_index_status(feature)
        key = (record_key, biological_id)
        if (
            not record_key_valid
            or not biological_id_valid
            or not stable_id_valid
            or not source_index_valid
            or not all(key)
            or record_key not in record_indexes
            or key in biological_identities
        ):
            raise GbdrawError(
                "Interactive feature catalog contains invalid or duplicate "
                "biological feature identity."
            )
        biological_identities[key] = (stable_id or biological_id, source_index)

    catalog_features: dict[
        str,
        tuple[int, str, int | None],
    ] = {}
    for feature in rendered_value:
        if not isinstance(feature, Mapping):
            raise GbdrawError(
                "Interactive feature catalog rendered features must be objects."
            )
        svg_id, svg_id_valid = _consistent_text_alias(
            feature,
            ("svgId", "svg_id", *_FEATURE_RENDERED_ID_KEYS),
        )
        record_key, record_key_valid = _consistent_text_alias(
            feature,
            ("recordKey", "record_key"),
        )
        biological_id, biological_id_valid = _consistent_text_alias(
            feature,
            ("biologicalFeatureId", "biological_feature_id"),
        )
        reference = (record_key, biological_id)
        biological_identity = biological_identities.get(reference)
        if (
            not svg_id_valid
            or not record_key_valid
            or not biological_id_valid
            or not svg_id
            or svg_id != _normalize_feature_id(svg_id)
            or svg_id in catalog_features
            or biological_identity is None
        ):
            raise GbdrawError(
                "Interactive feature catalog contains an invalid or duplicate "
                "rendered feature identity."
            )
        stable_id, source_index = biological_identity
        catalog_features[svg_id] = (
            record_indexes[record_key],
            stable_id,
            source_index,
        )

    rendered_entries = _collect_rendered_features(root)
    if catalog_features.keys() != rendered_entries.keys():
        raise GbdrawError(
            "Interactive feature catalog rendered feature IDs do not match the SVG."
        )
    for rendered_id, (record_index, stable_id, source_index) in catalog_features.items():
        if not _rendered_entry_agrees_with_identity(
            rendered_entries[rendered_id],
            rendered_id=rendered_id,
            record_index=record_index,
            stable_id=stable_id,
            source_index=source_index,
        ):
            raise GbdrawError(
                "Interactive feature catalog identity does not agree with rendered "
                f"SVG ID {rendered_id!r}."
            )
    return rendered_entries


def enrich_svg(
    svg_source: str,
    context: InteractiveSvgContext | None = None,
    *,
    result_index: int = 0,
    result_name: str = "interactive.svg",
    feature_catalog: Mapping[str, object] | None = None,
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
    from gbdraw.web_support.feature_catalog import (
        build_feature_catalog,
        build_feature_catalog_item,
        select_feature_catalog_item,
    )

    if feature_catalog is not None:
        catalog_item = select_feature_catalog_item(
            feature_catalog,
            result_index=result_index,
            result_name=result_name,
        )
        _validate_catalog_feature_bindings(root, catalog_item)
    else:
        if not context.biological_features:
            metadata_free_features = not context.features and not context.record_keys
            fallback_features = [
                dict(feature)
                for feature in (
                    context.features
                    if context.features
                    else _feature_payloads(root, context)
                )
            ]
            synthesize_local_record = (
                metadata_free_features
                and bool(fallback_features)
                and all(
                    _feature_record_index(feature) is None
                    for feature in fallback_features
                )
            )
            if synthesize_local_record:
                for feature in fallback_features:
                    feature["record_idx"] = 0
                    feature["rendered_feature_svg_id"] = _first_text(
                        feature.get("svg_id")
                    )
            context = replace(
                context,
                features=fallback_features,
                biological_features=fallback_features,
                record_keys=("record-1",) if synthesize_local_record else context.record_keys,
            )
        catalog_item = build_feature_catalog_item(
            ET.tostring(root, encoding="unicode"),
            context,
            result_index=result_index,
            result_name=result_name,
        )
    catalog = build_feature_catalog([catalog_item])

    biological_by_key = {
        (
            str(feature.get("recordKey") or ""),
            str(feature.get("biologicalFeatureId") or ""),
        ): feature
        for feature in catalog_item["biologicalFeatures"]
        if isinstance(feature, Mapping)
    }
    rendered_features_by_id = {
        str(feature.get("svgId") or ""): feature
        for feature in catalog_item["features"]
        if isinstance(feature, Mapping)
    }
    for element in root.iter():
        if not _is_feature_candidate(element):
            continue
        rendered = rendered_features_by_id.get(_element_feature_id(element))
        if rendered is None:
            continue
        biological = biological_by_key.get(
            (
                str(rendered.get("recordKey") or ""),
                str(rendered.get("biologicalFeatureId") or ""),
            )
        )
        stable_id = (
            str(
                biological.get("stableFeatureId")
                or biological.get("biologicalFeatureId")
                or ""
            )
            if biological is not None
            else ""
        )
        if stable_id:
            element.set("data-gbdraw-stable-feature-id", stable_id)
        element.set("data-gbdraw-interactive-feature", "true")
        _add_class_token(element, "gbdraw-interactive-feature")

    matches = [
        match
        for match in catalog_item["comparisonMatches"]
        if isinstance(match, Mapping)
    ]
    match_elements = [element for element in root.iter() if _is_match_candidate(element)]
    matches_by_id: dict[str, Mapping[str, object]] = {}
    for match in matches:
        match_id, valid = _catalog_match_id_status(match)
        if not valid or not match_id or match_id in matches_by_id:
            raise GbdrawError(
                "Interactive feature catalog contains invalid or duplicate match IDs."
            )
        matches_by_id[match_id] = match
    elements_by_match_id: dict[str, ET.Element] = {}
    for element in match_elements:
        match_id, valid = _element_match_id_status(element)
        if not valid or not match_id or match_id in elements_by_match_id:
            raise GbdrawError(
                "Rendered SVG contains invalid or duplicate match IDs."
            )
        elements_by_match_id[match_id] = element
    if matches_by_id.keys() != elements_by_match_id.keys():
        raise GbdrawError(
            "Interactive feature catalog match IDs do not match the rendered SVG."
        )
    for match_id, element in elements_by_match_id.items():
        element.set("data-gbdraw-match-id", match_id)
        element.set("data-gbdraw-pairwise-match-id", match_id)
        element.set("data-gbdraw-interactive-match", "true")
        _add_class_token(element, "gbdraw-interactive-pairwise-match")

    try:
        metadata_payload = json.dumps(
            catalog,
            ensure_ascii=False,
            separators=(",", ":"),
        )
    except (TypeError, ValueError) as exc:
        raise GbdrawError("Interactive SVG metadata payload could not be serialized.") from exc

    _ensure_glow_filters(root)

    metadata = ET.SubElement(root, _svg_tag("metadata"))
    metadata.set("id", INTERACTIVE_METADATA_ID)
    metadata.set("data-schema", str(INTERACTIVE_SCHEMA))
    metadata.set("data-popup-mode", popup_mode)
    metadata.set("data-result-index", str(result_index))
    metadata.set("data-result-name", result_name)
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
