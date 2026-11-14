"""Compact feature-catalog normalization for browser render responses."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, MutableMapping, Sequence
import copy
from dataclasses import dataclass, fields as dataclass_fields, replace
import json
import math
import re
from types import MappingProxyType
import xml.etree.ElementTree as ET
from typing import TYPE_CHECKING, Any, Iterable, Literal, TypedDict, TypeVar, cast, final, get_args

from svgwrite import Drawing

from gbdraw.core.numeric import JS_MAX_SAFE_INTEGER, ascii_float, parse_ascii_nonnegative_integer
from gbdraw.render.document import DocumentParts
from gbdraw.render.drawers.annotation_marks import ANNOTATION_ID_ATTRIBUTE, ANNOTATION_SET_ID_ATTRIBUTE
from gbdraw.render.group_content import GroupContent, GroupKind, held_group, leaf_groups
from gbdraw.render_stages import render_stage
from gbdraw.exceptions import GbdrawError
from gbdraw.features.ids import disambiguate_feature_ids
from gbdraw.render.interactive_svg import (
    CatalogSequenceRecords,
    InteractiveSvgContext,
    catalog_record_units,
    _RenderedFeatureEntry,
    _compact_wire_value,
    _feature_record_index_status,
    _feature_source_index_status,
    _feature_stable_id,
    _feature_stable_id_status,
    _first_text,
    _is_feature_candidate,
    _is_match_candidate,
    _match_id,
    _merge_rendered_feature_entry,
    _normalize_string_array,
    _rendered_feature_entry,
    _resolve_rendered_features,
    _validate_match_fragments,
)
from gbdraw.session_io import (
    CURRENT_FEATURE_CATALOG_SCHEMA as FEATURE_CATALOG_SCHEMA,
    SEQUENCE_FEATURE_CATALOG_SCHEMA,
)

if TYPE_CHECKING:
    from Bio.SeqRecord import SeqRecord

    from gbdraw.web_support.feature_metadata import FeatureSequences

# Schema 4 added source-anchor profiles; schema 5 adds each rendered feature's
# drawn selector values (``drawnSelector``) for live rule matching.
ANCHOR_FEATURE_CATALOG_SCHEMA = 4
LEGACY_FEATURE_CATALOG_SCHEMA = 3
_DRAWN_SELECTOR_FIELDS = frozenset({"hash", "location", "recordLocation"})

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
_BIOLOGICAL_REMOVED_KEYS = frozenset(
    _BIOLOGICAL_ALIAS_KEYS
    | _ORTHOGROUP_FEATURE_KEYS
    | {
        "selector",
        "drawn_selector",
        "nucleotide_fasta",
        "nucleotideFasta",
        "amino_acid_fasta",
        "aminoAcidFasta",
        "record_idx",
        "recordIndex",
        "record_index",
        "feature_index",
        "featureIndex",
        "source_feature_index",
        "sourceFeatureIndex",
    }
)
_SOURCE_IDENTITY_CONFLICT_KEYS = frozenset(
    {
        "record_index",
        "recordIndex",
        "fileIdx",
        "file_idx",
        "featureIndex",
        "source_feature_index",
        "sourceFeatureIndex",
        "stableFeatureSvgId",
        "stable_feature_svg_id",
        "stableFeatureId",
        "stableSvgId",
        "feature_svg_id",
        "featureSvgId",
    }
)
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
    if type(value) is list or type(value) is tuple:
        return list(value)
    if value is not None and isinstance(value, Sequence) and not isinstance(
        value, (str, bytes)
    ):
        return list(value)
    return []


def _text(value: object | None) -> str:
    if type(value) is str:
        return value.strip()
    if value is None:
        return ""
    return _first_text(value)


def _compact_wire_mapping(value: Mapping[str, object]) -> dict[str, object]:
    # _compact_wire_value returns a mapping (or None when empty) for a mapping.
    compacted = cast("Mapping[str, object] | None", _compact_wire_value(value))
    return dict(compacted or {})


def _legacy_feature_strand(
    feature: Mapping[str, object],
    parts: Sequence[Mapping[str, object]],
) -> str:
    strands: set[str] = set()
    for part in parts or (feature,):
        value = part.get("strand")
        if value in (1, "1", "+"):
            strands.add("+")
        elif value in (-1, "-1", "-"):
            strands.add("-")
    if len(strands) > 1:
        return "mixed"
    if strands:
        return next(iter(strands))
    return "unstranded"


def _legacy_exact_coordinate(value: Any) -> int | None:
    # Legacy payload values are arbitrary JSON; non-numeric ones raise TypeError.
    if isinstance(value, bool):
        return None
    try:
        number = int(value)
    except (TypeError, ValueError, OverflowError):
        return None
    if isinstance(value, float) and not value.is_integer():
        return None
    if isinstance(value, str) and value.strip() != str(number):
        return None
    return number if abs(number) <= JS_MAX_SAFE_INTEGER else None


_MISSING = object()


def _valid_drawn_selector(value: object) -> bool:
    """A drawn selector is null (unknown) or the three drawn values."""

    if value is None:
        return True
    return (
        isinstance(value, Mapping)
        and set(value) == _DRAWN_SELECTOR_FIELDS
        and all(item is None or (isinstance(item, str) and item) for item in value.values())
    )


def _drawn_selector(feature: Mapping[str, object] | None) -> dict[str, str | None] | None:
    raw = feature.get("drawn_selector") if isinstance(feature, Mapping) else None
    if not isinstance(raw, Mapping):
        return None
    return {field: (_text(raw.get(field)) or None) for field in ("hash", "location", "recordLocation")}


# What a catalog of schema 5 or older stored of the sequences (D5): a schema 6 row or item holds none of it.
STORED_SEQUENCE_KEYS = frozenset({
    "nucleotide_sequence", "nucleotideSequence", "amino_acid_sequence", "aminoAcidSequence",
    "sequence_warnings", "sequenceWarnings", "sequenceSourceIndex", "translationFromAminoAcidSequence",
})


def _without_stored_sequences(item: dict[str, object]) -> None:
    """Drop what an older catalog item stored of the sequences, giving back the ``/translation`` qualifier that a
    ``translationFromAminoAcidSequence`` mark stood for (Python reads every sequence again, D5)."""

    item.pop("sequenceSources", None)
    item.pop("matchSpanReads", None)
    features = item.get("biologicalFeatures")
    for feature in features if isinstance(features, list) else ():
        if not isinstance(feature, dict):
            continue
        amino_acid_sequence = _text(feature.get("amino_acid_sequence") or feature.get("aminoAcidSequence"))
        if feature.get("translationFromAminoAcidSequence") is True and amino_acid_sequence:
            qualifiers = feature.get("qualifiers")
            feature["qualifiers"] = {
                **(qualifiers if isinstance(qualifiers, Mapping) else {}), "translation": [amino_acid_sequence],
            }
        for key in STORED_SEQUENCE_KEYS:
            feature.pop(key, None)


def promote_legacy_feature_catalog(
    catalog: Mapping[str, object],
    *,
    schema: int = FEATURE_CATALOG_SCHEMA,
) -> dict[str, object]:
    """Promote schema 3, 4 or 5 without inventing values, to ``schema``: the current one, or 5 for the Session 46-50
    an older Session moves through (its move to drawings promotes it again).

    Every older schema stored the sequences; the current one holds none (``_without_stored_sequences``).

    Schema 3 gains only safely inferred source-anchor metadata and source
    feature indexes: its writer named a duplicated stable identity
    ``<stableFeatureId>~<source feature index>`` and some writers stored no
    index (OV-269). Rendered features of schema 3 and 4 gain
    ``drawnSelector: None``: their drawn selector values are unknown until the
    next Generate.
    """

    source_schema = catalog.get("schema")
    if source_schema not in (
        LEGACY_FEATURE_CATALOG_SCHEMA, ANCHOR_FEATURE_CATALOG_SCHEMA, SEQUENCE_FEATURE_CATALOG_SCHEMA,
    ):
        raise GbdrawError(
            f"Feature catalog must use schema {LEGACY_FEATURE_CATALOG_SCHEMA}, "
            f"{ANCHOR_FEATURE_CATALOG_SCHEMA} or {SEQUENCE_FEATURE_CATALOG_SCHEMA}."
        )
    migrated = copy.deepcopy(dict(catalog))
    items = migrated.get("items")
    if not isinstance(items, list):
        raise GbdrawError("Feature catalog items must be an array.")
    if schema not in (SEQUENCE_FEATURE_CATALOG_SCHEMA, FEATURE_CATALOG_SCHEMA) or source_schema > schema:
        raise GbdrawError(f"Feature catalog schema {source_schema} does not promote to schema {schema}.")
    migrated["schema"] = schema
    for item in items:
        if not isinstance(item, dict):
            raise GbdrawError("Feature catalog items must contain objects.")
        # Item validation stays with select_feature_catalog_item.
        if schema == FEATURE_CATALOG_SCHEMA:
            _without_stored_sequences(item)
        if source_schema == SEQUENCE_FEATURE_CATALOG_SCHEMA:
            continue
        rendered = item.get("features")
        for feature in rendered if isinstance(rendered, list) else ():
            if isinstance(feature, dict):
                feature["drawnSelector"] = None
        if source_schema == ANCHOR_FEATURE_CATALOG_SCHEMA:
            continue
        features = item.get("biologicalFeatures")
        if not isinstance(features, list):
            raise GbdrawError(
                "Feature catalog item is missing biologicalFeatures."
            )
        for feature in features:
            if not isinstance(feature, dict):
                raise GbdrawError(
                    "Feature catalog biologicalFeatures must contain objects."
                )
            stable_id = _text(feature.get("stableFeatureId"))
            biological_id = _text(feature.get("biologicalFeatureId"))
            suffix = (
                biological_id[len(stable_id) + 1:]
                if stable_id and biological_id.startswith(f"{stable_id}~")
                else ""
            )
            suffix_index = parse_ascii_nonnegative_integer(suffix)
            if suffix_index is not None and _feature_source_index_status(feature) == (None, True):
                feature["sourceFeatureIndex"] = suffix_index
            raw_parts = feature.get("location_parts", feature.get("locationParts"))
            parts = (
                raw_parts
                if isinstance(raw_parts, list) and raw_parts
                else [feature]
            )
            mapping_parts = [part for part in parts if isinstance(part, Mapping)]
            exact_single = len(mapping_parts) == len(parts) == 1
            start = (
                _legacy_exact_coordinate(mapping_parts[0].get("start"))
                if exact_single
                else None
            )
            end = (
                _legacy_exact_coordinate(mapping_parts[0].get("end"))
                if exact_single
                else None
            )
            exact_single = exact_single and start is not None and end is not None
            strand = _legacy_feature_strand(feature, mapping_parts)
            feature["anchorProfile"] = (
                {
                    "precision": "exact",
                    "operator": "single",
                    "partOrder": (
                        "source-forward" if strand == "unstranded" else "biological"
                    ),
                    "strand": strand,
                }
                if exact_single
                else {
                    "precision": "unavailable",
                    "operator": "unknown",
                    "partOrder": "ambiguous",
                    "strand": strand,
                }
            )
            if (
                exact_single
                and "location_parts" not in feature
                and "locationParts" not in feature
            ):
                feature["location_parts"] = [
                    {"start": start, "end": end, "strand": strand}
                ]
    return migrated


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

_MATCH_SNAPSHOT_ATTRIBUTES = tuple(
    dict.fromkeys(
        (
            "id",
            "data-gbdraw-match-fragment",
            "data-gbdraw-match-geometry",
            "data-gbdraw-match-id",
            "data-match-kind",
            "data-collinearity-block-id",
            "data-orthogroup-id",
            "fill",
            "data-query-record-id",
            "data-query",
            "data-subject-record-id",
            "data-subject",
            "data-collinearity-orientation",
            "data-orientation",
            *_MATCH_ATTRIBUTES.values(),
        )
    )
)
_ANNOTATION_SNAPSHOT_ATTRIBUTES = (
    "id",
    ANNOTATION_ID_ATTRIBUTE,
    "data-gbdraw-record-index",
    "data-gbdraw-annotation-track-id",
    ANNOTATION_SET_ID_ATTRIBUTE,
    "data-gbdraw-record-id",
    "data-gbdraw-annotation-mark",
    "data-gbdraw-annotation-label",
)


@dataclass(frozen=True)
class _SvgCatalogCandidate:
    attributes: Mapping[str, str]


@dataclass(frozen=True)
class _RenderedSvgCatalogIndex:
    """Read-only catalog state of the trees the catalog reads, in document order."""

    rendered_features: Mapping[str, _RenderedFeatureEntry]
    match_candidates: tuple[_SvgCatalogCandidate, ...]
    annotation_candidates: tuple[_SvgCatalogCandidate, ...]
    record_source_spans: Mapping[int, tuple[int, int, int]]


@dataclass(frozen=True)
class _CatalogElementReads:
    """What the catalog reads from the elements of one tree, each list in document order.

    Equal reads give an equal catalog item: the item reads nothing else of the
    SVG (no position, path, box or stroke).
    """

    record_source_spans: tuple[tuple[int, tuple[int, int, int]], ...]
    feature_entries: tuple[_RenderedFeatureEntry, ...]
    match_attributes: tuple[Mapping[str, str], ...]
    annotation_attributes: tuple[Mapping[str, str], ...]


@dataclass(frozen=True)
class _CatalogTraversal:
    """One traversal of one parsed tree: its reads and what it counted."""

    reads: _CatalogElementReads
    dom_element_count: int
    feature_candidate_count: int


@dataclass(frozen=True)
class _SourceFeatureCandidate:
    feature: Mapping[str, object]
    record_index: int
    stable_id: str
    feature_index: int | None
    record_key: str
    identity_base: str


@dataclass(frozen=True)
class _IndexedBiologicalFeature:
    normalized: Mapping[str, object]
    # The drawn units that show the feature: one, or each window of its record that shows it.
    record_indexes: tuple[int, ...]
    stable_id: str
    feature_index: int | None


def _snapshot_attributes(
    element: ET.Element,
    names: Sequence[str],
) -> Mapping[str, str]:
    return MappingProxyType(
        {
            name: value
            for name in names
            if (value := element.get(name)) is not None
        }
    )


def _catalog_read_metrics(traversal: _CatalogTraversal) -> dict[str, int]:
    """The metrics of one parse and traversal (one per ``catalog_snapshot`` build)."""

    return {
        "featureCatalogSvgParseCount": 1,
        "featureCatalogFullDomTraversalCount": 1,
        "featureCatalogDomElementCount": traversal.dom_element_count,
        "featureCatalogFeatureCandidateCount": traversal.feature_candidate_count,
    }


def _catalog_index_metrics(index: _RenderedSvgCatalogIndex) -> dict[str, int]:
    """The metrics of one index (one per ``catalog_item`` build)."""

    return {
        "featureCatalogRenderedFeatureCollectionCount": 1,
        "featureCatalogCatalogOnlyDomMutationCount": 0,
        "featureCatalogUniqueRenderedFeatureCount": len(
            index.rendered_features
        ),
        "featureCatalogMatchCandidateCount": len(index.match_candidates),
        "featureCatalogAnnotationCandidateCount": len(
            index.annotation_candidates
        ),
    }


def _record_catalog_metrics(
    diagnostics: MutableMapping[str, Any] | None,
    values: Mapping[str, int],
) -> None:
    if diagnostics is None:
        return
    metrics = diagnostics.setdefault("metrics", {})
    if not isinstance(metrics, MutableMapping):
        raise GbdrawError("Feature catalog diagnostics metrics must be a mapping.")
    for name, value in values.items():
        metrics[name] = int(metrics.get(name, 0)) + value


def _parse_catalog_svg(svg_source: str | Iterable[str]) -> ET.Element:
    """The SVG tree of ``svg_source``, a text or texts in order; texts are fed in turn, never joined (OV-3027)."""

    try:
        if isinstance(svg_source, str):
            root = ET.fromstring(svg_source)
        else:
            parser = ET.XMLParser()
            for text in svg_source:
                parser.feed(text)
            root = parser.close()
    except ET.ParseError as exc:
        raise GbdrawError(
            f"Malformed SVG source for feature catalog: {exc}"
        ) from exc
    if root.tag.rsplit("}", 1)[-1] != "svg":
        raise GbdrawError("Feature catalog generation expected an SVG root.")
    return root


def _catalog_element_reads(root: ET.Element) -> _CatalogTraversal:
    """Read every element of one SVG tree once, in document order."""

    spans: list[tuple[int, tuple[int, int, int]]] = []
    entries: list[_RenderedFeatureEntry] = []
    matches: list[Mapping[str, str]] = []
    annotations: list[Mapping[str, str]] = []
    dom_element_count = 0
    feature_candidate_count = 0
    for element in root.iter():
        dom_element_count += 1
        span = _record_source_span(element)
        if span is not None:
            spans.append(span)
        if _is_feature_candidate(element):
            feature_candidate_count += 1
            entry = _rendered_feature_entry(element)
            if entry is not None:
                entries.append(entry)
        if _is_match_candidate(element):
            matches.append(_snapshot_attributes(element, _MATCH_SNAPSHOT_ATTRIBUTES))
        if element.get(ANNOTATION_ID_ATTRIBUTE):
            annotations.append(_snapshot_attributes(element, _ANNOTATION_SNAPSHOT_ATTRIBUTES))
    return _CatalogTraversal(
        _CatalogElementReads(tuple(spans), tuple(entries), tuple(matches), tuple(annotations)),
        dom_element_count=dom_element_count,
        feature_candidate_count=feature_candidate_count,
    )


def _rendered_svg_catalog_index(reads: Sequence[_CatalogElementReads]) -> _RenderedSvgCatalogIndex:
    """The catalog state of trees read in document order.

    Spans and candidates join in order; features merge in order (the first
    entry of an id wins), so reading sibling trees one by one gives the index
    of one traversal over all of them.
    """

    rendered_features: dict[str, _RenderedFeatureEntry] = {}
    record_source_spans: dict[int, tuple[int, int, int]] = {}
    for part in reads:
        record_source_spans.update(part.record_source_spans)
        for entry in part.feature_entries:
            _merge_rendered_feature_entry(rendered_features, entry)
    return _RenderedSvgCatalogIndex(
        rendered_features=MappingProxyType(rendered_features),
        match_candidates=tuple(
            _SvgCatalogCandidate(attributes) for part in reads for attributes in part.match_attributes
        ),
        annotation_candidates=tuple(
            _SvgCatalogCandidate(attributes) for part in reads for attributes in part.annotation_attributes
        ),
        record_source_spans=MappingProxyType(record_source_spans),
    )


# One ``_RenderedFeatureEntry`` as ``_catalog_reads_text`` writes it: its field values in this order.
_ENTRY_FIELDS = tuple(field.name for field in dataclass_fields(_RenderedFeatureEntry))


def _catalog_reads_text(reads: _CatalogElementReads) -> str:
    """The reads of one tree as canonical JSON: equal text, equal reads.

    ``catalog_item`` is keyed by these texts by value (``_catalog_reads_from_text``
    reverses it exactly).
    """

    return json.dumps(
        [
            [[index, list(span)] for index, span in reads.record_source_spans],
            [[getattr(entry, name) for name in _ENTRY_FIELDS] for entry in reads.feature_entries],
            [dict(attributes) for attributes in reads.match_attributes],
            [dict(attributes) for attributes in reads.annotation_attributes],
        ],
        ensure_ascii=False,
        separators=(",", ":"),
    )


# The text of a content the catalog reads nothing from.
_NO_READS_TEXT = _catalog_reads_text(_CatalogElementReads((), (), (), ()))


def _catalog_reads_from_text(text: str) -> _CatalogElementReads:
    # Only ``_catalog_reads_text`` writes the text, in this shape.
    spans, entries, matches, annotations = cast(
        "tuple[list[tuple[int, tuple[int, int, int]]], list[list[object]], list[dict[str, str]], list[dict[str, str]]]",
        json.loads(text),
    )
    return _CatalogElementReads(
        tuple((int(index), (int(span[0]), int(span[1]), int(span[2]))) for index, span in spans),
        # By field name: ``_RenderedFeatureEntry`` alone owns the field types.
        tuple(
            _RenderedFeatureEntry(**cast("dict[str, Any]", dict(zip(_ENTRY_FIELDS, values, strict=True))))
            for values in entries
        ),
        tuple(MappingProxyType(attributes) for attributes in matches),
        tuple(MappingProxyType(attributes) for attributes in annotations),
    )


def _match_kind(attributes: Mapping[str, object]) -> str:
    value = _text(attributes.get("data-match-kind")).lower()
    if value in {"pairwise", "orthogroup", "collinear", "homology"}:
        return value
    if attributes.get("data-collinearity-block-id"):
        return "collinear"
    if attributes.get("data-orthogroup-id"):
        return "orthogroup"
    return "pairwise"


def _metadata_values(value: object | None) -> list[str]:
    return list(
        dict.fromkeys(
            part.strip() for part in _text(value).split(";") if part.strip()
        )
    )


def _match_payload(
    attributes: Mapping[str, object],
    index: int,
) -> dict[str, object]:
    match_id = _match_id(attributes) or f"match_{index + 1}"
    payload: dict[str, object] = {
        "id": match_id,
        "match_kind": _match_kind(attributes),
        "orthogroup_ids": _metadata_values(
            attributes.get("data-orthogroup-id")
        ),
        "collinearity_block_id": _text(
            attributes.get("data-collinearity-block-id")
        ),
        "fill": _first_text(attributes.get("fill"), "#94a3b8"),
        "query_record_id": _first_text(
            attributes.get("data-query-record-id"),
            attributes.get("data-query"),
        ),
        "subject_record_id": _first_text(
            attributes.get("data-subject-record-id"),
            attributes.get("data-subject"),
        ),
        "orientation": _first_text(
            attributes.get("data-collinearity-orientation"),
            attributes.get("data-orientation"),
        ),
    }
    payload.update(
        {
            field: _text(attributes.get(attribute))
            for field, attribute in _MATCH_ATTRIBUTES.items()
        }
    )
    return _compact_wire_mapping(payload)


def _record_source_span(element: ET.Element) -> tuple[int, tuple[int, int, int]] | None:
    """Read one record group's input-file span (PD-OI-076), as the Web reads it."""

    raw = (
        element.get("data-gbdraw-record-index"),
        element.get("data-gbdraw-record-source-start"),
        element.get("data-gbdraw-record-source-end"),
        element.get("data-gbdraw-record-source-step"),
    )
    if raw[3] is None or raw[0] is None:
        return None
    try:
        index, start, end, step = (int(str(value)) for value in raw)
    except ValueError:
        return None
    if index < 0 or start < 1 or end < start or step not in (1, -1):
        return None
    return index, (start, end, step)


def _record_source_interval_fields(
    payload: Mapping[str, object],
    record_source_spans: Mapping[int, tuple[int, int, int]],
) -> dict[str, str]:
    """Source coordinates of a match on cropped or reversed records (PD-OI-076).

    Mirrors ``recordSourceInterval`` in ``services/record-source-coordinates.js``,
    which fills the same fields for Interactive SVGs built without a catalog.
    """

    if _text(payload.get("match_kind")) == "homology":
        return {}
    fields: dict[str, str] = {}
    for role, prefix in (("query", "q"), ("subject", "s")):
        try:
            span = record_source_spans.get(int(str(payload.get(f"{role}_record_index"))))
            local = [float(str(payload.get(f"{prefix}{name}"))) for name in ("start", "end")]
        except (TypeError, ValueError):
            continue
        if span is None or not all(value.is_integer() and value >= 1 for value in local):
            continue
        low, high, step = span
        start, end = (low + int(value) - 1 if step == 1 else high - int(value) + 1 for value in local)
        fields[f"{prefix}source_start"] = str(start)
        fields[f"{prefix}source_end"] = str(end)
        if low != 1:
            fields[f"{prefix}table_interval"] = f"{start - low + 1}..{end - low + 1}"
    return fields


def _match_payloads(
    candidates: Sequence[_SvgCatalogCandidate],
    rendered_features: Mapping[str, Mapping[str, object]],
    record_source_spans: Mapping[int, tuple[int, int, int]],
) -> list[dict[str, object]]:
    payloads: list[dict[str, object]] = []
    grouped: dict[str, list[_SvgCatalogCandidate]] = {}
    for candidate in candidates:
        payload = _match_payload(candidate.attributes, len(grouped))
        match_id = _text(payload.get("id"))
        if not match_id:
            raise GbdrawError("Rendered SVG contains invalid or duplicate match IDs.")
        grouped.setdefault(match_id, []).append(candidate)
    for group in grouped.values():
        _validate_match_fragments([candidate.attributes for candidate in group])
        payload = _match_payload(group[0].attributes, len(payloads))
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
        payload.update(_record_source_interval_fields(payload, record_source_spans))
        payloads.append(_compact_wire_mapping(payload))
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
    resolved = ""
    for key in keys:
        if key not in payload:
            continue
        value = _text(payload.get(key))
        if not value:
            continue
        if resolved and value != resolved:
            return "", False
        resolved = value
    return resolved, True


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
    if value is None:
        return ""
    value_type = type(value)
    values: Iterable[object]
    if value_type is list or value_type is tuple:
        values = cast("Iterable[object]", value)
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        values = value
    else:
        values = (value,)
    for entry in values:
        raw = str(entry)
        if raw and raw.strip():
            return raw
    return ""


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
    return stable_id or ""


def _feature_index(feature: Mapping[str, object]) -> int | None:
    feature_index, valid = _feature_source_index_status(feature)
    if not valid:
        raise GbdrawError(
            "Feature metadata contains conflicting or invalid source feature-index aliases."
        )
    return feature_index


def _record_key_alias(feature: Mapping[str, object]) -> tuple[str, bool]:
    return _text_alias_status(feature, ("record_key", "recordKey"))


def _biological_feature_id_alias(
    feature: Mapping[str, object],
) -> tuple[str, bool]:
    return _text_alias_status(
        feature,
        ("biological_feature_id", "biologicalFeatureId"),
    )


def _validated_source_identity(
    feature: Mapping[str, object],
    record_keys: Sequence[str],
    *,
    description: str,
) -> tuple[int, str, int | None]:
    if type(feature) is dict:
        record_index = feature.get("record_idx")
        feature_index = feature.get("feature_index")
        stable_feature_id = feature.get("stable_feature_id")
        record_key = feature.get("recordKey")
        if (
            type(record_index) is int
            and 0 <= record_index < len(record_keys)
            and (
                feature_index is None
                or type(feature_index) is int
                and feature_index >= 0
            )
            and type(stable_feature_id) is str
            and stable_feature_id
            and feature.get("stable_svg_id") == stable_feature_id
            and feature.get("record_key") is None
            and (
                record_key is None
                or type(record_key) is str
                and record_key == record_keys[record_index]
            )
            and feature.keys().isdisjoint(_SOURCE_IDENTITY_CONFLICT_KEYS)
        ):
            stable_id = _stable_feature_id(feature)
            if stable_id:
                return record_index, stable_id, feature_index
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


def _record_identity_keys(
    context: InteractiveSvgContext,
    record_keys: Sequence[str],
) -> list[str]:
    """The record each drawn unit draws: the unit's own key, or for a window the key of its record (D-01)."""

    keys = [_text(value) for value in context.record_identity_keys]
    if not keys:
        return list(record_keys)
    if len(keys) != len(record_keys) or any(not key for key in keys):
        raise GbdrawError("Feature catalog record identity keys must name every drawn record.")
    return keys


def _windows_of_records(
    record_keys: Sequence[str],
    identity_keys: Sequence[str],
) -> dict[str, str]:
    """The catalog's ``windowOf``: each drawn unit that is a window of a record, by unit key, to that record's key."""

    return {
        unit_key: identity_key
        for unit_key, identity_key in zip(record_keys, identity_keys, strict=True)
        if unit_key != identity_key
    }


def _record_key(record_keys: Sequence[str], record_index: int | None) -> str:
    if record_index is None or not 0 <= record_index < len(record_keys):
        raise GbdrawError(
            "Feature catalog source identity has an invalid record index."
        )
    return str(record_keys[record_index])


def _source_feature_candidates(
    context: InteractiveSvgContext,
    record_keys: Sequence[str],
    identity_keys: Sequence[str],
) -> tuple[list[_SourceFeatureCandidate], dict[int, list[int]]]:
    """The biological features of the drawn records, and the drawn units that show each one.

    A feature that two windows of one record show is one biological feature of that record (OV-2717): the first
    unit's row stands for it, and every unit that shows it resolves to that row.
    """

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
    # A canceled Worker request stops between features: this build's state is local (D7).
    from gbdraw.api.prepared import cancel_checkpoint

    candidates: list[_SourceFeatureCandidate] = []
    first_by_identity: dict[tuple[str, str, int | None], int] = {}
    units_by_candidate: dict[int, list[int]] = {}
    for feature in source_features:
        cancel_checkpoint()
        if not isinstance(feature, Mapping):
            continue
        if _text(feature.get("type")).lower() == "source" and (
            _record_index(feature),
            _stable_feature_id(feature),
        ) not in rendered_source_identities:
            continue
        record_index, stable_id, feature_index = _validated_source_identity(
            feature,
            record_keys,
            description="Feature catalog source feature",
        )
        _record_key(record_keys, record_index)
        record_key = identity_keys[record_index]
        identity_base = stable_id or f"feature-{feature_index}"
        identity = (record_key, identity_base, feature_index)
        first = first_by_identity.setdefault(identity, len(candidates))
        if first < len(candidates) and candidates[first].record_index != record_index:
            units_by_candidate[first].append(record_index)
            continue
        units_by_candidate[len(candidates)] = [record_index]
        candidate = _SourceFeatureCandidate(
            feature=feature,
            record_index=record_index,
            stable_id=stable_id,
            feature_index=feature_index,
            record_key=record_key,
            identity_base=identity_base,
        )
        candidates.append(candidate)
    return candidates, units_by_candidate


def _normalized_biological_features(
    context: InteractiveSvgContext,
    record_keys: Sequence[str],
    identity_keys: Sequence[str],
) -> tuple[
    list[dict[str, object]],
    list[_IndexedBiologicalFeature],
]:
    # A canceled Worker request stops between features: this build's state is local (D7).
    from gbdraw.api.prepared import cancel_checkpoint

    candidates, units_by_candidate = _source_feature_candidates(context, record_keys, identity_keys)
    biological_ids = disambiguate_feature_ids(
        ((item.record_key, item.identity_base, item.feature_index) for item in candidates),
        description="Feature catalog",
    )
    normalized: list[dict[str, object]] = []
    indexed: list[_IndexedBiologicalFeature] = []

    for position, (candidate, biological_feature_id) in enumerate(zip(candidates, biological_ids, strict=True)):
        cancel_checkpoint()
        feature = candidate.feature
        record_key = candidate.record_key
        stable_id = candidate.stable_id
        feature_index = candidate.feature_index
        payload = {
            key: value
            for key, value in feature.items()
            if key not in _BIOLOGICAL_REMOVED_KEYS
        }
        payload["biologicalFeatureId"] = biological_feature_id
        payload["recordKey"] = record_key
        if feature_index is not None:
            payload["sourceFeatureIndex"] = feature_index
        if stable_id and stable_id != biological_feature_id:
            payload["stableFeatureId"] = stable_id

        qualifiers = payload.get("qualifiers")
        if isinstance(qualifiers, Mapping):
            normalized_qualifiers = {}
            for key, value in qualifiers.items():
                value_type = type(value)
                normalized_qualifiers[str(key)] = (
                    list(value)
                    if value_type is list
                    or value_type is tuple
                    or (
                        isinstance(value, Sequence)
                        and not isinstance(value, (str, bytes))
                    )
                    else value
                )
            payload["qualifiers"] = normalized_qualifiers
            for field in (
                "protein_id",
                "locus_tag",
                "gene_id",
                "old_locus_tag",
                "gene",
                "product",
            ):
                field_value = payload.get(field)
                if (
                    field_value is None
                    or field_value == ""
                    or field_value is False
                    or field_value
                    == _first_qualifier_value(
                        normalized_qualifiers.get(field)
                    )
                ):
                    payload.pop(field, None)
            note = payload.get("note")
            if (
                note is None
                or note == ""
                or note is False
                or note
                == _first_qualifier_value(
                    normalized_qualifiers.get("note")
                )[:50]
            ):
                payload.pop("note", None)

        source_protein_id = _text(
            payload.get("source_protein_id")
            or payload.get("sourceProteinId")
        )
        if source_protein_id:
            protein_id = _text(
                payload.get("protein_id")
                or payload.get("proteinId")
                or (
                    payload_qualifiers.get("protein_id")
                    if isinstance(
                        payload_qualifiers := payload.get("qualifiers"), Mapping
                    )
                    else None
                )
            )
            if protein_id and source_protein_id == protein_id:
                payload.pop("source_protein_id", None)
                payload.pop("sourceProteinId", None)

        if type(payload.get("qualifiers")) is dict:
            compact_qualifiers = _compact_wire_value(payload["qualifiers"])
            if compact_qualifiers is None:
                payload.pop("qualifiers")
            else:
                payload["qualifiers"] = compact_qualifiers

        raw_parts = payload.get("location_parts") or payload.get(
            "locationParts"
        )
        parts = _sequence(raw_parts) if raw_parts else []
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

        compact: dict[str, object] = {}
        for key, value in payload.items():
            value_type = type(value)
            if key == "qualifiers" and value_type is dict:
                normalized_value = value or None
            elif value_type is str or value_type is bool:
                normalized_value = value or None
            elif value is None:
                normalized_value = None
            elif value_type is int or value_type is float or value_type is bytes:
                normalized_value = value
            else:
                normalized_value = _compact_wire_value(value)
            if normalized_value is not None:
                compact[str(key)] = normalized_value
        normalized.append(compact)
        indexed.append(
            _IndexedBiologicalFeature(
                normalized=compact,
                record_indexes=tuple(units_by_candidate[position]),
                stable_id=stable_id,
                feature_index=feature_index,
            )
        )
    return normalized, indexed


_IndexKey = TypeVar("_IndexKey")


class _BiologicalFeatureIndex:
    def __init__(
        self,
        indexed_features: Sequence[_IndexedBiologicalFeature],
    ) -> None:
        self._by_record_feature_index: dict[
            tuple[int, int], Mapping[str, object] | None
        ] = {}
        self._by_record_stable: dict[
            tuple[int, str],
            Mapping[str, object] | list[Mapping[str, object]],
        ] = {}
        self._by_canonical: dict[
            tuple[str, str], Mapping[str, object] | None
        ] = {}
        for indexed in indexed_features:
            normalized = indexed.normalized
            stable_id = indexed.stable_id
            feature_index = indexed.feature_index
            for record_index in indexed.record_indexes:
                if feature_index is not None:
                    self._add_unique(
                        self._by_record_feature_index,
                        (record_index, feature_index),
                        normalized,
                    )
                if stable_id:
                    self._add_stable(
                        self._by_record_stable,
                        (record_index, stable_id),
                        normalized,
                    )
            canonical = _feature_reference(normalized)
            if all(canonical):
                self._add_unique(self._by_canonical, canonical, normalized)

    @staticmethod
    def _add_unique(
        index: dict[_IndexKey, Mapping[str, object] | None],
        key: _IndexKey,
        feature: Mapping[str, object],
    ) -> None:
        index[key] = None if key in index else feature

    @staticmethod
    def _add_stable(
        index: dict[
            tuple[int, str],
            Mapping[str, object] | list[Mapping[str, object]],
        ],
        key: tuple[int, str],
        feature: Mapping[str, object],
    ) -> None:
        existing = index.get(key)
        if existing is None:
            index[key] = feature
        elif isinstance(existing, list):
            existing.append(feature)
        else:
            index[key] = [existing, feature]

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
            stable_candidate = self._by_record_stable.get(
                (record_index, stable_id)
            )
            if feature_index is None:
                return (
                    stable_candidate
                    if isinstance(stable_candidate, Mapping)
                    else None
                )
            indexed = self._by_record_feature_index.get(
                (record_index, feature_index)
            )
            if indexed is None:
                return None
            if isinstance(stable_candidate, list):
                return (
                    indexed
                    if any(
                        candidate is indexed
                        for candidate in stable_candidate
                    )
                    else None
                )
            return indexed if indexed is stable_candidate else None
        if record_index is not None and feature_index is not None:
            return self._by_record_feature_index.get(
                (record_index, feature_index)
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
        return self._by_canonical.get(
            (record_key, biological_feature_id)
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
    index: _RenderedSvgCatalogIndex,
    context: InteractiveSvgContext,
    record_keys: Sequence[str],
    biological_index: _BiologicalFeatureIndex,
) -> tuple[
    list[dict[str, object]],
    dict[str, Mapping[str, object]],
    dict[str, dict[str, object]],
]:
    normalized: list[dict[str, object]] = []
    references_by_svg_id: dict[str, Mapping[str, object]] = {}
    match_identities_by_svg_id: dict[str, dict[str, object]] = {}
    for resolved in _resolve_rendered_features(
        index.rendered_features,
        context.features,
    ):
        svg_id = resolved.svg_id
        stable_id = _first_text(
            (
                resolved.stable_feature_id
                if resolved.feature is not None
                else ""
            ),
            resolved.entry.payload_stable_id,
        )
        record_index = resolved.record_index
        feature_index = resolved.feature_index
        description = f"Rendered feature {svg_id!r}"
        if (
            record_index is None
            or not 0 <= record_index < len(record_keys)
            or (not stable_id and feature_index is None)
        ):
            raise GbdrawError(
                f"{description} must identify one biological feature by a "
                "record-scoped stable feature ID or source feature index and a "
                "valid record index/key."
            )
        reference = biological_index.resolve(
            record_index=record_index,
            stable_id=stable_id,
            feature_index=feature_index,
        )
        if reference is None:
            raise GbdrawError(
                f"{description} identity fields do not resolve to the same "
                "biological feature by a record-scoped stable feature ID or "
                "source feature index."
            )
        record_key, biological_feature_id = _feature_reference(reference)
        if not record_key or not biological_feature_id:
            raise GbdrawError(
                "Rendered feature metadata does not resolve to one biological "
                f"feature: {svg_id}."
            )
        if svg_id in references_by_svg_id:
            raise GbdrawError(
                f"Rendered feature SVG ID {svg_id!r} is duplicated."
            )
        row: dict[str, object] = {"svgId": svg_id, "recordKey": record_key}
        if record_key != record_keys[record_index]:
            # A window of a record: the row names its record and the drawn unit it lies in (OV-2717).
            row["recordIndex"] = record_index
        row.update(
            {
                "biologicalFeatureId": biological_feature_id,
                "fillColor": resolved.entry.fill,
                "drawnSelector": _drawn_selector(resolved.feature),
            }
        )
        normalized.append(row)
        references_by_svg_id[svg_id] = reference
        match_identity: dict[str, object] = {
            "svg_id": svg_id,
            "record_idx": record_index,
            "stable_svg_id": stable_id,
        }
        if feature_index is not None:
            match_identity["feature_index"] = feature_index
        match_identities_by_svg_id[svg_id] = match_identity
    return normalized, references_by_svg_id, match_identities_by_svg_id


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
        # The drawn unit of each member: a feature that two windows of its record show is one member (OV-2717).
        known_member_references: dict[tuple[str, str], int | None] = {}
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
            member_unit = _record_index(member)
            if canonical_reference in known_member_references:
                if known_member_references[canonical_reference] in (None, member_unit):
                    raise GbdrawError(
                        f"Orthogroup {group_id!r} contains duplicate biological members."
                    )
                continue
            known_member_references[canonical_reference] = member_unit
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
            members.append(_compact_wire_mapping(compact_member))
        compact_group = _compact_wire_mapping(payload)
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
    indexes = [index for index in map(parse_ascii_nonnegative_integer, values) if index is not None]
    if len(indexes) != len(values):
        raise GbdrawError(
            f"{description} contains an invalid source feature index."
        )
    return indexes


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
        normalized.append(_compact_wire_mapping(payload))
    return normalized


def _normalized_annotations(
    candidates: Sequence[_SvgCatalogCandidate],
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
    for candidate in candidates:
        attributes = candidate.attributes
        annotation_id = _text(attributes.get(ANNOTATION_ID_ATTRIBUTE))
        key = (
            _text(attributes.get("data-gbdraw-record-index")),
            _text(attributes.get("data-gbdraw-annotation-track-id")),
            _text(attributes.get(ANNOTATION_SET_ID_ATTRIBUTE)),
            annotation_id,
        )
        payload = context_by_key.get(
            key,
            context_by_key.get((key[0], "", key[2], key[3]), {}),
        ).copy()
        record_index = parse_ascii_nonnegative_integer(key[0])
        payload.update(
            {
                "dom_id": _text(attributes.get("id")),
                "id": annotation_id,
                "set_id": key[2],
                "track_id": key[1],
                "record_id": _text(
                    attributes.get("data-gbdraw-record-id")
                ),
                "record_index": key[0] if record_index is None else record_index,
                "mark": _text(
                    attributes.get("data-gbdraw-annotation-mark")
                    or payload.get("mark")
                ),
                "label": _text(
                    attributes.get("data-gbdraw-annotation-label")
                    or payload.get("label")
                ),
            }
        )
        annotations.append(_compact_wire_mapping(payload))
    return annotations


def _deduplicated_sequence_sources(
    matches: Sequence[Mapping[str, object]],
    sources: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    selected = _sequence_sources_for_matches(matches, sources)
    normalized: list[dict[str, object]] = []
    by_key: dict[str, dict[str, object]] = {}
    for source in selected:
        if not isinstance(source, Mapping):
            continue
        payload = _compact_wire_mapping(dict(source))
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


def _integer_or_none(value: Any) -> int | None:
    try:
        number = ascii_float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number) or not number.is_integer():
        return None
    return int(number)


class _ComplementTable(dict[int, int]):
    """Translate table that maps every unlisted character to ``N``."""

    def __missing__(self, ordinal: int) -> int:
        return ord("N")


_COMPLEMENT_TABLE = _ComplementTable(
    str.maketrans("ACGTURYSWKMBDHVN-", "TGCAAYRSWMKVHDBN-")
)


def _reverse_complement_sequence(sequence: str) -> str:
    return "".join(sequence.split()).upper()[::-1].translate(_COMPLEMENT_TABLE)


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
            reverse = strand == "-" or ascii_float(strand) < 0
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
    """Reference a record source only after exact sequence reconstruction.

    ``record_keys`` names the record of each drawn unit; a record drawn in windows is read from its first unit, whose
    row stands for the feature (``_normalized_biological_features``).
    """

    record_index_by_key: dict[str, int] = {}
    for unit_index, record_key in enumerate(record_keys):
        record_index_by_key.setdefault(str(record_key), unit_index)
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


def _source_record_sequence(record: SeqRecord) -> str:
    """A match source record's bases. A record with none (a CONTIG line) stops with the input's reason
    (``SEQUENCE_UNDEFINED``), not as an invalid catalog."""

    from gbdraw.io.genome import record_sequence_text, require_record_sequence

    annotations = getattr(record, "annotations", None) or {}
    require_record_sequence(record, str(annotations.get("gbdraw_source_file") or "the record input"))
    return record_sequence_text(record) or ""


def match_sequence_sources(records: CatalogSequenceRecords) -> list[dict[str, object]]:
    """The whole-record sequences that match spans are read from (D5): each drawn unit (``linear:record:N`` or
    ``circular:record:N``, N its place in ``records``) and every comparison FASTA record
    (``homology:comparison:S:ID``). The one producer of a record's bases for a match: the Worker reads spans from
    it, and an exported Interactive SVG embeds the sources its matches read."""

    from gbdraw.api.prepared import cancel_checkpoint

    sources: list[dict[str, object]] = []
    for record_index, record in enumerate(records.records):
        cancel_checkpoint()
        aliases = {str(record.id or "").strip(), str(record.name or "").strip()}
        accessions = (getattr(record, "annotations", {}) or {}).get("accessions", [])
        if isinstance(accessions, str):
            aliases.add(accessions.strip())
        else:
            aliases.update(str(value or "").strip() for value in accessions)
        sources.append({
            "key": f"{records.mode}:record:{record_index}",
            "recordId": str(record.id),
            "aliases": sorted(value for value in aliases if value),
            "sequence": _source_record_sequence(record),
            "origin": "linear-record" if records.mode == "linear" else "circular-reference",
            "recordIndex": record_index,
        })
    for source_index, source_records in enumerate(records.comparison_records):
        for record in source_records:
            sources.append({
                "key": f"homology:comparison:{source_index}:{record.id}",
                "recordId": str(record.id),
                "aliases": [str(record.id), str(record.name)],
                "sequence": _source_record_sequence(record),
                "origin": "homology-comparison",
                "sourceIndex": source_index,
            })
    return sources


def read_match_span(sequence: str, start: int, end: int) -> str | None:
    """The bases of a match span, 1-based and inclusive on the displayed record; a start after the end reads the
    reverse strand (reverse complement, IUPAC letters kept). None when the span leaves the record."""

    length = len(sequence)
    if not (1 <= start <= length and 1 <= end <= length):
        return None
    low, high = min(start, end), max(start, end)
    bases = sequence[low - 1:high]
    return _reverse_complement_sequence(bases) if start > end else bases


SpanCoordinate = int | None
"""A span coordinate as the drawing gives it: a whole number, or ``None`` when the drawing's value is not a whole
number (the page reads the value as the match popup did before D5, and Python says so after the record)."""


@final
class RecordMatchSpan(TypedDict):
    """A span of a drawn record, by the record key of its drawn unit, with the record ID or alias the match names
    ('' when it names none), which the record must carry: displayed coordinates, 1-based and inclusive; a start after
    the end is the reverse strand."""

    kind: Literal["record"]
    recordKey: str
    recordId: str
    start: SpanCoordinate
    end: SpanCoordinate


@final
class ReferenceMatchSpan(TypedDict):
    """A span of a drawn record named by its record ID or an alias (a Circular comparison band names its reference
    by the ID its hit carries), searched among every drawn unit, or only the unit of ``recordKey`` (a band on a
    cropped or reversed reference also carries its record index); coordinates as in ``RecordMatchSpan``."""

    kind: Literal["reference"]
    recordId: str
    recordKey: str | None
    start: SpanCoordinate
    end: SpanCoordinate


@final
class ComparisonMatchSpan(TypedDict):
    """A span of a comparison FASTA record, by the index of its comparison source (None: the match names no source,
    and every source is searched) and its record ID or alias."""

    kind: Literal["comparison"]
    sourceIndex: int | None
    recordId: str
    start: SpanCoordinate
    end: SpanCoordinate


MatchSpanInput = RecordMatchSpan | ReferenceMatchSpan | ComparisonMatchSpan
MatchSpanReason = Literal[
    "NEVER_SUPPLIED", "SOURCE_UNAVAILABLE", "SOURCE_KEY_CONFLICT", "RECORD_UNKNOWN", "RECORD_AMBIGUOUS",
    "RECORD_ALIAS_AMBIGUOUS", "COORDINATE_NOT_WHOLE", "COORDINATE_BELOW_ONE", "OUT_OF_RANGE",
]
# Why a span cannot be read: the Web reads this list from here (tools/generate_web_contracts.py) and shows one text
# per reason, the text the match popup showed for that case before D5.
MATCH_SPAN_REASONS: tuple[MatchSpanReason, ...] = get_args(MatchSpanReason)
# Why a match names no span to read at all, known from the match alone: an index that is not a whole number, or no
# record named. The Web's match popup says so when it opens; an exported Interactive SVG writes it for that span.
MatchIdentityReason = Literal["IDENTITY_INVALID", "RECORD_UNNAMED"]
EmbeddedMatchSpanReason = MatchSpanReason | MatchIdentityReason
# The sentence of each reason, in English: the match popup of the Web shows it translated (a catalog key) and an
# exported Interactive SVG shows it as it is ({id}: the record the match names; {length}: the length of the record
# the span reads). The Web reads the table from here (tools/generate_web_contracts.py writes it into
# drawn-attributes.generated.js, not the persisted contracts, so its wording is not a contract); the export writes
# it into the file's script (``STANDALONE_MATCH_SPAN_REASON_TEXTS_SLOT``).
MATCH_SPAN_REASON_TEXTS: dict[EmbeddedMatchSpanReason, str] = {
    "NEVER_SUPPLIED": "Comparison sequence was not supplied for this BLAST source.",
    "SOURCE_UNAVAILABLE": "The sequence source for this match is unavailable.",
    "SOURCE_KEY_CONFLICT": "The sequence source key conflicts with the match identity.",
    "RECORD_UNKNOWN": "No sequence source matched record \u201c{id}\u201d.",
    "RECORD_AMBIGUOUS": "Record ID \u201c{id}\u201d is ambiguous in the available sequence sources.",
    "RECORD_ALIAS_AMBIGUOUS": "Record alias \u201c{id}\u201d is ambiguous in the available sequence sources.",
    "COORDINATE_NOT_WHOLE": "Coordinates must be whole 1-based values.",
    "COORDINATE_BELOW_ONE": "Coordinates must be greater than or equal to 1.",
    "OUT_OF_RANGE": "Coordinates exceed the sequence length ({length} bp).",
    "IDENTITY_INVALID": "The sequence identifier metadata for this match is invalid.",
    "RECORD_UNNAMED": "The match does not identify a sequence record.",
}



def _aliases(source: Mapping[str, object]) -> set[str]:
    """A source's record ID and aliases, each also by its first word (the match registry's rule before D5)."""

    names = {str(source["recordId"]).strip(), *(str(alias).strip() for alias in cast("list[object]", source["aliases"]))}
    return {name for value in names if value for name in (value, value.split()[0])}


def registered_match_sources(sources: Sequence[Mapping[str, object]]) -> list[Mapping[str, object]]:
    """The sources the match registry held before D5: none with an empty sequence, and of the sources that share a
    key (two comparison records with one ID), the one when they agree, else none."""

    by_key: dict[object, list[Mapping[str, object]]] = {}
    for source in sources:
        if source["sequence"]:
            by_key.setdefault(source["key"], []).append(source)
    return [
        group[0] for group in by_key.values()
        if all(
            (other["recordId"], other["sequence"], _aliases(other)) == (group[0]["recordId"], group[0]["sequence"], _aliases(group[0]))
            for other in group[1:]
        )
    ]


def match_span_source(
    records: CatalogSequenceRecords,
    sources: Sequence[Mapping[str, object]],
    span: MatchSpanInput,
) -> tuple[Mapping[str, object] | None, MatchSpanReason | None]:
    """The source a span reads among the ``registered_match_sources``, found in the match popup's order before D5: a drawn
    unit by its record key, which must carry the record ID or alias the match names; a drawn record by its record ID,
    else an alias (the record identity ``match_sequence_sources`` gives: ID, name, accessions), among every unit or
    the one a reference span names, which it reads when the span names no ID; or a comparison record of a source (of
    every source when the match names none) by its record ID, else an alias. Only a comparison span whose source has
    no comparison FASTA (no source has one, when the match names none) says its FASTA was not supplied; one whose
    FASTA lacks the record gives the record's reason, as a drawn record does (OV-2755). Ambiguity and absence are
    reasons, never a guess."""

    wanted = span["recordId"].strip()
    drawn = [source for source in sources if "sourceIndex" not in source]
    if span["kind"] == "record":
        keyed = _unit_sources(records, drawn, span["recordKey"])
        if not keyed:
            return None, "SOURCE_UNAVAILABLE"
        if wanted and wanted not in _aliases(keyed[0]):
            return None, "SOURCE_KEY_CONFLICT"
        return keyed[0], None
    if span["kind"] == "reference":
        candidates = drawn if span["recordKey"] is None else _unit_sources(records, drawn, span["recordKey"])
        if span["recordKey"] is not None and not wanted and candidates:
            return candidates[0], None
    else:
        index = span["sourceIndex"]
        supplied = records.comparison_records if index is None else records.comparison_records[index:index + 1]
        if not any(supplied):
            return None, "NEVER_SUPPLIED"
        candidates = [source for source in sources if "sourceIndex" in source and index in (None, source["sourceIndex"])]
    return _by_record_id(candidates, wanted)


def _unit_sources(
    records: CatalogSequenceRecords, drawn: Sequence[Mapping[str, object]], record_key: str
) -> list[Mapping[str, object]]:
    """The source of the drawn unit of ``record_key`` (none when no unit has that key)."""

    if record_key not in records.unit_keys:
        return []
    unit = records.unit_keys.index(record_key)
    return [source for source in drawn if source.get("recordIndex") == unit]


def _by_record_id(
    candidates: Sequence[Mapping[str, object]], wanted: str
) -> tuple[Mapping[str, object] | None, MatchSpanReason | None]:
    """The one candidate whose record ID is ``wanted``, else the one with that alias."""

    if not wanted:
        return None, "RECORD_UNKNOWN"
    reasons: tuple[MatchSpanReason, MatchSpanReason] = ("RECORD_AMBIGUOUS", "RECORD_ALIAS_AMBIGUOUS")
    for found, ambiguous in zip((
        [source for source in candidates if source["recordId"] == wanted],
        [source for source in candidates if wanted in _aliases(source)],
    ), reasons, strict=True):
        if len(found) == 1:
            return found[0], None
        if found:
            return None, ambiguous
    return None, "RECORD_UNKNOWN"


def span_reason(sequence: str, start: SpanCoordinate, end: SpanCoordinate) -> MatchSpanReason | None:
    """Why the coordinates of a span cannot be read on ``sequence``, checked in the match popup's order before D5."""

    if start is None or end is None:
        return "COORDINATE_NOT_WHOLE"
    if start < 1 or end < 1:
        return "COORDINATE_BELOW_ONE"
    if start > len(sequence) or end > len(sequence):
        return "OUT_OF_RANGE"
    return None


@final
class ReadableMatchSpan(TypedDict):
    """A span an exported Interactive SVG reads: the place, among the item's ``sequenceSources``, of the record whose
    bases it reads (on the coordinates and strand of its match)."""

    source: int


@final
class UnreadableMatchSpan(TypedDict):
    """A span an exported Interactive SVG cannot read: why, and the length of the record it names (0 when no record
    is found), as the Worker's read answers (``MatchSpanSequence``)."""

    reason: EmbeddedMatchSpanReason
    sourceLength: int


EmbeddedMatchSpan = ReadableMatchSpan | UnreadableMatchSpan


@final
class EmbeddedMatchSpans(TypedDict):
    """The two spans of one match in an exported Interactive SVG (``matchSpanReads``, by match ID)."""

    query: EmbeddedMatchSpan
    subject: EmbeddedMatchSpan


def _span_coordinate_of(value: object) -> SpanCoordinate:
    """A match's coordinate as the export reads it: a whole number, or ``None`` (missing or not whole)."""

    text = value if type(value) is str else str(value)
    # The common case, ASCII digits that a float holds exactly, needs no float.
    if text.isascii() and text.isdigit() and len(text) < 16:
        return int(text)
    try:
        number = float(text)
    except (TypeError, ValueError):
        return None
    return int(number) if math.isfinite(number) and number.is_integer() else None


def _match_index(value: object) -> tuple[bool, int | None]:
    """An index a match names (record or source) as the Web's match popup reads it: absent (valid, None), a whole
    number a JavaScript number holds exactly, or invalid."""

    text = _text(value)
    if not text:
        return True, None
    index = parse_ascii_nonnegative_integer(text)
    return (False, None) if index is None or index > JS_MAX_SAFE_INTEGER else (True, index)


# The name of the record a span reads, without its coordinates: (kind, record key, record ID, source index).
_SpanName = tuple[str, str | None, str, int | None]


def _embedded_match_span(
    match: Mapping[str, object], role: Literal["query", "subject"], unit_keys: Sequence[str]
) -> tuple[_SpanName, SpanCoordinate, SpanCoordinate] | MatchSpanReason | MatchIdentityReason:
    """The span of a match's ``role`` that the export reads, named as the Web's match popup names it
    (``matchSpanRead``): a comparison record of a ring's source for the side that is not the reference, a drawn unit
    by its record index for a Linear match, else a drawn record by ID or alias (in the unit of a record index, if
    any). A reason when the match names no record that can be read: an index that is not a whole number (for either
    span, as the popup reads it), a record index past the drawn units, or no record named."""

    record_valid, record_index = _match_index(match.get(f"{role}_record_index"))
    source_valid, source_index = _match_index(match.get("source_index"))
    if not (record_valid and source_valid):
        return "IDENTITY_INVALID"
    start = _span_coordinate_of(match.get(f"{role[0]}start"))
    end = _span_coordinate_of(match.get(f"{role[0]}end"))
    record_id = _text(match.get(f"{role}_record_id"))
    homology = _text(match.get("match_kind")) == "homology"
    if homology and role != _text(match.get("reference_side")).lower():
        return ("comparison", None, record_id, source_index), start, end
    record_key = unit_keys[record_index] if record_index is not None and record_index < len(unit_keys) else None
    if not homology and record_index is not None:
        if record_key is None:
            return "SOURCE_UNAVAILABLE"
        return ("record", record_key, record_id, None), start, end
    if not record_id and record_key is None:
        return "RECORD_UNNAMED"
    if record_index is not None and record_key is None:
        return "RECORD_UNKNOWN"
    return ("reference", record_key, record_id, None), start, end


def _span_input(name: _SpanName, start: SpanCoordinate, end: SpanCoordinate) -> MatchSpanInput:
    """The ``MatchSpanInput`` of a span name, as the Worker's read takes it."""

    kind, record_key, record_id, source_index = name
    if kind == "comparison":
        return {"kind": "comparison", "sourceIndex": source_index, "recordId": record_id, "start": start, "end": end}
    if kind == "record":
        return {"kind": "record", "recordKey": cast("str", record_key), "recordId": record_id, "start": start, "end": end}
    return {"kind": "reference", "recordId": record_id, "recordKey": record_key, "start": start, "end": end}


def embedded_match_span_reads(
    matches: Sequence[Mapping[str, object]],
    records: CatalogSequenceRecords,
    sources: Sequence[Mapping[str, object]],
    embedded: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], dict[str, EmbeddedMatchSpans]]:
    """What an exported Interactive SVG reads for the spans of ``matches`` (D5, OV-3542): each span's record found
    by ``match_span_source`` among the ``registered_match_sources`` of ``sources`` and judged by ``span_reason``, as
    the Worker reads a span for the Web's match popup. ``embedded`` are the sources the file embeds as the export
    wrote them before (``_deduplicated_sequence_sources``; the features of those records read their bases from
    them, ``sequenceSourceIndex``). Returns the sources a readable span reads that ``embedded`` lacks, which the file
    appends after them, and each match's two spans by match ID. A similarity group ribbon names no span. Each record
    a span names is searched once, however many matches name it."""

    registered = registered_match_sources(sources)
    found: dict[_SpanName, tuple[Mapping[str, object] | None, MatchSpanReason | None]] = {}
    Resolved = tuple[Mapping[str, object] | None, EmbeddedMatchSpanReason | None, int]

    def resolved(
        span: tuple[_SpanName, SpanCoordinate, SpanCoordinate] | MatchSpanReason | MatchIdentityReason
    ) -> Resolved:
        """The source a span reads (None when it cannot be read), why not, and the length of the record found."""

        if isinstance(span, str):
            return None, span, 0
        name, start, end = span
        if name not in found:
            found[name] = match_span_source(records, registered, _span_input(name, start, end))
        source, reason = found[name]
        sequence = "" if source is None else cast("str", source["sequence"])
        reason = reason or span_reason(sequence, start, end)
        return (None, reason or "SOURCE_UNAVAILABLE", len(sequence)) if source is None or reason else (source, None, 0)

    spans: dict[str, tuple[Resolved, Resolved]] = {
        _text(match.get("id")): (
            resolved(_embedded_match_span(match, "query", records.unit_keys)),
            resolved(_embedded_match_span(match, "subject", records.unit_keys)),
        )
        for match in matches
        if _text(match.get("match_kind")) != "orthogroup"
    }
    places = {source["key"]: index for index, source in enumerate(embedded)}
    read_keys = {source["key"] for pair in spans.values() for source, _reason, _length in pair if source is not None}
    appended = [
        _compact_wire_mapping(dict(source))
        for source in registered
        if source["key"] in read_keys and source["key"] not in places
    ]
    places.update((source["key"], len(embedded) + index) for index, source in enumerate(appended))

    def embedded_span(read: Resolved) -> EmbeddedMatchSpan:
        source, reason, length = read
        if source is not None:
            return {"source": places[source["key"]]}
        return {"reason": cast("EmbeddedMatchSpanReason", reason), "sourceLength": length}

    return appended, {
        match_id: {"query": embedded_span(query), "subject": embedded_span(subject)}
        for match_id, (query, subject) in spans.items()
    }


# The export profile of the catalog (D5-02): an exported Interactive SVG embeds the catalog with its sequences, in
# this schema and its compact encoding; the stored catalog (``FEATURE_CATALOG_SCHEMA``) holds no sequence.
INTERACTIVE_EMBED_SCHEMA = SEQUENCE_FEATURE_CATALOG_SCHEMA
# The keys a biological row reads from its feature before its sequences, in order (``feature_metadata``): the
# embedded sequences follow the last of them, where the catalog wrote them before D5.
_KEYS_BEFORE_SEQUENCES = (
    "record_id", "record_idx", "feature_index", "organism", "type", "start", "end", "strand", "protein_id",
    "source_protein_id", "locus_tag", "gene_id", "old_locus_tag", "gene", "product", "note", "qualifiers",
    "selector", "location_parts", "anchorProfile",
)


def _row_with_sequences(row: Mapping[str, object], sequences: FeatureSequences | None) -> dict[str, object]:
    """A stored biological row with its sequences in the export encoding: the three keys after the keys read from
    the feature, and, when its ``/translation`` qualifier is its amino-acid sequence, the qualifier dropped for the
    ``translationFromAminoAcidSequence`` mark."""

    if sequences is None:
        return dict(row)
    fields = {
        key: value
        for key, value in (
            ("nucleotide_sequence", sequences["nucleotideSequence"] or None),
            ("amino_acid_sequence", sequences["aminoAcidSequence"] or None),
            ("sequence_warnings", _compact_wire_value(list(sequences["sequenceWarnings"]))),
        )
        if value is not None
    }
    last = max((index for index, key in enumerate(row) if key in _KEYS_BEFORE_SEQUENCES), default=-1)
    embedded: dict[str, object] = {}
    for index, (key, value) in enumerate(row.items()):
        if index == 0 and last < 0:
            embedded.update(fields)
        embedded[key] = value
        if index == last:
            embedded.update(fields)
    amino_acid_sequence = sequences["aminoAcidSequence"]
    qualifiers = embedded.get("qualifiers")
    translation = qualifiers.get("translation") if isinstance(qualifiers, Mapping) else None
    if (
        amino_acid_sequence
        and isinstance(qualifiers, Mapping)
        and isinstance(translation, list)
        and len(translation) == 1
        and translation[0] == amino_acid_sequence
    ):
        kept = {key: value for key, value in qualifiers.items() if key != "translation"}
        if kept:
            embedded["qualifiers"] = kept
        else:
            embedded.pop("qualifiers")
        embedded["translationFromAminoAcidSequence"] = True
    return embedded


def embed_json_text(value: object) -> str:
    """The text of the export encoding (D5-02): compact JSON with characters as they are. The CLI export writes its
    metadata with it, and the Web export writes the parts Python returns as this text."""

    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


class EmbeddedSequenceSource(TypedDict):
    """One whole-record source of the export, with its place in the item's ``sequenceSources``."""

    index: int
    source: dict[str, object]


class EmbedCatalogPart(TypedDict):
    """A part of the export encoding of one catalog item (D5-02): the biological rows asked, with their sequences,
    and the whole-record sources the part owns, each with its place among the item's ``sourceCount`` sources; the
    part that holds the comparison sources also holds what the item's match spans read (``matchSpanReads``, None in
    the other parts and when no match names a span) and the sources only a span reads, after the others, so its
    ``sourceCount`` may be larger than the other parts'."""

    biologicalFeatures: list[dict[str, object]]
    sequenceSources: list[EmbeddedSequenceSource]
    sourceCount: int
    matchSpanReads: dict[str, EmbeddedMatchSpans] | None


def embed_catalog_part(
    item: Mapping[str, object],
    records: CatalogSequenceRecords,
    *,
    source_record_keys: Sequence[str] | None,
    comparison_sources: bool,
) -> EmbedCatalogPart:
    """The export encoder (D5-02), for a part of a stored catalog item: ``item["biologicalFeatures"]`` are the rows
    of the part; ``item`` names the drawn records (``recordKeys``, ``windowOf``) and the matches of the whole item.

    Each row gets its feature's sequences from ``feature_sequence_payload`` (the one producer), found by its record
    key and source feature index in the drawn records; a nucleotide sequence that a whole-record source of the item's
    matches reconstructs is written as ``sequenceSourceIndex``. The part returns the sources of the records the
    item's matches name that belong to the records ``source_record_keys`` name (each drawn unit of such a record;
    ``None``: every record) and, with ``comparison_sources``, those of the comparison FASTA records, what each match
    span reads (``matchSpanReads``, ``embedded_match_span_reads``, read once per item) and the sources only a span
    reads, each source with its index among the item's sources. The
    CLI and Python export read the whole item as one part (``embed_catalog_item``); the Web export reads parts of a
    bounded number of bases, so no reply of a large item holds all of it.
    """

    from gbdraw.features.source import bind_source_feature_indexes
    from gbdraw.web_support.feature_metadata import feature_sequence_payload

    record_keys = [str(key) for key in cast("Sequence[object]", item["recordKeys"])]
    if list(records.unit_keys) != record_keys:
        raise GbdrawError("Interactive SVG export records must be the catalog item's drawn records.")
    windows = item.get("windowOf", {})
    identity_keys = [str(windows.get(key, key)) if isinstance(windows, Mapping) else key for key in record_keys]
    rows = cast("Sequence[Mapping[str, object]]", item["biologicalFeatures"])

    def wanted(row: Mapping[str, object]) -> tuple[str, int] | None:
        key, source_index = _text(row.get("recordKey")), row.get("sourceFeatureIndex")
        return (key, source_index) if key in identity_keys and type(source_index) is int else None

    # Each row binds as the popup read binds its identity: the first drawn unit of its record that has the feature
    # (a record drawn in windows draws a feature in the window that holds it).
    bindings = bind_source_feature_indexes(
        records=records.records, record_keys=identity_keys,
        features=[found for found in map(wanted, rows) if found is not None],
    )

    def sequences_of(row: Mapping[str, object]) -> FeatureSequences | None:
        found = wanted(row)
        binding = None if found is None else bindings[found][0]
        if binding is None or binding.feature is None or binding.status != "present":
            return None
        return feature_sequence_payload(binding.feature, records.records[binding.record_index])

    biological = [_row_with_sequences(row, sequences_of(row)) for row in rows]
    matches = cast("Sequence[Mapping[str, object]]", item.get("comparisonMatches") or ())
    sources = match_sequence_sources(records)
    embedded = _deduplicated_sequence_sources(matches, sources)
    _reference_reconstructible_nucleotide_sequences(biological, identity_keys, embedded)
    # The spans are read once per item, by the part that holds the comparison sources: it appends the sources only
    # a span reads, after the others.
    appended, reads = embedded_match_span_reads(matches, records, sources, embedded) if comparison_sources else ([], {})
    owners = set(identity_keys if source_record_keys is None else source_record_keys)

    def owned(source: Mapping[str, object]) -> bool:
        unit = source.get("recordIndex")
        return comparison_sources if type(unit) is not int else identity_keys[unit] in owners

    owned_sources: list[EmbeddedSequenceSource] = [
        {"index": index, "source": source} for index, source in enumerate(embedded) if owned(source)
    ]
    owned_sources.extend({"index": len(embedded) + index, "source": source} for index, source in enumerate(appended))
    return {
        "biologicalFeatures": biological,
        "sequenceSources": owned_sources,
        "sourceCount": len(embedded) + len(appended),
        "matchSpanReads": reads or None,
    }


def embed_catalog_item(item: Mapping[str, object], records: CatalogSequenceRecords) -> dict[str, object]:
    """The export encoder (D5-02): a stored catalog item with every sequence an exported Interactive SVG embeds, its
    biological rows in the export encoding, the whole-record sources its matches read (``sequenceSources``) and what
    each match span reads (``matchSpanReads``), as one part (``embed_catalog_part``). The CLI and Python export write
    this, in ``INTERACTIVE_EMBED_SCHEMA``."""

    whole = embed_catalog_part(item, records, source_record_keys=None, comparison_sources=True)
    sources = [entry["source"] for entry in whole["sequenceSources"]]
    embedded: dict[str, object] = {**item, "biologicalFeatures": whole["biologicalFeatures"]}
    if sources:
        embedded["sequenceSources"] = sources
    if whole["matchSpanReads"] is not None:
        embedded["matchSpanReads"] = whole["matchSpanReads"]
    return embedded


def interactive_embed_catalog(
    item: Mapping[str, object], records: CatalogSequenceRecords | None
) -> dict[str, object]:
    """The catalog an exported Interactive SVG embeds (D5-02): the file has no Python, so its one item carries every
    sequence its popups, search and match FASTA read, from ``records`` (none when the caller holds no records)."""

    return {
        "schema": INTERACTIVE_EMBED_SCHEMA,
        "items": [dict(item) if records is None else embed_catalog_item(item, records)],
    }


# The group kinds whose elements the catalog reads: features (record and axis
# groups, with the record source span), comparison ribbons and annotations.
# The Circular axis circle (AXIS_CIRCLE) holds neither.
# A Circular conservation ring carries the hits the catalog binds to features.
# RECORD: a Multi-Record Circular record whose parts are not all placed
# contents; otherwise each part is read by its own kind, as at the top level
# (a definition, tick or GC part reads nothing, OV-2108).
_CATALOG_GROUP_KINDS = frozenset(
    {
        GroupKind.AXIS,
        GroupKind.FEATURES,
        GroupKind.COMPARISON,
        GroupKind.ANNOTATIONS,
        GroupKind.CONSERVATION,
        GroupKind.RECORD,
    }
)
# The namespaces of an svgwrite root, so the contents parse as in the document.
_CATALOG_CONTENTS_ROOT = (
    '<svg xmlns="http://www.w3.org/2000/svg" xmlns:ev="http://www.w3.org/2001/xml-events"'
    ' xmlns:xlink="http://www.w3.org/1999/xlink">'
)


def _catalog_contents(drawing: Drawing) -> tuple[GroupContent, ...] | None:
    """The contents whose elements the catalog reads, in document order.

    A record group adds no reads of its own: its parts (``leaf_groups``) are
    read, each when its kind is one the catalog reads. ``None`` when a
    top-level element other than ``defs`` is not a placed content: the catalog
    then reads the SVG.
    """

    contents: list[GroupContent] = []
    for element in drawing.elements:
        if element is drawing.defs:
            continue
        group = held_group(element)
        if group is None:
            return None
        contents.extend(leaf.content for leaf in leaf_groups(group) if leaf.key.kind in _CATALOG_GROUP_KINDS)
    return tuple(contents)


def catalog_content_snapshot(content: GroupContent) -> str:
    """The ``catalog_snapshot`` stage: what the catalog reads from one content.

    A content rebuilt with other geometry or strokes gives the same text, so
    the ``catalog_item`` keyed by the texts is reused.
    """

    from gbdraw.api.prepared import record_prepared_input_metric

    traversal = _catalog_element_reads(
        _parse_catalog_svg((_CATALOG_CONTENTS_ROOT, *content.texts(), "</svg>"))
    )
    for name, value in _catalog_read_metrics(traversal).items():
        record_prepared_input_metric(name, value)
    return _catalog_reads_text(traversal.reads)


def feature_catalog_item_from_snapshots(
    snapshots: tuple[str, ...],
    context: InteractiveSvgContext,
    *,
    result_index: int,
    result_name: str,
) -> dict[str, Any]:
    """The ``catalog_item`` stage: one catalog item from the reads of its contents.

    Keyed by the reads (``catalog_snapshot``), not by the contents: a change
    that moves or restyles a feature, match or annotation without changing
    what the catalog reads (identity attributes and fill) reuses the item and
    its revision. No other group of the drawing holds a feature, match,
    annotation or record source span.
    """

    from gbdraw.api.prepared import record_prepared_input_metric

    svg_index = _rendered_svg_catalog_index(
        [_catalog_reads_from_text(text) for text in snapshots]
    )
    for name, value in _catalog_index_metrics(svg_index).items():
        record_prepared_input_metric(name, value)
    return _feature_catalog_item(
        svg_index, context, result_index=result_index, result_name=result_name
    )


def build_feature_catalog_item(
    svg_source: str | DocumentParts,
    context: InteractiveSvgContext | None,
    *,
    result_index: int,
    result_name: str,
    drawing: Drawing | None = None,
    _diagnostics: MutableMapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Normalize one prepared render context without reading source files.

    In a cache transaction, a ``drawing`` whose groups are all placed contents
    gives one ``catalog_snapshot`` stage per content it reads and the
    ``catalog_item`` stage keyed by their texts; otherwise the catalog reads
    ``svg_source``, a text or a drawing's pieces.
    """

    from gbdraw.api.prepared import prepared_input_cache_active

    if not isinstance(context, InteractiveSvgContext):
        raise GbdrawError("Feature catalog generation requires prepared context.")
    contents = (
        _catalog_contents(drawing)
        if drawing is not None and prepared_input_cache_active()
        else None
    )
    if contents is not None:
        # A content the catalog reads nothing from adds nothing to the item, so it does not key it.
        snapshots = tuple(
            snapshot
            for snapshot in (
                render_stage("catalog_snapshot", catalog_content_snapshot, content)
                for content in contents
            )
            if snapshot != _NO_READS_TEXT
        )
        return render_stage(
            "catalog_item",
            feature_catalog_item_from_snapshots,
            snapshots,
            context,
            result_index=int(result_index),
            result_name=str(result_name),
        )
    traversal = _catalog_element_reads(
        _parse_catalog_svg(svg_source.texts() if isinstance(svg_source, DocumentParts) else svg_source)
    )
    svg_index = _rendered_svg_catalog_index((traversal.reads,))
    _record_catalog_metrics(
        _diagnostics, {**_catalog_read_metrics(traversal), **_catalog_index_metrics(svg_index)}
    )
    return _feature_catalog_item(
        svg_index, context, result_index=result_index, result_name=result_name
    )


def _feature_catalog_item(
    svg_index: _RenderedSvgCatalogIndex,
    context: InteractiveSvgContext,
    *,
    result_index: int,
    result_name: str,
) -> dict[str, Any]:
    # A canceled Worker request stops between the phases of this build, whose state is local (D7).
    from gbdraw.api.prepared import cancel_checkpoint

    record_keys = _record_keys(context)
    identity_keys = _record_identity_keys(context, record_keys)
    stored_features, stored_index = render_stage(
        "catalog_normalization",
        _normalized_biological_features,
        context,
        tuple(record_keys),
        tuple(identity_keys),
    )
    # Per-request copies (M11): the stored ones stay as built.
    biological_features = [dict(feature) for feature in stored_features]
    indexed_features = [
        replace(item, normalized=feature)
        for item, feature in zip(stored_index, biological_features, strict=True)
    ]
    biological_index = _BiologicalFeatureIndex(indexed_features)
    cancel_checkpoint()
    (
        rendered_features,
        rendered_references,
        rendered_match_identities,
    ) = _normalized_rendered_features(
        svg_index,
        context,
        record_keys,
        biological_index,
    )
    cancel_checkpoint()
    raw_matches = _match_payloads(
        svg_index.match_candidates,
        rendered_match_identities,
        svg_index.record_source_spans,
    )
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
        for listed_group_id in _sequence(match.get("orthogroup_ids")):
            normalized_group_id = _text(listed_group_id)
            if not normalized_group_id or normalized_group_id in known_group_ids:
                continue
            orthogroups.append(
                {"id": normalized_group_id, "scope": group_scope, "members": []}
            )
            known_group_ids.add(normalized_group_id)
    cancel_checkpoint()
    presentation_by_id = _orthogroup_presentation_by_id(raw_matches)
    cancel_checkpoint()
    normalized_orthogroups = _normalized_orthogroups(
        orthogroups,
        record_keys,
        biological_index,
        presentation_by_id,
        _presentation_from_search_scope(
            context.collinearity_search_scope
        ),
    )
    cancel_checkpoint()
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
        **({"windowOf": windows} if (windows := _windows_of_records(record_keys, identity_keys)) else {}),
        "features": rendered_features,
        "biologicalFeatures": biological_features,
        "orthogroups": normalized_orthogroups,
        "annotations": _normalized_annotations(
            svg_index.annotation_candidates,
            context,
        ),
        "comparisonMatches": normalized_matches,
    }
    return item


def build_feature_catalog(
    items: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Wrap normalized logical-result items in the current schema envelope."""

    return {
        "schema": FEATURE_CATALOG_SCHEMA,
        "items": [dict(item) for item in items],
    }


def _keyed_sequence_sources(value: object) -> Sequence[object]:
    """A catalog item's sequence sources (none when it has no array), each with its own non-empty key."""

    sequence_sources: Sequence[object] = (
        value if isinstance(value, Sequence) and not isinstance(value, (str, bytes)) else ()
    )
    known_sequence_source_keys: set[str] = set()
    for source in sequence_sources:
        key = _text(source.get("key")) if isinstance(source, Mapping) else ""
        if not key or key in known_sequence_source_keys:
            raise GbdrawError(
                "Feature catalog sequence sources require unique non-empty keys."
            )
        known_sequence_source_keys.add(key)
    return sequence_sources


def _require_catalog_sequence_fields(item: Mapping[str, object], catalog_schema: object) -> None:
    """The item's sequence fields as its schema allows them: ``sequenceSources`` is an array where there is one, and
    the stored schema holds no sequence at all (Python reads them on demand), nor what an export's match spans read
    (``matchSpanReads``)."""

    sources = item.get("sequenceSources")
    if "sequenceSources" in item and (not isinstance(sources, Sequence) or isinstance(sources, (str, bytes))):
        raise GbdrawError("Feature catalog sequenceSources must be an array.")
    features = item["biologicalFeatures"]
    if catalog_schema == FEATURE_CATALOG_SCHEMA and (
        "sequenceSources" in item
        or "matchSpanReads" in item
        or (
            isinstance(features, Sequence)
            and any(
                isinstance(feature, Mapping) and not STORED_SEQUENCE_KEYS.isdisjoint(feature)
                for feature in features
            )
        )
    ):
        raise GbdrawError(
            f"Feature catalog schema {FEATURE_CATALOG_SCHEMA} stores no sequences: Python reads them on demand."
        )


def select_feature_catalog_item(
    catalog: Mapping[str, object],
    *,
    result_index: int,
    result_name: str,
    expected_schema: int | None = None,
) -> dict[str, object]:
    """Return one validated catalog item matched to a logical Result."""

    catalog_schema = catalog.get("schema") if isinstance(catalog, Mapping) else None
    allowed_schemas = (
        {expected_schema}
        if expected_schema is not None
        else {
            LEGACY_FEATURE_CATALOG_SCHEMA,
            ANCHOR_FEATURE_CATALOG_SCHEMA,
            SEQUENCE_FEATURE_CATALOG_SCHEMA,
            FEATURE_CATALOG_SCHEMA,
        }
    )
    if (
        not isinstance(catalog, Mapping)
        or catalog_schema not in allowed_schemas
    ):
        required = (
            str(expected_schema)
            if expected_schema is not None
            else f"3, 4, 5 or {FEATURE_CATALOG_SCHEMA}"
        )
        raise GbdrawError(
            f"Feature catalog must use schema {required}."
        )
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
    _require_catalog_sequence_fields(item, catalog_schema)

    record_keys = [_text(value) for value in item["recordKeys"]]
    if any(not key for key in record_keys) or len(set(record_keys)) != len(
        record_keys
    ):
        raise GbdrawError(
            "Feature catalog record keys must be non-empty and unique."
        )
    record_units = catalog_record_units(item, record_keys)
    known_features: set[tuple[str, str]] = set()
    biological_source_indexes: dict[
        tuple[str, str],
        list[int | None],
    ] = defaultdict(list)
    sequence_sources = _keyed_sequence_sources(item.get("sequenceSources"))
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
            or reference[0] not in record_units
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
        if catalog_schema >= 4:
            profile = feature.get("anchorProfile")
            if (
                not isinstance(profile, Mapping)
                or set(profile) != {
                    "precision",
                    "operator",
                    "partOrder",
                    "strand",
                }
                or profile.get("precision")
                not in {"exact", "fuzzy", "unavailable"}
                or profile.get("operator")
                not in {"single", "join", "order", "unknown"}
                or profile.get("partOrder")
                not in {"biological", "source-forward", "ambiguous"}
                or profile.get("strand")
                not in {"+", "-", "unstranded", "mixed"}
            ):
                raise GbdrawError(
                    "Feature catalog contains an invalid source anchor profile."
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
            in_range = type(source_index) is int and 0 <= source_index < len(sequence_sources)
            indexed_source = sequence_sources[cast(int, source_index)] if in_range else None
            source = indexed_source if isinstance(indexed_source, Mapping) else None
            if (
                "nucleotide_sequence" in feature
                or "nucleotideSequence" in feature
                or source is None
                or _text(source.get("origin"))
                not in {"linear-record", "circular-reference"}
                or type(source.get("recordIndex")) is not int
                or source.get("recordIndex") != record_units[reference[0]][0]
                or not _reconstruct_nucleotide_sequence(
                    feature,
                    source,
                    # source is not None only when source_index is an int.
                    source_sequence=canonical_sources[cast(int, source_index)],
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
        # A row of a record drawn in windows names its unit; any other row's unit is its record's.
        windowed = reference[0] not in record_keys
        unit = feature.get("recordIndex", _MISSING)
        if (
            not svg_id
            or svg_id in rendered_ids
            or reference not in known_features
            or windowed != (unit is not _MISSING)
            or (windowed and (type(unit) is not int or unit not in record_units[reference[0]]))
        ):
            raise GbdrawError(
                "Feature catalog contains an invalid rendered feature "
                "reference."
            )
        if catalog_schema >= SEQUENCE_FEATURE_CATALOG_SCHEMA and not _valid_drawn_selector(
            feature.get("drawnSelector", _MISSING)
        ):
            raise GbdrawError(
                "Feature catalog contains an invalid rendered feature "
                "drawn selector."
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
    "interactive_embed_catalog",
    "build_feature_catalog_item",
    "canonical_catalog_sequence_sources",
    "materialize_catalog_nucleotide_sequence",
    "STORED_SEQUENCE_KEYS",
    "INTERACTIVE_EMBED_SCHEMA",
    "EmbedCatalogPart",
    "EmbeddedSequenceSource",
    "embed_catalog_item",
    "embed_catalog_part",
    "embed_json_text",
    "match_sequence_sources",
    "promote_legacy_feature_catalog",
    "read_match_span",
    "select_feature_catalog_item",
]
