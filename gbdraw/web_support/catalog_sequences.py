"""Worker helpers that read the sequences the feature catalog does not store (D5, PD-OI for
``web.feature-sequence.on-demand``).

The page sends the request the shown Result was made from (the committed request, with its staged resources), as for
Copy FASTA of Select region (``record_spans.py``). Python resolves the records as the render does and reads:

- ``read_feature_sequences``: the sequences of features named by identity (record key and biological feature ID),
  from ``feature_sequence_payload``, the one producer: the Sequence tab, member Copy and Download, and the sequence
  fields of Feature search.
- ``read_match_sequences``: the bases of match spans, on the displayed strand of a drawn record (named by its record
  key) or of a comparison FASTA record (named by its source and record ID): the Copy and FASTA of a match popup.
- ``read_embed_catalog_item``: a stored catalog item with every sequence an exported Interactive SVG embeds, in the
  export encoding (``embed_catalog_item``); ``read_embed_catalog_part`` reads one part of it (the rows of some records
  and the whole-record sources they own, ``embed_catalog_part``), which the Web export reads record by record.

Nothing is found by caption, label or display text.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from typing import Literal, TypedDict, cast, final

from gbdraw.api.record_planning import ResolvedRecordCollection
from gbdraw.api.request_render import load_comparison_sequence_records, resolve_request_records
from gbdraw.api.requests import DiagramRequest, LinearDiagramRequest
from gbdraw.exceptions import ValidationError
from gbdraw.features.source import FeatureIdentity, resolve_feature_identities
from gbdraw.render.interactive_svg import CatalogSequenceRecords
from gbdraw.session_request_codec import decode_canonical_request
from gbdraw.web_support.feature_catalog import (
    EmbedCatalogPart,
    embed_catalog_item,
    embed_catalog_part,
    match_sequence_sources,
    read_match_span,
)
from gbdraw.web_support.feature_metadata import feature_sequence_payload


@final
class FeatureSequenceInput(TypedDict):
    """One feature, by identity: the record key and the biological feature ID of the catalog row."""

    recordKey: str
    biologicalFeatureId: str


class FeatureSequencesOf(TypedDict):
    """One feature's sequences (``FeatureSequences`` of ``feature_metadata``) with the identity asked."""

    recordKey: str
    biologicalFeatureId: str
    nucleotideSequence: str
    aminoAcidSequence: str
    sequenceWarnings: list[str]


class FeatureSequencesReply(TypedDict):
    """The features in the order asked."""

    features: list[FeatureSequencesOf]


@final
class RecordMatchSpan(TypedDict):
    """A span of a drawn record, by the record key of its drawn unit: displayed coordinates, 1-based and inclusive;
    a start after the end is the reverse strand."""

    kind: Literal["record"]
    recordKey: str
    start: int
    end: int


@final
class ComparisonMatchSpan(TypedDict):
    """A span of a comparison FASTA record, by the index of its comparison source and its record ID or alias."""

    kind: Literal["comparison"]
    sourceIndex: int
    recordId: str
    start: int
    end: int


MatchSpanInput = RecordMatchSpan | ComparisonMatchSpan
MatchSpanReason = Literal["NEVER_SUPPLIED", "RECORD_UNKNOWN", "RECORD_AMBIGUOUS", "OUT_OF_RANGE"]


class MatchSpanSequence(TypedDict):
    """The bases of one span, on its strand (``orientation``), or why they cannot be read (``reason``)."""

    status: Literal["ok", "unavailable"]
    reason: MatchSpanReason | None
    recordId: str
    orientation: Literal["+", "-"]
    length: int
    sequence: str


class MatchSequencesReply(TypedDict):
    """The spans in the order asked."""

    spans: list[MatchSpanSequence]


def _protocol_error(message: str) -> ValidationError:
    return ValidationError(message, diagnostic={"code": "HELPER_PROTOCOL"})


def _feature_inputs(value: object) -> list[FeatureSequenceInput]:
    if not isinstance(value, list) or not value:
        raise _protocol_error("A feature sequence read needs a list of features.")
    features: list[FeatureSequenceInput] = []
    for item in value:
        if not (
            isinstance(item, Mapping)
            and set(item) == {"recordKey", "biologicalFeatureId"}
            and all(isinstance(item[key], str) and item[key] for key in ("recordKey", "biologicalFeatureId"))
        ):
            raise _protocol_error("A feature needs a record key and a biological feature ID.")
        features.append({"recordKey": item["recordKey"], "biologicalFeatureId": item["biologicalFeatureId"]})
    if len({(item["recordKey"], item["biologicalFeatureId"]) for item in features}) != len(features):
        raise _protocol_error("A feature sequence read names each feature once.")
    return features


def _span_inputs(value: object) -> list[MatchSpanInput]:
    if not isinstance(value, list) or not value:
        raise _protocol_error("A match sequence read needs a list of spans.")
    spans: list[MatchSpanInput] = []
    for item in value:
        if not isinstance(item, Mapping) or not all(type(item.get(key)) is int for key in ("start", "end")):
            raise _protocol_error("A match span needs a start and an end.")
        if set(item) == {"kind", "recordKey", "start", "end"} and item["kind"] == "record" and isinstance(
            item["recordKey"], str
        ):
            spans.append({"kind": "record", "recordKey": item["recordKey"], "start": item["start"], "end": item["end"]})
        elif (
            set(item) == {"kind", "sourceIndex", "recordId", "start", "end"}
            and item["kind"] == "comparison"
            and type(item["sourceIndex"]) is int
            and isinstance(item["recordId"], str)
        ):
            spans.append({
                "kind": "comparison", "sourceIndex": item["sourceIndex"], "recordId": item["recordId"],
                "start": item["start"], "end": item["end"],
            })
        else:
            raise _protocol_error("A match span names a drawn record by key or a comparison record by source and ID.")
    return spans


def _committed_records(
    request_payload: Mapping[str, object],
    resource_paths: Mapping[str, str],
    output_directory: str,
) -> tuple[DiagramRequest, ResolvedRecordCollection]:
    request = decode_canonical_request(request_payload, resource_paths=resource_paths, output_directory=output_directory)
    return request, resolve_request_records(request)


def _sequence_records(
    request: DiagramRequest,
    collection: ResolvedRecordCollection,
    unit_keys: Sequence[str] | None = None,
) -> CatalogSequenceRecords:
    """The drawn units (all, or those of ``unit_keys`` in that order) and comparison records of a request."""

    units = {item.record_key: index for index, item in enumerate(collection.provenance)}
    if unit_keys is None:
        indexes = list(range(len(collection.records)))
    else:
        missing = [key for key in unit_keys if key not in units]
        if missing:
            raise ValidationError(
                f"The drawing has no record {missing[0]!r}.",
                diagnostic={"code": "INPUT_INVALID", "field": "record_selector"},
            )
        indexes = [units[key] for key in unit_keys]
    return CatalogSequenceRecords(
        mode="linear" if isinstance(request, LinearDiagramRequest) else "circular",
        records=tuple(collection.records[index] for index in indexes),
        unit_keys=tuple(collection.provenance[index].record_key for index in indexes),
        comparison_records=load_comparison_sequence_records(request),
    )


def _aliases(source: Mapping[str, object]) -> set[str]:
    """A source's record ID and aliases, each also by its first word (the match registry's rule before D5)."""

    names = {str(source["recordId"]).strip(), *(str(alias).strip() for alias in cast("list[object]", source["aliases"]))}
    return {name for value in names if value for name in (value, value.split()[0])}


def _match_source(
    records: CatalogSequenceRecords,
    sources: Sequence[Mapping[str, object]],
    span: MatchSpanInput,
) -> tuple[Mapping[str, object] | None, MatchSpanReason | None]:
    """The source a span reads: a drawn unit by its record key, or a comparison record by its source and its record
    ID, else an alias; ambiguity and absence are reasons, never a guess."""

    if span["kind"] == "record":
        if span["recordKey"] not in records.unit_keys:
            return None, "RECORD_UNKNOWN"
        unit = records.unit_keys.index(span["recordKey"])
        return next(source for source in sources if source.get("recordIndex") == unit and "sourceIndex" not in source), None
    index = span["sourceIndex"]
    if not 0 <= index < len(records.comparison_records) or not records.comparison_records[index]:
        return None, "NEVER_SUPPLIED"
    candidates = [source for source in sources if source.get("sourceIndex") == index]
    for found in (
        [source for source in candidates if source["recordId"] == span["recordId"]],
        [source for source in candidates if span["recordId"] in _aliases(source)],
    ):
        if len(found) == 1:
            return found[0], None
        if found:
            return None, "RECORD_AMBIGUOUS"
    return None, "RECORD_UNKNOWN"


def read_feature_sequences(
    request_payload: Mapping[str, object],
    resource_paths: Mapping[str, str],
    output_directory: str,
    features: Sequence[FeatureSequenceInput],
) -> FeatureSequencesReply:
    """The sequences of ``features`` of the records that the canonical request ``request_payload`` draws."""

    asked = _feature_inputs(list(features))
    _request, collection = _committed_records(request_payload, resource_paths, output_directory)
    identities = [FeatureIdentity(record_key=item["recordKey"], biological_feature_id=item["biologicalFeatureId"]) for item in asked]
    bindings = resolve_feature_identities(
        records=collection.records,
        record_keys=tuple(item.record_identity_key for item in collection.provenance),
        source_catalogs=tuple(item.source_feature_catalog for item in collection.provenance),
        identities=identities,
    )
    read: list[FeatureSequencesOf] = []
    for item, identity in zip(asked, identities, strict=True):
        # A feature drawn in several windows of its record is one feature: the first unit that has it (D-18).
        binding = bindings[identity][0]
        if binding.status != "present" or binding.feature is None:
            # The code of the binding owner (``resolve_feature_identities``) for a record key outside the request.
            raise ValidationError(
                f"The drawing has no feature {identity.biological_feature_id!r} in {identity.record_key!r}"
                f" ({binding.status.replace('_', ' ')}).",
                diagnostic={"code": "FEATURE_IDENTITY"},
            )
        read.append({**item, **feature_sequence_payload(binding.feature, collection.records[binding.record_index])})
    return {"features": read}


def read_match_sequences(
    request_payload: Mapping[str, object],
    resource_paths: Mapping[str, str],
    output_directory: str,
    spans: Sequence[MatchSpanInput],
) -> MatchSequencesReply:
    """The bases of match ``spans`` of the drawn or comparison records of the canonical request ``request_payload``."""

    asked = _span_inputs(list(spans))
    request, collection = _committed_records(request_payload, resource_paths, output_directory)
    records = _sequence_records(request, collection)
    sources = match_sequence_sources(records)
    read: list[MatchSpanSequence] = []
    for span in asked:
        orientation: Literal["+", "-"] = "+" if span["start"] <= span["end"] else "-"
        source, reason = _match_source(records, sources, span)
        bases = None if source is None else read_match_span(cast("str", source["sequence"]), span["start"], span["end"])
        if source is None or bases is None:
            read.append({
                "status": "unavailable", "reason": reason or "OUT_OF_RANGE",
                "recordId": "" if source is None else str(source["recordId"]),
                "orientation": orientation, "length": 0, "sequence": "",
            })
            continue
        read.append({
            "status": "ok", "reason": None, "recordId": str(source["recordId"]),
            "orientation": orientation, "length": len(bases), "sequence": bases,
        })
    return {"spans": read}


def read_embed_catalog_item(
    request_payload: Mapping[str, object],
    resource_paths: Mapping[str, str],
    output_directory: str,
    item: Mapping[str, object],
) -> dict[str, object]:
    """The stored catalog ``item`` of a Result of the canonical request, with the sequences an exported Interactive
    SVG embeds, in the export encoding (``embed_catalog_item``, the encoder the CLI export uses)."""

    if not isinstance(item, Mapping) or not isinstance(item.get("recordKeys"), list):
        raise _protocol_error("An embed read needs a stored catalog item.")
    request, collection = _committed_records(request_payload, resource_paths, output_directory)
    unit_keys = [str(key) for key in cast("list[object]", item["recordKeys"])]
    return embed_catalog_item(item, _sequence_records(request, collection, unit_keys))


@final
class EmbedPartInput(TypedDict):
    """Which whole-record sources a part of the export returns: those of the records named by record key
    (``sourcesOf``) and, with ``comparisonSources``, those of the comparison FASTA records."""

    sourcesOf: list[str]
    comparisonSources: bool


def _embed_part_input(value: object) -> EmbedPartInput:
    if not (
        isinstance(value, Mapping)
        and set(value) == {"sourcesOf", "comparisonSources"}
        and isinstance(value["sourcesOf"], list)
        and all(isinstance(key, str) and key for key in value["sourcesOf"])
        and isinstance(value["comparisonSources"], bool)
    ):
        raise _protocol_error("An embed part names the records whose sources it returns.")
    return {"sourcesOf": list(value["sourcesOf"]), "comparisonSources": value["comparisonSources"]}


def read_embed_catalog_part(
    request_payload: Mapping[str, object],
    resource_paths: Mapping[str, str],
    output_directory: str,
    item: Mapping[str, object],
    part: Mapping[str, object],
) -> EmbedCatalogPart:
    """A part of the export encoding of a stored catalog item (``embed_catalog_part``): ``item`` holds the item's
    ``recordKeys``, ``windowOf`` and ``comparisonMatches`` and the ``biologicalFeatures`` of the part; ``part`` names the
    records whose whole-record sources it returns."""

    if not (
        isinstance(item, Mapping)
        and isinstance(item.get("recordKeys"), list)
        and isinstance(item.get("biologicalFeatures"), list)
    ):
        raise _protocol_error("An embed read needs the rows of a stored catalog item.")
    asked = _embed_part_input(part)
    request, collection = _committed_records(request_payload, resource_paths, output_directory)
    unit_keys = [str(key) for key in cast("list[object]", item["recordKeys"])]
    return embed_catalog_part(
        item, _sequence_records(request, collection, unit_keys),
        source_record_keys=asked["sourcesOf"], comparison_sources=asked["comparisonSources"],
    )


def read_feature_sequences_json(
    request_json: str, resource_paths_json: str, output_directory: str, features_json: str
) -> str:
    """JSON boundary of the diagram Worker."""

    return json.dumps(read_feature_sequences(
        json.loads(str(request_json)), json.loads(str(resource_paths_json)), str(output_directory),
        json.loads(str(features_json)),
    ))


def read_match_sequences_json(request_json: str, resource_paths_json: str, output_directory: str, spans_json: str) -> str:
    """JSON boundary of the diagram Worker."""

    return json.dumps(read_match_sequences(
        json.loads(str(request_json)), json.loads(str(resource_paths_json)), str(output_directory),
        json.loads(str(spans_json)),
    ))


def read_embed_sequences_json(
    request_json: str, resource_paths_json: str, output_directory: str, item_json: str, part_json: str
) -> str:
    """JSON boundary of the diagram Worker: one part of the export encoding."""

    return json.dumps(read_embed_catalog_part(
        json.loads(str(request_json)), json.loads(str(resource_paths_json)), str(output_directory),
        json.loads(str(item_json)), json.loads(str(part_json)),
    ), ensure_ascii=False, separators=(",", ":"))


__all__ = [
    "ComparisonMatchSpan",
    "EmbedPartInput",
    "FeatureSequenceInput",
    "FeatureSequencesOf",
    "FeatureSequencesReply",
    "MatchSequencesReply",
    "MatchSpanSequence",
    "RecordMatchSpan",
    "read_embed_catalog_item",
    "read_embed_catalog_part",
    "read_embed_sequences_json",
    "read_feature_sequences",
    "read_feature_sequences_json",
    "read_match_sequences",
    "read_match_sequences_json",
]
