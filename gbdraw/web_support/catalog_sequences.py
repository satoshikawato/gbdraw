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
  and the whole-record sources they own, ``embed_catalog_part``), which the Web export reads in parts of bounded size, as
  the export's text (``read_embed_sequences_json``).

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
    ComparisonMatchSpan,
    EmbedCatalogPart,
    MatchSpanInput,
    MatchSpanReason,
    RecordMatchSpan,
    ReferenceMatchSpan,
    SpanCoordinate,
    embed_catalog_item,
    embed_catalog_part,
    embed_json_text,
    match_sequence_sources,
    match_span_source,
    read_match_span,
    registered_match_sources,
    span_reason,
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


class MatchSpanSequence(TypedDict):
    """The bases of one span, on its strand (``orientation``), or why they cannot be read (``reason``);
    ``sourceLength`` is the length of the record the span reads (0 when no record is found)."""

    status: Literal["ok", "unavailable"]
    reason: MatchSpanReason | None
    recordId: str
    orientation: Literal["+", "-"]
    length: int
    sequence: str
    sourceLength: int


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


def _span_coordinate(value: object) -> SpanCoordinate:
    """A whole number (JSON writes one of 1e21 or more in exponent form, which reads as a float), or ``None``."""

    if value is None or type(value) is int:
        return value
    if type(value) is float and value.is_integer():
        return int(value)
    raise _protocol_error("A match span needs a start and an end.")


def _span_inputs(value: object) -> list[MatchSpanInput]:
    if not isinstance(value, list) or not value:
        raise _protocol_error("A match sequence read needs a list of spans.")
    spans: list[MatchSpanInput] = []
    for item in value:
        if not isinstance(item, Mapping) or not {"start", "end"} <= set(item):
            raise _protocol_error("A match span needs a start and an end.")
        start, end = _span_coordinate(item["start"]), _span_coordinate(item["end"])
        if (
            set(item) == {"kind", "recordKey", "recordId", "start", "end"}
            and item["kind"] == "record"
            and isinstance(item["recordKey"], str)
            and isinstance(item["recordId"], str)
        ):
            spans.append({
                "kind": "record", "recordKey": item["recordKey"], "recordId": item["recordId"], "start": start, "end": end,
            })
        elif (
            set(item) == {"kind", "recordId", "recordKey", "start", "end"}
            and item["kind"] == "reference"
            and isinstance(item["recordId"], str)
            and (item["recordKey"] is None or isinstance(item["recordKey"], str))
        ):
            spans.append({
                "kind": "reference", "recordId": item["recordId"], "recordKey": item["recordKey"], "start": start, "end": end,
            })
        elif (
            set(item) == {"kind", "sourceIndex", "recordId", "start", "end"}
            and item["kind"] == "comparison"
            and (item["sourceIndex"] is None or type(item["sourceIndex"]) is int)
            and isinstance(item["recordId"], str)
        ):
            spans.append({
                "kind": "comparison", "sourceIndex": item["sourceIndex"], "recordId": item["recordId"],
                "start": start, "end": end,
            })
        else:
            raise _protocol_error(
                "A match span names a drawn record by key or ID, or a comparison record by source and ID."
            )
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
    sources = registered_match_sources(match_sequence_sources(records))
    read: list[MatchSpanSequence] = []
    for span in asked:
        start, end = span["start"], span["end"]
        # A coordinate that is not whole has no strand; the page shows its own for that span.
        orientation: Literal["+", "-"] = "-" if start is not None and end is not None and start > end else "+"
        source, reason = match_span_source(records, sources, span)
        sequence = "" if source is None else cast("str", source["sequence"])
        reason = reason or span_reason(sequence, start, end)
        bases = (
            None if source is None or reason or start is None or end is None
            else read_match_span(sequence, start, end)
        )
        if source is None or bases is None:
            read.append({
                "status": "unavailable", "reason": reason or "OUT_OF_RANGE",
                "recordId": "" if source is None else str(source["recordId"]),
                "orientation": orientation, "length": 0, "sequence": "", "sourceLength": len(sequence),
            })
            continue
        read.append({
            "status": "ok", "reason": None, "recordId": str(source["recordId"]),
            "orientation": orientation, "length": len(bases), "sequence": bases, "sourceLength": len(sequence),
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


@final
class EmbeddedSourceText(TypedDict):
    """A whole-record source of an export part as the export's text, with its place among the item's sources."""

    index: int
    source: str


@final
class EmbedPartText(TypedDict):
    """The reply of one export part (``read_embed_sequences_json``): its rows and sources as the text the file
    embeds (``embed_json_text``), so neither the Worker nor the page parses them; ``biologicalFeatures`` is the rows'
    text joined by commas, ``rowCount`` their number; ``matchSpanReads`` is the text of what the item's match spans
    read, in the part that holds the comparison sources (None in the others and when no match names a span), whose
    ``sourceCount`` also counts the sources only a span reads."""

    biologicalFeatures: str
    rowCount: int
    sequenceSources: list[EmbeddedSourceText]
    sourceCount: int
    matchSpanReads: str | None


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
    request_json: str, resource_paths_json: str, output_directory: str, item_path: str, part_json: str
) -> str:
    """JSON boundary of the diagram Worker: one part of the export encoding, as its text (``EmbedPartText``). The
    part's item comes as the UTF-8 JSON file the Worker staged (``embedItem``)."""

    with open(str(item_path), encoding="utf-8") as file:
        item = json.load(file)
    part = read_embed_catalog_part(
        json.loads(str(request_json)), json.loads(str(resource_paths_json)), str(output_directory),
        item, json.loads(str(part_json)),
    )
    reply: EmbedPartText = {
        "biologicalFeatures": ",".join(embed_json_text(row) for row in part["biologicalFeatures"]),
        "rowCount": len(part["biologicalFeatures"]),
        "sequenceSources": [
            {"index": entry["index"], "source": embed_json_text(entry["source"])} for entry in part["sequenceSources"]
        ],
        "sourceCount": part["sourceCount"],
        "matchSpanReads": None if part["matchSpanReads"] is None else embed_json_text(part["matchSpanReads"]),
    }
    return embed_json_text(reply)


__all__ = [
    "ComparisonMatchSpan",
    "EmbedPartInput",
    "EmbedPartText",
    "EmbeddedSourceText",
    "FeatureSequenceInput",
    "FeatureSequencesOf",
    "FeatureSequencesReply",
    "MatchSequencesReply",
    "MatchSpanSequence",
    "RecordMatchSpan",
    "ReferenceMatchSpan",
    "read_embed_catalog_item",
    "read_embed_catalog_part",
    "read_embed_sequences_json",
    "read_feature_sequences",
    "read_feature_sequences_json",
    "read_match_sequences",
    "read_match_sequences_json",
]
