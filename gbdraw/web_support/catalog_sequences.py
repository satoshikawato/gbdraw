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
from typing import Literal, TypedDict, cast, final, get_args

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
    embed_json_text,
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


def _aliases(source: Mapping[str, object]) -> set[str]:
    """A source's record ID and aliases, each also by its first word (the match registry's rule before D5)."""

    names = {str(source["recordId"]).strip(), *(str(alias).strip() for alias in cast("list[object]", source["aliases"]))}
    return {name for value in names if value for name in (value, value.split()[0])}


def _registered(sources: Sequence[Mapping[str, object]]) -> list[Mapping[str, object]]:
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


def _match_source(
    records: CatalogSequenceRecords,
    sources: Sequence[Mapping[str, object]],
    span: MatchSpanInput,
) -> tuple[Mapping[str, object] | None, MatchSpanReason | None]:
    """The source a span reads among the ``_registered`` sources, found as the match popup found it before D5: a drawn
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


def _span_reason(sequence: str, start: SpanCoordinate, end: SpanCoordinate) -> MatchSpanReason | None:
    """Why the coordinates of a span cannot be read on ``sequence``, checked in the match popup's order before D5."""

    if start is None or end is None:
        return "COORDINATE_NOT_WHOLE"
    if start < 1 or end < 1:
        return "COORDINATE_BELOW_ONE"
    if start > len(sequence) or end > len(sequence):
        return "OUT_OF_RANGE"
    return None


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
    sources = _registered(match_sequence_sources(records))
    read: list[MatchSpanSequence] = []
    for span in asked:
        start, end = span["start"], span["end"]
        # A coordinate that is not whole has no strand; the page shows its own for that span.
        orientation: Literal["+", "-"] = "-" if start is not None and end is not None and start > end else "+"
        source, reason = _match_source(records, sources, span)
        sequence = "" if source is None else cast("str", source["sequence"])
        reason = reason or _span_reason(sequence, start, end)
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
    text joined by commas, ``rowCount`` their number."""

    biologicalFeatures: str
    rowCount: int
    sequenceSources: list[EmbeddedSourceText]
    sourceCount: int


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
    "MATCH_SPAN_REASONS",
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
