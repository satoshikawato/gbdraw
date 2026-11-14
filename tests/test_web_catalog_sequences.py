"""The sequences the feature catalog does not store (D5): read on demand from Python, by identity.

The stored catalog (schema 6) holds no sequence. ``gbdraw.web_support.catalog_sequences`` reads a feature's
sequences, the bases of a match span and the export embed from the records of the committed request, through the
one producer (``feature_sequence_payload``) and the whole-record sources (``match_sequence_sources``). These tests
pin that the values equal what dev stored before D5, that a feature is found by identity and never by caption, and
the reasons a read refuses.
"""

from __future__ import annotations

import json
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api import InMemoryRecordSource, LinearDiagramRequest, RecordInput, RenderOutputRequest
from gbdraw.exceptions import ValidationError
from gbdraw.render.interactive_svg import CatalogSequenceRecords
from gbdraw.session import build_session_document, load_session_document, materialize_session
from gbdraw.web_support.catalog_sequences import (
    read_embed_catalog_item,
    read_embed_catalog_part,
    read_embed_sequences_json,
    read_feature_sequences,
    read_match_sequences,
)
from gbdraw.web_support.feature_catalog import (
    FEATURE_CATALOG_SCHEMA,
    STORED_SEQUENCE_KEYS,
    embed_catalog_item,
    embed_catalog_part,
    match_sequence_sources,
    match_span_source,
    materialize_catalog_nucleotide_sequence,
    promote_legacy_feature_catalog,
    read_match_span,
    registered_match_sources,
    span_reason,
)
from gbdraw.web_support.request_render import render_embedded_canonical_web_request
from tests.utils.session_drawing_view import drawing_view
from tests.utils.session_files import embedded_session

SESSIONS = Path(__file__).resolve().parent / "fixtures" / "sessions"
SPAN_VECTORS = Path(__file__).resolve().parent / "fixtures" / "match-span-vectors.json"
# Sessions 40 whose catalogs (schema 3) store every sequence as dev wrote it: inline rows, translation marks, and
# (BGC) rows that reference whole-record sources, with comparison matches.
STORED_SEQUENCE_SESSIONS = ("HmmtDNA_ATskew", "lambda_basic_linear", "BGC0000708-BGC0000713")


@contextmanager
def _committed(path: Path, tmp_path: Path) -> Iterator[tuple[dict[str, object], dict[str, str], str]]:
    """The committed request of a Session's only drawing, with its resources as files."""

    document = load_session_document(path)
    with materialize_session(document, output_directory=tmp_path / "outputs") as materialized:
        payload = document.to_dict()
        request = drawing_view(payload)["renderRequest"] if "drawings" in payload else payload["renderRequest"]
        paths = {str(key): str(value) for key, value in materialized.resource_paths.items()}
        yield request, paths, str(materialized.output_directory)


def _stored_sequences(row: Mapping[str, object], sources: list[dict[str, object]]) -> dict[str, object]:
    """What dev stored for a row: the nucleotide sequence inline or by reference, the amino-acid sequence (kept also
    when the translation mark replaced the /translation qualifier) and the warnings."""

    return {
        "nucleotideSequence": materialize_catalog_nucleotide_sequence(row, sources),
        "aminoAcidSequence": row.get("amino_acid_sequence") or "",
        "sequenceWarnings": list(row.get("sequence_warnings") or []),  # type: ignore[call-overload]
    }


@pytest.mark.parametrize("name", STORED_SEQUENCE_SESSIONS)
def test_feature_sequences_equal_the_values_dev_stored(name: str, tmp_path: Path) -> None:
    path = SESSIONS / f"{name}.v40-schema5.json"
    [item] = json.loads(path.read_text(encoding="utf-8"))["editorState"]["featureCatalog"]["items"]
    rows = item["biologicalFeatures"]
    sources = list(item.get("sequenceSources") or [])
    asked = [{"recordKey": row["recordKey"], "biologicalFeatureId": row["biologicalFeatureId"]} for row in rows]

    with _committed(path, tmp_path) as (request, paths, output):
        read = read_feature_sequences(request, paths, output, asked)["features"]

    assert [(row["recordKey"], row["biologicalFeatureId"]) for row in read] == [
        (row["recordKey"], row["biologicalFeatureId"]) for row in asked
    ]
    differences = [
        row["biologicalFeatureId"]
        for row, got in zip(rows, read, strict=True)
        if _stored_sequences(row, sources)
        != {key: got[key] for key in ("nucleotideSequence", "aminoAcidSequence", "sequenceWarnings")}  # type: ignore[literal-required]
    ]
    assert differences == []
    assert all(row["nucleotideSequence"] for row in read)


def test_embed_of_the_stored_item_equals_the_item_dev_stored(tmp_path: Path) -> None:
    """The export embed (schema 5) of the promoted, sequence-free item is the item dev stored.

    Compared as objects: this catalog is schema 3, whose promotion appends keys after the sequence keys; the export
    bytes of current catalogs are pinned by the CLI export test."""

    path = SESSIONS / "BGC0000708-BGC0000713.v40-schema5.json"
    catalog = json.loads(path.read_text(encoding="utf-8"))["editorState"]["featureCatalog"]
    stored_five = promote_legacy_feature_catalog(catalog, schema=5)["items"][0]
    [bare] = promote_legacy_feature_catalog(catalog)["items"]
    assert not any(STORED_SEQUENCE_KEYS & set(row) for row in bare["biologicalFeatures"])
    assert "sequenceSources" not in bare

    with _committed(path, tmp_path) as (request, paths, output):
        embedded = read_embed_catalog_item(request, paths, output, bare)

    assert embedded == stored_five
    assert sum("sequenceSourceIndex" in row for row in embedded["biologicalFeatures"]) == 342


def _assembled_from_parts(item: Mapping[str, object], read_part: object) -> dict[str, object]:
    """The export item as the Web export writes it: one part per record of the item, in the order of its drawn
    records, the comparison sources with the last part; rows in their order, sources placed by their index."""

    windows = item.get("windowOf") or {}
    order = list(dict.fromkeys(str(windows.get(key, key)) for key in item["recordKeys"]))  # type: ignore[union-attr, attr-defined]
    rows: list[dict[str, object]] = []
    sources: dict[int, dict[str, object]] = {}
    reads = None
    counts = set()
    for position, key in enumerate(order):
        part_rows = [row for row in item["biologicalFeatures"] if row["recordKey"] == key]  # type: ignore[attr-defined]
        context = {name: item[name] for name in ("recordKeys", "windowOf", "comparisonMatches") if name in item}
        reply = read_part(  # type: ignore[operator]
            {**context, "biologicalFeatures": part_rows},
            {"sourcesOf": [key], "comparisonSources": position == len(order) - 1},
        )
        assert [(row["recordKey"], row["biologicalFeatureId"]) for row in reply["biologicalFeatures"]] == [
            (row["recordKey"], row["biologicalFeatureId"]) for row in part_rows
        ]
        rows.extend(reply["biologicalFeatures"])
        counts.add(reply["sourceCount"])
        for entry in reply["sequenceSources"]:
            assert entry["index"] not in sources
            sources[entry["index"]] = entry["source"]
        # Only the last part, which holds the comparison sources, holds what the match spans read.
        assert reply["matchSpanReads"] is None or position == len(order) - 1
        reads = reply["matchSpanReads"]
    [count] = counts
    assert sorted(sources) == list(range(count))
    embedded: dict[str, object] = {**item, "biologicalFeatures": rows}
    if sources:
        embedded["sequenceSources"] = [sources[index] for index in range(count)]
    if reads is not None:
        embedded["matchSpanReads"] = reads
    return embedded


def _compact(value: object) -> str:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":"))


def test_the_export_read_record_by_record_equals_the_whole_item_byte_for_byte(tmp_path: Path) -> None:
    """The Web export reads the embed one record at a time (no reply holds the whole item); assembled, the parts are
    the item the one encoder writes for the CLI export, with the same source references."""

    path = SESSIONS / "BGC0000708-BGC0000713.v40-schema5.json"
    catalog = json.loads(path.read_text(encoding="utf-8"))["editorState"]["featureCatalog"]
    [bare] = promote_legacy_feature_catalog(catalog)["items"]

    with _committed(path, tmp_path) as (request, paths, output):
        whole = read_embed_catalog_item(request, paths, output, bare)
        assembled = _assembled_from_parts(
            bare, lambda item, part: read_embed_catalog_part(request, paths, output, item, part)
        )

    assert len(whole["sequenceSources"]) == 5  # type: ignore[arg-type]
    assert _compact(assembled) == _compact(whole)


def test_an_embed_part_names_the_records_whose_sources_it_returns(tmp_path: Path) -> None:
    path = SESSIONS / "BGC0000708-BGC0000713.v40-schema5.json"
    catalog = json.loads(path.read_text(encoding="utf-8"))["editorState"]["featureCatalog"]
    [bare] = promote_legacy_feature_catalog(catalog)["items"]
    item = {"recordKeys": bare["recordKeys"], "comparisonMatches": bare["comparisonMatches"], "biologicalFeatures": []}

    with _committed(path, tmp_path) as (request, paths, output):
        reply = read_embed_catalog_part(
            request, paths, output, item, {"sourcesOf": [bare["recordKeys"][1]], "comparisonSources": False}
        )
        with pytest.raises(ValidationError, match="names the records"):
            read_embed_catalog_part(request, paths, output, item, {"sourcesOf": "record-1", "comparisonSources": False})

    assert reply["biologicalFeatures"] == []
    assert [entry["index"] for entry in reply["sequenceSources"]] == [1]
    assert reply["sequenceSources"][0]["source"]["key"] == "linear:record:1"
    assert reply["sourceCount"] == 5


def test_the_worker_reply_of_an_embed_part_is_the_text_the_file_embeds(tmp_path: Path) -> None:
    """The Worker and the page pass a part on as text: the reply holds each row and source as the encoder writes it,
    so the page parses none of them."""

    path = SESSIONS / "BGC0000708-BGC0000713.v40-schema5.json"
    catalog = json.loads(path.read_text(encoding="utf-8"))["editorState"]["featureCatalog"]
    [bare] = promote_legacy_feature_catalog(catalog)["items"]
    key = bare["recordKeys"][1]
    item = {"recordKeys": bare["recordKeys"], "comparisonMatches": bare["comparisonMatches"],
            "biologicalFeatures": [row for row in bare["biologicalFeatures"] if row["recordKey"] == key]}
    part = {"sourcesOf": [key], "comparisonSources": True}

    item_file = tmp_path / "embed-item.json"
    item_file.write_text(json.dumps(item), encoding="utf-8")

    with _committed(path, tmp_path) as (request, paths, output):
        objects = read_embed_catalog_part(request, paths, output, item, part)
        reply = json.loads(read_embed_sequences_json(
            json.dumps(request), json.dumps(paths), output, str(item_file), json.dumps(part)
        ))

    assert objects["biologicalFeatures"] and objects["sequenceSources"]
    assert reply == {
        "biologicalFeatures": ",".join(_compact(row) for row in objects["biologicalFeatures"]),
        "rowCount": len(objects["biologicalFeatures"]),
        "sequenceSources": [
            {"index": entry["index"], "source": _compact(entry["source"])} for entry in objects["sequenceSources"]
        ],
        "sourceCount": 5,
        # A similarity group ribbon names no span.
        "matchSpanReads": None,
    }


def test_the_export_writes_what_each_match_span_reads_with_the_comparison_sources() -> None:
    """The export finds each span's record as the Worker's read does (OV-3542): a band that names no source reads its
    record from any source, and a span that names no record reads none, even when only one record could be read.
    Only the part that holds the comparison sources (the last) carries the reads; the whole item carries the same."""

    records = CatalogSequenceRecords(
        mode="circular",
        records=(SeqRecord(Seq("AACCGGTT"), id="ref", name="REF"),),
        unit_keys=("ref-key",),
        comparison_records=((), (SeqRecord(Seq("TTGCAACC"), id="cmp", name="CMP"),)),
    )
    band = {"match_kind": "homology", "reference_side": "query", "qstart": "1", "qend": "2", "sstart": "4", "send": "1"}
    item = {"recordKeys": ["ref-key"], "biologicalFeatures": [], "comparisonMatches": [
        {**band, "id": "any-source", "query_record_id": "ref", "subject_record_id": "cmp"},
        {**band, "id": "no-record", "source_index": "1"},
    ]}

    first = embed_catalog_part(item, records, source_record_keys=["ref-key"], comparison_sources=False)
    last = embed_catalog_part(item, records, source_record_keys=[], comparison_sources=True)
    whole = embed_catalog_item(item, records)

    assert first["matchSpanReads"] is None
    assert last["matchSpanReads"] == whole["matchSpanReads"] == {
        "any-source": {"query": {"source": 0}, "subject": {"source": 1}},
        "no-record": {
            "query": {"reason": "RECORD_UNKNOWN", "sourceLength": 0},
            "subject": {"reason": "RECORD_UNKNOWN", "sourceLength": 0},
        },
    }
    assert [source["key"] for source in whole["sequenceSources"]] == [  # type: ignore[attr-defined]
        "circular:record:0", "homology:comparison:1:cmp",
    ]
    assert read_match_span(whole["sequenceSources"][1]["sequence"], 4, 1) == "GCAA"  # type: ignore[index]


def test_whole_record_sources_equal_stored_sequence_sources(tmp_path: Path) -> None:
    path = SESSIONS / "BGC0000708-BGC0000713.v40-schema5.json"
    [item] = json.loads(path.read_text(encoding="utf-8"))["editorState"]["featureCatalog"]["items"]
    stored = {source["key"]: source["sequence"] for source in item["sequenceSources"]}

    with _committed(path, tmp_path) as (request, paths, output):
        spans = [
            {"kind": "record", "recordKey": key, "recordId": "", "start": 1, "end": len(stored[f"linear:record:{index}"])}
            for index, key in enumerate(item["recordKeys"])
            if f"linear:record:{index}" in stored
        ]
        read = read_match_sequences(request, paths, output, spans)["spans"]

    assert len(spans) == len(stored) == 5
    assert [span["sequence"] for span in read] == list(stored.values())


@pytest.mark.parametrize(
    "case",
    json.loads(SPAN_VECTORS.read_text(encoding="utf-8"))["vectors"],
    ids=lambda case: str(case["name"]),
)
def test_match_spans_equal_the_dev_registry_bases(case: dict[str, object]) -> None:
    """The bases of a span equal what dev's page read from the stored whole-record source (vectors from dev's JS)."""

    start, end = case["start"], case["end"]
    assert isinstance(start, int) and isinstance(end, int)
    # A whole-record source is upper case (``match_sequence_sources``), as the source dev's page read.
    assert read_match_span(str(case["sequence"]).upper(), start, end) == case["bases"]


def _cds(start: int, end: int, strand: int, **qualifiers: list[str]) -> SeqFeature:
    return SeqFeature(FeatureLocation(start, end, strand=strand), type="CDS", qualifiers=dict(qualifiers))


@contextmanager
def _rendered(
    tmp_path: Path, *records: SeqRecord
) -> Iterator[tuple[dict[str, object], dict[str, str], str, dict[str, object]]]:
    """A Linear request of ``records`` (record keys ``key-N``), rendered as the Web does: its request, resources,
    output directory and stored catalog item."""

    for record in records:
        record.annotations["molecule_type"] = "DNA"
    document = build_session_document(LinearDiagramRequest(
        records=tuple(
            RecordInput(source=InMemoryRecordSource(record), record_key=f"key-{index}")
            for index, record in enumerate(records)
        ),
        output=RenderOutputRequest(output_prefix="on-demand", formats=("svg",)),
    ))
    resources = tmp_path / "resources"
    resources.mkdir()
    paths: dict[str, str] = {}
    for resource_id, source in document.resources.sources.items():
        with source.open_bytes() as stream:
            (resources / resource_id).write_bytes(stream.read())
        paths[resource_id] = str(resources / resource_id)
    view = drawing_view(embedded_session(document.to_dict(), document.resources.sources))
    request = view["renderRequest"]
    reply = render_embedded_canonical_web_request(request, resources=view["resources"], workspace=tmp_path / "render")
    output = tmp_path / "outputs"
    output.mkdir()
    [item] = reply["featureCatalog"]["items"]
    assert reply["featureCatalog"]["schema"] == FEATURE_CATALOG_SCHEMA
    yield request, paths, str(output), item


def test_read_is_by_identity_not_by_caption(tmp_path: Path) -> None:
    # Two CDSs with one product and one gene name: only their identities tell them apart.
    record = SeqRecord(
        Seq("ATGAAATAAGGGATGCCCTAA"),
        id="caption-record",
        features=[
            _cds(0, 9, 1, gene=["same"], product=["same protein"]),
            _cds(12, 21, 1, gene=["same"], product=["same protein"]),
        ],
    )
    with _rendered(tmp_path, record) as (request, paths, output, item):
        rows = item["biologicalFeatures"]
        assert len(rows) == 2 and len({row["biologicalFeatureId"] for row in rows}) == 2
        assert not any(STORED_SEQUENCE_KEYS & set(row) for row in rows)
        read = read_feature_sequences(request, paths, output, [
            {"recordKey": row["recordKey"], "biologicalFeatureId": row["biologicalFeatureId"]}
            for row in reversed(rows)
        ])["features"]

    assert [row["nucleotideSequence"] for row in read] == ["ATGCCCTAA", "ATGAAATAA"]
    assert [row["aminoAcidSequence"] for row in read] == ["MP", "MK"]


def test_identical_features_bind_to_their_own_rows(tmp_path: Path) -> None:
    """Two features with the same type and location (``<hash>~n``) each bind to the drawn feature of their row."""

    record = SeqRecord(
        Seq("ATGAAATAAGGG"),
        id="duplicate-record",
        features=[_cds(0, 9, 1, note=["first"]), _cds(0, 9, 1, note=["second"])],
    )
    with _rendered(tmp_path, record) as (request, paths, output, item):
        rows = item["biologicalFeatures"]
        read = read_feature_sequences(request, paths, output, [
            {"recordKey": row["recordKey"], "biologicalFeatureId": row["biologicalFeatureId"]} for row in rows
        ])["features"]

    assert len({row["biologicalFeatureId"] for row in rows}) == 2
    assert [row["nucleotideSequence"] for row in read] == ["ATGAAATAA", "ATGAAATAA"]


def test_unknown_identity_and_record_key_are_refused_with_a_reason(tmp_path: Path) -> None:
    record = SeqRecord(Seq("ATGAAATAAGGG"), id="refusal-record", features=[_cds(0, 9, 1)])
    with _rendered(tmp_path, record) as (request, paths, output, item):
        [row] = item["biologicalFeatures"]
        with pytest.raises(ValidationError) as unknown:
            read_feature_sequences(
                request, paths, output, [{"recordKey": row["recordKey"], "biologicalFeatureId": "fmissing"}]
            )
        with pytest.raises(ValidationError) as other_record:
            read_feature_sequences(
                request, paths, output, [{"recordKey": "key-9", "biologicalFeatureId": row["biologicalFeatureId"]}]
            )
        for malformed in ([], [{"recordKey": "key-0"}], [dict(row, extra=1)], [
            {"recordKey": row["recordKey"], "biologicalFeatureId": row["biologicalFeatureId"]}
        ] * 2):
            with pytest.raises(ValidationError) as protocol:
                read_feature_sequences(request, paths, output, malformed)  # type: ignore[arg-type]
            assert protocol.value.diagnostic["code"] == "HELPER_PROTOCOL"

    assert unknown.value.diagnostic == {"code": "FEATURE_IDENTITY"}
    assert "(unresolved)" in str(unknown.value)
    # A record key outside the request is the binding owner's error (``resolve_feature_identities``).
    assert other_record.value.diagnostic["code"] == "FEATURE_IDENTITY"


def test_match_span_identity_and_reasons(tmp_path: Path) -> None:
    """A drawn record by its key (not its ID or place), the strand by the order of start and end, and a reason for a
    span that cannot be read: an unknown record, a span outside the record, a comparison never supplied."""

    query = SeqRecord(Seq("AACCGGTT"), id="same-id", features=[_cds(0, 3, 1)])
    subject = SeqRecord(Seq("TTGGCCAA"), id="same-id", features=[_cds(0, 3, 1)])
    with _rendered(tmp_path, query, subject) as (request, paths, output, _item):
        read = read_match_sequences(request, paths, output, [
            {"kind": "record", "recordKey": "key-1", "recordId": "", "start": 2, "end": 4},
            {"kind": "record", "recordKey": "key-0", "recordId": "", "start": 4, "end": 2},
            {"kind": "record", "recordKey": "key-7", "recordId": "", "start": 1, "end": 2},
            {"kind": "record", "recordKey": "key-0", "recordId": "", "start": 1, "end": 9},
            {"kind": "comparison", "sourceIndex": 0, "recordId": "same-id", "start": 1, "end": 2},
            {"kind": "record", "recordKey": "key-0", "recordId": "", "start": 0, "end": 3},
            {"kind": "record", "recordKey": "key-0", "recordId": "", "start": 3, "end": -2},
        ])["spans"]

    assert [(span["status"], span["reason"], span["orientation"], span["sequence"]) for span in read] == [
        ("ok", None, "+", "TGG"),
        ("ok", None, "-", "GGT"),
        ("unavailable", "SOURCE_UNAVAILABLE", "+", ""),
        ("unavailable", "OUT_OF_RANGE", "+", ""),
        ("unavailable", "NEVER_SUPPLIED", "+", ""),
        ("unavailable", "COORDINATE_BELOW_ONE", "+", ""),
        ("unavailable", "COORDINATE_BELOW_ONE", "-", ""),
    ]
    assert read[0]["recordId"] == "same-id" and read[0]["length"] == 3
    # The length of the record a span reads, for the reason a span past its end shows (base: "({length} bp)").
    assert [span["sourceLength"] for span in read] == [8, 8, 0, 8, 0, 8, 8]


def test_a_span_of_an_empty_source_and_a_coordinate_below_one_have_their_reasons() -> None:
    """The coordinate reasons in the order the match popup checked them before D5: not a whole number, below 1, past
    the record's end (a record without bases is no source, as the match registry held none)."""

    assert span_reason("ACGT", None, 2) == "COORDINATE_NOT_WHOLE"
    assert span_reason("", 0, None) == "COORDINATE_NOT_WHOLE"
    assert span_reason("ACGT", 0, 2) == "COORDINATE_BELOW_ONE"
    assert span_reason("", 0, 2) == "COORDINATE_BELOW_ONE"
    assert span_reason("ACGT", 2, 5) == "OUT_OF_RANGE"
    assert span_reason("ACGT", 4, 1) is None


def test_a_coordinate_that_is_not_whole_is_judged_after_the_record(tmp_path: Path) -> None:
    """The page sends a coordinate that is not a whole number as ``None``; Python judges it in the match popup's order
    before D5: the record (unknown, ambiguous ID or alias, never supplied) first, then "not whole", then below 1,
    an empty source and past the end. A whole number JSON writes in exponent form is still whole."""

    query = SeqRecord(Seq("AACCGGTT"), id="query-id", name="shared", features=[_cds(0, 3, 1)])
    subject = SeqRecord(Seq("TTGGCCAA"), id="subject-id", name="shared", features=[_cds(0, 3, 1)])
    with _rendered(tmp_path, query, subject) as (request, paths, output, _item):
        read = read_match_sequences(request, paths, output, [
            {"kind": "record", "recordKey": "key-7", "recordId": "", "start": None, "end": 4},
            {"kind": "reference", "recordId": "missing", "recordKey": None, "start": 1, "end": None},
            {"kind": "reference", "recordId": "shared", "recordKey": None, "start": None, "end": 4},
            {"kind": "comparison", "sourceIndex": 0, "recordId": "query-id", "start": None, "end": None},
            {"kind": "record", "recordKey": "key-0", "recordId": "", "start": None, "end": -2},
            {"kind": "reference", "recordId": "subject-id", "recordKey": None, "start": 2, "end": None},
            {"kind": "record", "recordKey": "key-0", "recordId": "", "start": 1, "end": 1e21},
        ])["spans"]
        with pytest.raises(ValidationError, match="start and an end"):
            read_match_sequences(request, paths, output, [{"kind": "record", "recordKey": "key-0", "recordId": "", "start": 1.5, "end": 2}])

    assert [(span["status"], span["reason"], span["sourceLength"]) for span in read] == [
        ("unavailable", "SOURCE_UNAVAILABLE", 0),
        ("unavailable", "RECORD_UNKNOWN", 0),
        ("unavailable", "RECORD_ALIAS_AMBIGUOUS", 0),
        ("unavailable", "NEVER_SUPPLIED", 0),
        ("unavailable", "COORDINATE_NOT_WHOLE", 8),
        ("unavailable", "COORDINATE_NOT_WHOLE", 8),
        ("unavailable", "OUT_OF_RANGE", 8),
    ]


def test_a_drawn_reference_named_by_an_alias_is_read_without_catalog_rows(tmp_path: Path) -> None:
    """A Circular comparison band names its reference by the record ID its BLAST hit carries; the record is found
    among the drawn records by ID, then by alias (name, accession, first word), as the match registry found it
    before D5, also when the record has no feature rows in the catalog."""

    reference = SeqRecord(Seq("AACCGGTTAA"), id="NC_000001.1", name="REFLOCUS", features=[])
    reference.annotations["accessions"] = ["NC_000001"]
    other = SeqRecord(Seq("GGGGCCCCAA"), id="NC_000002.1", name="OTHERLOCUS", features=[])
    with _rendered(tmp_path, reference, other) as (request, paths, output, item):
        assert item["biologicalFeatures"] == []
        read = read_match_sequences(request, paths, output, [
            {"kind": "reference", "recordId": "NC_000001.1", "recordKey": None, "start": 1, "end": 4},
            {"kind": "reference", "recordId": "NC_000001", "recordKey": None, "start": 4, "end": 1},
            {"kind": "reference", "recordId": "REFLOCUS", "recordKey": None, "start": 2, "end": 3},
            {"kind": "reference", "recordId": "OTHERLOCUS", "recordKey": None, "start": 1, "end": 2},
            {"kind": "reference", "recordId": "NC_000009", "recordKey": None, "start": 1, "end": 2},
        ])["spans"]

    assert [(span["status"], span["reason"], span["recordId"], span["sequence"]) for span in read] == [
        ("ok", None, "NC_000001.1", "AACC"),
        ("ok", None, "NC_000001.1", "GGTT"),
        ("ok", None, "NC_000001.1", "AC"),
        ("ok", None, "NC_000002.1", "GG"),
        ("unavailable", "RECORD_UNKNOWN", "", ""),
    ]


def test_a_reference_span_that_names_its_unit_reads_that_unit_only(tmp_path: Path) -> None:
    """A Circular band on a cropped or reversed reference carries the reference's record index: the record is
    searched in that unit only (another record's ID is not found), and with no ID the unit is read, as the match
    popup found it before D5 (base c75daa1de, r4fix/probe-base-reference.out)."""

    reference = SeqRecord(Seq("AACCGGTT"), id="reference", name="REFLOCUS", features=[])
    other = SeqRecord(Seq("GGGGCCCC"), id="other", name="OTHER", features=[])
    with _rendered(tmp_path, reference, other) as (request, paths, output, _item):
        read = read_match_sequences(request, paths, output, [
            {"kind": "reference", "recordId": record_id, "recordKey": key, "start": 1, "end": 2}
            for record_id, key in (
                ("reference", "key-0"), ("REFLOCUS", "key-0"), ("other", "key-0"), ("", "key-0"), ("", "key-1"),
                ("other", None),
            )
        ])["spans"]

    assert [span["reason"] or span["sequence"] for span in read] == ["AA", "AA", "RECORD_UNKNOWN", "AA", "GG", "GG"]


def test_a_drawn_reference_id_shared_by_two_records_is_ambiguous() -> None:
    records = CatalogSequenceRecords(
        mode="circular",
        records=(
            SeqRecord(Seq("AAAA"), id="same", name="one"),
            SeqRecord(Seq("CCCC"), id="same", name="two"),
            SeqRecord(Seq("GGGG"), id="third described", name="three"),
            SeqRecord(Seq("TTTT"), id="fourth", name="shared"),
            SeqRecord(Seq("ACAC"), id="fifth", name="shared"),
        ),
        unit_keys=("key-0", "key-1", "key-2", "key-3", "key-4"),
    )
    sources = registered_match_sources(match_sequence_sources(records))

    def found(record_id: str) -> tuple[object, object]:
        source, reason = match_span_source(records, sources, {
            "kind": "reference", "recordId": record_id, "recordKey": None, "start": 1, "end": 2,
        })
        return (None if source is None else source["sequence"]), reason

    assert found("same") == (None, "RECORD_AMBIGUOUS")
    # Two records with one alias and no shared ID: the alias is ambiguous, as the match popup said before D5.
    assert found("shared") == (None, "RECORD_ALIAS_AMBIGUOUS")
    assert found("two") == ("CCCC", None)
    assert found("third") == ("GGGG", None)
    assert found("") == (None, "RECORD_UNKNOWN")


def test_a_comparison_span_found_by_id_then_alias_with_its_own_reason() -> None:
    """A comparison record is found by ID, else alias, in its source or, when the match names no source, in every
    source; two records with one ID (one registry key) count once when they agree and not at all otherwise, and a
    record without bases is none. Only a source without a comparison FASTA (or no source at all, when the match names
    none) says its FASTA was not supplied; a source whose FASTA lacks the record gives the record's own reason, as a
    drawn record does: unknown, ambiguous ID or alias (OV-2755)."""

    def record(record_id: str, name: str, sequence: str) -> SeqRecord:
        return SeqRecord(Seq(sequence), id=record_id, name=name)

    records = CatalogSequenceRecords(
        mode="circular",
        records=(record("reference", "reference", "AAAA"),),
        unit_keys=("reference-key",),
        comparison_records=(
            (
                record("contig_1", "c1", "CCCC"), record("dup", "x", "GGGG"), record("dup", "y", "TTTT"),
                record("twin", "z", "AAAA"), record("twin", "z", "AAAA"), record("p", "shared", "ACAC"),
                record("q", "shared", "GTGT"), record("e", "e", ""),
            ),
            (record("only2", "o", "CCGG"),),
            (),
        ),
    )
    sources = registered_match_sources(match_sequence_sources(records))

    def found(record_id: str, source_index: int | None = 0) -> object:
        source, reason = match_span_source(records, sources, {
            "kind": "comparison", "sourceIndex": source_index, "recordId": record_id, "start": 1, "end": 2,
        })
        return reason if source is None else source["sequence"]

    assert [found(name) for name in ("contig_1", "c1", "twin", "z")] == ["CCCC", "CCCC", "AAAA", "AAAA"]
    for name in ("dup", "x", "y", "absent", "", "e", "only2"):
        assert found(name) == "RECORD_UNKNOWN", name
    assert found("shared") == "RECORD_ALIAS_AMBIGUOUS"
    assert [found("only2", None), found("c1", None), found("absent", None), found("dup", None)] == [
        "CCGG", "CCCC", "RECORD_UNKNOWN", "RECORD_UNKNOWN",
    ]
    # A source without a FASTA, and a source the request does not have, were never supplied.
    assert [found("c1", 2), found("absent", 2), found("c1", 5)] == ["NEVER_SUPPLIED"] * 3
    without_fasta = CatalogSequenceRecords(
        mode="circular", records=records.records, unit_keys=records.unit_keys, comparison_records=((), ()),
    )
    for source_index in (None, 0, 1):
        assert match_span_source(without_fasta, registered_match_sources(match_sequence_sources(without_fasta)), {
            "kind": "comparison", "sourceIndex": source_index, "recordId": "c1", "start": 1, "end": 2,
        }) == (None, "NEVER_SUPPLIED"), source_index


def test_a_drawn_unit_by_key_carries_the_record_the_match_names(tmp_path: Path) -> None:
    """A drawn unit by its record key reads when the match names no record, or the unit's ID or alias; another name
    is a key that conflicts with the match identity, and a key without a unit is an unavailable source, each before
    the coordinates (base c75daa1de, r4fix/probe-base.out). A drawn record always has bases."""

    query = SeqRecord(Seq("AACCGGTT"), id="same-id", name="LOCUS0", features=[_cds(0, 3, 1)])
    with _rendered(tmp_path, query) as (request, paths, output, _item):
        read = read_match_sequences(request, paths, output, [
            {"kind": "record", "recordKey": "key-0", "recordId": record_id, "start": start, "end": 2}
            for record_id, start in (("", 1), ("same-id", 1), ("LOCUS0", 1), ("Other", 1), ("Other", None))
        ] + [
            {"kind": "record", "recordKey": "key-7", "recordId": "X", "start": None, "end": 10},
        ])["spans"]

    assert [span["reason"] or span["sequence"] for span in read] == [
        "AA", "AA", "AA", "SOURCE_KEY_CONFLICT", "SOURCE_KEY_CONFLICT", "SOURCE_UNAVAILABLE",
    ]


def test_the_stored_item_never_gains_references_and_the_embed_writes_them(tmp_path: Path) -> None:
    """The stored item of a Generate holds no sequence; only the export encoder writes sequences and references."""

    query = SeqRecord(Seq("ATGAAATAAGGG"), id="query", features=[_cds(0, 9, 1, translation=["MK"])])
    with _rendered(tmp_path, query) as (request, paths, output, item):
        embedded = read_embed_catalog_item(request, paths, output, item)

    assert not any(STORED_SEQUENCE_KEYS & set(row) for row in item["biologicalFeatures"])
    assert item["biologicalFeatures"][0]["qualifiers"]["translation"] == ["MK"]
    [row] = embedded["biologicalFeatures"]
    assert row["nucleotide_sequence"] == "ATGAAATAA" and row["amino_acid_sequence"] == "MK"
    assert row["translationFromAminoAcidSequence"] is True and "translation" not in row.get("qualifiers", {})
    with pytest.raises(Exception, match="drawn records"):
        embed_catalog_item(item, CatalogSequenceRecords(mode="linear", records=(query,), unit_keys=("other",)))


def test_the_export_of_a_record_drawn_in_windows_embeds_the_sequence_the_popup_reads(tmp_path: Path) -> None:
    """A feature drawn only in a later window of its record is bound in that window, as the popup read binds it (the
    first unit of the record that has the feature): the export embeds every row's sequence, equal to the read."""

    from gbdraw.api import GenBankInputSource
    from gbdraw.io.regions import RecordWindow

    mito = Path(__file__).resolve().parent / "test_inputs" / "HmmtDNA.gbk"
    document = build_session_document(LinearDiagramRequest(
        records=(RecordInput(
            GenBankInputSource(mito), windows=(RecordWindow(100, 900), RecordWindow(5000, 7600)), record_key="mito",
        ),),
        output=RenderOutputRequest(output_prefix="windows", formats=("svg",)),
    ))
    resources = tmp_path / "resources"
    resources.mkdir()
    paths: dict[str, str] = {}
    for resource_id, source in document.resources.sources.items():
        with source.open_bytes() as stream:
            (resources / resource_id).write_bytes(stream.read())
        paths[resource_id] = str(resources / resource_id)
    view = drawing_view(embedded_session(document.to_dict(), document.resources.sources))
    request = view["renderRequest"]
    reply = render_embedded_canonical_web_request(request, resources=view["resources"], workspace=tmp_path / "render")
    output = tmp_path / "outputs"
    output.mkdir()
    [item] = reply["featureCatalog"]["items"]
    rows = item["biologicalFeatures"]
    assert len(item["recordKeys"]) == 2 and {row["recordKey"] for row in rows} == {"mito"}

    read = read_feature_sequences(request, paths, output, [
        {"recordKey": row["recordKey"], "biologicalFeatureId": row["biologicalFeatureId"]} for row in rows
    ])["features"]
    embedded = read_embed_catalog_item(request, paths, output, item)
    sources = list(embedded.get("sequenceSources") or [])
    differences = [
        row["biologicalFeatureId"]
        for row, got in zip(embedded["biologicalFeatures"], read, strict=True)
        if _stored_sequences(row, sources)
        != {key: got[key] for key in ("nucleotideSequence", "aminoAcidSequence", "sequenceWarnings")}  # type: ignore[literal-required]
    ]
    assert differences == []
    assert len(rows) == 25 and all(row["nucleotideSequence"] for row in read)


def test_a_match_source_record_without_sequence_stops_with_the_input_reason():
    # A record with no bases (a CONTIG line) gives the input's reason, as str(record.seq) did, not "invalid catalog".
    from gbdraw.web_support.error_adapter import serialize_web_error

    contig = SeqRecord(Seq(None, 10), id="con1", name="con1")
    for records in (
        CatalogSequenceRecords(mode="linear", records=(contig,), unit_keys=("linear:record:0",)),
        CatalogSequenceRecords(mode="circular", records=(), comparison_records=((contig,),)),
    ):
        with pytest.raises(ValidationError) as caught:
            match_sequence_sources(records)
        model = serialize_web_error(caught.value, operation="readMatchSequences", stage="render")
        assert (model["code"], model["context"]) == ("INPUT_UNREADABLE", {"reason": "SEQUENCE_UNDEFINED"})
