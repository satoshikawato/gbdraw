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
    match_sequence_sources,
    materialize_catalog_nucleotide_sequence,
    promote_legacy_feature_catalog,
    read_match_span,
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
    [count] = counts
    assert sorted(sources) == list(range(count))
    embedded: dict[str, object] = {**item, "biologicalFeatures": rows}
    if sources:
        embedded["sequenceSources"] = [sources[index] for index in range(count)]
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
    }


def test_whole_record_sources_equal_stored_sequence_sources(tmp_path: Path) -> None:
    path = SESSIONS / "BGC0000708-BGC0000713.v40-schema5.json"
    [item] = json.loads(path.read_text(encoding="utf-8"))["editorState"]["featureCatalog"]["items"]
    stored = {source["key"]: source["sequence"] for source in item["sequenceSources"]}

    with _committed(path, tmp_path) as (request, paths, output):
        spans = [
            {"kind": "record", "recordKey": key, "start": 1, "end": len(stored[f"linear:record:{index}"])}
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
            {"kind": "record", "recordKey": "key-1", "start": 2, "end": 4},
            {"kind": "record", "recordKey": "key-0", "start": 4, "end": 2},
            {"kind": "record", "recordKey": "key-7", "start": 1, "end": 2},
            {"kind": "record", "recordKey": "key-0", "start": 1, "end": 9},
            {"kind": "comparison", "sourceIndex": 0, "recordId": "same-id", "start": 1, "end": 2},
            {"kind": "record", "recordKey": "key-0", "start": 0, "end": 3},
            {"kind": "record", "recordKey": "key-0", "start": 3, "end": -2},
        ])["spans"]

    assert [(span["status"], span["reason"], span["orientation"], span["sequence"]) for span in read] == [
        ("ok", None, "+", "TGG"),
        ("ok", None, "-", "GGT"),
        ("unavailable", "RECORD_UNKNOWN", "+", ""),
        ("unavailable", "OUT_OF_RANGE", "+", ""),
        ("unavailable", "NEVER_SUPPLIED", "+", ""),
        ("unavailable", "COORDINATE_BELOW_ONE", "+", ""),
        ("unavailable", "COORDINATE_BELOW_ONE", "-", ""),
    ]
    assert read[0]["recordId"] == "same-id" and read[0]["length"] == 3
    # The length of the record a span reads, for the reason a span past its end shows (base: "({length} bp)").
    assert [span["sourceLength"] for span in read] == [8, 8, 0, 8, 0, 8, 8]


def test_a_span_of_an_empty_source_and_a_coordinate_below_one_have_their_reasons() -> None:
    """The coordinate reasons in the order the match popup checked them before D5: below 1, an empty source, past
    the record's end."""

    from gbdraw.web_support.catalog_sequences import _span_reason

    assert _span_reason("ACGT", 0, 2) == "COORDINATE_BELOW_ONE"
    assert _span_reason("", 0, 2) == "COORDINATE_BELOW_ONE"
    assert _span_reason("", 1, 2) == "SOURCE_EMPTY"
    assert _span_reason("ACGT", 2, 5) == "OUT_OF_RANGE"
    assert _span_reason("ACGT", 4, 1) is None


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
            {"kind": "reference", "recordId": "NC_000001.1", "start": 1, "end": 4},
            {"kind": "reference", "recordId": "NC_000001", "start": 4, "end": 1},
            {"kind": "reference", "recordId": "REFLOCUS", "start": 2, "end": 3},
            {"kind": "reference", "recordId": "OTHERLOCUS", "start": 1, "end": 2},
            {"kind": "reference", "recordId": "NC_000009", "start": 1, "end": 2},
        ])["spans"]

    assert [(span["status"], span["reason"], span["recordId"], span["sequence"]) for span in read] == [
        ("ok", None, "NC_000001.1", "AACC"),
        ("ok", None, "NC_000001.1", "GGTT"),
        ("ok", None, "NC_000001.1", "AC"),
        ("ok", None, "NC_000002.1", "GG"),
        ("unavailable", "RECORD_UNKNOWN", "", ""),
    ]


def test_a_drawn_reference_id_shared_by_two_records_is_ambiguous() -> None:
    from gbdraw.web_support.catalog_sequences import _match_source

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
    sources = match_sequence_sources(records)

    def found(record_id: str) -> tuple[object, object]:
        source, reason = _match_source(records, sources, {
            "kind": "reference", "recordId": record_id, "start": 1, "end": 2,
        })
        return (None if source is None else source["sequence"]), reason

    assert found("same") == (None, "RECORD_AMBIGUOUS")
    # Two records with one alias and no shared ID: the alias is ambiguous, as the match popup said before D5.
    assert found("shared") == (None, "RECORD_ALIAS_AMBIGUOUS")
    assert found("two") == ("CCCC", None)
    assert found("third") == ("GGGG", None)
    assert found("") == (None, "RECORD_UNKNOWN")


def test_comparison_records_are_found_by_id_then_alias_and_ambiguity_is_a_reason() -> None:
    from gbdraw.web_support.catalog_sequences import _match_source

    def record(record_id: str, name: str, sequence: str) -> SeqRecord:
        return SeqRecord(Seq(sequence), id=record_id, name=name)

    records = CatalogSequenceRecords(
        mode="circular",
        records=(record("reference", "reference", "AAAA"),),
        unit_keys=("reference-key",),
        comparison_records=(
            (record("contig_1 described", "c1", "CCCC"), record("dup", "x", "GGGG"), record("dup", "y", "TTTT")),
        ),
    )
    sources = match_sequence_sources(records)

    def found(record_id: str) -> tuple[object, object]:
        source, reason = _match_source(records, sources, {
            "kind": "comparison", "sourceIndex": 0, "recordId": record_id, "start": 1, "end": 2,
        })
        return (None if source is None else source["sequence"]), reason

    assert found("contig_1 described") == ("CCCC", None)
    assert found("contig_1") == ("CCCC", None)
    assert found("c1") == ("CCCC", None)
    assert found("dup") == (None, "RECORD_AMBIGUOUS")
    assert found("x") == ("GGGG", None)
    assert found("absent") == (None, "RECORD_UNKNOWN")


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
