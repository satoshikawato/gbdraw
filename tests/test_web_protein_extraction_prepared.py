"""The blastp protein extraction reads the GenBank records the last Generate parsed (OV-1555).

The Worker stages the record's source under the resource id and token of the
last Generate, so ``extract_cds_protein_fasta`` finds that Generate's parsed
records in the prepared-input cache: it parses no GenBank file, publishes,
prunes and evicts nothing, changes no cached record, and returns exactly what
its own reader returns (R4). On a miss it keeps its own reader.
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import pytest

from gbdraw.api import GenBankInputSource, LinearDiagramOptions, LinearDiagramRequest, RecordInput
from gbdraw.api.requests import RecordCardinality
from gbdraw.api.options import LinearOutputOptions
from gbdraw.api.prepared import PreparedBiologicalInputCache
from gbdraw.api.request_render import plan_request
from gbdraw.features.overrides import FeatureOverride
from gbdraw.session_request_codec import encode_canonical_request
from gbdraw.web_support.request_render import prepared_input_lookup
from tests.test_web_feature_override_table import FIXTURE
from tests.test_web_helper_prepared_lookup import (  # noqa: F401 - the helpers fixture
    _cache_snapshot,
    _count_parses,
    _encoded_resources,
    _Stage,
    helpers,
)

VISIBILITY_HEADER = "record_id\tfeature_type\tqualifier\tvalue\taction\n"


@pytest.fixture
def drawn(tmp_path: Path) -> dict[str, object]:
    """A Linear request drawing both records of one GenBank file, and its Worker staging."""

    source = tmp_path / "two.gb"
    source.write_text(FIXTURE.read_text(encoding="utf-8"), encoding="utf-8")
    request = LinearDiagramRequest(
        records=(RecordInput(GenBankInputSource(source), record_key="r", cardinality=RecordCardinality.ALL),),
        options=LinearDiagramOptions(output=LinearOutputOptions(legend="none")),
    )
    encoded = encode_canonical_request(request)
    resources = _encoded_resources(encoded)
    cds = next(
        entry.biological_feature_id
        for entry in plan_request(request).provenance[0].source_feature_catalog
        if entry.feature_type == "CDS"
    )
    return {
        "source": source,
        "payload": encoded.payload,
        "stage": _Stage(tmp_path / "worker", resources),
        "override": FeatureOverride("r1", cds, feature_visibility="off").to_mapping(),
        "visibility": tmp_path / "visibility.tsv",
    }


# (region, selector, reverse flag, record index, overrides?, visibility rules?)
CASES = [
    (None, None, "0", 0, False, False),
    (None, "#2", "0", 1, False, False),
    (None, "#1", "1", 0, False, False),
    ("#1:200-3000", None, "0", 0, False, False),
    ("TESTB:200-3000:rc", "#2", "0", 1, False, False),
    (None, "#1", "0", 0, True, False),
    (None, "#1", "0", 0, False, True),
    ("#1:100-3500", None, "1", 0, True, True),
]


def _arguments(drawn: dict[str, object], path: str, case: tuple[object, ...]) -> list[object]:
    region, selector, reverse, index, overrides, visibility = case
    table = drawn["visibility"]
    assert isinstance(table, Path)
    if visibility:
        table.write_text(VISIBILITY_HEADER + "*\tCDS\tlocus_tag\tTESTA_0003\toff\n", encoding="utf-8")
    return [
        path, "genbank", None, region, selector, reverse, index, f"r{int(str(index)) + 1}",
        str(table) if visibility else None,
        json.dumps([drawn["override"]]) if overrides else None,
    ]


def _parsed_records(cache: PreparedBiologicalInputCache) -> bytes:
    """The cached records and what each catalog captured; catalog entries are built on first read."""
    return pickle.dumps([
        (entry.value.records, tuple(catalog.captured for catalog in entry.value.source_features))
        for entry in cache._layers["parsed"].values()
    ])


@pytest.mark.parametrize("case", CASES)
def test_extraction_after_a_generate_reads_its_records_and_equals_the_own_reader(
    helpers: dict[str, object], drawn: dict[str, object], case: tuple[object, ...], monkeypatch: pytest.MonkeyPatch,
) -> None:
    extract = helpers["extract_cds_protein_fasta"]
    assert callable(extract)
    expected = json.loads(extract(*_arguments(drawn, str(drawn["source"]), case)))
    assert "error" not in expected, expected

    cache = helpers["_WEB_PREPARED_INPUT_CACHE"]
    assert isinstance(cache, PreparedBiologicalInputCache)
    stage = drawn["stage"]
    assert isinstance(stage, _Stage)
    stage.render(cache, drawn["payload"])  # type: ignore[arg-type]
    before, records_before = _cache_snapshot(cache), _parsed_records(cache)

    parses = _count_parses(monkeypatch)
    _workspace, paths = stage.workspace("helper")
    (path,) = paths.values()
    with prepared_input_lookup(cache, resource_paths=paths, resource_identities=stage.identities):
        first = json.loads(extract(*_arguments(drawn, path, case)))
        again = json.loads(extract(*_arguments(drawn, path, case)))
    assert parses == []
    assert first == again == expected
    # Read-only: no entry published, pruned, evicted or reordered, and no cached record changed.
    assert _cache_snapshot(cache) == before
    assert _parsed_records(cache) == records_before


def test_the_worker_protocol_passes_the_staged_resources_after_the_arguments(
    helpers: dict[str, object], drawn: dict[str, object], monkeypatch: pytest.MonkeyPatch,
) -> None:
    cache = helpers["_WEB_PREPARED_INPUT_CACHE"]
    assert isinstance(cache, PreparedBiologicalInputCache)
    stage = drawn["stage"]
    assert isinstance(stage, _Stage)
    stage.render(cache, drawn["payload"])  # type: ignore[arg-type]
    call = helpers["call_web_json_helper"]
    assert callable(call)
    parses = _count_parses(monkeypatch)
    _workspace, paths = stage.workspace("helper")
    (path,) = paths.values()
    staged = json.loads(call(
        "extract_cds_protein_fasta", *_arguments(drawn, path, CASES[1]), json.dumps(paths), json.dumps(stage.identities),
    ))
    assert parses == []
    # Sent as files (no record source in the drawing): the helper's own reader, one parse.
    sent = json.loads(call(
        "extract_cds_protein_fasta", *_arguments(drawn, str(drawn["source"]), CASES[1]), "{}", "{}",
    ))
    assert len(parses) == 1
    assert staged == sent
    assert "error" not in staged, staged


def test_a_cold_lookup_that_does_not_publish_parses_once_and_stores_nothing(
    helpers: dict[str, object], drawn: dict[str, object], monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Every helper lookup but the protein extraction's keeps a miss for the call only."""
    cache = helpers["_WEB_PREPARED_INPUT_CACHE"]
    assert isinstance(cache, PreparedBiologicalInputCache)
    stage = drawn["stage"]
    assert isinstance(stage, _Stage)
    parses = _count_parses(monkeypatch)
    _workspace, paths = stage.workspace("helper")
    (path,) = paths.values()
    extract = helpers["extract_cds_protein_fasta"]
    assert callable(extract)
    with prepared_input_lookup(cache, resource_paths=paths, resource_identities=stage.identities):
        result = json.loads(extract(*_arguments(drawn, path, CASES[0])))
    assert "error" not in result, result
    assert len(parses) == 1
    assert all(not entries for entries in cache._layers.values())


def test_extractions_before_the_first_generate_parse_the_upload_once_for_all_cards_and_the_render(
    helpers: dict[str, object], drawn: dict[str, object], monkeypatch: pytest.MonkeyPatch,
) -> None:
    """OV-1595: the page stages a new upload once, under the id the next Generate names.

    Each card's extraction runs through the Worker's helper entry on a cold
    cache. The first one parses the file and publishes the parsed source (and
    nothing else); the other cards and the Generate that follows read it.
    """

    extract = helpers["extract_cds_protein_fasta"]
    call = helpers["call_web_json_helper"]
    assert callable(extract) and callable(call)
    cards = [CASES[0], CASES[1], CASES[3]]
    expected = [json.loads(extract(*_arguments(drawn, str(drawn["source"]), case))) for case in cards]
    assert all("error" not in result for result in expected), expected

    cache = helpers["_WEB_PREPARED_INPUT_CACHE"]
    assert isinstance(cache, PreparedBiologicalInputCache)
    stage = drawn["stage"]
    assert isinstance(stage, _Stage)
    assert all(not entries for entries in cache._layers.values())
    parses = _count_parses(monkeypatch)
    results = []
    for case in cards:
        _workspace, paths = stage.workspace("helper")
        (path,) = paths.values()
        results.append(json.loads(call(
            "extract_cds_protein_fasta", *_arguments(drawn, path, case), json.dumps(paths), json.dumps(stage.identities),
        )))
    assert results == expected
    assert len(parses) == 1
    assert len(cache._layers["parsed"]) == 1
    assert not cache._layers["resolved"] and not cache._layers["render_stage"]

    rendered = stage.render(cache, drawn["payload"])  # type: ignore[arg-type]
    assert "error" not in rendered["reply"], rendered["reply"]
    assert len(parses) == 1
    assert rendered["metrics"].get("parsedSourceParseCount", 0) == 0  # type: ignore[attr-defined]


def test_an_extraction_lookup_does_not_publish_a_parsed_source_it_changed(
    helpers: dict[str, object], drawn: dict[str, object],
) -> None:
    """The extraction's lookup publishes what it built unless it changed it: the lookup's end checks it in full."""

    from gbdraw.api.record_planning import prepared_parsed_source
    from gbdraw.api.request_render import load_render_source_records
    from gbdraw.exceptions import ValidationError

    cache = helpers["_WEB_PREPARED_INPUT_CACHE"]
    assert isinstance(cache, PreparedBiologicalInputCache)
    stage = drawn["stage"]
    assert isinstance(stage, _Stage)
    _workspace, paths = stage.workspace("helper")
    (path,) = paths.values()
    with pytest.raises(ValidationError, match="ownership was mutated"):
        with prepared_input_lookup(
            cache, resource_paths=paths, resource_identities=stage.identities, publish_inputs=True,
        ):
            parsed = prepared_parsed_source(GenBankInputSource(Path(path)), load_render_source_records)
            assert parsed is not None
            parsed.records[0].features[1].qualifiers["note"] = ["changed by a helper"]
    assert all(not entries for entries in cache._layers.values())


def test_an_extraction_that_replies_with_an_error_still_publishes_the_parse(
    helpers: dict[str, object], drawn: dict[str, object], monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A selector the file does not have is the extraction's own error; the file it parsed is kept for the Generate."""

    cache = helpers["_WEB_PREPARED_INPUT_CACHE"]
    assert isinstance(cache, PreparedBiologicalInputCache)
    stage = drawn["stage"]
    assert isinstance(stage, _Stage)
    call = helpers["call_web_json_helper"]
    assert callable(call)
    parses = _count_parses(monkeypatch)
    _workspace, paths = stage.workspace("helper")
    (path,) = paths.values()
    arguments = _arguments(drawn, path, CASES[0])
    arguments[4] = "#9"
    reply = json.loads(call("extract_cds_protein_fasta", *arguments, json.dumps(paths), json.dumps(stage.identities)))
    assert "error" in reply, reply
    assert len(cache._layers["parsed"]) == 1
    rendered = stage.render(cache, drawn["payload"])  # type: ignore[arg-type]
    assert "error" not in rendered["reply"], rendered["reply"]
    assert rendered["metrics"].get("parsedSourceParseCount", 0) == 0  # type: ignore[attr-defined]
    assert len(parses) == 1


def test_an_extraction_of_a_file_that_does_not_load_publishes_nothing(
    helpers: dict[str, object], tmp_path: Path, drawn: dict[str, object],
) -> None:
    """A failed load is never kept: the helper's own reader reports the error, and the cache stays empty."""

    cache = helpers["_WEB_PREPARED_INPUT_CACHE"]
    assert isinstance(cache, PreparedBiologicalInputCache)
    stage = _Stage(tmp_path / "broken", {"record-1-genbank": b"LOCUS       broken\nnot a flat file\n"})
    call = helpers["call_web_json_helper"]
    assert callable(call)
    _workspace, paths = stage.workspace("helper")
    (path,) = paths.values()
    reply = json.loads(call(
        "extract_cds_protein_fasta", *_arguments(drawn, path, CASES[0]), json.dumps(paths), json.dumps(stage.identities),
    ))
    assert "error" in reply, reply
    assert all(not entries for entries in cache._layers.values())
