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
    return pickle.dumps([entry.value for entry in cache._layers["parsed"].values()])


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


def test_a_cold_extraction_parses_once_and_stores_nothing(
    helpers: dict[str, object], drawn: dict[str, object], monkeypatch: pytest.MonkeyPatch,
) -> None:
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
