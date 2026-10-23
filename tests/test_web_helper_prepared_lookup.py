"""Worker helpers that plan a staged request read the prepared-input cache.

Align / Apply / Reset (``resolve_similarity_alignment_json`` with a
projection), Load Feature Edits TSV and the stored protein frames helper stage
the request's resources through the same transport as Generate. They read the
sources the last Generate parsed under the same resource identities, parse a
source at most once per call when the cache is cold, and leave the cache as
the Generate left it (OV-1550). An unchanged comparison table is parsed once
per Worker resource, not once per render (OV-1557).
"""

from __future__ import annotations

import io
import json
import shutil
from pathlib import Path
from typing import Callable

import pandas as pd
import pytest
from Bio import SeqIO

import gbdraw.analysis.conservation as conservation_module
import gbdraw.api.record_planning as record_planning_module
import gbdraw.diagrams.linear.assemble as assemble_module
import gbdraw.io.comparisons as comparisons_module
import gbdraw.session_request_codec as codec_module
from gbdraw.api import (
    CircularDiagramOptions,
    CircularDiagramRequest,
    GenBankInputSource,
    LinearDiagramOptions,
    LinearDiagramRequest,
    RecordInput,
)
from gbdraw.api.options import LinearOutputOptions
from gbdraw.api.prepared import PreparedBiologicalInputCache, prepared_resource_identity
from gbdraw.session_request_codec import encode_canonical_request
from gbdraw.web_support.request_render import (
    _render_staged_canonical_web_request_with_prepared_inputs as render_staged,
)
from tests.test_alignment_direction_projection import build_fixture
from tests.test_web_feature_override_table import FIXTURE as OVERRIDE_FIXTURE
from tests.test_web_feature_override_table import HEADER as OVERRIDE_HEADER

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def helpers() -> dict[str, object]:
    """The Worker's Python helpers, each test with its own prepared-input cache."""

    source = (ROOT / "gbdraw/web/js/app/python-helpers.js").read_text(encoding="utf-8")
    namespace: dict[str, object] = {}
    exec(source.split("`", 1)[1].rsplit("`", 1)[0].replace("\\\\", "\\"), namespace)
    return namespace


class _Stage:
    """Stage resources as the Worker does: one cache file per token, one symlink per request."""

    def __init__(self, root: Path, resources: dict[str, Path | bytes]) -> None:
        self.root = root
        self.cache_directory = root / "gbdraw-web-render-resource-cache"
        self.cache_directory.mkdir(parents=True)
        self.identities: dict[str, dict[str, object]] = {}
        self.cache_paths: dict[str, Path] = {}
        for index, (resource_id, source_path) in enumerate(sorted(resources.items()), start=1):
            token = f"render-resource-{index}"
            target = self.cache_directory / f"{token}.bin"
            if isinstance(source_path, bytes):
                target.write_bytes(source_path)
            else:
                shutil.copyfile(source_path, target)
            self.cache_paths[resource_id] = target
            self.identities[resource_id] = {"cacheToken": token, "size": target.stat().st_size}
        self.requests = 0

    def workspace(self, kind: str) -> tuple[Path, dict[str, str]]:
        self.requests += 1
        workspace = self.root / f"gbdraw-web-{kind}-{self.requests}"
        resource_directory = workspace / "resources"
        resource_directory.mkdir(parents=True)
        (workspace / ".gbdraw-worker-render-workspace").touch()
        paths: dict[str, str] = {}
        for index, (resource_id, target) in enumerate(sorted(self.cache_paths.items()), start=1):
            link = resource_directory / f"{index:04d}.bin"
            link.symlink_to(target)
            paths[resource_id] = str(link)
        return workspace, paths

    def render(self, cache: PreparedBiologicalInputCache, payload: dict[str, object]) -> dict[str, object]:
        workspace, paths = self.workspace("render")
        diagnostics: dict[str, object] = {"timingsMs": {}, "metrics": {}}
        reply = render_staged(
            payload,
            resource_paths=paths,
            workspace=workspace,
            _diagnostics=diagnostics,
            _prepared_input_cache=cache,
            _resource_identities=self.identities,
        )
        return {"reply": reply, "metrics": diagnostics["metrics"]}


def _encoded_resources(encoded: object) -> dict[str, Path | bytes]:
    return {
        item.resource_id: item.content if item.source_path is None else Path(item.source_path)
        for item in encoded.resources  # type: ignore[attr-defined]
    }


def _cache_snapshot(cache: PreparedBiologicalInputCache) -> dict[str, list[tuple[object, int]]]:
    return {layer: [(key, id(entry)) for key, entry in entries.items()] for layer, entries in cache._layers.items()}


def _count_parses(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    calls: list[str] = []
    original = SeqIO.parse

    def counting(handle, *args, **kwargs):  # type: ignore[no-untyped-def]
        calls.append(str(handle))
        return original(handle, *args, **kwargs)

    monkeypatch.setattr(SeqIO, "parse", counting)
    return calls


@pytest.fixture
def alignment(tmp_path: Path) -> dict[str, object]:
    (tmp_path / "sources").mkdir()
    _request, helper_request, context, paths_by_id, _response, _paths = build_fixture(tmp_path / "sources")
    return {
        "payload": context["canonicalRequest"],
        "helper": helper_request,
        "context": context,
        "stage": _Stage(tmp_path / "worker", {key: Path(path) for key, path in paths_by_id.items()}),
        "sources": len({path for path in paths_by_id.values()}),
    }


def _align(helpers: dict[str, object], alignment: dict[str, object]) -> dict[str, object]:
    stage = alignment["stage"]
    assert isinstance(stage, _Stage)
    workspace, paths = stage.workspace("helper")
    call = helpers["call_web_json_helper"]
    assert callable(call)
    result = json.loads(call(
        "resolve_similarity_alignment_json",
        json.dumps(alignment["helper"]),
        json.dumps(alignment["context"]),
        json.dumps(paths),
        str(workspace / "output"),
        json.dumps(stage.identities),
    ))
    assert "error" not in result, result
    return result


def test_alignment_projection_after_a_generate_parses_no_source_and_keeps_the_cache(
    helpers: dict[str, object], alignment: dict[str, object], monkeypatch: pytest.MonkeyPatch,
) -> None:
    cache = helpers["_WEB_PREPARED_INPUT_CACHE"]
    assert isinstance(cache, PreparedBiologicalInputCache)
    stage = alignment["stage"]
    assert isinstance(stage, _Stage)
    first = stage.render(cache, alignment["payload"])
    assert first["metrics"]["parsedSourceParseCount"] == alignment["sources"]
    before = _cache_snapshot(cache)
    assert before["render_stage"], "the Generate stored its render stages"

    parses = _count_parses(monkeypatch)
    _align(helpers, alignment)
    assert parses == []
    # Read-only: no entry published, pruned, evicted or reordered.
    assert _cache_snapshot(cache) == before
    assert cache._previous_transaction_failed is False

    # The next identical Generate reuses everything the first one stored.
    second = stage.render(cache, alignment["payload"])
    assert second["metrics"]["parsedSourceParseCount"] == 0
    assert second["metrics"]["renderStageBuildCount"] == 0
    assert parses == []


def test_a_cold_alignment_projection_parses_each_source_once_and_stores_nothing(
    helpers: dict[str, object], alignment: dict[str, object], monkeypatch: pytest.MonkeyPatch,
) -> None:
    cache = helpers["_WEB_PREPARED_INPUT_CACHE"]
    assert isinstance(cache, PreparedBiologicalInputCache)
    parses = _count_parses(monkeypatch)
    _align(helpers, alignment)
    # The plans (before, after, flipped) share one parse of each source.
    assert len(parses) == alignment["sources"]
    assert all(not entries for entries in cache._layers.values())


def test_an_alignment_plans_before_and_after_and_only_resolves_the_flipped_records(
    helpers: dict[str, object], alignment: dict[str, object], monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The flipped records give anchor facts only; they are never planned (OV-1585).

    With no orientations (the stale-target repair) the requested records are
    the flipped ones, so the after plan supplies them.
    """

    import gbdraw.web_support.similarity_alignment as alignment_module

    plans: list[object] = []
    original = alignment_module.plan_linear_request

    def counting(request):  # type: ignore[no-untyped-def]
        plans.append(request)
        return original(request)

    monkeypatch.setattr(alignment_module, "plan_linear_request", counting)
    reset = _align(helpers, alignment)
    assert len(plans) == 2
    # The Align click sends the current orientations: still the before and after plans only.
    plans.clear()
    current = {
        key: not value for key, value in reset["projection"]["geometryOrientations"].items()
    }
    _align(helpers, {**alignment, "context": {**alignment["context"], "orientations": current}})  # type: ignore[dict-item]
    assert len(plans) == 2


def test_alignment_identities_must_match_the_staged_manifest(
    helpers: dict[str, object], alignment: dict[str, object],
) -> None:
    stage = alignment["stage"]
    assert isinstance(stage, _Stage)
    workspace, paths = stage.workspace("helper")
    wrong = {key: {**value, "size": int(value["size"]) + 1} for key, value in stage.identities.items()}
    call = helpers["call_web_json_helper"]
    assert callable(call)
    result = json.loads(call(
        "resolve_similarity_alignment_json", json.dumps(alignment["helper"]),
        json.dumps(alignment["context"]), json.dumps(paths), str(workspace / "output"), json.dumps(wrong),
    ))
    assert result["error"]["operation"] == "resolveSimilarityAlignment"


def test_feature_edits_table_after_a_generate_parses_no_source(
    helpers: dict[str, object], tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    texts = [f"{chunk.strip()}\n//\n" for chunk in OVERRIDE_FIXTURE.read_text().split("\n//") if chunk.strip()]
    sources = []
    for name, text in zip(("testa.gb", "testb.gb"), texts, strict=True):
        sources.append(tmp_path / name)
        sources[-1].write_text(text)
    request = LinearDiagramRequest(
        records=tuple(RecordInput(GenBankInputSource(path), record_key=f"r{index}") for index, path in enumerate(sources)),
        options=LinearDiagramOptions(output=LinearOutputOptions(legend="none")),
    )
    encoded = encode_canonical_request(request)
    stage = _Stage(tmp_path / "worker", _encoded_resources(encoded))
    cache = helpers["_WEB_PREPARED_INPUT_CACHE"]
    assert isinstance(cache, PreparedBiologicalInputCache)
    stage.render(cache, encoded.payload)
    before = _cache_snapshot(cache)
    table = tmp_path / "overrides.tsv"
    table.write_text(OVERRIDE_HEADER + "#1\tlocus_tag=TESTA_0003\toff\t\t\n", encoding="utf-8")

    parses = _count_parses(monkeypatch)
    workspace, paths = stage.workspace("helper")
    call = helpers["call_web_json_helper"]
    assert callable(call)
    result = json.loads(call(
        "read_feature_override_table_json", str(table), json.dumps(encoded.payload),
        json.dumps(paths), str(workspace / "output"), json.dumps(stage.identities),
    ))
    assert "error" not in result, result
    assert len(result["rows"]) == 1
    assert parses == []
    assert _cache_snapshot(cache) == before


def test_stored_protein_frames_after_a_generate_parse_no_source(
    helpers: dict[str, object], alignment: dict[str, object], monkeypatch: pytest.MonkeyPatch,
) -> None:
    cache = helpers["_WEB_PREPARED_INPUT_CACHE"]
    assert isinstance(cache, PreparedBiologicalInputCache)
    stage = alignment["stage"]
    assert isinstance(stage, _Stage)
    stage.render(cache, alignment["payload"])
    before = _cache_snapshot(cache)

    parses = _count_parses(monkeypatch)
    workspace, paths = stage.workspace("helper")
    call = helpers["call_web_json_helper"]
    assert callable(call)
    from gbdraw.api.session_compat import DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION

    result = json.loads(call(
        "stored_protein_results_to_search_frame_json",
        json.dumps({"version": DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION, "renderRequest": alignment["payload"]}),
        json.dumps(paths), str(workspace / "output"), json.dumps(stage.identities),
    ))
    assert result == {"resources": {}}
    assert parses == []
    assert _cache_snapshot(cache) == before


# OV-1557: comparison tables


def _count_table_reads(monkeypatch: pytest.MonkeyPatch) -> list[pd.DataFrame]:
    """Each comparison table parse, as a copy of the frame it returned."""

    calls: list[pd.DataFrame] = []
    original = comparisons_module.read_comparison_table

    def counting(source, *args, **kwargs):  # type: ignore[no-untyped-def]
        frame = original(source, *args, **kwargs)
        calls.append(frame.copy(deep=True))
        return frame

    for module in (comparisons_module, conservation_module, assemble_module, record_planning_module, codec_module):
        if hasattr(module, "read_comparison_table"):
            monkeypatch.setattr(module, "read_comparison_table", counting)
    # A precomputed protein comparison is a canonical TSV.
    canonical = codec_module._read_canonical_table

    def counting_canonical(source, *args, **kwargs):  # type: ignore[no-untyped-def]
        frame = canonical(source, *args, **kwargs)
        calls.append(frame.copy(deep=True))
        return frame

    monkeypatch.setattr(codec_module, "_read_canonical_table", counting_canonical)
    return calls


def _two_record_sources(tmp_path: Path) -> tuple[Path, Path, Path]:
    texts = [f"{chunk.strip()}\n//\n" for chunk in OVERRIDE_FIXTURE.read_text().split("\n//") if chunk.strip()]
    sources = []
    for name, text in zip(("testa.gb", "testb.gb"), texts, strict=True):
        sources.append(tmp_path / name)
        sources[-1].write_text(text)
    records = [SeqIO.read(path, "genbank") for path in sources]
    rows = "".join(
        f"{records[0].id}\t{records[1].id}\t{90 + index % 9}.0\t{200}\t3\t0\t{1 + index * 10}\t{200 + index * 10}"
        f"\t{1 + index * 12}\t{200 + index * 12}\t1e-{20 + index % 30}\t{300 + index}.5\n"
        for index in range(200)
    )
    table = tmp_path / "pair.tsv"
    table.write_text(rows, encoding="utf-8")
    return sources[0], sources[1], table


def _linear_comparison_request(tmp_path: Path) -> object:
    first, second, table = _two_record_sources(tmp_path)
    return LinearDiagramRequest(
        records=(RecordInput(GenBankInputSource(first), record_key="a"), RecordInput(GenBankInputSource(second), record_key="b")),
        options=LinearDiagramOptions(blast_files=(str(table),), output=LinearOutputOptions(legend="none")),
    )


def _circular_ring_request(tmp_path: Path) -> object:
    first, _second, table = _two_record_sources(tmp_path)
    return CircularDiagramRequest(
        records=(RecordInput(GenBankInputSource(first)),),
        options=CircularDiagramOptions(conservation_blast_files=(str(table),)),
    )


def _linear_protein_request(tmp_path: Path) -> object:
    first, second, table = _two_record_sources(tmp_path)
    return LinearDiagramRequest(
        records=(RecordInput(GenBankInputSource(first), record_key="a"), RecordInput(GenBankInputSource(second), record_key="b")),
        options=LinearDiagramOptions(
            protein_comparisons=(comparisons_module.read_comparison_table(table),),
            output=LinearOutputOptions(legend="none"),
        ),
    )


@pytest.mark.parametrize(
    "build",
    [_linear_comparison_request, _circular_ring_request, _linear_protein_request],
    ids=["linear", "circular-ring", "linear-protein"],
)
def test_a_second_identical_generate_reads_no_comparison_table(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, build: Callable[[Path], object],
) -> None:
    encoded = encode_canonical_request(build(tmp_path))
    stage = _Stage(tmp_path / "worker", _encoded_resources(encoded))
    cache = PreparedBiologicalInputCache()
    reads = _count_table_reads(monkeypatch)
    first = stage.render(cache, encoded.payload)
    second = stage.render(cache, encoded.payload)
    assert len(reads) == 1
    assert second["metrics"]["decodedResourceCacheHitCount"] >= 1
    # After two renders the shared frame is still the frame the parse
    # returned (no reader changed it), and the drawing is the same.
    (stored,) = [
        entry.value if isinstance(entry.value, pd.DataFrame) else entry.value[1]
        for entry in cache._layers["decoded"].values()
        if isinstance(entry.value, (pd.DataFrame, tuple))
    ]
    pd.testing.assert_frame_equal(stored, reads[0])
    assert second["reply"]["results"][0]["content"] == first["reply"]["results"][0]["content"]


def test_a_changed_comparison_table_token_is_read_again(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    encoded = encode_canonical_request(_linear_comparison_request(tmp_path))
    stage = _Stage(tmp_path / "worker", _encoded_resources(encoded))
    cache = PreparedBiologicalInputCache()
    reads = _count_table_reads(monkeypatch)
    stage.render(cache, encoded.payload)
    table_id = next(item.resource_id for item in encoded.resources if item.source_path and item.source_path.name == "pair.tsv")
    stage.identities[table_id] = {**stage.identities[table_id], "cacheToken": "render-resource-99"}
    stage.render(cache, encoded.payload)
    assert len(reads) == 2


def test_comparison_resource_reader_outside_a_cache_reads_the_file(tmp_path: Path) -> None:
    _first, _second, table = _two_record_sources(tmp_path)
    assert prepared_resource_identity(table) is None
    frame = comparisons_module.read_comparison_resource(table)
    pd.testing.assert_frame_equal(frame, comparisons_module.read_comparison_table(io.StringIO(table.read_text())))
