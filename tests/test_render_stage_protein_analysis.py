"""The protein comparison analysis is one render stage (D-14, M4/P2).

A repeat render with the same stored raw results reuses the parse, the
orthogroups and the collinearity blocks. The display marks the analysis writes
into the request's raw cache come back with the stage result, so the reply's
LOSAT entries are the same on a build, on a hit and without a transaction.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import gbdraw.analysis.protein_colinearity as protein_colinearity_module
from gbdraw.analysis.protein_colinearity import (
    LosatpCacheManager,
    build_pairwise_protein_blastp_comparisons,
    extract_web_stable_cds_proteins,
    parse_losat_fasta_ids,
    parse_losatp_outfmt6,
)
from gbdraw.api.options import LinearDiagramOptions, LosatSearchOptions
from gbdraw.api.prepared import (
    PreparedBiologicalInputCache,
    PreparedResourceIdentity,
    register_prepared_resource_value,
)
from gbdraw.api.request_render import CurrentRequestArtifacts, build_request_diagram
from gbdraw.api.requests import GenBankInputSource, LinearDiagramRequest, RecordInput

_PROTEINS = ("MKLV" * 10, "MATG" * 12, "MSRE" * 9)


def _record(index: int) -> SeqRecord:
    sequence = "ATG" + "GCT" * 400 + "TAA"
    record = SeqRecord(Seq(sequence), id=f"rec{index}", name=f"rec{index}", description="")
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(
            FeatureLocation(start, start + 120, strand=1 if (start // 300) % 2 == 0 else -1),
            type="CDS",
            qualifiers={
                "translation": [protein],
                "protein_id": [f"rec{index}_p{number}"],
                "product": [f"protein {number}"],
            },
        )
        for number, (start, protein) in enumerate(zip((0, 300, 600), _PROTEINS))
    ]
    return record


def _request(mode: str, directory: Path) -> tuple[LinearDiagramRequest, dict]:
    """Three GenBank inputs, with the resource identities a Worker gives them."""

    paths = []
    for index in range(3):
        path = directory / f"rec{index}.gbk"
        SeqIO.write(_record(index), path, "genbank")
        paths.append(path)
    request = LinearDiagramRequest(
        records=tuple(
            RecordInput(source=GenBankInputSource(path), record_key=f"record-{index}")
            for index, path in enumerate(paths)
        ),
        # Collinear group inference on, as the recorded identities were built
        # before its default became off (OV-2303).
        options=LinearDiagramOptions(
            losat_search=LosatSearchOptions(program="losatp", losatp_mode=mode),
            collinear_infer_orthogroups=True,
        ),
    )
    identities = {
        path: PreparedResourceIdentity(f"rec{index}", f"token-{index}", path.stat().st_size)
        for index, path in enumerate(paths)
    }
    return request, identities


def _fake_losatp(query_fasta: str, subject_fasta: str, **kwargs) -> object:
    """Every query protein hits the subject protein at the same position (rank by header order)."""

    rows = [
        f"{query}\t{subject}\t{90 - rank}\t40\t0\t0\t1\t40\t1\t40\t1e-{20 + rank}\t{200 - rank}\n"
        for rank, (query, subject) in enumerate(
            zip(parse_losat_fasta_ids(query_fasta), parse_losat_fasta_ids(subject_fasta))
        )
    ]
    text = "".join(rows)
    kwargs["raw_output_callback"](text)
    kwargs["runtime_callback"]({"kind": "losat", "version": "test", "program": "blastp"})
    return parse_losatp_outfmt6(text)


def _no_search(*_args, **_kwargs):
    raise AssertionError("stored raw results answer every pair")


def _render(request, artifacts, cache=None, identities=None):
    if cache is None:
        return build_request_diagram(request, artifacts=artifacts), {}
    diagnostics: dict[str, object] = {"metrics": {}}
    with cache.transaction(resource_paths=identities, diagnostics=diagnostics):
        prepared = build_request_diagram(request, artifacts=artifacts)
    metrics = diagnostics["metrics"]
    assert isinstance(metrics, dict)
    return prepared, metrics


def _reply(prepared) -> tuple[object, ...]:
    return (
        prepared.drawing.tostring(),
        prepared.losat_cache_entries,
        prepared.losat_derived_cache_entries,
    )


@pytest.mark.parametrize("mode", ["pairwise", "similarity_groups", "collinear"])
def test_a_repeat_render_reuses_the_protein_analysis_and_its_display_marks(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    mode: str,
) -> None:
    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", _fake_losatp)
    request, identities = _request(mode, tmp_path)
    searched, _ = _render(request, None)
    raw = tuple(
        entry for entry in searched.losat_cache_entries if entry.get("kind") == "raw-losat"
    )
    assert raw and any(entry["display"] for entry in raw)
    # The stored entries come back undisplayed, so only the analysis marks them.
    artifacts = CurrentRequestArtifacts(
        losat_cache_entries=tuple({**entry, "display": False, "filename": ""} for entry in raw),
        protein_identity_manifest=searched.protein_identity_manifest,
    )
    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", _no_search)

    expected, _ = _render(request, artifacts)
    cache = PreparedBiologicalInputCache()
    cold, cold_metrics = _render(request, artifacts, cache, identities)
    warm, warm_metrics = _render(request, artifacts, cache, identities)

    assert cold_metrics["renderStageBuild:protein_analysis"] == 1
    assert warm_metrics["renderStageHit:protein_analysis"] == 1
    assert "renderStageBuild:protein_analysis" not in warm_metrics
    # The analysis gives its comparison frames an identity, so the stages that
    # draw them are reused too, and nothing in the warm render is uncacheable.
    assert warm_metrics["renderStageHit:comparison_rows"] == 1
    assert warm_metrics["renderStageHit:comparison_content"] >= 1
    assert not {name for name in warm_metrics if name.startswith("renderStageBuild:")}
    for metrics in (cold_metrics, warm_metrics):
        assert metrics.get("renderStageUncacheableCount", 0) == 0
        assert metrics.get("renderStageUndeclaredReadCount", 0) == 0
    assert [entry["display"] for entry in warm.losat_cache_entries if entry.get("kind") == "raw-losat"] == [
        entry["display"] for entry in expected.losat_cache_entries if entry.get("kind") == "raw-losat"
    ]
    assert any(entry["display"] for entry in expected.losat_cache_entries)
    assert _reply(cold) == _reply(expected)
    assert _reply(warm) == _reply(expected)


def test_a_stored_collinearity_result_reuses_its_projection_and_the_comparison_stages(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A Session's stored blocks are projected by one stage keyed by their resource."""

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", _fake_losatp)
    request, identities = _request("collinear", tmp_path)
    searched, _ = _render(request, None)
    assert searched.linear_metadata is not None
    stored = searched.linear_metadata.collinearity_result
    assert stored is not None
    given = replace(request, options=replace(request.options, losat_search=None, collinearity_blocks=stored))
    blocks_path = tmp_path / "blocks.json"
    blocks_path.write_text("stored blocks", encoding="utf-8")
    identities = {**identities, blocks_path: PreparedResourceIdentity("blocks", "token-blocks", 13)}
    cache = PreparedBiologicalInputCache()

    def render() -> tuple[object, dict]:
        diagnostics: dict[str, object] = {"metrics": {}}
        with cache.transaction(resource_paths=identities, diagnostics=diagnostics):
            # As the Session decoder registers a decoded collinearity result.
            register_prepared_resource_value(stored, blocks_path)
            prepared = build_request_diagram(given)
        metrics = diagnostics["metrics"]
        assert isinstance(metrics, dict)
        return prepared, metrics

    cold, cold_metrics = render()
    warm, warm_metrics = render()

    assert cold_metrics["renderStageBuild:collinearity_projection"] == 1
    assert warm_metrics["renderStageHit:collinearity_projection"] == 1
    assert warm_metrics["renderStageHit:comparison_rows"] == 1
    assert warm_metrics["renderStageHit:comparison_content"] >= 1
    assert not {name for name in warm_metrics if name.startswith("renderStageBuild:")}
    for metrics in (cold_metrics, warm_metrics):
        assert metrics.get("renderStageUncacheableCount", 0) == 0
    expected = build_request_diagram(given).drawing.tostring()
    assert cold.drawing.tostring() == expected
    assert warm.drawing.tostring() == expected


def _records_and_extraction():
    records = [_record(index) for index in range(3)]
    extraction = extract_web_stable_cds_proteins(
        records, record_instance_keys=("r0", "r1", "r2")
    )
    return records, extraction


def test_writes_replayed_on_an_equal_cache_give_the_same_session_entries(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A fork records what a search and the stored lookups write; replay is exact."""

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", _fake_losatp)
    records, extraction = _records_and_extraction()
    names = ("a.b.losatp.tsv", "b.c.losatp.tsv")
    seeded = LosatpCacheManager(identity_manifest=extraction.identity_manifest)
    build_pairwise_protein_blastp_comparisons(
        records, losatp_cache=seeded, protein_extraction=extraction, cache_filenames=names
    )
    (first, *rest) = seeded.session_entries()
    # One stored entry under another display name, the others missing: the
    # analysis marks the first under its new name and searches the rest.
    loaded = ({**first, "filename": "old.tsv", "display": True},)

    direct = LosatpCacheManager(loaded, identity_manifest=extraction.identity_manifest)
    build_pairwise_protein_blastp_comparisons(
        records, losatp_cache=direct, protein_extraction=extraction, cache_filenames=names
    )

    caller = LosatpCacheManager(loaded, identity_manifest=extraction.identity_manifest)
    fork = LosatpCacheManager.from_state(caller.state())
    build_pairwise_protein_blastp_comparisons(
        records, losatp_cache=fork, protein_extraction=extraction, cache_filenames=names
    )
    assert caller.session_entries() == LosatpCacheManager(
        loaded, identity_manifest=extraction.identity_manifest
    ).session_entries()
    caller.apply_writes(fork.writes)

    assert {write.kind for write in fork.writes} == {"entry", "display"}
    assert caller.session_entries() == direct.session_entries()
    assert [entry["filename"] for entry in caller.session_entries()] == list(names)
    assert len(rest) == 1
