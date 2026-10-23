from __future__ import annotations

import base64
import copy
import hashlib
import json
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
import pandas as pd
import pytest

from gbdraw.analysis.protein_colinearity import OrthogroupMember, OrthogroupResult
from gbdraw.api import (
    CircularBatchRequest,
    CircularDiagramOptions,
    CircularDiagramRequest,
    CircularOutputOptions,
    GenBankInputSource,
    GffFastaInputSource,
    LinearDiagramOptions,
    LinearDiagramRequest,
    RecordCardinality,
    RecordInput,
    RecordPresentation,
    RenderOutputRequest,
)
from gbdraw.api.prepared import (
    PreparedBiologicalInputCache,
    PreparedResourceIdentity,
    get_or_build_decoded_resource,
    get_or_build_parsed_source,
    register_prepared_value_identity,
)
from gbdraw.render_stages import render_stage
from gbdraw.annotations import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    RegionAnnotation,
)
from gbdraw.exceptions import ValidationError
import gbdraw.api.prepared as prepared_module
import gbdraw.api.request_render as request_render_module
from gbdraw.api.request_render import normalize_request_records, plan_request
from gbdraw.features.overrides import FeatureOverride
from gbdraw.io.regions import parse_region_spec
from gbdraw.io.record_select import parse_record_selector
from gbdraw.session import build_session_document
from gbdraw.web_support.request_render import (
    _render_staged_canonical_web_request_with_prepared_inputs as
    render_staged_canonical_web_request,
)


def _record(record_id: str, *, product: str = "cache target") -> SeqRecord:
    record = SeqRecord(
        Seq("ATG" * 240),
        id=record_id,
        name=record_id,
        description=f"{record_id} preparation cache fixture",
    )
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(
            FeatureLocation(30, 180, strand=1),
            type="CDS",
            qualifiers={
                "locus_tag": [f"{record_id}_0001"],
                "product": [product],
                "translation": ["M" * 49],
            },
        ),
        SeqFeature(
            FeatureLocation(240, 330, strand=-1),
            type="rRNA",
            qualifiers={"locus_tag": [f"{record_id}_r01"]},
        ),
    ]
    return record


def _write_genbank(path: Path, records: tuple[SeqRecord, ...]) -> None:
    SeqIO.write(records, path, "genbank")


def _document(
    source: Path,
    *,
    title: str | None = None,
    selected_features: tuple[str, ...] | None = None,
    region: str | None = None,
    reverse: bool = False,
    annotations: AnnotationOptions | None = None,
    feature_visibility_table: pd.DataFrame | None = None,
) -> dict[str, object]:
    request = CircularDiagramRequest(
        records=(
            RecordInput(
                source=GenBankInputSource(source),
                region=parse_region_spec(region) if region else None,
                presentation=RecordPresentation(reverse_complement=reverse),
                record_key="record-1",
            ),
        ),
        options=CircularDiagramOptions(
            plot_title=title,
            selected_features_set=selected_features,
            output=CircularOutputOptions(plot_title_position="top"),
            annotations=annotations,
            feature_visibility_table=feature_visibility_table,
        ),
        output=RenderOutputRequest(
            output_prefix="prepared-cache",
            formats=("svg",),
        ),
    )
    return build_session_document(request).to_dict()


def _render(
    tmp_path: Path,
    cache: PreparedBiologicalInputCache | None,
    document: dict[str, object],
    *,
    request_id: int,
    token_generation: int = 1,
    token_numbers: dict[str, int] | None = None,
    held_revision: str | None = None,
    held_segments: frozenset[str] = frozenset(),
    retention: dict[str, object] | None = None,
) -> tuple[dict[str, object], dict[str, object]]:
    resources = document["resources"]
    assert isinstance(resources, dict)
    workspace = tmp_path / f"gbdraw-web-render-{request_id}"
    resource_directory = workspace / "resources"
    resource_directory.mkdir(parents=True)
    (workspace / ".gbdraw-worker-render-workspace").touch()
    cache_directory = tmp_path / "gbdraw-web-render-resource-cache"
    cache_directory.mkdir(exist_ok=True)
    resource_paths: dict[str, Path] = {}
    identities: dict[str, dict[str, object]] = {}
    for index, (resource_id, entry) in enumerate(resources.items(), start=1):
        assert isinstance(entry, dict)
        content = base64.b64decode(str(entry["data"]), validate=True)
        token_number = (token_numbers or {}).get(
            resource_id,
            token_generation * 100 + index,
        )
        # A distinct stand-in for the sha256 hex token the Worker sends.
        cache_token = f"{token_number:064x}"
        cache_path = cache_directory / f"{cache_token}.bin"
        if cache_path.exists():
            assert cache_path.read_bytes() == content
        else:
            cache_path.write_bytes(content)
        request_path = resource_directory / f"{index:04d}.bin"
        request_path.symlink_to(cache_path)
        resource_paths[resource_id] = request_path
        identities[resource_id] = {
            "cacheToken": cache_token,
            "size": len(content),
        }
    diagnostics: dict[str, object] = {"timingsMs": {}, "metrics": {}}
    result = render_staged_canonical_web_request(
        document["renderRequest"],
        resource_paths=resource_paths,
        workspace=workspace,
        _diagnostics=diagnostics,
        _prepared_input_cache=cache,
        _resource_identities=identities if cache is not None else None,
        _held_feature_catalog_revision=held_revision,
        _held_segments=held_segments,
        _prepared_input_retention=retention,
    )
    return dict(result), diagnostics


def _metrics(diagnostics: dict[str, object]) -> dict[str, int]:
    metrics = diagnostics["metrics"]
    assert isinstance(metrics, dict)
    return {str(key): int(value) for key, value in metrics.items()}


def test_decoded_resource_cache_reuses_only_the_same_worker_identity(
    tmp_path: Path,
) -> None:
    path = tmp_path / "typed.json"
    path.write_text('{"schema":1}', encoding="utf-8")
    cache = PreparedBiologicalInputCache()
    identity = PreparedResourceIdentity("typed", "render-resource-1", path.stat().st_size)
    value = {"decoded": True}
    with cache.transaction(resource_paths={path: identity}, diagnostics=None):
        assert get_or_build_decoded_resource(
            path,
            ("canonical-typed-json", "result"),
            lambda: value,
        ) is value

    warm_diagnostics: dict[str, object] = {"metrics": {}}
    with cache.transaction(resource_paths={path: identity}, diagnostics=warm_diagnostics):
        assert get_or_build_decoded_resource(
            path,
            ("canonical-typed-json", "result"),
            lambda: pytest.fail("the unchanged typed resource was decoded twice"),
        ) is value
    assert _metrics(warm_diagnostics)["decodedResourceCacheHitCount"] == 1

    changed = PreparedResourceIdentity("typed", "render-resource-2", path.stat().st_size)
    replacement = {"decoded": "replacement"}
    changed_diagnostics: dict[str, object] = {"metrics": {}}
    with cache.transaction(resource_paths={path: changed}, diagnostics=changed_diagnostics):
        assert get_or_build_decoded_resource(
            path,
            ("canonical-typed-json", "result"),
            lambda: replacement,
        ) is replacement
    changed_metrics = _metrics(changed_diagnostics)
    assert changed_metrics["decodedResourceCacheMissCount"] == 1
    assert changed_metrics["decodedResourceBuildCount"] == 1
    assert changed_metrics["preparedInputCacheEvictionCount"] == 1


def test_warm_and_render_only_generates_reuse_biological_preparation(
    tmp_path: Path,
) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("cache-record"),))
    cache = PreparedBiologicalInputCache()
    cold_document = _document(source, title="Cold title")

    cold, cold_diagnostics = _render(
        tmp_path,
        cache,
        cold_document,
        request_id=1,
    )
    warm, warm_diagnostics = _render(
        tmp_path,
        cache,
        cold_document,
        request_id=2,
    )
    render_only, render_only_diagnostics = _render(
        tmp_path,
        cache,
        _document(source, title="Render-only title change"),
        request_id=3,
    )

    cold_metrics = _metrics(cold_diagnostics)
    assert cold_metrics["parsedSourceCacheMissCount"] == 1
    assert cold_metrics["parsedSourceParseCount"] == 1
    assert cold_metrics["resolvedRecordCacheMissCount"] == 1
    assert cold_metrics["resolvedRecordBuildCount"] == 1
    assert cold_metrics["renderStageBuild:interactive_context"] == 1
    assert cold_metrics["interactiveFeatureTraversalCount"] == 1
    assert cold_metrics["preparedInputCacheRetainedBytes"] > 0
    assert cold_metrics["preparedInputCacheMutationViolationCount"] == 0
    assert cold_metrics["preparedInputCacheFullCheckCount"] > 0

    for diagnostics in (warm_diagnostics, render_only_diagnostics):
        metrics = _metrics(diagnostics)
        assert metrics["parsedSourceCacheHitCount"] == 1
        assert metrics["parsedSourceParseCount"] == 0
        assert metrics["resolvedRecordCacheHitCount"] == 1
        assert metrics["resolvedRecordBuildCount"] == 0
        # Every counter of these stages, so a build would show as a second key.
        # The catalog item is a hit, so the normalization it reads is not asked for.
        for stage in ("interactive_context", "catalog_item"):
            assert {
                name: count for name, count in metrics.items() if name.endswith(f":{stage}")
            } == {f"renderStageHit:{stage}": 1}
        assert metrics["interactiveFeatureTraversalCount"] == 0
        assert metrics["preparedInputCacheRetainedBytes"] > 0
        assert metrics["preparedInputCacheMutationViolationCount"] == 0
        # O-2: neither request builds a stage that takes a record, so the
        # violation count above comes from the structure check only.
        assert metrics["preparedInputCacheFullCheckCount"] == 0

    assert warm["results"] == cold["results"]
    assert warm["metadata"] == cold["metadata"]
    assert warm["featureCatalog"] == cold["featureCatalog"]
    assert render_only["results"] != cold["results"]
    # Every top-level Circular group is a content (the plot title too), so the
    # catalog item is reused and no Generate after the first reads the SVG.
    assert _metrics(cold_diagnostics)["featureCatalogSvgParseCount"] == 1
    for diagnostics in (warm_diagnostics, render_only_diagnostics):
        metrics = _metrics(diagnostics)
        assert metrics.get("featureCatalogSvgParseCount", 0) == 0
        assert metrics.get("featureCatalogFullDomTraversalCount", 0) == 0
    assert not [name for name in _metrics(warm_diagnostics) if name.startswith("renderStageBuild:")]
    # A title change builds its own content and the composition plan, nothing else.
    assert {
        name for name in _metrics(render_only_diagnostics) if name.startswith("renderStageBuild:")
    } == {"renderStageBuild:circular_plot_title_content", "renderStageBuild:composition_plan"}
    assert render_only["featureCatalog"] == cold["featureCatalog"]


def _linear_document(source: Path, *, title: str) -> dict[str, object]:
    request = LinearDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(source), record_key="record-1"),),
        options=LinearDiagramOptions(plot_title=title),
        output=RenderOutputRequest(output_prefix="prepared-linear", formats=("svg",)),
    )
    return build_session_document(request).to_dict()


def test_a_linear_title_change_reuses_the_catalog_item(tmp_path: Path) -> None:
    """The catalog item reads the feature, comparison and annotation contents only."""

    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("linear-record"),))
    cache = PreparedBiologicalInputCache()
    document = _linear_document(source, title="Cold title")
    cold, cold_diagnostics = _render(tmp_path, cache, document, request_id=1)
    warm, warm_diagnostics = _render(tmp_path, cache, document, request_id=2)
    retitled, retitled_diagnostics = _render(
        tmp_path, cache, _linear_document(source, title="Render-only title change"), request_id=3
    )

    cold_metrics = _metrics(cold_diagnostics)
    assert cold_metrics["renderStageBuild:catalog_item"] == 1
    assert cold_metrics["featureCatalogSvgParseCount"] == 1
    assert cold_metrics["featureCatalogFeatureCandidateCount"] > 0
    for diagnostics in (warm_diagnostics, retitled_diagnostics):
        metrics = _metrics(diagnostics)
        assert metrics["renderStageHit:catalog_item"] == 1
        assert metrics.get("featureCatalogSvgParseCount", 0) == 0
        assert metrics.get("featureCatalogFullDomTraversalCount", 0) == 0
    assert not [name for name in _metrics(warm_diagnostics) if name.startswith("renderStageBuild:")]
    # A title change builds its own content and the composition plan, nothing else.
    assert {
        name for name in _metrics(retitled_diagnostics) if name.startswith("renderStageBuild:")
    } == {"renderStageBuild:plot_title_content", "renderStageBuild:composition_plan"}

    assert warm["results"] == cold["results"]
    assert warm["metadata"] == cold["metadata"]
    assert warm["featureCatalog"] == cold["featureCatalog"]
    assert retitled["results"] != cold["results"]
    assert retitled["featureCatalog"] == cold["featureCatalog"]


def _multi_record_document(source: Path, *, title: str) -> dict[str, object]:
    request = CircularDiagramRequest(
        records=(
            RecordInput(
                source=GenBankInputSource(source), record_key="records-1", cardinality=RecordCardinality.ALL
            ),
        ),
        options=CircularDiagramOptions(
            plot_title=title,
            output=CircularOutputOptions(plot_title_position="top"),
        ),
        output=RenderOutputRequest(output_prefix="prepared-multi-record", formats=("svg",)),
    )
    return build_session_document(request).to_dict()


def test_a_multi_record_circular_title_change_reuses_the_merged_records(tmp_path: Path) -> None:
    """The merged records and the catalog item read no title: a title change builds neither."""

    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("multi-a"), _record("multi-b")))
    cache = PreparedBiologicalInputCache()
    cold, cold_diagnostics = _render(tmp_path, cache, _multi_record_document(source, title="Cold title"), request_id=1)
    retitled, retitled_diagnostics = _render(
        tmp_path, cache, _multi_record_document(source, title="Render-only title change"), request_id=2
    )

    cold_metrics = _metrics(cold_diagnostics)
    assert cold_metrics["renderStageBuild:record_merge"] == 2
    assert cold_metrics["renderStageBuild:catalog_item"] == 1
    retitled_metrics = _metrics(retitled_diagnostics)
    assert retitled_metrics["renderStageHit:record_merge"] == 2
    assert retitled_metrics["renderStageHit:catalog_item"] == 1
    assert retitled_metrics.get("featureCatalogSvgParseCount", 0) == 0
    assert {name for name in retitled_metrics if name.startswith("renderStageBuild:")} == {
        "renderStageBuild:circular_plot_title_content",
        "renderStageBuild:composition_plan",
    }
    assert retitled["results"] != cold["results"]
    assert retitled["featureCatalog"] == cold["featureCatalog"]


def test_a_render_without_a_cache_always_returns_the_catalog_without_a_revision(tmp_path: Path) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("uncached-catalog"),))
    document = _document(source, title="Uncached")
    uncached, _ = _render(tmp_path, None, document, request_id=1, held_revision="e:1")
    cached, _ = _render(tmp_path, PreparedBiologicalInputCache(), document, request_id=2)

    assert uncached["featureCatalogRevision"] is None
    assert uncached["featureCatalog"] == cached["featureCatalog"]
    assert isinstance(cached["featureCatalogRevision"], str)
    assert "featureCatalog" not in uncached["metadata"]


def test_public_render_request_remains_uncached_by_default(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("public-record"),))
    request = CircularDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(source)),),
        output=RenderOutputRequest(
            output_prefix="public-uncached",
            output_directory=tmp_path,
            formats=("svg",),
            overwrite=True,
        ),
    )
    original_loader = request_render_module.load_gbks
    load_count = 0

    def counting_loader(paths: list[str]) -> list[SeqRecord]:
        nonlocal load_count
        load_count += 1
        return original_loader(paths)

    monkeypatch.setattr(request_render_module, "load_gbks", counting_loader)

    first = request_render_module.render_request(request)
    second = request_render_module.render_request(request)

    assert load_count == 2
    assert second is not first


def test_gallery_genbank_preparation_owner_remains_read_only(tmp_path: Path) -> None:
    session_path = (
        Path(__file__).parents[1]
        / "gbdraw"
        / "web"
        / "gallery"
        / "sessions"
        / "HmmtDNA_basic_circular.gbdraw-session.json"
    )
    document = json.loads(session_path.read_text(encoding="utf-8"))
    cache = PreparedBiologicalInputCache()

    first, first_diagnostics = _render(
        tmp_path,
        cache,
        document,
        request_id=1,
    )
    second, second_diagnostics = _render(
        tmp_path,
        cache,
        document,
        request_id=2,
    )

    assert second["results"] == first["results"]
    assert _metrics(first_diagnostics)["preparedInputCacheMutationViolationCount"] == 0
    assert _metrics(second_diagnostics)["parsedSourceCacheHitCount"] == 1


def test_annotation_change_invalidates_interactive_context_only(
    tmp_path: Path,
) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("cache-record"),))
    first_annotations = AnnotationOptions(
        sets=(
            AnnotationSet(
                "regions",
                (
                    RegionAnnotation(
                        "focus",
                        CoordinateSpan(None, 20, 80),
                        label="First",
                    ),
                ),
            ),
        )
    )
    changed_annotations = AnnotationOptions(
        sets=(
            AnnotationSet(
                "regions",
                (
                    RegionAnnotation(
                        "focus",
                        CoordinateSpan(None, 40, 100),
                        label="Changed",
                    ),
                ),
            ),
        )
    )
    cache = PreparedBiologicalInputCache()
    document = _document(source, annotations=first_annotations)
    first, _first_diagnostics = _render(tmp_path, cache, document, request_id=1)
    _warm, warm_diagnostics = _render(
        tmp_path,
        cache,
        document,
        request_id=2,
    )
    changed, changed_diagnostics = _render(
        tmp_path,
        cache,
        _document(source, annotations=changed_annotations),
        request_id=3,
    )

    assert _metrics(warm_diagnostics)["renderStageHit:interactive_context"] == 1
    changed_metrics = _metrics(changed_diagnostics)
    assert changed_metrics["parsedSourceCacheHitCount"] == 1
    assert changed_metrics["resolvedRecordCacheHitCount"] == 1
    assert changed_metrics["renderStageBuild:interactive_context"] == 1
    assert changed["featureCatalog"] != first["featureCatalog"]


def test_visibility_resource_invalidates_context_but_not_genbank_records(
    tmp_path: Path,
) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("cache-record"),))
    cache = PreparedBiologicalInputCache()
    first, _first_diagnostics = _render(
        tmp_path,
        cache,
        _document(source),
        request_id=1,
    )
    visibility = pd.DataFrame(
        [["*", "CDS", "product", ".*", "off"]],
        columns=["record_id", "feature_type", "qualifier", "value", "action"],
    )

    result, diagnostics = _render(
        tmp_path,
        cache,
        _document(source, feature_visibility_table=visibility),
        request_id=2,
    )
    metrics = _metrics(diagnostics)
    assert metrics["parsedSourceCacheHitCount"] == 1
    assert metrics["resolvedRecordCacheHitCount"] == 1
    assert metrics["renderStageBuild:interactive_context"] == 1
    assert result["results"] != first["results"]


def _orthogroup_document(source: Path, group_id: str) -> dict[str, object]:
    protein_id = "cache-record_0001"
    member = OrthogroupMember(
        orthogroup_id=group_id,
        protein_id=protein_id,
        record_index=0,
        feature_index=0,
        record_id="cache-record",
        label=group_id,
        start=30,
        end=180,
        strand=1,
        feature_svg_id=None,
        source_protein_id=protein_id,
    )
    orthogroups = OrthogroupResult(
        orthogroups={group_id: [member]},
        member_by_protein_id={protein_id: member},
    )
    request = LinearDiagramRequest(
        records=(
            RecordInput(
                source=GenBankInputSource(source),
                record_key="record-1",
            ),
        ),
        options=LinearDiagramOptions(orthogroups=orthogroups),
        output=RenderOutputRequest(
            output_prefix="prepared-orthogroups",
            formats=("svg",),
        ),
    )
    return build_session_document(request).to_dict()


def test_orthogroup_resource_invalidates_interactive_context_only(
    tmp_path: Path,
) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("cache-record"),))
    first_document = _orthogroup_document(source, "OG-A")
    changed_document = _orthogroup_document(source, "OG-B")
    changed_resources = changed_document["resources"]
    assert isinstance(changed_resources, dict)
    orthogroup_resource_id = next(
        resource_id
        for resource_id in changed_resources
        if "orthogroup" in resource_id
    )
    cache = PreparedBiologicalInputCache()

    first, _first_diagnostics = _render(
        tmp_path,
        cache,
        first_document,
        request_id=1,
    )
    changed, changed_diagnostics = _render(
        tmp_path,
        cache,
        changed_document,
        request_id=2,
        token_numbers={orthogroup_resource_id: 999},
    )

    metrics = _metrics(changed_diagnostics)
    assert metrics["parsedSourceCacheHitCount"] == 1
    assert metrics["parsedSourceParseCount"] == 0
    assert metrics["resolvedRecordCacheHitCount"] == 1
    assert metrics["resolvedRecordBuildCount"] == 0
    assert metrics["renderStageBuild:interactive_context"] == 1
    assert changed["featureCatalog"] != first["featureCatalog"]
    changed_payload = json.dumps(changed, sort_keys=True)
    assert "OG-B" in changed_payload
    assert "OG-A" not in changed_payload


def test_gff_fasta_type_filter_is_part_of_the_resolved_key_not_the_parsed_key(tmp_path: Path) -> None:
    fixture_root = Path(__file__).parent / "test_inputs"
    gff_path = fixture_root / "NC_013668.gff3"
    fasta_path = fixture_root / "NC_013668.fasta"
    identities = {
        gff_path: PreparedResourceIdentity(
            "record-1-gff3",
            "render-resource-1",
            gff_path.stat().st_size,
        ),
        fasta_path: PreparedResourceIdentity(
            "record-1-fasta",
            "render-resource-2",
            fasta_path.stat().st_size,
        ),
    }
    cache = PreparedBiologicalInputCache()

    def request(selected_features: tuple[str, ...], **options) -> CircularDiagramRequest:
        return CircularDiagramRequest(
            records=(
                RecordInput(
                    source=GffFastaInputSource(gff_path, fasta_path),
                    record_key="record-1",
                ),
            ),
            options=CircularDiagramOptions(
                selected_features_set=selected_features,
                **options,
            ),
        )

    cold_diagnostics: dict[str, object] = {"metrics": {}}
    with cache.transaction(resource_paths=identities, diagnostics=cold_diagnostics):
        assert len(normalize_request_records(request(("CDS",)))) == 1
    warm_diagnostics: dict[str, object] = {"metrics": {}}
    with cache.transaction(resource_paths=identities, diagnostics=warm_diagnostics):
        assert len(normalize_request_records(request(("CDS",)))) == 1
    changed_diagnostics: dict[str, object] = {"metrics": {}}
    with cache.transaction(resource_paths=identities, diagnostics=changed_diagnostics):
        assert len(normalize_request_records(request(("rRNA",)))) == 1

    assert _metrics(cold_diagnostics)["parsedSourceParseCount"] == 1
    assert _metrics(warm_diagnostics)["parsedSourceCacheHitCount"] == 1
    assert _metrics(warm_diagnostics)["resolvedRecordCacheHitCount"] == 1
    assert _metrics(changed_diagnostics)["parsedSourceCacheHitCount"] == 1
    assert _metrics(changed_diagnostics)["parsedSourceParseCount"] == 0
    assert _metrics(changed_diagnostics)["resolvedRecordBuildCount"] == 1

    # An override row that shows a type the filter dropped reloads from the
    # same parse: the type filter is applied in memory.
    gene = next(
        entry.biological_feature_id
        for entry in plan_request(request(("CDS",))).provenance[0].source_feature_catalog
        if entry.feature_type == "gene"
    )
    shown_diagnostics: dict[str, object] = {"metrics": {}}
    with cache.transaction(resource_paths=identities, diagnostics=shown_diagnostics):
        [record] = normalize_request_records(request(
            ("CDS",), feature_overrides=(FeatureOverride("record-1", gene, feature_visibility="on"),),
        ))
    assert "gene" in {feature.type for feature in record.features}
    assert _metrics(shown_diagnostics)["parsedSourceCacheHitCount"] == 1
    assert _metrics(shown_diagnostics)["parsedSourceParseCount"] == 0


def test_comparison_fasta_cache_reuses_only_nonempty_sources(tmp_path: Path) -> None:
    fasta_path = tmp_path / "comparison.fasta"
    fasta_path.write_text(">comparison\nATGATGATG\n", encoding="utf-8")
    identity = PreparedResourceIdentity(
        "comparison-1-fasta",
        "render-resource-1",
        fasta_path.stat().st_size,
    )
    cache = PreparedBiologicalInputCache()
    diagnostics: dict[str, object] = {"metrics": {}}
    with cache.transaction(
        resource_paths={fasta_path: identity},
        diagnostics=diagnostics,
    ):
        assert request_render_module._ComparisonSequenceSources(
            (str(fasta_path),)
        ).load()[0][0].id == "comparison"
    warm_diagnostics: dict[str, object] = {"metrics": {}}
    with cache.transaction(
        resource_paths={fasta_path: identity},
        diagnostics=warm_diagnostics,
    ):
        assert request_render_module._ComparisonSequenceSources(
            (str(fasta_path),)
        ).load()[0][0].id == "comparison"
    assert _metrics(warm_diagnostics)["parsedSourceCacheHitCount"] == 1
    assert _metrics(warm_diagnostics)["parsedSourceParseCount"] == 0

    empty_path = tmp_path / "empty.fasta"
    empty_path.touch()
    empty_identity = PreparedResourceIdentity(
        "comparison-1-fasta",
        "render-resource-2",
        0,
    )
    for _request_index in range(2):
        empty_diagnostics: dict[str, object] = {"metrics": {}}
        with cache.transaction(
            resource_paths={empty_path: empty_identity},
            diagnostics=empty_diagnostics,
        ):
            assert request_render_module._ComparisonSequenceSources(
                (str(empty_path),)
            ).load() == ((),)
        assert _metrics(empty_diagnostics)["parsedSourceCacheMissCount"] == 1
        assert _metrics(empty_diagnostics)["parsedSourceParseCount"] == 1


def test_selector_and_cardinality_changes_rebuild_resolved_records_only(
    tmp_path: Path,
) -> None:
    source = tmp_path / "multi-record.gbk"
    _write_genbank(source, (_record("first"), _record("second")))
    identity = PreparedResourceIdentity(
        "record-1-genbank",
        "render-resource-1",
        source.stat().st_size,
    )
    cache = PreparedBiologicalInputCache()

    def request(
        cardinality: RecordCardinality,
        selector: str | None = None,
    ) -> CircularDiagramRequest:
        return CircularDiagramRequest(
            records=(
                RecordInput(
                    source=GenBankInputSource(source),
                    cardinality=cardinality,
                    selector=parse_record_selector(selector),
                    record_key="record-1",
                ),
            ),
        )

    with cache.transaction(resource_paths={source: identity}, diagnostics=None):
        assert len(normalize_request_records(request(RecordCardinality.FIRST))) == 1
    all_diagnostics: dict[str, object] = {"metrics": {}}
    with cache.transaction(
        resource_paths={source: identity},
        diagnostics=all_diagnostics,
    ):
        assert len(normalize_request_records(request(RecordCardinality.ALL))) == 2
    selected_diagnostics: dict[str, object] = {"metrics": {}}
    with cache.transaction(
        resource_paths={source: identity},
        diagnostics=selected_diagnostics,
    ):
        records = normalize_request_records(
            request(RecordCardinality.EXACTLY_ONE, "#2")
        )
        assert [record.id for record in records] == ["second"]

    for diagnostics in (all_diagnostics, selected_diagnostics):
        metrics = _metrics(diagnostics)
        assert metrics["parsedSourceCacheHitCount"] == 1
        assert metrics["parsedSourceParseCount"] == 0
        assert metrics["resolvedRecordCacheMissCount"] == 1
        assert metrics["resolvedRecordBuildCount"] == 1


@pytest.mark.parametrize(
    ("kwargs", "expected_parse_hits"),
    [
        ({"region": "20-400"}, 1),
        ({"reverse": True}, 1),
        ({"selected_features": ("rRNA",)}, 1),
    ],
)
def test_semantic_changes_invalidate_only_the_dependent_layers(
    tmp_path: Path,
    kwargs: dict[str, object],
    expected_parse_hits: int,
) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("cache-record"),))
    cache = PreparedBiologicalInputCache()
    document = _document(source)
    _render(tmp_path, cache, document, request_id=1)
    changed_document = copy.deepcopy(document)
    render_request = changed_document["renderRequest"]
    assert isinstance(render_request, dict)
    records = render_request["records"]
    assert isinstance(records, list)
    record = records[0]
    assert isinstance(record, dict)
    if "region" in kwargs:
        record["region"] = {
            "selector": None,
            "start": 20,
            "end": 400,
            "reverseComplement": False,
        }
    if "reverse" in kwargs:
        presentation = record["presentation"]
        assert isinstance(presentation, dict)
        presentation["reverseComplement"] = True
    if "selected_features" in kwargs:
        options = render_request["diagramOptions"]
        assert isinstance(options, dict)
        options["selectedFeaturesSet"] = ["rRNA"]

    result, diagnostics = _render(
        tmp_path,
        cache,
        changed_document,
        request_id=2,
    )
    metrics = _metrics(diagnostics)
    assert metrics["parsedSourceCacheHitCount"] == expected_parse_hits
    assert metrics["parsedSourceParseCount"] == 0
    if "selected_features" in kwargs:
        assert metrics["resolvedRecordCacheHitCount"] == 1
        assert metrics["resolvedRecordBuildCount"] == 0
    else:
        assert metrics["resolvedRecordCacheMissCount"] == 1
        assert metrics["resolvedRecordBuildCount"] == 1
    assert metrics["renderStageBuild:interactive_context"] == 1
    assert metrics["interactiveFeatureTraversalCount"] == 1


def test_source_token_replacement_invalidates_every_layer(tmp_path: Path) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("cache-record"),))
    cache = PreparedBiologicalInputCache()
    original, _original_diagnostics = _render(
        tmp_path,
        cache,
        _document(source),
        request_id=1,
    )
    _write_genbank(source, (_record("cache-record", product="replacement"),))

    result, diagnostics = _render(
        tmp_path,
        cache,
        _document(source),
        request_id=2,
        token_generation=2,
    )
    metrics = _metrics(diagnostics)
    assert metrics["parsedSourceCacheMissCount"] == 1
    assert metrics["parsedSourceParseCount"] == 1
    assert metrics["resolvedRecordCacheMissCount"] == 1
    assert metrics["resolvedRecordBuildCount"] == 1
    assert metrics["renderStageBuild:interactive_context"] == 1
    assert metrics["preparedInputCacheEvictionCount"] >= 3
    assert result != original
    replacement_payload = json.dumps(result, sort_keys=True)
    assert "replacement" in replacement_payload
    assert "cache target" not in replacement_payload


def _token(number: int) -> str:
    return f"{number:064x}"


def test_each_drawing_keeps_the_prepared_inputs_of_its_last_request(tmp_path: Path) -> None:
    # A4b: two drawings alternate in one Worker. Each keeps the parsed and
    # resolved inputs of its own last request while the page retains its
    # resources; the render stages stay the last successful render's.
    first = tmp_path / "first.gbk"
    second = tmp_path / "second.gbk"
    _write_genbank(first, (_record("first-record"),))
    _write_genbank(second, (_record("second-record"),))
    cache = PreparedBiologicalInputCache()
    first_token, second_token = _token(101), _token(201)

    def render(
        document: dict[str, object],
        request_id: int,
        drawing: str,
        retained: tuple[str, ...],
    ) -> dict[str, int]:
        _result, diagnostics = _render(
            tmp_path,
            cache,
            document,
            request_id=request_id,
            # The first drawing's file has token 101, the second's 201.
            token_generation=1 if drawing == "circular" else 2,
            retention={"drawing": drawing, "resourceTokens": list(retained)},
        )
        return _metrics(diagnostics)

    first_document = _document(first)
    render(_document(first, region="20-400"), 1, "circular", ())
    render(first_document, 2, "circular", ())
    second_metrics = render(_document(second), 3, "linear", (first_token,))
    assert second_metrics["parsedSourceParseCount"] == 1
    returned = render(first_document, 4, "circular", (second_token,))
    assert returned["parsedSourceCacheHitCount"] == 1
    assert returned["parsedSourceParseCount"] == 0
    assert returned["resolvedRecordCacheHitCount"] == 1
    assert returned["resolvedRecordBuildCount"] == 0
    # The render stages of the other drawing's request were not kept.
    assert returned["renderStageBuild:interactive_context"] == 1
    # A drawing keeps only its last request's inputs: its earlier region
    # variant was dropped when its next request committed.
    variant = render(_document(first, region="20-400"), 5, "circular", (second_token,))
    assert variant["parsedSourceCacheHitCount"] == 1
    assert variant["resolvedRecordBuildCount"] == 1
    # A token no request retains any more releases its drawing's inputs.
    render(_document(second), 6, "linear", ())
    released = render(first_document, 7, "circular", ())
    assert released["parsedSourceParseCount"] == 1


@pytest.mark.parametrize(
    "retention",
    [
        {"drawing": "circular"},
        {"drawing": "Circular!", "resourceTokens": []},
        {"drawing": "circular", "resourceTokens": ["render-resource-1"]},
        {"drawing": "circular", "resourceTokens": "f" * 64},
        {"drawing": "circular", "resourceTokens": [], "extra": True},
    ],
)
def test_a_malformed_retention_is_rejected(
    tmp_path: Path,
    retention: dict[str, object],
) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("cache-record"),))
    with pytest.raises(ValidationError, match="retention"):
        _render(
            tmp_path,
            PreparedBiologicalInputCache(),
            _document(source),
            request_id=1,
            retention=retention,
        )


def test_identity_validation_rejects_paths_and_size_mismatches(tmp_path: Path) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("cache-record"),))
    document = _document(source)
    resources = document["resources"]
    assert isinstance(resources, dict)
    workspace = tmp_path / "gbdraw-web-render-1"
    resource_directory = workspace / "resources"
    resource_directory.mkdir(parents=True)
    (workspace / ".gbdraw-worker-render-workspace").touch()
    resource_paths: dict[str, Path] = {}
    identities: dict[str, dict[str, object]] = {}
    for index, (resource_id, entry) in enumerate(resources.items(), start=1):
        assert isinstance(entry, dict)
        content = base64.b64decode(str(entry["data"]), validate=True)
        path = resource_directory / f"{index:04d}.bin"
        path.write_bytes(content)
        resource_paths[resource_id] = path
        identities[resource_id] = {
            "cacheToken": f"{index:064x}",
            "size": len(content) + 1,
            "path": str(path),
        }

    with pytest.raises(ValidationError, match="Invalid prepared resource identity"):
        render_staged_canonical_web_request(
            document["renderRequest"],
            resource_paths=resource_paths,
            workspace=workspace,
            _prepared_input_cache=PreparedBiologicalInputCache(),
            _resource_identities=identities,
        )
    assert not workspace.exists()


def test_failed_fill_is_not_published_and_prior_entry_survives(
    tmp_path: Path,
) -> None:
    path = tmp_path / "resource.bin"
    path.write_bytes(b"resource")
    identity = PreparedResourceIdentity(
        resource_id="record-1",
        cache_token="render-resource-1",
        size=8,
    )
    identities = frozenset({identity})
    cache = PreparedBiologicalInputCache()
    original = (_record("original"),)
    with cache.transaction(resource_paths={path: identity}, diagnostics=None):
        assert get_or_build_parsed_source(
            ("source", identity),
            identities,
            lambda: original,
        ) is original

    with pytest.raises(RuntimeError, match="fill failed"):
        with cache.transaction(resource_paths={path: identity}, diagnostics=None):
            assert get_or_build_parsed_source(
                ("source", identity),
                identities,
                lambda: pytest.fail("the committed source should be reused"),
            ) is original
            get_or_build_parsed_source(
                ("replacement", identity),
                identities,
                lambda: (_ for _ in ()).throw(RuntimeError("fill failed")),
            )

    diagnostics: dict[str, object] = {"metrics": {}}
    with cache.transaction(resource_paths={path: identity}, diagnostics=diagnostics):
        assert get_or_build_parsed_source(
            ("source", identity),
            identities,
            lambda: pytest.fail("the prior committed source was poisoned"),
        ) is original
    assert _metrics(diagnostics)["parsedSourceCacheHitCount"] == 1


def test_mutation_guard_rejects_changed_cached_owners(tmp_path: Path) -> None:
    path = tmp_path / "resource.bin"
    path.write_bytes(b"resource")
    identity = PreparedResourceIdentity(
        resource_id="record-1",
        cache_token="render-resource-1",
        size=8,
    )
    identities = frozenset({identity})
    cache = PreparedBiologicalInputCache()
    records = (_record("mutable"),)
    with cache.transaction(resource_paths={path: identity}, diagnostics=None):
        get_or_build_parsed_source(("source", identity), identities, lambda: records)
    records[0].features.append(
        SeqFeature(FeatureLocation(400, 450), type="misc_feature")
    )

    diagnostics: dict[str, object] = {"metrics": {}}
    with pytest.raises(ValidationError, match="ownership was mutated"):
        with cache.transaction(
            resource_paths={path: identity},
            diagnostics=diagnostics,
        ):
            get_or_build_parsed_source(
                ("source", identity),
                identities,
                lambda: pytest.fail("a mutated owner must not be reused"),
            )
    assert _metrics(diagnostics)["preparedInputCacheMutationViolationCount"] == 1


def _parsed_source_cache(
    tmp_path: Path,
    sources: dict[str, tuple[SeqRecord, ...]],
) -> tuple[PreparedBiologicalInputCache, Path, PreparedResourceIdentity]:
    path = tmp_path / "resource.bin"
    path.write_bytes(b"resource")
    identity = PreparedResourceIdentity(
        resource_id="record-1",
        cache_token="render-resource-1",
        size=8,
    )
    cache = PreparedBiologicalInputCache()
    with cache.transaction(resource_paths={path: identity}, diagnostics=None):
        for name, records in sources.items():
            get_or_build_parsed_source(
                (name, identity),
                frozenset({identity}),
                lambda records=records: records,
            )
    return cache, path, identity


def _feature_count(record: SeqRecord, _salt: str) -> int:
    return len(record.features)


def _read_in_a_record_stage(record: SeqRecord, salt: str) -> None:
    """Build a render stage that takes ``record`` (a new ``salt`` builds again): the commit checks owners in full."""

    register_prepared_value_identity(record, (prepared_module.RECORD_IDENTITY_TAG, "test", salt))
    render_stage("gc_skew_data", _feature_count, record, salt)


def _replace_first_feature(records: tuple[SeqRecord, ...]) -> None:
    """Change a cached record without changing its feature count (D-05)."""

    records[0].features[0] = SeqFeature(
        FeatureLocation(33, 180, strand=1),
        type="CDS",
        qualifiers={"locus_tag": ["replacement"]},
    )


def test_a_read_only_lookup_hits_and_never_publishes_prunes_or_reorders(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    # A Worker helper's lookup (OV-1550): the entries a Generate stored are
    # read under the same keys; a miss builds once for the call and is not
    # kept; an entry the lookup does not name, and its order, stay.
    first = (_record("first"),)
    second = (_record("second"),)
    cache, path, identity = _parsed_source_cache(tmp_path, {"first": first, "second": second})
    stage_key = ("render-stage", "test", identity)
    with cache.transaction(resource_paths={path: identity}, diagnostics=None):
        get_or_build_parsed_source(("first", identity), frozenset({identity}), lambda: first)
        get_or_build_parsed_source(("second", identity), frozenset({identity}), lambda: second)
        transaction = prepared_module.active_render_stage_transaction()
        assert transaction is not None
        prepared_module.store_render_stage(transaction, stage_key, object(), frozenset())
    before = {layer: list(entries.items()) for layer, entries in cache._layers.items()}

    other_path = tmp_path / "other.bin"
    other_path.write_bytes(b"other")
    other = PreparedResourceIdentity("record-2", "render-resource-2", 5)
    builds: list[str] = []
    with cache.transaction(
        resource_paths={other_path: other, path: identity}, diagnostics=None, read_only=True,
    ):
        # The first entry was stored before the second; a hit does not reorder.
        assert get_or_build_parsed_source(
            ("first", identity), frozenset({identity}), lambda: pytest.fail("a stored source was parsed again"),
        ) is first
        miss = get_or_build_parsed_source(
            ("new", other), frozenset({other}), lambda: builds.append("new") or (_record("new"),),
        )
        assert get_or_build_parsed_source(("new", other), frozenset({other}), lambda: pytest.fail("parsed twice")) is miss
        transaction = prepared_module.active_render_stage_transaction()
        assert transaction is not None
        # A stage result is not kept, so its key is not walked for identities.
        with monkeypatch.context() as patch:
            patch.setattr(prepared_module, "_key_resource_identities", lambda key: pytest.fail("key walked"))
            prepared_module.store_render_stage(transaction, ("render-stage", "helper", other), object(), frozenset())
    assert builds == ["new"]
    assert {layer: list(entries.items()) for layer, entries in cache._layers.items()} == before


def test_a_read_only_lookup_checks_the_owners_it_read(tmp_path: Path) -> None:
    records = (_record("guarded"),)
    cache, path, identity = _parsed_source_cache(tmp_path, {"source": records})
    with pytest.raises(ValidationError, match="ownership was mutated"):
        with cache.transaction(resource_paths={path: identity}, diagnostics=None, read_only=True):
            cached = get_or_build_parsed_source(
                ("source", identity), frozenset({identity}), lambda: pytest.fail("parsed again"),
            )
            # A record built in the lookup makes the end check walk every feature.
            get_or_build_parsed_source(("built", identity), frozenset({identity}), lambda: (_record("built"),))
            _replace_first_feature(cached)
    # The changed entry is dropped, and the next Generate checks in full.
    assert ("source", identity) not in cache._layers["parsed"]
    assert cache._previous_transaction_failed is True


def test_structural_stamp_rejects_a_resized_feature_list_at_entry(
    tmp_path: Path,
) -> None:
    records = (_record("resized"),)
    cache, path, identity = _parsed_source_cache(tmp_path, {"source": records})
    records[0].features.append(SeqFeature(FeatureLocation(400, 450), type="misc_feature"))

    diagnostics: dict[str, object] = {"metrics": {}}
    reached_render = False
    with pytest.raises(ValidationError, match="ownership was mutated .*feature count"):
        with cache.transaction(resource_paths={path: identity}, diagnostics=diagnostics):
            get_or_build_parsed_source(
                ("source", identity),
                frozenset({identity}),
                lambda: pytest.fail("a mutated owner must not be reused"),
            )
            reached_render = True
    assert reached_render is False
    assert _metrics(diagnostics)["preparedInputCacheMutationViolationCount"] == 1


def test_feature_replaced_at_the_same_length_fails_at_commit_of_the_same_request(
    tmp_path: Path,
) -> None:
    # D-05 (R2-8): the structural check at entry passes, the render runs, and
    # the full check at commit (a record stage built) fails the same request.
    records = (_record("replaced"),)
    cache, path, identity = _parsed_source_cache(tmp_path, {"source": records})
    _replace_first_feature(records)

    diagnostics: dict[str, object] = {"metrics": {}}
    reached_render = False
    with pytest.raises(ValidationError, match="ownership was mutated .*feature 0"):
        with cache.transaction(resource_paths={path: identity}, diagnostics=diagnostics):
            assert get_or_build_parsed_source(
                ("source", identity),
                frozenset({identity}),
                lambda: pytest.fail("the entry check is structural"),
            ) is records
            _read_in_a_record_stage(records[0], "same-request")
            reached_render = True
    assert reached_render is True
    assert _metrics(diagnostics)["preparedInputCacheMutationViolationCount"] == 1
    assert _metrics(diagnostics)["preparedInputCacheFullCheckCount"] == 1

    rebuilt = (_record("rebuilt"),)
    with cache.transaction(resource_paths={path: identity}, diagnostics=None):
        assert get_or_build_parsed_source(
            ("source", identity),
            frozenset({identity}),
            lambda: rebuilt,
        ) is rebuilt


def test_a_request_without_a_record_stage_skips_the_walk_and_the_next_record_stage_finds_the_write(
    tmp_path: Path,
) -> None:
    # O-2 (D-07): a request that builds no record entry and no stage that
    # takes a record checks structure only; the change is found at the commit
    # of the next request that runs a record stage.
    records = (_record("skipped"),)
    cache, path, identity = _parsed_source_cache(tmp_path, {"source": records})
    _replace_first_feature(records)

    skipped: dict[str, object] = {"metrics": {}}
    with cache.transaction(resource_paths={path: identity}, diagnostics=skipped):
        assert get_or_build_parsed_source(
            ("source", identity),
            frozenset({identity}),
            lambda: pytest.fail("the entry check is structural"),
        ) is records
    assert _metrics(skipped)["preparedInputCacheFullCheckCount"] == 0
    assert _metrics(skipped)["preparedInputCacheMutationViolationCount"] == 0

    diagnostics: dict[str, object] = {"metrics": {}}
    with pytest.raises(ValidationError, match="ownership was mutated .*feature 0"):
        with cache.transaction(resource_paths={path: identity}, diagnostics=diagnostics):
            get_or_build_parsed_source(
                ("source", identity),
                frozenset({identity}),
                lambda: pytest.fail("the entry check is structural"),
            )
            _read_in_a_record_stage(records[0], "next-request")
    assert _metrics(diagnostics)["preparedInputCacheFullCheckCount"] == 1
    assert _metrics(diagnostics)["preparedInputCacheMutationViolationCount"] == 1


def test_a_request_that_builds_a_record_entry_checks_it_in_full_without_a_record_stage(
    tmp_path: Path,
) -> None:
    """O-2 (D-07): a new entry that holds records is walked at its commit, with no stage run."""

    cache, path, identity = _parsed_source_cache(tmp_path, {})
    diagnostics: dict[str, object] = {"metrics": {}}
    with cache.transaction(resource_paths={path: identity}, diagnostics=diagnostics):
        get_or_build_parsed_source(
            ("source", identity),
            frozenset({identity}),
            lambda: (_record("built"),),
        )
    assert _metrics(diagnostics)["preparedInputCacheFullCheckCount"] == 1
    assert _metrics(diagnostics).get("renderStageBuildCount", 0) == 0


def test_an_uncacheable_stage_checks_owners_in_full(tmp_path: Path) -> None:
    """A stage without a key may take a record without an identity: the commit walks."""

    records = (_record("uncacheable"),)
    cache, path, identity = _parsed_source_cache(tmp_path, {"source": records})
    _replace_first_feature(records)

    diagnostics: dict[str, object] = {"metrics": {}}
    with pytest.raises(ValidationError, match="ownership was mutated .*feature 0"):
        with cache.transaction(resource_paths={path: identity}, diagnostics=diagnostics):
            get_or_build_parsed_source(
                ("source", identity),
                frozenset({identity}),
                lambda: pytest.fail("the entry check is structural"),
            )
            # No identity is registered for this record: the stage is uncacheable.
            render_stage("gc_skew_data", _feature_count, records[0], "no-identity")
    assert _metrics(diagnostics)["renderStageUncacheableCount"] == 1
    assert _metrics(diagnostics)["preparedInputCacheFullCheckCount"] == 1


def test_request_after_a_mutation_error_checks_every_owner_in_full_at_entry(
    tmp_path: Path,
) -> None:
    first = (_record("first"),)
    second = (_record("second"),)
    cache, path, identity = _parsed_source_cache(
        tmp_path,
        {"first": first, "second": second},
    )
    _replace_first_feature(first)
    _replace_first_feature(second)

    with pytest.raises(ValidationError, match="ownership was mutated"):
        with cache.transaction(resource_paths={path: identity}, diagnostics=None):
            for name in ("first", "second"):
                get_or_build_parsed_source(
                    (name, identity),
                    frozenset({identity}),
                    lambda: pytest.fail("the entry check is structural"),
                )
            _read_in_a_record_stage(first[0], "mutation-error")

    # The commit stopped at the first changed owner; the other is still cached
    # and is now found before the render.
    reached_render = False
    with pytest.raises(ValidationError, match="ownership was mutated .*feature 0"):
        with cache.transaction(resource_paths={path: identity}, diagnostics=None):
            for name in ("first", "second"):
                # The evicted owner is rebuilt; the other one fails at entry.
                get_or_build_parsed_source(
                    (name, identity),
                    frozenset({identity}),
                    lambda: (_record("rebuilt"),),
                )
            reached_render = True
    assert reached_render is False


def test_request_after_a_failed_render_checks_owners_in_full_at_entry(
    tmp_path: Path,
) -> None:
    # D-05 (R2-8): when the render fails for another reason first, that request
    # shows the render's error, and the next request finds the change at entry.
    records = (_record("failed-render"),)
    cache, path, identity = _parsed_source_cache(tmp_path, {"source": records})
    _replace_first_feature(records)

    failed_diagnostics: dict[str, object] = {"metrics": {}}
    with pytest.raises(RuntimeError, match="render failed"):
        with cache.transaction(
            resource_paths={path: identity},
            diagnostics=failed_diagnostics,
        ):
            get_or_build_parsed_source(
                ("source", identity),
                frozenset({identity}),
                lambda: pytest.fail("the entry check is structural"),
            )
            raise RuntimeError("render failed")
    assert _metrics(failed_diagnostics)["preparedInputCacheMutationViolationCount"] == 0

    diagnostics: dict[str, object] = {"metrics": {}}
    reached_render = False
    with pytest.raises(ValidationError, match="ownership was mutated .*feature 0"):
        with cache.transaction(resource_paths={path: identity}, diagnostics=diagnostics):
            get_or_build_parsed_source(
                ("source", identity),
                frozenset({identity}),
                lambda: pytest.fail("a mutated owner must not be reused"),
            )
            reached_render = True
    assert reached_render is False
    metrics = _metrics(diagnostics)
    assert metrics["preparedInputCacheMutationViolationCount"] == 1
    assert metrics["preparedInputCacheEvictionCount"] == 1

    # After a successful transaction the entry check is structural again.
    rebuilt = (_record("rebuilt"),)
    with cache.transaction(resource_paths={path: identity}, diagnostics=None):
        get_or_build_parsed_source(
            ("source", identity),
            frozenset({identity}),
            lambda: rebuilt,
        )
    _replace_first_feature(rebuilt)
    reached_render = False
    with pytest.raises(ValidationError, match="ownership was mutated .*feature 0"):
        with cache.transaction(resource_paths={path: identity}, diagnostics=None):
            get_or_build_parsed_source(
                ("source", identity),
                frozenset({identity}),
                lambda: pytest.fail("the entry check is structural"),
            )
            _read_in_a_record_stage(rebuilt[0], "after-success")
            reached_render = True
    assert reached_render is True


def test_a_request_walks_each_cached_feature_once_only_when_a_record_stage_builds(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("walk-record"),))
    cache = PreparedBiologicalInputCache()
    document = _document(source, title="Walk count")
    _render(tmp_path, cache, document, request_id=1)

    walked: list[int] = []
    original = prepared_module._feature_stamp

    def counting(feature: object) -> object:
        walked.append(id(feature))
        return original(feature)

    monkeypatch.setattr(prepared_module, "_feature_stamp", counting)
    # O-2: a warm repeat builds no record entry and no stage that takes a record.
    _result, diagnostics = _render(tmp_path, cache, document, request_id=2)
    metrics = _metrics(diagnostics)
    assert metrics["parsedSourceCacheHitCount"] == 1
    assert metrics["resolvedRecordCacheHitCount"] == 1
    assert metrics["preparedInputCacheMutationViolationCount"] == 0
    assert metrics["preparedInputCacheFullCheckCount"] == 0
    assert walked == []

    # A feature selection change builds the feature layers, which take the
    # record: one full walk per cached entry that holds records, at commit.
    _result, diagnostics = _render(
        tmp_path, cache, _document(source, title="Walk count", selected_features=("CDS",)), request_id=3
    )
    metrics = _metrics(diagnostics)
    assert metrics["renderStageBuild:feature_layers"] == 1
    assert metrics["preparedInputCacheMutationViolationCount"] == 0
    record_entries = [
        entry
        for entries in cache._layers.values()
        for entry in entries.values()
        if prepared_module._records_from_value(entry.value) is not None
    ]
    assert metrics["preparedInputCacheFullCheckCount"] == len(record_entries) > 0
    cached_features = sum(
        len(record.features)
        for entry in record_entries
        for record in prepared_module._records_from_value(entry.value) or ()
    )
    assert cached_features > 0
    assert len(walked) == cached_features


# Step 3a: the feature catalog revision (gbdraw/web/CLAUDE.md, Runtime data
# flow). A held revision equal to this render's gives featureCatalog None; the
# page then keeps the catalog it admitted for that revision.


def _rewrite_genbank(path: Path, records: tuple[SeqRecord, ...]) -> None:
    path.unlink()
    _write_genbank(path, records)


def test_a_held_revision_is_answered_unchanged_while_its_catalog_builds_are_reused(tmp_path: Path) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("held-linear"),))
    cache = PreparedBiologicalInputCache()
    cold, _ = _render(tmp_path, cache, _linear_document(source, title="Cold"), request_id=1)
    revision = cold["featureCatalogRevision"]
    assert isinstance(revision, str) and cold["featureCatalog"] is not None
    warm, _ = _render(tmp_path, cache, _linear_document(source, title="Cold"), request_id=2, held_revision=revision)
    retitled, _ = _render(
        tmp_path, cache, _linear_document(source, title="Render-only title"), request_id=3, held_revision=revision
    )
    unheld, _ = _render(tmp_path, cache, _linear_document(source, title="Render-only title"), request_id=4)

    for reply in (warm, retitled):
        assert reply["featureCatalog"] is None
        assert reply["featureCatalogRevision"] == revision
    assert retitled["results"] != cold["results"]
    assert unheld["featureCatalogRevision"] == revision
    assert unheld["featureCatalog"] == cold["featureCatalog"]


def test_a_feature_change_sends_a_new_revision_and_a_rebuilt_equal_input_is_sent_again(tmp_path: Path) -> None:
    """Inputs A, B (a feature change), A: three revisions; A's first revision held at the third is not reused."""

    source = tmp_path / "records.gbk"
    cache = PreparedBiologicalInputCache()
    _write_genbank(source, (_record("change-a"),))
    first, _ = _render(tmp_path, cache, _document(source), request_id=1, token_generation=1)
    _rewrite_genbank(source, (_record("change-a", product="changed product"),))
    changed, _ = _render(
        tmp_path, cache, _document(source), request_id=2, token_generation=2,
        held_revision=str(first["featureCatalogRevision"]),
    )
    _rewrite_genbank(source, (_record("change-a"),))
    again, _ = _render(
        tmp_path, cache, _document(source), request_id=3, token_generation=1,
        held_revision=str(first["featureCatalogRevision"]),
    )

    revisions = [reply["featureCatalogRevision"] for reply in (first, changed, again)]
    assert None not in revisions and len(set(revisions)) == 3
    assert changed["featureCatalog"] is not None and changed["featureCatalog"] != first["featureCatalog"]
    assert again["featureCatalog"] == first["featureCatalog"]


def test_two_caches_never_share_a_revision(tmp_path: Path) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("two-caches"),))
    document = _linear_document(source, title="Same")
    one, _ = _render(tmp_path, PreparedBiologicalInputCache(), document, request_id=1)
    two, _ = _render(
        tmp_path, PreparedBiologicalInputCache(), document, request_id=2,
        held_revision=str(one["featureCatalogRevision"]),
    )
    assert two["featureCatalogRevision"] != one["featureCatalogRevision"]
    assert two["featureCatalog"] == one["featureCatalog"]


def test_a_cache_without_a_random_source_always_sends_the_catalog(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def no_random_source(_size: int) -> str:
        raise NotImplementedError("no random source")

    monkeypatch.setattr(prepared_module.secrets, "token_hex", no_random_source)
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("no-epoch"),))
    cache = PreparedBiologicalInputCache()
    document = _linear_document(source, title="Same")
    cold, _ = _render(tmp_path, cache, document, request_id=1)
    warm, _ = _render(tmp_path, cache, document, request_id=2, held_revision="x:1")
    for reply in (cold, warm):
        assert reply["featureCatalogRevision"] is None
        assert reply["featureCatalog"] is not None


def test_a_linear_revision_held_for_a_circular_request_is_not_reused(tmp_path: Path) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("other-mode"),))
    cache = PreparedBiologicalInputCache()
    linear, _ = _render(tmp_path, cache, _linear_document(source, title="Linear"), request_id=1)
    circular, _ = _render(
        tmp_path, cache, _document(source), request_id=2, held_revision=str(linear["featureCatalogRevision"])
    )
    assert circular["featureCatalog"] is not None
    assert circular["featureCatalogRevision"] not in (None, linear["featureCatalogRevision"])


def test_a_failed_render_publishes_no_build_and_keeps_the_held_revision_current(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import gbdraw.web_support.request_render as web_request_render

    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("failed-render"),))
    cache = PreparedBiologicalInputCache()
    document = _linear_document(source, title="Same")
    cold, _ = _render(tmp_path, cache, document, request_id=1)
    revision = str(cold["featureCatalogRevision"])
    original = web_request_render.collect_legend_row_facts

    def fail_once(*args: object, **kwargs: object) -> object:
        monkeypatch.setattr(web_request_render, "collect_legend_row_facts", original)
        raise RuntimeError("render failed after the catalog stage")

    monkeypatch.setattr(web_request_render, "collect_legend_row_facts", fail_once)
    with pytest.raises(RuntimeError, match="after the catalog stage"):
        _render(tmp_path, cache, document, request_id=2, held_revision=revision)
    after, _ = _render(tmp_path, cache, document, request_id=3, held_revision=revision)
    assert after["featureCatalog"] is None
    assert after["featureCatalogRevision"] == revision


def _circular_batch_document(first: Path, second: Path) -> dict[str, object]:
    request = CircularBatchRequest(
        records=tuple(
            RecordInput(source=GenBankInputSource(path), record_key=f"batch-{index}")
            for index, path in enumerate((first, second), start=1)
        ),
        outputs=tuple(
            RenderOutputRequest(output_prefix=f"batch-{index}", formats=("svg",)) for index in (1, 2)
        ),
    )
    return build_session_document(request).to_dict()


def test_a_circular_batch_revision_names_each_result_and_one_changed_record_sends_the_whole_catalog(
    tmp_path: Path,
) -> None:
    first = tmp_path / "first.gbk"
    second = tmp_path / "second.gbk"
    _write_genbank(first, (_record("batch-a"),))
    _write_genbank(second, (_record("batch-b"),))
    cache = PreparedBiologicalInputCache()
    cold, _ = _render(tmp_path, cache, _circular_batch_document(first, second), request_id=1)
    revision = str(cold["featureCatalogRevision"])
    assert len(revision.split(":")[1].split(".")) == 2
    resource_ids = list(_circular_batch_document(first, second)["resources"])
    _rewrite_genbank(second, (_record("batch-b", product="changed product"),))
    changed, _ = _render(
        tmp_path, cache, _circular_batch_document(first, second), request_id=2,
        token_numbers={resource_ids[0]: 101, resource_ids[1]: 202}, held_revision=revision,
    )

    assert changed["featureCatalogRevision"] not in (None, revision)
    assert str(changed["featureCatalogRevision"]).split(":")[1].split(".")[0] == revision.split(":")[1].split(".")[0]
    catalog = changed["featureCatalog"]
    cold_catalog = cold["featureCatalog"]
    assert isinstance(catalog, dict) and isinstance(cold_catalog, dict)
    assert [item["resultName"] for item in catalog["items"]] == ["batch-1.svg", "batch-2.svg"]
    assert catalog["items"][0] == cold_catalog["items"][0]
    assert catalog["items"][1] != cold_catalog["items"][1]


def _cw04_document(source: Path, overrides: tuple[FeatureOverride, ...]) -> dict[str, object]:
    request = LinearDiagramRequest(
        records=tuple(
            RecordInput(source=GenBankInputSource(source), record_key=f"record-{index}") for index in (1, 2)
        ),
        # The Web default: Show Labels None, so only a label set On is drawn.
        options=LinearDiagramOptions(
            config_overrides={"labels.linear.scope": "none"}, feature_overrides=overrides
        ),
        output=RenderOutputRequest(output_prefix="cw04", formats=("svg",)),
    )
    return build_session_document(request).to_dict()


def test_cw04_each_generate_leaves_the_catalog_a_fresh_cache_gives(tmp_path: Path) -> None:
    """Py-9, the oracle of the CW-04 loop (tests/web/computation-ownership.playwright.spec.js).

    The three requests of CW-04 on two Linear copies of one GenBank file, with
    the browser's Show Labels None and the rows it commits (observed on this
    branch, 2026-10-15): no row; then ``label_visibility="on"`` on the
    tRNA-Phe feature of ``record-1`` (the feature popup); then only
    ``label_text="CW_BULK_LABEL"`` on that feature (the Label TSV import
    replaces the label edits, so the visibility row goes, and its
    ``* * label`` row reaches the one label shown). Outcome on this code: each
    edit builds the catalog item again (a new revision; the catalog is sent),
    and every catalog equals a fresh cache's. In the browser each edit's
    automatic label reflow is that render, so the Generate after it reuses the
    build; the CW-04 loop pins that.
    """

    source = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    plan = plan_request(LinearDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(source), record_key="record-1"),),
        output=RenderOutputRequest(output_prefix="cw04", formats=("svg",)),
    ))
    trna = next(
        entry.biological_feature_id
        for entry in plan.provenance[0].source_feature_catalog or ()
        if entry.feature_type == "tRNA"
    )
    steps: tuple[tuple[FeatureOverride, ...], ...] = (
        (),
        (FeatureOverride("record-1", trna, label_visibility="on"),),
        (FeatureOverride("record-1", trna, label_text="CW_BULK_LABEL"),),
    )
    cache = PreparedBiologicalInputCache()
    replies = []
    for index, rows in enumerate(steps, start=1):
        # No held revision: every reply carries its catalog, so content is compared.
        reply, _ = _render(tmp_path, cache, _cw04_document(source, rows), request_id=index)
        fresh, _ = _render(
            tmp_path, PreparedBiologicalInputCache(), _cw04_document(source, rows), request_id=10 + index
        )
        assert reply["featureCatalog"] is not None
        assert reply["featureCatalog"] == fresh["featureCatalog"]
        replies.append(reply)
    # The rule (CW-04): a reused revision means an equal catalog.
    for earlier, later in zip(replies, replies[1:]):
        if later["featureCatalogRevision"] == earlier["featureCatalogRevision"]:
            assert later["featureCatalog"] == earlier["featureCatalog"]
    # The outcome recorded above, which the CW-04 loop pins.
    revisions = [reply["featureCatalogRevision"] for reply in replies]
    assert None not in revisions and len(set(revisions)) == 3
    # Holding the latest revision, a repeat of the last request keeps the catalog.
    repeat, _ = _render(
        tmp_path, cache, _cw04_document(source, steps[-1]), request_id=4, held_revision=revisions[-1]
    )
    assert repeat["featureCatalog"] is None
    assert repeat["featureCatalogRevision"] == revisions[-1]


# Document pieces (STEP3B stage S): the reply names every children text and
# carries it only when the page does not hold its revision.


def _content_pieces(reply: dict[str, object]) -> list[dict[str, object]]:
    results = reply["results"]
    assert isinstance(results, list)
    return [
        piece
        for result in results
        for piece in result["document"]["pieces"]
        if piece["kind"] == "content"
    ]


def _joined_document(result: dict[str, object]) -> str:
    document = result["document"]
    assert isinstance(document, dict)
    texts = [
        piece["start"] + piece["children"] + "</g>" if piece["kind"] == "content" else piece["text"]
        for piece in document["pieces"]
    ]
    return document["head"] + "".join(texts) + document["tail"]


def test_the_reply_document_joins_to_the_written_svg(tmp_path: Path) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("pieces-record"),))
    reply, _ = _render(tmp_path, PreparedBiologicalInputCache(), _document(source, title="Pieces"), request_id=1)

    (result,) = reply["results"]
    assert result["content"] == '<?xml version="1.0" encoding="utf-8" ?>\n' + _joined_document(result)
    pieces = _content_pieces(reply)
    assert pieces
    for piece in pieces:
        data = piece["children"].encode("utf-8")
        assert piece["revision"] == hashlib.sha256(data).hexdigest()
        assert piece["bytes"] == len(data)
    assert {piece["kind"] for piece in result["document"]["pieces"]} == {"content", "raw"}


def test_held_revisions_are_named_without_their_text(tmp_path: Path) -> None:
    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("held-record"),))
    cache = PreparedBiologicalInputCache()
    cold, _ = _render(tmp_path, cache, _document(source, title="Held"), request_id=1)
    empty_held, _ = _render(tmp_path, cache, _document(source, title="Held"), request_id=2)
    assert empty_held["results"] == cold["results"]

    revisions = [piece["revision"] for piece in _content_pieces(cold)]
    held = frozenset(revisions[:1])
    partial, _ = _render(tmp_path, cache, _document(source, title="Held"), request_id=3, held_segments=held)
    for cold_piece, piece in zip(_content_pieces(cold), _content_pieces(partial)):
        assert {key: value for key, value in piece.items() if key != "children"} == {
            key: value for key, value in cold_piece.items() if key != "children"
        }
        assert ("children" in piece) == (piece["revision"] not in held)
    assert partial["results"][0]["content"] == cold["results"][0]["content"]


def test_a_title_position_change_sends_no_children_text(tmp_path: Path) -> None:
    """Moving the title changes start tags only, so every children text is held."""

    source = tmp_path / "records.gbk"
    _write_genbank(source, (_record("moved-title"),))
    cache = PreparedBiologicalInputCache()
    document = _document(source, title="Moved title")
    cold, _ = _render(tmp_path, cache, document, request_id=1)
    held = frozenset(piece["revision"] for piece in _content_pieces(cold))
    moved_document = copy.deepcopy(document)
    output = moved_document["renderRequest"]["diagramOptions"]["output"]
    assert output["plotTitlePosition"] == "top"
    output["plotTitlePosition"] = "bottom"
    moved, _ = _render(tmp_path, cache, moved_document, request_id=2, held_segments=held)

    assert moved["results"][0]["content"] != cold["results"][0]["content"]
    assert _content_pieces(moved)
    assert not [piece for piece in _content_pieces(moved) if "children" in piece]


def _gc_percent_calls(monkeypatch: pytest.MonkeyPatch) -> list[int]:
    import gbdraw.render.groups.circular.definition as definition_module

    calls: list[int] = []
    real = definition_module.calculate_gc_percent

    def counted(seq: object) -> float:
        calls.append(1)
        return real(seq)

    monkeypatch.setattr(definition_module, "calculate_gc_percent", counted)
    return calls


def test_a_circular_title_change_reuses_the_gc_percent_of_every_record(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The GC percent reads only the record: a warm Generate keys it and does not recompute it."""

    calls = _gc_percent_calls(monkeypatch)
    source = tmp_path / "records.gbk"
    records = (_record("gc-a"), _record("gc-b"), _record("gc-c"))
    records[2].seq = Seq("ATGC" * 180)
    _write_genbank(source, records)
    cache = PreparedBiologicalInputCache()

    cold, cold_diagnostics = _render(tmp_path, cache, _multi_record_document(source, title="Cold"), request_id=1)
    cold_metrics = _metrics(cold_diagnostics)
    assert cold_metrics["renderStageBuild:record_gc_percent"] == len(records)
    assert len(calls) == len(records)

    calls.clear()
    retitled, retitled_diagnostics = _render(
        tmp_path, cache, _multi_record_document(source, title="Warm"), request_id=2
    )
    warm_metrics = _metrics(retitled_diagnostics)
    assert calls == []
    assert warm_metrics["renderStageHit:record_gc_percent"] >= len(records)
    assert "renderStageBuild:record_gc_percent" not in warm_metrics
    assert warm_metrics.get("renderStageUncacheableCount", 0) == 0
    assert retitled["featureCatalog"] == cold["featureCatalog"]

    # Replaced records are other records: each GC percent is computed again.
    replaced, replaced_diagnostics = _render(
        tmp_path,
        cache,
        _multi_record_document(source, title="Warm"),
        request_id=3,
        token_generation=2,
    )
    assert _metrics(replaced_diagnostics)["renderStageBuild:record_gc_percent"] == len(records)
    assert len(calls) == len(records)
    assert replaced["results"] == retitled["results"]


def test_a_single_record_title_change_reuses_the_gc_percent(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    calls = _gc_percent_calls(monkeypatch)
    source = tmp_path / "record.gbk"
    _write_genbank(source, (_record("gc-single"),))
    cache = PreparedBiologicalInputCache()
    _render(tmp_path, cache, _document(source, title="Cold"), request_id=1)
    assert len(calls) == 1

    calls.clear()
    _warm, diagnostics = _render(tmp_path, cache, _document(source, title="Warm"), request_id=2)
    metrics = _metrics(diagnostics)
    assert calls == []
    assert metrics["renderStageHit:record_gc_percent"] >= 1
    assert metrics.get("renderStageUncacheableCount", 0) == 0
