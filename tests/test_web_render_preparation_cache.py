from __future__ import annotations

import base64
import copy
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
)
from gbdraw.annotations import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    RegionAnnotation,
)
from gbdraw.exceptions import ValidationError
import gbdraw.api.request_render as request_render_module
from gbdraw.api.request_render import normalize_request_records
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
    cache: PreparedBiologicalInputCache,
    document: dict[str, object],
    *,
    request_id: int,
    token_generation: int = 1,
    token_numbers: dict[str, int] | None = None,
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
        cache_token = f"render-resource-{token_number}"
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
        _resource_identities=identities,
    )
    return result, diagnostics


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
    assert cold_metrics["interactiveContextCacheMissCount"] == 1
    assert cold_metrics["interactiveContextBuildCount"] == 1
    assert cold_metrics["interactiveFeatureTraversalCount"] == 1
    assert cold_metrics["selectorSafetyScopeBuildCount"] == 0
    assert cold_metrics["preparedInputCacheRetainedBytes"] > 0

    for diagnostics in (warm_diagnostics, render_only_diagnostics):
        metrics = _metrics(diagnostics)
        assert metrics["parsedSourceCacheHitCount"] == 1
        assert metrics["parsedSourceParseCount"] == 0
        assert metrics["resolvedRecordCacheHitCount"] == 1
        assert metrics["resolvedRecordBuildCount"] == 0
        assert metrics["interactiveContextCacheHitCount"] == 1
        assert metrics["interactiveContextBuildCount"] == 0
        assert metrics["interactiveFeatureTraversalCount"] == 0
        assert metrics["selectorSafetyScopeBuildCount"] == 0
        assert metrics["preparedInputCacheRetainedBytes"] > 0
        assert metrics["preparedInputCacheMutationViolationCount"] == 0

    assert warm["results"] == cold["results"]
    assert render_only["results"] != cold["results"]
    assert _metrics(warm_diagnostics)["featureCatalogSvgParseCount"] == 1
    assert _metrics(render_only_diagnostics)["featureCatalogSvgParseCount"] == 1
    assert float(render_only_diagnostics["timingsMs"]["drawing"]) > 0
    assert float(render_only_diagnostics["timingsMs"]["svgWrite"]) > 0


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

    assert _metrics(warm_diagnostics)["interactiveContextCacheHitCount"] == 1
    changed_metrics = _metrics(changed_diagnostics)
    assert changed_metrics["parsedSourceCacheHitCount"] == 1
    assert changed_metrics["resolvedRecordCacheHitCount"] == 1
    assert changed_metrics["interactiveContextCacheMissCount"] == 1
    assert changed_metrics["interactiveContextBuildCount"] == 1
    assert changed["metadata"] != first["metadata"]


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
    assert metrics["interactiveContextCacheMissCount"] == 1
    assert metrics["interactiveContextBuildCount"] == 1
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
    assert metrics["interactiveContextCacheMissCount"] == 1
    assert metrics["interactiveContextBuildCount"] == 1
    assert changed["metadata"] != first["metadata"]
    changed_payload = json.dumps(changed, sort_keys=True)
    assert "OG-B" in changed_payload
    assert "OG-A" not in changed_payload


def test_gff_fasta_policy_is_part_of_the_parsed_source_key(tmp_path: Path) -> None:
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

    def request(selected_features: tuple[str, ...]) -> CircularDiagramRequest:
        return CircularDiagramRequest(
            records=(
                RecordInput(
                    source=GffFastaInputSource(gff_path, fasta_path),
                    record_key="record-1",
                ),
            ),
            options=CircularDiagramOptions(
                selected_features_set=selected_features,
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
    assert _metrics(changed_diagnostics)["parsedSourceCacheMissCount"] == 1
    assert _metrics(changed_diagnostics)["parsedSourceParseCount"] == 1
    assert _metrics(changed_diagnostics)["resolvedRecordBuildCount"] == 1


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
    assert metrics["interactiveContextCacheMissCount"] == 1
    assert metrics["interactiveContextBuildCount"] == 1
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
    assert metrics["interactiveContextCacheMissCount"] == 1
    assert metrics["interactiveContextBuildCount"] == 1
    assert metrics["preparedInputCacheEvictionCount"] >= 3
    assert result != original
    replacement_payload = json.dumps(result, sort_keys=True)
    assert "replacement" in replacement_payload
    assert "cache target" not in replacement_payload


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
            "cacheToken": f"render-resource-{index}",
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
