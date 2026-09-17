"""Production PATH-B boundaries: compact default, exact legacy retrieval/wire."""

from dataclasses import replace
import importlib.util
import json
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from gbdraw.analysis import protein_colinearity as pc
from gbdraw.analysis.ortholog_paths import OrthologPathCollection
from gbdraw.exceptions import ValidationError
from gbdraw.session_request_codec import _read_typed_json_resource, encode_canonical_typed_resource
from gbdraw.web_support.orthogroup_metadata import serialize_orthogroups_payload

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("s06_runner", ROOT / "tools/benchmark_protein_comparison.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


def edge(pm, u, v, kind="rbh", prior=None):
    return pc.OrthologEdge("og_1", "og_1", "og_1", u, v,
        pm[u].record_index, pm[v].record_index, kind, "block_anchor", prior, 90., 1e-30, 200., 100)


def complete(size, extra=False):
    pm = {f"p{i}": runner.protein(pc, f"p{i}", i) for i in range(size + 2 * extra)}
    edges = [edge(pm, f"p{i}", f"p{j}") for i in range(size) for j in range(i+1, size)]
    if extra:
        edges.append(edge(pm, f"p{size}", f"p{size+1}"))
    updated, indexes = pc._build_ortholog_path_indexes({"og_1": edges}, pm)
    return pc.OrthogroupGraphResult({}, {}, ortholog_edges_by_orthogroup_id=updated,
                                   path_indexes_by_orthogroup_id=indexes)


def roundtrip(result, tmp_path):
    payload = encode_canonical_typed_resource("orthogroupResult", result)
    path = tmp_path / "paths.json"
    path.write_bytes(payload)
    decoded = _read_typed_json_resource("p", value_kind="orthogroupResult",
        expected=pc.OrthogroupResult | pc.OrthogroupGraphResult,
        path="comparison", resource_paths={"p": path})
    return payload, decoded


@pytest.mark.parametrize("size", [8, 12, 16, 24])
def test_default_inference_and_consumers_never_expand(size, tmp_path):
    pm, tables = runner.synthetic(pc, f"path-{size}", runner.SEED)
    with patch.object(OrthologPathCollection, "iter_paths", side_effect=AssertionError("enumerated")) as iterator, \
         patch.object(pc, "materialize_ortholog_paths", side_effect=AssertionError("materialized")) as adapter, \
         patch.object(OrthologPathCollection, "_path", side_effect=AssertionError("path object")) as constructor:
        result = pc.select_rbh_orthogroup_edges_from_directional_hits(tables, pm, record_count=size)
        assert type(result.orthogroups) is pc.OrthogroupGraphResult
        assert not hasattr(result.orthogroups, "ortholog_paths_by_orthogroup_id")
        index = result.orthogroups.path_indexes_by_orthogroup_id["og_1"]
        assert index.count == 2**(size-2)
        groups = serialize_orthogroups_payload(result.orthogroups)
        assert groups[0]["orthologPathCount"] == str(index.count)
        assert "orthologPaths" not in groups[0]
        payload, decoded = roundtrip(result.orthogroups, tmp_path)
        assert decoded == result.orthogroups
        assert b'orthologPathsByOrthogroupId' not in payload
        iterator.assert_not_called()
        adapter.assert_not_called()
        constructor.assert_not_called()


@pytest.mark.parametrize("size,extra", [(24, False), (32, False), (56, True)])
def test_production_big_integer_roundtrip_and_random_access(size, extra, tmp_path):
    result = complete(size, extra)
    wire, decoded = roundtrip(result, tmp_path)
    index = decoded.path_indexes_by_orthogroup_id["og_1"]
    assert index.count == 2**(size-2) + extra
    assert f'"count":"{index.count}"'.encode() in wire
    for rank in (1, 2, index.count//2, index.count-1, index.count):
        path = index.path_at(rank)
        assert index.rank_of(path.protein_ids) == rank
        assert index.path_by_id(path.path_id) == path
    for e in decoded.ortholog_edges_by_orthogroup_id["og_1"]:
        assert pc._edge_id(e) in index.path_by_id(e.path_id).edge_ids
    assert index.containing_count("p0") == 2**(size-2)
    if extra:
        assert index.count == 18014398509481985


def test_arbitrary_legacy_corpus_keeps_content_and_never_completes_closure(tmp_path):
    pm = {f"p{i}": runner.protein(pc, f"p{i}", i) for i in range(4)}
    edges = (edge(pm, "p0", "p1", prior="external"), edge(pm, "p1", "p0"), edge(pm, "p0", "p2"))
    paths = (pc.OrthologPath("other-group", "custom", ("p1", "p0", "p1"), ("e1", "e2"), ("p1",)),
             pc.OrthologPath("og_1", "custom-2", ("p0", "p2"), ("e3",), ()))
    old = pc.OrthogroupResult({}, {}, ortholog_edges_by_orthogroup_id={"og_1": edges},
        ortholog_paths_by_orthogroup_id={"og_1": paths, "empty": ()})
    _, decoded = roundtrip(old, tmp_path)
    assert type(decoded) is pc.OrthogroupGraphResult
    assert pc.materialize_ortholog_paths(decoded) == old
    index = decoded.path_indexes_by_orthogroup_id["og_1"]
    assert index.count == 2
    assert index.path_by_id("custom") == paths[0]
    assert index.rank_of(paths[0].protein_ids) == 1
    assert index.first_path_id(edges[0]) == "external"
    assert pc.materialize_ortholog_paths(old) is old


def test_dag_payload_integrity_rejects_cycles_wrong_counts_and_order():
    original = complete(4).path_indexes_by_orthogroup_id["og_1"].to_payload()
    for key, value in (("count", 4), ("count", "04"), ("count", "5"), ("ordering", "other")):
        with pytest.raises(ValidationError):
            OrthologPathCollection.from_payload(original | {key: value})
    with pytest.raises(ValidationError):
        OrthologPathCollection.from_payload(original | {"nodes": original["nodes"][::-1]})
    pm = {f"p{i}": runner.protein(pc, f"p{i}", i) for i in range(4)}
    for edges in ([edge(pm, "p0", "p0")],
                  [edge(pm, "p0", "p1"), edge(pm, "p1", "p0"), edge(pm, "p2", "p3")]):
        with pytest.raises(ValidationError, match="cycle"):
            OrthologPathCollection.from_edges("og_1", edges, pm)
    result = complete(4)
    edges = result.ortholog_edges_by_orthogroup_id["og_1"]
    with pytest.raises(ValidationError, match="first-path"):
        replace(result, ortholog_edges_by_orthogroup_id={"og_1": (replace(edges[0], path_id="wrong"), *edges[1:])})


def test_access_contract_errors_and_ambiguous_explicit_values():
    index = complete(4).path_indexes_by_orthogroup_id["og_1"]
    for rank in (True, 1., "1", None):
        with pytest.raises(TypeError):
            index.path_at(rank)
    for rank in (0, -1, 5):
        with pytest.raises(IndexError):
            index.path_at(rank)
    for path_id in ("og_2.path_1", "og_1.path_01", "og_1.path_0"):
        with pytest.raises(ValidationError):
            index.path_by_id(path_id)
    with pytest.raises(KeyError):
        index.containing_count("absent")
    path = index.path_at(1)
    explicit = OrthologPathCollection("og_1", "explicit", paths=(path, path))
    assert tuple(explicit.iter_paths()) == (path, path)
    with pytest.raises(ValidationError, match="Ambiguous"):
        explicit.path_by_id(path.path_id)
    with pytest.raises(ValidationError, match="ambiguous"):
        explicit.rank_of(path.protein_ids)
    with pytest.raises(KeyError):
        explicit.path_by_id("missing")


def test_published_typed_v1_v2_converge_once(tmp_path):
    # Actual v1/v2 tagged legacy format; current writer intentionally emits v3.
    old = pc.OrthogroupResult({}, {}, ortholog_paths_by_orthogroup_id={
        "x": (pc.OrthologPath("x", "noncanonical", ("a", "b"), ("e",)),)})
    def encode(value):
        from dataclasses import fields, is_dataclass
        from gbdraw.session_request_codec import _camel
        if is_dataclass(value):
            return {"type": type(value).__name__, "fields": {_camel(f.name): encode(getattr(value, f.name)) for f in fields(value)}}
        if isinstance(value, dict):
            return {k: encode(v) for k, v in value.items()}
        if isinstance(value, (tuple, list)):
            return [encode(v) for v in value]
        return value
    for schema in (1, 2):
        path = tmp_path / f"v{schema}.json"
        path.write_text(json.dumps({"schema": schema, "kind": "orthogroupResult", "value": encode(old)}))
        decoded = _read_typed_json_resource("p", value_kind="orthogroupResult",
            expected=pc.OrthogroupResult | pc.OrthogroupGraphResult,
            path="comparison", resource_paths={"p": path})
        assert pc.materialize_ortholog_paths(decoded) == old
        assert json.loads(encode_canonical_typed_resource("orthogroupResult", decoded))["schema"] == 3


def test_alias_keys_and_equal_protein_order_keys_preserve_legacy_order(tmp_path):
    import random
    from tests.prototypes.exhaustive_ortholog_paths import _build_ortholog_paths

    rng = random.Random(6006)
    for case in range(120):
        pm = {f"alias-{i}": replace(runner.protein(pc, f"internal-{i // 3}", i // 3),
                                     start=0, end=300) for i in range(8)}
        edges = [edge(pm, f"alias-{i}", f"alias-{j}") for i in range(8) for j in range(i+1, 8)
                 if rng.random() < .35]
        expected_edges, expected_paths = _build_ortholog_paths({"og_1": edges}, pm)
        updated, indexes = pc._build_ortholog_path_indexes({"og_1": edges}, pm)
        index = indexes["og_1"]
        assert updated == expected_edges
        assert tuple(index.iter_paths()) == expected_paths["og_1"]
        for rank, path in enumerate(expected_paths["og_1"], 1):
            assert index.path_at(rank) == path
            assert index.rank_of(path.protein_ids) == rank
        assert OrthologPathCollection.from_payload(index.to_payload()) == index

    # Exercise the public producer too; aliases are caller input, not rewritten
    # or normalized. Only the removed path phase is replaced by the old oracle.
    pm, tables = runner.synthetic(pc, "path-4", runner.SEED)
    pm = {pid: replace(p, protein_id="caller-internal") for pid, p in pm.items()}
    def old_path_boundary(edges, proteins):
        updated, paths = _build_ortholog_paths(edges, proteins)
        return updated, {g: OrthologPathCollection(g, "explicit", paths=ps) for g, ps in paths.items()}
    actual = pc.select_rbh_orthogroup_edges_from_directional_hits(
        tables, pm, record_count=4, path_representation="exhaustive")
    with patch.object(pc, "_build_ortholog_path_indexes", old_path_boundary):
        expected = pc.select_rbh_orthogroup_edges_from_directional_hits(
            tables, pm, record_count=4, path_representation="exhaustive")
    assert runner.canonical(actual) == runner.canonical(expected)


@pytest.mark.parametrize("name", ["gallery-collinear", "gallery-orthogroup"])
def test_gallery_science_matches_archived_s04_bytes(name):
    import gzip
    from gbdraw.analysis import collinearity as cc
    # Reproduce the archived pandas 2 inference setting, including its dtypes.
    with pd.option_context("future.infer_string", False):
        _, stages = runner.build_case(ROOT, pc, cc, name, runner.SEED)
        result = stages["post_search"]()
    legacy = replace(result, orthogroups=pc.materialize_ortholog_paths(result.orthogroups))
    saved = json.loads(gzip.decompress((ROOT / 'docs/internal/collinear_similarity_performance_plan_2026-09-15/results/data'
        / f's04-current-timing-cursor-{name}.json.gz').read_bytes()))
    expected = saved['cases'][name]['stages']['post_search']['semanticSha256']
    assert runner.digest(runner.json_bytes(runner.canonical(legacy))) == expected


@pytest.mark.parametrize("name", ["gallery-collinear", "gallery-orthogroup"])
def test_current_gallery_render_save_replay_do_not_enumerate(name, tmp_path):
    from gbdraw.analysis import collinearity as cc
    from gbdraw.api import (load_session_document, materialize_session, session_to_request,
                            save_session_document, render_request)
    _, stages = runner.build_case(ROOT, pc, cc, name, runner.SEED)
    result = stages["post_search"]()
    gallery = runner.GALLERIES[0 if name == "gallery-collinear" else 1]
    document = load_session_document(ROOT / f'gbdraw/web/gallery/sessions/{gallery}.gbdraw-session.json.gz')
    with materialize_session(document, output_directory=tmp_path) as materialized:
        request = session_to_request(materialized)
        options = replace(request.options, orthogroups=result.orthogroups,
                          collinearity_blocks=result if name == "gallery-collinear" else ())
        request = replace(request, options=options)
        with patch.object(OrthologPathCollection, 'iter_paths', side_effect=AssertionError('enumerated')), \
             patch.object(OrthologPathCollection, '_path', side_effect=AssertionError('path object')), \
             patch.object(pc, 'materialize_ortholog_paths', side_effect=AssertionError('materialized')):
            saved = tmp_path / 'compact-session.json'
            save_session_document(saved, request)
            with materialize_session(load_session_document(saved), output_directory=tmp_path / 'replay') as replay:
                decoded = session_to_request(replay)
                assert isinstance(decoded.options.orthogroups, pc.OrthogroupGraphResult)
                assert all(index.kind == 'dag' for index in decoded.options.orthogroups.path_indexes_by_orthogroup_id.values())
                rendered = render_request(decoded)
                assert rendered.output_paths and all(path.is_file() for path in rendered.output_paths)


def test_legacy_reference_promotion_rebuilds_only_declared_collection_data():
    from gbdraw.api.session_compat import rewrite_protein_artifact_references
    paths = (pc.OrthologPath('og_1', 'custom', ('p_r_old_0_9_1_deadbeefdead', 'p_r_other_10_19_1_cafebabecafe'),
                             ('p_r_old_0_9_1_deadbeefdead->p_r_other_10_19_1_cafebabecafe',), ('p_r_old_0_9_1_deadbeefdead',)),)
    value = pc.compact_ortholog_paths(pc.OrthogroupResult({}, {}, ortholog_paths_by_orthogroup_id={'og_1': paths}))
    mapping = {'p_r_old_0_9_1_deadbeefdead': 'h_first', 'p_r_other_10_19_1_cafebabecafe': 'h_second'}
    migrated = rewrite_protein_artifact_references(value, mapping)
    index = migrated.path_indexes_by_orthogroup_id['og_1']
    assert index.path_at(1) == pc.OrthologPath('og_1', 'custom', ('h_first', 'h_second'),
                                             ('h_first->h_second',), ('h_first',))
    assert index.containing_count('h_first') == 1
    assert 'p_r_old_0_9_1_deadbeefdead' not in index._containing
    # Current DAGs have current references; detached copies rebuild all DP state.
    current = complete(4)
    assert rewrite_protein_artifact_references(current, {}) == current
