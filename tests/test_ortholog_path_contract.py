"""S02 evidence and exhaustive small-input oracle; no production replacement.

Regenerate evidence: PYTHONPATH=. python tests/test_ortholog_path_contract.py --output FILE
The largest legacy enumeration is R=16 (16,384 paths). R>=24 is compact-only.
"""

from __future__ import annotations

from dataclasses import asdict, replace
import gzip
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import random
import subprocess
import sys

import pytest

from gbdraw.analysis import protein_colinearity as pc
from gbdraw.exceptions import ValidationError
from gbdraw.session_request_codec import (
    _read_typed_json_resource,
    encode_canonical_typed_resource,
)

ROOT = Path(__file__).resolve().parents[1]
legacy_spec = importlib.util.spec_from_file_location("s04_exhaustive", ROOT / "tests/prototypes/exhaustive_ortholog_paths.py")
legacy = importlib.util.module_from_spec(legacy_spec)
legacy_spec.loader.exec_module(legacy)
_build_ortholog_paths = legacy._build_ortholog_paths
prototype_spec = importlib.util.spec_from_file_location("s02_path_graph", ROOT / "tests/prototypes/ortholog_path_graph.py")
prototype = importlib.util.module_from_spec(prototype_spec)
prototype_spec.loader.exec_module(prototype)
DagPathIndex = prototype.DagPathIndex
spec = importlib.util.spec_from_file_location("s02_baseline", ROOT / "tools/benchmark_protein_comparison.py")
baseline = importlib.util.module_from_spec(spec)
spec.loader.exec_module(baseline)
SEED = 20260915


def proteins(n):
    return {f"p{i}": baseline.protein(pc, f"p{i}", i) for i in range(n)}


def edge(pm, u, v, kind="rbh", prior=None):
    return pc.OrthologEdge("og_1", "og_1", "og_1", u, v,
                          pm[u].record_index, pm[v].record_index,
                          kind, "block_anchor", prior, 90., 1e-30, 200., 100)


def compare(pm, edges):
    old_edges, old_paths = _build_ortholog_paths({"og_1": edges}, pm)
    production = pc.OrthologPathCollection.from_edges("og_1", edges, pm)
    assert production.count == len(old_paths["og_1"])
    assert tuple(production.iter_paths()) == old_paths["og_1"]
    assert tuple(replace(e, path_id=production.first_path_id(e)) for e in edges) == old_edges["og_1"]
    for rank, path in enumerate(old_paths["og_1"], 1):
        assert production.path_at(rank) == path
        assert production.rank_of(path.protein_ids) == rank
        assert production.path_by_id(path.path_id) == path
    index = DagPathIndex("og_1", edges, pm)
    assert index.materialized_paths == 0
    assert index.count == len(old_paths["og_1"])
    assert index.updated_edges() == old_edges["og_1"]
    # Equality covers dataclass type, path ID, protein/edge/shared tuples and order.
    actual = tuple(index.iter_paths())
    assert actual == old_paths["og_1"]
    for rank, path in enumerate(actual, 1):
        assert index.rank_of(path.protein_ids) == rank
        assert index.path_by_id(path.path_id) == path
    expected_counts = {pid: sum(pid in path.protein_ids for path in actual) for pid in index.children}
    assert index.containing_count == expected_counts
    oracle = {"edges": old_edges["og_1"], "paths": old_paths["og_1"]}
    return {
        "count": str(index.count),
        "oracleSha256": baseline.digest(baseline.json_bytes(baseline.canonical(oracle))),
        "inputSha256": baseline.digest(baseline.json_bytes(baseline.canonical((pm, tuple(edges))))),
        "index": index.summary(),
    }


def named_cases():
    pm = proteins(7)
    # Sort order deliberately differs from topological order and lexical IDs.
    pm = {pid: replace(p, record_index=0, start=(6-i)*100, end=(6-i)*100+90)
          for i, (pid, p) in enumerate(pm.items())}
    yield "empty", pm, []
    yield "excluded-and-missing", pm, [edge(pm, "p0", "p1", "same_record_inparalog", "keep"),
        replace(edge(pm, "p1", "p2"), subject_protein_id="missing", path_id="prior")]
    pairs = [("p0", "p2"), ("p1", "p2"), ("p2", "p3"), ("p2", "p4"), ("p0", "p4")]
    yield "multiple-starts-sinks-and-shared", pm, [edge(pm, u, v) for u, v in pairs]
    yield "parallel-protein-dedup", pm, [edge(pm, "p0", "p2", "rbh", "old-rbh"),
        edge(pm, "p0", "p2", "coortholog"), edge(pm, "p0", "p2", "coortholog"),
        edge(pm, "p2", "p3"), edge(pm, "p2", "p4"),
        edge(pm, "p5", "p6", "record_local_paralog", "keep-local")]
    yield "reverse-and-same-record-dag", pm, [edge(pm, "p5", "p2", "coortholog"),
        edge(pm, "p2", "p0"), edge(pm, "p0", "p4")]
    yield "disconnected", pm, [edge(pm, "p0", "p1"), edge(pm, "p2", "p3")]


def test_named_dag_metadata_equivalence():
    for _, pm, edges in named_cases():
        compare(pm, edges)
    _, pm, edges = next(c for c in named_cases() if c[0] == "parallel-protein-dedup")
    index = DagPathIndex("og_1", edges, pm)
    assert index.count == 2  # Six edge walks collapse to two protein sequences.
    assert index.updated_edges()[0].path_id == "old-rbh"  # Lexical coortholog edge wins.
    assert index.updated_edges()[1].path_id == index.updated_edges()[2].path_id


def all_small_dags():
    # All edge subsets of a five-node total order: 2^10 labeled DAGs.
    pm = proteins(5)
    pm = {pid: replace(p, record_index=(4-i) % 3, start=(i*71) % 113)
          for i, (pid, p) in enumerate(pm.items())}
    pairs = list(itertools.combinations(pm, 2))
    for mask in range(1 << len(pairs)):
        yield pm, [edge(pm, u, v, "coortholog" if i % 2 else "rbh")
                   for i, (u, v) in enumerate(pairs) if mask & (1 << i)]


def test_all_1024_small_dags_match_every_metadata_field():
    for pm, edges in all_small_dags():
        compare(pm, edges)


def randomized_dags():
    rng = random.Random(SEED)
    for _ in range(160):
        pm = proteins(8)
        pm = {pid: replace(p, record_index=rng.randrange(4), start=rng.randrange(3), end=900)
              for pid, p in pm.items()}
        order = list(pm)
        rng.shuffle(order)
        edges = []
        for u, v in itertools.combinations(order, 2):
            if rng.random() < .32:
                edges.append(edge(pm, u, v))
                if rng.random() < .4:
                    edges.extend([edge(pm, u, v, "coortholog")] * 2)
        rng.shuffle(edges)
        yield pm, edges


def test_seeded_dags_with_parallel_edges_and_coordinate_ties():
    for pm, edges in randomized_dags():
        compare(pm, edges)


def inference_cases():
    pm = proteins(6)
    # p4 -> p0 decreases record order; p1 -> p5 is incoming-only for p5.
    rows = [baseline.hit("p0", "p1"), baseline.hit("p1", "p0"),
            baseline.hit("p4", "p0", bitscore=160), baseline.hit("p1", "p5", bitscore=160)]
    yield "reverse-leaf-and-incoming-only", pm, rows
    pm = {pid: replace(p, record_index=i//2) for i, (pid, p) in enumerate(proteins(6).items())}
    rows = [baseline.hit("p0", "p2"), baseline.hit("p2", "p0"),
            baseline.hit("p0", "p1", bitscore=180), baseline.hit("p1", "p0", bitscore=180),
            baseline.hit("p4", "p5"), baseline.hit("p5", "p4")]
    yield "same-record-and-record-local", pm, rows
    rows = [baseline.hit(u, v) for u in ("p0", "p1") for v in ("p2", "p3")]
    rows += [baseline.hit(row["subject"], row["query"]) for row in rows]
    yield "near-reciprocal-core-and-duplicate-hsps", pm, rows + rows


def infer(pm, rows):
    by_pair = {}
    for row in rows:
        pair = (pm[row["query"]].record_index, pm[row["subject"]].record_index)
        by_pair.setdefault(pair, []).append(row)
    return pc.select_rbh_orthogroup_edges_from_directional_hits(
        {pair: baseline.frame(values) for pair, values in by_pair.items()}, pm,
        record_count=max(p.record_index for p in pm.values())+1,
        comparison_pairs=(), max_related_edges_per_orthogroup=2, path_representation="exhaustive",
    ).orthogroups


def check_inferred(pm, result):
    observed = []
    for group, edges in result.ortholog_edges_by_orthogroup_id.items():
        index = DagPathIndex(group, edges, pm)
        assert index.updated_edges() == edges
        assert tuple(index.iter_paths()) == result.ortholog_paths_by_orthogroup_id[group]
        path_edges = [e for e in edges if e.edge_kind in {"rbh", "coortholog"}]
        assert all(e.query_record_index != e.subject_record_index for e in path_edges)
        assert len({(e.query_protein_id, e.subject_protein_id) for e in path_edges}) == len(path_edges)
        observed.extend(asdict(e) for e in edges)
    return observed


def test_reachable_inference_graphs_and_nonmonotonic_record_order():
    observations = {}
    for name, pm, rows in inference_cases():
        result = infer(pm, rows)
        observations[name] = check_inferred(pm, result)
    reverse = observations["reverse-leaf-and-incoming-only"]
    assert any(e["query_record_index"] > e["subject_record_index"] and e["edge_kind"] == "coortholog" for e in reverse)
    assert any(e["query_protein_id"] == "p1" and e["subject_protein_id"] == "p5" for e in reverse)
    kinds = {e["edge_kind"] for e in observations["same-record-and-record-local"]}
    assert {"same_record_inparalog", "record_local_paralog"} <= kinds
    core = observations["near-reciprocal-core-and-duplicate-hsps"]
    assert any(e["edge_kind"] == "coortholog" for e in core)


def test_seeded_reachable_directional_evidence():
    rng = random.Random(SEED)
    for _ in range(60):
        pm = {pid: replace(p, record_index=i//2) for i, (pid, p) in enumerate(proteins(8).items())}
        rows = [baseline.hit(u, v, bitscore=rng.choice([90, 170, 200]))
                for u in pm for v in pm if rng.random() < .5]
        rng.shuffle(rows)
        check_inferred(pm, infer(pm, rows))


@pytest.mark.parametrize("size", [8, 12, 16])
def test_s01_full_selector_paths_and_archived_oracle(size):
    pm, tables = baseline.synthetic(pc, f"path-{size}", SEED)
    result = pc.select_rbh_orthogroup_edges_from_directional_hits(tables, pm, record_count=size, path_representation="exhaustive")
    check_inferred(pm, result.orthogroups)
    assert len(result.orthogroups.ortholog_paths_by_orthogroup_id["og_1"]) == 2**(size-2)
    archive = ROOT / f"docs/internal/collinear_similarity_performance_plan_2026-09-15/results/data/oracles/path-{size}.selector.json.gz"
    assert baseline.json_bytes(baseline.canonical(result)) == gzip.decompress(archive.read_bytes())


def complete_graph(size, extra=False):
    pm = proteins(size + (2 if extra else 0))
    edges = [edge(pm, f"p{i}", f"p{j}") for i in range(size) for j in range(i+1, size)]
    if extra:
        edges.append(edge(pm, f"p{size}", f"p{size+1}"))
    return DagPathIndex("og_1", edges, pm)


@pytest.mark.parametrize("size,extra", [(24, False), (32, False), (56, False), (56, True)])
def test_large_compact_count_rank_shared_without_legacy_enumeration(size, extra):
    index = complete_graph(size, extra)
    assert index.count == 2**(size-2) + int(extra)
    assert index.containing_count["p0"] == index.containing_count[f"p{size-1}"] == 2**(size-2)
    assert all(index.containing_count[f"p{i}"] == 2**(size-3) for i in range(1, size-1))
    assert index.materialized_paths == 0
    for rank in (1, 2, index.count//2, index.count-1, index.count):
        assert index.rank_of(index.path_at(rank).protein_ids) == rank
    # Edge first-path metadata is exact even when ranks exceed 2^53.
    for (u, v), eid in index.chosen.items():
        assert eid in index.path_at(index.first_edge_rank[eid]).edge_ids
    assert index.materialized_paths == 5 + len(index.chosen)


def test_cycle_observations_do_not_become_a_dag_fallback():
    pm = proteins(6)
    cycle = [edge(pm, "p0", "p1"), edge(pm, "p1", "p0")]
    _, paths = _build_ortholog_paths({"og_1": cycle}, pm)
    assert [p.protein_ids for p in paths["og_1"]] == [("p0", "p1"), ("p1", "p0")]
    # A cyclic component is not visited if another component has a source.
    _, paths = _build_ortholog_paths({"og_1": cycle + [edge(pm, "p2", "p3")]}, pm)
    assert [p.protein_ids for p in paths["og_1"]] == [("p2", "p3")]
    prefix_cycle = [edge(pm, "p4", "p0"), *cycle, edge(pm, "p1", "p5")]
    _, paths = _build_ortholog_paths({"og_1": prefix_cycle}, pm)
    assert {p.protein_ids for p in paths["og_1"]} == {("p4", "p0", "p1"), ("p4", "p0", "p1", "p5")}
    _, paths = _build_ortholog_paths({"og_1": [edge(pm, "p0", "p0")]}, pm)
    assert paths["og_1"] == ()
    for edges in (cycle, prefix_cycle, [edge(pm, "p0", "p0")]):
        with pytest.raises(ValidationError, match="cycle"):
            DagPathIndex("og_1", edges, pm)


def test_legacy_typed_tuple_accepts_explicit_cyclic_noncanonical_corpus(tmp_path):
    pm = proteins(2)
    edges = (edge(pm, "p0", "p1", prior="external"), edge(pm, "p1", "p0"))
    paths = (pc.OrthologPath("og_1", "user-path", ("p1", "p0"), (pc._edge_id(edges[1]),), ("p0",)),)
    legacy = pc.OrthogroupResult(orthogroups={}, member_by_protein_id={},
        ortholog_edges_by_orthogroup_id={"og_1": edges}, ortholog_paths_by_orthogroup_id={"og_1": paths})
    path = tmp_path / "typed.json"
    path.write_bytes(encode_canonical_typed_resource("orthogroupResult", legacy))
    decoded = _read_typed_json_resource("r", value_kind="orthogroupResult", expected=pc.OrthogroupResult | pc.OrthogroupGraphResult,
        path="comparison", resource_paths={"r": path})
    decoded = pc.materialize_ortholog_paths(decoded)
    assert decoded == legacy
    assert type(decoded.ortholog_paths_by_orthogroup_id["og_1"]) is tuple


def test_explicit_access_errors():
    index = complete_graph(4)
    for value in (True, 1., "1", None):
        with pytest.raises(TypeError):
            index.path_at(value)
    for rank in (-1, 0, index.count+1):
        with pytest.raises(IndexError):
            index.path_at(rank)
    for value in ("og_2.path_1", "og_1.path_01", "og_1.path_0", "og_1.path_1x"):
        with pytest.raises(ValidationError):
            index.path_by_id(value)
    for value in (("p0", "p1"), ("p1", "p3"), ("p0", "bad"), ()):
        with pytest.raises(ValidationError):
            index.rank_of(value)


def test_javascript_decimal_transport():
    values = [0, 1, 2**53-1, 2**53, 2**53+1, 2**54+1, 2**100+7]
    script = """
const fs = require('node:fs');
const assert = require('node:assert/strict');
const values = JSON.parse(fs.readFileSync(0, 'utf8'));
for (const v of values) {
  assert.match(v, /^(0|[1-9][0-9]*)$/);
  assert.equal(BigInt(v).toString(), v);
  assert.equal(JSON.parse(JSON.stringify({count: v})).count, v);
}
assert.notEqual(BigInt(Number(values[5])).toString(), values[5]);
process.stdout.write(JSON.stringify(values));
"""
    payload = json.dumps(list(map(str, values)))
    result = subprocess.run(["node", "-e", script], input=payload, text=True, capture_output=True, check=True)
    assert json.loads(result.stdout) == json.loads(payload)


def test_existing_catalog_and_popup_accept_exact_decimal_strings(tmp_path):
    from gbdraw.web_support.feature_catalog import select_feature_catalog_item

    source = ROOT / "gbdraw/web/gallery/sessions/hepatoplasmataceae_orthogroup.gbdraw-session.json.gz"
    session = json.loads(gzip.decompress(source.read_bytes()))
    catalog = session["editorState"]["featureCatalog"]
    exact = str(2**54+1)
    catalog["items"][0]["orthogroups"][0]["orthologPathCount"] = exact
    item = select_feature_catalog_item(catalog, result_index=0, result_name=session["results"][0]["name"])
    assert item["orthogroups"][0]["orthologPathCount"] == exact
    payload = tmp_path / "catalog.json"
    payload.write_text(json.dumps({"catalog": catalog, "results": session["results"]}), encoding="utf-8")
    script = """
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {pathToFileURL} from 'node:url';
const root = process.argv[1];
const {admitFeatureCatalog} = await import(pathToFileURL(root + '/gbdraw/web/js/services/feature-catalog.js'));
const {buildMatchPopupPayload} = await import(pathToFileURL(root + '/gbdraw/web/js/app/pairwise-match-popup.js'));
const data = JSON.parse(readFileSync(process.argv[2], 'utf8'));
const admission = admitFeatureCatalog(data.catalog, data.results, {mode: 'linear'});
const groups = admission.featureState.orthogroups;
const group = groups[0];
assert.equal(group.orthologPathCount, '18014398509481985');
for (const kind of ['orthogroup', 'collinear']) {
  const attrs = {'data-match-kind': kind, 'data-orthogroup-id': group.id};
  const popup = buildMatchPopupPayload({getAttribute: name => attrs[name] ?? ''}, {orthogroups: groups});
  assert.ok(JSON.stringify(popup).includes('18014398509481985'));
}
"""
    subprocess.run(["node", "--input-type=module", "-e", script, str(ROOT), str(payload)], check=True)


def write_evidence(output):
    records = {name: compare(pm, edges) for name, pm, edges in named_cases()}
    for name, pm, edges in named_cases():
        old_edges, old_paths = _build_ortholog_paths({"og_1": edges}, pm)
        records[name]["input"] = {"proteins": {k: asdict(v) for k, v in pm.items()}, "edges": [asdict(e) for e in edges]}
        records[name]["expected"] = {"edges": [asdict(e) for e in old_edges["og_1"]], "paths": [asdict(p) for p in old_paths["og_1"]]}
    exhaustive = [compare(pm, edges) for pm, edges in all_small_dags()]
    seeded = [compare(pm, edges) for pm, edges in randomized_dags()]
    inferred = {name: check_inferred(pm, infer(pm, rows)) for name, pm, rows in inference_cases()}
    large = {f"R{size}" + ("+1" if extra else ""): complete_graph(size, extra).summary()
             for size, extra in [(24, False), (32, False), (56, False), (56, True)]}
    sources = [Path(pc.__file__).resolve(), Path(__file__).resolve(), ROOT / "tests/prototypes/ortholog_path_graph.py"]
    report = {
        "kind": "S02-test-only-path-evidence", "seed": SEED,
        "sourceHead": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "python": sys.version, "node": subprocess.check_output(["node", "--version"], text=True).strip(),
        "sourceHashes": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        "named": records, "exhaustiveSmallDags": exhaustive, "seededDags": seeded,
        "reachableInference": inferred, "compactOnly": large,
        "limitations": "No production speedup claim. No legacy enumeration at R>=24. See pytest for additional cyclic/transport/S01 checks.",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(gzip.compress(json.dumps(report, sort_keys=True, separators=(",", ":")).encode(), mtime=0))


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    write_evidence(parser.parse_args().output)
