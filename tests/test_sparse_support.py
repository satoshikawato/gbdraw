"""S04 differential support, snapshot, and first-match metadata contracts."""
from dataclasses import replace
import importlib.util
from pathlib import Path
import random
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from gbdraw.analysis import protein_colinearity as pc, collinearity as cc
from gbdraw.exceptions import ParseError
from gbdraw.web_support.orthogroup_metadata import serialize_orthogroups_payload
from tests.prototypes import sparse_support as old

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("s04_benchmark", ROOT / "tools/benchmark_protein_comparison.py")
bench = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bench)


def row(q, s, score=1.0, coverage=1.0, **kwargs):
    return SimpleNamespace(**(bench.hit(q, s) | {"normalized_score": score, "min_coverage": coverage} | kwargs))


def candidate(pm, best, members, thresholds=None):
    ids = {pid: "og_1" for pid in members}
    indexed = pc._index_core_support_evidence(best, ids).get("u", {}).get("og_1", ())
    actual = pc._build_core_support_candidate("u", "og_1", indexed, thresholds or {}, pm)
    expected = old._build_core_support_candidate("u", "og_1", members, best, thresholds or {}, pm)
    assert actual == expected
    return actual


@pytest.mark.parametrize("direction", ["incoming", "outgoing", "both"])
@pytest.mark.parametrize("record", [0, 1])
def test_direction_and_membership_diagnostic_priority(direction, record):
    pm = {pid: bench.protein(pc, pid, r) for pid, r in [("u", 0), ("m", record), ("domain", record)]}
    pairs = [("m", "u")] if direction == "incoming" else [("u", "m")]
    if direction == "both":
        pairs.append(("m", "u"))
    best = {pair: row(*pair) for pair in pairs}
    best["domain", "u"] = row("domain", "u", 10, .1)
    # A self row and a duplicate member never contribute extra evidence.
    best["u", "u"] = row("u", "u", 100)
    result = candidate(pm, best, ["m", "domain", "m", "missing"])
    assert result.diagnostic_score == 10
    expected_pair = (("m", "u") if record == 0 else ("u", "m")) if direction == "both" else pairs[0]
    assert result.evidence_row is best[expected_pair]
    assert result.low_confidence_pass and not result.domain_only
    assert not result.high_confidence_pass


@pytest.mark.parametrize("seed", range(40))
def test_seeded_candidate_and_local_competition_equivalence(seed):
    rng = random.Random(20260915 + seed)
    pm = {pid: bench.protein(pc, pid, rng.randrange(3), i) for i, pid in enumerate(["u", "v"] + [f"m{i}" for i in range(18)])}
    groups = {f"og_{g}": {f"m{i}" for i in range(g * 6, (g + 1) * 6)} for g in range(3)}
    best = {}
    for u in ["u", "v"]:
        for m in pm:
            for q, s in [(u, m), (m, u)]:
                if rng.random() < .55:
                    best[q, s] = row(q, s, rng.choice([0., .4, 1., 1., 2.]), rng.choice([.1, .3, .7, 1.]), evalue=rng.choice([1e-20, 1e-30]))
    thresholds = {pid: pc._LocalThreshold(pid, rng.choice([.5, 1., 3.]), "fixture", 0, 0) for pid in pm if rng.random() < .6}
    mapping = {pid: gid for gid, members in groups.items() for pid in members}
    index = pc._index_core_support_evidence(best, mapping)
    for u in ["u", "v"]:
        for gid, members in groups.items():
            expected = old._build_core_support_candidate(u, gid, list(members), best, thresholds, pm)
            evidence = index.get(u, {}).get(gid, ())
            assert pc._build_core_support_candidate(u, gid, evidence, thresholds, pm) == expected
            # Full-rank ties include directional endpoint IDs; bucket insertion order is irrelevant.
            assert pc._build_core_support_candidate(u, gid, list(reversed(evidence)), thresholds, pm) == expected
    local = {"u": rng.choice([0., .5, 2., 5.]), "v": rng.choice([.5, 2., 5.])}
    assert pc._record_local_component_has_competing_core_support(["u", "v"], True, index, thresholds, pm, local) == old._record_local_component_has_competing_core_support(["u", "v"], groups, best, thresholds, pm, local)


def test_fixed_core_then_expanded_competition_and_local_addition():
    pm = {pid: bench.protein(pc, pid, r, i) for i, (pid, r) in enumerate([
        ("a", 0), ("b", 1), ("u", 0), ("v", 0), ("w", 0), ("x", 0), ("y", 0)])}
    best = {pair: row(*pair, score) for pair, score in [
        (("a", "b"), 1.), (("b", "a"), 1.), (("a", "u"), 1.), (("u", "b"), 1.),
        (("u", "v"), 1.), (("v", "w"), 1.), (("w", "v"), 1.),
        (("x", "y"), 1.), (("y", "x"), 1.)]}
    anchors = [pc._AnchorCoreEvidenceEdge("a", "b", best["a", "b"], 1., 1., "rbh")]
    snapshots = []
    original = pc._index_core_support_evidence
    def observe(evidence, mapping):
        snapshots.append(dict(mapping))
        return original(evidence, mapping)
    with patch.object(pc, "_index_core_support_evidence", observe):
        result = pc._build_anchor_core_orthogroups(best, anchors, pm, include_singletons=False, max_related_edges_per_orthogroup=2)
    assert "u" not in snapshots[0] and snapshots[1]["u"] == "og_1"
    assert result.member_by_protein_id["u"].role == "inparalog"
    # v cannot chain through u during assignment; u competes in the later local phase.
    assert "v" not in result.member_by_protein_id and "w" not in result.member_by_protein_id
    assert result.member_by_protein_id["x"].orthogroup_id == result.member_by_protein_id["y"].orthogroup_id == "og_2"
    assert result.scope_by_orthogroup_id == {"og_1": "cross_record", "og_2": "record_local"}
    assert all("x" not in snapshot for snapshot in snapshots)
    # The early snapshot would miss this competing core evidence.
    core = {"og_1": {"a", "b"}}
    expanded = {"og_1": {"a", "b", "u"}}
    thresholds = pc._derive_anchor_core_thresholds(best, anchors, pm)
    local = {"v": 1., "w": 1.}
    assert not old._record_local_component_has_competing_core_support(["v", "w"], core, best, thresholds, pm, local)
    assert old._record_local_component_has_competing_core_support(["v", "w"], expanded, best, thresholds, pm, local)


def test_sparse_visits_and_giant_groups():
    for name in ["sparse-200", "sparse-800", "giant-200", "giant-800"]:
        pm, tables = bench.synthetic(pc, name, bench.SEED)
        def fn():
            return pc.select_rbh_orthogroup_edges_from_directional_hits(tables, pm, record_count=3)
        result, counts, _ = bench.operation_probe(pc, cc, fn)
        assert counts.get("support.memberVisits", 0) == 0
        if name.startswith("sparse"):
            assert counts.get("_build_core_support_candidate.calls", 0) == 0
        else:
            assert len(result.orthogroups.orthogroups) == 2
            assert counts["support.evidenceVisits"] == 20
            assert counts["_build_core_support_candidate.calls"] == 20


def metadata_fixture():
    pm, tables = bench.synthetic(pc, "support-edges", bench.SEED)
    selection = pc.select_rbh_orthogroup_edges_from_directional_hits(tables, pm, record_count=3)
    # These metadata fixtures deliberately replace edges with arbitrary legacy data.
    return pm, pc.materialize_ortholog_paths(selection.orthogroups)


def test_endpoint_first_match_counts_and_unknown_errors():
    pm, groups = metadata_fixture()
    gid = next(iter(groups.orthogroups))
    edge = replace(groups.ortholog_edges_by_orthogroup_id[gid][0],
                   source_rbh_orthogroup_id="source", target_rbh_orthogroup_id="target")
    reverse = replace(edge, query_protein_id=edge.subject_protein_id, subject_protein_id=edge.query_protein_id, path_id="second")
    groups = replace(groups, ortholog_edges_by_orthogroup_id={gid: (edge, reverse, edge)}, related_edges_by_orthogroup_id={gid: (replace(reverse, path_id="related"),)})
    cache = {}
    for q, s in [(edge.query_protein_id, edge.subject_protein_id), (edge.subject_protein_id, edge.query_protein_id), ("missing", "unknown")]:
        assert pc._edge_metadata_for_protein_pair(groups, gid, q, s, cache) == old._edge_metadata_for_protein_pair(groups, gid, q, s)
    members = groups.orthogroups[gid]
    groups.orthogroups[gid] = [*members, members[0]]
    counts = pc._orthogroup_member_counts(groups.orthogroups[gid])
    for r in range(3):
        assert counts.get(r, 0) == sum(int(m.record_index) == r for m in groups.orthogroups[gid])
    table = bench.frame([bench.hit(edge.query_protein_id, edge.subject_protein_id)] * 3)
    links = pc.convert_pair_protein_hits_to_genomic_links(table, pm, pm, groups)
    assert links.ortholog_path_id.tolist() == [edge.path_id or ""] * 3
    assert links.rbh_orthogroup_id.tolist() == ["source;target"] * 3
    assert len(cache) == 1
    for bad in ["missing", "unknown"]:
        with pytest.raises(ParseError, match=bad):
            pc.convert_pair_protein_hits_to_genomic_links(bench.frame([bench.hit(bad, edge.subject_protein_id)]), pm, pm, groups)
    assert pc.convert_pair_protein_hits_to_genomic_links(bench.frame([]), pm, pm, groups).empty


def test_rbh_many_to_many_duplicate_text_ids_order_and_full_payload():
    _, groups = metadata_fixture()
    ids = list(groups.member_by_protein_id)
    rbh = {"last": (ids[-1], ids[0], ids[0]), 1: (ids[0],), "1": (ids[0], ids[1]), "empty": (), "unknown": ("absent",)}
    groups = replace(groups, rbh_orthogroups=rbh)
    actual = serialize_orthogroups_payload(groups)
    # All fields other than RBH IDs remain equal to a serialization with no RBH map.
    without = serialize_orthogroups_payload(replace(groups, rbh_orthogroups={}))
    for payload, plain, members in zip(actual, without, groups.orthogroups.values()):
        member_ids = {str(m.protein_id or "") for m in members}
        expected = [str(gid or "") for gid, pids in rbh.items() if member_ids.intersection(str(pid or "") for pid in pids)]
        assert payload["rbhOrthogroupIds"] == expected
        assert payload | {"rbhOrthogroupIds": []} == plain
    assert serialize_orthogroups_payload(None) == []
    assert serialize_orthogroups_payload(pc.OrthogroupResult({}, {}, rbh_orthogroups={"bad": None})) == []
    with pytest.raises(TypeError):
        serialize_orthogroups_payload(replace(groups, rbh_orthogroups={"bad": None}))


def test_projection_suppression_preserves_inference_and_all_edges():
    pm, tables = bench.synthetic(pc, "support-edges", bench.SEED)
    full = pc.select_rbh_orthogroup_edges_from_directional_hits(tables, pm, record_count=3)
    reduced = pc.select_rbh_orthogroup_edges_from_directional_hits(tables, pm, record_count=3, comparison_pairs=())
    assert full.orthogroups == reduced.orthogroups
    assert bench.canonical(full.all_edges_by_pair) == bench.canonical(reduced.all_edges_by_pair)
    assert reduced.adjacent_anchor_edges_by_pair == reduced.adjacent_display_edges_by_pair == {}


@pytest.mark.parametrize("cap", [1, 2])
def test_complete_ties_best_second_confidence_and_related_order(cap):
    pm = {pid: bench.protein(pc, pid, r, i) for i, (pid, r) in enumerate([
        ("a", 0), ("b", 1), ("c", 0), ("d", 1),
        ("high", 0), ("low", 0), ("tie0", 2), ("tie1", 2), ("tie2", 2)])}
    best = {pair: row(*pair) for pair in [("a", "b"), ("b", "a"), ("c", "d"), ("d", "c"), ("a", "high")]}
    best["a", "low"] = row("a", "low", .3)
    for pid in ["tie0", "tie1", "tie2"]:
        for member in ["a", "c"]:
            best[pid, member] = row(pid, member, .8)
    anchors = [pc._AnchorCoreEvidenceEdge(q, s, best[q, s], 1., 1., "rbh") for q, s in [("a", "b"), ("c", "d")]]
    result = pc._build_anchor_core_orthogroups(best, anchors, pm, include_singletons=False, max_related_edges_per_orthogroup=cap)
    assert (result.member_by_protein_id["high"].role, result.member_by_protein_id["high"].confidence) == ("inparalog", "high")
    assert (result.member_by_protein_id["low"].role, result.member_by_protein_id["low"].confidence) == ("low_confidence", "low")
    assert result.member_by_protein_id["low"].best_core_support == .15
    assert result.member_by_protein_id["low"].second_best_core_support == 0.
    assert all(pid not in result.member_by_protein_id for pid in ["tie0", "tie1", "tie2"])
    edges = result.related_edges_by_orthogroup_id["og_1"]
    assert len(edges) == cap
    assert [edge.query_protein_id for edge in edges] == [f"tie{i}" for i in range(cap)]
    assert all(edge.edge_kind == "ambiguous_paralog" for edge in edges)
    thresholds = pc._derive_anchor_core_thresholds(best, anchors, pm)
    scores = [old._build_core_support_candidate("tie0", gid, members, best, thresholds, pm) for gid, members in [("og_1", ["a", "b"]), ("og_2", ["c", "d"])]]
    assert [(s.support, s.diagnostic_score) for s in scores] == [(.8, .8), (.8, .8)]


@pytest.mark.parametrize("score", [0., -1., float("nan"), float("inf"), "bad"])
def test_empty_unknown_self_and_nonpositive_candidate_boundaries(score):
    pm = {"u": bench.protein(pc, "u", 0), "m": bench.protein(pc, "m", 1)}
    best = {("u", "m"): row("u", "m", score), ("u", "u"): row("u", "u", 20.),
            ("missing", "u"): row("missing", "u", 20.)}
    assert candidate(pm, best, ["m", "missing"]) is None
    assert candidate(pm, {}, []) is None
    assert not pc._record_local_component_has_competing_core_support(["u"], False, {}, {}, pm, {})
    assert pc._record_local_component_has_competing_core_support(["u"], True, {}, {}, pm, {})


def test_endpoint_index_consumes_each_edge_once_as_needed():
    _, groups = metadata_fixture()
    gid = next(iter(groups.orthogroups))
    template = groups.ortholog_edges_by_orthogroup_id[gid][0]
    edges = tuple(replace(template, query_protein_id=f"q{i}", subject_protein_id=f"s{i}", path_id=f"path-{i}") for i in range(10))
    visits = []
    class ObservedEdges(tuple):
        def __getitem__(self, index):
            edge = super().__getitem__(index)
            visits.append(edge.path_id)
            return edge
    plain = replace(groups, ortholog_edges_by_orthogroup_id={gid: edges}, related_edges_by_orthogroup_id={})
    observed = replace(plain, ortholog_edges_by_orthogroup_id={gid: ObservedEdges(edges)})
    indexes = {}
    for i, consumed in [(3, 4), (1, 4), (8, 9), (3, 9), (99, 10), (0, 10)]:
        q, s = f"q{i}", f"s{i}"
        assert pc._edge_metadata_for_protein_pair(observed, gid, s, q, indexes) == old._edge_metadata_for_protein_pair(plain, gid, s, q)
        assert len(visits) == consumed
    assert visits == [f"path-{i}" for i in range(10)]
