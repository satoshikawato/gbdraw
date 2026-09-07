"""Issue 464 Session 3: internal feature geometry, without later consumers."""
from __future__ import annotations

import copy
import re
from xml.etree import ElementTree as ET

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SimpleLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.prepared import PreparedBiologicalInputCache, PreparedResourceIdentity
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import (
    CircularDiagramRequest, LinearDiagramRequest, RecordInput, RecordDisplayOptions,
    RecordPresentation, InMemoryRecordSource, GenBankInputSource,
)
from gbdraw.features.coordinates import project_feature_parts


def _record(parts=None):
    record = SeqRecord(Seq("AAACCGTTAA"), id="duplicate", annotations={
        "topology": "circular", "molecule_type": "DNA",
    })
    parts = parts or [SimpleLocation(2, 7, strand=1)]
    record.features = [SeqFeature(
        parts[0] if len(parts) == 1 else CompoundLocation(parts), type="CDS",
        qualifiers={"product": ["test enzyme"], "locus_tag": ["cds1"]},
    )]
    return record


def _request(mode, source, start, reverse=False):
    options_type = CircularDiagramOptions if mode == "circular" else LinearDiagramOptions
    request_type = CircularDiagramRequest if mode == "circular" else LinearDiagramRequest
    cfg = apply_config_overrides(None, {
        "labels.circular.scope": "none", "labels.linear.scope": "none",
        "canvas.show_gc": False, "canvas.show_skew": False,
        "objects.features.arrow_geometry.head_length_ratio": 1.0,
    })
    return request_type(records=(RecordInput(
        source, display=RecordDisplayOptions(start_coordinate=start),
        presentation=RecordPresentation(reverse_complement=reverse),
    ),), options=options_type(config=cfg))


def _paths(plan):
    result = plan.build()
    svg = result.drawing if hasattr(result, "drawing") else result
    return [node for node in ET.fromstring(svg.tostring()).iter()
            if node.get("data-gbdraw-feature-part") == "block"]


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("reverse", [False, True])
def test_plan_factory_drawer_rotation_and_none_nonmutation(mode, reverse, monkeypatch):
    import gbdraw.diagrams.circular.assemble as circular
    import gbdraw.diagrams.linear.precalc as linear
    owner = circular if mode == "circular" else linear
    built = []
    real = owner.create_feature_layers
    def capture(*args, **kwargs):
        result = real(*args, **kwargs)
        built.append(result.foreground_features["gene_000000001"])
        return result
    monkeypatch.setattr(owner, "create_feature_layers", capture)
    record = _record()
    before = copy.deepcopy(record)
    plans = [plan_request(_request(mode, InMemoryRecordSource(record), start, reverse))
             for start in (None, 5, 1, None)]
    paths = [_paths(plan) for plan in plans]
    assert [p.attrib for p in paths[0]] == [p.attrib for p in paths[3]]
    assert [p.get("d") for p in paths[0]] != [p.get("d") for p in paths[1]]
    if reverse:
        assert [p.get("d") for p in paths[0]] != [p.get("d") for p in paths[2]]
    assert built[0].display_parts is None
    blocks = [p.fragment for p in built[1].display_parts if p.kind == "block"]
    expected = [(0, 3), (8, 10)] if reverse else [(8, 10), (0, 3)]
    assert [(p.display_start, p.display_end) for p in blocks] == expected
    assert sum(p.biological_end for p in blocks) == 1
    assert sum(p.source_end - p.source_start for p in blocks) == 5
    assert all(f.coordinates == built[0].coordinates for f in built)
    assert all(f.location == built[0].location for f in built)
    assert all(f.source_feature_index == 0 for f in built)
    stable = [{p.get("data-gbdraw-feature-id") for p in group} for group in paths]
    assert stable[0] == stable[1] == stable[2] == stable[3]
    for group in paths:
        assert len({p.get("id") for p in group}) == len(group)
    assert record.seq == before.seq
    assert record.annotations == before.annotations
    assert record.features == before.features


@pytest.mark.parametrize("strand", [1, -1, None, 0])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("bounds", [(0, 1), (9, 10), (2, 7), (0, 10), (5, 5)])
def test_feature_fragment_coverage_terminals_and_unknown_strand(strand, reverse, bounds):
    record = _record([SimpleLocation(*bounds, strand=strand)])
    plan = plan_request(_request("linear", InMemoryRecordSource(record), 5, reverse))
    fragments = project_feature_parts(plan.records[0].features[0].location.parts, plan.transforms[0])
    blocks = [p for p in fragments if p.kind == "block"]
    assert sum(p.fragment.source_end - p.fragment.source_start for p in blocks) == bounds[1] - bounds[0]
    assert sum(p.fragment.biological_end for p in blocks) == int(bounds[0] != bounds[1] and strand in (1, -1))
    assert all(0 <= p.fragment.display_start < p.fragment.display_end <= 10 for p in blocks)
    if strand not in (1, -1):
        assert all(p.strand == "undefined" for p in blocks)
    # Exercise both final path generators for the same true-boundary geometry.
    from gbdraw.features.factory import create_gene_object
    from gbdraw.render.drawers.circular.features import FeaturePathGenerator as Circular
    from gbdraw.render.drawers.linear.features import FeaturePathGenerator as Linear
    feature = create_gene_object("gene", plan.records[0].features[0], {}, {"CDS": "gray"}, 10, {}, None,
                                 glyph_kind="arrow", compute_label_text=False)
    feature.display_parts = fragments
    circular = Circular(100, 10, .1, .5, 0, "middle", False, head_length_ratio=1.0).generate_circular_gene_path(feature)
    linear = Linear(10, 100, 10, 1, feature.strand, False, 1).generate_linear_gene_path(feature)
    assert len(circular) == len(linear) == len(blocks)
    for part, cp, lp in zip(blocks, circular, linear, strict=True):
        fragment = part.fragment
        xs = [float(v) for v in re.findall(r"[ML]\s*([-+\d.eE]+),", lp[1])]
        assert min(xs) == pytest.approx(fragment.display_start * 10)
        assert max(xs) == pytest.approx(fragment.display_end * 10)
        points = [(float(x), float(y)) for x, y in re.findall(
            r"[ML]\s*([-+\d.eE]+),([-+\d.eE]+)", lp[1])]
        ys = {y for _, y in points}
        assert (len(ys) == 3) == (fragment.biological_end and strand in (1, -1))
        if fragment.biological_end and strand in (1, -1):
            middle = (min(ys) + max(ys)) / 2
            head_x = fragment.display_end if part.strand == "positive" else fragment.display_start
            assert (head_x * 10, middle) in points
        if fragment.artificial_start or fragment.artificial_end:
            assert len(cp) == len(lp) == 3
            assert not cp[2].rstrip().endswith("z")
            assert not lp[2].rstrip().endswith("z")
            previous = None
            for command, x, y in re.findall(r"([ML])\s*([-+\d.eE]+),([-+\d.eE]+)", lp[2]):
                point = (float(x), float(y))
                if previous is not None and command == "L" and previous[0] == point[0]:
                    assert not (fragment.artificial_start and point[0] == 10 * fragment.display_start)
                    assert not (fragment.artificial_end and point[0] == 10 * fragment.display_end)
                previous = point


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_prepared_cache_start_a_b_a_reuses_source_and_record_membership(tmp_path, monkeypatch, mode):
    import gbdraw.api.request_render as owner
    record = _record()
    path = tmp_path / "source.gb"
    SeqIO.write(record, path, "genbank")
    cache = PreparedBiologicalInputCache()
    identity = PreparedResourceIdentity("source", "render-resource-1", path.stat().st_size)
    calls = []
    real = owner.load_gbks
    def load(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)
    monkeypatch.setattr(owner, "load_gbks", load)
    plans, paths, snapshots = [], [], []
    for start in (3, 6, 3):
        with cache.transaction(resource_paths={path: identity}, diagnostics=None):
            plan = plan_request(_request(mode, GenBankInputSource(path), start))
            plans.append(plan)
            paths.append([p.attrib for p in _paths(plan)])
            snapshots.append(copy.deepcopy(plan.records[0]))
    assert len(calls) == 1
    assert plans[0].records[0] is plans[1].records[0] is plans[2].records[0]
    assert [p.transforms[0].start_coordinate for p in plans] == [3, 6, 3]
    assert paths[0] != paths[1]
    assert paths[0] == paths[2]
    assert all(s.seq == snapshots[0].seq and s.features == snapshots[0].features and
               s.annotations == snapshots[0].annotations for s in snapshots)


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("bounds", [[(1, 3), (7, 9)], [(8, 10), (0, 2)], [(8, 9), (1, 2)]])
def test_multipart_origin_and_connector_cut_keep_part_order(reverse, strand, bounds):
    exons = [SimpleLocation(start, end, strand=strand) for start, end in bounds]
    if strand == -1:
        exons.reverse()
    source = _record(exons)
    plan = plan_request(_request("linear", InMemoryRecordSource(source), 5, reverse))
    parts = plan.records[0].features[0].location.parts
    geometry = project_feature_parts(parts, plan.transforms[0])
    blocks = [p.fragment for p in geometry if p.kind == "block"]
    assert [p.part_index for p in blocks] == [0, 1]
    assert sum(p.biological_start for p in blocks) == 1
    assert sum(p.biological_end for p in blocks) == 1
    assert sum(p.source_end - p.source_start for p in blocks) == sum(end - start for start, end in bounds)
    lines = [p.fragment for p in geometry if p.kind == "line"]
    assert all(not p.biological_start and not p.biological_end for p in lines)
    expected_gap = 4 if bounds == [(1, 3), (7, 9)] else (0 if bounds == [(8, 10), (0, 2)] else 2)
    assert sum(p.display_end - p.display_start for p in lines) == expected_gap
    assert all(p.display_end - p.display_start < 10 for p in lines)
    if bounds == [(1, 3), (7, 9)]:
        assert len(lines) == 2  # The intron crosses the display cut; no whole-width line.
    trans_spliced = project_feature_parts(parts, plan.transforms[0], is_trans_spliced=True)
    assert all(p.kind == "block" for p in trans_spliced)


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_multi_record_builder_receives_each_instance_transform(mode, monkeypatch):
    from dataclasses import replace
    from gbdraw.api.options import CircularMultiRecordOptions, LinearMultiRecordOptions
    import gbdraw.diagrams.circular.assemble as circular
    import gbdraw.diagrams.linear.precalc as linear
    owner = circular if mode == "circular" else linear
    real = owner.create_feature_layers
    observed = []
    def capture(*args, **kwargs):
        observed.append(kwargs["record_transform"])
        return real(*args, **kwargs)
    monkeypatch.setattr(owner, "create_feature_layers", capture)
    source = InMemoryRecordSource(_record())
    first = _request(mode, source, 1)
    second = replace(first.records[0], record_key="second", display=RecordDisplayOptions(None, 10))
    layout = (CircularMultiRecordOptions() if mode == "circular" else
              LinearMultiRecordOptions(multi_record_positions=("#1@1", "#2@1")))
    kwargs = {"grouping": "grid"} if mode == "circular" else {}
    plan = plan_request(replace(first, records=(*first.records, second), layout=layout, **kwargs))
    plan.build()
    assert observed == list(plan.transforms)
    assert [t.start_coordinate for t in observed] == [1, 10]
    assert len(plan.provenance) == len(plan.displays) == len(plan.transforms) == 2
    assert plan.provenance[0].record_key != plan.provenance[1].record_key


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("start", [1, 10])
@pytest.mark.parametrize("bounds", [(0, 0), (0, 1), (9, 10), (0, 10)])
def test_full_builder_empty_and_boundary_features(mode, start, bounds):
    record = _record([SimpleLocation(*bounds, strand=1)])
    paths = _paths(plan_request(_request(mode, InMemoryRecordSource(record), start)))
    fills = [p for p in paths if not p.get("id", "").endswith("__outline")]
    if bounds == (0, 0):
        assert fills == []
    else:
        assert len(fills) == (2 if bounds == (0, 10) and start == 10 else 1)
        assert all(p.get("d") and "nan" not in p.get("d").lower() for p in paths)


def test_already_resolved_collection_rebinds_display_without_rc_or_crop():
    from dataclasses import replace
    from gbdraw.api.request_render import _coerce_resolved_collection
    from gbdraw.api.record_planning import resolve_record_inputs
    request = _request("linear", InMemoryRecordSource(_record()), 3, reverse=True)
    collection = resolve_record_inputs(request.records, gff_candidate_features=None, gff_keep_all_features=False)
    changed = replace(request, records=(replace(request.records[0], display=RecordDisplayOptions(None, 7)),))
    rebound = _coerce_resolved_collection(changed, collection)
    assert rebound.records is collection.records
    assert rebound.transforms[0].start_coordinate == 7
    assert rebound.transforms[0].source_base == 10
    assert rebound.transforms[0].source_step == -1
    assert collection.transforms[0].start_coordinate == 3


def test_empty_and_mixed_strand_parts_do_not_gain_connectors_or_new_terminals():
    from gbdraw.layout.record_coordinates import RecordDisplayTransform
    exons = [SimpleLocation(0, 0, strand=1), SimpleLocation(1, 3, strand=1),
             SimpleLocation(5, 5, strand=1), SimpleLocation(7, 9, strand=-1)]
    geometry = project_feature_parts(exons, RecordDisplayTransform(10, 1, 1, 5, True))
    assert [p.fragment.part_index for p in geometry] == [1, 3]
    assert all(p.kind == "block" for p in geometry)
    assert [p.strand for p in geometry] == ["positive", "negative"]
    assert [p.fragment.biological_end for p in geometry] == [False, True]
