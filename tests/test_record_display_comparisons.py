"""Issue 464 Session 5: source identities through annotation/comparison SVGs."""
from __future__ import annotations

import copy
import math
import re
from types import SimpleNamespace
from dataclasses import replace
from xml.etree import ElementTree as ET

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.annotations import (
    AnnotationOptions, AnnotationSet, CoordinateSpan, FeatureSelector, FeatureSpan,
    RegionAnnotation, RegionAnnotationStyle,
)
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import (
    CircularDiagramRequest, InMemoryRecordSource, LinearDiagramRequest,
    RecordDisplayOptions, RecordInput, RecordPresentation,
)
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.io.record_select import parse_record_selector
from gbdraw.linear_comparison import LinearComparison
from gbdraw.render.interactive_context import build_interactive_svg_context
from gbdraw.render.interactive_svg import enrich_svg
from gbdraw.exceptions import GbdrawError
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.linear_comparison import project_match_endpoints


def record(length=100, name="reference"):
    item = SeqRecord(Seq(("AAAACCGTGG" * length)[:length]), id=name,
                     annotations={"topology": "circular", "molecule_type": "DNA"})
    item.features = [
        SeqFeature(SimpleLocation(20, 60, strand=1), type="CDS",
                   qualifiers={"locus_tag": ["anchor"], "product": ["anchor enzyme"]}),
        SeqFeature(SimpleLocation(20, 60, strand=1), type="repeat_region",
                   qualifiers={"rpt_family": ["repeat"]}),
    ]
    return item


def request(mode, records, starts, reverses=None, **options):
    cls, opt = ((CircularDiagramRequest, CircularDiagramOptions) if mode == "circular"
                else (LinearDiagramRequest, LinearDiagramOptions))
    config = apply_config_overrides(None, {
        "canvas.show_gc": False, "canvas.show_skew": False,
        "labels.circular.scope": "none", "labels.linear.scope": "none",
    })
    return cls(records=tuple(RecordInput(
        InMemoryRecordSource(item), display=RecordDisplayOptions(start_coordinate=start),
        presentation=RecordPresentation(reverse_complement=reverse),
    ) for item, start, reverse in zip(records, starts, reverses or [False] * len(records), strict=True)),
        options=opt(config=config, **options))


def svg(req):
    plan = plan_request(req)
    result = plan.build()
    drawing = result.drawing if hasattr(result, "drawing") else result
    return plan, ET.fromstring(drawing.tostring())


def xy(path):
    return [(float(x), float(y)) for x, y in re.findall(
        r"[ML]\s*([-+\d.eE]+)[ ,]+([-+\d.eE]+)", path)]


def assert_source(item, before):
    assert item.seq == before.seq
    assert item.features == before.features
    assert item.annotations == before.annotations


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("target_kind", ["coordinate", "feature"])
def test_annotation_and_underlay_project_source_coverage(mode, reverse, target_kind):
    item = record()
    before = copy.deepcopy(item)
    target = (CoordinateSpan(None, 21, 60) if target_kind == "coordinate" else
              FeatureSpan(None, (FeatureSelector("anchor", key="locus_tag"),), envelope="segments"))
    annotations = AnnotationOptions(sets=(AnnotationSet("regions", (
        RegionAnnotation("region", target, mark="highlight", style=RegionAnnotationStyle(fill="#123456")),
    )),))
    plan, root = svg(request(mode, [item], [41], [reverse], annotations=annotations))
    annotation = next(n for n in root.iter() if n.get("data-gbdraw-annotation-id") == "region")
    shapes = [n for n in annotation if n.tag.endswith(("}path", "}rect")) and n.get("fill") != "none"]
    underlays = [n for n in root.iter() if n.get("data-gbdraw-auto-feature-underlay") == "true"]
    expected = [(0, 20), (80, 100)] if not reverse else [(0, 21), (81, 100)]
    assert len(shapes) == len(underlays) == 2
    for actual in (shapes, underlays):
        if mode == "linear":
            axis = next(n for n in root.iter() if n.get("data-gbdraw-record-id") == "reference"
                        and any(c.tag.endswith("}line") and c.get("y1") == c.get("y2") for c in n))
            bar = next(n for n in axis if n.tag.endswith("}line") and n.get("y1") == n.get("y2"))
            scale = float(bar.get("x2")) / 100
            observed = sorted((float(n.get("x")) / scale,
                               (float(n.get("x")) + float(n.get("width"))) / scale) for n in actual)
            for bounds, oracle in zip(observed, expected, strict=True):
                assert bounds == pytest.approx(oracle)
        else:
            angles = []
            for n in actual:
                x, y = xy(n.get("d"))[0]
                angles.append(((math.degrees(math.atan2(y, x)) + 90) % 360) / 3.6)
            assert sorted(angles) == pytest.approx([a for a, _ in expected])
    assert len({n.get("id") for n in underlays}) == 2
    assert len({n.get("data-gbdraw-feature-id") for n in underlays}) == 1
    assert plan.transforms[0].start_coordinate == 41
    assert_source(item, before)


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("reverse", [False, True])
def test_origin_wrap_annotation_uses_resolved_topology_in_both_modes(mode, reverse):
    item = record()
    annotation = RegionAnnotation("wrap", CoordinateSpan(None, 91, 10, wraps_origin=True), mark="band")
    options = AnnotationOptions(sets=(AnnotationSet("regions", (annotation,)),))
    _, root = svg(request(mode, [item], [41], [reverse], annotations=options))
    group = next(n for n in root.iter() if n.get("data-gbdraw-annotation-id") == "wrap")
    assert len(group) >= 1
    assert annotation.target.start == 91 and annotation.target.end == 10


def hit_frame(qstart=21, qend=60, sstart=31, send=90):
    return pd.DataFrame([("reference", "subject", 90., 60, 2, 1,
                          qstart, qend, sstart, send, 1e-20, 100.)], columns=COMPARISON_COLUMNS)


@pytest.mark.parametrize("starts,expected_fragments", [((None, None), 1), ((41, None), 2),
                                                     ((None, 51), 2), ((41, 51), 3)])
def test_linear_hsp_common_t_path_count_and_source_endpoints(starts, expected_fragments):
    frame = hit_frame()
    before = frame.copy(deep=True)
    _, root = svg(request("linear", [record(), record(120, "subject")], starts,
                          linear_comparisons=[LinearComparison(0, 1, frame)]))
    paths = [n for n in root.iter() if n.get("data-gbdraw-match-id")]
    assert len(paths) == expected_fragments
    assert len({n.get("data-gbdraw-match-id") for n in paths}) == 1
    assert len({n.get("id") for n in paths}) == len(paths)
    for n in paths:
        assert [n.get(k) for k in ("data-qstart", "data-qend", "data-sstart", "data-send")] == ["21", "60", "31", "90"]
        assert float(n.get("data-identity")) == 90.
    pd.testing.assert_frame_equal(frame, before)


@pytest.mark.parametrize("reverse", [False, True])
def test_circular_comparison_projects_reference_span_and_keeps_raw_endpoints(reverse):
    frame = hit_frame()
    before = frame.copy(deep=True)
    _, root = svg(request("circular", [record()], [41], [reverse],
                          conservation_dataframes=[frame], conservation_reference="query"))
    paths = [n for n in root.iter() if n.get("data-gbdraw-match-id")]
    assert len(paths) == 2
    assert len({n.get("data-gbdraw-match-id") for n in paths}) == 1
    assert len({n.get("id") for n in paths}) == 2
    for n in paths:
        assert [n.get(k) for k in ("data-qstart", "data-qend", "data-sstart", "data-send")] == ["21", "60", "31", "90"]
    pd.testing.assert_frame_equal(frame, before)


@pytest.mark.parametrize("mode", ["linear", "circular"])
def test_interactive_catalog_has_one_source_match_for_all_dom_fragments(mode):
    items = [record(), record(120, "subject")] if mode == "linear" else [record()]
    options = ({"linear_comparisons": [LinearComparison(0, 1, hit_frame())]} if mode == "linear" else
               {"conservation_dataframes": [hit_frame()], "conservation_reference": "query"})
    plan, root = svg(request(mode, items, [41, 51] if mode == "linear" else [41], **options))
    source = ET.tostring(root, encoding="unicode")
    context = build_interactive_svg_context(plan.records, mode=mode, linear_rendered_feature_ids=mode == "linear", record_transforms=plan.transforms)
    enriched = ET.fromstring(enrich_svg(source, context))
    paths = [n for n in enriched.iter() if n.get("data-gbdraw-interactive-match") == "true"]
    assert len(paths) == (3 if mode == "linear" else 2)
    import json
    catalog = json.loads(next(n.text for n in enriched.iter() if n.get("id") == "gbdraw-interactive-feature-metadata"))
    assert len(catalog["items"][0]["comparisonMatches"]) == 1
    assert len({n.get("id") for n in paths}) == len(paths)


def endpoint_oracle(endpoints, starts, lengths, reverses):
    """Independent D8/D11 oracle in existing record-local HSP coordinates."""
    if starts == (None, None):
        return [endpoints]
    spans = [(a - 1, b) if a <= b else (a, b - 1) for a, b in (endpoints[:2], endpoints[2:])]
    cuts = [(length - start if reverse else start - 1) if start is not None else 0
            for start, length, reverse in zip(starts, lengths, reverses, strict=True)]
    ts = {0., 1.}
    for (a, b), cut, start in zip(spans, cuts, starts, strict=True):
        if start is not None and min(a, b) < cut < max(a, b):
            ts.add((cut - a) / (b - a))
    result = []
    for left, right in zip(sorted(ts), sorted(ts)[1:]):
        values = []
        for (a, b), cut, length, start in zip(spans, cuts, lengths, starts, strict=True):
            offset = (length if (a + (b - a) * (left + right) / 2) < cut else 0) if start is not None else 0
            values.extend(a + (b - a) * t - cut + offset for t in (left, right))
        result.append(tuple(values))
    return result


@pytest.mark.parametrize("starts,reverses", [
    ((None, None), (False, False)), ((41, None), (False, False)),
    ((None, 51), (False, False)), ((41, 51), (False, False)),
    ((41, 51), (True, False)), ((41, 51), (False, True)),
    ((41, 51), (True, True)), ((1, 120), (False, True)),
    ((100, 1), (True, False)), ((21, 31), (False, False)),
    ((60, 90), (False, False)),
])
@pytest.mark.parametrize("directed", [(21, 60, 31, 90), (60, 21, 31, 90), (21, 60, 90, 31), (60, 21, 90, 31)])
def test_hsp_common_t_independent_oracle(starts, directed, reverses):
    transforms = tuple(RecordDisplayTransform(length, length if reverse else 1, -1 if reverse else 1, start, True)
                       for length, start, reverse in zip((100, 120), starts, reverses, strict=True))
    projected = project_match_endpoints(directed, transforms)
    expected = endpoint_oracle(directed, starts, (100, 120), reverses)
    assert len(projected) == len(expected)
    for actual, oracle in zip(projected, expected, strict=True):
        assert actual == pytest.approx(oracle)
    for side in (0, 2):
        assert sum(abs(part[side + 1] - part[side]) for part in projected) == pytest.approx(
            abs(directed[side + 1] - directed[side]) + (0 if starts == (None, None) else 1))


@pytest.mark.parametrize("ends", [(1, 1, 120, 120), (100, 100, 1, 1), (40, 41, 50, 51)])
def test_hsp_short_valid_spans_keep_coverage(ends):
    transforms = (RecordDisplayTransform(100, 1, 1, 41, True), RecordDisplayTransform(120, 1, 1, 51, True))
    actual = project_match_endpoints(ends, transforms)
    for part, expected in zip(actual, endpoint_oracle(ends, (41, 51), (100, 120), (False, False)), strict=True):
        assert part == pytest.approx(expected)


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("reverse", [False, True])
def test_annotation_artificial_seam_has_no_bracket_caps(mode, reverse):
    annotation = RegionAnnotation("bracket", CoordinateSpan(None, 21, 60), mark="bracket")
    _, root = svg(request(mode, [record()], [41], [reverse],
                          annotations=AnnotationOptions(sets=(AnnotationSet("a", (annotation,)),))))
    group = next(n for n in root.iter() if n.get("data-gbdraw-annotation-id") == "bracket")
    # Two span shapes and two real endpoint caps, independent of mode/orientation.
    assert len(group) == 4


@pytest.mark.parametrize("mutation", ["duplicate_index", "duplicate_dom", "missing_marker", "conflicting_endpoint", "missing_endpoint", "conflicting_alias"])
def test_interactive_rejects_real_fragment_identity_conflicts(mutation):
    plan, root = svg(request("linear", [record(), record(120, "subject")], [41, 51],
                             linear_comparisons=[LinearComparison(0, 1, hit_frame())]))
    paths = [n for n in root.iter() if n.get("data-gbdraw-match-id")]
    if mutation == "duplicate_index":
        paths[1].set("data-gbdraw-match-fragment", paths[0].get("data-gbdraw-match-fragment"))
    elif mutation == "duplicate_dom":
        paths[1].set("id", paths[0].get("id"))
    elif mutation == "missing_marker":
        del paths[1].attrib["data-gbdraw-match-fragment"]
    elif mutation == "conflicting_endpoint":
        paths[1].set("data-qstart", "22")
    elif mutation == "missing_endpoint":
        for path in paths:
            del path.attrib["data-qstart"]
    else:
        paths[1].set("data-gbdraw-pairwise-match-id", "other")
    context = build_interactive_svg_context(plan.records, mode="linear", linear_rendered_feature_ids=True)
    with pytest.raises(GbdrawError, match="match ID|fragment metadata"):
        enrich_svg(ET.tostring(root, encoding="unicode"), context)


def test_group_alignment_projects_member_centers_before_translation():
    from gbdraw.diagrams.linear.orthogroup_alignment import calculate_orthogroup_alignment_offsets
    frame = hit_frame().assign(orthogroup_id="og_1", query_record_index=0, subject_record_index=1,
                              query_feature_index=0, subject_feature_index=0,
                              query_feature_svg_id="anchor", subject_feature_svg_id="target")
    before = frame.copy(deep=True)
    config = SimpleNamespace(normalize_length=False, align_center=False, alignment_width=120., longest_genome=120.)
    offsets = calculate_orthogroup_alignment_offsets([record(), record(120, "subject")], [frame], config, "anchor",
        record_transforms=(RecordDisplayTransform(100, 1, 1, 41, True), RecordDisplayTransform(120, 1, 1, 51, True)))
    # Historical inclusive centers: 40.5 and 60.5; projected: 0.5 and 10.5.
    assert offsets == {0: 0., 1: -10.}
    pd.testing.assert_frame_equal(frame, before)


@pytest.mark.parametrize("reverse", [False, True])
def test_final_hsp_vertices_equal_independent_common_t_oracle(reverse):
    _, root = svg(request("linear", [record(), record(120, "subject")], [41, 51], [reverse, False],
                          linear_comparisons=[LinearComparison(0, 1, hit_frame(60, 21, 31, 90))]))
    paths = [n for n in root.iter() if n.get("data-gbdraw-match-id")]
    axes = [n for n in root.iter() if n.get("data-gbdraw-record-id") and any(c.tag.endswith("}line") for c in n)]
    scales = [float(next(c for c in axis if c.tag.endswith("}line")).get("x2")) / length
              for axis, length in zip(axes, (100, 120), strict=True)]
    for path, expected in zip(paths, endpoint_oracle((60, 21, 31, 90), (41, 51), (100, 120), (reverse, False)), strict=True):
        points = xy(path.get("d"))
        assert [p[0] for p in points] == pytest.approx([expected[0] * scales[0], expected[1] * scales[0], expected[3] * scales[1], expected[2] * scales[1]])
        assert points[0][1] == points[1][1] < points[2][1] == points[3][1]


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_comparison_annotation_prepared_cache_a_b_a(mode, monkeypatch, tmp_path):
    from Bio import SeqIO
    from gbdraw.api.prepared import PreparedBiologicalInputCache, PreparedResourceIdentity
    from gbdraw.api.requests import GenBankInputSource
    import gbdraw.api.request_render as owner
    path = tmp_path / "source.gb"
    SeqIO.write(record(), path, "genbank")
    source_bytes = path.read_bytes()
    cache = PreparedBiologicalInputCache()
    identity = PreparedResourceIdentity("source", "source-resource", len(source_bytes))
    load = owner.load_gbks
    calls = []
    def capture(*args, **kwargs):
        calls.append(1)
        return load(*args, **kwargs)
    monkeypatch.setattr(owner, "load_gbks", capture)
    frame = hit_frame()
    before = frame.copy(deep=True)
    annotations = AnnotationOptions(sets=(AnnotationSet("a", (RegionAnnotation("span", CoordinateSpan(parse_record_selector("#1"), 21, 60)),)),))
    options = ({"linear_comparisons": [LinearComparison(0, 1, frame)]} if mode == "linear" else
               {"conservation_dataframes": [frame], "conservation_reference": "query"})
    snapshots, records, contexts = [], [], []
    for start in (41, 11, 41):
        req = request(mode, [record(), record(120, "subject")] if mode == "linear" else [record()],
                      [start, 51] if mode == "linear" else [start], annotations=annotations, **options)
        req = replace(req, records=(replace(req.records[0], source=GenBankInputSource(path)), *req.records[1:]))
        with cache.transaction(resource_paths={path: identity}, diagnostics=None):
            prepared = owner.build_request_diagram(req)
            context = owner.build_prepared_interactive_context(prepared)
            contexts.append(context)
            records.append(prepared.records[0])
            snapshots.append(enrich_svg(prepared.drawing.tostring(), context))
    assert calls == [1]
    assert all(item.seq == records[0].seq and item.features == records[0].features for item in records)
    assert snapshots[0] != snapshots[1] and snapshots[0] == snapshots[2]
    assert contexts[0] == contexts[1] == contexts[2]
    assert path.read_bytes() == source_bytes
    pd.testing.assert_frame_equal(frame, before)


@pytest.mark.parametrize("starts", [(None, None), (201, 301)])
def test_rotated_mixed_groups_keep_semantic_membership_and_presentation_scope(starts):
    import json
    from tests.test_collinearity import _record_local_mixed_fixture
    from gbdraw.analysis.collinearity import build_orthogroup_collinearity_blocks_from_hits, convert_collinearity_blocks_to_comparisons
    items, extraction, hits = _record_local_mixed_fixture()
    for item in items:
        item.annotations["topology"] = "circular"
    source = copy.deepcopy(items)
    before_hits = {key: frame.copy(deep=True) for key, frame in hits.items()}
    result = build_orthogroup_collinearity_blocks_from_hits(hits, extraction, records=items)
    groups = result.orthogroups
    assert {gid: {m.protein_id for m in members} for gid, members in groups.orthogroups.items()} == {"og_1": {"a0", "a1", "b0"}, "og_2": {"a3", "a4"}}
    assert groups.scope_by_orthogroup_id == {"og_1": "cross_record", "og_2": "record_local"}
    assert groups.member_by_protein_id["a1"].role == "inparalog"
    assert {m.role for m in groups.orthogroups["og_2"]} == {"local_paralog"}
    assert {"a2", "b1"}.isdisjoint(groups.member_by_protein_id)
    frame = convert_collinearity_blocks_to_comparisons(result, records=items)[0]
    before_frame = frame.copy(deep=True)
    plan, root = svg(request("linear", items, starts, linear_comparisons=[LinearComparison(0, 1, frame)]))
    paths = [n for n in root.iter() if n.get("data-gbdraw-match-id")]
    assert len(paths) == (1 if starts[0] is None else 3)
    assert len({n.get("data-gbdraw-match-id") for n in paths}) == 1
    assert {n.get("data-orthogroup-id") for n in paths} == {"og_1"}
    context = build_interactive_svg_context(plan.records, mode="linear", linear_rendered_feature_ids=True, orthogroups=groups)
    enriched = ET.fromstring(enrich_svg(ET.tostring(root, encoding="unicode"), context))
    catalog = json.loads(next(n.text for n in enriched.iter() if n.get("id") == "gbdraw-interactive-feature-metadata"))["items"][0]
    catalog_groups = {g["id"]: g for g in catalog["orthogroups"]}
    assert catalog_groups["og_1"]["scope"] == "cross_record"
    assert catalog_groups["og_2"]["scope"] == "record_local"
    assert {m["orthogroup_ids"][0] for m in catalog["comparisonMatches"]} == {"og_1"}
    for item, original in zip(items, source, strict=True):
        assert_source(item, original)
    for key, frame in hits.items():
        pd.testing.assert_frame_equal(frame, before_hits[key])
    pd.testing.assert_frame_equal(before_frame, plan.request.options.linear_comparisons[0].matches)


@pytest.mark.parametrize("reference", ["query", "subject"])
def test_circular_multiple_rings_keep_coverage_and_unique_ids(reference):
    frame = hit_frame()
    _, root = svg(request("circular", [record(120 if reference == "subject" else 100, "subject" if reference == "subject" else "reference")], [41],
                          conservation_dataframes=[frame, frame.copy()], conservation_reference=reference))
    paths = [n for n in root.iter() if n.get("data-gbdraw-match-id")]
    assert len(paths) == 4
    assert len({n.get("id") for n in paths}) == 4
    assert len({n.get("data-gbdraw-match-id") for n in paths}) == 2
    length = 120 if reference == "subject" else 100
    expected_parts = [(110, 120), (0, 50)] if reference == "subject" else [(80, 100), (0, 20)]
    for pair in (paths[:2], paths[2:]):
        for path, (start, end) in zip(pair, expected_parts, strict=True):
            first = xy(path.get("d"))[0]
            angle = (math.degrees(math.atan2(first[1], first[0])) + 90) % 360
            assert angle == pytest.approx(start * 360 / length)
            assert float(path.get("data-identity")) == 90.


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_annotation_band_outlines_exclude_artificial_edges(mode):
    annotation = RegionAnnotation("band", CoordinateSpan(None, 21, 60), mark="band")
    _, root = svg(request(mode, [record()], [41], annotations=AnnotationOptions(sets=(AnnotationSet("a", (annotation,)),))))
    group = next(n for n in root.iter() if n.get("data-gbdraw-annotation-id") == "band")
    assert len(group) == 4
    outlines = [n for n in group if n.get("stroke-linecap") == "butt"]
    assert len(outlines) == 2
    for path in outlines:
        assert 'Z' not in path.get('d').upper()
        assert path.get('d').count('L') == (3 if mode == "linear" else 1)


@pytest.mark.parametrize("mode", ["linear", "circular"])
def test_duplicate_record_instances_use_their_own_comparison_transform(mode):
    from gbdraw.api.options import CircularMultiRecordOptions
    item = record()
    frame = hit_frame(21, 60, 21, 60).assign(sseqid="reference")
    options = ({"linear_comparisons": [LinearComparison(0, 1, frame)]} if mode == "linear" else
               {"conservation_dataframes": [frame], "conservation_reference": "query"})
    req = request(mode, [item, item], [41, 11], **options)
    if mode == "circular":
        req = replace(req, grouping="grid", layout=CircularMultiRecordOptions())
    plan, root = svg(req)
    assert [t.start_coordinate for t in plan.transforms] == [41, 11]
    paths = [n for n in root.iter() if n.get("data-gbdraw-match-id")]
    assert len(paths) == (2 if mode == "linear" else 3)
    assert len({n.get("id") for n in paths}) == len(paths)
    context = build_interactive_svg_context(plan.records, mode=mode, linear_rendered_feature_ids=mode == "linear", record_transforms=plan.transforms)
    enrich_svg(ET.tostring(root, encoding="unicode"), context)
