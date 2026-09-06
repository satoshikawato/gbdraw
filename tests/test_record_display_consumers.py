"""Issue 464 Session 4: typed plans through the final consumer SVG."""
from __future__ import annotations

import copy
import math
import re
from dataclasses import replace
from xml.etree import ElementTree as ET

import pytest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import (
    CircularDiagramRequest, InMemoryRecordSource, LinearDiagramRequest,
    RecordDisplayOptions, RecordInput, RecordPresentation,
)


def _record(length=10000):
    record = SeqRecord(Seq("AAAACCGTGG" * (length // 10)), id="duplicate",
                       annotations={"topology": "circular", "molecule_type": "DNA"})
    record.features = [SeqFeature(SimpleLocation(2000, 3000, strand=1), type="CDS",
                                  qualifiers={"product": ["anchor enzyme"], "locus_tag": ["anchor"]})]
    return record


def _request(mode, record, start, reverse=False, overrides=None):
    config = apply_config_overrides(None, {
        "labels.circular.scope": "outer", "labels.linear.scope": "all",
        "canvas.show_gc": False, "canvas.show_skew": False,
        **(overrides or {}),
    })
    request_type, options_type = ((CircularDiagramRequest, CircularDiagramOptions)
                                 if mode == "circular" else
                                 (LinearDiagramRequest, LinearDiagramOptions))
    return request_type(records=(RecordInput(
        InMemoryRecordSource(record), display=RecordDisplayOptions(start_coordinate=start),
        presentation=RecordPresentation(reverse_complement=reverse),
    ),), options=options_type(config=config))


def _svg(request):
    plan = plan_request(request)
    result = plan.build()
    drawing = result.drawing if hasattr(result, "drawing") else result
    return plan, ET.fromstring(drawing.tostring())


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_plan_rotation_moves_final_label_without_changing_source(mode):
    record = _record()
    before = copy.deepcopy(record)
    roots = [_svg(_request(mode, record, start))[1] for start in (None, 4001, None)]
    def label(root):
        return [ET.tostring(n) for n in root.iter()
                if n.tag.endswith("}text") and "anchor enzyme" in "".join(n.itertext())]
    assert len(label(roots[0])) == len(label(roots[1])) == 1
    assert label(roots[0]) != label(roots[1])
    assert label(roots[0]) == label(roots[2])
    assert record.seq == before.seq
    assert record.features == before.features
    assert record.annotations == before.annotations


def _xy(path):
    return [(float(x), float(y)) for x, y in re.findall(
        r"[ML]\s*([-+\d.eE]+)[ ,]+([-+\d.eE]+)", path)]


def _group(root, name):
    return next(node for node in root.iter() if node.get("id") == name)


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("start", [None, 1, 10000, 4001])
def test_source_tick_text_position_and_wrapped_definition(mode, reverse, start):
    plan, root = _svg(_request(mode, _record(), start, reverse, {
        "canvas.linear.ruler_on_axis": True, "objects.scale.interval": 1000,
        "objects.scale.style": "ruler",
        "canvas.linear.track_layout": "above",
    }))
    definition = next(n for n in root.iter() if n.get("data-gbdraw-role") == "record-definition")
    text = " ".join(definition.itertext())
    if start is None:
        assert "10,000 bp" in text and "[" not in text
    else:
        expected = {
            (False, 1): "[1..10,000] bp", (False, 10000): "[10,000..10,000], [1..9,999] bp",
            (False, 4001): "[4,001..10,000], [1..4,000] bp",
            (True, 1): "[1..1], [10,000..2] bp", (True, 10000): "[10,000..1] bp",
            (True, 4001): "[4,001..1], [10,000..4,002] bp",
        }[reverse, start]
        assert expected in text
    if mode == "linear":
        record_group = next(n for n in root.iter() if n.get("data-gbdraw-record-id") == "duplicate"
                            and n.get("data-gbdraw-role") is None)
        axis = next(n for n in record_group if n.tag.endswith("}line") and n.get("y1") == n.get("y2"))
        width = float(axis.get("x2"))
        tick = next(n for n in record_group if n.tag.endswith("}text") and "".join(n.itertext()) == "2 kbp")
        if start is not None:
            index = ((start - 2000) if reverse else (2000 - start)) % 10000
            assert float(tick.get("x")) == pytest.approx(width * index / 10000)
    else:
        ticks = _group(root, "tick")
        labels = [n for n in ticks.iter() if n.tag.endswith("}textPath")]
        assert any("".join(n.itertext()) == "2 kbp" for n in labels)
        if start is not None:
            points = _xy([n for n in ticks if n.tag.endswith("}path")][2].get("d"))
            angle = math.radians(360 * (((start - 2000) if reverse else (2000 - start)) % 10000) / 10000 - 90)
            radius = math.hypot(*points[0])
            assert points[0] == pytest.approx((radius * math.cos(angle), radius * math.sin(angle)))
    assert plan.transforms[0].start_coordinate == start


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("start", [1, 4001, 5001, 10000])
def test_plan_series_preserve_samples_and_open_svg_seam(mode, reverse, start, monkeypatch):
    import gbdraw.api.diagram as api_owner
    import gbdraw.diagrams.linear.assemble as linear_owner
    import gbdraw.svg.circular_tracks as circular_paths
    import gbdraw.svg.linear_tracks as linear_paths
    frame = pd.DataFrame({"GC content": [.1, .3, .8, .6], "GC skew": [.4, -.6, .2, -.2]},
                         index=[0, 2500, 5000, 7500])
    before = frame.copy(deep=True)
    monkeypatch.setattr(api_owner, "skew_df", lambda *a, **k: frame)
    monkeypatch.setattr(api_owner, "_build_circular_dinucleotide_content_dataframes", lambda *a, **k: {"GC": frame})
    monkeypatch.setattr(api_owner, "_build_circular_dinucleotide_skew_dataframes", lambda *a, **k: {"GC": frame})
    monkeypatch.setattr(linear_owner, "skew_df", lambda *a, **k: frame)
    owner = circular_paths if mode == "circular" else linear_paths
    project = owner.project_scalar_samples
    observed = []
    def capture(positions, values, transform, **kwargs):
        values = list(values)
        result = project(positions, values, transform, **kwargs)
        observed.append((list(positions), values, result, kwargs))
        return result
    monkeypatch.setattr(owner, "project_scalar_samples", capture)
    depth = pd.DataFrame({"reference_name": ["duplicate"] * 10000,
                          "position": range(1, 10001),
                          "depth": [10.] * 2500 + [30.] * 2500 + [80.] * 2500 + [60.] * 2500})
    depth_before = depth.copy(deep=True)
    request = _request(mode, _record(), start, reverse, {"canvas.show_gc": True, "canvas.show_skew": True})
    request = replace(request, options=replace(request.options, depth_table=depth, depth_window=2500, depth_step=2500))
    _, root = _svg(request)
    assert len(observed) == 3
    for positions, values, segments, kwargs in observed:
        segment, = segments
        assert positions == [0, 2500, 5000, 7500]
        assert segment.points[0].position == 0
        assert segment.points[-1].position == 10000
        assert segment.points[0].value == segment.points[-1].value
        assert all(a.position < b.position for a, b in zip(segment.points, segment.points[1:]))
        seen = set()
        for sample, index in zip(segment.points, segment.source_indices, strict=True):
            if index is not None:
                assert sample.value == values[index]
                seen.add(index)
        assert seen == {0, 1, 2, 3}
        # Independent boundary oracle: GC is already record-local after RC;
        # depth_df adapted the original source positions from 1-based input.
        local_cut = (10000 - start) if reverse else start - 1
        cut = (start if reverse else start - 1) if kwargs.get("source_positions") else local_cut
        sampled = cut % 2500 == 0
        assert segment.seam_sampled == sampled
        left = (cut // 2500) % 4
        fraction = (cut % 2500) / 2500
        assert segment.points[0].value == pytest.approx(values[left] * (1-fraction) + values[(left+1) % 4] * fraction)
        assert (segment.source_indices[0] is None) == (not sampled)
    for name in ("gc_content", "gc_skew", "depth"):
        group_name = "skew" if mode == "circular" and name == "gc_skew" else name
        paths = [n for n in _group(root, group_name) if n.tag.endswith("}path")]
        assert paths
        for path in paths:
            points = _xy(path.get("d"))
            assert len(points) >= 5
            cut = ((start if reverse else start - 1) if name == "depth"
                   else ((10000 - start) if reverse else start - 1))
            direction = -1 if name == "depth" and reverse else 1
            offsets = sorted({0, 10000, *[(direction * (p - cut)) % 10000 for p in (0, 2500, 5000, 7500)]})
            if mode == "linear":
                # Only baseline closure returns from the right edge to the left.
                assert all(b[0] >= a[0] or a[1] == pytest.approx(b[1])
                           for a, b in zip(points, points[1:]))
                width = max(x for x, _ in points)
                for offset in offsets:
                    assert any(x == pytest.approx(offset * width / 10000) for x, _ in points)
            elif name != "depth":
                # The outer ring meets itself at the duplicated seam value.
                outer = _xy(path.get("d").split("z")[0])
                assert outer[0] == pytest.approx(outer[-1], abs=1e-7)
                for offset in offsets:
                    angle = math.radians(360 * offset / 10000 - 90)
                    assert any((x / math.hypot(x, y), y / math.hypot(x, y)) == pytest.approx(
                        (math.cos(angle), math.sin(angle))) for x, y in outer)
    pd.testing.assert_frame_equal(frame, before)
    pd.testing.assert_frame_equal(depth, depth_before)


@pytest.mark.parametrize("mode,placement", [("circular", "horizontal"), ("circular", "radial"),
                                          ("linear", "auto"), ("linear", "above_feature")])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("shape", ["ordinary", "cut", "multipart", "origin", "tie", "empty"])
def test_plan_labels_select_projected_coverage_and_connect_leaders(mode, placement, reverse, shape, monkeypatch):
    from Bio.SeqFeature import CompoundLocation
    import gbdraw.diagrams.circular.assemble as circular
    import gbdraw.diagrams.linear.precalc as linear
    owner = circular if mode == "circular" else linear
    symbol = "prepare_label_list" if mode == "circular" else "prepare_label_list_linear"
    prepare = getattr(owner, symbol)
    captured = []
    def capture(*args, **kwargs):
        result = prepare(*args, **kwargs)
        captured.append((args[0], result))
        return result
    monkeypatch.setattr(owner, symbol, capture)
    intervals = {"ordinary": [(2000, 3000)], "cut": [(2000, 7000)],
                 "multipart": [(2000, 3000), (6000, 8000)], "origin": [(9000, 10000), (0, 1000)],
                 "tie": [(2000, 6000)], "empty": [(4000, 4000)]}[shape]
    record = _record()
    parts = [SimpleLocation(a, b, strand=1) for a, b in intervals]
    record.features[0].location = parts[0] if len(parts) == 1 else CompoundLocation(parts)
    before = copy.deepcopy(record)
    overrides = {f"labels.{mode}.placement": placement}
    if mode == "linear" and placement == "above_feature":
        overrides["labels.linear.rotation"] = 45
    elif shape == "ordinary":
        overrides["labels.rendering"] = "external_only"
    plan, root = _svg(_request(mode, record, 4001, reverse, overrides))
    features, labels = captured[-1]
    feature = next(iter(features.values()))
    blocks = [p.fragment for p in feature.display_parts if p.kind == "block"]
    if shape == "empty":
        assert not labels
        return
    assert len(labels) == 1
    label, = labels
    if mode == "linear":
        expected = {
            (False, "ordinary"): 8500, (True, "ordinary"): 1501,
            (False, "cut"): 1500, (True, "cut"): 8500.5,
            (False, "multipart"): 3000, (True, "multipart"): 7001,
            (False, "origin"): 5500, (True, "origin"): 3501,
            (False, "tie"): 1000, (True, "tie"): 1000.5,
        }[reverse, shape]
        record_group = next(n for n in root.iter() if n.get("data-gbdraw-record-id") == "duplicate"
                            and n.get("data-gbdraw-role") is None)
        axis = next(n for n in record_group if n.tag.endswith("}line") and n.get("y1") == n.get("y2"))
        scale = float(axis.get("x2")) / 10000
        assert label["feature_anchor_x"] == pytest.approx(expected * scale)
        assert any(f.display_start <= expected <= f.display_end for f in blocks)
        from gbdraw.labels.linear import calculate_label_bounds
        left, right, top, bottom = calculate_label_bounds(label)
        assert left < right and top < bottom
    else:
        expected = {"ordinary": 8500, "cut": 500, "multipart": 3000, "origin": 6000, "tie": 0}[shape]
        if reverse:
            expected = (1 - expected) % 10000
        assert label["middle"] == pytest.approx(expected)
        assert any(f.display_start <= expected <= f.display_end for f in blocks)
        if not label["is_embedded"]:
            assert any(n.tag.endswith("}line") or (n.tag.endswith("}path") and "L" in n.get("d", ""))
                       for n in root.iter())
    assert sum("anchor enzyme" in "".join(n.itertext()) for n in root.iter() if n.tag.endswith("}text")) == 1
    assert record.seq == before.seq and record.features == before.features and record.annotations == before.annotations
    assert plan.records[0].id == record.id


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_consumer_cache_a_b_a_preserves_parse_identity_and_entire_svg(mode, monkeypatch, tmp_path):
    from Bio import SeqIO
    from gbdraw.api.prepared import PreparedBiologicalInputCache, PreparedResourceIdentity
    from gbdraw.api.requests import GenBankInputSource
    import gbdraw.api.request_render as owner
    path = tmp_path / "source.gb"
    SeqIO.write(_record(), path, "genbank")
    cache = PreparedBiologicalInputCache()
    identity = PreparedResourceIdentity("source", "source-resource", path.stat().st_size)
    load = owner.load_gbks
    calls = []
    def capture(*args, **kwargs):
        calls.append(1)
        return load(*args, **kwargs)
    monkeypatch.setattr(owner, "load_gbks", capture)
    snapshots, plans = [], []
    for start in (2001, 4001, 2001):
        request = _request(mode, _record(), start, overrides={"canvas.show_gc": True, "canvas.show_skew": True,
            "canvas.linear.ruler_on_axis": True, "canvas.linear.track_layout": "above", "objects.scale.style": "ruler"})
        request = replace(request, records=(replace(request.records[0], source=GenBankInputSource(path)),))
        with cache.transaction(resource_paths={path: identity}, diagnostics=None):
            plan, root = _svg(request)
            plans.append(plan)
            snapshots.append(ET.tostring(root))
    assert calls == [1]
    assert plans[0].records[0] is plans[1].records[0] is plans[2].records[0]
    assert snapshots[0] != snapshots[1]
    assert snapshots[0] == snapshots[2]


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("percent", [False, True])
def test_rotation_preserves_physical_and_value_axes_and_definition_visibility(mode, percent):
    overrides = {"canvas.show_gc": True, "canvas.show_skew": True,
                 "objects.gc_content.mode": "percent" if percent else "deviation"}
    roots = [_svg(_request(mode, _record(), start, overrides=overrides))[1] for start in (1, 4001)]
    axes = lambda root: [[ET.tostring(child) for child in n] for n in root.iter()
                        if n.get("id", "").endswith("_axis") or n.get("id") == "length_bar"]
    assert axes(roots[0]) == axes(roots[1])
    if mode == "linear":
        _, root = _svg(_request(mode, _record(), 4001, overrides={"objects.definition.linear.show_length": False}))
        assert not any(n.get("data-definition-line-kind") == "length" for n in root.iter())


@pytest.mark.parametrize("mode,layout_kind", [("circular", "grid"), ("circular", "batch"),
                                            ("linear", "rows"), ("linear", "multi")])
def test_multi_record_consumers_keep_instance_order_and_sparse_depth(mode, layout_kind, monkeypatch, tmp_path):
    from gbdraw.api.options import CircularMultiRecordOptions, LinearMultiRecordOptions
    import gbdraw.svg.circular_tracks as circular_paths
    import gbdraw.svg.linear_tracks as linear_paths
    owner = circular_paths if mode == "circular" else linear_paths
    project = owner.project_scalar_samples
    observed = []
    def capture(positions, values, transform, **kwargs):
        observed.append((transform.length, transform.start_coordinate, transform.source_step, kwargs.get("source_positions", False)))
        return project(positions, values, transform, **kwargs)
    monkeypatch.setattr(owner, "project_scalar_samples", capture)
    first = _request(mode, _record(), 4001, overrides={"canvas.show_gc": True, "canvas.show_skew": True})
    second = replace(first.records[0], source=InMemoryRecordSource(_record(12000)), record_key="second",
                     display=RecordDisplayOptions(start_coordinate=12000), presentation=RecordPresentation(reverse_complement=True))
    tables = [pd.DataFrame({"reference_name": ["duplicate"], "position": [1], "depth": [value]}) for value in (20., 80.)]
    options = replace(first.options, depth_track_tables=((tables[0], None), (None, tables[1])),
                      depth_track_colors=("#123456", "#654321"), depth_track_labels=("first depth", "second depth"),
                      depth_window=1000, depth_step=1000)
    layout = (CircularMultiRecordOptions() if mode == "circular" else
              (LinearMultiRecordOptions(multi_record_positions=("#1@1", "#2@1")) if layout_kind == "multi" else None))
    kwargs = {"grouping": layout_kind} if mode == "circular" else {}
    if layout_kind == "batch":
        from gbdraw.api.requests import CircularBatchRequest, CircularBatchOutputPolicy
        request = CircularBatchRequest(records=(first.records[0], second), options=options,
                                       output_policy=CircularBatchOutputPolicy(output_prefix=tmp_path / "batch"))
    else:
        request = replace(first, records=(first.records[0], second), options=options, layout=layout, **kwargs)
    plan = plan_request(request)
    built = ([p.build() for p in plan.item_plans()] if layout_kind == "batch" else [plan.build()])
    roots = [ET.fromstring((item.drawing if hasattr(item, "drawing") else item).tostring()) for item in built]
    assert {(length, start, step) for length, start, step, _ in observed} == {(10000, 4001, 1), (12000, 12000, -1)}
    depth = [row for row in observed if row[3]]
    assert depth == [(10000, 4001, 1, True), (12000, 12000, -1, True)]
    defs = [n for root in roots for n in root.iter() if n.get("data-gbdraw-role") == "record-definition"]
    assert any("[4,001..10,000], [1..4,000] bp" in " ".join(n.itertext()) for n in defs)
    assert any("[12,000..1] bp" in " ".join(n.itertext()) for n in defs)
    paths = [n for root in roots for n in root.iter() if n.tag.endswith("}path")]
    assert any(n.get("fill") == "#123456" for n in paths)
    assert any(n.get("fill") == "#654321" for n in paths)
    assert len(plan.records) == len(plan.transforms) == 2


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_empty_gc_series_pass_through_real_plan_without_phantom_data(mode, monkeypatch):
    import gbdraw.api.diagram as api_owner
    import gbdraw.diagrams.linear.assemble as linear_owner
    empty = pd.DataFrame({"GC content": [], "GC skew": []})
    monkeypatch.setattr(api_owner, "_build_circular_dinucleotide_content_dataframes", lambda *a, **k: {"GC": empty})
    monkeypatch.setattr(api_owner, "_build_circular_dinucleotide_skew_dataframes", lambda *a, **k: {"GC": empty})
    monkeypatch.setattr(linear_owner, "skew_df", lambda *a, **k: empty)
    _, root = _svg(_request(mode, _record(), 4001, overrides={"canvas.show_gc": True, "canvas.show_skew": True}))
    assert all("nan" not in n.get("d", "").lower() for n in root.iter())
    assert all(not list(n) for n in root.iter() if n.get("id") in {"gc_content", "gc_skew", "skew"})


@pytest.mark.parametrize("mode,placement", [("circular", "horizontal"), ("circular", "radial"), ("linear", "auto")])
def test_dense_rotated_labels_keep_collision_bounds_and_actual_leader_endpoints(mode, placement, monkeypatch):
    import gbdraw.diagrams.circular.assemble as circular
    import gbdraw.diagrams.linear.precalc as linear
    from gbdraw.labels.circular import _count_label_overlaps
    from gbdraw.labels.linear import calculate_label_bounds
    owner = circular if mode == "circular" else linear
    symbol = "prepare_label_list" if mode == "circular" else "prepare_label_list_linear"
    real = getattr(owner, symbol)
    captured = []
    def capture(*args, **kwargs):
        result = real(*args, **kwargs)
        captured.append(result)
        return result
    monkeypatch.setattr(owner, symbol, capture)
    record = _record()
    record.features = [SeqFeature(SimpleLocation(300 + 500*i, 510 + 500*i, strand=1), type="CDS",
                                  qualifiers={"product": [f"enzyme {i:02d}"], "locus_tag": [f"cds{i}"]})
                       for i in range(18)]
    _, root = _svg(_request(mode, record, 4001, overrides={
        f"labels.{mode}.placement": placement, "labels.rendering": "external_only",
        "canvas.show_gc": True, "canvas.show_skew": True,
    }))
    labels = captured[-1]
    assert len(labels) == 18
    lines = [((float(n.get("x1")), float(n.get("y1"))), (float(n.get("x2")), float(n.get("y2"))))
             for n in root.iter() if n.tag.endswith("}line")]
    if mode == "circular":
        assert _count_label_overlaps(labels, 10000, use_min_gap=False) == 0
        # Leader groups stay in record-local coordinates before canvas translation.
        endpoints = [p for n in _group(root, "label_leaders").iter() if n.tag.endswith("}path")
                     for p in _xy(n.get("d", ""))] + [p for pair in lines for p in pair]
        for label in labels:
            expected = (label["feature_anchor_x"], label["feature_anchor_y"])
            assert any(point == pytest.approx(expected) for point in endpoints)
    else:
        boxes = [calculate_label_bounds(label) for label in labels]
        for i, (left, right, top, bottom) in enumerate(boxes):
            assert all(right <= a or left >= b or bottom <= c or top >= d
                       for a, b, c, d in boxes[i+1:])
        for label in labels:
            expected = (label["feature_anchor_x"], label["feature_middle_y"])
            assert any(point == pytest.approx(expected) for pair in lines for point in pair)


def test_circular_tick_radial_reservation_is_rotation_invariant():
    from gbdraw.svg.circular_ticks import resolve_circular_tick_label_geometry
    # The baseline/path direction change cancels radially. Existing radial
    # reservations therefore remain valid for every mapped source tick angle.
    bounds = []
    for tick in (0, 1, 2499, 2500, 5000, 7499, 7500, 9999):
        geometry = resolve_circular_tick_label_geometry(center_radius_px=300, total_len=10000,
            size="large", tick=tick, label_text="2 kbp", font_size=14, font_family="Arial",
            track_type="tuckin", strandedness=False, dpi=96)
        bounds.append((geometry.radial_inner_px, geometry.radial_outer_px))
    assert all(bound == pytest.approx(bounds[0]) for bound in bounds)
