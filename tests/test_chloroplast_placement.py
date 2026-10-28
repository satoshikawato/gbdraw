"""Regression coverage for editing the published chloroplast session."""

import copy
import json
import math
from pathlib import Path

import pytest

from gbdraw.api import materialize_session, session_to_request
from gbdraw.api.request_render import plan_request
from gbdraw.diagrams.circular import assemble
from gbdraw.features.placement import FeaturePlacementSlot
from gbdraw.labels.circular_leaders import leader_crossing_count
from gbdraw.labels.circular_radial import (
    _leader_text_collision_count,
    _text_collision_pairs,
)


SESSION = Path(__file__).parents[1] / "gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json"


@pytest.fixture
def chloroplast_session():
    return json.loads(SESSION.read_text())


def _build(document, tmp_path):
    # The published Session is the CLI command's one-record batch.
    with materialize_session(document, output_directory=tmp_path) as materialized:
        (item,) = plan_request(session_to_request(materialized)).item_plans()
        return item.build()


@pytest.mark.parametrize("separate", [False, True])
@pytest.mark.parametrize("scope", ["both", "outer"])
@pytest.mark.parametrize("reverse_placement", [False, True])
def test_chloroplast_placement_preserves_lanes_and_clears_labels_and_legend(
    chloroplast_session, separate, scope, reverse_placement, monkeypatch, tmp_path,
):
    document = copy.deepcopy(chloroplast_session)
    options = document["renderRequest"]["diagramOptions"]
    options["configOverrides"]["canvas.strandedness"] = separate
    options["configOverrides"]["labels.circular.scope"] = scope
    selected = {}
    for feature in document["editorState"]["featureCatalog"]["items"][0]["biologicalFeatures"]:
        if feature["type"] != "CDS" or feature["qualifiers"].get("gene", [""])[0] not in {
            "clpP", "rpl16", "rpoC1", "petB",
        }:
            continue
        assert len(feature["location_parts"]) > 1
        outward = (feature["strand"] == "+") != reverse_placement
        side = "outward" if outward else "inward"
        selected[feature["sourceFeatureIndex"]] = side
        options["featurePlacements"].append({
            "recordKey": feature["recordKey"],
            "biologicalFeatureId": feature["biologicalFeatureId"],
            "placement": {"kind": "lane", "side": side, "level": 1},
        })
    assert len(selected) == 4
    measured = []
    drawn = []
    original = assemble.prepare_label_list
    overlaps = assemble._region_labels_on_other_tracks

    def capture(features, *args, **kwargs):
        labels = original(features, *args, **kwargs)
        measured.append((features, labels))
        return labels

    def capture_layout(layout, labels):
        # The layout drawn after the labels: the one they were placed on (the labels take its parts).
        drawn.append(layout)
        return overlaps(layout, labels)

    monkeypatch.setattr(assemble, "prepare_label_list", capture)
    monkeypatch.setattr(assemble, "_region_labels_on_other_tracks", capture_layout)
    drawing = _build(document, tmp_path)
    result = drawing._gbdraw_circular_assembly_result
    features, labels = measured[-1]
    radial = drawn[-1]
    step = radial.features.width_px + radial.axis.radius_px * 0.01
    for feature in features.values():
        assignment = feature.placement
        lane = radial.features.lane_for_track_id(feature.feature_track_id)
        if feature.source_feature_index in selected:
            side = selected[feature.source_feature_index]
            assert (assignment.side, assignment.level) == (side, 1)
            expected = (1 if side == "outward" else -1) * (1.5 if separate else 1) * step
        else:
            assert assignment.level == 0
            expected = (0.5 if feature.strand == "positive" else -0.5) * step if separate else 0
        assert lane.center_px - radial.features.anchor_radius_px == pytest.approx(expected)

    for inner in (False, True):
        side_labels = [label for label in labels if not label["is_embedded"] and label.get("is_inner") == inner]
        assert not _text_collision_pairs(side_labels, 3)
        assert _leader_text_collision_count(side_labels) == 0
        assert leader_crossing_count(side_labels) == 0
        if inner and scope == "outer":
            assert side_labels == []
    for label in labels:
        lane = radial.features.lane_for_track_id(label["track_id"])
        assert math.hypot(label["feature_middle_x"], label["feature_middle_y"]) == pytest.approx(lane.center_px)
    legend = result.composition_plan.placement_for("legend").final_bounds
    assert all(not legend.intersects(obstacle, clearance=8 - 1e-6) for obstacle in result.overlay_obstacles)
    radius = radial.outer_content_radius_px
    # The body, GC track and definition stay clear even with no inner labels.
    assert any(obstacle.width >= 2 * radius and obstacle.height >= 2 * radius
               for obstacle in result.overlay_obstacles)
    assert {t.get("side") for t in FeaturePlacementSlot("circular", "split", separate).supported_targets()} == {
        None, "outward", "inward",
    }


def _render_layouts(document, monkeypatch, tmp_path):
    layouts = []
    original = assemble.resolve_circular_radial_layout

    def capture(*args, **kwargs):
        layout = original(*args, **kwargs)
        layouts.append(layout)
        return layout

    monkeypatch.setattr(assemble, "resolve_circular_radial_layout", capture)
    _build(document, tmp_path)
    return layouts


def _bands(layout):
    # An overlay span has no packing band; its drawn band is what it covers.
    bands = {}
    for slot in layout.slots:
        band = slot.packing_band_px or slot.draw_band_px
        bands[slot.id] = (slot.anchor_radius_px, band.inner_px, band.outer_px)
    return bands


def test_published_chloroplast_stack_keeps_its_bands(chloroplast_session, monkeypatch, tmp_path):
    # The bands that dev 98d4bc13 resolves for the published Session (after the label reflow), as moved by
    # laying labels and tick labels out in their drawn box (OV-165: the feature row +0.78 px, the GC row
    # 0.03 px wider on each side). The region tint spans from the GC track's inner edge to the feature row's
    # outer edge.
    layout = _render_layouts(chloroplast_session, monkeypatch, tmp_path)[-1]
    expected = {
        "features": (407.4970, 386.9346, 428.0595),
        "plastome_tint": (315.0798, 202.1001, 428.0595),
        "plastome_regions": (253.5, 243.5, 263.5),
        "gc_content": (218.4, 202.1001, 234.6999),
    }
    observed = _bands(layout)
    assert observed.keys() == expected.keys()
    for slot_id, bands in expected.items():
        assert observed[slot_id] == pytest.approx(bands, abs=1e-3), slot_id


def test_label_reflow_keeps_outside_rows_beyond_the_moved_feature_row(chloroplast_session, monkeypatch, tmp_path):
    # The Web's automatic stack (custom stack off): radial inner labels enlarge the
    # canvas and move the feature row outside the axis. The outside annotation row
    # stays an outside row beyond it instead of a frozen ring inside the axis.
    document = copy.deepcopy(chloroplast_session)
    options = document["renderRequest"]["diagramOptions"]

    def slot(slot_id, renderer, side, params, width=None):
        return {"kind": "circularTrackSlot", "id": slot_id, "renderer": renderer, "enabled": True, "side": side,
                "radius": None, "width": width, "z": 0, "params": params, "innerGapPx": None, "outerGapPx": None}

    options["tracks"]["circularTrackSlots"] = [
        slot("annotations_1", "annotations", "outside", {"set_id": "plastome_regions", "marks": ["bracket"]}),
        slot("features", "features", "inside", {"lane_direction": "inside"}, {"value": 16, "unit": "px"}),
        slot("ticks", "ticks", "inside", {"tick_label_layout": "label_in_tick_out"}),
        slot("gc_content", "dinucleotide_content", "inside", {"nt": "GC"}),
    ]
    options["tracks"]["circularTrackAxisIndex"] = 1
    for annotation in options["annotations"]["sets"][0]["annotations"]:
        annotation["target"] = {"kind": "featureSpan", "record": annotation["target"]["record"],
                                "selectors": [{"key": "gene", "value": "PRIVATE-MISSING"}],
                                "envelope": "outer_bounds", "circularPath": "shortest"}

    layouts = _render_layouts(document, monkeypatch, tmp_path)
    final = layouts[-1]
    by_id = {slot.id: slot for slot in final.slots}
    assert final.axis.radius_px > layouts[0].axis.radius_px  # the label reflow ran
    assert by_id["features"].side == "outside"
    assert by_id["features"].packing_band_px.inner_px >= final.axis.radius_px
    assert by_id["annotations_1"].packing_band_px.inner_px >= by_id["features"].packing_band_px.outer_px
