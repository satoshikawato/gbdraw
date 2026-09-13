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
from gbdraw.labels.circular_radial import (
    _leader_crossing_count,
    _leader_text_collision_count,
    _text_collision_pairs,
)


SESSION = Path(__file__).parents[1] / "gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json"


@pytest.fixture
def chloroplast_session():
    return json.loads(SESSION.read_text())


@pytest.mark.parametrize("separate", [False, True])
@pytest.mark.parametrize("scope", ["both", "outer"])
@pytest.mark.parametrize("reverse_placement", [False, True])
def test_chloroplast_placement_preserves_lanes_and_clears_labels_and_legend(
    chloroplast_session, separate, scope, reverse_placement, monkeypatch, tmp_path,
):
    document = copy.deepcopy(chloroplast_session)
    options = document["renderRequest"]["diagramOptions"]
    options["config"]["canvas"]["strandedness"] = separate
    options["config"]["labels"]["circular"]["scope"] = scope
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
    original = assemble.prepare_label_list

    def capture(features, *args, **kwargs):
        labels = original(features, *args, **kwargs)
        measured.append((features, kwargs["radial_layout"], labels))
        return labels

    monkeypatch.setattr(assemble, "prepare_label_list", capture)
    with materialize_session(document, output_directory=tmp_path) as materialized:
        drawing = plan_request(session_to_request(materialized)).build()
    result = drawing._gbdraw_circular_assembly_result
    features, radial, labels = measured[-1]
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
        assert _leader_crossing_count(side_labels) == 0
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
