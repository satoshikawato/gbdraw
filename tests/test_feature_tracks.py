from __future__ import annotations

import pytest

from gbdraw.exceptions import ValidationError
from gbdraw.features.objects import FeatureLocationPart, FeatureObject
from gbdraw.features.tracks import arrange_feature_tracks, get_feature_ends


def _make_feature(
    feature_id: str,
    strand: str,
    segments: list[tuple[int, int]],
) -> FeatureObject:
    location = [
        FeatureLocationPart("block", str(i).zfill(3), strand, start, end, i == len(segments))
        for i, (start, end) in enumerate(segments, start=1)
    ]
    return FeatureObject(
        feature_id=feature_id,
        location=location,
        is_directional=True,
        color="#54bcf8",
        note="",
        label_text="",
        coordinates=location,
        type="CDS",
        qualifiers={},
    )


def test_get_feature_ends_negative_multipart_is_not_origin_spanning() -> None:
    feature = _make_feature(
        feature_id="f1",
        strand="negative",
        # Descending order can occur on negative strand but this does not cross origin.
        segments=[(900, 920), (700, 730)],
    )

    start, end, strand = get_feature_ends(feature, genome_length=1000)

    assert strand == "negative"
    assert start == 700
    assert end == 920
    assert start < end


def test_get_feature_ends_detects_origin_spanning_when_touching_both_edges() -> None:
    feature = _make_feature(
        feature_id="f1",
        strand="negative",
        segments=[(990, 1000), (1, 20)],
    )

    start, end, strand = get_feature_ends(feature, genome_length=1000)

    assert strand == "negative"
    assert start == 990
    assert end == 20
    assert start > end


def test_arrange_feature_tracks_avoids_false_displacement_for_negative_multipart() -> None:
    feature_dict = {
        "a": _make_feature("a", "negative", [(900, 920), (700, 730)]),
        "b": _make_feature("b", "negative", [(10, 20)]),
        "c": _make_feature("c", "negative", [(40, 60)]),
    }

    arranged = arrange_feature_tracks(
        feature_dict=feature_dict,
        separate_strands=True,
        resolve_overlaps=True,
        genome_length=1000,
    )

    assert arranged["a"].feature_track_id == -1
    assert arranged["b"].feature_track_id == -1
    assert arranged["c"].feature_track_id == -1


def test_arrange_feature_tracks_non_stranded_resolve_keeps_legacy_indices_when_split_disabled() -> None:
    feature_dict = {
        "a": _make_feature("a", "positive", [(100, 300)]),
        "b": _make_feature("b", "negative", [(120, 280)]),
        "c": _make_feature("c", "negative", [(130, 260)]),
        "d": _make_feature("d", "positive", [(140, 240)]),
        "e": _make_feature("e", "negative", [(700, 730)]),
    }

    arranged = arrange_feature_tracks(
        feature_dict=feature_dict,
        separate_strands=False,
        resolve_overlaps=True,
        split_overlaps_by_strand=False,
        genome_length=1000,
    )

    assert arranged["a"].feature_track_id == 0
    assert arranged["b"].feature_track_id == 1
    assert arranged["c"].feature_track_id == 2
    assert arranged["d"].feature_track_id == 3
    assert arranged["e"].feature_track_id == 0


def test_arrange_feature_tracks_non_stranded_resolve_splits_inner_outer_indices() -> None:
    feature_dict = {
        "a": _make_feature("a", "positive", [(100, 300)]),
        "b": _make_feature("b", "negative", [(120, 280)]),
        "c": _make_feature("c", "negative", [(130, 260)]),
        "d": _make_feature("d", "positive", [(140, 240)]),
        "e": _make_feature("e", "negative", [(700, 730)]),
    }

    arranged = arrange_feature_tracks(
        feature_dict=feature_dict,
        separate_strands=False,
        resolve_overlaps=True,
        split_overlaps_by_strand=True,
        genome_length=1000,
    )

    # Center track remains shared across strands.
    assert arranged["a"].feature_track_id == 0
    assert arranged["e"].feature_track_id == 0

    # Inner/outer displacement tracks are resolved independently.
    assert arranged["b"].feature_track_id == -1
    assert arranged["c"].feature_track_id == -2
    assert arranged["d"].feature_track_id == 1


def test_arrange_feature_tracks_indexes_origin_spanning_features_conservatively() -> None:
    feature_dict = {
        "origin": _make_feature("origin", "positive", [(990, 1000), (1, 20)]),
        "middle": _make_feature("middle", "positive", [(500, 520)]),
        "left": _make_feature("left", "positive", [(10, 15)]),
    }

    arranged = arrange_feature_tracks(
        feature_dict=feature_dict,
        separate_strands=False,
        resolve_overlaps=True,
        genome_length=1000,
    )

    assert arranged["origin"].feature_track_id == 0
    assert arranged["middle"].feature_track_id == 0
    assert arranged["left"].feature_track_id == 1


def test_arrange_feature_tracks_treats_touching_intervals_as_non_overlapping() -> None:
    feature_dict = {
        "left": _make_feature("left", "positive", [(100, 200)]),
        "right": _make_feature("right", "positive", [(200, 300)]),
    }

    arranged = arrange_feature_tracks(
        feature_dict=feature_dict,
        separate_strands=False,
        resolve_overlaps=True,
        genome_length=1000,
    )

    assert arranged["left"].feature_track_id == 0
    assert arranged["right"].feature_track_id == 0


def test_arrange_feature_tracks_origin_span_can_touch_neighbor_without_overlap() -> None:
    feature_dict = {
        "origin": _make_feature("origin", "positive", [(990, 1000), (1, 20)]),
        "neighbor": _make_feature("neighbor", "positive", [(20, 30)]),
    }

    arranged = arrange_feature_tracks(
        feature_dict=feature_dict,
        separate_strands=False,
        resolve_overlaps=True,
        genome_length=1000,
    )

    assert arranged["origin"].feature_track_id == 0
    assert arranged["neighbor"].feature_track_id == 0


def test_arrange_feature_tracks_raises_when_all_feature_lanes_are_occupied() -> None:
    feature_dict = {
        f"feature_{index:03d}": _make_feature(
            f"feature_{index:03d}",
            "positive",
            [(100, 200)],
        )
        for index in range(101)
    }

    with pytest.raises(ValidationError, match=r"all 100 feature lanes are occupied"):
        arrange_feature_tracks(
            feature_dict=feature_dict,
            separate_strands=False,
            resolve_overlaps=True,
            genome_length=1000,
        )


def test_shared_overlap_length_and_tolerance_boundary():
    from gbdraw.features.tracks import feature_overlap_bp, features_conflict
    a = {"start": 990, "end": 20, "strand": "positive"}
    b = {"start": 995, "end": 10, "strand": "positive"}
    assert feature_overlap_bp(a, b, separate_strands=False, genome_length=1000) == 14
    assert not features_conflict(a, b, tolerance_bp=14, genome_length=1000)
    assert features_conflict(a, b, tolerance_bp=13, genome_length=1000)


@pytest.mark.parametrize("tolerance", [19, 20])
def test_undefined_negative_collision_uses_one_pool_in_both_orders(tolerance):
    from gbdraw.features.tracks import feature_overlap_bp, features_conflict, find_best_track
    for first, second in [("negative", "undefined"), ("undefined", "negative")]:
        a = {"start": 20, "end": 50, "strand": first, "id": "a"}
        b = {"start": 30, "end": 50, "strand": second, "id": "b"}
        assert feature_overlap_bp(a, b, separate_strands=True, genome_length=120) == 20
        assert features_conflict(a, b, separate_strands=True, tolerance_bp=tolerance,
                                 genome_length=120) == (tolerance < 20)
        assert find_best_track(b, {"track_1": [a]}, True, True, genome_length=120,
                               tolerance_bp=tolerance) == (-2 if tolerance < 20 else -1)
        positive = dict(b, strand="positive")
        assert not features_conflict(a, positive, separate_strands=True,
                                     genome_length=120, tolerance_bp=tolerance)


@pytest.mark.parametrize("overlap", [0, 1, 2, 3])
@pytest.mark.parametrize("tolerance", [0, 1, 2])
@pytest.mark.parametrize("split,separate,strand", [
    (False, False, "positive"), (True, False, "negative"),
    (True, False, "positive"), (False, True, "negative"),
    (False, True, "positive"),
])
def test_every_auto_path_uses_shared_tolerance(overlap, tolerance, split, separate, strand):
    from gbdraw.features.tracks import find_best_track
    features = {
        "a": _make_feature("a", strand, [(100, 200)]),
        "b": _make_feature("b", strand, [(200-overlap, 220)]),
    }
    arrange_feature_tracks(features, separate, True, split, 1000, tolerance_bp=tolerance)
    nominal = -1 if separate and strand == "negative" else 0
    step = -1 if strand == "negative" and (separate or split) else 1
    assert features["a"].feature_track_id == nominal
    assert features["b"].feature_track_id == nominal + (step if overlap > tolerance else 0)
    # The existing public/unindexed path must use the same final predicate too.
    a = {"start": 100, "end": 200, "strand": strand, "id": "a"}
    b = {"start": 200-overlap, "end": 220, "strand": strand, "id": "b"}
    expected = nominal + ((-1 if nominal == -1 else 1) if overlap > tolerance else 0)
    assert find_best_track(b, {f"track_{abs(nominal)}": [a]}, separate, True,
                           genome_length=1000, tolerance_bp=tolerance) == expected


def test_fixed_reservations_keep_index_candidate_filter(monkeypatch):
    import gbdraw.features.tracks as owner
    calls = []
    actual = owner.features_conflict

    def capture(a, b, **kwargs):
        calls.append((a["id"], b["id"]))
        return actual(a, b, **kwargs)

    monkeypatch.setattr(owner, "features_conflict", capture)
    features = {
        str(i): _make_feature(str(i), "positive", [(100+i*100, 110+i*100)])
        for i in range(1000)
    }
    owner.arrange_feature_tracks(features, False, True, genome_length=110000,
                                fixed_tracks={str(i): 0 for i in range(500)})
    assert len(calls) < 5000  # Far below 499,500 full pair comparisons.
    assert {f.feature_track_id for f in features.values()} == {0}


@pytest.mark.parametrize("separate", [False, True])
@pytest.mark.parametrize("split", [False, True])
def test_auto_retains_all_100_lanes_with_fixed_first(separate, split):
    features = {
        str(i): _make_feature(str(i), "positive", [(100, 200)])
        for i in range(100)
    }
    arrange_feature_tracks(features, separate, True, split, 1000, fixed_tracks={"99": 0})
    assert features["99"].feature_track_id == 0
    assert {f.feature_track_id for f in features.values()} == set(range(100))
    features["extra"] = _make_feature("extra", "positive", [(100, 200)])
    with pytest.raises(ValidationError, match="all 100 feature lanes"):
        arrange_feature_tracks(features, separate, True, split, 1000, fixed_tracks={"99": 0})
