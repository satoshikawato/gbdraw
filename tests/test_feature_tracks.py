from __future__ import annotations

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
