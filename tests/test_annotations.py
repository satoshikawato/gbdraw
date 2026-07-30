from __future__ import annotations

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.annotations import (
    AnnotationOptions,
    AnnotationSet,
    AnnotationTrackParams,
    CoordinateSpan,
    FeatureSelector,
    FeatureSpan,
    RegionAnnotation,
    RegionAnnotationStyle,
    annotation_sets_from_dataframe,
    annotation_track_params_from_mapping,
    assign_annotation_lanes,
    effective_annotation_style,
    layout_annotation_track,
    resolve_annotations,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import parse_record_selector


def _record(length: int = 100, record_id: str = "r1") -> SeqRecord:
    return SeqRecord(Seq("A" * length), id=record_id, name=record_id)


def test_coordinate_contract_is_one_based_inclusive() -> None:
    annotation = RegionAnnotation("a", CoordinateSpan(None, 2, 5))
    bundle = resolve_annotations((AnnotationSet("s", (annotation,)),), [_record()], mode="linear")
    assert bundle.annotations[0].segments == ((1, 5),)
    assert bundle.annotations[0].span_bp == 4


def test_highlight_is_a_supported_region_mark() -> None:
    annotation = RegionAnnotation(
        "highlighted",
        CoordinateSpan(None, 2, 5),
        mark="highlight",
    )

    assert annotation.mark == "highlight"


def test_highlight_track_flattens_overlaps_and_ignores_explicit_lanes() -> None:
    annotations = (
        RegionAnnotation("first", CoordinateSpan(None, 2, 10), mark="highlight", lane=0),
        RegionAnnotation("second", CoordinateSpan(None, 5, 15), mark="highlight", lane=0),
    )
    bundle = resolve_annotations(
        (AnnotationSet("s", annotations),),
        [_record()],
        mode="linear",
    )
    params = annotation_track_params_from_mapping(
        {"set_id": "s", "marks": "highlight", "cover_anchor": "true"}
    )
    track = layout_annotation_track(
        "highlight",
        "s",
        bundle.annotations,
        record_lengths={0: 100},
        params=params,
    )

    assert params.cover_anchor is True
    assert {placement.lane for placement in track.placements} == {0}
    assert track.required_extent_px == 18.0


def test_source_coordinates_follow_reverse_coordinate_map() -> None:
    record = _record(10)
    record.annotations["gbdraw_coord_base"] = 100
    record.annotations["gbdraw_coord_step"] = -1
    annotation = RegionAnnotation("a", CoordinateSpan(None, 95, 98))
    resolved = resolve_annotations((AnnotationSet("s", (annotation,)),), [record], mode="linear")
    assert resolved.annotations[0].segments == ((2, 6),)


def test_origin_span_splits_and_linear_rejects() -> None:
    annotation = RegionAnnotation("a", CoordinateSpan(None, 90, 10, wraps_origin=True))
    circular = resolve_annotations((AnnotationSet("s", (annotation,)),), [_record()], mode="circular")
    assert circular.annotations[0].segments == ((0, 10), (89, 100))
    with pytest.raises(ValidationError, match="linear"):
        resolve_annotations((AnnotationSet("s", (annotation,)),), [_record()], mode="linear")


def test_record_id_must_be_unique_but_index_disambiguates() -> None:
    records = [_record(record_id="same"), _record(record_id="same")]
    ambiguous = RegionAnnotation("a", CoordinateSpan(parse_record_selector("same"), 1, 2))
    with pytest.raises(ValidationError, match="multiple"):
        resolve_annotations((AnnotationSet("s", (ambiguous,)),), records, mode="linear")
    indexed = RegionAnnotation("a", CoordinateSpan(parse_record_selector("#2"), 1, 2))
    result = resolve_annotations((AnnotationSet("s", (indexed,)),), records, mode="linear")
    assert result.annotations[0].record_index == 1


def test_feature_target_uses_materialized_features_independent_of_visibility() -> None:
    record = _record()
    record.features = [
        SeqFeature(FeatureLocation(10, 20), type="CDS", qualifiers={"locus_tag": ["x"]}),
        SeqFeature(FeatureLocation(30, 35), type="CDS", qualifiers={"locus_tag": ["y"]}),
    ]
    target = FeatureSpan(
        None,
        (FeatureSelector("x", "locus_tag"), FeatureSelector("y", "locus_tag")),
        envelope="outer_bounds",
    )
    annotation = RegionAnnotation("a", target)
    result = resolve_annotations((AnnotationSet("s", (annotation,)),), [record], mode="linear")
    assert result.annotations[0].segments == ((10, 35),)


def test_table_groups_sets_in_first_appearance_order() -> None:
    table = pd.DataFrame(
        [
            {"set_id": "b", "id": "one", "mark": "line", "start": "1", "end": "5"},
            {"set_id": "a", "id": "two", "mark": "band", "start": "8", "end": "10"},
            {"set_id": "b", "id": "three", "mark": "bracket", "feature_selector": "locus_tag=x"},
        ]
    )
    sets = annotation_sets_from_dataframe(table)
    assert [item.id for item in sets] == ["b", "a"]
    assert [item.id for item in sets[0].annotations] == ["one", "three"]


def test_table_rejects_unknown_columns_with_context() -> None:
    table = pd.DataFrame([{"set_id": "s", "id": "a", "mark": "line", "start": 1, "end": 2, "typo": 1}])
    with pytest.raises(ValidationError, match="unknown columns.*typo"):
        annotation_sets_from_dataframe(table)


def test_lane_assignment_is_deterministic_and_separates_overlap() -> None:
    annotations = (
        RegionAnnotation("b", CoordinateSpan(None, 5, 15)),
        RegionAnnotation("a", CoordinateSpan(None, 1, 10)),
        RegionAnnotation("c", CoordinateSpan(None, 20, 25)),
    )
    bundle = resolve_annotations((AnnotationSet("s", annotations),), [_record()], mode="linear")
    lanes = assign_annotation_lanes(bundle.annotations, record_lengths={0: 100})
    assert {item.annotation.id: item.lane for item in lanes} == {"a": 0, "c": 0, "b": 1}


def test_annotation_options_reject_mixed_sources() -> None:
    with pytest.raises(ValidationError, match="not more than one"):
        AnnotationOptions(sets=(AnnotationSet("s", ()),), table=pd.DataFrame())


def test_annotation_style_has_priority_over_track_override() -> None:
    annotation = RegionAnnotation(
        "a",
        CoordinateSpan(None, 1, 5),
        style=RegionAnnotationStyle(stroke="#112233"),
    )
    bundle = resolve_annotations(
        (
            AnnotationSet(
                "s",
                (annotation,),
                default_style=RegionAnnotationStyle(stroke="#445566"),
            ),
        ),
        [_record()],
        mode="linear",
    )
    params = AnnotationTrackParams(
        set_id="s",
        style_override=RegionAnnotationStyle(stroke="#778899"),
    )
    assert effective_annotation_style(bundle.annotations[0], params).stroke == "#112233"
