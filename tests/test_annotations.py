from __future__ import annotations

import csv
import io
import json
import re
from pathlib import Path

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
    read_annotation_table,
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
        label_font_size=10.0,
    )

    assert params.cover_anchor is True
    assert {placement.lane for placement in track.placements} == {0}
    assert track.required_extent_px == 18.0


def test_annotation_track_params_normalize_both_slot_references() -> None:
    params = annotation_track_params_from_mapping(
        {"set_id": "s", "anchor_slot": " features ", "to_slot": " gc "}
    )

    assert (params.anchor_slot, params.to_slot) == ("features", "gc")
    assert annotation_track_params_from_mapping({"set_id": "s"}).to_slot is None


@pytest.mark.parametrize("name", ["anchor_slot", "to_slot"])
def test_annotation_track_params_reject_an_empty_slot_reference(name: str) -> None:
    with pytest.raises(ValidationError, match=f"^{name} cannot be empty"):
        annotation_track_params_from_mapping({"set_id": "s", name: " "})


def test_spanning_annotation_track_does_not_fit_lanes_to_its_own_extent() -> None:
    annotations = (
        RegionAnnotation("first", CoordinateSpan(None, 2, 10), mark="band"),
        RegionAnnotation("second", CoordinateSpan(None, 5, 15), mark="band"),
    )
    bundle = resolve_annotations(
        (AnnotationSet("s", annotations),),
        [_record()],
        mode="linear",
    )

    def layout(**extra: str):
        return layout_annotation_track(
            "regions",
            "s",
            bundle.annotations,
            record_lengths={0: 100},
            params=annotation_track_params_from_mapping(
                {"set_id": "s", "show_labels": False, "anchor_slot": "features", **extra}
            ),
            available_extent_px=14.0,
            label_font_size=10.0,
        )

    with pytest.raises(ValidationError, match="only 14.00px is available"):
        layout()
    spanning = layout(to_slot="gc")
    assert {placement.lane for placement in spanning.placements} == {0, 1}
    assert spanning.clipped is False


def test_annotation_marks_and_spacing_options_change_selected_layout_geometry() -> None:
    annotations = (
        RegionAnnotation("line", CoordinateSpan(None, 1, 20), mark="line"),
        RegionAnnotation("band", CoordinateSpan(None, 5, 25), mark="band"),
        RegionAnnotation("bracket", CoordinateSpan(None, 10, 30), mark="bracket"),
    )
    bundle = resolve_annotations(
        (AnnotationSet("s", annotations),),
        [_record()],
        mode="linear",
    )

    defaults = layout_annotation_track(
        "default",
        "s",
        bundle.annotations,
        record_lengths={0: 100},
        params=AnnotationTrackParams(set_id="s", marks=("line", "band")),
        label_font_size=10.0,
    )
    customized = layout_annotation_track(
        "custom",
        "s",
        bundle.annotations,
        record_lengths={0: 100},
        params=AnnotationTrackParams(
            set_id="s",
            marks=("line", "band"),
            lane_gap_px=7,
            padding_px=5,
        ),
        label_font_size=10.0,
    )

    assert {item.annotation.mark for item in defaults.placements} == {
        "line",
        "band",
    }
    assert defaults.lane_gap_px == pytest.approx(3.0)
    assert defaults.required_extent_px == pytest.approx(35.0)
    assert customized.lane_gap_px == pytest.approx(7.0)
    assert customized.required_extent_px == pytest.approx(45.0)


def test_source_coordinates_follow_reverse_coordinate_map() -> None:
    record = _record(10)
    record.annotations["gbdraw_coord_base"] = 100
    record.annotations["gbdraw_coord_step"] = -1
    annotation = RegionAnnotation("a", CoordinateSpan(None, 95, 98))
    resolved = resolve_annotations((AnnotationSet("s", (annotation,)),), [record], mode="linear")
    assert resolved.annotations[0].segments == ((2, 6),)


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_origin_span_uses_complete_record_topology(mode) -> None:
    annotation = RegionAnnotation("a", CoordinateSpan(None, 90, 10, wraps_origin=True))
    record = _record()
    record.annotations["topology"] = "circular"
    circular = resolve_annotations((AnnotationSet("s", (annotation,)),), [record], mode=mode)
    assert circular.annotations[0].segments == ((0, 10), (89, 100))
    for topology in ("linear", "unknown"):
        record.annotations["topology"] = topology
        with pytest.raises(ValidationError, match="complete circular"):
            resolve_annotations((AnnotationSet("s", (annotation,)),), [record], mode=mode)
    record.annotations.update(topology="circular", gbdraw_region_applied=True)
    with pytest.raises(ValidationError, match="complete circular"):
        resolve_annotations((AnnotationSet("s", (annotation,)),), [record], mode=mode)


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


TSV_IMPORT_CASES = json.loads(
    (Path(__file__).parent / "fixtures/annotations/tsv-import-cases.json").read_text()
)


def _literal_rows(text: str) -> list[list[str]]:
    """Split a TSV text the way the table readers do: ``"`` is ordinary and whole-line ``#`` comments are skipped."""
    return [
        line.split("\t")
        for line in text.split("\n")
        if line.strip() and not line.lstrip().startswith("#")
    ]


@pytest.mark.parametrize("case", TSV_IMPORT_CASES, ids=lambda case: case["name"])
def test_annotation_tsv_import_parity(case, tmp_path, caplog) -> None:
    path = tmp_path / "annotations.tsv"
    path.write_text(case["table"], encoding="utf-8")
    rows = _literal_rows(case["table"])
    if case["valid"]:
        controls = _literal_rows(case["control"])
        expected = annotation_sets_from_dataframe(pd.DataFrame(controls[1:], columns=controls[0]))
        assert read_annotation_table(path) == expected
        assert annotation_sets_from_dataframe(pd.DataFrame(rows[1:], columns=rows[0])) == expected
        notices = [record.message for record in caplog.records if "Ignored annotation table columns:" in record.message]
        assert len(notices) == (2 if case["ignored"] else 0)
        for notice in notices:
            assert all(name in notice for name in case["ignored"])
            assert "not saved in Sessions or TSV re-export" in notice
        assert "PRIVATE-CELL" not in caplog.text and "PRIVATE-GENE" not in caplog.text
        assert all(not annotation.metadata for item in expected for annotation in item.annotations)
    else:
        with pytest.raises(ValidationError):
            read_annotation_table(path)
        if all(len(row) == len(rows[0]) for row in rows[1:]):
            with pytest.raises(ValidationError):
                annotation_sets_from_dataframe(pd.DataFrame(rows[1:], columns=rows[0]))
        assert "Ignored annotation table columns:" not in caplog.text
        assert "PRIVATE-CELL" not in caplog.text


def test_annotation_file_keeps_quotes_literally_and_skips_blank_lines(tmp_path) -> None:
    path = tmp_path / "quotes.tsv"
    path.write_text(
        '\n\t\t\nset_id\tid\tmark\tstart\tend\tlabel\n\n'
        's\ta\tband\t1\t8\t"lead\n'
        's\tb\tband\t9\t12\t"quoted"\n'
        's\tc\tband\t13\t15\ta"b\n'
        's\td\tband\t16\t18\t"two words"\n'
    )
    labels = [item.label for item in read_annotation_table(path)[0].annotations]
    assert labels == ['"lead', '"quoted"', 'a"b', '"two words"']


def test_annotation_file_read_failure_has_context(tmp_path) -> None:
    with pytest.raises(ValidationError, match="Unable to read annotation table"):
        read_annotation_table(tmp_path / "missing.tsv")


def test_lane_assignment_is_deterministic_and_separates_overlap() -> None:
    annotations = (
        RegionAnnotation("b", CoordinateSpan(None, 5, 15)),
        RegionAnnotation("a", CoordinateSpan(None, 1, 10)),
        RegionAnnotation("c", CoordinateSpan(None, 20, 25)),
    )
    bundle = resolve_annotations((AnnotationSet("s", annotations),), [_record()], mode="linear")
    lanes = assign_annotation_lanes(bundle.annotations, record_lengths={0: 100})
    assert {item.annotation.id: item.lane for item in lanes} == {"a": 0, "c": 0, "b": 1}


def test_explicit_lane_conflict_names_its_cause_to_the_web() -> None:
    """OV-1262: the Web shows a classified message, not RENDER_FAILED."""
    from gbdraw.web_support.error_adapter import serialize_web_error

    annotations = (
        RegionAnnotation("first", CoordinateSpan(None, 1, 10), lane=0),
        RegionAnnotation("second", CoordinateSpan(None, 5, 15), lane=0),
    )
    bundle = resolve_annotations((AnnotationSet("s", annotations),), [_record()], mode="linear")
    with pytest.raises(ValidationError, match="Annotation s/second conflicts in explicit lane 0.") as raised:
        assign_annotation_lanes(bundle.annotations, record_lengths={0: 100})
    payload = serialize_web_error(raised.value, operation="generate", stage="render")
    assert (payload["code"], payload["stage"], payload["context"]) == (
        "TRACK_INVALID",
        "render",
        {"reason": "ANNOTATION_LANE_CONFLICT"},
    )


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


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("envelope", ["outer_bounds", "segments"])
@pytest.mark.parametrize("circular_path", ["forward", "reverse", "shortest"])
def test_selector_miss_skips_complete_row_with_one_private_warning(
    mode, envelope, circular_path
):
    record = _record()
    record.annotations["topology"] = "circular"
    record.features = [
        SeqFeature(
            FeatureLocation(10, 20), type="CDS", qualifiers={"gene": ["present"]}
        )
    ]
    matched = RegionAnnotation(
        "matched",
        FeatureSpan(
            None, ("gene=present",), envelope=envelope, circular_path=circular_path
        ),
    )
    partial = RegionAnnotation(
        "partial",
        FeatureSpan(
            None,
            (
                "gene=present",
                "gene=PRIVATE-MISSING",
                "gene=PRIVATE-MISSING",
                "gene=OTHER-MISSING",
            ),
            envelope=envelope,
            circular_path=circular_path,
        ),
    )
    absent = RegionAnnotation("absent", FeatureSpan(None, ("gene=PRIVATE-MISSING",)))
    control = resolve_annotations(
        (AnnotationSet("s", (matched,)),), [record], mode=mode
    )
    result = resolve_annotations(
        (AnnotationSet("s", (matched, partial, absent)),), [record], mode=mode
    )
    assert result.annotations == control.annotations
    assert [w.annotation_id for w in result.warnings] == ["partial", "absent"]
    assert [w.missing_count for w in result.warnings] == [2, 1]
    assert all(
        w.code == "feature_selector_unmatched"
        and w.set_id == "s"
        and w.record_id == "r1"
        and w.record_index == 0
        for w in result.warnings
    )
    assert "PRIVATE-MISSING" not in repr(result.warnings)
    assert "OTHER-MISSING" not in repr(result.warnings)


@pytest.mark.parametrize(
    "binding,record_ids,error",
    [
        (None, ["a", "b"], "ambiguous"),
        ("missing", ["a"], "did not match"),
        ("same", ["same", "same"], "multiple records"),
        ("#3", ["a", "b"], "out of range"),
    ],
)
def test_record_errors_remain_fatal_before_missing_feature(binding, record_ids, error):
    annotation = RegionAnnotation(
        "a",
        FeatureSpan(
            parse_record_selector(binding) if binding else None, ("gene=missing",)
        ),
    )
    with pytest.raises(ValidationError, match=error):
        resolve_annotations(
            (AnnotationSet("s", (annotation,)),),
            [_record(record_id=id) for id in record_ids],
            mode="linear",
        )


@pytest.mark.parametrize("selector", ["", "gene=", "=missing"])
def test_malformed_feature_selector_is_fatal(selector):
    with pytest.raises(ValidationError):
        FeatureSpan(None, (selector,))


def test_matched_feature_without_geometry_is_empty_span():
    record = _record()
    record.features = [
        SeqFeature(
            FeatureLocation(10, 10), type="CDS", qualifiers={"gene": ["present"]}
        )
    ]
    row = RegionAnnotation("a", FeatureSpan(None, ("gene=present",)))
    result = resolve_annotations((AnnotationSet("s", (row,)),), [record], mode="linear")
    assert not result.annotations
    assert [w.code for w in result.warnings] == ["empty_span"]
    assert result.warnings[0].missing_count == 0


@pytest.mark.parametrize(
    "policy,code",
    [
        ("clip", "out_of_bounds_clipped"),
        ("skip", "out_of_bounds_skipped"),
        ("error", None),
    ],
)
def test_coordinate_policy_remains_independent_of_selector_misses(policy, code):
    rows = (
        RegionAnnotation(
            "coordinate", CoordinateSpan(None, 90, 120, out_of_bounds=policy)
        ),
        RegionAnnotation("missing", FeatureSpan(None, ("gene=missing",))),
    )
    if policy == "error":
        with pytest.raises(ValidationError, match="extends outside"):
            resolve_annotations((AnnotationSet("s", rows),), [_record()], mode="linear")
    else:
        result = resolve_annotations(
            (AnnotationSet("s", rows),), [_record()], mode="linear"
        )
        assert [w.code for w in result.warnings] == [code, "feature_selector_unmatched"]
        assert [w.record_id for w in result.warnings] == ["r1", "r1"]
        assert [item.segments for item in result.annotations] == (
            [((89, 100),)] if policy == "clip" else []
        )


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("envelope", ["outer_bounds", "segments"])
@pytest.mark.parametrize("path", ["forward", "reverse", "shortest"])
def test_fully_matched_envelope_and_circular_path_are_preserved(mode, envelope, path):
    record = _record()
    record.annotations["topology"] = "circular"
    record.features = [
        SeqFeature(FeatureLocation(start, end), type="CDS", qualifiers={"gene": [gene]})
        for start, end, gene in [(10, 20, "first"), (80, 90, "last")]
    ]
    row = RegionAnnotation(
        "a",
        FeatureSpan(
            None, ("gene=first", "gene=last"), envelope=envelope, circular_path=path
        ),
    )
    result = resolve_annotations((AnnotationSet("s", (row,)),), [record], mode=mode)
    expected = (
        ((10, 20), (80, 90))
        if envelope == "segments"
        else (((80, 100), (0, 20)) if path != "forward" else ((10, 90),))
    )
    assert result.annotations[0].segments == expected
    assert result.warnings == ()


@pytest.mark.parametrize(
    ("field", "name"),
    [("fill", "annotation fill"), ("stroke", "annotation stroke"), ("label_color", "label color")],
)
def test_malformed_hex_annotation_color_names_its_column(field: str, name: str) -> None:
    # The hex check of resolve_color_to_hex raises its own ValidationError; the
    # annotation message keeps the column and the annotation-correct advice.
    with pytest.raises(ValidationError, match=f"{name} is not a valid color: '#zz'"):
        RegionAnnotationStyle(**{field: "#zz"})


def _drawn_label_fonts(
    mode: str,
    length: int,
    *,
    style: RegionAnnotationStyle | None = None,
    config_overrides: dict[str, object] | None = None,
) -> tuple[float, float]:
    """The drawn font size of a region label and of a feature label in one diagram of ``mode``."""

    from gbdraw.api import CircularDiagramOptions, LinearDiagramOptions
    from gbdraw.api.diagram import build_circular_diagram, build_linear_diagram

    record = SeqRecord(Seq("ACGT" * (length // 4)), id="r1", name="r1", annotations={"molecule_type": "DNA"})
    record.annotations["topology"] = "circular"
    record.features = [
        SeqFeature(FeatureLocation(100, 400, strand=1), type="CDS", qualifiers={"product": ["alpha protein"]})
    ]
    annotations = AnnotationOptions(
        sets=(AnnotationSet("regions", (RegionAnnotation("reg", CoordinateSpan(None, 500, 700), label="RegionX", style=style),)),)
    )
    overrides = {f"labels.{mode}.scope": "all" if mode == "linear" else "both", **(config_overrides or {})}
    if mode == "linear":
        drawing = build_linear_diagram(
            [record], options=LinearDiagramOptions(annotations=annotations, config_overrides=overrides)
        )
    else:
        drawing = build_circular_diagram(
            record, options=CircularDiagramOptions(annotations=annotations, config_overrides=overrides)
        )
    sizes = {
        text: float(size)
        for size, text in re.findall(
            r'<text[^>]*font-size="([^"]+)"[^>]*>(?:<textPath[^>]*>)?([^<]*)<', drawing.tostring()
        )
    }
    return sizes["RegionX"], sizes["alpha protein"]


@pytest.mark.parametrize("length", [1000, 200000], ids=["short", "long"])
@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_region_label_without_a_font_takes_the_feature_label_font(mode: str, length: int) -> None:
    region, feature = _drawn_label_fonts(mode, length)
    assert region == feature


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_region_label_follows_a_feature_label_font_override(mode: str) -> None:
    path = "labels.font_size.linear" if mode == "linear" else "labels.font_size"
    region, feature = _drawn_label_fonts(mode, 1000, config_overrides={f"{path}.short": 11.0})
    assert (region, feature) == (11.0, 11.0)


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_region_label_with_its_own_font_keeps_it(mode: str) -> None:
    region, feature = _drawn_label_fonts(mode, 1000, style=RegionAnnotationStyle(label_font_size=7.0))
    assert region == 7.0
    assert feature != 7.0


def test_lane_packing_measures_a_label_in_the_slot_style_font() -> None:
    annotations = (
        RegionAnnotation("first", CoordinateSpan(None, 40, 41), label="Wide label"),
        RegionAnnotation("second", CoordinateSpan(None, 60, 61), label="Next"),
    )
    bundle = resolve_annotations((AnnotationSet("s", annotations),), [_record()], mode="linear")

    def lanes(style_override: RegionAnnotationStyle | None) -> set[int]:
        track = layout_annotation_track(
            "regions",
            "s",
            bundle.annotations,
            record_lengths={0: 100},
            params=AnnotationTrackParams(set_id="s", style_override=style_override),
            label_font_size=1.0,
        )
        return {placement.lane for placement in track.placements}

    assert lanes(None) == {0}
    assert lanes(RegionAnnotationStyle(label_font_size=20.0)) == {0, 1}


def test_lane_packing_pads_a_bracket_by_the_slot_style_stroke() -> None:
    annotations = (
        RegionAnnotation("first", CoordinateSpan(None, 10, 40), mark="bracket"),
        RegionAnnotation("second", CoordinateSpan(None, 50, 80), mark="bracket"),
    )
    bundle = resolve_annotations((AnnotationSet("s", annotations),), [_record()], mode="linear")

    def lanes(style_override: RegionAnnotationStyle | None) -> set[int]:
        track = layout_annotation_track(
            "regions",
            "s",
            bundle.annotations,
            record_lengths={0: 100},
            params=AnnotationTrackParams(set_id="s", padding_px=0.0, style_override=style_override),
            label_font_size=10.0,
        )
        return {placement.lane for placement in track.placements}

    assert lanes(None) == {0}
    assert lanes(RegionAnnotationStyle(stroke_width=3.0)) == {0, 1}
