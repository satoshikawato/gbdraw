from __future__ import annotations

import copy
import re
import xml.etree.ElementTree as ET
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing
from svgwrite.container import Group

from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import (
    assemble_circular_diagram_from_record,
    assemble_circular_diagram_from_records,
    assemble_linear_diagram_from_records,
)
from gbdraw.annotations import AnnotationSet, CoordinateSpan, RegionAnnotation
from gbdraw.api.options import AnnotationOptions
from gbdraw.canvas import LinearCanvasConfigurator
from gbdraw.config.models import GbdrawConfig, LinearRenderProfile
from gbdraw.config.toml import load_config_toml
from gbdraw.configurators import LegendDrawingConfigurator
from gbdraw.features.ids import compute_feature_hash
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.linear_comparison import LinearComparison
from gbdraw.render.drawers.circular.gc_skew import SkewDrawer as CircularSkewDrawer
from gbdraw.render.drawers.linear.gc_skew import SkewDrawer as LinearSkewDrawer
from gbdraw.render.groups.linear.legend import LegendGroup
from gbdraw.render.interactive_svg import enrich_svg
from gbdraw.svg.ids import definition_group_svg_id, instance_svg_id
from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot

_URL_REFERENCE_RE = re.compile(r"url\(#([^)]+)\)")


def _root(svg: str) -> ET.Element:
    return ET.fromstring(svg)


def _ids(svg: str) -> list[str]:
    return [
        element.attrib["id"]
        for element in _root(svg).iter()
        if element.attrib.get("id")
    ]


def _assert_unique_ids_and_resolved_references(svg: str) -> None:
    root = _root(svg)
    ids = [
        element.attrib["id"]
        for element in root.iter()
        if element.attrib.get("id")
    ]
    assert len(ids) == len(set(ids))

    known_ids = set(ids)
    references: set[str] = set()
    for element in root.iter():
        for name, value in element.attrib.items():
            references.update(_URL_REFERENCE_RE.findall(value))
            if name.rsplit("}", 1)[-1] == "href" and value.startswith("#"):
                references.add(value[1:])
    assert references <= known_ids


def _skew_config() -> SimpleNamespace:
    return SimpleNamespace(
        high_fill_color="#ef4444",
        low_fill_color="#3b82f6",
        stroke_color="none",
        stroke_width=0.0,
        fill_opacity=1.0,
    )


def _skew_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {"GC skew": [0.25, -0.25, 0.1]},
        index=[0, 40, 80],
    )


def _duplicate_feature_record(record_id: str = "duplicate-feature") -> SeqRecord:
    record = SeqRecord(Seq("ATGC" * 80), id=record_id, name=record_id)
    feature = SeqFeature(
        SimpleLocation(30, 120, strand=1),
        type="CDS",
        qualifiers={"product": ["same"], "translation": ["M" * 30]},
    )
    record.features = [feature, copy.deepcopy(feature)]
    return record


def _linear_legend_svg() -> str:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    cfg = GbdrawConfig.from_dict(config_dict)
    profile = LinearRenderProfile(cfg)
    canvas_config = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=1_000,
        profile=profile,
        legend="bottom",
    )
    legend_table = {
        "CDS": {
            "type": "solid",
            "fill": "#54bcf8",
            "stroke": "none",
            "width": 0,
        },
        "Pairwise match identity": {
            "type": "gradient",
            "min_color": "#ffffff",
            "max_color": "#000000",
            "stroke": "none",
            "width": 0,
            "min_value": 0,
        },
    }
    configurator = LegendDrawingConfigurator(
        color_table=None,
        default_colors=None,
        selected_features_set=[],
        profile=profile,
        gc_config=None,
        skew_config=None,
        feature_config=None,
        canvas_config=canvas_config,
    )
    legend_measurement = configurator.measure_legend(
        legend_table,
        placement=canvas_config.legend_position,
        wrap_width=canvas_config.total_width,
    )
    drawing = Drawing(debug=False)
    drawing.add(
        LegendGroup(
            canvas_config,
            legend_measurement,
            legend_table,
            cfg=cfg,
        ).get_group()
    )
    return drawing.tostring()


def test_instance_svg_id_distinguishes_values_with_the_same_safe_fragment() -> None:
    spaced = instance_svg_id("feature", "a b")
    underscored = instance_svg_id("feature", "a_b")

    assert spaced != underscored
    assert spaced == instance_svg_id("feature", "a b")
    assert spaced.startswith("feature__instance_a_b_")
    assert underscored.startswith("feature__instance_a_b_")


def test_definition_group_svg_id_is_valid_collision_safe_and_deterministic() -> None:
    valid = definition_group_svg_id("safe.id-1", mode="linear")
    leading_digit = definition_group_svg_id("123", mode="linear")
    spaced = definition_group_svg_id("a b", mode="linear")
    slashed = definition_group_svg_id("a/b", mode="linear")
    underscored = definition_group_svg_id("a_b", mode="linear")

    assert valid == "safe.id-1_definition"
    assert len({spaced, slashed, underscored}) == 3
    assert leading_digit.startswith("record_123_")
    assert spaced == definition_group_svg_id("a b", mode="linear")
    assert spaced != definition_group_svg_id("a b", mode="circular")
    assert definition_group_svg_id(
        "same",
        mode="linear",
        record_index=0,
        record_count=2,
    ) != definition_group_svg_id(
        "same",
        mode="linear",
        record_index=1,
        record_count=2,
    )
    assert all(
        re.fullmatch(r"[A-Za-z_][A-Za-z0-9_.-]*", identifier)
        for identifier in (valid, leading_digit, spaced, slashed, underscored)
    )


def _record_definition_groups(svg: str) -> list[ET.Element]:
    return [
        element
        for element in _root(svg).iter()
        if element.tag.rsplit("}", 1)[-1] == "g"
        and element.attrib.get("data-gbdraw-role") == "record-definition"
    ]


@pytest.mark.parametrize(
    "mode",
    [
        pytest.param("linear", marks=pytest.mark.linear),
        pytest.param("circular", marks=pytest.mark.circular),
    ],
)
def test_single_record_definition_id_and_hooks_handle_unsafe_leading_id(
    mode: str,
) -> None:
    record = SeqRecord(Seq("ATGC" * 40), id="123/unsafe id")
    label_path = (
        "labels.circular.scope"
        if mode == "circular"
        else "labels.linear.scope"
    )
    label_value = "none"
    common = {
        "legend": "none",
        "cfg": apply_config_overrides(None, {
            label_path: label_value,
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
    }

    def render() -> str:
        canvas = (
            assemble_circular_diagram_from_record(record, **common)
            if mode == "circular"
            else assemble_linear_diagram_from_records([record], **common)
        )
        return canvas.tostring()

    first = render()
    second = render()
    assert _ids(first) == _ids(second)
    _assert_unique_ids_and_resolved_references(first)

    definitions = _record_definition_groups(first)
    assert len(definitions) == 1
    definition = definitions[0]
    assert definition.attrib["id"] == definition_group_svg_id(record.id, mode=mode)
    assert definition.attrib["data-gbdraw-definition-part"] == "main"
    assert definition.attrib["data-gbdraw-record-id"] == record.id
    assert definition.attrib["data-gbdraw-record-index"] == "0"
    assert re.fullmatch(r"[A-Za-z_][A-Za-z0-9_.-]*", definition.attrib["id"])


@pytest.mark.parametrize(
    "mode",
    [
        pytest.param("linear", marks=pytest.mark.linear),
        pytest.param("circular", marks=pytest.mark.circular),
    ],
)
def test_multi_record_definition_ids_disambiguate_duplicates_and_safe_aliases(
    mode: str,
) -> None:
    raw_record_ids = ["123/unsafe id", "123/unsafe id", "a b", "a_b"]
    records = [
        SeqRecord(Seq("ATGC" * 40), id=record_id)
        for record_id in raw_record_ids
    ]
    label_path = (
        "labels.circular.scope"
        if mode == "circular"
        else "labels.linear.scope"
    )
    label_value = "none"
    common = {
        "legend": "none",
        "cfg": apply_config_overrides(None, {
            label_path: label_value,
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
    }

    def render() -> str:
        canvas = (
            assemble_circular_diagram_from_records(records, **common)
            if mode == "circular"
            else assemble_linear_diagram_from_records(records, **common)
        )
        return canvas.tostring()

    first = render()
    second = render()
    assert _ids(first) == _ids(second)
    _assert_unique_ids_and_resolved_references(first)

    definitions = _record_definition_groups(first)
    assert len(definitions) == len(records)
    definitions_by_index = {
        int(definition.attrib["data-gbdraw-record-index"]): definition
        for definition in definitions
    }
    assert set(definitions_by_index) == set(range(len(records)))
    assert [
        definitions_by_index[index].attrib["data-gbdraw-record-id"]
        for index in range(len(records))
    ] == raw_record_ids

    definition_ids = [
        definitions_by_index[index].attrib["id"]
        for index in range(len(records))
    ]
    assert len(definition_ids) == len(set(definition_ids))
    assert all(
        re.fullmatch(r"[A-Za-z_][A-Za-z0-9_.-]*", identifier)
        for identifier in definition_ids
    )

    expected_ids = []
    for index, record_id in enumerate(raw_record_ids):
        record_count = (
            len(records)
            if mode == "linear" or raw_record_ids.count(record_id) > 1
            else 1
        )
        expected_ids.append(
            definition_group_svg_id(
                record_id,
                mode=mode,
                record_index=index,
                record_count=record_count,
            )
        )
    assert definition_ids == expected_ids


@pytest.mark.linear
def test_linear_skew_clip_id_is_deterministic_unique_and_referenced() -> None:
    def render() -> str:
        drawing = Drawing(debug=False)
        group = Group(id="gc_skew_record_1")
        drawing.add(
            LinearSkewDrawer(_skew_config()).draw(
                group,
                _skew_frame(),
                120,
                600.0,
                1.0,
                24.0,
                0.0,
                0.0,
                "GC",
                "gc_skew_record_1",
            )
        )
        return drawing.tostring()

    first = render()
    second = render()
    assert _ids(first) == _ids(second)
    _assert_unique_ids_and_resolved_references(first)


@pytest.mark.circular
def test_repeated_circular_skew_slots_get_distinct_deterministic_clip_ids() -> None:
    drawing = Drawing(debug=False)
    for group_id in ("gc_skew", "gc_skew_2"):
        group = Group(id=group_id)
        drawing.add(
            CircularSkewDrawer(_skew_config()).draw(
                100.0,
                group,
                _skew_frame(),
                120,
                20.0,
                0.8,
                "GC",
                record_identifier="record-a",
                group_identifier=group_id,
            )
        )

    svg = drawing.tostring()
    _assert_unique_ids_and_resolved_references(svg)
    clip_ids = [
        element.attrib["id"]
        for element in _root(svg).iter()
        if element.tag.rsplit("}", 1)[-1] == "clipPath"
    ]
    assert len(clip_ids) == 2
    assert len(set(clip_ids)) == 2


@pytest.mark.linear
def test_linear_dual_legend_ids_are_namespaced_deterministic_and_referenced() -> None:
    first = _linear_legend_svg()
    second = _linear_legend_svg()
    assert _ids(first) == _ids(second)
    _assert_unique_ids_and_resolved_references(first)

    root = _root(first)
    pairwise_groups = [
        element
        for element in root.iter()
        if element.attrib.get("data-gbdraw-role") == "comparison-legend"
    ]
    assert {group.attrib["id"] for group in pairwise_groups} == {
        "pairwise_legend_h",
        "pairwise_legend_v",
    }
    assert {group.attrib["data-gbdraw-orientation"] for group in pairwise_groups} == {
        "h",
        "v",
    }


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_identical_semantic_features_keep_stable_hooks_and_unique_dom_ids(
    mode: str,
) -> None:
    record = _duplicate_feature_record()
    stable_id = compute_feature_hash(record.features[0], record_id=record.id)
    label_path = (
        "labels.circular.scope"
        if mode == "circular"
        else "labels.linear.scope"
    )
    label_value = "none"
    common = {
        "legend": "none",
        "cfg": apply_config_overrides(None, {
            label_path: label_value,
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
    }
    canvas = (
        assemble_circular_diagram_from_record(record, **common)
        if mode == "circular"
        else assemble_linear_diagram_from_records([record], **common)
    )
    svg = canvas.tostring()
    _assert_unique_ids_and_resolved_references(svg)

    matching = [
        element
        for element in _root(svg).iter()
        if element.attrib.get("data-gbdraw-feature-id") == stable_id
        and element.attrib.get("data-gbdraw-feature-part") == "block"
    ]
    assert len(matching) == 2
    assert len({element.attrib["id"] for element in matching}) == 2
    assert all("__instance_" in element.attrib["id"] for element in matching)
    assert {
        element.attrib["data-gbdraw-source-feature-index"]
        for element in matching
    } == {"0", "1"}


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_identical_underlay_features_keep_stable_hooks_and_unique_dom_ids(
    mode: str,
) -> None:
    record = SeqRecord(Seq("ATGC" * 80), id="duplicate-underlay")
    feature = SeqFeature(
        SimpleLocation(30, 120, strand=1),
        type="repeat_region",
        qualifiers={"rpt_family": ["same"]},
    )
    record.features = [feature, copy.deepcopy(feature)]
    stable_id = compute_feature_hash(feature, record_id=record.id)
    label_path = (
        "labels.circular.scope"
        if mode == "circular"
        else "labels.linear.scope"
    )
    label_value = "none"
    common = {
        "selected_features_set": ["repeat_region"],
        "legend": "none",
        "cfg": apply_config_overrides(None, {label_path: label_value}),
    }
    canvas = (
        assemble_circular_diagram_from_record(record, **common)
        if mode == "circular"
        else assemble_linear_diagram_from_records([record], **common)
    )
    svg = canvas.tostring()
    _assert_unique_ids_and_resolved_references(svg)

    matching = [
        element
        for element in _root(svg).iter()
        if element.attrib.get("data-gbdraw-feature-id") == stable_id
        and element.attrib.get("data-gbdraw-auto-feature-underlay") == "true"
    ]
    assert len(matching) == 2
    assert len({element.attrib["id"] for element in matching}) == 2
    assert all("__instance_" in element.attrib["id"] for element in matching)


@pytest.mark.circular
def test_repeated_circular_record_grid_namespaces_colliding_ids_and_references() -> None:
    record = _duplicate_feature_record("same-record")
    canvas = assemble_circular_diagram_from_records(
        [record, copy.deepcopy(record)],
        legend="none",
        cfg=apply_config_overrides(None, {
            "labels.circular.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": True,
        }),
        window=40,
        step=40,
    )
    svg = canvas.tostring()
    _assert_unique_ids_and_resolved_references(svg)
    assert "__record_2" in svg


@pytest.mark.circular
def test_circular_multi_record_wrappers_expose_record_identity_hooks() -> None:
    records = [
        _duplicate_feature_record("wrapper-a"),
        _duplicate_feature_record("wrapper-b"),
    ]
    svg = assemble_circular_diagram_from_records(
        records,
        legend="none",
        cfg=apply_config_overrides(None, {
            "labels.circular.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
    ).tostring()

    wrappers = {
        element.attrib["data-gbdraw-record-index"]: element
        for element in list(_root(svg))
        if element.tag.rsplit("}", 1)[-1] == "g"
        and element.attrib.get("data-gbdraw-record-id")
    }
    assert set(wrappers) == {"0", "1"}
    assert wrappers["0"].attrib["data-gbdraw-record-id"] == "wrapper-a"
    assert wrappers["1"].attrib["data-gbdraw-record-id"] == "wrapper-b"
    assert wrappers["0"].attrib["id"] == "record_0"
    assert wrappers["1"].attrib["id"] == "record_1"


@pytest.mark.parametrize(
    "mode",
    [
        pytest.param("circular", marks=pytest.mark.circular),
        pytest.param("linear", marks=pytest.mark.linear),
    ],
)
def test_default_dinucleotide_tracks_expose_canonical_slot_hooks(
    mode: str,
) -> None:
    record = _duplicate_feature_record(f"default-{mode}-slots")
    label_path = (
        "labels.circular.scope"
        if mode == "circular"
        else "labels.linear.scope"
    )
    label_value = "none"
    common = {
        "legend": "none",
        "cfg": apply_config_overrides(None, {
            label_path: label_value,
            "canvas.show_gc": True,
            "canvas.show_skew": True,
        }),
        "window": 40,
        "step": 40,
    }
    canvas = (
        assemble_circular_diagram_from_record(record, **common)
        if mode == "circular"
        else assemble_linear_diagram_from_records([record], **common)
    )

    groups = {
        element.attrib["data-gbdraw-slot-id"]: element
        for element in _root(canvas.tostring()).iter()
        if element.attrib.get("data-gbdraw-slot-renderer")
        in {"dinucleotide_content", "dinucleotide_skew"}
    }
    assert {
        slot_id: group.attrib["data-gbdraw-slot-renderer"]
        for slot_id, group in groups.items()
    } == {
        "gc_content": "dinucleotide_content",
        "gc_skew": "dinucleotide_skew",
    }
    assert groups["gc_content"].attrib["id"] == "gc_content"
    assert groups["gc_skew"].attrib["id"] in {"skew", "gc_skew"}


@pytest.mark.circular
def test_default_circular_ticks_expose_slot_hooks_without_renaming_group() -> None:
    record = _duplicate_feature_record("default-circular-ticks")
    root = _root(
        assemble_circular_diagram_from_record(
            record,
            legend="none",
            cfg=apply_config_overrides(None, {
                "labels.circular.scope": "none",
                "canvas.show_gc": False,
                "canvas.show_skew": False,
            }),
        ).tostring()
    )

    tick_groups = [
        element
        for element in root.iter()
        if element.attrib.get("data-gbdraw-slot-renderer") == "ticks"
    ]
    assert len(tick_groups) == 1
    assert tick_groups[0].attrib["data-gbdraw-slot-id"] == "ticks"
    assert tick_groups[0].attrib["id"] == "tick"


@pytest.mark.circular
def test_circular_multi_reserved_record_ids_do_not_collide_with_shared_groups() -> None:
    records = [
        _duplicate_feature_record("legend"),
        _duplicate_feature_record("plot_title"),
    ]
    common = {
        "legend": "right",
        "plot_title": "Shared title",
        "plot_title_position": "top",
        "cfg": apply_config_overrides(None, {
            "labels.circular.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
    }

    first = assemble_circular_diagram_from_records(records, **common).tostring()
    second = assemble_circular_diagram_from_records(records, **common).tostring()
    assert _ids(first) == _ids(second)
    _assert_unique_ids_and_resolved_references(first)

    root = _root(first)
    assert len(
        [
            element
            for element in root.iter()
            if element.attrib.get("id") == "legend"
        ]
    ) == 1
    assert len(
        [
            element
            for element in root.iter()
            if element.attrib.get("id") == "plot_title"
        ]
    ) == 1
    assert {
        element.attrib.get("data-gbdraw-record-id")
        for element in root.iter()
        if element.attrib.get("data-gbdraw-record-id")
    } >= {"legend", "plot_title"}


@pytest.mark.parametrize(
    ("mode", "record_id"),
    [
        pytest.param("linear", "legend", marks=pytest.mark.linear),
        pytest.param("linear", "ordinary-record", marks=pytest.mark.linear),
        pytest.param("circular", "Axis", marks=pytest.mark.circular),
        pytest.param("circular", "tick", marks=pytest.mark.circular),
        pytest.param("circular", "ordinary-record", marks=pytest.mark.circular),
    ],
)
def test_record_group_ids_do_not_collide_with_renderer_owned_groups(
    mode: str,
    record_id: str,
) -> None:
    record = _duplicate_feature_record(record_id)
    label_path = (
        "labels.circular.scope"
        if mode == "circular"
        else "labels.linear.scope"
    )
    label_value = "none"
    common = {
        "legend": "right",
        "cfg": apply_config_overrides(None, {
            label_path: label_value,
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
    }

    def render() -> str:
        canvas = (
            assemble_circular_diagram_from_record(record, **common)
            if mode == "circular"
            else assemble_linear_diagram_from_records([record], **common)
        )
        return canvas.tostring()

    first = render()
    second = render()
    assert _ids(first) == _ids(second)
    _assert_unique_ids_and_resolved_references(first)

    record_groups = [
        element
        for element in _root(first).iter()
        if element.tag.rsplit("}", 1)[-1] == "g"
        and element.attrib.get("data-gbdraw-record-id") == record_id
        and element.attrib.get("id", "").startswith("record_group_")
    ]
    assert len(record_groups) == 1
    assert record_groups[0].attrib["id"] != record_id
    assert record_groups[0].attrib["id"].startswith("record_group_")
    assert record_groups[0].attrib["data-gbdraw-record-index"] == "0"


@pytest.mark.linear
def test_linear_record_id_does_not_collide_with_comparison_group() -> None:
    records = [
        _duplicate_feature_record("comparison1"),
        _duplicate_feature_record("comparison-subject"),
    ]
    row = ["q", "s", 90.0, 100, 0, 0, 10, 100, 20, 110, 1e-20, 200]
    comparison = LinearComparison(
        0,
        1,
        pd.DataFrame([row], columns=COMPARISON_COLUMNS),
    )

    canvas = assemble_linear_diagram_from_records(
        records,
        linear_comparisons=[comparison],
        legend="none",
        cfg=apply_config_overrides(None, {
            "labels.linear.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
    )
    svg = canvas.tostring()
    _assert_unique_ids_and_resolved_references(svg)

    root = _root(svg)
    comparison_groups = [
        element
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "g"
        and element.attrib.get("id") == "comparison1"
    ]
    record_group = next(
        element
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "g"
        and element.attrib.get("data-gbdraw-record-id") == "comparison1"
        and element.attrib.get("id", "").startswith("record_group_")
    )
    assert len(comparison_groups) == 1
    assert record_group.attrib["id"].startswith("record_group_")
    assert record_group.attrib["id"] != "comparison1"


@pytest.mark.linear
def test_annotation_track_prefix_record_id_is_semantic_only() -> None:
    record_id = "gbdraw-annotation-track-regions-1"
    record = _duplicate_feature_record(record_id)
    annotation = RegionAnnotation(
        "reviewed-region",
        CoordinateSpan(None, 10, 80),
        mark="band",
    )
    canvas = assemble_linear_diagram_from_records(
        [record],
        annotation_options=AnnotationOptions(
            sets=(AnnotationSet("regions", (annotation,)),),
        ),
        linear_track_slots=("regions:annotations@set_id=regions,h=20px",),
        legend="none",
        cfg=apply_config_overrides(None, {
            "labels.linear.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
    )
    svg = canvas.tostring()
    _assert_unique_ids_and_resolved_references(svg)

    root = _root(svg)
    annotation_group = next(
        element
        for element in root.iter()
        if element.attrib.get("data-gbdraw-slot-id") == "regions"
    )
    record_group = next(
        element
        for element in root.iter()
        if element.attrib.get("data-gbdraw-record-id") == record_id
        and element.tag.rsplit("}", 1)[-1] == "g"
        and element.attrib.get("id", "").startswith("record_group_")
        and element.attrib.get("id", "").startswith("record_group_")
    )
    assert annotation_group.attrib["id"].startswith("track_slot_")
    assert record_group.attrib["id"].startswith("record_group_")
    assert record_group.attrib["id"] != annotation_group.attrib["id"]


@pytest.mark.linear
@pytest.mark.parametrize(
    "record_id",
    [
        "gbdraw-interactive-feature-metadata",
        "gbdraw-interactive-feature-style",
        "gbdraw-interactive-feature-script",
        "gbdraw-interactive-feature-glow",
        "gbdraw-interactive-feature-match-glow",
        "gbdraw-viewport-controls",
    ],
)
def test_interactive_asset_ids_do_not_replace_record_groups(record_id: str) -> None:
    canvas = assemble_linear_diagram_from_records(
        [_duplicate_feature_record(record_id)],
        legend="none",
        cfg=apply_config_overrides(None, {
            "labels.linear.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
    )
    static_svg = canvas.tostring()
    static_record_group = next(
        element
        for element in _root(static_svg).iter()
        if element.attrib.get("data-gbdraw-record-id") == record_id
        and element.tag.rsplit("}", 1)[-1] == "g"
        and element.attrib.get("id", "").startswith("record_group_")
    )

    interactive_svg = enrich_svg(static_svg)
    _assert_unique_ids_and_resolved_references(interactive_svg)
    interactive_record_groups = [
        element
        for element in _root(interactive_svg).iter()
        if element.attrib.get("data-gbdraw-record-id") == record_id
        and element.tag.rsplit("}", 1)[-1] == "g"
        and element.attrib.get("id", "").startswith("record_group_")
    ]
    assert len(interactive_record_groups) == 1
    assert interactive_record_groups[0].attrib["id"] == static_record_group.attrib["id"]
    assert interactive_record_groups[0].attrib["id"].startswith("record_group_")


@pytest.mark.linear
def test_linear_user_slot_ids_are_namespaced_without_safe_fragment_aliases() -> None:
    record = _duplicate_feature_record("linear-slot-record")
    common = {
        "legend": "right",
        "cfg": apply_config_overrides(None, {
            "labels.linear.scope": "none",
            "canvas.show_gc": True,
            "canvas.show_skew": True,
        }),
        "linear_track_slots": [
            LinearTrackSlot("legend", "features", side="overlay"),
            LinearTrackSlot(
                "a b",
                "dinucleotide_skew",
                side="above",
                params={"nt": "GC"},
            ),
            LinearTrackSlot(
                "a_b",
                "dinucleotide_skew",
                side="below",
                params={"nt": "AT"},
            ),
        ],
        "window": 40,
        "step": 40,
    }

    first = assemble_linear_diagram_from_records([record], **common).tostring()
    second = assemble_linear_diagram_from_records([record], **common).tostring()
    assert _ids(first) == _ids(second)
    _assert_unique_ids_and_resolved_references(first)

    root = _root(first)
    slotted_groups = {
        element.attrib["data-gbdraw-slot-id"]: element
        for element in root.iter()
        if element.attrib.get("data-gbdraw-slot-id")
    }
    assert set(slotted_groups) == {"legend", "a b", "a_b"}
    assert all(
        group.attrib["id"].startswith("track_slot_")
        for group in slotted_groups.values()
    )
    assert slotted_groups["a b"].attrib["id"] != slotted_groups["a_b"].attrib["id"]
    assert slotted_groups["legend"].attrib["id"] != "legend"
    assert {
        slot_id: group.attrib["data-gbdraw-slot-renderer"]
        for slot_id, group in slotted_groups.items()
    } == {
        "legend": "features",
        "a b": "dinucleotide_skew",
        "a_b": "dinucleotide_skew",
    }
    assert sum(
        element.attrib.get("id") == "legend"
        for element in root.iter()
    ) == 1


@pytest.mark.linear
def test_linear_depth_alias_slots_namespace_outer_and_axis_ids() -> None:
    record = _duplicate_feature_record("linear-depth-slot-record")
    depth_table = pd.DataFrame(
        {
            "reference_name": [record.id, record.id, record.id],
            "position": [1, 120, 240],
            "depth": [5.0, 20.0, 10.0],
        }
    )
    canvas = assemble_linear_diagram_from_records(
        [record],
        depth_table=depth_table,
        linear_track_slots=[
            LinearTrackSlot(
                "a b",
                "depth",
                side="above",
                params={"track_index": 0},
            ),
            LinearTrackSlot(
                "a_b",
                "depth",
                side="below",
                params={"track_index": 0},
            ),
        ],
        legend="none",
        cfg=apply_config_overrides(None, {
            "labels.linear.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
        depth_window=40,
        depth_step=40,
    )
    svg = canvas.tostring()
    _assert_unique_ids_and_resolved_references(svg)

    slotted_groups = {
        element.attrib["data-gbdraw-slot-id"]: element
        for element in _root(svg).iter()
        if element.attrib.get("data-gbdraw-slot-id")
    }
    assert set(slotted_groups) == {"a b", "a_b"}
    assert len({group.attrib["id"] for group in slotted_groups.values()}) == 2
    assert all(
        group.attrib["id"].startswith("track_slot_")
        for group in slotted_groups.values()
    )
    assert {
        group.attrib["data-gbdraw-slot-renderer"]
        for group in slotted_groups.values()
    } == {"depth"}


@pytest.mark.linear
def test_linear_user_slot_ids_are_record_scoped_in_multi_record_svg() -> None:
    records = [
        _duplicate_feature_record("slot-record-a"),
        _duplicate_feature_record("slot-record-b"),
    ]
    canvas = assemble_linear_diagram_from_records(
        records,
        linear_track_slots=[
            LinearTrackSlot("features", "features", side="overlay"),
            LinearTrackSlot(
                "a b",
                "dinucleotide_skew",
                side="above",
                params={"nt": "GC"},
            ),
            LinearTrackSlot(
                "a_b",
                "dinucleotide_skew",
                side="below",
                params={"nt": "AT"},
            ),
        ],
        legend="none",
        cfg=apply_config_overrides(None, {
            "labels.linear.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": True,
        }),
        window=40,
        step=40,
    )
    svg = canvas.tostring()
    _assert_unique_ids_and_resolved_references(svg)

    slotted_groups = [
        element
        for element in _root(svg).iter()
        if element.attrib.get("data-gbdraw-slot-id")
    ]
    assert len(slotted_groups) == 6
    assert len({element.attrib["id"] for element in slotted_groups}) == 6
    assert {
        element.attrib["data-gbdraw-slot-id"]
        for element in slotted_groups
    } == {"features", "a b", "a_b"}


@pytest.mark.parametrize("mode", ["linear", "circular"])
def test_annotation_slot_ids_are_namespaced_without_safe_fragment_aliases(
    mode: str,
) -> None:
    record = _duplicate_feature_record(f"{mode}-annotation-slot-record")
    annotation_options = AnnotationOptions(
        sets=(
            AnnotationSet(
                "regions",
                (
                    RegionAnnotation(
                        "reviewed-region",
                        CoordinateSpan(None, 10, 80),
                        mark="band",
                    ),
                ),
            ),
        ),
    )
    if mode == "linear":
        canvas = assemble_linear_diagram_from_records(
            [record],
            cfg=apply_config_overrides(None, None),
            annotation_options=annotation_options,
            linear_track_slots=[
                LinearTrackSlot(
                    "a b",
                    "annotations",
                    side="above",
                    params={"set_id": "regions", "overflow": "clip"},
                ),
                LinearTrackSlot(
                    "a-b",
                    "annotations",
                    side="below",
                    params={"set_id": "regions", "overflow": "clip"},
                ),
            ],
            legend="none",
        )
    else:
        canvas = assemble_circular_diagram_from_record(
            record,
            cfg=apply_config_overrides(None, None),
            annotation_options=annotation_options,
            circular_track_slots=[
                CircularTrackSlot(
                    "a b",
                    "annotations",
                    side="outside",
                    params={"set_id": "regions", "overflow": "clip"},
                ),
                CircularTrackSlot(
                    "a-b",
                    "annotations",
                    side="inside",
                    params={"set_id": "regions", "overflow": "clip"},
                ),
            ],
            legend="none",
        )

    svg = canvas.tostring()
    _assert_unique_ids_and_resolved_references(svg)
    slotted_groups = {
        element.attrib["data-gbdraw-slot-id"]: element
        for element in _root(svg).iter()
        if element.attrib.get("data-gbdraw-slot-id")
    }
    assert set(slotted_groups) == {"a b", "a-b"}
    assert all(
        group.attrib["id"].startswith("track_slot_")
        for group in slotted_groups.values()
    )
    assert {
        group.attrib["data-gbdraw-slot-renderer"]
        for group in slotted_groups.values()
    } == {"annotations"}
    assert slotted_groups["a b"].attrib["id"] != slotted_groups["a-b"].attrib["id"]
    assert all(" " not in group.attrib["id"] for group in slotted_groups.values())


@pytest.mark.circular
def test_circular_text_paths_are_deterministic_and_slot_scoped() -> None:
    record = SeqRecord(
        Seq("A" * 40_000),
        id="circular-text-path-record",
        name="circular-text-path-record",
    )
    record.features = [
        SeqFeature(
            SimpleLocation(1_000, 16_000, strand=1),
            type="CDS",
            qualifiers={
                "product": ["alpha"],
                "translation": ["M" * 100],
            },
        ),
    ]
    common = {
        "legend": "none",
        "cfg": apply_config_overrides(None, {
            "labels.circular.scope": "outer",
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
        "circular_track_slots": [
            CircularTrackSlot("features a", "features", side="inside"),
            CircularTrackSlot("features_a", "features", side="outside"),
            CircularTrackSlot(
                "ticks a",
                "ticks",
                side="outside",
                params={"tick_label_layout": "label_out_tick_in"},
            ),
            CircularTrackSlot(
                "ticks_a",
                "ticks",
                side="inside",
                params={"tick_label_layout": "label_in_tick_out"},
            ),
        ],
    }

    first = assemble_circular_diagram_from_record(record, **common).tostring()
    second = assemble_circular_diagram_from_record(record, **common).tostring()
    assert _ids(first) == _ids(second)
    _assert_unique_ids_and_resolved_references(first)

    root = _root(first)
    all_ids = _ids(first)
    assert not any(re.fullmatch(r"id\d+", element_id) for element_id in all_ids)

    slotted_groups = {
        element.attrib["data-gbdraw-slot-id"]: element
        for element in root.iter()
        if element.attrib.get("data-gbdraw-slot-id")
    }
    expected_prefix_by_slot = {
        "features a": "circular_feature_label_path_",
        "features_a": "circular_feature_label_path_",
        "ticks a": "circular_tick_label_path_",
        "ticks_a": "circular_tick_label_path_",
    }
    path_ids_by_slot: dict[str, set[str]] = {}
    for slot_id, prefix in expected_prefix_by_slot.items():
        path_ids_by_slot[slot_id] = {
            element.attrib["id"]
            for element in slotted_groups[slot_id].iter()
            if element.attrib.get("id", "").startswith(prefix)
        }
        assert path_ids_by_slot[slot_id]
    assert {
        slot_id: group.attrib["data-gbdraw-slot-renderer"]
        for slot_id, group in slotted_groups.items()
    } == {
        "features a": "features",
        "features_a": "features",
        "ticks a": "ticks",
        "ticks_a": "ticks",
    }

    assert path_ids_by_slot["features a"].isdisjoint(
        path_ids_by_slot["features_a"]
    )
    assert path_ids_by_slot["ticks a"].isdisjoint(path_ids_by_slot["ticks_a"])


@pytest.mark.circular
def test_circular_numeric_alias_slots_namespace_outer_and_child_ids() -> None:
    record = _duplicate_feature_record("circular-skew-slot-record")
    common = {
        "legend": "none",
        "cfg": apply_config_overrides(None, {
            "labels.circular.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": True,
        }),
        "circular_track_slots": [
            CircularTrackSlot(
                "a b",
                "dinucleotide_skew",
                side="inside",
                params={"nt": "GC"},
            ),
            CircularTrackSlot(
                "a_b",
                "dinucleotide_skew",
                side="outside",
                params={"nt": "AT"},
            ),
        ],
        "window": 40,
        "step": 40,
    }

    first = assemble_circular_diagram_from_record(record, **common).tostring()
    second = assemble_circular_diagram_from_record(record, **common).tostring()
    assert _ids(first) == _ids(second)
    _assert_unique_ids_and_resolved_references(first)

    slotted_groups = {
        element.attrib["data-gbdraw-slot-id"]: element
        for element in _root(first).iter()
        if element.attrib.get("data-gbdraw-slot-id")
    }
    assert set(slotted_groups) == {"a b", "a_b"}
    assert all(
        group.attrib["id"].startswith("track_slot_")
        for group in slotted_groups.values()
    )
    assert {
        group.attrib["data-gbdraw-slot-renderer"]
        for group in slotted_groups.values()
    } == {"dinucleotide_skew"}
    assert slotted_groups["a b"].attrib["id"] != slotted_groups["a_b"].attrib["id"]

    child_ids = [
        element.attrib["id"]
        for group in slotted_groups.values()
        for element in group.iter()
        if element is not group and element.attrib.get("id")
    ]
    assert len(child_ids) == 2
    assert all(child_id.startswith("track_slot_child_") for child_id in child_ids)


@pytest.mark.circular
def test_repeated_circular_feature_and_tick_slots_have_unique_deterministic_ids() -> None:
    record = _duplicate_feature_record("slot-record")
    common = {
        "legend": "none",
        "cfg": apply_config_overrides(None, {
            "labels.circular.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
        "circular_track_slots": [
            "features_inner:features@side=inside,w=20px",
            "features_outer:features@side=outside,w=20px",
            "ticks_inner:ticks@side=inside,tick_label_layout=tick_only",
            "ticks_outer:ticks@side=outside,tick_label_layout=tick_only",
        ],
    }

    first = assemble_circular_diagram_from_record(record, **common).tostring()
    second = assemble_circular_diagram_from_record(record, **common).tostring()
    assert _ids(first) == _ids(second)
    _assert_unique_ids_and_resolved_references(first)

    root = _root(first)
    slotted_groups = {
        element.attrib["data-gbdraw-slot-id"]: element
        for element in root.iter()
        if element.attrib.get("data-gbdraw-slot-id")
    }
    assert set(slotted_groups) == {
        "features_inner",
        "features_outer",
        "ticks_inner",
        "ticks_outer",
    }
    assert {
        slotted_groups[slot_id].attrib["id"]
        for slot_id in ("features_inner", "features_outer")
    } != {"features_inner", "features_outer"}
    assert {
        slotted_groups[slot_id].attrib["id"]
        for slot_id in ("ticks_inner", "ticks_outer")
    } != {"ticks_inner", "ticks_outer"}
    assert all(
        group.attrib["id"].startswith("track_slot_")
        for group in slotted_groups.values()
    )
    assert all(
        slotted_groups[slot_id].attrib["data-gbdraw-record-id"] == record.id
        for slot_id in ("features_inner", "features_outer")
    )

    stable_id = compute_feature_hash(record.features[0], record_id=record.id)
    matching_features = [
        element
        for element in root.iter()
        if element.attrib.get("data-gbdraw-feature-id") == stable_id
        and element.attrib.get("data-gbdraw-feature-part") == "block"
    ]
    assert len(matching_features) == 4
    assert len({element.attrib["id"] for element in matching_features}) == 4
    assert all("__instance_slot_" in element.attrib["id"] for element in matching_features)


@pytest.mark.circular
def test_reserved_circular_slot_names_are_semantic_only_not_dom_ids() -> None:
    record = _duplicate_feature_record("reserved-slot-record")
    canvas = assemble_circular_diagram_from_record(
        record,
        legend="right",
        plot_title="Reserved slot names",
        plot_title_position="top",
        cfg=apply_config_overrides(None, {
            "labels.circular.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
        circular_track_slots=[
            "legend:features@side=inside,w=20px",
            "plot_title:features@side=outside,w=20px",
            "Axis:ticks@side=inside,tick_label_layout=tick_only",
            "tick:ticks@side=outside,tick_label_layout=tick_only",
            "label_text:dinucleotide_content@side=inside,w=20px",
        ],
    )
    svg = canvas.tostring()
    _assert_unique_ids_and_resolved_references(svg)

    root = _root(svg)
    slotted_groups = [
        element
        for element in root.iter()
        if element.attrib.get("data-gbdraw-slot-id")
    ]
    assert {
        element.attrib["data-gbdraw-slot-id"]
        for element in slotted_groups
    } == {"legend", "plot_title", "Axis", "tick", "label_text"}
    assert all(
        element.attrib["id"].startswith("track_slot_")
        for element in slotted_groups
    )
    assert all(
        element.attrib["id"] != element.attrib["data-gbdraw-slot-id"]
        for element in slotted_groups
    )
