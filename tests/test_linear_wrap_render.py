"""A wrapped record drawn: one record, one identity per feature, rows inside it, each with its own height."""

from __future__ import annotations

from collections import Counter
from pathlib import Path
import re
import subprocess
import sys
import time
from xml.etree import ElementTree

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
import pandas as pd
import pytest

import gbdraw.analysis.skew as skew_module
from gbdraw.api import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    RegionAnnotation,
    render_request,
)
from gbdraw.api.options import LinearDiagramOptions, LinearRequestTrackOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import (
    GenBankInputSource,
    LinearDiagramRequest,
    RecordDisplayOptions,
    RecordInput,
    RecordPresentation,
    RenderOutputRequest,
)
from gbdraw.auto_sizes import DrawnExtent, auto_setting_values
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.toml import load_config_toml
import gbdraw.diagrams.linear.assemble as assemble
import gbdraw.render.drawers.linear.depth as depth_drawer
import gbdraw.render.drawers.linear.gc_content as gc_drawer
from gbdraw.io.regions import parse_region_spec
from gbdraw.layout.linear import VerticalBand
from gbdraw.session import build_session_document
from gbdraw.web_support.request_render import render_embedded_canonical_web_request
from gbdraw.layout.linear_wrap import LinearWrapSummary, drawn_axis_length, stack_row_offsets

WRAP = "canvas.linear.wrap_length"
LAMBDA = Path(__file__).parent / "test_inputs" / "NC_001416.gb"
LAMBDA_LENGTH = 48_502
SVG = "{http://www.w3.org/2000/svg}"
WRAP_ROW = "data-gbdraw-wrap-row"
GC_ON = {"canvas.show_gc": True, "canvas.show_skew": True, "labels.linear.scope": "all"}


def _build(
    *,
    source: Path = LAMBDA,
    also: tuple[Path, ...] = (),
    wrap: object = 10_000,
    overrides: dict[str, object] | None = None,
    display: RecordDisplayOptions | None = None,
    presentation: RecordPresentation | None = None,
    region: str | None = None,
    **options: object,
):
    plan = plan_request(
        LinearDiagramRequest(
            records=tuple(
                RecordInput(
                    GenBankInputSource(str(path)),
                    record_key=f"r{index}",
                    region=parse_region_spec(region) if region else None,
                    display=display or RecordDisplayOptions(),
                    presentation=presentation or RecordPresentation(),
                )
                for index, path in enumerate((source, *also))
            ),
            options=LinearDiagramOptions(
                config_overrides={**({WRAP: wrap} if wrap is not None else {}), **(overrides or {})},
                **options,
            ),
        )
    )
    return plan.build()


def _root(built) -> ElementTree.Element:
    return ElementTree.fromstring(built.drawing.tostring())


def _ids(root: ElementTree.Element) -> Counter[str]:
    """Every id but a clip path's: each row's skew path has its own clip."""

    return Counter(
        element.attrib["id"] for element in root.iter() if "id" in element.attrib and not element.tag.endswith("clipPath")
    )


def _folded(ids: Counter[str]) -> set[str]:
    """Dom ids with a block's or connector's part number dropped: a row boundary cuts a block as a seam does."""

    return {re.sub(r"__(part|line)\d+(__outline)?$", "", identifier) for identifier in ids}


def _feature_ids(root: ElementTree.Element) -> set[str]:
    return {
        element.attrib["data-gbdraw-feature-id"]
        for element in root.iter()
        if "data-gbdraw-feature-id" in element.attrib
    }


def _rows_of(element: ElementTree.Element) -> list[ElementTree.Element]:
    return [item for item in element.iter() if "data-gbdraw-wrap-row" in item.attrib]


def _translate_y(element: ElementTree.Element) -> float:
    match = re.search(r"translate\(\s*[-0-9.e]+[ ,]+([-0-9.e]+)\s*\)", element.attrib.get("transform", ""))
    return float(match.group(1)) if match else 0.0


def _record_groups(root: ElementTree.Element) -> list[ElementTree.Element]:
    return [element for element in root.iter(f"{SVG}g") if element.attrib.get("id", "").startswith("record_group_")]


def _path_points(data: str) -> list[tuple[float, float]]:
    return [(float(x), float(y)) for x, y in re.findall(r"(-?\d*\.?\d+)[ ,]+(-?\d*\.?\d+)", data)]


def _element_ys(element: ElementTree.Element) -> list[float]:
    """The y coordinates an element draws at, in its own frame."""

    tag = element.tag.removeprefix(SVG)
    attributes = element.attrib
    if tag == "path":
        return [y for _, y in _path_points(attributes["d"])]
    if tag == "line":
        return [float(attributes["y1"]), float(attributes["y2"])]
    if tag == "rect":
        return [float(attributes["y"]), float(attributes["y"]) + float(attributes["height"])]
    if tag == "text":
        half = float(attributes.get("font-size", 0)) / 2.0
        return [float(attributes["y"]) - half, float(attributes["y"]) + half]
    return []


def _shift(element: ElementTree.Element) -> float:
    return sum(float(y) for y in re.findall(r"translate\(\s*[-0-9.e]+[ ,]+([-0-9.e]+)\s*\)", element.attrib.get("transform", "")))


def _ink_by_row(root: ElementTree.Element) -> tuple[dict[int, tuple[float, float]], list[tuple[float, float]]]:
    """Each row's drawn vertical extent on the canvas, read from the drawn elements, and each feature label's."""

    rows: dict[int, list[float]] = {}
    labels: list[tuple[float, float]] = []

    def walk(element: ElementTree.Element, offset: float, row: int | None) -> None:
        if element.tag == f"{SVG}defs" or element.tag == f"{SVG}clipPath":
            return
        offset += _shift(element)
        row = int(element.attrib[WRAP_ROW]) if WRAP_ROW in element.attrib else row
        ys = [offset + y for y in _element_ys(element)]
        if ys and "data-label-feature-id" in element.attrib:
            labels.append((min(ys), max(ys)))
        elif ys and row is not None:
            rows.setdefault(row, []).extend(ys)
        for child in element:
            walk(child, offset, row)

    walk(root, 0.0, None)
    return {row: (min(ys), max(ys)) for row, ys in rows.items()}, labels


def _assert_rows_do_not_overlap(root: ElementTree.Element, count: int) -> None:
    rows, labels = _ink_by_row(root)
    assert sorted(rows) == list(range(count))
    extents = [rows[row] for row in sorted(rows)]
    assert all(above[1] < below[0] for above, below in zip(extents, extents[1:]))
    # A label stands in its row's band or in the gap above it, so it meets the ink of at most one row.
    for top, bottom in labels:
        assert sum(1 for low, high in extents if top < high and low < bottom) <= 1


def _labels(root: ElementTree.Element) -> Counter[str]:
    return Counter(
        element.attrib["data-label-feature-id"]
        for element in root.iter(f"{SVG}text")
        if "data-label-feature-id" in element.attrib
    )


def _join_record(path: Path) -> Path:
    """A 30 kb record with a CDS and a join() across row cuts at 10 kb, and small genes on one row only."""

    record = SeqRecord(Seq("ATGC" * 7500), id="JOIN1", name="JOIN1", description="wrap fixture")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "circular"
    record.features.append(SeqFeature(SimpleLocation(8000, 14000, strand=1), type="CDS", qualifiers={"gene": ["cut"]}))
    record.features.append(SeqFeature(SimpleLocation(2000, 11000, strand=-1), type="CDS", qualifiers={"gene": ["back"]}))
    record.features.append(
        SeqFeature(
            CompoundLocation([SimpleLocation(19_000, 19_500, strand=1), SimpleLocation(20_500, 21_500, strand=1)]),
            type="CDS",
            qualifiers={"gene": ["spliced"]},
        )
    )
    # Row 2 only: many short genes, so its external labels stack and it grows taller than the others.
    for index in range(30):
        start = 22_000 + index * 40
        record.features.append(
            SeqFeature(SimpleLocation(start, start + 30, strand=1), type="CDS", qualifiers={"gene": [f"crowded_gene_{index}"]})
        )
    SeqIO.write(record, path, "genbank")
    return path


def _geometry(built, record_index: int = 0) -> dict:
    """The resolved slot geometry a record reports (the same payload the Web reads)."""

    return built.drawing._gbdraw_track_slot_geometry["records"][record_index]


def _drawn_labels(monkeypatch, **build_options: object) -> list:
    """The record's feature labels as the assembly placed them."""

    seen: list = []
    original = assemble.linear_label_fit

    def spy(by_record, *args, **kwargs):
        seen.extend(by_record[0][1])
        return original(by_record, *args, **kwargs)

    monkeypatch.setattr(assemble, "linear_label_fit", spy)
    _build(**build_options)
    return seen


# The owners of the row arithmetic.


def test_drawn_axis_length_is_the_row_width_only_on_a_wrapped_record() -> None:
    built = _build()
    record = plan_request(
        LinearDiagramRequest(records=(RecordInput(GenBankInputSource(str(LAMBDA))),), options=LinearDiagramOptions())
    )
    assert drawn_axis_length(record.records[0], None) == LAMBDA_LENGTH
    assert drawn_axis_length(record.records[0], record.transforms[0]) == LAMBDA_LENGTH
    assert drawn_axis_length(record.records[0], record.transforms[0].with_wrap(10_000, LAMBDA_LENGTH)) == 10_000
    cropped = plan_request(
        LinearDiagramRequest(
            records=(RecordInput(GenBankInputSource(str(LAMBDA)), region=parse_region_spec("NC_001416.1:5000-35000")),),
            options=LinearDiagramOptions(),
        )
    )
    assert drawn_axis_length(cropped.records[0], cropped.transforms[0]) == 30_001
    assert built.metadata.wrap == LinearWrapSummary(10_000, 10_000, (LAMBDA_LENGTH,))


def test_stacking_adds_each_rows_own_height() -> None:
    bands = [VerticalBand(-50.0, 10.0), VerticalBand(-120.0, 10.0), VerticalBand(-50.0, 10.0)]
    assert stack_row_offsets(bands, 5.0) == (0.0, 135.0, 200.0)


# One record after drawing.


@pytest.mark.parametrize("overrides", [{}, GC_ON], ids=["features", "gc_skew_labels"])
def test_one_record_group_and_one_identity_per_feature(overrides) -> None:
    wrapped, plain = _root(_build(overrides=overrides)), _root(_build(wrap=None, overrides={
        **overrides, "canvas.linear.track_layout": "above", "objects.scale.style": "ruler", "canvas.linear.ruler_on_axis": True
    }))
    assert _feature_ids(wrapped) == _feature_ids(plain)
    assert _folded(_ids(wrapped)) == _folded(_ids(plain))
    assert not [identifier for identifier in _ids(wrapped) if "_row_" in identifier]
    assert [identifier for identifier, count in _ids(wrapped).items() if count > 1] == []
    assert not [element for element in _rows_of(wrapped) if element.tag == f"{SVG}g" and "id" in element.attrib]
    assert len(_record_groups(wrapped)) == 1
    assert [group.attrib["id"] for group in _record_groups(wrapped)] == [group.attrib["id"] for group in _record_groups(plain)]


def test_two_cut_records_draw_two_record_groups() -> None:
    assert len(_record_groups(_root(_build(also=(LAMBDA,))))) == 2


def test_every_row_has_an_axis() -> None:
    root = _root(_build())
    axes = [
        element
        for element in root.iter(f"{SVG}g")
        if element.attrib.get("data-gbdraw-wrap-row") is not None
        and element.attrib.get("data-gbdraw-wrap-part") == "axis"
    ]
    assert [element.attrib["data-gbdraw-wrap-row"] for element in axes] == ["0", "1", "2", "3", "4"]
    offsets = [_translate_y(element) for element in axes]
    assert offsets == sorted(offsets) and offsets[0] == 0.0 and len(set(offsets)) == 5
    widths = [float(next(element.iter(f"{SVG}line")).attrib["x2"]) for element in axes]
    # The last row holds 8,502 bp at the scale of the others.
    assert widths[:4] == [widths[0]] * 4
    assert widths[4] == pytest.approx(widths[0] * 8_502 / 10_000)


def test_definition_is_drawn_once_with_the_whole_record_length() -> None:
    svg = _build().drawing.tostring()
    assert svg.count("48,502 bp") == 1


def test_gc_is_computed_once_per_record(monkeypatch) -> None:
    calls = []
    original = skew_module.skew_df

    def counting(*args, **kwargs):
        calls.append(args[0].id)
        return original(*args, **kwargs)

    monkeypatch.setattr(assemble, "skew_df", counting)
    _build(overrides={"canvas.show_gc": True, "canvas.show_skew": True})
    assert len(calls) == 1


def _depth_table() -> pd.DataFrame:
    positions = range(1, LAMBDA_LENGTH, 50)
    return pd.DataFrame({"reference_name": "NC_001416.1", "position": positions, "depth": [float(p % 40) for p in positions]})


def _counting(monkeypatch, module: object, name: str) -> list[int]:
    calls: list[int] = []
    original = getattr(module, name)

    def counting(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(module, name, counting)
    return calls


def test_depth_plot_data_is_prepared_once_per_record(monkeypatch) -> None:
    calls = _counting(monkeypatch, depth_drawer, "_prepare_depth_plot_dataframe")
    root = _root(_build(depth_table=_depth_table(), overrides={"canvas.show_depth": True}))
    assert len(_rows_of(root)) >= 5
    assert len(calls) == 1


def test_gc_percent_plot_data_is_prepared_once_per_record(monkeypatch) -> None:
    calls = _counting(monkeypatch, gc_drawer, "gc_content_percent_df")
    _build(overrides={"canvas.show_gc": True, "objects.gc_content.mode": "percent"})
    assert len(calls) == 1


def _half_gc_record(path: Path) -> Path:
    """20 kb: all GC, then all AT. The record's mean is 0.5, while each half's own mean is its only value."""

    record = SeqRecord(Seq("GC" * 5_000 + "AT" * 5_000), id="HALF1", name="HALF1", description="half")
    record.annotations["molecule_type"] = "DNA"
    SeqIO.write(record, path, "genbank")
    return path


def test_each_row_draws_gc_against_the_whole_records_mean(tmp_path: Path) -> None:
    root = _root(_build(source=_half_gc_record(tmp_path / "half.gb"), overrides={"canvas.show_gc": True}))
    slot = next(element for element in root.iter(f"{SVG}g") if element.attrib.get("data-gbdraw-slot-renderer") == "dinucleotide_content")
    deviations = {}
    for row in slot:
        if WRAP_ROW not in row.attrib:
            continue
        points = _path_points(next(row.iter(f"{SVG}path")).attrib["d"])
        baseline = points[0][1]
        # The path opens and closes on the baseline: every sample lies between. A row's own mean would draw it flat.
        deviations[row.attrib[WRAP_ROW]] = sorted(y - baseline for _, y in points[1:-2])
    assert deviations["0"] and deviations["1"]
    medians = {row: values[len(values) // 2] for row, values in deviations.items()}
    # Most of each row is one composition, a full half track from the record's mean, in opposite directions.
    # Against its own row's mean, nearly all of a row would lie close to the baseline.
    assert medians["0"] < 0.0 < medians["1"]
    assert medians["0"] == pytest.approx(-medians["1"])
    assert abs(medians["0"]) == pytest.approx(max(abs(value) for value in deviations["0"]))


def test_gc_and_skew_have_one_path_per_row() -> None:
    root = _root(_build(overrides=GC_ON))
    slots = {
        element.attrib["data-gbdraw-slot-id"]: element
        for element in root.iter(f"{SVG}g")
        if element.attrib.get("data-gbdraw-slot-renderer") in {"dinucleotide_content", "dinucleotide_skew"}
    }
    assert len(slots) == 2
    for slot in slots.values():
        rows = [element for element in slot if "data-gbdraw-wrap-row" in element.attrib]
        assert [element.attrib["data-gbdraw-wrap-row"] for element in rows] == ["0", "1", "2", "3", "4"]
        assert all(next(element.iter(f"{SVG}path")).attrib["d"] for element in rows)
        assert [_translate_y(element) for element in rows] == sorted(_translate_y(element) for element in rows)
        assert len({_translate_y(element) for element in rows}) == 5


def test_a_cut_feature_keeps_its_identity_over_its_rows(tmp_path: Path) -> None:
    source = _join_record(tmp_path / "join.gb")
    root = _root(_build(source=source, overrides={"labels.linear.scope": "all"}))
    paths = [element for element in root.iter(f"{SVG}path") if element.attrib.get("data-gbdraw-feature-part") == "block"]
    by_feature: dict[str, list[ElementTree.Element]] = {}
    for path in paths:
        by_feature.setdefault(path.attrib["data-gbdraw-feature-id"], []).append(path)
    parents = {child: parent for parent in root.iter() for child in parent}
    rows = sorted(
        len({parents[part].attrib["data-gbdraw-wrap-row"] for part in parts}) for parts in by_feature.values()
    )
    # "cut", "back" and "spliced" each draw on two rows; every other feature on one.
    assert rows[-3:] == [2, 2, 2] and rows[-4] == 1


def test_label_is_drawn_once_on_the_row_of_the_longest_fragment(tmp_path: Path, monkeypatch) -> None:
    source = _join_record(tmp_path / "join.gb")
    options = {"source": source, "overrides": {"labels.linear.scope": "all"}}
    texts = [element.text for element in _root(_build(**options)).iter(f"{SVG}text")]
    assert texts.count("cut") == 1
    assert texts.count("back") == 1
    labels = {label["label_text"]: label for label in _drawn_labels(monkeypatch, **options)}
    # "cut" holds 2 kb on row 0 and 4 kb on row 1; "back" 8 kb on row 0 and 1 kb on row 1.
    assert labels["cut"]["wrap_row"] == 1
    assert labels["back"]["wrap_row"] == 0


def test_a_tie_goes_to_the_earlier_row(tmp_path: Path, monkeypatch) -> None:
    record = SeqRecord(Seq("ATGC" * 5000), id="TIE1", name="TIE1", description="tie")
    record.annotations["molecule_type"] = "DNA"
    record.features.append(SeqFeature(SimpleLocation(9000, 11000, strand=1), type="CDS", qualifiers={"gene": ["tie"]}))
    path = tmp_path / "tie.gb"
    SeqIO.write(record, path, "genbank")
    labels = _drawn_labels(monkeypatch, source=path, overrides={"labels.linear.scope": "all"})
    assert labels[0]["wrap_row"] == 0


# Each row its own height.


def test_a_row_with_more_labels_is_taller_and_a_sparse_one_is_not_padded(tmp_path: Path) -> None:
    source = _join_record(tmp_path / "join.gb")
    wrap = _geometry(_build(source=source, overrides={"labels.linear.scope": "all"}))["wrap"]
    offsets = wrap["rowOffsets"]
    assert len(offsets) == 3
    # Row 2 holds the crowded genes: its labels stand above its axis, so the step down to it is the largest.
    assert offsets[2] - offsets[1] > offsets[1] - offsets[0]
    row_heights = [band["bottomPx"] - band["topPx"] for band in wrap["rowCanvasBands"]]
    assert row_heights[2] > row_heights[1]


def test_the_records_extent_is_the_sum_of_its_rows_and_gaps(tmp_path: Path) -> None:
    source = _join_record(tmp_path / "join.gb")
    built = _build(source=source, overrides={"labels.linear.scope": "all"})
    record = _geometry(built)
    bands = record["wrap"]["rowCanvasBands"]
    heights = [band["bottomPx"] - band["topPx"] for band in bands]
    extent = record["canvasBand"]["bottomPx"] - record["canvasBand"]["topPx"]
    gaps = [below["topPx"] - above["bottomPx"] for above, below in zip(bands, bands[1:])]
    assert all(gap > 0.0 for gap in gaps) and gaps[0] == pytest.approx(gaps[-1])
    assert extent == pytest.approx(sum(heights) + sum(gaps))
    # Not a uniform pitch: the tallest row times the rows would be more.
    assert extent < max(heights) * len(heights) + sum(gaps) - 1.0
    # The canvas holds that extent.
    assert float(built.drawing.attribs["height"].removesuffix("px")) >= extent




@pytest.mark.parametrize(
    "fields",
    [
        {"overrides": GC_ON},
        {"overrides": {"labels.linear.scope": "all"}},
        {"overrides": {**GC_ON, "canvas.linear.track_layout": "below"}},
        {"overrides": {"canvas.strandedness": True, "labels.linear.scope": "all"}},
        {
            "overrides": GC_ON,
            "tracks": LinearRequestTrackOptions(
                linear_track_slots=(
                    "gc:dinucleotide_content@side=above",
                    "features:features@side=below",
                    "skew:dinucleotide_skew@side=below",
                )
            ),
        },
    ],
    ids=["tracks", "labels", "below", "separate_strands", "custom_slots"],
)
def test_rows_do_not_overlap_in_the_drawn_output(fields) -> None:
    _assert_rows_do_not_overlap(_root(_build(**fields)), 5)


def test_labels_of_the_first_record_only_on_rows() -> None:
    root = _root(_build(also=(LAMBDA,), overrides={"labels.linear.scope": "first"}))
    labels_by_group = [
        sum(1 for element in group.iter(f"{SVG}text") if "data-label-feature-id" in element.attrib)
        for group in _record_groups(root)
    ]
    assert labels_by_group[0] > 0 and labels_by_group[1] == 0
    assert len([element for group in _record_groups(root) for element in group.iter() if element.attrib.get("data-gbdraw-wrap-part") == "axis"]) == 10


def test_separate_strands_keep_both_strands_on_each_row() -> None:
    root = _root(_build(overrides={"canvas.strandedness": True}))
    plain = _root(_build(wrap=None, overrides={"canvas.strandedness": True}))
    assert _feature_ids(root) == _feature_ids(plain)
    tops_by_row: dict[str, set[float]] = {}
    parents = {child: parent for parent in root.iter() for child in parent}
    for path in root.iter(f"{SVG}path"):
        if path.attrib.get("data-gbdraw-feature-part") == "block":
            row = parents[path].attrib[WRAP_ROW]
            tops_by_row.setdefault(row, set()).add(min(y for _, y in _path_points(path.attrib["d"])))
    # Each strand has its own lane on a row that holds both.
    assert any(len(tops) >= 2 for tops in tops_by_row.values())


def test_a_custom_slot_order_reaches_every_row() -> None:
    built = _build(
        overrides=GC_ON,
        tracks=LinearRequestTrackOptions(
            linear_track_slots=("features:features@side=above", "gc:dinucleotide_content@side=below"),
        ),
    )
    root = _root(built)
    slot = next(element for element in root.iter(f"{SVG}g") if element.attrib.get("data-gbdraw-slot-id") == "gc")
    assert [row.attrib[WRAP_ROW] for row in slot if WRAP_ROW in row.attrib] == ["0", "1", "2", "3", "4"]
    assert not [element for element in root.iter(f"{SVG}g") if element.attrib.get("data-gbdraw-slot-id") == "gc_skew"]


def test_the_slot_geometry_payload_carries_the_row_offsets(tmp_path: Path) -> None:
    source = _join_record(tmp_path / "join.gb")
    wrap = _geometry(_build(source=source, overrides={"labels.linear.scope": "all"}))["wrap"]
    assert wrap["rowOffsets"][0] == 0.0 and len(wrap["rowOffsets"]) == len(wrap["rowCanvasBands"]) == 3
    assert "rowPitch" not in wrap
    assert "wrap" not in _geometry(_build(wrap=None))


def test_a_label_wider_than_its_fragment_keeps_a_neighbour_out_of_its_lane(recwarn) -> None:
    """The label stands on the longest fragment: its width is measured against that, not the whole span."""

    root = _root(
        _build(
            annotations=_annotations(
                # 900 bp on row 0 and 1,000 bp on row 1, where the label stands, wider than those 1,000 bp.
                RegionAnnotation(id="across", target=CoordinateSpan(None, 9_100, 11_000), label="Across xxxxxxxxxxxxxxxxxxx", mark="band"),
                RegionAnnotation(id="neighbour", target=CoordinateSpan(None, 11_040, 11_500), label="Nbr", mark="band"),
            )
        )
    )
    top = {
        group.attrib["data-gbdraw-annotation-id"]: next(group.iter(f"{SVG}rect")).attrib["y"]
        for group in root.iter(f"{SVG}g")
        if "data-gbdraw-annotation-id" in group.attrib
    }
    assert top["across"] != top["neighbour"]
    assert not [warning for warning in recwarn if "do not fit" in str(warning.message)]


def test_a_clipped_annotation_track_clips_per_row_at_the_slots_height() -> None:
    built = _build(
        # Overlapping marks take more lanes than the slot's height holds.
        annotations=_annotations(*(_region(f"whole{index}", 2_000 + index, 45_000 + index) for index in range(8))),
        tracks=LinearRequestTrackOptions(
            linear_track_slots=("regions:annotations@set_id=regions,side=above,h=20px,overflow=clip", "features:features"),
        ),
    )
    root = _root(built)
    clipped = next(group for group in root.iter(f"{SVG}g") if group.attrib.get("data-gbdraw-annotation-clipped"))
    clip_id = clipped.attrib["clip-path"].removeprefix("url(#").removesuffix(")")
    rectangles = list(next(clip for clip in root.iter(f"{SVG}clipPath") if clip.attrib["id"] == clip_id))
    assert [float(rect.attrib["height"]) for rect in rectangles] == [20.0] * 5
    offsets = _geometry(built)["wrap"]["rowOffsets"]
    tops = [float(rect.attrib["y"]) for rect in rectangles]
    assert [top - tops[0] for top in tops] == pytest.approx(offsets)


# The row length sets the size class, and one row length serves every record.


def _sparse_record(path: Path, length: int) -> Path:
    record = SeqRecord(Seq("ATGC" * (length // 4)), id="LONG1", name="LONG1", description="long")
    record.annotations["molecule_type"] = "DNA"
    record.features.append(SeqFeature(SimpleLocation(2_000, 8_000, strand=1), type="CDS", qualifiers={"gene": ["a"]}))
    SeqIO.write(record, path, "genbank")
    return path


def _block_height(built) -> float:
    root = _root(built)
    block = next(path for path in root.iter(f"{SVG}path") if path.attrib.get("data-gbdraw-feature-part") == "block")
    ys = [y for _, y in _path_points(block.attrib["d"])]
    return max(ys) - min(ys)


def test_the_drawn_size_class_follows_the_row_length(tmp_path: Path) -> None:
    long_record = _sparse_record(tmp_path / "long.gb", 300_000)
    short_record = _sparse_record(tmp_path / "short.gb", 50_000)
    above = {"canvas.linear.track_layout": "above"}
    long_height = _block_height(_build(source=long_record, wrap=None, overrides=above))
    short_height = _block_height(_build(source=short_record, wrap=None, overrides=above))
    assert long_height != short_height
    assert _block_height(_build(source=long_record, wrap=50_000)) == pytest.approx(short_height)


def test_the_size_class_of_an_extent_reads_the_row_length_and_the_window_the_record() -> None:
    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    whole = DrawnExtent("linear", (300_000,))
    rows = DrawnExtent("linear", (300_000,), wrap_row_length=50_000)
    assert whole.size_class(cfg) == "long"
    assert rows.size_class(cfg) == "short"
    assert auto_setting_values(rows, cfg)["window"] == auto_setting_values(whole, cfg)["window"]


def test_the_figure_scale_comes_from_the_row_length(tmp_path: Path) -> None:
    short = SeqRecord(Seq("ATGC" * 1000), id="SHORT1", name="SHORT1", description="short")
    short.annotations["molecule_type"] = "DNA"
    path = tmp_path / "short.gb"
    SeqIO.write(short, path, "genbank")
    built = _build(also=(path,), wrap=10_000)
    assert len(_geometry(built)["wrap"]["rowOffsets"]) == 5
    assert "wrap" not in _geometry(built, 1)
    assert built.metadata.wrap.record_rows == (5, 1)
    assert _geometry(built, 1)["sequenceWidthPx"] == pytest.approx(_geometry(built)["sequenceWidthPx"] * 4_000 / 10_000)


@pytest.mark.parametrize("overrides", [{}, GC_ON], ids=["default", "tracks"])
def test_nothing_to_cut_draws_exactly_the_same_request_without_a_wrap_length(overrides) -> None:
    """Same layout and scale style as the plain request (the CLI's byte test is in ``test_linear_wrap_plan``)."""

    wrapped = _build(wrap=100_000, overrides=overrides)
    assert wrapped.drawing.tostring() == _build(wrap=None, overrides=overrides).drawing.tostring()
    assert wrapped.metadata.wrap is None


# Rotated, reversed, cropped, annotation tracks, depth.


@pytest.mark.parametrize(
    "fields",
    [
        {"display": RecordDisplayOptions(is_circular=True, start_coordinate=20_000)},
        {"presentation": RecordPresentation(reverse_complement=True)},
        {"region": "NC_001416.1:5000-35000"},
        {"region": "NC_001416.1:5000-35000", "presentation": RecordPresentation(reverse_complement=True)},
    ],
    ids=["rotated", "reversed", "cropped", "cropped_reversed"],
)
def test_transformed_records_wrap(fields) -> None:
    root = _root(_build(overrides=GC_ON, **fields))
    assert [identifier for identifier, count in _ids(root).items() if count > 1] == []
    assert {element.attrib["data-gbdraw-wrap-row"] for element in _rows_of(root)} >= {"0", "1", "2"}


def _annotations(*items: RegionAnnotation) -> AnnotationOptions:
    return AnnotationOptions(sets=(AnnotationSet("regions", tuple(items)),))


def _region(name: str, start: int, end: int, mark: str = "band") -> RegionAnnotation:
    return RegionAnnotation(id=name, target=CoordinateSpan(None, start, end), label=name, mark=mark)


@pytest.mark.parametrize("mark", ["band", "bracket", "highlight"])
def test_annotation_cut_at_the_seam_stays_one_item_labelled_once(mark: str) -> None:
    built = _build(annotations=_annotations(_region("across", 8_000, 13_000, mark)))
    root = _root(built)
    items = [element for element in root.iter(f"{SVG}g") if element.attrib.get("data-gbdraw-annotation-id") == "across"]
    assert len(items) == 1
    rows = {element.attrib["data-gbdraw-wrap-row"] for element in _rows_of(items[0])}
    assert rows == {"0", "1"}
    assert [element.text for element in root.iter(f"{SVG}text")].count("across") == 1


def test_annotation_lanes_are_kept_per_row() -> None:
    root = _root(
        _build(
            annotations=_annotations(
                _region("first", 1_000, 2_000),
                _region("second", 1_500, 2_500),
                # Same row-local span as "first" on the next row: it does not push a lane.
                _region("next_row", 11_000, 12_000),
            ),
        )
    )
    lane_y = {
        element.attrib["data-gbdraw-annotation-id"]: float(next(element.iter(f"{SVG}rect")).attrib["y"])
        for element in root.iter(f"{SVG}g")
        if "data-gbdraw-annotation-id" in element.attrib
    }
    assert lane_y["first"] == lane_y["next_row"] != lane_y["second"]


def test_annotation_label_stands_on_its_longest_row() -> None:
    built = _build(annotations=_annotations(_region("across", 8_000, 13_000)))
    record = _geometry(built)
    slot = next(slot for slot in record["slots"] if slot["renderer"] == "annotations")
    text = next(element for element in _root(built).iter(f"{SVG}text") if element.text == "across")
    # 2 kb on row 0, 3 kb on row 1: the label's y, from the slot's first-row origin, lies in row 1.
    y = float(text.attrib["y"]) + slot["resolvedOriginPx"]
    row = record["wrap"]["rowCanvasBands"][1]
    assert row["topPx"] <= y <= row["bottomPx"]


def test_depth_track_wraps(tmp_path: Path) -> None:
    depth = tmp_path / "depth.tsv"
    depth.write_text("".join(f"NC_001416.1\t{position}\t{position % 50}\n" for position in range(1, LAMBDA_LENGTH, 100)))
    root = _root(_build(depth_file=str(depth)))
    groups = [element for element in root.iter(f"{SVG}g") if str(element.attrib.get("id", "")).startswith("depth")]
    assert groups
    rows = [element for element in groups[0] if "data-gbdraw-wrap-row" in element.attrib]
    assert len(rows) == 5


def test_hide_scale_draws_no_ticks_and_still_wraps() -> None:
    root = _root(_build(overrides={"objects.scale.show": False}))
    assert {element.attrib["data-gbdraw-wrap-row"] for element in _rows_of(root)} == {"0", "1", "2", "3", "4"}


def test_ruler_ticks_are_placed_on_their_rows() -> None:
    root = _root(_build(overrides={"objects.scale.interval": 5_000}))
    axes = [
        element
        for element in root.iter(f"{SVG}g")
        if element.attrib.get("data-gbdraw-wrap-part") == "axis"
    ]
    labels = {
        element.attrib["data-gbdraw-wrap-row"]: [text.text for text in element.iter(f"{SVG}text")]
        for element in axes
    }
    assert labels["0"][0].replace(" ", "").lower().startswith("5")
    # Base 15,000 is on row 1, at index 4,999 in it.
    assert any(text.replace(" ", "").lower().startswith("15") for text in labels["1"])


# The surfaces.


def test_render_request_and_cli_draw_the_same_rows(tmp_path: Path) -> None:
    render_request(
        LinearDiagramRequest(
            records=(RecordInput(GenBankInputSource(str(LAMBDA))),),
            options=LinearDiagramOptions(config_overrides={WRAP: 10_000}),
            output=RenderOutputRequest(output_prefix="api", output_directory=tmp_path, formats=("svg",)),
        )
    )
    result = subprocess.run(
        [sys.executable, "-c", "from gbdraw.cli import main; main()", "linear", "--gbk", str(LAMBDA),
         "--wrap_length", "10kb", "-o", "cli", "-f", "svg"],
        cwd=tmp_path, capture_output=True, text=True,
        env={"PYTHONPATH": str(Path(__file__).parents[1]), "PATH": "/usr/bin:/bin"},
    )
    assert result.returncode == 0, result.stderr
    assert "Wrapping at 10 kb per row: 5 rows" in result.stderr + result.stdout
    api = ElementTree.parse(tmp_path / "api.svg").getroot()
    cli = ElementTree.parse(tmp_path / "cli.svg").getroot()
    assert len(_rows_of(api)) == len(_rows_of(cli)) > 0


def test_cli_log_names_auto(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, "-c", "from gbdraw.cli import main; main()", "linear", "--gbk", str(LAMBDA),
         "--wrap_length", "auto", "-o", "cli", "-f", "svg"],
        cwd=tmp_path, capture_output=True, text=True,
        env={"PYTHONPATH": str(Path(__file__).parents[1]), "PATH": "/usr/bin:/bin"},
    )
    assert result.returncode == 0, result.stderr
    assert "Wrapping at 10 kb per row (auto): 5 rows" in result.stderr + result.stdout


def test_unwrapped_metadata_says_so() -> None:
    assert _build(wrap=None).metadata.wrap is None


def test_the_web_reply_names_the_rows(tmp_path: Path) -> None:
    def reply(**overrides: object) -> dict:
        document = build_session_document(
            LinearDiagramRequest(
                records=(RecordInput(GenBankInputSource(str(LAMBDA))),),
                options=LinearDiagramOptions(config_overrides=overrides),
                output=RenderOutputRequest(output_prefix="web", formats=("svg",)),
            )
        ).to_dict()
        return render_embedded_canonical_web_request(
            document["renderRequest"], resources=document["resources"], workspace=tmp_path / "workspace"
        )["metadata"]

    assert reply(**{WRAP: "auto"})["linearWrap"] == [
        {
            "resultIndex": 0,
            "resultName": "web.svg",
            "rowLength": 10_000,
            "rowCount": 5,
            "requested": "auto",
            "recordRows": [5],
        }
    ]
    assert reply()["linearWrap"] == []


def test_a_hundred_rows_stay_fast_and_tile(tmp_path: Path) -> None:
    record = SeqRecord(Seq("ATGC" * 25_000), id="LONG1", name="LONG1", description="hundred rows")
    record.annotations["molecule_type"] = "DNA"
    for index in range(1000):
        record.features.append(
            SeqFeature(SimpleLocation(index * 100, index * 100 + 150 if index < 999 else 99_990, strand=1), type="CDS", qualifiers={"gene": [f"g{index}"]})
        )
    path = tmp_path / "long.gb"
    SeqIO.write(record, path, "genbank")
    started = time.perf_counter()
    built = _build(source=path, wrap=1_000, overrides=GC_ON)
    elapsed = time.perf_counter() - started
    assert len(_geometry(built)["wrap"]["rowOffsets"]) == 100
    # A guard against per-row work that is quadratic in the rows, not a benchmark.
    assert elapsed < 120.0
