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
import pytest

import gbdraw.analysis.skew as skew_module
from gbdraw.api import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    RegionAnnotation,
    render_request,
)
from gbdraw.api.options import LinearDiagramOptions
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
from gbdraw.io.regions import parse_region_spec
from gbdraw.layout.linear import VerticalBand
from gbdraw.session import build_session_document
from gbdraw.web_support.request_render import render_embedded_canonical_web_request
from gbdraw.layout.linear_wrap import LinearWrapSummary, drawn_axis_length, stack_row_offsets

WRAP = "canvas.linear.wrap_length"
LAMBDA = Path(__file__).parent / "test_inputs" / "NC_001416.gb"
LAMBDA_LENGTH = 48_502
SVG = "{http://www.w3.org/2000/svg}"
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
    assert built.metadata.wrap == LinearWrapSummary(10_000, 5, 10_000, (5,))


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
    groups = [element for element in wrapped.iter(f"{SVG}g") if "data-gbdraw-record-id" in element.attrib and element.attrib.get("id")]
    assert len([group for group in groups if group.attrib["id"].startswith("NC_001416")]) >= 1


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


def test_gc_and_skew_have_one_path_per_row_at_the_records_mean() -> None:
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


def test_the_records_bottom_extent_is_the_sum_of_its_rows(tmp_path: Path) -> None:
    source = _join_record(tmp_path / "join.gb")
    record = _geometry(_build(source=source, overrides={"labels.linear.scope": "all"}))
    assert record["canvasBand"]["bottomPx"] == pytest.approx(record["wrap"]["rowCanvasBands"][-1]["bottomPx"])
    assert record["canvasBand"]["topPx"] == pytest.approx(record["wrap"]["rowCanvasBands"][0]["topPx"])


def test_rows_do_not_overlap() -> None:
    bands = _geometry(_build(overrides=GC_ON))["wrap"]["rowCanvasBands"]
    assert all(below["topPx"] > above["bottomPx"] for above, below in zip(bands, bands[1:]))


def test_the_slot_geometry_payload_carries_the_row_offsets(tmp_path: Path) -> None:
    source = _join_record(tmp_path / "join.gb")
    wrap = _geometry(_build(source=source, overrides={"labels.linear.scope": "all"}))["wrap"]
    assert wrap["rowOffsets"][0] == 0.0 and len(wrap["rowOffsets"]) == len(wrap["rowCanvasBands"]) == 3
    assert "rowPitch" not in wrap
    assert "wrap" not in _geometry(_build(wrap=None))


# The row length sets the size class, and one row length serves every record.


def test_the_size_class_follows_the_row_length() -> None:
    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    whole = DrawnExtent("linear", (300_000,))
    rows = DrawnExtent("linear", (300_000,), wrap_row_length=50_000)
    assert whole.size_class(cfg) == "long"
    assert rows.size_class(cfg) == "short"
    # The window still reads the whole record.
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


def test_a_record_shorter_than_the_row_draws_exactly_as_unwrapped() -> None:
    wrapped = _build(wrap=100_000, overrides={"canvas.linear.track_layout": "above", "objects.scale.style": "ruler", "canvas.linear.ruler_on_axis": True})
    plain = _build(wrap=None, overrides={"canvas.linear.track_layout": "above", "objects.scale.style": "ruler", "canvas.linear.ruler_on_axis": True})
    assert wrapped.drawing.tostring() == plain.drawing.tostring()
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
