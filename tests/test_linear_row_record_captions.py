"""The record captions of a Linear row of whole records (OV-3521, OV-3727).

In a row of several whole records each record's accession and length stand above it. Where two of them would overlap,
the row leaves out the caption of the shorter record (the longer record keeps its caption first; two kept captions stay
half their font size apart, closer they read as one) and says once how many it left out: a CLI WARNING, and the same text in the Web's input summary. Every caption stands at least above one lane of
features on each strand a record of the row draws, so the caption of a record with no drawn feature, or with features only on the lower strand,
is level with its neighbours'; a row with nothing to leave out otherwise keeps the captions where they were.
"""

from __future__ import annotations

import logging
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api import (
    GenBankInputSource,
    LinearDiagramRequest,
    LinearMultiRecordOptions,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.core.text import measure_mixed_text_line
from gbdraw.session import build_session_document
from gbdraw.web_support.request_render import render_embedded_canonical_web_request
from tests.utils.session_files import embedded_session

SVG = "{http://www.w3.org/2000/svg}"
LONG = 100_000
SHORT = 400
SHORT_COUNT = 8


def _record(record_id: str, length: int, *, features: bool, strand: int = 1) -> SeqRecord:
    record = SeqRecord(Seq("ATGC" * (length // 4)), id=record_id, name=record_id, description=f"{record_id} test")
    record.annotations["molecule_type"] = "DNA"
    if features:
        step = max(length // 4, 120)
        for start in range(0, length - 99, step):
            record.features.append(
                SeqFeature(FeatureLocation(start, start + 99, strand=strand), type="CDS", qualifiers={"product": ["p"]})
            )
    return record


def _files(tmp_path: Path, records: list[SeqRecord]) -> list[Path]:
    paths = []
    for record in records:
        path = tmp_path / f"{record.id}.gb"
        SeqIO.write([record], path, "genbank")
        paths.append(path)
    return paths


def _one_long_and_many_short(tmp_path: Path) -> list[Path]:
    """One long record and eight short ones (every other one with no feature), as one row."""

    records = [_record("LONG0001.1", LONG, features=True)] + [
        _record(f"SHORT{index:03d}.1", SHORT, features=index % 2 == 0) for index in range(1, SHORT_COUNT + 1)
    ]
    return _files(tmp_path, records)


def _draw_one_row(tmp_path: Path, paths: list[Path], *extra: str) -> ET.Element:
    from gbdraw.linear import linear_main

    positions = [token for index in range(1, len(paths) + 1) for token in ("--multi_record_position", f"#{index}@1")]
    linear_main(["--gbk", *map(str, paths), *positions, *extra, "-o", str(tmp_path / "out"), "-f", "svg"])
    return ET.parse(tmp_path / "out.svg").getroot()


def _translate(element: ET.Element) -> tuple[float, float]:
    found = re.findall(r"translate\(\s*([-+0-9.eE]+)[,\s]+([-+0-9.eE]+)\s*\)", element.attrib.get("transform", ""))
    assert found
    return sum(float(x) for x, _y in found), sum(float(y) for _x, y in found)


def _captions(root: ET.Element) -> dict[str, tuple[float, float, float, float]]:
    """Each drawn record caption by record ID: its left, right, line y and font size, from its text as the SVG
    draws it."""

    boxes = {}
    for group in root.iter(f"{SVG}g"):
        if group.attrib.get("data-gbdraw-role") != "record-definition":
            continue
        texts = [text for text in group.iter(f"{SVG}text") if "".join(text.itertext()).strip()]
        if not texts:
            continue
        x, y = _translate(group)
        width = max(
            measure_mixed_text_line(
                [("".join(text.itertext()), "normal")],
                text.attrib["font-family"],
                float(text.attrib["font-size"]),
                font_weight=text.attrib.get("font-weight", "normal"),
            )[0]
            for text in texts
        )
        assert all(text.attrib.get("text-anchor") == "middle" for text in texts)
        size = max(float(text.attrib["font-size"]) for text in texts)
        boxes[group.attrib["data-gbdraw-record-id"]] = (x - 0.5 * width, x + 0.5 * width, y, size)
    return boxes


def _warnings(caplog: pytest.LogCaptureFixture) -> list[str]:
    return [record.getMessage() for record in caplog.records if record.levelno == logging.WARNING]


def test_a_row_leaves_out_captions_that_would_overlap_and_says_how_many(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    caplog.set_level(logging.WARNING)
    root = _draw_one_row(tmp_path, _one_long_and_many_short(tmp_path))

    captions = _captions(root)
    assert "LONG0001.1" in captions  # the longest record keeps its caption first
    drawn = sorted(captions.values())
    for (_left, right, _y, size), (next_left, _right, _next_y, _size) in zip(drawn, drawn[1:]):
        assert next_left - right >= 0.5 * size
    left_out = SHORT_COUNT + 1 - len(captions)
    assert 0 < left_out < SHORT_COUNT
    message = (
        f"{left_out} of {SHORT_COUNT + 1} record captions in row 1 were left out because they overlap; "
        "to name the row once instead, give it a name (record_label) and hide Accession and Length."
    )
    assert _warnings(caplog).count(message) == 1


@pytest.mark.parametrize(
    ("layout", "height"),
    [
        pytest.param((), "189.03125px", id="middle"),
        # Review r3: above the axis an empty features slot is packed at the axis band's top edge, not at the axis.
        pytest.param(("--track_layout", "above"), "201.03125px", id="above"),
    ],
)
def test_a_record_with_no_feature_has_its_caption_where_its_neighbours_have_theirs(
    tmp_path: Path, layout: tuple[str, ...], height: str
) -> None:
    """OV-3727: a record with no drawn feature has no feature lanes; its caption stands above the lanes it would draw,
    measured from its axis as its neighbours' lanes are, as their captions stand above theirs. The figure is as high as
    on dev."""

    paths = _files(tmp_path, [
        _record("FEAT0001.1", 20_000, features=True),
        _record("EMPTY001.1", 20_000, features=False),
        _record("FEAT0002.1", 20_000, features=True),
    ])
    root = _draw_one_row(tmp_path, paths, *layout)
    captions = _captions(root)

    assert set(captions) == {"FEAT0001.1", "EMPTY001.1", "FEAT0002.1"}
    assert len({round(y, 6) for _left, _right, y, _size in captions.values()}) == 1
    assert root.attrib["height"] == height


@pytest.mark.parametrize(
    ("neighbour", "height"),
    [
        # No record of the row draws on the upper strand: the captions stand where dev has them, and so does the row.
        pytest.param(("EMPTY001.1", False), "197.03125px", id="no-feature"),
        pytest.param(("PLUS0001.1", True), "245.03125px", id="upper-strand"),
    ],
)
def test_a_record_with_features_only_on_the_lower_strand_has_its_caption_level_with_its_neighbours(
    tmp_path: Path, neighbour: tuple[str, bool], height: str
) -> None:
    """OV-3727 with strands drawn apart (the Web's Linear default): a record whose features are all on the lower strand
    draws no upper lane; where another record of its row draws on the upper strand, its caption stands above the upper
    lane it would draw, level with its neighbour's caption. The figure is as high as on dev in both rows."""

    neighbour_id, neighbour_features = neighbour
    paths = _files(tmp_path, [
        _record(neighbour_id, 20_000, features=neighbour_features),
        _record("MINUS001.1", 20_000, features=True, strand=-1),
    ])
    root = _draw_one_row(tmp_path, paths, "--separate_strands")
    captions = _captions(root)

    assert set(captions) == {neighbour_id, "MINUS001.1"}
    assert len({round(y, 6) for _left, _right, y, _size in captions.values()}) == 1
    assert root.attrib["height"] == height


def _labelled_row(tmp_path: Path, *, with_labels: bool) -> dict[str, tuple[float, float, float, float]]:
    """Captions, as offsets from their record's axis, of a row whose middle record has labelled features (or none)."""

    def record(record_id: str, labelled: bool) -> SeqRecord:
        record = SeqRecord(Seq("ATGC" * 5000), id=record_id, name=record_id, description=f"{record_id} test")
        record.annotations["molecule_type"] = "DNA"
        for index, start in enumerate(range(0, 19000, 800)):
            qualifiers = {"gene": [f"longgenename{index:02d}"], "product": ["a fairly long product name"]} if labelled else {}
            record.features.append(SeqFeature(FeatureLocation(start, start + 600, strand=1), type="CDS", qualifiers=qualifiers))
        return record

    folder = tmp_path / ("labelled" if with_labels else "plain")
    folder.mkdir()
    records = [record("PLAIN001.1", False)]
    if with_labels:
        records.append(record("LABEL001.1", True))
    records.append(record("PLAIN002.1", False))
    root = _draw_one_row(folder, _files(folder, records), "--show_labels", "all")
    axes = {
        group.attrib["data-gbdraw-record-id"]: _translate(group)[1]
        for group in root.iter(f"{SVG}g")
        if group.attrib.get("data-gbdraw-slot-id") == "features" and "data-gbdraw-record-id" in group.attrib
    }
    return {
        record_id: (left, right, y - axes[record_id], size)
        for record_id, (left, right, y, size) in _captions(root).items()
    }


def test_a_caption_stays_above_its_own_record_beside_a_record_with_labels(tmp_path: Path) -> None:
    """Review r1: feature labels above one record do not lift the captions of the other records of its row."""

    beside_labels = _labelled_row(tmp_path, with_labels=True)
    alone = _labelled_row(tmp_path, with_labels=False)

    assert set(beside_labels) == {"PLAIN001.1", "LABEL001.1", "PLAIN002.1"}
    for record_id in ("PLAIN001.1", "PLAIN002.1"):
        assert beside_labels[record_id][2] == pytest.approx(alone[record_id][2])
    assert beside_labels["LABEL001.1"][2] < beside_labels["PLAIN001.1"][2] - 100


def test_a_row_whose_captions_fit_keeps_every_caption_and_warns_nothing(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    caplog.set_level(logging.WARNING)
    paths = _files(tmp_path, [_record(f"FIT{index:04d}.1", 20_000, features=True) for index in range(1, 4)])
    captions = _captions(_draw_one_row(tmp_path, paths))

    assert len(captions) == 3
    assert not [message for message in _warnings(caplog) if "record captions" in message]


def test_the_web_reply_carries_the_warning_the_cli_logs(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """CLI = Web: the Web's input summary says the same, from the same layout owner."""

    caplog.set_level(logging.WARNING)
    paths = _one_long_and_many_short(tmp_path)
    _draw_one_row(tmp_path, paths)
    cli_messages = [message for message in _warnings(caplog) if "record captions" in message]
    assert len(cli_messages) == 1

    request = LinearDiagramRequest(
        records=tuple(RecordInput(source=GenBankInputSource(path)) for path in paths),
        layout=LinearMultiRecordOptions(multi_record_positions=tuple(f"#{index}@1" for index in range(1, len(paths) + 1))),
        output=RenderOutputRequest(output_prefix="row", formats=("svg",)),
    )
    document = build_session_document(request)
    embedded = embedded_session(document.to_dict(), document.resources.sources)
    reply = render_embedded_canonical_web_request(
        embedded["drawings"][0]["renderRequest"], resources=embedded["resources"], workspace=tmp_path / "web",
    )

    warnings = reply["metadata"]["inputSummaries"][0]["warnings"]
    assert [item for item in warnings if item["code"] == "record_captions_overlap"] == [
        {"code": "record_captions_overlap", "message": cli_messages[0]}
    ]
