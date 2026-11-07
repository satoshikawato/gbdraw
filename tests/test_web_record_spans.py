"""The Worker helper behind Copy FASTA of Select region (PD-OI-115).

The page sends the request the shown Result was made from; the helper resolves
its records as the render does and slices the source sequence of a record by
source coordinates, on the source strand, so a reverse-complemented display or
a crop gives the bases of the input file.
"""
from __future__ import annotations

from pathlib import Path

import pytest
from Bio import SeqIO

from gbdraw.api import (
    CircularDiagramRequest,
    GenBankInputSource,
    LinearDiagramRequest,
    RecordDisplayOptions,
    RecordInput,
    RecordPresentation,
    resolve_request,
)
from gbdraw.api.requests import DiagramRequest
from gbdraw.exceptions import ValidationError
from gbdraw.io.regions import parse_region_spec
from gbdraw.session_request_codec import encode_canonical_request
from gbdraw.web_support.record_spans import read_record_span_sequences

_INPUTS = Path(__file__).resolve().parent / "test_inputs"
_MITO = _INPUTS / "HmmtDNA.gbk"
_BASES = str(SeqIO.read(_MITO, "genbank").seq)


def _committed(request: DiagramRequest, directory: Path) -> tuple[dict[str, object], dict[str, str]]:
    """The request as a Generate commits it: planned, encoded, its resources staged."""

    encoded = encode_canonical_request(resolve_request(request))
    staged = directory / "resources"
    staged.mkdir(parents=True, exist_ok=True)
    paths = {}
    for index, resource in enumerate(encoded.resources, start=1):
        target = staged / f"{index:04d}.bin"
        target.write_bytes(
            resource.content if resource.content is not None else Path(str(resource.source_path)).read_bytes()
        )
        paths[resource.resource_id] = str(target)
    return encoded.payload, paths


def _read(request: DiagramRequest, directory: Path, *spans: tuple[str, int, int]):
    payload, paths = _committed(request, directory)
    return read_record_span_sequences(
        payload, paths, str(directory / "out"),
        [{"recordKey": key, "start": start, "end": end} for key, start, end in spans],
    )


def test_a_span_of_a_forward_record_reads_its_source_bases_with_a_source_header(tmp_path: Path) -> None:
    request = CircularDiagramRequest(records=[RecordInput(GenBankInputSource(_MITO), record_key="mito")])
    reply = _read(request, tmp_path, ("mito", 4001, 4100))
    [span] = reply["spans"]
    assert span == {
        "recordKey": "mito", "recordId": "NC_012920.1", "start": 4001, "end": 4100, "length": 100,
        "wrapsOrigin": False, "originWindow": False, "recordLength": 16569, "header": "NC_012920.1:4001-4100",
        "sequence": _BASES[4000:4100],
    }
    # 60 columns, uppercase: the one FASTA writer of the package (losat_fasta_text).
    lines = reply["fasta"].splitlines()
    assert lines[0] == ">NC_012920.1:4001-4100"
    assert [len(line) for line in lines[1:]] == [60, 40]
    assert "".join(lines[1:]) == _BASES[4000:4100].upper()


def test_a_reverse_complemented_display_reads_the_source_strand(tmp_path: Path) -> None:
    forward = LinearDiagramRequest(records=(RecordInput(GenBankInputSource(_MITO), record_key="card"),))
    reversed_ = LinearDiagramRequest(records=(
        RecordInput(GenBankInputSource(_MITO), record_key="card", presentation=RecordPresentation(reverse_complement=True)),
    ))
    assert _read(reversed_, tmp_path / "rc", ("card", 101, 400)) == _read(forward, tmp_path / "fw", ("card", 101, 400))


def test_a_crop_reads_source_coordinates_of_the_whole_record(tmp_path: Path) -> None:
    cropped = LinearDiagramRequest(records=(
        RecordInput(GenBankInputSource(_MITO), record_key="card", region=parse_region_spec("#1:1001-9000:rc")),
    ))
    [span] = _read(cropped, tmp_path, ("card", 2001, 2010))["spans"]
    assert (span["header"], span["sequence"]) == ("NC_012920.1:2001-2010", _BASES[2000:2010])


def test_a_span_across_the_origin_reads_the_end_then_the_start(tmp_path: Path) -> None:
    request = CircularDiagramRequest(records=[RecordInput(GenBankInputSource(_MITO), record_key="mito")])
    [span] = _read(request, tmp_path, ("mito", 16001, 501))["spans"]
    assert span["wrapsOrigin"] is True
    assert span["header"] == "NC_012920.1:16001-501"
    assert span["sequence"] == _BASES[16000:] + _BASES[:501]
    assert span["length"] == 569 + 501


@pytest.mark.parametrize(
    ("span", "code"),
    [(("other", 1, 10), "INPUT_INVALID"), (("mito", 0, 10), "REGION_INVALID"), (("mito", 1, 16570), "REGION_INVALID")],
)
def test_a_span_outside_the_drawn_records_is_refused_with_a_typed_reason(tmp_path: Path, span, code: str) -> None:
    request = CircularDiagramRequest(records=[RecordInput(GenBankInputSource(_MITO), record_key="mito")])
    with pytest.raises(ValidationError) as raised:
        _read(request, tmp_path, span)
    assert raised.value.diagnostic["code"] == code


def test_a_span_across_the_end_of_a_linear_record_is_refused(tmp_path: Path) -> None:
    request = LinearDiagramRequest(records=(
        RecordInput(GenBankInputSource(_INPUTS.parent.parent / "examples" / "MjeNMV.gb"), record_key="card"),
    ))
    with pytest.raises(ValidationError) as raised:
        _read(request, tmp_path, ("card", 500, 100))
    assert raised.value.diagnostic["reason"] == "RECORD_BOUNDS"


_MJENMV = _INPUTS.parent.parent / "examples" / "MjeNMV.gb"


def _read_drawn(request: DiagramRequest, directory: Path, *keys: str, reach: int = 600):
    payload, paths = _committed(request, directory)
    return read_record_span_sequences(
        payload, paths, str(directory / "out"), [{"recordKey": key, "originWindow": reach} for key in keys]
    )


def _spans(reply) -> list[tuple[str, int, int, bool]]:
    return [(span["recordKey"], span["start"], span["end"], span["originWindow"]) for span in reply["spans"]]


# A sequence search asks for the records the shown drawing draws by key alone: Python reads the bases each draws, in
# source coordinates, from the transform the record is drawn by (the one the axis source map is written from), and
# adds the origin window of a record a span may pass the origin of (circular, no region, drawn whole; OV-2571, D-52).
@pytest.mark.parametrize(("reach", "width"), [(600, 600), (20000, 16569 // 2)])
def test_a_drawn_circular_record_reads_whole_then_its_origin_window(tmp_path: Path, reach: int, width: int) -> None:
    request = CircularDiagramRequest(records=[RecordInput(GenBankInputSource(_MITO), record_key="mito")])
    reply = _read_drawn(request, tmp_path, "mito", reach=reach)
    start = 16569 - width + 1
    assert _spans(reply) == [("mito", 1, 16569, False), ("mito", start, width, True)]
    whole, window = reply["spans"]
    assert whole["sequence"] == _BASES and whole["wrapsOrigin"] is False
    assert (window["wrapsOrigin"], window["recordLength"]) == (True, 16569)
    assert window["header"] == f"NC_012920.1:{start}-{width}"
    assert window["sequence"] == _BASES[start - 1:] + _BASES[:width]


def test_a_drawn_record_with_a_display_start_reads_whole_with_its_origin_window(tmp_path: Path) -> None:
    request = CircularDiagramRequest(records=[RecordInput(
        GenBankInputSource(_MITO), record_key="mito", display=RecordDisplayOptions(start_coordinate=16000),
    )])
    assert _spans(_read_drawn(request, tmp_path, "mito")) == [("mito", 1, 16569, False), ("mito", 15970, 600, True)]


def test_a_reverse_complemented_drawn_record_reads_the_whole_source_strand(tmp_path: Path) -> None:
    request = CircularDiagramRequest(records=[RecordInput(
        GenBankInputSource(_MITO), record_key="mito", presentation=RecordPresentation(reverse_complement=True),
    )])
    reply = _read_drawn(request, tmp_path, "mito")
    assert _spans(reply) == [("mito", 1, 16569, False), ("mito", 15970, 600, True)]
    assert reply["spans"][0]["sequence"] == _BASES


@pytest.mark.parametrize("spec", ["#1:1001-9000", "#1:1001-9000:rc"])
def test_a_drawn_crop_reads_the_bases_it_draws_and_has_no_origin_window(tmp_path: Path, spec: str) -> None:
    request = CircularDiagramRequest(records=[
        RecordInput(GenBankInputSource(_MITO), record_key="mito", region=parse_region_spec(spec)),
    ])
    reply = _read_drawn(request, tmp_path, "mito")
    assert _spans(reply) == [("mito", 1001, 9000, False)]
    assert reply["spans"][0]["sequence"] == _BASES[1000:9000]


def test_a_drawn_linear_record_has_no_origin_window(tmp_path: Path) -> None:
    request = LinearDiagramRequest(records=(RecordInput(GenBankInputSource(_MJENMV), record_key="card"),))
    assert _spans(_read_drawn(request, tmp_path, "card")) == [("card", 1, 306008, False)]


def test_a_record_drawn_twice_is_read_once_and_windows_follow_every_drawn_span(tmp_path: Path) -> None:
    request = LinearDiagramRequest(records=(
        RecordInput(GenBankInputSource(_MITO), record_key="mito"), RecordInput(GenBankInputSource(_MJENMV), record_key="card"),
    ))
    assert _spans(_read_drawn(request, tmp_path, "mito", "card", "mito", reach=30)) == [
        ("mito", 1, 16569, False), ("card", 1, 306008, False), ("mito", 16540, 30, True),
    ]


# Two drawn records with the same ID and length (two files that hold one record each): each searched sequence keeps
# a header of its own, which names its record key, so LOSAT's rows map back to the record they hit.
def test_spans_with_one_header_name_their_record_keys_and_map_back_to_their_records(tmp_path: Path) -> None:
    from gbdraw.web_support.sequence_search import sequence_search_annotation_set

    request = LinearDiagramRequest(records=(
        RecordInput(GenBankInputSource(_MITO), record_key="a"), RecordInput(GenBankInputSource(_MITO), record_key="b"),
    ))
    reply = _read(request, tmp_path, ("a", 1, 16569), ("b", 1, 16569), ("a", 101, 200))
    headers = [span["header"] for span in reply["spans"]]
    assert headers == ["NC_012920.1:1-16569@a", "NC_012920.1:1-16569@b", "NC_012920.1:101-200"]
    assert [line[1:] for line in reply["fasta"].splitlines() if line.startswith(">")] == headers
    rows = "".join(
        f"q\t{header}\t100.000\t100\t0\t0\t1\t100\t{start}\t{start + 99}\t1e-50\t180\n"
        for header, start in ((headers[0], 1001), (headers[1], 2001))
    )
    search = sequence_search_annotation_set({
        "table": rows, "setId": "Search 1",
        "subjects": [
            {"subjectId": span["header"], "record": {"kind": "recordIndex", "index": index},
             "recordKey": span["recordKey"], "sourceStart": span["start"], "originWindow": span["originWindow"],
             "recordLength": span["recordLength"]}
            for index, span in enumerate(reply["spans"][:2])
        ],
        "thresholds": {"evalue": 1e-2, "bitscore": 50, "identity": 70, "alignmentLength": 0},
        "program": "blastn", "queries": [{"name": "q", "length": 100}], "colorMode": "query", "queryColors": {},
    })
    assert [(hit["recordKey"], hit["start"]) for hit in search["hits"]] == [("a", 1001), ("b", 2001)]
