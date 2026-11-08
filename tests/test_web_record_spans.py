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
from gbdraw.api.request_render import resolve_request_records
from gbdraw.api.requests import DiagramRequest
from gbdraw.exceptions import ValidationError
from gbdraw.io.regions import RecordWindow, parse_region_spec
from gbdraw.render.track_slot_metadata import RecordAxisInput
from gbdraw.session_request_codec import decode_canonical_request, encode_canonical_request
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
# source coordinates (the crop start a crop records, ``_record_source_start``, for as many bases as the record has),
# and adds the origin window of a record a span may cross the origin of (``spans_may_cross_origin``, the axis source
# map's ``isCircular``; OV-2571, D-52).
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


def _drawn_mito(mode: str, **record: object) -> DiagramRequest:
    item = RecordInput(GenBankInputSource(_MITO), record_key="mito", **record)  # type: ignore[arg-type]
    return CircularDiagramRequest(records=[item]) if mode == "circular" else LinearDiagramRequest(records=(item,))


_RC = {"presentation": RecordPresentation(reverse_complement=True)}
_START = {"display": RecordDisplayOptions(start_coordinate=16000)}


@pytest.mark.parametrize(
    ("mode", "record", "span", "window"),
    [
        ("circular", _START, (1, 16569), True),
        ("circular", _RC, (1, 16569), True),
        ("circular", {**_RC, "display": RecordDisplayOptions(start_coordinate=500)}, (1, 16569), True),
        ("circular", {"region": parse_region_spec("#1:1001-9000")}, (1001, 9000), False),
        ("circular", {"region": parse_region_spec("#1:1001-9000:rc")}, (1001, 9000), False),
        ("circular", {"region": parse_region_spec("#1:1001-9000"), **_RC}, (1001, 9000), False),
        ("linear", {}, (1, 16569), True),
        ("linear", {"region": parse_region_spec("#1:1001-9000")}, (1001, 9000), False),
        ("linear", {"region": parse_region_spec("#1:1001-9000:rc")}, (1001, 9000), False),
        ("linear", {"display": RecordDisplayOptions(is_circular=False)}, (1, 16569), False),
    ],
    ids=[
        "display-start", "rc", "rc-display-start", "crop", "crop-rc", "crop-then-rc",
        "linear", "linear-crop", "linear-crop-rc", "linear-forced-linear",
    ],
)
def test_a_drawn_record_reads_the_span_and_window_its_axis_source_map_shows(
    tmp_path: Path, mode: str, record: dict[str, object], span: tuple[int, int], window: bool,
) -> None:
    request = _drawn_mito(mode, **record)
    reply = _read_drawn(request, tmp_path, "mito")
    start, end = span
    assert _spans(reply)[0] == ("mito", start, end, False)
    assert reply["spans"][0]["sequence"] == _BASES[start - 1:end]
    assert [item["originWindow"] for item in reply["spans"][1:]] == ([True] if window else [])
    # The same record's axis source map, which Select region reads (``sourceSpanOf`` over the whole drawn span):
    # display indexes 0 and drawnLength - 1 show the two ends, or the record is drawn whole; isCircular decides the window.
    payload, paths = _committed(request, tmp_path / "axis")
    collection = resolve_request_records(decode_canonical_request(payload, resource_paths=paths, output_directory=str(tmp_path)))
    source = RecordAxisInput(collection.records[0], "mito", collection.transforms[0]).source()
    shown = [(source["anchorBase"] - 1 + source["step"] * index) % source["sourceLength"] + 1
             for index in (0, source["drawnLength"] - 1)]
    whole = source["drawnLength"] == source["sourceLength"]
    assert (start, end) == ((1, source["sourceLength"]) if whole else (min(shown), max(shown)))
    assert source["isCircular"] is window


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
        "queries": [{"name": "q", "length": 100, "program": "blastn"}], "colorMode": "query", "queryColors": {},
    })
    assert [(hit["recordKey"], hit["start"]) for hit in search["hits"]] == [("a", 1001), ("b", 2001)]


@pytest.mark.parametrize(
    "spans",
    [
        [{"recordKey": "mito", "originWindow": 0}],
        [{"recordKey": "mito", "originWindow": -5}],
        [{"recordKey": "mito", "originWindow": 600}, {"recordKey": "mito", "start": 1, "end": 10}],
    ],
    ids=["window-0", "window-negative", "mixed-kinds"],
)
def test_a_read_with_a_window_below_one_or_both_kinds_is_a_protocol_error(tmp_path: Path, spans) -> None:
    payload, paths = _committed(CircularDiagramRequest(records=[RecordInput(GenBankInputSource(_MITO), record_key="mito")]), tmp_path)
    with pytest.raises(ValidationError) as raised:
        read_record_span_sequences(payload, paths, str(tmp_path / "out"), spans)
    assert raised.value.diagnostic["code"] == "HELPER_PROTOCOL"


# A record drawn through windows (D-01) is named by its own key: a span reads its source, a drawn read reads each
# window, and every reply names the record.
def _windows_card(*windows) -> LinearDiagramRequest:

    return LinearDiagramRequest(records=(
        RecordInput(GenBankInputSource(_MITO), record_key="card", windows=tuple(RecordWindow(*w) for w in windows)),
    ))


def test_a_span_of_a_record_drawn_in_windows_is_read_by_the_record_key(tmp_path: Path) -> None:
    [span] = _read(_windows_card((1001, 3000), (6001, 8000, True)), tmp_path, ("card", 2001, 2010))["spans"]
    assert (span["recordKey"], span["header"], span["sequence"]) == ("card", "NC_012920.1:2001-2010", _BASES[2000:2010])


def test_a_drawn_read_of_a_record_in_windows_reads_every_window(tmp_path: Path) -> None:
    reply = _read_drawn(_windows_card((1001, 3000), (6001, 8000, True)), tmp_path, "card")
    assert _spans(reply) == [("card", 1001, 3000, False), ("card", 6001, 8000, False)]
    assert [span["sequence"] for span in reply["spans"]] == [_BASES[1000:3000], _BASES[6000:8000]]


def test_a_drawn_read_of_a_window_across_the_origin_wraps(tmp_path: Path) -> None:
    request = LinearDiagramRequest(records=(
        RecordInput(
            GenBankInputSource(_MITO), record_key="card", display=RecordDisplayOptions(is_circular=True),
            windows=(RecordWindow(16001, 500),),
        ),
    ))
    [span, *_origin] = _read_drawn(request, tmp_path, "card")["spans"]
    assert (span["recordKey"], span["start"], span["end"], span["wrapsOrigin"]) == ("card", 16001, 500, True)
    assert span["sequence"] == _BASES[16000:] + _BASES[:500]
