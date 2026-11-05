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
        "wrapsOrigin": False, "header": "NC_012920.1:4001-4100", "sequence": _BASES[4000:4100],
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
