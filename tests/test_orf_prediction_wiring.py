"""ORF prediction through the one record loading path: the CLI options, the typed request field, its canonical
form, and the Python API (D-87: off unless asked; D-126 defaults). Predicted CDS are added only to records without
CDS of their own until their placement beside annotated CDS is decided."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest
from Bio import SeqIO

import gbdraw
from gbdraw.analysis.orf_prediction import PREDICTED_ORF_CAPTION, PREDICTED_ORF_COLOR, OrfPredictionParameters
from gbdraw.api.options import LinearDiagramOptions
from gbdraw.api.requests import FastaInputSource, LinearDiagramRequest, RecordInput
from gbdraw.exceptions import ValidationError
from gbdraw.session_request_codec import (
    CanonicalRequestDecodingError,
    decode_canonical_request,
    encode_canonical_request,
)
from tests.utils.cli_in_process import run_gbdraw

INPUTS = Path(__file__).resolve().parent / "test_inputs"


def _fasta(tmp_path: Path, name: str) -> Path:
    """A FASTA copy of a test record: its sequence without features."""

    record = SeqIO.read(INPUTS / name, "genbank")
    path = tmp_path / f"{Path(name).stem}.fasta"
    path.write_text(f">{record.id}\n{record.seq}\n")
    return path


def _draw(tmp_path: Path, mode: str, *args: str) -> tuple[int, str, str]:
    result = run_gbdraw([mode, *args, "-o", str(tmp_path / "out"), "-f", "svg", "--overwrite"], cwd=tmp_path)
    svg = tmp_path / "out.svg"
    return result.returncode, result.stdout + result.stderr, svg.read_text() if svg.exists() else ""


def test_prediction_is_off_unless_asked(tmp_path: Path) -> None:
    code, output, svg = _draw(tmp_path, "linear", "--fasta", str(_fasta(tmp_path, "NC_001416.gb")))
    assert code == 0, output
    assert "drew 0 features" in output
    assert PREDICTED_ORF_CAPTION not in svg


@pytest.mark.parametrize("mode", ["linear", "circular"])
def test_a_fasta_record_draws_its_predicted_cds_in_their_own_legend_row(tmp_path: Path, mode: str) -> None:
    code, output, svg = _draw(tmp_path, mode, "--fasta", str(_fasta(tmp_path, "NC_001416.gb")), "--predict_orfs")
    assert code == 0, output
    # D-126 defaults: no CDS states a table, so table 11; methionine starts; 50 residues (as orfind with -m 50).
    assert "drew 110 features (110 CDS)" in output
    assert PREDICTED_ORF_CAPTION in svg
    assert PREDICTED_ORF_COLOR in svg.lower()


def test_the_cli_options_reach_the_prediction(tmp_path: Path) -> None:
    code, output, _ = _draw(
        tmp_path, "linear", "--fasta", str(_fasta(tmp_path, "NC_001416.gb")), "--predict_orfs",
        "--orf_min_length", "100", "--orf_genetic_code", "1", "--orf_start_codons", "table", "--orf_keep_nested",
    )
    assert code == 0, output
    expected = gbdraw.api.OrfPredictionParameters(
        genetic_code=1, start_codons="table", min_protein_length=100, keep_nested=True
    )
    record = SeqIO.read(INPUTS / "NC_001416.gb", "genbank")
    count = len(gbdraw.analysis.orf_prediction.predict_orfs(record.seq, topology="linear", parameters=expected))
    assert f"drew {count} features ({count} CDS)" in output


def test_the_record_is_read_in_the_topology_it_is_drawn_in(tmp_path: Path) -> None:
    fasta = str(_fasta(tmp_path, "NZ_CP046414.1.gbk"))
    code, output, _ = _draw(tmp_path, "circular", "--fasta", fasta, "--predict_orfs")
    assert code == 0, output
    assert "drew 119 features (119 CDS)" in output
    # Read as circular, an ORF runs across the origin and holds a minus-strand ORF over it (prh HANDOFF).
    code, output, _ = _draw(tmp_path, "circular", "--fasta", fasta, "--predict_orfs", "--record_topology", "circular")
    assert code == 0, output
    assert "drew 118 features (118 CDS)" in output


def test_a_record_with_cds_of_its_own_is_refused_for_now(tmp_path: Path) -> None:
    code, output, _ = _draw(tmp_path, "linear", "--gbk", str(INPUTS / "NC_001416.gb"), "--predict_orfs")
    assert code != 0
    assert "ORF prediction adds CDS only to a record without CDS of its own" in output
    assert "has 73" in output


@pytest.mark.parametrize(
    ("args", "message"),
    [
        (["--orf_min_length", "60"], "--orf_min_length needs --predict_orfs"),
        (["--predict_orfs", "--orf_genetic_code", "7"], "Unknown genetic code 7"),
        (["--predict_orfs", "--orf_min_length", "0"], "at least 1"),
        (["--predict_orfs", "--orf_start_codons", "all"], "invalid choice"),
    ],
)
def test_orf_options_are_checked(tmp_path: Path, args: list[str], message: str) -> None:
    code, output, _ = _draw(tmp_path, "linear", "--fasta", str(_fasta(tmp_path, "NC_001416.gb")), *args)
    assert code != 0
    assert message in output


def _request(tmp_path: Path, orf_prediction: OrfPredictionParameters | None) -> LinearDiagramRequest:
    return LinearDiagramRequest(
        records=(RecordInput(source=FastaInputSource(_fasta(tmp_path, "NC_001416.gb")), record_key="lambda"),),
        options=LinearDiagramOptions(orf_prediction=orf_prediction),
    )


def _diagram_options(payload: object) -> dict[str, object]:
    assert isinstance(payload, dict)
    return payload["diagramOptions"]


@pytest.mark.parametrize(
    "value",
    [
        {"geneticCode": 7, "startCodons": "methionine", "minProteinLength": 50, "keepNested": False},
        {"geneticCode": "records", "startCodons": "methionine", "minProteinLength": 50, "keepNested": False},
        {"geneticCode": "record", "startCodons": "methionine", "minProteinLength": 0, "keepNested": False},
        {"geneticCode": "record", "startCodons": "methionine", "minProteinLength": 50},
    ],
)
def test_a_stored_field_outside_its_domain_is_refused(tmp_path: Path, value: dict[str, object]) -> None:
    payload = encode_canonical_request(_request(tmp_path, OrfPredictionParameters())).payload
    payload = {**payload, "diagramOptions": {**_diagram_options(payload), "orfPrediction": value}}
    resources = {resource.resource_id: resource.source_path for resource in encode_canonical_request(
        _request(tmp_path, None)).resources}
    with pytest.raises(CanonicalRequestDecodingError, match="orfPrediction"):
        decode_canonical_request(payload, resource_paths=resources, output_directory=tmp_path)


def test_the_python_api_option_predicts_on_records_given_in_memory(tmp_path: Path) -> None:
    records = gbdraw.read_fasta(_fasta(tmp_path, "HmmtDNA.gbk"))
    features_before = list(records[0].features)
    options = gbdraw.LinearOptions(
        features=gbdraw.FeatureOptions(orf_prediction=gbdraw.OrfPredictionParameters(genetic_code=2))
    )
    svg = gbdraw.draw_linear(records, options=options).to_svg()
    assert PREDICTED_ORF_CAPTION in svg
    # The caller's record is not changed.
    assert records[0].features == features_before
    with pytest.raises(ValidationError, match="orf_prediction"):
        replace(LinearDiagramOptions(), orf_prediction="record")  # type: ignore[arg-type]
