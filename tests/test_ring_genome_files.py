"""A similarity ring whose comparison genome is split over several files (PR-B).

One strain's genome is often several files (a chromosome and plasmids
downloaded separately). ``--conservation_genome FILE [FILE ...]`` (one
occurrence is one ring), a ``;``-separated ``comparison_sequence`` table cell
and ``CircularDiagramOptions.conservation_sequence_files`` with a list per ring
give one ring all of those files: each file whole, in the listed order, records
concatenated, as if the files were one file.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gbdraw.api import CircularDiagramOptions, read_conservation_table
from gbdraw.circular import circular_main
from gbdraw.comparisons.losat_runtime import resolve_losat_runtime
from gbdraw.exceptions import ExportError, ValidationError
from gbdraw.io.comparison_sequences import comparison_genome, read_comparison_sequence_file
from gbdraw.session_request_codec import (
    decode_canonical_request,
    encode_canonical_request,
)
from tests.utils.session_drawing_view import session_view

REPO_ROOT = Path(__file__).resolve().parents[1]
TUTORIAL = REPO_ROOT / "gbdraw" / "web" / "tutorial-data"
HUMAN = TUTORIAL / "human-mitochondrion" / "HmmtDNA.gbk"
FOUR = TUTORIAL / "metazoan-mitochondria-four"
DANIO = FOUR / "NC_002333.2.gb"
DROSOPHILA = FOUR / "NC_024511.2.gb"
DANIO_TSV = TUTORIAL / "metazoan-mitochondria-comparison" / "danio-human.tlosatx.tsv"
DANIO_FASTA = TUTORIAL / "metazoan-mitochondria-comparison" / "NC_002333.2.fna"


def _one_file_genome(tmp_path: Path, *parts: Path) -> Path:
    """``cat parts > <dir>/<first part's name>``: the same genome as one file of the same name."""

    joined = tmp_path / "joined" / parts[0].name
    joined.parent.mkdir(parents=True, exist_ok=True)
    joined.write_bytes(b"".join(part.read_bytes() for part in parts))
    return joined


def _svg(prefix: Path, suffix: str = ".svg") -> bytes:
    return prefix.with_suffix(suffix).read_bytes()


@pytest.fixture(scope="module")
def native_losat() -> None:
    try:
        runtime = resolve_losat_runtime("tlosatx")
    except ValidationError as error:
        pytest.skip(f"No native LOSAT runtime resolves: {error}")
    if runtime.kind != "losat":
        pytest.skip(f"Resolved runtime is {runtime.kind}, not native LOSAT.")


_LOSAT = ["--losat", "tlosatx", "--losat_gencode", "2", "--conservation_losat_gencode", "2"]


@pytest.mark.circular
def test_a_two_file_ring_draws_the_figure_of_one_file_with_both_records(
    native_losat, tmp_path: Path
) -> None:
    joined = _one_file_genome(tmp_path, DANIO, DROSOPHILA)
    one = tmp_path / "one"
    circular_main(["--gbk", str(HUMAN), *_LOSAT, "--conservation_sequence", str(joined), "-o", str(one), "-f", "svg"])

    two = tmp_path / "two"
    raw = tmp_path / "raw"
    circular_main([
        "--gbk", str(HUMAN), *_LOSAT, "--conservation_genome", str(DANIO), str(DROSOPHILA),
        "--losat_output_dir", str(raw), "-o", str(two), "-f", "svg",
    ])
    assert _svg(two) == _svg(one)

    # The table cell lists the files of the ring, separated by ';'.
    table = tmp_path / "rings.tsv"
    table.write_text(f"comparison_sequence\tlosat_gencode\n{DANIO};{DROSOPHILA}\t2\n", encoding="utf-8")
    from_table = tmp_path / "from-table"
    circular_main([
        "--gbk", str(HUMAN), "--losat", "tlosatx", "--losat_gencode", "2",
        "--conservation_table", str(table), "-o", str(from_table), "-f", "svg",
    ])
    assert _svg(from_table) == _svg(one)

    # --losat_output_dir writes the ring's files as one ';' cell, and the
    # manifest is a --conservation_table that draws the same figure.
    manifest = read_conservation_table(str(raw / "conservation.tsv"))
    assert manifest.comparison_sequence_files == [(str(DANIO.resolve()), str(DROSOPHILA.resolve()))]
    assert manifest.labels == ["Danio rerio mitochondrion, complete genome"]
    replay = tmp_path / "replay"
    circular_main([
        "--gbk", str(HUMAN), "--conservation_table", str(raw / "conservation.tsv"),
        "--conservation_reference", "subject", "-o", str(replay), "-f", "svg",
    ])
    assert _svg(replay) == _svg(one)


@pytest.mark.circular
def test_a_precomputed_ring_exports_spans_of_every_file_of_its_genome(tmp_path: Path) -> None:
    """Without LOSAT, the files of a ring are its comparison genome for span export."""

    joined = _one_file_genome(tmp_path, DANIO, DROSOPHILA)
    common = ["--gbk", str(HUMAN), "--conservation_blast", str(DANIO_TSV), "--conservation_reference", "subject"]
    # The interactive SVG names its output prefix: one name in two folders.
    one = tmp_path / "one" / "ring"
    circular_main([*common, "--conservation_sequence", str(joined), "-o", str(one), "-f", "interactive_svg"])
    two = tmp_path / "two" / "ring"
    circular_main([
        *common, "--conservation_genome", str(DANIO), str(DROSOPHILA), "-o", str(two), "-f", "interactive_svg",
    ])
    assert _svg(two, ".interactive.svg") == _svg(one, ".interactive.svg")


def test_two_files_of_one_genome_with_one_record_id_are_refused(tmp_path: Path) -> None:
    copy = tmp_path / "copy" / DANIO_FASTA.name
    copy.parent.mkdir()
    copy.write_bytes(DANIO_FASTA.read_bytes())
    with pytest.raises(ValidationError, match="NC_002333.2") as excinfo:
        comparison_genome((read_comparison_sequence_file(DANIO_FASTA), read_comparison_sequence_file(copy)))
    assert excinfo.value.diagnostic == {
        "code": "INPUT_INVALID",
        "reason": "UNIQUE_IDS",
        "field": "comparison_sequence",
        "recordId": "NC_002333.2",
    }
    # One file is read as before (OV-2954 records its repeated IDs).
    assert comparison_genome((read_comparison_sequence_file(DANIO_FASTA),)).label == "NC_002333.2"

    # Without LOSAT the genome is read for the span export of the interactive SVG.
    with pytest.raises(ExportError, match=f"record NC_002333.2 in {DANIO_FASTA} and {copy}"):
        circular_main([
            "--gbk", str(HUMAN), "--conservation_blast", str(DANIO_TSV), "--conservation_reference", "subject",
            "--conservation_genome", str(DANIO_FASTA), str(copy), "-o", str(tmp_path / "dup"),
            "-f", "interactive_svg",
        ])


def test_a_comparison_file_that_repeats_a_record_id_draws_as_before_and_warns(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    """OV-2954 (Owner 2026-11-09 「Q2=C」): hover and export name its later records by the first; gbdraw says so."""

    repeated = tmp_path / "repeated.fna"
    repeated.write_bytes(DANIO_FASTA.read_bytes() * 2)
    common = ["--gbk", str(HUMAN), "--conservation_blast", str(DANIO_TSV), "--conservation_reference", "subject"]
    prefix = tmp_path / "out" / "ring"
    with caplog.at_level("WARNING"):
        circular_main([*common, "--conservation_sequence", str(repeated), "-o", str(prefix), "-f", "interactive_svg"])
    warned = [record.getMessage() for record in caplog.records if record.levelname == "WARNING"]
    assert warned.count(
        "Comparison sequence file repeated.fna: record ID NC_002333.2 names more than one record. The ring draws "
        "every record, but the hover popup and the span export show the first record with that ID."
    ) == 1
    # Without the span export the file is not read, so nothing is named.
    caplog.clear()
    with caplog.at_level("WARNING"):
        circular_main([*common, "--conservation_sequence", str(repeated), "-o", str(tmp_path / "plain" / "ring")])
    assert not [record for record in caplog.records if "repeated.fna" in record.getMessage()]


def test_the_genome_label_is_the_rule_applied_to_its_first_file(tmp_path: Path) -> None:
    genome = comparison_genome((read_comparison_sequence_file(DROSOPHILA), read_comparison_sequence_file(DANIO)))
    assert genome.label == "Drosophila melanogaster mitochondrion, complete genome"
    assert [record.id for record in genome.records] == ["NC_024511.2", "NC_002333.2"]
    fasta_first = comparison_genome((read_comparison_sequence_file(DANIO_FASTA), read_comparison_sequence_file(DROSOPHILA)))
    assert fasta_first.label == "NC_002333.2"


@pytest.mark.parametrize(
    "extra",
    [
        ["--conservation_sequence", str(DANIO)],
        ["--conservation_table", str(DANIO_TSV)],
    ],
    ids=["with --conservation_sequence", "with --conservation_table"],
)
def test_conservation_genome_is_refused_with_the_other_ring_sequence_inputs(
    extra: list[str], tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        circular_main([
            "--gbk", str(HUMAN), "--conservation_blast", str(DANIO_TSV),
            "--conservation_genome", str(DANIO), *extra, "-o", str(tmp_path / "x"),
        ])
    assert excinfo.value.code == 2
    assert "--conservation_genome cannot be combined with" in capsys.readouterr().err


def test_conservation_genome_gives_one_ring_per_occurrence(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit):
        circular_main([
            "--gbk", str(HUMAN), "--conservation_blast", str(DANIO_TSV),
            "--conservation_genome", str(DANIO), "--conservation_genome", str(DROSOPHILA),
            "-o", str(tmp_path / "x"),
        ])
    assert "--conservation_genome must give one genome per --conservation_blast" in capsys.readouterr().err


def test_the_typed_field_takes_a_path_or_a_list_of_paths_per_ring() -> None:
    options = CircularDiagramOptions(
        conservation_blast_files=["a.tsv", "b.tsv", "c.tsv"],
        conservation_sequence_files=["a.gb", ["b1.gb", "b2.gb"], None],
    )
    assert options.conservation_sequence_files == (("a.gb",), ("b1.gb", "b2.gb"), ())
    assert options.conservation_genome_files == (("a.gb",), ("b1.gb", "b2.gb"), ())
    assert CircularDiagramOptions().conservation_genome_files == ()
    assert CircularDiagramOptions(conservation_sequence_files=["", []]).conservation_genome_files == ((), ())
    for bad in (["a.gb", ["b.gb", ""]], ["a.gb", [None]], [3]):
        with pytest.raises(ValidationError, match="conservation_sequence_files"):
            CircularDiagramOptions(conservation_sequence_files=bad)  # type: ignore[arg-type]


def test_a_table_cell_lists_the_files_of_its_ring(tmp_path: Path) -> None:
    table = tmp_path / "rings.tsv"
    table.write_text(
        f"blast\tcomparison_sequence\n{DANIO_TSV}\ta.gb ; sub/b.gb\n{DANIO_TSV}\t\n", encoding="utf-8"
    )
    parsed = read_conservation_table(str(table))
    assert parsed.comparison_sequence_files == [(str(tmp_path / "a.gb"), str(tmp_path / "sub" / "b.gb")), ()]
    assert [(dep.row_index, dep.column, dep.path) for dep in parsed.path_dependencies] == [
        (0, "blast", str(DANIO_TSV)),
        (0, "comparison_sequence", str(tmp_path / "a.gb")),
        (0, "comparison_sequence", str(tmp_path / "sub" / "b.gb")),
        (1, "blast", str(DANIO_TSV)),
    ]
    table.write_text(f"blast\tcomparison_sequence\n{DANIO_TSV}\ta.gb;;b.gb\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="comparison_sequence"):
        read_conservation_table(str(table))


def _ring_request(sequence_files: object, *, rings: int = 2) -> object:
    from gbdraw.api import CircularDiagramRequest, RecordInput, RenderOutputRequest
    from gbdraw.api.requests import GenBankInputSource

    return CircularDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(path=str(HUMAN))),),
        options=CircularDiagramOptions(
            conservation_blast_files=[str(DANIO_TSV)] * rings,
            conservation_sequence_files=sequence_files,  # type: ignore[arg-type]
            conservation_reference="subject",
        ),
        output=RenderOutputRequest(output_prefix="ring", formats=("svg",)),
    )


def _paths(encoded: object) -> dict[str, Path]:
    return {
        resource.resource_id: resource.source_path
        for resource in encoded.resources  # type: ignore[attr-defined]
        if resource.source_path is not None
    }


def test_the_request_writes_a_list_of_file_references_for_every_ring() -> None:
    encoded = encode_canonical_request(_ring_request([[str(DANIO), str(DROSOPHILA)], None]))
    payload = encoded.payload
    rings = payload["diagramOptions"]["conservationFastaFiles"]
    assert [[ref["representation"] for ref in ring] for ring in rings] == [["file", "file"], []]
    decoded = decode_canonical_request(payload, resource_paths=_paths(encoded), output_directory=".")
    assert decoded.options.conservation_sequence_files == ((str(DANIO), str(DROSOPHILA)), ())


def test_a_released_request_reads_each_ring_file_as_a_one_file_genome() -> None:
    """Schema 9 (v0.14.5) wrote one file reference or null per ring; the Web promotion reads the same vectors."""

    vectors = json.loads((REPO_ROOT / "tests" / "fixtures" / "ring-genome-files-vectors.json").read_text(encoding="utf-8"))
    for case in vectors["cases"]:
        rings = len(case["schema9"])
        encoded = encode_canonical_request(_ring_request([str(DANIO)] * rings, rings=rings))
        payload = json.loads(json.dumps(encoded.payload))
        payload["schema"] = 9
        payload["diagramOptions"]["conservationFastaFiles"] = case["schema9"]
        for record in payload["records"]:
            record["region"] = None
            record.pop("windows")
        paths = {**_paths(encoded), "ring-a": DANIO, "ring-c": DROSOPHILA}
        decoded = decode_canonical_request(payload, resource_paths=paths, output_directory=".")
        expected = tuple(tuple(str(paths[ref["resourceId"]]) for ref in ring) for ring in case["schema10"])
        assert decoded.options.conservation_sequence_files == expected, case["name"]


@pytest.mark.circular
def test_a_saved_session_stores_every_file_of_the_ring(tmp_path: Path) -> None:
    prefix = tmp_path / "saved"
    circular_main([
        "--gbk", str(HUMAN), "--conservation_blast", str(DANIO_TSV), "--conservation_reference", "subject",
        "--conservation_genome", str(DANIO), str(DROSOPHILA), "-o", str(prefix), "-f", "svg", "--save_session",
    ])
    session = session_view(json.loads(prefix.with_suffix(".gbdraw-session.json").read_text(encoding="utf-8")))
    [ring] = session["renderRequest"]["diagramOptions"]["conservationFastaFiles"]
    assert [ref["resourceId"] for ref in ring] == ["conservation-fasta-files-1-1", "conservation-fasta-files-1-2"]
    replay = tmp_path / "replay"
    circular_main(["--session", str(prefix.with_suffix(".gbdraw-session.json")), "-o", str(replay), "-f", "svg"])
    assert _svg(replay) == _svg(prefix)


def test_the_web_helper_reads_every_file_of_a_genome_as_the_cli_ring_does(tmp_path: Path) -> None:
    """The Web reads a ring's files in one helper call: the query FASTA and its raw key are the CLI's."""

    from gbdraw.comparisons.circular_losat import comparison_query_fasta
    from gbdraw.web_support.comparison_sequences import read_comparison_sequence_json

    result = json.loads(read_comparison_sequence_json(str(DANIO), str(DROSOPHILA)))
    genome = comparison_genome((read_comparison_sequence_file(DANIO), read_comparison_sequence_file(DROSOPHILA)))
    assert result["fasta"] == comparison_query_fasta(genome, ordinal=1)
    assert result["recordIds"] == ["NC_002333.2", "NC_024511.2"]
    assert result["recordLabel"] == read_comparison_sequence_file(DANIO).record_label
    joined = _one_file_genome(tmp_path, DANIO, DROSOPHILA)
    assert result["fasta"] == json.loads(read_comparison_sequence_json(str(joined)))["fasta"]
    copy = tmp_path / "copy" / DANIO.name
    copy.parent.mkdir()
    copy.write_bytes(DANIO.read_bytes())
    with pytest.raises(ValidationError) as excinfo:
        read_comparison_sequence_json(str(DANIO), str(copy))
    assert excinfo.value.diagnostic["reason"] == "UNIQUE_IDS"


def test_the_beginner_api_takes_the_files_of_a_ring_genome() -> None:
    """D-66: ``comparison_sequence_source`` takes one path or a sequence of paths, normalized once."""

    from gbdraw.interface import ComparisonRingOptions, ComparisonRingTrackOptions, _ring_losat_values

    rings = ComparisonRingOptions(
        losat="losatn",
        tracks=(
            ComparisonRingTrackOptions(comparison_sequence_source=(DANIO, str(DROSOPHILA))),
            ComparisonRingTrackOptions(comparison_sequence_source=DANIO_FASTA),
        ),
    )
    values = _ring_losat_values(rings)
    options = CircularDiagramOptions(conservation_sequence_files=values["conservation_sequence_files"])
    assert options.conservation_genome_files == ((str(DANIO), str(DROSOPHILA)), (str(DANIO_FASTA),))
    with pytest.raises(ValidationError, match="or a sequence of them"):
        _ring_losat_values(ComparisonRingOptions(
            losat="losatn", tracks=(ComparisonRingTrackOptions(comparison_sequence_source=()),)
        ))


def test_a_restored_cli_table_cell_keeps_every_file_of_its_ring(tmp_path: Path) -> None:
    """A legacy Session's table dependencies of one cell are its ``;`` parts, in order."""

    import base64

    from gbdraw.session_io import _restore_cli_table_paths

    def embedded(name: str, data: bytes) -> dict[str, str]:
        return {"name": name, "encoding": "base64", "data": base64.b64encode(data).decode("ascii")}

    table = b"blast\tcomparison_sequence\nring.tsv\ta.gb;b.gb\n"
    session = {
        "files": {
            "cliInputs": [
                embedded("table.tsv", table), embedded("a.gb", b"A"), embedded("b.gb", b"B"),
            ],
            "cliTables": [{
                "argIndex": 1,
                "slot": "files.cliInputs[0]",
                "dependencies": [
                    {"rowIndex": 0, "column": "comparison_sequence", "slot": "files.cliInputs[1]"},
                    {"rowIndex": 0, "column": "comparison_sequence", "slot": "files.cliInputs[2]"},
                ],
            }],
        }
    }
    run_args = ["--conservation_table", str(tmp_path / "missing.tsv")]
    _restore_cli_table_paths(session, run_args, temp_dir=tmp_path / "restored")
    restored = Path(run_args[1])
    [_header, row] = restored.read_text(encoding="utf-8").splitlines()
    parts = row.split("\t")[1].split(";")
    assert [(restored.parent / part).read_bytes() for part in parts] == [b"A", b"B"]
