"""``gbdraw locate``: the records and comparisons tables, the new drawing, and the command.

Searches run with a fake LOSAT that answers ``blastn``, ``tblastn`` and
``blastp`` by exact matching (``_FAKE_LOSAT``), on synthetic genomes whose
genes are copies of the reference genes. One case runs a real LOSAT with
tblastn when one resolves (``GBDRAW_TEST_LOSAT_BIN``).
"""

from __future__ import annotations

import json
import os
import random
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw import cli
from gbdraw.api.locate_drawing import LocateFile, derive_locate_request, locate_unit_pairs
from gbdraw.api.options import LinearDiagramOptions, LinearMultiRecordOptions, LosatSearchOptions
from gbdraw.api.record_planning import (
    _resolve_linear_comparison_selector,
    comparison_unit_pairs,
    record_input_manifest_from_table,
)
from gbdraw.api.request_render import resolve_request_records
from gbdraw.api.requests import GenBankInputSource, LinearDiagramRequest, RecordInput
from gbdraw.cli_utils.locate import locate_tables
from gbdraw.comparisons.locate import (
    LocateInputFile,
    LocateOptions,
    LocatePlan,
    LocateSearchResult,
    finish_locate,
    locate_drawing,
    locate_results,
    prepare_locate_search,
)
from gbdraw.comparisons.losat_runtime import require_losat_program, run_losat_search
from gbdraw.exceptions import ValidationError
from gbdraw.io.cli_tables import (
    ComparisonsTableEntry,
    RecordsTableEntry,
    comparisons_table_text,
    read_comparisons_table,
    read_records_table,
    records_table_text,
)
from gbdraw.io.regions import RecordWindow
from gbdraw.web_support.locate import LocateFileInput, LocateSearchInput, locate_search_plan_json, locate_windows_json

_FAKE_LOSAT = """#!{python}
import sys
from Bio import SeqIO
from Bio.Seq import Seq

args = sys.argv[1:]
if args == ["--version"]:
    print("losat 0.2.0")
    sys.exit(0)
if len(args) == 2 and args[1] == "--help":
    sys.exit(0 if args[0] in ("blastn", "tblastn", "blastp") else 2)
program = args[0]
query = list(SeqIO.parse(args[args.index("-query") + 1], "fasta"))
subject = list(SeqIO.parse(args[args.index("-subject") + 1], "fasta"))


def row(q, s, length, qend, sstart, send):
    print(f"{{q.id}}\\t{{s.id}}\\t100.00\\t{{length}}\\t0\\t0\\t1\\t{{qend}}\\t{{sstart}}\\t{{send}}\\t1e-50\\t{{2.0 * length:.1f}}")


def hits(text, word):
    start = text.find(word)
    while start >= 0:
        yield start
        start = text.find(word, start + 1)


for q in query:
    word = str(q.seq)
    for s in subject:
        sequence = str(s.seq)
        size = len(sequence)
        if program == "blastn":
            for start in hits(sequence, word):
                row(q, s, len(word), len(word), start + 1, start + len(word))
            for start in hits(sequence, str(Seq(word).reverse_complement())):
                row(q, s, len(word), len(word), start + len(word), start + 1)
        elif program == "tblastn":
            reverse = str(Seq(sequence).reverse_complement())
            for frame in range(3):
                for strand, text in (("+", sequence), ("-", reverse)):
                    usable = text[frame:]
                    protein = str(Seq(usable[: len(usable) // 3 * 3]).translate())
                    for start in hits(protein, word):
                        low = frame + 3 * start + 1
                        high = low + 3 * len(word) - 1
                        if strand == "+":
                            row(q, s, len(word), len(word), low, high)
                        else:
                            row(q, s, len(word), len(word), size - low + 1, size - high + 1)
        elif word == sequence:
            row(q, s, len(word), len(word), 1, len(word))
"""

_CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG", "TGA")]
GENE_COUNT = 6


def _fake_losat(tmp_path: Path) -> str:
    path = tmp_path / "losat"
    path.write_text(_FAKE_LOSAT.format(python=sys.executable), encoding="utf-8")
    path.chmod(0o755)
    return str(path)


def _orf(rng: random.Random, codons: int) -> str:
    return "ATG" + "".join(rng.choice(_CODONS) for _ in range(codons)) + "TAA"


def _filler(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _record(record_id: str, sequence: str) -> SeqRecord:
    return SeqRecord(Seq(sequence), id=record_id, name=record_id, description=record_id,
                     annotations={"molecule_type": "DNA", "topology": "linear"})


def _genomes(tmp_path: Path) -> dict[str, Path]:
    """REF holds genes g0..g7 at 1000 bp steps from 1,001; the window 1,501-8,000 holds g1..g6 whole and cuts g7.
    A holds the window reversed in one record; B holds g1..g3 and g4..g6 on two contigs; C holds none."""

    rng = random.Random(7)
    sequence = list(_filler(rng, 12_000))
    reference = _record("REF", "")
    features = []
    for index in range(8):
        start = 1000 * (index + 1)
        gene = _orf(rng, 99)
        sequence[start:start + len(gene)] = gene
        features.append(SeqFeature(
            FeatureLocation(start, start + len(gene), strand=1), type="CDS",
            qualifiers={"locus_tag": [f"g{index}"], "translation": [str(Seq(gene).translate(to_stop=True))]},
        ))
    reference = _record("REF", "".join(sequence))
    reference.features = features
    a = _record("A1", _filler(rng, 4000)) + reference[1500:8000].reverse_complement(
        id="A1", name="A1", description="A1", features=True, annotations=True
    ) + _record("A1", _filler(rng, 4000))
    b1 = reference[1900:4400] + _record("B1", _filler(rng, 3000))
    b2 = _record("B2", _filler(rng, 3000)) + reference[4900:7400] + _record("B2", _filler(rng, 3000))
    c1 = _record("C1", _filler(rng, 9000))
    paths = {}
    for name, records in (("ref", [reference]), ("a", [a]), ("b", [b1, b2]), ("c", [c1])):
        for record, record_id in zip(records, {"ref": ["REF"], "a": ["A1"], "b": ["B1", "B2"], "c": ["C1"]}[name]):
            record.id = record.name = record_id
            record.annotations = {"molecule_type": "DNA", "topology": "linear"}
        paths[name] = tmp_path / f"{name}.gb"
        SeqIO.write(records, paths[name], "genbank")
    return paths


def _run(argv: list[str], monkeypatch: pytest.MonkeyPatch) -> int:
    monkeypatch.setattr(sys, "argv", ["gbdraw", "locate", *argv])
    try:
        cli.main()
    except SystemExit as exc:
        return int(exc.code or 0)
    return 0


def _locate_args(paths: dict[str, Path], losat: str, prefix: Path, *targets: str) -> list[str]:
    return [
        "--reference_gbk", str(paths["ref"]), "--reference_region", "1501-8100",
        "--targets", *(str(paths[name]) for name in targets), "--losat_bin", losat, "-o", str(prefix),
    ]


def test_tables_written_by_their_owner_read_back_unchanged(tmp_path: Path) -> None:
    (tmp_path / "g.gb").write_text("", encoding="utf-8")
    records = tmp_path / "plan.records.tsv"
    records.write_text(records_table_text([
        RecordsTableEntry(gbk="g.gb", record_label="g.gb", record_id="R1", region="30-10:origin",
                          reverse_complement=True, row=1, column=1, record_group="1:R1"),
        RecordsTableEntry(gbk="g.gb", record_id="R1", region="40-60", row=1, column=3, record_group="1:R1"),
        RecordsTableEntry(gbk="g.gb", record_id="R2", row=1, column=2),
    ]), encoding="utf-8")
    comparisons = tmp_path / "plan.comparisons.tsv"
    comparisons.write_text(comparisons_table_text([ComparisonsTableEntry("R1", "#3")]), encoding="utf-8")

    table = read_records_table(str(records))
    assert [
        (row.gbk, row.record_label, row.record_id, row.region, row.reverse_complement, row.row, row.column, row.record_group)
        for row in table.rows
    ] == [
        (str(tmp_path / "g.gb"), "g.gb", "R1", "30-10:origin", True, 1, 1, "1:R1"),
        (str(tmp_path / "g.gb"), "", "R1", "40-60", False, 1, 3, "1:R1"),
        (str(tmp_path / "g.gb"), "", "R2", "", False, 1, 2, ""),
    ]
    assert [len(record.rows) for record in table.records] == [2, 1]
    assert [(row.blast, row.query, row.subject, row.source) for row in read_comparisons_table(str(comparisons)).rows] == [
        ("", "R1", "#3", "losat")
    ]


def _files(paths: dict[str, Path], *names: str) -> list[LocateInputFile]:
    return [LocateInputFile(str(paths[name]), str(paths[name]), None, None) for name in names]


def _plan(paths: dict[str, Path], losat: str, targets: tuple[str, ...]) -> LocatePlan:
    [reference] = _files(paths, "ref")
    search = prepare_locate_search(
        reference, _files(paths, *targets), record=None, windows=[RecordWindow(1501, 8100)], program="tlosatn"
    )
    outcomes = [
        LocateSearchResult(job.file_key, run_losat_search("tlosatn", job.query_fasta, job.subject_fasta, losat_bin=losat))
        for job in search.plan.jobs
    ]
    return finish_locate(search, locate_results(search, outcomes), LocateOptions.for_program())


def test_the_plan_finds_the_genes_reversed_and_split_over_contigs(tmp_path: Path) -> None:
    paths = _genomes(tmp_path)
    plan = _plan(paths, _fake_losat(tmp_path), ("a", "b", "c"))

    a, b, c = plan.files
    assert [(w.record_id, w.window.reverse_complement) for w in a.windows] == [("A1", True)]
    assert [(w.record_id, len(w.genes)) for w in b.windows] == [("B1", 3), ("B2", 3)]
    assert b.windows[0].touches_start
    assert [w.code for w in c.warnings] == ["TARGET_NOT_ANNOTATED", "FILE_NO_WINDOWS"]
    assert [cds.locus_tag for cds in plan.reference.cut_genes] == ["g7"]
    drawing = locate_drawing(plan)
    assert [[record.record_id for record in row.records] for row in drawing.rows] == [["REF"], ["A1"], ["B1", "B2"]]
    assert drawing.pairs == (((0, 0), (1, 0)), ((1, 0), (2, 0)), ((1, 0), (2, 1)))


def test_the_tables_and_the_new_drawing_give_the_same_records_and_pairs(tmp_path: Path) -> None:
    paths = _genomes(tmp_path)
    plan = _plan(paths, _fake_losat(tmp_path), ("a", "b"))
    drawing = locate_drawing(plan)
    inputs = {str(path): LocateInputFile(str(path), str(path), None, None) for path in paths.values()}
    records_text, comparisons_text = locate_tables(drawing, inputs, tmp_path)
    (tmp_path / "p.records.tsv").write_text(records_text, encoding="utf-8")
    (tmp_path / "p.comparisons.tsv").write_text(comparisons_text, encoding="utf-8")

    from_table = record_input_manifest_from_table(str(tmp_path / "p.records.tsv")).records
    source = LinearDiagramRequest(records=(RecordInput(GenBankInputSource(str(paths["ref"])), record_key="reference"),))
    files = {key: LocateFile(GenBankInputSource(key), Path(key).name) for key in inputs}
    derived = derive_locate_request(source, plan, files).request

    def without_key(record: RecordInput) -> RecordInput:
        return RecordInput(record.source, record.selector, presentation=record.presentation, cardinality=record.cardinality,
                           display=record.display, windows=tuple(record.windows))

    assert [without_key(record) for record in derived.records] == [without_key(record) for record in from_table]
    assert [record.presentation.label for record in derived.records] == ["ref.gb", "a.gb", "b.gb", None]

    table_request = LinearDiagramRequest(
        records=from_table,
        options=LinearDiagramOptions(losat_search=LosatSearchOptions("losatp", losatp_mode="pairwise")),
        layout=LinearMultiRecordOptions(),
    )
    units = resolve_request_records(table_request).records
    table_pairs = tuple(
        pair
        for row in read_comparisons_table(str(tmp_path / "p.comparisons.tsv")).rows
        for pair in comparison_unit_pairs(*(
            _resolve_linear_comparison_selector(units, selector, table_path="t", row_number=row.row_number, column="q")
            for selector in (row.query, row.subject)
        ))
    )
    assert derived.options.losat_search is not None
    assert tuple(derived.options.losat_search.pairs or ()) == table_pairs == locate_unit_pairs(drawing)
    assert derived.options.losat_search.losatp_mode == "pairwise"


def test_the_command_writes_the_tables_the_plan_and_a_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    paths = _genomes(tmp_path)
    prefix = tmp_path / "out" / "fen"
    prefix.parent.mkdir()
    caplog.set_level("INFO")
    code = _run([*_locate_args(paths, _fake_losat(tmp_path), prefix, "a", "b", "c"), "--losat_output_dir", str(tmp_path / "raw")], monkeypatch)

    assert code == 0
    table = read_records_table(str(prefix) + ".records.tsv")
    assert [(row.record_id, row.row, row.column, row.reverse_complement, bool(row.record_label)) for row in table.rows] == [
        ("REF", 1, 1, False, True), ("A1", 2, 1, True, True), ("B1", 3, 1, False, True), ("B2", 3, 2, False, False),
    ]
    assert table.rows[0].gbk == str(paths["ref"])
    assert [(row.query, row.subject) for row in read_comparisons_table(str(prefix) + ".comparisons.tsv").rows] == [
        ("REF", "A1"), ("A1", "B1"), ("A1", "B2"),
    ]
    plan = json.loads(Path(str(prefix) + ".plan.json").read_text(encoding="utf-8"))
    assert plan["options"]["program"] == "tlosatn" and plan["runtime"]["version"] == "0.2.0"
    assert [len(file["windows"]) for file in plan["files"]] == [1, 2, 0]
    assert sorted(path.name for path in (tmp_path / "raw").iterdir()) == ["1.a.gb.tlosatn.tsv", "2.b.gb.tlosatn.tsv", "3.c.gb.tlosatn.tsv"]
    messages = [record.getMessage() for record in caplog.records]
    assert "1 reference CDS cut by a window edge are not searched: g7." in messages
    assert f"  1 A1 {plan['files'][0]['windows'][0]['window']['start']:,}-{plan['files'][0]['windows'][0]['window']['end']:,} (−): 6 gene(s), g1 to g6." in messages
    assert any(message.startswith(f"{paths['c']}: no window; the thresholds kept 0 of 0 rows") for message in messages)

    # An existing output stops the command unless --overwrite is given.
    assert _run(_locate_args(paths, _fake_losat(tmp_path), prefix, "a"), monkeypatch) == 1
    assert _run([*_locate_args(paths, _fake_losat(tmp_path), prefix, "a"), "--overwrite"], monkeypatch) == 0


_NO_TBLASTN = """#!{python}
import sys
if sys.argv[1:] == ["--version"]:
    print("losat 0.1.0")
    sys.exit(0)
sys.exit(2 if sys.argv[1:2] == ["tblastn"] else 0)
"""


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda args, paths: args + ["--reference_record", "NOPE"], "has no record 'NOPE'"),
        (lambda args, paths: args + ["--reference_region", "7000-9000"], "Reference windows overlap."),
        (lambda args, paths: [*args[:3], "9500-9900", *args[4:]], "holds no whole CDS"),
        (lambda args, paths: args + ["--targets", str(paths["ref"])], "is the reference file"),
        (lambda args, paths: args + ["--ncbi_blast_bin", "/opt/tblastn"], "not both"),
        (lambda args, paths: [*args[:-2], "-o", args[-1], "--losat_bin", "NO_TBLASTN"], "(0.1.0) has no tblastn"),
        (lambda args, paths: [*args[:5], str(paths["c"]), *args[6:]], "No target file has a window"),
    ],
    ids=["unknown-record", "overlap", "no-cds", "reference-target", "two-runtimes", "no-tblastn", "no-window"],
)
def test_the_command_stops_with_a_message(
    change, message: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    paths = _genomes(tmp_path)
    no_tblastn = tmp_path / "losat-0.1.0"
    no_tblastn.write_text(_NO_TBLASTN.format(python=sys.executable), encoding="utf-8")
    no_tblastn.chmod(0o755)
    args = [str(no_tblastn) if item == "NO_TBLASTN" else item
            for item in change(_locate_args(paths, _fake_losat(tmp_path), tmp_path / "p", "a"), paths)]

    assert _run(args, monkeypatch) == 1
    error = capsys.readouterr().err
    assert message in error
    if "tblastn" in message:
        assert "--losat (losatn)" in error and "--ncbi_blast_bin" in error
    assert not (tmp_path / "p.records.tsv").exists()


def _tblastn_runtime() -> str:
    losat_bin = os.environ.get("GBDRAW_TEST_LOSAT_BIN") or None
    try:
        require_losat_program("tlosatn", losat_bin=losat_bin)
    except ValidationError as error:
        pytest.skip(
            "needs a LOSAT runtime with tblastn (the pinned LOSAT 0.1.0 has none and LOSAT v0.2.0 is unreleased); "
            f"set GBDRAW_TEST_LOSAT_BIN: {error}"
        )
    return losat_bin or "losat"


def test_the_command_finds_the_windows_with_a_real_tblastn(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    losat = _tblastn_runtime()
    paths = _genomes(tmp_path)
    prefix = tmp_path / "real"
    args = _locate_args(paths, losat, prefix, "a", "b")
    if losat == "losat":
        args = args[: args.index("--losat_bin")] + args[args.index("--losat_bin") + 2:]

    assert _run(args, monkeypatch) == 0
    fake = _plan(paths, _fake_losat(tmp_path), ("a", "b"))
    plan = json.loads(Path(str(prefix) + ".plan.json").read_text(encoding="utf-8"))
    assert [[(w["record_id"], w["window"]) for w in file["windows"]] for file in plan["files"]] == [
        [(w.record_id, json.loads(json.dumps(fake.to_json()["files"][index]["windows"][n]["window"]))) for n, w in enumerate(file_plan.windows)]
        for index, file_plan in enumerate(fake.files)
    ]


def test_the_worker_payload_gives_the_plan_the_command_writes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    paths = _genomes(tmp_path)
    losat = _fake_losat(tmp_path)
    prefix = tmp_path / "cli"
    assert _run([*_locate_args(paths, losat, prefix, "a", "b", "c"), "--flank", "500", "--max_windows", "4"], monkeypatch) == 0

    def file(name: str) -> LocateFileInput:
        return {"key": str(paths[name]), "gbk": str(paths[name]), "gff": None, "fasta": None}

    search: LocateSearchInput = {
        "reference": file("ref"), "recordId": None, "windows": [{"start": 1501, "end": 8100, "reverseComplement": False}],
        "targets": [file("a"), file("b"), file("c")], "program": "tlosatn",
    }
    jobs = locate_search_plan_json(search)
    assert [job["fileKey"] for job in jobs["jobs"]] == [str(paths[name]) for name in ("a", "b", "c")]
    reply = locate_windows_json({
        **search,
        "outcomes": [
            {"fileKey": job["fileKey"], "text": run_losat_search("tlosatn", job["queryFasta"], job["subjectFasta"], losat_bin=losat), "error": None}
            for job in jobs["jobs"]
        ],
        "options": {"evalue": 1e-2, "bitscore": 50.0, "identity": 0.0, "alignmentLength": 0, "flank": 500, "mergeDistance": 10_000, "maxWindows": 4},
    })

    written = json.loads(Path(str(prefix) + ".plan.json").read_text(encoding="utf-8"))
    assert written.pop("runtime")["program"] == "tblastn"
    assert reply["plan"].pop("runtime") is None
    assert json.loads(json.dumps(reply["plan"])) == written
