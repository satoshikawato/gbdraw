"""``gbdraw locate``: the records and comparisons tables, the new drawing, and the command.

Searches run with a fake LOSAT that answers ``blastn``, ``tblastn`` and
``blastp`` by exact matching (``_FAKE_LOSAT``), on synthetic genomes whose
genes are copies of the reference genes. One case runs a real LOSAT with
tblastn when one resolves (``GBDRAW_TEST_LOSAT_BIN``).
"""

from __future__ import annotations

import json
import re
import os
from dataclasses import replace
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
    LocateWindowChoice,
    choose_locate_windows,
    finish_locate,
    locate_drawing,
    locate_file_counts,
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
from tests.utils.session_drawing_view import drawing_view

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


def _run(argv: list[str], monkeypatch: pytest.MonkeyPatch, command: str = "locate") -> int:
    monkeypatch.setattr(sys, "argv", ["gbdraw", command, *argv])
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


# The reference span as one window, and as two windows (g1..g3 and g4..g6).
ONE_WINDOW = (RecordWindow(1501, 8100),)
TWO_WINDOWS = (RecordWindow(1501, 4500), RecordWindow(4501, 8100, column=2))


# The key `gbdraw linear` and the record planner give the record of a drawing's first input.
SOURCE_KEY = "record-1"


def _plan(
    paths: dict[str, Path], losat: str, targets: tuple[str, ...], windows: tuple[RecordWindow, ...] = ONE_WINDOW
) -> LocatePlan:
    [reference] = _files(paths, "ref")
    search = prepare_locate_search(
        reference, _files(paths, *targets), record=None, record_key=SOURCE_KEY, windows=list(windows), program="tlosatn"
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
    assert drawing.pairs == (((0, 0, 0), (1, 0, 0)), ((1, 0, 0), (2, 0, 0)), ((1, 0, 0), (2, 1, 0)))


def test_the_drawing_pairs_only_the_windows_that_share_a_reference_gene(tmp_path: Path) -> None:
    paths = _genomes(tmp_path)
    drawing = locate_drawing(_plan(paths, _fake_losat(tmp_path), ("b", "a"), TWO_WINDOWS))

    assert [[len(record.windows) for record in row.records] for row in drawing.rows] == [[2], [1, 1], [1]]
    # REF 1501-4500 holds g1..g3 (B1), REF 4501-8100 holds g4..g6 (B2); A1 holds all six.
    assert drawing.pairs == (
        ((0, 0, 0), (1, 0, 0)), ((0, 0, 1), (1, 1, 0)), ((1, 0, 0), (2, 0, 0)), ((1, 1, 0), (2, 0, 0)),
    )
    assert locate_unit_pairs(drawing) == ((0, 2), (1, 3), (2, 4), (3, 4))


@pytest.mark.parametrize(
    ("targets", "windows", "comparisons"),
    [
        (("a", "b"), ONE_WINDOW, [("REF", "A1"), ("A1", "B1"), ("A1", "B2")]),
        (("b", "a"), TWO_WINDOWS, [("REF:1501-4500", "B1"), ("REF:4501-8100", "B2"), ("B1", "A1"), ("B2", "A1")]),
    ],
    ids=["one-window", "two-windows"],
)
def test_the_tables_and_the_new_drawing_give_the_same_records_and_pairs(
    targets: tuple[str, ...], windows: tuple[RecordWindow, ...], comparisons: list[tuple[str, str]], tmp_path: Path
) -> None:
    paths = _genomes(tmp_path)
    plan = _plan(paths, _fake_losat(tmp_path), targets, windows)
    drawing = locate_drawing(plan)
    inputs = {str(path): LocateInputFile(str(path), str(path), None, None) for path in paths.values()}
    records_text, comparisons_text = locate_tables(drawing, inputs, tmp_path)
    (tmp_path / "p.records.tsv").write_text(records_text, encoding="utf-8")
    (tmp_path / "p.comparisons.tsv").write_text(comparisons_text, encoding="utf-8")

    from_table = record_input_manifest_from_table(str(tmp_path / "p.records.tsv")).records
    source = LinearDiagramRequest(records=(RecordInput(GenBankInputSource(str(paths["ref"]))),))
    files = {key: LocateFile(GenBankInputSource(key), Path(key).name) for key in inputs}
    derived = derive_locate_request(source, plan, files).request

    def without_key(record: RecordInput) -> RecordInput:
        return RecordInput(record.source, record.selector, presentation=record.presentation, cardinality=record.cardinality,
                           display=record.display, windows=tuple(record.windows))

    assert [without_key(record) for record in derived.records] == [without_key(record) for record in from_table]
    assert [record.presentation.label for record in derived.records][0] == "ref.gb"
    rows = read_comparisons_table(str(tmp_path / "p.comparisons.tsv")).rows
    assert [(row.query, row.subject) for row in rows] == comparisons

    table_request = LinearDiagramRequest(
        records=from_table,
        options=LinearDiagramOptions(losat_search=LosatSearchOptions("losatp", losatp_mode="pairwise")),
        layout=LinearMultiRecordOptions(),
    )
    units = resolve_request_records(table_request).records
    table_pairs = tuple(
        pair
        for row in rows
        for pair in comparison_unit_pairs(*(
            _resolve_linear_comparison_selector(units, selector, table_path="t", row_number=row.row_number, column="q")
            for selector in (row.query, row.subject)
        ))
    )
    assert derived.options.losat_search is not None
    assert tuple(derived.options.losat_search.pairs or ()) == table_pairs == locate_unit_pairs(drawing)
    assert derived.options.losat_search.losatp_mode == "pairwise"


@pytest.mark.parametrize("whole", [False, True], ids=["part-of-a-record", "whole-record"])
def test_the_new_drawing_writes_a_reversed_window_in_the_form_the_web_cards_write(tmp_path: Path, whole: bool) -> None:
    # A record drawn in one window carries its orientation on the window, a record drawn whole on its presentation
    # (``RecordInput.oriented``), as the Web's File cards write it back; the drawn units are reversed. The
    # part-of-a-record case failed on dev before OV-3365's fix (Python wrote the presentation form). The whole-record
    # case passed on dev already: it pins the form a (-) whole record keeps, which the Web then reads (OV-3376 was
    # the Web reading another drawing's row, tests/web/record-draw-request.test.mjs).
    paths = _genomes(tmp_path)
    plan = _plan(paths, _fake_losat(tmp_path), ("a",))
    [planned] = plan.files[0].windows
    assert planned.window.reverse_complement
    if whole:
        plan = choose_locate_windows(plan, [LocateWindowChoice(
            plan.files[0].file_key, planned.window.column or 0, 1, planned.record_length, True)])
    source = LinearDiagramRequest(records=(RecordInput(GenBankInputSource(str(paths["ref"]))),))
    files = {str(path): LocateFile(GenBankInputSource(str(path)), path.name) for path in paths.values()}
    derived = derive_locate_request(source, plan, files).request
    target = derived.records[1]

    assert [record == record.oriented(record.is_reversed) for record in derived.records] == [True, True]
    assert target.is_reversed
    assert [window.reverse_complement for window in target.windows] == ([] if whole else [True])
    assert target.presentation.reverse_complement is whole
    units = resolve_request_records(derived)
    assert [str(record.id) for record in units.records][1] == "A1"
    assert len(units.records[1]) == (planned.record_length if whole else planned.window.end - planned.window.start + 1)


def test_the_new_drawing_carries_the_edits_of_the_reference_record(tmp_path: Path) -> None:
    from gbdraw.annotations.models import AnnotationOptions, AnnotationSet, CoordinateSpan, RegionAnnotation
    from gbdraw.features.overrides import FeatureOverride

    paths = _genomes(tmp_path)
    plan = _plan(paths, _fake_losat(tmp_path), ("a",))
    plain = LinearDiagramRequest(records=(RecordInput(GenBankInputSource(str(paths["ref"]))),))
    [provenance] = resolve_request_records(plain).provenance
    assert provenance.record_identity_key == SOURCE_KEY
    catalog = {feature.location_parts[0][0]: feature for feature in provenance.source_feature_catalog}
    inside, outside = catalog[3000], catalog[1000]  # g2 lies inside 1501-8100, g0 outside it (0-based starts).
    source = replace(plain, options=LinearDiagramOptions(
        feature_overrides=(
            FeatureOverride(SOURCE_KEY, inside.biological_feature_id, label_text="kept"),
            FeatureOverride(SOURCE_KEY, outside.biological_feature_id, label_text="dropped"),
        ),
        annotations=AnnotationOptions(sets=(AnnotationSet(id="marks", annotations=(
            RegionAnnotation("in", CoordinateSpan(None, 3001, 3300)),
            RegionAnnotation("out", CoordinateSpan(None, 9001, 9300)),
        )),)),
    ))
    files = {str(path): LocateFile(GenBankInputSource(str(path)), path.name) for path in paths.values()}
    derived = derive_locate_request(source, plan, files)

    assert [(row.record_key, row.label_text) for row in derived.request.options.feature_overrides] == [(SOURCE_KEY, "kept")]
    assert derived.request.records[0].record_key == SOURCE_KEY
    assert derived.request.options.annotations is not None
    assert [item.id for item in derived.request.options.annotations.sets[0].annotations] == ["in"]
    assert [(item.kind, item.count) for item in derived.not_carried] == [("ANNOTATIONS_OUTSIDE", 1)]


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

# A LOSAT with tblastn whose every search fails.
_FAILS = """#!{python}
import sys
if sys.argv[1:] == ["--version"]:
    print("losat 0.2.0")
    sys.exit(0)
if sys.argv[2:] == ["--help"]:
    sys.exit(0)
print("out of memory", file=sys.stderr)
sys.exit(3)
"""


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda args, paths: args + ["--reference_record", "NOPE"], "has no record 'NOPE'"),
        (lambda args, paths: args + ["--reference_region", "7000-9000"], "Reference windows overlap."),
        (lambda args, paths: args + ["--reference_region", "11001-13000"], "runs past the record end (12000)"),
        (lambda args, paths: [*args[:3], "9500-9900", *args[4:]], "holds no whole CDS"),
        (lambda args, paths: args + ["--targets", str(paths["ref"])], "is the reference file"),
        (lambda args, paths: args + ["--ncbi_blast_bin", "/opt/tblastn"], "not both"),
        (lambda args, paths: [*args[:-2], "-o", args[-1], "--losat_bin", "NO_TBLASTN"], "(0.1.0) has no tblastn"),
        (lambda args, paths: [*args[:5], str(paths["c"]), *args[6:]], "No target file has a window"),
        (lambda args, paths: [*args[:-2], "-o", args[-1], "--losat_bin", "FAILS"], "Every search failed: "),
    ],
    ids=[
        "unknown-record", "overlap", "outside", "no-cds", "reference-target", "two-runtimes", "no-tblastn", "no-window",
        "searches-failed",
    ],
)
def test_the_command_stops_with_a_message(
    change, message: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    paths = _genomes(tmp_path)
    scripts = {}
    for name, text in (("NO_TBLASTN", _NO_TBLASTN), ("FAILS", _FAILS)):
        scripts[name] = tmp_path / name.lower()
        scripts[name].write_text(text.format(python=sys.executable), encoding="utf-8")
        scripts[name].chmod(0o755)
    args = [str(scripts[item]) if item in scripts else item
            for item in change(_locate_args(paths, _fake_losat(tmp_path), tmp_path / "p", "a"), paths)]

    assert _run(args, monkeypatch) == 1
    error = capsys.readouterr().err
    assert message in error
    if "tblastn" in message:
        assert "--losat (losatn)" in error and "--ncbi_blast_bin" in error
    if message.startswith("Every search failed"):
        assert "the search failed" in error and "Lower the thresholds" not in error
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


def _losatp_runtime(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """``gbdraw linear`` resolves the stand-in LOSAT (or ``GBDRAW_TEST_LOSAT_BIN``) as if it were the installed
    one, so the receipt's command runs as written, without ``--losat_bin``."""

    from gbdraw.comparisons import losat_runtime

    losat = os.environ.get("GBDRAW_TEST_LOSAT_BIN") or _fake_losat(tmp_path)
    monkeypatch.setattr(
        losat_runtime, "resolve_losat_runtime", lambda *_args, **_kwargs: losat_runtime.LosatRuntime("losat", losat, "explicit")
    )


def _locate_two_windows(paths: dict[str, Path], tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    prefix = tmp_path / "P"
    args = _locate_args(paths, _fake_losat(tmp_path), prefix, "b", "a")
    args[args.index("1501-8100")] = "1501-4500"
    assert _run([*args, "--reference_region", "4501-8100"], monkeypatch) == 0
    return prefix


def test_the_receipts_command_draws_the_plan_with_only_the_listed_window_pairs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    paths = _genomes(tmp_path)
    _locate_two_windows(paths, tmp_path, monkeypatch)
    _losatp_runtime(monkeypatch, tmp_path)
    monkeypatch.chdir(tmp_path)
    caplog.set_level("INFO")

    # PD-OI-123: `gbdraw linear --records_table PREFIX.records.tsv --comparisons_table PREFIX.comparisons.tsv
    # --losat losatp -o FIG` draws the plan.
    command = "--records_table P.records.tsv --comparisons_table P.comparisons.tsv --losat losatp -o FIG".split()
    assert _run(command, monkeypatch, "linear") == 0
    assert (tmp_path / "FIG.svg").is_file()
    assert not [record.getMessage() for record in caplog.records if "no protein match" in record.getMessage()]

    # The Session stores the tables of the four window pairs, and replays them.
    assert _run([*command[:-1], "FIG2", "--save_session"], monkeypatch, "linear") == 0
    session = drawing_view(json.loads((tmp_path / "FIG2.gbdraw-session.json").read_text(encoding="utf-8")))
    assert _stored_pairs(session) == [(0, 2), (1, 3), (2, 4), (3, 4)]
    assert _run(["--session", "FIG2.gbdraw-session.json", "-o", "FIG3"], monkeypatch, "linear") == 0
    assert (tmp_path / "FIG3.svg").read_bytes() == (tmp_path / "FIG2.svg").read_bytes()


def test_a_comparisons_table_with_another_losatp_mode_names_the_pairwise_mode(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    paths = _genomes(tmp_path)
    prefix = _locate_two_windows(paths, tmp_path, monkeypatch)
    _losatp_runtime(monkeypatch, tmp_path)
    capsys.readouterr()

    assert _run([
        "--records_table", f"{prefix}.records.tsv", "--comparisons_table", f"{prefix}.comparisons.tsv",
        "--losat", "losatp", "--losatp_mode", "similarity_groups", "-o", str(tmp_path / "FIG"),
    ], monkeypatch, "linear") == 1
    error = capsys.readouterr().err
    assert "source=losat rows list the pairs to compare" in error and "use --losatp_mode pairwise" in error


def test_a_comparisons_table_names_a_window_by_its_record_and_span(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from gbdraw.api.record_planning import resolve_linear_options

    paths = _genomes(tmp_path)
    prefix = _locate_two_windows(paths, tmp_path, monkeypatch)
    units = resolve_request_records(LinearDiagramRequest(
        records=record_input_manifest_from_table(f"{prefix}.records.tsv").records, layout=LinearMultiRecordOptions(),
    )).records

    def resolve(selector: str) -> tuple[int, ...]:
        return _resolve_linear_comparison_selector(units, selector, table_path="t", row_number=1, column="query")

    assert [resolve(selector) for selector in ("REF", "#1", "REF:1501-4500", "#1:4501-8100", "B2", "#3:2101-6403")] == [
        (0, 1), (0, 1), (0,), (1,), (3,), (3,),
    ]
    with pytest.raises(ValidationError, match=r"record REF has no window 1-100; its windows are 1501-4500, 4501-8100\."):
        resolve("REF:1-100")
    with pytest.raises(ValidationError, match="unresolved record selector 'REF:4501-8100:rc'"):
        resolve("REF:4501-8100:rc")
    # A record ID that reads as a region still names its record.
    named = [_record("chr1:1-100", "ACGT"), _record("chr2", "ACGT")]
    assert _resolve_linear_comparison_selector(named, "chr1:1-100", table_path="t", row_number=1, column="query") == (0,)
    # A comparison file compares whole drawn records, also when a selector names one window.
    blast = tmp_path / "pair.tsv"
    blast.write_text("q\ts\t90\t100\t0\t0\t1\t100\t1\t100\t1e-30\t200\n", encoding="utf-8")
    table = tmp_path / "file.tsv"
    table.write_text(f"blast\tquery\tsubject\n{blast}\tREF:1501-4500\tB1\n", encoding="utf-8")
    layout = LinearMultiRecordOptions(multi_record_positions=("#1@1", "#2@1", "#3@2", "#4@2", "#5@3"))
    with pytest.raises(ValidationError, match="drawn in several windows"):
        resolve_linear_options(LinearDiagramOptions(comparison_table_file=str(table)), records=units, layout=layout)


def _write_genbank(path: Path, record: SeqRecord, record_id: str) -> Path:
    record.id = record.name = record_id
    record.annotations = {"molecule_type": "DNA", "topology": "linear"}
    SeqIO.write([record], path, "genbank")
    return path


def _stored_pairs(session: dict) -> list[tuple[int, int]]:
    """The unit pairs whose tables a Pairwise Session stores; its pipeline marker searches nothing (OV-3037)."""

    comparisons = session["renderRequest"]["comparisons"]
    [protein] = [item for item in comparisons if item["kind"] == "generatedProteinComparison"]
    assert (protein["mode"], protein["pairs"]) == ("none", [])
    return [
        (item["queryRecordIndex"], item["subjectRecordIndex"])
        for item in comparisons if item["kind"] == "precomputedProteinComparison"
    ]


@pytest.mark.parametrize("case", ["id-with-colon", "id-of-another-window"])
def test_a_window_whose_record_id_would_not_read_back_is_named_by_its_number(
    case: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # A record drawn in two windows (g1..g3 and g4..g6, 20 kb apart). "ctg:1:START-END" would read as file "ctg",
    # record "1"; "chr1:START-END" would name the other drawn record whose ID is literally that text.
    paths = _genomes(tmp_path)
    reference = SeqIO.read(paths["ref"], "genbank")
    rng = random.Random(11)
    split = reference[1900:4400] + _record("x", _filler(rng, 20_000)) + reference[4900:7400] + _record("x", _filler(rng, 2000))
    split_id = "ctg:1" if case == "id-with-colon" else "chr1"
    paths["x"] = _write_genbank(tmp_path / "x.gb", split, split_id)
    targets = ["x"]
    if case == "id-of-another-window":
        first = next(w for w in _plan(paths, _fake_losat(tmp_path), ("x",)).files[0].windows if w.window.column == 1)
        paths["y"] = _write_genbank(tmp_path / "y.gb", reference[1500:8000], f"chr1:{first.window.start}-{first.window.end}")
        targets.append("y")
    prefix = tmp_path / "P"
    assert _run(_locate_args(paths, _fake_losat(tmp_path), prefix, *targets), monkeypatch) == 0

    rows = [(row.query, row.subject) for row in read_comparisons_table(f"{prefix}.comparisons.tsv").rows]
    assert [query for query, _subject in rows[:2]] == ["REF", "REF"]
    if case == "id-with-colon":
        assert all(subject.startswith("#2:") for _query, subject in rows)
    else:
        # The first window takes its number; the second keeps its ID; the other record keeps its literal ID.
        assert rows[0][1] == f"#2:{first.window.start}-{first.window.end}"
        assert rows[1][1].startswith("chr1:") and rows[2] == (rows[0][1], f"chr1:{first.window.start}-{first.window.end}")
    _losatp_runtime(monkeypatch, tmp_path)
    monkeypatch.chdir(tmp_path)
    command = "--records_table P.records.tsv --comparisons_table P.comparisons.tsv --losat losatp -o FIG --save_session"
    assert _run(command.split(), monkeypatch, "linear") == 0
    session = drawing_view(json.loads((tmp_path / "FIG.gbdraw-session.json").read_text(encoding="utf-8")))
    expected = [(0, 1), (0, 2)] + ([(1, 3), (2, 3)] if case == "id-of-another-window" else [])
    assert _stored_pairs(session) == expected


def _interleaved(tmp_path: Path) -> dict[str, Path]:
    """D holds g1-g2 at the start of contig D1, g3-g4 on contig D2 and g5-g6 far down D1: a File whose windows run
    D1, D2, D1 in the order of the reference, two of one record around one of another."""

    paths = _genomes(tmp_path)
    reference = SeqIO.read(paths["ref"], "genbank")
    rng = random.Random(11)
    d1 = reference[1900:3400] + _record("D1", _filler(rng, 14_000)) + reference[5900:7400] + _record("D1", _filler(rng, 2000))
    d2 = _record("D2", _filler(rng, 2000)) + reference[3900:5400] + _record("D2", _filler(rng, 2000))
    for record, record_id in ((d1, "D1"), (d2, "D2")):
        record.id = record.name = record_id
        record.annotations = {"molecule_type": "DNA", "topology": "linear"}
    paths["d"] = tmp_path / "d.gb"
    SeqIO.write([d1, d2], paths["d"], "genbank")
    return paths


def test_interleaved_windows_of_a_file_are_numbered_and_drawn_in_the_order_of_the_reference(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    # PD-OI-122: a row orders its windows by column, ties by input order, and numbers them in that order; PD-OI-123:
    # the report and the proposal give each window the number it has in the drawing. A record whose windows stand
    # around another record's keeps that order through the records table and `gbdraw linear`.
    paths = _interleaved(tmp_path)
    caplog.set_level("INFO")
    assert _run(_locate_args(paths, _fake_losat(tmp_path), tmp_path / "P", "d"), monkeypatch) == 0
    plan = LocatePlan.from_json(json.loads((tmp_path / "P.plan.json").read_text(encoding="utf-8")))
    windows = plan.files[0].windows
    assert [(item.record_id, item.window.column) for item in windows] == [("D1", 1), ("D2", 2), ("D1", 3)]
    spans = [f"{item.window.start:,}-{item.window.end:,}" for item in windows]

    # The report: each window under the number the drawing gives it.
    listed = [message for message in (record.getMessage() for record in caplog.records) if message.startswith("  ")]
    assert [line.split(":")[0] for line in listed[:3]] == [f"  1 D1 {spans[0]}", f"  2 D2 {spans[1]}", f"  3 D1 {spans[2]}"]

    # The records table: the D1 rows share a record_group, and each row keeps its window's column.
    table = read_records_table(str(tmp_path / "P.records.tsv"))
    rows = sorted(table.rows[1:], key=lambda row: row.column or 0)
    assert [(row.column, row.record_id, row.row) for row in rows] == [(1, "D1", 2), (2, "D2", 2), (3, "D1", 2)]
    assert rows[0].record_group == rows[2].record_group != rows[1].record_group

    # `gbdraw linear` draws the row's windows in that order, numbered 1 to 3 in its window list.
    monkeypatch.chdir(tmp_path)
    assert _run("--records_table P.records.tsv -o FIG".split(), monkeypatch, "linear") == 0
    drawn = re.findall(r">(D[12]:[^<]*)<", (tmp_path / "FIG.svg").read_text(encoding="utf-8"))
    assert drawn == [f"D1:{spans[0]}", f"D2:{spans[1]}", f"D1:{spans[2]}"]


def _upstream_copy(tmp_path: Path) -> dict[str, Path]:
    """E holds a copy of g1-g2 more than 10 kb before a whole copy of g1..g6 on one contig: the first is a window
    the plan does not select (column 1), the second the window it keeps (column 2)."""

    paths = _genomes(tmp_path)
    reference = SeqIO.read(paths["ref"], "genbank")
    rng = random.Random(3)
    e = (_record("E1", _filler(rng, 1000)) + reference[1900:3400] + _record("E1", _filler(rng, 15_000))
         + reference[1500:8000] + _record("E1", _filler(rng, 2000)))
    paths["e"] = _write_genbank(tmp_path / "e.gb", e, "E1")
    return paths


def test_a_kept_window_after_one_not_selected_is_numbered_as_the_drawing_numbers_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    # A window's number is its place in the drawn row, found by its column in the plan (its identity there, windows
    # not selected counted): the report and the drawing give the kept window the number 1.
    paths = _upstream_copy(tmp_path)
    caplog.set_level("INFO")
    assert _run(_locate_args(paths, _fake_losat(tmp_path), tmp_path / "P", "e"), monkeypatch) == 0
    plan = LocatePlan.from_json(json.loads((tmp_path / "P.plan.json").read_text(encoding="utf-8")))
    [kept] = plan.files[0].windows
    assert kept.window.column == 2 and [item.window.column for item in plan.files[0].candidates] == [1]
    [numbers] = [counts.window_numbers for counts in locate_file_counts(plan)]
    assert numbers == ((2, 1),)
    span = f"{kept.window.start:,}-{kept.window.end:,}"
    listed = [message for message in (record.getMessage() for record in caplog.records) if message.startswith("  ")]
    assert listed[0].split(":")[0] == f"  1 E1 {span}"
