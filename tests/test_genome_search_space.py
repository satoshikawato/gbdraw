"""Comparison E-values are computed per genome (OV-3143).

A genome is the records of the files of one Linear row, every record of those files counted, drawn or not. Every
search between two rows passes LOSAT one search space, the query genome's length times the subject genome's, so
a hit has one E-value whether its records are drawn whole, cropped, in windows, or from a genome split over files.
A search between two records that are each the whole of a one-record file passes none: LOSAT's own space is then
the genome pair's, and those figures stay as they were.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import pytest
from Bio import SeqIO

import gbdraw.comparisons.losat_runtime as runtime_module
from gbdraw.comparisons.losat_jobs import genome_database_sizes, genome_search_spaces
from gbdraw.linear import linear_main
from tests.utils.session_drawing_view import drawing_view
from tests.utils.session_files import read_session_payload

REPO_ROOT = Path(__file__).resolve().parents[1]
LAMBDA = REPO_ROOT / "tests" / "test_inputs" / "NC_001416.gb"
DE3 = REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "de3" / "NC_042057.1.gb"
LAMBDA_LENGTH = 48502
DE3_LENGTH = 42925
GENOMES = LAMBDA_LENGTH * DE3_LENGTH
THRESHOLDS = ["--bitscore", "0", "--evalue", "10", "--identity", "0", "--alignment_length", "0", "-f", "svg"]
_TABLE_HEADER = "gbk\trecord_id\trow\tcolumn\tregion\treverse_complement\trecord_group"

_FAKE = """#!{python}
import json, sys
args = sys.argv[1:]
with open({log!r}, "a", encoding="utf-8") as handle:
    handle.write(json.dumps(args) + "\\n")
if args in (["--version"], ["-version"]):
    print("losat {version}")
    sys.exit(0)
if len(args) == 2 and args[1] == "--help":
    print("      -task <TASK>")
    print("          [default: megablast] [possible values: megablast, blastn]")
    if {search_space}:
        print("      -searchsp <INT_VALUE>")
        print("      -dbsize <LETTERS>")
    sys.exit(0)
def ids(path):
    with open(path, encoding="utf-8") as handle:
        return [line[1:].split()[0] for line in handle if line.startswith(">")]
for query in ids(args[args.index("-query") + 1])[:3]:
    for subject in ids(args[args.index("-subject") + 1])[:3]:
        print("\\t".join([query, subject, "95.000", "30", "1", "0", "1", "30", "1", "30", "1e-5", "50"]))
"""


@pytest.fixture(autouse=True)
def _fresh_probe_caches(monkeypatch: pytest.MonkeyPatch) -> None:
    for name in ("_RUNTIME_VERSIONS", "_TASK_VALUES", "_HELP_TEXTS"):
        monkeypatch.setattr(runtime_module, name, {}, raising=False)


def _fake_losat(tmp_path: Path, *, search_space: bool = True) -> tuple[str, Path]:
    log = tmp_path / "fake.log"
    path = tmp_path / "bin" / "losat"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        _FAKE.format(
            python=sys.executable, log=str(log), search_space=search_space,
            version="0.2.0" if search_space else "0.1.0",
        ),
        encoding="utf-8",
    )
    path.chmod(0o755)
    return str(path), log


def _search_spaces(log: Path, flag: str = "-searchsp") -> list[int | None]:
    """The ``flag`` value of each search the fake ran, in order (None for a search without it)."""

    calls = [json.loads(line) for line in log.read_text(encoding="utf-8").splitlines()]
    return [
        int(call[call.index(flag) + 1]) if flag in call else None
        for call in calls if "-query" in call
    ]


def _table(tmp_path: Path, rows: list[str]) -> Path:
    table = tmp_path / "records.tsv"
    table.write_text(_TABLE_HEADER + "\n" + "".join(f"{row}\n" for row in rows), encoding="utf-8")
    return table


def _split_de3(tmp_path: Path) -> tuple[Path, Path]:
    """DE3 written as two files of one record each: one genome split over two files."""

    record = SeqIO.read(DE3, "genbank")
    left, right = record[:20000], record[20000:]
    for part, suffix in ((left, "a"), (right, "b")):
        part.id = part.name = f"DE3{suffix}"
        part.annotations["molecule_type"] = "DNA"
    paths = (tmp_path / "de3a.gb", tmp_path / "de3b.gb")
    SeqIO.write(left, paths[0], "genbank")
    SeqIO.write(right, paths[1], "genbank")
    return paths


def _two_record_de3(tmp_path: Path) -> Path:
    """DE3 written as one file of two records: two replicons of one genome."""

    record = SeqIO.read(DE3, "genbank")
    parts = [record[:20000], record[20000:]]
    for part, suffix in zip(parts, "ab"):
        part.id = part.name = f"DE3{suffix}"
        part.annotations["molecule_type"] = "DNA"
    path = tmp_path / "de3-two.gb"
    SeqIO.write(parts, path, "genbank")
    return path


def _run(tmp_path: Path, argv: list[str], *, program: str = "losatn", search_space: bool = True) -> tuple[list[int | None], dict]:
    losat, log = _fake_losat(tmp_path, search_space=search_space)
    prefix = tmp_path / "out"
    linear_main([*argv, *THRESHOLDS, "--losat", program, "--losat_bin", losat, "--save_session", "-o", str(prefix)])
    session = drawing_view(read_session_payload(prefix.with_suffix(".gbdraw")))
    return _search_spaces(log), session


# The owner.


def test_the_search_space_of_a_record_pair_is_its_genomes_lengths_multiplied() -> None:
    # Rows 0 and 1; row 1 is one genome split over files "b" and "c", and file "c" has a record that is not drawn.
    space = genome_search_spaces(
        rows=[0, 1, 1],
        source_ids=["a", "b", "c"],
        source_lengths=[100, 40, 70],
        searched_lengths=[100, 40, 30],
        letters=(1, 1),
    )
    assert space(0, 1) == space(0, 2) == space(1, 0) == 100 * 110
    translated = genome_search_spaces(
        rows=[0, 1, 1], source_ids=["a", "b", "c"], source_lengths=[100, 40, 70],
        searched_lengths=[100, 40, 30], letters=(3, 3),
    )
    assert translated(0, 1) == (100 // 3) * (110 // 3)


def test_two_whole_records_each_the_whole_of_its_genome_need_no_search_space() -> None:
    space = genome_search_spaces(
        rows=[0, 1, 1], source_ids=["a", "b", "b"], source_lengths=[100, 40, 40],
        searched_lengths=[100, 40, 40], letters=(1, 1),
    )
    # Record 2 draws file "b" again: still the whole genome on both sides.
    assert space(0, 1) is None and space(0, 2) is None
    cropped = genome_search_spaces(
        rows=[0, 1], source_ids=["a", "b"], source_lengths=[100, 40], searched_lengths=[100, 39], letters=(1, 1),
    )
    assert cropped(0, 1) == 100 * 40


def test_a_source_named_with_two_lengths_is_refused() -> None:
    with pytest.raises(ValueError, match="two lengths"):
        genome_search_spaces(
            rows=[0, 1], source_ids=["a", "a"], source_lengths=[100, 90], searched_lengths=[50, 40],
            letters=(1, 1),
        )


def test_a_losatp_search_counts_the_proteins_of_its_subject_genome_as_its_database() -> None:
    # LOSATP's scale is the subject genome's proteins (every record of its files, uncropped), whichever records,
    # crops or windows a search reads; the query side counts each protein as before.
    size = genome_database_sizes(
        rows=[0, 1, 1], source_ids=["a", "b", "c"], source_lengths=[100, 40, 70],
        searched_lengths=[100, 40, 30], source_proteins=[30, 12, 20],
    )
    assert size(0, 1) == size(0, 2) == 12 + 20
    assert size(1, 0) == size(2, 0) == 30
    whole = genome_database_sizes(
        rows=[0, 1], source_ids=["a", "b"], source_lengths=[100, 40], searched_lengths=[100, 40],
        source_proteins=[30, 12],
    )
    assert whole(0, 1) is None and whole(1, 0) is None
    cropped = genome_database_sizes(
        rows=[0, 1], source_ids=["a", "b"], source_lengths=[100, 40], searched_lengths=[100, 39],
        source_proteins=[30, 12],
    )
    assert (cropped(0, 1), cropped(1, 0)) == (12, 30)
    with pytest.raises(ValueError, match="two protein counts"):
        genome_database_sizes(
            rows=[0, 1], source_ids=["a", "a"], source_lengths=[100, 100], searched_lengths=[50, 40],
            source_proteins=[30, 29],
        )


def test_a_source_counts_the_letters_of_every_protein_losatp_would_search_in_it() -> None:
    from gbdraw.analysis.protein_colinearity import (
        extract_web_stable_cds_proteins,
        proteins_to_fasta,
        source_protein_letters,
    )

    de3 = SeqIO.read(DE3, "genbank")
    lambda_record = SeqIO.read(LAMBDA, "genbank")

    def searched_letters(record: SeqIO.SeqRecord) -> int:
        fasta = proteins_to_fasta(extract_web_stable_cds_proteins([record]).proteins_by_record[0])
        return sum(len(line) for line in fasta.splitlines() if not line.startswith(">"))

    assert source_protein_letters([de3]) == searched_letters(de3) > 12000
    # A file of two records counts both, drawn or not.
    assert source_protein_letters([de3, lambda_record]) == searched_letters(de3) + searched_letters(lambda_record)


# The command line.


@pytest.mark.linear
def test_a_whole_record_pair_passes_no_search_space_and_keeps_its_raw_key(tmp_path: Path) -> None:
    spaces, session = _run(tmp_path, ["--gbk", str(LAMBDA), str(DE3)])
    assert spaces == [None]
    entry = session["losatCache"]["entries"][0]
    assert entry["args"] == ["--task", "megablast"]


@pytest.mark.linear
def test_a_cropped_record_searches_in_its_genomes_space(tmp_path: Path) -> None:
    table = _table(tmp_path, [
        f"{LAMBDA}\tNC_001416.1\t1\t1\t1001-30000\tfalse\t",
        f"{DE3}\tNC_042057.1\t2\t1\t\tfalse\t",
    ])
    spaces, session = _run(tmp_path, ["--records_table", str(table)])
    assert spaces == [GENOMES]
    assert session["losatCache"]["entries"][0]["args"] == ["--task", "megablast", "--searchsp", str(GENOMES)]


@pytest.mark.linear
def test_each_window_of_a_record_searches_in_its_genomes_space(tmp_path: Path) -> None:
    table = _table(tmp_path, [
        f"{LAMBDA}\tNC_001416.1\t1\t1\t1001-10000\tfalse\tlambda",
        f"{LAMBDA}\tNC_001416.1\t1\t2\t30001-40000\tfalse\tlambda",
        f"{DE3}\tNC_042057.1\t2\t1\t\tfalse\t",
    ])
    spaces, _session = _run(tmp_path, ["--records_table", str(table)])
    # Both windows are queries of one search of their file.
    assert spaces == [GENOMES]


@pytest.mark.linear
def test_a_genome_split_over_two_files_in_one_row_searches_in_one_space(tmp_path: Path) -> None:
    left, right = _split_de3(tmp_path)
    table = _table(tmp_path, [
        f"{LAMBDA}\tNC_001416.1\t1\t1\t\tfalse\t",
        f"{left}\tDE3a\t2\t1\t\tfalse\t",
        f"{right}\tDE3b\t2\t2\t\tfalse\t",
    ])
    spaces, _session = _run(tmp_path, ["--records_table", str(table)])
    assert spaces == [GENOMES, GENOMES]


@pytest.mark.linear
def test_the_records_of_one_file_in_one_row_search_in_their_files_space(tmp_path: Path) -> None:
    two = _two_record_de3(tmp_path)
    # Only the first record of the file is drawn: the second still counts in its genome.
    table = _table(tmp_path, [
        f"{LAMBDA}\tNC_001416.1\t1\t1\t\tfalse\t",
        f"{two}\tDE3a\t2\t1\t\tfalse\t",
    ])
    spaces, _session = _run(tmp_path, ["--records_table", str(table)])
    assert spaces == [GENOMES]


@pytest.mark.linear
def test_tlosatx_counts_its_search_space_in_translated_letters(tmp_path: Path) -> None:
    table = _table(tmp_path, [
        f"{LAMBDA}\tNC_001416.1\t1\t1\t1001-30000\tfalse\t",
        f"{DE3}\tNC_042057.1\t2\t1\t\tfalse\t",
    ])
    spaces, _session = _run(tmp_path, ["--records_table", str(table)], program="tlosatx")
    assert spaces == [(LAMBDA_LENGTH // 3) * (DE3_LENGTH // 3)]


@pytest.mark.linear
def test_a_losat_without_a_search_space_option_searches_as_before_and_says_so(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    table = _table(tmp_path, [
        f"{LAMBDA}\tNC_001416.1\t1\t1\t1001-30000\tfalse\t",
        f"{DE3}\tNC_042057.1\t2\t1\t\tfalse\t",
    ])
    with caplog.at_level("WARNING"):
        spaces, session = _run(tmp_path, ["--records_table", str(table)], search_space=False)
    assert spaces == [None]
    assert session["losatCache"]["entries"][0]["args"] == ["--task", "megablast"]
    assert "-searchsp" in caplog.text and "-dbsize" in caplog.text


# One E-value whatever the view (NCBI BLAST+, when installed).


def _evalues(tmp_path: Path, name: str, rows: list[str]) -> dict[tuple[str, ...], str]:
    """The E-value of each hit, by its subject-frame coordinates, in the raw rows of the searches."""

    work = tmp_path / name
    work.mkdir()
    table = work / "records.tsv"
    table.write_text(_TABLE_HEADER + "\n" + "".join(f"{row}\n" for row in rows), encoding="utf-8")
    raw = work / "raw"
    linear_main([
        "--records_table", str(table), *THRESHOLDS, "--losat", "losatn",
        "--ncbi_blast_bin", str(shutil.which("blastn")), "--losat_output_dir", str(raw), "-o", str(work / "out"),
    ])
    hits: dict[tuple[str, ...], str] = {}
    for path in sorted(raw.glob("*.tsv")):
        if path.name == "comparisons.tsv":
            continue
        for line in path.read_text(encoding="utf-8").splitlines():
            columns = line.split("\t")
            # Bit score and length name a hit in every view; the E-value is what the view must not change.
            hits[(columns[3], columns[11])] = columns[10]
    return hits


@pytest.mark.linear
@pytest.mark.skipif(shutil.which("blastn") is None, reason="NCBI BLAST+ blastn is not installed.")
def test_a_hit_has_one_evalue_whether_its_genome_is_cropped_windowed_or_split(tmp_path: Path) -> None:
    left, right = _split_de3(tmp_path)
    lam = f"{LAMBDA}\tNC_001416.1\t1\t1\t\tfalse\t"
    cropped = _evalues(tmp_path, "cropped", [lam, f"{DE3}\tNC_042057.1\t2\t1\t1-30000\tfalse\t"])
    windows = _evalues(tmp_path, "windows", [
        lam,
        f"{DE3}\tNC_042057.1\t2\t1\t1-15000\tfalse\tde3",
        f"{DE3}\tNC_042057.1\t2\t2\t15001-30000\tfalse\tde3",
    ])
    split = _evalues(tmp_path, "split", [
        lam, f"{left}\tDE3a\t2\t1\t\tfalse\t", f"{right}\tDE3b\t2\t2\t\tfalse\t",
    ])
    shared = set(cropped) & set(windows) & set(split)
    assert len(shared) >= 3
    for hit in shared:
        assert cropped[hit] == windows[hit] == split[hit], hit


# A genome is the files linked through the rows they share.


@pytest.mark.linear
def test_files_linked_through_shared_rows_are_one_genome_with_one_search(tmp_path: Path) -> None:
    """File F has one record beside file H in row 1 and one alone in row 3; file G is row 2. F and H are one genome,
    so both F records search G in one space, and one job searches F against G, as it did before genomes."""

    two = _two_record_de3(tmp_path)
    g_file, _right = _split_de3(tmp_path)
    table = _table(tmp_path, [
        f"{two}\tDE3a\t3\t1\t\tfalse\t",
        f"{two}\tDE3b\t1\t1\t\tfalse\t",
        f"{LAMBDA}\tNC_001416.1\t1\t2\t\tfalse\t",
        f"{g_file}\tDE3a\t2\t1\t\tfalse\t",
    ])
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("blast\tquery\tsubject\tsource\n\t#1\t#4\tlosat\n\t#2\t#4\tlosat\n", encoding="utf-8")
    spaces, _session = _run(tmp_path, ["--records_table", str(table), "--comparisons_table", str(pairs)])
    assert spaces == [(DE3_LENGTH + LAMBDA_LENGTH) * 20000]


# The length of a GFF3/FASTA source is every FASTA entry.


def _gff_naming_one_of_two_entries(tmp_path: Path) -> tuple[Path, Path]:
    record = SeqIO.read(DE3, "genbank")
    fasta = tmp_path / "two.fa"
    fasta.write_text(f">seqA\n{record.seq[:30000]}\n>seqB\n{record.seq[30000:]}\n", encoding="utf-8")
    gff = tmp_path / "one.gff3"
    gff.write_text(
        "##gff-version 3\n##sequence-region seqA 1 30000\n"
        "seqA\ttest\tCDS\t101\t400\t.\t+\t0\tID=cds1;locus_tag=A1\n",
        encoding="utf-8",
    )
    return gff, fasta


@pytest.mark.linear
def test_a_gff3_fasta_source_counts_every_fasta_entry_in_its_genome(tmp_path: Path) -> None:
    gff, fasta = _gff_naming_one_of_two_entries(tmp_path)
    other_gff = REPO_ROOT / "tests" / "test_inputs" / "NC_013668.gff3"
    other_fasta = REPO_ROOT / "tests" / "test_inputs" / "NC_013668.fasta"
    other_length = len(SeqIO.read(other_fasta, "fasta"))
    spaces, _session = _run(tmp_path, ["--gff", str(gff), str(other_gff), "--fasta", str(fasta), str(other_fasta)])
    # The GFF3 file names seqA only; its genome is seqA and seqB, the whole of DE3.
    assert spaces == [DE3_LENGTH * other_length]


# LOSATP counts the proteins of the subject genome as its database (-dbsize).


def _letters(*paths: Path) -> int:
    """The protein letters LOSATP searches in every record of ``paths`` (the owner, with no visibility rules)."""

    from gbdraw.analysis.protein_colinearity import source_protein_letters

    return sum(source_protein_letters(list(SeqIO.parse(path, "genbank"))) for path in paths)


def _run_losatp(
    tmp_path: Path, argv: list[str], *, search_space: bool = True, mode: str = "pairwise"
) -> tuple[list[int | None], dict]:
    losat, log = _fake_losat(tmp_path, search_space=search_space)
    prefix = tmp_path / "out"
    linear_main([
        *argv, *THRESHOLDS, "--losat", "losatp", "--losatp_mode", mode, "--losat_bin", losat,
        "--save_session", "-o", str(prefix),
    ])
    session = drawing_view(read_session_payload(prefix.with_suffix(".gbdraw")))
    return _search_spaces(log, "-dbsize"), session


@pytest.mark.linear
def test_a_whole_losatp_pair_passes_no_database_size_and_keeps_its_raw_key(tmp_path: Path) -> None:
    sizes, session = _run_losatp(tmp_path, ["--gbk", str(LAMBDA), str(DE3)])
    assert sizes == [None]
    assert all("--dbsize" not in entry["args"] for entry in session["losatCache"]["entries"])


@pytest.mark.linear
def test_a_cropped_losatp_subject_counts_every_protein_of_its_whole_genome(tmp_path: Path) -> None:
    table = _table(tmp_path, [
        f"{LAMBDA}\tNC_001416.1\t1\t1\t\tfalse\t",
        f"{DE3}\tNC_042057.1\t2\t1\t1-20000\tfalse\t",
    ])
    sizes, session = _run_losatp(tmp_path, ["--records_table", str(table)])
    assert sizes == [_letters(DE3)]
    args = session["losatCache"]["entries"][0]["args"]
    assert args[args.index("--dbsize") + 1] == str(_letters(DE3))


@pytest.mark.linear
@pytest.mark.parametrize("mode", ["similarity_groups", "collinear"])
def test_similarity_group_and_collinear_searches_count_their_subject_genomes_proteins(
    tmp_path: Path, mode: str
) -> None:
    table = _table(tmp_path, [
        f"{LAMBDA}\tNC_001416.1\t1\t1\t\tfalse\t",
        f"{DE3}\tNC_042057.1\t2\t1\t1-20000\tfalse\t",
    ])
    sizes, _session = _run_losatp(tmp_path, ["--records_table", str(table)], mode=mode)
    # Lambda against the cropped DE3 counts DE3's proteins, and DE3 against lambda lambda's (a self search of
    # the whole lambda, when the mode makes one, passes none).
    assert {_letters(DE3), _letters(LAMBDA)} <= set(sizes) <= {_letters(DE3), _letters(LAMBDA), None}


@pytest.mark.linear
def test_a_losatp_subject_split_over_two_files_in_one_row_is_one_database(tmp_path: Path) -> None:
    left, right = _split_de3(tmp_path)
    table = _table(tmp_path, [
        f"{LAMBDA}\tNC_001416.1\t1\t1\t\tfalse\t",
        f"{left}\tDE3a\t2\t1\t\tfalse\t",
        f"{right}\tDE3b\t2\t2\t\tfalse\t",
    ])
    sizes, _session = _run_losatp(tmp_path, ["--records_table", str(table)])
    assert sizes == [_letters(left, right)] * 2


@pytest.mark.linear
def test_an_undrawn_record_of_a_losatp_subject_file_counts_in_its_database(tmp_path: Path) -> None:
    two = _two_record_de3(tmp_path)
    table = _table(tmp_path, [
        f"{LAMBDA}\tNC_001416.1\t1\t1\t\tfalse\t",
        f"{two}\tDE3a\t2\t1\t\tfalse\t",
    ])
    sizes, _session = _run_losatp(tmp_path, ["--records_table", str(table)])
    assert sizes == [_letters(two)]
    # The drawn record DE3a alone has fewer.
    assert _letters(two) > _letters_of_first(two)


def _letters_of_first(path: Path) -> int:
    from gbdraw.analysis.protein_colinearity import source_protein_letters

    return source_protein_letters([next(SeqIO.parse(path, "genbank"))])


@pytest.mark.linear
def test_a_losat_without_a_database_size_option_searches_losatp_as_before_and_says_so(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    table = _table(tmp_path, [
        f"{LAMBDA}\tNC_001416.1\t1\t1\t\tfalse\t",
        f"{DE3}\tNC_042057.1\t2\t1\t1-20000\tfalse\t",
    ])
    with caplog.at_level("WARNING"):
        sizes, session = _run_losatp(tmp_path, ["--records_table", str(table)], search_space=False)
    assert sizes == [None]
    assert all("--dbsize" not in entry["args"] for entry in session["losatCache"]["entries"])
    assert "-dbsize" in caplog.text and "LOSATP" in caplog.text


def _protein_evalues(tmp_path: Path, name: str, rows: list[str]) -> dict[tuple[str, ...], str]:
    """The E-value of each protein hit of ``losatp.raw.tsv``, by its proteins, length and bit score."""

    work = tmp_path / name
    work.mkdir()
    table = work / "records.tsv"
    table.write_text(_TABLE_HEADER + "\n" + "".join(f"{row}\n" for row in rows), encoding="utf-8")
    raw = work / "raw"
    linear_main([
        "--records_table", str(table), *THRESHOLDS, "--losat", "losatp", "--losatp_mode", "pairwise",
        "--ncbi_blast_bin", str(shutil.which("blastp")), "--losat_output_dir", str(raw), "-o", str(work / "out"),
    ])
    hits: dict[tuple[str, ...], str] = {}
    for line in (raw / "losatp.raw.tsv").read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        columns = line.split("\t")
        hits[(columns[0], columns[1], columns[3], columns[11])] = columns[10]
    return hits


@pytest.mark.linear
@pytest.mark.skipif(shutil.which("blastp") is None, reason="NCBI BLAST+ blastp is not installed.")
def test_a_protein_hit_has_one_evalue_whether_its_subject_is_cropped_windowed_or_split(tmp_path: Path) -> None:
    left, right = _split_de3(tmp_path)
    lam = f"{LAMBDA}\tNC_001416.1\t1\t1\t\tfalse\t"
    cropped = _protein_evalues(tmp_path, "cropped", [lam, f"{DE3}\tNC_042057.1\t2\t1\t1-30000\tfalse\t"])
    windows = _protein_evalues(tmp_path, "windows", [
        lam,
        f"{DE3}\tNC_042057.1\t2\t1\t1-15000\tfalse\tde3",
        f"{DE3}\tNC_042057.1\t2\t2\t15001-30000\tfalse\tde3",
    ])
    split = _protein_evalues(tmp_path, "split", [
        lam, f"{left}\tDE3a\t2\t1\t\tfalse\t", f"{right}\tDE3b\t2\t2\t\tfalse\t",
    ])
    shared = set(cropped) & set(windows) & set(split)
    assert len(shared) >= 3
    for hit in shared:
        assert cropped[hit] == windows[hit] == split[hit], hit


def test_a_losatp_cache_that_answers_every_pair_needs_no_runtime(tmp_path: Path) -> None:
    """The Web's Generate reads LOSATP results from its cache with no LOSAT to probe: a sized search whose every
    pair the cache holds must not reach the runtime."""

    from gbdraw.analysis.protein_colinearity import (
        LosatpCacheManager,
        _search_losatp_record_pairs,
        extract_web_stable_cds_proteins,
    )
    from gbdraw.comparisons.losat_jobs import GenomeDatabaseSizes

    records = [SeqIO.read(LAMBDA, "genbank"), SeqIO.read(DE3, "genbank")]
    extraction = extract_web_stable_cds_proteins(records)
    sizes = GenomeDatabaseSizes(proteins=(100, 200), whole=(False, False))

    def search(cache: LosatpCacheManager, losat_bin: str) -> dict:
        return _search_losatp_record_pairs(
            records, extraction, [(0, 1)], losatp_bin=losat_bin, ncbi_blastp_bin=None, losatp_threads=None,
            candidate_limit=None, max_hsps_per_subject=1, runner=None, losatp_cache=cache,
        )

    losat, log = _fake_losat(tmp_path)
    first = LosatpCacheManager(identity_manifest=extraction.identity_manifest, database_sizes=sizes)
    searched = search(first, losat)
    assert _search_spaces(log, "-dbsize") == [200]
    replay = LosatpCacheManager(
        first.session_entries(), identity_manifest=extraction.identity_manifest, database_sizes=sizes
    )
    replayed = search(replay, str(tmp_path / "no-such-losat"))
    assert replayed[(0, 1)].equals(searched[(0, 1)])


# The Web passes the scale Python plans for the drawing's own request (PD-OI-018 revision 6, receipt f):
# ``plan_genome_scales`` gives each drawn-unit pair the -searchsp or -dbsize the command line's search passes.


def _captured_cli_request(monkeypatch: pytest.MonkeyPatch) -> list:
    """The requests the command line plans, in order: the first is the figure's own."""

    import gbdraw.api.request_render as request_render

    captured: list = []
    planner = request_render.plan_linear_request

    def capture(request):  # type: ignore[no-untyped-def]
        captured.append(request)
        return planner(request)

    monkeypatch.setattr(request_render, "plan_linear_request", capture)
    return captured


def web_genome_scales(request, directory: Path) -> dict:
    """The reply the Web's Generate gets for ``request``: the request encoded and its resources staged, as the page
    sends a drawing's draft request. A request holds no LOSATN or TLOSATX search (the Web resolves them before it
    renders), so the page names the program."""

    from dataclasses import replace

    from gbdraw.session_request_codec import encode_canonical_request
    from gbdraw.web_support.linear_units import plan_genome_scales

    program = request.options.losat_search.program
    if program != "losatp":
        request = replace(request, options=replace(request.options, losat_search=None))
    encoded = encode_canonical_request(request)
    staged = directory / "staged"
    staged.mkdir(parents=True, exist_ok=True)
    paths = {}
    for index, resource in enumerate(encoded.resources, start=1):
        target = staged / f"{index:04d}.bin"
        target.write_bytes(
            resource.content if resource.content is not None else Path(str(resource.source_path)).read_bytes()
        )
        paths[resource.resource_id] = str(target)
    names = {resource.resource_id: resource.name for resource in encoded.resources}
    return plan_genome_scales(encoded.payload, paths, str(directory / "web-out"), program, names)


_LAM = f"{LAMBDA}\tNC_001416.1\t1\t1\t\tfalse\t"


_DE3_WINDOWS = [
    f"{DE3}\tNC_042057.1\t2\t1\t1-15000\tfalse\tde3", f"{DE3}\tNC_042057.1\t2\t2\t15001-30000\tfalse\tde3",
]
_DE3_WHOLE_AND_WINDOW = [
    f"{DE3}\tNC_042057.1\t2\t1\t1-42925\tfalse\tde3", f"{DE3}\tNC_042057.1\t2\t2\t1-15000\tfalse\tde3",
]
_DE3_CROP = f"{DE3}\tNC_042057.1\t2\t1\t1-20000\tfalse\t"


@pytest.mark.linear
@pytest.mark.parametrize(
    ("program", "mode", "rows", "jobs"),
    [
        ("losatn", None, [_LAM, f"{DE3}\tNC_042057.1\t2\t1\t\tfalse\t"], [[(0, 1)]]),
        ("losatn", None, [_LAM, _DE3_CROP], [[(0, 1)]]),
        ("losatn", None, [_LAM, *_DE3_WINDOWS], [[(0, 1), (0, 2)]]),
        ("losatn", None, [_LAM, *_DE3_WHOLE_AND_WINDOW], [[(0, 1), (0, 2)]]),
        ("tlosatx", None, [_LAM, _DE3_CROP], [[(0, 1)]]),
        ("losatp", "pairwise", [_LAM, f"{DE3}\tNC_042057.1\t2\t1\t\tfalse\t"], [[(0, 1)]]),
        ("losatp", "pairwise", [_LAM, _DE3_CROP], [[(0, 1)]]),
        ("losatp", "pairwise", [_LAM, *_DE3_WINDOWS], [[(0, 1), (0, 2)]]),
        ("losatp", "pairwise", [_LAM, *_DE3_WHOLE_AND_WINDOW], [[(0, 1), (0, 2)]]),
        ("losatp", "similarity_groups", [_LAM, _DE3_CROP], [[(0, 0)], [(0, 1)], [(1, 0)], [(1, 1)]]),
        ("losatp", "collinear", [_LAM, _DE3_CROP], [[(0, 1)], [(1, 0)]]),
    ],
    ids=[
        "losatn-whole", "losatn-crop", "losatn-windows", "losatn-whole-and-window", "tlosatx-crop", "losatp-whole",
        "losatp-crop", "losatp-windows", "losatp-whole-and-window", "losatp-similarity-groups", "losatp-collinear",
    ],
)
def test_python_plans_the_web_the_scale_the_command_line_searches_with(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, program: str, mode: str | None, rows: list[str],
    jobs: list[list[tuple[int, int]]],
) -> None:
    """Each search the command line runs, in its order, passes the scale the Web's reply gives each of its unit
    pairs, as the decimal string the Web passes on."""

    captured = _captured_cli_request(monkeypatch)
    table = ["--records_table", str(_table(tmp_path, rows))]
    if program == "losatp":
        values, _session = _run_losatp(tmp_path, table, mode=str(mode))
    else:
        values, _session = _run(tmp_path, table, program=program)
    reply = web_genome_scales(captured[0], tmp_path)
    scales = reply["scales"]
    assert len(scales) == len(reply["recordKeys"]) == len(rows)
    assert all(len(row) == len(rows) for row in scales)
    assert [{scales[query][subject] for query, subject in job} for job in jobs] == [
        {None if value is None else str(value)} for value in values
    ]


@pytest.mark.linear
def test_the_web_counts_every_fasta_entry_of_a_gff3_fasta_genome(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """OV-3250: the Web's scale counts the FASTA entry the GFF3 file does not name, as the command line does."""

    captured = _captured_cli_request(monkeypatch)
    gff, fasta = _gff_naming_one_of_two_entries(tmp_path)
    other_gff = REPO_ROOT / "tests" / "test_inputs" / "NC_013668.gff3"
    other_fasta = REPO_ROOT / "tests" / "test_inputs" / "NC_013668.fasta"
    spaces, _session = _run(tmp_path, ["--gff", str(gff), str(other_gff), "--fasta", str(fasta), str(other_fasta)])
    reply = web_genome_scales(captured[0], tmp_path)
    assert reply["scales"][0][1] == str(spaces[0]) == str(DE3_LENGTH * len(SeqIO.read(other_fasta, "fasta")))


def _cds_in(record, start: int, end: int):  # type: ignore[no-untyped-def]
    return next(
        feature for feature in record.features
        if feature.type == "CDS" and start <= int(feature.location.start) and int(feature.location.end) <= end
    )


@pytest.mark.linear
def test_an_edit_in_the_second_window_of_a_losatp_subject_counts_in_its_database() -> None:
    """OV-3264: the edits of every window of a record count in its genome's proteins, not the first window's only."""

    from gbdraw.analysis.protein_colinearity import extract_cds_proteins, source_protein_letters
    from gbdraw.api import (
        GenBankInputSource, LinearDiagramOptions, LinearDiagramRequest, LinearMultiRecordOptions, RecordInput,
        RecordPresentation, plan_linear_request,
    )
    from gbdraw.api.options import LosatSearchOptions
    from gbdraw.features.ids import compute_source_feature_hash
    from gbdraw.features.overrides import FeatureOverride
    from gbdraw.io.regions import RecordWindow

    de3 = SeqIO.read(DE3, "genbank")
    hidden = _cds_in(de3, 25000, 42925)
    request = LinearDiagramRequest(
        records=(
            RecordInput(GenBankInputSource(LAMBDA), record_key="lambda", presentation=RecordPresentation(grid_row=1)),
            RecordInput(
                GenBankInputSource(DE3), record_key="de3",
                windows=(RecordWindow(1, 15000), RecordWindow(25001, 42925)),
                presentation=RecordPresentation(grid_row=2),
            ),
        ),
        layout=LinearMultiRecordOptions(),
        options=LinearDiagramOptions(
            losat_search=LosatSearchOptions(program="losatp", losatp_mode="pairwise"),
            feature_overrides=(
                FeatureOverride("de3", str(compute_source_feature_hash(hidden, de3)), feature_visibility="off"),
            ),
        ),
    )
    sizes = plan_linear_request(request).losatp_database_sizes
    assert sizes is not None
    [proteins] = extract_cds_proteins([de3], prefer_source_ids=False).proteins_by_record
    hidden_letters = sum(len(protein.sequence) for protein in proteins if protein.feature_index == de3.features.index(hidden))
    assert hidden_letters > 0
    assert sizes(0, 1) == sizes(0, 2) == source_protein_letters([de3]) - hidden_letters


def _de3_request(windows, *, hidden=None, de3_source=None):  # type: ignore[no-untyped-def]
    """Lambda in row 1 and DE3 (``windows``) in row 2, LOSATP Pairwise, optionally hiding the CDS ``hidden``."""

    from gbdraw.api import (
        GenBankInputSource, LinearDiagramOptions, LinearDiagramRequest, LinearMultiRecordOptions, RecordInput,
        RecordPresentation,
    )
    from gbdraw.api.options import LosatSearchOptions
    from gbdraw.features.ids import compute_source_feature_hash
    from gbdraw.features.overrides import FeatureOverride

    de3 = SeqIO.read(DE3, "genbank")
    edits = () if hidden is None else (
        FeatureOverride("de3", str(compute_source_feature_hash(hidden, de3)), feature_visibility="off"),
    )
    return LinearDiagramRequest(
        records=(
            RecordInput(GenBankInputSource(LAMBDA), record_key="lambda", presentation=RecordPresentation(grid_row=1)),
            RecordInput(
                de3_source or GenBankInputSource(DE3), record_key="de3", windows=windows,
                presentation=RecordPresentation(grid_row=2),
            ),
        ),
        layout=LinearMultiRecordOptions(),
        options=LinearDiagramOptions(
            losat_search=LosatSearchOptions(program="losatp", losatp_mode="pairwise"), feature_overrides=edits,
        ),
    )


@pytest.mark.linear
def test_an_edit_outside_the_crop_of_a_losatp_subject_counts_in_its_database() -> None:
    """A hidden CDS the crop leaves out still leaves the genome's proteins: the crop does not change the scale."""

    from gbdraw.analysis.protein_colinearity import extract_cds_proteins, source_protein_letters
    from gbdraw.api import plan_linear_request
    from gbdraw.io.regions import RecordWindow

    de3 = SeqIO.read(DE3, "genbank")
    hidden = _cds_in(de3, 30000, 42925)
    [proteins] = extract_cds_proteins([de3], prefer_source_ids=False).proteins_by_record
    expected = source_protein_letters([de3]) - sum(
        len(protein.sequence) for protein in proteins if protein.feature_index == de3.features.index(hidden)
    )
    assert expected < source_protein_letters([de3])
    outside = plan_linear_request(_de3_request((RecordWindow(1, 20000),), hidden=hidden)).losatp_database_sizes
    inside = plan_linear_request(_de3_request((RecordWindow(1, 40000),), hidden=hidden)).losatp_database_sizes
    assert outside is not None and inside is not None
    assert outside(0, 1) == inside(0, 1) == expected


@pytest.mark.linear
def test_a_losatp_subject_cropped_in_memory_counts_the_proteins_of_the_record_it_was_cut_from() -> None:
    """An in-memory crop counts its uncropped record's proteins, as its genome length is the uncropped record's."""

    from gbdraw.analysis.protein_colinearity import source_protein_letters
    from gbdraw.api import InMemoryRecordSource, plan_linear_request
    from gbdraw.io.regions import crop_record_to_region, parse_region_spec

    de3 = SeqIO.read(DE3, "genbank")
    crop = crop_record_to_region(de3, parse_region_spec(f"{de3.id}:1-20000"))
    sizes = plan_linear_request(_de3_request((), de3_source=InMemoryRecordSource(crop))).losatp_database_sizes
    assert sizes is not None
    assert sizes(0, 1) == source_protein_letters([de3])


@pytest.mark.linear
@pytest.mark.parametrize("program", ["losatn", "losatp"])
def test_a_record_drawn_whole_and_as_a_window_is_one_search_with_its_genome_scale(
    tmp_path: Path, program: str
) -> None:
    """A genome is searched whole only when every drawn unit of it reads all of it: one window that reads less gives
    every search of the genome the scale, so the source batch stays one job, as before."""

    table = ["--records_table", str(_table(tmp_path, [
        _LAM, f"{DE3}\tNC_042057.1\t2\t1\t1-42925\tfalse\tde3", f"{DE3}\tNC_042057.1\t2\t2\t1-15000\tfalse\tde3",
    ]))]
    if program == "losatp":
        values, _session = _run_losatp(tmp_path, table)
        assert values == [_letters(DE3)]
    else:
        values, _session = _run(tmp_path, table, program=program)
        assert values == [GENOMES]


# A genome pair past what LOSAT reads (Owner decision 2026-11-10): searched per file, with a warning.

_LIMIT = 2**63 - 1
# 2**63 - 1 = 7 * 7 * 73 * 127 * 337 * 92737 * 649657.
_AT_LIMIT = (7 * 7 * 73 * 127 * 337, 92737 * 649657)


@pytest.mark.parametrize(("lengths", "fits"), [(_AT_LIMIT, True), ((2**31, 2**32), False)])
def test_a_genome_pair_past_the_largest_search_space_losat_reads_is_searched_per_file(
    lengths: tuple[int, int], fits: bool
) -> None:
    from gbdraw.comparisons.losat_jobs import LOSAT_SEARCH_SPACE_LIMIT

    assert LOSAT_SEARCH_SPACE_LIMIT == _LIMIT
    space = genome_search_spaces(
        rows=[0, 1], source_ids=["a", "b"], source_lengths=list(lengths), searched_lengths=[10, 10], letters=(1, 1),
    )
    product = lengths[0] * lengths[1]
    assert product == (_LIMIT if fits else _LIMIT + 1)
    warning = space.per_file_warning(0, 1, ["a.gb", "b.gb"])
    if fits:
        assert space(0, 1) == product and warning is None
    else:
        assert space(0, 1) is None and space(1, 0) is None
        assert warning == space.per_file_warning(1, 0, ["a.gb", "b.gb"])
        assert warning is not None and warning.message() == (
            "The search space of genomes a.gb and b.gb (9,223,372,036,854,775,808) is past the largest LOSAT reads "
            "(9,223,372,036,854,775,807), so their E-values are computed per file, not per genome."
        )


def test_the_command_line_warns_once_for_a_genome_pair_searched_per_file(caplog: pytest.LogCaptureFixture) -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from gbdraw.api.options import LosatSearchOptions
    from gbdraw.comparisons.linear_losat import _genome_search_space

    with caplog.at_level("WARNING"):
        space = _genome_search_space(
            LosatSearchOptions(program="losatn"),
            records=[SeqRecord(Seq("ACGT"), id="a"), SeqRecord(Seq("ACGT"), id="b")],
            rows_by_record=[0, 1],
            source_ids=[("/inputs/big-a.gb",), ("/inputs/big-b.gb",)],
            source_lengths=[2**31, 2**32],
            pairs=[(0, 1), (1, 0)],
        )
    assert space(0, 1) is None
    warnings = [record.getMessage() for record in caplog.records if record.levelname == "WARNING"]
    assert len(warnings) == 1
    assert "big-a.gb" in warnings[0] and "big-b.gb" in warnings[0] and "per file, not per genome" in warnings[0]


@pytest.mark.linear
def test_the_web_gets_the_per_file_warning_the_command_line_logs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    from gbdraw.api import request_render

    captured = _captured_cli_request(monkeypatch)
    monkeypatch.setattr(request_render, "_source_genome_lengths", lambda collection: [2**31, 2**32])
    with caplog.at_level("WARNING"):
        values, _session = _run(tmp_path, ["--records_table", str(_table(tmp_path, [_LAM, _DE3_CROP]))])
    assert values == [None]
    logged = [record.getMessage() for record in caplog.records if "per file, not per genome" in record.getMessage()]
    reply = web_genome_scales(captured[0], tmp_path)
    assert reply["scales"][0][1] is None
    # The reply names every unit pair past the limit (DE3 against itself too); the page shows those it searches.
    # The page shows `tm(text, params)`: the English text with its parameters is byte for byte the logged one.
    from gbdraw.comparisons.losat_jobs import PER_FILE_WARNING_TEXT

    assert {item["text"] for item in reply["warnings"]} == {PER_FILE_WARNING_TEXT}
    searched = {
        (item["queryUnit"], item["subjectUnit"]): item["text"].format(**item["params"]) for item in reply["warnings"]
    }
    assert {(0, 1), (1, 0), (1, 1)} == set(searched)
    assert {searched[(0, 1)], searched[(1, 0)]} == set(logged)
    assert "NC_001416.gb" in logged[0] and "NC_042057.1.gb" in logged[0]
