"""CLI and Python LOSATP search the Web source-file databases (design D7, PR-5).

A source file is one genome: records of one file are searched together, a
record never searches itself unless that search is requested, and raw keys
carry the Web ``searchContext`` when a side holds several records.
"""

from __future__ import annotations

from tests.utils.session_files import read_session_payload
import hashlib
import json
from pathlib import Path
import shutil
import sys

import pytest

import gbdraw.analysis.protein_colinearity as protein_colinearity
import gbdraw.comparisons.losat_runtime as runtime_module
from gbdraw.analysis.protein_colinearity import validate_protein_raw_entry_references
from gbdraw.linear import linear_main
from tests.utils.session_drawing_view import drawing_view

REPO_ROOT = Path(__file__).resolve().parents[1]
INPUTS = REPO_ROOT / "tests" / "test_inputs"

_FAKE = """#!{python}
import json, sys
args = sys.argv[1:]
if args in (["--version"], ["-version"]):
    print("losat 0.1.0")
    sys.exit(0)
if len(args) == 2 and args[1] == "--help":
    sys.exit(0)
def read(path):
    with open(path, encoding="utf-8") as handle:
        return handle.read()
query = read(args[args.index("-query") + 1])
subject = read(args[args.index("-subject") + 1])
with open({log!r}, "a", encoding="utf-8") as handle:
    handle.write(json.dumps({{"args": args, "query": query, "subject": subject}}) + "\\n")
def ids(text):
    return [line[1:].split()[0] for line in text.splitlines() if line.startswith(">")]
for query_id in ids(query):
    for subject_id in ids(subject):
        print("\\t".join([query_id, subject_id, "80.000", "10", "2", "0", "1", "10",
                         "1", "10", "1e-20", "90"]))
"""


@pytest.fixture(autouse=True)
def _fresh_probe_caches(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(runtime_module, "_RUNTIME_VERSIONS", {})
    monkeypatch.setattr(runtime_module, "_TASK_VALUES", {}, raising=False)


def _fake_losat(tmp_path: Path) -> tuple[str, Path]:
    log = tmp_path / "searches.jsonl"
    path = tmp_path / "bin" / "losat"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_FAKE.format(python=sys.executable, log=str(log)), encoding="utf-8")
    path.chmod(0o755)
    return str(path), log


def _searches(log: Path) -> list[dict]:
    return [json.loads(line) for line in log.read_text(encoding="utf-8").splitlines()]


def _fasta_ids(text: str) -> list[str]:
    return [line[1:].split()[0] for line in text.splitlines() if line.startswith(">")]


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _inputs(tmp_path: Path) -> Path:
    """Two BGC records packaged in one file, and a third in its own file."""

    (tmp_path / "two.gbk").write_text(
        (INPUTS / "BGC0000708.gbk").read_text(encoding="utf-8")
        + (INPUTS / "BGC0000709.gbk").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    shutil.copy(INPUTS / "BGC0000711.gbk", tmp_path / "BGC0000711.gbk")
    for name in ("BGC0000708.gbk", "BGC0000709.gbk"):
        shutil.copy(INPUTS / name, tmp_path / name)
    table = tmp_path / "records.tsv"
    table.write_text(
        "gbk\trecord_id\n"
        f"{tmp_path / 'two.gbk'}\tBGC0000708\n"
        f"{tmp_path / 'two.gbk'}\tBGC0000709\n"
        f"{tmp_path / 'BGC0000711.gbk'}\tBGC0000711\n",
        encoding="utf-8",
    )
    return table


def _run(tmp_path: Path, name: str, inputs: list[str], mode: str) -> tuple[dict, list[dict]]:
    losat, log = _fake_losat(tmp_path / name)
    prefix = tmp_path / name / "out"
    linear_main([
        *inputs, "--losat", "losatp", "--losatp_mode", mode, "--losat_bin", losat,
        "--losat_threads", "1", "-o", str(prefix), "-f", "svg", "--save_session",
    ])
    session = drawing_view(read_session_payload(prefix.with_suffix(".gbdraw")))
    return session, _searches(log)


def _entries_by_pair(session: dict) -> dict[tuple[str, str], dict]:
    return {
        (entry["queryRecordInstanceKey"], entry["subjectRecordInstanceKey"]): entry
        for entry in session["losatCache"]["entries"]
    }


@pytest.mark.linear
def test_records_of_one_file_are_one_database(tmp_path: Path) -> None:
    table = _inputs(tmp_path)
    session, searches = _run(tmp_path, "groups", ["--records_table", str(table)], "similarity_groups")

    # 3 self searches, 2 within the file, 2 between the files (Web plan),
    # not 9 record-pair searches.
    assert len(searches) == 7
    pairs = [
        (frozenset(_fasta_ids(search["query"])), frozenset(_fasta_ids(search["subject"])))
        for search in searches
    ]
    proteins = [query for query, subject in pairs if query == subject]
    assert len(proteins) == 3
    first, second, third = proteins
    # A record never searches itself within its file; between the files, each
    # file is one database.
    assert sorted(
        (sorted(query), sorted(subject)) for query, subject in pairs if query != subject
    ) == sorted(
        (sorted(query), sorted(subject))
        for query, subject in (
            (first, second), (second, first),
            (first | second, third), (third, first | second),
        )
    )

    entries = _entries_by_pair(session)
    assert len(entries) == 9
    with_context = {pair for pair, entry in entries.items() if entry.get("searchContext")}
    assert with_context == {
        ("record-1", "record-3"), ("record-2", "record-3"),
        ("record-3", "record-1"), ("record-3", "record-2"),
    }
    manifest = session["proteinIdentityManifest"]
    for entry in entries.values():
        assert validate_protein_raw_entry_references(entry, manifest)
    # searchContext is the Web hash of the searched query and subject FASTA.
    contexts = {
        _sha(json.dumps([_sha(search["query"]), _sha(search["subject"])], separators=(",", ":")))
        for search in searches
    }
    assert {entries[pair]["searchContext"] for pair in with_context} <= contexts


@pytest.mark.linear
def test_pairwise_searches_the_query_file_against_the_subject_file(tmp_path: Path) -> None:
    table = _inputs(tmp_path)
    session, searches = _run(tmp_path, "pairwise", ["--records_table", str(table)], "pairwise")

    # record-1 -> record-2 searches two.gbk without record-1; record-2 -> record-3
    # searches every record of two.gbk against BGC0000711.gbk.
    assert len(searches) == 2
    assert [len(_fasta_ids(search["query"])) > len(_fasta_ids(searches[0]["query"])) for search in searches] == [False, True]
    entries = _entries_by_pair(session)
    assert set(entries) == {("record-1", "record-2"), ("record-2", "record-3")}
    assert "searchContext" not in entries[("record-1", "record-2")]
    assert entries[("record-2", "record-3")]["searchContext"]
    assert [entry["display"] for entry in session["losatCache"]["entries"]] == [True, True]


@pytest.mark.linear
def test_unlisted_pairwise_pairs_search_each_upper_row_record_against_the_row_below(tmp_path: Path) -> None:
    """OV-3870: records #1 and #2 stand in the lower row, record #3 above them. Pairwise without listed pairs searches
    each record of a row against each record of the next row, so record #3 is the query of both searches, as in a
    nucleotide search and in the Web, whatever the record order."""

    _inputs(tmp_path)
    table = tmp_path / "rows.tsv"
    table.write_text(
        "gbk\trow\n"
        f"{tmp_path / 'BGC0000708.gbk'}\t2\n"
        f"{tmp_path / 'BGC0000709.gbk'}\t2\n"
        f"{tmp_path / 'BGC0000711.gbk'}\t1\n",
        encoding="utf-8",
    )
    session, searches = _run(tmp_path, "rows", ["--records_table", str(table)], "pairwise")

    # One query (record #3) against two subjects.
    assert len(searches) == 2
    assert len({frozenset(_fasta_ids(search["query"])) for search in searches}) == 1
    assert len({frozenset(_fasta_ids(search["subject"])) for search in searches}) == 2
    assert [
        (entry["queryRecordInstanceKey"], entry["subjectRecordInstanceKey"])
        for entry in session["losatCache"]["entries"]
    ] == [("record-3", "record-1"), ("record-3", "record-2")]


@pytest.mark.linear
def test_single_record_files_search_each_pair(tmp_path: Path) -> None:
    _inputs(tmp_path)
    files = [str(tmp_path / name) for name in ("BGC0000708.gbk", "BGC0000709.gbk", "BGC0000711.gbk")]
    session, searches = _run(tmp_path, "single", ["--gbk", *files], "similarity_groups")

    assert len(searches) == 9
    entries = _entries_by_pair(session)
    assert len(entries) == 9
    assert not any("searchContext" in entry for entry in entries.values())


def _replay(session_path: Path, prefix: Path) -> dict:
    linear_main(["--session", str(session_path), "-o", str(prefix), "-f", "svg", "--save_session"])
    return drawing_view(read_session_payload(prefix.with_suffix(".gbdraw")))


@pytest.mark.linear
@pytest.mark.parametrize("mode", ["similarity_groups", "pairwise"])
def test_cli_session_keeps_one_file_as_one_source_on_replay(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str
) -> None:
    """PR5-B1: a CLI Session stores the records of one file in one resource."""

    table = _inputs(tmp_path)
    session, _ = _run(tmp_path, mode, ["--records_table", str(table)], mode)
    assert any(entry.get("searchContext") for entry in session["losatCache"]["entries"])

    def search_attempted(*_args: object, **_kwargs: object) -> object:
        raise AssertionError("SEARCH ATTEMPTED")

    monkeypatch.setattr(protein_colinearity, "run_losatp_blastp", search_attempted)
    source = tmp_path / mode / "out"
    replay = tmp_path / mode / "replay"
    replayed = _replay(source.with_suffix(".gbdraw"), replay)
    assert replay.with_suffix(".svg").read_bytes() == source.with_suffix(".svg").read_bytes()

    # The Web shape: one GenBank resource per source file, named after it,
    # record-ID selectors for the records of a multi-record file, a
    # single-record file unchanged.
    records = session["renderRequest"]["records"]
    assert [(record["source"], record["selector"]) for record in records] == [
        ({"kind": "genbank", "resourceId": "record-1-genbank"}, {"kind": "recordId", "value": "BGC0000708"}),
        ({"kind": "genbank", "resourceId": "record-1-genbank"}, {"kind": "recordId", "value": "BGC0000709"}),
        ({"kind": "genbank", "resourceId": "record-3-genbank"}, None),
    ]
    # The Session stores the tables it drew, and Similarity groups its typed
    # result, so the Session draws them without LOSAT (OV-1068, OV-3037).
    assert sorted(session["resources"]) == [
        *(["comparison-orthogroups"] if mode == "similarity_groups" else []),
        "comparison-protein-1",
        "comparison-protein-2",
        "record-1-genbank",
        "record-3-genbank",
    ]
    assert session["resources"]["record-1-genbank"]["name"] == "two.gbk"
    # The search is stored as done: its pipeline marker keeps the settings.
    assert [
        (comparison["kind"], comparison.get("mode"), comparison.get("queryRecordIndex"))
        for comparison in session["renderRequest"]["comparisons"]
        if comparison["kind"] != "orthogroupResult"
    ] == [
        ("precomputedProteinComparison", None, 0),
        ("precomputedProteinComparison", None, 1),
        ("generatedProteinComparison", "none", None),
    ]
    # A replay writes the same layout and keeps every raw entry.
    assert replayed["renderRequest"]["records"] == records
    assert replayed["resources"] == session["resources"]
    assert replayed["losatCache"] == session["losatCache"]


def _table_beside_search(tmp_path: Path, rows: str) -> tuple[list[str], Path]:
    """Three BGC records, a BLAST table of #1 and #2, and a comparisons table of ``rows`` (OV-3688)."""

    for name in ("BGC0000708.gbk", "BGC0000709.gbk", "BGC0000711.gbk"):
        shutil.copy(INPUTS / name, tmp_path / name)
    (tmp_path / "pair.tsv").write_text(
        "BGC0000708\tBGC0000709\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n", encoding="utf-8"
    )
    (tmp_path / "comparisons.tsv").write_text("blast\tquery\tsubject\tsource\n" + rows, encoding="utf-8")
    inputs = ["--gbk", *(str(tmp_path / f"BGC000070{n}.gbk") for n in (8, 9)), str(tmp_path / "BGC0000711.gbk")]
    return [*inputs, "--comparisons_table", str(tmp_path / "comparisons.tsv")], tmp_path / "out"


def _comparison_groups(svg: str) -> dict[tuple[str, str], list[str]]:
    """The match paths of each comparison group of ``svg``, by its query and subject record indexes."""

    import re

    groups: dict[tuple[str, str], list[str]] = {}
    for group in re.finditer(r'<g [^>]*id="comparison\d+"[^>]*>(.*?)</g>', svg, flags=re.S):
        head = group.group(0)[: group.group(0).index(">")]
        query = re.search(r'data-query-record-index="(\d+)"', head)
        subject = re.search(r'data-subject-record-index="(\d+)"', head)
        assert query and subject, head
        groups[(query.group(1), subject.group(1))] = re.findall(r"<path [^>]*>", group.group(1))
    return groups


@pytest.mark.linear
def test_a_blast_table_draws_beside_a_protein_search_of_other_pairs(tmp_path: Path) -> None:
    """OV-3688: a comparisons table draws its ``source=table`` pair and searches its ``source=losat`` pair with
    LOSATP, as the Web draws that Selected plan."""

    losat, log = _fake_losat(tmp_path / "bin-dir")
    argv, prefix = _table_beside_search(tmp_path, "pair.tsv\t#1\t#2\ttable\n\t#2\t#3\tlosat\n")
    linear_main([
        *argv, "--losat", "losatp", "--losatp_mode", "pairwise", "--losat_bin", losat,
        "--losat_threads", "1", "-o", str(prefix), "-f", "svg", "--save_session",
    ])
    # Only the listed pair is searched: BGC0000709's proteins against BGC0000711's.
    [search] = _searches(log)
    described = [
        [line.split()[1] for line in search[side].splitlines() if line.startswith(">")][0] for side in ("query", "subject")
    ]
    assert described == ["CAF33300.1", "CAH58684.1"]
    groups = _comparison_groups(prefix.with_suffix(".svg").read_text(encoding="utf-8"))
    assert sorted(groups) == [("0", "1"), ("1", "2")]
    # The table's one row is a nucleotide ribbon; the search draws protein matches.
    assert len(groups[("0", "1")]) == 1 and "data-query-protein-id" not in groups[("0", "1")][0]
    assert groups[("1", "2")] and all("data-query-protein-id" in path for path in groups[("1", "2")])
    session = drawing_view(read_session_payload(prefix.with_suffix(".gbdraw")))
    assert [
        (comparison["kind"], comparison.get("queryRecordIndex"), comparison.get("mode"))
        for comparison in session["renderRequest"]["comparisons"]
    ] == [
        ("nucleotideBlast", 0, None),
        ("precomputedProteinComparison", 1, None),
        ("generatedProteinComparison", None, "none"),
    ]


@pytest.mark.linear
def test_a_pair_that_a_table_and_a_protein_search_both_name_is_refused(tmp_path: Path) -> None:
    """OV-3688: one record pair drawn from a BLAST table and searched by LOSATP is refused in words."""

    from gbdraw.exceptions import ValidationError

    losat, log = _fake_losat(tmp_path / "bin-dir")
    argv, prefix = _table_beside_search(tmp_path, "pair.tsv\t#1\t#2\ttable\n\t#2\t#1\tlosat\n")
    with pytest.raises(ValidationError, match="record pair #1 and #2 .*both"):
        linear_main([
            *argv, "--losat", "losatp", "--losatp_mode", "pairwise", "--losat_bin", losat,
            "-o", str(prefix), "-f", "svg",
        ])
    assert not log.exists()


@pytest.mark.linear
@pytest.mark.parametrize("given", [(0, 1), (1, 0)])
def test_a_given_comparison_of_a_pair_that_pairwise_searches_unlisted_is_refused(
    tmp_path: Path, given: tuple[int, int]
) -> None:
    """OV-3728: a Pairwise LOSATP search that lists no pairs searches every adjacent pair, so a given comparison of one
    is refused in words (either direction), not merged into the protein matches."""

    import pandas as pd
    from Bio import SeqIO

    from gbdraw.api import (
        InMemoryRecordSource, LinearDiagramOptions, LinearDiagramRequest, LosatSearchOptions, RecordInput,
        render_request,
    )
    from gbdraw.api.options import LosatRuntimeOptions
    from gbdraw.api.requests import RenderOutputRequest
    from gbdraw.exceptions import ValidationError
    from gbdraw.io.comparisons import COMPARISON_COLUMNS
    from gbdraw.linear_comparison import LinearComparison

    losat, log = _fake_losat(tmp_path / "bin-dir")
    names = ("BGC0000708", "BGC0000709", "BGC0000711")
    row = [names[given[0]], names[given[1]], 90.0, 100, 0, 0, 1, 100, 1, 100, 1e-20, 200.0]
    request = LinearDiagramRequest(
        output=RenderOutputRequest(output_prefix="out", output_directory=tmp_path, overwrite=True),
        records=tuple(
            RecordInput(source=InMemoryRecordSource(SeqIO.read(INPUTS / f"{name}.gbk", "genbank")), record_key=f"r{i}")
            for i, name in enumerate(names)
        ),
        options=LinearDiagramOptions(
            linear_comparisons=(LinearComparison(*given, pd.DataFrame([row], columns=COMPARISON_COLUMNS)),),
            losat_search=LosatSearchOptions(
                program="losatp", losatp_mode="pairwise",
                runtime=LosatRuntimeOptions(losat_executable=losat, threads=1),
            ),
        ),
    )
    with pytest.raises(ValidationError, match="record pair #1 and #2 .*both"):
        render_request(request)
    assert not log.exists()
