"""CLI and Python LOSATP search the Web source-file databases (design D7, PR-5).

A source file is one genome: records of one file are searched together, a
record never searches itself unless that search is requested, and raw keys
carry the Web ``searchContext`` when a side holds several records.
"""

from __future__ import annotations

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
    monkeypatch.setattr(runtime_module, "_CLI_DIALECTS", {})
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
    session = json.loads(prefix.with_suffix(".gbdraw-session.json").read_text(encoding="utf-8"))
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
    return json.loads(prefix.with_suffix(".gbdraw-session.json").read_text(encoding="utf-8"))


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
    replayed = _replay(source.with_suffix(".gbdraw-session.json"), replay)
    assert replay.with_suffix(".svg").read_bytes() == source.with_suffix(".svg").read_bytes()

    # The Web shape: one GenBank resource per source file, `#k` selectors for
    # the records of a multi-record file, a single-record file unchanged.
    records = session["renderRequest"]["records"]
    assert [(record["source"], record["selector"]) for record in records] == [
        ({"kind": "genbank", "resourceId": "record-1-genbank"}, {"kind": "recordIndex", "index": 0}),
        ({"kind": "genbank", "resourceId": "record-1-genbank"}, {"kind": "recordIndex", "index": 1}),
        ({"kind": "genbank", "resourceId": "record-3-genbank"}, None),
    ]
    assert sorted(session["resources"]) == ["record-1-genbank", "record-3-genbank"]
    # A replay writes the same layout and keeps every raw entry.
    assert replayed["renderRequest"]["records"] == records
    assert replayed["resources"] == session["resources"]
    assert replayed["losatCache"] == session["losatCache"]
