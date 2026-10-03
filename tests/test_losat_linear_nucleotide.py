"""Linear LOSATN / TLOSATX from the CLI and the Python APIs (design 3.2-3.7, PR-3)."""

from __future__ import annotations

from contextlib import ExitStack
import json
from pathlib import Path
import sys

import pytest

import gbdraw.comparisons.losat_runtime as runtime_module
from gbdraw import LinearComparisonOptions, LinearOptions, draw_linear, read_genbank
from gbdraw.api import LosatRuntimeOptions, LosatSearchOptions
from gbdraw.comparisons.losat_runtime import resolve_losat_runtime
from gbdraw.exceptions import ValidationError
from gbdraw.linear import linear_main

REPO_ROOT = Path(__file__).resolve().parents[1]
LAMBDA = REPO_ROOT / "tests" / "test_inputs" / "NC_001416.gb"
DE3 = REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "de3" / "NC_042057.1.gb"
TUTORIAL_TSV = (
    REPO_ROOT / "gbdraw" / "web" / "tutorial-data" / "lambda-de3-comparison" / "lambda-de3.losatn.tsv"
)
T_CLI_07_SVG = REPO_ROOT / "docs" / "images" / "t-cli-07" / "lambda-de3-losatn.svg"
# T-CLI-07 thresholds and canvas options.
T_CLI_07_OPTIONS = [
    "--record_id", "NC_001416.1", "--record_id", "NC_042057.1",
    "--bitscore", "50", "--evalue", "0.01", "--identity", "0",
    "--alignment_length", "0", "--comparison_height", "120", "-f", "svg",
]

_FAKE = """#!{python}
import json, os, sys
args = sys.argv[1:]
with open({log!r}, "a", encoding="utf-8") as handle:
    handle.write(json.dumps(args) + "\\n")
if args in (["--version"], ["-version"]):
    print("losat 0.1.0")
    sys.exit(0)
if len(args) == 2 and args[1] == "--help":
    print("      -query_gencode <QUERY_GENCODE>")
    print("      -task <TASK>")
    print("          [default: megablast] [possible values: megablast, blastn]")
    sys.exit(0)
if {exit_code}:
    sys.stderr.write("fake failure\\n")
    sys.exit({exit_code})
def ids(path):
    with open(path, encoding="utf-8") as handle:
        return [line[1:].split()[0] for line in handle if line.startswith(">")]
gencode = args[args.index("-query_gencode") + 1] if "-query_gencode" in args else "1"
for query in ids(args[args.index("-query") + 1]):
    for subject in ids(args[args.index("-subject") + 1]):
        print("\\t".join([query, subject, "95.000", "300", "1", "0", "101", "400",
                         "201", "500", "1e-30", str(400 + int(gencode))]))
"""


@pytest.fixture(autouse=True)
def _fresh_probe_caches(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(runtime_module, "_CLI_DIALECTS", {})
    monkeypatch.setattr(runtime_module, "_RUNTIME_VERSIONS", {})
    monkeypatch.setattr(runtime_module, "_TASK_VALUES", {}, raising=False)


def _fake_losat(tmp_path: Path, *, exit_code: int = 0) -> tuple[str, Path]:
    log = tmp_path / "fake.log"
    path = tmp_path / "bin" / "losat"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        _FAKE.format(python=sys.executable, log=str(log), exit_code=exit_code),
        encoding="utf-8",
    )
    path.chmod(0o755)
    return str(path), log


def _searches(log: Path) -> list[list[str]]:
    if not log.exists():
        return []
    calls = [json.loads(line) for line in log.read_text(encoding="utf-8").splitlines()]
    return [call for call in calls if "-query" in call]


def _session(prefix: Path) -> dict:
    return json.loads(prefix.with_suffix(".gbdraw-session.json").read_text(encoding="utf-8"))


def _diagnostic(excinfo: pytest.ExceptionInfo) -> dict:
    return dict(excinfo.value.diagnostic or {})


# Acceptance with the resolved native LOSAT (design 4 PR-3).


@pytest.fixture(scope="module")
def native_losat() -> None:
    with ExitStack() as stack:
        try:
            runtime = resolve_losat_runtime("losatn", stack=stack)
        except ValidationError as error:
            pytest.skip(f"No native LOSAT runtime resolves: {error}")
        if runtime.kind != "losat":
            pytest.skip(f"Resolved runtime is {runtime.kind}, not native LOSAT.")


@pytest.mark.linear
def test_cli_losatn_reproduces_the_t_cli_07_figure(native_losat, tmp_path: Path) -> None:
    prefix = tmp_path / "lambda-de3-losatn"
    raw = tmp_path / "raw"
    linear_main([
        "--gbk", str(LAMBDA), str(DE3), *T_CLI_07_OPTIONS,
        "--losat", "losatn", "--losat_output_dir", str(raw), "--save_session",
        "-o", str(prefix),
    ])

    assert prefix.with_suffix(".svg").read_bytes() == T_CLI_07_SVG.read_bytes()
    edge = raw / "NC_001416.1.NC_042057.1.losatn.tsv"
    assert edge.read_bytes() == TUTORIAL_TSV.read_bytes()
    manifest = (raw / "comparisons.tsv").read_text(encoding="utf-8")
    assert manifest == "blast\tquery\tsubject\nNC_001416.1.NC_042057.1.losatn.tsv\t#1\t#2\n"

    # The manifest is a --comparisons_table without LOSAT.
    again = tmp_path / "from-table"
    linear_main([
        "--gbk", str(LAMBDA), str(DE3), *T_CLI_07_OPTIONS,
        "--comparisons_table", str(raw / "comparisons.tsv"), "-o", str(again),
    ])
    assert again.with_suffix(".svg").read_bytes() == T_CLI_07_SVG.read_bytes()

    session = _session(prefix)
    entry = session["losatCache"]["entries"][0]
    assert entry["key"] == "113e10ccc562a0c207554d16fe73e6df91c6828cdba0d13544fef95358d5f5e3"
    assert entry["runtime"]["program"] == "blastn"
    assert [item["kind"] for item in session["renderRequest"]["comparisons"]][0] == "nucleotideBlast"


# Session and replay.


@pytest.mark.linear
def test_saved_session_replays_without_losat(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    losat, _log = _fake_losat(tmp_path)
    prefix = tmp_path / "saved"
    linear_main([
        "--gbk", str(LAMBDA), str(DE3), *T_CLI_07_OPTIONS,
        "--losat", "losatn", "--losat_bin", losat, "--save_session", "-o", str(prefix),
    ])
    session = _session(prefix)
    comparison = session["renderRequest"]["comparisons"][0]
    assert comparison["kind"] == "nucleotideBlast"
    assert (comparison["queryRecordIndex"], comparison["subjectRecordIndex"]) == (0, 1)
    entry = session["losatCache"]["entries"][0]
    assert entry["schema"] == 2 and entry["identityKind"] == "nucleotide"
    assert entry["args"] == ["--task", "megablast"]
    assert entry["runtime"]["source"] == "explicit"
    assert entry["display"] is True and entry["filename"].endswith(".losatn.tsv")

    def no_runtime(*_args, **_kwargs):
        raise AssertionError("Session replay resolved a LOSAT runtime.")

    monkeypatch.setattr(runtime_module, "resolve_losat_runtime", no_runtime)
    replay = tmp_path / "replay"
    linear_main([
        "--session", str(prefix.with_suffix(".gbdraw-session.json")),
        "-o", str(replay), "-f", "svg",
    ])
    assert replay.with_suffix(".svg").read_bytes() == prefix.with_suffix(".svg").read_bytes()


# TLOSATX translation tables.


@pytest.mark.linear
def test_tlosatx_record_gencodes_reach_the_argv_and_the_raw_key(tmp_path: Path) -> None:
    losat, log = _fake_losat(tmp_path)
    prefix = tmp_path / "coded"
    linear_main([
        "--gbk", str(LAMBDA), str(DE3), *T_CLI_07_OPTIONS,
        "--losat", "tlosatx", "--losat_gencode", "4", "11", "--losat_bin", losat,
        "--save_session", "-o", str(prefix),
    ])
    (search,) = _searches(log)
    assert search[0] == "tblastx"
    assert search[search.index("-query_gencode") + 1] == "4"
    assert search[search.index("-db_gencode") + 1] == "11"
    entry = _session(prefix)["losatCache"]["entries"][0]
    assert entry["args"] == ["--query-gencode", "4", "--db-gencode", "11"]
    # The fake runtime scores by query gencode: the table changed the result.
    assert entry["text"].rstrip("\n").endswith("\t404")

    default = tmp_path / "default"
    linear_main([
        "--gbk", str(LAMBDA), str(DE3), *T_CLI_07_OPTIONS,
        "--losat", "tlosatx", "--losat_bin", losat, "--save_session", "-o", str(default),
    ])
    assert "-query_gencode" not in _searches(log)[-1]
    default_entry = _session(default)["losatCache"]["entries"][0]
    assert default_entry["args"] == []
    assert default_entry["key"] != entry["key"]


def test_records_table_losat_gencode_column(tmp_path: Path) -> None:
    losat, log = _fake_losat(tmp_path)
    table = tmp_path / "records.tsv"
    table.write_text(
        f"gbk\trecord_id\tlosat_gencode\n{LAMBDA}\tNC_001416.1\t5\n{DE3}\tNC_042057.1\t\n",
        encoding="utf-8",
    )
    linear_main([
        "--records_table", str(table), "--losat", "tlosatx", "--losat_bin", losat,
        "-o", str(tmp_path / "table"), "-f", "svg",
    ])
    (search,) = _searches(log)
    assert search[search.index("-query_gencode") + 1] == "5"
    assert "-db_gencode" not in search
    with pytest.raises(ValidationError, match="not both"):
        linear_main([
            "--records_table", str(table), "--losat", "tlosatx", "--losat_gencode", "11",
            "--losat_bin", losat, "-o", str(tmp_path / "both"), "-f", "svg",
        ])


# Mixed LOSAT and table edges.


@pytest.mark.linear
def test_comparisons_table_mixes_losat_and_table_edges(tmp_path: Path) -> None:
    losat, log = _fake_losat(tmp_path)
    table = tmp_path / "pairs.tsv"
    table.write_text(
        "source\tblast\tquery\tsubject\n"
        "losat\t\t#1\t#2\n"
        f"table\t{TUTORIAL_TSV}\t#2\t#3\n",
        encoding="utf-8",
    )
    prefix = tmp_path / "mixed"
    linear_main([
        "--gbk", str(DE3), str(LAMBDA), str(DE3),
        "--comparisons_table", str(table), "--losat", "losatn", "--losat_bin", losat,
        "--save_session", "-o", str(prefix), "-f", "svg",
    ])
    assert len(_searches(log)) == 1
    comparisons = _session(prefix)["renderRequest"]["comparisons"]
    endpoints = {
        (item["kind"], item["queryRecordIndex"], item["subjectRecordIndex"])
        for item in comparisons
        if "queryRecordIndex" in item
    }
    assert ("nucleotideBlast", 0, 1) in endpoints
    assert ("precomputedProteinComparison", 1, 2) in endpoints


# Diagnostics (design 3.7).


@pytest.mark.parametrize(
    ("kwargs", "field"),
    [
        ({"program": "losatn", "losatp_mode": "pairwise"}, "losatp_mode"),
        ({"program": "losatn", "record_gencodes": (11,)}, "record_gencodes"),
        ({"program": "tlosatx", "losatn_task": "blastn"}, "losatn_task"),
        ({"program": "losatp", "losatp_mode": "pairwise", "losatn_task": "blastn"}, "losatn_task"),
        ({"program": "tlosatx", "losatp_max_hits": 3}, "losatp_max_hits"),
    ],
)
def test_options_of_another_program_are_rejected(kwargs: dict, field: str) -> None:
    with pytest.raises(ValidationError) as excinfo:
        LosatSearchOptions(**kwargs)
    assert _diagnostic(excinfo) == {
        "code": "COMPARISON_INPUT",
        "reason": "LOSAT_OPTION_PROGRAM",
        "field": field,
        "program": kwargs["program"],
    }


@pytest.mark.parametrize(
    ("argv", "message", "field"),
    [
        (["--losat", "losatn", "--losat_gencode", "11"], "--losat_gencode", "record_gencodes"),
        (["--losat", "tlosatx", "--losatn_task", "blastn"], "--losatn_task", "losatn_task"),
        (["--losat", "losatn", "--losatp_mode", "pairwise"], "--losatp_mode", "losatp_mode"),
        (["--losat", "tlosatx", "--collinear_min_anchors", "2"], "--collinear_min_anchors", "collinear_min_anchors"),
    ],
)
def test_cli_option_program_diagnostics(tmp_path: Path, argv: list[str], message: str, field: str) -> None:
    with pytest.raises(ValidationError, match=message) as excinfo:
        linear_main(["--gbk", str(LAMBDA), str(DE3), *argv, "-o", str(tmp_path / "x")])
    assert _diagnostic(excinfo)["reason"] == "LOSAT_OPTION_PROGRAM"
    assert _diagnostic(excinfo)["field"] == field


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("case", "match"),
    [
        ("losat-row-without-losat", "source=losat"),
        ("losat-without-losat-row", "no comparisons table row"),
        ("single-record", "at least two records"),
        ("losatp-selected-not-pairwise", "pairwise"),
    ],
)
def test_losat_plan_diagnostics(tmp_path: Path, case: str, match: str) -> None:
    losat, log = _fake_losat(tmp_path)
    losat_rows = _write(tmp_path, "losat.tsv", "source\tblast\tquery\tsubject\nlosat\t\t#1\t#2\n")
    table_rows = _write(tmp_path, "table.tsv", f"blast\tquery\tsubject\n{TUTORIAL_TSV}\t#1\t#2\n")
    argv = {
        "losat-row-without-losat": ["--gbk", str(LAMBDA), str(DE3), "--comparisons_table", str(losat_rows)],
        "losat-without-losat-row": [
            "--gbk", str(LAMBDA), str(DE3), "--comparisons_table", str(table_rows), "--losat", "losatn",
        ],
        "single-record": ["--gbk", str(LAMBDA), "--losat", "losatn"],
        "losatp-selected-not-pairwise": [
            "--gbk", str(LAMBDA), str(DE3), "--comparisons_table", str(losat_rows),
            "--losat", "losatp", "--losatp_mode", "similarity_groups",
        ],
    }[case]
    with pytest.raises(ValidationError, match=match) as excinfo:
        linear_main([*argv, "--losat_bin", losat, "-o", str(tmp_path / "x"), "-f", "svg"])
    assert _diagnostic(excinfo)["code"] == "COMPARISON_INPUT"
    assert _diagnostic(excinfo)["reason"] == "LOSAT_PLAN"
    assert _searches(log) == []


def test_losat_with_blast_files_is_a_plan_error() -> None:
    from gbdraw.api import LinearDiagramOptions
    from gbdraw.comparisons.linear_losat import resolve_linear_nucleotide_losat

    records = read_genbank([LAMBDA, DE3])
    with pytest.raises(ValidationError) as excinfo:
        resolve_linear_nucleotide_losat(
            LinearDiagramOptions(
                blast_files=(str(TUTORIAL_TSV),),
                losat_search=LosatSearchOptions(program="losatn"),
            ),
            records=records,
            rows_by_record=(0, 1),
            source_ids=("a", "b"),
            record_keys=("record-1", "record-2"),
            record_labels=("a", "b"),
            input_indexes=(0, 1),
        )
    assert _diagnostic(excinfo)["reason"] == "LOSAT_PLAN"


def test_unsupported_losatn_task_fails_before_searching(tmp_path: Path) -> None:
    losat, log = _fake_losat(tmp_path)
    with pytest.raises(ValidationError, match="dc-megablast") as excinfo:
        linear_main([
            "--gbk", str(LAMBDA), str(DE3), "--losat", "losatn", "--losatn_task", "dc-megablast",
            "--losat_bin", losat, "-o", str(tmp_path / "x"), "-f", "svg",
        ])
    assert _diagnostic(excinfo)["code"] == "COMPARISON_INPUT"
    assert _diagnostic(excinfo)["reason"] == "LOSAT_TASK"
    assert "0.1.0" in str(excinfo.value)
    assert _searches(log) == []


def test_runtime_failure_has_a_losat_runtime_diagnostic(tmp_path: Path) -> None:
    losat, _log = _fake_losat(tmp_path, exit_code=3)
    with pytest.raises(ValidationError, match="exit code 3") as excinfo:
        linear_main([
            "--gbk", str(LAMBDA), str(DE3), "--losat", "losatn", "--losat_bin", losat,
            "-o", str(tmp_path / "x"), "-f", "svg",
        ])
    assert _diagnostic(excinfo) == {
        "code": "LOSAT_RUNTIME", "reason": "FAILED", "program": "losatn", "exitCode": 3,
    }


def test_missing_runtime_is_unavailable(tmp_path: Path) -> None:
    with pytest.raises(ValidationError) as excinfo:
        linear_main([
            "--gbk", str(LAMBDA), str(DE3), "--losat", "losatn",
            "--losat_bin", str(tmp_path / "missing" / "losat"),
            "-o", str(tmp_path / "x"), "-f", "svg",
        ])
    assert _diagnostic(excinfo) == {"code": "LOSAT_RUNTIME", "reason": "UNAVAILABLE"}


# Typed and introductory APIs.


@pytest.mark.linear
def test_typed_and_introductory_apis_run_tlosatx(tmp_path: Path) -> None:
    losat, log = _fake_losat(tmp_path)
    records = read_genbank([LAMBDA, DE3])
    draw_linear(
        records,
        options=LinearOptions(
            comparisons=LinearComparisonOptions(
                losat="tlosatx", gencodes=[4, None], pairs=[(0, 1)], losat_executable=losat,
            )
        ),
    )
    (search,) = _searches(log)
    assert search[search.index("-query_gencode") + 1] == "4"
    assert "-db_gencode" not in search

    search_options = LosatSearchOptions(
        program="tlosatx",
        record_gencodes=(11,),
        pairs=((0, 1),),
        runtime=LosatRuntimeOptions(losat_executable=losat),
    )
    assert search_options.record_gencodes == (11,)
    assert LosatSearchOptions(program="losatn").losatn_task == "megablast"


def test_unresolved_nucleotide_intent_is_never_encoded() -> None:
    from gbdraw.api import LinearDiagramOptions, LinearDiagramRequest, RecordInput
    from gbdraw.api.requests import GenBankInputSource
    from gbdraw.session_request_codec import (
        CanonicalRequestEncodingError,
        encode_canonical_request,
    )

    request = LinearDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(str(LAMBDA))),),
        options=LinearDiagramOptions(losat_search=LosatSearchOptions(program="losatn")),
    )
    with pytest.raises(CanonicalRequestEncodingError, match="losatn search"):
        encode_canonical_request(request)


def test_edge_filenames_follow_the_shared_vectors() -> None:
    from gbdraw.comparisons.linear_losat import losat_edge_filename

    vectors = Path(__file__).parent / "fixtures" / "losat_edge_filename_cases.json"
    for case in json.loads(vectors.read_text(encoding="utf-8"))["cases"]:
        assert losat_edge_filename(case["left"], case["right"], case["suffix"]) == case["expected"]
