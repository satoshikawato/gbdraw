"""One LOSAT runtime owner for LOSATN, TLOSATX, and LOSATP.

Resolution order, CLI dialect, argv, raw-cache args, and the runtime record
are program data in ``gbdraw.comparisons.losat_runtime``. The integration
tests at the end run whichever native LOSAT resolves (the bundled build in a
source checkout) and reproduce the tracked tutorial evidence. Set
``GBDRAW_TEST_LOSAT_BIN`` to run them against another LOSAT executable.
"""

from __future__ import annotations

from collections import Counter
from contextlib import ExitStack
import json
import os
from pathlib import Path
import sys

import pandas as pd
import pytest
from Bio import SeqIO

import gbdraw.analysis.protein_colinearity as protein_colinearity_module
import gbdraw.comparisons.losat_runtime as runtime_module
import gbdraw.losat_setup as losat_setup_module
from gbdraw.analysis.protein_colinearity import (
    PROTEIN_LOSAT_CACHE_SCHEMA,
    LosatpCacheManager,
    build_pairwise_protein_blastp_comparisons,
    build_protein_losat_cache_key,
    build_protein_losat_pair_identity,
    extract_web_stable_cds_proteins,
    parse_losatp_outfmt6,
    validate_protein_raw_entry_references,
)
from gbdraw.comparisons.losat_runtime import (
    LOSAT_PROGRAMS,
    LosatRuntime,
    LosatSearchArgs,
    build_losat_command,
    detect_losat_cli_dialect,
    losat_cache_args,
    losat_runtime_record,
    resolve_losat_runtime,
    run_losat_search,
)
from gbdraw.exceptions import ValidationError
from tests.test_protein_colinearity import _cds, _record

REPO_ROOT = Path(__file__).resolve().parents[1]
TUTORIAL_DATA = REPO_ROOT / "gbdraw" / "web" / "tutorial-data"
PROGRAMS = ("losatn", "tlosatx", "losatp")
NCBI_EXECUTABLES = {"losatn": "blastn", "tlosatx": "tblastx", "losatp": "blastp"}
OPTIONS = {
    "losatn": LosatSearchArgs(task="megablast"),
    "tlosatx": LosatSearchArgs(query_gencode=5, db_gencode=2),
    "losatp": LosatSearchArgs(max_hsps=1, max_target_seqs=5),
}

_FAKE_RUNTIME = """#!{python}
import json, os, sys
args = sys.argv[1:]
log = os.environ.get("GBDRAW_FAKE_RUNTIME_LOG")
if log:
    with open(log, "a", encoding="utf-8") as handle:
        handle.write(json.dumps([os.path.basename(sys.argv[0]), *args]) + "\\n")
DIALECT = {dialect!r}
if args in (["--version"], ["-version"]):
    print({version!r})
    sys.exit(0)
if len(args) == 2 and args[1] == "--help":
    print("      --query-gencode <QUERY_GENCODE>" if DIALECT == "v1" else "      -query_gencode <QUERY_GENCODE>")
    sys.exit(0)
rejected = "-query_gencode" if DIALECT == "v1" else "--query-gencode"
if rejected in args:
    sys.stderr.write("error: unexpected argument " + rejected + "\\n")
    sys.exit(2)
def first_id(path):
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                return line[1:].split()[0]
query = first_id(args[args.index("-query") + 1])
subject = first_id(args[args.index("-subject") + 1])
print(query + "\\t" + subject + "\\t100.000\\t10\\t0\\t0\\t1\\t10\\t1\\t10\\t1e-10\\t50.0")
"""


def _fake_runtime(
    path: Path,
    *,
    dialect: str = "v2",
    version: str = "losat 0.1.0",
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        _FAKE_RUNTIME.format(python=sys.executable, dialect=dialect, version=version),
        encoding="utf-8",
    )
    path.chmod(0o755)
    return path.absolute()


def _logged_calls(log: Path) -> list[list[str]]:
    if not log.exists():
        return []
    return [json.loads(line) for line in log.read_text(encoding="utf-8").splitlines()]


@pytest.fixture(autouse=True)
def _fresh_probe_caches(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(runtime_module, "_CLI_DIALECTS", {})
    monkeypatch.setattr(runtime_module, "_RUNTIME_VERSIONS", {})


# Resolution order (design 3.6): explicit, conda, managed, bundled, PATH losat,
# PATH NCBI <program>, else unavailable. Identical for every program.
RESOLUTION_ORDER = ("conda", "managed", "bundled", "path-losat", "path-ncbi")


@pytest.fixture
def runtime_sources(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Install fake candidates for each automatic source, enabled on request."""

    candidates = {
        "conda": tmp_path / "conda" / "bin" / "losat",
        "managed": tmp_path / "managed" / "LOSAT",
        "bundled": tmp_path / "bundled" / "losat",
        "path-losat": tmp_path / "path" / "losat",
    }
    path_dir = tmp_path / "path"
    path_dir.mkdir()
    venv_prefix = tmp_path / "venv"
    venv_prefix.mkdir()
    monkeypatch.setenv("PATH", str(path_dir))
    monkeypatch.setattr(runtime_module.sys, "prefix", str(venv_prefix))
    monkeypatch.setattr(losat_setup_module, "managed_losat", lambda: None)
    monkeypatch.setattr(runtime_module, "_bundled_losat_resource", lambda: None)

    def enable(source: str, program: str) -> str:
        if source == "conda":
            (tmp_path / "conda" / "conda-meta").mkdir(parents=True, exist_ok=True)
            monkeypatch.setattr(runtime_module.sys, "prefix", str(tmp_path / "conda"))
            return str(_fake_runtime(candidates["conda"]))
        if source == "managed":
            binary = _fake_runtime(candidates["managed"])
            monkeypatch.setattr(losat_setup_module, "managed_losat", lambda: binary)
            return str(binary)
        if source == "bundled":
            binary = _fake_runtime(candidates["bundled"])
            monkeypatch.setattr(runtime_module, "_bundled_losat_resource", lambda: binary)
            return str(binary)
        if source == "path-losat":
            return str(_fake_runtime(candidates["path-losat"]))
        binary = _fake_runtime(path_dir / NCBI_EXECUTABLES[program], version="2.16.0+")
        return str(binary)

    return enable


@pytest.mark.parametrize("program", PROGRAMS)
@pytest.mark.parametrize("winner_index", range(len(RESOLUTION_ORDER)))
def test_every_program_resolves_runtimes_in_one_order(
    runtime_sources,
    program: str,
    winner_index: int,
) -> None:
    executables = {
        source: runtime_sources(source, program)
        for source in RESOLUTION_ORDER[winner_index:]
    }
    winner = RESOLUTION_ORDER[winner_index]

    with ExitStack() as stack:
        runtime = resolve_losat_runtime(program, stack=stack)

    expected_kind = "ncbi-blast" if winner == "path-ncbi" else "losat"
    expected_source = {
        "conda": "conda",
        "managed": "managed",
        "bundled": "bundled",
        "path-losat": "path",
        "path-ncbi": "path",
    }[winner]
    assert (runtime.kind, runtime.source, runtime.executable) == (
        expected_kind,
        expected_source,
        executables[winner],
    )


@pytest.mark.parametrize("program", PROGRAMS)
def test_explicit_runtimes_bypass_discovery_for_every_program(
    monkeypatch: pytest.MonkeyPatch,
    program: str,
) -> None:
    def fail(*_args, **_kwargs):
        raise AssertionError("automatic discovery ran")

    monkeypatch.setattr(runtime_module, "_conda_losat_runtime", fail)
    monkeypatch.setattr(runtime_module, "_bundled_losat_resource", fail)
    monkeypatch.setattr(runtime_module, "_path_executable", fail)
    monkeypatch.setattr(losat_setup_module, "managed_losat", fail)

    with ExitStack() as stack:
        explicit_losat = resolve_losat_runtime(program, losat_bin="/opt/losat", stack=stack)
        explicit_ncbi = resolve_losat_runtime(
            program,
            ncbi_blast_bin=f"/opt/ncbi/{NCBI_EXECUTABLES[program]}",
            stack=stack,
        )
        with pytest.raises(ValidationError, match="not both"):
            resolve_losat_runtime(
                program,
                losat_bin="/opt/losat",
                ncbi_blast_bin="/opt/ncbi/blast",
                stack=stack,
            )

    assert explicit_losat == LosatRuntime("losat", "/opt/losat", "explicit")
    assert explicit_ncbi == LosatRuntime(
        "ncbi-blast", f"/opt/ncbi/{NCBI_EXECUTABLES[program]}", "explicit"
    )


@pytest.mark.parametrize("program", PROGRAMS)
def test_unavailable_runtime_names_the_program_specific_ncbi_fallback(
    runtime_sources,
    tmp_path: Path,
    program: str,
) -> None:
    # Another program's NCBI executable on PATH is not a fallback.
    other = next(name for key, name in NCBI_EXECUTABLES.items() if key != program)
    _fake_runtime(tmp_path / "path" / other, version="2.16.0+")

    with ExitStack() as stack, pytest.raises(ValidationError) as exc_info:
        resolve_losat_runtime(program, stack=stack)

    message = str(exc_info.value)
    assert "needs LOSAT or NCBI BLAST+" in message
    assert "`losat` was not found on PATH" in message
    assert f"`{NCBI_EXECUTABLES[program]}` was not found on PATH" in message
    assert "--losat_bin" in message and "--ncbi_blast_bin" in message


# Argv goldens. LOSATP argv is byte-identical to the pre-owner builders for
# both LOSAT dialects; only TLOSATX gencode flags differ between CLI v1 and v2.
_Q = Path("q.fasta")
_S = Path("s.fasta")
_IO = ["-query", "q.fasta", "-subject", "s.fasta", "-outfmt", "6"]
ARGV_GOLDENS = {
    ("losatn", "v2"): ["/opt/losat", "blastn", *_IO, "-task", "megablast", "-num_threads", "4"],
    ("losatn", "v1"): ["/opt/losat", "blastn", *_IO, "-task", "megablast", "-num_threads", "4"],
    ("losatn", "ncbi"): ["/opt/ncbi/blastn", *_IO, "-task", "megablast", "-num_threads", "4"],
    ("tlosatx", "v2"): [
        "/opt/losat", "tblastx", *_IO,
        "-query_gencode", "5", "-db_gencode", "2", "-num_threads", "4",
    ],
    ("tlosatx", "v1"): [
        "/opt/losat", "tblastx", *_IO,
        "--query-gencode", "5", "--db-gencode", "2", "-num_threads", "4",
    ],
    ("tlosatx", "ncbi"): [
        "/opt/ncbi/tblastx", *_IO,
        "-query_gencode", "5", "-db_gencode", "2", "-num_threads", "4",
    ],
    ("losatp", "v2"): [
        "/opt/losat", "blastp", *_IO,
        "-max_hsps", "1", "-max_target_seqs", "5", "-num_threads", "4",
    ],
    ("losatp", "v1"): [
        "/opt/losat", "blastp", *_IO,
        "-max_hsps", "1", "-max_target_seqs", "5", "-num_threads", "4",
    ],
    ("losatp", "ncbi"): [
        "/opt/ncbi/blastp", *_IO,
        "-max_hsps", "1", "-max_target_seqs", "5", "-num_threads", "4",
    ],
}


@pytest.mark.parametrize(("program", "dialect"), sorted(ARGV_GOLDENS))
def test_argv_goldens_for_losat_dialects_and_ncbi(program: str, dialect: str) -> None:
    if dialect == "ncbi":
        runtime = LosatRuntime("ncbi-blast", f"/opt/ncbi/{NCBI_EXECUTABLES[program]}", "path")
        command = build_losat_command(
            runtime, program, query_path=_Q, subject_path=_S,
            options=OPTIONS[program], threads=4,
        )
    else:
        runtime = LosatRuntime("losat", "/opt/losat", "explicit")
        command = build_losat_command(
            runtime, program, query_path=_Q, subject_path=_S,
            options=OPTIONS[program], threads=4, dialect=dialect,
        )

    assert command == ARGV_GOLDENS[(program, dialect)]


def test_losatp_argv_omits_unset_options_like_the_previous_builder() -> None:
    losat = LosatRuntime("losat", "losat", "explicit")
    ncbi = LosatRuntime("ncbi-blast", "blastp", "explicit")
    options = LosatSearchArgs(max_hsps=None, max_target_seqs=None)

    assert build_losat_command(
        losat, "losatp", query_path=_Q, subject_path=_S, options=options, dialect="v1"
    ) == ["losat", "blastp", *_IO]
    assert build_losat_command(
        ncbi, "losatp", query_path=_Q, subject_path=_S, options=options
    ) == ["blastp", *_IO]


def test_cli_dialect_is_probed_once_per_executable_and_only_when_flags_differ(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    log = tmp_path / "calls.jsonl"
    monkeypatch.setenv("GBDRAW_FAKE_RUNTIME_LOG", str(log))
    v1 = _fake_runtime(tmp_path / "v1" / "losat", dialect="v1")
    v2 = _fake_runtime(tmp_path / "v2" / "losat", dialect="v2")

    for program in ("losatn", "losatp"):
        build_losat_command(
            LosatRuntime("losat", str(v1), "explicit"), program,
            query_path=_Q, subject_path=_S, options=OPTIONS[program],
        )
    assert _logged_calls(log) == []

    v1_command = build_losat_command(
        LosatRuntime("losat", str(v1), "explicit"), "tlosatx",
        query_path=_Q, subject_path=_S, options=OPTIONS["tlosatx"],
    )
    v2_command = build_losat_command(
        LosatRuntime("losat", str(v2), "explicit"), "tlosatx",
        query_path=_Q, subject_path=_S, options=OPTIONS["tlosatx"],
    )
    assert detect_losat_cli_dialect(str(v1)) == "v1"
    assert detect_losat_cli_dialect(str(v2)) == "v2"

    assert "--query-gencode" in v1_command and "-query_gencode" not in v1_command
    assert "-query_gencode" in v2_command and "--query-gencode" not in v2_command
    assert _logged_calls(log) == [["losat", "tblastx", "--help"]] * 2


CACHE_ARG_GOLDENS = {
    ("losatn", LosatSearchArgs(task="megablast")): ["--task", "megablast"],
    ("losatn", LosatSearchArgs(task="dc-megablast")): ["--task", "dc-megablast"],
    ("tlosatx", LosatSearchArgs(query_gencode=5, db_gencode=2)): [
        "--query-gencode", "5", "--db-gencode", "2",
    ],
    # The Web passes a translation table only when the record sets one.
    ("tlosatx", LosatSearchArgs(query_gencode=11)): ["--query-gencode", "11"],
    ("tlosatx", LosatSearchArgs()): [],
    ("losatp", LosatSearchArgs(max_hsps=1, max_target_seqs=5)): [
        "--max-hsps-per-subject", "1", "--max-target-seqs", "5",
    ],
    ("losatp", LosatSearchArgs(max_target_seqs=5)): ["--max-target-seqs", "5"],
    ("losatp", LosatSearchArgs()): [],
}


@pytest.mark.parametrize(("program", "options"), list(CACHE_ARG_GOLDENS))
def test_cache_args_keep_the_web_v1_form(program: str, options: LosatSearchArgs) -> None:
    assert losat_cache_args(program, options) == CACHE_ARG_GOLDENS[(program, options)]


@pytest.mark.parametrize(
    ("program", "options", "message"),
    [
        ("losatn", LosatSearchArgs(), "requires"),
        ("losatp", LosatSearchArgs(task="megablast"), "does not apply"),
        ("losatn", LosatSearchArgs(task="megablast", query_gencode=1), "does not apply"),
    ],
)
def test_program_options_are_validated_as_data(
    program: str,
    options: LosatSearchArgs,
    message: str,
) -> None:
    with pytest.raises(ValidationError, match=message):
        losat_cache_args(program, options)


def test_program_table_names_search_programs() -> None:
    assert {name: spec.search for name, spec in LOSAT_PROGRAMS.items()} == NCBI_EXECUTABLES


def test_moved_runtime_paths_leave_no_shims_in_protein_colinearity() -> None:
    for name in (
        "ProteinBlastpRuntime",
        "_resolve_protein_blastp_runtime",
        "_bundled_losatp_resource",
        "_conda_losatp_runtime",
        "_build_losat_blastp_command",
        "_build_ncbi_blastp_command",
        "_build_protein_blastp_command",
        "_run_protein_blastp_subprocess",
        "_losatp_cache_args",
    ):
        assert not hasattr(protein_colinearity_module, name), name


@pytest.mark.parametrize("dialect", ["v1", "v2"])
def test_tlosatx_search_runs_in_the_detected_dialect_and_records_the_runtime(
    tmp_path: Path,
    dialect: str,
) -> None:
    binary = _fake_runtime(tmp_path / dialect / "losat", dialect=dialect)
    records: list[dict[str, object]] = []

    text = run_losat_search(
        "tlosatx",
        ">query_a\nATGAAA\n",
        ">subject_b\nATGAAA\n",
        options=LosatSearchArgs(query_gencode=11, db_gencode=4),
        losat_bin=str(binary),
        threads=2,
        runtime_callback=records.append,
    )

    assert text.startswith("query_a\tsubject_b\t")
    assert records == [
        {
            "kind": "losat",
            "version": "0.1.0",
            "source": "explicit",
            "path": str(binary),
            "program": "tblastx",
            "cli": dialect,
        }
    ]


def test_ncbi_runtime_record_has_version_and_no_losat_dialect(tmp_path: Path) -> None:
    binary = _fake_runtime(
        tmp_path / "ncbi" / "blastn",
        version="blastn: 2.16.0+\n Package: blast 2.16.0, build Jun 25 2024",
    )
    records: list[dict[str, object]] = []

    run_losat_search(
        "losatn",
        ">q\nACGT\n",
        ">s\nACGT\n",
        options=LosatSearchArgs(task="megablast"),
        ncbi_blast_bin=str(binary),
        runtime_callback=records.append,
    )

    assert records == [
        {
            "kind": "ncbi-blast",
            "version": "2.16.0+",
            "source": "explicit",
            "path": str(binary),
            "program": "blastn",
        }
    ]


def test_bundled_runtime_record_uses_the_package_relative_path(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    binary = _fake_runtime(tmp_path / "checkout" / "losat", dialect="v1")
    monkeypatch.setattr(runtime_module, "_bundled_platform_dir", lambda: "linux-x86_64")

    record = losat_runtime_record(LosatRuntime("losat", str(binary), "bundled"), "losatp")

    assert record == {
        "kind": "losat",
        "version": "0.1.0",
        "source": "bundled",
        "path": "gbdraw/bin/linux-x86_64/losat",
        "program": "blastp",
        "cli": "v1",
    }


@pytest.mark.linear
def test_losatp_cache_records_runtime_outside_the_key_and_round_trips(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
    ]
    extraction = extract_web_stable_cds_proteins(
        records,
        record_instance_keys=("r_left", "r_right"),
    )
    runtime_record = {
        "kind": "losat",
        "version": "0.1.0",
        "source": "bundled",
        "path": "gbdraw/bin/linux-x86_64/losat",
        "program": "blastp",
        "cli": "v1",
    }

    def fake_losatp(query_fasta: str, subject_fasta: str, **kwargs) -> pd.DataFrame:
        query_id = query_fasta.splitlines()[0][1:].split()[0]
        subject_id = subject_fasta.splitlines()[0][1:].split()[0]
        raw_text = f"{query_id}\t{subject_id}\t90\t3\t0\t0\t1\t3\t1\t3\t1e-20\t200\n"
        kwargs["raw_output_callback"](raw_text)
        kwargs["runtime_callback"](dict(runtime_record))
        return parse_losatp_outfmt6(raw_text)

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", fake_losatp)
    cache = LosatpCacheManager(identity_manifest=extraction.identity_manifest)
    build_pairwise_protein_blastp_comparisons(
        records,
        losatp_cache=cache,
        protein_extraction=extraction,
        cache_filenames=("record_a.record_b.losatp.tsv",),
    )

    (entry,) = cache.session_entries()
    pair_identity = build_protein_losat_pair_identity(
        extraction.identity_manifest,
        query_record_instance_key="r_left",
        subject_record_instance_key="r_right",
    )
    assert entry["schema"] == PROTEIN_LOSAT_CACHE_SCHEMA
    assert entry["runtime"] == runtime_record
    assert entry["key"] == build_protein_losat_cache_key(
        pair_identity,
        args=losat_cache_args("losatp", LosatSearchArgs(max_hsps=1)),
    )
    assert validate_protein_raw_entry_references(entry, extraction.identity_manifest)

    def fail_run(*_args, **_kwargs):
        raise AssertionError("LOSATP should not run on a cache hit")

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", fail_run)
    replay = LosatpCacheManager([entry], identity_manifest=extraction.identity_manifest)
    build_pairwise_protein_blastp_comparisons(
        records,
        losatp_cache=replay,
        protein_extraction=extraction,
        cache_filenames=("record_a.record_b.losatp.tsv",),
    )
    assert replay.session_entries() == (entry,)


# Integration: the resolved native LOSAT reproduces tracked evidence.


def _native_losat_bin() -> str | None:
    return os.environ.get("GBDRAW_TEST_LOSAT_BIN") or None


@pytest.fixture(scope="module")
def native_losat() -> dict[str, object]:
    with ExitStack() as stack:
        try:
            runtime = resolve_losat_runtime(
                "tlosatx",
                losat_bin=_native_losat_bin(),
                stack=stack,
            )
        except ValidationError as error:
            pytest.skip(f"No native LOSAT runtime resolves: {error}")
        if runtime.kind != "losat":
            pytest.skip(f"Resolved runtime is {runtime.kind}, not native LOSAT.")
        return losat_runtime_record(runtime, "tlosatx")


def _fasta(path: Path, *, width: int | None = None) -> str:
    record = SeqIO.read(path, "genbank")
    sequence = str(record.seq).upper()
    if width is None:
        return f">{record.id}\n{sequence}\n"
    lines = [f">{record.id}"]
    lines.extend(sequence[index : index + width] for index in range(0, len(sequence), width))
    return "\n".join(lines) + "\n"


def _lambda_de3() -> tuple[str, str]:
    return (
        _fasta(REPO_ROOT / "tests" / "test_inputs" / "NC_001416.gb"),
        _fasta(TUTORIAL_DATA / "de3" / "NC_042057.1.gb"),
    )


def _native_search(program: str, query: str, subject: str, options, *, threads: int = 1) -> str:
    return run_losat_search(
        program,
        query,
        subject,
        options=options,
        losat_bin=_native_losat_bin(),
        threads=threads,
    )


@pytest.mark.linear
def test_native_losat_records_its_runtime_identity(native_losat: dict[str, object]) -> None:
    assert native_losat["kind"] == "losat"
    assert native_losat["version"] == "0.1.0"
    assert native_losat["program"] == "tblastx"
    assert native_losat["cli"] in {"v1", "v2"}


@pytest.mark.linear
def test_native_losatn_reproduces_lambda_de3_bytes(native_losat) -> None:
    query, subject = _lambda_de3()
    expected = (TUTORIAL_DATA / "lambda-de3-comparison" / "lambda-de3.losatn.tsv").read_text(
        encoding="utf-8"
    )

    assert _native_search("losatn", query, subject, LosatSearchArgs(task="megablast")) == expected


@pytest.mark.linear
def test_native_tlosatx_reproduces_lambda_de3_rows_for_any_thread_count(native_losat) -> None:
    query, subject = _lambda_de3()
    options = LosatSearchArgs(query_gencode=1, db_gencode=1)
    expected = (TUTORIAL_DATA / "lambda-de3-comparison" / "lambda-de3.tlosatx.tsv").read_text(
        encoding="utf-8"
    )

    one_thread = _native_search("tlosatx", query, subject, options, threads=1)
    four_threads = _native_search("tlosatx", query, subject, options, threads=4)

    assert len(one_thread.splitlines()) == 397
    assert Counter(one_thread.splitlines()) == Counter(expected.splitlines())
    assert four_threads == one_thread


@pytest.mark.linear
@pytest.mark.parametrize(
    ("query_fasta", "query_gencode", "expected_name"),
    [
        ("NC_002333.2.fna", 2, "danio-human.tlosatx.tsv"),
        ("NC_024511.2.fna", 5, "drosophila-human.tlosatx.tsv"),
        ("NC_001328.1.fna", 5, "caenorhabditis-human.tlosatx.tsv"),
    ],
)
def test_native_tlosatx_reproduces_metazoan_mtdna_bytes(
    native_losat,
    query_fasta: str,
    query_gencode: int,
    expected_name: str,
) -> None:
    fixture = TUTORIAL_DATA / "metazoan-mitochondria-comparison"
    query = (fixture / query_fasta).read_text(encoding="utf-8")
    subject = _fasta(TUTORIAL_DATA / "human-mitochondrion" / "HmmtDNA.gbk", width=60)
    options = LosatSearchArgs(query_gencode=query_gencode, db_gencode=2)

    assert _native_search("tlosatx", query, subject, options) == (
        fixture / expected_name
    ).read_text(encoding="utf-8")
