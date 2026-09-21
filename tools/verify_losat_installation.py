"""Verify a public LOSAT install and offline searches from an installed gbdraw wheel."""

from __future__ import annotations

import argparse
from contextlib import ExitStack
import ctypes
import errno
import hashlib
import json
from pathlib import Path
import platform
import re
import shutil
import socket
import subprocess
import sys
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import gbdraw
from gbdraw import losat_setup as setup
from gbdraw.analysis import collinearity, protein_colinearity as protein


def runtime_identity(expected_mode: str) -> dict:
    translated = False
    if platform.system() == "Darwin":
        # Apple: About the Rosetta translation environment, processIsTranslated.
        # Query this Python process; a separately spawned universal tool may run natively.
        sysctl = ctypes.CDLL(None, use_errno=True).sysctlbyname
        sysctl.argtypes = [ctypes.c_char_p, ctypes.c_void_p, ctypes.POINTER(ctypes.c_size_t),
                          ctypes.c_void_p, ctypes.c_size_t]
        sysctl.restype = ctypes.c_int
        value = ctypes.c_int()
        size = ctypes.c_size_t(ctypes.sizeof(value))
        if sysctl(b"sysctl.proc_translated", ctypes.byref(value), ctypes.byref(size), None, 0) == -1:
            error = ctypes.get_errno()
            if error != errno.ENOENT:
                raise OSError(error, "Cannot determine Rosetta execution")
        else:
            if value.value not in (0, 1):
                raise ValueError("Unexpected Rosetta process identity")
            translated = value.value == 1
    mode = "rosetta" if translated else "native"
    if mode != expected_mode or (translated and platform.machine() != "x86_64"):
        raise ValueError(f"Expected {expected_mode} execution; observed {mode} / {platform.machine()}")
    return {"execution_mode": mode, "python_machine": platform.machine(), "python_executable": sys.executable}


def macos_binary_metadata(binary: Path) -> dict | None:
    if platform.system() != "Darwin":
        return None
    # Inspect the setup result without changing quarantine, signatures or OS policy.
    names = subprocess.check_output(["/usr/bin/xattr", str(binary)], text=True).splitlines()
    attributes = {name: subprocess.check_output(
        ["/usr/bin/xattr", "-p", "-x", name, str(binary)], text=True
    ).replace(" ", "").replace("\n", "").lower() for name in names}
    signed = subprocess.run(["/usr/bin/codesign", "--display", "--verbose=4", str(binary)],
                            capture_output=True, text=True, check=False)
    return {"extended_attributes_hex": attributes,
            "quarantine_present": "com.apple.quarantine" in attributes,
            "codesign": {"returncode": signed.returncode, "stdout": signed.stdout, "stderr": signed.stderr}}


def source_identity(source: Path, declared_sha: str | None) -> str:
    if (source / ".git").exists():
        return subprocess.check_output(
            ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
        ).strip()
    sha = str(declared_sha or "").strip()
    if not re.fullmatch(r"[0-9a-f]{40}", sha):
        raise ValueError(
            "An extracted LOSAT source requires --losat-source-sha with its 40-character commit SHA"
        )
    return sha


def cache_snapshot() -> dict:
    root = setup.cache_root()
    if not root.exists():
        return {"root": str(root), "exists": False, "entries": []}
    entries = sorted(str(path.relative_to(root)) for path in root.rglob("*"))
    return {"root": str(root), "exists": True, "entries": entries}


def conda_subdir() -> str:
    system = platform.system()
    machine = platform.machine().lower()
    if system == "Linux" and machine == "x86_64":
        return "linux-64"
    if system == "Darwin" and machine == "x86_64":
        return "osx-64"
    if system == "Darwin" and machine in {"arm64", "aarch64"}:
        return "osx-arm64"
    raise ValueError(f"Unsupported initial conda target: {system} / {machine}")


def prepare_runtime(installation_mode: str, launcher: list[str]) -> protein.ProteinBlastpRuntime:
    if installation_mode == "managed":
        if setup.managed_losat() is not None:
            raise ValueError("Acceptance requires a fresh target cache for the first download")
        subprocess.run([*launcher, "-m", "gbdraw.cli", "setup-losat"], check=True)

    with ExitStack() as stack:
        runtime = protein._resolve_protein_blastp_runtime("losat", None, stack)

    if installation_mode == "managed":
        installed = setup.managed_losat()
        if installed is None or runtime.source != "managed" or Path(runtime.executable) != installed:
            raise ValueError("Managed setup did not resolve to its installed LOSAT executable")
        return runtime

    expected = (Path(sys.prefix) / "bin" / "losat").absolute()
    if runtime.source != "conda" or Path(runtime.executable) != expected:
        raise ValueError(
            f"Conda acceptance requires source=conda at {expected}; "
            f"observed source={runtime.source} at {runtime.executable}"
        )
    return runtime


def fixture_records(query_path: Path) -> list[SeqRecord]:
    # Six real fixture proteins in two annotated records exercise each consumer.
    proteins = list(SeqIO.parse(query_path, "fasta"))[:6]
    records = []
    for record_index in range(2):
        record = SeqRecord(
            Seq("N" * sum(len(item.seq) * 3 + 30 for item in proteins)),
            id=f"genome{record_index}",
        )
        record.annotations["molecule_type"] = "DNA"
        offset = 0
        for index, item in enumerate(proteins):
            length = len(item.seq) * 3
            record.features.append(
                SeqFeature(
                    FeatureLocation(offset, offset + length, strand=1),
                    type="CDS",
                    qualifiers={
                        "translation": [str(item.seq)],
                        "protein_id": [f"protein{index}"],
                        "locus_tag": [f"gene{index}"],
                    },
                )
            )
            offset += length + 30
        records.append(record)
    return records


def run_cli_smoke(records: list[SeqRecord], output_dir: Path) -> dict:
    cli_dir = output_dir / "cli-pairwise"
    cli_dir.mkdir()
    inputs = []
    for index, record in enumerate(records):
        path = cli_dir / f"record-{index}.gb"
        SeqIO.write(record, path, "genbank")
        inputs.append(path)
    output_prefix = cli_dir / "diagram"
    raw_output = cli_dir / "protein-search.tsv"
    command = [
        sys.executable,
        "-m",
        "gbdraw.cli",
        "linear",
        "--gbk",
        *(str(path) for path in inputs),
        "--protein_blastp_mode",
        "pairwise",
        "--losatp_threads",
        "1",
        "--protein_blastp_output",
        str(raw_output),
        "--format",
        "svg",
        "--output",
        str(output_prefix),
    ]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    (cli_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (cli_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise ValueError(
            f"Installed-package CLI smoke failed with exit code {completed.returncode}: "
            f"{completed.stderr.strip()}"
        )
    rows = [
        line for line in raw_output.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    if not rows or not output_prefix.with_suffix(".svg").is_file():
        raise ValueError("Installed-package CLI smoke did not produce comparison rows and SVG output")
    return {
        "argv": command,
        "returncode": completed.returncode,
        "comparison_rows": len(rows),
        "svg": str(output_prefix.with_suffix(".svg")),
        "raw_output": str(raw_output),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--losat-source", type=Path, required=True)
    parser.add_argument("--losat-source-sha")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--execution-mode", choices=("native", "rosetta"), default="native")
    parser.add_argument(
        "--installation-mode",
        choices=("managed", "conda"),
        default="managed",
    )
    parser.add_argument("--expected-binary-sha256")
    parser.add_argument(
        "--offline-context",
        default="Python socket and managed-download entry points blocked",
    )
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    runtime = runtime_identity(args.execution_mode)
    lock = setup.read_release_lock()
    if args.installation_mode == "managed" and not lock["artifacts"]:
        raise ValueError("Public release pin is required before acceptance")
    source = args.losat_source.resolve()
    source_sha = source_identity(source, args.losat_source_sha)
    if source_sha != lock["candidate_sha"]:
        raise ValueError("Smoke fixture source does not match the pinned release candidate")
    contract = json.loads((source / "docs/release/v0.1.0_rc_contract.json").read_text())
    smoke = contract["smoke"]
    query_path, subject_path = source / smoke["query"], source / smoke["subject"]
    launcher = (["/usr/bin/arch", "-x86_64", sys.executable]
                if runtime["execution_mode"] == "rosetta" else [sys.executable])
    cache_before = cache_snapshot()
    resolved_runtime = prepare_runtime(args.installation_mode, launcher)
    installed = Path(resolved_runtime.executable)
    binary_sha256 = setup.file_sha256(installed)
    expected_binary_sha256 = str(args.expected_binary_sha256 or "").strip()
    if expected_binary_sha256 and binary_sha256 != expected_binary_sha256:
        raise ValueError(
            f"LOSAT binary SHA-256 mismatch: expected {expected_binary_sha256}, observed {binary_sha256}"
        )
    version = subprocess.run(
        [str(installed), "--version"], capture_output=True, text=True, check=True
    )
    target = (
        setup.runtime_target()
        if args.installation_mode == "managed"
        else conda_subdir()
    )
    report = {"status": "PASS", "candidate_sha": source_sha, "target": target,
              "runtime": runtime, "macos_binary_metadata": macos_binary_metadata(installed),
              "platform": platform.platform(), "python": sys.version,
              "installation_mode": args.installation_mode,
              "gbdraw_file": str(Path(gbdraw.__file__).resolve()),
              "sys_executable": sys.executable, "sys_prefix": sys.prefix,
              "resolver": {"kind": resolved_runtime.kind, "source": resolved_runtime.source,
                           "path": resolved_runtime.executable},
              "losat": {"path": str(installed), "sha256": binary_sha256,
                        "version_stdout": version.stdout.strip(),
                        "version_stderr": version.stderr.strip()},
              "offline": {"context": args.offline_context,
                          "python_socket_blocked": True,
                          "managed_download_blocked": True},
              "fixture_hashes": {p.name: setup.file_sha256(p) for p in (query_path, subject_path)}, "cases": []}
    if args.installation_mode == "managed":
        report["installation"] = json.loads((installed.parent / "INSTALL.json").read_text())
    searches = []
    original_search = protein._run_protein_blastp_subprocess

    def recorded_search(command, **kwargs):
        folder = args.output.parent / f"search-{len(searches) + 1:03d}"
        folder.mkdir()
        invocation = {"argv": command, "stdout_representation": "gbdraw text-mode stdout encoded as UTF-8"}
        for flag, name in [("-query", "query.faa"), ("-subject", "subject.faa")]:
            shutil.copyfile(command[command.index(flag) + 1], folder / name)
        (folder / "invocation.json").write_text(json.dumps(invocation, indent=2) + "\n")
        completed = original_search(command, **kwargs)
        (folder / "stdout.txt").write_bytes(completed.stdout.encode("utf-8"))
        (folder / "stderr.txt").write_bytes(completed.stderr.encode("utf-8"))
        searches.append({"directory": folder.name, "argv": command,
                         "returncode": completed.returncode,
                         "stdout_sha256": setup.file_sha256(folder / "stdout.txt")})
        return completed

    records = fixture_records(query_path)
    with patch.object(setup, "urlopen", side_effect=AssertionError("offline download attempted")), patch.object(
        socket, "create_connection", side_effect=AssertionError("offline connection attempted")
    ), patch.object(protein, "_run_protein_blastp_subprocess", side_effect=recorded_search):
        if args.installation_mode == "managed":
            assert setup.setup_losat() == installed
        raw = []
        hits = protein.run_losatp_blastp(query_path.read_text(), subject_path.read_text(), threads=1, raw_output_callback=raw.append)
        digest = hashlib.sha256(raw[0].encode()).hexdigest()
        if digest != smoke["expected_output_sha256"]:
            raise ValueError(f"Public binary BLASTP output differs: {digest}")
        report["cases"].append({"case": smoke["case_id"], "rows": len(hits), "output_sha256": digest})
        for name, builder in [("pairwise", protein.build_pairwise_protein_blastp_comparisons),
                              ("similarity_groups", protein.build_rbh_orthogroup_protein_blastp_comparisons),
                              ("collinear", collinearity.build_orthogroup_collinearity_blocks)]:
            automatic = builder(records, losatp_threads=1)
            explicit = builder(records, losatp_threads=1, losatp_bin=str(installed))
            actual = collinearity.convert_collinearity_blocks_to_comparisons(automatic, records=records) if name == "collinear" else automatic.comparisons
            expected = collinearity.convert_collinearity_blocks_to_comparisons(explicit, records=records) if name == "collinear" else explicit.comparisons
            assert len(actual) == len(expected) and all(a.equals(e) for a, e in zip(actual, expected)), name
            count = sum(len(frame) for frame in actual)
            assert count > 0, name
            for index, (actual_frame, expected_frame) in enumerate(zip(actual, expected)):
                for kind, frame in [("automatic", actual_frame), ("explicit", expected_frame)]:
                    frame.to_csv(args.output.parent / f"{name}-{kind}-{index}.tsv", sep="\t", index=False)
            report["cases"].append({"case": name, "rows": count, "automatic_equals_explicit": True})
        report["cli"] = run_cli_smoke(records, args.output.parent)
    report["searches"] = searches
    cache_after = cache_snapshot()
    report["managed_cache"] = {"before": cache_before, "after": cache_after,
                               "unchanged": cache_before == cache_after}
    if args.installation_mode == "conda" and cache_before != cache_after:
        raise ValueError("Conda acceptance changed the managed LOSAT cache")
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
