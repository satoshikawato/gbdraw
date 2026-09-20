"""Verify a public LOSAT install and offline searches from an installed gbdraw wheel."""

from __future__ import annotations

import argparse
import ctypes
import errno
import hashlib
import json
from pathlib import Path
import platform
import shutil
import socket
import subprocess
import sys
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--losat-source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--execution-mode", choices=("native", "rosetta"), default="native")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    runtime = runtime_identity(args.execution_mode)
    lock = setup.read_release_lock()
    if not lock["artifacts"]:
        raise ValueError("Public release pin is required before acceptance")
    source = args.losat_source.resolve()
    source_sha = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    if source_sha != lock["candidate_sha"]:
        raise ValueError("Smoke fixture source does not match the pinned release candidate")
    contract = json.loads((source / "docs/release/v0.1.0_rc_contract.json").read_text())
    smoke = contract["smoke"]
    query_path, subject_path = source / smoke["query"], source / smoke["subject"]
    if setup.managed_losat() is not None:
        raise ValueError("Acceptance requires a fresh target cache for the first download")
    launcher = (["/usr/bin/arch", "-x86_64", sys.executable]
                if runtime["execution_mode"] == "rosetta" else [sys.executable])
    subprocess.run([*launcher, "-m", "gbdraw.cli", "setup-losat"], check=True)
    installed = setup.managed_losat()
    assert installed is not None
    report = {"status": "PASS", "candidate_sha": source_sha, "target": setup.runtime_target(),
              "runtime": runtime, "macos_binary_metadata": macos_binary_metadata(installed),
              "platform": platform.platform(), "python": sys.version,
              "gbdraw_package": setup.__file__, "installation": json.loads((installed.parent / "INSTALL.json").read_text()),
              "fixture_hashes": {p.name: setup.file_sha256(p) for p in (query_path, subject_path)}, "cases": []}
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
        searches.append({"directory": folder.name, "returncode": completed.returncode,
                         "stdout_sha256": setup.file_sha256(folder / "stdout.txt")})
        return completed

    # Six real fixture proteins in two annotated records exercise each consumer.
    proteins = list(SeqIO.parse(query_path, "fasta"))[:6]
    records = []
    for record_index in range(2):
        record = SeqRecord(Seq("N" * sum(len(p.seq) * 3 + 30 for p in proteins)), id=f"genome{record_index}")
        offset = 0
        for index, item in enumerate(proteins):
            length = len(item.seq) * 3
            record.features.append(SeqFeature(FeatureLocation(offset, offset + length, strand=1), type="CDS",
                qualifiers={"translation": [str(item.seq)], "protein_id": [f"protein{index}"], "locus_tag": [f"gene{index}"]}))
            offset += length + 30
        records.append(record)
    with patch.object(setup, "urlopen", side_effect=AssertionError("offline download attempted")), patch.object(
        socket, "create_connection", side_effect=AssertionError("offline connection attempted")
    ), patch.object(protein, "_run_protein_blastp_subprocess", side_effect=recorded_search):
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
            managed = builder(records, losatp_threads=1)
            explicit = builder(records, losatp_threads=1, losatp_bin=str(installed))
            actual = collinearity.convert_collinearity_blocks_to_comparisons(managed, records=records) if name == "collinear" else managed.comparisons
            expected = collinearity.convert_collinearity_blocks_to_comparisons(explicit, records=records) if name == "collinear" else explicit.comparisons
            assert len(actual) == len(expected) and all(a.equals(e) for a, e in zip(actual, expected)), name
            count = sum(len(frame) for frame in actual)
            assert count > 0, name
            for index, (actual_frame, expected_frame) in enumerate(zip(actual, expected)):
                for kind, frame in [("managed", actual_frame), ("explicit", expected_frame)]:
                    frame.to_csv(args.output.parent / f"{name}-{kind}-{index}.tsv", sep="\t", index=False)
            report["cases"].append({"case": name, "rows": count, "managed_equals_explicit": True})
    report["searches"] = searches
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
