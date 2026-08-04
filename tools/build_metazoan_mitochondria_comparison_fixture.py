#!/usr/bin/env python3
"""Build frozen TLOSATX evidence for the complete metazoan mtDNA examples."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from Bio import SeqIO


REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_ROOT = REPO_ROOT / "gbdraw" / "web" / "tutorial-data"
OUTPUT_ROOT = FIXTURE_ROOT / "metazoan-mitochondria-comparison"
LOSAT_PATH = REPO_ROOT / "gbdraw" / "bin" / "linux-x86_64" / "losat"
LOSAT_SHA256 = "3d357ae9062d5ed4b2fa591e0316c2dd67fb7b3b986a079902bfa078ee18d3dd"
THRESHOLDS = {
    "evalue": 1e-5,
    "bitscore": 50.0,
    "identity": 40.0,
    "alignment_length": 50,
}


@dataclass(frozen=True)
class RecordSpec:
    path: Path
    record_id: str
    length: int
    sha256: str
    genetic_code: int


@dataclass(frozen=True)
class ComparisonSpec:
    query: RecordSpec
    filename: str
    size: int
    sha256: str
    raw_rows: int
    retained_rows: int


HUMAN = RecordSpec(
    FIXTURE_ROOT / "human-mitochondrion" / "HmmtDNA.gbk",
    "NC_012920.1",
    16_569,
    "7530d659e7174272372814edfecb2ece1f87a444395a861fcdf1b977c4aa5c1f",
    2,
)
COMPARISONS = (
    ComparisonSpec(
        RecordSpec(
            FIXTURE_ROOT / "metazoan-mitochondria-four" / "NC_002333.2.gb",
            "NC_002333.2",
            16_596,
            "94ab35da6f81abc2595fcd425c23585ed78d9396b5143918d9d1025d8a4d2140",
            2,
        ),
        "danio-human.tlosatx.tsv",
        19_284,
        "4135025ba9fc6f346551757ad3d842b30c2f813bb2a412882d571997ff95c597",
        276,
        68,
    ),
    ComparisonSpec(
        RecordSpec(
            FIXTURE_ROOT / "metazoan-mitochondria-four" / "NC_024511.2.gb",
            "NC_024511.2",
            19_524,
            "79fa36199682961919c4a11f3a8fc50c9e598e68b867eb25e847bce1aa1c4229",
            5,
        ),
        "drosophila-human.tlosatx.tsv",
        6_587,
        "1d9287499564e24bc52063625be977c71976c2e5f2bad236ec08733584358117",
        93,
        24,
    ),
    ComparisonSpec(
        RecordSpec(
            FIXTURE_ROOT / "metazoan-mitochondria-four" / "NC_001328.1.gb",
            "NC_001328.1",
            13_794,
            "8de5f7cf3686f493ee5b8068dba39b31c5d02e70d997063f65ed19d0fa859a59",
            5,
        ),
        "caenorhabditis-human.tlosatx.tsv",
        4_681,
        "b7244a0154965935d60694549a8d2e779c8f0788737b2a1f4f1f96ebc57bd416",
        66,
        14,
    ),
)
COMPARISON_FASTA_FILENAMES = {
    "NC_002333.2": "NC_002333.2.fna",
    "NC_024511.2": "NC_024511.2.fna",
    "NC_001328.1": "NC_001328.1.fna",
}


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _load_record(spec: RecordSpec):
    payload = spec.path.read_bytes()
    if _sha256(payload) != spec.sha256:
        raise ValueError(f"source checksum mismatch: {spec.path}")
    record = SeqIO.read(spec.path, "genbank")
    if record.id != spec.record_id or len(record) != spec.length:
        raise ValueError(f"source identity mismatch: {spec.record_id}")
    if record.annotations.get("topology") != "circular":
        raise ValueError(f"source is not circular: {spec.record_id}")
    if "complete genome" not in record.description:
        raise ValueError(f"source is not complete: {spec.record_id}")
    return record


def _render_fasta(spec: RecordSpec) -> bytes:
    record = _load_record(spec)
    sequence = str(record.seq).upper()
    lines = [f">{record.id}"]
    lines.extend(sequence[index : index + 60] for index in range(0, len(sequence), 60))
    return ("\n".join(lines) + "\n").encode("ascii")


def _write_fasta(spec: RecordSpec, path: Path) -> None:
    path.write_bytes(_render_fasta(spec))


def _retained(row: tuple[str, ...]) -> bool:
    return (
        float(row[10]) <= THRESHOLDS["evalue"]
        and float(row[11]) >= THRESHOLDS["bitscore"]
        and float(row[2]) >= THRESHOLDS["identity"]
        and int(row[3]) >= THRESHOLDS["alignment_length"]
    )


def _summarize(data: bytes, spec: ComparisonSpec) -> dict[str, Any]:
    text = data.decode("utf-8")
    if not text.endswith("\n"):
        raise ValueError(f"output is not newline terminated: {spec.filename}")
    rows = [tuple(line.split("\t")) for line in text.splitlines()]
    if any(len(row) != 12 for row in rows):
        raise ValueError(f"output is not 12-column outfmt 6: {spec.filename}")
    orientations: Counter[str] = Counter()
    for row in rows:
        if row[:2] != (spec.query.record_id, HUMAN.record_id):
            raise ValueError(f"endpoint mismatch: {spec.filename}")
        qstart, qend, sstart, send = (int(value) for value in row[6:10])
        if not 1 <= min(qstart, qend) <= max(qstart, qend) <= spec.query.length:
            raise ValueError(f"query coordinate outside complete record: {spec.filename}")
        if not 1 <= min(sstart, send) <= max(sstart, send) <= HUMAN.length:
            raise ValueError(f"subject coordinate outside complete record: {spec.filename}")
        orientations[("+" if qend >= qstart else "-") + ("+" if send >= sstart else "-")] += 1
    summary = {
        "sizeBytes": len(data),
        "sha256": _sha256(data),
        "rawRowCount": len(rows),
        "retainedRowCount": sum(_retained(row) for row in rows),
        "orientationCounts": dict(sorted(orientations.items())),
    }
    expected = {
        "sizeBytes": spec.size,
        "sha256": spec.sha256,
        "rawRowCount": spec.raw_rows,
        "retainedRowCount": spec.retained_rows,
    }
    for key, value in expected.items():
        if summary[key] != value:
            raise ValueError(
                f"{spec.filename} {key} changed: {summary[key]} != {value}"
            )
    return summary


def _run_once(spec: ComparisonSpec, directory: Path) -> bytes:
    directory.mkdir(parents=True, exist_ok=True)
    query_path = directory / "query.fna"
    subject_path = directory / "subject.fna"
    output_path = directory / spec.filename
    _write_fasta(spec.query, query_path)
    _write_fasta(HUMAN, subject_path)
    subprocess.run(
        [
            str(LOSAT_PATH),
            "tblastx",
            "--query",
            str(query_path),
            "--subject",
            str(subject_path),
            "--query-gencode",
            str(spec.query.genetic_code),
            "--db-gencode",
            str(HUMAN.genetic_code),
            "--num-threads",
            "1",
            "--out",
            str(output_path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return output_path.read_bytes()


def build(*, write: bool) -> dict[str, dict[str, Any]]:
    if _sha256(LOSAT_PATH.read_bytes()) != LOSAT_SHA256:
        raise ValueError("pinned LOSAT executable checksum changed")
    version = subprocess.run(
        [str(LOSAT_PATH), "--version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if version != "losat 0.1.0":
        raise ValueError(f"pinned LOSAT version changed: {version}")

    summaries: dict[str, dict[str, Any]] = {}
    with TemporaryDirectory(prefix="gbdraw-mtdna-tlosatx-") as temp_name:
        temp_root = Path(temp_name)
        for spec in COMPARISONS:
            first = _run_once(spec, temp_root / f"{spec.query.record_id}-1")
            second = _run_once(spec, temp_root / f"{spec.query.record_id}-2")
            if first != second:
                raise ValueError(f"non-deterministic LOSAT output: {spec.filename}")
            summaries[spec.filename] = _summarize(first, spec)
            destination = OUTPUT_ROOT / spec.filename
            if write:
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(first)
            elif destination.read_bytes() != first:
                raise ValueError(f"frozen output differs: {spec.filename}")
    for spec in COMPARISONS:
        filename = COMPARISON_FASTA_FILENAMES[spec.query.record_id]
        payload = _render_fasta(spec.query)
        destination = OUTPUT_ROOT / filename
        if write:
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(payload)
        elif destination.read_bytes() != payload:
            raise ValueError(f"frozen FASTA differs: {filename}")
        summaries[filename] = {
            "sizeBytes": len(payload),
            "sha256": _sha256(payload),
            "recordId": spec.query.record_id,
            "sequenceLength": spec.query.length,
        }
    return summaries


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    try:
        summaries = build(write=args.write)
    except (FileNotFoundError, OSError, subprocess.SubprocessError, ValueError) as error:
        parser.exit(1, f"build failed: {error}\n")
    print(json.dumps(summaries, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
