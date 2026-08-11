#!/usr/bin/env python3
"""Qualify frozen browser LOSAT evidence for whole Lambda and DE3 genomes."""

from __future__ import annotations

import argparse
import hashlib
import json
import threading
from collections import Counter
from dataclasses import dataclass
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

from Bio import SeqIO


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
FIXTURE_ROOT = WEB_ROOT / "tutorial-data"
OUTPUT_ROOT = FIXTURE_ROOT / "lambda-de3-comparison"
WASM_PATH = WEB_ROOT / "wasm" / "losat" / "losat.wasm"
RUNTIME_VERSION = "0.1.0"
RUNTIME_SHA256 = "1c597c9e366785051af698cb69f315bec10ac6c55e7a3b9c2e9325dbbd46f312"
QUERY_ID = "NC_001416.1"
SUBJECT_ID = "NC_042057.1"
OUTFMT6_COLUMN_COUNT = 12
DISPLAY_THRESHOLDS = {
    "evalueMaximum": 0.01,
    "bitscoreMinimum": 50.0,
    "identityMinimum": 0.0,
    "alignmentLengthMinimum": 0,
}


@dataclass(frozen=True)
class RecordSpec:
    record_id: str
    path: Path
    length: int
    file_sha256: str
    sequence_sha256: str
    browser_fasta_sha256: str


@dataclass(frozen=True)
class ProgramSpec:
    name: str
    runtime_program: str
    extra_args: tuple[str, ...]
    expected_size: int | None
    expected_sha256: str | None
    expected_rows: int | None
    expected_retained_rows: int | None

    @property
    def output_path(self) -> Path:
        return OUTPUT_ROOT / f"lambda-de3.{self.name}.tsv"


RECORDS = {
    QUERY_ID: RecordSpec(
        record_id=QUERY_ID,
        path=FIXTURE_ROOT / "lambda" / "NC_001416.gb",
        length=48_502,
        file_sha256="4b76b8bacc8026aac3f19a4a915f4ac772ad61e7ec18f0e2cc859229f95a66e7",
        sequence_sha256="36432a40f602258d19ae7c8152ddbc30390b559f2859c01d7047c77b048c71b3",
        browser_fasta_sha256="1471a46a2b559a7c0cd6d018c442044e92aff8b4b1c4f6f1ab888efd44d20329",
    ),
    SUBJECT_ID: RecordSpec(
        record_id=SUBJECT_ID,
        path=FIXTURE_ROOT / "de3" / "NC_042057.1.gb",
        length=42_925,
        file_sha256="288eb87480f8fe6eab6246fe1fc4af78c85a6cb7e591c47fd7d7f0170c932e09",
        sequence_sha256="533c360de52e4d5ac1f5dc31e0b927a518ff04951676d287b4b6d0238758aa45",
        browser_fasta_sha256="c89d0dbeeb4552822e7c79f166ae65a66abb578378d2a8189febb98ac1134117",
    ),
}
PROGRAMS = (
    ProgramSpec(
        name="losatn",
        runtime_program="blastn",
        extra_args=("--task", "megablast"),
        expected_size=436,
        expected_sha256="703b0ac749669152a8e6d5fa6fb246cf5973cb1a3e0ca9db2f346ee33628317c",
        expected_rows=6,
        expected_retained_rows=6,
    ),
    ProgramSpec(
        name="tlosatx",
        runtime_program="tblastx",
        extra_args=("--query-gencode", "1", "--db-gencode", "1"),
        expected_size=28_303,
        expected_sha256="483e98f8b3dce172523cf00c82eb9d47e4faee13e4781b06f3d58ea4fb63532d",
        expected_rows=397,
        expected_retained_rows=266,
    ),
)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _load_record(spec: RecordSpec):
    payload = spec.path.read_bytes()
    if _sha256(payload) != spec.file_sha256:
        raise ValueError(f"source checksum mismatch: {spec.path.relative_to(REPO_ROOT)}")
    record = SeqIO.read(spec.path, "genbank")
    sequence = str(record.seq).upper().encode("ascii")
    if record.id != spec.record_id or len(record) != spec.length:
        raise ValueError(f"source record identity mismatch: {spec.record_id}")
    if _sha256(sequence) != spec.sequence_sha256:
        raise ValueError(f"source sequence checksum mismatch: {spec.record_id}")
    return record


def _browser_fasta(spec: RecordSpec) -> str:
    record = _load_record(spec)
    sequence = str(record.seq).upper()
    lines = [f">{record.id}"]
    lines.extend(sequence[index : index + 60] for index in range(0, len(sequence), 60))
    fasta = "\n".join(lines) + "\n"
    if _sha256(fasta.encode("ascii")) != spec.browser_fasta_sha256:
        raise ValueError(f"browser FASTA checksum mismatch: {spec.record_id}")
    return fasta


def _parse_rows(data: bytes) -> list[tuple[str, ...]]:
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("LOSAT evidence is not UTF-8") from error
    if not text or not text.endswith("\n"):
        raise ValueError("LOSAT evidence must be non-empty and newline terminated")
    rows: list[tuple[str, ...]] = []
    for line_number, line in enumerate(text.splitlines(), start=1):
        columns = tuple(line.split("\t"))
        if len(columns) != OUTFMT6_COLUMN_COUNT:
            raise ValueError(f"LOSAT row {line_number} does not have 12 columns")
        rows.append(columns)
    return rows


def _retained(row: tuple[str, ...]) -> bool:
    return (
        float(row[10]) <= DISPLAY_THRESHOLDS["evalueMaximum"]
        and float(row[11]) >= DISPLAY_THRESHOLDS["bitscoreMinimum"]
        and float(row[2]) >= DISPLAY_THRESHOLDS["identityMinimum"]
        and int(row[3]) >= DISPLAY_THRESHOLDS["alignmentLengthMinimum"]
    )


def _summarize(data: bytes) -> dict[str, Any]:
    rows = _parse_rows(data)
    retained = [row for row in rows if _retained(row)]
    orientations: Counter[str] = Counter()
    for line_number, row in enumerate(rows, start=1):
        if row[0] != QUERY_ID or row[1] != SUBJECT_ID:
            raise ValueError(f"LOSAT row {line_number} has the wrong endpoints")
        identity = float(row[2])
        alignment_length = int(row[3])
        mismatches = int(row[4])
        gap_opens = int(row[5])
        qstart, qend, sstart, send = (int(value) for value in row[6:10])
        evalue, bitscore = float(row[10]), float(row[11])
        if not 0 <= identity <= 100 or alignment_length <= 0:
            raise ValueError(f"LOSAT row {line_number} has invalid alignment values")
        if mismatches < 0 or gap_opens < 0 or evalue < 0 or bitscore < 0:
            raise ValueError(f"LOSAT row {line_number} has invalid score values")
        if not 1 <= qstart <= RECORDS[QUERY_ID].length or not 1 <= qend <= RECORDS[QUERY_ID].length:
            raise ValueError(f"LOSAT row {line_number} is outside the complete Lambda record")
        if not 1 <= sstart <= RECORDS[SUBJECT_ID].length or not 1 <= send <= RECORDS[SUBJECT_ID].length:
            raise ValueError(f"LOSAT row {line_number} is outside the complete DE3 record")
        orientations[("+" if qend >= qstart else "-") + ("+" if send >= sstart else "-")] += 1
    return {
        "sizeBytes": len(data),
        "sha256": _sha256(data),
        "rowCount": len(rows),
        "retainedRowCount": len(retained),
        "queryCoordinateRange": [
            min(min(int(row[6]), int(row[7])) for row in rows),
            max(max(int(row[6]), int(row[7])) for row in rows),
        ],
        "subjectCoordinateRange": [
            min(min(int(row[8]), int(row[9])) for row in rows),
            max(max(int(row[8]), int(row[9])) for row in rows),
        ],
        "orientationCounts": dict(sorted(orientations.items())),
    }


def _validate(spec: ProgramSpec, data: bytes, *, frozen: bool) -> dict[str, Any]:
    summary = _summarize(data)
    expected = {
        "sizeBytes": spec.expected_size,
        "sha256": spec.expected_sha256,
        "rowCount": spec.expected_rows,
        "retainedRowCount": spec.expected_retained_rows,
    }
    if frozen:
        missing = [key for key, value in expected.items() if value is None]
        if missing:
            raise ValueError(f"{spec.name} has no frozen {', '.join(missing)}")
        for key, value in expected.items():
            if summary[key] != value:
                raise ValueError(
                    f"{spec.name} {key} changed: {summary[key]} != {value}"
                )
    if summary["retainedRowCount"] == 0:
        raise ValueError(f"{spec.name} produced no displayed evidence")
    return summary


class _Handler(SimpleHTTPRequestHandler):
    def do_GET(self) -> None:  # noqa: N802 - inherited HTTP handler API
        if self.path == "/__lambda_de3_qualification__.html":
            body = b"<!doctype html><meta charset=utf-8><title>Lambda DE3 qualification</title>"
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        super().do_GET()

    def log_message(self, _format: str, *_args: object) -> None:
        return


def _check_runtime() -> None:
    payload = WASM_PATH.read_bytes()
    if _sha256(payload) != RUNTIME_SHA256 or RUNTIME_VERSION.encode("ascii") not in payload:
        raise ValueError("pinned browser LOSAT runtime changed")


def _run_browser(runs: int) -> tuple[dict[str, list[bytes]], dict[str, list[float]]]:
    try:
        from playwright.sync_api import sync_playwright
    except ImportError as error:
        raise RuntimeError("Python Playwright with Chromium is required") from error

    query_fasta = _browser_fasta(RECORDS[QUERY_ID])
    subject_fasta = _browser_fasta(RECORDS[SUBJECT_ID])
    handler = partial(_Handler, directory=str(WEB_ROOT))
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    origin = f"http://127.0.0.1:{server.server_port}"
    outputs = {spec.name: [] for spec in PROGRAMS}
    timings = {spec.name: [] for spec in PROGRAMS}
    try:
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(headless=True)
            try:
                page = browser.new_page()
                page.goto(
                    f"{origin}/__lambda_de3_qualification__.html",
                    wait_until="domcontentloaded",
                )
                for spec in PROGRAMS:
                    for _ in range(runs):
                        result = page.evaluate(
                            """async ({ program, queryFasta, subjectFasta, extraArgs }) => {
                                const { runLosatPair } = await import('/js/services/losat.js');
                                const startedAt = performance.now();
                                const text = await runLosatPair({
                                    program,
                                    queryFasta,
                                    subjectFasta,
                                    outfmt: '6',
                                    extraArgs
                                });
                                return { text, elapsedMs: performance.now() - startedAt };
                            }""",
                            {
                                "program": spec.runtime_program,
                                "queryFasta": query_fasta,
                                "subjectFasta": subject_fasta,
                                "extraArgs": list(spec.extra_args),
                            },
                        )
                        outputs[spec.name].append(result["text"].encode("utf-8"))
                        timings[spec.name].append(float(result["elapsedMs"]) / 1_000)
            finally:
                browser.close()
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
    return outputs, timings


def _check() -> dict[str, dict[str, Any]]:
    _check_runtime()
    for spec in RECORDS.values():
        _load_record(spec)
    return {
        spec.name: _validate(spec, spec.output_path.read_bytes(), frozen=True)
        for spec in PROGRAMS
    }


def _qualify(*, runs: int, write: bool) -> dict[str, dict[str, Any]]:
    if runs < 2:
        raise ValueError("qualification requires at least two serial runs")
    _check_runtime()
    outputs, timings = _run_browser(runs)
    summaries: dict[str, dict[str, Any]] = {}
    for spec in PROGRAMS:
        candidates = outputs[spec.name]
        if any(candidate != candidates[0] for candidate in candidates[1:]):
            raise ValueError(f"{spec.name} was not byte-deterministic")
        summary = _validate(spec, candidates[0], frozen=not write)
        summary["serialRunSeconds"] = [round(value, 3) for value in timings[spec.name]]
        summaries[spec.name] = summary
        if write:
            OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
            spec.output_path.write_bytes(candidates[0])
        elif spec.output_path.read_bytes() != candidates[0]:
            raise ValueError(f"frozen {spec.name} evidence differs from the browser runtime")
    return summaries


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--qualify", action="store_true")
    action.add_argument("--write", action="store_true")
    parser.add_argument("--runs", type=int, default=2)
    args = parser.parse_args(argv)
    try:
        summaries = (
            _qualify(runs=args.runs, write=args.write)
            if args.qualify or args.write
            else _check()
        )
    except (FileNotFoundError, OSError, RuntimeError, ValueError) as error:
        parser.exit(1, f"qualification failed: {error}\n")
    print(json.dumps(summaries, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
