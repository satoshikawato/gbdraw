#!/usr/bin/env python3
"""Qualify frozen browser LOSAT evidence for whole Lambda and DE3 genomes.

The tutorial fixture manifest (``gbdraw/web/tutorial-data/manifest.json``) owns the frozen values: each evidence
entry whose derivation names this script holds its search (``derivation.arguments``), the runtime that qualified it
(``derivation.toolVersion``, ``derivation.toolSha256``) and what its TSV holds (``sizeBytes``, ``sha256``, the row
counts, coordinate ranges and orientations).

- No option: check, without a browser, that each TSV matches its entry and that its entry names the shipped serial
  LOSAT file (``losat-build.json``).
- ``--qualify``: run each search in Chromium with the shipped serial file and require the frozen bytes.
- ``--write``: run each search in Chromium and write each TSV and its manifest entry together. Run it after the
  LOSAT Wasm files change (``gbdraw/web/wasm/losat/README.md``).
"""

from __future__ import annotations

import argparse
import copy
import dataclasses
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
WASM_PATH = WEB_ROOT / "wasm" / "losat" / "losat.wasm"
# The shipped serial file's build (tools/record_losat_wasm_build.py) and the runtime each frozen output was
# qualified with (its derivation in the tutorial fixture manifest).
BUILD_MANIFEST_PATH = WASM_PATH.with_name("losat-build.json")
FIXTURE_MANIFEST_PATH = FIXTURE_ROOT / "manifest.json"
# The one value of a source record the manifest does not hold: the FASTA text the browser search reads, which this
# script builds from the record (60 columns, upper case).
BROWSER_FASTA_SHA256 = {
    "NC_001416.1": "1471a46a2b559a7c0cd6d018c442044e92aff8b4b1c4f6f1ab888efd44d20329",
    "NC_042057.1": "c89d0dbeeb4552822e7c79f166ae65a66abb578378d2a8189febb98ac1134117",
}
OUTFMT6_COLUMN_COUNT = 12
DISPLAY_THRESHOLDS = {
    "evalueMaximum": 0.01,
    "bitscoreMinimum": 50.0,
    "identityMinimum": 0.0,
    "alignmentLengthMinimum": 0,
}


@dataclass(frozen=True)
class RecordSpec:
    """A source record of the searches, from its tutorial fixture manifest entry (``derivation.sourceFileIds``)."""

    record_id: str
    path: Path
    length: int
    file_sha256: str
    sequence_sha256: str

    @classmethod
    def from_entry(cls, entry: dict[str, Any]) -> RecordSpec:
        (record,) = entry["records"]
        return cls(
            record_id=str(record["id"]),
            path=FIXTURE_ROOT / entry["relativePath"],
            length=int(record["length"]),
            file_sha256=str(entry["sha256"]),
            sequence_sha256=str(record["sequenceSha256"]),
        )


@dataclass(frozen=True)
class Evidence:
    """One frozen TSV and its tutorial fixture manifest entry."""

    file_id: str
    entry: dict[str, Any]
    # The query and the subject record, in the order of the entry's ``derivation.sourceFileIds``.
    sources: tuple[RecordSpec, RecordSpec]

    @property
    def output_path(self) -> Path:
        return FIXTURE_ROOT / self.entry["relativePath"]

    @property
    def runtime_program(self) -> str:
        return str(self.entry["derivation"]["arguments"]["program"])

    @property
    def extra_args(self) -> list[str]:
        """The search's options in the Web's raw-cache spelling, from its recorded arguments."""

        arguments = self.entry["derivation"]["arguments"]
        if self.runtime_program == "blastn":
            return ["--task", str(arguments["task"])]
        if self.runtime_program == "tblastx":
            return [
                "--query-gencode", str(arguments["queryGeneticCode"]),
                "--db-gencode", str(arguments["subjectGeneticCode"]),
            ]
        raise ValueError(f"{self.file_id} records program {self.runtime_program}, which this script does not run")


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
    if _sha256(fasta.encode("ascii")) != BROWSER_FASTA_SHA256[spec.record_id]:
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


def _summarize(data: bytes, evidence: Evidence) -> dict[str, Any]:
    query, subject = evidence.sources
    rows = _parse_rows(data)
    retained = [row for row in rows if _retained(row)]
    orientations: Counter[str] = Counter()
    for line_number, row in enumerate(rows, start=1):
        if row[0] != query.record_id or row[1] != subject.record_id:
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
        if not 1 <= qstart <= query.length or not 1 <= qend <= query.length:
            raise ValueError(f"LOSAT row {line_number} is outside the complete {query.record_id} record")
        if not 1 <= sstart <= subject.length or not 1 <= send <= subject.length:
            raise ValueError(f"LOSAT row {line_number} is outside the complete {subject.record_id} record")
        orientations[("+" if qend >= qstart else "-") + ("+" if send >= sstart else "-")] += 1
    summary: dict[str, Any] = {
        "sizeBytes": len(data),
        "sha256": _sha256(data),
        "rawRowCount": len(rows),
        "defaultRetainedRowCount": len(retained),
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
    minimum = evidence.entry["expectedSemantics"].get("tutorialMinimumAlignmentLength")
    if minimum is not None:
        summary["tutorialRetainedRowCount"] = sum(int(row[3]) >= minimum for row in rows)
    return summary


def _frozen(entry: dict[str, Any]) -> dict[str, Any]:
    """What the manifest entry records of its TSV, keyed as `_summarize` names it."""

    semantics = entry["expectedSemantics"]
    frozen = {"sizeBytes": entry["sizeBytes"], "sha256": entry["sha256"]}
    for key in ("rawRowCount", "defaultRetainedRowCount", "queryCoordinateRange", "subjectCoordinateRange",
                "orientationCounts", "tutorialRetainedRowCount"):
        if key in semantics:
            frozen[key] = semantics[key]
    return frozen


def _validate(evidence: Evidence, data: bytes, *, frozen: bool) -> dict[str, Any]:
    summary = _summarize(data, evidence)
    if frozen:
        for key, value in _frozen(evidence.entry).items():
            if summary[key] != value:
                raise ValueError(f"{evidence.file_id} {key} changed: {summary[key]} != {value}")
    if summary["defaultRetainedRowCount"] == 0:
        raise ValueError(f"{evidence.file_id} produced no displayed evidence")
    return summary


def _value_spans(text: str) -> dict[tuple[str | int, ...], tuple[int, int]]:
    """Where each value of a JSON text starts and ends, by its path of keys and indexes."""

    decoder = json.JSONDecoder()
    spans: dict[tuple[str | int, ...], tuple[int, int]] = {}

    def skip(index: int) -> int:
        while text[index] in " \t\r\n":
            index += 1
        return index

    def value(index: int, path: tuple[str | int, ...]) -> int:
        start = index = skip(index)
        if text[index] in "{[":
            closing, position = ("}" if text[index] == "{" else "]"), 0
            index = skip(index + 1)
            while text[index] != closing:
                key: str | int = position
                if closing == "}":
                    key, index = decoder.raw_decode(text, index)
                    index = skip(index) + 1  # the colon
                index = skip(value(index, (*path, key)))
                position += 1
                if text[index] == ",":
                    index = skip(index + 1)
            index += 1
        else:
            index = decoder.raw_decode(text, index)[1]
        spans[path] = (start, index)
        return index

    value(0, ())
    return spans


def _changed_values(old: Any, new: Any, path: tuple[str | int, ...] = ()) -> list[tuple[tuple[str | int, ...], Any]]:
    """The values that differ between two JSON documents of one shape, each at its deepest path that still exists in
    both."""

    if isinstance(old, dict) and isinstance(new, dict) and old.keys() == new.keys():
        return [change for key in old for change in _changed_values(old[key], new[key], (*path, key))]
    if isinstance(old, list) and isinstance(new, list) and len(old) == len(new):
        return [change for index, item in enumerate(old) for change in _changed_values(item, new[index], (*path, index))]
    return [] if old == new else [(path, new)]


def _write_manifest(evidence: list[Evidence], build: dict[str, Any], runs: int) -> None:
    """Record each evidence TSV as written into its manifest entry, changing only the values that differ, so the
    manifest keeps its own layout."""

    text = FIXTURE_MANIFEST_PATH.read_text(encoding="utf-8")
    old = json.loads(text)
    new = copy.deepcopy(old)
    for item in evidence:
        written = dataclasses.replace(item, entry=new["files"][item.file_id])
        _record(written, _summarize(item.output_path.read_bytes(), written), build, runs)
    spans = _value_spans(text)
    for path, value in sorted(_changed_values(old, new), key=lambda change: spans[change[0]][0], reverse=True):
        start, end = spans[path]
        text = text[:start] + json.dumps(value) + text[end:]
    if json.loads(text) != new:
        raise ValueError(f"{FIXTURE_MANIFEST_PATH.name} could not be updated in place")
    FIXTURE_MANIFEST_PATH.write_text(text, encoding="utf-8")


def _sources(manifest: dict[str, Any], file_id: str, derivation: dict[str, Any]) -> tuple[RecordSpec, RecordSpec]:
    """The query and the subject record an evidence entry names (``derivation.sourceFileIds``)."""

    query_id, subject_id = derivation["sourceFileIds"]
    entries = [manifest["files"][query_id], manifest["files"][subject_id]]
    if derivation["sourceChecksums"] != [entry["sha256"] for entry in entries]:
        raise ValueError(f"{file_id} records source checksums that are not its source files' sha256")
    query, subject = (RecordSpec.from_entry(entry) for entry in entries)
    return query, subject


def _evidence(manifest: dict[str, Any]) -> list[Evidence]:
    """The manifest entries this script derives."""

    script = Path(__file__).resolve().relative_to(REPO_ROOT).as_posix()
    found = [
        Evidence(file_id, entry, _sources(manifest, file_id, entry["derivation"]))
        for file_id, entry in manifest["files"].items()
        if (entry.get("derivation") or {}).get("script") == script
    ]
    if not found:
        raise ValueError(f"{FIXTURE_MANIFEST_PATH.name} names no file derived by {script}")
    return found


def _record(evidence: Evidence, summary: dict[str, Any], build: dict[str, Any], runs: int) -> None:
    """Write a qualified summary into its manifest entry: the TSV's bytes, its rows and the runtime."""

    entry = evidence.entry
    entry["sizeBytes"], entry["sha256"] = summary["sizeBytes"], summary["sha256"]
    semantics = entry["expectedSemantics"]
    for key in _frozen(entry):
        if key in semantics:
            semantics[key] = summary[key]
    derivation = entry["derivation"]
    derivation["toolVersion"] = build["version"]
    derivation["toolSha256"] = build["files"][WASM_PATH.name]["sha256"]
    derivation["arguments"]["runs"] = runs
    statistics = derivation["outputStatistics"]
    for key in ("rawRowCount", "defaultRetainedRowCount", "tutorialRetainedRowCount"):
        if key in statistics:
            statistics[key] = summary[key]
    statistics["byteIdenticalRuns"] = True


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


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _check_runtime(build: dict[str, Any], evidence: list[Evidence] | None) -> None:
    """The serial file must be the recorded build and, unless the evidence is being rewritten (`evidence` None), the
    runtime that qualified each frozen output."""

    shipped = build["files"][WASM_PATH.name]
    if _sha256(WASM_PATH.read_bytes()) != shipped["sha256"]:
        raise ValueError(f"{WASM_PATH.name} does not match {BUILD_MANIFEST_PATH.name}")
    for item in evidence or []:
        derivation = item.entry["derivation"]
        if (derivation["toolVersion"], derivation["toolSha256"]) != (build["version"], shipped["sha256"]):
            raise ValueError(
                f"frozen {item.file_id} evidence was qualified with browser LOSAT {derivation['toolVersion']} "
                f"(sha256 {derivation['toolSha256']}); the shipped file is LOSAT {build['version']} "
                f"(sha256 {shipped['sha256']}); qualify it again with --write"
            )


def _run_browser(evidence: list[Evidence], runs: int) -> tuple[dict[str, list[bytes]], dict[str, list[float]]]:
    try:
        from playwright.sync_api import sync_playwright
    except ImportError as error:
        raise RuntimeError("Python Playwright with Chromium is required") from error

    fasta = {spec.record_id: _browser_fasta(spec) for item in evidence for spec in item.sources}
    handler = partial(_Handler, directory=str(WEB_ROOT))
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    origin = f"http://127.0.0.1:{server.server_port}"
    outputs: dict[str, list[bytes]] = {item.file_id: [] for item in evidence}
    timings: dict[str, list[float]] = {item.file_id: [] for item in evidence}
    try:
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(headless=True)
            try:
                page = browser.new_page()
                page.goto(
                    f"{origin}/__lambda_de3_qualification__.html",
                    wait_until="domcontentloaded",
                )
                for item in evidence:
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
                                "program": item.runtime_program,
                                "queryFasta": fasta[item.sources[0].record_id],
                                "subjectFasta": fasta[item.sources[1].record_id],
                                "extraArgs": item.extra_args,
                            },
                        )
                        outputs[item.file_id].append(result["text"].encode("utf-8"))
                        timings[item.file_id].append(float(result["elapsedMs"]) / 1_000)
            finally:
                browser.close()
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
    return outputs, timings


def _check() -> dict[str, dict[str, Any]]:
    """Check the frozen TSVs against their manifest entries and the shipped serial file, without a browser."""

    evidence = _evidence(_load_json(FIXTURE_MANIFEST_PATH))
    _check_runtime(_load_json(BUILD_MANIFEST_PATH), evidence)
    for spec in {spec for item in evidence for spec in item.sources}:
        _browser_fasta(spec)
    return {item.file_id: _validate(item, item.output_path.read_bytes(), frozen=True) for item in evidence}


def _qualify(*, runs: int, write: bool) -> dict[str, dict[str, Any]]:
    if runs < 2:
        raise ValueError("qualification requires at least two serial runs")
    manifest = _load_json(FIXTURE_MANIFEST_PATH)
    build = _load_json(BUILD_MANIFEST_PATH)
    evidence = _evidence(manifest)
    _check_runtime(build, None if write else evidence)
    outputs, timings = _run_browser(evidence, runs)
    summaries: dict[str, dict[str, Any]] = {}
    for item in evidence:
        candidates = outputs[item.file_id]
        if any(candidate != candidates[0] for candidate in candidates[1:]):
            raise ValueError(f"{item.file_id} was not byte-deterministic")
        summary = _validate(item, candidates[0], frozen=not write)
        if write:
            item.output_path.parent.mkdir(parents=True, exist_ok=True)
            item.output_path.write_bytes(candidates[0])
        elif item.output_path.read_bytes() != candidates[0]:
            raise ValueError(f"frozen {item.file_id} evidence differs from the browser runtime")
        summaries[item.file_id] = {**summary, "serialRunSeconds": [round(value, 3) for value in timings[item.file_id]]}
    if write:
        _write_manifest(evidence, build, runs)
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
