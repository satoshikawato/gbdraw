#!/usr/bin/env python3
"""Run the required LOSAT cache migration browser acceptance without skipping."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import threading
from contextlib import contextmanager
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Iterator


REPO_ROOT = Path(__file__).resolve().parents[1]
NODE_SPEC = REPO_ROOT / "tests" / "web" / "losat-cache-migration.playwright.spec.js"
FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "sessions"
FIXTURE_PATH = FIXTURE_DIR / "BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz"
EXPECTED_PATH = FIXTURE_DIR / "BGC0000708-BGC0000713.schema-v2.expected.json"


class AcceptanceChecks:
    def __init__(self) -> None:
        self.count = 0

    def require(self, condition: object, message: str) -> None:
        self.count += 1
        if not condition:
            raise AssertionError(message)


class _QuietHandler(SimpleHTTPRequestHandler):
    extensions_map = {
        **SimpleHTTPRequestHandler.extensions_map,
        ".css": "text/css; charset=utf-8",
        ".data": "application/octet-stream",
        ".js": "text/javascript; charset=utf-8",
        ".mjs": "text/javascript; charset=utf-8",
        ".wasm": "application/wasm",
        ".whl": "application/octet-stream",
        ".woff2": "font/woff2",
    }

    def log_message(self, format: str, *args: object) -> None:
        del format, args


@contextmanager
def _serve_repo() -> Iterator[str]:
    handler = lambda *args, **kwargs: _QuietHandler(  # noqa: E731
        *args, directory=str(REPO_ROOT), **kwargs
    )
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        yield f"http://127.0.0.1:{server.server_address[1]}"
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=10)


def _read_session(path: Path) -> dict[str, Any]:
    raw = path.read_bytes()
    if raw.startswith(b"\x1f\x8b"):
        raw = gzip.decompress(raw)
    value = json.loads(raw.decode("utf-8"))
    if not isinstance(value, dict):
        raise AssertionError(f"Session is not an object: {path}")
    return value


def _load_expected() -> dict[str, Any]:
    value = json.loads(EXPECTED_PATH.read_text(encoding="utf-8"))
    required = {
        "storedRawEntries",
        "totalPairs",
        "cacheHits",
        "cacheMisses",
        "uniqueJobs",
        "workerCalls",
    }
    missing = sorted(required.difference(value))
    if missing:
        raise AssertionError(f"Acceptance oracle is missing: {', '.join(missing)}")
    return value


def _node_playwright_cli() -> tuple[str, str] | None:
    node = shutil.which("node")
    if not node:
        return None
    resolved = subprocess.run(
        [node, "-p", "require.resolve('@playwright/test/cli')"],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if resolved.returncode != 0:
        return None
    cli = resolved.stdout.strip()
    return (node, cli) if cli else None


def _run_node_adapter(node: str, cli: str) -> int:
    env = os.environ.copy()
    env["GBDRAW_REPO"] = str(REPO_ROOT)
    command = [
        node,
        cli,
        "test",
        str(NODE_SPEC.relative_to(REPO_ROOT)),
        "--workers=1",
        "--reporter=line",
    ]
    print("Running LOSAT cache browser acceptance with Node Playwright.", flush=True)
    return subprocess.run(command, cwd=REPO_ROOT, env=env, check=False).returncode


def _import_session(page: Any, path: Path, checks: AcceptanceChecks) -> None:
    with page.expect_event("dialog", timeout=120_000) as dialog_info:
        page.locator('input[accept^=".json,"]').first.set_input_files(str(path))
    dialog = dialog_info.value
    checks.require(
        dialog.message == "Session loaded successfully!",
        f"Unexpected import dialog: {dialog.message!r}",
    )
    dialog.accept()
    page.wait_for_function(
        """() => {
          const app = window.__GBDRAW_APP__;
          return app?.mode === 'linear' && app?.linearSeqs?.length === 5;
        }""",
        timeout=120_000,
    )


def _save_session(page: Any, checks: AcceptanceChecks) -> Path:
    with page.expect_download(timeout=120_000) as download_info:
        page.evaluate("async () => window.__GBDRAW_APP__.saveSessionWithTitle()")
    path = download_info.value.path()
    checks.require(bool(path), "Session download has no local path.")
    return Path(path)


def _assert_legacy_preserved(
    session: dict[str, Any], expected: dict[str, Any], checks: AcceptanceChecks
) -> None:
    entries = session.get("losatCache", {}).get("entries", [])
    candidates = (
        session.get("legacyArtifacts", {})
        .get("proteinRawCandidates", {})
        .get("entries", [])
    )
    checks.require(entries == [], "Legacy protein entries leaked into the current cache.")
    checks.require(
        len(candidates) == expected["storedRawEntries"],
        "Load -> Save lost legacy protein candidates.",
    )
    checks.require(
        all(
            candidate.get("state") == "pending"
            and candidate.get("originalEntry", {}).get("schema") == 2
            and candidate.get("originalEntry", {}).get("program") == "blastp"
            for candidate in candidates
        ),
        "Saved legacy candidate envelope is not lossless schema 2.",
    )


def _assert_renamed_resources(session: dict[str, Any], checks: AcceptanceChecks) -> None:
    resources = [
        value
        for value in session.get("resources", {}).values()
        if value.get("kind") == "genbank"
    ]
    expected_names = {
        f"renamed-cache-input-{index}.gbk" for index in range(1, 6)
    }
    checks.require(len(resources) == 5, "Expected five embedded GenBank resources.")
    checks.require(
        all(resource.get("lastModified") == 0 for resource in resources),
        "A zero file mtime was replaced while saving.",
    )
    checks.require(
        all(any(name in resource.get("name", "") for name in expected_names) for resource in resources),
        "Renamed input filenames were not serialized.",
    )
    original_names = set(
        session.get("webFiles", {}).get("resourceOriginalNames", {}).values()
    )
    checks.require(
        expected_names.issubset(original_names),
        "Original-name metadata does not contain every renamed input.",
    )


def _assert_current_artifacts(
    session: dict[str, Any], expected: dict[str, Any], checks: AcceptanceChecks
) -> None:
    manifest = session.get("proteinIdentityManifest", {})
    entries = session.get("losatCache", {}).get("entries", [])
    checks.require(manifest.get("schema") == 1, "Protein identity manifest is not schema 1.")
    checks.require(
        len(entries) == expected["totalPairs"],
        "The migrated current cache does not contain every expected pair.",
    )
    checks.require(
        all(
            entry.get("schema") == 3
            and entry.get("kind") == "raw-losat"
            and entry.get("identityKind") == "protein"
            and entry.get("program") == "blastp"
            for entry in entries
        ),
        "The current protein cache contains a non-schema-3 entry.",
    )

    record_analyses = manifest.get("recordAnalyses", {})
    instances = manifest.get("recordInstances", {})
    owners: dict[str, tuple[str, str]] = {}
    for instance_key, instance in instances.items():
        checks.require(
            instance.get("recordAnalysisId") in record_analyses,
            f"Missing record analysis for {instance_key}.",
        )
        for feature_id, transport_id in instance.get("transportIds", {}).items():
            checks.require(feature_id.startswith("f_"), "Manifest feature ID is not stable.")
            checks.require(not transport_id.startswith("p_r_"), "Legacy p_r_ transport ID remains.")
            checks.require(not any(char.isspace() for char in transport_id), "Whitespace in transport ID.")
            checks.require(transport_id not in owners, "Transport ID is not globally unique.")
            owners[transport_id] = (instance_key, feature_id)

    resolved = 0
    for entry in entries:
        query_key = entry.get("queryRecordInstanceKey")
        subject_key = entry.get("subjectRecordInstanceKey")
        query = instances.get(query_key, {})
        subject = instances.get(subject_key, {})
        checks.require(
            query.get("bindingHash") == entry.get("queryBindingHash"),
            "Query binding hash does not resolve through the manifest.",
        )
        checks.require(
            subject.get("bindingHash") == entry.get("subjectBindingHash"),
            "Subject binding hash does not resolve through the manifest.",
        )
        query_ids = set(query.get("transportIds", {}).values())
        subject_ids = set(subject.get("transportIds", {}).values())
        for raw_line in entry.get("text", "").splitlines():
            if not raw_line.strip() or raw_line.lstrip().startswith("#"):
                continue
            query_id, subject_id, *_ = raw_line.split("\t")
            checks.require(query_id in query_ids, "QUERY does not resolve in its binding.")
            checks.require(subject_id in subject_ids, "SUBJECT does not resolve in its binding.")
            checks.require(owners.get(query_id, (None,))[0] == query_key, "QUERY owner is ambiguous.")
            checks.require(
                owners.get(subject_id, (None,))[0] == subject_key,
                "SUBJECT owner is ambiguous.",
            )
            resolved += 2
    checks.require(resolved > 0, "No QUERY/SUBJECT references were asserted.")


def _generate(page: Any) -> dict[str, Any]:
    return page.evaluate(
        """async () => {
          const app = window.__GBDRAW_APP__;
          const result = await app.runAnalysis();
          return {
            result,
            errorSummary: String(app.errorLog?.summary || ''),
            errorDetails: Array.isArray(app.errorLog?.details)
              ? app.errorLog.details.map((detail) => String(detail))
              : [],
            telemetry: app.lastRunInfo?.losatTelemetry ||
              window.__GBDRAW_LAST_LOSAT_TELEMETRY__ || null,
            executorCalls: Number(window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ || 0)
          };
        }"""
    )


def _assert_telemetry(
    run: dict[str, Any], expected: dict[str, Any], checks: AcceptanceChecks
) -> None:
    checks.require(run.get("result") == {"status": "ok"}, f"Generation failed: {run}")
    checks.require(not run.get("errorSummary"), f"Generation error: {run}")
    checks.require(not run.get("errorDetails"), f"Generation details: {run}")
    telemetry = run.get("telemetry") or {}
    for field in ("totalPairs", "cacheHits", "cacheMisses", "uniqueJobs", "workerCalls"):
        checks.require(
            telemetry.get(field) == expected[field],
            f"Telemetry {field}: expected {expected[field]!r}, got {telemetry.get(field)!r}",
        )
    checks.require(
        run.get("executorCalls") == expected["workerCalls"],
        "Counting LOSAT executor was called.",
    )


def _run_python_adapter() -> int:
    try:
        from playwright.sync_api import sync_playwright
    except ImportError as error:
        print(f"Neither Node nor Python Playwright is available: {error}", file=sys.stderr)
        return 2

    checks = AcceptanceChecks()
    expected = _load_expected()
    fixture_bytes = gzip.decompress(FIXTURE_PATH.read_bytes())
    fixture = json.loads(fixture_bytes.decode("utf-8"))
    checks.require(
        hashlib.sha256(fixture_bytes).hexdigest() == expected.get("sourceSha256"),
        "Legacy fixture SHA-256 does not match its acceptance oracle.",
    )
    checks.require(
        fixture.get("version") == expected.get("sessionVersion"),
        "Legacy fixture session version does not match its acceptance oracle.",
    )
    checks.require(
        fixture.get("renderRequest", {}).get("schema") == expected.get("renderRequestSchema"),
        "Legacy fixture renderRequest schema does not match its acceptance oracle.",
    )
    checks.require(
        len(fixture.get("losatCache", {}).get("entries", []))
        == expected["storedRawEntries"],
        "Legacy fixture raw-entry count does not match its acceptance oracle.",
    )
    print("Running LOSAT cache browser acceptance with Python Playwright.", flush=True)
    try:
        with _serve_repo() as base_url, sync_playwright() as playwright:
            browser = playwright.chromium.launch()
            try:
                context = browser.new_context(accept_downloads=True)
                page = context.new_page()
                page.set_default_timeout(120_000)
                page.add_init_script(
                    """window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ = 0;
                    window.__GBDRAW_LOSAT_EXECUTOR__ = async () => {
                      window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ += 1;
                      throw new Error('LOSAT executor must not run during verified cache migration.');
                    };"""
                )
                page.goto(f"{base_url}/gbdraw/web/index.html", wait_until="domcontentloaded")
                page.wait_for_function("() => window.__GBDRAW_APP__")
                _import_session(page, FIXTURE_PATH, checks)
                renamed = page.evaluate(
                    """async () => {
                      const app = window.__GBDRAW_APP__;
                      for (let index = 0; index < app.linearSeqs.length; index += 1) {
                        const source = app.linearSeqs[index].gb;
                        app.linearSeqs[index].gb = new File(
                          [await source.arrayBuffer()],
                          `renamed-cache-input-${index + 1}.gbk`,
                          { type: source.type || 'text/plain', lastModified: 0 }
                        );
                      }
                      return app.linearSeqs.map((record) => ({
                        name: record.gb.name,
                        lastModified: record.gb.lastModified
                      }));
                    }"""
                )
                checks.require(
                    renamed
                    == [
                        {"name": f"renamed-cache-input-{index}.gbk", "lastModified": 0}
                        for index in range(1, 6)
                    ],
                    "Browser File rename/zero-mtime setup failed.",
                )

                before_generate_path = _save_session(page, checks)
                before_generate = _read_session(before_generate_path)
                _assert_legacy_preserved(before_generate, expected, checks)
                _assert_renamed_resources(before_generate, checks)

                page.reload(wait_until="domcontentloaded")
                page.wait_for_function("() => window.__GBDRAW_APP__")
                _import_session(page, before_generate_path, checks)
                page.wait_for_function(
                    "() => window.__GBDRAW_APP__?.pyodideReady === true", timeout=240_000
                )
                _assert_telemetry(_generate(page), expected, checks)

                migrated_path = _save_session(page, checks)
                migrated = _read_session(migrated_path)
                _assert_current_artifacts(migrated, expected, checks)
                _assert_renamed_resources(migrated, checks)

                page.reload(wait_until="domcontentloaded")
                page.wait_for_function("() => window.__GBDRAW_APP__")
                _import_session(page, migrated_path, checks)
                page.wait_for_function(
                    "() => window.__GBDRAW_APP__?.pyodideReady === true", timeout=240_000
                )
                _assert_telemetry(_generate(page), expected, checks)
                context.close()
            finally:
                browser.close()
    except Exception as error:
        print(f"Python Playwright acceptance failed: {error}", file=sys.stderr)
        return 1

    if checks.count == 0:
        print("Python Playwright acceptance executed no assertions.", file=sys.stderr)
        return 1
    print(f"LOSAT cache browser acceptance passed ({checks.count} assertions).")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--python",
        action="store_true",
        help="Force the Python Playwright adapter (useful for validating the fallback).",
    )
    args = parser.parse_args(argv)

    for required in (NODE_SPEC, FIXTURE_PATH, EXPECTED_PATH):
        if not required.is_file():
            print(f"Required acceptance input is missing: {required}", file=sys.stderr)
            return 2

    node_adapter = None if args.python else _node_playwright_cli()
    if node_adapter:
        return _run_node_adapter(*node_adapter)
    return _run_python_adapter()


if __name__ == "__main__":
    raise SystemExit(main())
