#!/usr/bin/env python3
"""Run the required LOSAT cache migration browser acceptance without skipping."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import traceback
from contextlib import contextmanager
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Iterator


REPO_ROOT = Path(__file__).resolve().parents[1]
NODE_SPEC = REPO_ROOT / "tests" / "web" / "losat-cache-migration.playwright.spec.js"
FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "sessions"
FIXTURE_PATH = FIXTURE_DIR / "BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz"
EXPECTED_PATH = FIXTURE_DIR / "BGC0000708-BGC0000713.schema-v2.expected.json"
CURRENT_SESSION_VERSION = 38
CURRENT_RENDER_REQUEST_SCHEMA = 5
CURRENT_PROTEIN_RAW_SCHEMA = 4
CURRENT_PROTEIN_DERIVED_SCHEMA = 3


_LAYOUT_INSPECTION_SCRIPT = """async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const geometry = JSON.parse(JSON.stringify(state.trackSlotResolvedGeometry.value || null));
  const svg = state.svgContainer.value?.querySelector('svg');
  if (!geometry || !svg) {
    return { error: 'Resolved Linear geometry or the preview SVG is unavailable.' };
  }

  const epsilon = 1e-6;
  const records = (Array.isArray(geometry.records) ? geometry.records : [])
    .slice()
    .sort((left, right) => Number(left.recordIndex) - Number(right.recordIndex));
  const directGroups = Array.from(svg.children).filter(
    (element) => element.tagName?.toLowerCase() === 'g'
  );
  const comparisonGroups = directGroups
    .filter((element) => (
      element.hasAttribute('data-query-row') &&
      element.hasAttribute('data-subject-row')
    ))
    .sort((left, right) => (
      Number(left.getAttribute('data-query-row')) -
      Number(right.getAttribute('data-query-row'))
    ));
  const activeBoundaries = new Set(comparisonGroups.map((element) => (
    `${element.getAttribute('data-query-row')}:${element.getAttribute('data-subject-row')}`
  )));

  const xOverlaps = (left, right) => (
    Math.min(Number(left.xEndPx), Number(right.xEndPx)) -
    Math.max(Number(left.xStartPx), Number(right.xStartPx))
  ) > epsilon;
  const pairIsEligible = (left, right, boundaryActive) => (
    (left.kind === 'body' && right.kind === 'body') ||
    (boundaryActive && left.kind === 'comparison' && right.kind === 'comparison') ||
    (left.kind === 'definition' && right.kind === 'definition') ||
    (left.kind === 'definition' && right.kind === 'body') ||
    (left.kind === 'body' && right.kind === 'definition')
  );

  const collisionChecks = [];
  const axisGaps = [];
  for (let boundary = 0; boundary < records.length - 1; boundary += 1) {
    const nextRow = boundary + 1;
    const current = records[boundary];
    const following = records[nextRow];
    const boundaryActive = activeBoundaries.has(`${boundary}:${nextRow}`);
    axisGaps.push(Number(following?.axisYpx) - Number(current?.axisYpx));
    for (const currentBand of current?.collisionBands || []) {
      for (const nextBand of following?.collisionBands || []) {
        if (
          !xOverlaps(currentBand, nextBand) ||
          !pairIsEligible(currentBand, nextBand, boundaryActive)
        ) continue;
        collisionChecks.push({
          boundaryRow: boundary,
          kindPair: `${currentBand.kind}:${nextBand.kind}`,
          gapPx: Number(nextBand.absoluteTopPx) - Number(currentBand.absoluteBottomPx)
        });
      }
    }
  }

  const indexedSemanticGroups = (predicate) => directGroups
    .filter(predicate)
    .map((element) => {
      const rawIndex = element.getAttribute('data-gbdraw-record-index');
      const index = rawIndex === null || String(rawIndex).trim() === ''
        ? Number.NaN
        : Number(rawIndex);
      return Number.isInteger(index) && index >= 0 ? { element, index } : null;
    })
    .filter(Boolean)
    .sort((left, right) => left.index - right.index);
  const definitionRoles = new Set(['record-definition', 'record-definition-row']);
  const recordGroups = indexedSemanticGroups((element) => (
    element.hasAttribute('data-gbdraw-record-id') &&
    !definitionRoles.has(String(element.getAttribute('data-gbdraw-role') || '').trim())
  ));
  const definitionGroups = indexedSemanticGroups((element) => (
    String(element.getAttribute('data-gbdraw-role') || '').trim() === 'record-definition'
  ));
  const screenBox = (element) => {
    const box = element.getBoundingClientRect();
    return {
      id: String(element.id || ''),
      left: Number(box.left),
      right: Number(box.right),
      top: Number(box.top),
      bottom: Number(box.bottom)
    };
  };
  const overlap = (left, right) => ({
    overlapX: Math.min(left.right, right.right) - Math.max(left.left, right.left),
    overlapY: Math.min(left.bottom, right.bottom) - Math.max(left.top, right.top)
  });
  const domChecks = [];
  const addDomCheck = (label, leftElement, rightElement) => {
    if (!leftElement || !rightElement) {
      domChecks.push({ label, missing: true });
      return;
    }
    const left = screenBox(leftElement);
    const right = screenBox(rightElement);
    domChecks.push({ label, missing: false, ...overlap(left, right) });
  };
  for (let index = 0; index < records.length; index += 1) {
    addDomCheck(
      `definition-body:${index}`,
      definitionGroups[index]?.element,
      recordGroups[index]?.element
    );
    if (index >= records.length - 1) continue;
    addDomCheck(
      `body-body:${index}`,
      recordGroups[index]?.element,
      recordGroups[index + 1]?.element
    );
    addDomCheck(
      `definition-definition:${index}`,
      definitionGroups[index]?.element,
      definitionGroups[index + 1]?.element
    );
    addDomCheck(
      `comparison-upper-body:${index}`,
      comparisonGroups[index],
      recordGroups[index]?.element
    );
    addDomCheck(
      `comparison-lower-body:${index}`,
      comparisonGroups[index],
      recordGroups[index + 1]?.element
    );
  }

  return {
    error: '',
    schema: geometry.schema,
    mode: geometry.mode,
    recordCount: records.length,
    comparisonCount: comparisonGroups.length,
    recordGroupCount: recordGroups.length,
    definitionGroupCount: definitionGroups.length,
    axisGaps,
    collisionChecks,
    collisionDomains: Array.from(new Set(
      collisionChecks.map((check) => check.kindPair)
    )).sort(),
    domChecks,
    geometrySignature: JSON.stringify(geometry)
  };
}"""


_SVG_GEOMETRY_PARITY_SCRIPT = """async (exportedSvg) => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const preview = state.svgContainer.value?.querySelector('svg');
  const parsed = new DOMParser().parseFromString(String(exportedSvg || ''), 'image/svg+xml');
  const exported = parsed.documentElement;
  if (
    !preview ||
    !exported ||
    exported.tagName?.toLowerCase() !== 'svg' ||
    parsed.querySelector('parsererror')
  ) {
    return { error: 'The preview or exported SVG could not be parsed.' };
  }
  const geometryAttributes = [
    'viewBox', 'width', 'height', 'data-horizontal-viewbox', 'data-vertical-viewbox',
    'transform', 'd', 'points', 'x', 'y', 'x1', 'x2', 'y1', 'y2',
    'cx', 'cy', 'r', 'rx', 'ry', 'dx', 'dy', 'font-family', 'font-size',
    'font-weight', 'font-style', 'text-anchor', 'dominant-baseline',
    'stroke-width', 'display', 'href', 'xlink:href'
  ];
  const signature = (root) => [root, ...root.querySelectorAll('*')].map((element) => {
    const tag = element.tagName.toLowerCase();
    const attributes = geometryAttributes
      .filter((name) => element.hasAttribute(name))
      .map((name) => [name, element.getAttribute(name)]);
    const text = tag === 'text' || tag === 'tspan' ? element.textContent : '';
    return [tag, String(element.id || ''), attributes, text];
  });
  const previewSignature = signature(preview);
  const exportedSignature = signature(exported);
  const length = Math.max(previewSignature.length, exportedSignature.length);
  let mismatchIndex = -1;
  for (let index = 0; index < length; index += 1) {
    if (
      JSON.stringify(previewSignature[index] ?? null) !==
      JSON.stringify(exportedSignature[index] ?? null)
    ) {
      mismatchIndex = index;
      break;
    }
  }
  return {
    error: '',
    equal: mismatchIndex === -1,
    previewElementCount: previewSignature.length,
    exportedElementCount: exportedSignature.length,
    mismatchIndex,
    previewMismatch: mismatchIndex < 0 ? null : previewSignature[mismatchIndex],
    exportedMismatch: mismatchIndex < 0 ? null : exportedSignature[mismatchIndex]
  };
}"""


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
        "downloads",
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


def _download_svg(page: Any, checks: AcceptanceChecks) -> Path:
    with page.expect_download(timeout=120_000) as download_info:
        page.evaluate("() => window.__GBDRAW_APP__.downloadSVG()")
    path = download_info.value.path()
    checks.require(bool(path), "SVG download has no local path.")
    return Path(path)


def _download_losat_pair(
    page: Any,
    pair_index: int,
    checks: AcceptanceChecks,
) -> Path:
    with page.expect_download(timeout=120_000) as download_info:
        page.evaluate(
            "async (index) => window.__GBDRAW_APP__.downloadLosatPair(index)",
            pair_index,
        )
    path = download_info.value.path()
    checks.require(bool(path), "Pair LOSAT download has no local path.")
    return Path(path)


def _download_all_losat_pairs(
    page: Any,
    expected_count: int,
    checks: AcceptanceChecks,
) -> list[Path]:
    downloads: list[Any] = []
    complete = threading.Event()

    def capture(download: Any) -> None:
        downloads.append(download)
        if len(downloads) >= expected_count:
            complete.set()

    page.on("download", capture)
    try:
        page.evaluate("async () => window.__GBDRAW_APP__.downloadLosatCache()")
        complete.wait(timeout=30)
        checks.require(
            len(downloads) == expected_count,
            f"Bulk LOSAT export downloaded {len(downloads)} files; "
            f"expected {expected_count}.",
        )
        paths = [download.path() for download in downloads]
        checks.require(
            all(paths),
            "A bulk LOSAT download has no local path.",
        )
        return [Path(path) for path in paths]
    finally:
        page.remove_listener("download", capture)


def _losat_data_rows(text: str) -> list[list[str]]:
    return [
        line.split("\t")
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def _assert_hydrated_download(
    text: str,
    expected: dict[str, Any],
    checks: AcceptanceChecks,
) -> None:
    rows = _losat_data_rows(text)
    checks.require(bool(rows), "Hydrated LOSAT download contains no data rows.")
    checks.require(
        all(len(columns) == 12 for columns in rows),
        "Hydrated LOSAT download contains a non-12-column row.",
    )
    for pattern in expected["downloads"]["forbiddenPatterns"]:
        checks.require(
            re.search(pattern, text, flags=re.IGNORECASE) is None,
            f"Hydrated LOSAT download leaked internal ID pattern {pattern!r}.",
        )
    alias_pattern = re.compile(r"^[A-Za-z0-9._%~-]+$")
    checks.require(
        all(
            alias_pattern.fullmatch(columns[0])
            and alias_pattern.fullmatch(columns[1])
            for columns in rows
        ),
        "Hydrated LOSAT QUERY/SUBJECT is not a normal export alias.",
    )


def _assert_hydrated_downloads(
    page: Any,
    expected: dict[str, Any],
    checks: AcceptanceChecks,
) -> None:
    pair_path = _download_losat_pair(
        page,
        int(expected["downloads"]["pairIndex"]),
        checks,
    )
    pair_text = pair_path.read_text(encoding="utf-8")
    _assert_hydrated_download(pair_text, expected, checks)
    checks.require(
        _losat_data_rows(pair_text)[0][:2]
        == expected["downloads"]["firstDataIds"],
        "Pair download did not expose the expected ordinary aliases.",
    )

    bulk_paths = _download_all_losat_pairs(
        page,
        int(expected["downloads"]["bulkEntryCount"]),
        checks,
    )
    for path in bulk_paths:
        _assert_hydrated_download(path.read_text(encoding="utf-8"), expected, checks)


def _assert_hydrated_utf8_prompt_boundary(
    page: Any,
    expected: dict[str, Any],
    checks: AcceptanceChecks,
) -> None:
    probe = page.evaluate(
        """async () => {
          const {
            LOSAT_EXPORT_CONFIRM_THRESHOLD_BYTES,
            confirmHydratedLosatExport,
            totalHydratedLosatExportBytes
          } = await import('/gbdraw/web/js/app/run-analysis.js');
          const hydrated = { text: '界界' };
          const utf8Bytes = totalHydratedLosatExportBytes([hydrated]);
          let prompts = 0;
          const accepted = confirmHydratedLosatExport(
            utf8Bytes,
            () => {
              prompts += 1;
              return true;
            },
            5
          );
          return {
            threshold: LOSAT_EXPORT_CONFIRM_THRESHOLD_BYTES,
            codeUnits: hydrated.text.length,
            utf8Bytes,
            prompts,
            accepted
          };
        }"""
    )
    checks.require(
        probe
        == {
            "threshold": expected["downloads"]["hydratedPromptThresholdBytes"],
            "codeUnits": 2,
            "utf8Bytes": 6,
            "prompts": 1,
            "accepted": True,
        },
        f"Hydrated UTF-8 prompt boundary is incorrect: {probe}",
    )


def _install_user_uploaded_tsv(
    page: Any,
    expected: dict[str, Any],
) -> list[int]:
    return page.evaluate(
        """async ({ filename, text }) => {
          const bytes = new TextEncoder().encode(text);
          const app = window.__GBDRAW_APP__;
          app.blastSource = 'upload';
          app.losatProgram = 'blastn';
          app.linearSeqs[0].blast = new File(
            [bytes],
            filename,
            { type: 'text/tab-separated-values', lastModified: 0 }
          );
          return Array.from(
            new Uint8Array(await app.linearSeqs[0].blast.arrayBuffer())
          );
        }""",
        {
            "filename": expected["downloads"]["userUploadedFilename"],
            "text": expected["downloads"]["userUploadedTsv"],
        },
    )


def _read_first_uploaded_tsv(page: Any) -> dict[str, Any]:
    return page.evaluate(
        """async () => {
          const file = window.__GBDRAW_APP__.linearSeqs[0].blast;
          return {
            name: file?.name || '',
            bytes: file
              ? Array.from(new Uint8Array(await file.arrayBuffer()))
              : []
          };
        }"""
    )


def _assert_current_session_boundary(
    session: dict[str, Any],
    checks: AcceptanceChecks,
    *,
    require_derived: bool,
) -> None:
    raw_entries = session.get("losatCache", {}).get("entries", [])
    derived_entries = session.get("losatDerivedCache", {}).get("entries", [])
    checks.require(
        session.get("format") == "gbdraw-session",
        "Saved payload is not a gbdraw session.",
    )
    checks.require(
        session.get("version") == CURRENT_SESSION_VERSION,
        f"Saved session is not version {CURRENT_SESSION_VERSION}.",
    )
    checks.require(
        session.get("renderRequest", {}).get("schema")
        == CURRENT_RENDER_REQUEST_SCHEMA,
        f"Saved renderRequest is not schema {CURRENT_RENDER_REQUEST_SCHEMA}.",
    )
    checks.require(
        all(
            entry.get("identityKind") != "protein"
            or entry.get("schema") == CURRENT_PROTEIN_RAW_SCHEMA
            for entry in raw_entries
        ),
        "A legacy protein entry remains in the current raw cache.",
    )
    checks.require(
        all(
            entry.get("schema") == CURRENT_PROTEIN_DERIVED_SCHEMA
            and entry.get("idEncoding") == "runtime-handle-v1"
            for entry in derived_entries
        ),
        "A legacy entry remains in the current derived cache.",
    )
    if require_derived:
        checks.require(
            bool(derived_entries),
            "Generation did not persist a current derived cache entry.",
        )


def _assert_legacy_preserved(
    session: dict[str, Any],
    source_session: dict[str, Any],
    expected: dict[str, Any],
    checks: AcceptanceChecks,
) -> None:
    _assert_current_session_boundary(session, checks, require_derived=False)
    entries = session.get("losatCache", {}).get("entries", [])
    legacy_artifacts = session.get("legacyArtifacts", {})
    raw_envelope = legacy_artifacts.get("proteinRawCandidates", {})
    candidates = raw_envelope.get("entries", [])
    derived_entries = session.get("losatDerivedCache", {}).get("entries", [])
    derived_envelope = legacy_artifacts.get("proteinDerivedEvidence", {})
    checks.require(entries == [], "Legacy protein entries leaked into the current cache.")
    checks.require(
        derived_entries == [],
        "Legacy derived entries leaked into the current derived cache.",
    )
    checks.require(
        raw_envelope.get("schema") == 1,
        "Legacy protein candidate envelope is not schema 1.",
    )
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
    checks.require(
        [candidate.get("originalEntry") for candidate in candidates]
        == source_session.get("losatCache", {}).get("entries", []),
        "Load -> Save changed a legacy raw candidate.",
    )
    checks.require(
        derived_envelope.get("schema") == 1
        and derived_envelope.get("entries", [])
        == source_session.get("losatDerivedCache", {}).get("entries", []),
        "Load -> Save changed the quarantined legacy derived evidence.",
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
    session: dict[str, Any],
    source_session: dict[str, Any],
    expected: dict[str, Any],
    checks: AcceptanceChecks,
) -> None:
    _assert_current_session_boundary(session, checks, require_derived=True)
    manifest = session.get("proteinIdentityManifest", {})
    entries = session.get("losatCache", {}).get("entries", [])
    derived_entries = session.get("losatDerivedCache", {}).get("entries", [])
    legacy_artifacts = session.get("legacyArtifacts", {})
    checks.require(manifest.get("schema") == 2, "Protein identity manifest is not schema 2.")
    checks.require(
        len(entries) == expected["totalPairs"],
        "The migrated current cache does not contain every expected pair.",
    )
    checks.require(
        all(
            entry.get("schema") == 4
            and entry.get("kind") == "raw-losat"
            and entry.get("identityKind") == "protein"
            and entry.get("idEncoding") == "runtime-handle-v1"
            and entry.get("program") == "blastp"
            for entry in entries
        ),
        "The current protein cache contains a non-schema-4 entry.",
    )
    checks.require(
        all(
            entry.get("schema") == CURRENT_PROTEIN_DERIVED_SCHEMA
            and entry.get("kind") == "derived-losatp-payload"
            and entry.get("idEncoding") == "runtime-handle-v1"
            for entry in derived_entries
        ),
        "The current derived cache contains a non-schema-3 protein payload.",
    )
    legacy_candidates = (
        legacy_artifacts.get("proteinRawCandidates", {}).get("entries", [])
    )
    legacy_derived = (
        legacy_artifacts.get("proteinDerivedEvidence", {}).get("entries", [])
    )
    checks.require(
        len(legacy_candidates)
        == expected["storedRawEntries"] - expected["totalPairs"],
        "The saved session did not remove exactly the promoted legacy candidates.",
    )
    checks.require(
        all(
            candidate.get("state") == "pending"
            and candidate.get("originalEntry")
            in source_session.get("losatCache", {}).get("entries", [])
            for candidate in legacy_candidates
        ),
        "A promoted or mutated legacy candidate remains in the saved envelope.",
    )
    checks.require(
        legacy_derived
        == source_session.get("losatDerivedCache", {}).get("entries", []),
        "Pending-candidate derived evidence was not preserved losslessly.",
    )

    record_analyses = manifest.get("recordAnalyses", {})
    instances = manifest.get("recordInstances", {})
    owners: dict[str, tuple[str, str]] = {}
    for instance_key, instance in instances.items():
        checks.require(
            instance.get("recordAnalysisId") in record_analyses,
            f"Missing record analysis for {instance_key}.",
        )
        checks.require(
            bool(instance.get("runtimeBindingHash"))
            and bool(instance.get("displayBindingHash")),
            f"Missing compact binding hashes for {instance_key}.",
        )
        checks.require(
            set(instance.get("runtimeIds", {}))
            == set(instance.get("featureMetadata", {})),
            f"Runtime and display maps differ for {instance_key}.",
        )
        for feature_id, runtime_handle in instance.get("runtimeIds", {}).items():
            checks.require(feature_id.startswith("f_"), "Manifest feature ID is not stable.")
            checks.require(
                re.fullmatch(r"h_[a-z2-7]{26}", runtime_handle) is not None,
                "Runtime handle has the wrong encoding.",
            )
            checks.require(runtime_handle not in owners, "Runtime handle is not globally unique.")
            owners[runtime_handle] = (instance_key, feature_id)

    resolved = 0
    for entry in entries:
        query_key = entry.get("queryRecordInstanceKey")
        subject_key = entry.get("subjectRecordInstanceKey")
        query = instances.get(query_key, {})
        subject = instances.get(subject_key, {})
        checks.require(
            query.get("runtimeBindingHash") == entry.get("queryRuntimeBindingHash"),
            "Query binding hash does not resolve through the manifest.",
        )
        checks.require(
            subject.get("runtimeBindingHash") == entry.get("subjectRuntimeBindingHash"),
            "Subject binding hash does not resolve through the manifest.",
        )
        query_ids = set(query.get("runtimeIds", {}).values())
        subject_ids = set(subject.get("runtimeIds", {}).values())
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

    reference_keys = {
        "query_protein_id",
        "subject_protein_id",
        "queryProteinId",
        "subjectProteinId",
        "proteinId",
    }
    reference_list_keys = {"proteinIds", "sharedProteinIds"}
    derived_references: list[str] = []

    def collect_references(value: object) -> None:
        if isinstance(value, list):
            for item in value:
                collect_references(item)
            return
        if not isinstance(value, dict):
            return
        for key, item in value.items():
            if key in reference_keys and isinstance(item, str):
                derived_references.append(item)
            elif key in reference_list_keys and isinstance(item, list):
                derived_references.extend(
                    reference for reference in item if isinstance(reference, str)
                )
            collect_references(item)

    collect_references(derived_entries)
    checks.require(
        bool(derived_references),
        "No derived protein references were asserted.",
    )
    unresolved = sorted(set(derived_references).difference(owners))
    checks.require(
        unresolved == [],
        f"Derived protein references do not resolve through the manifest: {unresolved[:3]}",
    )
    checks.require(
        "p_r_" not in json.dumps(derived_entries, ensure_ascii=False),
        "A legacy metadata-derived protein ID remains in the derived cache.",
    )
    checks.require(
        re.search(
            r"@.+\|.+~f_[0-9a-f]{64}",
            json.dumps(derived_entries, ensure_ascii=False),
        )
        is None,
        "A version-35 readable transport ID remains in the derived cache.",
    )


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


def _migration_ui_snapshot(page: Any) -> dict[str, Any]:
    return page.evaluate(
        """async () => {
          const { state } = await import('/gbdraw/web/js/state.js');
          return JSON.parse(JSON.stringify({
            orthogroups: state.orthogroups.value,
            selectedOrthogroupAlignmentFeature:
              state.selectedOrthogroupAlignmentFeature.value,
            selectedOrthogroupId: state.selectedOrthogroupId.value,
            extractedFeatures: state.extractedFeatures.value,
            biologicalFeatures: state.biologicalFeatures.value
          }));
        }"""
    )


def _cancel_during_render(page: Any) -> dict[str, Any]:
    return page.evaluate(
        """async () => {
          const app = window.__GBDRAW_APP__;
          const { state } = await import('/gbdraw/web/js/state.js');
          const before = {
            proteinIdentityManifest: state.proteinIdentityManifest.value,
            legacyProteinRawCandidates: state.legacyProteinRawCandidates.value,
            legacyProteinDerivedEvidence:
              state.legacyProteinDerivedEvidence.value,
            losatCache: Array.from(state.losatCache.value.entries()),
            losatDerivedCache: Array.from(state.losatDerivedCache.value.entries()),
            losatCacheInfo: state.losatCacheInfo.value
          };
          let sawRendering = false;
          let cancelInvoked = false;
          const cancelPoll = setInterval(() => {
            if (
              sawRendering ||
              String(app.processingStatus || '') !== 'Rendering SVG...'
            ) return;
            sawRendering = true;
            app.cancelGeneration();
            cancelInvoked = true;
          }, 0);
          try {
            const result = await app.runAnalysis();
            const sameMapEntries = (entries, current) => (
              entries.length === current.size &&
              entries.every(([key, value]) => current.get(key) === value)
            );
            return {
              result,
              sawRendering,
              cancelInvoked,
              errorSummary: String(app.errorLog?.summary || ''),
              executorCalls: Number(window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ || 0),
              authorityRestored: (
                state.proteinIdentityManifest.value === before.proteinIdentityManifest &&
                state.legacyProteinRawCandidates.value ===
                  before.legacyProteinRawCandidates &&
                state.legacyProteinDerivedEvidence.value ===
                  before.legacyProteinDerivedEvidence &&
                sameMapEntries(before.losatCache, state.losatCache.value) &&
                sameMapEntries(
                  before.losatDerivedCache,
                  state.losatDerivedCache.value
                ) &&
                state.losatCacheInfo.value === before.losatCacheInfo
              )
            };
          } finally {
            clearInterval(cancelPoll);
          }
        }"""
    )


def _fail_renderer_after_migration(page: Any) -> dict[str, Any]:
    return page.evaluate(
        """async () => {
          const app = window.__GBDRAW_APP__;
          const { state } = await import('/gbdraw/web/js/state.js');
          const before = {
            proteinIdentityManifest: state.proteinIdentityManifest.value,
            legacyProteinRawCandidates: state.legacyProteinRawCandidates.value,
            legacyProteinDerivedEvidence:
              state.legacyProteinDerivedEvidence.value,
            losatCache: Array.from(state.losatCache.value.entries()),
            losatDerivedCache: Array.from(state.losatDerivedCache.value.entries()),
            losatCacheInfo: state.losatCacheInfo.value
          };
          const previousLegendFontSize = app.adv.legend_font_size;
          app.adv.legend_font_size = 'not-a-number';
          try {
            const result = await app.runAnalysis();
            const sameMapEntries = (entries, current) => (
              entries.length === current.size &&
              entries.every(([key, value]) => current.get(key) === value)
            );
            return {
              result,
              errorSummary: String(app.errorLog?.summary || ''),
              executorCalls: Number(window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ || 0),
              authorityRestored: (
                state.proteinIdentityManifest.value === before.proteinIdentityManifest &&
                state.legacyProteinRawCandidates.value ===
                  before.legacyProteinRawCandidates &&
                state.legacyProteinDerivedEvidence.value ===
                  before.legacyProteinDerivedEvidence &&
                sameMapEntries(before.losatCache, state.losatCache.value) &&
                sameMapEntries(
                  before.losatDerivedCache,
                  state.losatDerivedCache.value
                ) &&
                state.losatCacheInfo.value === before.losatCacheInfo
              )
            };
          } finally {
            app.adv.legend_font_size = previousLegendFontSize;
          }
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


def _inspect_layout(page: Any) -> dict[str, Any]:
    page.wait_for_function(
        """() => {
          const app = window.__GBDRAW_APP__;
          return app?.svgContainer?.querySelector('svg')
            ?.querySelectorAll('g[data-query-row][data-subject-row]').length === 4;
        }""",
        timeout=120_000,
    )
    return page.evaluate(_LAYOUT_INSPECTION_SCRIPT)


def _assert_layout(
    layout: dict[str, Any],
    checks: AcceptanceChecks,
) -> None:
    checks.require(not layout.get("error"), f"Layout inspection failed: {layout}")
    # Browser run metadata publishes collision bands in the schema-1 wrapper.
    # The canvas-private schema-2 axisGapConstraints are covered by Python unit tests.
    checks.require(
        layout.get("schema") == 1 and layout.get("mode") == "linear",
        "Published Linear run geometry is not schema 1: "
        f"schema={layout.get('schema')!r}, mode={layout.get('mode')!r}.",
    )
    checks.require(layout.get("recordCount") == 5, "Expected five resolved BGC records.")
    checks.require(
        layout.get("comparisonCount") == 4,
        "Expected four rendered BGC comparison corridors.",
    )
    checks.require(
        layout.get("recordGroupCount") == 5,
        "Expected five rendered BGC record bodies.",
    )
    checks.require(
        layout.get("definitionGroupCount") == 5,
        "Expected five rendered BGC definition groups.",
    )

    axis_gaps = layout.get("axisGaps", [])
    checks.require(
        len(axis_gaps) == 4 and all(float(gap) > 0 for gap in axis_gaps),
        f"Resolved BGC record axes are not strictly separated: {axis_gaps}",
    )

    collision_checks = layout.get("collisionChecks", [])
    checks.require(
        bool(collision_checks),
        "No eligible collision-domain pairs were inspected.",
    )
    collisions = [
        check
        for check in collision_checks
        if float(check.get("gapPx", float("-inf"))) < -1e-6
    ]
    checks.require(
        collisions == [],
        f"Resolved collision bands overlap: {collisions[:3]}",
    )
    domains = set(layout.get("collisionDomains", []))
    checks.require(
        {"body:body", "comparison:comparison", "definition:definition"}
        <= domains,
        f"BGC collision coverage is incomplete: {sorted(domains)}",
    )

    dom_checks = layout.get("domChecks", [])
    checks.require(bool(dom_checks), "No rendered SVG bbox pairs were inspected.")
    dom_collisions = [
        check
        for check in dom_checks
        if check.get("missing")
        or (
            float(check.get("overlapX", 0.0)) > 1e-6
            and float(check.get("overlapY", 0.0)) > 1e-6
        )
    ]
    checks.require(
        dom_collisions == [],
        f"Rendered definition/body/comparison bboxes overlap: {dom_collisions[:3]}",
    )


def _assert_svg_geometry_parity(
    page: Any,
    exported_path: Path,
    checks: AcceptanceChecks,
) -> None:
    parity = page.evaluate(
        _SVG_GEOMETRY_PARITY_SCRIPT,
        exported_path.read_text(encoding="utf-8"),
    )
    checks.require(not parity.get("error"), f"SVG parity inspection failed: {parity}")
    checks.require(
        parity.get("previewElementCount", 0) > 0,
        "Preview SVG geometry signature is empty.",
    )
    checks.require(
        parity.get("previewElementCount") == parity.get("exportedElementCount"),
        f"Preview/export SVG element counts differ: {parity}",
    )
    checks.require(
        parity.get("equal") is True,
        f"Preview/export SVG geometry differs: {parity}",
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
                page.on(
                    "pageerror",
                    lambda error: print(
                        f"Browser page error: {error}",
                        file=sys.stderr,
                    ),
                )
                page.on(
                    "console",
                    lambda message: (
                        print(
                            f"Browser console error: {message.text}",
                            file=sys.stderr,
                        )
                        if message.type == "error"
                        else None
                    ),
                )
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
                _assert_legacy_preserved(
                    before_generate,
                    fixture,
                    expected,
                    checks,
                )
                _assert_renamed_resources(before_generate, checks)

                page.reload(wait_until="domcontentloaded")
                page.wait_for_function("() => window.__GBDRAW_APP__")
                _import_session(page, before_generate_path, checks)
                page.wait_for_function(
                    "() => window.__GBDRAW_APP__?.pyodideReady === true", timeout=240_000
                )
                legacy_ui_before_cancel = _migration_ui_snapshot(page)
                canceled_legacy_run = _cancel_during_render(page)
                checks.require(
                    canceled_legacy_run.get("sawRendering")
                    and canceled_legacy_run.get("cancelInvoked")
                    and canceled_legacy_run.get("authorityRestored")
                    and canceled_legacy_run.get("result") == {"status": "canceled"},
                    "Legacy render cancellation did not roll migration back: "
                    f"{canceled_legacy_run}",
                )
                checks.require(
                    not canceled_legacy_run.get("errorSummary"),
                    f"Legacy render cancellation surfaced an error: {canceled_legacy_run}",
                )
                checks.require(
                    _migration_ui_snapshot(page) == legacy_ui_before_cancel,
                    "Legacy UI references changed after the canceled migration.",
                )
                legacy_ui_before_render_error = _migration_ui_snapshot(page)
                failed_legacy_run = _fail_renderer_after_migration(page)
                checks.require(
                    failed_legacy_run.get("result") == {"status": "error"}
                    and failed_legacy_run.get("authorityRestored")
                    and bool(failed_legacy_run.get("errorSummary"))
                    and failed_legacy_run.get("executorCalls") == 0,
                    "Legacy renderer failure did not roll migration back: "
                    f"{failed_legacy_run}",
                )
                checks.require(
                    _migration_ui_snapshot(page) == legacy_ui_before_render_error,
                    "Legacy UI references changed after the failed render.",
                )
                _assert_telemetry(_generate(page), expected, checks)
                first_layout = _inspect_layout(page)
                _assert_layout(first_layout, checks)
                _assert_svg_geometry_parity(
                    page,
                    _download_svg(page, checks),
                    checks,
                )
                _assert_hydrated_downloads(page, expected, checks)
                _assert_hydrated_utf8_prompt_boundary(page, expected, checks)

                migrated_path = _save_session(page, checks)
                migrated = _read_session(migrated_path)
                _assert_current_artifacts(migrated, fixture, expected, checks)
                _assert_renamed_resources(migrated, checks)

                page.reload(wait_until="domcontentloaded")
                page.wait_for_function("() => window.__GBDRAW_APP__")
                _import_session(page, migrated_path, checks)
                page.wait_for_function(
                    "() => window.__GBDRAW_APP__?.pyodideReady === true", timeout=240_000
                )
                _assert_telemetry(_generate(page), expected, checks)
                second_layout = _inspect_layout(page)
                _assert_layout(second_layout, checks)
                checks.require(
                    second_layout.get("geometrySignature")
                    == first_layout.get("geometrySignature"),
                    "Resolved BGC geometry changed after save/reload/regeneration.",
                )

                expected_upload_bytes = _install_user_uploaded_tsv(page, expected)
                uploaded_session_path = _save_session(page, checks)
                page.reload(wait_until="domcontentloaded")
                page.wait_for_function("() => window.__GBDRAW_APP__")
                _import_session(page, uploaded_session_path, checks)
                restored_upload = _read_first_uploaded_tsv(page)
                checks.require(
                    restored_upload
                    == {
                        "name": expected["downloads"]["userUploadedFilename"],
                        "bytes": expected_upload_bytes,
                    },
                    "User-uploaded TSV bytes changed across save/reload.",
                )
                context.close()
            finally:
                browser.close()
    except Exception as error:
        print(
            f"Python Playwright acceptance failed: {error}\n{traceback.format_exc()}",
            file=sys.stderr,
        )
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
