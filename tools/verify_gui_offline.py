#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

from playwright.sync_api import TimeoutError as PlaywrightTimeoutError
from playwright.sync_api import sync_playwright


REPO_ROOT = Path(__file__).resolve().parent.parent
ALLOWED_HTTP_HOSTS = {"127.0.0.1", "localhost"}


@dataclass
class ScenarioResult:
    started_at_utc: str
    server_url: str | None = None
    pyodide_ready: bool = False
    circular_svg_rendered: bool = False
    circular_png_exported: bool = False
    circular_pdf_exported: bool = False
    linear_losat_generated: bool = False
    linear_losat_tsv_exported: bool = False
    external_http_requests: list[str] = field(default_factory=list)
    blocked_external_http_requests: list[str] = field(default_factory=list)
    request_failures: list[dict[str, Any]] = field(default_factory=list)
    downloaded_files: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    inputs: dict[str, str] = field(default_factory=dict)
    server_stdout: list[str] = field(default_factory=list)
    success: bool = False


def _now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _require_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing required {label}: {path}")


def _start_gui_server() -> tuple[subprocess.Popen[str], str, list[str]]:
    env = os.environ.copy()
    env["BROWSER"] = "true"
    env["PYTHONUNBUFFERED"] = "1"
    cmd = [sys.executable, "-u", "-m", "gbdraw.cli", "gui"]
    proc = subprocess.Popen(
        cmd,
        cwd=str(REPO_ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env,
    )

    collected: list[str] = []
    deadline = time.time() + 30
    while time.time() < deadline:
        line = proc.stdout.readline() if proc.stdout is not None else ""
        if not line:
            if proc.poll() is not None:
                break
            time.sleep(0.1)
            continue
        stripped = line.rstrip()
        collected.append(stripped)
        match = re.search(r"http://localhost:(\d+)", stripped)
        if match:
            return proc, f"http://127.0.0.1:{match.group(1)}", collected

    raise RuntimeError(
        "Failed to start gbdraw GUI local server.\n"
        + "\n".join(collected[-30:])
    )


def _stop_server(proc: subprocess.Popen[str]) -> None:
    if proc.poll() is not None:
        return
    proc.terminate()
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait(timeout=5)


def _wait_for_pyodide_ready(page, timeout_ms: int) -> None:
    page.wait_for_function(
        """
        () => {
          const startupError = document.body.innerText.includes('Startup Error:');
          if (startupError) return false;
          const button = Array.from(document.querySelectorAll('button'))
            .find((btn) => btn.innerText.includes('Generate Diagram'));
          return Boolean(button && !button.disabled);
        }
        """,
        timeout=timeout_ms,
    )


def _click_generate(page, timeout_ms: int) -> None:
    button = page.locator("button:has-text('Generate Diagram')").first
    button.click()
    page.wait_for_selector("button:has-text('PNG')", timeout=timeout_ms)


def _upload_genbank_files(page, *paths: Path) -> None:
    input_card = page.locator("div.card.border-l-4.border-l-blue-500").first
    file_inputs = input_card.locator("input[type='file']")
    count = file_inputs.count()
    if count < len(paths):
        raise RuntimeError(
            f"Not enough file inputs in Input Genomes card: needed {len(paths)}, found {count}"
        )
    for idx, path in enumerate(paths):
        file_inputs.nth(idx).set_input_files(str(path))


def _export_button_download(page, button_text: str, timeout_ms: int) -> str:
    button = page.locator(f"button:has-text('{button_text}')").first
    with page.expect_download(timeout=timeout_ms) as download_info:
        button.click()
    download = download_info.value
    return download.suggested_filename


def _run_browser_scenario(
    server_url: str,
    circular_gb: Path,
    linear_gb_1: Path,
    linear_gb_2: Path,
    timeout_ms: int,
    result: ScenarioResult,
) -> None:
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(headless=True)
        context = browser.new_context(accept_downloads=True)
        page = context.new_page()

        all_http_requests: list[str] = []
        blocked_external_http_requests: list[str] = []
        request_failures: list[dict[str, Any]] = []

        def route_handler(route, request) -> None:
            parsed = urlparse(request.url)
            if parsed.scheme in {"http", "https"}:
                all_http_requests.append(request.url)
                if parsed.hostname not in ALLOWED_HTTP_HOSTS:
                    blocked_external_http_requests.append(request.url)
                    route.abort()
                    return
            route.continue_()

        page.route("**/*", route_handler)

        def on_request_failed(request) -> None:
            request_failures.append(
                {
                    "url": request.url,
                    "failure": request.failure,
                }
            )

        page.on("requestfailed", on_request_failed)

        try:
            page.goto(server_url, wait_until="domcontentloaded", timeout=timeout_ms)
            _wait_for_pyodide_ready(page, timeout_ms=timeout_ms)
            result.pyodide_ready = True

            _upload_genbank_files(page, circular_gb)
            _click_generate(page, timeout_ms=timeout_ms)
            page.wait_for_selector("text=Result Preview", timeout=timeout_ms)
            page.wait_for_selector("main svg", timeout=timeout_ms)
            result.circular_svg_rendered = True

            png_name = _export_button_download(page, "PNG", timeout_ms=timeout_ms)
            pdf_name = _export_button_download(page, "PDF", timeout_ms=timeout_ms)
            result.downloaded_files.extend([png_name, pdf_name])
            result.circular_png_exported = True
            result.circular_pdf_exported = True

            page.locator("button:has-text('Linear')").first.click()
            page.locator("input[type='radio'][value='losat']").first.check()
            page.locator("button:has-text('Add Seq')").first.click()
            _upload_genbank_files(page, linear_gb_1, linear_gb_2)
            _click_generate(page, timeout_ms=timeout_ms)
            page.wait_for_selector("text=Result Preview", timeout=timeout_ms)
            losat_save = page.locator("button:has-text('Save LOSAT TSV')").first
            losat_save.wait_for(state="visible", timeout=timeout_ms)
            page.wait_for_function(
                """
                () => {
                  const btn = Array.from(document.querySelectorAll('button'))
                    .find((x) => x.innerText.includes('Save LOSAT TSV'));
                  return Boolean(btn && !btn.disabled);
                }
                """,
                timeout=timeout_ms,
            )
            result.linear_losat_generated = True

            losat_name = _export_button_download(page, "Save LOSAT TSV", timeout_ms=timeout_ms)
            result.downloaded_files.append(losat_name)
            result.linear_losat_tsv_exported = True
        finally:
            external_only = [
                url
                for url in all_http_requests
                if urlparse(url).hostname not in ALLOWED_HTTP_HOSTS
            ]
            result.external_http_requests = sorted(set(external_only))
            result.blocked_external_http_requests = sorted(set(blocked_external_http_requests))
            result.request_failures = request_failures
            context.close()
            browser.close()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run offline GUI verification for gbdraw with headless Chromium."
    )
    parser.add_argument(
        "--circular-gb",
        default=str(REPO_ROOT / "tests" / "test_inputs" / "NC_010162.gb"),
        help="GenBank file for circular generation test.",
    )
    parser.add_argument(
        "--linear-gb-1",
        default=str(REPO_ROOT / "tests" / "test_inputs" / "AP027078.gb"),
        help="First GenBank file for linear LOSAT test.",
    )
    parser.add_argument(
        "--linear-gb-2",
        default=str(REPO_ROOT / "tests" / "test_inputs" / "AP027131.gb"),
        help="Second GenBank file for linear LOSAT test.",
    )
    parser.add_argument(
        "--timeout-seconds",
        type=int,
        default=240,
        help="Timeout for major browser waits.",
    )
    parser.add_argument(
        "--report",
        default=str(REPO_ROOT / "gbdraw" / "web" / "offline_verification_report.json"),
        help="Path to write JSON verification report.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    circular_gb = Path(args.circular_gb).resolve()
    linear_gb_1 = Path(args.linear_gb_1).resolve()
    linear_gb_2 = Path(args.linear_gb_2).resolve()
    report_path = Path(args.report).resolve()
    timeout_ms = max(1, args.timeout_seconds) * 1000

    result = ScenarioResult(
        started_at_utc=_now_utc(),
        inputs={
            "circular_gb": str(circular_gb),
            "linear_gb_1": str(linear_gb_1),
            "linear_gb_2": str(linear_gb_2),
        },
    )

    for p, label in (
        (circular_gb, "circular input"),
        (linear_gb_1, "linear input #1"),
        (linear_gb_2, "linear input #2"),
    ):
        try:
            _require_file(p, label)
        except FileNotFoundError as exc:
            result.errors.append(str(exc))
            report_path.parent.mkdir(parents=True, exist_ok=True)
            report_path.write_text(json.dumps(asdict(result), indent=2), encoding="utf-8")
            print(f"[FAIL] {exc}")
            print(f"[INFO] Wrote report: {report_path}")
            return 1

    proc: subprocess.Popen[str] | None = None
    try:
        proc, server_url, server_stdout = _start_gui_server()
        result.server_url = server_url
        result.server_stdout = server_stdout
        _run_browser_scenario(
            server_url=server_url,
            circular_gb=circular_gb,
            linear_gb_1=linear_gb_1,
            linear_gb_2=linear_gb_2,
            timeout_ms=timeout_ms,
            result=result,
        )
    except PlaywrightTimeoutError as exc:
        result.errors.append(f"Playwright timeout: {exc}")
    except Exception as exc:  # noqa: BLE001
        result.errors.append(str(exc))
    finally:
        if proc is not None:
            _stop_server(proc)

    if result.external_http_requests:
        result.errors.append(
            "External HTTP(S) requests were observed: "
            + ", ".join(result.external_http_requests)
        )

    result.success = (
        result.pyodide_ready
        and result.circular_svg_rendered
        and result.circular_png_exported
        and result.circular_pdf_exported
        and result.linear_losat_generated
        and result.linear_losat_tsv_exported
        and not result.external_http_requests
        and not result.errors
    )

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(asdict(result), indent=2), encoding="utf-8")

    if result.success:
        print("[OK] Offline GUI verification succeeded")
        print(f"[INFO] Wrote report: {report_path}")
        return 0

    print("[FAIL] Offline GUI verification failed")
    for err in result.errors:
        print(f"  - {err}")
    print(f"[INFO] Wrote report: {report_path}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
