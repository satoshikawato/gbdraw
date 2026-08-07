"""Visible-UI capture for all-record Hepatoplasmataceae Collinear blocks."""

from __future__ import annotations

import gzip
import json
import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from Bio import SeqIO
from Bio.Seq import Seq
from playwright.sync_api import BrowserType, Page, expect

from assertions.svg_semantics import (
    assert_static_svg_safety,
    inspect_gui_bgc_losatp_svg,
    inspect_svg_file,
)
from config import ACTION_TIMEOUT_MS, GUI_HEPATOPLASMATACEAE_FIXTURES
from flows.web_capture import (
    assert_fixture_identity,
    capture_screenshot,
    fit_complete_linear_preview,
    generate_and_inspect,
    open_browser_capture,
    set_feature_search_visible,
    wait_for_worker,
)


RECORD_IDS = tuple(fixture[3] for fixture in GUI_HEPATOPLASMATACEAE_FIXTURES)
RECORD_LENGTH_TEXT = {
    record_id: f"{length:,} bp"
    for _, _, _, record_id, length, _ in GUI_HEPATOPLASMATACEAE_FIXTURES
}
ADJACENT_PAIRS = set(zip(RECORD_IDS[:-1], RECORD_IDS[1:], strict=True))
EXPECTED_FEATURE_ELEMENTS = 2_994
ALL_VS_ALL_SESSION = (
    Path(__file__).resolve().parents[3]
    / "gbdraw"
    / "web"
    / "gallery"
    / "sessions"
    / "hepatoplasmataceae_orthogroup.gbdraw-session.json.gz"
)
ALL_VS_ALL_SESSION_SIZE = 14_973_140
ALL_VS_ALL_SESSION_SHA256 = (
    "d9aeaa9ab70cede297f6ae2eaa51a8ed0a6ffc40216b70cb97570d2518ead500"
)
GALLERY_COLLINEAR_SESSION = (
    Path(__file__).resolve().parents[3]
    / "gbdraw"
    / "web"
    / "gallery"
    / "sessions"
    / "hepatoplasmataceae_collinear.gbdraw-session.json.gz"
)
GALLERY_COLLINEAR_SESSION_SIZE = 11_398_634
GALLERY_COLLINEAR_SESSION_SHA256 = (
    "4a054b6239b83b98f22105dd5bf409592a94b8a75d8fa12710ab118ddbaf825c"
)


@dataclass(frozen=True)
class HepatoplasmataceaeCollinearResult:
    screenshot_bytes: dict[str, int]
    final_svg_semantics: dict[str, Any]
    download: dict[str, Any]
    fasta_download: dict[str, Any]
    popup: dict[str, Any]
    block_count: int
    source_records: tuple[dict[str, Any], ...]


def _validate_fixtures() -> tuple[dict[str, Any], ...]:
    reports = []
    for path, size, digest, record_id, length, organism in (
        GUI_HEPATOPLASMATACEAE_FIXTURES
    ):
        assert_fixture_identity(path, expected_size=size, expected_sha256=digest)
        records = list(SeqIO.parse(path, "genbank"))
        if len(records) != 1:
            raise AssertionError(f"Expected one complete record in {path.name}")
        record = records[0]
        if (
            record.id != record_id
            or len(record) != length
            or record.annotations.get("organism") != organism
            or str(record.annotations.get("topology", "")).lower() != "circular"
            or "complete" not in record.description.lower()
        ):
            raise AssertionError(f"Unexpected Hepatoplasmataceae record: {path}")
        cds_count = sum(feature.type == "CDS" for feature in record.features)
        if cds_count < 500:
            raise AssertionError(f"{record_id} has too few annotated CDS features")
        reports.append(
            {
                "path": path,
                "record_id": record.id,
                "length": len(record),
                "organism": organism,
                "cds_count": cds_count,
            }
        )
    return tuple(reports)


def _load_all_vs_all_session(page: Page) -> frozenset[str]:
    assert_fixture_identity(
        ALL_VS_ALL_SESSION,
        expected_size=ALL_VS_ALL_SESSION_SIZE,
        expected_sha256=ALL_VS_ALL_SESSION_SHA256,
    )
    with gzip.open(ALL_VS_ALL_SESSION, "rt", encoding="utf-8") as handle:
        payload = json.load(handle)
    entries = payload.get("losatCache", {}).get("entries", [])
    session_keys = frozenset(str(entry.get("key", "")) for entry in entries)
    if len(entries) != 25 or len(session_keys) != 25 or "" in session_keys:
        raise AssertionError("Starting session is not the fixed 5x5 LOSATP cache")
    page.once("dialog", lambda dialog: dialog.accept())
    with page.expect_event("dialog", timeout=ACTION_TIMEOUT_MS) as dialog_info:
        with page.expect_file_chooser(timeout=ACTION_TIMEOUT_MS) as chooser_info:
            page.get_by_role("button", name="Load Session", exact=True).click()
        chooser_info.value.set_files(ALL_VS_ALL_SESSION)
    if dialog_info.value.message != "Session loaded successfully!":
        raise AssertionError(f"Session load failed: {dialog_info.value.message}")
    expect(page.get_by_role("button", name="Linear", exact=True)).to_have_attribute(
        "aria-pressed", "true"
    )
    expect(page.get_by_role("group", name="GenBank File selection", exact=True)).to_have_count(5)
    return session_keys


def _load_gallery_collinear_session(page: Page) -> frozenset[str]:
    assert_fixture_identity(
        GALLERY_COLLINEAR_SESSION,
        expected_size=GALLERY_COLLINEAR_SESSION_SIZE,
        expected_sha256=GALLERY_COLLINEAR_SESSION_SHA256,
    )
    with gzip.open(GALLERY_COLLINEAR_SESSION, "rt", encoding="utf-8") as handle:
        payload = json.load(handle)
    entries = payload.get("losatCache", {}).get("entries", [])
    session_keys = frozenset(str(entry.get("key", "")) for entry in entries)
    blastp = payload.get("config", {}).get("losat", {}).get("blastp", {})
    if (
        len(entries) != 13
        or len(session_keys) != 13
        or "" in session_keys
        or blastp.get("mode") != "collinear"
        or blastp.get("collinearSearchScope") != "adjacent"
    ):
        raise AssertionError("Gallery session is not the fixed adjacent Collinear project")
    page.once("dialog", lambda dialog: dialog.accept())
    with page.expect_event("dialog", timeout=ACTION_TIMEOUT_MS) as dialog_info:
        with page.expect_file_chooser(timeout=ACTION_TIMEOUT_MS) as chooser_info:
            page.get_by_role("button", name="Load Session", exact=True).click()
        chooser_info.value.set_files(GALLERY_COLLINEAR_SESSION)
    if dialog_info.value.message != "Session loaded successfully!":
        raise AssertionError(f"Session load failed: {dialog_info.value.message}")
    expect(page.get_by_role("button", name="Linear", exact=True)).to_have_attribute(
        "aria-pressed", "true"
    )
    expect(
        page.get_by_role("group", name="GenBank File selection", exact=True)
    ).to_have_count(5)
    return session_keys


def _set_presentation(page: Page, *, title: str) -> None:
    page.get_by_label("Track Layout", exact=True).select_option("middle")
    for checkbox_name in (
        "Separate Strands",
        "Align to Center",
        "Show GC Content",
        "Show GC Skew",
    ):
        checkbox = page.get_by_role("checkbox", name=checkbox_name, exact=True)
        checkbox.check()
        expect(checkbox).to_be_checked()

    axis = page.get_by_label("Axis & Scale", exact=True)
    axis.click()
    page.get_by_label("Show Coordinate Scale (Linear)", exact=True).check()
    page.get_by_label("Linear scale style", exact=True).select_option("ruler")
    axis.click()

    colors = page.get_by_label("Colors", exact=True)
    colors.click()
    page.get_by_label("Palette", exact=True).select_option("ajisai")
    colors.click()

    title_panel = page.get_by_label("Title & Legend", exact=True)
    title_panel.click()
    page.get_by_label("Plot Title", exact=True).fill(title)
    page.get_by_label("Plot Title Position", exact=True).select_option("top")
    page.get_by_label("Legend Position", exact=True).select_option("right")
    title_panel.click()


def _configure_all_record_collinear(
    page: Page, *, output_prefix: str, evidence_scope: str = "all"
) -> None:
    global_source = page.locator('[data-capture="linear-blast-source"]')
    expect(global_source).to_have_count(1)
    global_source.get_by_role("radio", name="Run LOSAT", exact=True).check()
    page.get_by_role("radio", name="LOSATP", exact=True).check()
    page.get_by_label("LOSAT execution", exact=True).select_option(
        "threaded" if evidence_scope == "all" else "auto"
    )
    total_threads = page.get_by_label("LOSAT total threads", exact=True)
    total_threads.select_option("32" if evidence_scope == "all" else "safe")
    expect(total_threads).to_have_value("32" if evidence_scope == "all" else "safe")
    threads = page.get_by_label("LOSAT threads per run", exact=True)
    if threads.is_enabled():
        threads.select_option("32" if evidence_scope == "all" else "auto")
    parallel_runs = page.get_by_label("LOSAT parallel runs", exact=True)
    if evidence_scope == "all":
        parallel_runs.select_option("1")
    else:
        parallel_runs.select_option(index=0)
    page.get_by_label("LOSATP blastp mode", exact=True).select_option("collinear")
    page.get_by_label("Collinear max unit gap", exact=True).fill("0")
    page.get_by_label("Collinear minimum block genes", exact=True).fill("1")
    page.get_by_label("Collinear color mode", exact=True).select_option(
        "orientation_identity"
    )
    page.get_by_label("Collinear evidence scope", exact=True).select_option(
        evidence_scope
    )
    page.get_by_label("Collinear diagonal drift", exact=True).fill("0")
    page.get_by_label("Collinear merge conflicts", exact=True).fill("1")
    page.get_by_label("LOSATP minimum bitscore", exact=True).fill("50")
    page.get_by_label("LOSATP maximum e-value", exact=True).fill("0.01")
    page.get_by_label("LOSATP minimum identity", exact=True).fill("0")
    page.get_by_label("LOSATP minimum alignment length", exact=True).fill("0")
    pairwise = page.get_by_label("Pairwise Match", exact=True)
    pairwise.click()
    page.get_by_label("Pairwise Match Style", exact=True).select_option("curve")
    pairwise.click()
    page.get_by_label("Output Prefix", exact=True).fill(output_prefix)


def _assert_records(report: dict[str, Any]) -> None:
    if report.get("svgCount") != 1:
        raise AssertionError("Expected one Hepatoplasmataceae SVG")
    if set(report.get("recordIds", [])) != set(RECORD_IDS):
        raise AssertionError("The five complete Hepatoplasmataceae records changed")
    if report.get("featureElementCount") != EXPECTED_FEATURE_ELEMENTS:
        raise AssertionError(
            "Expected 2,994 rendered feature elements, found "
            f"{report.get('featureElementCount')}"
        )
    texts = set(report.get("texts", []))
    for record_id, length_text in RECORD_LENGTH_TEXT.items():
        if record_id not in texts or length_text not in texts:
            raise AssertionError(f"Missing complete record text: {record_id}")
    assert_static_svg_safety(report)


def _assert_all_vs_all_groups(report: dict[str, Any]) -> None:
    _assert_records(report)
    matches = report.get("comparisonMatches", [])
    if not matches or {match.get("kind") for match in matches} != {"orthogroup"}:
        raise AssertionError("The starting project has no Similarity-group links")


def _assert_collinear(
    report: dict[str, Any], *, title: str, evidence_scope: str
) -> None:
    _assert_records(report)
    texts = set(report.get("texts", []))
    expected_texts = ["GC content", "GC skew (+)", "GC skew (-)"]
    if title:
        expected_texts.insert(0, title)
    for expected in expected_texts:
        if expected not in texts:
            raise AssertionError(f"Hepatoplasmataceae figure is missing {expected!r}")
    if not report.get("coordinateTicks"):
        raise AssertionError("Hepatoplasmataceae figure has no ruler ticks")

    matches = report.get("comparisonMatches", [])
    if not matches or {match.get("kind") for match in matches} != {"collinear"}:
        raise AssertionError("No LOSATP Collinear blocks were rendered")
    endpoints = {(match.get("query"), match.get("subject")) for match in matches}
    if endpoints != ADJACENT_PAIRS:
        raise AssertionError("Collinear ribbons do not connect every adjacent row")
    orientations = Counter(match.get("orientation") for match in matches)
    if not orientations["plus"] or not orientations["minus"]:
        raise AssertionError(f"Collinear orientations changed: {orientations!r}")
    if {match.get("colorMode") for match in matches} != {
        "orientation_identity"
    }:
        raise AssertionError("Collinear blocks lost orientation-and-identity colors")
    group_scopes = {match.get("groupScope") for match in matches}
    if evidence_scope == "all":
        if group_scopes != {"global_collinear"}:
            raise AssertionError(
                "Collinear blocks were not reduced from all-record evidence"
            )
    elif group_scopes - {"", "adjacent_collinear", "adjacent_local"}:
        raise AssertionError(
            f"Adjacent Collinear blocks have an unexpected scope: {group_scopes!r}"
        )


def _cache_state(page: Page) -> dict[str, Any]:
    return page.evaluate(
        """
        () => ({
          mode: window.__GBDRAW_APP__?.losat?.blastp?.mode,
          scope: window.__GBDRAW_APP__?.losat?.blastp?.collinearSearchScope,
          telemetry: window.__GBDRAW_APP__?.lastRunInfo?.losatTelemetry ||
            globalThis.__GBDRAW_LAST_LOSAT_TELEMETRY__ || null,
          cache: (window.__GBDRAW_APP__?.losatCacheInfo || []).map((entry) => ({
            key: String(entry?.key || ''),
            edgeKey: String(entry?.edgeKey || '')
          }))
        })
        """
    )


def _assert_loaded_all_vs_all(
    page: Page, *, session_keys: frozenset[str]
) -> None:
    state = _cache_state(page)
    cache = state.get("cache", [])
    visible_keys = {entry.get("key", "") for entry in cache}
    if (
        state.get("mode") != "orthogroup"
        or len(cache) != 4
        or not visible_keys <= session_keys
    ):
        raise AssertionError(f"Starting project did not restore its cache: {state!r}")


def _assert_all_record_search(
    page: Page, *, session_keys: frozenset[str]
) -> None:
    state = _cache_state(page)
    if state.get("mode") != "collinear" or state.get("scope") != "all":
        raise AssertionError(f"Collinear evidence scope is not all: {state!r}")
    cache = state.get("cache", [])
    keys = {entry.get("key", "") for entry in cache}
    telemetry = state.get("telemetry") or {}
    if (
        len(cache) != 4
        or not keys <= session_keys
        or telemetry.get("totalPairs") != 25
        or telemetry.get("cacheHits") != 25
        or telemetry.get("cacheMisses") != 0
        or telemetry.get("uniqueJobs") != 0
    ):
        raise AssertionError(f"Collinear conversion reran the all-vs-all search: {state!r}")


def _assert_gallery_adjacent_search(
    page: Page, *, session_keys: frozenset[str]
) -> None:
    state = _cache_state(page)
    cache = state.get("cache", [])
    keys = {entry.get("key", "") for entry in cache}
    telemetry = state.get("telemetry") or {}
    if (
        state.get("mode") != "collinear"
        or state.get("scope") != "adjacent"
        or len(cache) != 4
        or not keys <= session_keys
        or telemetry.get("cacheMisses") not in (None, 0)
        or telemetry.get("uniqueJobs") not in (None, 0)
    ):
        raise AssertionError(
            f"Gallery Collinear project did not reuse adjacent evidence: {state!r}"
        )


def _download_button(page: Page, button: Any, download_dir: Path, name: str) -> Path:
    with page.expect_download(timeout=ACTION_TIMEOUT_MS) as info:
        button.click()
    download = info.value
    if download.failure() is not None:
        raise AssertionError(f"Browser download failed: {download.failure()}")
    path = download_dir / name
    download.save_as(path)
    return path


def _download_collinear_envelopes(
    page: Page, popup: Any, download_dir: Path
) -> dict[str, Any]:
    buttons = popup.get_by_role("button", name=re.compile(r"FASTA"))
    expect(buttons).to_have_count(3)
    path = _download_button(page, buttons.last, download_dir, "collinear_members.fasta")
    records = list(SeqIO.parse(path, "fasta"))
    source_sequences = {
        record.id: str(record.seq).upper()
        for fixture in GUI_HEPATOPLASMATACEAE_FIXTURES
        for record in SeqIO.parse(fixture[0], "genbank")
    }
    if len(records) != 2:
        raise AssertionError("Collinear envelope export does not contain two records")
    for record in records:
        source_id = next(
            (record_id for record_id in source_sequences if record_id in record.id),
            None,
        )
        sequence = str(record.seq).upper()
        if source_id is None or (
            sequence not in source_sequences[source_id]
            and str(Seq(sequence).reverse_complement())
            not in source_sequences[source_id]
        ):
            raise AssertionError("A Collinear envelope differs from its source record")
    return {
        "filename": path.name,
        "bytes": path.stat().st_size,
        "records": len(records),
        "record_ids": [record.id for record in records],
    }


def capture_hepatoplasmataceae_collinear(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
    *,
    output_prefix: str,
    screenshot_names: Mapping[str, str],
    starting_project: str = "all-vs-all",
) -> HepatoplasmataceaeCollinearResult:
    """Run one cached LOSATP-to-Collinear project journey."""

    source_records = _validate_fixtures()
    for path in output_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    download_dir.mkdir(parents=True, exist_ok=True)

    if starting_project not in {"all-vs-all", "gallery-collinear"}:
        raise ValueError(f"Unsupported Collinear starting project: {starting_project}")
    evidence_scope = "all" if starting_project == "all-vs-all" else "adjacent"
    title = (
        "All-vs-all LOSATP Collinear blocks across Hepatoplasmataceae"
        if evidence_scope == "all"
        else "LOSATP Collinear blocks across Hepatoplasmataceae"
    )
    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshots: dict[str, int] = {}
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        session_cache_keys = (
            _load_all_vs_all_session(page)
            if starting_project == "all-vs-all"
            else _load_gallery_collinear_session(page)
        )
        result_region = page.get_by_role("region", name="Result Preview", exact=True)
        expect(result_region).to_be_visible(timeout=ACTION_TIMEOUT_MS)
        restored_report = inspect_gui_bgc_losatp_svg(result_region)
        if starting_project == "all-vs-all":
            _assert_all_vs_all_groups(restored_report)
            _assert_loaded_all_vs_all(page, session_keys=session_cache_keys)
        else:
            _assert_collinear(restored_report, title="", evidence_scope="adjacent")
            restored_state = _cache_state(page)
            if (
                restored_state.get("mode") != "collinear"
                or restored_state.get("scope") != "adjacent"
            ):
                raise AssertionError(
                    f"Gallery session did not restore Collinear mode: {restored_state!r}"
                )
        set_feature_search_visible(page, visible=False)
        if "input" in screenshot_names:
            screenshots[screenshot_names["input"]] = capture_screenshot(
                page, output_paths[screenshot_names["input"]], "Linear"
            )

        if "plain" in screenshot_names:
            fit_complete_linear_preview(page, target_zoom="30%")
            screenshots[screenshot_names["plain"]] = capture_screenshot(
                page, output_paths[screenshot_names["plain"]], "Linear"
            )

        _set_presentation(page, title=title)
        _configure_all_record_collinear(
            page,
            output_prefix=output_prefix,
            evidence_scope=evidence_scope,
        )
        page.get_by_label("Collinear evidence scope", exact=True).scroll_into_view_if_needed()
        screenshots[screenshot_names["settings"]] = capture_screenshot(
            page, output_paths[screenshot_names["settings"]], "Linear"
        )

        final_report = generate_and_inspect(
            page,
            inspect_gui_bgc_losatp_svg,
            lambda report: _assert_collinear(
                report, title=title, evidence_scope=evidence_scope
            ),
            timeout_ms=600_000,
        )
        if evidence_scope == "all":
            _assert_all_record_search(page, session_keys=session_cache_keys)
        else:
            _assert_gallery_adjacent_search(page, session_keys=session_cache_keys)
        fit_complete_linear_preview(
            page, target_zoom="40%", pan_left=True
        )
        screenshots[screenshot_names["result"]] = capture_screenshot(
            page, output_paths[screenshot_names["result"]], "Linear"
        )

        first_match = page.get_by_role(
            "region", name="Result Preview", exact=True
        ).locator('[data-match-kind="collinear"]').first
        expect(first_match).to_be_visible()
        first_match.click()
        popup = page.get_by_role("dialog", name="Pairwise match details", exact=True)
        expect(popup).to_be_visible()
        popup_text = popup.inner_text()
        for token in ("Block ID", "Anchor", "Orientation", "Collinear block spans"):
            if token.casefold() not in popup_text.casefold():
                raise AssertionError(f"Collinear popup is missing {token!r}")
        screenshots[screenshot_names["popup"]] = capture_screenshot(
            page, output_paths[screenshot_names["popup"]], "Linear"
        )

        fasta_report = _download_collinear_envelopes(page, popup, download_dir)
        svg_path = _download_button(
            page,
            page.get_by_role("button", name="SVG", exact=True),
            download_dir,
            f"{output_prefix}.svg",
        )
        download_report = {
            "filename": svg_path.name,
            "bytes": svg_path.stat().st_size,
            "semantics": inspect_svg_file(svg_path),
        }
        _assert_collinear(
            download_report["semantics"],
            title=title,
            evidence_scope=evidence_scope,
        )
        capture.assert_clean()
    finally:
        capture.close()

    return HepatoplasmataceaeCollinearResult(
        screenshot_bytes=screenshots,
        final_svg_semantics=final_report,
        download=download_report,
        fasta_download=fasta_report,
        popup={"mode": "collinear", "text": popup_text},
        block_count=len(final_report["comparisonMatches"]),
        source_records=source_records,
    )


__all__ = [
    "HepatoplasmataceaeCollinearResult",
    "capture_hepatoplasmataceae_collinear",
]
