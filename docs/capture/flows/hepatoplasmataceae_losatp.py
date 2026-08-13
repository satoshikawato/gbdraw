"""Visible-UI capture from raw Hepatoplasmataceae records to Collinear blocks."""

from __future__ import annotations

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
from config import (
    ACTION_TIMEOUT_MS,
    GUI_HEPATOPLASMATACEAE_FIXTURES,
    VIEWPORT_HEIGHT,
    VIEWPORT_WIDTH,
    PYTHON_OPERATION_TIMEOUT_MS,
)
from flows.web_capture import (
    assert_fixture_identity,
    capture_screenshot,
    fit_complete_linear_preview,
    generate_and_inspect,
    open_browser_capture,
    open_linear_comparison_disclosure,
    select_linear_losat_mode,
    set_feature_search_visible,
    wait_for_app_shell,
)


RECORD_IDS = tuple(fixture[3] for fixture in GUI_HEPATOPLASMATACEAE_FIXTURES)
RECORD_LENGTH_TEXT = {
    record_id: f"{length:,} bp"
    for _, _, _, record_id, length, _ in GUI_HEPATOPLASMATACEAE_FIXTURES
}
ADJACENT_PAIRS = set(zip(RECORD_IDS[:-1], RECORD_IDS[1:], strict=True))
EXPECTED_FEATURE_ELEMENTS = 2_994
EXPECTED_COLLINEAR_MATCH_ELEMENTS = 500
ALL_RECORD_GENERATION_TIMEOUT_MS = 1_200_000
POPUP_MARGIN_PX = 12
PREVIEW_INSET_PX = 12
PREVIEW_TARGET_GUTTER_PX = 24


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


def _set_source_inputs(page: Page) -> None:
    linear = page.get_by_role("button", name="Linear", exact=True)
    linear.click()
    expect(linear).to_have_attribute("aria-pressed", "true")
    page.get_by_role("radio", name="GenBank", exact=True).check()
    expect(page.get_by_role("status").filter(has_text="Current:")).to_contain_text(
        "Current: No comparison"
    )

    first_fixture = GUI_HEPATOPLASMATACEAE_FIXTURES[0]
    page.get_by_test_id("linear-genbank-1").set_input_files(first_fixture[0])
    add_sequence = page.get_by_role("button", name="Add sequence", exact=True)
    expect(add_sequence).to_have_count(2)
    for index, fixture in enumerate(
        GUI_HEPATOPLASMATACEAE_FIXTURES[1:], start=2
    ):
        add_sequence.first.click()
        page.get_by_test_id(f"linear-genbank-{index}").set_input_files(fixture[0])

    selected_files = page.get_by_role(
        "group", name="GenBank File selection", exact=True
    )
    expect(selected_files).to_have_count(5)
    for index, fixture in enumerate(GUI_HEPATOPLASMATACEAE_FIXTURES):
        expect(selected_files.nth(index)).to_contain_text(fixture[0].name)

    for index, fixture in enumerate(GUI_HEPATOPLASMATACEAE_FIXTURES, start=1):
        record_options = page.get_by_role(
            "button",
            name=f"Record options for sequence {index}",
            exact=True,
        )
        record_options.click()
        selector = page.get_by_label(
            f"Record selector for sequence {index}", exact=True
        )
        expect(selector).to_be_enabled(timeout=PYTHON_OPERATION_TIMEOUT_MS)
        expect(selector).to_contain_text(
            fixture[3], timeout=PYTHON_OPERATION_TIMEOUT_MS
        )
        expect(selector).to_have_value("")
        expect(
            page.get_by_role(
                "group", name=f"Linear sequence {index}", exact=True
            )
        ).not_to_contain_text("Loading records...")
        expect(
            page.get_by_label(f"Region start for sequence {index}", exact=True)
        ).to_have_value("")
        expect(
            page.get_by_label(f"Region end for sequence {index}", exact=True)
        ).to_have_value("")
        expect(
            page.get_by_label(
                f"Reverse complement for sequence {index}", exact=True
            )
        ).not_to_be_checked()
        record_options.click()


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
) -> tuple[Any, Any]:
    commands = page.get_by_role(
        "group", name="Set all adjacent comparisons", exact=True
    )
    commands.get_by_role(
        "button", name="Run LOSAT for all adjacent pairs", exact=True
    ).click()
    expect(page.get_by_role("status").filter(has_text="Current:")).to_contain_text(
        "Current: Run LOSAT for all adjacent pairs"
    )

    selected_pairs = open_linear_comparison_disclosure(
        page,
        "selected-pairs",
        "Selected pairs",
    )
    final_boundary = selected_pairs.get_by_role(
        "region",
        name="Comparison boundary from display row 4 to 5",
        exact=True,
    )
    expect(final_boundary).to_be_visible()
    expect(final_boundary).to_contain_text("#4 → #5")
    selected_pairs.get_by_role(
        "button", name="Selected pairs", exact=True
    ).click()

    settings = open_linear_comparison_disclosure(
        page,
        "settings",
        "Comparison Settings",
    )
    select_linear_losat_mode(
        settings,
        label="LOSATP",
        mode_key="blastp",
    )
    losatp_mode = settings.get_by_role(
        "combobox", name="LOSATP mode", exact=True
    )
    losatp_mode.select_option("pairwise")
    expect(losatp_mode).to_have_value("pairwise")
    match_style = settings.get_by_label("Pairwise Match Style", exact=True)
    match_style.select_option("curve")
    expect(match_style).to_have_value("curve")
    losatp_mode.select_option("collinear")
    expect(losatp_mode).to_have_value("collinear")

    settings.get_by_label("Collinear max unit gap", exact=True).fill("0")
    settings.get_by_label(
        "Collinear minimum block genes", exact=True
    ).fill("1")
    settings.get_by_label("Collinear color mode", exact=True).select_option(
        "orientation_identity"
    )
    settings.get_by_label("Collinear evidence scope", exact=True).select_option(
        evidence_scope
    )
    settings.get_by_label(
        "Linear comparison minimum bitscore", exact=True
    ).fill("50")
    settings.get_by_label(
        "Linear comparison maximum e-value", exact=True
    ).fill("0.01")
    settings.get_by_label(
        "Linear comparison minimum identity", exact=True
    ).fill("0")
    settings.get_by_label(
        "Linear comparison minimum alignment length", exact=True
    ).fill("0")

    advanced = open_linear_comparison_disclosure(
        page,
        "advanced",
        "Advanced comparison and layout",
    )
    advanced.get_by_label("LOSAT execution", exact=True).select_option(
        "threaded" if evidence_scope == "all" else "auto"
    )
    total_threads = advanced.get_by_label("LOSAT total threads", exact=True)
    total_threads.select_option("32" if evidence_scope == "all" else "safe")
    expect(total_threads).to_have_value("32" if evidence_scope == "all" else "safe")
    threads = advanced.get_by_label("LOSAT threads per run", exact=True)
    if threads.is_enabled():
        threads.select_option("8" if evidence_scope == "all" else "auto")
    parallel_runs = advanced.get_by_label("LOSAT parallel runs", exact=True)
    if evidence_scope == "all":
        parallel_runs.select_option("4")
    else:
        parallel_runs.select_option(index=0)
    advanced.get_by_label("Collinear diagonal drift", exact=True).fill("0")
    advanced.get_by_label("Collinear merge conflicts", exact=True).fill("1")
    page.get_by_label("Output Prefix", exact=True).fill(output_prefix)
    return settings, advanced


def _assert_records(report: dict[str, Any]) -> None:
    if report.get("svgCount") != 1:
        raise AssertionError("Expected one Hepatoplasmataceae SVG")
    ordered_record_ids = tuple(dict.fromkeys(report.get("recordIds", [])))
    if ordered_record_ids != RECORD_IDS:
        raise AssertionError(
            "The five complete Hepatoplasmataceae records changed order or identity"
        )
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


def _assert_plain_records(report: dict[str, Any]) -> None:
    _assert_records(report)
    if report.get("comparisonMatches"):
        raise AssertionError("The five-record baseline unexpectedly has comparison links")


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
    if len(matches) != EXPECTED_COLLINEAR_MATCH_ELEMENTS:
        raise AssertionError(
            "Expected 500 rendered Collinear match elements, found "
            f"{len(matches)}"
        )
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
          threading: window.__GBDRAW_APP__?.losatThreadingStatus || null,
          scheduling: {
            totalThreads: String(window.__GBDRAW_APP__?.losat?.totalThreadBudget || ''),
            parallelRuns: String(window.__GBDRAW_APP__?.losat?.parallelWorkers || ''),
            threadsPerRun: String(window.__GBDRAW_APP__?.losat?.threadsPerJob || '')
          },
          cache: (window.__GBDRAW_APP__?.losatCacheInfo || []).map((entry) => ({
            key: String(entry?.key || ''),
            edgeKey: String(entry?.edgeKey || '')
          }))
        })
        """
    )


def _assert_empty_cache(page: Page) -> None:
    state = _cache_state(page)
    if state.get("cache") or state.get("telemetry") not in (None, {}):
        raise AssertionError(
            f"Raw-input Collinear journey did not start with an empty cache: {state!r}"
        )


def _assert_fresh_search(page: Page, *, evidence_scope: str) -> None:
    state = _cache_state(page)
    if state.get("mode") != "collinear" or state.get("scope") != evidence_scope:
        raise AssertionError(
            f"Collinear evidence scope is not {evidence_scope}: {state!r}"
        )
    cache = state.get("cache", [])
    telemetry = state.get("telemetry") or {}
    expected_jobs = 25 if evidence_scope == "all" else 13
    if (
        len(cache) != 4
        or telemetry.get("totalPairs") != expected_jobs
        or telemetry.get("cacheHits") != 0
        or telemetry.get("cacheMisses") != expected_jobs
        or telemetry.get("uniqueJobs") != expected_jobs
    ):
        raise AssertionError(f"Collinear run did not execute a fresh search: {state!r}")
    if evidence_scope == "all":
        threading = state.get("threading") or {}
        scheduling = state.get("scheduling") or {}
        if (
            threading.get("state") != "available"
            or threading.get("mode") != "threaded"
            or threading.get("pairWorkers") != 4
            or threading.get("threadsPerJob") != 8
            or scheduling
            != {
                "totalThreads": "32",
                "parallelRuns": "4",
                "threadsPerRun": "8",
            }
        ):
            raise AssertionError(
                "The all-record Collinear run used an unexpected runtime plan: "
                f"{state!r}"
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


def _require_bounding_box(locator: Any, name: str) -> dict[str, float]:
    box = locator.bounding_box()
    if box is None:
        raise AssertionError(f"Could not resolve {name} bounds")
    return box


def _box_is_inside(
    inner: Mapping[str, float], outer: Mapping[str, float], *, tolerance: float = 1
) -> bool:
    return (
        inner["x"] >= outer["x"] - tolerance
        and inner["y"] >= outer["y"] - tolerance
        and inner["x"] + inner["width"]
        <= outer["x"] + outer["width"] + tolerance
        and inner["y"] + inner["height"]
        <= outer["y"] + outer["height"] + tolerance
    )


def _boxes_overlap(left: Mapping[str, float], right: Mapping[str, float]) -> bool:
    return not (
        left["x"] + left["width"] <= right["x"]
        or right["x"] + right["width"] <= left["x"]
        or left["y"] + left["height"] <= right["y"]
        or right["y"] + right["height"] <= left["y"]
    )


def _assert_input_capture_framing(page: Page) -> None:
    selected_files = page.get_by_role(
        "group", name="GenBank File selection", exact=True
    )
    expect(selected_files).to_have_count(5)
    selected_files.nth(4).evaluate(
        "(element) => element.scrollIntoView({ block: 'center' })"
    )

    fourth_box = _require_bounding_box(selected_files.nth(3), "fourth input")
    fifth_box = _require_bounding_box(selected_files.nth(4), "fifth input")
    viewport = {
        "x": 0,
        "y": 0,
        "width": VIEWPORT_WIDTH,
        "height": VIEWPORT_HEIGHT,
    }
    for name, box in (("fourth input", fourth_box), ("fifth input", fifth_box)):
        if not _box_is_inside(box, viewport, tolerance=0):
            raise AssertionError(f"The input capture clips the {name}: {box!r}")
    if fourth_box["y"] + fourth_box["height"] > fifth_box["y"]:
        raise AssertionError(
            "The fourth and fifth input uploaders are not visibly ordered in "
            "the capture frame"
        )


def _result_preview_content_box(result_preview: Any) -> dict[str, float]:
    geometry = result_preview.evaluate(
        """
        (region, inset) => {
          const svgs = Array.from(region.getElementsByTagName('svg'));
          const svg = svgs.find((candidate) => (
            Array.from(candidate.getElementsByTagName('g')).some((group) => (
              group.dataset?.gbdrawRole === 'record-definition'
            ))
          ));
          const wrapper = svg?.parentElement;
          const canvas = wrapper?.parentElement;
          if (!svg || !wrapper || !canvas || !region.contains(canvas)) return null;
          const regionRect = region.getBoundingClientRect();
          const canvasRect = canvas.getBoundingClientRect();
          const left = Math.max(regionRect.left, canvasRect.left) + inset;
          const right = Math.min(regionRect.right, canvasRect.right) - inset;
          const top = Math.max(regionRect.top, canvasRect.top) + inset;
          const bottom = Math.min(regionRect.bottom, canvasRect.bottom) - inset;
          return {
            x: left,
            y: top,
            width: right - left,
            height: bottom - top
          };
        }
        """,
        PREVIEW_INSET_PX,
    )
    if geometry is None or geometry["width"] <= 0 or geometry["height"] <= 0:
        raise AssertionError("Could not resolve the inset Result Preview content box")
    return geometry


def _record_definition_geometry(result_preview: Any) -> tuple[dict[str, Any], ...]:
    definitions = result_preview.evaluate(
        """
        (region) => Array.from(region.getElementsByTagName('g'))
          .filter((group) => group.dataset?.gbdrawRole === 'record-definition')
          .map((group) => {
            const bounds = group.getBoundingClientRect();
            return {
              index: Number(group.dataset?.gbdrawRecordIndex),
              recordId: String(group.dataset?.gbdrawRecordId || ''),
              text: String(group.textContent || '').replace(/\\s+/g, ' ').trim(),
              bounds: {
                x: Number(bounds.x),
                y: Number(bounds.y),
                width: Number(bounds.width),
                height: Number(bounds.height)
              }
            };
          })
          .sort((left, right) => left.index - right.index)
        """
    )
    return tuple(definitions)


def _definition_union(
    definitions: tuple[dict[str, Any], ...],
) -> dict[str, float]:
    if not definitions:
        raise AssertionError("The Linear preview has no record-definition groups")
    left = min(item["bounds"]["x"] for item in definitions)
    top = min(item["bounds"]["y"] for item in definitions)
    right = max(
        item["bounds"]["x"] + item["bounds"]["width"]
        for item in definitions
    )
    bottom = max(
        item["bounds"]["y"] + item["bounds"]["height"]
        for item in definitions
    )
    return {
        "x": left,
        "y": top,
        "width": right - left,
        "height": bottom - top,
    }


def _blank_preview_drag_point(
    result_preview: Any, *, delta_x: float
) -> dict[str, float]:
    point = result_preview.evaluate(
        """
        (region, options) => {
          const svgs = Array.from(region.getElementsByTagName('svg'));
          const svg = svgs.find((candidate) => (
            Array.from(candidate.getElementsByTagName('g')).some((group) => (
              group.dataset?.gbdrawRole === 'record-definition'
            ))
          ));
          const canvas = svg?.parentElement?.parentElement;
          if (!svg || !canvas || !region.contains(canvas)) return null;

          const canvasRect = canvas.getBoundingClientRect();
          const regionRect = region.getBoundingClientRect();
          const left = Math.max(canvasRect.left, regionRect.left) + options.inset;
          const right = Math.min(canvasRect.right, regionRect.right) - options.inset;
          const top = Math.max(canvasRect.top, regionRect.top) + options.inset;
          const bottom = Math.min(canvasRect.bottom, regionRect.bottom) - options.inset;

          const isBlankCanvasTarget = (target) => {
            let current = target;
            while (current && current !== canvas) {
              const tag = String(current.tagName || '').toLowerCase();
              const dataset = current.dataset || {};
              if (
                ['button', 'input', 'textarea', 'select', 'a', 'label'].includes(tag)
                || current.getAttribute?.('role') === 'button'
                || String(current.id || '').startsWith('f')
                || dataset.nodrag === 'true'
                || dataset.gbdrawFeatureId
                || dataset.gbdrawPairwiseMatchId
                || dataset.matchKind
                || dataset.pairwiseMatchStyle
                || dataset.collinearityBlockId
                || dataset.collinearGroupScope
              ) return false;
              current = current.parentElement;
            }
            return current === canvas;
          };

          const xFractions = [0.12, 0.25, 0.38, 0.5, 0.62, 0.75, 0.88];
          const yFractions = [0.15, 0.28, 0.42, 0.58, 0.72, 0.85];
          for (const yFraction of yFractions) {
            const y = top + ((bottom - top) * yFraction);
            for (const xFraction of xFractions) {
              const x = left + ((right - left) * xFraction);
              const endX = x + options.deltaX;
              if (endX < left || endX > right) continue;
              const startTarget = document.elementFromPoint(x, y);
              const endTarget = document.elementFromPoint(endX, y);
              if (
                isBlankCanvasTarget(startTarget)
                && endTarget
                && canvas.contains(endTarget)
              ) return { x, y, endX };
            }
          }
          return null;
        }
        """,
        {"deltaX": delta_x, "inset": PREVIEW_INSET_PX},
    )
    if point is None:
        raise AssertionError(
            "Could not find a visible blank Result Preview drag path for "
            f"{delta_x:.1f}px"
        )
    return point


def _drag_preview_horizontally(
    page: Page, result_preview: Any, *, delta_x: float, label: str
) -> None:
    content = _result_preview_content_box(result_preview)
    max_delta = content["width"] * 0.7
    if abs(delta_x) > max_delta:
        raise AssertionError(
            f"{label} requires an unsafe Result Preview drag: {delta_x:.1f}px"
        )
    point = _blank_preview_drag_point(result_preview, delta_x=delta_x)
    page.mouse.move(point["x"], point["y"])
    page.mouse.down()
    page.mouse.move(point["endX"], point["y"], steps=12)
    page.mouse.up()
    page.wait_for_timeout(250)
    page.evaluate("() => window.getSelection()?.removeAllRanges()")
    if page.evaluate("() => window.getSelection()?.rangeCount ?? 0") != 0:
        raise AssertionError(f"{label} left a text selection after preview panning")


def _assert_record_definitions_inside_preview(result_preview: Any) -> None:
    definitions = _record_definition_geometry(result_preview)
    if tuple(item["recordId"] for item in definitions) != RECORD_IDS:
        raise AssertionError(
            "The 80% baseline does not contain the five ordered record definitions"
        )
    content = _result_preview_content_box(result_preview)
    for item in definitions:
        record_id = item["recordId"]
        if (
            record_id not in item["text"]
            or RECORD_LENGTH_TEXT[record_id] not in item["text"]
        ):
            raise AssertionError(
                f"The 80% baseline definition is missing ID or length: {item!r}"
            )
        if not _box_is_inside(item["bounds"], content, tolerance=0):
            raise AssertionError(
                "A record definition is clipped by the inset Result Preview: "
                f"{item!r}; content={content!r}"
            )


def _prepare_plain_definition_detail(page: Page) -> None:
    reset_zoom = page.get_by_role("button", name="Reset zoom", exact=True)
    zoom_in = page.get_by_role("button", name="Zoom in", exact=True)
    expect(reset_zoom).to_contain_text("40%")
    for _ in range(4):
        zoom_in.click()
    expect(reset_zoom).to_contain_text("80%")
    page.wait_for_timeout(250)

    result_preview = page.get_by_role("region", name="Result Preview", exact=True)
    definitions = _record_definition_geometry(result_preview)
    if tuple(item["recordId"] for item in definitions) != RECORD_IDS:
        raise AssertionError("Could not resolve all five baseline record definitions")
    definition_bounds = _definition_union(definitions)
    content = _result_preview_content_box(result_preview)
    if definition_bounds["width"] > content["width"]:
        raise AssertionError("The 80% record-definition column cannot fit the preview")
    target_x = content["x"] + PREVIEW_TARGET_GUTTER_PX
    _drag_preview_horizontally(
        page,
        result_preview,
        delta_x=target_x - definition_bounds["x"],
        label="80% record-definition detail",
    )
    _assert_record_definitions_inside_preview(result_preview)


def _return_to_overview_zoom(page: Page) -> None:
    reset_zoom = page.get_by_role("button", name="Reset zoom", exact=True)
    zoom_out = page.get_by_role("button", name="Zoom out", exact=True)
    expect(reset_zoom).to_contain_text("80%")
    for _ in range(4):
        zoom_out.click()
    expect(reset_zoom).to_contain_text("40%")
    page.wait_for_timeout(250)


def _assert_match_inside_result_preview(
    result_preview: Any, first_match: Any
) -> None:
    content = _result_preview_content_box(result_preview)
    match_box = _require_bounding_box(first_match, "Pairwise match 1")
    if not _box_is_inside(match_box, content, tolerance=0):
        raise AssertionError(
            "Pairwise match 1 is clipped by the inset Result Preview: "
            f"{match_box!r}; content={content!r}"
        )


def _prepare_collinear_detail(page: Page) -> Any:
    reset_zoom = page.get_by_role("button", name="Reset zoom", exact=True)
    zoom_in = page.get_by_role("button", name="Zoom in", exact=True)
    expect(reset_zoom).to_contain_text("40%")
    for _ in range(4):
        zoom_in.click()
    expect(reset_zoom).to_contain_text("80%")
    page.wait_for_timeout(250)

    result_preview = page.get_by_role("region", name="Result Preview", exact=True)
    first_match = result_preview.get_by_role(
        "button", name="Pairwise match 1", exact=True
    )
    expect(first_match).to_be_visible()
    match_box = _require_bounding_box(first_match, "Pairwise match 1")
    content = _result_preview_content_box(result_preview)
    if match_box["width"] > content["width"]:
        raise AssertionError("Pairwise match 1 cannot fit the Result Preview")
    target_x = content["x"] + PREVIEW_TARGET_GUTTER_PX
    _drag_preview_horizontally(
        page,
        result_preview,
        delta_x=target_x - match_box["x"],
        label="80% Pairwise match 1 detail",
    )
    _assert_match_inside_result_preview(result_preview, first_match)
    return first_match


def _move_popup_opposite_match(page: Page, popup: Any, first_match: Any) -> None:
    match_box = _require_bounding_box(first_match, "Pairwise match 1")
    popup_box = _require_bounding_box(popup, "Pairwise match popup")
    match_center_x = match_box["x"] + (match_box["width"] / 2)
    target_x = (
        VIEWPORT_WIDTH - popup_box["width"] - POPUP_MARGIN_PX
        if match_center_x < VIEWPORT_WIDTH / 2
        else POPUP_MARGIN_PX
    )
    target_y = POPUP_MARGIN_PX

    drag_x = popup_box["x"] + 24
    drag_y = popup_box["y"] + 24
    page.mouse.move(drag_x, drag_y)
    page.mouse.down()
    page.mouse.move(
        drag_x + target_x - popup_box["x"],
        drag_y + target_y - popup_box["y"],
        steps=12,
    )
    page.mouse.up()
    page.wait_for_timeout(100)

    moved_popup_box = _require_bounding_box(popup, "moved Pairwise match popup")
    viewport = {
        "x": 0,
        "y": 0,
        "width": VIEWPORT_WIDTH,
        "height": VIEWPORT_HEIGHT,
    }
    if not _box_is_inside(moved_popup_box, viewport):
        raise AssertionError(
            "The Pairwise match popup is not fully inside the pinned viewport: "
            f"{moved_popup_box!r}"
        )
    if abs(moved_popup_box["x"] - target_x) > 2 or abs(
        moved_popup_box["y"] - target_y
    ) > 2:
        raise AssertionError(
            "The Pairwise match popup did not reach the opposite top corner: "
            f"{moved_popup_box!r}"
        )

    match_box = _require_bounding_box(first_match, "selected Pairwise match 1")
    result_preview = page.get_by_role("region", name="Result Preview", exact=True)
    _assert_match_inside_result_preview(result_preview, first_match)
    if _boxes_overlap(moved_popup_box, match_box):
        raise AssertionError(
            "The Pairwise match popup geometrically covers the selected match"
        )

    fasta_buttons = popup.get_by_role("button", name=re.compile(r"FASTA"))
    expect(fasta_buttons).to_have_count(3)
    for index in range(3):
        button = fasta_buttons.nth(index)
        expect(button).to_be_visible()
        button_box = _require_bounding_box(button, f"FASTA button {index + 1}")
        if not _box_is_inside(button_box, moved_popup_box):
            raise AssertionError(
                f"FASTA button {index + 1} is clipped from the popup top view"
            )


def _scroll_popup_to_collinearity_details(page: Page, popup: Any) -> None:
    content = popup.get_by_role(
        "region", name="Pairwise match details content", exact=True
    )
    expect(content).to_be_visible()
    before = content.evaluate(
        """
        (node) => ({
          scrollTop: Number(node.scrollTop),
          maxScroll: Number(node.scrollHeight - node.clientHeight)
        })
        """
    )
    if before["maxScroll"] <= 0:
        raise AssertionError("The Pairwise match details have no scrollable content")
    content.hover()
    page.mouse.wheel(0, before["maxScroll"])
    page.wait_for_timeout(150)
    after = content.evaluate(
        """
        (node) => ({
          scrollTop: Number(node.scrollTop),
          maxScroll: Number(node.scrollHeight - node.clientHeight)
        })
        """
    )
    if (
        after["scrollTop"] <= before["scrollTop"]
        or after["scrollTop"] < after["maxScroll"] - 2
    ):
        raise AssertionError(
            "The visible popup scroll did not reach the Collinearity details: "
            f"before={before!r}; after={after!r}"
        )


def _download_collinear_envelopes(
    page: Page, popup: Any, download_dir: Path
) -> dict[str, Any]:
    buttons = popup.get_by_role("button", name=re.compile(r"FASTA"))
    expect(buttons).to_have_count(3)
    with page.expect_download(timeout=ACTION_TIMEOUT_MS) as info:
        buttons.last.click()
    download = info.value
    if download.failure() is not None:
        raise AssertionError(f"Browser download failed: {download.failure()}")
    suggested_name = download.suggested_filename
    if (
        re.fullmatch(
            r"comparison[0-9]+_match[0-9]+_both\.fna", suggested_name
        )
        is None
    ):
        raise AssertionError(
            f"Unexpected Collinear envelope filename: {suggested_name!r}"
        )
    path = download_dir / "collinear_members.fasta"
    download.save_as(path)
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
        "suggested_filename": suggested_name,
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
    evidence_scope: str = "all",
) -> HepatoplasmataceaeCollinearResult:
    """Upload five raw records and run one fresh LOSATP Collinear journey."""

    source_records = _validate_fixtures()
    for path in output_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    download_dir.mkdir(parents=True, exist_ok=True)

    if evidence_scope not in {"all", "adjacent"}:
        raise ValueError(f"Unsupported Collinear evidence scope: {evidence_scope}")
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
        wait_for_app_shell(page)
        _set_source_inputs(page)
        _assert_empty_cache(page)
        page.mouse.move(1_400, 880)
        if "input" in screenshot_names:
            _assert_input_capture_framing(page)
            screenshots[screenshot_names["input"]] = capture_screenshot(
                page, output_paths[screenshot_names["input"]], "Linear"
            )

        if "plain" in screenshot_names:
            generate_and_inspect(
                page,
                inspect_gui_bgc_losatp_svg,
                _assert_plain_records,
                timeout_ms=600_000,
            )
            set_feature_search_visible(page, visible=False)
            fit_complete_linear_preview(page, target_zoom="40%")
            _prepare_plain_definition_detail(page)
            page.mouse.move(1_400, 880)
            screenshots[screenshot_names["plain"]] = capture_screenshot(
                page, output_paths[screenshot_names["plain"]], "Linear"
            )
            _return_to_overview_zoom(page)

        _set_presentation(page, title=title)
        settings, advanced = _configure_all_record_collinear(
            page,
            output_prefix=output_prefix,
            evidence_scope=evidence_scope,
        )
        settings.get_by_role(
            "group", name="LOSAT Mode", exact=True
        ).scroll_into_view_if_needed()
        screenshots[screenshot_names["settings"]] = capture_screenshot(
            page, output_paths[screenshot_names["settings"]], "Linear"
        )

        final_report = generate_and_inspect(
            page,
            inspect_gui_bgc_losatp_svg,
            lambda report: _assert_collinear(
                report, title=title, evidence_scope=evidence_scope
            ),
            timeout_ms=(
                ALL_RECORD_GENERATION_TIMEOUT_MS
                if evidence_scope == "all"
                else 600_000
            ),
        )
        _assert_fresh_search(page, evidence_scope=evidence_scope)
        set_feature_search_visible(page, visible=False)
        fit_complete_linear_preview(page, target_zoom="40%")
        screenshots[screenshot_names["result"]] = capture_screenshot(
            page, output_paths[screenshot_names["result"]], "Linear"
        )

        first_match = _prepare_collinear_detail(page)
        page.mouse.move(VIEWPORT_WIDTH - 40, VIEWPORT_HEIGHT - 30)
        screenshots[screenshot_names["detail"]] = capture_screenshot(
            page, output_paths[screenshot_names["detail"]], "Linear"
        )

        first_match.focus()
        first_match.press("Enter")
        popup = page.get_by_role("dialog", name="Pairwise match details", exact=True)
        expect(popup).to_be_visible()
        expect(first_match).to_have_class(
            re.compile(r"(?:^|\s)gbdraw-match-selected(?:\s|$)")
        )
        popup_text = popup.inner_text()
        for token in ("Block ID", "Anchor", "Orientation", "Collinear block spans"):
            if token.casefold() not in popup_text.casefold():
                raise AssertionError(f"Collinear popup is missing {token!r}")
        _move_popup_opposite_match(page, popup, first_match)
        screenshots[screenshot_names["popup"]] = capture_screenshot(
            page, output_paths[screenshot_names["popup"]], "Linear"
        )

        fasta_report = _download_collinear_envelopes(page, popup, download_dir)
        _scroll_popup_to_collinearity_details(page, popup)
        page.get_by_role(
            "button", name="Close match popup", exact=True
        ).click()
        expect(popup).not_to_be_visible()
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
