"""Capture the four mitochondrial presentation states through public controls."""

from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

from Bio import SeqIO
from playwright.sync_api import expect, sync_playwright

from config import ACTION_TIMEOUT_MS
from flows.how_to.presentation import (
    GENE_PRIORITY_RULE_PATH,
    HUMAN_MITOCHONDRION_PATH,
    _inspect_presentation_svg,
    _open_human_circular,
)
from flows.tutorials.gui_feature_highlight import COLOR_TABLE, REGION_TABLE
from flows.web_capture import generate_and_wait_for_result, open_browser_capture, wait_for_app_shell
from video.capture import rasterize_svg
from video.model import sha256
from web_server import CaptureWebServer


STATES = (
    "human.labels-product", "human.labels-gene", "human.functional-colors",
    "human.dloop-bracket",
)


def _base_tracks(page) -> None:
    page.get_by_role("button", name="Custom Track Slots").click()
    page.get_by_role("checkbox", name="Use custom stack", exact=True).check()
    for slot_id in ("gc_content", "gc_skew"):
        slot = page.get_by_role("group", name=f"Circular track slot {slot_id}", exact=True)
        if slot.count():
            slot.get_by_title("Remove", exact=True).click()
    ticks = page.get_by_role("group", name="Circular track slot ticks", exact=True)
    move = ticks.get_by_title("Move outside Axis", exact=True)
    if move.is_enabled():
        move.click()
    ticks.locator("select").last.select_option("label_out_tick_in")
    page.get_by_role("group", name="Circular track slot features", exact=True).locator("select").last.select_option("split")
    page.get_by_role("button", name="Custom Track Slots").click()


def _add_dloop(page, regions_path: Path) -> None:
    panel = page.get_by_label("Region Annotations", exact=True)
    panel.click()
    page.get_by_label("Import TSV", exact=True).set_input_files(regions_path)
    expect(page.get_by_label("Annotation set id", exact=True)).to_have_value("mitochondrial_regions")
    panel.click()
    page.get_by_role("button", name="Custom Track Slots").click()
    page.get_by_label("New circular track renderer", exact=True).select_option("annotations")
    page.get_by_role("button", name="Add track").click()
    slot = page.get_by_role("group", name="Circular track slot annotations", exact=True)
    slot.get_by_label("Annotation set", exact=True).select_option("mitochondrial_regions")
    slot.get_by_label("Annotation placement", exact=True).select_option("inside")
    slot.get_by_label("Circular track slot id annotations", exact=True).fill("mitochondrial_regions")
    slot = page.get_by_role("group", name="Circular track slot mitochondrial_regions", exact=True)
    slot.get_by_title("Width", exact=True).fill("24px")
    slot.get_by_label("Show annotation labels", exact=True).check()
    slot.locator('select:has(option[value="compress"])').select_option("compress")
    slot.locator('input[type="number"]').last.fill("1")
    page.get_by_role("button", name="Custom Track Slots").click()


def _save_svg(page, path: Path) -> dict:
    report = _inspect_presentation_svg(page.get_by_role("region", name="Result Preview", exact=True))
    if report.get("featureElementCount") != 37 or "NC_012920.1" not in report.get("recordIds", []):
        raise AssertionError("The generated result lacks the 37 source features")
    with page.expect_download(timeout=ACTION_TIMEOUT_MS) as info:
        page.get_by_role("button", name="SVG", exact=True).click()
    download = info.value
    if download.failure() is not None:
        raise AssertionError(f"SVG download failed: {download.failure()}")
    path.parent.mkdir(parents=True, exist_ok=True)
    download.save_as(path)
    return report


def _check_cds(report: dict, labels: set[str], colors: dict[str, str] | None) -> None:
    rows = [row for row in report["features"] if row["type"] == "CDS"]
    if len(rows) != 13:
        raise AssertionError(f"Expected 13 CDS, found {len(rows)}")
    if labels - set(report["texts"]):
        raise AssertionError(f"Missing labels: {sorted(labels - set(report['texts']))}")
    if colors is None:
        if len({row["fill"] for row in rows}) != 1:
            raise AssertionError("CDS base color is not uniform")
    elif any(row["fill"] != colors.get(row["gene"]) for row in rows):
        raise AssertionError("Functional CDS colors do not match the four explicit rules")


def capture_human_edits(run: Path, assets: dict) -> None:
    source = SeqIO.read(HUMAN_MITOCHONDRION_PATH, "genbank")
    cds = [feature for feature in source.features if feature.type == "CDS"]
    if source.id != "NC_012920.1" or len(source) != 16569 or len(cds) != 13:
        raise AssertionError("The fixed mitochondrial source contract changed")
    products = {feature.qualifiers["product"][0] for feature in cds}
    genes = {feature.qualifiers["gene"][0] for feature in cds}
    if len(products) != 13 or len(genes) != 13:
        raise AssertionError("Expected 13 distinct CDS products and genes")
    colors = {
        **{gene: "#3b82f6" for gene in genes if gene.startswith("ND")},
        **{gene: "#ef4444" for gene in genes if gene.startswith("COX")},
        **{gene: "#f59e0b" for gene in genes if gene.startswith("ATP")},
        "CYTB": "#8b5cf6",
    }
    if set(colors) != genes or tuple(sum(value == color for value in colors.values()) for color in
                                    ("#3b82f6", "#ef4444", "#f59e0b", "#8b5cf6")) != (7, 3, 2, 1):
        raise AssertionError("Unexpected CDS functional classification")

    evidence_dir = run / "evidence" / "human-edits"
    tables = evidence_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    product_priority = tables / "cds_product_priority.tsv"
    product_priority.write_text("CDS\tproduct\n", encoding="utf-8")
    colors_path = tables / "presentation_colors.tsv"
    colors_path.write_text("\n".join(COLOR_TABLE.splitlines()[:4]) + "\n", encoding="utf-8")
    regions_path = tables / "mitochondrial_regions.tsv"
    regions_path.write_text(REGION_TABLE, encoding="utf-8")
    recording_dir = run / "raw" / "browser"
    reports = {}
    with CaptureWebServer() as server, sync_playwright() as playwright:
        capture = open_browser_capture(playwright.chromium, server.base_url,
                                       viewport=(1920, 1080), record_video_dir=recording_dir)
        page = capture.page
        recording = page.video
        try:
            page.goto(server.base_url, wait_until="domcontentloaded")
            wait_for_app_shell(page)
            _open_human_circular(page, "meet_gbdraw_human")
            page.get_by_label("Separate Strands", exact=True).uncheck()
            page.get_by_label("Hide GC Content", exact=True).check()
            page.get_by_label("Hide GC Skew", exact=True).check()
            _base_tracks(page)
            labels = page.get_by_label("Labels", exact=True)
            labels.click()
            page.get_by_label("Label Mode", exact=True).select_option("both")
            page.get_by_label("Priority File (TSV)", exact=True).set_input_files(product_priority)
            labels.click()
            title = page.get_by_label("Titles and Record Labels", exact=True)
            title.click()
            page.get_by_role("textbox", name="Plot Title", exact=True).fill("Human mitochondrial genome")
            title.click()
            legend = page.get_by_label("Legend settings", exact=True)
            legend.click()
            page.get_by_label("Legend position", exact=True).select_option("right")
            legend.click()

            for index, asset_id in enumerate(STATES):
                if index == 1:
                    labels.click()
                    page.get_by_label("Priority File (TSV)", exact=True).set_input_files(GENE_PRIORITY_RULE_PATH)
                    labels.click()
                elif index == 2:
                    panel = page.get_by_label("Colors", exact=True)
                    panel.click()
                    page.get_by_label("Specific Table (-t)", exact=True).set_input_files(colors_path)
                    panel.click()
                elif index == 3:
                    _add_dloop(page, regions_path)
                generate_and_wait_for_result(page)
                svg = run / "figures" / f"{asset_id}.svg"
                report = _save_svg(page, svg)
                _check_cds(report, products if index == 0 else genes, colors if index >= 2 else None)
                if ("D-loop" in report["texts"]) != (index == 3):
                    raise AssertionError(f"Incorrect D-loop visibility in {asset_id}")
                reports[asset_id] = {"svg_sha256": sha256(svg), "semantics": report}

            # Turn only the annotation slot off and back on. The 37 source
            # features, labels, and four CDS color groups must survive both.
            slot_panel = page.get_by_role("button", name="Custom Track Slots")
            slot_panel.click()
            slot = page.get_by_role("group", name="Circular track slot mitochondrial_regions", exact=True)
            enabled = slot.get_by_role("checkbox", name="Enable circular track slot mitochondrial_regions", exact=True)
            enabled.uncheck()
            slot_panel.click()
            generate_and_wait_for_result(page)
            before_svg = evidence_dir / "roundtrip-p2.svg"
            before_report = _save_svg(page, before_svg)
            _check_cds(before_report, genes, colors)
            if "D-loop" in before_report["texts"]:
                raise AssertionError("D-loop remained visible after disabling its slot")
            slot_panel.click()
            enabled.check()
            slot_panel.click()
            generate_and_wait_for_result(page)
            after_svg = evidence_dir / "roundtrip-p3.svg"
            after_report = _save_svg(page, after_svg)
            _check_cds(after_report, genes, colors)
            if "D-loop" not in after_report["texts"] or sha256(after_svg) != reports[STATES[-1]]["svg_sha256"]:
                raise AssertionError("P3 was not restored after re-enabling the D-loop slot")
            reports["roundtrip"] = {
                "p2_svg_sha256": sha256(before_svg), "p3_svg_sha256": sha256(after_svg),
                "p2_dloop_visible": False, "p3_matches_original": True,
            }

            button = page.get_by_role("button", name="SVG", exact=True)
            button.scroll_into_view_if_needed()
            page.wait_for_timeout(1000)
            with page.expect_download(timeout=ACTION_TIMEOUT_MS) as info:
                button.click()
            download = info.value
            if download.failure() is not None:
                raise AssertionError(f"Final SVG export failed: {download.failure()}")
            exported = run / "raw" / "final-export.svg"
            exported.parent.mkdir(parents=True, exist_ok=True)
            download.save_as(exported)
            if sha256(exported) != reports[STATES[-1]]["svg_sha256"]:
                raise AssertionError("Final SVG export differs from the P3 figure")
            page.wait_for_timeout(1500)
            capture.assert_clean()
        finally:
            capture.close()
        if recording is None:
            raise AssertionError("Playwright did not provide a video recording")
        raw = run / "raw" / "human-edits-full.webm"
        shutil.copyfile(recording.path(), raw)
        probe = subprocess.run(["ffprobe", "-v", "error", "-show_entries", "format=duration",
                                "-of", "default=noprint_wrappers=1:nokey=1", str(raw)],
                               capture_output=True, text=True, check=True)
        duration = float(probe.stdout.strip())
        clip = run / "raw" / "human.svg-export.mp4"
        subprocess.run(["ffmpeg", "-hide_banner", "-loglevel", "error", "-y", "-ss", str(max(0, duration - 4)),
                        "-i", str(raw), "-vf", "fps=30,scale=1920:1080,format=yuv420p", "-frames:v", "120",
                        "-an", "-c:v", "libx264", "-preset", "veryfast", "-crf", "20", str(clip)], check=True)
        evidence = evidence_dir / "transitions.json"
        evidence.write_text(json.dumps({"source_sha256": sha256(HUMAN_MITOCHONDRION_PATH),
                                        "gene_priority_sha256": sha256(GENE_PRIORITY_RULE_PATH),
                                        "product_priority_sha256": sha256(product_priority),
                                        "color_table_sha256": sha256(colors_path),
                                        "region_table_sha256": sha256(regions_path),
                                        "recording_sha256": sha256(raw), "states": reports}, indent=2), encoding="utf-8")
        browser = playwright.chromium.launch(headless=True)
        try:
            for asset_id in STATES:
                svg = run / "figures" / f"{asset_id}.svg"
                png = run / "figures" / f"{asset_id}.png"
                width, height = rasterize_svg(browser, svg, png)
                assets[asset_id] = {
                    "kind": "image", "path": str(png.relative_to(run)), "sha256": sha256(png),
                    "width": width, "height": height,
                    "svg": {"path": str(svg.relative_to(run)), "sha256": sha256(svg)},
                    "evidence": {"path": str(evidence.relative_to(run)), "sha256": sha256(evidence)},
                }
            assets["human.svg-export"] = {
                "kind": "video", "path": str(clip.relative_to(run)), "sha256": sha256(clip),
                "width": 1920, "height": 1080,
                "evidence": {"path": str(evidence.relative_to(run)), "sha256": sha256(evidence)},
            }
        finally:
            browser.close()
    (run / "assets.json").write_text(json.dumps({"schema_version": 1, "assets": assets}, indent=2), encoding="utf-8")
