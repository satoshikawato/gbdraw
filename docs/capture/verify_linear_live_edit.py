#!/usr/bin/env python3
"""Verify Web reference editing checkpoints using the existing capture owners."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path

from Bio import SeqIO
from playwright.sync_api import expect, sync_playwright

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "docs" / "capture"))
from config import FIRST_LINEAR_LABEL_RULE_PATH, GENERATION_TIMEOUT_MS  # noqa: E402
from flows.web_capture import (  # noqa: E402
    assert_fixture_identity,
    fit_complete_linear_preview,
    generate_and_wait_for_result,
    open_browser_capture,
    set_feature_search_visible,
    wait_for_app_shell,
)
from flows.how_to.interactive_sessions import (  # noqa: E402
    _load_current_session,
    _save_current_session,
)
from web_server import CaptureWebServer  # noqa: E402

SOURCE_RECEIPT = Path(__file__).with_name("linear-live-edit-source-verification.json")
PUBLIC_IMAGE = ROOT / "docs/assets/web-app/linear-current-result.png"
EDIT_LABEL = "portal"


def checkpoint(page, name, output):
    """Observe state without a second request builder or artifact owner."""
    page.evaluate("""async () => { await window.Vue.nextTick();
      await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))); }""")
    page.wait_for_function(
        "() => !window.__GBDRAW_APP__.processing && !window.__GBDRAW_APP__.labelReflowPending && !window.__GBDRAW_APP__.featureExtractionPending && !window.__GBDRAW_HISTORY__.capturing.value && !window.__GBDRAW_HISTORY__.restoring.value",
        timeout=GENERATION_TIMEOUT_MS,
    )
    result = page.evaluate("""async () => {
      const { state: s } = await import('./js/state.js');
      const { buildConfigData, getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
      const { serializeCleanSvg } = await import('./js/services/svg-serialization.js');
      const root = s.svgContainer.value?.querySelector('svg');
      return {
        result: s.results.value[s.selectedResultIndex.value]?.content || '',
        mounted: root ? serializeCleanSvg(root) : '',
        request: getCommittedCanonicalRenderRequest(), draft: await buildConfigData(s.activeDrawing()),
        identity: s.extractedFeatures.value.map(f => ({id:f.id, biologicalFeatureId:f.biologicalFeatureId,
          record:f.record_id, start:f.start, end:f.end, strand:f.strand})),
        history: [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()],
        generation: s.resultGenerationKey.value,
        intent: window.__GBDRAW_HISTORY__.getCurrentIntent(),
        historyLabel: window.__GBDRAW_HISTORY__.undoLabel(), diagnostics: window.__GBDRAW_HISTORY__.getDiagnostics()
      };
    }""")
    result["resultSha256"] = hashlib.sha256(result["result"].encode()).hexdigest()
    (output / f"{name}.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        name,
        result["history"],
        result["resultSha256"],
        flush=True,
    )
    return result


def same_artifact(left, right):
    for key in ("result", "mounted", "request", "identity", "generation"):
        assert left[key] == right[key], key


def source_fixtures(output):
    receipts = json.loads(SOURCE_RECEIPT.read_text())
    for receipt, length in zip(receipts, (48502, 42925), strict=True):
        path = ROOT / receipt["mirrorPath"]
        assert receipt["equal"]
        assert_fixture_identity(
            path,
            expected_size=receipt["mirrorBytes"],
            expected_sha256=receipt["mirrorSha256"],
        )
        original = path.read_bytes()
        if path.suffix == ".gz":
            original = gzip.decompress(original)
        assert len(original) == receipt["sourceBytes"]
        assert hashlib.sha256(original).hexdigest() == receipt["sourceSha256"]
        inputs = output / "source-inputs"
        inputs.mkdir(exist_ok=True)
        path = inputs / f"{receipt['accession']}.gb"
        path.write_bytes(original)
        (record,) = SeqIO.parse(path, "genbank")
        assert record.id == receipt["accession"] and len(record) == length
        assert record.annotations["topology"] == "linear"
        yield path


def export_svg(page, output, name):
    with page.expect_download() as pending:
        page.get_by_role("button", name="SVG", exact=True).click()
    path = output / f"{name}.svg"
    pending.value.save_as(path)
    return path.read_text()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=Path("/tmp/gbdraw-linear-live-edit-evidence")
    )
    parser.add_argument(
        "--capture",
        action="store_true",
        help="Replace the reference page image after validation.",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    paths = tuple(source_fixtures(args.output))
    with CaptureWebServer() as server, sync_playwright() as pw:
        capture = open_browser_capture(pw.chromium, server.base_url)
        fresh_capture = None
        page = capture.page
        try:
            page.goto(server.base_url, wait_until="domcontentloaded")
            wait_for_app_shell(page)
            page.get_by_role("button", name="Linear", exact=True).click()
            lock = page.get_by_label("Lock Definition Column", exact=True)
            expect(lock).to_be_checked()
            titles = page.get_by_label("Titles and Record Labels", exact=True)
            titles.click()
            accession = page.get_by_label("Accession visibility", exact=True)
            length = page.get_by_label("Length / Coordinates visibility", exact=True)
            expect(accession).to_have_value("auto")
            expect(length).to_have_value("auto")
            page.get_by_label("Default font size", exact=True).fill("24")
            titles.click()
            page.get_by_role("button", name="Add sequence", exact=True).click()
            for i, path in enumerate(paths, 1):
                page.get_by_test_id(f"linear-genbank-{i}").set_input_files(path)
            page.get_by_label("Advanced comparison and layout", exact=True).click()
            rows = page.get_by_label("Arrange linear records in rows", exact=True)
            rows.check()
            row2 = page.get_by_label("Linear record row for sequence 2", exact=True)
            row2.fill("1")
            row2.press("Tab")
            # CSS finding: this existing status region has no accessible name/test ID.
            reason = page.locator("[data-linear-label-auto-layout]")
            expect(reason).to_contain_text(
                "throughout the diagram on the next successful Generate"
            )
            before = checkpoint(page, "01-shared-auto", args.output)
            page.get_by_role(
                "button", name="Record Labels: Accession", exact=True
            ).click()
            expect(accession).to_be_focused()
            after = checkpoint(page, "02-navigation", args.output)
            same_artifact(before, after)
            assert before["history"] == after["history"]
            expect(accession).to_have_value("auto")
            # Keep metadata visible in this finished comparison; Auto itself stays the default.
            accession.select_option("show")
            length.select_option("show")
            row2.fill("2")
            row2.press("Tab")
            lock.uncheck()
            expect(lock).not_to_be_checked()
            checkpoint(page, "03-lock-off", args.output)
            lock.check()
            expect(lock).to_be_checked()
            page.get_by_label("Labels", exact=True).click()
            page.get_by_label("Show Labels", exact=True).select_option("all")
            page.get_by_label("Label Font Size", exact=True).fill("24")
            page.get_by_label("Priority File (TSV)", exact=True).set_input_files(
                FIRST_LINEAR_LABEL_RULE_PATH
            )
            page.get_by_label("Labels", exact=True).click()
            page.get_by_role(
                "button", name="Run LOSAT for all adjacent pairs", exact=True
            ).click()
            page.get_by_label("Output Prefix", exact=True).fill("lambda_de3_live_edit")
            generate_and_wait_for_result(page)
            generated = checkpoint(page, "04-generated", args.output)
            assert len(generated["identity"]) >= 130
            assert (
                "NC_001416.1" in generated["mounted"]
                and "NC_042057.1" in generated["mounted"]
            )
            assert "comparison" in generated["mounted"].lower()
            # Use the same live Feature editor through its public controls.
            page.get_by_role("button", name="Editor", exact=True).click()
            page.get_by_label("Auto Reflow", exact=True).uncheck()
            # CSS/placeholder findings: list search and popup label lack associated labels.
            page.get_by_placeholder("Search by feature or annotation...").fill(
                "portal protein"
            )
            page.get_by_role("button", name="Edit", exact=True).first.click()
            page.get_by_placeholder("Edit label text", exact=True).fill(EDIT_LABEL)
            page.get_by_role("button", name="Apply Label", exact=True).click()
            page.get_by_role("button", name="Close feature popup", exact=True).click()
            page.get_by_title("Close editor", exact=True).click()
            live = checkpoint(page, "05-live", args.output)
            assert EDIT_LABEL in live["mounted"]
            assert live["request"] == generated["request"]
            assert live["identity"] == generated["identity"]
            assert live["history"][0] > generated["history"][0]
            page.get_by_label("Axis & Scale", exact=True).click()
            scale = page.get_by_label("Linear scale font size", exact=True)
            scale.fill("19")
            scale.press("Tab")
            edited = checkpoint(page, "06-pending", args.output)
            # The scale edit is draft only; the committed request and Result are unchanged.
            assert edited["draft"]["adv"]["scale_font_size"] == 19
            same_artifact(live, edited)
            fit_complete_linear_preview(page, target_zoom="40%")
            set_feature_search_visible(page, visible=False)
            camera = checkpoint(page, "06b-camera", args.output)
            same_artifact(edited, camera)
            assert edited["history"] == camera["history"]
            preview = page.get_by_role("region", name="Result Preview", exact=True)
            preview.scroll_into_view_if_needed()
            top = preview.bounding_box()
            # CSS finding: generated SVG has no independent accessible capture name.
            figure = preview.locator("svg").bounding_box()
            assert top and figure
            page.screenshot(
                path=str(args.output / "linear-current-result.png"),
                clip={
                    "x": top["x"],
                    "y": top["y"],
                    "width": top["width"],
                    "height": figure["y"] + figure["height"] - top["y"] + 24,
                },
            )
            saved_path = _save_current_session(
                page, args.output, title="lambda_de3_pending"
            )
            saved = json.loads(gzip.decompress(saved_path.read_bytes()))
            after_save = checkpoint(page, "07-saved", args.output)
            same_artifact(edited, after_save)
            # Naming the saved Session changes only its existing document metadata.
            assert after_save["history"] == [edited["history"][0] + 1, 0]
            named_intent = json.loads(json.dumps(edited["intent"]))
            named_intent["ui"]["title"] = "lambda_de3_pending"
            assert after_save["intent"] == named_intent
            for receipt in json.loads(SOURCE_RECEIPT.read_text()):
                import base64

                resource_hashes = [
                    hashlib.sha256(base64.b64decode(r["data"])).hexdigest()
                    for r in saved["resources"].values()
                    if "data" in r
                ]
                assert receipt["sourceSha256"] in resource_hashes
            fresh_capture = open_browser_capture(pw.chromium, server.base_url)
            fresh = fresh_capture.page
            fresh.goto(server.base_url, wait_until="domcontentloaded")
            wait_for_app_shell(fresh)
            _load_current_session(fresh, saved_path)
            loaded = checkpoint(fresh, "08-loaded", args.output)
            admitted = fresh.evaluate(
                """async text => (await import('./js/services/svg-sanitization.js')).sanitizeSvgContent(text)""",
                saved["results"][0]["content"],
            )
            assert loaded["result"] == admitted
            assert loaded["request"] == edited["request"]
            assert loaded["identity"] == edited["identity"]
            assert loaded["draft"]["adv"]["scale_font_size"] == 19
            exported = export_svg(fresh, args.output, "current-result")
            assert exported == loaded["mounted"]
            same_artifact(loaded, checkpoint(fresh, "09-exported", args.output))
            generate_and_wait_for_result(fresh)
            final = checkpoint(fresh, "10-regenerated", args.output)
            assert final["request"] != loaded["request"]
            assert EDIT_LABEL in final["mounted"]
            assert final["identity"] == loaded["identity"]
            fresh.get_by_role("button", name="Undo", exact=True).click()
            undone = checkpoint(fresh, "11-undo", args.output)
            assert (
                undone["result"] == loaded["result"]
                and undone["request"] == loaded["request"]
            )
            fresh.get_by_role("button", name="Redo", exact=True).click()
            redone = checkpoint(fresh, "12-redo", args.output)
            assert (
                redone["result"] == final["result"]
                and redone["request"] == final["request"]
            )
            assert (
                export_svg(fresh, args.output, "regenerated-result")
                == redone["mounted"]
            )
            capture.assert_clean()
            fresh_capture.assert_clean()
            if args.capture:
                PUBLIC_IMAGE.write_bytes(
                    (args.output / "linear-current-result.png").read_bytes()
                )
            print(
                "PASS: raw sources, Auto navigation, Lock, live edits, draft, Save/fresh Load, Export, Generate, Undo/Redo",
                flush=True,
            )
        finally:
            capture.close()
            if fresh_capture:
                fresh_capture.close()


if __name__ == "__main__":
    main()
