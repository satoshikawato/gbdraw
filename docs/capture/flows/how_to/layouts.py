"""Capture documentation journeys for multi-record layout controls."""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from Bio import SeqIO
from playwright.sync_api import BrowserType, Page, expect

from assertions.downloads import (
    EXPECTED_GUI_CIRCULAR_LAYOUT_SVG,
    assert_gui_circular_layout_svg_download,
)
from assertions.svg_semantics import (
    assert_gui_circular_layout_svg,
    assert_static_svg_safety,
    inspect_first_circular_svg,
    inspect_first_linear_svg,
    parse_translate_chain,
)
from config import (
    ACTION_TIMEOUT_MS,
    GUI_CIRCULAR_LAYOUT_COMBINED_FILENAME,
    GUI_CIRCULAR_LAYOUT_FIXTURES,
    GUI_CIRCULAR_LAYOUT_SCREENSHOT_NAMES,
    GUI_LINEAR_LAYOUT_SCREENSHOT_NAMES,
    GUI_LOSATN_DE3_FIXTURE_PATH,
    GUI_LOSATN_DE3_FIXTURE_SHA256,
    GUI_LOSATN_DE3_FIXTURE_SIZE,
    FIRST_LINEAR_FIXTURE_PATH,
    FIRST_LINEAR_FIXTURE_SHA256,
    FIRST_LINEAR_FIXTURE_SIZE,
    FIRST_LINEAR_LABEL_RULE_PATH,
    FIRST_LINEAR_LABEL_RULE_SHA256,
    FIRST_LINEAR_LABEL_RULE_SIZE,
    WORKER_READY_TIMEOUT_MS,
)
from flows.web_capture import (
    assert_fixture_identity,
    assert_output_paths,
    capture_screenshot,
    fit_complete_linear_preview,
    generate_and_inspect,
    open_browser_capture,
    set_feature_search_visible,
    wait_for_worker,
)


@dataclass(frozen=True)
class CaptureResult:
    screenshot_bytes: dict[str, int]
    final_svg_semantics: dict[str, Any]
    download: dict[str, Any]
    source_records: tuple[dict[str, Any], ...]


EXPECTED_LINEAR_LAYOUT_SVG = "linear_regions_orientation.svg"


@dataclass(frozen=True)
class LinearLayoutCaptureResult:
    screenshot_bytes: dict[str, int]
    final_svg_semantics: dict[str, Any]
    download: dict[str, Any]
    source_records: tuple[dict[str, Any], ...]
    rendered_state: dict[str, Any]


def _write_complete_mitochondrial_container(path: Path) -> tuple[dict[str, Any], ...]:
    """Place four unchanged complete GenBank records in one upload container."""

    source_records = []
    source_bytes = []
    for fixture_path, size, digest, record_id, length, organism in (
        GUI_CIRCULAR_LAYOUT_FIXTURES
    ):
        assert_fixture_identity(
            fixture_path,
            expected_size=size,
            expected_sha256=digest,
        )
        records = list(SeqIO.parse(fixture_path, "genbank"))
        if len(records) != 1:
            raise AssertionError(f"Expected one record in {fixture_path.name}")
        record = records[0]
        if record.id != record_id or len(record) != length:
            raise AssertionError(
                f"Unexpected complete-record identity in {fixture_path.name}: "
                f"{record.id} ({len(record)} bp)"
            )
        if str(record.annotations.get("topology", "")).lower() != "circular":
            raise AssertionError(f"{record_id} is not annotated as circular")
        if record.annotations.get("organism") != organism:
            raise AssertionError(
                f"Unexpected organism for {record_id}: {record.annotations.get('organism')}"
            )
        if "complete" not in record.description.lower():
            raise AssertionError(f"{record_id} is not described as a complete genome")
        source_records.append(record)
        source_bytes.append(fixture_path.read_bytes())

    combined_bytes = b"".join(source_bytes)
    path.write_bytes(combined_bytes)
    if path.read_bytes() != combined_bytes:
        raise AssertionError("The multi-record upload container changed source bytes")

    combined_records = list(SeqIO.parse(path, "genbank"))
    if len(combined_records) != len(source_records):
        raise AssertionError("The upload container did not preserve all four records")
    report = []
    for source, combined in zip(source_records, combined_records, strict=True):
        if (
            combined.id != source.id
            or str(combined.seq) != str(source.seq)
            or len(combined) != len(source)
            or combined.annotations.get("topology") != source.annotations.get("topology")
        ):
            raise AssertionError(
                f"The upload container changed complete record {source.id}"
            )
        report.append(
            {
                "record_id": combined.id,
                "length": len(combined),
                "topology": combined.annotations.get("topology"),
                "organism": combined.annotations.get("organism"),
            }
        )
    return tuple(report)


def _fit_complete_grid_preview(page: Page, target_zoom: str = "30%") -> None:
    """Fit the complete 2 by 2 canvas in the public result viewport."""

    reset_zoom = page.get_by_role("button", name="Reset zoom", exact=True)
    zoom_out = page.get_by_role("button", name="Zoom out", exact=True)
    for _ in range(10):
        if target_zoom in reset_zoom.inner_text():
            break
        zoom_out.click()
    else:
        raise AssertionError(
            f"Could not reach the documented Circular preview zoom: {target_zoom}"
        )

    result_region = page.get_by_role("region", name="Result Preview", exact=True)
    box = result_region.bounding_box()
    if box is None:
        raise AssertionError("Could not resolve the Circular result preview bounds")
    page.mouse.move(
        box["x"] + (box["width"] * 0.82),
        box["y"] + (box["height"] * 0.42),
    )
    page.mouse.down()
    page.mouse.move(
        box["x"] + (box["width"] * 0.225),
        box["y"] + (box["height"] * 0.50),
        steps=12,
    )
    page.mouse.up()
    page.evaluate("() => window.getSelection()?.removeAllRanges()")
    selection_range_count = page.evaluate(
        "() => window.getSelection()?.rangeCount ?? 0"
    )
    if selection_range_count != 0:
        raise AssertionError("Circular preview retained a text selection after panning")


def capture_gui_circular_layout(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> CaptureResult:
    """Run H-GUI-02 using four complete, naturally circular RefSeq records."""

    assert_output_paths(
        output_paths,
        GUI_CIRCULAR_LAYOUT_SCREENSHOT_NAMES,
        "H-GUI-02",
    )
    download_dir.mkdir(parents=True, exist_ok=True)
    upload_path = download_dir / GUI_CIRCULAR_LAYOUT_COMBINED_FILENAME
    source_records = _write_complete_mitochondrial_container(upload_path)
    assert_fixture_identity(
        FIRST_LINEAR_LABEL_RULE_PATH,
        expected_size=FIRST_LINEAR_LABEL_RULE_SIZE,
        expected_sha256=FIRST_LINEAR_LABEL_RULE_SHA256,
    )

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}

    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)

        circular = page.get_by_role("button", name="Circular", exact=True)
        circular.click()
        expect(circular).to_have_attribute("aria-pressed", "true")
        genbank = page.get_by_role("radio", name="GenBank", exact=True)
        genbank.check()
        expect(genbank).to_be_checked()
        page.get_by_label("GenBank/DDBJ File", exact=True).set_input_files(upload_path)
        expect(
            page.get_by_role(
                "group", name="GenBank/DDBJ File selection", exact=True
            )
        ).to_contain_text(GUI_CIRCULAR_LAYOUT_COMBINED_FILENAME)

        prefix = page.get_by_label("Output Prefix", exact=True)
        prefix.fill("multi_record_circular")
        expect(prefix).to_have_value("multi_record_circular")
        multi_record = page.get_by_label("Multi-Record Canvas", exact=True)
        multi_record.check()
        expect(multi_record).to_be_checked()

        expected_rows = (
            ("#1", "NC_012920.1", "1"),
            ("#2", "NC_002333.2", "1"),
            ("#3", "NC_024511.2", "2"),
            ("#4", "NC_001328.1", "2"),
        )
        row_controls = []
        for selector, record_id, row in expected_rows:
            control = page.get_by_label(
                f"Row for {selector} ({record_id})",
                exact=True,
            )
            expect(control).to_be_visible(timeout=WORKER_READY_TIMEOUT_MS)
            control.select_option(row)
            expect(control).to_have_value(row)
            row_controls.append(control)

        size_mode = page.get_by_label("Record Size Mode", exact=True)
        size_mode.select_option("equal")
        expect(size_mode).to_have_value("equal")
        min_radius = page.get_by_label("Min Radius Ratio", exact=True)
        min_radius.fill("0.75")
        expect(min_radius).to_have_value("0.75")
        column_gap = page.get_by_label("Column Gap Ratio", exact=True)
        column_gap.fill("0.40")
        expect(column_gap).to_have_value("0.40")
        row_gap = page.get_by_label("Row Gap Ratio", exact=True)
        row_gap.fill("0.08")
        expect(row_gap).to_have_value("0.08")

        title_and_legend = page.get_by_label("Title & Legend", exact=True)
        title_and_legend.click()
        plot_title = page.get_by_label("Plot Title", exact=True)
        plot_title.fill("Complete metazoan mitochondrial genomes")
        expect(plot_title).to_have_value("Complete metazoan mitochondrial genomes")
        title_position = page.get_by_label("Plot Title Position", exact=True)
        title_position.select_option("top")
        expect(title_position).to_have_value("top")
        definition_font_size = page.get_by_label("Definition Font Size", exact=True)
        definition_font_size.fill("20")
        expect(definition_font_size).to_have_value("20")
        keep_definitions = page.get_by_label(
            "Keep Full Definition with Plot Title",
            exact=True,
        )
        keep_definitions.check()
        expect(keep_definitions).to_be_checked()
        title_and_legend.click()

        labels = page.get_by_label("Labels", exact=True)
        labels.click()
        label_mode = page.get_by_label("Label Mode", exact=True)
        label_mode.select_option("out")
        expect(label_mode).to_have_value("out")
        page.get_by_label("Priority File (TSV)", exact=True).set_input_files(
            FIRST_LINEAR_LABEL_RULE_PATH
        )
        label_font_size = page.get_by_label("Label Font Size", exact=True)
        label_font_size.fill("16")
        expect(label_font_size).to_have_value("16")
        labels.click()

        size_mode.scroll_into_view_if_needed()
        expect(row_controls[-1]).to_be_visible()
        screenshot_bytes["grid-settings.png"] = capture_screenshot(
            page,
            output_paths["grid-settings.png"],
            "Circular",
        )

        final_report = generate_and_inspect(
            page,
            inspect_first_circular_svg,
            assert_gui_circular_layout_svg,
        )
        set_feature_search_visible(page, visible=False)
        _fit_complete_grid_preview(page)
        screenshot_bytes["grid-result.png"] = capture_screenshot(
            page,
            output_paths["grid-result.png"],
            "Circular",
        )

        svg_button = page.get_by_role("button", name="SVG", exact=True)
        with page.expect_download(timeout=ACTION_TIMEOUT_MS) as download_info:
            svg_button.click()
        download = download_info.value
        if download.failure() is not None:
            raise AssertionError(f"SVG download failed: {download.failure()}")
        if download.suggested_filename != EXPECTED_GUI_CIRCULAR_LAYOUT_SVG:
            raise AssertionError(
                "Unexpected SVG download name: "
                f"expected {EXPECTED_GUI_CIRCULAR_LAYOUT_SVG}, "
                f"found {download.suggested_filename}"
            )
        downloaded_svg = download_dir / download.suggested_filename
        download.save_as(downloaded_svg)
        download_report = assert_gui_circular_layout_svg_download(downloaded_svg)

        capture.assert_clean()
    finally:
        capture.close()

    return CaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        source_records=source_records,
    )


def _complete_linear_source_records() -> tuple[dict[str, Any], ...]:
    records = []
    for path, size, digest, record_id, length in (
        (
            FIRST_LINEAR_FIXTURE_PATH,
            FIRST_LINEAR_FIXTURE_SIZE,
            FIRST_LINEAR_FIXTURE_SHA256,
            "NC_001416.1",
            48_502,
        ),
        (
            GUI_LOSATN_DE3_FIXTURE_PATH,
            GUI_LOSATN_DE3_FIXTURE_SIZE,
            GUI_LOSATN_DE3_FIXTURE_SHA256,
            "NC_042057.1",
            42_925,
        ),
    ):
        assert_fixture_identity(path, expected_size=size, expected_sha256=digest)
        parsed = list(SeqIO.parse(path, "genbank"))
        if len(parsed) != 1:
            raise AssertionError(f"Expected one complete record in {path.name}")
        record = parsed[0]
        if record.id != record_id or len(record) != length:
            raise AssertionError(
                f"Unexpected source record in {path.name}: {record.id} ({len(record)} bp)"
            )
        if str(record.annotations.get("topology", "")).lower() != "linear":
            raise AssertionError(f"{record.id} is not annotated as linear")
        if "complete" not in record.description.lower():
            raise AssertionError(f"{record.id} is not described as complete")
        records.append(
            {
                "record_id": record.id,
                "length": len(record),
                "topology": record.annotations.get("topology"),
                "description": record.description,
            }
        )
    return tuple(records)


def _inspect_complete_linear_layout(result_region: Any) -> dict[str, Any]:
    return inspect_first_linear_svg(result_region)


def _assert_complete_linear_layout(report: dict[str, Any]) -> None:
    if report.get("svgCount") != 1:
        raise AssertionError(f"Expected one Linear SVG, found {report.get('svgCount')}")
    if set(report.get("recordIds", [])) != {"NC_001416.1", "NC_042057.1"}:
        raise AssertionError("The result does not contain the complete Lambda and DE3 records")
    if report.get("featureElementCount") != 130:
        raise AssertionError(
            "Expected 73 Lambda and 57 DE3 CDS elements, found "
            f"{report.get('featureElementCount')}"
        )
    texts = set(report.get("texts", []))
    for expected in (
        "NC_001416.1",
        "1-48502",
        "NC_042057.1",
        "42925-1",
        "Enterobacteria phage lambda",
        "Enterobacteria phage DE3",
    ):
        if expected not in texts:
            raise AssertionError(
                f"Linear layout SVG is missing {expected!r}: {sorted(texts)!r}"
            )
    record_placements = {
        str(item.get("recordId")): item
        for item in report.get("recordPlacements", [])
        if item.get("recordId")
    }
    if set(record_placements) != {"NC_001416.1", "NC_042057.1"}:
        raise AssertionError(
            f"Unexpected Linear record placements: {record_placements!r}"
        )
    lambda_bounds = record_placements["NC_001416.1"].get("bounds") or {}
    de3_bounds = record_placements["NC_042057.1"].get("bounds") or {}
    if float(lambda_bounds.get("bottom", 0)) >= float(de3_bounds.get("y", 0)):
        raise AssertionError(
            "The two documented Linear rows overlap: "
            f"{record_placements!r}"
        )
    if len(report.get("coordinateTicks", [])) < 4:
        raise AssertionError("The Linear ruler did not emit enough coordinate labels")
    if report.get("matches"):
        raise AssertionError("The layout-only diagram must not include comparison links")
    assert_static_svg_safety(report)


def _assert_linear_layout_download(path: Path) -> dict[str, Any]:
    contents = path.read_bytes()
    if len(contents) < 20_000:
        raise AssertionError(f"Downloaded Linear SVG is unexpectedly small: {len(contents)}")
    root = ET.fromstring(contents)
    if root.tag.rsplit("}", 1)[-1] != "svg":
        raise AssertionError("Downloaded Linear output is not an SVG")
    elements = list(root.iter())
    if any(element.tag.rsplit("}", 1)[-1] == "script" for element in elements):
        raise AssertionError("Static Linear SVG contains a script")
    record_ids = {
        element.attrib.get("data-gbdraw-record-id", "") for element in elements
    }
    if not {"NC_001416.1", "NC_042057.1"} <= record_ids:
        raise AssertionError("Downloaded SVG lost a complete input record")
    record_groups = {
        element.attrib.get("data-gbdraw-record-id", ""): element
        for element in elements
        if element.tag.rsplit("}", 1)[-1] == "g"
        and element.attrib.get("id", "").startswith("record_group_")
        and element.attrib.get("data-gbdraw-record-id")
    }
    translations = {}
    for record_id, group in record_groups.items():
        translation = parse_translate_chain(group.attrib.get("transform", ""))
        if translation is not None:
            translations[record_id] = translation
    if set(translations) != {"NC_001416.1", "NC_042057.1"}:
        raise AssertionError(f"Downloaded SVG record groups changed: {translations!r}")
    if translations["NC_001416.1"][1] >= translations["NC_042057.1"][1]:
        raise AssertionError(f"Downloaded SVG row order changed: {translations!r}")
    text = " ".join("".join(element.itertext()) for element in elements)
    for expected in ("1-48502", "42925-1"):
        if expected not in text:
            raise AssertionError(f"Downloaded SVG is missing {expected}")
    return {
        "filename": path.name,
        "bytes": len(contents),
        "record_ids": sorted(record_ids - {""}),
        "record_translations": translations,
    }


def capture_gui_linear_layout(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> LinearLayoutCaptureResult:
    """Run H-GUI-03 with two complete linear phage genomes."""

    assert_output_paths(
        output_paths,
        GUI_LINEAR_LAYOUT_SCREENSHOT_NAMES,
        "H-GUI-03",
    )
    source_records = _complete_linear_source_records()
    download_dir.mkdir(parents=True, exist_ok=True)

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}

    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)

        linear = page.get_by_role("button", name="Linear", exact=True)
        linear.click()
        expect(linear).to_have_attribute("aria-pressed", "true")
        page.get_by_role("radio", name="GenBank", exact=True).check()
        page.get_by_role("radio", name="No comparison", exact=True).first.check()
        page.get_by_role("button", name="Add sequence", exact=True).click()
        page.get_by_test_id("linear-genbank-1").set_input_files(
            FIRST_LINEAR_FIXTURE_PATH
        )
        page.get_by_test_id("linear-genbank-2").set_input_files(
            GUI_LOSATN_DE3_FIXTURE_PATH
        )

        definitions = (
            (1, "Enterobacteria phage lambda", 48_502, False),
            (2, "Enterobacteria phage DE3", 42_925, True),
        )
        for index, definition, length, reverse in definitions:
            definition_input = page.get_by_label(
                f"Definition for sequence {index}", exact=True
            )
            definition_input.fill(definition)
            expect(definition_input).to_have_value(definition)
            start = page.get_by_label(
                f"Region start for sequence {index}", exact=True
            )
            end = page.get_by_label(f"Region end for sequence {index}", exact=True)
            start.fill("1")
            end.fill(str(length))
            expect(start).to_have_value("1")
            expect(end).to_have_value(str(length))
            reverse_control = page.get_by_label(
                f"Reverse complement for sequence {index}", exact=True
            )
            if reverse:
                reverse_control.check()
                expect(reverse_control).to_be_checked()
            else:
                expect(reverse_control).not_to_be_checked()

        arrange_rows = page.get_by_label(
            "Arrange linear records in rows", exact=True
        )
        arrange_rows.check()
        expect(arrange_rows).to_be_checked()
        for index, row in ((1, "1"), (2, "2")):
            row_control = page.get_by_label(
                f"Linear record row for sequence {index}", exact=True
            )
            row_control.fill(row)
            row_control.press("Tab")
            expect(row_control).to_have_value(row)
        record_gap = page.get_by_label("Linear record gap", exact=True)
        record_gap.fill("24")
        expect(record_gap).to_have_value("24")

        track_layout = page.get_by_label("Track Layout", exact=True)
        track_layout.select_option("above")
        expect(track_layout).to_have_value("above")
        page.get_by_label("Axis & Scale", exact=True).click()
        show_scale = page.get_by_label(
            "Show Coordinate Scale (Linear)", exact=True
        )
        show_scale.check()
        expect(show_scale).to_be_checked()
        scale_style = page.get_by_label("Linear scale style", exact=True)
        scale_style.select_option("ruler")
        expect(scale_style).to_have_value("ruler")
        ruler_on_axis = page.get_by_label("Ruler on Axis", exact=True)
        ruler_on_axis.check()
        expect(ruler_on_axis).to_be_checked()
        prefix = page.get_by_label("Output Prefix", exact=True)
        prefix.fill("linear_regions_orientation")
        expect(prefix).to_have_value("linear_regions_orientation")

        arrange_rows.scroll_into_view_if_needed()
        screenshot_bytes["record-layout.png"] = capture_screenshot(
            page, output_paths["record-layout.png"], "Linear"
        )

        final_report = generate_and_inspect(
            page,
            _inspect_complete_linear_layout,
            _assert_complete_linear_layout,
        )
        rendered_state = page.evaluate(
            r"""
            () => ({
              sequences: (window.__GBDRAW_APP__?.linearSeqs || []).map((seq) => ({
                definition: String(seq?.definition || ''),
                start: Number(seq?.region_start),
                end: Number(seq?.region_end),
                reverse: Boolean(seq?.region_reverse)
              })),
              layoutEnabled: Boolean(window.__GBDRAW_APP__?.linearRecordLayoutEnabled),
              gap: Number(window.__GBDRAW_APP__?.linearRecordGap),
              scaleStyle: String(window.__GBDRAW_APP__?.form?.scale_style || ''),
              rulerOnAxis: Boolean(window.__GBDRAW_APP__?.form?.linear_ruler_on_axis)
            })
            """
        )
        rendered_state["rows"] = [
            int(
                page.get_by_label(
                    f"Linear record row for sequence {index}", exact=True
                ).input_value()
            )
            for index in (1, 2)
        ]
        expected_sequences = [
            {
                "definition": "Enterobacteria phage lambda",
                "start": 1,
                "end": 48_502,
                "reverse": False,
            },
            {
                "definition": "Enterobacteria phage DE3",
                "start": 1,
                "end": 42_925,
                "reverse": True,
            },
        ]
        if rendered_state.get("sequences") != expected_sequences:
            raise AssertionError(
                f"The generated complete-record orientation state changed: {rendered_state!r}"
            )
        if rendered_state.get("gap") != 24:
            raise AssertionError(f"Unexpected record gap: {rendered_state!r}")
        if not rendered_state.get("layoutEnabled") or rendered_state.get("rows") != [
            1,
            2,
        ]:
            raise AssertionError(f"The documented row controls changed: {rendered_state!r}")
        if rendered_state.get("scaleStyle") != "ruler" or not rendered_state.get(
            "rulerOnAxis"
        ):
            raise AssertionError(f"The documented ruler is not active: {rendered_state!r}")

        fit_complete_linear_preview(page, target_zoom="30%", pan_left=True)
        screenshot_bytes["orientation-result.png"] = capture_screenshot(
            page, output_paths["orientation-result.png"], "Linear"
        )

        svg_button = page.get_by_role("button", name="SVG", exact=True)
        with page.expect_download(timeout=ACTION_TIMEOUT_MS) as download_info:
            svg_button.click()
        download = download_info.value
        if download.failure() is not None:
            raise AssertionError(f"Linear layout SVG download failed: {download.failure()}")
        if download.suggested_filename != EXPECTED_LINEAR_LAYOUT_SVG:
            raise AssertionError(
                "Unexpected Linear layout SVG filename: "
                f"{download.suggested_filename}"
            )
        downloaded_svg = download_dir / download.suggested_filename
        download.save_as(downloaded_svg)
        download_report = _assert_linear_layout_download(downloaded_svg)

        capture.assert_clean()
    finally:
        capture.close()

    return LinearLayoutCaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        source_records=source_records,
        rendered_state=rendered_state,
    )
