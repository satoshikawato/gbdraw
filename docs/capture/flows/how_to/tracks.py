"""Capture quantitative and annotation-track GUI how-to journeys."""

from __future__ import annotations

import csv
import hashlib
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Any, Mapping

from Bio import SeqIO
from playwright.sync_api import BrowserType, Page, expect

from assertions.svg_semantics import inspect_first_circular_svg
from config import ACTION_TIMEOUT_MS
from flows.web_capture import (
    assert_fixture_identity,
    assert_output_paths,
    capture_screenshot,
    generate_and_inspect,
    open_browser_capture,
    set_feature_search_visible,
    wait_for_app_shell,
)


REPO_ROOT = Path(__file__).resolve().parents[4]
TUTORIAL_DATA_ROOT = REPO_ROOT / "gbdraw" / "web" / "tutorial-data"

GUI_QUANTITATIVE_SCENARIO_ID = "H-GUI-09"
GUI_QUANTITATIVE_SCREENSHOT_NAMES = (
    "track-settings.png",
    "track-result.png",
)
GUI_QUANTITATIVE_GENBANK_PATH = (
    TUTORIAL_DATA_ROOT / "depth-1kb" / "AP027133.gb"
)
GUI_QUANTITATIVE_GENBANK_SIZE = 1_344_094
GUI_QUANTITATIVE_GENBANK_SHA256 = (
    "913af50dd9d37cc2107be5e46484b885c5d586fb414b4b501380fc8f17a659d6"
)
GUI_QUANTITATIVE_DEPTH_PATH = (
    TUTORIAL_DATA_ROOT
    / "depth-1kb"
    / "AP027133.DRR394922.depth-1kb.tsv"
)
GUI_QUANTITATIVE_DEPTH_SIZE = 16_913
GUI_QUANTITATIVE_DEPTH_SHA256 = (
    "6f57cfd89a165ad97a162aa2f0b1f3b3ad21fb5638f4f9ac5cbd069badd6aab7"
)
EXPECTED_GUI_QUANTITATIVE_SVG = "quantitative_tracks.svg"

GUI_ANNOTATION_SCENARIO_ID = "H-GUI-10"
GUI_ANNOTATION_SCREENSHOT_NAMES = (
    "slot-settings.png",
    "annotation-result.png",
)
GUI_ANNOTATION_GENBANK_PATH = (
    TUTORIAL_DATA_ROOT / "tobacco-plastome-regions" / "NC_001879.gbk"
)
GUI_ANNOTATION_GENBANK_SIZE = 331_860
GUI_ANNOTATION_GENBANK_SHA256 = (
    "25c5b39fd25d702c0a390fe5e7480eda0ccc1e4d6d7c388445b4686049412a24"
)
GUI_ANNOTATION_TABLE_PATH = (
    TUTORIAL_DATA_ROOT
    / "tobacco-plastome-regions"
    / "nicotiana-tabacum-regions.tsv"
)
GUI_ANNOTATION_TABLE_SIZE = 594
GUI_ANNOTATION_TABLE_SHA256 = (
    "3a85aed5145c88f93b4478d1901fab53714b9d47afc754d32cc9e5c0b8412b88"
)
EXPECTED_GUI_ANNOTATION_SVG = "region_annotations_and_slots.svg"

EXPECTED_QUANTITATIVE_SLOT_ORDER = (
    ("features", "features"),
    ("ticks", "ticks"),
    ("depth", "depth"),
    ("gc_content", "dinucleotide_content"),
    ("gc_skew", "dinucleotide_skew"),
)
EXPECTED_ANNOTATION_SLOT_ORDER = (
    ("annotations", "annotations"),
    ("features", "features"),
    ("ticks", "ticks"),
    ("gc_content", "dinucleotide_content"),
    ("gc_skew", "dinucleotide_skew"),
)


@dataclass(frozen=True)
class CaptureResult:
    screenshot_bytes: dict[str, int]
    final_svg_semantics: dict[str, Any]
    download: dict[str, Any]
    fixture_report: dict[str, Any]
    track_slots: tuple[dict[str, Any], ...]


def _local_name(name: str) -> str:
    return name.rsplit("}", 1)[-1].lower()


def _normalized_text(element: ET.Element) -> str:
    return " ".join("".join(element.itertext()).split())


def _validate_complete_record(
    path: Path,
    *,
    expected_size: int,
    expected_sha256: str,
    expected_id: str,
    expected_length: int,
) -> dict[str, Any]:
    assert_fixture_identity(
        path,
        expected_size=expected_size,
        expected_sha256=expected_sha256,
    )
    records = list(SeqIO.parse(path, "genbank"))
    if len(records) != 1:
        raise AssertionError(f"Expected one complete record in {path.name}")
    record = records[0]
    if record.id != expected_id or len(record) != expected_length:
        raise AssertionError(
            f"Unexpected record identity in {path.name}: {record.id} ({len(record)} bp)"
        )
    topology = str(record.annotations.get("topology", "")).lower()
    if topology != "circular":
        raise AssertionError(f"{expected_id} is not annotated as circular")
    if "complete" not in record.description.lower():
        raise AssertionError(f"{expected_id} is not described as complete")
    return {
        "recordId": record.id,
        "recordLength": len(record),
        "topology": topology,
        "description": record.description,
    }


def _validate_depth_fixture() -> dict[str, Any]:
    assert_fixture_identity(
        GUI_QUANTITATIVE_DEPTH_PATH,
        expected_size=GUI_QUANTITATIVE_DEPTH_SIZE,
        expected_sha256=GUI_QUANTITATIVE_DEPTH_SHA256,
    )
    with GUI_QUANTITATIVE_DEPTH_PATH.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["reference_name", "position", "depth"]:
            raise AssertionError("The 1 kbp depth fixture header changed")
        rows = list(reader)
    if len(rows) != 607:
        raise AssertionError(f"Expected 607 depth bins, found {len(rows)}")
    positions = [int(row["position"]) for row in rows]
    expected_positions = [1 + (index * 1_000) for index in range(607)]
    if positions != expected_positions:
        raise AssertionError("Depth bins are not consecutive 1 kbp starts")
    if {row["reference_name"] for row in rows} != {"AP027133.1"}:
        raise AssertionError("Depth rows do not all target AP027133.1")
    values = [Decimal(row["depth"]) for row in rows]
    expected_edge_values = (Decimal("21.901000"), Decimal("15.938144"))
    if (values[0], values[-1]) != expected_edge_values:
        raise AssertionError("The frozen first or final depth-bin mean changed")
    if min(values) != Decimal("12.446000") or max(values) != Decimal("74.546000"):
        raise AssertionError("The frozen depth range changed")
    return {
        "binCount": len(rows),
        "binSizeBp": 1_000,
        "firstPosition": positions[0],
        "lastPosition": positions[-1],
        "minimum": str(min(values)),
        "maximum": str(max(values)),
        "firstValue": str(values[0]),
        "lastValue": str(values[-1]),
    }


def _validate_annotation_fixture() -> dict[str, Any]:
    assert_fixture_identity(
        GUI_ANNOTATION_TABLE_PATH,
        expected_size=GUI_ANNOTATION_TABLE_SIZE,
        expected_sha256=GUI_ANNOTATION_TABLE_SHA256,
    )
    with GUI_ANNOTATION_TABLE_PATH.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = (
        ("lsc", "bracket", 1, 86_686, "LSC", "0"),
        ("irb", "bracket", 86_687, 112_029, "IRb", "0"),
        ("ssc", "bracket", 112_030, 130_600, "SSC", "0"),
        ("ira", "bracket", 130_601, 155_943, "IRa", "0"),
    )
    actual = tuple(
        (
            row["id"],
            row["mark"],
            int(row["start"]),
            int(row["end"]),
            row["label"],
            row["lane"],
        )
        for row in rows
    )
    if actual != expected:
        raise AssertionError("The frozen four-region annotation table changed")
    if {row["set_id"] for row in rows} != {"plastome_regions"}:
        raise AssertionError("The tobacco annotation set id changed")
    if {row["record"] for row in rows} != {"NC_001879.2"}:
        raise AssertionError("The tobacco annotations target a different record")
    return {
        "setId": "plastome_regions",
        "annotationCount": len(rows),
        "annotationIds": [row[0] for row in actual],
        "labels": [row[4] for row in actual],
        "lanes": [int(row[5]) for row in actual],
    }


_TRACK_SVG_INSPECTION_SCRIPT = r"""
root => {
  const svg = root.getElementsByTagName('svg')[0];
  if (!svg) return { svgCount: 0 };
  const elements = [svg, ...Array.from(svg.getElementsByTagName('*'))];
  const normalizedText = (node) => String(node?.textContent || '')
    .replace(/\s+/g, ' ')
    .trim();
  const texts = Array.from(svg.getElementsByTagName('text'))
    .map(normalizedText)
    .filter(Boolean);
  const slots = elements
    .filter((element) => Boolean(element.getAttribute('data-gbdraw-slot-renderer')))
    .map((element) => {
      const descendants = [element, ...Array.from(element.getElementsByTagName('*'))];
      const slotTexts = descendants
        .filter((item) => String(item.tagName || '').toLowerCase() === 'text')
        .map(normalizedText)
        .filter(Boolean);
      return {
        id: String(element.getAttribute('id') || ''),
        slotId: String(element.getAttribute('data-gbdraw-slot-id') || ''),
        renderer: String(element.getAttribute('data-gbdraw-slot-renderer') || ''),
        pathCount: descendants.filter((item) => String(item.tagName || '').toLowerCase() === 'path').length,
        tickTexts: slotTexts.filter((text) => /^(?:\d+(?:\.\d+)?x|\d+(?:\.\d+)?%)$/.test(text)),
        texts: slotTexts
      };
    });
  const annotations = elements
    .filter((element) => Boolean(element.getAttribute('data-gbdraw-annotation-id')))
    .map((element) => ({
      id: String(element.getAttribute('data-gbdraw-annotation-id') || ''),
      setId: String(element.getAttribute('data-gbdraw-annotation-set-id') || ''),
      trackId: String(element.getAttribute('data-gbdraw-annotation-track-id') || ''),
      recordId: String(element.getAttribute('data-gbdraw-record-id') || ''),
      mark: String(element.getAttribute('data-gbdraw-annotation-mark') || ''),
      label: String(element.getAttribute('data-gbdraw-annotation-label') || '')
    }));
  const eventAttributes = [];
  const unsafeHrefs = [];
  for (const element of elements) {
    for (const attribute of Array.from(element.attributes || [])) {
      const name = String(attribute.name || '').toLowerCase();
      const value = String(attribute.value || '').trim();
      if (name.startsWith('on')) eventAttributes.push(`${name}=${value}`);
      if ((name === 'href' || name === 'xlink:href') && /^(?:javascript:|https?:|\/\/|data:text\/html)/i.test(value)) {
        unsafeHrefs.push(value);
      }
    }
  }
  return {
    svgCount: 1,
    recordIds: elements
      .map((element) => String(element.getAttribute('data-gbdraw-record-id') || ''))
      .filter(Boolean),
    logicalFeatureCount: Array.isArray(window.__GBDRAW_APP__?.extractedFeatures)
      ? window.__GBDRAW_APP__.extractedFeatures.length
      : null,
    featureElementCount: elements.filter((element) => Boolean(element.getAttribute('data-gbdraw-feature-id'))).length,
    texts,
    slots,
    annotations,
    scriptCount: svg.getElementsByTagName('script').length,
    eventAttributes,
    unsafeHrefs
  };
}
"""


def _inspect_tracks_svg(result_region: Any) -> dict[str, Any]:
    report = result_region.evaluate(_TRACK_SVG_INSPECTION_SCRIPT)
    base = inspect_first_circular_svg(result_region)
    report["interactiveSvg"] = base["interactiveSvg"]
    return report


def _inspect_tracks_svg_file(path: Path) -> dict[str, Any]:
    root = ET.parse(path).getroot()
    elements = list(root.iter())
    texts = [
        _normalized_text(element)
        for element in elements
        if _local_name(element.tag) == "text" and _normalized_text(element)
    ]
    slots: list[dict[str, Any]] = []
    for element in elements:
        renderer = element.attrib.get("data-gbdraw-slot-renderer")
        if not renderer:
            continue
        descendants = list(element.iter())
        slot_texts = [
            _normalized_text(descendant)
            for descendant in descendants
            if _local_name(descendant.tag) == "text" and _normalized_text(descendant)
        ]
        slots.append(
            {
                "id": element.attrib.get("id", ""),
                "slotId": element.attrib.get("data-gbdraw-slot-id", ""),
                "renderer": renderer,
                "pathCount": sum(
                    1 for descendant in descendants if _local_name(descendant.tag) == "path"
                ),
                "tickTexts": [
                    text
                    for text in slot_texts
                    if re.fullmatch(r"(?:\d+(?:\.\d+)?x|\d+(?:\.\d+)?%)", text)
                ],
                "texts": slot_texts,
            }
        )
    annotations = [
        {
            "id": element.attrib.get("data-gbdraw-annotation-id", ""),
            "setId": element.attrib.get("data-gbdraw-annotation-set-id", ""),
            "trackId": element.attrib.get("data-gbdraw-annotation-track-id", ""),
            "recordId": element.attrib.get("data-gbdraw-record-id", ""),
            "mark": element.attrib.get("data-gbdraw-annotation-mark", ""),
            "label": element.attrib.get("data-gbdraw-annotation-label", ""),
        }
        for element in elements
        if element.attrib.get("data-gbdraw-annotation-id")
    ]
    event_attributes: list[str] = []
    unsafe_hrefs: list[str] = []
    for element in elements:
        for raw_name, value in element.attrib.items():
            name = _local_name(raw_name)
            if name.startswith("on"):
                event_attributes.append(f"{name}={value}")
            if name == "href" and re.match(
                r"^(?:javascript:|https?:|//|data:text/html)", value, re.I
            ):
                unsafe_hrefs.append(value)
    return {
        "svgCount": 1 if _local_name(root.tag) == "svg" else 0,
        "recordIds": [
            element.attrib["data-gbdraw-record-id"]
            for element in elements
            if element.attrib.get("data-gbdraw-record-id")
        ],
        "logicalFeatureCount": None,
        "featureElementCount": sum(
            1 for element in elements if element.attrib.get("data-gbdraw-feature-id")
        ),
        "texts": texts,
        "slots": slots,
        "annotations": annotations,
        "scriptCount": sum(1 for element in elements if _local_name(element.tag) == "script"),
        "eventAttributes": event_attributes,
        "unsafeHrefs": unsafe_hrefs,
    }


def _slot_pairs(report: Mapping[str, Any]) -> tuple[tuple[str, str], ...]:
    return tuple(
        (slot["slotId"], slot["renderer"])
        for slot in report["slots"]
        if not str(slot["slotId"]).startswith("__gbdraw_auto_")
    )


def _assert_safe_svg(report: Mapping[str, Any]) -> None:
    if report.get("svgCount") != 1:
        raise AssertionError("Expected one generated SVG")
    if report.get("scriptCount") != 0:
        raise AssertionError("Generated SVG contains a script")
    if report.get("eventAttributes"):
        raise AssertionError("Generated SVG contains event-handler attributes")
    if report.get("unsafeHrefs"):
        raise AssertionError("Generated SVG contains an unsafe external reference")


def _assert_quantitative_svg(report: Mapping[str, Any]) -> None:
    _assert_safe_svg(report)
    if "AP027133.1" not in report["recordIds"]:
        raise AssertionError("The quantitative diagram does not identify AP027133.1")
    if report["featureElementCount"] != 576:
        raise AssertionError(
            "Expected 576 displayed AP027133.1 features, found "
            f"{report['featureElementCount']}"
        )
    if report.get("logicalFeatureCount") not in (None, 576):
        raise AssertionError(
            "Expected 576 logical AP027133.1 features, found "
            f"{report['logicalFeatureCount']}"
        )
    if _slot_pairs(report) != EXPECTED_QUANTITATIVE_SLOT_ORDER:
        raise AssertionError(f"Unexpected quantitative slot order: {_slot_pairs(report)}")
    by_renderer = {slot["renderer"]: slot for slot in report["slots"]}
    if by_renderer["depth"]["pathCount"] < 1:
        raise AssertionError("The depth slot has no plotted path")
    if set(by_renderer["depth"]["tickTexts"]) != {
        "0x",
        "20x",
        "40x",
        "60x",
        "80x",
    }:
        raise AssertionError(
            f"Unexpected depth axis ticks: {by_renderer['depth']['tickTexts']}"
        )
    if set(by_renderer["dinucleotide_content"]["tickTexts"]) != {
        "0%",
        "20%",
        "40%",
        "60%",
        "80%",
        "100%",
    }:
        raise AssertionError(
            "Unexpected GC percentage ticks: "
            f"{by_renderer['dinucleotide_content']['tickTexts']}"
        )
    required_texts = {
        "AP027133.1",
        "606,194 bp",
        "Complete AP027133.1 genome with quantitative tracks",
        "DRR394922 mean depth",
        "GC content",
        "GC skew (+)",
        "GC skew (-)",
    }
    missing = required_texts.difference(report["texts"])
    if missing:
        raise AssertionError(f"Missing quantitative SVG text: {sorted(missing)}")


def _assert_annotation_svg(report: Mapping[str, Any]) -> None:
    _assert_safe_svg(report)
    if "NC_001879.2" not in report["recordIds"]:
        raise AssertionError("The annotation diagram does not identify NC_001879.2")
    if report["featureElementCount"] != 195:
        raise AssertionError(
            "Expected 195 rendered tobacco-plastome feature elements, found "
            f"{report['featureElementCount']}"
        )
    if report.get("logicalFeatureCount") not in (None, 145):
        raise AssertionError(
            "Expected 145 logical tobacco-plastome features, found "
            f"{report['logicalFeatureCount']}"
        )
    if _slot_pairs(report) != EXPECTED_ANNOTATION_SLOT_ORDER:
        raise AssertionError(f"Unexpected annotation slot order: {_slot_pairs(report)}")
    expected_annotations = {
        ("lsc", "plastome_regions", "annotations", "NC_001879.2", "bracket", "LSC"),
        ("irb", "plastome_regions", "annotations", "NC_001879.2", "bracket", "IRb"),
        ("ssc", "plastome_regions", "annotations", "NC_001879.2", "bracket", "SSC"),
        ("ira", "plastome_regions", "annotations", "NC_001879.2", "bracket", "IRa"),
    }
    actual_annotations = {
        (
            item["id"],
            item["setId"],
            item["trackId"],
            item["recordId"],
            item["mark"],
            item["label"],
        )
        for item in report["annotations"]
    }
    if actual_annotations != expected_annotations:
        raise AssertionError(f"Unexpected rendered annotations: {actual_annotations}")
    required_texts = {
        "NC_001879.2",
        "155,943 bp",
        "Complete Nicotiana tabacum plastome regions",
        "LSC",
        "IRb",
        "SSC",
        "IRa",
        "AT skew (+)",
        "AT skew (-)",
    }
    missing = required_texts.difference(report["texts"])
    if missing:
        raise AssertionError(f"Missing annotation SVG text: {sorted(missing)}")


def _resize_sidebar(page: Page) -> None:
    resize_handle = page.get_by_title("Drag to resize", exact=True)
    box = resize_handle.bounding_box()
    if box is None:
        raise AssertionError("Could not resolve the settings-pane resize handle")
    page.mouse.move(box["x"] + (box["width"] / 2), box["y"] + (box["height"] / 2))
    page.mouse.down()
    page.mouse.move(515, box["y"] + (box["height"] / 2), steps=8)
    page.mouse.up()


def _fit_circular_preview(
    page: Page,
    target_zoom: str = "70%",
    *,
    pan_left_ratio: float = 0.24,
    pan_up_ratio: float = 0.0,
) -> None:
    reset_zoom = page.get_by_role("button", name="Reset zoom", exact=True)
    zoom_out = page.get_by_role("button", name="Zoom out", exact=True)
    if "100%" not in reset_zoom.inner_text():
        reset_zoom.click()
        expect(reset_zoom).to_contain_text("100%")
    for _ in range(12):
        if target_zoom in reset_zoom.inner_text():
            break
        zoom_out.click()
    else:
        raise AssertionError(f"Could not reach Circular preview zoom {target_zoom}")

    if pan_left_ratio > 0.0:
        result_region = page.get_by_role("region", name="Result Preview", exact=True)
        box = result_region.bounding_box()
        if box is None:
            raise AssertionError("Could not resolve the Circular result preview bounds")
        y = box["y"] + (box["height"] * 0.52)
        start_ratio = 0.70
        page.mouse.move(box["x"] + (box["width"] * start_ratio), y)
        page.mouse.down()
        page.mouse.move(
            box["x"] + (box["width"] * (start_ratio - pan_left_ratio)),
            y - (box["height"] * pan_up_ratio),
            steps=10,
        )
        page.mouse.up()
    page.evaluate("() => window.getSelection()?.removeAllRanges()")


def _configure_title(page: Page, title: str) -> None:
    page.get_by_label("Title & Legend", exact=True).click()
    plot_title = page.get_by_label("Plot Title", exact=True)
    plot_title.fill(title)
    expect(plot_title).to_have_value(title)
    title_position = page.get_by_label("Plot Title Position", exact=True)
    title_position.select_option("top")
    expect(title_position).to_have_value("top")
    definition_size = page.get_by_label("Definition Font Size", exact=True)
    definition_size.fill("17")
    expect(definition_size).to_have_value("17")
    legend_position = page.get_by_label("Legend Position", exact=True)
    legend_position.select_option("right")
    expect(legend_position).to_have_value("right")


def _track_slot_snapshot(page: Page) -> tuple[dict[str, Any], ...]:
    snapshot = page.evaluate(
        """
        () => (window.__GBDRAW_APP__?.adv?.circular_track_slots || []).map((slot) => ({
          id: String(slot?.id || ''),
          renderer: String(slot?.renderer || ''),
          enabled: slot?.enabled !== false,
          side: slot?.side == null ? null : String(slot.side),
          params: JSON.parse(JSON.stringify(slot?.params || {}))
        }))
        """
    )
    return tuple(snapshot)


def _assert_slot_snapshot(
    snapshot: tuple[dict[str, Any], ...],
    expected: tuple[tuple[str, str], ...],
) -> None:
    actual = tuple(
        (slot["id"], slot["renderer"])
        for slot in snapshot
        if slot["enabled"]
    )
    if actual != expected:
        raise AssertionError(f"Unexpected browser track-slot state: {actual}")


def _download_svg(
    page: Page,
    download_dir: Path,
    *,
    expected_filename: str,
    assert_svg: Any,
) -> dict[str, Any]:
    svg_button = page.get_by_role("button", name="SVG", exact=True)
    with page.expect_download(timeout=ACTION_TIMEOUT_MS) as download_info:
        svg_button.click()
    download = download_info.value
    if download.failure() is not None:
        raise AssertionError(f"SVG download failed: {download.failure()}")
    if download.suggested_filename != expected_filename:
        raise AssertionError(
            f"Expected download {expected_filename}, found {download.suggested_filename}"
        )
    path = download_dir / download.suggested_filename
    download.save_as(path)
    report = _inspect_tracks_svg_file(path)
    assert_svg(report)
    return {
        "filename": download.suggested_filename,
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "path": str(path),
        "semantics": report,
    }


def capture_gui_quantitative_tracks(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> CaptureResult:
    """Run H-GUI-09 through visible controls in one fresh browser context."""

    genome_report = _validate_complete_record(
        GUI_QUANTITATIVE_GENBANK_PATH,
        expected_size=GUI_QUANTITATIVE_GENBANK_SIZE,
        expected_sha256=GUI_QUANTITATIVE_GENBANK_SHA256,
        expected_id="AP027133.1",
        expected_length=606_194,
    )
    depth_report = _validate_depth_fixture()
    assert_output_paths(
        output_paths,
        GUI_QUANTITATIVE_SCREENSHOT_NAMES,
        GUI_QUANTITATIVE_SCENARIO_ID,
    )
    download_dir.mkdir(parents=True, exist_ok=True)

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_app_shell(page)
        _resize_sidebar(page)

        circular = page.get_by_role("button", name="Circular", exact=True)
        circular.click()
        expect(circular).to_have_attribute("aria-pressed", "true")
        genbank = page.get_by_role("radio", name="GenBank", exact=True)
        genbank.check()
        page.get_by_label("GenBank/DDBJ File", exact=True).set_input_files(
            GUI_QUANTITATIVE_GENBANK_PATH
        )
        expect(
            page.get_by_role(
                "group", name="GenBank/DDBJ File selection", exact=True
            )
        ).to_contain_text("AP027133.gb")

        prefix = page.get_by_label("Output Prefix", exact=True)
        prefix.fill("quantitative_tracks")
        expect(prefix).to_have_value("quantitative_tracks")

        depth_section = page.get_by_label("Depth TSV tracks", exact=True)
        depth_section.click()
        page.get_by_label("Depth TSV", exact=True).set_input_files(
            GUI_QUANTITATIVE_DEPTH_PATH
        )
        expect(
            page.get_by_role("group", name="Depth TSV selection", exact=True)
        ).to_contain_text("AP027133.DRR394922.depth-1kb.tsv")

        depth_window = page.get_by_label("Depth Window", exact=True)
        depth_window.fill("1")
        expect(depth_window).to_have_value("1")
        depth_step = page.get_by_label("Depth Step", exact=True)
        depth_step.fill("1000")
        expect(depth_step).to_have_value("1000")
        depth_min = page.get_by_label("Depth Min", exact=True)
        depth_min.fill("0")
        expect(depth_min).to_have_value("0")
        depth_max = page.get_by_label("Depth Max", exact=True)
        depth_max.fill("80")
        expect(depth_max).to_have_value("80")

        log_scale = page.get_by_role("checkbox", name="Log Scale", exact=True)
        log_scale.uncheck()
        depth_axis = page.get_by_role("checkbox", name="Depth Axis", exact=True)
        depth_axis.check()
        depth_ticks = page.get_by_role("checkbox", name="Depth Ticks", exact=True)
        depth_ticks.check()

        depth_title = page.get_by_label("Depth legend title", exact=True)
        depth_title.fill("DRR394922 mean depth")
        expect(depth_title).to_have_value("DRR394922 mean depth")
        large_tick = page.get_by_label("Depth large tick", exact=True)
        large_tick.fill("20")
        expect(large_tick).to_have_value("20")
        small_tick = page.get_by_label("Depth small tick", exact=True)
        small_tick.fill("10")
        expect(small_tick).to_have_value("10")

        dinucleotide_section = page.get_by_label(
            "Dinucleotide content/skew", exact=True
        )
        dinucleotide_section.click()
        gc_mode = page.get_by_role("combobox").filter(
            has=page.get_by_role("option", name="Percent", exact=True)
        )
        gc_mode.select_option("percent")
        expect(gc_mode).to_have_value("percent")
        page.get_by_role("checkbox", name="Percent Axis", exact=True).check()
        page.get_by_role("checkbox", name="Percent Ticks", exact=True).check()

        custom_slots = page.get_by_role(
            "button", name=re.compile(r"Custom Track Slots$"), exact=False
        )
        custom_slots.click()
        use_custom = page.get_by_role(
            "checkbox", name="Use custom stack", exact=True
        )
        use_custom.check()
        expect(use_custom).to_be_checked()
        slots = _track_slot_snapshot(page)
        _assert_slot_snapshot(slots, EXPECTED_QUANTITATIVE_SLOT_ORDER)

        _configure_title(page, "Complete AP027133.1 genome with quantitative tracks")

        depth_section.scroll_into_view_if_needed()
        screenshot_bytes["track-settings.png"] = capture_screenshot(
            page,
            output_paths["track-settings.png"],
            "Circular",
        )

        final_report = generate_and_inspect(
            page,
            _inspect_tracks_svg,
            _assert_quantitative_svg,
        )
        set_feature_search_visible(page, visible=False)
        _fit_circular_preview(page)
        depth_section.scroll_into_view_if_needed()
        screenshot_bytes["track-result.png"] = capture_screenshot(
            page,
            output_paths["track-result.png"],
            "Circular",
        )
        download_report = _download_svg(
            page,
            download_dir,
            expected_filename=EXPECTED_GUI_QUANTITATIVE_SVG,
            assert_svg=_assert_quantitative_svg,
        )
        capture.assert_clean()
    finally:
        capture.close()

    return CaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        fixture_report={**genome_report, "depth": depth_report},
        track_slots=slots,
    )


def capture_gui_annotation_tracks(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> CaptureResult:
    """Run H-GUI-10 through visible controls in one fresh browser context."""

    genome_report = _validate_complete_record(
        GUI_ANNOTATION_GENBANK_PATH,
        expected_size=GUI_ANNOTATION_GENBANK_SIZE,
        expected_sha256=GUI_ANNOTATION_GENBANK_SHA256,
        expected_id="NC_001879.2",
        expected_length=155_943,
    )
    annotation_report = _validate_annotation_fixture()
    assert_output_paths(
        output_paths,
        GUI_ANNOTATION_SCREENSHOT_NAMES,
        GUI_ANNOTATION_SCENARIO_ID,
    )
    download_dir.mkdir(parents=True, exist_ok=True)

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_app_shell(page)
        _resize_sidebar(page)

        circular = page.get_by_role("button", name="Circular", exact=True)
        circular.click()
        expect(circular).to_have_attribute("aria-pressed", "true")
        genbank = page.get_by_role("radio", name="GenBank", exact=True)
        genbank.check()
        page.get_by_label("GenBank/DDBJ File", exact=True).set_input_files(
            GUI_ANNOTATION_GENBANK_PATH
        )
        expect(
            page.get_by_role(
                "group", name="GenBank/DDBJ File selection", exact=True
            )
        ).to_contain_text("NC_001879.gbk")

        prefix = page.get_by_label("Output Prefix", exact=True)
        prefix.fill("region_annotations_and_slots")
        expect(prefix).to_have_value("region_annotations_and_slots")
        track_preset = page.get_by_label("Track Preset", exact=True)
        track_preset.select_option("middle")
        expect(track_preset).to_have_value("middle")

        annotations_section = page.get_by_label("Region Annotations", exact=True)
        annotations_section.click()
        page.get_by_label("Import TSV", exact=True).set_input_files(
            GUI_ANNOTATION_TABLE_PATH
        )
        annotation_set_id = page.get_by_label("Annotation set id", exact=True)
        expect(annotation_set_id).to_have_value("plastome_regions")
        legend_label = page.get_by_placeholder("Set legend label (optional)")
        legend_label.fill("Plastome structural regions")
        expect(legend_label).to_have_value("Plastome structural regions")
        expect(page.get_by_label("Annotation target record", exact=True)).to_have_count(4)
        lane_controls = page.get_by_label("Annotation lane", exact=True)
        expect(lane_controls).to_have_count(4)
        for index, lane in enumerate((0, 1, 0, 1)):
            lane_control = lane_controls.nth(index)
            lane_control.fill(str(lane))
            expect(lane_control).to_have_value(str(lane))
        rendered_lanes = tuple(
            page.evaluate(
                """
                () => (window.__GBDRAW_APP__?.annotationSets?.[0]?.annotations || [])
                  .map((item) => Number(item?.lane))
                """
            )
        )
        if rendered_lanes != (0, 1, 0, 1):
            raise AssertionError(
                f"Unexpected visible annotation-lane assignment: {rendered_lanes}"
            )

        custom_slots = page.get_by_role(
            "button", name=re.compile(r"Custom Track Slots$"), exact=False
        )
        custom_slots.click()
        use_custom = page.get_by_role(
            "checkbox", name="Use custom stack", exact=True
        )
        use_custom.check()
        expect(use_custom).to_be_checked()

        skew_slot = page.get_by_role(
            "group", name="Circular track slot gc_skew", exact=True
        )
        track_nt = skew_slot.get_by_label("Track dinucleotide", exact=True)
        track_nt.fill("AT")
        expect(track_nt).to_have_value("AT")
        track_legend = skew_slot.get_by_label("Track legend label", exact=True)
        track_legend.fill("AT skew")
        expect(track_legend).to_have_value("AT skew")

        new_renderer = page.get_by_label("New circular track renderer", exact=True)
        new_renderer.select_option("annotations")
        expect(new_renderer).to_have_value("annotations")
        page.get_by_role(
            "button", name=re.compile(r"Add track$"), exact=False
        ).click()
        annotation_slot = page.get_by_role(
            "group", name="Circular track slot annotations", exact=True
        )
        expect(annotation_slot).to_be_visible()
        annotation_set = annotation_slot.get_by_label("Annotation set", exact=True)
        annotation_set.select_option("plastome_regions")
        expect(annotation_set).to_have_value("plastome_regions")
        move_outside = annotation_slot.get_by_title("Move outside Axis", exact=True)
        expect(move_outside).to_be_enabled()
        move_outside.click()
        annotation_slot = page.get_by_role(
            "group", name="Circular track slot annotations", exact=True
        )
        placement = annotation_slot.get_by_label("Annotation placement", exact=True)
        expect(placement).to_have_value("outside")
        show_labels = annotation_slot.get_by_label(
            "Show annotation labels", exact=True
        )
        show_labels.check()
        expect(show_labels).to_be_checked()

        slots = _track_slot_snapshot(page)
        _assert_slot_snapshot(slots, EXPECTED_ANNOTATION_SLOT_ORDER)
        skew_snapshot = next(slot for slot in slots if slot["id"] == "gc_skew")
        annotation_snapshot = next(
            slot for slot in slots if slot["id"] == "annotations"
        )
        if skew_snapshot["params"].get("nt") != "AT":
            raise AssertionError("The visible custom skew control did not select AT")
        if annotation_snapshot["side"] != "outside":
            raise AssertionError("The annotation slot is not outside the circular axis")
        if annotation_snapshot["params"].get("set_id") != "plastome_regions":
            raise AssertionError("The annotation slot is not bound to plastome_regions")

        _configure_title(page, "Complete Nicotiana tabacum plastome regions")

        annotation_slot.scroll_into_view_if_needed()
        screenshot_bytes["slot-settings.png"] = capture_screenshot(
            page,
            output_paths["slot-settings.png"],
            "Circular",
        )

        final_report = generate_and_inspect(
            page,
            _inspect_tracks_svg,
            _assert_annotation_svg,
        )
        set_feature_search_visible(page, visible=False)
        _fit_circular_preview(page, pan_up_ratio=0.01)
        annotation_slot.scroll_into_view_if_needed()
        screenshot_bytes["annotation-result.png"] = capture_screenshot(
            page,
            output_paths["annotation-result.png"],
            "Circular",
        )
        download_report = _download_svg(
            page,
            download_dir,
            expected_filename=EXPECTED_GUI_ANNOTATION_SVG,
            assert_svg=_assert_annotation_svg,
        )
        capture.assert_clean()
    finally:
        capture.close()

    return CaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        fixture_report={
            **genome_report,
            "annotations": {**annotation_report, "renderedLanes": rendered_lanes},
        },
        track_slots=slots,
    )


__all__ = [
    "GUI_ANNOTATION_SCREENSHOT_NAMES",
    "GUI_QUANTITATIVE_SCREENSHOT_NAMES",
    "capture_gui_annotation_tracks",
    "capture_gui_quantitative_tracks",
]
