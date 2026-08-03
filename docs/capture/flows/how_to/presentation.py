"""Capture the H-GUI-11 and H-GUI-12 presentation-control journeys."""

from __future__ import annotations

import math
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from Bio import SeqIO
from playwright.sync_api import BrowserType, Page, expect

from assertions.svg_semantics import assert_static_svg_safety, inspect_svg_file
from config import ACTION_TIMEOUT_MS
from flows.web_capture import (
    assert_fixture_identity,
    assert_output_paths,
    capture_screenshot,
    generate_and_inspect,
    open_browser_capture,
    wait_for_worker,
)


REPO_ROOT = Path(__file__).resolve().parents[4]
HUMAN_MITOCHONDRION_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "human-mitochondrion"
    / "HmmtDNA.gbk"
)
HUMAN_MITOCHONDRION_SIZE = 64_640
HUMAN_MITOCHONDRION_SHA256 = (
    "7530d659e7174272372814edfecb2ece1f87a444395a861fcdf1b977c4aa5c1f"
)
DEFAULT_COLOR_RULE_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "tobacco-plastome-regions"
    / "modified_default_colors.tsv"
)
DEFAULT_COLOR_RULE_SIZE = 95
DEFAULT_COLOR_RULE_SHA256 = (
    "e48654dfc5225c8c1eec251f773fc07892228dee906cb1e105e4d24cb5ae8bc1"
)
GENE_PRIORITY_RULE_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "shared"
    / "cds_gene_qualifier_priority.tsv"
)
GENE_PRIORITY_RULE_SIZE = 9
GENE_PRIORITY_RULE_SHA256 = (
    "1d3c787c5768b52191cc7e8a325cd00fdc4b1e9b495d37e580c3a69dec21b87a"
)
VISIBILITY_RULE_PATH = (
    REPO_ROOT
    / "gbdraw"
    / "web"
    / "tutorial-data"
    / "human-mitochondrion"
    / "HmmtDNA_feature_visibility.tsv"
)
VISIBILITY_RULE_SIZE = 212
VISIBILITY_RULE_SHA256 = (
    "aca42671c08de16dfbb2317d8b8f53801d1601e956c04aab6de5deaed2a5f348"
)

STYLE_SCREENSHOT_NAMES = ("style-settings.png", "style-result.png")
PRESENTATION_SCREENSHOT_NAMES = (
    "presentation-settings.png",
    "presentation-result.png",
)
STYLE_OUTPUT_PREFIX = "styled_features_labels_legend"
PRESENTATION_OUTPUT_PREFIX = "feature_visibility_shapes"
EXPECTED_STYLE_SVG = f"{STYLE_OUTPUT_PREFIX}.svg"
EXPECTED_PRESENTATION_SVG = f"{PRESENTATION_OUTPUT_PREFIX}.svg"
STYLE_TITLE = "Human mitochondrial genome: selected gene labels"
PRESENTATION_TITLE = "Human mitochondrial feature presentation"
STYLE_LABEL_PATTERN = r"^(COX[123]|ATP[68]|CYTB|ND4L|ND4)$"
STYLE_LABELS = frozenset(
    {"COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB", "ND4L", "ND4"}
)
EXPECTED_LEGEND_ORDER = (
    "tRNA",
    "rRNA",
    "CDS",
    "GC content",
    "GC skew (+)",
    "GC skew (-)",
)
EXPECTED_VISIBILITY_RULES = (
    {
        "recordId": "NC_012920.1",
        "featureType": "D-loop",
        "qualifier": "location",
        "value": r"^0\.\.16569$",
        "action": "show",
    },
    {
        "recordId": "NC_012920.1",
        "featureType": "CDS",
        "qualifier": "product",
        "value": r"^cytochrome c oxidase subunit I$",
        "action": "off",
    },
    {
        "recordId": "*",
        "featureType": "CDS",
        "qualifier": "product",
        "value": r"^ATP synthase F0 subunit 6$",
        "action": "exclude_matching",
    },
)


_PRESENTATION_INSPECTION_SCRIPT = r"""
root => {
  const svg = root.getElementsByTagName('svg')[0];
  if (!svg) return { svgCount: 0 };
  const elements = [svg, ...Array.from(svg.getElementsByTagName('*'))];
  const features = elements.filter((element) => (
    Boolean(element.getAttribute('data-gbdraw-feature-id')) &&
    String(element.getAttribute('data-gbdraw-feature-part') || 'block') === 'block'
  ));
  const metadata = Array.isArray(window.__GBDRAW_APP__?.extractedFeatures)
    ? window.__GBDRAW_APP__.extractedFeatures
    : [];
  const firstText = (value) => Array.isArray(value)
    ? String(value[0] || '')
    : String(value || '');
  const metadataRows = metadata.map((feature) => ({
    type: String(feature?.type || ''),
    recordId: String(feature?.record_id || ''),
    featureId: String(
      feature?.rendered_svg_id ||
      feature?.rendered_feature_svg_id ||
      feature?.svg_id ||
      ''
    ),
    gene: firstText(feature?.gene || feature?.qualifiers?.gene),
    product: firstText(feature?.product || feature?.qualifiers?.product),
    start: Number(feature?.start ?? feature?.location?.start ?? NaN),
    end: Number(feature?.end ?? feature?.location?.end ?? NaN),
    strand: String(feature?.strand || '')
  }));
  const metadataById = new Map(metadataRows.map((row) => [row.featureId, row]));
  const featureRows = features.map((element) => {
    const id = String(element.getAttribute('data-gbdraw-feature-id') || '');
    const path = String(element.getAttribute('d') || '');
    let radialMidpoint = null;
    try {
      const length = element.getTotalLength();
      const point = element.getPointAtLength(length * 0.5);
      radialMidpoint = Math.hypot(point.x, point.y);
    } catch (_error) {
      radialMidpoint = null;
    }
    return {
      id,
      type: String(metadataById.get(id)?.type || ''),
      gene: String(metadataById.get(id)?.gene || ''),
      product: String(metadataById.get(id)?.product || ''),
      tag: String(element.localName || element.tagName || '').toLowerCase(),
      fill: String(element.getAttribute('fill') || '').toLowerCase(),
      stroke: String(element.getAttribute('stroke') || '').toLowerCase(),
      strokeWidth: String(element.getAttribute('stroke-width') || ''),
      underlay: element.getAttribute('data-gbdraw-auto-feature-underlay') === 'true',
      path,
      lineCommandCount: (path.match(/[Ll]/g) || []).length,
      radialMidpoint
    };
  });
  const legendRows = elements
    .filter((element) => Boolean(element.getAttribute('data-legend-key')))
    .map((element) => {
      const swatch = element.querySelector('path, rect, circle');
      return {
        key: String(element.getAttribute('data-legend-key') || ''),
        fill: String(swatch?.getAttribute('fill') || '').toLowerCase(),
        stroke: String(swatch?.getAttribute('stroke') || '').toLowerCase(),
        strokeWidth: String(swatch?.getAttribute('stroke-width') || '')
      };
    });
  const texts = Array.from(svg.getElementsByTagName('text'))
    .map((node) => String(node.textContent || '').replace(/\s+/g, ' ').trim())
    .filter(Boolean);
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
  const legend = elements.find((element) => element.getAttribute('id') === 'legend');
  const recordGroup = elements.find((element) => (
    String(element.localName || element.tagName || '').toLowerCase() === 'g' &&
    element.getAttribute('data-gbdraw-record-id') === 'NC_012920.1'
  ));
  const parseTranslate = (element) => {
    const match = String(element?.getAttribute('transform') || '')
      .match(/translate\(\s*([-+\d.eE]+)[,\s]+([-+\d.eE]+)\s*\)/);
    return match ? [Number(match[1]), Number(match[2])] : null;
  };
  const app = window.__GBDRAW_APP__;
  const fileName = (file) => String(file?.name || '');
  return {
    svgCount: 1,
    recordIds: elements
      .map((element) => String(element.getAttribute('data-gbdraw-record-id') || ''))
      .filter(Boolean),
    featureElementCount: featureRows.length,
    features: featureRows,
    metadata: metadataRows,
    texts,
    legend: legendRows,
    recordTranslate: parseTranslate(recordGroup),
    legendTranslate: parseTranslate(legend),
    state: {
      prefix: String(app?.form?.prefix || ''),
      title: String(app?.form?.plot_title || ''),
      titlePosition: String(app?.adv?.plot_title_position || ''),
      legendPosition: String(app?.form?.legend || ''),
      labelsMode: String(app?.form?.labels_mode || ''),
      separateStrands: Boolean(app?.form?.separate_strands),
      trackPreset: String(app?.form?.track_type || ''),
      palette: String(app?.selectedPalette || ''),
      filterMode: String(app?.filterMode || ''),
      whitelist: Array.isArray(app?.manualWhitelist)
        ? app.manualWhitelist.map((rule) => ({
            feat: String(rule?.feat || ''),
            qual: String(rule?.qual || ''),
            key: String(rule?.key || '')
          }))
        : [],
      defaultColorFile: fileName(app?.files?.d_color),
      priorityFile: fileName(app?.files?.qualifier_priority),
      featureShapes: { ...(app?.adv?.feature_shapes || {}) },
      arrowHeadLengthRatio: String(app?.adv?.arrow_head_length_ratio ?? ''),
      arrowShaftWidthRatio: Number(app?.adv?.arrow_shaft_width_ratio),
      blockStrokeColor: String(app?.adv?.block_stroke_color || '').toLowerCase(),
      blockStrokeWidth: Number(app?.adv?.block_stroke_width),
      lineStrokeColor: String(app?.adv?.line_stroke_color || '').toLowerCase(),
      lineStrokeWidth: Number(app?.adv?.line_stroke_width),
      resolveOverlaps: Boolean(app?.adv?.resolve_overlaps),
      visibilityRules: Array.isArray(app?.featureVisibilityManualRules)
        ? app.featureVisibilityManualRules.map((rule) => ({
            recordId: String(rule?.recordId || ''),
            featureType: String(rule?.featureType || ''),
            qualifier: String(rule?.qualifier || ''),
            value: String(rule?.value || ''),
            action: String(rule?.action || '')
          }))
        : []
    },
    interactiveSvg: svg.getAttribute('data-gbdraw-interactive-svg') === 'true',
    scriptCount: svg.getElementsByTagName('script').length,
    eventAttributes,
    unsafeHrefs
  };
}
"""


@dataclass(frozen=True)
class CaptureResult:
    screenshot_bytes: dict[str, int]
    final_svg_semantics: dict[str, Any]
    download: dict[str, Any]
    source_record: dict[str, Any]


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1].lower()


def _inspect_downloaded_presentation_svg(path: Path) -> dict[str, Any]:
    report = inspect_svg_file(path)
    root = ET.parse(path).getroot()
    elements = list(root.iter())
    feature_elements = [
        element
        for element in elements
        if element.attrib.get("data-gbdraw-feature-id")
        and element.attrib.get("data-gbdraw-feature-part", "block") == "block"
    ]
    legend_rows = []
    for element in elements:
        key = element.attrib.get("data-legend-key")
        if not key:
            continue
        swatch = next(
            (
                child
                for child in element.iter()
                if _local_name(child.tag) in {"path", "rect", "circle"}
            ),
            None,
        )
        legend_rows.append(
            {
                "key": key,
                "fill": (swatch.attrib.get("fill", "") if swatch is not None else "").lower(),
                "stroke": (
                    swatch.attrib.get("stroke", "") if swatch is not None else ""
                ).lower(),
                "strokeWidth": (
                    swatch.attrib.get("stroke-width", "")
                    if swatch is not None
                    else ""
                ),
            }
        )
    report.update(
        {
            "featureElementCount": len(feature_elements),
            "features": [
                {
                    "id": element.attrib.get("data-gbdraw-feature-id", ""),
                    "fill": element.attrib.get("fill", "").lower(),
                    "stroke": element.attrib.get("stroke", "").lower(),
                    "strokeWidth": element.attrib.get("stroke-width", ""),
                    "underlay": element.attrib.get(
                        "data-gbdraw-auto-feature-underlay"
                    )
                    == "true",
                    "path": element.attrib.get("d", ""),
                    "lineCommandCount": sum(
                        element.attrib.get("d", "").count(command)
                        for command in ("L", "l")
                    ),
                }
                for element in feature_elements
            ],
            "legend": legend_rows,
        }
    )
    return report


def _inspect_presentation_svg(result_region: Any) -> dict[str, Any]:
    report = result_region.evaluate(_PRESENTATION_INSPECTION_SCRIPT)
    if not isinstance(report, dict):
        raise AssertionError("Presentation SVG inspection did not return a report")
    return report


def _assert_whole_human_mitochondrion(report: dict[str, Any]) -> None:
    if report.get("svgCount") != 1:
        raise AssertionError(f"Expected one SVG, found {report.get('svgCount')}")
    if "NC_012920.1" not in set(report.get("recordIds", [])):
        raise AssertionError("The result does not contain complete NC_012920.1")
    texts = set(report.get("texts", []))
    if not {"NC_012920.1", "16,569 bp"}.issubset(texts):
        raise AssertionError("The complete human mitochondrial identity is missing")
    assert_static_svg_safety(report)


def _assert_style_svg(report: dict[str, Any]) -> None:
    _assert_whole_human_mitochondrion(report)
    if report.get("featureElementCount") != 37:
        raise AssertionError(
            "Styling must preserve all 37 CDS, rRNA, and tRNA features; found "
            f"{report.get('featureElementCount')}"
        )
    texts = set(report.get("texts", []))
    if STYLE_TITLE not in texts:
        raise AssertionError("The styled SVG is missing its plot title")
    missing_labels = STYLE_LABELS - texts
    if missing_labels:
        raise AssertionError(f"Missing gene-qualifier labels: {sorted(missing_labels)}")
    product_labels = {
        "cytochrome c oxidase subunit I",
        "cytochrome c oxidase subunit II",
        "cytochrome c oxidase subunit III",
        "ATP synthase F0 subunit 6",
    }
    if texts & product_labels:
        raise AssertionError("CDS labels came from product instead of gene")

    fills = {feature.get("fill") for feature in report.get("features", [])}
    for expected_fill in ("#d3d3d3", "#009e73", "#e69f00"):
        if expected_fill not in fills:
            raise AssertionError(f"Missing default-color override fill {expected_fill}")
    legend = report.get("legend", [])
    legend_order = tuple(row.get("key") for row in legend)
    if legend_order != EXPECTED_LEGEND_ORDER:
        raise AssertionError(f"Unexpected legend order: {legend_order!r}")
    record_translate = report.get("recordTranslate")
    legend_translate = report.get("legendTranslate")
    if (
        record_translate is not None
        and legend_translate is not None
        and legend_translate[0] <= record_translate[0]
    ):
        raise AssertionError("The legend is not positioned to the right of the record")
    legend_fills = {row.get("key"): row.get("fill") for row in legend}
    if legend_fills.get("GC content") != "#fdd7aa":
        raise AssertionError("The soft_pastels GC-content color was not applied")
    if legend_fills.get("GC skew (+)") != "#f9e2ae":
        raise AssertionError("The soft_pastels positive-skew color was not applied")
    if legend_fills.get("GC skew (-)") != "#b5d8eb":
        raise AssertionError("The soft_pastels negative-skew color was not applied")

    state = report.get("state")
    if state is not None:
        expected_state = {
            "prefix": STYLE_OUTPUT_PREFIX,
            "title": STYLE_TITLE,
            "titlePosition": "top",
            "legendPosition": "right",
            "labelsMode": "out",
            "separateStrands": True,
            "trackPreset": "middle",
            "palette": "soft_pastels",
            "filterMode": "Whitelist",
            "whitelist": [
                {"feat": "CDS", "qual": "gene", "key": STYLE_LABEL_PATTERN}
            ],
            "defaultColorFile": DEFAULT_COLOR_RULE_PATH.name,
            "priorityFile": GENE_PRIORITY_RULE_PATH.name,
        }
        for key, value in expected_state.items():
            if state.get(key) != value:
                raise AssertionError(
                    f"Unexpected H-GUI-11 control state for {key}: {state.get(key)!r}"
                )


def _normalized_visibility_rules(report: dict[str, Any]) -> tuple[dict[str, str], ...]:
    return tuple(
        {
            key: str(rule.get(key, ""))
            for key in ("recordId", "featureType", "qualifier", "value", "action")
        }
        for rule in report.get("state", {}).get("visibilityRules", [])
    )


def _assert_feature_presentation_svg(report: dict[str, Any]) -> None:
    _assert_whole_human_mitochondrion(report)
    if report.get("featureElementCount") != 37:
        raise AssertionError(
            "The show/off rules must replace COX1 with D-loop and keep 37 visible "
            f"features; found {report.get('featureElementCount')}; "
            f"rules={_normalized_visibility_rules(report)!r}; "
            f"types={sorted({row.get('type') for row in report.get('metadata', [])})!r}"
        )
    if PRESENTATION_TITLE not in set(report.get("texts", [])):
        raise AssertionError("The feature-presentation SVG is missing its title")

    state = report.get("state")
    if state is not None:
        shapes = state.get("featureShapes", {})
        expected_shapes = {"CDS": "arrow", "rRNA": "rectangle", "tRNA": "underlay"}
        for feature_type, expected_shape in expected_shapes.items():
            if shapes.get(feature_type) != expected_shape:
                raise AssertionError(
                    f"Expected {feature_type}={expected_shape}, found "
                    f"{shapes.get(feature_type)!r}"
                )
        expected_state = {
            "prefix": PRESENTATION_OUTPUT_PREFIX,
            "title": PRESENTATION_TITLE,
            "titlePosition": "top",
            "legendPosition": "right",
            "separateStrands": False,
            "trackPreset": "middle",
            "arrowHeadLengthRatio": "1",
            "arrowShaftWidthRatio": 0.55,
            "blockStrokeColor": "#263238",
            "blockStrokeWidth": 2.5,
            "lineStrokeColor": "#455a64",
            "lineStrokeWidth": 2,
            "resolveOverlaps": True,
        }
        for key, value in expected_state.items():
            if state.get(key) != value:
                raise AssertionError(
                    f"Unexpected H-GUI-12 control state for {key}: {state.get(key)!r}"
                )
        if _normalized_visibility_rules(report) != EXPECTED_VISIBILITY_RULES:
            raise AssertionError(
                "Visible feature rules differ from the repository-managed table: "
                f"{_normalized_visibility_rules(report)!r}"
            )

    features = report.get("features", [])
    underlays = [feature for feature in features if feature.get("underlay")]
    if len(underlays) != 22:
        raise AssertionError(f"Expected 22 tRNA underlays, found {len(underlays)}")
    typed_underlays = [feature for feature in underlays if feature.get("type")]
    if typed_underlays and {feature.get("type") for feature in typed_underlays} != {
        "tRNA"
    }:
        raise AssertionError("A non-tRNA feature was rendered as the automatic underlay")

    metadata = report.get("metadata", [])
    if metadata:
        if any(
            row.get("product") == "cytochrome c oxidase subunit I"
            for row in metadata
        ):
            raise AssertionError("The off rule did not remove COX1")
        if not any(row.get("type") == "D-loop" for row in metadata):
            raise AssertionError("The show rule did not reveal the D-loop")
        if not any(
            row.get("product") == "ATP synthase F0 subunit 6" for row in metadata
        ):
            raise AssertionError("Exclude from matching incorrectly hid ATP6")

    foreground = [feature for feature in features if not feature.get("underlay")]
    if not foreground:
        raise AssertionError("The SVG has no foreground feature glyphs")
    for feature in foreground:
        if feature.get("stroke") != "#263238" or not math.isclose(
            float(feature.get("strokeWidth") or 0), 2.5
        ):
            raise AssertionError(
                "Foreground feature strokes do not match #263238 at 2.5 px"
            )

    typed_foreground = [feature for feature in foreground if feature.get("type")]
    if typed_foreground:
        typed_counts = {
            feature_type: sum(
                feature.get("type") == feature_type for feature in typed_foreground
            )
            for feature_type in ("D-loop", "rRNA", "CDS")
        }
        if typed_counts != {"D-loop": 1, "rRNA": 2, "CDS": 12}:
            raise AssertionError(f"Unexpected foreground shape counts: {typed_counts!r}")
        rectangles = [
            feature
            for feature in typed_foreground
            if feature.get("type") in {"D-loop", "rRNA"}
        ]
        if any(feature.get("lineCommandCount") != 1 for feature in rectangles):
            raise AssertionError("D-loop or rRNA did not render as a rectangle")
        arrows = [
            feature for feature in typed_foreground if feature.get("type") == "CDS"
        ]
        if any(feature.get("lineCommandCount", 0) < 2 for feature in arrows):
            raise AssertionError("A CDS did not render with arrow-head geometry")

    radial_values = [
        float(feature["radialMidpoint"])
        for feature in foreground
        if feature.get("radialMidpoint") is not None
    ]
    if radial_values:
        radial_bands = {round(value / 4) * 4 for value in radial_values}
        if len(radial_bands) < 2:
            raise AssertionError(
                "Resolve Overlaps did not produce more than one visible feature lane"
            )


def _verify_fixture_contracts(*, include_style_tables: bool) -> dict[str, Any]:
    assert_fixture_identity(
        HUMAN_MITOCHONDRION_PATH,
        expected_size=HUMAN_MITOCHONDRION_SIZE,
        expected_sha256=HUMAN_MITOCHONDRION_SHA256,
    )
    if include_style_tables:
        assert_fixture_identity(
            DEFAULT_COLOR_RULE_PATH,
            expected_size=DEFAULT_COLOR_RULE_SIZE,
            expected_sha256=DEFAULT_COLOR_RULE_SHA256,
        )
        assert_fixture_identity(
            GENE_PRIORITY_RULE_PATH,
            expected_size=GENE_PRIORITY_RULE_SIZE,
            expected_sha256=GENE_PRIORITY_RULE_SHA256,
        )
    else:
        assert_fixture_identity(
            VISIBILITY_RULE_PATH,
            expected_size=VISIBILITY_RULE_SIZE,
            expected_sha256=VISIBILITY_RULE_SHA256,
        )
    records = list(SeqIO.parse(HUMAN_MITOCHONDRION_PATH, "genbank"))
    if len(records) != 1:
        raise AssertionError("HmmtDNA.gbk must contain exactly one source record")
    source = records[0]
    if (
        source.id != "NC_012920.1"
        or len(source) != 16_569
        or str(source.annotations.get("topology", "")).lower() != "circular"
        or "complete" not in source.description.lower()
    ):
        raise AssertionError("H-GUI-11/12 require the complete circular NC_012920.1")
    return {
        "record_id": source.id,
        "length": len(source),
        "topology": source.annotations.get("topology"),
        "description": source.description,
    }


def _open_human_circular(page: Page, prefix: str) -> None:
    circular = page.get_by_role("button", name="Circular", exact=True)
    circular.click()
    expect(circular).to_have_attribute("aria-pressed", "true")
    genbank = page.get_by_role("radio", name="GenBank", exact=True)
    genbank.check()
    expect(genbank).to_be_checked()
    page.get_by_label("GenBank/DDBJ File", exact=True).set_input_files(
        HUMAN_MITOCHONDRION_PATH
    )
    expect(
        page.get_by_role(
            "group", name="GenBank/DDBJ File selection", exact=True
        )
    ).to_contain_text(HUMAN_MITOCHONDRION_PATH.name)
    output_prefix = page.get_by_label("Output Prefix", exact=True)
    output_prefix.fill(prefix)
    expect(output_prefix).to_have_value(prefix)
    track_preset = page.get_by_label("Track Preset", exact=True)
    track_preset.select_option("middle")
    expect(track_preset).to_have_value("middle")


def _set_plot_title(page: Page, title: str) -> None:
    title_panel = page.get_by_label("Title & Legend", exact=True)
    title_panel.click()
    plot_title = page.get_by_label("Plot Title", exact=True)
    plot_title.fill(title)
    expect(plot_title).to_have_value(title)
    title_position = page.get_by_label("Plot Title Position", exact=True)
    title_position.select_option("top")
    expect(title_position).to_have_value("top")
    legend_position = page.get_by_label("Legend Position", exact=True)
    legend_position.select_option("right")
    expect(legend_position).to_have_value("right")
    keep_definition = page.get_by_label(
        "Keep Full Definition with Plot Title", exact=True
    )
    keep_definition.check()
    expect(keep_definition).to_be_checked()
    title_panel.click()


def _fit_finished_circular_preview(page: Page) -> None:
    reset_zoom = page.get_by_role("button", name="Reset zoom", exact=True)
    zoom_out = page.get_by_role("button", name="Zoom out", exact=True)
    for _ in range(10):
        if "60%" in reset_zoom.inner_text():
            break
        zoom_out.click()
    else:
        raise AssertionError("Could not reach the documented 60% preview zoom")


def _park_feature_search(page: Page) -> None:
    """Drag the public search palette above the clipped preview viewport."""

    search = page.get_by_role("searchbox", name="Search features", exact=True)
    box = search.bounding_box()
    if box is None:
        raise AssertionError("Could not resolve the feature-search palette")
    x = box["x"] - 4
    y = box["y"] + box["height"] + 3
    page.mouse.move(x, y)
    page.mouse.down()
    page.mouse.move(x, 2, steps=10)
    page.mouse.up()


def _pan_preview_left(page: Page, distance_ratio: float = 0.27) -> None:
    """Use the public canvas gesture to reveal the right-side legend."""

    result_region = page.get_by_role("region", name="Result Preview", exact=True)
    box = result_region.bounding_box()
    if box is None:
        raise AssertionError("Could not resolve the result preview for panning")
    y = box["y"] + (box["height"] * 0.75)
    page.mouse.move(box["x"] + (box["width"] * 0.88), y)
    page.mouse.down()
    page.mouse.move(
        box["x"] + (box["width"] * (0.88 - distance_ratio)),
        y,
        steps=10,
    )
    page.mouse.up()
    page.keyboard.press("Escape")


def _wait_for_preview_transform(page: Page) -> None:
    """Let the public 200 ms pan/zoom transition reach its final pixels."""

    page.wait_for_timeout(250)


def _download_svg(
    page: Page,
    download_dir: Path,
    expected_name: str,
    assert_svg: Any,
) -> dict[str, Any]:
    svg_button = page.get_by_role("button", name="SVG", exact=True)
    with page.expect_download(timeout=ACTION_TIMEOUT_MS) as download_info:
        svg_button.click()
    download = download_info.value
    if download.failure() is not None:
        raise AssertionError(f"SVG download failed: {download.failure()}")
    if download.suggested_filename != expected_name:
        raise AssertionError(
            f"Expected {expected_name}, downloaded {download.suggested_filename}"
        )
    path = download_dir / expected_name
    download.save_as(path)
    if path.stat().st_size < 10_000:
        raise AssertionError(f"Downloaded SVG is unexpectedly small: {path.stat().st_size}")
    report = _inspect_downloaded_presentation_svg(path)
    assert_svg(report)
    return {
        "filename": path.name,
        "bytes": path.stat().st_size,
        "semantics": report,
    }


def capture_gui_styling(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> CaptureResult:
    """Run H-GUI-11 from raw GenBank input in one fresh browser context."""

    source_record = _verify_fixture_contracts(include_style_tables=True)
    assert_output_paths(output_paths, STYLE_SCREENSHOT_NAMES, "H-GUI-11")
    download_dir.mkdir(parents=True, exist_ok=True)
    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}

    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        _open_human_circular(page, STYLE_OUTPUT_PREFIX)
        species = page.get_by_label("Species", exact=True)
        species.fill("<i>Homo sapiens</i>")
        expect(species).to_have_value("<i>Homo sapiens</i>")
        separate_strands = page.get_by_label("Separate Strands", exact=True)
        separate_strands.check()
        expect(separate_strands).to_be_checked()

        colors_panel = page.get_by_label("Colors", exact=True)
        colors_panel.click()
        palette = page.get_by_label("Palette", exact=True)
        palette.select_option("soft_pastels")
        expect(palette).to_have_value("soft_pastels")
        page.get_by_label("Override File (-d)", exact=True).set_input_files(
            DEFAULT_COLOR_RULE_PATH
        )
        expect(
            page.get_by_role(
                "group", name="Override File (-d) selection", exact=True
            )
        ).to_contain_text(DEFAULT_COLOR_RULE_PATH.name)
        colors_panel.click()

        labels_panel = page.get_by_label("Labels", exact=True)
        labels_panel.click()
        label_mode = page.get_by_label("Label Mode", exact=True)
        label_mode.select_option("out")
        expect(label_mode).to_have_value("out")
        page.get_by_role("button", name="Whitelist", exact=True).click()
        page.get_by_role("button", name="+ Add Rule", exact=True).click()
        page.get_by_role("textbox", name="Feat", exact=True).fill("CDS")
        page.get_by_role("textbox", name="Qual", exact=True).fill("gene")
        page.get_by_role("textbox", name="Pattern", exact=True).fill(
            STYLE_LABEL_PATTERN
        )
        page.get_by_label("Priority File (TSV)", exact=True).set_input_files(
            GENE_PRIORITY_RULE_PATH
        )
        expect(
            page.get_by_role(
                "group", name="Priority File (TSV) selection", exact=True
            )
        ).to_contain_text(GENE_PRIORITY_RULE_PATH.name)
        labels_panel.click()
        _set_plot_title(page, STYLE_TITLE)

        final_report = generate_and_inspect(
            page, _inspect_presentation_svg, _assert_style_svg
        )
        _park_feature_search(page)
        zoom_out = page.get_by_role("button", name="Zoom out", exact=True)
        zoom_out.click()
        expect(
            page.get_by_role("button", name="Reset zoom", exact=True)
        ).to_contain_text("90%")
        _wait_for_preview_transform(page)
        labels_panel.click()
        page.get_by_role(
            "group", name="Priority File (TSV) selection", exact=True
        ).scroll_into_view_if_needed()
        page.get_by_role(
            "button", name="Whitelist", exact=True
        ).scroll_into_view_if_needed()
        screenshot_bytes["style-settings.png"] = capture_screenshot(
            page, output_paths["style-settings.png"], "Circular"
        )
        labels_panel.click()
        _fit_finished_circular_preview(page)
        _pan_preview_left(page)
        _wait_for_preview_transform(page)
        screenshot_bytes["style-result.png"] = capture_screenshot(
            page, output_paths["style-result.png"], "Circular"
        )
        download_report = _download_svg(
            page, download_dir, EXPECTED_STYLE_SVG, _assert_style_svg
        )
        capture.assert_clean()
    finally:
        capture.close()

    return CaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        source_record=source_record,
    )


def _set_visibility_rule(page: Page, index: int, rule: Mapping[str, str]) -> None:
    rule_number = index + 1
    record_id = page.get_by_label(
        f"Visibility rule {rule_number} record ID", exact=True
    )
    record_id.fill(rule["recordId"])
    record_id.blur()
    expect(record_id).to_have_value(rule["recordId"])
    page.get_by_label(
        f"Visibility rule {rule_number} feature type", exact=True
    ).select_option(rule["featureType"])
    qualifier = page.get_by_label(
        f"Visibility rule {rule_number} qualifier", exact=True
    )
    qualifier.fill(rule["qualifier"])
    qualifier.blur()
    expect(qualifier).to_have_value(rule["qualifier"])
    page.get_by_label(
        f"Visibility rule {rule_number} action", exact=True
    ).select_option(rule["action"])
    value_regex = page.get_by_label(
        f"Visibility rule {rule_number} value regex", exact=True
    )
    value_regex.fill(rule["value"])
    value_regex.blur()
    expect(value_regex).to_have_value(rule["value"])


def capture_gui_feature_presentation(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> CaptureResult:
    """Run H-GUI-12 from raw GenBank input in one fresh browser context."""

    source_record = _verify_fixture_contracts(include_style_tables=False)
    assert_output_paths(
        output_paths, PRESENTATION_SCREENSHOT_NAMES, "H-GUI-12"
    )
    download_dir.mkdir(parents=True, exist_ok=True)
    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}

    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        _open_human_circular(page, PRESENTATION_OUTPUT_PREFIX)
        separate_strands = page.get_by_label("Separate Strands", exact=True)
        separate_strands.uncheck()
        expect(separate_strands).not_to_be_checked()
        resolve_overlaps = page.get_by_label("Resolve Overlaps", exact=True)
        resolve_overlaps.check()
        expect(resolve_overlaps).to_be_checked()

        features_panel = page.get_by_label("Features", exact=True)
        features_panel.click()
        page.get_by_label("Rendering for CDS", exact=True).select_option("arrow")
        page.get_by_label("Rendering for rRNA", exact=True).select_option(
            "rectangle"
        )
        page.get_by_label("Rendering for tRNA", exact=True).select_option(
            "underlay"
        )
        head_ratio = page.get_by_label("Arrow head length ratio", exact=True)
        head_ratio.fill("1")
        expect(head_ratio).to_have_value("1")
        shaft_ratio = page.get_by_label("Arrow shaft width ratio", exact=True)
        shaft_ratio.fill("0.55")
        expect(shaft_ratio).to_have_value("0.55")

        add_visibility_rule = page.get_by_role(
            "button", name="Add feature visibility rule", exact=True
        )
        for index, rule in enumerate(EXPECTED_VISIBILITY_RULES):
            add_visibility_rule.click()
            _set_visibility_rule(page, index, rule)

        block_color_mode = page.get_by_label(
            "Block stroke color mode", exact=True
        )
        block_color_mode.select_option("color")
        block_color = page.get_by_label("Block stroke color", exact=True)
        block_color.fill("#263238")
        expect(block_color).to_have_value("#263238")
        block_width = page.get_by_label("Block Stroke Width", exact=True)
        block_width.fill("2.5")
        expect(block_width).to_have_value("2.5")
        line_color_mode = page.get_by_label("Line stroke color mode", exact=True)
        line_color_mode.select_option("color")
        line_color = page.get_by_label("Line stroke color", exact=True)
        line_color.fill("#455a64")
        expect(line_color).to_have_value("#455a64")
        line_width = page.get_by_label("Line Stroke Width", exact=True)
        line_width.fill("2")
        expect(line_width).to_have_value("2")

        _set_plot_title(page, PRESENTATION_TITLE)
        final_report = generate_and_inspect(
            page, _inspect_presentation_svg, _assert_feature_presentation_svg
        )
        _park_feature_search(page)
        page.get_by_label(
            "Visibility rule 1 record ID", exact=True
        ).scroll_into_view_if_needed()
        screenshot_bytes["presentation-settings.png"] = capture_screenshot(
            page, output_paths["presentation-settings.png"], "Circular"
        )
        features_panel.click()
        _fit_finished_circular_preview(page)
        _pan_preview_left(page, 0.06)
        _wait_for_preview_transform(page)
        screenshot_bytes["presentation-result.png"] = capture_screenshot(
            page, output_paths["presentation-result.png"], "Circular"
        )
        download_report = _download_svg(
            page,
            download_dir,
            EXPECTED_PRESENTATION_SVG,
            _assert_feature_presentation_svg,
        )
        capture.assert_clean()
    finally:
        capture.close()

    return CaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        source_record=source_record,
    )


__all__ = [
    "CaptureResult",
    "PRESENTATION_SCREENSHOT_NAMES",
    "STYLE_SCREENSHOT_NAMES",
    "capture_gui_feature_presentation",
    "capture_gui_styling",
]
