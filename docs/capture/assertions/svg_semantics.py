"""Read-only semantic checks for documentation-owned result SVGs."""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from decimal import Decimal
from pathlib import Path
from typing import Any


EXPECTED_GUI_LOSATN_MATCHES = (
    ("NC_001416.1", "NC_042057.1", "99.981", "21232", "2", "2", "1", "21231", "20081", "41311", "0.0", "39185"),
    ("NC_001416.1", "NC_042057.1", "99.922", "6412", "4", "1", "27971", "34381", "4697", "11108", "0.0", "11812"),
    ("NC_001416.1", "NC_042057.1", "99.866", "5205", "7", "0", "43298", "48502", "14876", "20080", "0.0", "9574"),
    ("NC_001416.1", "NC_042057.1", "99.634", "1914", "4", "3", "38591", "40501", "12972", "14885", "0.0", "3493"),
    ("NC_001416.1", "NC_042057.1", "100.000", "1620", "0", "0", "26104", "27723", "41306", "42925", "0.0", "2992"),
    ("NC_001416.1", "NC_042057.1", "100.000", "254", "0", "0", "27724", "27977", "1", "254", "6.07e-133", "470"),
)


_SVG_INSPECTION_SCRIPT = r"""
root => {
  const svgs = Array.from(root.getElementsByTagName('svg'));
  if (svgs.length !== 1) {
    return { svgCount: svgs.length };
  }

  const svg = svgs[0];
  const elements = [svg, ...Array.from(svg.getElementsByTagName('*'))];
  const texts = Array.from(svg.getElementsByTagName('text'))
    .map((node) => String(node.textContent || '').replace(/\s+/g, ' ').trim())
    .filter(Boolean);
  const ids = elements
    .map((element) => String(element.getAttribute('id') || ''))
    .filter(Boolean);
  const recordIds = elements
    .map((element) => String(element.getAttribute('data-gbdraw-record-id') || ''))
    .filter(Boolean);
  const featureElements = elements.filter((element) => (
    Boolean(element.getAttribute('data-gbdraw-feature-id'))
  ));
  const featureGeometry = featureElements.map((element) => {
    const bounds = element.getBBox();
    return {
      height: Number(bounds.height),
      strokeWidth: String(element.getAttribute('stroke-width') || ''),
    };
  });
  const matchElements = elements.filter((element) => (
    Boolean(element.getAttribute('data-gbdraw-pairwise-match-id')) &&
    String(element.getAttribute('data-match-kind') || '').toLowerCase() === 'pairwise'
  ));
  const comparisonElements = elements.filter((element) => (
    Boolean(element.getAttribute('data-gbdraw-pairwise-match-id'))
  ));
  const elementById = (id) => elements.find((element) => element.getAttribute('id') === id);
  const childCountById = (id) => Array.from(elementById(id)?.children || []).length;
  const labelTextsWithin = (element) => [
    ...Array.from(element.querySelectorAll(
      'g[id^="label_text"] text, text[text-anchor="start"][transform^="rotate("]'
    )),
    ...Array.from(element.getElementsByTagName('textPath')).filter((node) => {
      const href = String(
        node.getAttribute('href') || node.getAttribute('xlink:href') || ''
      );
      return href.startsWith('#circular_feature_label_path_');
    })
  ]
    .map((node) => String(node.textContent || '').replace(/\s+/g, ' ').trim())
    .filter(Boolean);
  const labelTexts = labelTextsWithin(svg);
  const parseTranslate = (element) => {
    const match = String(element?.getAttribute('transform') || '')
      .match(/translate\(\s*([-+\d.eE]+)[,\s]+([-+\d.eE]+)\s*\)/);
    return match ? [Number(match[1]), Number(match[2])] : null;
  };
  const recordPlacementCandidates = elements
    .filter((element) => (
      String(element.tagName || '').toLowerCase() === 'g' &&
      /^(?:record_\d+|record_group_[\w-]+)$/.test(String(element.getAttribute('id') || '')) &&
      Boolean(element.getAttribute('data-gbdraw-record-id'))
    ));
  const indexedRecordPlacements = recordPlacementCandidates.filter((element) => (
    /^record_\d+$/.test(String(element.getAttribute('id') || ''))
  ));
  const recordPlacements = (
    indexedRecordPlacements.length > 0
      ? indexedRecordPlacements
      : recordPlacementCandidates.filter((element) => (
          /^record_group_[\w-]+$/.test(String(element.getAttribute('id') || ''))
        ))
  )
    .map((element) => {
      const bounds = element.getBoundingClientRect();
      return {
        id: String(element.getAttribute('id') || ''),
        recordId: String(element.getAttribute('data-gbdraw-record-id') || ''),
        translate: parseTranslate(element),
        bounds: {
          x: Number(bounds.x),
          y: Number(bounds.y),
          width: Number(bounds.width),
          height: Number(bounds.height),
          right: Number(bounds.right),
          bottom: Number(bounds.bottom),
        },
        labelTexts: labelTextsWithin(element),
      };
    });
  const recordGroup = elements.find((element) => (
    String(element.tagName || '').toLowerCase() === 'g' &&
    element.getAttribute('data-gbdraw-record-id') === 'NC_012920.1'
  ));
  const linearRecordGroup = elements.find((element) => (
    String(element.tagName || '').toLowerCase() === 'g' &&
    element.getAttribute('data-gbdraw-record-id') === 'NC_001416.1'
  ));
  const eventAttributes = [];
  const unsafeHrefs = [];

  for (const element of elements) {
    for (const attribute of Array.from(element.attributes || [])) {
      const name = String(attribute.name || '').toLowerCase();
      const value = String(attribute.value || '').trim();
      if (name.startsWith('on')) {
        eventAttributes.push(`${name}=${value}`);
      }
      if ((name === 'href' || name === 'xlink:href') && /^(?:javascript:|https?:|\/\/|data:text\/html)/i.test(value)) {
        unsafeHrefs.push(value);
      }
    }
  }

  return {
    svgCount: svgs.length,
    ids,
    recordIds,
    featureElementCount: featureElements.length,
    featureHeights: featureGeometry.map((entry) => entry.height),
    featureStrokeWidths: featureGeometry.map((entry) => entry.strokeWidth),
    featureIds: featureElements.map((element) => (
      String(element.getAttribute('data-gbdraw-feature-id') || '')
    )),
    matches: matchElements.map((element) => ({
      query: String(element.getAttribute('data-query-record-id') || ''),
      subject: String(element.getAttribute('data-subject-record-id') || ''),
      identity: String(element.getAttribute('data-identity') || ''),
      alignmentLength: String(element.getAttribute('data-alignment-length') || ''),
      mismatches: String(element.getAttribute('data-mismatches') || ''),
      gapOpens: String(element.getAttribute('data-gap-opens') || ''),
      qstart: String(element.getAttribute('data-qstart') || ''),
      qend: String(element.getAttribute('data-qend') || ''),
      sstart: String(element.getAttribute('data-sstart') || ''),
      send: String(element.getAttribute('data-send') || ''),
      evalue: String(element.getAttribute('data-evalue') || ''),
      bitscore: String(element.getAttribute('data-bitscore') || ''),
    })),
    comparisonMatches: comparisonElements.map((element) => ({
      id: String(element.getAttribute('data-gbdraw-pairwise-match-id') || ''),
      kind: String(element.getAttribute('data-match-kind') || ''),
      query: String(element.getAttribute('data-query-record-id') || ''),
      subject: String(element.getAttribute('data-subject-record-id') || ''),
      orthogroupId: String(element.getAttribute('data-orthogroup-id') || ''),
      blockId: String(element.getAttribute('data-collinearity-block-id') || ''),
      anchorCount: String(element.getAttribute('data-collinearity-anchor-count') || ''),
      orientation: String(element.getAttribute('data-collinearity-orientation') || ''),
      colorMode: String(element.getAttribute('data-collinearity-color-mode') || ''),
      groupScope: String(element.getAttribute('data-collinear-group-scope') || ''),
    })),
    texts,
    labelTexts,
    italicTexts: elements
      .filter((element) => (
        String(element.getAttribute('font-style') || '').toLowerCase() === 'italic' ||
        /font-style\s*:\s*italic/i.test(String(element.getAttribute('style') || ''))
      ))
      .map((element) => String(element.textContent || '').replace(/\s+/g, ' ').trim())
      .filter(Boolean),
    coordinateTicks: texts.filter((text) => /^\d+\s+kbp$/.test(text)),
    groupChildCounts: {
      label_leaders: childCountById('label_leaders'),
      label_text: childCountById('label_text'),
      legend: childCountById('legend'),
    },
    recordTranslate: parseTranslate(recordGroup),
    linearRecordTranslate: parseTranslate(linearRecordGroup),
    recordPlacements,
    legendTranslate: parseTranslate(elementById('legend')),
    rootTag: String(svg.localName || svg.tagName || '').toLowerCase(),
    interactiveSvg: svg.getAttribute('data-gbdraw-interactive-svg') === 'true',
    scriptCount: svg.getElementsByTagName('script').length,
    eventAttributes,
    unsafeHrefs,
  };
}
"""

_RENDERED_FEATURE_STATE_SCRIPT = r"""
() => {
  const features = Array.isArray(window.__GBDRAW_APP__?.extractedFeatures)
    ? window.__GBDRAW_APP__.extractedFeatures
    : [];
  return features.map((feature) => ({
    type: String(feature?.type || ''),
    strand: String(feature?.strand || ''),
    recordId: String(feature?.record_id || ''),
    featureId: String(
      feature?.rendered_feature_svg_id || feature?.svg_id || ''
    )
  }));
}
"""


def _local_name(name: str) -> str:
    return name.rsplit("}", 1)[-1].lower()


def _normalized_text(element: ET.Element) -> str:
    return " ".join("".join(element.itertext()).split())


def _translate(element: ET.Element | None) -> list[float] | None:
    if element is None:
        return None
    match = re.search(
        r"translate\(\s*([-+\d.eE]+)[,\s]+([-+\d.eE]+)\s*\)",
        element.attrib.get("transform", ""),
    )
    return [float(match.group(1)), float(match.group(2))] if match else None


def _label_texts_within(element: ET.Element) -> list[str]:
    grouped_labels = [
        _normalized_text(text)
        for group in element.iter()
        if _local_name(group.tag) == "g"
        and group.attrib.get("id", "").startswith("label_text")
        for text in group.iter()
        if _local_name(text.tag) == "text" and _normalized_text(text)
    ]
    rotated_linear_labels = [
        _normalized_text(text)
        for text in element.iter()
        if _local_name(text.tag) == "text"
        and text.attrib.get("text-anchor") == "start"
        and text.attrib.get("transform", "").startswith("rotate(")
        and _normalized_text(text)
    ]
    curved_circular_labels = [
        _normalized_text(text_path)
        for text_path in element.iter()
        if _local_name(text_path.tag) == "textpath"
        and any(
            _local_name(name) == "href"
            and value.startswith("#circular_feature_label_path_")
            for name, value in text_path.attrib.items()
        )
        and _normalized_text(text_path)
    ]
    return list(
        dict.fromkeys(
            [
                *grouped_labels,
                *rotated_linear_labels,
                *curved_circular_labels,
            ]
        )
    )


def inspect_svg_file(path: Path) -> dict[str, Any]:
    """Collect the same semantic report from a downloaded SVG file."""

    root = ET.parse(path).getroot()
    elements = list(root.iter())
    texts = [
        _normalized_text(element)
        for element in elements
        if _local_name(element.tag) == "text" and _normalized_text(element)
    ]
    ids = [element.attrib["id"] for element in elements if element.attrib.get("id")]
    by_id = {element.attrib["id"]: element for element in elements if element.attrib.get("id")}
    feature_elements = [
        element for element in elements if element.attrib.get("data-gbdraw-feature-id")
    ]
    match_elements = [
        element
        for element in elements
        if element.attrib.get("data-gbdraw-pairwise-match-id")
        and element.attrib.get("data-match-kind", "").lower() == "pairwise"
    ]
    comparison_elements = [
        element
        for element in elements
        if element.attrib.get("data-gbdraw-pairwise-match-id")
    ]
    record_group = next(
        (
            element
            for element in elements
            if _local_name(element.tag) == "g"
            and element.attrib.get("data-gbdraw-record-id") == "NC_012920.1"
        ),
        None,
    )
    linear_record_group = next(
        (
            element
            for element in elements
            if _local_name(element.tag) == "g"
            and element.attrib.get("data-gbdraw-record-id") == "NC_001416.1"
        ),
        None,
    )
    record_container_candidates = [
        element
        for element in elements
        if _local_name(element.tag) == "g"
        and re.fullmatch(
            r"(?:record_\d+|record_group_[\w-]+)",
            element.attrib.get("id", ""),
        )
        and element.attrib.get("data-gbdraw-record-id")
    ]
    indexed_record_containers = [
        element
        for element in record_container_candidates
        if re.fullmatch(r"record_\d+", element.attrib.get("id", ""))
    ]
    record_containers = indexed_record_containers or [
        element
        for element in record_container_candidates
        if re.fullmatch(r"record_group_[\w-]+", element.attrib.get("id", ""))
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
        "ids": ids,
        "recordIds": [
            element.attrib["data-gbdraw-record-id"]
            for element in elements
            if element.attrib.get("data-gbdraw-record-id")
        ],
        "featureElementCount": len(feature_elements),
        "featureHeights": [
            float(element.attrib.get("data-gbdraw-feature-height", "nan"))
            for element in feature_elements
            if element.attrib.get("data-gbdraw-feature-height")
        ],
        "featureStrokeWidths": [
            element.attrib.get("stroke-width", "") for element in feature_elements
        ],
        "featureIds": [
            element.attrib["data-gbdraw-feature-id"] for element in feature_elements
        ],
        "matches": [
            {
                "query": element.attrib.get("data-query-record-id", ""),
                "subject": element.attrib.get("data-subject-record-id", ""),
                "identity": element.attrib.get("data-identity", ""),
                "alignmentLength": element.attrib.get("data-alignment-length", ""),
                "mismatches": element.attrib.get("data-mismatches", ""),
                "gapOpens": element.attrib.get("data-gap-opens", ""),
                "qstart": element.attrib.get("data-qstart", ""),
                "qend": element.attrib.get("data-qend", ""),
                "sstart": element.attrib.get("data-sstart", ""),
                "send": element.attrib.get("data-send", ""),
                "evalue": element.attrib.get("data-evalue", ""),
                "bitscore": element.attrib.get("data-bitscore", ""),
            }
            for element in match_elements
        ],
        "comparisonMatches": [
            {
                "id": element.attrib.get("data-gbdraw-pairwise-match-id", ""),
                "kind": element.attrib.get("data-match-kind", ""),
                "query": element.attrib.get("data-query-record-id", ""),
                "subject": element.attrib.get("data-subject-record-id", ""),
                "orthogroupId": element.attrib.get("data-orthogroup-id", ""),
                "blockId": element.attrib.get("data-collinearity-block-id", ""),
                "anchorCount": element.attrib.get(
                    "data-collinearity-anchor-count", ""
                ),
                "orientation": element.attrib.get(
                    "data-collinearity-orientation", ""
                ),
                "colorMode": element.attrib.get(
                    "data-collinearity-color-mode", ""
                ),
                "groupScope": element.attrib.get("data-collinear-group-scope", ""),
            }
            for element in comparison_elements
        ],
        "texts": texts,
        "labelTexts": _label_texts_within(root),
        "italicTexts": [
            _normalized_text(element)
            for element in elements
            if (
                element.attrib.get("font-style", "").lower() == "italic"
                or re.search(r"font-style\s*:\s*italic", element.attrib.get("style", ""), re.I)
            )
            and _normalized_text(element)
        ],
        "coordinateTicks": [text for text in texts if re.fullmatch(r"\d+\s+kbp", text)],
        "groupChildCounts": {
            group_id: len(list(by_id[group_id])) if group_id in by_id else 0
            for group_id in ("label_leaders", "label_text", "legend")
        },
        "recordTranslate": _translate(record_group),
        "linearRecordTranslate": _translate(linear_record_group),
        "recordPlacements": [
            {
                "id": element.attrib.get("id", ""),
                "recordId": element.attrib.get("data-gbdraw-record-id", ""),
                "translate": _translate(element),
                "bounds": None,
                "labelTexts": _label_texts_within(element),
            }
            for element in record_containers
        ],
        "legendTranslate": _translate(by_id.get("legend")),
        "rootTag": _local_name(root.tag),
        "interactiveSvg": root.attrib.get("data-gbdraw-interactive-svg") == "true",
        "scriptCount": sum(_local_name(element.tag) == "script" for element in elements),
        "eventAttributes": event_attributes,
        "unsafeHrefs": unsafe_hrefs,
    }


def inspect_first_circular_svg(result_region: Any) -> dict[str, Any]:
    """Inspect the rendered SVG without changing page or application state."""

    report = result_region.evaluate(_SVG_INSPECTION_SCRIPT)
    if not isinstance(report, dict):
        raise AssertionError("SVG inspection did not return a report")
    return report


def inspect_first_linear_svg(result_region: Any) -> dict[str, Any]:
    """Inspect the rendered Linear SVG without changing page or application state."""

    report = result_region.evaluate(_SVG_INSPECTION_SCRIPT)
    if not isinstance(report, dict):
        raise AssertionError("SVG inspection did not return a report")
    return report


def inspect_lambda_input_svg(result_region: Any, page: Any) -> dict[str, Any]:
    """Inspect the rendered Lambda SVG and its generated feature metadata."""

    report = inspect_first_linear_svg(result_region)
    feature_metadata = page.evaluate(_RENDERED_FEATURE_STATE_SCRIPT)
    if not isinstance(feature_metadata, list):
        raise AssertionError("Rendered feature inspection did not return a list")

    cds_features = [
        feature
        for feature in feature_metadata
        if isinstance(feature, dict) and feature.get("type") == "CDS"
    ]
    report["metadataFeatureIds"] = [
        str(feature.get("featureId") or "") for feature in cds_features
    ]
    report["metadataRecordIds"] = [
        str(feature.get("recordId") or "") for feature in cds_features
    ]
    report["cdsStrandCounts"] = {
        strand: sum(feature.get("strand") == strand for feature in cds_features)
        for strand in ("+", "-", "undefined")
    }
    return report


def inspect_gui_losatn_svg(result_region: Any) -> dict[str, Any]:
    """Inspect the two-record Lambda/DE3 Linear comparison."""

    return inspect_first_linear_svg(result_region)


def inspect_gui_bgc_losatp_svg(result_region: Any) -> dict[str, Any]:
    """Inspect the five-record BGC LOSATP result."""

    return inspect_first_linear_svg(result_region)


def _match_tuple(match: dict[str, Any]) -> tuple[str, ...]:
    return tuple(
        str(match.get(key, ""))
        for key in (
            "query",
            "subject",
            "identity",
            "alignmentLength",
            "mismatches",
            "gapOpens",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        )
    )


def _numeric_match_signature(match: tuple[str, ...]) -> tuple[Any, ...]:
    return (*match[:2], *(Decimal(value) for value in match[2:]))


def assert_plain_gui_losatn_svg(report: dict[str, Any]) -> None:
    """Verify both complete input records before any comparison is enabled."""

    if report.get("svgCount") != 1:
        raise AssertionError(f"Expected one result SVG, found {report.get('svgCount')}")
    if set(report.get("recordIds", [])) != {"NC_001416.1", "NC_042057.1"}:
        raise AssertionError(
            "Plain Linear result does not contain the complete Lambda and DE3 records"
        )
    if report.get("featureElementCount") != 130:
        raise AssertionError(
            "Expected 73 Lambda and 57 DE3 CDS feature elements, found "
            f"{report.get('featureElementCount')}"
        )
    for expected_text in (
        "NC_001416.1",
        "48,502 bp",
        "NC_042057.1",
        "42,925 bp",
    ):
        if expected_text not in set(report.get("texts", [])):
            raise AssertionError(f"Missing whole-record SVG text: {expected_text}")
    if report.get("matches"):
        raise AssertionError("The Step 2 Linear map must not contain comparison links")
    assert_static_svg_safety(report)


def assert_gui_losatn_svg(report: dict[str, Any]) -> None:
    """Verify all six qualified LOSATN links and their whole-record endpoints."""

    assert_plain_gui_losatn_svg({**report, "matches": []})
    matches = tuple(_match_tuple(match) for match in report.get("matches", []))
    actual_signatures = sorted(_numeric_match_signature(match) for match in matches)
    expected_signatures = sorted(
        _numeric_match_signature(match) for match in EXPECTED_GUI_LOSATN_MATCHES
    )
    if actual_signatures != expected_signatures:
        raise AssertionError(
            "Rendered LOSATN link metadata does not match the qualified six-row TSV: "
            f"{matches!r}"
        )

    for match in matches:
        qstart, qend = int(match[6]), int(match[7])
        sstart, send = int(match[8]), int(match[9])
        if not (1 <= min(qstart, qend) <= max(qstart, qend) <= 48_502):
            raise AssertionError(f"LOSATN query endpoint is outside whole Lambda: {match}")
        if not (1 <= min(sstart, send) <= max(sstart, send) <= 42_925):
            raise AssertionError(f"LOSATN subject endpoint is outside whole DE3: {match}")


_BGC_RECORD_LENGTHS = {
    "BGC0000708": "40,579 bp",
    "BGC0000709": "50,466 bp",
    "BGC0000711": "30,837 bp",
    "BGC0000712": "48,169 bp",
    "BGC0000713": "31,892 bp",
}
_BGC_ADJACENT_PAIRS = {
    ("BGC0000708", "BGC0000709"),
    ("BGC0000709", "BGC0000711"),
    ("BGC0000711", "BGC0000712"),
    ("BGC0000712", "BGC0000713"),
}
_BGC_GALLERY_GENE_LABELS = {
    "livA",
    "livB",
    "livC",
    "livD",
    "livE",
    "livF",
    "livG",
    "livH",
    "livI",
    "livK",
    "livL",
    "livM",
    "livN",
    "livO",
    "livP",
    "livQ",
    "livS",
    "livT",
    "livU",
    "livV",
    "livW",
    "livX",
    "livY",
    "livZ",
}


def _assert_bgc_gallery_presentation(report: dict[str, Any]) -> None:
    """Require the public BGC recipe's readable feature and label geometry."""

    feature_heights = [float(value) for value in report.get("featureHeights", [])]
    if len(feature_heights) != 155 or any(
        abs(value - 37.5) > 0.01 for value in feature_heights
    ):
        raise AssertionError(
            "BGC feature glyphs do not use the Gallery 75 px feature-height recipe: "
            f"{sorted(set(feature_heights))!r}"
        )
    if set(report.get("featureStrokeWidths", [])) != {"2.0"}:
        raise AssertionError(
            "BGC feature outlines do not use the Gallery 2 px stroke: "
            f"{sorted(set(report.get('featureStrokeWidths', [])))}"
        )

    label_texts = set(report.get("labelTexts", []))
    if label_texts != _BGC_GALLERY_GENE_LABELS:
        raise AssertionError(
            "BGC labels must match the Gallery's first-record gene labels exactly: "
            f"{sorted(label_texts)!r}"
        )
    placements = {
        str(placement.get("recordId") or ""): set(
            placement.get("labelTexts", [])
        )
        for placement in report.get("recordPlacements", [])
    }
    if placements.get("BGC0000708") != _BGC_GALLERY_GENE_LABELS or any(
        placements.get(record_id, set())
        for record_id in ("BGC0000709", "BGC0000711", "BGC0000712", "BGC0000713")
    ):
        raise AssertionError(
            "BGC labels must appear on the first record only: "
            f"{ {record_id: sorted(labels) for record_id, labels in placements.items()}!r}"
        )

    texts = set(report.get("texts", []))
    expected_record_text = {
        "Streptomyces lividus CBS 844.73",
        "Streptomyces fradiae ATCC 10745",
        "Streptomyces fradiae MCIMB 8233",
        "Streptomyces rimosus subsp. paromomycinus NRRL 2455",
        "Streptomyces ribosidificus ATCC 21294",
    }
    if not expected_record_text.issubset(texts):
        raise AssertionError("BGC record labels do not match the Gallery presentation")
    if not {"10 kbp", "20 kbp", "30 kbp", "40 kbp", "50 kbp"}.issubset(
        set(report.get("coordinateTicks", []))
    ):
        raise AssertionError("The Gallery-style BGC coordinate ruler is incomplete")


def assert_plain_gui_bgc_svg(report: dict[str, Any]) -> None:
    """Verify the five complete linear BGC records before LOSATP runs."""

    if report.get("svgCount") != 1:
        raise AssertionError(f"Expected one result SVG, found {report.get('svgCount')}")
    record_ids = set(report.get("recordIds", []))
    if record_ids != set(_BGC_RECORD_LENGTHS):
        raise AssertionError(
            "BGC result does not contain exactly the five complete source records: "
            f"{sorted(record_ids)}"
        )
    if report.get("featureElementCount") != 155:
        raise AssertionError(
            "Expected 155 CDS feature elements across the five BGC records, found "
            f"{report.get('featureElementCount')}"
        )
    texts = set(report.get("texts", []))
    for record_id, length_text in _BGC_RECORD_LENGTHS.items():
        if record_id not in texts or length_text not in texts:
            raise AssertionError(
                f"Missing complete BGC record text: {record_id} / {length_text}"
            )
    if report.get("comparisonMatches"):
        raise AssertionError("The plain five-record BGC map contains comparison links")
    assert_static_svg_safety(report)


def assert_gui_bgc_similarity_groups_svg(report: dict[str, Any]) -> None:
    """Verify the canonical five-record LOSATP Similarity groups result."""

    assert_plain_gui_bgc_svg({**report, "comparisonMatches": []})
    _assert_bgc_gallery_presentation(report)
    matches = report.get("comparisonMatches", [])
    if len(matches) != 77 or {match.get("kind") for match in matches} != {
        "orthogroup"
    }:
        raise AssertionError(
            "Expected 77 LOSATP Similarity-group links, found "
            f"{len(matches)} links with kinds "
            f"{sorted({str(match.get('kind')) for match in matches})}"
        )
    endpoints = {(match.get("query"), match.get("subject")) for match in matches}
    if endpoints != _BGC_ADJACENT_PAIRS:
        raise AssertionError(f"Unexpected Similarity-group endpoints: {endpoints!r}")
    if ("BGC0000708", "BGC0000713") in endpoints:
        raise AssertionError("The canonical BGC result contains a direct 0708 to 0713 edge")
    linked_group_ids = {str(match.get("orthogroupId")) for match in matches}
    expected_group_ids = {
        *(f"og_{index}" for index in range(1, 14)),
        "og_15",
        "og_16",
        "og_17",
        "og_18",
        "og_19",
        "og_21",
        "og_22",
        "og_23",
    }
    if linked_group_ids != expected_group_ids:
        raise AssertionError(
            "Similarity-group IDs changed: "
            f"expected {sorted(expected_group_ids)}, found {sorted(linked_group_ids)}"
        )


def assert_gui_bgc_collinear_svg(report: dict[str, Any]) -> None:
    """Verify the separate five-record LOSATP Collinear block result."""

    assert_plain_gui_bgc_svg({**report, "comparisonMatches": []})
    _assert_bgc_gallery_presentation(report)
    matches = report.get("comparisonMatches", [])
    if len(matches) != 7 or {match.get("kind") for match in matches} != {
        "collinear"
    }:
        raise AssertionError(
            f"Expected seven Collinear blocks, found {len(matches)}: {matches!r}"
        )
    endpoints = {(match.get("query"), match.get("subject")) for match in matches}
    if endpoints != _BGC_ADJACENT_PAIRS:
        raise AssertionError(f"Unexpected Collinear endpoints: {endpoints!r}")
    expected = {
        "block_0001": ("13", "plus"),
        "block_0002": ("3", "minus"),
        "block_0003": ("21", "plus"),
        "block_0004": ("2", "plus"),
        "block_0005": ("15", "plus"),
        # BGC0000713 is reverse-complemented for the Gallery layout, so the
        # two raw minus-orientation blocks are plus in displayed coordinates.
        "block_0006": ("13", "plus"),
        "block_0007": ("2", "plus"),
    }
    actual = {
        str(match.get("blockId")): (
            str(match.get("anchorCount")),
            str(match.get("orientation")),
        )
        for match in matches
    }
    if actual != expected:
        raise AssertionError(
            f"Collinear block IDs, anchors, or orientations changed: {actual!r}"
        )
    if {match.get("colorMode") for match in matches} != {"orientation_identity"}:
        raise AssertionError("Collinear blocks are not using Orientation + identity")


def assert_gui_circular_layout_svg(report: dict[str, Any]) -> None:
    """Verify the four complete mitochondrial records and their 2 by 2 grid."""

    if report.get("svgCount") != 1:
        raise AssertionError(f"Expected one result SVG, found {report.get('svgCount')}")

    expected_records = {
        "NC_012920.1": ("16,569 bp", "Homo sapiens"),
        "NC_002333.2": ("16,596 bp", "Danio rerio"),
        "NC_024511.2": ("19,524 bp", "Drosophila melanogaster"),
        "NC_001328.1": ("13,794 bp", "Caenorhabditis elegans"),
    }
    record_ids = set(report.get("recordIds", []))
    if record_ids != set(expected_records):
        raise AssertionError(
            "Multi-record Circular result does not contain exactly the four complete "
            f"mitochondrial RefSeq records: {sorted(record_ids)}"
        )
    if report.get("featureElementCount") != 147:
        raise AssertionError(
            "Expected 147 mitochondrial CDS, tRNA, and rRNA feature elements, "
            f"found {report.get('featureElementCount')}"
        )

    texts = set(report.get("texts", []))
    for record_id, (length_text, organism) in expected_records.items():
        for expected_text in (record_id, length_text, organism):
            if expected_text not in texts:
                raise AssertionError(
                    f"Missing multi-record SVG text: {expected_text}; "
                    f"found {sorted(texts)!r}"
                )
    if not {"GC content", "GC skew (+)", "GC skew (-)"}.issubset(texts):
        raise AssertionError("The multi-record Circular SVG is missing quantitative tracks")
    if "Complete metazoan mitochondrial genomes" not in texts:
        raise AssertionError("The multi-record Circular SVG is missing its shared title")

    placements = report.get("recordPlacements", [])
    if len(placements) != 4 or {
        placement.get("recordId") for placement in placements
    } != set(expected_records):
        raise AssertionError(f"Unexpected Circular record placements: {placements!r}")
    translations = [placement.get("translate") for placement in placements]
    if any(not translate or len(translate) != 2 for translate in translations):
        raise AssertionError(f"Missing Circular record translations: {translations!r}")

    shared_cds_genes = {
        "ATP6",
        "COX1",
        "COX2",
        "COX3",
        "CYTB",
        "ND1",
        "ND2",
        "ND3",
        "ND4",
        "ND4L",
        "ND5",
        "ND6",
    }
    expected_cds_genes = {
        "NC_012920.1": shared_cds_genes | {"ATP8"},
        "NC_002333.2": shared_cds_genes | {"ATP8"},
        "NC_024511.2": shared_cds_genes | {"ATP8"},
        "NC_001328.1": shared_cds_genes,
    }
    forbidden_cds_products = {
        "ATP synthase F0 subunit 6",
        "ATP synthase F0 subunit 8",
        "NADH dehydrogenase subunit 1",
        "NADH dehydrogenase subunit 2",
        "NADH dehydrogenase subunit 3",
        "NADH dehydrogenase subunit 4",
        "NADH dehydrogenase subunit 4L",
        "NADH dehydrogenase subunit 5",
        "NADH dehydrogenase subunit 6",
        "cytochrome b",
        "cytochrome c oxidase subunit I",
        "cytochrome c oxidase subunit II",
        "cytochrome c oxidase subunit III",
    }
    for placement in placements:
        record_id = str(placement.get("recordId") or "")
        label_texts = set(placement.get("labelTexts", []))
        missing = expected_cds_genes[record_id] - label_texts
        if missing:
            raise AssertionError(
                f"{record_id} is missing CDS gene labels: {sorted(missing)!r}"
            )
        products = forbidden_cds_products & label_texts
        if products:
            raise AssertionError(
                f"{record_id} uses CDS product text instead of gene labels: "
                f"{sorted(products)!r}"
            )
    x_coordinates = sorted(float(translate[0]) for translate in translations)
    y_coordinates = sorted(float(translate[1]) for translate in translations)
    if (
        len(x_coordinates) != 4
        or abs(x_coordinates[0] - x_coordinates[1]) > 25
        or abs(x_coordinates[2] - x_coordinates[3]) > 25
        or abs(x_coordinates[1] - x_coordinates[2]) < 100
        or len(y_coordinates) != 4
        or abs(y_coordinates[0] - y_coordinates[1]) > 1
        or abs(y_coordinates[2] - y_coordinates[3]) > 1
        or abs(y_coordinates[1] - y_coordinates[2]) < 100
    ):
        raise AssertionError(
            "Expected a 2 by 2 Circular grid; found translations "
            f"{translations!r}"
        )

    bounds = [placement.get("bounds") for placement in placements]
    if all(isinstance(bound, dict) for bound in bounds):
        for index, left in enumerate(bounds):
            for right in bounds[index + 1 :]:
                overlap_width = min(float(left["right"]), float(right["right"])) - max(
                    float(left["x"]), float(right["x"])
                )
                overlap_height = min(
                    float(left["bottom"]), float(right["bottom"])
                ) - max(float(left["y"]), float(right["y"]))
                if overlap_width > 0.5 and overlap_height > 0.5:
                    raise AssertionError(
                        "Circular record bounds overlap on the shared grid: "
                        f"{left!r} and {right!r}"
                    )

    assert_static_svg_safety(report)


def assert_static_svg_safety(report: dict[str, Any]) -> None:
    """Reject active or externally linked content in a standard SVG."""

    if report.get("interactiveSvg"):
        raise AssertionError("Standard SVG is marked as interactive")
    if report.get("scriptCount") != 0:
        raise AssertionError("Static preview SVG contains a script element")
    if report.get("eventAttributes"):
        raise AssertionError(
            f"Static preview SVG contains event handlers: {report['eventAttributes']}"
        )
    if report.get("unsafeHrefs"):
        raise AssertionError(
            f"Static preview SVG contains unsafe links: {report['unsafeHrefs']}"
        )


def assert_first_circular_svg(report: dict[str, Any]) -> None:
    """Assert the biological, track, coordinate, and safety contract."""

    if report.get("svgCount") != 1:
        raise AssertionError(f"Expected one result SVG, found {report.get('svgCount')}")

    record_ids = set(report.get("recordIds", []))
    if "NC_012920.1" not in record_ids:
        raise AssertionError(f"Missing NC_012920.1 record semantics: {sorted(record_ids)}")

    feature_count = report.get("featureElementCount")
    if feature_count != 37:
        raise AssertionError(f"Expected 37 mitochondrial feature elements, found {feature_count}")

    ids = set(report.get("ids", []))
    for expected_id in ("gc_content", "skew"):
        if expected_id not in ids:
            raise AssertionError(f"Missing SVG track group: {expected_id}")
    if not ({"tick", "ticks"} & ids):
        raise AssertionError("Missing the Circular coordinate-tick group")

    texts = set(report.get("texts", []))
    for expected_text in (
        "NC_012920.1",
        "16,569 bp",
        "GC content",
        "GC skew (+)",
        "GC skew (-)",
    ):
        if expected_text not in texts:
            raise AssertionError(f"Missing SVG text: {expected_text}")

    coordinate_ticks = set(report.get("coordinateTicks", []))
    if not {"1 kbp", "16 kbp"}.issubset(coordinate_ticks):
        raise AssertionError(
            "Expected coordinate ticks from 1 kbp through 16 kbp; "
            f"found {sorted(coordinate_ticks)}"
        )

    assert_static_svg_safety(report)


def assert_first_linear_svg(report: dict[str, Any]) -> None:
    """Assert the whole-record Lambda identity, feature count, and SVG safety."""

    if report.get("svgCount") != 1:
        raise AssertionError(f"Expected one result SVG, found {report.get('svgCount')}")

    record_ids = set(report.get("recordIds", []))
    if "NC_001416.1" not in record_ids:
        raise AssertionError(f"Missing NC_001416.1 record semantics: {sorted(record_ids)}")

    feature_count = report.get("featureElementCount")
    if feature_count != 73:
        raise AssertionError(f"Expected 73 Lambda CDS feature elements, found {feature_count}")

    texts = set(report.get("texts", []))
    for expected_text in ("NC_001416.1", "48,502 bp", "CDS"):
        if expected_text not in texts:
            raise AssertionError(f"Missing SVG text: {expected_text}")

    assert_static_svg_safety(report)


def assert_lambda_input_svg(report: dict[str, Any]) -> None:
    """Assert whole-record Lambda identity, CDS strands, and rendered IDs."""

    assert_first_linear_svg(report)
    if set(report.get("metadataRecordIds", [])) != {"NC_001416.1"}:
        raise AssertionError(
            "Generated feature metadata does not belong to NC_001416.1"
        )
    strand_counts = report.get("cdsStrandCounts", {})
    if strand_counts != {"+": 47, "-": 26, "undefined": 0}:
        raise AssertionError(
            "Expected 47 positive-strand and 26 negative-strand Lambda CDS "
            f"features, found {strand_counts}"
        )

    rendered_ids = sorted(str(value) for value in report.get("featureIds", []))
    metadata_ids = sorted(
        str(value) for value in report.get("metadataFeatureIds", [])
    )
    if len(metadata_ids) != 73 or metadata_ids != rendered_ids:
        raise AssertionError(
            "Rendered Lambda feature IDs do not match generated CDS metadata"
        )


def assert_matching_lambda_input_semantics(
    genbank_report: dict[str, Any],
    gff3_report: dict[str, Any],
) -> None:
    """Require GenBank and GFF3 + FASTA to render the same Lambda CDS set."""

    assert_lambda_input_svg(genbank_report)
    assert_lambda_input_svg(gff3_report)
    for key in ("recordIds", "featureIds", "metadataFeatureIds", "metadataRecordIds"):
        genbank_value = sorted(str(value) for value in genbank_report.get(key, []))
        gff3_value = sorted(str(value) for value in gff3_report.get(key, []))
        if genbank_value != gff3_value:
            raise AssertionError(
                f"GenBank and GFF3 + FASTA differ for {key}: "
                f"{genbank_value!r} != {gff3_value!r}"
            )

    for key in ("cdsStrandCounts",):
        genbank_value = genbank_report.get(key)
        gff3_value = gff3_report.get(key)
        if genbank_value != gff3_value:
            raise AssertionError(
                f"GenBank and GFF3 + FASTA differ for {key}: "
                f"{genbank_value!r} != {gff3_value!r}"
            )


def assert_finished_linear_svg(report: dict[str, Any]) -> None:
    """Verify concise labels, the coordinate ruler, and the left legend."""

    assert_first_linear_svg(report)
    texts = set(report.get("texts", []))

    coordinate_ticks = set(report.get("coordinateTicks", []))
    if not {"5 kbp", "45 kbp"}.issubset(coordinate_ticks):
        raise AssertionError(
            "Expected the whole-record ruler from 5 kbp through 45 kbp; "
            f"found {sorted(coordinate_ticks)}"
        )
    if not {"A", "B", "J", "int"}.issubset(texts):
        raise AssertionError("The final Linear SVG is missing concise Lambda gene labels")

    ids = set(report.get("ids", []))
    if "legend" not in ids or report.get("groupChildCounts", {}).get("legend", 0) < 1:
        raise AssertionError("The final Linear SVG is missing its CDS legend")

    record_translate = report.get("linearRecordTranslate")
    legend_translate = report.get("legendTranslate")
    if not record_translate or not legend_translate:
        raise AssertionError("Could not resolve Linear record and legend positions")
    if float(legend_translate[0]) >= float(record_translate[0]):
        raise AssertionError("The final legend is not positioned left of the Linear record")


def assert_species_markup(report: dict[str, Any]) -> None:
    """Verify that the species markup became italic SVG text."""

    if "Homo sapiens" not in set(report.get("texts", [])):
        raise AssertionError("Missing the Homo sapiens center label")
    if "Homo sapiens" not in set(report.get("italicTexts", [])):
        raise AssertionError("The Homo sapiens center label is not italicized")


def assert_finished_circular_svg(report: dict[str, Any]) -> None:
    """Verify the final labels and right-side legend in addition to base semantics."""

    assert_first_circular_svg(report)
    assert_species_markup(report)

    ids = set(report.get("ids", []))
    for expected_id in ("label_leaders", "label_text", "legend"):
        if expected_id not in ids:
            raise AssertionError(f"Missing final SVG group: {expected_id}")

    counts = report.get("groupChildCounts", {})
    if counts.get("label_leaders", 0) < 10 or counts.get("label_text", 0) < 10:
        raise AssertionError(f"External feature labels are incomplete: {counts}")
    if counts.get("legend", 0) < 1:
        raise AssertionError("The final SVG legend is empty")

    texts = set(report.get("texts", []))
    if not {"tRNA", "rRNA", "CDS"}.issubset(texts):
        raise AssertionError("The final SVG is missing feature legend entries")

    record_translate = report.get("recordTranslate")
    legend_translate = report.get("legendTranslate")
    if not record_translate or not legend_translate:
        raise AssertionError("Could not resolve record and legend positions")
    if float(legend_translate[0]) <= float(record_translate[0]):
        raise AssertionError(
            "The final legend is not positioned to the right of the Circular record"
        )
