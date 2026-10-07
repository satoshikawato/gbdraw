#!/usr/bin/env python3
"""Generate the shared Legend layout vectors (Python's results, exact).

``tests/fixtures/legend_layout_vectors.json`` holds inputs and Python's own
results for:

- text measurement (``calculate_bbox_dimensions``);
- the Legend of every Gallery Session Result, and edited copies of it (rows
  renamed, deleted, added, reordered; Circular Legends moved to each side);
- synthetic Legends for the branches the Gallery does not reach (conservation
  gradients, wrapping, an empty table).

``tests/test_legend_layout_vectors.py`` recomputes the expected values with
Python and ``tests/web/legend-layout-vectors.test.mjs`` with the JavaScript
port (``gbdraw/web/js/services/legend-layout.js``); both require equality.

Usage:
  python tools/generate_legend_layout_vectors.py            # recompute expected values
  python tools/generate_legend_layout_vectors.py --check    # fail when stale
  python tools/generate_legend_layout_vectors.py --from-gallery  # rebuild inputs, then expected
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
TARGET = REPO_ROOT / "tests" / "fixtures" / "legend_layout_vectors.json"
GALLERY_SESSIONS = REPO_ROOT / "gbdraw" / "web" / "gallery" / "sessions"
sys.path.insert(0, str(REPO_ROOT))

from gbdraw.configurators.legend import (  # noqa: E402
    _circular_legend_local_bounds,
    _linear_legend_local_bounds,
)
from gbdraw.core.text import _resolve_font_path, calculate_bbox_dimensions  # noqa: E402
from gbdraw.layout.composition import (  # noqa: E402
    OVERLAY_CANDIDATE_SCORE_ORDER,
    OVERLAY_CANVAS_GROWTH_CANDIDATE_ORDER,
    OVERLAY_CANVAS_GROWTH_SCORE_ORDER,
    OVERLAY_QUADRANT_BOUNDARY_RATIO,
    CompositionItem,
    CompositionRequest,
    CompositionSpacing,
    plan_composition,
)
from gbdraw.layout.spatial import Aabb  # noqa: E402
from gbdraw.legend.circular_layout import build_circular_legend_layout  # noqa: E402
from gbdraw.legend.linear_layout import build_linear_legend_layout  # noqa: E402

SVG = "{http://www.w3.org/2000/svg}"
NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
DPI = 96
LONG_CAPTION = "A much longer Legend caption for this row"
OVERLAY_POLICY = {
    "candidateScoreOrder": list(OVERLAY_CANDIDATE_SCORE_ORDER),
    "canvasGrowthCandidateOrder": list(OVERLAY_CANVAS_GROWTH_CANDIDATE_ORDER),
    "canvasGrowthScoreOrder": list(OVERLAY_CANVAS_GROWTH_SCORE_ORDER),
    "quadrantBoundaryRatio": OVERLAY_QUADRANT_BOUNDARY_RATIO,
}
MEASUREMENT_TEXTS = (
    "", "CDS", "tRNA", "GC content", "GC skew (+)", "other proteins", "AVATAR", "To Ty Ye",
    "Wolf", "Pairwise match identity", "rRNA [#ff0000]", "Core biosynthetic genes",
    "Café crème", "Δ-proteobacteria", "日本語", "emoji 🧬 row", "ﬁ ligature", "  spaced  ",
    "a b", "ab", "a中b", "中", "😀", "אל", "שלום עולם", "e\u0301 combining", "Ünïcödé Ωmega",
)
MEASUREMENT_FAMILIES = (
    "'Liberation Sans', 'Arial', 'Helvetica', 'Nimbus Sans L', sans-serif",
    "Times New Roman, serif",
    "Courier New",
    "Unknown Family",
)


# ---- inputs ----

def _font_file(family: str) -> str:
    path = _resolve_font_path(family)
    if not path:
        raise SystemExit(f"No bundled font for {family!r}")
    return Path(path).stem


def _translation(element: ET.Element | None) -> list[float]:
    x = y = 0.0
    if element is None:
        return [x, y]
    for first, second in re.findall(
        rf"translate\(\s*({NUMBER})\s*[, ]\s*({NUMBER})?\s*\)", element.get("transform", "")
    ):
        x += float(first)
        y += float(second or 0)
    return [x, y]


def _children(element: ET.Element, tag: str) -> list[ET.Element]:
    return [child for child in element if child.tag == SVG + tag]


def _by_id(root: ET.Element, ident: str) -> ET.Element | None:
    return next((element for element in root.iter() if element.get("id") == ident), None)


def _solid_rows(group: ET.Element) -> list[dict[str, Any]]:
    rows = []
    for entry in _children(group, "g"):
        key = entry.get("data-legend-key")
        if key is None:
            continue
        swatch = _children(entry, "path")[0]
        label = _children(entry, "text")[0]
        rows.append({
            "key": key, "type": "solid", "stroke": swatch.get("stroke") or "none",
            "strokeWidth": float(swatch.get("stroke-width") or 0),
            "_font": (label.get("font-family"), float(label.get("font-size"))),
        })
    return rows


def _gradient_rows(group: ET.Element | None) -> tuple[list[dict[str, Any]], str | None]:
    if group is None:
        return [], None
    rows = []
    labels = [label.text or "" for label in _children(group, "text")]
    for entry in _children(group, "g"):
        key = entry.get("data-legend-key")
        if key is None:
            continue
        bar = _children(entry, "path")[0]
        labels.extend(label.text or "" for label in _children(entry, "text")[1:])
        rows.append({
            "key": key, "type": "gradient", "stroke": bar.get("stroke") or "none",
            "strokeWidth": float(bar.get("stroke-width") or 0),
        })
    min_label = next((label for label in labels if label.endswith("%") and label != "100%"), None)
    return rows, min_label


def _box_from_payload(payload: dict[str, float]) -> list[float]:
    return [payload["x"], payload["y"], payload["x"] + payload["width"], payload["y"] + payload["height"]]


def _gallery_cases() -> list[dict[str, Any]]:
    cases = []
    for path in sorted(GALLERY_SESSIONS.iterdir()):
        if not path.name.endswith((".gbdraw-session.json", ".gbdraw-session.json.gz")):
            continue
        opener = gzip.open if path.name.endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8") as handle:
            session = json.load(handle)
        for index, result in enumerate(session.get("results") or []):
            root = ET.fromstring(result["content"].encode())
            legend = _by_id(root, "legend")
            if legend is None:
                continue
            meta = json.loads(root.get("data-gbdraw-composition"))
            linear = _by_id(root, "legend_horizontal") is not None
            if linear:
                horizontal = _by_id(root, "legend_horizontal")
                solids = _solid_rows(_by_id(horizontal, "feature_legend_h"))
                gradients, min_label = _gradient_rows(_by_id(horizontal, "pairwise_legend_h"))
            else:
                feature = _by_id(legend, "feature_legend")
                solids = _solid_rows(feature if feature is not None else legend)
                gradients, min_label = _gradient_rows(_by_id(legend, "conservation_identity_legend"))
            family, size = solids[0].pop("_font")
            for row in solids[1:]:
                row.pop("_font")
            if gradients and min_label:
                for row in gradients:
                    row["minValue"] = float(min_label.rstrip("%"))
            auto = meta["primary"]["automaticTranslation"]
            final = meta["primary"]["finalBounds"]
            reflow = meta["legendReflow"]
            if "wrapWidth" in reflow:
                # The inputs Python laid the Legend out with (`legendReflow`).
                if reflow["fontFile"] != _font_file(family):
                    raise SystemExit(f"{path.name}#{index}: legendReflow.fontFile {reflow['fontFile']!r} "
                                     f"is not the face of {family!r}")
                recorded = reflow["primaryLocalBounds"]
                primary = [recorded["minX"], recorded["minY"], recorded["maxX"], recorded["maxY"]]
                wrap_width, size, dpi = reflow["wrapWidth"], reflow["fontSize"], reflow["dpi"]
            else:
                # A Result written before Python recorded them: the Web's fallbacks.
                primary = [final["x"] - auto[0], final["y"] - auto[1]]
                primary += [primary[0] + final["width"], primary[1] + final["height"]]
                wrap_width, dpi = final["width"], DPI
            title = meta.get("title")
            cases.append({
                "name": f"{path.name.split('.gbdraw')[0]}#{index}",
                "mode": "linear" if linear else "circular",
                "edit": "none",
                "rows": solids + gradients,
                "options": {
                    "side": meta["legendSide"], "wrapWidth": wrap_width, "fontFile": _font_file(family),
                    "fontFamily": family, "fontSize": size, "dpi": dpi,
                    "colorRectSize": reflow["colorRectSize"],
                },
                "composition": {
                    "primary": primary,
                    "title": _box_from_payload(title["localBounds"]) if title else None,
                    "titleSide": meta.get("titleSide", "none") if title else "none",
                    "overlayObstacles": [
                        [o["x"] - auto[0], o["y"] - auto[1], o["x"] - auto[0] + o["width"], o["y"] - auto[1] + o["height"]]
                        for o in meta.get("overlayObstacles", [])
                    ],
                    "spacing": meta["spacing"],
                },
            })
    return cases


def _edited(case: dict[str, Any]) -> list[dict[str, Any]]:
    rows = case["rows"]
    solids = [row for row in rows if row["type"] == "solid"]
    variants: list[tuple[str, list[dict[str, Any]], str]] = []
    side = case["options"]["side"]

    def renamed(old: str, new: str) -> list[dict[str, Any]]:
        return [{**row, "key": new} if row["key"] == old else row for row in rows]

    if solids:
        variants.append(("rename the first row to a long caption", renamed(solids[0]["key"], LONG_CAPTION), side))
        variants.append(("rename the last row", renamed(solids[-1]["key"], "GC percent"), side))
    if len(solids) >= 2:
        variants.append(("delete the second row", [row for row in rows if row is not solids[1]], side))
        variants.append(("delete the first row", [row for row in rows if row is not solids[0]], side))
        variants.append(("reverse the row order", list(reversed(solids)) + [row for row in rows if row["type"] != "solid"], side))
    if solids:
        first = solids[0]
        added = {"key": "Manual row", "type": "solid", "stroke": first["stroke"], "strokeWidth": first["strokeWidth"]}
        variants.append(("add a row", rows + [added], side))
        another = {**added, "key": "Another manually added row with a long name"}
        variants.append(("add two rows", rows + [added, another], side))
    if case["mode"] == "circular":
        variants += [(f"move the Legend to {other}", rows, other) for other in ("top", "bottom", "left", "right") if other != side]
    return [
        {**case, "edit": label, "rows": edited_rows, "options": {**case["options"], "side": edited_side}}
        for label, edited_rows, edited_side in variants
    ]


def _synthetic(template: dict[str, Any]) -> list[dict[str, Any]]:
    options = template["options"]
    composition = template["composition"]
    solid = [{"key": caption, "type": "solid", "stroke": "gray", "strokeWidth": 2.0} for caption in (
        "CDS", "tRNA", "rRNA", "repeat_region", "GC content", "GC skew (+)", "GC skew (-)",
        "hypothetical chloroplast reading frames (ycf)", "Café", "日本語",
    )]
    conservation = [{"key": "Conservation identity", "type": "gradient", "stroke": "none", "strokeWidth": 0.0, "minValue": 70}]
    compact = [
        {"key": "Sequence A identity", "type": "gradient", "stroke": "none", "strokeWidth": 0.0, "minValue": 52.5},
        {"key": "B", "type": "gradient", "stroke": "black", "strokeWidth": 1.0, "minValue": 52.5},
    ]
    cases = []
    for mode in ("circular", "linear"):
        for side in ("left", "top", "upper_right", "lower_left"):
            for label, rows in (
                ("solid rows", solid),
                ("solid rows and one gradient", solid + conservation),
                ("solid rows and two gradients", solid + compact),
                ("one gradient only", conservation),
                ("one very wide row", [{"key": "W" * 160, "type": "solid", "stroke": "none", "strokeWidth": 0.0}] + solid[:2]),
            ):
                cases.append({
                    "name": f"synthetic {mode}", "mode": mode, "edit": f"{label}, {side}", "rows": rows,
                    "options": {**options, "side": side, "wrapWidth": 900.0},
                    "composition": composition,
                })
        cases.append({
            "name": f"synthetic {mode}", "mode": mode, "edit": "empty table", "rows": [],
            "options": {**options, "side": "left"}, "composition": composition,
        })
    return cases


def build_inputs() -> dict[str, Any]:
    gallery = _gallery_cases()
    layout = []
    for case in gallery:
        layout.append(case)
        layout.extend(_edited(case))
    circular_template = next(case for case in gallery if case["mode"] == "circular")
    layout.extend(_synthetic(circular_template))
    measurement = []
    for family in MEASUREMENT_FAMILIES:
        for size in (14.0, 20.0, 16.5):
            for caption in MEASUREMENT_TEXTS:
                measurement.append({
                    "text": caption, "fontFamily": family, "fontFile": _font_file(family),
                    "fontSize": size, "dpi": DPI,
                })
    # A Result laid out at another dpi (`legendReflow.dpi`) measures its rows,
    # an added row too, at that dpi.
    for dpi in (72, 150):
        for caption in ("Added legend entry", "GC skew (+)", "To Ty Ye"):
            family = MEASUREMENT_FAMILIES[0]
            measurement.append({
                "text": caption, "fontFamily": family, "fontFile": _font_file(family),
                "fontSize": 20.0, "dpi": dpi,
            })
    return {"overlayPolicy": OVERLAY_POLICY, "measurement": measurement, "layout": layout}


# ---- expected values (Python's own functions) ----

def _table(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    table: dict[str, dict[str, Any]] = {}
    for row in rows:
        properties: dict[str, Any] = {
            "type": row["type"], "fill": "#000000", "stroke": row["stroke"], "width": row["strokeWidth"],
        }
        if row["type"] == "gradient":
            properties.update({"min_value": row.get("minValue", 0), "min_color": "#000000", "max_color": "#ffffff"})
        table[row["key"]] = properties
    return table


def _box(box: Aabb) -> list[float]:
    return [box.min_x, box.min_y, box.max_x, box.max_y]


def _aabb(values: list[float]) -> Aabb:
    return Aabb(*values)


def _entries(entries) -> list[dict[str, Any]]:
    return [
        {"key": entry.key, "rectX": entry.rect_x, "rectY": entry.rect_y, "textX": entry.text_x, "textY": entry.text_y}
        for entry in entries
    ]


def _linear_expected(case: dict[str, Any]) -> tuple[dict[str, Any], Aabb]:
    options = case["options"]
    layout = build_linear_legend_layout(
        _table(case["rows"]), legend_position=options["side"], canvas_width=options["wrapWidth"],
        font_family=options["fontFamily"], font_size=options["fontSize"], dpi=options["dpi"],
        color_rect_size=options["colorRectSize"],
    )

    def orientation(value) -> dict[str, Any]:
        gradient = value.gradient
        return {
            "feature": {
                "entries": _entries(value.feature.entries), "width": value.feature.width,
                "height": value.feature.height, "numLines": value.feature.num_lines,
            },
            "gradient": None if gradient is None else {
                "compact": gradient.compact,
                "entries": [
                    {"key": entry.key, "titleX": entry.title_x, "titleY": entry.title_y, "barX": entry.bar_x, "barY": entry.bar_y}
                    for entry in gradient.entries
                ],
                "width": gradient.width, "height": gradient.height, "barWidth": gradient.bar_width,
                "minLabelText": gradient.min_label_text, "minLabelX": gradient.min_label_x,
                "maxLabelX": gradient.max_label_x, "scaleLabelY": gradient.scale_label_y,
            },
            "featureX": value.feature_x, "featureY": value.feature_y,
            "gradientX": value.gradient_x, "gradientY": value.gradient_y,
            "width": value.width, "height": value.height,
        }

    bounds = _linear_legend_local_bounds(layout, color_rect_size=options["colorRectSize"])
    return {
        "horizontal": orientation(layout.horizontal), "vertical": orientation(layout.vertical),
        "activeOrientation": layout.active_orientation,
    }, bounds


def _circular_expected(case: dict[str, Any]) -> tuple[dict[str, Any], Aabb]:
    options = case["options"]
    layout = build_circular_legend_layout(
        _table(case["rows"]), legend_position=options["side"], canvas_width=options["wrapWidth"],
        font_family=options["fontFamily"], font_size=options["fontSize"], dpi=options["dpi"],
        color_rect_size=options["colorRectSize"],
    )
    gradient = layout.gradient
    bounds = _circular_legend_local_bounds(layout, color_rect_size=options["colorRectSize"])
    return {
        "horizontal": layout.horizontal, "width": layout.width, "height": layout.height,
        "featureWidth": layout.feature_width, "featureHeight": layout.feature_height,
        "pairwiseLegendWidth": layout.pairwise_legend_width, "lineMargin": layout.line_margin,
        "xMargin": layout.x_margin, "numLines": layout.num_lines, "numColumns": layout.num_columns,
        "numItemsPerLine": layout.num_items_per_line, "entries": _entries(layout.solid_entries),
        "gradient": None if gradient is None else {
            "compact": gradient.compact, "width": gradient.width, "height": gradient.height,
            "barWidth": gradient.bar_width, "barX": gradient.bar_x, "minLabelText": gradient.min_label_text,
            "scaleY": gradient.scale_y,
            "compactEntries": [
                {"key": entry.key, "labelY": entry.label_y, "barY": entry.bar_y} for entry in gradient.compact_entries
            ],
            "singleEntries": [
                {
                    "key": entry.key, "titleX": entry.title_x, "titleY": entry.title_y, "barX": entry.bar_x,
                    "barY": entry.bar_y, "minLabelX": entry.min_label_x, "maxLabelX": entry.max_label_x,
                    "scaleLabelY": entry.scale_label_y,
                }
                for entry in gradient.single_entries
            ],
        },
        "gradientX": layout.gradient_x, "gradientY": layout.gradient_y,
    }, bounds


def _composition_expected(case: dict[str, Any], legend: Aabb) -> dict[str, Any]:
    inputs = case["composition"]
    spacing = inputs["spacing"]
    side = case["options"]["side"]
    plan = plan_composition(CompositionRequest(
        primary=CompositionItem("primary", _aabb(inputs["primary"])),
        legend=CompositionItem("legend", legend),
        title=CompositionItem("title", _aabb(inputs["title"])) if inputs["title"] else None,
        legend_placement=side,
        title_placement=inputs["titleSide"] if inputs["title"] else "none",
        overlay_obstacles=tuple(_aabb(obstacle) for obstacle in inputs["overlayObstacles"]),
        spacing=CompositionSpacing(
            edge_padding_px=spacing["edgePaddingPx"], dock_gap_px=spacing["dockGapPx"],
            title_gap_px=spacing["titleGapPx"], stack_gap_px=spacing["stackGapPx"],
            overlay_clearance_px=spacing["overlayClearancePx"],
        ),
    ))
    return {
        "canvas": _box(plan.canvas_bounds),
        "placements": [
            {"role": placement.role, "translation": list(placement.translation), "finalBounds": _box(placement.final_bounds)}
            for placement in plan.placements
        ],
        "overlayObstacles": [_box(obstacle) for obstacle in plan.overlay_obstacles],
        "overlayConflictIndices": list(plan.overlay_conflict_indices),
        "overlayResolution": plan.overlay_resolution.value if plan.overlay_resolution else None,
    }


def expected_values(inputs: dict[str, Any]) -> dict[str, Any]:
    measurement = []
    for case in inputs["measurement"]:
        width, height = calculate_bbox_dimensions(case["text"], case["fontFamily"], case["fontSize"], case["dpi"])
        measurement.append({**case, "expected": {"width": width, "height": height}})
    layout = []
    for case in inputs["layout"]:
        build = _linear_expected if case["mode"] == "linear" else _circular_expected
        result, bounds = build(case)
        result["localBounds"] = _box(bounds)
        result["composition"] = _composition_expected(case, bounds)
        layout.append({**case, "expected": result})
    return {
        "description": (
            "Python's Legend measurement, layout, local bounds and composition for shared inputs. "
            "Generated by tools/generate_legend_layout_vectors.py; boxes are [minX, minY, maxX, maxY]."
        ),
        "overlayPolicy": inputs["overlayPolicy"],
        "measurement": measurement,
        "layout": layout,
    }


def _stored_inputs() -> dict[str, Any]:
    stored = json.loads(TARGET.read_text(encoding="utf-8"))
    strip = lambda cases: [{key: value for key, value in case.items() if key != "expected"} for case in cases]  # noqa: E731
    return {
        "overlayPolicy": stored["overlayPolicy"],
        "measurement": strip(stored["measurement"]),
        "layout": strip(stored["layout"]),
    }


def render(inputs: dict[str, Any]) -> str:
    # One case per line keeps a regeneration's diff readable case by case.
    values = expected_values(inputs)
    compact = lambda value: json.dumps(value, ensure_ascii=False, separators=(",", ":"))  # noqa: E731
    lines = ["{"]
    lines.append(f' "description": {compact(values["description"])},')
    lines.append(f' "overlayPolicy": {compact(values["overlayPolicy"])},')
    for key in ("measurement", "layout"):
        cases = values[key]
        lines.append(f' "{key}": [')
        lines.extend(f"  {compact(case)}{',' if index < len(cases) - 1 else ''}" for index, case in enumerate(cases))
        lines.append(" ]," if key == "measurement" else " ]")
    lines.append("}")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="Fail if the expected values are stale.")
    parser.add_argument("--from-gallery", action="store_true", help="Rebuild the inputs from the Gallery Sessions.")
    args = parser.parse_args()
    inputs = build_inputs() if args.from_gallery or not TARGET.is_file() else _stored_inputs()
    expected = render(inputs)
    if args.check:
        if not TARGET.is_file() or TARGET.read_text(encoding="utf-8") != expected:
            print("Legend layout vectors are stale. Run: python tools/generate_legend_layout_vectors.py")
            return 1
        return 0
    TARGET.write_text(expected, encoding="utf-8", newline="\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
