#!/usr/bin/env python
"""Audit visible SVG composition geometry in a local Chromium browser."""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from pathlib import Path
from typing import Any, Iterable, Sequence


_BROWSER_MEASUREMENT = r"""
() => {
  const svg = document.documentElement;
  if (!svg || svg.localName !== 'svg') {
    throw new Error('document root is not an SVG element');
  }

  const finite = (value) => Number.isFinite(Number(value));
  const rootScreenMatrix = svg.getScreenCTM();
  if (!rootScreenMatrix) {
    throw new Error('SVG root has no screen transformation matrix');
  }
  const screenToRoot = rootScreenMatrix.inverse();

  const graphicsSelector = 'path,rect,circle,ellipse,line,polyline,polygon,text,use,image';
  const leafBBoxInRoot = (element) => {
    const style = getComputedStyle(element);
    if (style.display === 'none' || style.visibility === 'hidden' || style.visibility === 'collapse') {
      return null;
    }
    const bbox = element.getBBox();
    const elementScreenMatrix = element.getScreenCTM();
    if (!elementScreenMatrix || ![bbox.x, bbox.y, bbox.width, bbox.height].every(finite)) {
      return null;
    }
    const elementToRoot = screenToRoot.multiply(elementScreenMatrix);
    const corners = [
      new DOMPoint(bbox.x, bbox.y),
      new DOMPoint(bbox.x + bbox.width, bbox.y),
      new DOMPoint(bbox.x, bbox.y + bbox.height),
      new DOMPoint(bbox.x + bbox.width, bbox.y + bbox.height),
    ].map((point) => point.matrixTransform(elementToRoot));
    const xs = corners.map((point) => point.x);
    const ys = corners.map((point) => point.y);
    const strokeVisible = style.stroke !== 'none'
      && style.stroke !== 'transparent'
      && Number.parseFloat(style.strokeOpacity || '1') > 0;
    const halfStroke = strokeVisible
      ? 0.5 * Math.max(0, Number.parseFloat(style.strokeWidth) || 0)
      : 0;
    const expandX = halfStroke * Math.hypot(elementToRoot.a, elementToRoot.c);
    const expandY = halfStroke * Math.hypot(elementToRoot.b, elementToRoot.d);
    return {
      minX: Math.min(...xs) - expandX,
      minY: Math.min(...ys) - expandY,
      maxX: Math.max(...xs) + expandX,
      maxY: Math.max(...ys) + expandY,
    };
  };

  const bboxInRoot = (element) => {
    const leaves = [];
    if (element.matches(graphicsSelector)) leaves.push(element);
    leaves.push(...element.querySelectorAll(graphicsSelector));
    const bounds = leaves.map(leafBBoxInRoot).filter(Boolean);
    if (!bounds.length) return null;
    return {
      minX: Math.min(...bounds.map((item) => item.minX)),
      minY: Math.min(...bounds.map((item) => item.minY)),
      maxX: Math.max(...bounds.map((item) => item.maxX)),
      maxY: Math.max(...bounds.map((item) => item.maxY)),
    };
  };

  const union = (bounds) => {
    const present = bounds.filter(Boolean);
    if (!present.length) return null;
    return {
      minX: Math.min(...present.map((item) => item.minX)),
      minY: Math.min(...present.map((item) => item.minY)),
      maxX: Math.max(...present.map((item) => item.maxX)),
      maxY: Math.max(...present.map((item) => item.maxY)),
    };
  };

  const decorationSelector = [
    '#legend',
    '[data-gbdraw-role="comparison-legend"]',
    '#plot_title',
    '[data-gbdraw-role="plot-title"]',
  ].join(',');
  const nonPaintTags = new Set(['defs', 'style', 'title', 'desc', 'metadata', 'script']);
  const topLevelPaint = Array.from(svg.children).filter((element) => {
    if (nonPaintTags.has(element.localName)) return false;
    if (element.matches(decorationSelector)) return false;
    return true;
  });

  const primaryBounds = union(topLevelPaint.map(bboxInRoot));
  const legendElements = Array.from(
    svg.querySelectorAll('#legend, [data-gbdraw-role="comparison-legend"]')
  ).filter((element, index, all) =>
    !all.some((candidate, candidateIndex) =>
      candidateIndex !== index && candidate.contains(element)
    )
  );
  const titleElements = Array.from(
    svg.querySelectorAll('#plot_title, [data-gbdraw-role="plot-title"]')
  ).filter((element, index, all) =>
    !all.some((candidate, candidateIndex) =>
      candidateIndex !== index && candidate.contains(element)
    )
  );
  const legendBounds = union(legendElements.map(bboxInRoot));
  const titleBounds = union(titleElements.map(bboxInRoot));
  const paintedBounds = union([primaryBounds, legendBounds, titleBounds]);

  const viewBox = svg.viewBox.baseVal;
  const view = {
    x: Number(viewBox.x),
    y: Number(viewBox.y),
    width: Number(viewBox.width),
    height: Number(viewBox.height),
  };
  const rootWidth = Number(svg.width.baseVal.value);
  const rootHeight = Number(svg.height.baseVal.value);

  let dock = 'none';
  let dockGap = null;
  if (primaryBounds && legendBounds) {
    if (legendBounds.minX >= primaryBounds.maxX) {
      dock = 'right';
      dockGap = legendBounds.minX - primaryBounds.maxX;
    } else if (legendBounds.maxX <= primaryBounds.minX) {
      dock = 'left';
      dockGap = primaryBounds.minX - legendBounds.maxX;
    } else if (legendBounds.minY >= primaryBounds.maxY) {
      dock = 'bottom';
      dockGap = legendBounds.minY - primaryBounds.maxY;
    } else if (legendBounds.maxY <= primaryBounds.minY) {
      dock = 'top';
      dockGap = primaryBounds.minY - legendBounds.maxY;
    } else {
      dock = 'overlay';
      const overlapX = Math.min(primaryBounds.maxX, legendBounds.maxX)
        - Math.max(primaryBounds.minX, legendBounds.minX);
      const overlapY = Math.min(primaryBounds.maxY, legendBounds.maxY)
        - Math.max(primaryBounds.minY, legendBounds.minY);
      dockGap = -Math.min(overlapX, overlapY);
    }
  }

  const margins = paintedBounds ? {
    left: paintedBounds.minX - view.x,
    top: paintedBounds.minY - view.y,
    right: view.x + view.width - paintedBounds.maxX,
    bottom: view.y + view.height - paintedBounds.maxY,
  } : null;
  const clipped = margins
    ? Object.values(margins).some((margin) => margin < -0.01)
    : false;
  const titleCenterOffsetX = primaryBounds && titleBounds
    ? 0.5 * (titleBounds.minX + titleBounds.maxX)
      - 0.5 * (primaryBounds.minX + primaryBounds.maxX)
    : null;

  const axis = svg.querySelector('#Axis, [id^="Axis_"]');
  let circularAxisCenter = null;
  if (axis) {
    const axisScreenMatrix = axis.getScreenCTM();
    if (axisScreenMatrix) {
      const point = new DOMPoint(0, 0).matrixTransform(
        screenToRoot.multiply(axisScreenMatrix)
      );
      circularAxisCenter = { x: point.x, y: point.y };
    }
  }

  return {
    mode: axis ? 'circular' : 'linear',
    viewBox: view,
    rootSize: { width: rootWidth, height: rootHeight },
    primaryBounds,
    legendBounds,
    titleBounds,
    paintedBounds,
    dock,
    dockGap,
    margins,
    clipped,
    titleCenterOffsetX,
    circularAxisCenter,
  };
}
"""


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Measure SVG composition bounds with local Chromium.",
    )
    parser.add_argument(
        "path",
        type=Path,
        help="SVG file or directory to scan recursively.",
    )
    parser.add_argument(
        "--json-out",
        type=Path,
        help="Write the complete audit result as JSON.",
    )
    return parser.parse_args(argv)


def _svg_files(path: Path) -> list[Path]:
    resolved = path.resolve()
    if resolved.is_file():
        if resolved.suffix.lower() != ".svg":
            raise ValueError(f"not an SVG file: {path}")
        return [resolved]
    if not resolved.is_dir():
        raise ValueError(f"path does not exist: {path}")
    return sorted(candidate for candidate in resolved.rglob("*.svg") if candidate.is_file())


def _round_floats(value: Any) -> Any:
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("browser returned a non-finite measurement")
        return round(value, 6)
    if isinstance(value, list):
        return [_round_floats(item) for item in value]
    if isinstance(value, dict):
        return {key: _round_floats(item) for key, item in value.items()}
    return value


def _measure_files(files: Iterable[Path], *, display_root: Path) -> list[dict[str, Any]]:
    try:
        from playwright.sync_api import Error as PlaywrightError
        from playwright.sync_api import sync_playwright
    except ImportError as exc:  # pragma: no cover - depends on the local toolchain
        raise RuntimeError(
            "Python Playwright is required; install it and its Chromium browser first."
        ) from exc

    measurements: list[dict[str, Any]] = []
    try:
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(headless=True)
            try:
                page = browser.new_page(viewport={"width": 1600, "height": 1200})
                for svg_path in files:
                    page.goto(svg_path.as_uri(), wait_until="load")
                    measured = page.evaluate(_BROWSER_MEASUREMENT)
                    measured["path"] = str(svg_path.relative_to(display_root))
                    measurements.append(_round_floats(measured))
            finally:
                browser.close()
    except PlaywrightError as exc:  # pragma: no cover - host sandbox dependent
        raise RuntimeError(f"Chromium audit failed: {exc}") from exc
    return measurements


def _stats(values: Sequence[float]) -> dict[str, float | int | None]:
    if not values:
        return {"count": 0, "median": None, "minimum": None, "maximum": None}
    return {
        "count": len(values),
        "median": round(float(statistics.median(values)), 6),
        "minimum": round(float(min(values)), 6),
        "maximum": round(float(max(values)), 6),
    }


def _summarize(measurements: Sequence[dict[str, Any]]) -> dict[str, Any]:
    modes: dict[str, Any] = {}
    for mode in ("circular", "linear"):
        mode_rows = [row for row in measurements if row["mode"] == mode]
        side_gaps = [
            float(row["dockGap"])
            for row in mode_rows
            if row["dock"] in {"left", "right"} and row["dockGap"] is not None
        ]
        modes[mode] = {
            "files": len(mode_rows),
            "sideLegendGap": _stats(side_gaps),
            "sideLegendGapOver150": sum(gap > 150.0 for gap in side_gaps),
            "clipped": sum(bool(row["clipped"]) for row in mode_rows),
        }
    return {"files": len(measurements), "modes": modes}


def _format_number(value: object) -> str:
    if value is None:
        return "—"
    return f"{float(value):.1f}"


def _print_table(measurements: Sequence[dict[str, Any]]) -> None:
    header = f"{'MODE':8} {'DOCK':8} {'GAP':>8} {'CLIP':>5} {'VIEWBOX':>18}  PATH"
    print(header)
    print("-" * len(header))
    for row in measurements:
        view = row["viewBox"]
        view_size = f"{view['width']:.1f}x{view['height']:.1f}"
        print(
            f"{row['mode']:8} {row['dock']:8} {_format_number(row['dockGap']):>8} "
            f"{('yes' if row['clipped'] else 'no'):>5} {view_size:>18}  {row['path']}"
        )


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    try:
        files = _svg_files(args.path)
        if not files:
            raise ValueError(f"no SVG files found under: {args.path}")
        display_root = args.path.resolve().parent if args.path.is_file() else args.path.resolve()
        measurements = _measure_files(files, display_root=display_root)
    except (RuntimeError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    payload = {
        "schema": 1,
        "root": str(args.path.resolve()),
        "summary": _summarize(measurements),
        "files": measurements,
    }
    _print_table(measurements)
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
