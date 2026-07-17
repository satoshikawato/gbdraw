"""Reusable SVG pattern definitions."""

from __future__ import annotations

import hashlib

from svgwrite import Drawing  # type: ignore[reportMissingImports]
from svgwrite.shapes import Line  # type: ignore[reportMissingImports]

from gbdraw.annotations.models import HatchStyle


def hatch_pattern_id(hatch: HatchStyle) -> str:
    """Return the stable SVG ID shared by marks and legend previews."""

    key = (
        f"{hatch.angle:g}|{hatch.spacing:g}|{hatch.color}|"
        f"{hatch.width:g}|{int(hatch.cross)}"
    )
    return f"gbdraw-hatch-{hashlib.sha1(key.encode('utf-8')).hexdigest()[:12]}"


def hatch_pattern_paint(hatch: HatchStyle) -> str:
    """Return a paint reference for a hatch registered on the diagram."""

    return f"url(#{hatch_pattern_id(hatch)})"


def ensure_hatch_pattern(drawing: Drawing, hatch: HatchStyle) -> str:
    """Return a deduplicated ``url(#...)`` paint for one hatch style."""

    pattern_id = hatch_pattern_id(hatch)
    existing = {
        str(getattr(item, "attribs", {}).get("id", ""))
        for item in getattr(drawing.defs, "elements", ())
    }
    if pattern_id not in existing:
        pattern = drawing.pattern(
            id=pattern_id,
            size=(hatch.spacing, hatch.spacing),
            patternUnits="userSpaceOnUse",
            patternTransform=f"rotate({hatch.angle:g})",
        )
        pattern.add(
            Line(
                start=(0, 0),
                end=(0, hatch.spacing),
                stroke=hatch.color,
                stroke_width=hatch.width,
            )
        )
        if hatch.cross:
            pattern.add(
                Line(
                    start=(0, 0),
                    end=(hatch.spacing, 0),
                    stroke=hatch.color,
                    stroke_width=hatch.width,
                )
            )
        drawing.defs.add(pattern)
    return hatch_pattern_paint(hatch)


__all__ = ["ensure_hatch_pattern", "hatch_pattern_id", "hatch_pattern_paint"]
