"""Annotation entries for the shared diagram legend."""

from __future__ import annotations

from collections.abc import Sequence

from gbdraw.annotations.models import (
    ResolvedAnnotationBundle,
    annotation_track_params_from_mapping,
    effective_annotation_style,
)
from gbdraw.legend.table import LegendTable
from gbdraw.render.patterns import hatch_pattern_paint


def sync_annotation_legend_entries(
    out: LegendTable,
    bundle: ResolvedAnnotationBundle,
    slots: Sequence[object] | None,
) -> LegendTable:
    """Append explicit annotation legend labels, one source per caption and style."""

    if not slots or not bundle.annotations:
        return out
    seen: set[tuple[object, ...]] = set()
    for slot in slots:
        if not getattr(slot, "enabled", True) or str(getattr(slot, "renderer", "")) != "annotations":
            continue
        params = annotation_track_params_from_mapping(getattr(slot, "params", {}) or {})
        for annotation in bundle.annotations:
            if (
                annotation.set_id != params.set_id
                or (params.marks is not None and annotation.mark not in params.marks)
                or not annotation.legend_label
            ):
                continue
            style = effective_annotation_style(annotation, params)
            hatch = style.hatch
            signature = (
                annotation.legend_label,
                style.stroke,
                style.stroke_width,
                style.fill,
                style.fill_opacity,
                hatch.angle if hatch else None,
                hatch.spacing if hatch else None,
                hatch.color if hatch else None,
                hatch.width if hatch else None,
                hatch.cross if hatch else None,
            )
            if signature in seen:
                continue
            seen.add(signature)
            out.add_source(str(annotation.legend_label), {
                "type": "solid",
                "fill": hatch_pattern_paint(hatch) if hatch else (
                    style.fill
                    or ("#94a3b8" if annotation.mark == "highlight" else "none")
                ),
                "stroke": style.stroke,
                "width": style.stroke_width,
                "annotation_hatch": hatch,
            })
    return out


__all__ = ["sync_annotation_legend_entries"]
