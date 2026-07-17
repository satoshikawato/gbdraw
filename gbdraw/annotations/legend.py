"""Annotation entries for the shared diagram legend."""

from __future__ import annotations

from collections.abc import Sequence

from gbdraw.annotations.models import (
    ResolvedAnnotationBundle,
    annotation_track_params_from_mapping,
    effective_annotation_style,
)
from gbdraw.legend.table import _unique_legend_key
from gbdraw.render.patterns import hatch_pattern_paint


def sync_annotation_legend_entries(
    legend_table: dict,
    bundle: ResolvedAnnotationBundle,
    slots: Sequence[object] | None,
) -> dict:
    """Append explicit annotation legend labels, deduplicated by caption and style."""

    if not slots or not bundle.annotations:
        return legend_table
    out = dict(legend_table)
    seen: set[tuple[object, ...]] = set()
    for slot in slots:
        if not getattr(slot, "enabled", True) or str(getattr(slot, "renderer", "")) != "annotations":
            continue
        params = annotation_track_params_from_mapping(getattr(slot, "params", {}) or {})
        for annotation in bundle.annotations:
            if annotation.set_id != params.set_id or not annotation.legend_label:
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
            caption = str(annotation.legend_label)
            out[_unique_legend_key(out, caption)] = {
                "type": "solid",
                "fill": hatch_pattern_paint(hatch) if hatch else (style.fill or "none"),
                "stroke": style.stroke,
                "width": style.stroke_width,
                "annotation_hatch": hatch,
            }
    return out


__all__ = ["sync_annotation_legend_entries"]
