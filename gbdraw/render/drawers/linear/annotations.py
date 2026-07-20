"""Linear annotation SVG geometry."""

from __future__ import annotations

import hashlib
import re

from svgwrite import Drawing  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]
from svgwrite.path import Path  # type: ignore[reportMissingImports]
from svgwrite.shapes import Line, Rect  # type: ignore[reportMissingImports]
from svgwrite.text import Text  # type: ignore[reportMissingImports]

from gbdraw.annotations import (
    AnnotationTrackParams,
    ResolvedAnnotationTrack,
    apply_feature_underlay_dom_attributes,
    effective_annotation_style,
    is_auto_feature_underlay,
)
from gbdraw.render.patterns import ensure_hatch_pattern


def _dom_token(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "-", str(value)).strip("-") or "annotation"
    return cleaned[:48]


def annotation_dom_id(
    *,
    record_index: int,
    slot_id: str,
    set_id: str,
    annotation_id: str,
) -> str:
    raw = f"{record_index}|{slot_id}|{set_id}|{annotation_id}"
    digest = hashlib.sha1(raw.encode("utf-8")).hexdigest()[:10]
    return f"gbdraw-annotation-{_dom_token(slot_id)}-{_dom_token(annotation_id)}-{digest}"


def draw_linear_annotation_track(
    drawing: Drawing,
    track: ResolvedAnnotationTrack,
    *,
    record_id: str,
    record_index: int,
    record_length: int,
    bar_length_px: float,
    y_offset_px: float,
    side: str,
    height_px: float,
    font_family: str,
    params: AnnotationTrackParams,
) -> Group:
    """Create a positioned annotation group in record-local coordinates."""

    group = Group(id=f"gbdraw-annotation-track-{_dom_token(track.slot_id)}-{record_index + 1}", debug=False)
    placements = [item for item in track.placements if item.annotation.record_index == record_index]
    lane_count = max((item.lane for item in placements), default=-1) + 1
    usable_height = max(1.0, float(height_px) - 2.0 * params.padding_px)
    lane_height = 14.0 if track.clipped else usable_height / max(1, lane_count)
    scale = float(bar_length_px) / max(1, int(record_length))
    label_direction = (
        -1.0
        if side == "above" or (side == "overlay" and params.layer == "underlay")
        else 1.0
    )

    for placement in placements:
        annotation = placement.annotation
        auto_feature_underlay = is_auto_feature_underlay(annotation)
        style = effective_annotation_style(annotation, params)
        lane_center = float(y_offset_px) + (
            params.padding_px + (placement.lane + 0.5) * lane_height
        )
        item_id = annotation_dom_id(
            record_index=record_index,
            slot_id=track.slot_id,
            set_id=track.set_id,
            annotation_id=annotation.id,
        )
        item_group = Group(id=item_id, debug=False)
        if not auto_feature_underlay:
            item_group.attribs.update(
                {
                    "data-gbdraw-annotation-id": annotation.id,
                    "data-gbdraw-annotation-set-id": track.set_id,
                    "data-gbdraw-annotation-track-id": track.slot_id,
                    "data-gbdraw-record-id": record_id,
                    "data-gbdraw-record-index": str(record_index),
                    "data-gbdraw-annotation-mark": annotation.mark,
                    "data-gbdraw-annotation-label": annotation.label,
                }
            )
        dash = ",".join(f"{value:g}" for value in style.stroke_dasharray) or None
        part_count = len(annotation.segments)
        for part_index, (start, end) in enumerate(annotation.segments, start=1):
            x1, x2 = float(start) * scale, float(end) * scale
            if annotation.mark in {"band", "highlight"}:
                fill = ensure_hatch_pattern(drawing, style.hatch) if style.hatch else (
                    style.fill or ("#94a3b8" if annotation.mark == "highlight" else "none")
                )
                height_factor = 1.0 if annotation.mark == "highlight" else 0.84
                rect = Rect(
                    insert=(x1, lane_center - 0.5 * height_factor * lane_height),
                    size=(max(0.1, x2 - x1), height_factor * lane_height),
                    fill=fill,
                    fill_opacity=style.fill_opacity,
                    stroke="none" if auto_feature_underlay else style.stroke,
                    stroke_width=0 if auto_feature_underlay else style.stroke_width,
                    debug=False,
                )
                if dash:
                    rect.attribs["stroke-dasharray"] = dash
                if auto_feature_underlay:
                    apply_feature_underlay_dom_attributes(
                        rect,
                        annotation,
                        record_id=record_id,
                        record_index=record_index,
                        part_index=part_index,
                        part_count=part_count,
                        include_stable_id=True,
                    )
                item_group.add(rect)
                continue
            line = Line(
                start=(x1, lane_center),
                end=(x2, lane_center),
                fill="none",
                stroke=style.stroke,
                stroke_width=style.stroke_width,
                stroke_linecap="round",
            )
            if dash:
                line.attribs["stroke-dasharray"] = dash
            item_group.add(line)
            if annotation.mark == "bracket" and style.line_cap != "none":
                cap = min(5.0, 0.35 * lane_height)
                if style.line_cap == "arrow":
                    path = Path(
                        d=f"M {x2 - cap:g},{lane_center - cap:g} L {x2:g},{lane_center:g} "
                        f"L {x2 - cap:g},{lane_center + cap:g}",
                        fill="none",
                        stroke=style.stroke,
                        stroke_width=style.stroke_width,
                    )
                    item_group.add(path)
                    item_group.add(Line(start=(x1, lane_center - cap), end=(x1, lane_center + cap), stroke=style.stroke, stroke_width=style.stroke_width))
                else:
                    item_group.add(Line(start=(x1, lane_center - cap), end=(x1, lane_center + cap), stroke=style.stroke, stroke_width=style.stroke_width))
                    item_group.add(Line(start=(x2, lane_center - cap), end=(x2, lane_center + cap), stroke=style.stroke, stroke_width=style.stroke_width))

        if params.show_labels and annotation.label:
            if style.label_position == "start":
                label_bp, anchor = annotation.segments[0][0], "start"
            elif style.label_position == "end":
                label_bp, anchor = annotation.segments[-1][1], "end"
            else:
                label_bp, anchor = annotation.midpoint_bp, "middle"
            label_y = lane_center + label_direction * (0.5 * lane_height + style.label_offset)
            item_group.add(
                Text(
                    annotation.label,
                    insert=(float(label_bp) * scale, label_y),
                    fill=style.label_color,
                    stroke="none",
                    font_size=style.label_font_size or 10.0,
                    font_family=font_family,
                    text_anchor=anchor,
                    dominant_baseline="central",
                )
            )
        group.add(item_group)
    if track.clipped:
        group.attribs["data-gbdraw-annotation-clipped"] = "true"
        clip_id = f"gbdraw-annotation-clip-{_dom_token(track.slot_id)}-{record_index + 1}"
        clip_path = drawing.clipPath(id=clip_id)
        clip_path.add(
            Rect(
                insert=(0, float(y_offset_px)),
                size=(float(bar_length_px), float(height_px)),
            )
        )
        drawing.defs.add(clip_path)
        group.attribs["clip-path"] = f"url(#{clip_id})"
    return group


__all__ = ["annotation_dom_id", "draw_linear_annotation_track"]
