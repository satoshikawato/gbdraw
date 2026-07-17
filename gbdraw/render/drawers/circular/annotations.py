"""Circular annotation SVG geometry."""

from __future__ import annotations

import math

from svgwrite import Drawing  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]
from svgwrite.path import Path  # type: ignore[reportMissingImports]
from svgwrite.text import Text  # type: ignore[reportMissingImports]

from gbdraw.annotations import (
    AnnotationTrackParams,
    ResolvedAnnotationTrack,
    effective_annotation_style,
)
from gbdraw.render.drawers.linear.annotations import annotation_dom_id
from gbdraw.render.patterns import ensure_hatch_pattern


def _point(position: float, total_length: int, radius: float) -> tuple[float, float]:
    angle = 2.0 * math.pi * (float(position) / max(1, total_length)) - 0.5 * math.pi
    return radius * math.cos(angle), radius * math.sin(angle)


def _arc_path(start: float, end: float, total_length: int, radius: float) -> str:
    x1, y1 = _point(start, total_length, radius)
    x2, y2 = _point(end, total_length, radius)
    span = max(0.0, end - start)
    large = 1 if span > 0.5 * total_length else 0
    return f"M {x1:g},{y1:g} A {radius:g},{radius:g} 0 {large} 1 {x2:g},{y2:g}"


def _band_path(start: float, end: float, total_length: int, inner: float, outer: float) -> str:
    ox1, oy1 = _point(start, total_length, outer)
    ox2, oy2 = _point(end, total_length, outer)
    ix2, iy2 = _point(end, total_length, inner)
    ix1, iy1 = _point(start, total_length, inner)
    span = max(0.0, end - start)
    large = 1 if span > 0.5 * total_length else 0
    return (
        f"M {ox1:g},{oy1:g} A {outer:g},{outer:g} 0 {large} 1 {ox2:g},{oy2:g} "
        f"L {ix2:g},{iy2:g} A {inner:g},{inner:g} 0 {large} 0 {ix1:g},{iy1:g} Z"
    )


def _safe_segments(start: int, end: int, total_length: int):
    if end - start >= total_length:
        middle = start + total_length / 2.0
        return ((float(start), middle), (middle, float(end)))
    return ((float(start), float(end)),)


def draw_circular_annotation_track(
    drawing: Drawing,
    track: ResolvedAnnotationTrack,
    *,
    record_id: str,
    record_index: int,
    record_length: int,
    inner_radius_px: float,
    outer_radius_px: float,
    side: str,
    font_family: str,
    params: AnnotationTrackParams,
) -> Group:
    group = Group(id=f"gbdraw-annotation-track-{track.slot_id}-{record_index + 1}", debug=False)
    placements = [item for item in track.placements if item.annotation.record_index == record_index]
    lane_count = max((item.lane for item in placements), default=-1) + 1
    usable = max(1.0, outer_radius_px - inner_radius_px - 2.0 * params.padding_px)
    lane_width = 14.0 if track.clipped else usable / max(1, lane_count)

    for placement in placements:
        annotation = placement.annotation
        style = effective_annotation_style(annotation, params)
        if side == "inside":
            lane_outer = outer_radius_px - params.padding_px - placement.lane * lane_width
            lane_inner = lane_outer - lane_width
        else:
            lane_inner = inner_radius_px + params.padding_px + placement.lane * lane_width
            lane_outer = lane_inner + lane_width
        radius = 0.5 * (lane_inner + lane_outer)
        item_group = Group(
            id=annotation_dom_id(
                record_index=record_index,
                slot_id=track.slot_id,
                set_id=track.set_id,
                annotation_id=annotation.id,
            ),
            debug=False,
        )
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
        for raw_start, raw_end in annotation.segments:
            for start, end in _safe_segments(raw_start, raw_end, record_length):
                if annotation.mark == "band":
                    fill = ensure_hatch_pattern(drawing, style.hatch) if style.hatch else (style.fill or "none")
                    path = Path(
                        d=_band_path(start, end, record_length, lane_inner, lane_outer),
                        fill=fill,
                        fill_opacity=style.fill_opacity,
                        stroke=style.stroke,
                        stroke_width=style.stroke_width,
                    )
                else:
                    path = Path(
                        d=_arc_path(start, end, record_length, radius),
                        fill="none",
                        stroke=style.stroke,
                        stroke_width=style.stroke_width,
                        stroke_linecap="round",
                    )
                if dash:
                    path.attribs["stroke-dasharray"] = dash
                item_group.add(path)
                if annotation.mark == "bracket" and style.line_cap != "none":
                    cap = min(5.0, 0.35 * lane_width)
                    for position in (start, end):
                        x1, y1 = _point(position, record_length, radius - cap)
                        x2, y2 = _point(position, record_length, radius + cap)
                        cap_path = Path(
                            d=f"M {x1:g},{y1:g} L {x2:g},{y2:g}",
                            fill="none",
                            stroke=style.stroke,
                            stroke_width=style.stroke_width,
                        )
                        item_group.add(cap_path)

        if params.show_labels and annotation.label:
            label_radius = lane_outer + style.label_offset if side != "inside" else lane_inner - style.label_offset
            x, y = _point(annotation.midpoint_bp, record_length, label_radius)
            angle = 360.0 * annotation.midpoint_bp / max(1, record_length) - 90.0
            orientation = style.label_orientation
            rotation = 0.0
            if orientation in {"auto", "tangent", "arc"}:
                rotation = angle + 90.0
                if 90.0 < rotation % 360.0 < 270.0:
                    rotation += 180.0
            elif orientation == "radial":
                rotation = angle
            text = Text(
                annotation.label,
                insert=(x, y),
                fill=style.label_color,
                stroke="none",
                font_size=style.label_font_size or 10.0,
                font_family=font_family,
                text_anchor="middle",
                dominant_baseline="central",
            )
            if rotation:
                text.rotate(rotation, center=(x, y))
            item_group.add(text)
        group.add(item_group)
    if track.clipped:
        group.attribs["data-gbdraw-annotation-clipped"] = "true"
        clip_id = f"gbdraw-annotation-clip-{track.slot_id}-{record_index + 1}"
        ring = Path(
            d=(
                f"M {outer_radius_px:g},0 A {outer_radius_px:g},{outer_radius_px:g} 0 1 0 "
                f"{-outer_radius_px:g},0 A {outer_radius_px:g},{outer_radius_px:g} 0 1 0 "
                f"{outer_radius_px:g},0 Z M {inner_radius_px:g},0 "
                f"A {inner_radius_px:g},{inner_radius_px:g} 0 1 1 {-inner_radius_px:g},0 "
                f"A {inner_radius_px:g},{inner_radius_px:g} 0 1 1 {inner_radius_px:g},0 Z"
            )
        )
        ring.attribs["clip-rule"] = "evenodd"
        clip_path = drawing.clipPath(id=clip_id)
        clip_path.add(ring)
        drawing.defs.add(clip_path)
        group.attribs["clip-path"] = f"url(#{clip_id})"
    return group


__all__ = ["draw_circular_annotation_track"]
