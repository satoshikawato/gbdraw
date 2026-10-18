"""Circular annotation SVG geometry."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

from svgwrite.base import BaseElement

from gbdraw.svg.elements import ClipPath, Group, Path, Text

from gbdraw.annotations import (
    AnnotationTrackParams,
    ResolvedAnnotationTrack,
    apply_feature_underlay_dom_attributes,
    effective_annotation_style,
    is_auto_feature_underlay,
)
from gbdraw.annotations.layout import PlacedRegionAnnotation
from gbdraw.layout.circular import RadialBand
from gbdraw.layout.spatial import Aabb
from gbdraw.render.drawers.annotation_marks import (
    annotation_dom_id,
    cap_size,
    label_font_size,
    label_ink_offsets,
    mark_reach,
)
from gbdraw.render.group_content import GroupContent
from gbdraw.render.patterns import hatch_pattern, hatch_pattern_id, hatch_pattern_paint
from gbdraw.svg.ids import stable_svg_id


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


@dataclass(frozen=True)
class CircularAnnotationLabel:
    """One Circular annotation label as it is drawn (``text-anchor=middle``, ``dominant-baseline=central``).

    ``x`` and ``y`` anchor the text around the circle centre, ``rotation`` turns
    it about that point (degrees, clockwise), and ``corners`` are the corners of
    the box of its painted letters after the turn, in the same frame.
    """

    annotation_id: str
    text: str
    x: float
    y: float
    rotation: float
    font_size: float
    color: str
    corners: tuple[tuple[float, float], ...]

    @property
    def ink_bounds(self) -> Aabb:
        xs = [x for x, _y in self.corners]
        ys = [y for _x, y in self.corners]
        return Aabb(min(xs), min(ys), max(xs), max(ys))


def _record_lanes(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    inner_radius_px: float,
    outer_radius_px: float,
    params: AnnotationTrackParams,
) -> tuple[list[PlacedRegionAnnotation], float]:
    """The record's placements and the width of one lane."""

    placements = [item for item in track.placements if item.annotation.record_index == record_index]
    lane_count = max((item.lane for item in placements), default=-1) + 1
    usable = max(1.0, outer_radius_px - inner_radius_px - 2.0 * params.padding_px)
    return placements, 14.0 if track.clipped else usable / max(1, lane_count)


def _lane_radii(
    lane: int,
    lane_width: float,
    *,
    side: str,
    inner_radius_px: float,
    outer_radius_px: float,
    params: AnnotationTrackParams,
) -> tuple[float, float]:
    """The inner and outer radius of a lane: ``inside`` slots stack their lanes inward."""

    if side == "inside":
        lane_outer = outer_radius_px - params.padding_px - lane * lane_width
        return lane_outer - lane_width, lane_outer
    lane_inner = inner_radius_px + params.padding_px + lane * lane_width
    return lane_inner, lane_inner + lane_width


def _label_rotation(orientation: str, angle: float) -> float:
    """The turn of a label at ``angle`` (degrees, clockwise from 3 o'clock): along the ring kept upright, along the
    ray, or none."""

    if orientation in {"auto", "tangent", "arc"}:
        rotation = angle + 90.0
        return rotation + 180.0 if 90.0 < rotation % 360.0 < 270.0 else rotation
    return angle if orientation == "radial" else 0.0


def _ink_extent(
    radius: float, direction: tuple[float, float], offsets: tuple[float, float, float, float]
) -> tuple[float, float]:
    """How near the centre and how far from it a label's ink box comes with its anchor at ``radius``.

    ``direction`` is the unit ray in the text's own frame; the box excludes the centre.
    """

    left, right, top, bottom = offsets
    wx, wy = direction
    # The centre lies at -radius * direction from the anchor.
    ox, oy = -radius * wx, -radius * wy
    near = math.hypot(max(left - ox, 0.0, ox - right), max(top - oy, 0.0, oy - bottom))
    far = max(math.hypot(x - ox, y - oy) for x in (left, right) for y in (top, bottom))
    return near, far


def _label_radius(
    target: float,
    outward: bool,
    direction: tuple[float, float],
    offsets: tuple[float, float, float, float],
) -> float:
    """The anchor radius at which a label's ink comes nearest the centre at ``target`` (``outward``) or reaches
    farthest from it at ``target``."""

    left, right, top, bottom = offsets
    projections = [x * direction[0] + y * direction[1] for x in (left, right) for y in (top, bottom)]
    # The projection on the ray bounds the reach; each step closes the rest, as the reach moves almost one for one.
    radius = target - (min(projections) if outward else max(projections))
    for _ in range(32):
        near, far = _ink_extent(radius, direction, offsets)
        miss = target - (near if outward else far)
        if abs(miss) < 1e-9 or radius + miss <= 0.0:
            break
        radius += miss
    return max(radius, 0.0)


def circular_annotation_labels(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    record_length: int,
    inner_radius_px: float,
    outer_radius_px: float,
    side: str,
    font_family: str,
    params: AnnotationTrackParams,
    rings: Sequence[RadialBand] = (),
) -> tuple[CircularAnnotationLabel, ...]:
    """The labels of one record's annotation track, around the centre of the circle.

    Each label's painted letters keep ``label_offset`` from the painted edge of
    its mark (stroke included), outward or, on an ``inside`` slot, inward, at
    the label's own turn (OV-849, the Linear rule of OV-844). A painted ring of
    ``rings`` (the axis circle) that the letters would meet, or come within
    ``label_offset`` of, is cleared the same way beyond its paint.
    """

    if not params.show_labels:
        return ()
    placements, lane_width = _record_lanes(
        track, record_index=record_index, inner_radius_px=inner_radius_px, outer_radius_px=outer_radius_px, params=params
    )
    outward = side != "inside"
    labels: list[CircularAnnotationLabel] = []
    for placement in placements:
        annotation = placement.annotation
        if not annotation.label:
            continue
        style = effective_annotation_style(annotation, params)
        lane_inner, lane_outer = _lane_radii(
            placement.lane,
            lane_width,
            side=side,
            inner_radius_px=inner_radius_px,
            outer_radius_px=outer_radius_px,
            params=params,
        )
        # A Circular band fills its whole lane; an arrow cap is drawn as a straight cap.
        reach = mark_reach(annotation, style, lane_width, fill_share=1.0, arrow_heads=False).from_centre
        centre = 0.5 * (lane_inner + lane_outer)
        edge = centre + reach if outward else centre - reach
        gap = float(style.label_offset)
        size = label_font_size(style)
        offsets = label_ink_offsets(annotation.label, font_family, size, "middle")
        angle = 360.0 * annotation.geometry_midpoint_bp / max(1, record_length) - 90.0
        rotation = _label_rotation(style.label_orientation, angle)
        theta, turn = math.radians(angle), math.radians(rotation)
        # The ray's direction in the text's frame: the screen direction turned back by the label's rotation.
        direction = (math.cos(theta - turn), math.sin(theta - turn))
        radius = _label_radius(edge + gap if outward else edge - gap, outward, direction, offsets)
        for ring in sorted(rings, key=lambda band: band.inner_px if outward else -band.outer_px):
            near, far = _ink_extent(radius, direction, offsets)
            if outward and edge < ring.outer_px and ring.inner_px < far + gap:
                edge = ring.outer_px
            elif not outward and ring.inner_px < edge and near - gap < ring.outer_px:
                edge = ring.inner_px
            else:
                continue
            radius = _label_radius(edge + gap if outward else edge - gap, outward, direction, offsets)
        x, y = radius * math.cos(theta), radius * math.sin(theta)
        left, right, top, bottom = offsets
        cos, sin = math.cos(turn), math.sin(turn)
        labels.append(
            CircularAnnotationLabel(
                annotation_id=annotation.id,
                text=annotation.label,
                x=x,
                y=y,
                rotation=rotation,
                font_size=size,
                color=style.label_color,
                corners=tuple(
                    (x + cx * cos - cy * sin, y + cx * sin + cy * cos)
                    for cx, cy in ((left, top), (right, top), (right, bottom), (left, bottom))
                ),
            )
        )
    return tuple(labels)


def circular_annotation_track_content(
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
    labels: Sequence[CircularAnnotationLabel] = (),
    dom_group_id: str | None = None,
    semantic_slot_id: str | None = None,
) -> GroupContent:
    """Build an annotation track around the circle centre (render stage ``annotation_content``).

    ``labels`` are the track's labels (``circular_annotation_labels``), drawn
    with their annotations. The hatch patterns the marks use, in order of first
    use, and the clip path of a clipped track are the content's definitions.
    """

    group_id = dom_group_id or f"gbdraw-annotation-track-{track.slot_id}-{record_index + 1}"
    group = Group(id=group_id)
    definitions: dict[str, BaseElement] = {}
    if semantic_slot_id:
        group.attribs["data-gbdraw-slot-id"] = str(semantic_slot_id)
        group.attribs["data-gbdraw-slot-renderer"] = "annotations"
    placements, lane_width = _record_lanes(
        track, record_index=record_index, inner_radius_px=inner_radius_px, outer_radius_px=outer_radius_px, params=params
    )
    labels_by_id = {label.annotation_id: label for label in labels}

    for placement in placements:
        annotation = placement.annotation
        auto_feature_underlay = is_auto_feature_underlay(annotation)
        style = effective_annotation_style(annotation, params)
        lane_inner, lane_outer = _lane_radii(
            placement.lane,
            lane_width,
            side=side,
            inner_radius_px=inner_radius_px,
            outer_radius_px=outer_radius_px,
            params=params,
        )
        radius = 0.5 * (lane_inner + lane_outer)
        item_group = Group(
            id=annotation_dom_id(
                record_index=record_index,
                slot_id=track.slot_id,
                set_id=track.set_id,
                annotation_id=annotation.id,
            ),
        )
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
        part_count = sum(
            len(_safe_segments(start, end, record_length))
            for start, end in annotation.geometry_segments
        )
        part_index = 0
        for segment_index, (raw_start, raw_end) in enumerate(annotation.geometry_segments):
            fragment = annotation.display_parts[segment_index] if annotation.display_parts is not None else None
            for start, end in _safe_segments(raw_start, raw_end, record_length):
                part_index += 1
                if annotation.mark in {"band", "highlight"}:
                    if style.hatch:
                        definitions.setdefault(hatch_pattern_id(style.hatch), hatch_pattern(style.hatch))
                    fill = hatch_pattern_paint(style.hatch) if style.hatch else (
                        style.fill or ("#94a3b8" if annotation.mark == "highlight" else "none")
                    )
                    path = Path(
                        d=_band_path(start, end, record_length, lane_inner, lane_outer),
                        fill=fill,
                        fill_opacity=style.fill_opacity,
                        stroke="none" if auto_feature_underlay else style.stroke,
                        stroke_width=0 if auto_feature_underlay else style.stroke_width,
                    )
                else:
                    path = Path(
                        d=_arc_path(start, end, record_length, radius),
                        fill="none",
                        stroke=style.stroke,
                        stroke_width=style.stroke_width,
                        stroke_linecap=("butt" if fragment is not None and (fragment.artificial_start or fragment.artificial_end) else "round"),
                    )
                if dash:
                    path.attribs["stroke-dasharray"] = dash
                if auto_feature_underlay:
                    apply_feature_underlay_dom_attributes(
                        path,
                        annotation,
                        record_id=record_id,
                        record_index=record_index,
                        part_index=part_index,
                        part_count=part_count,
                        include_stable_id=False,
                    )
                item_group.add(path)
                if annotation.mark in {"band", "highlight"} and not auto_feature_underlay and fragment is not None and (fragment.artificial_start or fragment.artificial_end):
                    path.attribs["stroke"] = "none"
                    outline = _arc_path(start, end, record_length, lane_inner) + " " + _arc_path(start, end, record_length, lane_outer)
                    for position, capped in ((start, not fragment.artificial_start), (end, not fragment.artificial_end)):
                        if capped:
                            x1, y1 = _point(position, record_length, lane_inner)
                            x2, y2 = _point(position, record_length, lane_outer)
                            outline += f" M {x1:g},{y1:g} L {x2:g},{y2:g}"
                    border = Path(d=outline, fill="none", stroke=style.stroke, stroke_width=style.stroke_width, stroke_linecap="butt")
                    if dash:
                        border.attribs["stroke-dasharray"] = dash
                    item_group.add(border)
                if annotation.mark == "bracket" and style.line_cap != "none":
                    cap = cap_size(lane_width)
                    for position in (start, end):
                        if fragment is not None and (
                            (position == start and (start != raw_start or fragment.artificial_start))
                            or (position == end and (end != raw_end or fragment.artificial_end))
                        ):
                            continue
                        x1, y1 = _point(position, record_length, radius - cap)
                        x2, y2 = _point(position, record_length, radius + cap)
                        cap_path = Path(
                            d=f"M {x1:g},{y1:g} L {x2:g},{y2:g}",
                            fill="none",
                            stroke=style.stroke,
                            stroke_width=style.stroke_width,
                        )
                        item_group.add(cap_path)

        label = labels_by_id.get(annotation.id)
        if label is not None:
            text = Text(
                label.text,
                insert=(label.x, label.y),
                fill=label.color,
                stroke="none",
                font_size=label.font_size,
                font_family=font_family,
                text_anchor="middle",
                dominant_baseline="central",
            )
            if label.rotation:
                text.rotate(label.rotation, center=(label.x, label.y))
            item_group.add(text)
        group.add(item_group)
    clip_paths: list[BaseElement] = []
    if track.clipped:
        group.attribs["data-gbdraw-annotation-clipped"] = "true"
        clip_id = (
            stable_svg_id("annotation_clip", "circular", group_id)
            if dom_group_id
            else f"gbdraw-annotation-clip-{track.slot_id}-{record_index + 1}"
        )
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
        clip_path = ClipPath(id=clip_id)
        clip_path.add(ring)
        clip_paths.append(clip_path)
        group.attribs["clip-path"] = f"url(#{clip_id})"
    return GroupContent.of(group, (*definitions.values(), *clip_paths))


__all__ = ["CircularAnnotationLabel", "circular_annotation_labels", "circular_annotation_track_content"]
