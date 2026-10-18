"""Linear annotation SVG geometry."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace

from svgwrite.base import BaseElement

from gbdraw.svg.elements import ClipPath, Group, Line, Path, Rect, Text

from gbdraw.annotations import (
    AnnotationTrackParams,
    RegionAnnotationStyle,
    ResolvedAnnotationTrack,
    ResolvedRegionAnnotation,
    annotation_label_font_size,
    apply_feature_underlay_dom_attributes,
    effective_annotation_style,
    is_auto_feature_underlay,
    linear_annotation_labels_above,
)
from gbdraw.annotations.layout import PlacedRegionAnnotation
from gbdraw.layout.spatial import Aabb
from gbdraw.layout.text_geometry import TextAnchor
from gbdraw.render.drawers.annotation_marks import (
    annotation_dom_id,
    cap_size,
    capped_ends,
    dom_token,
    label_ink_offsets,
    mark_reach,
    past_span_edge,
)
from gbdraw.render.group_content import GroupContent
from gbdraw.render.patterns import hatch_pattern, hatch_pattern_id, hatch_pattern_paint
from gbdraw.svg.ids import stable_svg_id


@dataclass(frozen=True)
class LinearAnnotationLabel:
    """One Linear annotation label as it is drawn (``dominant-baseline=central``).

    ``x`` and ``y`` anchor the text in record-local x and slot-local y;
    ``ink`` is the box of its painted letters in the same frame, and
    ``segments`` the x ranges of its marks: the label's centre may move along
    them, never into a gap between them.
    """

    annotation_id: str
    text: str
    x: float
    y: float
    anchor: TextAnchor
    font_size: float
    color: str
    ink: Aabb
    segments: tuple[tuple[float, float], ...]

    def moved_to(self, x: float) -> LinearAnnotationLabel:
        return replace(self, x=x, ink=self.ink.translated(x - self.x, 0.0))


def _band_height_factor(mark: str) -> float:
    """The share of its lane a ``highlight`` or ``band`` rectangle fills."""

    return 1.0 if mark == "highlight" else 0.84


def _label_anchor(style: RegionAnnotationStyle) -> TextAnchor:
    return "start" if style.label_position == "start" else "end" if style.label_position == "end" else "middle"


def _record_lanes(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    height_px: float,
    params: AnnotationTrackParams,
) -> tuple[list[PlacedRegionAnnotation], float]:
    """The record's placements and the height of one lane."""

    placements = [item for item in track.placements if item.annotation.record_index == record_index]
    lane_count = max((item.lane for item in placements), default=-1) + 1
    usable_height = max(1.0, float(height_px) - 2.0 * params.padding_px)
    lane_height = 14.0 if track.clipped else usable_height / max(1, lane_count)
    return placements, lane_height


def _label_rows(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    side: str,
    covers_slots: bool,
    height_px: float,
    font_family: str,
    params: AnnotationTrackParams,
) -> list[tuple[ResolvedRegionAnnotation, RegionAnnotationStyle, float, tuple[float, float, float, float]]]:
    """Each labelled region of the record: its style, the label's slot-local ``y`` and its ink offsets.

    The letters keep ``label_offset`` from the painted edge of the mark on the
    side ``linear_annotation_labels_above`` picks (OV-844). A slot that
    ``covers_slots`` (the plan's ``spans_slots``) keeps them outside every slot
    it covers too: ``label_offset`` beyond the slot's edge, or beyond the
    mark's paint where it passes that edge, so the place does not depend on
    the overlay's own height (OV-1156, review 1176). The paint passes the edge
    where ``MarkReach.past_lane_edge`` exceeds ``padding_px``: with lane
    height ``h``, a stroke wider than ``2 * padding_px`` (highlight),
    ``2 * padding_px + 0.16 * h`` (band) or ``h + 2 * padding_px`` (bracket,
    line), or an arrow cap whose ``min(5, 0.35 * h) + stroke_width * sqrt(0.5) / 2``
    exceeds ``h / 2 + padding_px``.
    """

    if not params.show_labels:
        return []
    placements, lane_height = _record_lanes(track, record_index=record_index, height_px=height_px, params=params)
    above = linear_annotation_labels_above(side, params)
    rows = []
    for placement in placements:
        annotation = placement.annotation
        if not annotation.label:
            continue
        style = effective_annotation_style(annotation, params)
        size = annotation_label_font_size(style, track.label_font_size)
        offsets = label_ink_offsets(annotation.label, font_family, size, _label_anchor(style))
        _left, _right, top, bottom = offsets
        # A Linear band fills _band_height_factor of its lane; an arrow cap is a head of two diagonals.
        reach = mark_reach(
            annotation, style, lane_height, fill_share=_band_height_factor(annotation.mark), arrow_heads=True
        )
        if covers_slots:
            # Beyond the slot edge, or beyond the mark's paint where it passes that edge.
            inset = (
                params.padding_px + placement.lane * lane_height
                if above
                else float(height_px) - params.padding_px - (placement.lane + 1) * lane_height
            )
            near = past_span_edge(reach, inset) + float(style.label_offset)
            y = -near - bottom if above else float(height_px) + near - top
        else:
            lane_center = params.padding_px + (placement.lane + 0.5) * lane_height
            near = reach.from_centre + float(style.label_offset)
            y = lane_center - near - bottom if above else lane_center + near - top
        rows.append((annotation, style, y, offsets))
    return rows


def linear_annotation_label_room(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    side: str,
    covers_slots: bool,
    height_px: float,
    font_family: str,
    params: AnnotationTrackParams,
) -> float:
    """How far the record's labels reach beyond the slot's ``height_px`` on their side.

    The room holds the letters and ``label_offset`` beyond them too, so the
    axis line, the next slot or the record edge keeps that gap (OV-844).
    """

    above = linear_annotation_labels_above(side, params)
    room = 0.0
    for _annotation, style, y, (_left, _right, top, bottom) in _label_rows(
        track,
        record_index=record_index,
        side=side,
        covers_slots=covers_slots,
        height_px=height_px,
        font_family=font_family,
        params=params,
    ):
        far = -(y + top) if above else y + bottom - float(height_px)
        room = max(room, far + float(style.label_offset))
    return room


def linear_annotation_labels(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    record_length: int,
    bar_length_px: float,
    side: str,
    covers_slots: bool,
    height_px: float,
    font_family: str,
    params: AnnotationTrackParams,
) -> tuple[LinearAnnotationLabel, ...]:
    """The labels of one record's annotation track at their ``label_position``."""

    scale = float(bar_length_px) / max(1, int(record_length))
    labels: list[LinearAnnotationLabel] = []
    for annotation, style, y, (left, right, top, bottom) in _label_rows(
        track,
        record_index=record_index,
        side=side,
        covers_slots=covers_slots,
        height_px=height_px,
        font_family=font_family,
        params=params,
    ):
        segments = annotation.geometry_segments
        anchor = _label_anchor(style)
        if anchor == "start":
            label_bp = float(segments[0][0])
        elif anchor == "end":
            label_bp = float(segments[-1][1])
        else:
            label_bp = float(annotation.geometry_midpoint_bp)
        x = label_bp * scale
        labels.append(
            LinearAnnotationLabel(
                annotation_id=annotation.id,
                text=annotation.label,
                x=x,
                y=y,
                anchor=anchor,
                font_size=annotation_label_font_size(style, track.label_font_size),
                color=style.label_color,
                ink=Aabb(x + left, y + top, x + right, y + bottom),
                segments=tuple((float(start) * scale, float(end) * scale) for start, end in segments),
            )
        )
    return tuple(labels)


def linear_annotation_track_content(
    track: ResolvedAnnotationTrack,
    *,
    record_id: str,
    record_index: int,
    record_length: int,
    bar_length_px: float,
    y_offset_px: float,
    height_px: float,
    font_family: str,
    params: AnnotationTrackParams,
    labels: Sequence[LinearAnnotationLabel] = (),
    dom_group_id: str | None = None,
    semantic_slot_id: str | None = None,
) -> GroupContent:
    """Build an annotation track in record-local coordinates (render stage ``annotation_content``).

    ``labels`` are the track's placed labels (``linear_annotation_labels``,
    moved by the record's label placement), drawn with their annotations.
    The hatch patterns the marks use, in order of first use, and the clip path
    of a clipped track are the content's definitions.
    """

    group_id = dom_group_id or (
        f"gbdraw-annotation-track-{dom_token(track.slot_id)}-{record_index + 1}"
    )
    group = Group(id=group_id)
    definitions: dict[str, BaseElement] = {}
    if semantic_slot_id:
        group.attribs["data-gbdraw-slot-id"] = str(semantic_slot_id)
        group.attribs["data-gbdraw-slot-renderer"] = "annotations"
    placements, lane_height = _record_lanes(track, record_index=record_index, height_px=height_px, params=params)
    scale = float(bar_length_px) / max(1, int(record_length))
    labels_by_id = {label.annotation_id: label for label in labels}

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
        item_group = Group(id=item_id)
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
        part_count = len(annotation.geometry_segments)
        for part_index, (start, end) in enumerate(annotation.geometry_segments, start=1):
            x1, x2 = float(start) * scale, float(end) * scale
            fragment = annotation.display_parts[part_index - 1] if annotation.display_parts is not None else None
            cap_start, cap_end = capped_ends(annotation, part_index - 1)
            if annotation.mark in {"band", "highlight"}:
                if style.hatch:
                    definitions.setdefault(hatch_pattern_id(style.hatch), hatch_pattern(style.hatch))
                fill = hatch_pattern_paint(style.hatch) if style.hatch else (
                    style.fill or ("#94a3b8" if annotation.mark == "highlight" else "none")
                )
                height_factor = _band_height_factor(annotation.mark)
                rect = Rect(
                    insert=(x1, lane_center - 0.5 * height_factor * lane_height),
                    size=(max(0.1, x2 - x1), height_factor * lane_height),
                    fill=fill,
                    fill_opacity=style.fill_opacity,
                    stroke="none" if auto_feature_underlay else style.stroke,
                    stroke_width=0 if auto_feature_underlay else style.stroke_width,
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
                if not auto_feature_underlay and fragment is not None and (fragment.artificial_start or fragment.artificial_end):
                    rect.attribs['stroke'] = 'none'
                    top = lane_center - 0.5 * height_factor * lane_height
                    bottom = lane_center + 0.5 * height_factor * lane_height
                    outline = f"M {x1:g},{top:g} L {x2:g},{top:g} M {x1:g},{bottom:g} L {x2:g},{bottom:g}"
                    for x, capped in ((x1, cap_start), (x2, cap_end)):
                        if capped:
                            outline += f" M {x:g},{top:g} L {x:g},{bottom:g}"
                    border = Path(d=outline, fill="none", stroke=style.stroke, stroke_width=style.stroke_width, stroke_linecap="butt")
                    if dash:
                        border.attribs["stroke-dasharray"] = dash
                    item_group.add(rect)
                    item_group.add(border)
                else:
                    item_group.add(rect)
                continue
            line = Line(
                start=(x1, lane_center),
                end=(x2, lane_center),
                fill="none",
                stroke=style.stroke,
                stroke_width=style.stroke_width,
                stroke_linecap=("butt" if fragment is not None and (fragment.artificial_start or fragment.artificial_end) else "round"),
            )
            if dash:
                line.attribs["stroke-dasharray"] = dash
            item_group.add(line)
            if annotation.mark == "bracket" and style.line_cap != "none":
                cap = cap_size(lane_height)
                if style.line_cap == "arrow":
                    path = Path(
                        d=f"M {x2 - cap:g},{lane_center - cap:g} L {x2:g},{lane_center:g} "
                        f"L {x2 - cap:g},{lane_center + cap:g}",
                        fill="none",
                        stroke=style.stroke,
                        stroke_width=style.stroke_width,
                    )
                    if cap_end:
                        item_group.add(path)
                    if cap_start:
                        item_group.add(Line(start=(x1, lane_center - cap), end=(x1, lane_center + cap), stroke=style.stroke, stroke_width=style.stroke_width))
                else:
                    if cap_start:
                        item_group.add(Line(start=(x1, lane_center - cap), end=(x1, lane_center + cap), stroke=style.stroke, stroke_width=style.stroke_width))
                    if cap_end:
                        item_group.add(Line(start=(x2, lane_center - cap), end=(x2, lane_center + cap), stroke=style.stroke, stroke_width=style.stroke_width))

        label = labels_by_id.get(annotation.id)
        if label is not None:
            item_group.add(
                Text(
                    label.text,
                    insert=(label.x, float(y_offset_px) + label.y),
                    fill=label.color,
                    stroke="none",
                    font_size=label.font_size,
                    font_family=font_family,
                    text_anchor=label.anchor,
                    dominant_baseline="central",
                )
            )
        group.add(item_group)
    clip_paths: list[BaseElement] = []
    if track.clipped:
        group.attribs["data-gbdraw-annotation-clipped"] = "true"
        clip_id = (
            stable_svg_id("annotation_clip", "linear", group_id)
            if dom_group_id
            else f"gbdraw-annotation-clip-{dom_token(track.slot_id)}-{record_index + 1}"
        )
        clip_path = ClipPath(id=clip_id)
        clip_path.add(
            Rect(
                insert=(0, float(y_offset_px)),
                size=(float(bar_length_px), float(height_px)),
            )
        )
        clip_paths.append(clip_path)
        group.attribs["clip-path"] = f"url(#{clip_id})"
    return GroupContent.of(group, (*definitions.values(), *clip_paths))


__all__ = [
    "LinearAnnotationLabel",
    "linear_annotation_label_room",
    "linear_annotation_labels",
    "linear_annotation_track_content",
]
