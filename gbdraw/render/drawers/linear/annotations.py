"""Linear annotation SVG geometry."""

from __future__ import annotations

import hashlib
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace

from svgwrite.base import BaseElement

from gbdraw.svg.elements import ClipPath, Group, Line, Path, Rect, Text

from gbdraw.annotations import (
    AnnotationTrackParams,
    ResolvedAnnotationTrack,
    ResolvedRegionAnnotation,
    apply_feature_underlay_dom_attributes,
    effective_annotation_style,
    is_auto_feature_underlay,
    linear_annotation_labels_above,
)
from gbdraw.annotations.layout import PlacedRegionAnnotation
from gbdraw.core.text import dominant_baseline_shift_em, get_font_vertical_metrics, text_ink
from gbdraw.layout.spatial import Aabb
from gbdraw.layout.text_geometry import TextAnchor
from gbdraw.render.group_content import GroupContent
from gbdraw.render.patterns import hatch_pattern, hatch_pattern_id, hatch_pattern_paint
from gbdraw.svg.ids import stable_svg_id


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


_LABEL_ANCHOR_SHARE: Mapping[TextAnchor, float] = {"start": 0.0, "middle": 0.5, "end": 1.0}
_DEFAULT_LABEL_FONT_SIZE = 10.0


@dataclass(frozen=True)
class LinearAnnotationLabel:
    """One Linear annotation label as it is drawn (``dominant-baseline=central``).

    ``x`` and ``y`` anchor the text in record-local x and slot-local y;
    ``ink`` is the box of its painted letters in the same frame, and ``span``
    the x range of its mark, along which the label may move.
    """

    annotation_id: str
    text: str
    x: float
    y: float
    anchor: TextAnchor
    font_size: float
    color: str
    ink: Aabb
    span: tuple[float, float]

    def moved_to(self, x: float) -> LinearAnnotationLabel:
        return replace(self, x=x, ink=self.ink.translated(x - self.x, 0.0))


def _label_ink_offsets(
    text: str, font_family: str, font_size: float
) -> tuple[float, float, float, float, float]:
    """The advance of ``text`` and its ink box from the start of the advance and its central baseline:
    advance, left, right, top, bottom."""

    ink = text_ink(text, font_family)
    metrics = get_font_vertical_metrics(font_family)
    shift = None if metrics is None else dominant_baseline_shift_em("central", metrics)
    if ink is None or shift is None:
        # A font that cannot be read: an em box with the fallback width of calculate_bbox_dimensions.
        width = 0.6 * len(text) * font_size
        return width, 0.0, width, -0.5 * font_size, 0.5 * font_size
    return (
        ink.advance * font_size,
        ink.left * font_size,
        ink.right * font_size,
        (shift - ink.top) * font_size,
        (shift - ink.bottom) * font_size,
    )


def _label_clearance(annotation: ResolvedRegionAnnotation, params: AnnotationTrackParams) -> float:
    """How far a label's letters keep from its lane edge: half the mark's stroke plus ``label_offset``."""

    style = effective_annotation_style(annotation, params)
    half_stroke = 0.0 if is_auto_feature_underlay(annotation) else 0.5 * float(style.stroke_width)
    return half_stroke + float(style.label_offset)


def linear_annotation_label_room(
    annotation: ResolvedRegionAnnotation,
    params: AnnotationTrackParams,
    font_family: str,
) -> float:
    """How far a label's painted letters reach beyond its lane edge (0 without a drawn label)."""

    if not (params.show_labels and annotation.label):
        return 0.0
    style = effective_annotation_style(annotation, params)
    size = float(style.label_font_size or _DEFAULT_LABEL_FONT_SIZE)
    _advance, _left, _right, top, bottom = _label_ink_offsets(annotation.label, font_family, size)
    return _label_clearance(annotation, params) + (bottom - top)


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


def linear_annotation_labels(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    record_length: int,
    bar_length_px: float,
    side: str,
    height_px: float,
    font_family: str,
    params: AnnotationTrackParams,
) -> tuple[LinearAnnotationLabel, ...]:
    """The labels of one record's annotation track at their ``label_position`` (OV-844).

    A label lies beyond its lane on the side ``linear_annotation_labels_above``
    picks, with its painted letters ``label_offset`` clear of the mark's
    stroke (``linear_annotation_label_room``).
    """

    if not params.show_labels:
        return ()
    placements, lane_height = _record_lanes(track, record_index=record_index, height_px=height_px, params=params)
    scale = float(bar_length_px) / max(1, int(record_length))
    above = linear_annotation_labels_above(side, params)
    labels: list[LinearAnnotationLabel] = []
    for placement in placements:
        annotation = placement.annotation
        if not annotation.label:
            continue
        style = effective_annotation_style(annotation, params)
        size = float(style.label_font_size or _DEFAULT_LABEL_FONT_SIZE)
        advance, left, right, top, bottom = _label_ink_offsets(annotation.label, font_family, size)
        lane_top = params.padding_px + placement.lane * lane_height
        clearance = _label_clearance(annotation, params)
        y = lane_top - clearance - bottom if above else lane_top + lane_height + clearance - top
        segments = annotation.geometry_segments
        anchor: TextAnchor
        if style.label_position == "start":
            label_bp, anchor = float(segments[0][0]), "start"
        elif style.label_position == "end":
            label_bp, anchor = float(segments[-1][1]), "end"
        else:
            label_bp, anchor = float(annotation.geometry_midpoint_bp), "middle"
        x = label_bp * scale
        origin = x - _LABEL_ANCHOR_SHARE[anchor] * advance
        labels.append(
            LinearAnnotationLabel(
                annotation_id=annotation.id,
                text=annotation.label,
                x=x,
                y=y,
                anchor=anchor,
                font_size=size,
                color=style.label_color,
                ink=Aabb(origin + left, y + top, origin + right, y + bottom),
                span=(
                    min(float(start) for start, _end in segments) * scale,
                    max(float(end) for _start, end in segments) * scale,
                ),
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
        f"gbdraw-annotation-track-{_dom_token(track.slot_id)}-{record_index + 1}"
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
            cap_start = fragment is None or not fragment.artificial_start
            cap_end = fragment is None or not fragment.artificial_end
            if annotation.mark in {"band", "highlight"}:
                if style.hatch:
                    definitions.setdefault(hatch_pattern_id(style.hatch), hatch_pattern(style.hatch))
                fill = hatch_pattern_paint(style.hatch) if style.hatch else (
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
                cap = min(5.0, 0.35 * lane_height)
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
            else f"gbdraw-annotation-clip-{_dom_token(track.slot_id)}-{record_index + 1}"
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
    "annotation_dom_id",
    "linear_annotation_label_room",
    "linear_annotation_labels",
    "linear_annotation_track_content",
]
