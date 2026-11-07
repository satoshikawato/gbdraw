"""Linear annotation SVG geometry."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from functools import partial
from typing import Literal, NamedTuple

from svgwrite.base import BaseElement

from gbdraw.svg.elements import ClipPath, Group, Line, Path, Rect, Text

from gbdraw.annotations import (
    AnnotationTrackParams,
    MarkPart,
    RegionAnnotationStyle,
    ResolvedAnnotationTrack,
    ResolvedRegionAnnotation,
    annotation_label_font_size,
    apply_feature_underlay_dom_attributes,
    effective_annotation_style,
    is_auto_feature_underlay,
    linear_annotation_labels_above,
    mark_parts,
)
from gbdraw.annotations.layout import PlacedRegionAnnotation
from gbdraw.layout.linear import VerticalBand
from gbdraw.layout.spatial import Aabb
from gbdraw.layout.text_geometry import TextAnchor
from gbdraw.render.drawers.annotation_marks import (
    ANNOTATION_ID_ATTRIBUTE,
    ANNOTATION_SET_ID_ATTRIBUTE,
    annotation_dom_id,
    cap_size,
    dom_token,
    gradient_fill_group,
    gradient_steps,
    label_ink_offsets,
    mark_reach,
    past_span_edge,
)
from gbdraw.render.drawers.label_backing import backing_box, backing_rect
from gbdraw.render.drawers.linear.wrap_rows import wrap_row_group
from gbdraw.render.group_content import GroupContent
from gbdraw.render.patterns import hatch_pattern, hatch_pattern_id, hatch_pattern_paint
from gbdraw.svg.ids import stable_svg_id
from gbdraw.svg.path_data import path_number


@dataclass(frozen=True)
class LinearAnnotationLabel:
    """One Linear annotation label as it is drawn (``dominant-baseline=central``).

    ``x`` and ``y`` anchor the text in record-local x and slot-local y;
    ``ink`` is the box of its painted letters in the same frame, and
    ``segments`` the x ranges of its marks: the label's centre may move along
    them, never into a gap between them. ``mark_edge`` is the y, on the
    letters' side, where the paint of its mark (or the slot it keeps clear of)
    begins, ``label_offset`` from the letters; the ``backing`` never crosses it
    (OV-1290). ``boxes`` are the filled rectangles of a band or highlight
    mark, where the label may sit instead (LP-D75). A ``backed`` label is
    drawn in front of the record's leader lines on its ``backing`` (LP-D82).
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
    mark_edge: float
    boxes: tuple[Aabb, ...] = ()
    backed: bool = False

    @property
    def backing(self) -> Aabb:
        """The box a backed label is drawn on: the one geometry its drawer and its count read. It
        stops at ``mark_edge``, so its padding never reaches past ``label_offset`` into its own
        mark or the feature glyphs under it."""

        return backing_box(self.ink, self.font_size, edge_y=self.mark_edge)

    @property
    def footprint(self) -> Aabb:
        """What the drawn label covers: its backing when it has one, else its letters."""

        return self.backing if self.backed else self.ink

    def translated(self, dx: float, dy: float) -> LinearAnnotationLabel:
        return replace(self, x=self.x + dx, y=self.y + dy, ink=self.ink.translated(dx, dy))

    def moved_to(self, x: float) -> LinearAnnotationLabel:
        return replace(self, x=x, ink=self.ink.translated(x - self.x, 0.0))


class AnnotationRowFrame(NamedTuple):
    """A slot on one row of a wrapped record: how far it lies below the slot's first-row frame, and its height there."""

    shift_y: float
    height_px: float


@dataclass(frozen=True)
class AnnotationSlotFrame:
    """An annotation slot as its marks are drawn, slot-local: its height on the first row, its frame on each row of a
    wrapped record (none on an unwrapped one), and each strand's lane of the feature track it covers (none when it
    covers no track that separates the strands; the same on every row)."""

    height_px: float
    row_frames: tuple[AnnotationRowFrame, ...] = ()
    strand_bands: Mapping[Literal["+", "-"], VerticalBand] = field(default_factory=dict)

    def clip_path(self, clip_id: str, *, y_offset_px: float, bar_length_px: float) -> ClipPath:
        """The clip of a clipped track: one rectangle per row, the slot's height on each of them."""

        clip_path = ClipPath(id=clip_id)
        for frame in self.row_frames or (AnnotationRowFrame(0.0, float(self.height_px)),):
            clip_path.add(
                Rect(
                    insert=(0, float(y_offset_px) + frame.shift_y),
                    size=(float(bar_length_px), frame.height_px),
                )
            )
        return clip_path


def _band_height_factor(mark: str) -> float:
    """The share of its lane a ``highlight`` or ``band`` rectangle fills."""

    return 1.0 if mark == "highlight" else 0.84


def _lane_center(placement: PlacedRegionAnnotation, lane_height: float, padding_px: float) -> float:
    """The slot-local y of the middle of a placement's lane."""

    return padding_px + (placement.lane + 0.5) * lane_height


def _band_rect(mark: str, lane_center: float, lane_height: float, x1: float, x2: float) -> tuple[float, float, float, float]:
    """The ``(x, y, width, height)`` of one part of a ``band`` or ``highlight`` mark."""

    height_factor = _band_height_factor(mark)
    return x1, lane_center - 0.5 * height_factor * lane_height, max(0.1, x2 - x1), height_factor * lane_height


def _solid_fill(annotation: ResolvedRegionAnnotation, style: RegionAnnotationStyle) -> str:
    """The fill colour of a ``band`` or ``highlight`` mark without hatching (``none`` for none)."""

    return style.fill or ("#94a3b8" if annotation.mark == "highlight" else "none")


def _band_paint(
    target: Group,
    definitions: dict[str, BaseElement],
    annotation: ResolvedRegionAnnotation,
    style: RegionAnnotationStyle,
    part: MarkPart,
    start: float,
    end: float,
    record_length: int,
    scale: float,
    band: tuple[float, float],
) -> str:
    """The fill of one part of a band or highlight mark: its hatch pattern (added to ``definitions``), its solid
    fill, or ``none`` over the steps of its gradient (added to ``target`` before the mark's outline). ``band``
    is the part's top and height."""

    if style.hatch:
        definitions.setdefault(hatch_pattern_id(style.hatch), hatch_pattern(style.hatch))
        return hatch_pattern_paint(style.hatch)
    if style.gradient is None:
        return _solid_fill(annotation, style)
    top, height = band
    x2 = float(end) * scale
    steps = gradient_steps(
        annotation, style.gradient, part, start, end, record_length=record_length, length_px=x2 - start * scale
    )
    target.add(gradient_fill_group(steps, style.fill_opacity, lambda step, color: Rect(
        insert=(float(step) * scale, top), size=(max(0.1, x2 - float(step) * scale), height), fill=color, stroke="none",
    )))
    return "none"


def _label_boxes(
    placement: PlacedRegionAnnotation,
    style: RegionAnnotationStyle,
    *,
    lane_height: float,
    padding_px: float,
    scale: float,
    segments: Sequence[tuple[int, int]],
) -> tuple[Aabb, ...]:
    """The filled rectangles of a mark a label may sit inside (LP-D75): a band or highlight
    with a solid fill; a hatched, unfilled, bracket or line mark has none."""

    annotation = placement.annotation
    if annotation.mark not in {"band", "highlight"} or style.hatch or _solid_fill(annotation, style) == "none":
        return ()
    lane_center = _lane_center(placement, lane_height, padding_px)
    boxes: list[Aabb] = []
    for start, end in segments:
        x, y, width, height = _band_rect(annotation.mark, lane_center, lane_height, float(start) * scale, float(end) * scale)
        boxes.append(Aabb(x, y, x + width, y + height))
    return tuple(boxes)


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
    return placements, _lane_height(track, placements, height_px, params)


def _lane_height(
    track: ResolvedAnnotationTrack,
    placements: Sequence[PlacedRegionAnnotation],
    height_px: float,
    params: AnnotationTrackParams,
) -> float:
    """The height of one lane of the record's ``placements`` in a slot ``height_px`` high."""

    lane_count = max((item.lane for item in placements), default=-1) + 1
    usable_height = max(1.0, float(height_px) - 2.0 * params.padding_px)
    return 14.0 if track.clipped else usable_height / max(1, lane_count)


def _mark_lane(
    track: ResolvedAnnotationTrack,
    placement: PlacedRegionAnnotation,
    placements: Sequence[PlacedRegionAnnotation],
    *,
    y_offset_px: float,
    frame: AnnotationSlotFrame,
    lane_height: float,
    params: AnnotationTrackParams,
) -> tuple[float, float]:
    """The middle y of a placement's lane and the lane's height: in its strand's lane of ``frame.strand_bands`` (a
    slot that covers a feature track that separates the strands), else in the slot's band."""

    strand = placement.annotation.display_strand
    band = frame.strand_bands.get(strand) if strand is not None else None
    if band is not None:
        y_offset_px = float(y_offset_px) + band.top_y
        lane_height = _lane_height(track, placements, band.height, params)
    return float(y_offset_px) + _lane_center(placement, lane_height, params.padding_px), lane_height


def _label_rows(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    side: str,
    covers_slots: bool,
    height_px: float,
    font_family: str,
    params: AnnotationTrackParams,
    wrap_row: int | None = None,
) -> list[tuple[PlacedRegionAnnotation, RegionAnnotationStyle, float, tuple[float, float, float, float]]]:
    """Each labelled region of the record (of one row when ``wrap_row`` is given): its style, the label's
    slot-local ``y`` and its ink offsets.

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
        if not annotation.label or (wrap_row is not None and annotation.label_row(join_origin=False) != wrap_row):
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
            lane_center = _lane_center(placement, lane_height, params.padding_px)
            near = reach.from_centre + float(style.label_offset)
            y = lane_center - near - bottom if above else lane_center + near - top
        rows.append((placement, style, y, offsets))
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
    wrap_row: int | None = None,
) -> float:
    """How far the record's labels (one row's, when ``wrap_row`` is given) reach beyond the slot's ``height_px``
    on their side.

    The room holds the letters and ``label_offset`` beyond them too, so the
    axis line, the next slot or the record edge keeps that gap (OV-844).
    """

    above = linear_annotation_labels_above(side, params)
    room = 0.0
    for _placement, style, y, (_left, _right, top, bottom) in _label_rows(
        track,
        record_index=record_index,
        side=side,
        covers_slots=covers_slots,
        height_px=height_px,
        font_family=font_family,
        params=params,
        wrap_row=wrap_row,
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
    frame: AnnotationSlotFrame,
    font_family: str,
    params: AnnotationTrackParams,
) -> tuple[LinearAnnotationLabel, ...]:
    """The labels of one record's annotation track at their ``label_position``, in the slot's ``frame``, the one
    its marks are drawn in (``linear_annotation_track_content``).

    On a wrapped record (``frame.row_frames``) a label stands over its marks on the row of its longest one, its y
    moved from that row's frame into the slot's first-row frame.
    """

    height_px = frame.height_px
    row_frames = frame.row_frames

    scale = float(bar_length_px) / max(1, int(record_length))
    placements, _lane = _record_lanes(track, record_index=record_index, height_px=height_px, params=params)
    above = linear_annotation_labels_above(side, params)
    labels: list[LinearAnnotationLabel] = []
    candidates = (
        [
            (row_item, shift)
            for row, (shift, row_height) in enumerate(row_frames)
            for row_item in _label_rows(
                track,
                record_index=record_index,
                side=side,
                covers_slots=covers_slots,
                height_px=row_height,
                font_family=font_family,
                params=params,
                wrap_row=row,
            )
        ]
        if row_frames
        else [
            (row_item, 0.0)
            for row_item in _label_rows(
                track,
                record_index=record_index,
                side=side,
                covers_slots=covers_slots,
                height_px=height_px,
                font_family=font_family,
                params=params,
            )
        ]
    )
    for (placement, style, row_y, (left, right, top, bottom)), shift_y in candidates:
        annotation = placement.annotation
        row = annotation.label_row(join_origin=False)
        y = row_y + shift_y
        lane_height = _lane_height(
            track, placements, row_frames[row].height_px if row_frames else height_px, params
        )
        segments = annotation.row_segments(row)
        anchor = _label_anchor(style)
        # The row's parts lie in the order of the region's segments, which a display start can turn round.
        if anchor == "start":
            label_bp = float(min(start for start, _end in segments))
        elif anchor == "end":
            label_bp = float(max(end for _start, end in segments))
        else:
            label_bp = annotation.label_midpoint_bp(join_origin=False)
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
                mark_edge=y + bottom + style.label_offset if above else y + top - style.label_offset,
                # A clipped track may hide part of its marks.
                boxes=()
                if track.clipped
                else tuple(
                    box.translated(0.0, shift_y)
                    for box in _label_boxes(
                        placement, style, lane_height=lane_height, padding_px=params.padding_px, scale=scale, segments=segments
                    )
                ),
            )
        )
    return tuple(labels)


def _track_group_id(track: ResolvedAnnotationTrack, record_index: int, dom_group_id: str | None) -> str:
    return dom_group_id or f"gbdraw-annotation-track-{dom_token(track.slot_id)}-{record_index + 1}"


def _label_text(label: LinearAnnotationLabel, font_family: str, y_offset_px: float = 0.0) -> Text:
    return Text(
        label.text,
        insert=(label.x, y_offset_px + label.y),
        fill=label.color,
        stroke="none",
        font_size=label.font_size,
        font_family=font_family,
        text_anchor=label.anchor,
        dominant_baseline="central",
    )


def linear_annotation_track_content(
    track: ResolvedAnnotationTrack,
    *,
    record_id: str,
    record_index: int,
    record_length: int,
    bar_length_px: float,
    y_offset_px: float,
    frame: AnnotationSlotFrame,
    font_family: str,
    params: AnnotationTrackParams,
    labels: Sequence[LinearAnnotationLabel] = (),
    dom_group_id: str | None = None,
    semantic_slot_id: str | None = None,
) -> GroupContent:
    """Build an annotation track in record-local coordinates (render stage ``annotation_content``).

    On a wrapped record (``frame.row_frames``) each part of a mark goes into an id-less group of its row in its
    item. A region with a strand is drawn in its strand's lane of ``frame.strand_bands``, any other in the slot's
    band.

    ``labels`` are the track's placed labels (``linear_annotation_labels``,
    moved by the record's label placement), drawn with their annotations.
    The hatch patterns the marks use, in order of first use, and the clip path
    of a clipped track are the content's definitions.
    """

    group_id = _track_group_id(track, record_index, dom_group_id)
    group = Group(id=group_id)
    definitions: dict[str, BaseElement] = {}
    if semantic_slot_id:
        group.attribs["data-gbdraw-slot-id"] = str(semantic_slot_id)
        group.attribs["data-gbdraw-slot-renderer"] = "annotations"
    placements, slot_lane_height = _record_lanes(
        track, record_index=record_index, height_px=frame.height_px, params=params
    )
    scale = float(bar_length_px) / max(1, int(record_length))
    labels_by_id = {label.annotation_id: label for label in labels}

    for placement in placements:
        annotation = placement.annotation
        auto_feature_underlay = is_auto_feature_underlay(annotation)
        style = effective_annotation_style(annotation, params)
        lane = partial(_mark_lane, track, placement, placements, y_offset_px=y_offset_px, frame=frame, params=params)
        lane_center, lane_height = lane(lane_height=slot_lane_height)
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
                    ANNOTATION_ID_ATTRIBUTE: annotation.id,
                    ANNOTATION_SET_ID_ATTRIBUTE: track.set_id,
                    "data-gbdraw-annotation-track-id": track.slot_id,
                    "data-gbdraw-record-id": record_id,
                    "data-gbdraw-record-index": str(record_index),
                    "data-gbdraw-annotation-mark": annotation.mark,
                    "data-gbdraw-annotation-label": annotation.label,
                }
            )
        dash = ",".join(f"{value:g}" for value in style.stroke_dasharray) or None
        parts = mark_parts(annotation, join_origin=False)
        part_count = len(parts)
        row_targets: dict[int, Group] = {}
        for part_index, part in enumerate(parts, start=1):
            x1, x2 = float(part.start) * scale, float(part.end) * scale
            cap_start, cap_end = part.capped_start, part.capped_end
            part_target = item_group
            if frame.row_frames and annotation.display_parts is not None:
                row = frame.row_frames[part.row]
                part_target = row_targets.get(part.row) or wrap_row_group(part.row, row.shift_y)
                if part.row not in row_targets:
                    row_targets[part.row] = part_target
                    item_group.add(part_target)
                lane_center, lane_height = lane(lane_height=_lane_height(track, placements, row.height_px, params))
            if annotation.mark in {"band", "highlight"}:
                height_factor = _band_height_factor(annotation.mark)
                rect_x, rect_y, rect_width, rect_height = _band_rect(annotation.mark, lane_center, lane_height, x1, x2)
                fill = _band_paint(part_target, definitions, annotation, style, part, part.start, part.end,
                                   record_length, scale, (rect_y, rect_height))
                rect = Rect(
                    insert=(rect_x, rect_y),
                    size=(rect_width, rect_height),
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
                if not auto_feature_underlay and not (cap_start and cap_end):
                    rect.attribs['stroke'] = 'none'
                    top = lane_center - 0.5 * height_factor * lane_height
                    bottom = lane_center + 0.5 * height_factor * lane_height
                    outline = f"M {path_number(x1)},{path_number(top)} L {path_number(x2)},{path_number(top)} M {path_number(x1)},{path_number(bottom)} L {path_number(x2)},{path_number(bottom)}"
                    for x, capped in ((x1, cap_start), (x2, cap_end)):
                        if capped:
                            outline += f" M {path_number(x)},{path_number(top)} L {path_number(x)},{path_number(bottom)}"
                    border = Path(d=outline, fill="none", stroke=style.stroke, stroke_width=style.stroke_width, stroke_linecap="butt")
                    if dash:
                        border.attribs["stroke-dasharray"] = dash
                    part_target.add(rect)
                    part_target.add(border)
                else:
                    part_target.add(rect)
                continue
            line = Line(
                start=(x1, lane_center),
                end=(x2, lane_center),
                fill="none",
                stroke=style.stroke,
                stroke_width=style.stroke_width,
                stroke_linecap="round" if cap_start and cap_end else "butt",
            )
            if dash:
                line.attribs["stroke-dasharray"] = dash
            part_target.add(line)
            if annotation.mark == "bracket" and style.line_cap != "none":
                cap = cap_size(lane_height)
                if style.line_cap == "arrow":
                    path = Path(
                        d=f"M {path_number(x2 - cap)},{path_number(lane_center - cap)} L {path_number(x2)},{path_number(lane_center)} "
                        f"L {path_number(x2 - cap)},{path_number(lane_center + cap)}",
                        fill="none",
                        stroke=style.stroke,
                        stroke_width=style.stroke_width,
                    )
                    if cap_end:
                        part_target.add(path)
                    if cap_start:
                        part_target.add(Line(start=(x1, lane_center - cap), end=(x1, lane_center + cap), stroke=style.stroke, stroke_width=style.stroke_width))
                else:
                    if cap_start:
                        part_target.add(Line(start=(x1, lane_center - cap), end=(x1, lane_center + cap), stroke=style.stroke, stroke_width=style.stroke_width))
                    if cap_end:
                        part_target.add(Line(start=(x2, lane_center - cap), end=(x2, lane_center + cap), stroke=style.stroke, stroke_width=style.stroke_width))

        label = labels_by_id.get(annotation.id)
        # A backed label is drawn in front of the leader lines (linear_region_label_front_content).
        if label is not None and not label.backed:
            item_group.add(_label_text(label, font_family, float(y_offset_px)))
        group.add(item_group)
    clip_paths: list[BaseElement] = []
    if track.clipped:
        group.attribs["data-gbdraw-annotation-clipped"] = "true"
        clip_id = (
            stable_svg_id("annotation_clip", "linear", group_id)
            if dom_group_id
            else f"gbdraw-annotation-clip-{dom_token(track.slot_id)}-{record_index + 1}"
        )
        clip_paths.append(frame.clip_path(clip_id, y_offset_px=y_offset_px, bar_length_px=bar_length_px))
        group.attribs["clip-path"] = f"url(#{clip_id})"
    return GroupContent.of(group, (*definitions.values(), *clip_paths))


def linear_region_label_front_content(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    font_family: str,
    labels: Sequence[LinearAnnotationLabel],
    dom_group_id: str | None = None,
) -> GroupContent | None:
    """The track's backed labels in front of the record's leader lines (LP-D82 ``RGN-FRONT``), in
    the track's slot-local frame; ``None`` when it has none.

    Each label is drawn on its ``backing``, a white box with rounded corners at
    70 % opacity. Its group names the annotation's item group
    (``data-gbdraw-region-label-for``), where the Web selects the annotation.
    """

    backed = [label for label in labels if label.backed]
    if not backed:
        return None
    group = Group(id=f"{_track_group_id(track, record_index, dom_group_id)}-front")
    for label in backed:
        item = Group()
        item.attribs["data-gbdraw-region-label-for"] = annotation_dom_id(
            record_index=record_index,
            slot_id=track.slot_id,
            set_id=track.set_id,
            annotation_id=label.annotation_id,
        )
        item.add(backing_rect(label.backing, label.font_size))
        item.add(_label_text(label, font_family))
        group.add(item)
    return GroupContent.of(group)


__all__ = [
    "AnnotationRowFrame",
    "AnnotationSlotFrame",
    "LinearAnnotationLabel",
    "linear_annotation_label_room",
    "linear_annotation_labels",
    "linear_annotation_track_content",
    "linear_region_label_front_content",
]
