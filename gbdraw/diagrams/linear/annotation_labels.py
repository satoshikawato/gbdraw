"""Where a Linear record's annotation labels go along their spans (OV-845).

An annotation label keeps the place ``linear_annotation_labels`` gives it
when its painted letters are clear of the record's feature labels, their
leaders, the record names and definition text, the other records of its row
and the annotation labels placed before it. Otherwise it moves, with its
centre over one of its marks (never into a gap between them), to the nearest
clear place; with none, to the place with the fewest conflicts, nearest
first. Its letters may hang past either end of the bar into empty space (the
figure's bounds follow them).

When that place still overlaps something, the label goes into the empty area
inside its own filled mark (LP-D75) when that clears every conflict (LP-D77).
Otherwise, at that place, the record's other labels make way (LP-D82), only
when that clears every conflict of the label and makes no new overlap: each
feature label whose letters it meets moves to another row of its side, the
nearest row first, then one new row beyond the others (``RGN-YIELD-ROW``);
and a label that then meets only the record's leader lines is drawn in front
of them on a white backing (``RGN-FRONT``), where a leader behind it is no
conflict. A move to a new row needs room the record's plan does not have:
the figure proposes it (``PlacedAnnotationLabels.new_rows``) and is assembled
again with it. The font never changes, and the label keeps one line. A
conflict that remains is drawn and counted: one pair per label and each thing
it overlaps (a feature label with its leader, a text band, a bar, an
annotation label), which the label fit warns about (OV-1155, LP-D68).
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections.abc import Mapping, Sequence
from dataclasses import replace
from typing import NamedTuple, TypeAlias, cast

from gbdraw.annotations import ResolvedAnnotationTrack, annotation_track_params_from_mapping
from gbdraw.annotations.label_search import nearest_fewest_conflicts
from gbdraw.core.text import text_ink_offsets
from gbdraw.features.objects import FeatureObject
from gbdraw.labels.fit import LabelFit
from gbdraw.labels.linear import (
    FeatureLabelRow,
    RecordRegionLabels,
    calculate_label_bounds,
    linear_label_in_row,
    linear_label_row,
    linear_label_side,
    linear_leader_polyline,
    linear_row_label_fit,
)
from gbdraw.layout.linear import LinearFeatureLaneGeometry
from gbdraw.layout.spatial import Aabb
from gbdraw.render.drawers.linear.annotations import LinearAnnotationLabel, linear_annotation_labels

from .track_slots import LinearRecordVerticalPlan

AnnotationLabelsBySlot = Mapping[str, tuple[LinearAnnotationLabel, ...]]
# What drew an obstacle: (record index, kind, index); a feature label's box and leader share one.
ObstacleOwner: TypeAlias = tuple[int, str, int]


class LabelObstacle(NamedTuple):
    box: Aabb
    owner: ObstacleOwner
    # A piece of a feature label's leader line, which a backed region label may cover (LP-D82).
    leader: bool = False


class PlacedAnnotationLabels(NamedTuple):
    """One record's annotation labels by slot, its region labels as the label fit counts them,
    the feature labels moved to another row for them (``rows``), and the moves that need a new
    row (``new_rows``, with the region label pairs they would clear, ``new_row_overlaps``, and
    the region labels that asked for them by placement order, ``proposers``)."""

    by_slot: dict[str, tuple[LinearAnnotationLabel, ...]]
    region_labels: RecordRegionLabels
    rows: tuple[FeatureLabelRow, ...] = ()
    new_rows: tuple[FeatureLabelRow, ...] = ()
    new_row_overlaps: int = 0
    proposers: tuple[int, ...] = ()


def _feature_label_box(label: Mapping[str, object]) -> Aabb:
    """A feature label's painted letters as ``LabelDrawer`` draws them (central baseline); its layout box
    when it is rotated or the font cannot be measured."""

    offsets = (
        text_ink_offsets(
            str(label["label_text"]),
            str(label["font_family"]),
            float(cast(float, label["font_size"])),
            dominant_baseline="central",
            text_anchor=str(label.get("text_anchor", "middle")),
        )
        if float(cast(float, label.get("rotation_deg", 0.0))) == 0.0
        else None
    )
    if offsets is None:
        left, right, top, bottom = calculate_label_bounds(dict(label))
        return Aabb(left, top, right, bottom)
    left, right, top, bottom = offsets
    x = float(cast(float, label["middle_x"]))
    y = float(cast(float, label["middle_y"]))
    return Aabb(x + left, y + top, x + right, y + bottom)


def _feature_label_obstacles(
    label: Mapping[str, object], owner: ObstacleOwner, origin_y: float, leader_half_width: float
) -> list[LabelObstacle]:
    """One feature label's letters and leader line, from its record's bar start and axis."""

    obstacles = [LabelObstacle(_feature_label_box(label).translated(0.0, origin_y), owner)]
    leader = linear_leader_polyline(dict(label))
    for (x1, y1), (x2, y2) in zip(leader, leader[1:]):
        obstacles.append(
            LabelObstacle(Aabb(x1, origin_y + y1, x2, origin_y + y2).expanded(leader_half_width), owner, True)
        )
    return obstacles


def _features_origin_y(plan: LinearRecordVerticalPlan) -> float:
    return next((float(slot.origin_y) for slot in plan.slots if slot.renderer == "features"), 0.0)


def linear_feature_label_obstacles(
    feature_labels: Sequence[dict[str, object]],
    plan: LinearRecordVerticalPlan,
    leader_stroke_width: float,
    record_index: int,
) -> list[LabelObstacle]:
    """A record's feature label boxes and leader lines, from its bar start and axis."""

    origin_y = _features_origin_y(plan)
    return [
        obstacle
        for index, label in enumerate(feature_labels)
        for obstacle in _feature_label_obstacles(
            label, (record_index, "feature", index), origin_y, 0.5 * float(leader_stroke_width)
        )
    ]


def linear_feature_glyph_boxes(
    features: Mapping[str, FeatureObject],
    geometry: LinearFeatureLaneGeometry,
    plan: LinearRecordVerticalPlan,
    *,
    separate_strands: bool,
    scale: float,
    stroke_width: float,
) -> list[Aabb]:
    """Each feature part's lane across its x range and stroke, from the record's bar start and
    axis: the ink a label inside a filled mark keeps clear of (LP-D75)."""

    origin_y = next((float(slot.origin_y) for slot in plan.slots if slot.renderer == "features"), None)
    if origin_y is None:
        return []
    half_stroke = 0.5 * stroke_width
    boxes: list[Aabb] = []
    for feature in features.values():
        band = geometry.lane_for(
            strand=feature.strand, track_id=feature.feature_track_id, separate_strands=separate_strands
        ).band.translate(origin_y)
        parts = (
            feature.location if feature.display_parts is None else [part.location for part in feature.display_parts]
        )
        boxes.extend(
            Aabb(part.start * scale - half_stroke, band.top_y, part.end * scale + half_stroke, band.bottom_y)
            for part in parts
        )
    return boxes


def _place_along_span(label: LinearAnnotationLabel, origin_y: float, obstacles: Sequence[LabelObstacle], gap: float) -> float:
    """The anchor x for ``label`` (in a slot at ``origin_y``).

    Its own x when what it covers (its letters, or its backing) is clear;
    else, with the label's centre over one of its marks, the x
    ``nearest_fewest_conflicts`` picks. A sorted sweep keeps the conflict
    count O(log B) in the B obstacles in the label's row.
    """

    ink = label.footprint.translated(0.0, origin_y)
    left, right = ink.min_x - label.x, ink.max_x - label.x
    centre = 0.5 * (left + right)
    # Open intervals of anchor x that put the letters within ``gap`` of an obstacle in their row.
    blocked = [
        (box.min_x - gap - right, box.max_x + gap - left)
        for box, *_rest in obstacles
        if box.min_y < ink.max_y + gap and box.max_y > ink.min_y - gap
    ]
    starts = sorted(start for start, _end in blocked)
    ends = sorted(end for _start, end in blocked)

    def conflicts(x: float) -> int:
        # Intervals with start < x, less those already closed (end <= x): start < x < end.
        return bisect_left(starts, x) - bisect_right(ends, x)

    allowed = [(low - centre, high - centre) for low, high in label.segments]
    candidates: set[float] = set()
    for low, high in allowed:
        candidates.update((low, high))
        candidates.update(value for value in (*starts, *ends) if low <= value <= high)
    return nearest_fewest_conflicts(label.x, candidates, conflicts)


def _place_inside(
    label: LinearAnnotationLabel,
    origin_y: float,
    obstacles: Sequence[LabelObstacle],
    hard: Sequence[Aabb],
    gap: float,
) -> LinearAnnotationLabel | None:
    """``label`` within one of its filled boxes, ``gap`` clear of the box edges, the ``hard`` ink
    and every obstacle (LP-D75), or ``None`` where no such place exists.

    Of the clear places, the one nearest the label's own x, then centred in
    the empty height it occupies, nearest the box's middle."""

    ink = label.ink.translated(0.0, origin_y)
    best: tuple[tuple[float, float], float, float] | None = None
    blockers = [*hard, *(box for box, *_rest in obstacles)]
    for box in (item.translated(0.0, origin_y) for item in label.boxes):
        low_x, high_x = box.min_x + gap, box.max_x - gap - ink.width
        low_y, high_y = box.min_y + gap, box.max_y - gap
        if high_x < low_x or high_y - low_y < ink.height:
            continue
        near = [item for item in blockers if item.intersects(box, gap)]
        want = min(max(ink.min_x, low_x), high_x)
        lefts = {low_x, high_x, want}
        lefts.update(
            value
            for item in near
            for value in (item.max_x + gap, item.min_x - gap - ink.width)
            if low_x <= value <= high_x
        )
        middle = 0.5 * (box.min_y + box.max_y)
        for left in sorted(lefts):
            column = sorted(
                (item.min_y - gap, item.max_y + gap)
                for item in near
                if item.min_x < left + ink.width + gap and item.max_x > left - gap
            )
            cursor, spaces = low_y, []
            for top, bottom in column:
                if top > cursor:
                    spaces.append((cursor, min(top, high_y)))
                cursor = max(cursor, bottom)
            if cursor < high_y:
                spaces.append((cursor, high_y))
            for top, bottom in spaces:
                if bottom - top < ink.height:
                    continue
                centre_y = 0.5 * (top + bottom)
                key = (abs(left - want), abs(centre_y - middle))
                if best is None or key < best[0]:
                    best = (key, left - ink.min_x, centre_y - 0.5 * (ink.min_y + ink.max_y))
    return None if best is None else label.translated(best[1], best[2])


def _hits(box: Aabb, obstacles: Sequence[LabelObstacle]) -> dict[ObstacleOwner, bool]:
    """The owners of the obstacles ``box`` overlaps, each with whether more than its leader line is hit."""

    hits: dict[ObstacleOwner, bool] = {}
    for obstacle in obstacles:
        if obstacle.box.intersects(box):
            hits[obstacle.owner] = hits.get(obstacle.owner, False) or not obstacle.leader
    return hits


def _conflicts(label: LinearAnnotationLabel, origin_y: float, obstacles: Sequence[LabelObstacle]) -> int:
    """The pairs ``label`` makes: one per obstacle owner its letters overlap."""

    return len(_hits(label.ink.translated(0.0, origin_y), obstacles))


def _track_inks(plan: LinearRecordVerticalPlan, slot_id: str, glyphs: Sequence[Aabb], width: float) -> list[Aabb]:
    """What a label placed off its one-line room keeps clear of: the feature glyphs, the axis
    and the paint of every other slot that is not an annotation track, across the bar."""

    bands = [plan.axis_band] + [
        slot.paint_band
        for slot in plan.slots
        if slot.id != slot_id and slot.paint_band is not None and slot.renderer not in {"features", "spacer", "annotations"}
    ]
    return [*glyphs, *(Aabb(0.0, band.top_y, width, band.bottom_y) for band in bands if band.height > 0.0)]


def _choose_candidate(
    label: LinearAnnotationLabel,
    *,
    origin_y: float,
    obstacles: Sequence[LabelObstacle],
    hard: Sequence[Aabb],
    gap: float,
) -> LinearAnnotationLabel:
    """The normal place along the span, unless it keeps a conflict and the place inside the
    filled box clears every conflict of the label (LP-D75, LP-D77)."""

    normal = label.moved_to(_place_along_span(label, origin_y, obstacles, gap))
    if _conflicts(normal, origin_y, obstacles) == 0:
        return normal
    inside = _place_inside(label, origin_y, obstacles, hard, gap)
    return inside if inside is not None and _conflicts(inside, origin_y, obstacles) == 0 else normal


class _Record(NamedTuple):
    """What the record's feature labels move within: their dicts, their origin and leader width,
    the ink they must not newly meet (feature glyphs, the axis and other tracks), and where the
    label fit counts them: the x of the record's start in its row and the feature labels of the
    row's other records, each with its x (``LinearLabelFrame``)."""

    record_index: int
    labels: list[dict[str, object]]
    origin_y: float
    leader_half_width: float
    inks: tuple[LabelObstacle, ...]
    label_x: float = 0.0
    row_labels: tuple[tuple[Sequence[dict[str, object]], float], ...] = ()

    def label_fit(self, labels: Sequence[dict[str, object]]) -> LabelFit:
        """The label fit's count of the row's feature labels with the record's labels ``labels``."""

        return linear_row_label_fit([*self.row_labels, (labels, self.label_x)])


class _Way(NamedTuple):
    """The other labels making way for a region label (LP-D82)."""

    label: LinearAnnotationLabel
    labels: list[dict[str, object]]
    obstacles: list[LabelObstacle]
    rows: tuple[FeatureLabelRow, ...]
    new_row: bool


def _owners_met(boxes: Sequence[LabelObstacle], owner: ObstacleOwner, others: Sequence[LabelObstacle]) -> set[ObstacleOwner]:
    return {other.owner for other in others if other.owner != owner and any(other.box.intersects(item.box) for item in boxes)}


def _another_row(
    owner: ObstacleOwner,
    region: Aabb,
    record: _Record,
    labels: Sequence[dict[str, object]],
    obstacles: Sequence[LabelObstacle],
    *,
    gap: float,
    new_rows: bool,
) -> tuple[int, bool, dict[str, object], list[LabelObstacle]] | None:
    """The row the feature label ``owner`` moves to, off the region label's letters ``region``:
    another row of its side, the nearest first (the outer one on a tie), then, when ``new_rows``,
    one new row beyond them. A move into an existing row makes no new overlap with the ink and
    obstacles; no move makes the label fit's count of the row's feature labels larger, counted
    where and as the label fit counts them (``_Record.label_fit``)."""

    label = labels[owner[2]]
    side = linear_label_side(label)
    if side is None:
        return None
    current = linear_label_row(label)
    top = max(linear_label_row(item) for item in labels if linear_label_side(item) == side)
    rows = sorted((row for row in range(1, top + 1) if row != current), key=lambda row: (abs(row - current), -row))
    others = [*(item for item in obstacles if item.owner != owner), *record.inks]
    before_boxes = [item for item in obstacles if item.owner == owner]
    before = _owners_met(before_boxes, owner, others)
    fit = record.label_fit(labels)
    for row in [*rows, *([top + 1] if new_rows else [])]:
        moved = linear_label_in_row(dict(label), side, row)
        boxes = _feature_label_obstacles(moved, owner, record.origin_y, record.leader_half_width)
        room = row > top
        # A new row is placed by the next assembly, whose plan makes its room: that one checks it.
        if not room and (boxes[0].box.intersects(region, gap) or not _owners_met(boxes, owner, others) <= before):
            continue
        trial = [moved if index == owner[2] else item for index, item in enumerate(labels)]
        after = record.label_fit(trial)
        if after.overlap_pairs > fit.overlap_pairs or after.crossing_pairs > fit.crossing_pairs:
            continue
        return row, room, moved, boxes
    return None


def _make_way(
    label: LinearAnnotationLabel,
    *,
    origin_y: float,
    obstacles: Sequence[LabelObstacle],
    hard: Sequence[Aabb],
    record: _Record,
    gap: float,
    new_rows: bool,
    may_back: bool,
) -> _Way | None:
    """The record's other labels making way for ``label``, which keeps conflicts (LP-D82), or
    ``None`` when they cannot clear every conflict of it (LP-D77).

    Each feature label of the record whose letters ``label`` meets moves to
    another row (``RGN-YIELD-ROW``); then ``label`` may meet only the record's
    leader lines, and is backed (``RGN-FRONT``) where its backing covers
    nothing but those leaders and no other ink. A move to a new row stops the
    search there: the next assembly decides.
    """

    ink = label.ink.translated(0.0, origin_y)
    labels = list(record.labels)
    obstacles = list(obstacles)
    rows: list[FeatureLabelRow] = []
    new_row = False
    for owner, letters in sorted(_hits(ink, obstacles).items()):
        if not letters:
            continue
        if owner[0] != record.record_index or owner[1] != "feature":
            return None
        found = _another_row(owner, ink, record, labels, obstacles, gap=gap, new_rows=new_rows)
        if found is None:
            return None
        row, room, moved, boxes = found
        labels[owner[2]] = moved
        obstacles = [*(item for item in obstacles if item.owner != owner), *boxes]
        rows.append(FeatureLabelRow(record.record_index, cast(int, moved["source_feature_index"]), row))
        new_row = new_row or room
    if new_row or not _hits(ink, obstacles):
        return _Way(label, labels, obstacles, tuple(rows), new_row)
    if not may_back:
        return None

    def behind(obstacle: LabelObstacle) -> bool:
        return obstacle.leader and obstacle.owner[:2] == (record.record_index, "feature")

    # The backing is wider than the letters: it keeps clear of what is not behind it, as the letters do.
    backed = replace(label, backed=True)
    backed = backed.moved_to(_place_along_span(backed, origin_y, [item for item in obstacles if not behind(item)], gap))
    backing = backed.backing.translated(0.0, origin_y)
    met = [obstacle for obstacle in obstacles if obstacle.box.intersects(backing)]
    if not all(behind(obstacle) for obstacle in met) or any(box.intersects(backing) for box in hard):
        return None
    return _Way(backed, labels, obstacles, tuple(rows), False)


def place_linear_annotation_labels(
    *,
    record_index: int,
    record_length: int,
    sequence_width: float,
    plan: LinearRecordVerticalPlan,
    annotation_track_layouts: Mapping[str, ResolvedAnnotationTrack],
    obstacles: Sequence[LabelObstacle],
    feature_labels: Sequence[dict[str, object]] = (),
    glyphs: Sequence[Aabb] = (),
    font_family: str,
    leader_stroke_width: float = 0.0,
    gap: float,
    new_rows: bool = True,
    label_x: float = 0.0,
    row_labels: Sequence[tuple[Sequence[dict[str, object]], float]] = (),
) -> PlacedAnnotationLabels:
    """Each annotation slot's labels on one record, placed along their spans (slot-local y).

    ``obstacles`` is the ink placed before these labels (the record's feature
    labels and leaders, the record names and definition text, the other
    records of its row), from this record's bar start and axis; ``glyphs``
    the record's feature glyph boxes (``linear_feature_glyph_boxes``). The
    record's ``feature_labels`` may make way (LP-D82); ``new_rows`` allows a
    move to a new row. ``label_x`` and ``row_labels`` are where the label fit
    counts the row's feature labels (the record's start and the other records'
    labels, each with its start); no move makes that count larger. A placed
    label still over an obstacle counts one pair per obstacle owner.
    """

    obstacles = list(obstacles)
    record = _Record(
        record_index,
        list(feature_labels),
        _features_origin_y(plan),
        0.5 * float(leader_stroke_width),
        tuple(
            LabelObstacle(box, (record_index, "ink", index))
            for index, box in enumerate(_track_inks(plan, "", glyphs, sequence_width))
        ),
        float(label_x),
        tuple(row_labels),
    )
    placed: dict[str, tuple[LinearAnnotationLabel, ...]] = {}
    font_sizes: list[float] = []
    rows: list[FeatureLabelRow] = []
    proposed: list[FeatureLabelRow] = []
    proposers: list[int] = []
    conflicted: list[int] = []
    overlaps = proposed_overlaps = 0
    for slot in plan.slots:
        track = annotation_track_layouts.get(slot.id)
        if slot.renderer != "annotations" or track is None:
            continue
        params = annotation_track_params_from_mapping(slot.params)
        labels = linear_annotation_labels(
            track,
            record_index=record_index,
            record_length=record_length,
            bar_length_px=sequence_width,
            side=slot.side,
            covers_slots=slot.spans_slots,
            height_px=float(slot.height),
            font_family=font_family,
            params=params,
        )
        slot_labels: list[LinearAnnotationLabel] = []
        origin_y = float(slot.origin_y)
        hard = _track_inks(plan, slot.id, glyphs, sequence_width)
        for label in labels:
            moved = _choose_candidate(label, origin_y=origin_y, obstacles=obstacles, hard=hard, gap=gap)
            conflicts = _conflicts(moved, origin_y, obstacles)
            way = (
                _make_way(
                    moved,
                    origin_y=origin_y,
                    obstacles=obstacles,
                    hard=hard,
                    record=record,
                    gap=gap,
                    new_rows=new_rows,
                    # A clipped track's letters are clipped with it, so they cannot come in front.
                    may_back=not track.clipped,
                )
                if conflicts
                else None
            )
            if way is not None and way.new_row:
                proposed.extend(way.rows)
                proposed_overlaps += conflicts
                proposers.append(len(font_sizes))
            elif way is not None:
                moved, conflicts = way.label, 0
                record = record._replace(labels=way.labels)
                obstacles = way.obstacles
                rows.extend(way.rows)
            overlaps += conflicts
            if conflicts:
                conflicted.append(len(font_sizes))
            obstacles.append(
                LabelObstacle(moved.footprint.translated(0.0, origin_y), (record_index, "region", len(font_sizes)))
            )
            slot_labels.append(moved)
            font_sizes.append(moved.font_size)
        placed[slot.id] = tuple(slot_labels)
    return PlacedAnnotationLabels(
        placed,
        RecordRegionLabels(tuple(font_sizes), overlaps, tuple(conflicted)),
        tuple(rows),
        tuple(proposed),
        proposed_overlaps,
        tuple(proposers),
    )


def linear_annotation_label_inks(
    labels: AnnotationLabelsBySlot, plan: LinearRecordVerticalPlan, record_index: int
) -> list[LabelObstacle]:
    """What a record's placed annotation labels cover, from its bar start and axis."""

    placed = [
        (label, float(plan.slot_by_id(slot_id).origin_y))
        for slot_id, slot_labels in labels.items()
        for label in slot_labels
    ]
    return [
        LabelObstacle(label.footprint.translated(0.0, origin_y), (record_index, "region", index))
        for index, (label, origin_y) in enumerate(placed)
    ]


__all__ = [
    "AnnotationLabelsBySlot",
    "LabelObstacle",
    "ObstacleOwner",
    "PlacedAnnotationLabels",
    "linear_annotation_label_inks",
    "linear_feature_glyph_boxes",
    "linear_feature_label_obstacles",
    "place_linear_annotation_labels",
]
