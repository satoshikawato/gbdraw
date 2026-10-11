"""New drawing from a region of a source drawing.

A region drawing draws the same records of the same files as its source, cut to
regions given in source coordinates (1-based, inclusive). It keeps the source's
look for those records: settings, colors and rules, label tables and filters,
per-feature edits and placements, annotations and Depth inside the regions. A
hand-set size is kept only where its Auto value is the same at the new length
and mode (``gbdraw.auto_sizes``). Comparison tables are not carried: they are
in the coordinates of the uncut records.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, fields, replace
from numbers import Real
from typing import Any, cast

from gbdraw.annotations.models import (
    AnnotationOptions,
    CoordinateSpan,
    FeatureIdentitySpan,
    RegionAnnotation,
)
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import (
    CircularDiagramOptions,
    CircularRequestTrackOptions,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    _ModeDiagramOptions,
    config_override_mode,
)
from gbdraw.api.record_planning import ResolvedRecordCollection, ResolvedRecordProvenance
from gbdraw.api.request_render import resolve_request_records
from gbdraw.api.requests import (
    CircularDiagramRequest,
    DiagramRequest,
    LinearDiagramRequest,
    RecordDisplayOptions,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.auto_sizes import AutoValue, DrawnExtent, SettingAdaptation, adapt_explicit_settings
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import RecordSelector
from gbdraw.io.regions import RegionSpec
from gbdraw.session import MaterializedSession, SessionDrawingSpec, session_to_request
from gbdraw.session_drawings import DrawingMode
from gbdraw.tracks.circular import CircularTrackSlot

_SIZE_OPTIONS = ("window", "step", "depth_window", "depth_step")
_COMPARISON_THRESHOLDS = ("evalue", "bitscore", "identity", "alignment_length")
_SHARED_FIELDS = tuple(option.name for option in fields(_ModeDiagramOptions) if option.init)
_LINEAR_ONLY_FIELDS = tuple(
    option.name for option in fields(LinearDiagramOptions) if option.init and option.name not in _SHARED_FIELDS
)
_CIRCULAR_ONLY_FIELDS = tuple(
    option.name for option in fields(CircularDiagramOptions) if option.init and option.name not in _SHARED_FIELDS
)
_LINEAR_COMPARISON_INPUTS = (
    "blast_files",
    "linear_comparisons",
    "protein_comparisons",
    "orthogroups",
    "collinearity_blocks",
    "comparison_table_file",
)
# Ring results are in the coordinates of the whole record; the plan is ring-aligned.
_RING_RESULTS = (
    "conservation_blast_files",
    "conservation_dataframes",
    "conservation_table_file",
    "conservation_search_results",
)
_RING_PLAN = (
    "conservation_sequence_files",
    "conservation_labels",
    "conservation_colors",
    "conservation_losat_gencodes",
)


@dataclass(frozen=True)
class RegionSelection:
    """A region of one record of the source drawing.

    ``record_key`` is the source drawing's record key; ``start`` and ``end``
    are source coordinates, 1-based and inclusive.
    """

    record_key: str
    start: int
    end: int


@dataclass(frozen=True)
class RegionDrawing:
    """A derived region drawing, the sizes that returned to Auto, and what was not carried."""

    drawing: SessionDrawingSpec
    adaptation: SettingAdaptation
    dropped: tuple[str, ...]


@dataclass(frozen=True)
class _Region:
    item: ResolvedRecordProvenance
    index: int
    start: int
    end: int
    whole: bool

    def contains(self, location_parts: Sequence[tuple[int, int, int | None]]) -> bool:
        """Whether a feature (0-based, end-exclusive parts) intersects the region."""
        return self.whole or any(start < self.end and end >= self.start for start, end, _strand in location_parts)


def derive_region_drawing(
    materialized: MaterializedSession,
    regions: Sequence[RegionSelection],
    *,
    drawing: str | None = None,
    mode: DrawingMode | None = None,
    margin: int = 0,
    reverse_complement: bool | None = None,
    adapt_sizes: bool = True,
    name: str | None = None,
) -> RegionDrawing:
    """Derive a drawing of ``regions`` from a drawing of ``materialized``.

    ``drawing`` selects the source drawing by ID or name, as in
    :func:`gbdraw.session.session_to_request`; without it, the Session's only
    drawing. One region per record; the new drawing lists the records in source order.
    ``margin`` bp are added on each side and clamped to the record ends. The
    mode defaults to Linear, or to the source's mode for one whole record; a
    Circular drawing takes exactly one record. Each record keeps the source's
    orientation unless ``reverse_complement`` is given. With ``adapt_sizes``
    off, every hand-set size is kept. A region that crosses the origin
    (``start > end``) is refused. Use the result while ``materialized`` is
    active: its request reads the materialized files.
    """
    # The one read of the source drawing's settings and edits.
    source = session_to_request(materialized, drawing=drawing)
    collection = resolve_request_records(source)
    selected = _resolve_regions(collection, regions, margin=margin)
    source_mode: DrawingMode = "linear" if isinstance(source, LinearDiagramRequest) else "circular"
    target_mode: DrawingMode = mode or (source_mode if len(selected) == 1 and selected[0].whole else "linear")
    if target_mode not in {"circular", "linear"}:
        raise ValidationError(f"Unknown drawing mode: {target_mode!r}.", diagnostic={"code": "INPUT_INVALID"})
    if target_mode == "circular" and len(selected) != 1:
        raise ValidationError(
            "A Circular region drawing draws exactly one record.", diagnostic={"code": "INPUT_INVALID"}
        )

    options = source.options
    cfg = apply_config_overrides(options.config, None)
    source_extent = DrawnExtent(source_mode, tuple(len(record) for record in collection.records))
    target_extent = DrawnExtent(target_mode, tuple(region.end - region.start + 1 for region in selected))
    adaptation = adapt_explicit_settings(
        _explicit_values(options, source_extent.size_class(cfg)),
        source=source_extent,
        target=target_extent,
        cfg=cfg,
    )
    if not adapt_sizes:
        kept = {**adaptation.kept, **{reset.setting: reset.explicit for reset in adaptation.reset}}
        adaptation = SettingAdaptation(kept=kept, reset=())
    returned = {reset.setting for reset in adaptation.reset}

    dropped: list[str] = []
    config_overrides: dict[str, object] = {}
    other_mode_paths: list[str] = []
    for path, value in (options.config_overrides or {}).items():
        if _setting_of(path) in returned:
            continue
        if config_override_mode(path) not in {None, target_mode}:
            other_mode_paths.append(path)
            continue
        config_overrides[path] = value
    if other_mode_paths:
        dropped.append(f"{source_mode.title()}-only settings: {', '.join(other_mode_paths)}.")
    by_key = {region.item.record_key: region for region in selected}
    values: dict[str, Any] = {name: getattr(options, name) for name in _SHARED_FIELDS}
    values.update(
        {name: None if name in returned else getattr(options, name) for name in _SIZE_OPTIONS},
        config_overrides=config_overrides or None,
        feature_overrides=tuple(
            row for row in options.feature_overrides if _keeps_feature(by_key, row.record_key, row.biological_feature_id)
        ),
        feature_override_table=None,
        feature_override_table_file=None,
        feature_placements=tuple(
            row
            for row in options.feature_placements
            if target_mode == source_mode and _keeps_feature(by_key, row.record_key, row.biological_feature_id)
        ),
        annotations=_region_annotations(options.annotations, collection, selected, dropped),
        depth_tracks=_region_depth_tracks(options, selected, target_mode),
    )
    if target_mode != source_mode:
        # Comparison thresholds filter rings in Circular and ribbons in Linear:
        # the new mode starts from its own defaults.
        values.update(dict.fromkeys(_COMPARISON_THRESHOLDS))
        if options.feature_placements:
            dropped.append(f"{len(options.feature_placements)} feature placement(s) of the {source_mode.title()} drawing.")
    if options.feature_placement_table is not None or options.feature_placement_table_file is not None:
        dropped.append("The feature placement table.")
        values.update(feature_placement_table=None, feature_placement_table_file=None)

    records = tuple(_record_input(source, collection, region, reverse_complement) for region in selected)
    output = getattr(source, "output", None) or RenderOutputRequest()
    request: DiagramRequest
    if target_mode == "linear":
        layout = None
        if isinstance(source, LinearDiagramRequest):
            assert isinstance(options, LinearDiagramOptions)
            values.update(
                {name: getattr(options, name) for name in _LINEAR_ONLY_FIELDS},
                **_linear_comparisons(options, selected, dropped),
            )
            if source.layout is not None:
                layout = LinearMultiRecordOptions(record_gap_px=source.layout.record_gap_px)
            if source.similarity_alignment is not None:
                dropped.append("The similarity alignment of the source records.")
        elif _has_any(options, _RING_RESULTS + _RING_PLAN):
            dropped.append("The conservation rings of the Circular drawing.")
        request = LinearDiagramRequest(
            records=records, options=LinearDiagramOptions(**values), layout=layout, output=output
        )
    else:
        if isinstance(options, CircularDiagramOptions):
            values.update({name: getattr(options, name) for name in _CIRCULAR_ONLY_FIELDS})
            if options.tracks is not None:
                values["tracks"] = _circular_tracks(options.tracks, returned)
            if _has_any(options, _RING_RESULTS):
                values.update(dict.fromkeys(_RING_RESULTS))
                if options.losat_search is None:
                    values.update(dict.fromkeys(_RING_PLAN))
                    dropped.append("The conservation ring tables (they use the coordinates of the whole record).")
        elif _has_any(options, _LINEAR_COMPARISON_INPUTS) or getattr(options, "losat_search", None) is not None:
            dropped.append("The comparisons of the Linear drawing.")
        request = CircularDiagramRequest(
            records=records, options=CircularDiagramOptions(**values), output=output, grouping="single"
        )
    return RegionDrawing(
        drawing=SessionDrawingSpec(request=request, mode=target_mode, name=name),
        adaptation=adaptation,
        dropped=tuple(dropped),
    )


def _resolve_regions(
    collection: ResolvedRecordCollection,
    selections: Sequence[RegionSelection],
    *,
    margin: int,
) -> tuple[_Region, ...]:
    if not selections:
        raise ValidationError(
            "A region drawing needs at least one region.", diagnostic={"code": "INPUT_INVALID", "reason": "REQUIRED"}
        )
    if isinstance(margin, bool) or not isinstance(margin, int) or margin < 0:
        raise ValidationError(
            "margin must be a non-negative integer.",
            diagnostic={"code": "INPUT_INVALID", "reason": "NONNEGATIVE_INTEGER"},
        )
    items = {item.record_key: item for item in collection.provenance}
    resolved: dict[str, tuple[ResolvedRecordProvenance, int, int, bool]] = {}
    for selection in selections:
        item = items.get(selection.record_key)
        if item is None:
            raise ValidationError(
                f"The source drawing has no record {selection.record_key!r}.", diagnostic={"code": "INPUT_INVALID"}
            )
        if item.record_key in resolved:
            raise ValidationError(
                f"Record {selection.record_key!r} has more than one region.", diagnostic={"code": "INPUT_INVALID"}
            )
        length = int(item.source_length or len(collection.records[item.resolved_index]))
        start, end = selection.start, selection.end
        if start > end:
            raise ValidationError(
                f"Region {start}..{end} of {selection.record_key!r} ends before it starts; "
                "regions that cross the origin of a circular record are not supported.",
                diagnostic={"code": "INPUT_INVALID"},
            )
        if start < 1 or end > length:
            raise ValidationError(
                f"Region {start}..{end} of {selection.record_key!r} is outside 1..{length}.",
                diagnostic={"code": "INPUT_INVALID"},
            )
        start, end = max(1, start - margin), min(length, end + margin)
        resolved[item.record_key] = (item, start, end, (start, end) == (1, length))
    ordered = sorted(resolved.values(), key=lambda value: value[0].resolved_index)
    return tuple(_Region(item, index, start, end, whole) for index, (item, start, end, whole) in enumerate(ordered))


def _record_input(
    source: DiagramRequest,
    collection: ResolvedRecordCollection,
    region: _Region,
    reverse_complement: bool | None,
) -> RecordInput:
    item = region.item
    selector = item.selector
    if selector is None and item.source_record_count > 1:
        index = item.source_record_index
        selector = RecordSelector(raw=f"#{index + 1}", record_id=None, record_index=index)
    reverse = (
        item.presentation.reverse_complement or bool(item.region and item.region.reverse_complement)
        if reverse_complement is None
        else reverse_complement
    )
    annotations = collection.records[item.resolved_index].annotations
    label, subtitle = annotations.get("gbdraw_record_label"), annotations.get("gbdraw_record_subtitle")
    crop = None
    if not region.whole:
        crop = RegionSpec(
            raw=f"{region.start}-{region.end}" + (":rc" if reverse else ""),
            file_selector=None,
            record_id=None,
            record_index=None,
            start=region.start,
            end=region.end,
            reverse_complement=reverse,
        )
    return RecordInput(
        source=source.records[item.input_index].source,
        selector=selector,
        region=crop,
        presentation=replace(
            item.presentation,
            label=None if label is None else str(label),
            subtitle=None if subtitle is None else str(subtitle),
            reverse_complement=reverse and region.whole,
            grid_row=None,
            grid_column=None,
        ),
        record_key=item.record_key,
        # A display start cannot be combined with a crop.
        display=item.display if region.whole else RecordDisplayOptions(),
    )


def _setting_of(path: str) -> str:
    head, _, leaf = path.rpartition(".")
    return head if leaf in {"short", "long"} else path


def _explicit_values(options: _ModeDiagramOptions, size_class: str) -> dict[str, AutoValue]:
    """The source's hand-set numbers, keyed as ``auto_setting_values`` is.

    A short/long pair counts once, with the value of the source's size class.
    """
    explicit: dict[str, AutoValue] = {}
    overrides = options.config_overrides or {}
    for path, value in overrides.items():
        if isinstance(value, Real) and not isinstance(value, bool):
            setting = _setting_of(path)
            chosen = overrides.get(f"{setting}.{size_class}", value) if setting != path else value
            explicit[setting] = cast(AutoValue, chosen)
    for option in _SIZE_OPTIONS:
        value = getattr(options, option)
        if value is not None:
            explicit[option] = value
    tracks = getattr(options, "tracks", None)
    if isinstance(tracks, CircularRequestTrackOptions):
        for slot in tracks.circular_track_slots or ():
            if isinstance(slot, CircularTrackSlot) and slot.width is not None and slot.width.unit == "px":
                explicit[f"circularSlot.{slot.id}.width"] = slot.width.value
    return explicit


def _keeps_feature(by_key: Mapping[str, _Region], record_key: str, feature_id: str) -> bool:
    region = by_key.get(record_key)
    if region is None:
        return False
    for feature in region.item.source_feature_catalog or ():
        if feature.biological_feature_id == feature_id:
            return region.contains(feature.location_parts)
    # Without a catalog entry the edit is kept and reports itself as absent.
    return True


def _region_annotations(
    annotations: AnnotationOptions | None,
    collection: ResolvedRecordCollection,
    regions: Sequence[_Region],
    dropped: list[str],
) -> AnnotationOptions | None:
    """Annotation items inside or overlapping the regions, retargeted to the new records."""
    if annotations is None:
        return None
    if annotations.table is not None or annotations.table_file is not None:
        dropped.append("The annotation table.")
        return None
    by_source = {region.item.resolved_index: region for region in regions}
    by_key = {region.item.record_key: region for region in regions}
    source_ids = [str(record.id) for record in collection.records]
    new_ids = [source_ids[region.item.resolved_index] for region in regions]
    outside = 0

    def retarget(annotation: RegionAnnotation) -> RegionAnnotation | None:
        nonlocal outside
        target = annotation.target
        if isinstance(target, FeatureIdentitySpan):
            if _keeps_feature(by_key, target.record_key, target.biological_feature_id):
                return annotation
            outside += 1
            return None
        region = by_source.get(_bound_record(source_ids, target.record))
        if region is None:
            outside += 1
            return None
        if isinstance(target, CoordinateSpan):
            if target.wraps_origin or target.coordinate_space == "local":
                reason = "crosses the origin" if target.wraps_origin else "uses drawn coordinates"
                dropped.append(f"Annotation {annotation.id!r} ({reason}).")
                return None
            inside = region.start <= target.start and target.end <= region.end
            if target.end < region.start or target.start > region.end or (target.out_of_bounds == "error" and not inside):
                outside += 1
                return None
        selector = target.record
        if selector is None or selector.record_id is None or new_ids.count(selector.record_id) != 1:
            selector = (
                None
                if len(regions) == 1
                else RecordSelector(raw=f"#{region.index + 1}", record_id=None, record_index=region.index)
            )
        return replace(annotation, target=replace(target, record=selector))

    sets = tuple(
        replace(
            annotation_set,
            annotations=tuple(
                kept for kept in (retarget(annotation) for annotation in annotation_set.annotations) if kept is not None
            ),
        )
        for annotation_set in annotations.sets
    )
    if outside:
        dropped.append(f"{outside} annotation(s) outside the regions.")
    return AnnotationOptions(sets=sets)


def _bound_record(record_ids: Sequence[str], selector: RecordSelector | None) -> int:
    """The source record an annotation selector names, as the resolver binds it; -1 for none."""
    if selector is None:
        return 0 if len(record_ids) == 1 else -1
    if selector.record_index is not None:
        return selector.record_index
    matches = [index for index, record_id in enumerate(record_ids) if record_id == selector.record_id]
    return matches[0] if len(matches) == 1 else -1


def _region_depth_tracks(
    options: _ModeDiagramOptions,
    regions: Sequence[_Region],
    target_mode: DrawingMode,
) -> tuple[Any, ...] | None:
    """Depth tracks of the drawn records; a per-record source follows its record."""
    tracks = []
    for track in options.depth_tracks or ():
        source = track.source
        if isinstance(source, tuple):
            source = tuple(source[region.item.resolved_index] for region in regions)
            if all(item is None for item in source):
                continue
        tracks.append(replace(track, source=source, height=track.height if target_mode == "linear" else None))
    return tuple(tracks) or None


def _linear_comparisons(
    options: LinearDiagramOptions,
    regions: Sequence[_Region],
    dropped: list[str],
) -> dict[str, Any]:
    """Comparison results are dropped; a LOSAT search runs again on the drawn records."""
    changes: dict[str, Any] = dict.fromkeys(_LINEAR_COMPARISON_INPUTS)
    search = options.losat_search
    if search is None:
        if _has_any(options, _LINEAR_COMPARISON_INPUTS):
            dropped.append("The comparison tables (they use the coordinates of the whole records).")
        return changes
    new_index = {region.item.resolved_index: region.index for region in regions}
    pairs = None
    if search.pairs is not None:
        pairs = [
            (new_index[query], new_index[subject])
            for query, subject in search.pairs
            if query in new_index and subject in new_index
        ]
    gencodes = tuple(search.record_gencodes)
    if len(gencodes) > 1:
        gencodes = tuple(gencodes[region.item.input_index] for region in regions)
    if len(regions) < 2 or pairs == []:
        changes["losat_search"] = None
        dropped.append("The LOSAT comparison (no pair of drawn records to compare).")
    else:
        changes["losat_search"] = replace(search, pairs=pairs, record_gencodes=gencodes)
    return changes


def _circular_tracks(tracks: CircularRequestTrackOptions, returned: set[str]) -> CircularRequestTrackOptions:
    """Track widths whose Auto value changed return to Auto."""
    if tracks.circular_track_slots is None:
        return tracks
    return replace(
        tracks,
        circular_track_slots=[
            replace(slot, width=None)
            if isinstance(slot, CircularTrackSlot) and f"circularSlot.{slot.id}.width" in returned
            else slot
            for slot in tracks.circular_track_slots
        ],
    )


def _has_any(options: object, names: Sequence[str]) -> bool:
    """Whether any of the option fields holds a value (an empty sequence holds none)."""
    for name in names:
        value = getattr(options, name, None)
        if value is not None and not (isinstance(value, (tuple, list)) and not value):
            return True
    return False


__all__ = ["RegionDrawing", "RegionSelection", "derive_region_drawing"]
