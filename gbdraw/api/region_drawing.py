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

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, fields, replace
from numbers import Real
from typing import Any, Literal, cast

from gbdraw.analysis.origin_arcs import location_origin_arc
from gbdraw.annotations.models import (
    AnnotationOptions,
    CoordinateSpan,
    FeatureIdentitySpan,
    RegionAnnotation,
)
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import (
    DEPTH_COMPATIBILITY_FIELDS,
    CircularDiagramOptions,
    CircularRequestTrackOptions,
    DepthTrackInput,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    _ModeDiagramOptions,
    config_override_mode,
)
from gbdraw.api.record_planning import ResolvedRecordCollection, ResolvedRecordProvenance
from gbdraw.api.request_render import resolve_request_records
from gbdraw.api.requests import (
    CircularBatchRequest,
    CircularDiagramRequest,
    DiagramRequest,
    LinearDiagramRequest,
    RecordDisplayOptions,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.analysis.depth_tracks import DepthTrackStyle, depth_track_config
from gbdraw.auto_sizes import (
    ADAPTABLE_SETTINGS,
    AutoValue,
    DrawnDepthRow,
    DrawnExtent,
    SettingAdaptation,
    SIZE_CLASSES,
    SizeClass,
    adapt_explicit_settings,
)
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.models.render_profiles import CircularRenderProfile
from gbdraw.configurators.depth import DepthConfigurator
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import RecordSelector
from gbdraw.io.regions import RecordWindow, window_drawn_reverse
from gbdraw.session import MaterializedSession, SessionDrawingSpec, session_to_request
from gbdraw.session_drawings import DrawingMode
from gbdraw.session_request_codec import canonical_depth_tracks
from gbdraw.tracks.circular import CircularTrackSlot, default_circular_track_slots, parse_circular_track_slots
from gbdraw.tracks.parsing import depth_slot_track_index, feature_density_slot_track_id

_SIZE_OPTIONS = ("window", "step", "depth_window", "depth_step")
_SIZE_SETTINGS = frozenset(setting for settings in ADAPTABLE_SETTINGS.values() for setting in settings)
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
# Ring results are in the search frame of the source drawing's records; the plan is ring-aligned.
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
class FeatureRegionSelection:
    """The span of one feature of a record of the source drawing.

    ``feature_id`` is the feature's biological feature ID in the source
    drawing's feature catalog; the region runs from the first to the last base
    of its location, in source coordinates. A feature across the origin is
    refused as an origin-crossing region.
    """

    record_key: str
    feature_id: str


RegionChoice = RegionSelection | FeatureRegionSelection


NotCarriedKind = Literal[
    "OTHER_MODE_SETTINGS",
    "FEATURE_PLACEMENTS",
    "FEATURE_PLACEMENT_TABLE",
    "SIMILARITY_ALIGNMENT",
    "CONSERVATION_RINGS",
    "CONSERVATION_RING_TABLES",
    "LINEAR_COMPARISONS",
    "COMPARISON_TABLES",
    "LOSAT_NO_PAIR",
    "ANNOTATION_TABLE",
    "ANNOTATION_CROSSES_ORIGIN",
    "ANNOTATION_DRAWN_COORDINATES",
    "ANNOTATIONS_OUTSIDE",
]


@dataclass(frozen=True)
class NotCarried:
    """One part of the source drawing that the region drawing does not carry.

    ``kind`` names the part; ``count`` and ``items`` (setting paths, an
    annotation ID) are its locators. The CLI and API text is ``message``; the
    Web words the kind in its own wording owner.
    """

    kind: NotCarriedKind
    source_mode: DrawingMode
    count: int = 0
    items: tuple[str, ...] = ()

    @property
    def message(self) -> str:
        mode = self.source_mode.title()
        return {
            "OTHER_MODE_SETTINGS": lambda: f"{mode}-only settings: {', '.join(self.items)}.",
            "FEATURE_PLACEMENTS": lambda: f"{self.count} feature placement(s) of the {mode} drawing.",
            "FEATURE_PLACEMENT_TABLE": lambda: "The feature placement table.",
            "SIMILARITY_ALIGNMENT": lambda: "The similarity alignment of the source records.",
            "CONSERVATION_RINGS": lambda: "The conservation rings of the Circular drawing.",
            "CONSERVATION_RING_TABLES": lambda: (
                "The conservation ring tables (they use the coordinates of the records of the source drawing)."
            ),
            "LINEAR_COMPARISONS": lambda: "The comparisons of the Linear drawing.",
            "COMPARISON_TABLES": lambda: "The comparison tables (they use the coordinates of the whole records).",
            "LOSAT_NO_PAIR": lambda: "The LOSAT comparison (no pair of drawn records to compare).",
            "ANNOTATION_TABLE": lambda: "The annotation table.",
            "ANNOTATION_CROSSES_ORIGIN": lambda: f"Annotation {self.items[0]!r} (crosses the origin).",
            "ANNOTATION_DRAWN_COORDINATES": lambda: f"Annotation {self.items[0]!r} (uses drawn coordinates).",
            "ANNOTATIONS_OUTSIDE": lambda: f"{self.count} annotation(s) outside the regions.",
        }[self.kind]()


StartAfterEndReason = Literal["START_AFTER_END"]


class StartAfterEndRegion(ValidationError):
    """A region of a linear record whose start is after its end: a linear record has no origin to cross.

    On a circular record such a region runs across the origin and is drawn as
    one window (D-28, D-34).
    """

    def __init__(self, record_key: str, start: int, end: int) -> None:
        reason: StartAfterEndReason = "START_AFTER_END"
        super().__init__(
            f"Region {start}..{end} of {record_key!r} ends before it starts; choose a start before the end.",
            diagnostic={"code": "REGION_INVALID", "field": "region", "reason": reason},
        )
        self.record_key = record_key
        self.start = start
        self.end = end
        self.reason = reason


@dataclass(frozen=True)
class RegionRequest:
    """A derived region request, its mode, the sizes that returned to Auto, and what was not carried."""

    request: DiagramRequest
    mode: DrawingMode
    adaptation: SettingAdaptation
    not_carried: tuple[NotCarried, ...]


@dataclass(frozen=True)
class RegionDrawing:
    """A derived region drawing, the sizes that returned to Auto, and what was not carried."""

    drawing: SessionDrawingSpec
    adaptation: SettingAdaptation
    not_carried: tuple[NotCarried, ...]

    @property
    def dropped(self) -> tuple[str, ...]:
        """What was not carried, as the CLI and API word it."""
        return tuple(item.message for item in self.not_carried)


@dataclass(frozen=True)
class SourceRegion:
    """The spans of one source record a derived drawing draws: one window each, in source order (D-01).

    ``index`` is the record's place in the new drawing; ``first_unit`` the drawn index of its first window there.
    A span whose start is after its end runs across the origin of the circular record of ``length`` bases.
    """

    item: ResolvedRecordProvenance
    index: int
    spans: tuple[tuple[int, int], ...]
    whole: bool
    first_unit: int
    length: int

    @property
    def pieces(self) -> tuple[tuple[int, int], ...]:
        """The spans as runs of ascending source bases: a span across the origin is two."""
        return tuple(piece for span in self.spans for piece in _span_pieces(span, self.length))

    def contains(self, location_parts: Sequence[tuple[int, int, int | None]]) -> bool:
        """Whether a feature (0-based, end-exclusive parts) intersects the region."""
        return self.whole or any(
            start < span_end and end >= span_start
            for span_start, span_end in self.pieces
            for start, end, _strand in location_parts
        )

    @property
    def drawn_lengths(self) -> tuple[int, ...]:
        return tuple(sum(end - start + 1 for start, end in _span_pieces(span, self.length)) for span in self.spans)

    @property
    def units(self) -> range:
        """The drawn indexes of the record's windows in the new drawing."""
        return range(self.first_unit, self.first_unit + len(self.spans))


def _span_pieces(span: tuple[int, int], length: int) -> tuple[tuple[int, int], ...]:
    """A 1-based inclusive span as ascending runs: one, or two for a span across the origin (start after end)."""
    start, end = span
    return ((start, length), (1, end)) if start > end else (span,)


def derive_region_drawing(
    materialized: MaterializedSession,
    regions: Sequence[RegionChoice],
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
    off, every hand-set size is kept. A region whose start is after its end
    runs across the origin of a circular record and is drawn as one window
    (D-34); on a linear record it is refused. Use the result while
    ``materialized`` is active: its request reads the materialized files.
    """
    # The one read of the source drawing's settings and edits.
    derived = derive_region_request(
        session_to_request(materialized, drawing=drawing),
        regions,
        mode=mode,
        margin=margin,
        reverse_complement=reverse_complement,
        adapt_sizes=adapt_sizes,
    )
    return RegionDrawing(
        drawing=SessionDrawingSpec(request=derived.request, mode=derived.mode, name=name),
        adaptation=derived.adaptation,
        not_carried=derived.not_carried,
    )


def derive_region_request(
    source: DiagramRequest,
    regions: Sequence[RegionChoice],
    *,
    mode: DrawingMode | None = None,
    margin: int = 0,
    reverse_complement: bool | None = None,
    adapt_sizes: bool = True,
) -> RegionRequest:
    """Derive the request of a drawing of ``regions`` from the request of its source drawing.

    The one derivation behind :func:`derive_region_drawing` and the Web's New
    drawing from a region; the arguments are as there.
    """
    collection = resolve_request_records(source)
    selected = _resolve_regions(collection, regions, margin=margin)
    source_mode: DrawingMode = "linear" if isinstance(source, LinearDiagramRequest) else "circular"
    target_mode: DrawingMode = mode or (source_mode if len(selected) == 1 and selected[0].whole else "linear")
    if target_mode not in {"circular", "linear"}:
        raise ValidationError(f"Unknown drawing mode: {target_mode!r}.", diagnostic={"code": "INPUT_INVALID"})
    if target_mode == "circular" and (len(selected) != 1 or len(selected[0].spans) != 1):
        raise ValidationError(
            "A Circular region drawing draws exactly one record.", diagnostic={"code": "INPUT_INVALID"}
        )
    source_ids = [str(record.id) for record in collection.records]
    derived = DerivedRecords(
        mode=target_mode,
        records=tuple(_record_input(source, region, reverse_complement) for region in selected),
        record_ids=tuple(source_ids[region.item.resolved_index] for region in selected),
        lengths=tuple(length for region in selected for length in region.drawn_lengths),
        regions=selected,
    )
    return derive_from_source(
        source,
        collection,
        derived,
        adapt_sizes=adapt_sizes,
        linear_comparisons=lambda options, dropped: _linear_comparisons(options, selected, dropped, source_mode),
    )


@dataclass(frozen=True)
class DerivedRecords:
    """The records of a drawing derived from a source drawing.

    ``records`` are its record inputs and ``record_ids`` their record IDs, in
    drawing order; ``lengths`` the drawn length of each of its units;
    ``regions`` the records that come from records of the source drawing, each
    placed by its ``index`` in ``records``. A record that is not in the source
    drawing carries none of its edits. ``layout`` places records that share a
    Linear row; without it each record takes a row of its own.
    """

    mode: DrawingMode
    records: tuple[RecordInput, ...]
    record_ids: tuple[str, ...]
    lengths: tuple[int, ...]
    regions: tuple[SourceRegion, ...]
    layout: LinearMultiRecordOptions | None = None


# The comparison fields of a derived Linear drawing, from the source's options; it adds what it does not carry.
LinearComparisonRule = Callable[[LinearDiagramOptions, list[NotCarried]], Mapping[str, object]]


def derive_from_source(
    source: DiagramRequest,
    collection: ResolvedRecordCollection,
    derived: DerivedRecords,
    *,
    adapt_sizes: bool,
    linear_comparisons: LinearComparisonRule,
) -> RegionRequest:
    """The request of ``derived`` with the settings and edits of ``source`` that apply to it.

    The one reader of a source drawing's settings for a derived drawing (New
    drawing from a region, Find in other files): a hand-set size whose Auto
    value changes returns to Auto; per-feature edits, annotations and Depth
    follow their records; what is not carried is listed.
    """
    source_mode: DrawingMode = "linear" if isinstance(source, LinearDiagramRequest) else "circular"
    target_mode = derived.mode
    selected = derived.regions
    options = source.options
    by_index = {region.index: region for region in selected}
    source_tracks = canonical_depth_tracks(options, record_count=len(collection.records))
    depth_tracks = _region_depth_tracks(
        source_tracks, [by_index.get(index) for index in range(len(derived.records))], target_mode
    )
    # The source's settings without its sizes, as its drawing reads them; the depth rows each drawing draws.
    cfg = apply_config_overrides(
        options.config,
        {path: value for path, value in (options.config_overrides or {}).items() if _setting_of(path) not in _SIZE_SETTINGS},
    )
    source_extent = DrawnExtent(
        source_mode,
        tuple(len(record) for record in collection.records),
        depth_rows=_drawn_depth_rows(options, source_tracks, cfg, source_mode),
        density_window_settings=_drawn_density_window_settings(options, source_mode),
        separate_records=isinstance(source, CircularBatchRequest),
    )
    target_extent = DrawnExtent(
        target_mode,
        derived.lengths,
        depth_rows=_drawn_depth_rows(options, depth_tracks, cfg, target_mode),
        density_window_settings=_drawn_density_window_settings(options, target_mode),
    )
    # The source record each drawn unit comes from, by its index there: each record compares its own Auto windows.
    record_sources: list[int | None] = [None] * len(derived.lengths)
    for region in selected:
        for unit in region.units:
            record_sources[unit] = region.item.resolved_index
    adaptation = adapt_explicit_settings(
        _explicit_values(options),
        source=source_extent,
        target=target_extent,
        cfg=cfg,
        record_sources=record_sources,
    )
    if not adapt_sizes:
        kept = {**adaptation.kept, **{reset.setting: reset.explicit for reset in adaptation.reset}}
        adaptation = SettingAdaptation(kept=kept, reset=())
    returned = {reset.setting for reset in adaptation.reset}

    dropped: list[NotCarried] = []
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
        dropped.append(NotCarried("OTHER_MODE_SETTINGS", source_mode, items=tuple(other_mode_paths)))
    by_key = {region.item.record_identity_key: region for region in selected}
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
        annotations=_region_annotations(
            options.annotations, collection, selected, dropped, source_mode, derived.record_ids
        ),
        # The depth tracks in their canonical form, which replaces the compatibility inputs.
        **dict.fromkeys(name for name in DEPTH_COMPATIBILITY_FIELDS if name in _SHARED_FIELDS),
        depth_tracks=depth_tracks,
    )
    if target_mode != source_mode:
        # Comparison thresholds filter rings in Circular and ribbons in Linear:
        # the new mode starts from its own defaults.
        values.update(dict.fromkeys(_COMPARISON_THRESHOLDS))
        if options.feature_placements:
            dropped.append(NotCarried("FEATURE_PLACEMENTS", source_mode, count=len(options.feature_placements)))
    if options.feature_placement_table is not None or options.feature_placement_table_file is not None:
        dropped.append(NotCarried("FEATURE_PLACEMENT_TABLE", source_mode))
        values.update(feature_placement_table=None, feature_placement_table_file=None)

    records = derived.records
    output = getattr(source, "output", None) or RenderOutputRequest()
    request: DiagramRequest
    if target_mode == "linear":
        # Each derived record takes a row of its own (no grid row), so the
        # record gap, which spaces the records of one row, draws nothing: the
        # drawing states no layout, as the Web writes it.
        if isinstance(source, LinearDiagramRequest):
            assert isinstance(options, LinearDiagramOptions)
            values.update(
                {name: getattr(options, name) for name in _LINEAR_ONLY_FIELDS if name not in DEPTH_COMPATIBILITY_FIELDS},
                **linear_comparisons(options, dropped),
            )
            if source.similarity_alignment is not None:
                dropped.append(NotCarried("SIMILARITY_ALIGNMENT", source_mode))
        elif _has_any(options, _RING_RESULTS + _RING_PLAN):
            dropped.append(NotCarried("CONSERVATION_RINGS", source_mode))
        request = LinearDiagramRequest(
            records=records, options=LinearDiagramOptions(**values), output=output, layout=derived.layout
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
                    dropped.append(NotCarried("CONSERVATION_RING_TABLES", source_mode))
        elif _has_any(options, _LINEAR_COMPARISON_INPUTS) or getattr(options, "losat_search", None) is not None:
            dropped.append(NotCarried("LINEAR_COMPARISONS", source_mode))
        request = CircularDiagramRequest(
            records=records, options=CircularDiagramOptions(**values), output=output, grouping="single"
        )
    return RegionRequest(request=request, mode=target_mode, adaptation=adaptation, not_carried=tuple(dropped))


def _resolve_regions(
    collection: ResolvedRecordCollection,
    selections: Sequence[RegionChoice],
    *,
    margin: int,
) -> tuple[SourceRegion, ...]:
    if not selections:
        raise ValidationError(
            "A region drawing needs at least one region.", diagnostic={"code": "INPUT_INVALID", "reason": "REQUIRED"}
        )
    if isinstance(margin, bool) or not isinstance(margin, int) or margin < 0:
        raise ValidationError(
            "margin must be a non-negative integer.",
            diagnostic={"code": "INPUT_INVALID", "reason": "NONNEGATIVE_INTEGER"},
        )
    # A record drawn in windows is selected by its own key; its first window stands for it.
    items: dict[str, ResolvedRecordProvenance] = {}
    for unit in collection.provenance:
        items.setdefault(unit.record_identity_key, unit)
    resolved: dict[str, tuple[ResolvedRecordProvenance, list[tuple[int, int]], int]] = {}
    for selection in selections:
        item = items.get(selection.record_key)
        if item is None:
            raise ValidationError(
                f"The source drawing has no record {selection.record_key!r}.", diagnostic={"code": "INPUT_INVALID"}
            )
        length = int(item.source_length or len(collection.records[item.resolved_index]))
        start, end = _selected_span(selection, item, length)
        if start > end and not item.is_circular:
            raise StartAfterEndRegion(selection.record_key, start, end)
        if min(start, end) < 1 or max(start, end) > length:
            raise ValidationError(
                f"Region {start}..{end} of {selection.record_key!r} is outside 1..{length}.",
                diagnostic={"code": "REGION_INVALID", "field": "region", "reason": "RECORD_BOUNDS"},
            )
        resolved.setdefault(selection.record_key, (item, [], length))[1].append(
            _with_margin(start, end, margin, length)
        )
    ordered = sorted(resolved.values(), key=lambda value: value[0].resolved_index)
    regions: list[SourceRegion] = []
    first_unit = 0
    for index, (item, ranges, length) in enumerate(ordered):
        spans = _merged_spans(ranges, length)
        regions.append(SourceRegion(item, index, spans, spans == ((1, length),), first_unit, length))
        first_unit += len(spans)
    return tuple(regions)


def _with_margin(start: int, end: int, margin: int, length: int) -> tuple[int, int]:
    """A span grown by ``margin`` bp on each side: clamped to the record ends, or, across the origin of a circular
    record, grown round the circle up to one whole turn from its new start."""
    if start <= end:
        return max(1, start - margin), min(length, end + margin)
    if (length - start + 1) + end + 2 * margin < length:
        return start - margin, end + margin
    first = (start - margin - 1) % length + 1
    return (1, length) if first == 1 else (first, first - 1)


def _merged_spans(ranges: Sequence[tuple[int, int]], length: int) -> tuple[tuple[int, int], ...]:
    """Regions of one record in source order, overlapping or touching ones joined: each base is drawn once.

    When a region runs across the origin (start after end, D-34), the runs that
    meet at the origin are joined into that one window across it.
    """
    merged: list[tuple[int, int]] = []
    for start, end in sorted(piece for span in ranges for piece in _span_pieces(span, length)):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    crosses = any(start > end for start, end in ranges)
    if crosses and len(merged) > 1 and merged[0][0] == 1 and merged[-1][1] == length:
        merged = [(merged[-1][0], merged[0][1]), *merged[1:-1]]
    return tuple(merged)


def _selected_span(selection: RegionChoice, item: ResolvedRecordProvenance, length: int) -> tuple[int, int]:
    """A selection's 1-based inclusive span; an origin-crossing feature spans its arc (start after end)."""
    if isinstance(selection, RegionSelection):
        return selection.start, selection.end
    for feature in item.source_feature_catalog:
        if feature.biological_feature_id == selection.feature_id:
            parts = feature.location_parts
            arc = location_origin_arc(parts, length)
            if arc is not None:
                return arc[0] + 1, arc[1]
            return min(part[0] for part in parts) + 1, max(part[1] for part in parts)
    raise ValidationError(
        f"Record {selection.record_key!r} has no feature {selection.feature_id!r}.",
        diagnostic={"code": "INPUT_INVALID"},
    )


def _record_input(
    source: DiagramRequest,
    region: SourceRegion,
    reverse_complement: bool | None,
) -> RecordInput:
    item = region.item
    selector = item.selector
    if selector is None and item.source_record_count > 1:
        index = item.source_record_index
        selector = RecordSelector(raw=f"#{index + 1}", record_id=None, record_index=index)
    reverse = (
        window_drawn_reverse(bool(item.region and item.region.reverse_complement), item.presentation.reverse_complement)
        if reverse_complement is None
        else reverse_complement
    )
    # The windows of one record share its row in source order (D-01; the planner places them).
    windows = () if region.whole else tuple(RecordWindow(start, end, reverse) for start, end in region.spans)
    return RecordInput(
        source=source.records[item.input_index].source,
        selector=selector,
        windows=windows,
        presentation=replace(
            item.presentation,
            reverse_complement=reverse and region.whole,
            grid_row=None,
            grid_column=None,
        ),
        record_key=item.record_identity_key,
        # A display start cannot be combined with a crop.
        # A window across the origin keeps the record circular, as the source drew it (D-34).
        display=(
            item.display
            if region.whole
            else RecordDisplayOptions(is_circular=True if any(start > end for start, end in region.spans) else None)
        ),
    )


def _setting_of(path: str) -> str:
    head, _, leaf = path.rpartition(".")
    return head if leaf in {"short", "long"} else path


def _explicit_values(options: _ModeDiagramOptions) -> dict[str, dict[SizeClass, AutoValue]]:
    """The source's hand-set numbers, keyed as ``auto_setting_values`` is, with the value a diagram of each size
    class draws.

    A short/long pair counts once and gives a class only its own leaf: a class without one draws its Auto value.
    """
    explicit: dict[str, dict[SizeClass, AutoValue]] = {}
    for path, value in (options.config_overrides or {}).items():
        if isinstance(value, Real) and not isinstance(value, bool):
            setting = _setting_of(path)
            classes = SIZE_CLASSES if setting == path else tuple(size for size in SIZE_CLASSES if path.endswith(f".{size}"))
            explicit.setdefault(setting, {}).update(dict.fromkeys(classes, cast(AutoValue, value)))
    for option in _SIZE_OPTIONS:
        value = getattr(options, option)
        if value is not None:
            explicit[option] = dict.fromkeys(SIZE_CLASSES, value)
    for slot in _stated_slots(getattr(options, "tracks", None)) or ():
        if slot.width is not None and slot.width.unit == "px":
            explicit[f"circularSlot.{slot.id}.width"] = dict.fromkeys(SIZE_CLASSES, slot.width.value)
    return explicit


def _stated_slots(tracks: object) -> list[CircularTrackSlot] | None:
    """The Circular track slots a request states, as the drawing reads them: a slot given as a string (the CLI
    spelling) is the slot parsed from it. ``None`` where the request states none."""
    slots = tracks.circular_track_slots if isinstance(tracks, CircularRequestTrackOptions) else None
    return None if slots is None else parse_circular_track_slots(slots)


def _drawn_depth_rows(
    options: _ModeDiagramOptions,
    depth_tracks: Sequence[DepthTrackInput] | None,
    cfg: GbdrawConfig,
    mode: DrawingMode,
) -> tuple[DrawnDepthRow, ...]:
    """The depth rows a Circular drawing of ``options`` with ``depth_tracks`` draws, in stack order: its stack's
    enabled depth slots, else the default stack's; each with the axis of the series it draws (that series' own
    tick font over the Depth settings of ``cfg``, as the drawing configures it) and its width as stated."""
    if mode != "circular" or not depth_tracks:
        return ()
    stated = _stated_slots(getattr(options, "tracks", None))
    slots = (
        stated
        if stated is not None
        else default_circular_track_slots(show_depth=True, depth_track_count=len(depth_tracks))
    )
    base = DepthConfigurator(1, 1, CircularRenderProfile(cfg))
    axes = {
        index: depth_track_config(DepthTrackStyle(None, None, None, track.tick_font_size), base, 1, 1, None)
        for index, track in enumerate(depth_tracks)
    }
    rows: list[DrawnDepthRow] = []
    for slot in slots:
        track_index = depth_slot_track_index(slot.id, slot.params) if slot.renderer == "depth" else -1
        axis = axes.get(track_index) if slot.enabled else None
        if axis is None:
            continue
        # The Depth row's hand-set width is the size setting circularSlot.depth.width (_explicit_values).
        depth_width = slot.id == "depth" and (slot.width is None or slot.width.unit == "px")
        width = None if slot.width is None or depth_width else slot.width.resolve(cfg.canvas.circular.radius)
        rows.append(DrawnDepthRow(axis, track_index, width, depth_width))
    return tuple(rows)


def _drawn_density_window_settings(options: _ModeDiagramOptions, mode: DrawingMode) -> tuple[str, ...]:
    """The window settings (``window``, ``step``) that a feature density track a drawing of ``options`` in ``mode``
    draws takes, being one that sets none of its own.

    A stack the request does not state draws every track; a stated Circular stack, which a new Circular drawing
    carries from a Circular source, only those its enabled ``feature_density`` slots name. A stated Linear stack is
    read as drawing every track: it is in a Linear source, whose new drawing is Linear or one Circular record, and
    in each of these drawings the density window is the GC content window (``feature_density_windows``,
    ``record_sliding_windows``), so the density window compares as the window does.
    """
    tracks = options.feature_density_tracks or ()
    stated = _stated_slots(getattr(options, "tracks", None)) if mode == "circular" else None
    if stated is not None:
        named = {
            feature_density_slot_track_id(slot.id, slot.params)
            for slot in stated
            if slot.enabled and slot.renderer == "feature_density"
        }
        tracks = tuple(track for track in tracks if track.id in named)
    return tuple(setting for setting in ("window", "step") if any(getattr(track, setting) is None for track in tracks))


def _keeps_feature(by_key: Mapping[str, SourceRegion], record_key: str, feature_id: str) -> bool:
    region = by_key.get(record_key)
    if region is None:
        return False
    for feature in region.item.source_feature_catalog:
        if feature.biological_feature_id == feature_id:
            return region.contains(feature.location_parts)
    # Without a catalog entry the edit is kept and reports itself as absent.
    return True


def _region_annotations(
    annotations: AnnotationOptions | None,
    collection: ResolvedRecordCollection,
    regions: Sequence[SourceRegion],
    dropped: list[NotCarried],
    source_mode: DrawingMode,
    new_ids: Sequence[str],
) -> AnnotationOptions | None:
    """Annotation items inside or overlapping the regions, retargeted to the new records (``new_ids`` in order)."""
    if annotations is None:
        return None
    if annotations.table is not None or annotations.table_file is not None:
        dropped.append(NotCarried("ANNOTATION_TABLE", source_mode))
        return None
    by_source = {region.item.resolved_index: region for region in regions}
    by_key = {region.item.record_identity_key: region for region in regions}
    source_ids = [str(record.id) for record in collection.records]
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
                kind: NotCarriedKind = (
                    "ANNOTATION_CROSSES_ORIGIN" if target.wraps_origin else "ANNOTATION_DRAWN_COORDINATES"
                )
                dropped.append(NotCarried(kind, source_mode, items=(annotation.id,)))
                return None
            inside = any(start <= target.start and target.end <= end for start, end in region.pieces)
            overlaps = any(target.end >= start and target.start <= end for start, end in region.pieces)
            if not overlaps or (target.out_of_bounds == "error" and not inside):
                outside += 1
                return None
        selector = target.record
        if selector is None or selector.record_id is None or new_ids.count(selector.record_id) != 1:
            selector = (
                None
                if len(new_ids) == 1
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
        dropped.append(NotCarried("ANNOTATIONS_OUTSIDE", source_mode, count=outside))
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
    source_tracks: Sequence[DepthTrackInput] | None,
    slots: Sequence[SourceRegion | None],
    target_mode: DrawingMode,
) -> tuple[DepthTrackInput, ...] | None:
    """Depth tracks of the drawn records, one slot per new record; a per-record source follows its record, and a
    record not in the source drawing has none."""
    tracks = []
    for track in source_tracks or ():
        source = track.source
        if isinstance(source, tuple):
            source = tuple(None if slot is None else source[slot.item.resolved_index] for slot in slots)
            if all(item is None for item in source):
                continue
        tracks.append(replace(track, source=source, height=track.height if target_mode == "linear" else None))
    return tuple(tracks) or None


def _linear_comparisons(
    options: LinearDiagramOptions,
    regions: Sequence[SourceRegion],
    dropped: list[NotCarried],
    source_mode: DrawingMode,
) -> dict[str, Any]:
    """Comparison results are dropped; a LOSAT search runs again on the drawn records."""
    changes: dict[str, Any] = dict.fromkeys(_LINEAR_COMPARISON_INPUTS)
    search = options.losat_search
    if search is None:
        if _has_any(options, _LINEAR_COMPARISON_INPUTS):
            dropped.append(NotCarried("COMPARISON_TABLES", source_mode))
        return changes
    new_units = {region.item.resolved_index: region.units for region in regions}
    pairs = None
    if search.pairs is not None:
        pairs = [
            (new_query, new_subject)
            for query, subject in search.pairs
            if query in new_units and subject in new_units
            for new_query in new_units[query]
            for new_subject in new_units[subject]
        ]
    gencodes = tuple(search.record_gencodes)
    if len(gencodes) > 1:
        gencodes = tuple(gencodes[region.item.input_index] for region in regions)
    if sum(len(region.spans) for region in regions) < 2 or pairs == []:
        changes["losat_search"] = None
        dropped.append(NotCarried("LOSAT_NO_PAIR", source_mode))
    else:
        changes["losat_search"] = replace(search, pairs=pairs, record_gencodes=gencodes)
    return changes


def _circular_tracks(tracks: CircularRequestTrackOptions, returned: set[str]) -> CircularRequestTrackOptions:
    """Track widths whose Auto value changed return to Auto."""
    slots = _stated_slots(tracks)
    if slots is None:
        return tracks
    return replace(
        tracks,
        circular_track_slots=[
            replace(slot, width=None) if f"circularSlot.{slot.id}.width" in returned else slot for slot in slots
        ],
    )


def _has_any(options: object, names: Sequence[str]) -> bool:
    """Whether any of the option fields holds a value (an empty sequence holds none)."""
    for name in names:
        value = getattr(options, name, None)
        if value is not None and not (isinstance(value, (tuple, list)) and not value):
            return True
    return False


__all__ = [
    "DerivedRecords",
    "FeatureRegionSelection",
    "LinearComparisonRule",
    "NotCarried",
    "NotCarriedKind",
    "StartAfterEndReason",
    "StartAfterEndRegion",
    "RegionChoice",
    "RegionDrawing",
    "RegionRequest",
    "RegionSelection",
    "SourceRegion",
    "derive_from_source",
    "derive_region_drawing",
    "derive_region_request",
]
