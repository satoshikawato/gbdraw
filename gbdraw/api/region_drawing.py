"""Copy records or regions of a source drawing to a new drawing or into an existing one.

A region drawing draws the same records of the same files as its source, cut to
regions given in source coordinates (1-based, inclusive). It keeps the source's
look for those records: settings, colors and rules, label tables and filters,
per-feature edits and placements, annotations and Depth inside the regions. A
hand-set size is kept only where its Auto value is the same at the new length
and mode (``gbdraw.auto_sizes``). Comparison tables are not carried: they are
in the coordinates of the uncut records. With ``settings="defaults"`` the new
drawing draws the records with the mode's defaults instead.

A paste (:func:`paste_region_request`) appends the records or regions to an
existing drawing, which keeps its own settings; the records' per-feature edits
come along unless the paste takes the records only.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, fields, replace
from numbers import Real
from types import SimpleNamespace
from typing import Any, Literal, cast

import pandas as pd
from pandas import DataFrame

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
    ColorOptions,
    CircularMultiRecordOptions,
    CircularRequestTrackOptions,
    DepthTrackInput,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    LinearRecordTranslation,
    LosatSearchOptions,
    _ModeDiagramOptions,
    config_override_mode,
)
from gbdraw.api.record_planning import (
    ResolvedRecordCollection,
    ResolvedRecordProvenance,
    SourceFileContents,
    resolve_implicit_record_output_prefix,
)
from gbdraw.api.request_render import resolve_request_records
from gbdraw.api.requests import (
    CircularBatchRequest,
    CircularDiagramRequest,
    EXPANDED_RECORD_KEY_SEPARATOR,
    DiagramRequest,
    LinearDiagramRequest,
    RecordDisplayOptions,
    RecordInput,
    RecordPresentation,
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
from gbdraw.core.record_metadata import RecordSourceKind
from gbdraw.config.models.render_profiles import CircularRenderProfile
from gbdraw.configurators.depth import DepthConfigurator
from gbdraw.exceptions import ValidationError
from gbdraw.features.colors import build_specific_color_rules
from gbdraw.features.ids import exact_value_pattern, feature_instance
from gbdraw.features.overrides import FeatureOverride
from gbdraw.features.selector_values import find_specific_color_rule
from gbdraw.io.colors import read_color_table
from gbdraw.io.record_select import RecordSelector
from gbdraw.io.regions import RecordWindow, window_drawn_reverse
from gbdraw.layout.record_placement import parse_record_row_position
from gbdraw.layout.similarity_alignment import (
    AlignmentDecisionStatus,
    AlignmentRecordDecision,
    AlignmentResolutionRationale,
    SimilarityAlignmentPlan,
    SimilarityAlignmentReference,
)
from gbdraw.session import (
    MaterializedSession,
    SessionDocument,
    SessionDrawingSpec,
    replace_drawing_request,
    session_to_request,
)
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


# Settings of a new drawing: the source drawing's settings and edits, or the mode's defaults (PD-OI-090 (1)).
SettingsChoice = Literal["keep", "defaults"]
# What comes along into an existing drawing: the records with their per-feature edits, or the records only.
PasteEdits = Literal["edits", "records"]
# A selection the existing drawing draws already: left out, or added again as a record of its own.
DuplicatesChoice = Literal["skip", "add"]
# Why a copy is refused:
# CIRCULAR_WHOLE_RECORDS: a Circular drawing takes several whole records, or one record or range alone;
# CIRCULAR_RANGE_TARGET: the existing Circular drawing draws a range, so it takes nothing more (PD-OI-090 (4));
# SAME_DRAWING: a copy goes to another drawing, and the source drawing does not change.
CopyRefusalReason = Literal["CIRCULAR_WHOLE_RECORDS", "CIRCULAR_RANGE_TARGET", "SAME_DRAWING"]


class CopyRefused(ValidationError):
    """A copy the drawing cannot take, by ``reason`` (the Web words it)."""

    _MESSAGES: Mapping[CopyRefusalReason, str] = {
        "CIRCULAR_WHOLE_RECORDS": "A Circular drawing draws several whole records, or one record or range alone.",
        "CIRCULAR_RANGE_TARGET": "The Circular drawing draws a range; it takes no other record or range.",
        "SAME_DRAWING": "A copy goes to another drawing; choose a drawing other than the source.",
    }

    def __init__(self, reason: CopyRefusalReason) -> None:
        super().__init__(self._MESSAGES[reason], diagnostic={"code": "INPUT_INVALID"})
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
    settings: SettingsChoice = "keep",
) -> RegionDrawing:
    """Derive a drawing of ``regions`` from a drawing of ``materialized``.

    ``drawing`` selects the source drawing by ID or name, as in
    :func:`gbdraw.session.session_to_request`; without it, the Session's only
    drawing. One region per record; the new drawing lists the records in source order.
    ``margin`` bp are added on each side and clamped to the record ends. The
    mode defaults to the source's mode when every region is a whole record, and
    to Linear when any is a range; a Circular drawing of several whole records
    is a Multi-Record Canvas, and a Circular drawing of a range draws that one
    range (else :class:`CopyRefused`). Each record keeps the source's
    orientation unless ``reverse_complement`` is given. With ``adapt_sizes``
    off, every hand-set size is kept. With ``settings="defaults"`` the drawing
    takes the mode's defaults and none of the source's settings or edits. A
    region whose start is after its end runs across the origin of a circular
    record and is drawn as one window (D-34); on a linear record it is refused.
    Use the result while ``materialized`` is active: its request reads the
    materialized files.
    """
    # The one read of the source drawing's settings and edits.
    derived = derive_region_request(
        session_to_request(materialized, drawing=drawing),
        regions,
        mode=mode,
        margin=margin,
        reverse_complement=reverse_complement,
        adapt_sizes=adapt_sizes,
        settings=settings,
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
    settings: SettingsChoice = "keep",
    circular_canvas: bool = True,
) -> RegionRequest:
    """Derive the request of a drawing of ``regions`` from the request of its source drawing.

    The one derivation behind :func:`derive_region_drawing` and the Web's New
    drawing from a region; the arguments are as there. Without
    ``circular_canvas`` a Circular drawing takes one record and the mode's
    default for several whole records is Linear: the Web's New drawing from a
    region asks so until its dialog offers the Multi-Record Canvas.
    """
    collection = resolve_request_records(source)
    selected = _resolve_regions(collection, regions, margin=margin)
    return _derive_selected(
        source,
        collection,
        selected,
        mode=mode,
        reverse_complement=reverse_complement,
        adapt_sizes=adapt_sizes,
        settings=settings,
        circular_canvas=circular_canvas,
    )


def _derive_selected(
    source: DiagramRequest,
    collection: ResolvedRecordCollection,
    selected: tuple[SourceRegion, ...],
    *,
    mode: DrawingMode | None,
    reverse_complement: bool | None,
    adapt_sizes: bool,
    settings: SettingsChoice,
    circular_canvas: bool,
) -> RegionRequest:
    """The request of a drawing of the resolved regions ``selected`` of ``source``."""
    if settings not in {"keep", "defaults"}:
        raise ValidationError(f"Unknown settings choice: {settings!r}.", diagnostic={"code": "INPUT_INVALID"})
    source_mode: DrawingMode = "linear" if isinstance(source, LinearDiagramRequest) else "circular"
    whole = all(region.whole for region in selected)
    target_mode: DrawingMode = mode or (
        source_mode if whole and (len(selected) == 1 or circular_canvas) else "linear"
    )
    if target_mode not in {"circular", "linear"}:
        raise ValidationError(f"Unknown drawing mode: {target_mode!r}.", diagnostic={"code": "INPUT_INVALID"})
    one_span = len(selected) == 1 and len(selected[0].spans) == 1
    if target_mode == "circular" and not (one_span or (circular_canvas and whole)):
        raise CopyRefused("CIRCULAR_WHOLE_RECORDS")
    defaults = settings == "defaults"
    records = tuple(_record_input(source, region, reverse_complement, defaults=defaults) for region in selected)
    # A Circular drawing of several records is a Multi-Record Canvas: the source canvas's layout without its
    # placements (each record of the new drawing takes its place in order), else the canvas defaults.
    layout = (
        replace(source.layout, multi_record_positions=None)
        if not defaults and isinstance(source, CircularDiagramRequest) and source.layout is not None
        else CircularMultiRecordOptions()
    ) if target_mode == "circular" and len(records) > 1 else None
    if defaults:
        # The mode's defaults draw the records: nothing of the source to adapt or to leave behind.
        request: DiagramRequest = (
            LinearDiagramRequest(records=records)
            if target_mode == "linear"
            else CircularDiagramRequest(records=records, layout=layout, grouping="single" if layout is None else "grid")
        )
        # The uploaded comparison and ring tables are not copied either, and are listed (PD-OI-090 (4)).
        tables = tuple(NotCarried(kind, source_mode) for kind, _count in _uploaded_tables(source.options))
        return RegionRequest(
            request=request, mode=target_mode, adaptation=SettingAdaptation(kept={}, reset=()), not_carried=tables
        )
    source_ids = [str(record.id) for record in collection.records]
    derived = DerivedRecords(
        mode=target_mode,
        records=records,
        record_ids=tuple(source_ids[region.item.resolved_index] for region in selected),
        lengths=tuple(length for region in selected for length in region.drawn_lengths),
        regions=selected,
        layout=layout,
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
    Linear row (without it each record takes a row of its own), or lays out a
    Circular Multi-Record Canvas (without it the Circular drawing draws one
    record).
    """

    mode: DrawingMode
    records: tuple[RecordInput, ...]
    record_ids: tuple[str, ...]
    lengths: tuple[int, ...]
    regions: tuple[SourceRegion, ...]
    layout: LinearMultiRecordOptions | CircularMultiRecordOptions | None = None


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
            records=records,
            options=LinearDiagramOptions(**values),
            output=output,
            layout=derived.layout if isinstance(derived.layout, LinearMultiRecordOptions) else None,
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
        # A Circular drawing of several records is a grid (the request's grouping rule).
        request = CircularDiagramRequest(
            records=records,
            options=CircularDiagramOptions(**values),
            output=output,
            layout=derived.layout if isinstance(derived.layout, CircularMultiRecordOptions) else None,
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
    *,
    defaults: bool = False,
) -> RecordInput:
    """The record input of ``region``: as the source draws the record, or, with ``defaults``, as a new drawing of
    the record's file draws it (forward, with the record's own caption and topology)."""
    item = region.item
    selector = item.selector
    if selector is None and item.source_record_count > 1:
        index = item.source_record_index
        selector = RecordSelector(raw=f"#{index + 1}", record_id=None, record_index=index)
    if reverse_complement is not None:
        reverse = reverse_complement
    elif defaults:
        reverse = False
    else:
        reverse = window_drawn_reverse(
            bool(item.region and item.region.reverse_complement), item.presentation.reverse_complement
        )
    # The windows of one record share its row in source order (D-01; the planner places them).
    windows = () if region.whole else tuple(RecordWindow(start, end, reverse) for start, end in region.spans)
    return RecordInput(
        source=source.records[item.input_index].source,
        selector=selector,
        windows=windows,
        presentation=replace(
            RecordPresentation() if defaults else item.presentation,
            reverse_complement=reverse and region.whole,
            grid_row=None,
            grid_column=None,
        ),
        record_key=item.record_identity_key,
        # A display start cannot be combined with a crop.
        # A window across the origin keeps the record circular, as the source drew it (D-34).
        display=(
            (RecordDisplayOptions() if defaults else item.display)
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



PasteNotCarriedKind = Literal[
    "ANNOTATIONS",
    "FEATURE_PLACEMENTS",
    "DEPTH",
    "COMPARISON_TABLES",
    "CONSERVATION_RING_TABLES",
    "NO_COMPARISON_ROWS",
]


@dataclass(frozen=True)
class PasteNotCarried:
    """One part a paste does not bring, or the pasted records lack, by kind, with its ``count``.

    ANNOTATIONS, FEATURE_PLACEMENTS and DEPTH are the source drawing's
    annotations, feature placements and Depth tracks of the copied records and
    regions; COMPARISON_TABLES and CONSERVATION_RING_TABLES its uploaded
    comparison and ring tables (counted by file or table). A paste carries
    none of them, whatever it brings (PD-OI-090 (4)). NO_COMPARISON_ROWS
    counts the pasted records the existing drawing's comparisons have no rows
    or pairs for. The CLI and API text is ``message``; the Web words the kind
    in its own wording owner.
    """

    kind: PasteNotCarriedKind
    count: int

    @property
    def message(self) -> str:
        return {
            "ANNOTATIONS": f"{self.count} annotation(s) of the source drawing.",
            "FEATURE_PLACEMENTS": f"{self.count} feature placement(s) of the source drawing.",
            "DEPTH": f"{self.count} Depth track(s) of the source drawing.",
            "COMPARISON_TABLES": f"{self.count} uploaded comparison table(s) of the source drawing.",
            "CONSERVATION_RING_TABLES": f"{self.count} uploaded ring table(s) of the source drawing.",
            "NO_COMPARISON_ROWS": f"The comparisons of the drawing have no rows for {self.count} pasted record(s).",
        }[self.kind]


@dataclass(frozen=True)
class AlreadyPresent:
    """A selected span the existing drawing draws already.

    ``record_key`` is the source drawing's record key; ``start`` and ``end``
    the span in source coordinates (``1`` to the record's length for a whole
    record, then ``whole``); ``target_record_key`` the record of the existing
    drawing that draws it.
    """

    record_key: str
    start: int
    end: int
    whole: bool
    target_record_key: str


@dataclass(frozen=True)
class LeftOut:
    """A selected span the existing drawing cannot take, with the reason, while it takes the others: a range into a
    Circular drawing of whole records (``CIRCULAR_WHOLE_RECORDS``). Source coordinates, as in :class:`AlreadyPresent`.
    """

    record_key: str
    start: int
    end: int
    reason: CopyRefusalReason


@dataclass(frozen=True)
class PasteRequest:
    """The existing drawing's request with the pasted records, their keys there, and the lists for the dialog."""

    request: DiagramRequest
    mode: DrawingMode
    pasted_record_keys: tuple[str, ...]
    already_present: tuple[AlreadyPresent, ...]
    left_out: tuple[LeftOut, ...]
    not_carried: tuple[PasteNotCarried, ...]


@dataclass(frozen=True)
class PastedDrawing:
    """The existing drawing with the pasted records, the Session with it in its place, and the lists.

    ``drawing`` is the target drawing with its ID, name and Web state (its
    settings draft, Legend and feature edits, ``ui`` and editor state), ready
    to stand in its place; ``document`` is the Session with it there, every
    other drawing as it was. When records were pasted the drawing's saved
    Result drew the old request, so it is dropped (``result_dropped``) and the
    next Generate draws the drawing again.
    """

    drawing: SessionDrawingSpec
    document: SessionDocument
    result_dropped: bool
    pasted_record_keys: tuple[str, ...]
    already_present: tuple[AlreadyPresent, ...]
    left_out: tuple[LeftOut, ...]
    not_carried: tuple[PasteNotCarried, ...]

    @property
    def dropped(self) -> tuple[str, ...]:
        """What was not carried, as the CLI and API word it."""
        return tuple(item.message for item in self.not_carried)


def paste_region_drawing(
    materialized: MaterializedSession,
    regions: Sequence[RegionChoice],
    *,
    target: str,
    drawing: str | None = None,
    margin: int = 0,
    reverse_complement: bool | None = None,
    edits: PasteEdits = "edits",
    duplicates: DuplicatesChoice = "skip",
) -> PastedDrawing:
    """Paste ``regions`` of a drawing of ``materialized`` into another of its drawings, ``target`` (by ID or name).

    ``drawing`` selects the source drawing as in :func:`derive_region_drawing`;
    the target is another drawing (:class:`CopyRefused`, ``SAME_DRAWING``).
    The records are appended after the target's own
    (:func:`paste_region_request`). The result keeps everything else of the
    target and of the Session (:class:`PastedDrawing`); save its ``document``.
    Use it while ``materialized`` is active: its request reads the
    materialized files.
    """
    document = materialized.document
    source_id = document.drawing(drawing).id
    target_id = document.drawing(target).id
    if source_id == target_id:
        raise CopyRefused("SAME_DRAWING")
    pasted = paste_region_request(
        session_to_request(materialized, drawing=source_id),
        regions,
        session_to_request(materialized, drawing=target_id),
        margin=margin,
        reverse_complement=reverse_complement,
        edits=edits,
        duplicates=duplicates,
    )
    replaced = replace_drawing_request(
        document, target_id, pasted.request, keep_result=not pasted.pasted_record_keys
    )
    return PastedDrawing(
        drawing=replaced.drawing,
        document=replaced.document,
        result_dropped=replaced.result_dropped,
        pasted_record_keys=pasted.pasted_record_keys,
        already_present=pasted.already_present,
        left_out=pasted.left_out,
        not_carried=pasted.not_carried,
    )


def paste_region_request(
    source: DiagramRequest,
    regions: Sequence[RegionChoice],
    target: DiagramRequest,
    *,
    margin: int = 0,
    reverse_complement: bool | None = None,
    edits: PasteEdits = "edits",
    duplicates: DuplicatesChoice = "skip",
) -> PasteRequest:
    """The request of the existing drawing ``target`` with ``regions`` of the drawing ``source`` appended.

    The target keeps its settings, palette, rules, tables and comparisons; the
    pasted records follow its own records, each on a row of its own (Linear) or
    on one new row of the canvas (Circular), and take keys of their own there.
    With ``edits="edits"`` each record comes as the source draws it, with the
    per-feature edits (color, label, hidden) inside its regions rebound to its
    new key. A per-feature color is a ``hash`` row of the color table: the
    rows that color a feature inside the regions are added after the target's
    rows, naming the feature of the pasted record. With ``edits="records"`` the
    records come as a new drawing of their files draws them. Either way the
    source's annotations, feature placements, Depth and uploaded tables stay,
    and are listed (``not_carried``). A selected span the target draws already
    (the same record of the same file, by file content, and the same span) is
    listed; ``duplicates="skip"`` leaves it out, ``"add"`` adds it as a record
    of its own. An aligned Linear target keeps its alignment and draws the
    pasted records unmoved; separate Circular diagrams take each record as a
    diagram of its own. A Circular target takes whole records, and only while
    it draws whole records: a range is left out and listed (``left_out``), and
    the paste is refused (:class:`CopyRefused`) when nothing whole is left.
    Neither request changes.
    """
    if edits not in {"edits", "records"}:
        raise ValidationError(f"Unknown edits choice: {edits!r}.", diagnostic={"code": "INPUT_INVALID"})
    if duplicates not in {"skip", "add"}:
        raise ValidationError(f"Unknown duplicates choice: {duplicates!r}.", diagnostic={"code": "INPUT_INVALID"})
    target_mode: DrawingMode = "linear" if isinstance(target, LinearDiagramRequest) else "circular"
    collection = resolve_request_records(source)
    selected = _resolve_regions(collection, regions, margin=margin)
    target_collection = resolve_request_records(target)
    left_out: tuple[LeftOut, ...] = ()
    if target_mode == "circular":
        if any(item.window is not None or item.region is not None for item in target_collection.provenance):
            raise CopyRefused("CIRCULAR_RANGE_TARGET")
        if not any(region.whole for region in selected):
            raise CopyRefused("CIRCULAR_WHOLE_RECORDS")
        left_out = tuple(
            LeftOut(region.item.record_identity_key, start, end, "CIRCULAR_WHOLE_RECORDS")
            for region in selected if not region.whole for start, end in region.spans
        )
        selected = _without_spans(selected, {(item.record_key, (item.start, item.end)) for item in left_out})
    present = _already_present(source, selected, target, target_collection)
    kept = selected if duplicates == "add" else _without_spans(
        selected, {(item.record_key, (item.start, item.end)) for item in present}
    )
    if not kept:
        return PasteRequest(target, target_mode, (), present, left_out, ())
    carried = _derive_selected(
        source,
        collection,
        kept,
        mode=target_mode,
        reverse_complement=reverse_complement,
        adapt_sizes=False,
        settings="keep" if edits == "edits" else "defaults",
        circular_canvas=True,
    ).request
    keys = _fresh_record_keys(
        [record.record_key or "" for record in carried.records],
        inputs={record.record_key for record in target.records if record.record_key is not None},
        units={item.record_identity_key for item in target_collection.provenance},
    )
    added = tuple(
        replace(record, record_key=key) for record, key in zip(carried.records, keys, strict=True)
    )
    rename = dict(zip((record.record_key for record in carried.records), keys, strict=True))
    rows = tuple(replace(row, record_key=rename[row.record_key]) for row in carried.options.feature_overrides)
    colors = _pasted_colors(
        source.options.colors if edits == "edits" else None,
        target.options.colors,
        [(region, rename[region.item.record_identity_key]) for region in kept],
    )
    units = sum(len(region.spans) for region in kept)
    pasted_options = _pasted_options(target.options, rows, colors, len(target_collection.records), units, len(added))
    request: DiagramRequest
    if isinstance(target, CircularBatchRequest):
        assert isinstance(pasted_options, CircularDiagramOptions)
        request = _pasted_diagrams(target, added, [collection.records[r.item.resolved_index].id for r in kept])
        request = replace(request, options=pasted_options)
    elif isinstance(target, LinearDiagramRequest):
        records, placed = _placed(target, added, kept, unit_count=len(target_collection.records))
        assert isinstance(placed, LinearMultiRecordOptions) or placed is None
        assert isinstance(pasted_options, LinearDiagramOptions)
        layout, alignment = _unmoved(placed, target.similarity_alignment, keys)
        request = replace(
            target, records=records, layout=layout, options=pasted_options, similarity_alignment=alignment
        )
    else:
        records, canvas = _placed(target, added, kept, unit_count=len(target_collection.records))
        assert isinstance(canvas, CircularMultiRecordOptions) and isinstance(pasted_options, CircularDiagramOptions)
        request = replace(target, records=records, layout=canvas, grouping="grid", options=pasted_options)
    not_carried = _paste_not_carried(source, collection, kept, target_mode, compares_by_rows=_compares_by_rows(target))
    return PasteRequest(request, target_mode, keys, present, left_out, not_carried)


UploadedTableKind = Literal["COMPARISON_TABLES", "CONSERVATION_RING_TABLES"]


def _uploaded_tables(options: object) -> tuple[tuple[UploadedTableKind, int], ...]:
    """The uploaded comparison and ring tables of a drawing, by kind and count (files or tables): its comparison and
    ring inputs when it has no LOSAT search, whose results it would otherwise hold. No copy carries them: a new
    drawing that keeps the settings lists them in its own words (``_linear_comparisons``, the ring tables)."""
    if getattr(options, "losat_search", None) is not None:
        return ()

    def count(names: Sequence[str]) -> int:
        values = [getattr(options, name, None) for name in names]
        return sum(len(value) if isinstance(value, (list, tuple)) else 1 for value in values if value is not None)

    found: tuple[tuple[UploadedTableKind, int], ...] = (
        ("COMPARISON_TABLES", count(_LINEAR_COMPARISON_INPUTS) if isinstance(options, LinearDiagramOptions) else 0),
        ("CONSERVATION_RING_TABLES", count(_RING_RESULTS) if isinstance(options, CircularDiagramOptions) else 0),
    )
    return tuple(item for item in found if item[1])


def _paste_not_carried(
    source: DiagramRequest,
    collection: ResolvedRecordCollection,
    kept: Sequence[SourceRegion],
    target_mode: DrawingMode,
    *,
    compares_by_rows: bool,
) -> tuple[PasteNotCarried, ...]:
    """What a paste of ``kept`` leaves behind, read from the source drawing whatever the paste brings and whatever
    the target's mode: the annotations, feature placements and Depth tracks of the copied records and regions, and
    the uploaded tables; and the pasted records the target's comparisons have no rows for."""
    options = source.options
    source_mode: DrawingMode = "linear" if isinstance(source, LinearDiagramRequest) else "circular"
    by_key = {region.item.record_identity_key: region for region in kept}
    source_ids = [str(record.id) for record in collection.records]
    annotations = _region_annotations(
        options.annotations, collection, kept, [], source_mode,
        tuple(source_ids[region.item.resolved_index] for region in kept),
    )
    depth = _region_depth_tracks(
        canonical_depth_tracks(options, record_count=len(collection.records)), kept, target_mode
    )
    counts: tuple[tuple[PasteNotCarriedKind, int], ...] = (
        ("ANNOTATIONS", sum(len(item.annotations) for item in annotations.sets) if annotations else 0),
        ("FEATURE_PLACEMENTS", sum(
            _keeps_feature(by_key, row.record_key, row.biological_feature_id) for row in options.feature_placements
        )),
        ("DEPTH", len(depth or ())),
        *_uploaded_tables(options),
        ("NO_COMPARISON_ROWS", len(kept) if compares_by_rows else 0),
    )
    return tuple(PasteNotCarried(kind, count) for kind, count in counts if count)


def _compares_by_rows(target: DiagramRequest) -> bool:
    """Whether the drawing's comparisons are results or pairs of its own records, which a pasted record has none of;
    a search with no pairs compares adjacent rows, the pasted ones too."""
    options = target.options
    search = getattr(options, "losat_search", None)
    return _has_any(options, _LINEAR_COMPARISON_INPUTS + _RING_RESULTS) or (
        isinstance(target, LinearDiagramRequest) and search is not None and search.pairs is not None
    )


def _already_present(
    source: DiagramRequest,
    selected: Sequence[SourceRegion],
    target: DiagramRequest,
    target_collection: ResolvedRecordCollection,
) -> tuple[AlreadyPresent, ...]:
    """The selected spans whose record and span (whole, or the same start and end) the target draws already.

    A record is the same when it is the same record (by its place) of a file with the same bytes
    (``SourceFileContents``): each drawing stages its own copy of a file. The files are read from the requests'
    record sources, and a file is read only when its size matches.
    """
    contents = SourceFileContents()
    drawn: list[tuple[RecordSourceKind, tuple[str, ...], int, tuple[int, int], str]] = []
    for item, record in zip(target_collection.provenance, target_collection.records, strict=True):
        kind, paths = contents.files(target.records[item.input_index].source)
        if not paths:
            continue
        length = int(item.source_length or len(record))
        span = (
            (item.window.span.start, item.window.span.end) if item.window is not None
            else (item.region.start, item.region.end) if item.region is not None
            else (1, length)
        )
        drawn.append((kind, paths, item.source_record_index, span, item.record_identity_key))
    present: list[AlreadyPresent] = []
    for region in selected:
        kind, paths = contents.files(source.records[region.item.input_index].source)
        for span in region.spans:
            target_key = next((
                key for drawn_kind, drawn_paths, index, drawn_span, key in drawn
                if (drawn_kind, index, drawn_span) == (kind, region.item.source_record_index, span)
                and paths and contents.same(paths, drawn_paths)
            ), None)
            if target_key is not None:
                key = region.item.record_identity_key
                present.append(AlreadyPresent(key, span[0], span[1], span == (1, region.length), target_key))
    return tuple(present)


def _without_spans(
    selected: Sequence[SourceRegion], spans: set[tuple[str, tuple[int, int]]]
) -> tuple[SourceRegion, ...]:
    """``selected`` without the given spans of their records, renumbered; a record left with none is left out."""
    kept: list[SourceRegion] = []
    first_unit = 0
    for region in selected:
        left = tuple(span for span in region.spans if (region.item.record_identity_key, span) not in spans)
        if left:
            kept.append(replace(
                region, index=len(kept), spans=left, whole=left == ((1, region.length),), first_unit=first_unit
            ))
            first_unit += len(left)
    return tuple(kept)


def _fresh_record_keys(keys: Sequence[str], *, inputs: set[str], units: set[str]) -> tuple[str, ...]:
    """A key for each pasted record that names no record of the drawing: ``<key>``, else ``<key>-2``, ``<key>-3``, ...

    A key is taken by a record the drawing draws (``units``) and by any key the
    Web reads as one of its records: a member ``<input key>:<n>`` of a record
    input that draws a whole file, or an input key ``<key>:<n>`` of which the
    pasted key would be the File card.
    """
    def taken(key: str) -> bool:
        return key in units or key in inputs or any(
            _expands(key, other) or _expands(other, key) for other in inputs
        )

    fresh: list[str] = []
    for key in keys:
        candidate, number = key, 1
        while taken(candidate) or candidate in fresh:
            number += 1
            candidate = f"{key}-{number}"
        fresh.append(candidate)
    return tuple(fresh)


def _expands(key: str, card: str) -> bool:
    """Whether ``key`` is the key of a member of the File card ``card`` (``expanded_record_key``)."""
    prefix = f"{card}{EXPANDED_RECORD_KEY_SEPARATOR}"
    return key.startswith(prefix) and key[len(prefix):].isdigit() and not key[len(prefix):].startswith("0")


def _pasted_options(
    options: CircularDiagramOptions | LinearDiagramOptions,
    rows: Sequence[FeatureOverride],
    colors: ColorOptions | None,
    record_count: int,
    added_units: int,
    added_inputs: int,
) -> CircularDiagramOptions | LinearDiagramOptions:
    """The existing drawing's options with the pasted records' feature edits ``rows`` after its own, and its color
    table with their color rows (``colors``, None for none); a per-record Depth track has none for the
    ``added_units`` drawn records, and a genetic code per record input is the runtime default for the
    ``added_inputs`` records."""
    search: LosatSearchOptions | None = options.losat_search
    if search is not None and len(search.record_gencodes) > 1:
        search = replace(search, record_gencodes=(*search.record_gencodes, *(None,) * added_inputs))
    pasted = replace(
        options,
        feature_overrides=(*options.feature_overrides, *rows),
        losat_search=search,
        colors=options.colors if colors is None else colors,
    )
    depth: list[DepthTrackInput] = []
    for track in canonical_depth_tracks(options, record_count=record_count) or ():
        source = track.source
        depth.append(replace(track, source=(*source, *(None,) * added_units)) if isinstance(source, tuple) else track)
    if not any(isinstance(track.source, tuple) for track in depth):
        return pasted
    # The depth tracks in their canonical form, which replaces the compatibility inputs.
    pasted = replace(
        pasted,
        depth_tracks=tuple(depth),
        depth_table=None,
        depth_file=None,
        depth_tables=None,
        depth_files=None,
        depth_track_tables=None,
        depth_track_files=None,
        depth_track_labels=None,
        depth_track_colors=None,
        depth_track_large_tick_intervals=None,
        depth_track_small_tick_intervals=None,
        depth_track_tick_font_sizes=None,
    )
    return replace(pasted, depth_track_heights=None) if isinstance(pasted, LinearDiagramOptions) else pasted


_COLOR_TABLE_COLUMNS = ("feature_type", "qualifier_key", "value", "color", "caption")


def _color_rows(colors: ColorOptions | None) -> DataFrame | None:
    """The rows of a drawing's color table, read as its render reads them."""
    if colors is None:
        return None
    if colors.color_table is not None:
        return colors.color_table
    return read_color_table(colors.color_table_file) if colors.color_table_file else None


def _pasted_colors(
    source: ColorOptions | None,
    target: ColorOptions | None,
    pasted: Sequence[tuple[SourceRegion, str]],
) -> ColorOptions | None:
    """The target's colors with a ``hash`` row for each feature inside a pasted region that a ``hash`` row of the
    source colors: the row names the feature in the pasted record (its key there), with the source row's color and
    caption, after the target's own rows. None when no such row comes along: the target's colors stay as they are.

    The source's rows are matched as the renderer matches them (``find_specific_color_rule``): a row naming the
    feature's instance in the source record, or an older row naming its stable hash, colors it.
    """
    table = _color_rows(source)
    if table is None or "qualifier_key" not in table.columns:
        return None
    hash_rules = build_specific_color_rules(table[table["qualifier_key"] == "hash"])
    added: list[tuple[str, str, str, str, str]] = []
    for region, key in pasted:
        for entry in region.item.source_feature_catalog:
            if not region.contains(entry.location_parts):
                continue
            matched = find_specific_color_rule(
                SimpleNamespace(type=entry.feature_type, qualifiers={}),
                hash_rules,
                selector={"hash": entry.stable_feature_id},
                instance=feature_instance(entry.biological_feature_id, region.item.record_identity_key),
            )
            if matched is not None:
                value = exact_value_pattern(feature_instance(entry.biological_feature_id, key))
                added.append((entry.feature_type, "hash", value, str(matched[0]), str(matched[1] or "")))
    if not added:
        return None
    rows = DataFrame(added, columns=list(_COLOR_TABLE_COLUMNS))
    own = _color_rows(target)
    merged = rows if own is None else pd.concat([own, rows], ignore_index=True)
    return replace(target or ColorOptions(), color_table=merged, color_table_file=None)


def _unmoved(
    layout: LinearMultiRecordOptions | None,
    plan: SimilarityAlignmentPlan | SimilarityAlignmentReference | None,
    keys: Sequence[str],
) -> tuple[LinearMultiRecordOptions | None, SimilarityAlignmentPlan | SimilarityAlignmentReference | None]:
    """An aligned drawing's translations and Similarity alignment plan, which name every record it draws, naming
    the pasted records ``keys`` too as records drawn where they are: no translation of their own, and the plan's
    decision for a record with no candidate (``resolve_similarity_alignment``)."""
    if layout is not None and layout.record_translations:
        layout = replace(layout, record_translations=(
            *layout.record_translations, *(LinearRecordTranslation(key) for key in keys)
        ))
    if not isinstance(plan, SimilarityAlignmentPlan):
        return layout, plan
    skipped = tuple(
        AlignmentRecordDecision(
            record_key=key,
            status=AlignmentDecisionStatus.SKIPPED,
            rationale=AlignmentResolutionRationale.SKIPPED_NO_CANDIDATE,
        )
        for key in keys
    )
    return layout, replace(plan, records=(*plan.records, *skipped))


def _pasted_diagrams(
    target: CircularBatchRequest, added: Sequence[RecordInput], record_ids: Sequence[str]
) -> CircularBatchRequest:
    """Separate Circular diagrams with a diagram of each pasted record after them: no grid place, and an output
    named by its record ID, as ``resolve_circular_batch_outputs`` names one, where the diagrams name each output."""
    records = (*target.records, *(
        replace(record, presentation=replace(record.presentation, grid_row=None, grid_column=None)) for record in added
    ))
    if not target.outputs:
        return replace(target, records=records)
    used = {(output.output_directory, output.output_prefix) for output in target.outputs}
    outputs = list(target.outputs)
    for record_id in record_ids:
        base = resolve_implicit_record_output_prefix(record_id)
        candidate, suffix = base, 2
        while (outputs[-1].output_directory, candidate) in used:
            candidate, suffix = f"{base}_{suffix}", suffix + 1
        used.add((outputs[-1].output_directory, candidate))
        outputs.append(replace(outputs[-1], output_prefix=candidate))
    return replace(target, records=records, outputs=tuple(outputs))


def _placed(
    target: CircularDiagramRequest | LinearDiagramRequest,
    added: Sequence[RecordInput],
    regions: Sequence[SourceRegion],
    *,
    unit_count: int,
) -> tuple[tuple[RecordInput, ...], LinearMultiRecordOptions | CircularMultiRecordOptions | None]:
    """The target's records with ``added`` after them, and its layout: where the target places its records by row,
    a Linear record takes a new row of its own and the Circular records one new row together; otherwise they follow
    in order. A Circular target becomes, or stays, a Multi-Record Canvas."""
    linear = isinstance(target, LinearDiagramRequest)
    layout = target.layout if linear or target.layout is not None else CircularMultiRecordOptions()
    rows = [record.presentation.grid_row for record in target.records]
    positions = layout.multi_record_positions if layout is not None else None
    new_rows = [index + 1 if linear else 1 for index in range(len(added))]
    if any(row is not None for row in rows):
        top = max(row for row in rows if row is not None)
        columns = any(
            record.presentation.grid_column is not None or any(window.column is not None for window in record.windows)
            for record in target.records
        )
        placed = []
        for index, (record, row) in enumerate(zip(added, new_rows, strict=True)):
            column = None if not columns else 1 if linear else index + 1
            if columns and len(record.windows) > 1:
                record = replace(
                    record, windows=tuple(replace(window, column=number) for number, window in enumerate(record.windows, 1))
                )
                column = None
            placed.append(replace(record, presentation=replace(record.presentation, grid_row=top + row, grid_column=column)))
        return (*target.records, *placed), layout
    if positions:
        top = max(parse_record_row_position(str(token))[1] for token in positions)
        tokens = list(positions)
        unit = unit_count
        for region, row in zip(regions, new_rows, strict=True):
            for _span in region.spans:
                unit += 1
                tokens.append(f"#{unit}@{top + row}")
        assert layout is not None
        layout = replace(layout, multi_record_positions=tuple(tokens))
    return (*target.records, *added), layout

__all__ = [
    "AlreadyPresent",
    "CopyRefusalReason",
    "CopyRefused",
    "DerivedRecords",
    "DuplicatesChoice",
    "FeatureRegionSelection",
    "LeftOut",
    "LinearComparisonRule",
    "NotCarried",
    "NotCarriedKind",
    "PasteEdits",
    "PasteNotCarried",
    "PasteNotCarriedKind",
    "PasteRequest",
    "PastedDrawing",
    "SettingsChoice",
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
    "paste_region_drawing",
    "paste_region_request",
]
