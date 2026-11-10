"""Resolve typed record inputs once while retaining their source provenance."""

from __future__ import annotations

import copy
import logging
from collections import Counter
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Callable, Hashable, Iterable, Literal, Mapping, MutableMapping, NamedTuple, Sequence, TypeAlias, TypeVar, cast

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.collinearity import (
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityResult,
)
from gbdraw.analysis.origin_arcs import location_center
from gbdraw.analysis.protein_colinearity import (
    OrthogroupGraphResult,
    OrthogroupMember,
    OrthogroupResult,
)
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.models.canvas import LINEAR_WRAP_LENGTH_PATH, LinearWrapLength, linear_wrap_length_from_config
from gbdraw.config.modify import modify_config_dict, raw_config_value, with_config_overrides
from gbdraw.exceptions import ValidationError
from gbdraw.analysis.windows import WINDOW_SOURCE_ANNOTATION, WindowSource, cut_content
from gbdraw.core.numeric import ascii_float, parse_ascii_integer
from gbdraw.core.record_metadata import (
    WINDOW_OF_ANNOTATION,
    WINDOW_SPAN_ANNOTATION,
    RECORD_ORIGIN_ANNOTATION,
    RecordOrigin,
    RecordSourceKind,
    held_record_origin,
    RecordContent,
    record_window_span,
    _feature_source_index_map,
    _iter_source_features,
    _SourceFeatureAnchorProfile,
    _read_coord_map,
    _source_period,
)
from gbdraw.io.cli_tables import (
    ConservationTable,
    RecordsTableRecord,
    read_comparisons_table,
    read_conservation_table,
    read_records_table,
)
from gbdraw.io.comparisons import read_comparison_resource
from gbdraw.io.genome import (
    feature_types_left_out, filter_features_by_type, load_fasta, load_gbks, parse_gff_fasta, require_record_sequence,
)
from gbdraw.io.record_sequences import detached_records
from gbdraw.io.record_select import (
    RecordSelector,
    drawn_record_units,
    in_draw_order,
    parse_record_selector,
    record_index_number,
    reverse_records,
    select_record_units,
)
from gbdraw.io.regions import (
    RecordWindow,
    RegionSpec,
    assign_region_specs,
    crop_record_windows,
    parse_region_spec,
    region_text,
)
from gbdraw.layout.linear_wrap import refuse_undrawable_wrap, resolve_wrap_length, wrap_blocking_inputs
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.record_transform import detected_topology, drawn_circular
from gbdraw.render_stages import render_stage
from gbdraw.layout.record_placement import FeatureAnchorPlacement, FeatureAnchorStatus, resolve_record_row_positions
from gbdraw.layout.similarity_alignment import (
    AlignmentAnchorIdentity,
    AlignmentDecisionStatus,
    AlignmentEvidenceEdge,
    SimilarityAlignmentCandidate,
    SimilarityAlignmentPlan,
)
from gbdraw.linear_comparison import LinearComparison, reverse_unbound_endpoint_span

from .config import apply_config_overrides
from .options import (
    LinearWrapDefaults,
    CircularDiagramOptions,
    DepthTrackInput,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
)
from gbdraw.features.source import SourceFeatureCatalog, SourceFeatureIdentity, capture_source_features
from gbdraw.annotations.models import FeatureSelector
from gbdraw.features.ids import compute_feature_hash_from_location_parts

from .prepared import (
    ParsedRecordInputs,
    PreparedResourceIdentity,
    get_or_build_parsed_source,
    prepared_resource_identity,
)
from .requests import (
    CircularBatchOutputPolicy,
    InputPlace,
    expanded_record_key,
    window_unit_key,
    FastaInputSource,
    GenBankInputSource,
    GffFastaInputSource,
    InMemoryRecordSource,
    RecordCardinality,
    RecordCollectionOptions,
    RecordDisplayOptions,
    RecordInput,
    RecordInputSource,
    RecordPresentation,
    RenderOutputRequest,
)

logger = logging.getLogger(__name__)

GenBankLoader = Callable[..., list[SeqRecord]]
GffFastaLoader = Callable[..., list[SeqRecord]]
FastaLoader = Callable[..., list[SeqRecord]]
# What a record was read from: the planner's name for each record source type.


@dataclass(frozen=True)
class ResolvedRecordDisplay:
    """Transient display intent resolved from complete source facts."""

    detected_topology: Literal["circular", "linear", "unknown"]
    is_circular: bool
    start_coordinate: int | None
    current_start_coordinate: int
    orientation_step: Literal[-1, 1]


def resolve_record_display(
    options: RecordDisplayOptions,
    *,
    source_length: int,
    detected_topology: Literal["circular", "linear", "unknown"],
    source_base: int,
    source_step: Literal[-1, 1],
    has_input_region: bool = False,
    has_collection_region: bool = False,
    is_cropped: bool = False,
) -> ResolvedRecordDisplay:
    """Validate supplied source facts without loading or transforming records.

    The caller supplies the full source length/topology, the existing affine
    map from core.record_metadata, and crop applicability for this record.
    Neither source facts nor crop state are inferred from a materialized length.
    current_start_coordinate is the existing affine start before extra rotation.
    """
    if not isinstance(options, RecordDisplayOptions):
        raise ValidationError("Record display has an unsupported type.")
    if detected_topology not in ("circular", "linear", "unknown"):
        raise ValidationError("Detected topology must be circular, linear, or unknown.")
    crop_flags = (has_input_region, has_collection_region, is_cropped)
    if not all(isinstance(value, bool) for value in crop_flags):
        raise ValidationError("Resolved crop flags must be booleans.")
    if options.start_coordinate is not None and any(crop_flags):
        raise ValidationError("An explicit display start cannot be combined with a crop.")
    is_circular = drawn_circular(options.is_circular, detected_topology)
    # The transform owns length, affine-domain, orientation and anchor validation.
    RecordDisplayTransform(
        source_length, source_base, source_step, options.start_coordinate, is_circular,
    )
    return ResolvedRecordDisplay(
        detected_topology=detected_topology,
        is_circular=is_circular,
        start_coordinate=options.start_coordinate,
        current_start_coordinate=source_base,
        orientation_step=source_step,
    )


@dataclass(frozen=True)
class ResolvedWindow:
    """The window of its record a drawn unit shows (D-01).

    ``record_key`` is the record's own key, shared by its units; ``index`` is
    the window's place among the record's ``count`` windows; ``span`` is the
    window. A record with one window draws one unit keyed like the record.
    """

    record_key: str
    index: int
    count: int
    span: RecordWindow

    @property
    def is_one_of_several(self) -> bool:
        """Whether the record draws more than one unit (each then has its own unit key)."""
        return self.count > 1

    @property
    def unit_key(self) -> str:
        """The key of the drawn unit of this window (``window_unit_key``)."""
        return window_unit_key(self.record_key, self.index, self.count)

    @property
    def window_of(self) -> str | None:
        """The key of the record this unit is one of several windows of; ``None`` when the record draws one unit."""
        return self.record_key if self.is_one_of_several else None

    @property
    def search_region(self) -> str:
        """The region a comparison search of this unit reads: the window's crop before its reverse complement
        (``region_text``). The Web's extraction crops by it (``LinearUnitWindow.searchRegion``), and a protein record
        analysis names it on the CLI and on the Web alike."""
        return region_text(self.span.start, self.span.end, False)


def record_drawn_windows(record_key: str, windows: Sequence[RecordWindow]) -> tuple[ResolvedWindow | None, ...]:
    """The drawn units of the record ``record_key`` drawn through ``windows``, in drawing order (D-01).

    One unit per window; a record without windows draws one unit, the record itself (``None``). The one owner of a
    record's units: the record resolver draws them, and the Web plans its comparison searches over them
    (``gbdraw.web_support.linear_units``).
    """

    return tuple(ResolvedWindow(record_key, index, len(windows), window) for index, window in enumerate(windows)) or (
        None,
    )


def refuse_file_comparison_of_windows(
    query_units: Sequence[int],
    subject_units: Sequence[int],
    *,
    records: tuple[str, str],
    where: str,
    row: int | None = None,
) -> None:
    """Refuse a comparison file between records one of which is drawn in several windows: the file gives the
    coordinates of one drawn unit, so such a pair is searched (``comparison_unit_pairs``)."""

    if len(query_units) > 1 or len(subject_units) > 1:
        windowed = records[0] if len(query_units) > 1 else records[1]
        raise ValidationError(
            f"{where}: record {windowed!r} is drawn in several windows, and a comparison file gives the coordinates "
            "of one drawn record; use source=losat for it.",
            diagnostic={
                "code": "COMPARISON_INPUT",
                "reason": "WINDOWED_FILE_COMPARISON",
                **({} if row is None else {"row": row}),
            },
        )


def comparison_unit_pairs(query_units: Sequence[int], subject_units: Sequence[int]) -> tuple[tuple[int, int], ...]:
    """The drawn-unit pairs a comparison of two records searches: every window of one against every window of the
    other (D-01, D-37), in the order of the query's units, then the subject's."""

    return tuple((query, subject) for query in query_units for subject in subject_units)


@dataclass(frozen=True)
class ResolvedRecordProvenance:
    """Stable source coordinates for one displayed record (one drawn unit of a record with windows)."""

    resolved_index: int
    input_index: int
    source_record_index: int
    source_record_count: int
    source_record_id: str
    source_kind: RecordSourceKind
    source_paths: tuple[str, ...]
    record_key: str
    cardinality: RecordCardinality
    selector: RecordSelector | None
    region: RegionSpec | None
    presentation: RecordPresentation
    display: RecordDisplayOptions = field(default_factory=RecordDisplayOptions)
    source_length: int | None = None
    # The length of every record of the source, each uncropped: the genome a source file is (OV-3143).
    source_file_length: int | None = None
    detected_topology: Literal["circular", "linear", "unknown"] = "unknown"
    is_cropped: bool = False
    has_collection_region: bool = False
    # No other record of the source has this record's ID, so the ID selects it.
    source_record_id_unique: bool = True
    resolved_display: ResolvedRecordDisplay | None = None
    source_features: SourceFeatureCatalog = field(kw_only=True, repr=False)
    # GFF3 features whose type the request does not draw, left out when the record
    # was planned, by type (the GFF3 type filter); the input summary counts them.
    feature_types_left_out_at_load: tuple[tuple[str, int], ...] = field(default=(), kw_only=True)
    # The window this unit draws, when the record is drawn through windows; ``record_key`` is then the unit's key.
    window: ResolvedWindow | None = field(default=None, kw_only=True)
    # The record whose window an earlier plan drew this record as, when this plan draws it whole (``RecordOrigin``).
    origin_window_of: str | None = field(default=None, kw_only=True)
    # The record's place among the resolved records before its windows expand: labels and subtitles apply by it.
    # ``None`` reads as ``resolved_index`` (``place``).
    record_position: int | None = field(default=None, kw_only=True)

    @property
    def place(self) -> int:
        """The record's place among the resolved records, by which collection labels and subtitles apply."""
        return self.resolved_index if self.record_position is None else self.record_position

    @property
    def record_identity_key(self) -> str:
        """The key of the record this unit draws: per-record edits, colours and evidence bind to it."""
        return self.window.record_key if self.window is not None else self.record_key

    @property
    def window_of(self) -> str | None:
        """The key of the record whose window this unit draws (D-01), from this plan or the one before."""
        return self.window.window_of if self.window is not None else self.origin_window_of

    def origin(self, record: SeqRecord, content: RecordContent | None = None) -> RecordOrigin:
        """The origin of ``record``, this unit's planned record, as the projection carries it; ``content`` is the
        record's ``RecordContent`` when the caller holds it already."""
        return RecordOrigin(
            content=content if content is not None else RecordContent(record),
            source_kind=self.source_kind,
            source_paths=self.source_paths,
            source_record_index=self.source_record_index,
            source_record_count=self.source_record_count,
            source_record_id=self.source_record_id,
            source_record_id_unique=self.source_record_id_unique,
            window_of=self.window_of,
            source_file_length=self.source_file_length,
        )

    @property
    def source_name(self) -> str:
        """The file name of the record's source, as the input summary names it."""
        return Path(self.source_paths[0]).name if self.source_paths else "in-memory record"

    @property
    def is_circular(self) -> bool:
        """Whether the record is drawn as circular: its display option, else its detected topology."""
        display = self.resolved_display
        return display.is_circular if display is not None else self.detected_topology == "circular"

    @property
    def source_feature_catalog(self) -> tuple[SourceFeatureIdentity, ...]:
        """The identities of the original source's features, built on first read."""
        return self.source_features.entries


@dataclass(frozen=True)
class ResolvedRecordCollection:
    """Flattened records and one provenance entry per displayed record."""

    records: tuple[SeqRecord, ...]
    provenance: tuple[ResolvedRecordProvenance, ...]
    displays: tuple[ResolvedRecordDisplay, ...] = field(init=False)
    transforms: tuple[RecordDisplayTransform, ...] = field(init=False)

    def __post_init__(self) -> None:
        if not self.records:
            raise ValidationError("A resolved record collection cannot be empty.")
        if len(self.records) != len(self.provenance):
            raise ValidationError(
                "Resolved record provenance must align with displayed records."
            )
        provenance = []
        displays = []
        transforms = []
        for record, item in zip(self.records, self.provenance, strict=True):
            base, step = cast("tuple[int, Literal[-1, 1]]", _read_coord_map(record))
            cropped = item.is_cropped or bool(record.annotations.get("gbdraw_region_applied"))
            length = item.source_length
            if length is None and not cropped:
                length = len(record)
            try:
                display = resolve_record_display(
                    item.display,
                    # A None length is rejected by RecordDisplayTransform.
                    source_length=cast(int, length),
                    detected_topology=item.detected_topology,
                    source_base=base,
                    source_step=step,
                    has_input_region=item.region is not None,
                    has_collection_region=item.has_collection_region,
                    is_cropped=cropped,
                )
                transform = RecordDisplayTransform(
                    cast(int, length), base, step, display.start_coordinate, display.is_circular,
                    crosses_origin=_source_period(record) is not None,
                )
            except ValidationError as exc:
                raise ValidationError(
                    f"Record {item.record_key!r} (input {item.input_index + 1}, "
                    f"source {item.source_paths or item.source_kind}, "
                    f"selector {item.selector}, record {item.source_record_id!r}): {exc}"
                ) from exc
            displays.append(display)
            transforms.append(transform)
            provenance.append(replace(item, resolved_display=display))
        object.__setattr__(self, "provenance", tuple(provenance))
        object.__setattr__(self, "displays", tuple(displays))
        object.__setattr__(self, "transforms", tuple(transforms))


# The record key of the n-th record input when it names none (`record-<n>`),
# and the separator before the source record number of each record an ALL
# input expands to (`<key>:<k>`). The Web Source recipe names a CLI record by
# them (gbdraw/web/js/contracts/schema.generated.js).
CLI_RECORD_KEY_PREFIX = "record-"

_SOURCE_BINDING_COLUMNS = frozenset(
    f"{role}_{name}" for role in ("query", "subject") for name in ("feature_index", "feature_svg_id")
)


SourceFeatureSpan: TypeAlias = tuple[int, int, int | None]


@dataclass(frozen=True)
class ComparisonRecordBinding:
    """A drawn record's source feature IDs, view feature IDs and spans, by source feature index.

    ``current`` holds the view IDs of the record as drawn (in the search frame
    when the binding was made for it), ``opposite`` those of the other strand;
    ``current_spans`` and ``opposite_spans`` hold the 0-based [start, end) and
    strand of each feature in the same two frames. ``reverse`` is whether the
    drawn record is reverse-complemented; it decides items without tags.
    """

    source: Mapping[int, str]
    current: Mapping[int, str]
    opposite: Mapping[int, str]
    length: int
    current_spans: Mapping[int, SourceFeatureSpan]
    opposite_spans: Mapping[int, SourceFeatureSpan]
    reverse: bool


def comparison_record_binding(
    record: SeqRecord,
    source_features: SourceFeatureCatalog,
    search_frame: bool,
) -> ComparisonRecordBinding:
    """The feature binding of one record (render stage ``comparison_rows``)."""

    current: dict[int, str] = {}
    opposite: dict[int, str] = {}
    current_spans: dict[int, SourceFeatureSpan] = {}
    opposite_spans: dict[int, SourceFeatureSpan] = {}
    indexes = _feature_source_index_map(record.features)
    for feature in _iter_source_features(record.features):
        if feature.location is None:
            continue
        index = indexes[id(feature)]
        for target, spans, location in ((current, current_spans, feature.location),
                                        (opposite, opposite_spans, feature.location._flip(len(record)))):
            target[index] = compute_feature_hash_from_location_parts(
                feature.type,
                [(int(part.start), int(part.end), part.strand) for part in location.parts],
                record_id=record.id,
            )
            spans[index] = (int(location.start), int(location.end),
                            location.strand if location.strand in (-1, 1) else None)
    source = source_features.source_hashes
    reverse = record_is_reverse_complemented(record)
    if search_frame and reverse:
        current, opposite = opposite, current
        current_spans, opposite_spans = opposite_spans, current_spans
    return ComparisonRecordBinding(
        source, current, opposite, len(record), current_spans, opposite_spans, reverse,
    )


def record_is_reverse_complemented(record: SeqRecord) -> bool:
    """Whether a resolved record is drawn reverse-complemented (its search frame differs)."""

    return _read_coord_map(record)[1] == -1


def comparison_record_bindings(
    records: Sequence[SeqRecord],
    provenance: Sequence[ResolvedRecordProvenance],
    *,
    search_frame: bool,
) -> tuple[ComparisonRecordBinding, ...]:
    """One binding per record: drawn views, or each record's search frame."""

    return tuple(
        render_stage(
            "comparison_rows",
            comparison_record_binding,
            record,
            item.source_features,
            search_frame,
        )
        for record, item in zip(records, provenance, strict=True)
    )


def project_comparison_frame(
    frame: pd.DataFrame,
    query: ComparisonRecordBinding,
    subject: ComparisonRecordBinding,
    drawn_feature_ids: bool = False,
) -> pd.DataFrame:
    """Rows of ``frame`` in the views of the bindings (render stage ``comparison_rows``).

    With ``drawn_feature_ids`` (``LinearComparison.drawn_feature_ids``), a side
    without a ``*_view_feature_svg_id`` column reads its ``*_feature_svg_id``
    as the view ID, takes the source ID of its feature index, and gets both
    written. Returns ``frame`` itself when no row changes.
    """

    if frame.empty or not _SOURCE_BINDING_COLUMNS <= set(frame.columns):
        return frame
    updated = None
    for row_index, row in frame.iterrows():
        flips = []
        for role, binding, prefix in (("query", query, "q"), ("subject", subject, "s")):
            parsed_indices = [
                parse_ascii_integer(value.strip()) for value in str(row[f"{role}_feature_index"]).split(";")
            ]
            indices = [index for index in parsed_indices if index is not None]
            if len(indices) != len(parsed_indices):
                raise ValidationError("Comparison source feature indexes must be integers.")
            view_column = f"{role}_view_feature_svg_id"
            drawn_ids = drawn_feature_ids and view_column not in frame.columns
            view_ids = str(row.get(view_column, row[f"{role}_feature_svg_id"])).split(";")
            source_ids = ([str(binding.source.get(index, "")) for index in indices] if drawn_ids
                          else str(row[f"{role}_feature_svg_id"]).split(";"))
            if not len(indices) == len(source_ids) == len(view_ids):
                raise ValidationError("Comparison source feature binding has inconsistent coverage.")
            if any(binding.source.get(index) != source_id
                   for index, source_id in zip(indices, source_ids, strict=True)):
                raise ValidationError("Comparison source feature index conflicts with its source feature ID.")
            current_ids, opposite_ids = ([views.get(index) for index in indices]
                                         for views in (binding.current, binding.opposite))
            flip = view_ids != current_ids
            if flip and view_ids != opposite_ids:
                raise ValidationError("Comparison view feature IDs do not match the current source/crop binding.")
            flips.append(flip)
            if not flip and not drawn_ids:
                continue
            if updated is None:
                updated = frame.copy(deep=True)
            if drawn_ids:
                updated.at[row_index, f"{role}_feature_svg_id"] = ";".join(source_ids)
            if not flip:
                updated.at[row_index, view_column] = ";".join(view_ids)
                continue
            for coordinate_field in (f"{prefix}start", f"{prefix}end"):
                coordinate = row[coordinate_field]
                if not 1 <= coordinate <= binding.length or int(coordinate) != coordinate:
                    raise ValidationError("Comparison view coordinate must be a finite genomic base within its record.")
                updated.at[row_index, coordinate_field] = binding.length + 1 - int(coordinate)
            # Every ID is present: view_ids matched opposite_ids, which has the same keys.
            updated.at[row_index, view_column] = ";".join(
                cast("list[str]", current_ids)
            )
        if updated is not None and flips[0] != flips[1] and "collinearity_orientation" in frame.columns:
            orientation = row["collinearity_orientation"]
            if orientation in ("plus", "minus"):
                updated.at[row_index, "collinearity_orientation"] = "minus" if orientation == "plus" else "plus"
    return frame if updated is None else updated


class UnmatchedFrameItem(NamedTuple):
    """A stored item whose tag matches neither frame of its record's binding."""

    kind: Literal["member", "anchor", "block"]
    record_index: int
    feature_index: int | None
    protein_id: str


@dataclass(frozen=True)
class FrameProjection:
    """Options with their source-bound results projected, and what stayed unmatched."""

    options: LinearDiagramOptions
    unmatched: tuple[UnmatchedFrameItem, ...] = ()


def has_source_bound_results(
    comparisons: Sequence[LinearComparison] | None,
    protein_comparisons: Sequence[pd.DataFrame] | None,
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None,
    collinearity: CollinearityResult | Sequence[CollinearityBlock] | None,
) -> bool:
    """Whether any table row or typed result is bound to source features."""

    frames = [comparison.matches for comparison in comparisons or ()]
    frames.extend(protein_comparisons or ())
    return (
        orthogroups is not None
        or bool(collinearity)
        or any(not frame.empty and _SOURCE_BINDING_COLUMNS <= set(frame.columns) for frame in frames)
    )


def _binding(bindings: Sequence[ComparisonRecordBinding], index: int) -> ComparisonRecordBinding:
    if not 0 <= index < len(bindings):
        raise ValidationError("Comparison record endpoint is outside the resolved collection.")
    return bindings[index]


class _AnchorSide(NamedTuple):
    record_index: int
    feature_index: int | None
    view_id: str
    start: int
    end: int
    strand: int | None
    protein_id: str


def _anchor_side(anchor: CollinearityAnchor, query: bool) -> _AnchorSide:
    if query:
        return _AnchorSide(anchor.query_record_index, anchor.query_feature_index,
                           anchor.query_view_feature_svg_id, anchor.query_start, anchor.query_end,
                           anchor.query_strand, anchor.query_protein_id)
    return _AnchorSide(anchor.subject_record_index, anchor.subject_feature_index,
                       anchor.subject_view_feature_svg_id, anchor.subject_start, anchor.subject_end,
                       anchor.subject_strand, anchor.subject_protein_id)


class _ResultProjector:
    """Projects the typed results of one request, item by item by tag (OV-1145).

    A member is tagged by its span against its feature's spans, an anchor
    side by its view ID; an item matching ``opposite`` moves to ``current``.
    An anchor side without tags (Session 39) follows the record-level rule of
    :func:`gbdraw.linear_comparison.reverse_unbound_endpoint_span`. An item
    matching neither frame, or a block whose sides do not each agree, is kept
    unchanged and named in ``unmatched``. Ranks and IDs never change.
    """

    def __init__(self, bindings: Sequence[ComparisonRecordBinding]) -> None:
        self.bindings = bindings
        self.unmatched: list[UnmatchedFrameItem] = []
        self._orthogroups: dict[int, OrthogroupResult | OrthogroupGraphResult] = {}

    def member(self, member: OrthogroupMember) -> OrthogroupMember:
        binding = _binding(self.bindings, int(member.record_index))
        index = int(member.feature_index)
        span = (int(member.start), int(member.end), member.strand if member.strand in (-1, 1) else None)
        if span == binding.current_spans.get(index):
            return member
        if span != binding.opposite_spans.get(index):
            self.unmatched.append(UnmatchedFrameItem("member", member.record_index, index, member.protein_id))
            return member
        start, end, strand = binding.current_spans[index]
        return replace(member, start=start, end=end, strand=strand)

    def orthogroups(
        self, result: OrthogroupResult | OrthogroupGraphResult,
    ) -> OrthogroupResult | OrthogroupGraphResult:
        known = self._orthogroups.get(id(result))
        if known is not None:
            return known
        projected: dict[int, OrthogroupMember] = {}

        def project(member: OrthogroupMember) -> OrthogroupMember:
            if id(member) not in projected:
                projected[id(member)] = self.member(member)
            return projected[id(member)]

        groups = {group_id: [project(member) for member in members]
                  for group_id, members in result.orthogroups.items()}
        by_protein = {protein_id: project(member) for protein_id, member in result.member_by_protein_id.items()}
        changed = any(
            after is not before
            for before, after in zip(
                (member for members in result.orthogroups.values() for member in members),
                (member for members in groups.values() for member in members),
                strict=True,
            )
        ) or any(by_protein[key] is not member for key, member in result.member_by_protein_id.items())
        output = replace(result, orthogroups=groups, member_by_protein_id=by_protein) if changed else result
        self._orthogroups[id(result)] = output
        return output

    def _side_flip(self, side: _AnchorSide) -> bool | None:
        """Whether an anchor side moves to the other frame; None when unmatched."""

        binding = _binding(self.bindings, int(side.record_index))
        if not side.view_id and side.feature_index is None:
            return binding.reverse
        if side.view_id and side.feature_index is not None:
            if side.view_id == binding.current.get(int(side.feature_index)):
                return False
            if side.view_id == binding.opposite.get(int(side.feature_index)):
                return True
        self.unmatched.append(UnmatchedFrameItem(
            "anchor", int(side.record_index), side.feature_index, str(side.protein_id),
        ))
        return None

    def _flipped(self, side: _AnchorSide) -> tuple[int, int, int | None, str]:
        binding = self.bindings[int(side.record_index)]
        start, end, strand = reverse_unbound_endpoint_span(
            int(side.start), int(side.end), side.strand, (binding.length, True),
        )
        view_id = side.view_id
        if view_id and side.feature_index is not None:
            view_id = binding.current[int(side.feature_index)]
        return start, end, strand, view_id

    def anchor(self, anchor: CollinearityAnchor) -> tuple[CollinearityAnchor, tuple[bool, bool] | None]:
        query, subject = _anchor_side(anchor, True), _anchor_side(anchor, False)
        query_flip, subject_flip = self._side_flip(query), self._side_flip(subject)
        if query_flip is None or subject_flip is None:
            return anchor, None
        projected = anchor
        if query_flip:
            start, end, strand, view_id = self._flipped(query)
            projected = replace(projected, query_start=start, query_end=end, query_strand=strand,
                                query_view_feature_svg_id=view_id)
        if subject_flip:
            start, end, strand, view_id = self._flipped(subject)
            projected = replace(projected, subject_start=start, subject_end=end, subject_strand=strand,
                                subject_view_feature_svg_id=view_id)
        return projected, (query_flip, subject_flip)

    def block(self, block: CollinearityBlock) -> CollinearityBlock:
        before = len(self.unmatched)
        projected = [self.anchor(anchor) for anchor in block.anchors]
        flips = {sides for _, sides in projected if sides is not None}
        if len(self.unmatched) != before:
            return block
        if len(flips) > 1:
            first = block.anchors[0]
            self.unmatched.append(UnmatchedFrameItem(
                "block", int(block.query_record_index), first.query_feature_index, str(first.query_protein_id),
            ))
            return block
        if not flips or flips == {(False, False)}:
            return block
        ((query_flip, subject_flip),) = flips
        orientation = block.orientation
        if query_flip != subject_flip:
            orientation = "minus" if orientation == "plus" else "plus"
        return replace(block, anchors=tuple(anchor for anchor, _ in projected), orientation=orientation)

    def collinearity(
        self, value: CollinearityResult | Sequence[CollinearityBlock],
    ) -> CollinearityResult | Sequence[CollinearityBlock]:
        blocks = value.blocks if isinstance(value, CollinearityResult) else tuple(value)
        projected_blocks = tuple(self.block(block) for block in blocks)
        changed = any(before is not after for before, after in zip(blocks, projected_blocks, strict=True))
        if not isinstance(value, CollinearityResult):
            return projected_blocks if changed else value
        unblocked = tuple(self.anchor(anchor)[0] for anchor in value.unblocked_anchors)
        orthogroups = None if value.orthogroups is None else self.orthogroups(value.orthogroups)
        changed = (changed or orthogroups is not value.orthogroups
                   or any(before is not after for before, after in
                          zip(value.unblocked_anchors, unblocked, strict=True)))
        if not changed:
            return value
        return replace(value, blocks=projected_blocks, unblocked_anchors=unblocked, orthogroups=orthogroups)


@dataclass(frozen=True)
class _TypedResultProjection:
    """The typed results of one request in the frames of its bindings.

    ``None`` keeps the given result: it was absent or already in that frame.
    """

    orthogroups: OrthogroupResult | OrthogroupGraphResult | None
    collinearity: CollinearityResult | Sequence[CollinearityBlock] | None
    unmatched: tuple[UnmatchedFrameItem, ...]


def _project_typed_results(
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None,
    collinearity: CollinearityResult | Sequence[CollinearityBlock] | None,
    bindings: tuple[ComparisonRecordBinding, ...],
) -> _TypedResultProjection:
    """Typed results in the frames of ``bindings`` (render stage ``comparison_rows``).

    As a stage, a projected result has the identity of this stage's key (the
    results' and bindings' identities), so a stage that takes it keys that
    identity instead of walking every member on each render (OV-1582).
    """

    projector = _ResultProjector(bindings)
    projected_orthogroups = None if orthogroups is None else projector.orthogroups(orthogroups)
    projected_collinearity = None if collinearity is None else projector.collinearity(collinearity)
    return _TypedResultProjection(
        None if projected_orthogroups is orthogroups else projected_orthogroups,
        None if projected_collinearity is collinearity else projected_collinearity,
        tuple(projector.unmatched),
    )


def project_source_bound_results(
    options: LinearDiagramOptions,
    bindings: Sequence[ComparisonRecordBinding],
) -> FrameProjection:
    """Project source-bound tables and typed results into the frames of ``bindings``.

    Drawn bindings give the drawn views (the planner); search-frame bindings
    give each record's search frame (every Session writer, R5). Each table row
    and typed item tells its frame by its tags, so either direction is exact
    and an item already in the target frame stays. Source identities, scores,
    ranks and block membership remain unchanged. Coordinates are record-local;
    the existing renderer alone applies a circular display cut. Returns
    ``options`` itself when nothing changes. Unmatched table rows raise; an
    unmatched typed item stays as given and is listed (OV-1145).
    """

    def project(frame: pd.DataFrame, query_index: int, subject_index: int,
                drawn_feature_ids: bool = False) -> pd.DataFrame:
        if frame.empty or not _SOURCE_BINDING_COLUMNS <= set(frame.columns):
            return frame
        return render_stage(
            "comparison_rows",
            project_comparison_frame,
            frame,
            _binding(bindings, query_index),
            _binding(bindings, subject_index),
            drawn_feature_ids,
        )

    # A projected table holds source and view IDs, so it is no longer read
    # with the drawn-frame IDs of a Session up to 40.
    explicit = tuple(replace(comparison, matches=project(comparison.matches,
        comparison.query_record_index, comparison.subject_record_index, comparison.drawn_feature_ids),
        drawn_feature_ids=False)
        for comparison in options.linear_comparisons or ())
    protein = None if options.protein_comparisons is None else tuple(
        project(frame, index, index + 1) for index, frame in enumerate(options.protein_comparisons))
    typed = (
        _TypedResultProjection(None, None, ())
        if options.orthogroups is None and options.collinearity_blocks is None
        else render_stage("comparison_rows", _project_typed_results, options.orthogroups,
                          options.collinearity_blocks, tuple(bindings))
    )
    orthogroups = options.orthogroups if typed.orthogroups is None else typed.orthogroups
    collinearity = options.collinearity_blocks if typed.collinearity is None else typed.collinearity
    changed = (
        any(before.matches is not after.matches or before.drawn_feature_ids for before, after in
            zip(options.linear_comparisons or (), explicit, strict=True))
        or any(before is not after for before, after in
               zip(options.protein_comparisons or (), protein or (), strict=True))
        or orthogroups is not options.orthogroups
        or collinearity is not options.collinearity_blocks
    )
    projected = replace(
        options,
        linear_comparisons=explicit if options.linear_comparisons is not None else None,
        protein_comparisons=protein,
        orthogroups=orthogroups,
        collinearity_blocks=collinearity,
    ) if changed else options
    return FrameProjection(projected, typed.unmatched)


def alignment_record_keys(collection: ResolvedRecordCollection) -> tuple[str, ...]:
    """The records a Similarity alignment plan covers, in drawing order: the windows of one record are one (D-50)."""

    return tuple(dict.fromkeys(item.record_identity_key for item in collection.provenance))


def project_similarity_alignment_centers(
    collection: ResolvedRecordCollection,
    plan: SimilarityAlignmentPlan | None,
) -> tuple[float | None, ...]:
    """Project selected anchor centers from already resolved record orientations.

    One centre per drawn unit: a record drawn in several windows has it in the window that holds its anchor
    (``alignment_anchor_unit``) and ``None`` in its other windows, which move with it (D-50).
    """

    centers: list[float | None] = [None] * len(collection.records)
    if plan is None:
        return tuple(centers)
    plan.validate_record_coverage(alignment_record_keys(collection))
    for decision in plan.records:
        if decision.status is AlignmentDecisionStatus.SKIPPED:
            continue
        if decision.anchor is None:
            raise ValidationError(
                f"Similarity alignment record {decision.record_key!r} has no anchor.",
                diagnostic={"code": "INPUT_INVALID"},
            )
        index, center = alignment_anchor_unit(collection, decision.anchor)
        if center is None:
            raise ValidationError(
                "Similarity alignment anchor is outside the current crop for record "
                f"{decision.record_key!r}.",
                diagnostic={"code": "INPUT_INVALID"},
            )
        centers[index] = center
    return tuple(centers)


def alignment_anchor_unit(
    collection: ResolvedRecordCollection, anchor: AlignmentAnchorIdentity,
) -> tuple[int, float | None]:
    """The drawn unit of the anchor's record that holds the anchor, and the anchor's centre on it (D-50).

    A record drawn whole or in one window has one unit. A record drawn in several windows is anchored in the first
    window, in drawing order, that holds the anchor's centre; when none holds it the centre is ``None`` (the anchor
    is unusable) and the unit is the record's first.
    """

    units = [index for index, item in enumerate(collection.provenance) if item.record_identity_key == anchor.record_key]
    if not units:
        raise ValidationError("Alignment anchor record coverage is invalid.", diagnostic={"code": "INPUT_INVALID"})
    position, center = first_unit_holding(_anchor_center_on_unit(collection, anchor, index) for index in units)
    return units[position], center


def first_unit_holding(offsets: Iterable[float | None]) -> tuple[int, float | None]:
    """The first drawn unit of a record, in drawing order, that holds a source position, and the position's display
    offset on it (D-50).

    ``offsets`` are the position's display offsets on the record's units in drawing order, ``None`` on a unit that
    does not hold it. When no unit holds it the result is the first unit and ``None``. The renderer
    (``alignment_anchor_unit``) anchors a record drawn in several windows by this rule, and the Web helper reads its
    placements from the renderer (``project_similarity_alignment_anchor_fact``).
    """

    for position, offset in enumerate(offsets):
        if offset is not None:
            return position, offset
    return 0, None


def _anchor_center_on_unit(
    collection: ResolvedRecordCollection, anchor: AlignmentAnchorIdentity, index: int,
) -> float | None:
    """The anchor's display centre on drawn unit ``index``, ``None`` when the unit does not hold it."""

    feature = _alignment_source_feature(anchor, collection.provenance[index])
    projected = _source_feature_display_center(feature, collection.transforms[index])
    return projected if 0.0 <= projected <= float(len(collection.records[index])) else None


def _source_feature_display_center(
    feature: SourceFeatureIdentity, transform: RecordDisplayTransform,
) -> float:
    """The display offset of a source feature's middle, the one alignment centre."""

    return transform.source_position_to_display_offset(
        location_center(feature.location_parts, transform.length)
    )


def _alignment_source_feature(
    anchor: AlignmentAnchorIdentity | None,
    provenance: ResolvedRecordProvenance,
) -> SourceFeatureIdentity:
    if anchor is None:
        raise ValidationError(
            f"Similarity alignment record {provenance.record_key!r} has no anchor."
        )
    catalog = provenance.source_feature_catalog
    if anchor.source_feature_index is not None:
        matches = [
            feature
            for feature in catalog
            if feature.source_feature_index == anchor.source_feature_index
        ]
    elif anchor.stable_feature_svg_id is not None:
        matches = [
            feature
            for feature in catalog
            if feature.stable_feature_id == anchor.stable_feature_svg_id
        ]
    else:
        matches = [
            feature
            for feature in catalog
            if anchor.biological_feature_id
            in {feature.biological_feature_id, feature.stable_feature_id}
        ]
    if (
        len(matches) == 1
        and anchor.stable_feature_svg_id is not None
        and matches[0].stable_feature_id != anchor.stable_feature_svg_id
    ):
        matches = []
    if len(matches) != 1:
        raise ValidationError(
            "Similarity alignment anchor must resolve to exactly one source feature "
            f"for record {provenance.record_key!r}."
        )
    return matches[0]


def _selector_text(selector: FeatureSelector) -> str:
    return selector.value if selector.key is None else f"{selector.key}={selector.value}"


def feature_drawn_start(feature: SourceFeatureIdentity, transform: RecordDisplayTransform) -> float:
    """The display offset of a source feature's first base as drawn.

    The one feature-start rule of ``anchor_feature`` and ``origin_feature``,
    and the Web's Rotate to feature "Start of the record" (``resolveFeatureAnchor``).
    The feature is read along its parts in biological order (the order of
    ``location_parts``); the first base as drawn is where that path begins,
    its 5' base, when the feature reads forward as drawn, and where it ends,
    its 3' base, otherwise. An unstranded or mixed-strand feature is read along
    its parts as listed, as a plus-strand one. The source base is then mapped
    through the transform, so a feature across the origin or a display start
    begins on the piece its path begins on.
    """

    parts = feature.location_parts
    minus = all(strand == -1 for _start, _end, strand in parts)
    five_prime, three_prime = (
        (parts[0][1], parts[-1][0] + 1) if minus else (parts[0][0] + 1, parts[-1][1])
    )
    reads_forward = minus == (transform.source_step == -1)
    return float(transform.source_base_to_display_index(five_prime if reads_forward else three_prime))


_Start = TypeVar("_Start", int, float)


def _match_feature_starts(
    collection: ResolvedRecordCollection,
    selector: FeatureSelector,
    option: str,
    kept: str,
    start_of: Callable[[SourceFeatureIdentity, ResolvedRecordProvenance, RecordDisplayTransform, int], _Start | None],
) -> tuple[_Start | None, ...]:
    """Where ``selector``'s feature starts in each record, by ``start_of``.

    The one owner of how a feature-start option (``anchor_feature``,
    ``origin_feature``) matches features. The selector matches a record's
    source features exactly (:class:`FeatureSelector`); ``start_of`` gives a
    match's start, or ``None`` when it starts outside the drawn region. Matches
    with the same start (a gene and its CDS) are one match, and matches at
    different starts are refused with the candidates. A record without a
    match has ``None`` and keeps its ``kept`` with a WARNING; a selector that
    no record matches is refused. A refusal ``start_of`` raises names the
    option and the record.
    """

    text = _selector_text(selector)
    starts: list[_Start | None] = []
    for index, (record, provenance, transform) in enumerate(zip(
        collection.records, collection.provenance, collection.transforms, strict=True,
    )):
        by_start: dict[_Start, list[SourceFeatureIdentity]] = {}
        outside = 0
        name = f"record #{index + 1} {str(record.id)!r}"
        for feature in provenance.source_feature_catalog:
            if not feature.matches(key=selector.key, value=selector.value, record_id=provenance.source_record_id):
                continue
            try:
                start = start_of(feature, provenance, transform, len(record))
            except ValidationError as exc:
                raise ValidationError(
                    f"{option} {text} in {name}: {exc}", diagnostic={"code": "INPUT_INVALID"},
                ) from exc
            if start is None:
                outside += 1
            else:
                by_start.setdefault(start, []).append(feature)
        if len(by_start) > 1:
            candidates = "; ".join(
                f"{feature.feature_type} {_feature_location_text(feature)}"
                for matches in by_start.values()
                for feature in matches
            )
            raise ValidationError(
                f"{option} {text} matches {sum(map(len, by_start.values()))} features at "
                f"different starts in {name}: {candidates}. Use a selector that names one "
                "feature, such as locus_tag=<tag>.",
                diagnostic={"code": "INPUT_INVALID", "field": option, "reason": "SELECTOR_AMBIGUOUS"},
            )
        if not by_start:
            logger.warning(
                "%s %s %s; the record keeps its %s.",
                option,
                text,
                f"matched no feature in {name}"
                if not outside
                else f"matched {outside} feature(s) in {name} whose start lies outside its drawn region",
                kept,
            )
        starts.append(next(iter(by_start), None))
    if all(start is None for start in starts):
        raise ValidationError(
            f"{option} {text} matched no feature in any record.",
            diagnostic={"code": "INPUT_INVALID", "field": option, "reason": "SELECTOR_NO_MATCH"},
        )
    return tuple(starts)


def _is_cropped(record: SeqRecord, provenance: ResolvedRecordProvenance, transform: RecordDisplayTransform) -> bool:
    """Whether a record draws part of its source: a region, a collection region or a crop."""

    return (
        len(record) != transform.length
        or provenance.is_cropped
        or provenance.region is not None
        or provenance.has_collection_region
    )


def _drawn_start_source_base(transform: RecordDisplayTransform, drawn_start: float) -> int:
    """The 1-based source base a complete circular record starts at so that ``drawn_start`` is its left end."""

    return transform.display_index_to_source_base(int(drawn_start))


def resolve_anchor_feature(
    collection: ResolvedRecordCollection,
    selector: FeatureSelector,
    rows_by_record: Sequence[int],
) -> tuple[ResolvedRecordCollection, FeatureAnchorPlacement]:
    """``collection`` with the records ``anchor_feature`` rotates, and where it puts each record.

    The one owner of what ``anchor_feature`` does to each record and of what it cannot align. Matching and its
    ambiguity rule are :func:`_match_feature_starts`; the start is :func:`feature_drawn_start` in the record's
    drawn transform, so it is read after any crop, reverse complement or display start, and only matches that start
    inside the drawn region count. A complete circular record with the feature is rotated: its display start becomes
    the source base of that start (replacing any start it had), so the feature starts at its left end. Any other
    record with the feature (linear, or cropped) is aligned: the build moves it sideways. The reference is the first
    rotated record, else the first record with the feature. Records that share a row (``rows_by_record``), which one
    x per record cannot align, are refused.
    """

    if len(set(rows_by_record)) != len(rows_by_record):
        raise ValidationError(
            f"anchor_feature {_selector_text(selector)} cannot align records that share a row; "
            "give each record its own row.",
            diagnostic={"code": "INPUT_INVALID", "field": "anchor_feature", "reason": "ANCHOR_SHARED_ROW"},
        )

    def drawn_start(
        feature: SourceFeatureIdentity,
        _provenance: ResolvedRecordProvenance,
        transform: RecordDisplayTransform,
        drawn_length: int,
    ) -> float | None:
        start = feature_drawn_start(feature, transform)
        return start if 0.0 <= start < float(drawn_length) else None

    matched = _match_feature_starts(collection, selector, "anchor_feature", "place", drawn_start)
    rotates = tuple(
        start is not None and transform.is_circular and not _is_cropped(record, provenance, transform)
        for start, record, provenance, transform in zip(
            matched, collection.records, collection.provenance, collection.transforms, strict=True,
        )
    )
    reference = next(
        (index for index, rotated in enumerate(rotates) if rotated),
        next(index for index, start in enumerate(matched) if start is not None),
    )
    statuses: tuple[FeatureAnchorStatus, ...] = tuple(
        "rotated" if rotated
        else "unmatched" if start is None
        else "reference" if index == reference
        else "aligned"
        for index, (start, rotated) in enumerate(zip(matched, rotates, strict=True))
    )
    placement = FeatureAnchorPlacement(
        starts=tuple(0.0 if rotated else start for start, rotated in zip(matched, rotates, strict=True)),
        statuses=statuses,
        reference_index=reference,
    )
    if not placement.rotates:
        return collection, placement
    return ResolvedRecordCollection(
        collection.records,
        tuple(
            replace(
                item,
                display=replace(
                    item.display, start_coordinate=_drawn_start_source_base(transform, cast(float, start)),
                ),
            )
            if rotated
            else item
            for item, transform, start, rotated in zip(
                collection.provenance, collection.transforms, matched, rotates, strict=True,
            )
        ),
    ), placement


# The features without a start as drawn, by the Web capability code
# that refuses each (``resolveFeatureAnchor``, the twin of
# :func:`resolve_origin_feature_start`), in the order the Web checks them.
FEATURE_START_REFUSALS: Mapping[str, str] = {
    "location-parts-overlap": "its parts overlap",
    "location-fuzzy": "it has inexact (fuzzy) coordinates",
    "location-order-ambiguous": "its parts are joined with order()",
    "location-operator-unknown": "its parts are joined in an unknown way",
    "location-traversal-ambiguous": "its parts lie on both strands or in no known order",
    "feature-start-direction-unavailable": "it has no strand",
}


def _feature_start_refusal(feature: SourceFeatureIdentity, profile: _SourceFeatureAnchorProfile) -> str | None:
    bounds = sorted(part[:2] for part in feature.location_parts)
    if any(start < previous_end for (_start, previous_end), (start, _end) in zip(bounds, bounds[1:])):
        return "location-parts-overlap"
    if profile.precision != "exact":
        return "location-fuzzy"
    if profile.operator == "order":
        return "location-order-ambiguous"
    if profile.operator not in {"single", "join"}:
        return "location-operator-unknown"
    if profile.part_order == "ambiguous" or profile.strand == "mixed":
        return "location-traversal-ambiguous"
    if profile.strand not in {"+", "-"}:
        return "feature-start-direction-unavailable"
    return None


def resolve_origin_feature_start(
    feature: SourceFeatureIdentity,
    transform: RecordDisplayTransform,
    profile: _SourceFeatureAnchorProfile,
) -> int:
    """The 1-based source base at which a complete circular record starts so it opens with ``feature``.

    It is the feature's start as drawn (:func:`feature_drawn_start` on the
    uncropped record), converted back to
    the source base ``RecordDisplayOptions.start_coordinate`` takes: the Web
    Rotate to feature "Start of the record" (``resolveFeatureAnchor``, the twin
    in tools/twins.json; tests/fixtures/feature-start-anchor-vectors.json).
    ``profile`` is the feature's source-location profile, the one the Web
    feature metadata carries. The first base needs exact, non-overlapping
    parts on one strand, joined in a known order; any other feature is
    refused, as on the Web.
    """

    if not transform.is_circular:
        raise ValidationError(
            "Record rotation requires an effectively circular record.",
            diagnostic={"code": "INPUT_INVALID"},
        )
    refusal = _feature_start_refusal(feature, profile)
    if refusal is not None:
        raise ValidationError(
            f"{feature.feature_type} {_feature_location_text(feature)} has no start as drawn: "
            f"{FEATURE_START_REFUSALS[refusal]}.",
            diagnostic={"code": "INPUT_INVALID"},
        )
    return _drawn_start_source_base(transform, feature_drawn_start(feature, transform))


def project_origin_feature_starts(
    collection: ResolvedRecordCollection,
    selector: FeatureSelector,
) -> tuple[int | None, ...]:
    """The display start of each record for ``origin_feature``, a 1-based source base or ``None``.

    The one owner of the ``origin_feature`` starts and of the records it
    cannot rotate: a record that is cropped, not circular, or already given a
    display start is refused. Matching and its ambiguity rule are
    :func:`_match_feature_starts`, shared with ``anchor_feature``; the start
    is :func:`resolve_origin_feature_start`.
    """

    text = _selector_text(selector)
    for index, (record, provenance, transform) in enumerate(zip(
        collection.records, collection.provenance, collection.transforms, strict=True,
    )):
        name = f"record #{index + 1} {str(record.id)!r}"
        refusal = (
            f"cannot be combined with a display start (the display start of {name})"
            if provenance.display.start_coordinate is not None
            else f"cannot rotate the cropped {name}" if _is_cropped(record, provenance, transform)
            else f"cannot rotate the linear {name}; give it a circular topology"
            if not transform.is_circular
            else None
        )
        if refusal is not None:
            raise ValidationError(f"origin_feature {text} {refusal}.", diagnostic={"code": "INPUT_INVALID"})

    def origin_start(
        feature: SourceFeatureIdentity,
        provenance: ResolvedRecordProvenance,
        transform: RecordDisplayTransform,
        _drawn_length: int,
    ) -> int:
        return resolve_origin_feature_start(
            feature, transform, provenance.source_features.anchor_profile(feature.source_feature_index)
        )

    return _match_feature_starts(collection, selector, "origin_feature", "start", origin_start)


def with_origin_feature_starts(
    collection: ResolvedRecordCollection,
    selector: FeatureSelector | None,
) -> ResolvedRecordCollection:
    """``collection`` with each matched record's display start set for ``origin_feature``.

    Returns ``collection`` itself when ``selector`` is ``None``.
    """

    if selector is None:
        return collection
    starts = project_origin_feature_starts(collection, selector)
    return ResolvedRecordCollection(
        collection.records,
        tuple(
            item if start is None else replace(item, display=replace(item.display, start_coordinate=start))
            for item, start in zip(collection.provenance, starts, strict=True)
        ),
    )


def _feature_location_text(feature: SourceFeatureIdentity) -> str:
    """The feature's span in 1-based source coordinates, with ``complement`` for the minus strand."""

    parts = ",".join(f"{start + 1}..{end}" for start, end, _strand in feature.location_parts)
    text = parts if len(feature.location_parts) == 1 else f"join({parts})"
    return f"complement({text})" if feature.location_parts[0][2] == -1 else text


@dataclass(frozen=True)
class AlignmentAnchorDisplayFact:
    """Source-bound anchor geometry in one already resolved display transform.

    ``unit_index`` is the drawn unit that holds the anchor (``alignment_anchor_unit``), whose placement draws it.
    """

    source_start: int
    source_end: int
    source_strand: int | None
    displayed_strand: int | None
    display_center: float | None
    unit_index: int


def project_similarity_alignment_anchor_fact(
    collection: ResolvedRecordCollection,
    anchor: AlignmentAnchorIdentity,
) -> AlignmentAnchorDisplayFact:
    """Use the same source identity and center projection as the renderer."""

    # An exact source anchor outside every window of its record is unusable (centre None).
    index, center = alignment_anchor_unit(collection, anchor)
    feature = _alignment_source_feature(anchor, collection.provenance[index])
    if anchor.biological_feature_id not in {feature.biological_feature_id, feature.stable_feature_id}:
        raise ValidationError("Alignment anchor biological identity changed.")
    strands = {part[2] for part in feature.location_parts}
    strand = next(iter(strands)) if len(strands) == 1 and strands <= {-1, 1} else None
    return AlignmentAnchorDisplayFact(
        min(part[0] for part in feature.location_parts),
        max(part[1] for part in feature.location_parts),
        strand, None if strand is None else strand * collection.transforms[index].source_step, center, index,
    )

def similarity_alignment_evidence_edges(
    collection: ResolvedRecordCollection,
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None,
    group_id: str,
) -> tuple[AlignmentEvidenceEdge, ...]:
    """Bind one group's direct ortholog edges to the current source anchors.

    An endpoint is bound only when its record index names the member's source
    record and its feature resolves uniquely there. Edges with an unbound
    endpoint are omitted, so missing or stale graph evidence is never guessed.
    """

    if orthogroups is None:
        return ()
    identity_by_protein: dict[str, tuple[int, AlignmentAnchorIdentity]] = {}
    for member in orthogroups.orthogroups.get(group_id, ()):
        if not 0 <= member.record_index < len(collection.records):
            continue
        provenance = collection.provenance[member.record_index]
        if member.record_id and member.record_id != provenance.source_record_id:
            continue
        candidate = _similarity_candidate_from_orthogroup_member(collection, member)
        if candidate.identity_is_unique:
            identity_by_protein[member.protein_id] = (member.record_index, candidate.anchor)
    edges = []
    for edge in orthogroups.ortholog_edges_by_orthogroup_id.get(group_id, ()):
        query = identity_by_protein.get(edge.query_protein_id)
        subject = identity_by_protein.get(edge.subject_protein_id)
        if query is None or subject is None or (query[0], subject[0]) != (
            edge.query_record_index, edge.subject_record_index,
        ):
            continue
        edges.append(AlignmentEvidenceEdge(
            group_id=group_id, query=query[1], subject=subject[1], edge_kind=str(edge.edge_kind),
        ))
    return tuple(edges)


def _similarity_candidate_from_orthogroup_member(
    collection: ResolvedRecordCollection,
    member: OrthogroupMember,
) -> SimilarityAlignmentCandidate:
    provenance = collection.provenance[member.record_index]
    catalog = provenance.source_feature_catalog
    matches = [
        feature
        for feature in catalog
        if member.feature_svg_id is not None
        and feature.stable_feature_id == member.feature_svg_id
    ]
    if not matches:
        matches = [
            feature
            for feature in catalog
            if feature.source_feature_index == member.feature_index
        ]
    identity_is_unique = len(matches) == 1
    feature = matches[0] if identity_is_unique else None
    anchor = AlignmentAnchorIdentity(
        # A plan names a record drawn in several windows by the record's key (D-50).
        record_key=provenance.record_identity_key,
        biological_feature_id=(
            feature.biological_feature_id
            if feature is not None
            else str(member.feature_svg_id or member.protein_id)
        ),
        source_feature_index=(
            feature.source_feature_index if feature is not None else member.feature_index
        ),
        stable_feature_svg_id=(
            feature.stable_feature_id
            if feature is not None
            else member.feature_svg_id
        ),
    )
    display_center: float | None = None
    if feature is not None:
        projected = _source_feature_display_center(feature, collection.transforms[member.record_index])
        if 0.0 <= projected <= float(len(collection.records[member.record_index])):
            display_center = projected
    return SimilarityAlignmentCandidate(
        group_id=member.orthogroup_id,
        anchor=anchor,
        displayed_strand=member.strand if member.strand in (-1, 1) else None,
        center_mappable=display_center is not None,
        display_center=display_center,
        identity_is_unique=identity_is_unique,
        representative=member.representative,
        role=str(member.role),
        source_start=member.start,
        source_end=member.end,
        display_name=member.gene or member.product or member.label or member.protein_id,
    )


@dataclass(frozen=True)
class RecordInputManifest:
    """Unresolved typed record inputs derived from CLI path syntax."""

    records: tuple[RecordInput, ...]
    record_options: RecordCollectionOptions
    source_paths: tuple[str, ...]
    # Records-table ``losat_gencode`` per record input; empty without values.
    losat_gencodes: tuple[int | None, ...] = ()


def _expanded_cli_track_values(
    values: Sequence[object] | None,
    *,
    track_count: int,
    field_name: str,
) -> list[object | None]:
    items = list(values or ())
    if not items:
        return [None] * track_count
    if len(items) == 1:
        return items * track_count
    if len(items) != track_count:
        raise ValidationError(
            f"{field_name} count must be one or equal to the number of "
            f"depth tracks ({track_count})."
        )
    return items


def _optional_positive_float(value: object | None, *, field_name: str) -> float | None:
    if value is None or str(value).strip().lower() in {
        "",
        "auto",
        "none",
        "null",
        "-",
    }:
        return None
    try:
        numeric = ascii_float(cast(Any, value))
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{field_name} values must be numeric or auto.") from exc
    if numeric <= 0:
        raise ValidationError(f"{field_name} values must be > 0.")
    return numeric


def depth_track_inputs_from_cli(
    groups: Sequence[Sequence[str]] | None,
    *,
    labels: Sequence[str] | None = None,
    colors: Sequence[str] | None = None,
    heights: Sequence[object] | None = None,
    large_tick_intervals: Sequence[object] | None = None,
    small_tick_intervals: Sequence[object] | None = None,
    tick_font_sizes: Sequence[object] | None = None,
) -> tuple[DepthTrackInput, ...] | None:
    """Translate logical CLI depth groups without knowing record cardinality."""

    if not groups:
        return None
    track_count = len(groups)
    label_values = _expanded_cli_track_values(
        labels,
        track_count=track_count,
        field_name="depth_track_labels",
    )
    color_values = _expanded_cli_track_values(
        colors,
        track_count=track_count,
        field_name="depth_track_colors",
    )
    height_values = _expanded_cli_track_values(
        heights,
        track_count=track_count,
        field_name="depth_track_heights",
    )
    large_values = _expanded_cli_track_values(
        large_tick_intervals,
        track_count=track_count,
        field_name="depth_track_large_tick_intervals",
    )
    small_values = _expanded_cli_track_values(
        small_tick_intervals,
        track_count=track_count,
        field_name="depth_track_small_tick_intervals",
    )
    font_values = _expanded_cli_track_values(
        tick_font_sizes,
        track_count=track_count,
        field_name="depth_track_tick_font_sizes",
    )
    inputs: list[DepthTrackInput] = []
    for index, group in enumerate(groups):
        sources = tuple(
            None
            if str(path).strip().lower() in {"", "-", "none", "null"}
            else str(path)
            for path in group
        )
        if not sources or not any(source is not None for source in sources):
            raise ValidationError(
                f"--depth_track #{index + 1} must include at least one file."
            )
        source: object = sources[0] if len(sources) == 1 else sources
        inputs.append(
            DepthTrackInput(
                source=source,
                label=(
                    str(label_values[index]).strip()
                    if label_values[index] is not None
                    else None
                ),
                color=(
                    str(color_values[index]).strip()
                    if color_values[index] is not None
                    else None
                ),
                height=_optional_positive_float(
                    height_values[index],
                    field_name="depth_track_heights",
                ),
                large_tick_interval=_optional_positive_float(
                    large_values[index],
                    field_name="depth_track_large_tick_intervals",
                ),
                small_tick_interval=_optional_positive_float(
                    small_values[index],
                    field_name="depth_track_small_tick_intervals",
                ),
                tick_font_size=_optional_positive_float(
                    font_values[index],
                    field_name="depth_track_tick_font_sizes",
                ),
            )
        )
    return tuple(inputs)


def _source_key(source: RecordInputSource) -> tuple[object, ...]:
    if isinstance(source, InMemoryRecordSource):
        return ("memory", id(source.record), id(source.source_features))
    kind, paths = _source_details(source)
    return (kind, *paths)


def _source_details(
    source: RecordInputSource,
) -> tuple[RecordSourceKind, tuple[str, ...]]:
    """The one reading of a record source: its kind and the paths of its files, in resource order."""
    if isinstance(source, GenBankInputSource):
        return "genbank", (str(source.path),)
    if isinstance(source, GffFastaInputSource):
        return "gff_fasta", (str(source.gff_path), str(source.fasta_path))
    if isinstance(source, FastaInputSource):
        return "fasta", (str(source.fasta_path),)
    if isinstance(source, InMemoryRecordSource):
        return "memory", ()
    raise ValidationError("Unsupported record input source.")


def _file_source(kind: RecordSourceKind, paths: Sequence[str]) -> RecordInputSource:
    """The file source :func:`_source_details` reads as ``kind`` and ``paths``: a record held from a file plan
    (``held_record_origin``) names the file it was read from this way."""
    if kind == "genbank":
        return GenBankInputSource(paths[0])
    if kind == "gff_fasta":
        return GffFastaInputSource(paths[0], paths[1])
    if kind == "fasta":
        return FastaInputSource(paths[0])
    raise ValidationError(f"A {kind} record source names no file.", diagnostic={"code": "INPUT_INVALID"})


def _gff_feature_types(
    gff_candidate_features: Sequence[str] | None, gff_keep_all_features: bool,
) -> frozenset[str] | None:
    """The GFF3 feature types a request draws; ``None`` keeps every type."""
    return None if gff_keep_all_features or gff_candidate_features is None else frozenset(gff_candidate_features)


def _load_source_records(
    source: RecordInputSource,
    *,
    genbank_loader: GenBankLoader,
    gff_loader: GffFastaLoader,
    fasta_loader: FastaLoader = load_fasta,
) -> ParsedRecordInputs:
    def load() -> ParsedRecordInputs:
        gff_catalogs: list[SourceFeatureCatalog] = []
        if isinstance(source, GenBankInputSource):
            records = genbank_loader([str(source.path)])
        elif isinstance(source, FastaInputSource):
            records = fasta_loader([str(source.fasta_path)])
        elif isinstance(source, GffFastaInputSource):
            # One parse serves every type filter: ``resolve_record_inputs`` applies it in memory.
            records = gff_loader(
                str(source.gff_path),
                str(source.fasta_path),
                source_feature_catalogs=gff_catalogs,
            )
        elif isinstance(source, InMemoryRecordSource):
            # The caller keeps its record: a plan must not change with it.
            records = [in_draw_order(copy.deepcopy(source.record))]
        else:  # pragma: no cover - RecordInput validates this union.
            raise ValidationError("Unsupported record input source.")
        if not records:
            raise ValidationError("A record input source resolved to no records.")
        catalogs = (
            (source.source_features,)
            if isinstance(source, InMemoryRecordSource) and source.source_features is not None
            else tuple(gff_catalogs) if isinstance(source, GffFastaInputSource)
            else tuple(capture_source_features(record) for record in records)
        )
        return ParsedRecordInputs(tuple(records), catalogs)

    cache_spec = _prepared_source_cache_spec(source)
    if cache_spec is None:
        return load()
    key, identities = cache_spec
    return get_or_build_parsed_source(
        key,
        identities,
        load,
        publish=lambda value: bool(value),
    )


class _NotPrepared(Exception):
    """A parsed-source lookup that found no entry: it builds nothing."""


def prepared_parsed_source(
    source: RecordInputSource,
    load: Callable[[RecordInputSource], ParsedRecordInputs],
) -> ParsedRecordInputs | None:
    """The records a Generate parses from this staged Web source, or ``None``.

    A Worker helper inside ``prepared_input_lookup`` reads the entry stored
    under the same resource identities. On a miss it parses with ``load``, the
    render's own source load (``request_render.load_render_source_records``,
    which the caller passes): the protein extraction's lookup publishes that
    entry for the next Generate (OV-1595), any other lookup keeps it for the
    call. It gets ``None`` for a source the Worker did not stage, or when that
    load fails: it then keeps its own reader, so its errors are the ones it
    has without the cache. The records are shared: a reader must not change
    them (the lookup checks their owner stamps).
    """

    cache_spec = _prepared_source_cache_spec(source)
    if cache_spec is None:
        return None

    def not_prepared() -> ParsedRecordInputs:
        raise _NotPrepared

    try:
        return cast(ParsedRecordInputs, get_or_build_parsed_source(*cache_spec, not_prepared))
    except _NotPrepared:
        pass
    try:
        return load(source)
    except Exception:
        return None


def _prepared_source_cache_spec(
    source: RecordInputSource,
) -> tuple[Hashable, frozenset[PreparedResourceIdentity]] | None:
    """Build a path-independent parsed-source cache key for a Web resource."""

    if isinstance(source, InMemoryRecordSource):
        return None
    kind, paths = _source_details(source)
    identities = tuple(prepared_resource_identity(path) for path in paths)
    if any(identity is None for identity in identities):
        return None
    return ("parsed-source-v1", kind, *identities), frozenset(
        identity for identity in identities if identity is not None
    )


def _selected_source_indexes(
    records: Sequence[SeqRecord],
    selector: RecordSelector | None,
) -> list[int]:
    if selector is None:
        return list(range(len(records)))
    if selector.record_index is not None:
        index = selector.record_index
        if index < 0 or index >= len(records):
            raise ValidationError(
                f"Record selector {selector.label()} is out of range "
                f"(loaded {len(records)} record(s))."
            )
        return [index]
    record_id = selector.record_id or ""
    matches = [
        index for index, record in enumerate(records) if record.id == record_id
    ]
    if not matches:
        raise ValidationError(
            f"Record selector '{record_id}' did not match any record ID."
        )
    if len(matches) > 1:
        raise ValidationError(
            f"Record selector '{record_id}' matched multiple records. "
            "Use #index to disambiguate."
        )
    return matches


def _cardinality_indexes(
    indexes: Sequence[int],
    *,
    cardinality: RecordCardinality,
    input_index: int,
) -> tuple[int, ...]:
    selected = tuple(indexes)
    if cardinality is RecordCardinality.EXACTLY_ONE:
        if len(selected) != 1:
            raise ValidationError(
                f"RecordInput #{input_index + 1} requires exactly one record; "
                f"resolved {len(selected)}. Add a selector, use "
                "RecordCardinality.FIRST, or explicitly use RecordCardinality.ALL. "
                "For per-record CLI settings, use --records_table."
            )
        return selected
    if not selected:
        raise ValidationError(
            f"RecordInput #{input_index + 1} resolved no records."
        )
    if cardinality is RecordCardinality.FIRST:
        return selected[:1]
    if cardinality is RecordCardinality.ALL:
        return selected
    raise ValidationError("Unsupported record cardinality.")


def _apply_record_key(record: SeqRecord, *, record_key: str) -> None:
    if getattr(record, "annotations", None) is None:
        record.annotations = {}
    record.annotations["gbdraw_record_key"] = record_key


def effective_record_presentation(
    presentation: RecordPresentation,
    options: RecordCollectionOptions,
    resolved_index: int,
) -> RecordPresentation:
    """The presentation a resolved record is drawn with.

    A non-empty collection label or subtitle at the record's place replaces the
    one of its input. The record itself carries neither, so a record is the
    same object, with the same render-stage identity, whatever its label.
    """

    label = options.labels[resolved_index] if resolved_index < len(options.labels) else None
    subtitle = options.subtitles[resolved_index] if resolved_index < len(options.subtitles) else None
    if not label and not subtitle:
        return presentation
    return replace(
        presentation,
        label=label or presentation.label,
        subtitle=subtitle or presentation.subtitle,
    )


def unit_presentation(
    presentation: RecordPresentation,
    options: RecordCollectionOptions,
    item: ResolvedRecordProvenance,
) -> RecordPresentation:
    """The presentation a drawn unit is drawn with: its record's, placed by its window's column (D-08)."""

    effective = effective_record_presentation(presentation, options, item.place)
    return replace(effective, grid_column=unit_grid_column(effective.grid_column, item.window))


def unit_grid_column(record_column: int | None, window: ResolvedWindow | None) -> int | None:
    """The column a drawn unit is placed by: its window's own column, else its record's (D-08, D-14)."""

    column = window.span.column if window is not None else None
    return record_column if column is None else column


def _apply_provenance_annotations(
    record: SeqRecord,
    provenance: ResolvedRecordProvenance,
) -> None:
    if getattr(record, "annotations", None) is None:
        record.annotations = {}
    # The file the record was read from: a record given in memory names it only while its origin holds
    # (``held_record_origin``), so an edited record names no file.
    record.annotations.update(
        {
            "gbdraw_input_index": provenance.input_index,
            "gbdraw_source_record_index": provenance.source_record_index,
            "gbdraw_source_record_count": provenance.source_record_count,
            "gbdraw_source_record_id": provenance.source_record_id,
            "gbdraw_source_record_id_unique": provenance.source_record_id_unique,
            "gbdraw_record_key": provenance.record_key,
            "gbdraw_record_cardinality": provenance.cardinality.value,
            "gbdraw_source_kind": provenance.source_kind,
            "gbdraw_source_paths": provenance.source_paths,
        }
    )
    if provenance.source_length is not None:
        record.annotations["gbdraw_source_length"] = provenance.source_length
    # The window this unit draws and the record it is a window of (D-01): this plan's window, or an earlier plan's
    # while the record is unchanged since (``origin_window_of``). An edited unit given again is no window, so
    # per-record readers (``record_identity_key``, Linear rows) and the Session writer agree.
    if provenance.window is not None:
        record.annotations[WINDOW_SPAN_ANNOTATION] = (provenance.window.span.start, provenance.window.span.end)
    elif provenance.origin_window_of is None:
        record.annotations.pop(WINDOW_SPAN_ANNOTATION, None)
    if provenance.window_of is not None:
        record.annotations[WINDOW_OF_ANNOTATION] = provenance.window_of
    else:
        record.annotations.pop(WINDOW_OF_ANNOTATION, None)
    # The file stamps stay when this plan has no paths: the file loader writes them on records it reads
    # (gbdraw/io/genome.py), and file-qualified region selectors of an in-memory record read them.
    if provenance.source_paths:
        record.annotations["gbdraw_source_file"] = provenance.source_paths[0]
        record.annotations["gbdraw_source_basename"] = Path(provenance.source_paths[0]).name


def numbered_places(places: Sequence[InputPlace]) -> tuple[InputPlace, ...]:
    """``places`` after each placed input that draws several records whole (a File card, ``RecordCardinality.ALL``)
    takes one column per record: consecutive columns from its own, in file order (``InputPlace.record_column``), with
    every later column of its row moved after them (OV-3072). Places without such an input stay as given.

    The one owner of that numbering: the record resolver places the drawn records by it, and the Web's card plan
    reads its rows by it (``plan_card_windows`` of ``gbdraw.web_support.linear_units``).
    """

    widened = [(place.row, place.column, place.records - 1) for place in places if place.column is not None and place.records > 1]
    if not widened:
        return tuple(places)

    def moved(row: int | None, column: int | None) -> int | None:
        return None if column is None else column + sum(
            extra for widened_row, widened_column, extra in widened if widened_row == row and widened_column < column
        )

    return tuple(
        place._replace(
            column=moved(place.row, place.column),
            window_columns=tuple(moved(place.row, column) for column in place.window_columns),
        )
        for place in places
    )


def _placed_input(record_input: RecordInput, place: InputPlace) -> tuple[RecordPresentation, tuple[RecordWindow, ...]]:
    """The presentation and windows of ``record_input`` at ``place`` (``numbered_places``): the input's own while its
    place is unchanged."""

    presentation, windows = record_input.presentation, tuple(record_input.windows)
    if place.window_columns and place.window_columns != tuple(window.column for window in windows):
        windows = tuple(replace(window, column=column) for window, column in zip(windows, place.window_columns, strict=True))
    if place.column != presentation.grid_column:
        presentation = replace(presentation, grid_column=place.column)
    return presentation, windows


@dataclass
class _PendingRecord:
    """One resolved record before its windows expand: its drawn records so far and how to finish them."""

    units: tuple[SeqRecord, ...]
    provenance: ResolvedRecordProvenance
    windows: tuple[RecordWindow, ...]
    pending_reverse: bool
    # Whether ``units`` were cut from the record by its windows in this plan (each keeps its cut's digest).
    cut: bool = False


def _expand_record_windows(
    pending: Sequence[_PendingRecord],
    options: RecordCollectionOptions,
) -> ResolvedRecordCollection:
    """Crop by collection regions, orient, and give each window of a record its own drawn unit (D-01).

    A record's presentation reverse complement applies after its crop, the
    order a region's own ``:rc`` uses, so both spellings draw one window.
    Several collection regions that select one record are its windows. A unit
    keeps the record's provenance with its own key (``window_unit_key``), its
    window, and the record's place, by which labels and subtitles apply.
    """

    try:
        collection_specs = (
            assign_region_specs([item.units[0] for item in pending], options.regions) if options.regions else {}
        )
    except ValueError as exc:
        raise ValidationError(str(exc)) from exc
    # Labels and subtitles apply by place through effective_record_presentation.
    for field_name in ("labels", "subtitles"):
        if len(getattr(options, field_name)) > len(pending):
            logger.warning(
                "More record %s were provided than records resolved; "
                "extra values will be ignored.",
                field_name,
            )
    records: list[SeqRecord] = []
    provenance: list[ResolvedRecordProvenance] = []
    for position, item in enumerate(pending):
        specs = collection_specs.get(position, ())
        windows = item.windows
        units = item.units
        cut = item.cut
        if specs:
            if windows:
                raise ValidationError(
                    f"Record {item.provenance.record_key!r} has windows in its input and a collection region; "
                    "give its windows in one place.",
                    diagnostic={"code": "INPUT_INVALID"},
                )
            windows = tuple(RecordWindow.from_region(spec) for spec in specs)
            try:
                units = crop_record_windows(
                    units[0],
                    specs,
                    item.pending_reverse,
                    log=logger,
                    circular=drawn_circular(item.provenance.display.is_circular, item.provenance.detected_topology),
                )
            except ValueError as exc:
                raise ValidationError(str(exc), diagnostic={"code": "INPUT_INVALID"}) from exc
            cut = True
        elif item.pending_reverse:
            units = crop_record_windows(units[0], (), True, log=logger)
        # The crop draws one unit per window, as ``record_drawn_windows`` plans them.
        for unit, window in zip(units, record_drawn_windows(item.provenance.record_key, windows), strict=True):
            unit_item = replace(
                item.provenance,
                resolved_index=len(records),
                record_position=position,
                record_key=window.unit_key if window is not None else item.provenance.record_key,
                # A unit of an input window draws that window's span as its region.
                region=window.span.region() if window is not None and not specs else None,
                has_collection_region=bool(specs),
                window=window,
            )
            unit_item = replace(
                unit_item, presentation=unit_presentation(item.provenance.presentation, options, unit_item)
            )
            _apply_record_key(unit, record_key=unit_item.record_key)
            _apply_provenance_annotations(unit, unit_item)
            # Since its cut, a unit changed only in gbdraw's annotations, which its digest leaves out.
            _keep_record_origin(unit, unit_item, cut_content(unit) if cut else None)
            records.append(unit)
            provenance.append(unit_item)
    return ResolvedRecordCollection(tuple(records), tuple(provenance))


def reverse_resolved_records(
    collection: ResolvedRecordCollection,
    indexes: frozenset[int],
) -> ResolvedRecordCollection:
    """Reverse-complement the resolved records at ``indexes`` after their crop.

    The order :func:`crop_then_reverse_records` applies to a record's own
    reverse flag, so each record draws as an explicit reverse of it draws.
    """

    records = list(collection.records)
    for index in sorted(indexes):
        (records[index],) = reverse_records((records[index],), True, log=logger)
        _apply_provenance_annotations(records[index], collection.provenance[index])
        # A reversed crop keeps the digest of its reverse complement as cut (``reversed_window_source``).
        _keep_record_origin(records[index], collection.provenance[index], cut_content(records[index]))
    return ResolvedRecordCollection(tuple(records), collection.provenance)


def _keep_record_origin(
    record: SeqRecord, provenance: ResolvedRecordProvenance, content: RecordContent | None = None
) -> None:
    """Keep on ``record``, a planned record now complete, where it was read from (OV-2923).

    A later plan of the record, or a Session writer, reads it while the record is unchanged (``held_record_origin``).
    ``content`` is the record's content when its caller holds it (a crop's, from its cut).
    """

    record.annotations[RECORD_ORIGIN_ANNOTATION] = provenance.origin(record, content)


def records_holding_their_sequences(records: Iterable[SeqRecord]) -> list[SeqRecord]:
    """Planned records handed to a caller: each holds its sequence (``detached_records``), and a crop the record it
    was cut from, so the caller may rewrite or delete the input files (OV-3353).

    The contents a reader compares later (where a record was read from, what a crop was when cut) take their digests
    first, while the files still hold the bases. Crops of one record share one copy of it.
    """

    sources: dict[int, SeqRecord] = {}
    held: list[SeqRecord] = []
    for record in records:
        _take_digests(record)
        [copied] = detached_records([record])
        stored = record.annotations.get(WINDOW_SOURCE_ANNOTATION)
        if isinstance(stored, WindowSource):
            if id(stored.record) not in sources:
                _take_digests(stored.record)
                [sources[id(stored.record)]] = detached_records([stored.record])
            if sources[id(stored.record)] is not stored.record:
                copied = copy.copy(copied)
                copied.annotations = {
                    **record.annotations, WINDOW_SOURCE_ANNOTATION: replace(stored, record=sources[id(stored.record)]),
                }
        held.append(copied)
    return held


def _take_digests(record: SeqRecord) -> None:
    origin = record.annotations.get(RECORD_ORIGIN_ANNOTATION)
    if isinstance(origin, RecordOrigin):
        origin.content.digest
    crop_content = cut_content(record)
    if crop_content is not None:
        crop_content.digest


def _detached_record(record: SeqRecord) -> SeqRecord:
    """Return a record whose own fields can be written without changing ``record``.

    Planning writes the annotations of each resolved record; parsed sources are
    shared by every record resolved from them. The sequence and the feature
    objects stay shared: no later step writes them (reverse complement and
    region crops build new features). A caller's in-memory record is copied
    whole when its source loads.
    """

    detached = copy.copy(record)
    detached.annotations = dict(record.annotations)
    detached.features = list(record.features)
    detached.dbxrefs = list(record.dbxrefs)
    return detached


def _uncropped_length(record: SeqRecord) -> int:
    """The length of a loaded record before the crop its loader applied."""

    length = record.annotations.get("gbdraw_source_length") if record.annotations.get("gbdraw_region_applied") else None
    return length if isinstance(length, int) else len(record)


def resolve_record_inputs(
    record_inputs: Sequence[RecordInput],
    *,
    record_options: RecordCollectionOptions | None = None,
    gff_candidate_features: Sequence[str] | None,
    gff_keep_all_features: bool,
    genbank_loader: GenBankLoader = load_gbks,
    gff_loader: GffFastaLoader = parse_gff_fasta,
    fasta_loader: FastaLoader = load_fasta,
    parsed_sources: MutableMapping[tuple[object, ...], ParsedRecordInputs] | None = None,
) -> ResolvedRecordCollection:
    """Load each unique source once, then apply typed selection and transforms.

    ``parsed_sources`` keeps the parses across calls: a caller that resolves the
    same sources again with other GFF3 candidate types does not parse them again.
    """

    inputs = tuple(record_inputs)
    if not inputs:
        raise ValidationError("A request requires at least one RecordInput.")
    gff_types = _gff_feature_types(gff_candidate_features, gff_keep_all_features)
    cache = {} if parsed_sources is None else parsed_sources
    id_counts: dict[tuple[object, ...], Counter[str]] = {}
    pending: list[_PendingRecord] = []
    selections: list[tuple[tuple[object, ...], ParsedRecordInputs, tuple[int, ...]]] = []
    for input_index, record_input in enumerate(inputs):
        key = _source_key(record_input.source)
        parsed = cache.get(key)
        if parsed is None:
            parsed = cache[key] = _load_source_records(
                record_input.source,
                genbank_loader=genbank_loader,
                gff_loader=gff_loader,
                fasta_loader=fasta_loader,
            )
        if key not in id_counts:
            id_counts[key] = Counter(str(record.id) for record in parsed.records)
        selections.append((key, parsed, _cardinality_indexes(
            _selected_source_indexes(parsed.records, record_input.selector),
            cardinality=record_input.cardinality,
            input_index=input_index,
        )))
    # Every input's count of records is known: a placed File card takes a column per record (OV-3072).
    places = numbered_places([
        InputPlace.of(record_input, len(source_indexes)) for record_input, (*_, source_indexes) in zip(inputs, selections)
    ])
    for input_index, (record_input, (key, parsed, source_indexes), place) in enumerate(zip(inputs, selections, places)):
        raw_records = parsed.records
        file_length = sum(_uncropped_length(record) for record in raw_records)
        selector = record_input.selector
        presentation, input_windows = _placed_input(record_input, place)
        source_kind, source_paths = _source_details(record_input.source)
        # A record planned before and unchanged since is still the record of the file it was read from (OV-2923).
        origin = held_record_origin(record_input.source.record) if isinstance(record_input.source, InMemoryRecordSource) else None
        expands = len(source_indexes) > 1
        for position, source_record_index in enumerate(source_indexes):
            source_record = raw_records[source_record_index]
            require_record_sequence(source_record, source_paths[0] if source_paths else "the record input")
            types_left_out: tuple[tuple[str, int], ...] = ()
            if source_kind == "gff_fasta":
                types_left_out = feature_types_left_out(source_record.features, gff_types)
                source_record = filter_features_by_type(
                    source_record, gff_types, source_indexes=_feature_source_index_map(source_record.features),
                )
            source_cropped = bool(source_record.annotations.get("gbdraw_region_applied"))
            record = _detached_record(source_record)
            reverse = record_input.presentation.reverse_complement
            units: tuple[SeqRecord, ...] = (record,)
            if input_windows:
                try:
                    units = crop_record_windows(
                        record,
                        [window.region() for window in input_windows],
                        reverse,
                        log=logger,
                        circular=drawn_circular(record_input.display.is_circular, detected_topology(source_record, source_kind)),
                    )
                except ValueError as exc:
                    raise ValidationError(str(exc)) from exc
                reverse = False
            base_key = record_input.record_key or f"{CLI_RECORD_KEY_PREFIX}{input_index + 1}"
            record_key = expanded_record_key(base_key, source_record_index + 1) if expands else base_key
            item = ResolvedRecordProvenance(
                resolved_index=len(pending),
                input_index=input_index,
                source_record_index=origin.source_record_index if origin else source_record_index,
                source_record_count=origin.source_record_count if origin else len(raw_records),
                source_record_id=origin.source_record_id if origin else str(raw_records[source_record_index].id),
                source_record_id_unique=(
                    origin.source_record_id_unique if origin
                    else id_counts[key][str(raw_records[source_record_index].id)] == 1
                ),
                source_kind=origin.source_kind if origin else source_kind,
                source_paths=origin.source_paths if origin else source_paths,
                origin_window_of=origin.window_of if origin else None,
                record_key=record_key,
                cardinality=record_input.cardinality,
                selector=selector,
                region=None,
                presentation=(
                    replace(presentation, grid_column=place.record_column(position))
                    if expands and place.column is not None else presentation
                ),
                display=record_input.display,
                source_length=_uncropped_length(source_record),
                source_file_length=origin.source_file_length if origin else file_length,
                detected_topology=detected_topology(source_record, source_kind),
                is_cropped=source_cropped,
                source_features=parsed.source_features[source_record_index],
                feature_types_left_out_at_load=types_left_out,
            )
            for unit in units:
                _apply_record_key(unit, record_key=record_key)
                _apply_provenance_annotations(unit, item)
            pending.append(_PendingRecord(units, item, input_windows, reverse, bool(input_windows)))
    return _expand_record_windows(pending, record_options or RecordCollectionOptions())


def _normalized_cli_values(
    values: Sequence[object] | None,
    *,
    count: int,
    field_name: str,
    default: object,
) -> list[object]:
    normalized = list(values or ())
    if len(normalized) > count:
        raise ValidationError(
            f"Too many {field_name} values (expected at most {count})."
        )
    normalized.extend(default for _ in range(count - len(normalized)))
    return normalized


def parse_cli_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value or "").strip().lower()
    if text in {"1", "true", "yes", "y", "on"}:
        return True
    if text in {"0", "false", "no", "n", "off", "", "none", "null", "-"}:
        return False
    raise ValidationError(f"Invalid reverse_complement value: {value}")


def record_input_manifest_from_paths(
    *,
    gbk_paths: Sequence[str] | None = None,
    gff_paths: Sequence[str] | None = None,
    fasta_paths: Sequence[str] | None = None,
    cardinalities: Sequence[RecordCardinality] | RecordCardinality,
    selectors: Sequence[str] | None = None,
    reverse_flags: Sequence[object] | None = None,
    labels: Sequence[str] | None = None,
    subtitles: Sequence[str] | None = None,
    regions: Sequence[str] | None = None,
) -> RecordInputManifest:
    """Translate CLI path lists into unresolved typed record inputs."""

    gbks = tuple(gbk_paths or ())
    gffs = tuple(gff_paths or ())
    fastas = tuple(fasta_paths or ())
    if gffs and len(gffs) != len(fastas):
        raise ValidationError("GFF3 and FASTA path counts must match.")
    if bool(gbks) == bool(fastas):
        raise ValidationError(
            "Specify either GenBank paths, or FASTA paths with or without GFF3 paths."
        )
    count = len(gbks) if gbks else len(fastas)
    if isinstance(cardinalities, RecordCardinality):
        cardinality_values = [cardinalities] * count
    else:
        cardinality_values = list(cardinalities)
        if len(cardinality_values) != count:
            raise ValidationError(
                "Record cardinality count must match the number of input sources."
            )
    if not all(
        isinstance(cardinality, RecordCardinality)
        for cardinality in cardinality_values
    ):
        raise ValidationError(
            "Record cardinalities must be RecordCardinality values."
        )
    selector_values = _normalized_cli_values(
        selectors,
        count=count,
        field_name="record_id",
        default="",
    )
    reverse_values = _normalized_cli_values(
        reverse_flags,
        count=count,
        field_name="reverse_complement",
        default=False,
    )
    records: list[RecordInput] = []
    source_paths: list[str] = []
    for index in range(count):
        source: RecordInputSource = (
            GenBankInputSource(gbks[index]) if gbks
            else GffFastaInputSource(gffs[index], fastas[index]) if gffs
            else FastaInputSource(fastas[index])
        )
        source_paths.append(_source_details(source)[1][0])
        selector = parse_record_selector(str(selector_values[index] or ""))
        cardinality = cardinality_values[index]
        if selector is not None:
            cardinality = RecordCardinality.EXACTLY_ONE
        records.append(
            RecordInput(
                source=source,
                cardinality=cardinality,
                selector=selector,
                presentation=RecordPresentation(
                    reverse_complement=parse_cli_bool(reverse_values[index]),
                ),
                record_key=f"{CLI_RECORD_KEY_PREFIX}{index + 1}",
            )
        )
    return RecordInputManifest(
        records=tuple(records),
        record_options=RecordCollectionOptions(
            regions=tuple(parse_region_spec(value) for value in (regions or ())),
            labels=tuple(labels or ()),
            subtitles=tuple(subtitles or ()),
        ),
        source_paths=tuple(source_paths),
    )


def _records_table_input(record: RecordsTableRecord, index: int) -> RecordInput:
    """The record input of one records-table record: a row, or the rows that share a ``record_group`` (D-32)."""

    first = record.first
    source: RecordInputSource = (
        GenBankInputSource(first.gbk) if first.gbk
        else GffFastaInputSource(first.gff, first.fasta) if first.gff
        else FastaInputSource(first.fasta)
    )
    windowed = record.is_windowed and len(record.rows) > 1
    windows: tuple[RecordWindow, ...] = ()
    if windowed:
        # One table row per window, placed by its column (ties in table order).
        ordered = sorted(record.rows, key=lambda row: (row.column is None, row.column or 0))
        windows = tuple(
            replace(
                RecordWindow.from_region(spec := parse_region_spec(row.region), column=row.column),
                reverse_complement=spec.reverse_complement or row.reverse_complement,
            )
            for row in ordered
        )
    topology = record.topology
    drawn = RecordInput(
        source=source,
        cardinality=RecordCardinality.EXACTLY_ONE,
        selector=parse_record_selector(first.record_id),
        display=RecordDisplayOptions(
            is_circular=None if topology is None else topology == "circular",
            start_coordinate=first.display_start,
        ),
        # One region is the record's span; a group's windows are its spans (a region passed replaces windows).
        **({"windows": windows} if windowed or not first.region else {"region": parse_region_spec(first.region)}),
        presentation=RecordPresentation(
            label=record.record_label or None,
            subtitle=record.record_subtitle or None,
            reverse_complement=first.reverse_complement and not windowed,
            grid_row=first.row,
            grid_column=None if windowed else first.column,
        ),
        record_key=f"{CLI_RECORD_KEY_PREFIX}{index + 1}",
    )
    # A record drawn in one window carries its orientation on the window, as every request writer writes it.
    return drawn.oriented(drawn.is_reversed) if len(drawn.windows) == 1 else drawn


def record_input_manifest_from_table(path: str) -> RecordInputManifest:
    """Translate one records table into exact-one unresolved record inputs, one per record it draws."""

    table = read_records_table(path)
    records = [_records_table_input(record, index) for index, record in enumerate(table.records)]
    return RecordInputManifest(
        records=tuple(records),
        record_options=RecordCollectionOptions(),
        source_paths=tuple(_source_details(record.source)[1][0] for record in records),
        losat_gencodes=(
            tuple(record.losat_gencode for record in table.records)
            if any(row.losat_gencode is not None for row in table.rows)
            else ()
        ),
    )


def resolve_circular_options(
    options: CircularDiagramOptions,
) -> CircularDiagramOptions:
    """Resolve an optional Circular comparison manifest once."""

    if options.conservation_table_file is None:
        return options
    losat = options.losat_search is not None
    table = read_conservation_table(options.conservation_table_file, losat=losat)
    try:
        return _with_conservation_table(options, table, losat=losat)
    except ValidationError as exc:
        # The typed check names its field; the table's user knows the column (OV-3326).
        if (exc.diagnostic or {}).get("field") != "conservation_losat_gencodes":
            raise
        column = f"{options.conservation_table_file}: losat_gencode column"
        raise ValidationError(
            str(exc).replace("conservation_losat_gencodes", column, 1), diagnostic=exc.diagnostic
        ) from exc


def _with_conservation_table(
    options: CircularDiagramOptions, table: ConservationTable, *, losat: bool
) -> CircularDiagramOptions:
    return replace(
        options,
        conservation_table_file=None,
        conservation_blast_files=None if losat else tuple(table.conservation_blast_files),
        conservation_losat_gencodes=(
            tuple(table.losat_gencodes)
            if losat and table.losat_gencodes is not None
            else None
        ),
        conservation_sequence_files=(
            tuple(table.comparison_sequence_files)
            if table.comparison_sequence_files is not None
            else None
        ),
        conservation_labels=(
            tuple(table.labels) if table.labels is not None else None
        ),
        conservation_colors=(
            tuple(table.colors) if table.colors is not None else None
        ),
    )


def _resolve_linear_comparison_selector(
    records: Sequence[SeqRecord],
    selector: str,
    *,
    table_path: str,
    row_number: int,
    column: str,
) -> tuple[int, ...]:
    """The drawn units a comparisons-table selector names: ``#n`` counts records and a record ID names the record
    whose units carry it (``drawn_record_units``); a record without windows is one unit. A selector that names no
    record and reads as a region (``RECORD:START-END[:origin]``, as the records table writes a window) names the
    one window of that record with that span."""

    where = f"{table_path}: row {row_number}, column {column!r}"
    diagnostic = {"code": "COMPARISON_INPUT", "row": row_number}
    groups = drawn_record_units(records)
    number = record_index_number(selector) if selector.startswith("#") else None
    if number is not None:
        if 0 <= number - 1 < len(groups):
            return groups[number - 1]
    elif not selector.startswith("#"):
        matches = [group for group in groups if any(str(records[index].id) == selector for index in group)]
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ValidationError(
                f"{where}: selector {selector!r} matched multiple record IDs; use #index.", diagnostic=diagnostic
            )
    window = _comparison_window_selector(selector)
    if window is not None:
        record_selector, start, end = window
        try:
            units = select_record_units(records, record_selector)
        except ValueError as exc:
            raise ValidationError(
                f"{where}: unresolved record in window selector {selector!r}: {exc}", diagnostic=diagnostic
            ) from exc
        named = [unit for unit in units if record_window_span(records[unit]) == (start, end)]
        if len(named) == 1:
            return (named[0],)
        spans = ", ".join(
            region_text(*span, False) for unit in units if (span := record_window_span(records[unit])) is not None
        )
        raise ValidationError(
            f"{where}: record {record_selector.label()} has no window {region_text(start, end, False)}"
            + (f"; its windows are {spans}." if spans else "; it is drawn whole, so name it without a span."),
            diagnostic=diagnostic,
        )
    if selector.startswith("#") and number is None:
        raise ValidationError(f"{where}: invalid record selector {selector!r}.", diagnostic=diagnostic)
    raise ValidationError(f"{where}: unresolved record selector {selector!r}.", diagnostic=diagnostic)


class DrawnUnitName(NamedTuple):
    """What a comparisons-table selector reads of a drawn unit: its record ID, the key of the record it draws
    (shared by the windows of one record), and its window span (``None`` for a record drawn whole)."""

    record_id: str
    record_key: str
    window_span: tuple[int, int] | None


def comparison_selector_units(units: Sequence[DrawnUnitName], selector: str) -> tuple[int, ...] | None:
    """The units ``selector`` names among ``units`` as the comparisons table reader resolves it
    (:func:`_resolve_linear_comparison_selector`), or ``None`` when it names none: a writer checks a selector
    with it before writing one."""

    records = []
    for unit in units:
        annotations: dict[str, object] = {WINDOW_OF_ANNOTATION: unit.record_key}
        if unit.window_span is not None:
            annotations[WINDOW_SPAN_ANNOTATION] = unit.window_span
        records.append(SeqRecord(Seq(""), id=unit.record_id, annotations=annotations))
    try:
        return _resolve_linear_comparison_selector(records, selector, table_path="", row_number=0, column="")
    except ValidationError:
        return None


def _comparison_window_selector(selector: str) -> tuple[RecordSelector, int, int] | None:
    """The record and span of a window selector ``RECORD:START-END[:origin]``, or ``None`` when ``selector`` is
    not one. Orientation is not part of a window's name, so ``:rc`` and a reversed span are not window selectors."""

    try:
        spec = parse_region_spec(selector)
    except ValueError:
        return None
    record_selector = spec.record_selector()
    if record_selector is None or spec.file_selector is not None or spec.reverse_complement:
        return None
    return record_selector, spec.start, spec.end


def _requested_wrap_length(options: LinearDiagramOptions) -> LinearWrapLength:
    """A request's wrap length, read without building its config: a figure that does not wrap pays nothing."""

    overrides = options.config_overrides or {}
    if LINEAR_WRAP_LENGTH_PATH in overrides:
        return linear_wrap_length_from_config(overrides[LINEAR_WRAP_LENGTH_PATH])
    if isinstance(options.config, GbdrawConfig):
        return options.config.canvas.linear.wrap_length
    if options.config is None:
        return "none"
    return linear_wrap_length_from_config(raw_config_value(options.config, LINEAR_WRAP_LENGTH_PATH))


def _with_wrap_defaults(options: LinearDiagramOptions, defaults: LinearWrapDefaults | None) -> LinearDiagramOptions:
    """``options`` with the leaves and slots a request left to the wrap defaults.

    The defaults apply over what the request carries, a config or overrides alike: they hold only the leaves the
    request left unset (``linear_wrap_default_overrides``), so they replace only values its mode profile filled in,
    such as the plain axis colour a cut record's ruler takes the place of (OV-2495).
    """

    if defaults is None:
        return replace(options, wrap_defaults=None)
    config = options.config
    tracks = options.tracks if defaults.track_order is None else defaults.track_order.tracks(cut=True)
    if config is None:
        return replace(
            options,
            config_overrides={**(options.config_overrides or {}), **defaults.config_overrides},
            tracks=tracks,
            wrap_defaults=None,
        )
    return replace(
        options,
        config=(
            with_config_overrides(config, defaults.config_overrides)
            if isinstance(config, GbdrawConfig)
            else modify_config_dict(config, defaults.config_overrides)
        ),
        tracks=tracks,
        wrap_defaults=None,
    )


def decide_linear_wrap(
    options: LinearDiagramOptions,
    *,
    records: Sequence[SeqRecord],
    layout: LinearMultiRecordOptions | None,
    similarity_alignment: object = None,
) -> LinearDiagramOptions:
    """The one place a wrap length meets the record lengths: what the figure draws with, or its refusal.

    With nothing to cut the request is the same request without the option. When a record is cut, the wrap defaults
    apply and what rows cannot draw is refused, before any comparison search.
    """

    defaults = options.wrap_defaults
    row_length = resolve_wrap_length(_requested_wrap_length(options), [len(record.seq) for record in records])
    if row_length is None:
        if defaults is None:
            return options
        tracks = options.tracks if defaults.track_order is None else defaults.track_order.tracks(cut=False)
        return replace(options, tracks=tracks, wrap_defaults=None)
    decided = _with_wrap_defaults(options, defaults)
    refuse_undrawable_wrap(
        apply_config_overrides(decided.config, decided.config_overrides),
        resolve_record_row_positions(records, layout.multi_record_positions if layout is not None else None)[1],
        wrap_blocking_inputs(
            blast_files=decided.blast_files,
            comparison_table_file=decided.comparison_table_file,
            linear_comparisons=decided.linear_comparisons,
            protein_comparisons=decided.protein_comparisons,
            orthogroups=decided.orthogroups,
            collinearity_blocks=decided.collinearity_blocks,
            losat_search=decided.losat_search,
            similarity_alignment=similarity_alignment,
        ),
        row_length=row_length,
        records=records,
        # The x of the records an anchor places is decided with the anchor, after the records are known: the
        # build refuses what it cannot wrap then (``refuse_wrapped_anchor_moves``).
        translations=layout.record_translations if layout is not None and decided.anchor_feature is None else (),
    )
    return decided


def resolve_linear_options(
    options: LinearDiagramOptions,
    *,
    records: Sequence[SeqRecord],
    layout: LinearMultiRecordOptions | None,
) -> LinearDiagramOptions:
    """Resolve an optional comparison manifest against displayed records."""

    positions = layout.multi_record_positions if layout is not None else None
    _ordered, rows_by_record = resolve_record_row_positions(records, positions)
    multi_record_rows = len(set(rows_by_record)) < len(records)
    if multi_record_rows and options.blast_files:
        raise ValidationError(
            "-b/--blast is ambiguous when a Linear row contains multiple records; "
            "use a comparison table with explicit query and subject selectors."
        )
    adjacent_pairs = max(len(records) - 1, 0)
    if options.blast_files and len(options.blast_files) > adjacent_pairs:
        # A file past the last pair would name a record that is not loaded.
        raise ValidationError(
            f"Too many -b/--blast files (expected at most {adjacent_pairs}): each file "
            f"belongs to one adjacent pair of the {len(records)} loaded record(s).",
            diagnostic={"code": "COMPARISON_INPUT"},
        )
    protein_tables = options.protein_comparisons or ()
    if not multi_record_rows and len(protein_tables) > adjacent_pairs:
        # Ordered table n joins records n and n + 1 in every protein mode
        # (OV-976, OV-977); a multi-record row refuses ordered tables later.
        raise ValidationError(
            f"Protein comparison table {adjacent_pairs + 1} has no adjacent record pair: "
            f"the request has {len(records)} records.",
            diagnostic={"code": "COMPARISON_INPUT"},
        )
    if options.comparison_table_file is None:
        return options
    table = read_comparisons_table(options.comparison_table_file)
    record_units = {unit: group for group in drawn_record_units(records) for unit in group}
    comparisons: list[LinearComparison] = []
    losat_pairs: list[tuple[int, int]] = []
    for row in table.rows:
        query_units = _resolve_linear_comparison_selector(
            records,
            row.query,
            table_path=table.table_path,
            row_number=row.row_number,
            column="query",
        )
        subject_units = _resolve_linear_comparison_selector(
            records,
            row.subject,
            table_path=table.table_path,
            row_number=row.row_number,
            column="subject",
        )
        query_index, subject_index = query_units[0], subject_units[0]
        query_row = int(rows_by_record[query_index])
        subject_row = int(rows_by_record[subject_index])
        if query_units == subject_units:
            raise ValidationError(
                f"{table.table_path}: row {row.row_number}, column 'subject': "
                "query and subject resolved to the same record."
            )
        if query_row == subject_row or abs(query_row - subject_row) != 1:
            topology = (
                "different rows"
                if query_row == subject_row
                else "adjacent rows"
            )
            raise ValidationError(
                f"{table.table_path}: row {row.row_number}, column 'subject': "
                f"query row {query_row + 1} and subject row "
                f"{subject_row + 1} must be in {topology}."
            )
        if row.source == "losat":
            # A record drawn in windows is searched window by window against each window of the other (D-01).
            losat_pairs.extend(comparison_unit_pairs(query_units, subject_units))
            continue
        # A file compares whole drawn records, also when a selector names one window of a record.
        refuse_file_comparison_of_windows(
            record_units[query_index],
            record_units[subject_index],
            records=(row.query, row.subject),
            where=f"{table.table_path}: row {row.row_number}",
            row=row.row_number,
        )
        try:
            # The file's staged identity keys the stages that take its rows.
            matches = read_comparison_resource(row.blast)
        except ValidationError as exc:
            raise ValidationError(
                f"{table.table_path}: row {row.row_number}, column 'blast': {exc}"
            ) from exc
        comparisons.append(
            LinearComparison(query_index, subject_index, matches, source_name=Path(row.blast).name)
        )
    search = options.losat_search
    if losat_pairs and search is None:
        row_number = next(row.row_number for row in table.rows if row.source == "losat")
        raise ValidationError(
            f"{table.table_path}: row {row_number} has source=losat; pass --losat "
            "(Python: losat_search) to run the search.",
            diagnostic={"code": "COMPARISON_INPUT", "reason": "LOSAT_PLAN", "row": row_number},
        )
    if search is not None and not losat_pairs:
        raise ValidationError(
            f"{table.table_path}: --losat is set but no comparisons table row has "
            "source=losat.",
            diagnostic={"code": "COMPARISON_INPUT", "reason": "LOSAT_PLAN"},
        )
    if search is not None and search.pairs is not None:
        raise ValidationError(
            "Pass LOSAT record pairs either in the comparisons table or as "
            "losat_search pairs, not both.",
            diagnostic={"code": "COMPARISON_INPUT", "reason": "LOSAT_PLAN"},
        )
    if search is not None and search.program == "losatp" and search.losatp_mode != "pairwise":
        raise ValidationError(
            f"{table.table_path}: source=losat rows list the pairs to compare, and --losatp_mode "
            f"{search.losatp_mode} does not compare listed pairs; use --losatp_mode pairwise "
            "(Python: losatp_mode='pairwise').",
            diagnostic={"code": "COMPARISON_INPUT", "reason": "LOSAT_PLAN", "field": "losatp_mode"},
        )
    return replace(
        options,
        comparison_table_file=None,
        linear_comparisons=tuple(comparisons),
        losat_search=(
            replace(search, pairs=tuple(losat_pairs)) if search is not None else None
        ),
    )


def resolve_circular_batch_outputs(
    policy: CircularBatchOutputPolicy,
    records: Sequence[SeqRecord],
) -> tuple[RenderOutputRequest, ...]:
    """Resolve one output request per displayed record after expansion."""

    prefix = str(policy.output_prefix) if policy.output_prefix is not None else None
    if prefix is None:
        raw_prefixes: list[str] = []
        used: set[str] = set()
        for record in records:
            base = resolve_implicit_record_output_prefix(record.id)
            candidate = base
            suffix = 2
            while candidate in used:
                candidate = f"{base}_{suffix}"
                suffix += 1
            used.add(candidate)
            raw_prefixes.append(candidate)
    elif len(records) == 1:
        raw_prefixes = [prefix]
    else:
        raw_prefixes = [
            f"{prefix}_{index}" for index in range(1, len(records) + 1)
        ]
    outputs: list[RenderOutputRequest] = []
    for raw_prefix in raw_prefixes:
        path = Path(raw_prefix)
        outputs.append(
            RenderOutputRequest(
                output_prefix=path.name,
                output_directory=(
                    path.parent if path.parent != Path(".") else None
                ),
                formats=policy.formats,
                overwrite=policy.overwrite,
                interactive_metadata_policy=policy.interactive_metadata_policy,
                embed_fonts=policy.embed_fonts,
                dpi=policy.dpi,
            )
        )
    return tuple(outputs)


def resolve_implicit_record_output_prefix(record_id: object) -> str:
    """Require a record-derived output name to be one filename component."""

    raw_record_id = str(record_id)
    try:
        output = RenderOutputRequest(output_prefix=raw_record_id)
    except ValidationError as exc:
        raise ValidationError(
            f"Record ID {raw_record_id!r} cannot be used as an implicit output "
            "filename prefix. Specify an explicit output prefix."
        ) from exc
    return output.output_prefix


__all__ = [
    "DrawnUnitName",
    "RecordInputManifest",
    "ResolvedRecordCollection",
    "ResolvedRecordProvenance",
    "depth_track_inputs_from_cli",
    "parse_cli_bool",
    "prepared_parsed_source",
    "feature_drawn_start",
    "resolve_anchor_feature",
    "comparison_selector_units",
    "project_origin_feature_starts",
    "resolve_origin_feature_start",
    "with_origin_feature_starts",
    "record_input_manifest_from_paths",
    "record_input_manifest_from_table",
    "resolve_circular_batch_outputs",
    "resolve_circular_options",
    "resolve_implicit_record_output_prefix",
    "resolve_linear_options",
    "resolve_record_inputs",
    "records_holding_their_sequences",
    "reverse_resolved_records",
]
