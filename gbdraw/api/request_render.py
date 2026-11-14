"""Normalize and render typed diagram requests without CLI/session orchestration."""

from __future__ import annotations

import copy
import json
import logging
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field, replace
from pathlib import Path
from time import perf_counter
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Hashable,
    Iterator,
    Literal,
    Mapping,
    MutableMapping,
    NamedTuple,
    Sequence,
    TypeAlias,
    TypedDict,
    TypeVar,
    cast,
    overload,
)

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from svgwrite import Drawing

from gbdraw.analysis.collinearity import (
    CollinearityResult,
    LosslessCollinearityParameters,
)
from gbdraw.analysis.depth_tracks import (
    DepthTrackSpec,
    depth_track_count,
    normalize_depth_tracks,
)
from gbdraw.analysis.protein_artifacts import (
    CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA,
    DerivedLosatpPair,
    DerivedLosatpRecord,
    DerivedLosatpSettings,
    derived_losatp_cache_key,
    derived_losatp_identity,
    derived_losatp_provenance,
    derived_pair_payload,
    effective_collinearity_unit_kinds,
    is_current_derived_protein_artifact,
    normalize_derived_losatp_settings,
    validate_current_derived_protein_artifacts,
)
from gbdraw.analysis.protein_colinearity import (
    LosatpCacheManager,
    OrthogroupGraphResult,
    OrthogroupMember,
    OrthogroupResult,
    ProteinExtractionResult,
    ProteinIdentityManifest,
    extract_web_stable_cds_proteins,
    is_json_exact_protein_identity_manifest,
    is_protein_losat_cache_entry,
    source_protein_letters,
    validate_protein_identity_manifest,
    validate_protein_raw_entry_references,
)
from gbdraw.analysis.windows import window_source
from gbdraw.annotations import (
    FeatureIdentitySpan,
    FeatureSelector,
    ResolutionWarning,
    ResolvedAnnotationBundle,
    read_annotation_options,
    resolve_annotations,
)
from gbdraw.annotations.planning import annotation_sets_left_out
from gbdraw.comparisons.circular_losat import (
    resolve_circular_conservation_losat,
    ring_subject_sources,
)
from gbdraw.comparisons.linear_losat import (
    NucleotideLosatRecord,
    linear_search_spaces,
    resolve_linear_nucleotide_losat,
)
from gbdraw.comparisons.losat_jobs import (
    GenomeDatabaseSizes,
    PerFileWarning,
    genome_database_sizes,
    genome_file_names,
    losat_record_uids,
    losat_source_ids,
    losatp_job_specs,
    source_file_name,
)
from gbdraw.exceptions import ValidationError
from gbdraw.features.overrides import (
    FeatureIdentityNotice,
    FeatureOverride,
    ResolvedFeatureOverride,
    bind_feature_overrides,
    normalize_feature_overrides,
    read_feature_override_table,
)
from gbdraw.features.placement import (
    ResolvedRecordFeatureInputs,
    read_feature_placement_table,
    resolve_record_feature_inputs,
)
from gbdraw.features.source import (
    FeatureIdentity,
    IdentityBindings,
    SourceFeatureIdentity,
    capture_source_features,
    resolve_feature_identities,
)
from gbdraw.features.visibility import (
    read_feature_visibility_file,
    resolve_candidate_feature_types,
)
from gbdraw.io.colors import read_color_table, read_default_color_overrides
from gbdraw.io.comparison_sequences import (
    ComparisonGenome,
    ComparisonSequenceFile,
    RepeatedRecordIds,
    comparison_genome,
    read_comparison_sequence_file,
)
from gbdraw.io.comparisons import COMPARISON_COLUMNS, read_comparison_resource
from gbdraw.io.genome import load_fasta, parse_gff_fasta, source_genome_length
from gbdraw.io.record_select import select_record_units
from gbdraw.labels.filtering import (
    read_filter_list_file,
    read_label_override_file,
    read_qualifier_priority_file,
)
from gbdraw.labels.fit import WARN_SEARCH_SECONDS, LabelFitReport
from gbdraw.layout.linear_multi_record import drawn_unit_rows, grid_display_pairs
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.record_transform import detected_topology
from gbdraw.layout.record_units import left_off_circle, record_units
from gbdraw.linear_comparison import auto_orient_reversed_records, search_frame_linear_comparisons
from gbdraw.mode_profiles import ComparisonThresholds
from gbdraw.io.record_sequences import detached_records
from gbdraw.layout.record_placement import (
    FeatureAnchorPlacement,
    parse_record_row_position,
    record_order_selector,
    resolve_record_row_positions,
)
from gbdraw.layout.similarity_alignment import (
    SimilarityAlignmentCandidate,
    SimilarityAlignmentPlan,
    SimilarityAlignmentReference,
    SimilarityAlignmentReferenceError,
    resolve_similarity_alignment,
)
from gbdraw.render.formats import resolve_output_paths
from gbdraw.render.interactive_context import (
    interactive_annotation_rows,
    interactive_svg_context,
    require_interactive_svg_metadata,
)
from gbdraw.render.interactive_svg import InteractiveSvgContext
from gbdraw.render_stages import render_stage
from gbdraw.tracks import (
    CircularTrackSlot,
    LinearTrackSlot,
    parse_circular_track_slots,
    parse_linear_track_slots,
)
from gbdraw.tracks.parsing import slots_draw_features

from .diagram import (
    DEFAULT_SELECTED_FEATURES,
    LinearDiagramBuildResult,
    LinearDiagramMetadata,
    build_circular_diagram,
    build_circular_multi_diagram,
    build_linear_diagram_result,
    label_fit_report,
    resolve_feature_inputs,
)
from .input_summary import InputSummary, summarize_inputs
from .io import load_render_gbks as load_gbks  # the render's reader; tests replace it under this name
from .migration_report import SessionMigrationReport
from .options import (
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    ColorOptions,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    LinearRecordTranslation,
    LosatSearchOptions,
    losatp_analysis_mode,
)
from .prepared import (
    RECORD_IDENTITY_TAG,
    ParsedRecordInputs,
    ResolvedFeatureInputs,
    get_or_build_parsed_source,
    get_or_build_resolved_records,
    prepared_input_cache_active,
    prepared_resource_identity,
    register_prepared_resource_value,
    register_prepared_value_identity,
)
from .record_planning import (
    FastaLoader,
    GenBankLoader,
    GffFastaLoader,
    ResolvedRecordCollection,
    ResolvedRecordDisplay,
    ResolvedRecordProvenance,
    _file_source,
    _gff_feature_types,
    _load_source_records,
    _prepared_source_cache_spec,
    _similarity_candidate_from_orthogroup_member,
    _source_details,
    _source_key,
    alignment_record_keys,
    comparison_record_bindings,
    decide_linear_wrap,
    effective_record_presentation,
    has_source_bound_results,
    project_similarity_alignment_centers,
    project_source_bound_results,
    record_is_reverse_complemented,
    resolve_anchor_feature,
    resolve_circular_batch_outputs,
    resolve_circular_options,
    resolve_implicit_record_output_prefix,
    resolve_linear_options,
    resolve_record_inputs,
    reverse_resolved_records,
    similarity_alignment_evidence_edges,
    with_origin_feature_starts,
)
from .render import preflight_output_paths, save_figure_to
from .requests import (
    CircularBatchRequest,
    CircularDiagramRequest,
    DiagramRequest,
    FastaInputSource,
    GenBankInputSource,
    GffFastaInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordCardinality,
    RecordCollectionOptions,
    RecordInput,
    RecordInputSource,
    RecordPresentation,
    RenderOutputRequest,
    drawn_record_identity,
)

if TYPE_CHECKING:
    from gbdraw.linear_comparison import ComparisonTableTally, LinearComparison

    from .requests import CircularBatchOutputPolicy

logger = logging.getLogger(__name__)

# A helper that takes one of these returns the same type: mode never changes.
_RequestT = TypeVar("_RequestT", CircularDiagramRequest, CircularBatchRequest, LinearDiagramRequest)
_OptionsT = TypeVar("_OptionsT", CircularDiagramOptions, LinearDiagramOptions)


_REQUEST_RENDER_DIAGNOSTICS: ContextVar[MutableMapping[str, Any] | None] = (
    ContextVar("gbdraw_request_render_diagnostics", default=None)
)


@contextmanager
def _capture_request_render_diagnostics(
    diagnostics: MutableMapping[str, Any] | None,
) -> Iterator[None]:
    """Collect coarse Web-only render timings without changing public results."""

    token = _REQUEST_RENDER_DIAGNOSTICS.set(diagnostics)
    try:
        yield
    finally:
        _REQUEST_RENDER_DIAGNOSTICS.reset(token)


_DERIVED_LOSAT_ENTRIES_WANTED: ContextVar[bool] = ContextVar(
    "gbdraw_derived_losat_entries_wanted", default=True
)


@contextmanager
def _without_derived_losat_entries() -> Iterator[None]:
    """Leave the derived LOSAT cache entries of a linear build unbuilt.

    The Web wrapper and the CLI enter this: no Generate reply, Result or
    catalog reads the entries, and the Web Session builds its own. It changes
    one field of the result: ``losat_derived_cache_entries`` is ``()``. The
    checks on their inputs (identity manifest, canonical comparison columns,
    the canonical encoding of the computed result) still run, so a request
    fails where it failed before. The Python API never enters it.
    """

    token = _DERIVED_LOSAT_ENTRIES_WANTED.set(False)
    try:
        yield
    finally:
        _DERIVED_LOSAT_ENTRIES_WANTED.reset(token)


@contextmanager
def _request_render_diagnostic_phase(name: str) -> Iterator[None]:
    diagnostics = _REQUEST_RENDER_DIAGNOSTICS.get()
    if diagnostics is None:
        yield
        return
    started_at = perf_counter()
    try:
        yield
    finally:
        timings = diagnostics.setdefault("timingsMs", {})
        timings[name] = float(timings.get(name, 0.0)) + (
            perf_counter() - started_at
        ) * 1000.0


def _record_request_render_diagnostic_metric(name: str, value: int = 1) -> None:
    diagnostics = _REQUEST_RENDER_DIAGNOSTICS.get()
    if diagnostics is None:
        return
    metrics = diagnostics.setdefault("metrics", {})
    metrics[name] = int(metrics.get(name, 0)) + int(value)


def _artifact_text_characters(value: object) -> int:
    if isinstance(value, str):
        return len(value)
    if isinstance(value, Mapping):
        return sum(
            len(str(key)) + _artifact_text_characters(entry)
            for key, entry in value.items()
        )
    if isinstance(value, Sequence) and not isinstance(value, (bytes, bytearray)):
        return sum(_artifact_text_characters(entry) for entry in value)
    return 0


def _default_color_overrides(colors: ColorOptions) -> DataFrame | None:
    """The decoded override rows of the default-colors file, keyed by the file."""

    if not colors.default_colors_file:
        return None
    return register_prepared_resource_value(
        read_default_color_overrides(colors.default_colors_file),
        colors.default_colors_file,
    )


def _resolve_request_option_tables(
    options: _OptionsT,
    *,
    mode: Literal["circular", "linear"],
) -> _OptionsT:
    """Resolve file-backed option tables once before record and drawing plans."""

    colors = options.colors
    changed = False
    if colors is not None:
        color_table = colors.color_table
        if color_table is None and colors.color_table_file is not None:
            color_table = register_prepared_resource_value(
                read_color_table(colors.color_table_file),
                colors.color_table_file,
            )
        # A Default colors file stays a file unless a table shadows it: a
        # Session saves the rows the user gave, which apply over the palette
        # when the request is drawn (`_prepare_diagram_inputs`), so a later
        # palette change reaches the Session (OV-3302).
        default_colors_file = None if colors.default_colors is not None else colors.default_colors_file
        if colors.color_table_file is not None or default_colors_file != colors.default_colors_file:
            colors = replace(
                colors, color_table=color_table, color_table_file=None, default_colors_file=default_colors_file
            )
            changed = True

    feature_visibility_table = options.feature_visibility_table
    if options.feature_visibility_table_file is not None:
        if feature_visibility_table is None:
            feature_visibility_table = register_prepared_resource_value(
                read_feature_visibility_file(
                    options.feature_visibility_table_file
                ),
                options.feature_visibility_table_file,
            )
        changed = True
    label_whitelist_table = options.label_whitelist_table
    if options.label_whitelist_file is not None:
        if label_whitelist_table is not None:
            raise ValidationError(
                "Pass either label_whitelist_table or label_whitelist_file, not both."
            )
        if label_whitelist_table is None:
            label_whitelist_table = register_prepared_resource_value(
                read_filter_list_file(options.label_whitelist_file),
                options.label_whitelist_file,
            )
        changed = True
    qualifier_priority_table = options.qualifier_priority_table
    if options.qualifier_priority_file is not None:
        if qualifier_priority_table is not None:
            raise ValidationError(
                "Pass either qualifier_priority_table or "
                "qualifier_priority_file, not both."
            )
        if qualifier_priority_table is None:
            qualifier_priority_table = register_prepared_resource_value(
                read_qualifier_priority_file(options.qualifier_priority_file),
                options.qualifier_priority_file,
            )
        changed = True
    label_override_table = options.label_override_table
    if options.label_override_file is not None:
        if label_override_table is not None:
            raise ValidationError(
                "Pass either label_override_table or label_override_file, not both."
            )
        if label_override_table is None:
            label_override_table = register_prepared_resource_value(
                read_label_override_file(options.label_override_file),
                options.label_override_file,
            )
        changed = True
    annotations = options.annotations
    if annotations is not None and annotations.table_file is not None:
        annotations = register_prepared_resource_value(
            read_annotation_options(annotations, mode=mode),
            annotations.table_file,
        )
        changed = True
    if not changed:
        return options
    return replace(
        options,
        colors=colors,
        annotations=annotations,
        feature_visibility_table=feature_visibility_table,
        feature_visibility_table_file=None,
        label_whitelist_table=label_whitelist_table,
        label_whitelist_file=None,
        qualifier_priority_table=qualifier_priority_table,
        qualifier_priority_file=None,
        label_override_table=label_override_table,
        label_override_file=None,
    )


@dataclass
class _ComparisonSequenceSources:
    """Memoized Circular comparison-genome records shared by a request batch.

    ``genomes`` are the files of each ring (``conservation_genome_files``);
    each file is parsed once, through the prepared parsed-source cache.
    """

    genomes: tuple[tuple[str, ...], ...]
    _genomes: tuple[ComparisonGenome | None, ...] | None = None
    _records: tuple[tuple[SeqRecord, ...], ...] | None = None

    def load(self) -> tuple[tuple[SeqRecord, ...], ...]:
        if self._records is None:
            self._records = tuple(
                genome.records if genome is not None else () for genome in self.files()
            )
        return self._records

    @staticmethod
    def _read(path: str) -> ComparisonSequenceFile:
        identity = prepared_resource_identity(path)
        if identity is None:
            return read_comparison_sequence_file(path)
        return get_or_build_parsed_source(
            ("parsed-source-v1", "comparison-sequence", identity),
            frozenset({identity}),
            lambda: read_comparison_sequence_file(path),
            publish=lambda value: bool(value.records),
        )

    def files(self) -> tuple[ComparisonGenome | None, ...]:
        """Each comparison genome read once by the one reader (design D12); ``None`` for a ring without one."""

        if self._genomes is None:
            self._genomes = tuple(
                comparison_genome(tuple(self._read(path) for path in paths)) if paths else None
                for paths in self.genomes
            )
        return self._genomes

    def repeated_record_ids(self) -> tuple[RepeatedRecordIds, ...]:
        """The repeated record IDs of the genomes this render read (OV-2954), each file once.

        Only a genome the render read (for the LOSAT query or the interactive
        span export) is reported: a ring that reads no file names no record.
        """

        found: dict[RepeatedRecordIds, None] = {}
        for genome in self._genomes or ():
            if genome is not None:
                found.update(dict.fromkeys(genome.repeated_record_ids))
        return tuple(found)


@dataclass(frozen=True)
class PreparedDiagramInputs:
    """Resolved values reused by record loading, drawing, and metadata."""

    features: ResolvedFeatureInputs
    gff_candidate_features: tuple[str, ...]
    gff_keep_all_features: bool
    comparison_sequences: _ComparisonSequenceSources | None = None
    # Aligned with the planned records; empty without identity-addressed inputs.
    record_features: tuple[ResolvedRecordFeatureInputs, ...] = ()
    feature_identity_notices: tuple[FeatureIdentityNotice, ...] = ()
    # Annotation featureIdentity targets, keyed by (record key, feature ID).
    feature_bindings: Mapping[tuple[str, str], IdentityBindings] | None = None
    # Source parses shared by the loads of one ``_load_request_records`` call.
    parsed_sources: dict[tuple[object, ...], ParsedRecordInputs] | None = field(default=None, repr=False, compare=False)


def _is_current_nucleotide_losat_entry(entry: Mapping[str, Any]) -> bool:
    return (
        entry.get("schema") == 2
        and entry.get("kind") == "raw-losat"
        and entry.get("identityKind") in {None, "nucleotide"}
        and str(entry.get("program") or "").lower() != "blastp"
        and isinstance(entry.get("key"), str)
        and bool(entry.get("key"))
        and isinstance(entry.get("text"), str)
    )


def _detached_artifact_entries(
    entries: Sequence[Mapping[str, Any]],
    *,
    label: str,
) -> tuple[Mapping[str, Any], ...]:
    if isinstance(entries, (str, bytes, bytearray)):
        raise ValidationError(f"{label} must be a sequence of objects.")
    if not all(isinstance(entry, Mapping) for entry in entries):
        raise ValidationError(f"{label} must be a sequence of objects.")
    try:
        detached = tuple(copy.deepcopy(dict(entry)) for entry in entries)
    except (TypeError, ValueError) as exc:
        raise ValidationError(f"{label} must be a sequence of objects.") from exc
    if len(detached) != len(entries):
        raise ValidationError(f"{label} must be a sequence of objects.")
    return detached


@dataclass(frozen=True)
class CurrentRequestArtifacts:
    """Validated current artifacts accepted by the typed render boundary."""

    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()
    losat_derived_cache_entries: tuple[Mapping[str, Any], ...] = ()
    protein_identity_manifest: Mapping[str, Any] | None = None
    protein_source_mode: Literal[
        "none", "pairwise", "orthogroup", "collinear"
    ] | None = None

    def __post_init__(self) -> None:
        with _request_render_diagnostic_phase("artifactCopy"):
            raw_entries = _detached_artifact_entries(
                self.losat_cache_entries,
                label="losat_cache_entries",
            )
            derived_entries = _detached_artifact_entries(
                self.losat_derived_cache_entries,
                label="losat_derived_cache_entries",
            )
            manifest = self.protein_identity_manifest
            # A proof was checked when it was made and is read-only: it is
            # kept as it is. Raw data is detached here and checked below.
            authority = manifest if isinstance(manifest, ProteinIdentityManifest) else None
            stored_manifest: Mapping[str, object] | None = authority
            detached_manifest: Mapping[str, Any] | None = None
            if manifest is not None and authority is None:
                if not isinstance(manifest, Mapping):
                    raise ValidationError("protein_identity_manifest must be an object.")
                detached_manifest = copy.deepcopy(dict(manifest))
            if (
                _REQUEST_RENDER_DIAGNOSTICS.get() is not None
                and (raw_entries or derived_entries or detached_manifest is not None)
            ):
                _record_request_render_diagnostic_metric("artifactBundleCopyCount")
                _record_request_render_diagnostic_metric(
                    "artifactEntryCount",
                    len(raw_entries) + len(derived_entries),
                )
                _record_request_render_diagnostic_metric(
                    "artifactTextCharacters",
                    _artifact_text_characters(raw_entries)
                    + _artifact_text_characters(derived_entries)
                    + _artifact_text_characters(detached_manifest),
                )

        with _request_render_diagnostic_phase("artifactValidation"):
            protein_entries: list[Mapping[str, Any]] = []
            raw_keys: set[str] = set()
            for index, entry in enumerate(raw_entries):
                if is_protein_losat_cache_entry(entry):
                    protein_entries.append(entry)
                elif not _is_current_nucleotide_losat_entry(entry):
                    raise ValidationError(
                        "Unsupported current LOSAT artifact at "
                        f"losat_cache_entries[{index}]."
                    )
                key = str(entry["key"])
                if key in raw_keys:
                    raise ValidationError(
                        f"Duplicate current LOSAT artifact key: {key!r}."
                    )
                raw_keys.add(key)

            derived_keys: set[str] = set()
            for index, entry in enumerate(derived_entries):
                if not is_current_derived_protein_artifact(entry):
                    raise ValidationError(
                        "Unsupported current derived LOSATP artifact at "
                        f"losat_derived_cache_entries[{index}]."
                    )
                key = str(entry["key"])
                if key in derived_keys:
                    raise ValidationError(
                        f"Duplicate current derived LOSATP artifact key: {key!r}."
                    )
                derived_keys.add(key)

            if detached_manifest is not None:
                _record_request_render_diagnostic_metric(
                    "proteinManifestFullValidationCount"
                )
                authority = validate_protein_identity_manifest(detached_manifest)
                # The proof checks the entries either way; it is kept in place
                # of the detached copy only when it holds plain JSON data.
                stored_manifest = (
                    authority
                    if is_json_exact_protein_identity_manifest(authority)
                    else detached_manifest
                )
            if protein_entries and authority is None:
                raise ValidationError(
                    "Current protein LOSATP artifacts require protein_identity_manifest."
                )
            if authority is not None:
                for index, entry in enumerate(protein_entries):
                    _record_request_render_diagnostic_metric(
                        "proteinRawReferenceValidationCount"
                    )
                    if not validate_protein_raw_entry_references(
                        entry,
                        authority,
                    ):
                        raise ValidationError(
                            "Current protein LOSATP artifact does not resolve through "
                            f"protein_identity_manifest: losat_cache_entries[{index}]."
                        )
            if derived_entries:
                _record_request_render_diagnostic_metric(
                    "derivedArtifactFullValidationCount"
                )
                validate_current_derived_protein_artifacts(
                    derived_entries,
                    authority,
                )
            if self.protein_source_mode not in {
                None,
                "none",
                "pairwise",
                "orthogroup",
                "collinear",
            }:
                raise ValidationError(
                    f"Unsupported protein_source_mode: {self.protein_source_mode!r}."
                )

        object.__setattr__(self, "losat_cache_entries", raw_entries)
        object.__setattr__(self, "losat_derived_cache_entries", derived_entries)
        object.__setattr__(self, "protein_identity_manifest", stored_manifest)


@dataclass(frozen=True)
class PreparedDiagramRequest:
    """A validated request with normalized records and its SVG drawing."""

    mode: Literal["circular", "linear"]
    request: CircularDiagramRequest | LinearDiagramRequest
    records: tuple[SeqRecord, ...]
    drawing: Drawing
    inputs: PreparedDiagramInputs | None = None
    linear_metadata: LinearDiagramMetadata | None = None
    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()
    losat_derived_cache_entries: tuple[Mapping[str, Any], ...] = ()
    protein_identity_manifest: Mapping[str, Any] | None = None
    transforms: tuple[RecordDisplayTransform, ...] = ()
    resolved_annotations: ResolvedAnnotationBundle = ResolvedAnnotationBundle(())
    input_summary: InputSummary = InputSummary(())

    @property
    def annotation_warnings(self) -> tuple[ResolutionWarning, ...]:
        return self.resolved_annotations.warnings

    @property
    def feature_identity_notices(self) -> tuple[FeatureIdentityNotice, ...]:
        return self.inputs.feature_identity_notices if self.inputs is not None else ()

    @property
    def label_fit_report(self) -> LabelFitReport | None:
        """Crossing label leaders drawn anyway: the fit report."""
        return label_fit_report(self.drawing)


@dataclass(frozen=True)
class RequestRenderResult:
    """Files and normalized inputs produced by one request render."""

    mode: Literal["circular", "linear"]
    request: CircularDiagramRequest | LinearDiagramRequest
    records: tuple[SeqRecord, ...]
    drawing: Drawing
    output_paths: tuple[Path, ...]
    interactive_context: InteractiveSvgContext | None = None
    linear_metadata: LinearDiagramMetadata | None = None
    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()
    losat_derived_cache_entries: tuple[Mapping[str, Any], ...] = ()
    protein_identity_manifest: Mapping[str, Any] | None = None
    annotation_warnings: tuple[ResolutionWarning, ...] = ()
    feature_identity_notices: tuple[FeatureIdentityNotice, ...] = ()
    # What the diagram read from its inputs and drew of them (CLI log, Web status).
    input_summary: InputSummary = InputSummary(())
    migration_report: SessionMigrationReport = field(
        default_factory=SessionMigrationReport
    )

    @property
    def label_fit_report(self) -> LabelFitReport | None:
        """Crossing label leaders drawn anyway: the fit report."""
        return label_fit_report(self.drawing)


@dataclass(frozen=True)
class PreparedCircularBatchRequest:
    """A validated Circular batch with one prepared diagram per record."""

    request: CircularBatchRequest
    records: tuple[SeqRecord, ...]
    items: tuple[PreparedDiagramRequest, ...]
    inputs: PreparedDiagramInputs | None = None
    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()
    protein_identity_manifest: Mapping[str, Any] | None = None

    @property
    def mode(self) -> Literal["circular"]:
        return "circular"


@dataclass(frozen=True)
class CircularBatchRenderResult:
    """Saved outputs for an explicit separate-diagram Circular batch.

    A Session replay fills ``migration_report`` with the Session's evidence,
    as it does on ``RequestRenderResult`` for one diagram.
    """

    request: CircularBatchRequest
    records: tuple[SeqRecord, ...]
    items: tuple[RequestRenderResult, ...]
    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()
    protein_identity_manifest: Mapping[str, Any] | None = None
    migration_report: SessionMigrationReport = field(
        default_factory=SessionMigrationReport
    )

    @property
    def mode(self) -> Literal["circular"]:
        return "circular"

    @property
    def drawings(self) -> tuple[Drawing, ...]:
        return tuple(item.drawing for item in self.items)

    @property
    def output_paths(self) -> tuple[Path, ...]:
        return tuple(path for item in self.items for path in item.output_paths)

    @property
    def interactive_contexts(self) -> tuple[InteractiveSvgContext | None, ...]:
        return tuple(item.interactive_context for item in self.items)


def _validated_plan_records(
    value: Sequence[SeqRecord],
    *,
    expected_count: int,
) -> tuple[SeqRecord, ...]:
    if isinstance(value, (str, bytes)):
        raise ValidationError(
            "Request plan records must be a sequence of SeqRecord values."
        )
    try:
        records = tuple(value)
    except TypeError as exc:
        raise ValidationError(
            "Request plan records must be a sequence of SeqRecord values."
        ) from exc
    if not records or not all(isinstance(record, SeqRecord) for record in records):
        raise ValidationError(
            "Request plan records must be a non-empty sequence of SeqRecord values."
        )
    if len(records) != expected_count:
        raise ValidationError(
            "Request plan record count must match the request record count."
        )
    return records


def _initialize_plan_display_context(plan) -> None:
    """Keep direct constructors and planner-produced context record-aligned."""
    if not plan.transforms:
        collection = (
            ResolvedRecordCollection(plan.records, plan.provenance)
            if plan.provenance
            else _coerce_resolved_collection(plan.request, plan.records)
        )
        object.__setattr__(plan, "provenance", collection.provenance)
        object.__setattr__(plan, "displays", collection.displays)
        object.__setattr__(plan, "transforms", collection.transforms)
    if any(len(values) != len(plan.records) for values in (
        plan.provenance, plan.displays, plan.transforms,
    )):
        raise ValidationError("Plan display context must align with its records.")
    if plan.resolved_annotations is None:
        object.__setattr__(plan, "resolved_annotations", resolve_annotations(
            plan.request.options.annotations, plan.records,
            mode=plan.mode, record_transforms=plan.transforms,
            feature_bindings=plan.inputs.feature_bindings if plan.inputs is not None else None,
        ))


@dataclass(frozen=True)
class CircularRequestPlan:
    """Normalized records, layout, and builder choice for a Circular request."""

    request: CircularDiagramRequest
    records: tuple[SeqRecord, ...]
    layout: CircularMultiRecordOptions | None
    precomputed_depth_track_specs: tuple[DepthTrackSpec, ...] | None = None
    precomputed_depth_track_count: int | None = None
    inputs: PreparedDiagramInputs | None = None
    provenance: tuple[ResolvedRecordProvenance, ...] = ()
    displays: tuple[ResolvedRecordDisplay, ...] = ()
    transforms: tuple[RecordDisplayTransform, ...] = ()
    resolved_annotations: ResolvedAnnotationBundle | None = None
    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.request, CircularDiagramRequest):
            raise ValidationError(
                "CircularRequestPlan request must be CircularDiagramRequest."
            )
        object.__setattr__(
            self,
            "records",
            _validated_plan_records(
                self.records,
                expected_count=len(self.request.records),
            ),
        )
        if self.layout is not None and not isinstance(
            self.layout,
            CircularMultiRecordOptions,
        ):
            raise ValidationError(
                "CircularRequestPlan layout must be CircularMultiRecordOptions or None."
            )
        if (self.layout is None) != (self.request.layout is None):
            raise ValidationError(
                "CircularRequestPlan layout presence must match the request."
            )
        if self.provenance and len(self.provenance) != len(self.records):
            raise ValidationError(
                "CircularRequestPlan provenance must align with its records."
            )

        _initialize_plan_display_context(self)

    @property
    def mode(self) -> Literal["circular"]:
        return "circular"

    def preflight_outputs(self) -> None:
        """Validate every materialized output before diagram construction."""

        _preflight_render_output(self.request.output)

    def build(self, *, label_fit_warn_search_seconds: float = WARN_SEARCH_SECONDS) -> Drawing:
        """Build the drawing; crossing label leaders search rows for the given seconds."""

        shared_kwargs: dict[str, Any] = (
            {"_resolved_feature_inputs": self.inputs.features}
            if self.inputs is not None
            else {}
        )
        if self.inputs is not None and self.inputs.record_features:
            shared_kwargs["_resolved_record_features"] = self.inputs.record_features
        shared_kwargs["_resolved_annotations"] = self.resolved_annotations
        shared_kwargs["_label_fit_warn_search_seconds"] = label_fit_warn_search_seconds
        if self.layout is None:
            depth_kwargs: dict[str, Any] = dict(shared_kwargs)
            if self.precomputed_depth_track_specs is not None:
                depth_kwargs.update({
                    "_precomputed_depth_track_specs": self.precomputed_depth_track_specs,
                    "_precomputed_depth_track_count": self.precomputed_depth_track_count,
                })
            return build_circular_diagram(
                self.records[0],
                options=self.request.options,
                _record_transform=self.transforms[0],
                _record_presentation=self.provenance[0].presentation,
                **depth_kwargs,
            )
        return build_circular_multi_diagram(
            self.records,
            options=self.request.options,
            layout=self.layout,
            _record_transforms=self.transforms,
            _record_presentations=tuple(item.presentation for item in self.provenance),
            **shared_kwargs,
        )


def _without_depth_inputs(options: CircularDiagramOptions) -> CircularDiagramOptions:
    """Strip depth sources after a batch plan has normalized them once."""

    return replace(
        options,
        depth_tracks=None,
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


@dataclass(frozen=True)
class CircularBatchRequestPlan:
    """Normalized records and one single-diagram plan per Circular batch item."""

    request: CircularBatchRequest
    records: tuple[SeqRecord, ...]
    inputs: PreparedDiagramInputs | None = None
    provenance: tuple[ResolvedRecordProvenance, ...] = ()
    displays: tuple[ResolvedRecordDisplay, ...] = ()
    transforms: tuple[RecordDisplayTransform, ...] = ()
    resolved_annotations: ResolvedAnnotationBundle | None = None
    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.request, CircularBatchRequest):
            raise ValidationError(
                "CircularBatchRequestPlan request must be CircularBatchRequest."
            )
        object.__setattr__(
            self,
            "records",
            _validated_plan_records(
                self.records,
                expected_count=len(self.request.records),
            ),
        )
        if self.provenance and len(self.provenance) != len(self.records):
            raise ValidationError(
                "CircularBatchRequestPlan provenance must align with its records."
            )

        _initialize_plan_display_context(self)

    @property
    def mode(self) -> Literal["circular"]:
        return "circular"

    def preflight_outputs(self) -> None:
        """Validate every materialized batch output before diagram construction."""

        preflight_diagram_request_outputs((self.request,))

    def item_plans(self) -> tuple[CircularRequestPlan, ...]:
        assert self.resolved_annotations is not None  # filled by __post_init__
        plans: list[CircularRequestPlan] = []
        options = self.request.options
        normalized_depth = normalize_depth_tracks(
            self.records,
            depth_tracks=options.depth_tracks,
            depth_table=options.depth_table,
            depth_file=options.depth_file,
            depth_tables=options.depth_tables,
            depth_files=options.depth_files,
            depth_track_tables=options.depth_track_tables,
            depth_track_files=options.depth_track_files,
            depth_track_labels=options.depth_track_labels,
            depth_track_colors=options.depth_track_colors,
            depth_track_large_tick_intervals=options.depth_track_large_tick_intervals,
            depth_track_small_tick_intervals=options.depth_track_small_tick_intervals,
            depth_track_tick_font_sizes=options.depth_track_tick_font_sizes,
        )
        logical_depth_count = depth_track_count(normalized_depth)
        item_options = _without_depth_inputs(options) if normalized_depth is not None else options
        for index, (record, output) in enumerate(
            zip(self.records, self.request.outputs, strict=True)
        ):
            item_request = CircularDiagramRequest(
                records=(
                    RecordInput(
                        source=InMemoryRecordSource(record, source_features=self.provenance[index].source_features),
                        record_key=self.request.records[index].record_key,
                        display=self.request.records[index].display,
                    ),
                ),
                options=replace(
                    item_options,
                    feature_placements=tuple(
                        item for item in item_options.feature_placements
                        if item.record_key == self.provenance[index].record_key
                    ),
                    feature_overrides=tuple(
                        item for item in item_options.feature_overrides
                        if item.record_key == self.provenance[index].record_key
                    ),
                ) if item_options.feature_placements or item_options.feature_overrides
                else item_options,
                output=output,
                grouping="single",
            )
            plans.append(
                CircularRequestPlan(
                    request=item_request,
                    records=(record,),
                    layout=None,
                    resolved_annotations=replace(
                        self.resolved_annotations,
                        annotations=tuple(replace(item, record_index=0)
                            for item in self.resolved_annotations.annotations
                            if item.record_index == index),
                        warnings=tuple(item for item in self.resolved_annotations.warnings
                            if item.record_index == index),
                    ),
                    precomputed_depth_track_specs=(
                        tuple(normalized_depth[index])
                        if normalized_depth is not None
                        else None
                    ),
                    precomputed_depth_track_count=(
                        logical_depth_count if normalized_depth is not None else None
                    ),
                    inputs=(replace(
                        self.inputs,
                        record_features=self.inputs.record_features[index:index + 1],
                        feature_identity_notices=tuple(
                            notice for notice in self.inputs.feature_identity_notices
                            if notice.record_index == index
                        ),
                    ) if self.inputs is not None and (
                        self.inputs.record_features or self.inputs.feature_identity_notices
                    ) else self.inputs),
                    displays=(self.displays[index],),
                    transforms=(self.transforms[index],),
                    provenance=(
                        (self.provenance[index],)
                        if self.provenance
                        else ()
                    ),
                )
            )
        return tuple(plans)


@dataclass(frozen=True)
class LinearRequestPlan:
    """Normalized records, layout, and builder choice for a Linear request."""

    request: LinearDiagramRequest
    records: tuple[SeqRecord, ...]
    layout: LinearMultiRecordOptions | None
    inputs: PreparedDiagramInputs | None = None
    provenance: tuple[ResolvedRecordProvenance, ...] = ()
    displays: tuple[ResolvedRecordDisplay, ...] = ()
    transforms: tuple[RecordDisplayTransform, ...] = ()
    resolved_annotations: ResolvedAnnotationBundle | None = None
    # The Similarity alignment's anchor centres (one per record, or none).
    alignment_anchor_centers: tuple[float | None, ...] = ()
    # Where ``anchor_feature`` puts each record (the planner's decision), or None without one.
    feature_anchor: FeatureAnchorPlacement | None = None
    # Raw cache entries of the LOSATN / TLOSATX searches this plan resolved.
    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()
    # Each LOSATP search's database size (``_losatp_database_sizes``); None when no LOSATP search runs.
    losatp_database_sizes: GenomeDatabaseSizes | None = None
    # The request keeps its comparison rows in the search frame; the plan
    # projects them into the drawn record views once, which also checks their
    # source bindings, and the build draws these options (OV-399).
    _drawn_options: LinearDiagramOptions = field(
        init=False,
        repr=False,
        compare=False,
    )
    # Whether any record is drawn reverse-complemented; only then do typed
    # items move between the search and drawn frames.
    _reverses: bool = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.request, LinearDiagramRequest):
            raise ValidationError(
                "LinearRequestPlan request must be LinearDiagramRequest."
            )
        object.__setattr__(
            self,
            "records",
            _validated_plan_records(
                self.records,
                expected_count=len(self.request.records),
            ),
        )
        if self.layout is not None and not isinstance(
            self.layout,
            LinearMultiRecordOptions,
        ):
            raise ValidationError(
                "LinearRequestPlan layout must be LinearMultiRecordOptions or None."
            )
        if (self.layout is None) != (self.request.layout is None):
            raise ValidationError(
                "LinearRequestPlan layout presence must match the request."
            )
        if self.provenance and len(self.provenance) != len(self.records):
            raise ValidationError(
                "LinearRequestPlan provenance must align with its records."
            )

        _initialize_plan_display_context(self)
        if not self.alignment_anchor_centers:
            object.__setattr__(
                self,
                "alignment_anchor_centers",
                tuple(None for _ in self.records),
            )
        if len(self.alignment_anchor_centers) != len(self.records):
            raise ValidationError(
                "Linear alignment anchor centers must align with its records."
            )
        # The planner decides the anchor for the records it plans.
        assert self.feature_anchor is None or len(self.feature_anchor.statuses) == len(self.records)
        options = self.request.options
        # Rows are checked against their bindings; a typed item moves only on a
        # reversed record (a forward record draws its search frame).
        reverses = any(map(record_is_reverse_complemented, self.records))
        object.__setattr__(self, "_reverses", reverses)
        object.__setattr__(
            self,
            "_drawn_options",
            project_source_bound_results(
                options,
                comparison_record_bindings(self.records, self.provenance, search_frame=False),
            ).options
            if has_source_bound_results(
                options.linear_comparisons, options.protein_comparisons,
                options.orthogroups if reverses else None,
                options.collinearity_blocks if reverses else None,
            )
            else options,
        )

    @property
    def mode(self) -> Literal["linear"]:
        return "linear"

    def preflight_outputs(self) -> None:
        """Validate every materialized output before diagram construction."""

        _preflight_render_output(self.request.output)

    def build(
        self,
        *,
        losatp_cache: LosatpCacheManager | None = None,
        protein_extraction: ProteinExtractionResult | None = None,
    ) -> LinearDiagramBuildResult:
        options = self._drawn_options
        kwargs: dict[str, Any] = {"options": options}
        if self.layout is not None:
            kwargs["layout"] = self.layout
        if losatp_cache is not None:
            kwargs["losatp_cache"] = losatp_cache
        if protein_extraction is not None:
            kwargs["protein_extraction"] = protein_extraction
        if self.inputs is not None:
            kwargs["_resolved_feature_inputs"] = self.inputs.features
            if self.inputs.record_features:
                kwargs["_resolved_record_features"] = self.inputs.record_features
        kwargs["_resolved_annotations"] = self.resolved_annotations
        kwargs["_record_transforms"] = self.transforms
        kwargs["_record_presentations"] = tuple(item.presentation for item in self.provenance)
        kwargs["similarity_alignment"] = self.request.similarity_alignment
        kwargs["_alignment_anchor_centers"] = self.alignment_anchor_centers
        kwargs["_feature_anchor"] = self.feature_anchor
        built = build_linear_diagram_result(self.records, **kwargs)
        if not isinstance(built, LinearDiagramBuildResult):
            if not isinstance(built, Drawing):
                raise ValidationError(
                    "The Linear diagram builder returned an unsupported result."
                )
            built = LinearDiagramBuildResult(
                drawing=built,
                metadata=LinearDiagramMetadata(
                    protein_comparisons=(
                        tuple(options.protein_comparisons)
                        if options.protein_comparisons is not None
                        else None
                    ),
                    linear_comparisons=tuple(options.linear_comparisons or ()),
                    orthogroups=options.orthogroups,
                    collinearity_result=(
                        options.collinearity_blocks
                        if isinstance(options.collinearity_blocks, CollinearityResult)
                        else None
                    ),
                ),
            )
        metadata = built.metadata
        if not self._reverses or not has_source_bound_results(
            metadata.linear_comparisons, metadata.protein_comparisons,
            metadata.orthogroups, metadata.collinearity_result,
        ):
            return built
        # The drawn results go back to each record's search frame when a
        # Session stores them (OV-1145); the bindings are made once, here.
        # A forward record draws its search frame, so it needs none.
        return replace(built, metadata=replace(
            metadata,
            search_frame_bindings=comparison_record_bindings(
                self.records, self.provenance, search_frame=True,
            ),
        ))


DiagramRequestPlan: TypeAlias = (
    CircularRequestPlan | CircularBatchRequestPlan | LinearRequestPlan
)


@dataclass(frozen=True)
class _PreparedLinearArtifacts:
    cache: LosatpCacheManager | None
    extraction: ProteinExtractionResult | None
    nucleotide_entries: tuple[Mapping[str, Any], ...]
    passthrough_derived_entries: tuple[Mapping[str, Any], ...]
    source_mode: str


def _linear_request_uses_comparisons(request: LinearDiagramRequest) -> bool:
    options = request.options
    return bool(
        options.blast_files
        or options.linear_comparisons
        or options.comparison_table_file
        or options.protein_comparisons
        or options.collinearity_blocks
        or options.losat_search is not None
    )


def _prepare_diagram_inputs(request: DiagramRequest) -> PreparedDiagramInputs:
    if not isinstance(
        request,
        (CircularDiagramRequest, CircularBatchRequest, LinearDiagramRequest),
    ):
        raise ValidationError("Unsupported diagram request type.")
    options = request.options
    color_table = _color_table(options)
    feature_visibility_table = _visibility_table(options)
    colors = options.colors
    default_colors = colors.default_colors if colors is not None else None
    palette = colors.default_colors_palette if colors is not None else "default"
    load_comparison = isinstance(request, LinearDiagramRequest) and _linear_request_uses_comparisons(request)
    has_gff_source = any(
        isinstance(record_input.source, GffFastaInputSource)
        for record_input in request.records
    )
    candidate_features, keep_all_features = (
        resolve_candidate_feature_types(
            options.selected_features_set or DEFAULT_SELECTED_FEATURES,
            color_table=color_table,
            feature_visibility_table=feature_visibility_table,
        )
        if has_gff_source
        else (set(options.selected_features_set or DEFAULT_SELECTED_FEATURES), False)
    )
    comparison_sequences = (
        _ComparisonSequenceSources(options.conservation_genome_files)
        if isinstance(options, CircularDiagramOptions)
        and options.conservation_sequence_files
        else None
    )
    return PreparedDiagramInputs(
        features=resolve_feature_inputs(
            color_table=color_table,
            default_colors=default_colors,
            feature_visibility_table=feature_visibility_table,
            default_color_overrides=(
                _default_color_overrides(colors)
                if default_colors is None and colors is not None
                else None
            ),
            default_colors_palette=palette,
            load_comparison_colors=load_comparison,
        ),
        gff_candidate_features=tuple(sorted(candidate_features)),
        gff_keep_all_features=keep_all_features,
        comparison_sequences=comparison_sequences,
    )


def _with_prepared_colors(request: _RequestT, inputs: PreparedDiagramInputs) -> _RequestT:
    """Record the canonical table actually compiled for this request."""
    if inputs.features.color_table is None:
        return request
    colors = replace(
        request.options.colors or ColorOptions(),
        color_table=inputs.features.color_table,
        color_table_file=None,
    )
    return replace(request, options=replace(request.options, colors=colors))


def _gff_types_named_by_visibility_rows(
    rows: Sequence[FeatureOverride], collection: ResolvedRecordCollection,
) -> set[str]:
    """Return the types of the GFF3 source features that identity rows set the
    Feature visibility of: a row turned on is drawn, and a feature with its own
    visibility stays in the Web Features list so it can be shown again (R-5).
    """
    shown: dict[str, set[str]] = {}
    for row in rows:
        if row.feature_visibility is not None:
            shown.setdefault(row.record_key, set()).add(row.biological_feature_id)
    return {
        entry.feature_type
        for item in collection.provenance
        if item.source_kind == "gff_fasta" and item.record_key in shown
        for entry in item.source_feature_catalog
        if entry.biological_feature_id in shown[item.record_key]
    } if shown else set()


def _source_catalogs(
    collection: ResolvedRecordCollection,
) -> tuple[tuple[SourceFeatureIdentity, ...], ...]:
    return tuple(item.source_feature_catalog for item in collection.provenance)


class _IdentityTableContext(TypedDict):
    records: tuple[SeqRecord, ...]
    record_keys: tuple[str, ...]
    source_record_ids: tuple[str, ...]
    source_catalogs: tuple[tuple[SourceFeatureIdentity, ...], ...]


def _identity_table_context(collection: ResolvedRecordCollection) -> _IdentityTableContext:
    """The records an identity table names: ``#index`` and record IDs select among them."""
    return {
        "records": collection.records,
        "record_keys": tuple(item.record_identity_key for item in collection.provenance),
        "source_record_ids": tuple(item.source_record_id for item in collection.provenance),
        "source_catalogs": _source_catalogs(collection),
    }


def _materialize_identity_tables(
    request: _RequestT, collection: ResolvedRecordCollection,
) -> _RequestT:
    """Replace feature placement and feature override tables with their exact rows."""
    options = request.options
    placement_table = (options.feature_placement_table if options.feature_placement_table is not None
                       else options.feature_placement_table_file)
    override_table = (options.feature_override_table if options.feature_override_table is not None
                      else options.feature_override_table_file)
    if placement_table is None and override_table is None:
        return request
    context = _identity_table_context(collection)
    changes: dict[str, Any] = {}
    if placement_table is not None:
        changes.update(
            feature_placements=read_feature_placement_table(
                placement_table,
                mode="linear" if isinstance(request, LinearDiagramRequest) else "circular",
                **context,
            ),
            feature_placement_table=None, feature_placement_table_file=None,
        )
    if override_table is not None:
        changes.update(
            feature_overrides=read_feature_override_table(override_table, **context),
            feature_override_table=None, feature_override_table_file=None,
        )
    return replace(request, options=replace(options, **changes))


def read_request_feature_override_table(
    request: DiagramRequest,
    table: DataFrame | str | Path,
    *,
    unmatched: list[int],
) -> tuple[FeatureOverride, ...]:
    """Resolve a feature override table against the records ``request`` loads.

    The Web's Load Feature Edits TSV reads a table against its committed request
    (design Q4 6.4). Rows that select no record or no feature are left out and
    listed in ``unmatched`` (Owner Q3 = A); other defects fail as on the CLI.
    """
    collection = _coerce_resolved_collection(
        request, _normalize_request_records(request, _prepare_diagram_inputs(request)),
    )
    return read_feature_override_table(table, **_identity_table_context(collection), unmatched=unmatched)


def _load_request_records(
    request: _RequestT,
    inputs: PreparedDiagramInputs,
    parsed_sources: dict[tuple[object, ...], ParsedRecordInputs] | None = None,
) -> tuple[_RequestT, ResolvedRecordCollection]:
    """Load the records and turn identity tables into exact rows against them.

    Tables resolve before the GFF3 check below, so a table row and an exact row
    load the same features. ``parsed_sources`` receives the source parses for a
    caller that reads the sources again (the ring subject database).
    """
    inputs = replace(inputs, parsed_sources={} if parsed_sources is None else parsed_sources)
    collection = _coerce_resolved_collection(request, _normalize_request_records(request, inputs))
    request = _materialize_identity_tables(request, collection)
    shown_types: set[str] = _gff_types_named_by_visibility_rows(request.options.feature_overrides, collection)
    if inputs.gff_keep_all_features or shown_types <= set(inputs.gff_candidate_features):
        return request, collection
    # A row can set the visibility of a GFF3 feature whose type the type filter
    # dropped. Load that type as well, from the parse the first load made.
    return request, _coerce_resolved_collection(request, _normalize_request_records(request, replace(
        inputs,
        gff_candidate_features=tuple(sorted({*inputs.gff_candidate_features, *shown_types})),
    )))


def _render_source_loaders() -> tuple[GenBankLoader, GffFastaLoader, FastaLoader]:
    """The loaders a request render reads its record sources with (read when called).

    The prepared parsed-source key names no loader, so every build under it
    takes them from here: the render, its warm-up, and a Worker helper's
    lookup (:func:`load_render_source_records`).
    """

    return load_gbks, parse_gff_fasta, load_fasta


def load_render_source_records(source: RecordInputSource) -> ParsedRecordInputs:
    """Parse one record source as a request render does, under its prepared key when staged."""

    genbank_loader, gff_loader, fasta_loader = _render_source_loaders()
    return _load_source_records(
        source, genbank_loader=genbank_loader, gff_loader=gff_loader, fasta_loader=fasta_loader
    )


def _normalize_request_records(
    request: DiagramRequest,
    inputs: PreparedDiagramInputs,
) -> ResolvedRecordCollection:
    def resolve() -> ResolvedRecordCollection:
        genbank_loader, gff_loader, fasta_loader = _render_source_loaders()
        return resolve_record_inputs(
            request.records,
            record_options=request.record_options,
            gff_candidate_features=inputs.gff_candidate_features,
            gff_keep_all_features=inputs.gff_keep_all_features,
            genbank_loader=genbank_loader,
            gff_loader=gff_loader,
            fasta_loader=fasta_loader,
            parsed_sources=inputs.parsed_sources,
        )

    if not prepared_input_cache_active():
        return resolve()
    source_specs = tuple(
        _prepared_source_cache_spec(record_input.source)
        for record_input in request.records
    )
    if any(spec is None for spec in source_specs):
        return resolve()
    typed_source_specs = tuple(
        spec for spec in source_specs if spec is not None
    )
    seen_source_keys: set[Hashable] = set()
    for record_input, (source_key, _identities) in zip(
        request.records,
        typed_source_specs,
        strict=True,
    ):
        if source_key in seen_source_keys:
            continue
        seen_source_keys.add(source_key)
        load_render_source_records(record_input.source)
    resource_identities = frozenset(
        identity
        for _key, identities in typed_source_specs
        for identity in identities
    )
    # A GFF3 source is cached with every type; the type filter shapes what resolves from it.
    gff_types = (
        _gff_feature_types(inputs.gff_candidate_features, inputs.gff_keep_all_features)
        if any(isinstance(item.source, GffFastaInputSource) for item in request.records)
        else frozenset()
    )
    gff_type_key = None if gff_types is None else tuple(sorted(gff_types))
    input_keys = tuple(_record_input_preparation_key(item) for item in request.records)
    collection_key = _record_collection_preparation_key(request.record_options)
    key = (
        "resolved-records-v1",
        tuple(source_key for source_key, _identities in typed_source_specs),
        input_keys,
        collection_key,
        gff_type_key,
    )
    collection = get_or_build_resolved_records(key, resource_identities, resolve)
    # Per-record identity (design 2.1, M1): the elements of the key above that
    # shape one record, so a record's render stages survive a change to another
    # record. Collection regions apply by the place in the whole collection,
    # which every input shapes, so with collection regions a record has the
    # whole key and its place. Labels and subtitles shape no record.
    # The source feature catalog of a record is resolved with it.
    default_collection = not request.record_options.regions
    for record, provenance in zip(collection.records, collection.provenance, strict=True):
        input_index = provenance.input_index
        record_identity = (
            (
                RECORD_IDENTITY_TAG,
                typed_source_specs[input_index][0],
                input_keys[input_index],
                input_index,
                provenance.source_record_index,
                gff_type_key if isinstance(request.records[input_index].source, GffFastaInputSource) else None,
                # Each window of a record is its own drawn record (D-01).
                *(() if provenance.window is None else (provenance.window.index,)),
            )
            if default_collection
            else (RECORD_IDENTITY_TAG, key, provenance.resolved_index)
        )
        register_prepared_value_identity(record, record_identity)
        register_prepared_value_identity(
            provenance.source_features, ("source-feature-catalog", record_identity)
        )
    return collection


def _selector_preparation_key(value: Any) -> Hashable:
    if value is None:
        return None
    return (value.record_id, value.record_index)


def _region_preparation_key(value: Any) -> Hashable:
    if value is None:
        return None
    return (
        value.file_selector,
        value.record_id,
        value.record_index,
        int(value.start),
        int(value.end),
        bool(value.reverse_complement),
    )


def _presentation_preparation_key(value: RecordPresentation) -> Hashable:
    # Label and subtitle shape no record: the provenance carries them
    # (``effective_record_presentation``), so they are not in the key.
    return (
        bool(value.reverse_complement),
        value.grid_row,
        value.grid_column,
    )


def _windows_preparation_key(value: RecordInput) -> Hashable:
    # A record with one window keys as the region it was before windows (D-09), so its key does not change.
    if not value.windows:
        return None
    if len(value.windows) == 1:
        return _region_preparation_key(value.windows[0].region())
    return tuple((window.start, window.end, window.reverse_complement, window.column) for window in value.windows)


def _record_input_preparation_key(value: RecordInput) -> Hashable:
    return (
        _selector_preparation_key(value.selector),
        value.cardinality.value,
        _windows_preparation_key(value),
        _presentation_preparation_key(value.presentation),
        value.record_key,
    )


def _record_collection_preparation_key(
    value: RecordCollectionOptions,
) -> Hashable:
    return tuple(_region_preparation_key(region) for region in value.regions)


def _relabelled(presentation: RecordPresentation, labelled: RecordPresentation) -> RecordPresentation:
    """``presentation`` with the label and subtitle of ``labelled``."""

    if (presentation.label, presentation.subtitle) == (labelled.label, labelled.subtitle):
        return presentation
    return replace(presentation, label=labelled.label, subtitle=labelled.subtitle)


def _coerce_resolved_collection(
    request: DiagramRequest,
    value: ResolvedRecordCollection | Sequence[SeqRecord],
) -> ResolvedRecordCollection:
    """Keep private test seams compatible while planners consume provenance."""

    if isinstance(value, ResolvedRecordCollection):
        # Cached records represent selection/RC/crop, not display intent. Rebind
        # transient context (display options, label and subtitle) without
        # changing records, source membership, keys or places: the resolver
        # placed each record (a File card's records by ``numbered_places``).
        rebound = tuple(
            replace(
                item,
                display=request.records[item.input_index].display,
                presentation=_relabelled(
                    item.presentation,
                    effective_record_presentation(
                        request.records[item.input_index].presentation, request.record_options, item.place
                    ),
                ),
            )
            for item in value.provenance
        )
        if rebound == value.provenance:
            return value
        return ResolvedRecordCollection(value.records, rebound)
    records = tuple(value)
    if len(records) != len(request.records):
        raise ValidationError(
            "A record resolver without provenance must return one record per input."
        )
    provenance: list[ResolvedRecordProvenance] = []
    for index, (record, record_input) in enumerate(
        zip(records, request.records, strict=True)
    ):
        source_kind, source_paths = _source_details(record_input.source)
        provenance.append(
            ResolvedRecordProvenance(
                resolved_index=index,
                input_index=index,
                source_record_index=0,
                source_record_count=1,
                source_record_id=str(record.id),
                source_kind=source_kind,
                source_paths=source_paths,
                record_key=record_input.record_key or f"record-{index + 1}",
                cardinality=record_input.cardinality,
                selector=record_input.selector,
                region=None,
                presentation=effective_record_presentation(
                    record_input.presentation, request.record_options, index
                ),
                display=record_input.display,
                source_length=record.annotations.get("gbdraw_source_length"),
                detected_topology=detected_topology(record, source_kind),
                # The records are the caller's: capture without stamping them.
                source_features=capture_source_features(record, stamp=False),
            )
        )
    return ResolvedRecordCollection(records, tuple(provenance))


def resolve_request_records(request: DiagramRequest) -> ResolvedRecordCollection:
    """Load the records ``request`` draws, with the provenance of each."""

    # Every request mode loads its records the same way.
    return _load_request_records(cast(Any, request), _prepare_diagram_inputs(request))[1]


def normalize_request_records(request: DiagramRequest) -> tuple[SeqRecord, ...]:
    """Resolve typed record inputs according to their explicit cardinality."""

    return tuple(detached_records(resolve_request_records(request).records))


def _materialized_record_inputs(
    collection: ResolvedRecordCollection,
) -> tuple[RecordInput, ...]:
    """Project resolved records to canonical exact-one inputs."""

    materialized: list[RecordInput] = []
    for record, provenance in zip(
        collection.records,
        collection.provenance,
        strict=True,
    ):
        source_presentation = provenance.presentation
        materialized.append(
            RecordInput(
                source=InMemoryRecordSource(record, source_features=provenance.source_features),
                presentation=RecordPresentation(
                    label=source_presentation.label,
                    subtitle=source_presentation.subtitle,
                    reverse_complement=False,
                    grid_row=source_presentation.grid_row,
                    grid_column=source_presentation.grid_column,
                ),
                record_key=provenance.record_key,
                display=provenance.display,
            )
        )
    return tuple(materialized)


def _is_materialized_exact_one_request(request: DiagramRequest) -> bool:
    return (
        request.record_options == RecordCollectionOptions()
        and all(
            isinstance(record_input.source, InMemoryRecordSource)
            and record_input.cardinality is RecordCardinality.EXACTLY_ONE
            and record_input.selector is None
            and not record_input.windows
            and not record_input.presentation.reverse_complement
            for record_input in request.records
        )
    )


def _resolve_first_record_output(
    output: RenderOutputRequest,
    records: Sequence[SeqRecord],
) -> RenderOutputRequest:
    if not output.resolve_prefix_from_first_record:
        return output
    return replace(
        output,
        output_prefix=resolve_implicit_record_output_prefix(records[0].id),
        output_directory=None,
        resolve_prefix_from_first_record=False,
    )


def _warn_circular_topologies(records: Sequence[SeqRecord]) -> None:
    for record in records:
        topology = getattr(record, "annotations", {}).get("topology")
        if topology == "linear":
            logger.warning(
                "The annotation indicates that record %s is linear. "
                "Are you sure you want to visualize it as circular?",
                record.id,
            )
        elif topology is not None and topology != "circular":
            logger.warning(
                "Topology information not available for %s.",
                record.id,
            )


def grid_placement_order(
    positioned: Sequence[tuple[int, RecordPresentation]],
) -> list[tuple[int, RecordPresentation]]:
    """Records with a grid row (or, in a projected Linear request, drawn units), by index, in drawing order: by row,
    then by column (one without a column by its index), then by index. The one owner of that order: the layouts below
    place by it, and the Web's card numbers a row's windows by it (``plan_card_windows`` of
    ``gbdraw.web_support.linear_units``)."""

    return sorted(
        positioned,
        key=lambda item: (
            int(item[1].grid_row or 0),
            int(item[1].grid_column) if item[1].grid_column is not None else item[0],
            item[0],
        ),
    )


def _layout_with_record_placements(
    request: CircularDiagramRequest,
) -> CircularMultiRecordOptions | None:
    layout = request.layout
    if layout is None or layout.multi_record_positions:
        return layout
    positioned = [
        (index, record_input.presentation)
        for index, record_input in enumerate(request.records)
        if record_input.presentation.grid_row is not None
    ]
    if not positioned:
        return layout
    positioned = grid_placement_order(positioned)
    positions = tuple(
        f"#{index + 1}@{presentation.grid_row}"
        for index, presentation in positioned
    )
    return replace(layout, multi_record_positions=positions)


def _linear_layout_with_record_placements(
    request: LinearDiagramRequest,
) -> LinearMultiRecordOptions | None:
    layout = request.layout
    if layout is not None and layout.multi_record_positions:
        return _unit_positions(request, layout)
    positioned = [
        (index, record_input.presentation)
        for index, record_input in enumerate(request.records)
        if record_input.presentation.grid_row is not None
    ]
    if not positioned:
        return _windows_share_rows(request, layout)
    if layout is None:  # pragma: no cover - LinearDiagramRequest refuses grid rows without a layout.
        raise ValidationError("Grid placement requires LinearMultiRecordOptions.", diagnostic={"code": "INPUT_INVALID"})
    rows = _request_unit_rows(request)
    positions = tuple(f"#{index + 1}@{rows[index] + 1}" for index, _presentation in grid_placement_order(positioned))
    return replace(layout, multi_record_positions=positions)


def _request_unit_rows(request: LinearDiagramRequest) -> tuple[int, ...]:
    """The row of each record of a projected request that names no row positions (``drawn_unit_rows``): the one rule
    the Web's unit plan reads too (``plan_linear_units``)."""

    return drawn_unit_rows(
        [record_input.presentation.grid_row for record_input in request.records],
        [drawn_record_identity(record_input) for record_input in request.records],
    )


def _unit_positions(request: LinearDiagramRequest, layout: LinearMultiRecordOptions) -> LinearMultiRecordOptions:
    """Given row positions name records (``#n`` counts records, a record ID names its record, I-6); each becomes the
    positions of the record's drawn units, in window order, at that row.

    ``request`` is projected (one in-memory record per drawn unit). Without a record of several windows records and
    units are one, and the positions stay as given.
    """

    keys = [drawn_record_identity(record_input) for record_input in request.records]
    if len(set(keys)) == len(keys):
        return layout
    records = [cast(InMemoryRecordSource, record_input.source).record for record_input in request.records]
    positions: list[str] = []
    for token in layout.multi_record_positions or ():
        selector, row = parse_record_row_position(str(token))
        try:
            units = select_record_units(records, record_order_selector(selector, records))
        except ValueError as exc:
            raise ValidationError(f"multi_record_position selector '{selector}': {exc}") from exc
        positions.extend(f"#{unit + 1}@{row}" for unit in units)
    return replace(layout, multi_record_positions=tuple(positions))


def _with_window_layout(request: LinearDiagramRequest) -> LinearDiagramRequest:
    """A projected request whose windows share rows has a Linear layout, the default one when it set none (D-01)."""

    keys = [drawn_record_identity(record_input) for record_input in request.records]
    if request.layout is not None or len(set(keys)) == len(keys):
        return request
    return replace(request, layout=LinearMultiRecordOptions())


def _windows_share_rows(
    request: LinearDiagramRequest, layout: LinearMultiRecordOptions | None
) -> LinearMultiRecordOptions | None:
    """Without grid rows, each record takes a row and the windows of a record share it, in window order (D-01).

    ``request`` is projected: one in-memory record per drawn unit. A request without a record of several windows
    keeps its layout (one record per row).
    """

    keys = [drawn_record_identity(record_input) for record_input in request.records]
    if layout is None or len(set(keys)) == len(keys):
        return layout
    rows = _request_unit_rows(request)
    return replace(
        layout,
        multi_record_positions=tuple(f"#{index + 1}@{row + 1}" for index, row in enumerate(rows)),
    )


def _materialize_feature_inputs(
    request: DiagramRequest, collection: ResolvedRecordCollection, inputs: PreparedDiagramInputs,
) -> PreparedDiagramInputs:
    """Resolve every identity-addressed input once per Generate (design Q4, 3.1).

    ``_load_request_records`` has already turned identity tables into exact rows.
    """
    options = request.options
    targets = tuple(
        FeatureIdentity(target.record_key, target.biological_feature_id)
        for annotation_set in (options.annotations.sets if options.annotations else ())
        for annotation in annotation_set.annotations
        if isinstance(target := annotation.target, FeatureIdentitySpan)
    )
    if not options.feature_placements and not options.feature_overrides and not targets:
        return inputs
    resolution = resolve_record_feature_inputs(
        records=collection.records,
        record_keys=tuple(item.record_identity_key for item in collection.provenance),
        source_catalogs=_source_catalogs(collection), placements=options.feature_placements,
        mode="linear" if isinstance(request, LinearDiagramRequest) else "circular",
        feature_overrides=options.feature_overrides,
        target_identities=targets,
        selected_features=options.selected_features_set or DEFAULT_SELECTED_FEATURES,
        feature_visibility_rules=inputs.features.feature_visibility_rules,
        specific_color_rules=inputs.features.specific_color_rules,
        feature_shapes=options.feature_shapes,
    )
    return replace(
        inputs,
        record_features=(
            resolution.records if options.feature_placements or options.feature_overrides else ()
        ),
        feature_identity_notices=resolution.notices,
        feature_bindings={
            (identity.record_key, identity.biological_feature_id): resolution.bindings[identity]
            for identity in targets
        },
    )


def _resolve_ring_losat(
    options: CircularDiagramOptions,
    record_inputs: Sequence[RecordInput],
    parsed_sources: Mapping[tuple[object, ...], ParsedRecordInputs],
    inputs: PreparedDiagramInputs,
    drawn_records: Sequence[SeqRecord],
) -> tuple[CircularDiagramOptions, tuple[Mapping[str, Any], ...]]:
    """Replace ring LOSATN / TLOSATX intent with ring rows (design 3.3).

    The subject database is every subject-side source file of the request,
    whole (Owner D-25), from the parses the record load made; the ring rows are
    in the search frame of ``drawn_records`` (PD-OI-073). Each comparison
    genome is read once through the request's memoized reader.
    """

    if options.losat_search is None:
        return options, ()

    def file_records(source: GenBankInputSource | FastaInputSource) -> tuple[SeqRecord, ...]:
        parsed = parsed_sources.get(_source_key(source))
        # A request served from the prepared resolved-record cache made no parse
        # here: the prepared parsed source holds it.
        return (parsed if parsed is not None else load_render_source_records(source)).records

    sources = inputs.comparison_sequences or _ComparisonSequenceSources(
        options.conservation_genome_files
    )

    def load_sequences() -> tuple[ComparisonGenome, ...]:
        return tuple(genome for genome in sources.files() if genome is not None)

    return resolve_circular_conservation_losat(
        options,
        subject_sources=ring_subject_sources(
            (item.source for item in record_inputs), file_records=file_records
        ),
        drawn_records=drawn_records,
        load_sequences=load_sequences,
    )


_CircularRequestT = TypeVar("_CircularRequestT", CircularDiagramRequest, CircularBatchRequest)


def _with_origin_feature(
    request: _CircularRequestT,
    collection: ResolvedRecordCollection,
) -> tuple[_CircularRequestT, ResolvedRecordCollection, bool]:
    """Resolve ``origin_feature`` into the display start of each record it matches.

    The request then carries the starts in its records, not the option, so a
    Session stores the start the figure drew with. The flag says whether it
    rotated.
    """

    selector = request.options.origin_feature
    if not isinstance(selector, FeatureSelector):
        return request, collection, False
    return (
        replace(request, options=replace(request.options, origin_feature=None)),
        with_origin_feature_starts(collection, selector),
        True,
    )


def plan_circular_request(
    request: CircularDiagramRequest,
) -> CircularRequestPlan:
    """Normalize a Circular request into one explicit builder plan."""

    if not isinstance(request, CircularDiagramRequest):
        raise ValidationError("request must be CircularDiagramRequest.")
    with _request_render_diagnostic_phase("preparation"):
        resolved_options = _resolve_request_option_tables(
            resolve_circular_options(request.options),
            mode="circular",
        )
        unresolved_request = (
            request
            if resolved_options is request.options
            else replace(request, options=resolved_options)
        )
        inputs = _prepare_diagram_inputs(unresolved_request)
        unresolved_request = _with_prepared_colors(unresolved_request, inputs)
    parsed_sources: dict[tuple[object, ...], ParsedRecordInputs] = {}
    with _request_render_diagnostic_phase("recordLoad"):
        unresolved_request, collection = _load_request_records(unresolved_request, inputs, parsed_sources)
    with _request_render_diagnostic_phase("preparation"):
        unresolved_request, collection, rotated = _with_origin_feature(unresolved_request, collection)
        inputs = _materialize_feature_inputs(unresolved_request, collection, inputs)
        records = collection.records
        projected_request = replace(
            unresolved_request,
            records=_materialized_record_inputs(collection),
            output=_resolve_first_record_output(
                unresolved_request.output,
                records,
            ),
            record_options=RecordCollectionOptions(),
        )
        resolved_layout = _layout_with_record_placements(projected_request)
        materialized_request = (
            unresolved_request
            if not rotated
            and _is_materialized_exact_one_request(unresolved_request)
            and not unresolved_request.output.resolve_prefix_from_first_record
            and resolved_layout == unresolved_request.layout
            else projected_request
        )
        ring_options, losat_cache_entries = _resolve_ring_losat(
            materialized_request.options, unresolved_request.records, parsed_sources, inputs, records
        )
        if ring_options is not materialized_request.options:
            materialized_request = replace(materialized_request, options=ring_options)
        # Records on one circle are drawn as arcs of it: a linear record is what that layout is for.
        if materialized_request.layout is None or materialized_request.layout.record_layout != "one_circle":
            _warn_circular_topologies(records)
    return CircularRequestPlan(
        request=materialized_request,
        records=records,
        layout=resolved_layout,
        inputs=inputs,
        provenance=collection.provenance,
        displays=collection.displays,
        transforms=collection.transforms,
        losat_cache_entries=losat_cache_entries,
    )


def plan_circular_batch_request(
    request: CircularBatchRequest,
) -> CircularBatchRequestPlan:
    """Normalize an explicit separate-diagram Circular batch."""

    if not isinstance(request, CircularBatchRequest):
        raise ValidationError("request must be CircularBatchRequest.")
    with _request_render_diagnostic_phase("preparation"):
        resolved_options = _resolve_request_option_tables(
            resolve_circular_options(request.options),
            mode="circular",
        )
        unresolved_request = (
            request
            if resolved_options is request.options
            else replace(request, options=resolved_options)
        )
        inputs = _prepare_diagram_inputs(unresolved_request)
        unresolved_request = _with_prepared_colors(unresolved_request, inputs)
    parsed_sources: dict[tuple[object, ...], ParsedRecordInputs] = {}
    with _request_render_diagnostic_phase("recordLoad"):
        unresolved_request, collection = _load_request_records(unresolved_request, inputs, parsed_sources)
    with _request_render_diagnostic_phase("preparation"):
        unresolved_request, collection, rotated = _with_origin_feature(unresolved_request, collection)
        inputs = _materialize_feature_inputs(unresolved_request, collection, inputs)
        records = collection.records
        outputs = (
            unresolved_request.outputs
            if unresolved_request.outputs
            else resolve_circular_batch_outputs(
                # A batch request has explicit outputs or an output policy.
                cast("CircularBatchOutputPolicy", unresolved_request.output_policy),
                records,
            )
        )
        materialized_request = (
            unresolved_request
            if not rotated
            and _is_materialized_exact_one_request(unresolved_request)
            and unresolved_request.output_policy is None
            else replace(
                unresolved_request,
                records=_materialized_record_inputs(collection),
                outputs=outputs,
                output_policy=None,
                record_options=RecordCollectionOptions(),
            )
        )
        ring_options, losat_cache_entries = _resolve_ring_losat(
            materialized_request.options, unresolved_request.records, parsed_sources, inputs, records
        )
        if ring_options is not materialized_request.options:
            materialized_request = replace(materialized_request, options=ring_options)
        _warn_circular_topologies(records)
    return CircularBatchRequestPlan(
        request=materialized_request,
        records=records,
        inputs=inputs,
        provenance=collection.provenance,
        displays=collection.displays,
        transforms=collection.transforms,
        losat_cache_entries=losat_cache_entries,
    )


def plan_linear_request(
    request: LinearDiagramRequest,
) -> LinearRequestPlan:
    """Normalize a Linear request into one explicit builder plan.

    A :class:`SimilarityAlignmentReference` is resolved first: the requested
    orthogroup analysis runs once and the plan's request carries its result
    and the resolved alignment plan. Use :func:`render_request` or
    :func:`build_request_diagram` to keep that analysis's LOSATP artifacts.
    """

    if not isinstance(request, LinearDiagramRequest):
        raise ValidationError("request must be LinearDiagramRequest.")
    if isinstance(request.similarity_alignment, SimilarityAlignmentReference):
        # Resolving the reference keeps the request mode.
        request, _artifacts = cast(
            "tuple[LinearDiagramRequest, CurrentRequestArtifacts | None]",
            _resolve_similarity_alignment_reference(request, None),
        )
    loaded = _load_linear_request(request)
    unresolved_request, collection, prepared_inputs = loaded.request, loaded.collection, loaded.prepared_inputs
    inputs, projected_request, resolved_layout = loaded.inputs, loaded.projected_request, loaded.layout
    with _request_render_diagnostic_phase("preparation"):
        resolved_options = resolve_linear_options(
            decide_linear_wrap(
                projected_request.options,
                records=collection.records,
                layout=resolved_layout,
                similarity_alignment=projected_request.similarity_alignment,
            ),
            records=collection.records,
            layout=resolved_layout,
        )
        resolved_options, losat_cache_entries = _resolve_nucleotide_losat(
            resolved_options, collection, resolved_layout
        )
        oriented = _auto_oriented_collection(resolved_options, collection, resolved_layout)
        reoriented = oriented is not collection
        if resolved_options.auto_orient:
            # The request carries the decided reverse flags in its records, not
            # the option, so a Session stores the orientation the figure drew.
            resolved_options = replace(resolved_options, auto_orient=False)
        if reoriented:
            collection = oriented
        # One record-position alignment at most: the request refuses an
        # anchor_feature beside a similarity alignment. The anchor is read
        # after auto_orient, and the starts of the records it rotates are
        # stored in the request's records (a Session stores starts, never the option).
        anchor_feature = resolved_options.anchor_feature
        feature_anchor = None
        if isinstance(anchor_feature, FeatureSelector):
            collection, feature_anchor = resolve_anchor_feature(
                collection,
                anchor_feature,
                resolve_record_row_positions(
                    collection.records,
                    resolved_layout.multi_record_positions if resolved_layout is not None else None,
                )[1],
            )
        records_changed = reoriented or (feature_anchor is not None and feature_anchor.rotates)
        if records_changed:
            # The search-frame rows, the layout and the resolved options hold
            # for either orientation and start; the records and their features do not.
            inputs = _materialize_feature_inputs(unresolved_request, collection, prepared_inputs)
            projected_request = replace(
                projected_request, records=_materialized_record_inputs(collection)
            )
        materialized_request = (
            unresolved_request
            if not records_changed
            and _is_materialized_exact_one_request(unresolved_request)
            and resolved_layout == unresolved_request.layout
            and resolved_options is unresolved_request.options
            else replace(
                projected_request,
                options=resolved_options,
            )
        )
        # The reference was resolved to a plan at the top of this function.
        alignment_anchor_centers = project_similarity_alignment_centers(
            collection,
            cast("SimilarityAlignmentPlan | None", materialized_request.similarity_alignment),
        )
        # The plan checks the source bindings of the comparison rows when it is
        # built, so a Session writer rejects them; the request keeps the rows as
        # given (OV-399). Building it here keeps that projection in `preparation`.
        return LinearRequestPlan(
            request=materialized_request,
            records=collection.records,
            layout=resolved_layout,
            inputs=inputs,
            provenance=collection.provenance,
            displays=collection.displays,
            transforms=collection.transforms,
            alignment_anchor_centers=alignment_anchor_centers,
            feature_anchor=feature_anchor,
            losat_cache_entries=losat_cache_entries,
            losatp_database_sizes=(
                None if losatp_analysis_mode(resolved_options.losat_search) == "none"
                else _losatp_database_sizes(replace(loaded, collection=collection, inputs=inputs))
            ),
        )


@dataclass(frozen=True)
class _LoadedLinearRequest:
    """A Linear request with its records loaded, as its plan reads them before it resolves its options."""

    # The request with its option tables and colours resolved, and its record inputs as the record load wrote them.
    request: LinearDiagramRequest
    collection: ResolvedRecordCollection
    prepared_inputs: PreparedDiagramInputs
    # ``prepared_inputs`` with the identity-addressed inputs resolved on the collection.
    inputs: PreparedDiagramInputs
    projected_request: LinearDiagramRequest
    layout: LinearMultiRecordOptions | None
    # Each source file the record load parsed, by its source key.
    parsed_sources: Mapping[tuple[object, ...], ParsedRecordInputs]


def _load_linear_request(request: LinearDiagramRequest) -> _LoadedLinearRequest:
    """Load the records of a Linear request and resolve what its plan reads of them: its feature inputs and layout."""

    with _request_render_diagnostic_phase("preparation"):
        unresolved_request = replace(
            request,
            options=_resolve_request_option_tables(
                request.options,
                mode="linear",
            ),
        )
        inputs = _prepare_diagram_inputs(unresolved_request)
        unresolved_request = _with_prepared_colors(unresolved_request, inputs)
    parsed_sources: dict[tuple[object, ...], ParsedRecordInputs] = {}
    with _request_render_diagnostic_phase("recordLoad"):
        unresolved_request, collection = _load_request_records(unresolved_request, inputs, parsed_sources)
    with _request_render_diagnostic_phase("preparation"):
        projected_request = _with_window_layout(replace(
            unresolved_request,
            records=_materialized_record_inputs(collection),
            record_options=RecordCollectionOptions(),
        ))
        return _LoadedLinearRequest(
            request=cast(LinearDiagramRequest, unresolved_request),
            collection=collection,
            prepared_inputs=inputs,
            inputs=_materialize_feature_inputs(unresolved_request, collection, inputs),
            projected_request=projected_request,
            layout=_linear_layout_with_record_placements(projected_request),
            parsed_sources=parsed_sources,
        )


GenomeScaleProgram: TypeAlias = Literal["losatn", "tlosatx", "losatp"]
# The LOSAT programs whose searches pass an E-value scale of their genomes.
GENOME_SCALE_PROGRAMS: tuple[GenomeScaleProgram, ...] = ("losatn", "tlosatx", "losatp")


def linear_genome_scales(
    request: LinearDiagramRequest, program: GenomeScaleProgram, source_names: Mapping[str, str] | None = None
) -> tuple[tuple[str, ...], Callable[[int, int], int | None], Callable[[int, int], PerFileWarning | None]]:
    """The E-value scale each record-pair search of a Linear request's ``program`` search passes LOSAT, by drawn
    record indexes (a record drawn in windows is a record per window), with the key of the record each draws and
    the warning of a pair searched per file.

    LOSATN and TLOSATX pass the search space of the genome pair (``-searchsp``), LOSATP the proteins of the subject
    genome (``-dbsize``); ``None`` keeps LOSAT's own (PD-OI-018 revision 6). A genome pair whose search space is past
    what LOSAT reads is searched per file, with the warning the command line logs (Owner decision 2026-11-10);
    ``source_names`` names the files the warning names, by the name the request reads them under. They are the
    values the plan's own searches pass, from the same owners: the Web's Generate asks for them before it searches
    (``gbdraw.web_support.linear_units.plan_genome_scales``), so the CLI and the Web pass one value.
    """

    loaded = _load_linear_request(request)
    keys = tuple(item.record_identity_key for item in loaded.collection.provenance)
    if program == "losatp":
        return keys, _losatp_database_sizes(loaded), lambda _query, _subject: None
    rows, source_ids, source_lengths = _genome_inputs(loaded.collection, loaded.layout)
    spaces = linear_search_spaces(
        program, records=loaded.collection.records, rows_by_record=rows, source_ids=source_ids,
        source_lengths=source_lengths,
    )
    names = genome_file_names(
        rows, source_ids, lambda source: (source_names or {}).get(source_file_name(source), source_file_name(source))
    )
    return keys, spaces, lambda query, subject: spaces.per_file_warning(query, subject, names)


def _auto_oriented_collection(
    options: LinearDiagramOptions,
    collection: ResolvedRecordCollection,
    layout: LinearMultiRecordOptions | None,
) -> ResolvedRecordCollection:
    """Reverse-complement the records ``auto_orient`` decides from nucleotide hits.

    Reads the comparison files and the comparisons of ``options`` (LOSATN and
    TLOSATX results included) in the search frame, filtered by the match
    thresholds as the drawing filters them; protein comparisons are not read.
    :func:`auto_orient_reversed_records` decides; a record the request already
    reverses keeps its orientation. Returns ``collection`` itself when no
    record changes.
    """

    if not options.auto_orient:
        return collection
    records = collection.records
    _ordered, rows_by_record = resolve_record_row_positions(
        records, layout.multi_record_positions if layout is not None else None
    )
    comparisons = search_frame_linear_comparisons(
        tuple(read_comparison_resource(path) for path in options.blast_files or ()),
        None,
        tuple(options.linear_comparisons or ()),
        ComparisonThresholds(
            evalue=float(cast(float, options.evalue)),
            bitscore=float(cast(float, options.bitscore)),
            identity=float(cast(float, options.identity)),
            alignment_length=int(cast(int, options.alignment_length)),
        ),
        rows_by_record,
        len(set(rows_by_record)) < len(records),
    )
    reversed_indexes = auto_orient_reversed_records(
        comparisons, rows_by_record, tuple(map(record_is_reverse_complemented, records))
    )
    if not reversed_indexes:
        return collection
    for index in sorted(reversed_indexes):
        logger.info(
            "auto_orient reverse-complemented record #%d %r: its reverse hits "
            "against the row above weigh more than half of its hits.",
            index + 1, str(records[index].id),
        )
    return reverse_resolved_records(collection, reversed_indexes)


def _resolve_nucleotide_losat(
    options: LinearDiagramOptions,
    collection: ResolvedRecordCollection,
    layout: LinearMultiRecordOptions | None,
) -> tuple[LinearDiagramOptions, tuple[Mapping[str, Any], ...]]:
    """Replace LOSATN / TLOSATX intent with comparisons (design 3.4).

    A source file is one search job's side; an in-memory record is its own source
    (:func:`gbdraw.comparisons.losat_jobs.losat_source_ids`, as for LOSATP). A
    row's source files are one genome, whose E-value scale every search between
    two rows shares (OV-3143).
    """

    if options.losat_search is None or options.losat_search.program == "losatp":
        return options, ()
    rows_by_record, source_ids, source_lengths = _genome_inputs(collection, layout)
    records = [
        NucleotideLosatRecord(
            record=record,
            source_id=source_id,
            source_length=source_length,
            record_key=record_key,
            label=item.presentation.label or str(record.id),
            window=item.window.index + 1 if item.window is not None and item.window.is_one_of_several else None,
            input_index=item.input_index,
        )
        for item, record, source_id, source_length, record_key in zip(
            collection.provenance,
            collection.records,
            source_ids,
            source_lengths,
            losat_record_uids(collection.records),
            strict=True,
        )
    ]
    return resolve_linear_nucleotide_losat(options, records=records, rows_by_record=rows_by_record)


def _source_genome_lengths(collection: ResolvedRecordCollection) -> list[int]:
    """The genome length of each record's source (``source_genome_length``), read once per file (a GFF3/FASTA
    pair reads its FASTA file again)."""

    file_lengths: dict[tuple[str, ...], int] = {}

    def genome_length(item: ResolvedRecordProvenance, record: SeqRecord) -> int:
        loaded = item.source_file_length if item.source_file_length is not None else len(record)
        if item.source_kind == "memory":
            return loaded
        if item.source_paths not in file_lengths:
            file_lengths[item.source_paths] = source_genome_length(item.source_kind, item.source_paths, loaded)
        return file_lengths[item.source_paths]

    return [genome_length(item, record) for item, record in zip(collection.provenance, collection.records, strict=True)]


def _genome_inputs(
    collection: ResolvedRecordCollection, layout: LinearMultiRecordOptions | None
) -> tuple[list[int], list[Hashable], list[int]]:
    """What a genome of each drawn record is made of (OV-3143): its row, its source file, and the length of every
    record of that source (``_source_genome_lengths``); the genomes are the files linked through shared rows."""

    rows = resolve_record_row_positions(
        collection.records, layout.multi_record_positions if layout is not None else None
    )[1]
    return list(rows), list(losat_source_ids(collection.records)), _source_genome_lengths(collection)


def _losatp_database_sizes(loaded: _LoadedLinearRequest) -> GenomeDatabaseSizes:
    """Each LOSATP record-pair search's database size, its subject genome's proteins (OV-3143).

    A source counts every record of its file, uncropped, drawn or not (a record held from a file plan, its
    file's; an in-memory crop, the record it was cut from, as its genome length does), with the proteins LOSATP's
    extraction takes under the drawing's visibility rules and the request's feature edits of each drawn record,
    bound to the uncropped record, so an edit outside a crop or window counts too (an undrawn record has none):
    ``source_protein_letters``, read once per source.
    """

    collection, inputs = loaded.collection, loaded.inputs
    provenance = collection.provenance
    edits = normalize_feature_overrides(loaded.request.options.feature_overrides or ())

    def source_of(item: ResolvedRecordProvenance) -> Hashable:
        return item.source_paths or ("memory", item.input_index)

    letters: dict[Hashable, int] = {}
    for item in provenance:
        key = source_of(item)
        if key in letters:
            continue
        source = loaded.request.records[item.input_index].source
        if item.source_kind != "memory" and isinstance(source, InMemoryRecordSource):
            source = _file_source(item.source_kind, item.source_paths)
        parsed = loaded.parsed_sources.get(_source_key(source))
        source_records = [
            window_source(record).record if record.annotations.get("gbdraw_region_applied") else record
            for record in (parsed if parsed is not None else load_render_source_records(source)).records
        ]
        drawn: dict[int, ResolvedRecordProvenance] = {}
        for other in provenance:
            if source_of(other) == key:
                drawn.setdefault(other.source_record_index, other)
        letters[key] = source_protein_letters(
            source_records,
            feature_visibility_rules=inputs.features.feature_visibility_rules,
            feature_overrides=_source_record_edits(edits, source_records, drawn),
        )
    rows, source_ids, source_lengths = _genome_inputs(collection, loaded.layout)
    return genome_database_sizes(
        rows=rows,
        source_ids=source_ids,
        source_lengths=source_lengths,
        searched_lengths=[len(record) for record in collection.records],
        source_proteins=[letters[source_of(item)] for item in provenance],
    )


def _source_record_edits(
    edits: Sequence[FeatureOverride],
    source_records: Sequence[SeqRecord],
    drawn: Mapping[int, ResolvedRecordProvenance],
) -> list[Mapping[int, ResolvedFeatureOverride]]:
    """The feature edits of each record of one source, bound to the uncropped records by feature identity.

    ``drawn`` is a drawn unit of each drawn record, by source record index (its key and source catalog name the
    record); an undrawn record has no edits.
    """

    indexes = sorted(drawn)
    keys = [drawn[index].record_identity_key for index in indexes]
    wanted = [row for row in edits if row.identity.record_key in set(keys)]
    if not wanted:
        return [{} for _ in source_records]
    catalogs = [drawn[index].source_feature_catalog for index in indexes]
    records = [source_records[index] for index in indexes]
    bindings = resolve_feature_identities(
        records=records, record_keys=keys, source_catalogs=catalogs, identities={row.identity for row in wanted},
    )
    bound = dict(zip(indexes, bind_feature_overrides(wanted, bindings, catalogs), strict=True))
    return [bound.get(index, {}) for index in range(len(source_records))]


def resolve_similarity_alignment_plan(
    plan: LinearRequestPlan,
    orthogroups: OrthogroupResult | OrthogroupGraphResult | None,
    reference: SimilarityAlignmentReference,
) -> SimilarityAlignmentPlan:
    """Resolve an exact feature/protein reference into one complete plan.

    ``orthogroups`` is the Similarity group result of the orthogroup analysis
    run on ``plan``'s records (``linear_metadata.orthogroups`` of the built
    plan). The reference names one member by its exact protein ID or feature
    SVG ID; a Similarity Group ID is rejected. Each other record aligns on its
    sole usable member or its sole direct reciprocal-best-hit member, and a
    record without a candidate keeps its position. A record that still has
    several candidates raises :class:`SimilarityAlignmentReferenceError`
    naming the record and its exact candidate IDs; nothing prompts.
    """

    if not isinstance(plan, LinearRequestPlan):
        raise ValidationError(
            "plan must be LinearRequestPlan.",
            diagnostic={"code": "INPUT_INVALID"},
        )
    if not isinstance(reference, SimilarityAlignmentReference):
        raise ValidationError(
            "reference must be SimilarityAlignmentReference.",
            diagnostic={"code": "INPUT_INVALID"},
        )
    if orthogroups is None:
        raise SimilarityAlignmentReferenceError(
            reference,
            "requires orthogroup metadata from the requested analysis.",
        )
    target = reference.feature_id
    if target in orthogroups.orthogroups:
        raise SimilarityAlignmentReferenceError(
            reference,
            "accepts an exact feature/protein ID, not "
            f"Similarity Group ID {target!r}.",
        )
    collection = ResolvedRecordCollection(plan.records, plan.provenance)
    record_keys = alignment_record_keys(collection)
    candidate_rows: list[tuple[OrthogroupMember, SimilarityAlignmentCandidate]] = []
    exact_matches: list[tuple[OrthogroupMember, SimilarityAlignmentCandidate]] = []
    for members in orthogroups.orthogroups.values():
        for member in members:
            if member.record_index < 0 or member.record_index >= len(collection.records):
                continue
            candidate = _similarity_candidate_from_orthogroup_member(
                collection,
                member,
            )
            candidate_rows.append((member, candidate))
            aliases = {
                str(value)
                for value in (
                    member.protein_id,
                    member.source_protein_id,
                    member.feature_svg_id,
                    candidate.anchor.biological_feature_id,
                    candidate.anchor.stable_feature_svg_id,
                )
                if value
            }
            if target in aliases:
                exact_matches.append((member, candidate))
    if not exact_matches:
        raise SimilarityAlignmentReferenceError(
            reference,
            "did not match an exact feature/protein ID.",
        )
    if len(exact_matches) != 1:
        matches = ", ".join(
            sorted(
                f"{candidate.anchor.record_key}:{member.protein_id}"
                for member, candidate in exact_matches
            )
        )
        raise SimilarityAlignmentReferenceError(
            reference,
            f"is ambiguous; use one exact candidate ID from: {matches}.",
        )
    reference_member, reference_candidate = exact_matches[0]
    group_id = reference_member.orthogroup_id
    resolution = resolve_similarity_alignment(
        record_keys=record_keys,
        group_id=group_id,
        reference=reference_candidate.anchor,
        candidates=tuple(
            candidate
            for member, candidate in candidate_rows
            if member.orthogroup_id == group_id
        ),
        edges=similarity_alignment_evidence_edges(collection, orthogroups, group_id),
    )
    if resolution.ambiguities:
        protein_by_key = {
            candidate.anchor.canonical_key: member.protein_id
            for member, candidate in candidate_rows
        }
        details = "; ".join(
            f"record {ambiguity.record_key!r} candidates ["
            + ", ".join(
                protein_by_key.get(
                    candidate.anchor.canonical_key,
                    candidate.anchor.biological_feature_id,
                )
                for candidate in ambiguity.candidates
            )
            + "]"
            for ambiguity in resolution.ambiguities
        )
        raise SimilarityAlignmentReferenceError(
            reference,
            "cannot choose among multiple candidates; select an exact candidate "
            f"for each record: {details}.",
        )
    return resolution.require_plan()


def _reference_record_inputs(
    request: LinearDiagramRequest,
    plan: LinearRequestPlan,
    reference: SimilarityAlignmentReference,
) -> tuple[RecordInput, ...]:
    """Key each input by its one displayed record, as the plan names them."""

    keys_by_input: dict[int, list[str]] = {}
    for item in plan.provenance:
        keys_by_input.setdefault(item.input_index, []).append(item.record_key)
    records: list[RecordInput] = []
    for index, record in enumerate(request.records):
        keys = keys_by_input.get(index, [])
        if len(keys) != 1:
            raise SimilarityAlignmentReferenceError(
                reference,
                "requires every RecordInput to resolve to exactly one displayed record.",
            )
        records.append(
            record if record.record_key == keys[0] else replace(record, record_key=keys[0])
        )
    return tuple(records)


def _resolved_losat_search(
    search: LosatSearchOptions | None,
) -> LosatSearchOptions | None:
    """Keep the search settings of evidence that the request now carries."""

    if search is None:
        return None
    return replace(search, losatp_mode="none", pairs=None)


def stored_protein_result_request(
    request: DiagramRequest,
    metadata: LinearDiagramMetadata | None,
) -> DiagramRequest:
    """The request a Session stores after a render that ran a protein pipeline.

    A Linear request that generated a LOSATP result stores it as the Web
    Session writer does after a Generate (OV-1068): the pipeline marker of
    mode ``"none"`` reuses the stored result; Pairwise and Similarity groups
    store each drawn pair's table, Similarity groups also its typed result,
    Collinear only its typed result. So a Session draws its comparison
    without LOSAT, the Web's included. Any other request is stored as given.
    ``metadata`` is the render's; nothing is computed again. The drawn
    results are stored in each record's search frame (R5, OV-1145).
    """

    if not isinstance(request, LinearDiagramRequest):
        return request
    options = request.options
    mode = losatp_analysis_mode(options.losat_search)
    if mode not in {"pairwise", "orthogroup", "collinear"}:
        return request
    if metadata is None or not {
        "pairwise": metadata.protein_comparisons is not None or bool(metadata.linear_comparisons),
        "orthogroup": metadata.orthogroups is not None,
        "collinear": metadata.collinearity_result is not None,
    }[mode]:
        raise ValueError(f"computed no {mode} result")
    metadata = _search_frame_metadata(metadata)
    stored = replace(options, losat_search=_resolved_losat_search(options.losat_search))
    return replace(
        request,
        options=(
            replace(stored, collinearity_blocks=metadata.collinearity_result)
            if mode == "collinear"
            else replace(
                stored,
                orthogroups=metadata.orthogroups,
                protein_comparisons=tuple(metadata.protein_comparisons or ()),
                linear_comparisons=_stored_linear_comparisons(options, metadata),
            )
        ),
    )


def saved_protein_result_request(
    request: DiagramRequest,
    artifacts: CurrentRequestArtifacts,
) -> DiagramRequest:
    """``request`` as a Session stores it after a run (:func:`stored_protein_result_request`), its protein
    pipeline answered by the raw LOSAT rows in ``artifacts`` that a Session saved.

    A request without a protein pipeline is returned as given. A search the rows
    do not answer runs where a LOSAT runtime exists; without one (the browser) the
    build raises its ``LOSAT_RUNTIME`` error.
    """

    if not isinstance(request, LinearDiagramRequest) or losatp_analysis_mode(
        request.options.losat_search
    ) not in {"pairwise", "orthogroup", "collinear"}:
        return request
    built = build_request_plan_diagram(plan_linear_request(request), artifacts=artifacts)
    return stored_protein_result_request(built.request, built.linear_metadata)


def _stored_linear_comparisons(
    options: LinearDiagramOptions,
    metadata: LinearDiagramMetadata,
) -> tuple[LinearComparison, ...]:
    """The explicit comparisons a stored request keeps: the given ones as
    given, then the pairs the analysis drew (a row layout's adjacent-row
    pairs, OV-2297), from ``metadata`` in the search frame."""

    given = tuple(options.linear_comparisons or ())
    return (*given, *tuple(metadata.linear_comparisons or ())[len(given):])


def _search_frame_metadata(metadata: LinearDiagramMetadata) -> LinearDiagramMetadata:
    """The drawn results of ``metadata`` in each record's search frame.

    Every Python writer stores them so (R5, OV-1145): the request of a
    Session, the resolved alignment request and the derived payloads.
    Metadata without bindings (no build attached them) is returned as given.
    """

    if not metadata.search_frame_bindings or not has_source_bound_results(
        metadata.linear_comparisons, metadata.protein_comparisons,
        metadata.orthogroups, metadata.collinearity_result,
    ):
        return metadata
    stored = project_source_bound_results(
        LinearDiagramOptions(
            linear_comparisons=metadata.linear_comparisons,
            protein_comparisons=metadata.protein_comparisons,
            orthogroups=metadata.orthogroups,
            collinearity_blocks=metadata.collinearity_result,
        ),
        metadata.search_frame_bindings,
    ).options
    collinearity = stored.collinearity_blocks
    assert collinearity is None or isinstance(collinearity, CollinearityResult)
    return replace(
        metadata,
        linear_comparisons=tuple(stored.linear_comparisons or ()),
        protein_comparisons=(
            None if stored.protein_comparisons is None else tuple(stored.protein_comparisons)
        ),
        orthogroups=stored.orthogroups,
        collinearity_result=collinearity,
    )


def _resolve_similarity_alignment_reference(
    request: DiagramRequest,
    artifacts: CurrentRequestArtifacts | None,
) -> tuple[DiagramRequest, CurrentRequestArtifacts | None]:
    """Run the requested analysis once and replace a reference with its plan.

    The resolved request carries that analysis as precomputed comparisons, so
    rendering it starts no second search, and the returned artifacts carry
    its LOSATP entries. Any other request passes through unchanged.
    """

    if not isinstance(request, LinearDiagramRequest) or not isinstance(
        request.similarity_alignment,
        SimilarityAlignmentReference,
    ):
        return request, artifacts
    reference = request.similarity_alignment
    analysis_plan = plan_linear_request(replace(request, similarity_alignment=None))
    analysis = build_request_plan_diagram(analysis_plan, artifacts=artifacts)
    drawn = analysis.linear_metadata
    alignment = resolve_similarity_alignment_plan(
        analysis_plan,
        drawn.orthogroups if drawn is not None else None,
        reference,
    )
    assert drawn is not None  # a missing result has no orthogroups
    metadata = _search_frame_metadata(drawn)
    records = _reference_record_inputs(request, analysis_plan, reference)
    layout = request.layout or LinearMultiRecordOptions()
    # The build drew the comparisons; the request keeps the planned ones as
    # given and stores the analysis rows and typed results, all in the
    # search frame (OV-399, OV-1145).
    planned = analysis_plan.request.options
    resolved = replace(
        request,
        records=records,
        options=replace(
            request.options,
            losat_search=_resolved_losat_search(request.options.losat_search),
            protein_comparisons=(
                planned.protein_comparisons
                if planned.protein_comparisons is not None or metadata.protein_comparisons is None
                else tuple(metadata.protein_comparisons)
            ),
            linear_comparisons=_stored_linear_comparisons(planned, metadata),
            orthogroups=metadata.orthogroups,
            collinearity_blocks=metadata.collinearity_result,
        ),
        layout=replace(
            layout,
            record_translations=layout.record_translations or tuple(
                LinearRecordTranslation(record_key=str(record.record_key))
                for record in records
            ),
        ),
        similarity_alignment=alignment,
    )
    return resolved, CurrentRequestArtifacts(
        losat_cache_entries=analysis.losat_cache_entries,
        losat_derived_cache_entries=analysis.losat_derived_cache_entries,
        protein_identity_manifest=analysis.protein_identity_manifest,
        protein_source_mode="orthogroup",
    )


def plan_request(request: DiagramRequest) -> DiagramRequestPlan:
    """Resolve one typed request into a reusable, not-yet-built render plan."""

    if isinstance(request, CircularBatchRequest):
        return plan_circular_batch_request(request)
    if isinstance(request, CircularDiagramRequest):
        return plan_circular_request(request)
    if isinstance(request, LinearDiagramRequest):
        return plan_linear_request(request)
    raise ValidationError("Unsupported diagram request type.")


def resolve_request(request: DiagramRequest) -> DiagramRequest:
    """Resolve an unresolved typed request to an exact-one projection.

    A Linear ``anchor_feature`` resolves to the record translations it was
    drawn with (:func:`resolved_plan_request`).
    """

    return resolved_plan_request(plan_request(request))


def resolved_plan_request(plan: DiagramRequestPlan) -> DiagramRequest:
    """The request ``plan`` resolves to, as :func:`resolve_request` and every Session writer store it.

    A Linear ``anchor_feature`` is resolved when the figure is laid out, so
    such a plan is built once (:func:`with_anchor_feature_translations`).
    """

    if isinstance(plan, LinearRequestPlan) and plan.request.options.anchor_feature is not None:
        return build_request_plan_diagram(plan).request
    return plan.request


def with_anchor_feature_translations(
    request: LinearDiagramRequest,
    metadata: LinearDiagramMetadata,
) -> LinearDiagramRequest:
    """``request`` with its ``anchor_feature`` replaced by the record translations its drawing used.

    The one owner of what an anchored figure stores: every record's
    translation as ``_final_record_translations`` applied it (the records the
    anchor moved, and every other record's own translation), without the
    option, so a Session replays the figure without aligning again and the
    Web's Align at feature writes the same values. ``metadata`` is the build
    result of ``request``'s plan. Returns ``request`` itself when it sets no
    anchor.
    """

    if request.options.anchor_feature is None:
        return request
    return replace(
        request,
        options=replace(request.options, anchor_feature=None),
        layout=replace(
            request.layout if request.layout is not None else LinearMultiRecordOptions(),
            record_translations=metadata.record_translations,
        ),
    )


def _source_protein_mode(
    request: LinearDiagramRequest,
    artifacts: CurrentRequestArtifacts,
) -> str:
    requested = losatp_analysis_mode(request.options.losat_search)
    if requested != "none":
        return requested
    if artifacts.protein_source_mode in {"pairwise", "orthogroup", "collinear"}:
        return str(artifacts.protein_source_mode)
    for entry in artifacts.losat_derived_cache_entries:
        mode = str(entry.get("mode") or "")
        if mode in {"pairwise", "orthogroup", "collinear"}:
            return mode
    return "orthogroup" if request.options.orthogroups is not None else "pairwise"


def _record_analysis_selector(item: ResolvedRecordProvenance) -> str | None:
    """The record selector a protein record analysis of ``item`` names: its input's, else ``#n`` for the n-th record
    of a file its ALL input expanded to (``<key>:<n>``), the selector of the Web card each such record becomes.

    The planned request holds no selector (its records are materialized), so the plan's provenance gives it.
    """

    if item.selector is not None:
        return item.selector.raw
    if item.cardinality is RecordCardinality.ALL and item.source_record_count > 1:
        return f"#{item.source_record_index + 1}"
    return None


def _extract_linear_request_proteins(
    request: LinearDiagramRequest,
    records: tuple[SeqRecord, ...],
    inputs: PreparedDiagramInputs,
    provenance: Sequence[ResolvedRecordProvenance],
) -> ProteinExtractionResult:
    """The proteins of the drawn units of a planned request, with the record analyses the Web's extraction of each
    unit writes: each unit of a window names the region its search reads (``ResolvedWindow.search_region``), each
    unit the record selector of its card (:func:`_record_analysis_selector`), and the manifest names every protein
    in its record's search frame (``extract_protein_identity_manifest``)."""

    record_instance_keys = _request_record_instance_keys(request)
    extraction = render_stage(
        "protein_extraction",
        extract_web_stable_cds_proteins,
        records,
        record_instance_keys=record_instance_keys,
        record_source_ids=tuple(str(record.id) for record in records),
        record_selectors=tuple(_record_analysis_selector(item) for item in provenance),
        regions=tuple(
            item.window.search_region if item.window is not None else None
            for item, _record in zip(provenance, records, strict=True)
        ),
        feature_visibility_rules=inputs.features.feature_visibility_rules,
        feature_overrides=tuple(item.overrides for item in inputs.record_features) or None,
    )
    if extraction.identity_manifest is None:
        raise ValidationError("Protein extraction did not produce an identity manifest.")
    return extraction


def _prepare_linear_artifacts(
    request: LinearDiagramRequest,
    records: tuple[SeqRecord, ...],
    artifacts: CurrentRequestArtifacts,
    inputs: PreparedDiagramInputs,
    provenance: Sequence[ResolvedRecordProvenance],
    *,
    extraction: ProteinExtractionResult | None,
    database_sizes: GenomeDatabaseSizes | None = None,
) -> _PreparedLinearArtifacts:
    current_protein = tuple(
        entry
        for entry in artifacts.losat_cache_entries
        if is_protein_losat_cache_entry(entry)
    )
    nucleotide_entries = tuple(
        copy.deepcopy(dict(entry))
        for entry in artifacts.losat_cache_entries
        if _is_current_nucleotide_losat_entry(entry)
    )
    needs_protein_identity = bool(
        current_protein
        or losatp_analysis_mode(request.options.losat_search) != "none"
        or request.options.protein_comparisons is not None
        or request.options.orthogroups is not None
        or request.options.collinearity_blocks is not None
    )
    passthrough_derived = tuple(
        copy.deepcopy(dict(entry))
        for entry in artifacts.losat_derived_cache_entries
    )
    if not needs_protein_identity:
        return _PreparedLinearArtifacts(
            cache=None,
            extraction=None,
            nucleotide_entries=nucleotide_entries,
            passthrough_derived_entries=passthrough_derived,
            source_mode="none",
        )

    if extraction is None:
        extraction = _extract_linear_request_proteins(request, records, inputs, provenance)
    manifest = extraction.identity_manifest
    assert manifest is not None
    reusable_current = tuple(
        entry
        for entry in current_protein
        if validate_protein_raw_entry_references(entry, manifest)
    )
    cache = LosatpCacheManager(
        reusable_current,
        identity_manifest=manifest,
        threads_per_job=(
            request.options.losat_search.runtime.threads
            if request.options.losat_search is not None
            else None
        )
        or "auto",
        database_sizes=database_sizes,
    )
    return _PreparedLinearArtifacts(
        cache=cache,
        extraction=extraction,
        nucleotide_entries=nucleotide_entries,
        passthrough_derived_entries=(),
        source_mode=_source_protein_mode(request, artifacts),
    )


class _ProteinPairFrame(NamedTuple):
    frame: DataFrame
    pair_index: int
    query_index: int
    subject_index: int


def _require_comparison_columns(frame: DataFrame) -> None:
    missing = [column for column in COMPARISON_COLUMNS if column not in frame.columns]
    if missing:
        raise ValidationError(
            "Protein comparison is missing canonical column(s): "
            + ", ".join(missing)
        )


def _resolved_protein_pair_frames(
    metadata: LinearDiagramMetadata | None,
    request: LinearDiagramRequest,
) -> list[_ProteinPairFrame]:
    """The comparison frames of the derived pair payloads, columns checked."""

    direct: Sequence[Any] | None = (
        metadata.protein_comparisons if metadata is not None else None
    )
    if direct is None:
        direct = request.options.protein_comparisons
    explicit: Sequence[LinearComparison] | None = (
        metadata.linear_comparisons if metadata is not None else None
    )

    result: list[_ProteinPairFrame] = []
    # resolve_linear_options refused a table past the last record pair.
    for pair_index, frame in enumerate(direct or ()):
        if not isinstance(frame, DataFrame):
            continue
        _require_comparison_columns(frame)
        result.append(_ProteinPairFrame(frame, pair_index, pair_index, pair_index + 1))
    if result:
        return result

    explicit = explicit if explicit is not None else request.options.linear_comparisons
    for comparison in explicit or ():
        frame = getattr(comparison, "matches", None)
        if not isinstance(frame, DataFrame) or not {
            "query_protein_id",
            "subject_protein_id",
        }.issubset(frame.columns):
            continue
        _require_comparison_columns(frame)
        result.append(
            _ProteinPairFrame(
                frame,
                len(result),
                int(comparison.query_record_index),
                int(comparison.subject_record_index),
            )
        )
    return result


def _derived_compact_resource(
    mode: str,
    collinearity: object,
    orthogroups: object,
) -> dict[str, bytes]:
    """The canonical bytes of the computed result a derived entry embeds."""

    from gbdraw.session_request_codec import encode_canonical_typed_resource

    # The canonical bytes are the stage result, so no reader can change them.
    if mode == "collinear" and isinstance(collinearity, CollinearityResult):
        return {
            "collinearityResult": render_stage(
                "derived_losat_entries", encode_canonical_typed_resource, "result", collinearity
            )
        }
    if orthogroups is not None:
        return {
            "orthogroupResult": render_stage(
                "derived_losat_entries", encode_canonical_typed_resource, "orthogroupResult", orthogroups
            )
        }
    return {}


def _request_record_instance_keys(request: DiagramRequest) -> tuple[str, ...]:
    """The record instance key of each drawn record (the Web's record keys)."""

    return tuple(
        record_input.record_key or f"record-{index + 1}"
        for index, record_input in enumerate(request.records)
    )


def _request_derived_settings(
    request: LinearDiagramRequest,
    mode: str,
    rows_by_record: Sequence[int],
) -> DerivedLosatpSettings:
    """The converter settings a request's derived payload was computed with. ``rows_by_record`` are the rows the
    drawing places its records in: a drawing whose records stand in a grid displays its grid pairs
    (``grid_display_pairs``), as the Web's unit plan says."""

    options = request.options
    params = options.collinearity_params or LosslessCollinearityParameters()
    search = options.losat_search
    explicit_display_pairs = grid_display_pairs(rows_by_record) is not None
    return normalize_derived_losatp_settings(
        mode=mode,
        bitscore=options.bitscore,
        evalue=options.evalue,
        identity=options.identity,
        alignment_length=options.alignment_length,
        membership_mode=str(options.orthogroup_membership_mode),
        member_max_hits=search.losatp_member_max_hits if search is not None else None,
        min_anchors=params.min_anchors,
        max_unit_gap=params.max_unit_gap,
        max_diagonal_drift=params.max_diagonal_drift,
        max_conflicts=params.max_conflicts,
        merge_orientation=str(params.merge_orientation),
        unit_mode=str(options.collinearity_unit_mode),
        color_mode=str(options.collinearity_color_mode),
        anchor_mode=str(options.collinearity_anchor_mode),
        search_scope=str(options.collinearity_search_scope),
        max_paralog_links=options.collinear_max_paralog_links_per_orthogroup,
        infer_orthogroups=options.collinear_infer_orthogroups,
        explicit_display_pairs=explicit_display_pairs,
    )


def _request_losatp_job_specs(
    request: LinearDiagramRequest,
    mode: str,
    pairs: Sequence[tuple[int, int]] | None,
) -> tuple[tuple[int, int], ...]:
    """The directed record pairs that a LOSATP run of ``request`` in ``mode`` searches.

    ``pairs`` are its displayed comparison edges (``None``: the consecutive
    records); the registered job plan (``losatp_job_specs``) gives the searches.
    """

    options = request.options
    return losatp_job_specs(
        mode,
        record_count=len(request.records),
        pairs=pairs,
        infer_orthogroups=options.collinear_infer_orthogroups,
        search_scope=str(options.collinearity_search_scope),
    )


def _derived_identity_pairs(
    request: LinearDiagramRequest,
    settings: DerivedLosatpSettings,
    pair_frames: Sequence[_ProteinPairFrame],
    raw_entries: Sequence[Mapping[str, Any]],
) -> tuple[DerivedLosatpPair, ...] | None:
    """The searches of the run in the Web's search order, or None when one has
    no single raw entry (the payload was not computed from raw LOSAT results)."""

    instance_keys = _request_record_instance_keys(request)
    displayed = {
        (pair.query_index, pair.subject_index): pair.pair_index for pair in pair_frames
    }
    raw_keys: dict[tuple[str, str], list[str]] = {}
    for entry in raw_entries:
        if is_protein_losat_cache_entry(entry):
            raw_keys.setdefault(
                (
                    str(entry.get("queryRecordInstanceKey") or ""),
                    str(entry.get("subjectRecordInstanceKey") or ""),
                ),
                [],
            ).append(str(entry.get("key") or ""))
    pairs: list[DerivedLosatpPair] = []
    for query_index, subject_index in _request_losatp_job_specs(
        request, settings.mode, tuple(displayed)
    ):
        keys = raw_keys.get((instance_keys[query_index], instance_keys[subject_index]), [])
        if len(keys) != 1 or not keys[0]:
            return None
        number = displayed.get((query_index, subject_index))
        pairs.append(
            DerivedLosatpPair(
                query_index=query_index,
                subject_index=subject_index,
                raw_key=keys[0],
                display_pair=number is not None,
                pair_index=number,
            )
        )
    return tuple(pairs)


def _build_current_derived_entries(
    metadata: LinearDiagramMetadata | None,
    request: LinearDiagramRequest,
    records: tuple[SeqRecord, ...],
    layout: LinearMultiRecordOptions | None,
    artifacts: _PreparedLinearArtifacts,
    raw_entries: Sequence[Mapping[str, Any]],
) -> tuple[Mapping[str, Any], ...]:
    wanted = _DERIVED_LOSAT_ENTRIES_WANTED.get()
    if artifacts.extraction is None:
        return artifacts.passthrough_derived_entries if wanted else ()
    mode = artifacts.source_mode
    if mode not in {"orthogroup", "collinear"}:
        return ()
    manifest = artifacts.extraction.identity_manifest
    if manifest is None:
        raise ValidationError("Derived protein payload requires an identity manifest.")
    # A derived payload holds the stored form, as a Web payload does (OV-1145).
    metadata = None if metadata is None else _search_frame_metadata(metadata)
    pair_frames = _resolved_protein_pair_frames(metadata, request)
    orthogroups = (
        metadata.orthogroups
        if metadata is not None
        else request.options.orthogroups
    )
    collinearity = metadata.collinearity_result if metadata is not None else request.options.collinearity_blocks
    if not pair_frames and not (orthogroups is not None and getattr(orthogroups, "orthogroups", None)):
        return ()
    # The Web switch keeps the encoding check of the computed result.
    compact_resource = _derived_compact_resource(mode, collinearity, orthogroups)
    # As on the Web, a derived payload holds its computed result and comes
    # from raw LOSAT results: one owner gives its identity, key and pairs.
    if not wanted or not compact_resource:
        return ()
    settings = _request_derived_settings(
        request,
        mode,
        resolve_record_row_positions(records, layout.multi_record_positions if layout is not None else None)[1],
    )
    pairs = _derived_identity_pairs(request, settings, pair_frames, raw_entries)
    if pairs is None:
        return ()
    record_instances = manifest.record_instances
    derived_records = []
    for index, instance_key in enumerate(_request_record_instance_keys(request)):
        binding = record_instances.get(instance_key, {})
        # The drawn unit's orientation, as the Web's view transform of it states (the planned request's records
        # are materialized and hold none).
        derived_records.append(
            DerivedLosatpRecord(
                index=index,
                runtime_binding_hash=str(binding.get("runtimeBindingHash") or ""),
                display_binding_hash=str(binding.get("displayBindingHash") or ""),
                length=len(records[index].seq),
                reverse=record_is_reverse_complemented(records[index]),
            )
        )
    identity = derived_losatp_identity(settings, derived_records, pairs)
    entry: Mapping[str, Any] = {
        "schema": CURRENT_DERIVED_PROTEIN_ARTIFACT_SCHEMA,
        "kind": "derived-losatp-payload",
        "idEncoding": "runtime-handle-v1",
        "key": derived_losatp_cache_key(identity),
        "mode": mode,
        "payload": {
            "identity": identity,
            "provenance": derived_losatp_provenance(
                identity,
                settings,
                effective_unit_kinds=effective_collinearity_unit_kinds(
                    collinearity if isinstance(collinearity, CollinearityResult) else None
                ),
            ),
            "pairs": [
                derived_pair_payload(
                    pair.pair_index, pair.query_index, pair.subject_index, pair.frame
                )
                for pair in pair_frames
            ],
            **{
                name: json.loads(encoded)
                for name, encoded in compact_resource.items()
            },
        },
    }
    return (entry,)


def _merged_losat_entries(
    loaded: Sequence[Mapping[str, Any]],
    searched: Sequence[Mapping[str, Any]],
) -> tuple[Mapping[str, Any], ...]:
    """Session entries first, then entries this run searched (one per key)."""

    keys = {entry.get("key") for entry in loaded}
    return (*loaded, *(entry for entry in searched if entry.get("key") not in keys))


@overload
def build_request_plan_diagram(
    plan: CircularBatchRequestPlan,
    *,
    artifacts: CurrentRequestArtifacts | None = None,
    label_fit_warn_search_seconds: float = WARN_SEARCH_SECONDS,
) -> PreparedCircularBatchRequest: ...


@overload
def build_request_plan_diagram(
    plan: CircularRequestPlan | LinearRequestPlan,
    *,
    artifacts: CurrentRequestArtifacts | None = None,
    label_fit_warn_search_seconds: float = WARN_SEARCH_SECONDS,
) -> PreparedDiagramRequest: ...


@overload
def build_request_plan_diagram(
    plan: DiagramRequestPlan,
    *,
    artifacts: CurrentRequestArtifacts | None = None,
    label_fit_warn_search_seconds: float = WARN_SEARCH_SECONDS,
) -> PreparedDiagramRequest | PreparedCircularBatchRequest: ...


def build_request_plan_diagram(
    plan: DiagramRequestPlan,
    *,
    artifacts: CurrentRequestArtifacts | None = None,
    label_fit_warn_search_seconds: float = WARN_SEARCH_SECONDS,
) -> PreparedDiagramRequest | PreparedCircularBatchRequest:
    """Build a previously resolved plan without loading or planning it again.

    ``label_fit_warn_search_seconds`` bounds the label fit search of a
    Circular figure drawn with crossing leaders (0: counts only, every row
    ``not_computed``).
    """

    return _build_request_plan_diagram(
        plan,
        artifacts=artifacts,
        label_fit_warn_search_seconds=label_fit_warn_search_seconds,
    )


def _build_request_plan_diagram(
    plan: DiagramRequestPlan,
    *,
    artifacts: CurrentRequestArtifacts | None,
    protein_extraction: ProteinExtractionResult | None = None,
    label_fit_warn_search_seconds: float = WARN_SEARCH_SECONDS,
) -> PreparedDiagramRequest | PreparedCircularBatchRequest:
    """Build ``plan``; ``protein_extraction`` is one its caller already ran on it."""

    if not isinstance(
        plan,
        (CircularRequestPlan, CircularBatchRequestPlan, LinearRequestPlan),
    ):
        raise ValidationError("plan must be a DiagramRequestPlan.")
    current_artifacts = (
        CurrentRequestArtifacts()
        if artifacts is None
        else artifacts
    )
    if not isinstance(current_artifacts, CurrentRequestArtifacts):
        raise ValidationError("artifacts must be CurrentRequestArtifacts.")
    if isinstance(plan, CircularBatchRequestPlan):
        with _request_render_diagnostic_phase("drawing"):
            items = tuple(
                PreparedDiagramRequest(
                    mode="circular",
                    request=item_plan.request,
                    records=item_plan.records,
                    drawing=item_plan.build(
                        label_fit_warn_search_seconds=label_fit_warn_search_seconds
                    ),
                    inputs=item_plan.inputs,
                    transforms=item_plan.transforms,
                    # Filled by CircularRequestPlan.__post_init__.
                    resolved_annotations=cast(
                        ResolvedAnnotationBundle, item_plan.resolved_annotations
                    ),
                    input_summary=_plan_input_summary(item_plan),
                )
                for item_plan in plan.item_plans()
            )
        return PreparedCircularBatchRequest(
            request=plan.request,
            records=plan.records,
            items=items,
            inputs=plan.inputs,
            losat_cache_entries=_merged_losat_entries(
                current_artifacts.losat_cache_entries, plan.losat_cache_entries
            ),
            protein_identity_manifest=current_artifacts.protein_identity_manifest,
        )
    request = plan.request
    records = plan.records
    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()
    losat_derived_cache_entries: tuple[Mapping[str, Any], ...] = ()
    protein_identity_manifest: Mapping[str, Any] | None = None
    linear_metadata: LinearDiagramMetadata | None = None
    if isinstance(plan, CircularRequestPlan):
        with _request_render_diagnostic_phase("drawing"):
            drawing = plan.build(label_fit_warn_search_seconds=label_fit_warn_search_seconds)
        losat_cache_entries = _merged_losat_entries(
            current_artifacts.losat_cache_entries, plan.losat_cache_entries
        )
        losat_derived_cache_entries = current_artifacts.losat_derived_cache_entries
        protein_identity_manifest = current_artifacts.protein_identity_manifest
    else:
        if plan.inputs is None:  # pragma: no cover - planners always prepare inputs.
            raise ValidationError("Linear request plan has no prepared diagram inputs.")
        with _request_render_diagnostic_phase("comparisonPreparation"):
            linear_artifacts = _prepare_linear_artifacts(
                plan.request,
                records,
                current_artifacts,
                plan.inputs,
                plan.provenance,
                extraction=protein_extraction,
                database_sizes=plan.losatp_database_sizes,
            )
        with _request_render_diagnostic_phase("drawing"):
            linear_build = plan.build(
                losatp_cache=linear_artifacts.cache,
                protein_extraction=linear_artifacts.extraction,
            )
        drawing = linear_build.drawing
        linear_metadata = linear_build.metadata
        request = with_anchor_feature_translations(plan.request, linear_metadata)
        with _request_render_diagnostic_phase("comparisonPreparation"):
            protein_entries = (
                tuple(linear_artifacts.cache.session_entries())
                if linear_artifacts.cache is not None
                else ()
            )
            nucleotide_keys = {
                entry.get("key") for entry in linear_artifacts.nucleotide_entries
            }
            losat_cache_entries = (
                *protein_entries,
                *linear_artifacts.nucleotide_entries,
                *(
                    entry
                    for entry in plan.losat_cache_entries
                    if entry.get("key") not in nucleotide_keys
                ),
            )
            losat_derived_cache_entries = _build_current_derived_entries(
                linear_metadata,
                plan.request,
                records,
                plan.layout,
                linear_artifacts,
                losat_cache_entries,
            )
            protein_identity_manifest = (
                linear_artifacts.extraction.identity_manifest
                if linear_artifacts.extraction is not None
                and linear_artifacts.extraction.identity_manifest is not None
                else current_artifacts.protein_identity_manifest
            )
    return PreparedDiagramRequest(
        mode=plan.mode,
        transforms=plan.transforms,
        # Filled by the plan's __post_init__.
        resolved_annotations=cast(ResolvedAnnotationBundle, plan.resolved_annotations),
        request=request,
        records=records,
        drawing=drawing,
        inputs=plan.inputs,
        linear_metadata=linear_metadata,
        losat_cache_entries=losat_cache_entries,
        losat_derived_cache_entries=losat_derived_cache_entries,
        protein_identity_manifest=protein_identity_manifest,
        input_summary=_plan_input_summary(
            plan, linear_metadata.comparison_tallies if linear_metadata is not None else (),
        ),
    )


def _plan_input_summary(
    plan: CircularRequestPlan | LinearRequestPlan,
    comparisons: Sequence[ComparisonTableTally] = (),
) -> InputSummary:
    """The input summary of a built plan, from the feature inputs and the track stack its drawing read."""

    slots = _plan_track_slots(plan)
    left_out = annotation_sets_left_out(cast(ResolvedAnnotationBundle, plan.resolved_annotations), slots)
    if not plan.provenance:  # A plan built without record planning has no sources to name.
        return InputSummary((), tuple(comparisons), left_out)
    inputs = plan.inputs
    return summarize_inputs(
        plan.records,
        plan.provenance,
        selected_features=plan.request.options.selected_features_set or DEFAULT_SELECTED_FEATURES,
        feature_visibility_rules=inputs.features.feature_visibility_rules if inputs is not None else None,
        specific_color_rules=inputs.features.specific_color_rules if inputs is not None else None,
        record_features=inputs.record_features if inputs is not None else (),
        comparisons=comparisons,
        draws_features=slots_draw_features(slots),
        annotation_sets_left_out=left_out,
        left_off_circle=_plan_records_left_off_circle(plan),
    )


def _plan_records_left_off_circle(plan: CircularRequestPlan | LinearRequestPlan) -> dict[int, int]:
    """The records one circle leaves out (by index), each with the shortest record it draws."""

    if not isinstance(plan, CircularRequestPlan) or plan.layout is None or plan.layout.record_layout != "one_circle":
        return {}
    shortest = plan.layout.one_circle_min_length
    return {
        unit.record_index: shortest
        for unit in record_units(plan.records, plan.transforms or None)
        if left_off_circle(unit, shortest)
    }


def _plan_track_slots(
    plan: CircularRequestPlan | LinearRequestPlan,
) -> list[CircularTrackSlot] | list[LinearTrackSlot] | None:
    """The plan's explicit track stack, parsed as the drawing parses it; ``None`` without one."""

    if isinstance(plan, CircularRequestPlan):
        circular = plan.request.options.tracks
        if circular is not None and circular.circular_track_slots is not None:
            return parse_circular_track_slots(list(circular.circular_track_slots))
        return None
    linear = plan.request.options.tracks
    if linear is not None and linear.linear_track_slots is not None:
        return parse_linear_track_slots(list(linear.linear_track_slots))
    return None


def build_request_diagram(
    request: DiagramRequest,
    *,
    artifacts: CurrentRequestArtifacts | None = None,
) -> PreparedDiagramRequest | PreparedCircularBatchRequest:
    """Normalize inputs and build a drawing from current typed artifacts."""

    request, artifacts = _resolve_similarity_alignment_reference(request, artifacts)
    return build_request_plan_diagram(
        plan_request(request),
        artifacts=artifacts,
    )


def _visibility_table(
    options: CircularDiagramOptions | LinearDiagramOptions,
) -> DataFrame | None:
    table = options.feature_visibility_table
    file_path = options.feature_visibility_table_file
    if table is None and file_path is not None:
        return read_feature_visibility_file(file_path)
    return table


def _color_table(
    options: CircularDiagramOptions | LinearDiagramOptions,
) -> DataFrame | None:
    colors = options.colors
    if colors is None or colors.color_table is not None:
        return colors.color_table if colors is not None else None
    if colors.color_table_file is not None:
        return read_color_table(colors.color_table_file)
    return None


def build_prepared_interactive_context(
    prepared: PreparedDiagramRequest,
    *,
    comparison_sequence_records: Sequence[Sequence[SeqRecord]] = (),
) -> InteractiveSvgContext:
    """Build interactive metadata from the same inputs used for drawing."""

    def build() -> InteractiveSvgContext:
        options = prepared.request.options
        inputs = prepared.inputs or _prepare_diagram_inputs(prepared.request)
        computed_orthogroups = _prepared_orthogroups(prepared)
        collinearity_search_scope = None
        if (
            isinstance(options, LinearDiagramOptions)
            and prepared.linear_metadata is not None
            and prepared.linear_metadata.collinearity_result is not None
            and prepared.linear_metadata.collinearity_result.orthogroups
            is not None
        ):
            collinearity_search_scope = str(
                options.collinearity_search_scope
            )
        # Keyed by the annotation rows the context reads, by value: a Region
        # Annotation style change keeps the context and the catalog.
        return render_stage(
            "interactive_context",
            interactive_svg_context,
            tuple(prepared.records),
            selected_features_set=options.selected_features_set,
            feature_visibility_rules=inputs.features.feature_visibility_rules,
            record_features=inputs.record_features,
            specific_color_rules=inputs.features.specific_color_rules,
            orthogroups=computed_orthogroups,
            linear_rendered_feature_ids=prepared.mode == "linear",
            annotation_rows=interactive_annotation_rows(prepared.resolved_annotations),
            mode=prepared.mode,
            comparison_sequence_records=tuple(tuple(items) for items in comparison_sequence_records),
            collinearity_search_scope=collinearity_search_scope,
            record_transforms=prepared.transforms or None,
        )

    return require_interactive_svg_metadata(build)


def _prepared_orthogroups(prepared: PreparedDiagramRequest) -> object | None:
    computed = (
        prepared.linear_metadata.orthogroups
        if prepared.linear_metadata is not None
        else None
    )
    if computed is None and isinstance(prepared.request.options, LinearDiagramOptions):
        computed = prepared.request.options.orthogroups
    return computed


def _interactive_context(
    prepared: PreparedDiagramRequest,
    *,
    comparison_sequence_records: Sequence[Sequence[SeqRecord]] = (),
    include_feature_catalog: bool = False,
) -> InteractiveSvgContext | None:
    output = prepared.request.output
    needs_interactive_metadata = (
        "interactive_svg" in output.formats
        and output.interactive_metadata_policy != "omit"
    )
    if not needs_interactive_metadata and not include_feature_catalog:
        return None
    # The render stage ``interactive_context`` reuses the context by its arguments.
    return build_prepared_interactive_context(
        prepared,
        comparison_sequence_records=comparison_sequence_records,
    )


def render_request(
    request: DiagramRequest,
    *,
    artifacts: CurrentRequestArtifacts | None = None,
    include_feature_catalog: bool = False,
    label_fit_warn_search_seconds: float = WARN_SEARCH_SECONDS,
) -> RequestRenderResult | CircularBatchRenderResult:
    """Build and save one typed request from current typed artifacts.

    ``label_fit_warn_search_seconds`` bounds the label fit search of a
    Circular figure drawn with crossing leaders; the Web passes 0 and computes
    rows on demand (OD-F14).
    """

    if isinstance(request, LinearDiagramRequest) and isinstance(
        request.similarity_alignment,
        SimilarityAlignmentReference,
    ):
        _preflight_render_output(request.output)  # before the analysis runs
    request, artifacts = _resolve_similarity_alignment_reference(request, artifacts)
    plan = plan_request(request)
    batch_outputs_preflighted = isinstance(plan, CircularBatchRequestPlan)
    plan.preflight_outputs()
    prepared = build_request_plan_diagram(
        plan,
        artifacts=artifacts,
        label_fit_warn_search_seconds=label_fit_warn_search_seconds,
    )
    return render_prepared_request(
        prepared,
        batch_outputs_preflighted=batch_outputs_preflighted,
        include_feature_catalog=include_feature_catalog,
    )


def render_prepared_request(
    prepared: PreparedDiagramRequest | PreparedCircularBatchRequest,
    *,
    batch_outputs_preflighted: bool = False,
    include_feature_catalog: bool = False,
) -> RequestRenderResult | CircularBatchRenderResult:
    """Save an already-built request without planning or loading it again."""

    if not isinstance(
        prepared,
        (PreparedDiagramRequest, PreparedCircularBatchRequest),
    ):
        raise ValidationError("prepared must be a prepared diagram request.")
    if isinstance(prepared, PreparedCircularBatchRequest):
        if not batch_outputs_preflighted:
            preflight_diagram_request_outputs((prepared.request,))
    return _render_request_diagram(
        prepared,
        include_feature_catalog=include_feature_catalog,
    )


def _render_request_diagram(
    prepared: PreparedDiagramRequest | PreparedCircularBatchRequest,
    *,
    include_feature_catalog: bool = False,
) -> RequestRenderResult | CircularBatchRenderResult:
    with _request_render_diagnostic_phase("interactivePreparation"):
        comparison_sequence_records = _comparison_sequence_records(
            prepared,
            include_feature_catalog=include_feature_catalog,
        )
    sources = prepared.inputs.comparison_sequences if prepared.inputs is not None else None
    repeated_ids = sources.repeated_record_ids() if sources is not None else ()
    if isinstance(prepared, PreparedCircularBatchRequest):
        return CircularBatchRenderResult(
            request=prepared.request,
            records=prepared.records,
            items=tuple(
                _render_prepared_request(
                    item,
                    comparison_sequence_records=comparison_sequence_records,
                    include_feature_catalog=include_feature_catalog,
                    repeated_comparison_ids=repeated_ids,
                )
                for item in prepared.items
            ),
            losat_cache_entries=prepared.losat_cache_entries,
            protein_identity_manifest=prepared.protein_identity_manifest,
        )
    return _render_prepared_request(
        prepared,
        comparison_sequence_records=comparison_sequence_records,
        include_feature_catalog=include_feature_catalog,
        repeated_comparison_ids=repeated_ids,
    )


def _comparison_sequence_records(
    prepared: PreparedDiagramRequest | PreparedCircularBatchRequest,
    *,
    include_feature_catalog: bool = False,
) -> tuple[tuple[SeqRecord, ...], ...]:
    """Load optional Circular comparison FASTA sources once per render."""

    options = prepared.request.options
    if not isinstance(options, CircularDiagramOptions):
        return ()
    outputs = (
        prepared.request.outputs
        if isinstance(prepared, PreparedCircularBatchRequest)
        else (prepared.request.output,)
    )
    if not include_feature_catalog and not any(
        "interactive_svg" in output.formats
        and output.interactive_metadata_policy != "omit"
        for output in outputs
    ):
        return ()

    return require_interactive_svg_metadata(
        lambda: load_comparison_sequence_records(prepared.request, prepared.inputs)
    )


def load_comparison_sequence_records(
    request: DiagramRequest,
    inputs: PreparedDiagramInputs | None = None,
) -> tuple[tuple[SeqRecord, ...], ...]:
    """The records of each comparison FASTA a Circular request names (none for a source without one), by source.

    The one loader: an interactive render and the sequence reads of a shown Result (D5) call it; ``inputs`` are the
    request's prepared inputs when the caller holds them.
    """

    options = request.options
    if not isinstance(options, CircularDiagramOptions):
        return ()
    if inputs is not None and inputs.comparison_sequences is not None:
        return inputs.comparison_sequences.load()
    return _ComparisonSequenceSources(options.conservation_genome_files).load()


def _render_output_paths(output: RenderOutputRequest) -> tuple[Path, ...]:
    base = Path(output.output_directory or ".") / output.output_prefix
    return tuple(
        Path(path)
        for path in resolve_output_paths(
            str(base),
            output.formats,
            include_base_svg=True,
        )
    )


def _preflight_render_output(output: RenderOutputRequest) -> None:
    preflight_output_paths(
        _render_output_paths(output),
        overwrite=output.overwrite,
    )


def _diagram_request_outputs(
    request: DiagramRequest,
) -> tuple[RenderOutputRequest, ...]:
    return (
        tuple(request.outputs)
        if isinstance(request, CircularBatchRequest)
        else (request.output,)
    )


def diagram_request_output_paths(request: DiagramRequest) -> tuple[Path, ...]:
    """Return every file path one resolved typed request will write."""

    return tuple(
        path
        for output in _diagram_request_outputs(request)
        for path in _render_output_paths(output)
    )


def preflight_diagram_request_outputs(requests: Sequence[DiagramRequest]) -> None:
    """Check every output of the requests together before any is written.

    Each output keeps its own overwrite permission; two outputs that resolve
    to one file are rejected.
    """

    path_groups = tuple(
        (output, _render_output_paths(output))
        for request in requests
        for output in _diagram_request_outputs(request)
    )
    for output, paths in path_groups:
        preflight_output_paths(paths, overwrite=output.overwrite)
    try:
        path_identities = [
            path.resolve(strict=False)
            for _output, paths in path_groups
            for path in paths
        ]
    except (OSError, ValueError) as exc:
        raise ValidationError(
            "Could not resolve one or more diagram output paths."
        ) from exc
    repeated = sorted(
        {str(path) for path in path_identities if path_identities.count(path) > 1}
    )
    if repeated:
        raise ValidationError(
            "Diagram output requests resolve to duplicate file paths: "
            + ", ".join(repeated)
            + "."
        )


def _render_prepared_request(
    prepared: PreparedDiagramRequest,
    *,
    comparison_sequence_records: Sequence[Sequence[SeqRecord]] = (),
    include_feature_catalog: bool = False,
    repeated_comparison_ids: tuple[RepeatedRecordIds, ...] = (),
) -> RequestRenderResult:
    """Write one already-planned diagram using its canonical output request.

    Its input summary adds the repeated record IDs of the comparison genomes the
    render read (OV-2954), known only once they are read.
    """

    output = prepared.request.output
    with _request_render_diagnostic_phase("interactivePreparation"):
        interactive_context = _interactive_context(
            prepared,
            comparison_sequence_records=comparison_sequence_records,
            include_feature_catalog=include_feature_catalog,
        )
    export_interactive_context = (
        interactive_context
        if (
            "interactive_svg" in output.formats
            and output.interactive_metadata_policy != "omit"
        )
        else None
    )
    with _request_render_diagnostic_phase("svgWrite"):
        paths = save_figure_to(
            prepared.drawing,
            output.formats,
            output_dir=(
                str(output.output_directory)
                if output.output_directory is not None
                else None
            ),
            output_prefix=output.output_prefix,
            overwrite=output.overwrite,
            interactive_context=export_interactive_context,
            embed_fonts=output.embed_fonts,
            dpi=output.dpi,
        )
    return RequestRenderResult(
        mode=prepared.mode,
        request=prepared.request,
        records=prepared.records,
        drawing=prepared.drawing,
        output_paths=tuple(Path(path) for path in paths),
        annotation_warnings=prepared.annotation_warnings,
        feature_identity_notices=prepared.feature_identity_notices,
        input_summary=replace(prepared.input_summary, repeated_comparison_ids=repeated_comparison_ids),
        interactive_context=interactive_context,
        linear_metadata=prepared.linear_metadata,
        losat_cache_entries=prepared.losat_cache_entries,
        losat_derived_cache_entries=prepared.losat_derived_cache_entries,
        protein_identity_manifest=prepared.protein_identity_manifest,
    )


__all__ = [
    "CircularBatchRenderResult",
    "CircularBatchRequestPlan",
    "CircularRequestPlan",
    "CurrentRequestArtifacts",
    "DiagramRequestPlan",
    "LinearRequestPlan",
    "PreparedCircularBatchRequest",
    "PreparedDiagramInputs",
    "PreparedDiagramRequest",
    "RequestRenderResult",
    "build_request_plan_diagram",
    "build_request_diagram",
    "build_prepared_interactive_context",
    "diagram_request_output_paths",
    "grid_placement_order",
    "normalize_request_records",
    "plan_request",
    "plan_circular_batch_request",
    "plan_circular_request",
    "GENOME_SCALE_PROGRAMS",
    "GenomeScaleProgram",
    "linear_genome_scales",
    "plan_linear_request",
    "preflight_diagram_request_outputs",
    "render_request",
    "render_prepared_request",
    "resolve_request",
    "resolved_plan_request",
    "resolve_similarity_alignment_plan",
    "with_anchor_feature_translations",
]
