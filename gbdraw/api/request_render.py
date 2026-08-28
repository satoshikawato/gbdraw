"""Normalize and render typed diagram requests without CLI/session orchestration."""

from __future__ import annotations

import copy
from contextlib import contextmanager
from contextvars import ContextVar
import hashlib
import json
import logging
from dataclasses import dataclass, replace
from pathlib import Path
from time import perf_counter
from typing import (
    Any,
    Hashable,
    Iterator,
    Literal,
    Mapping,
    MutableMapping,
    Sequence,
    TypeAlias,
)

from Bio import SeqIO  # type: ignore[reportMissingImports]
from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]

from gbdraw.analysis.protein_artifacts import (
    is_current_derived_protein_artifact,
    validate_current_derived_protein_artifacts,
)
from gbdraw.analysis.collinearity import (
    CollinearityResult,
    LosslessCollinearityParameters,
)
from gbdraw.analysis.depth_tracks import (
    DepthTrackSpec,
    depth_track_count,
    normalize_depth_tracks,
)
from gbdraw.exceptions import ValidationError
from gbdraw.analysis.protein_colinearity import (
    LosatpCacheManager,
    ProteinExtractionResult,
    extract_web_stable_cds_proteins,
    is_protein_losat_cache_entry,
    validate_protein_identity_manifest,
    validate_protein_raw_entry_references,
)
from gbdraw.features.visibility import (
    read_feature_visibility_file,
    resolve_candidate_feature_types,
)
from gbdraw.annotations import AnnotationOptions, read_annotation_table
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.io.colors import load_default_colors, read_color_table
from gbdraw.labels.filtering import (
    read_filter_list_file,
    read_label_override_file,
    read_qualifier_priority_file,
)
from gbdraw.render.interactive_context import (
    build_interactive_svg_context,
    require_interactive_svg_metadata,
)
from gbdraw.render.interactive_svg import InteractiveSvgContext
from gbdraw.render.formats import resolve_output_paths
from gbdraw.web_support.orthogroup_metadata import serialize_orthogroups_payload
from gbdraw.version import __version__

from .diagram import (
    DEFAULT_SELECTED_FEATURES,
    LinearDiagramBuildResult,
    LinearDiagramMetadata,
    build_circular_diagram,
    build_circular_multi_diagram,
    build_linear_diagram_result,
)
from .io import load_gbks, load_gff_fasta
from .options import (
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
)
from .prepared import (
    PreparedResourceIdentity,
    ResolvedFeatureInputs,
    get_or_build_interactive_context,
    get_or_build_parsed_source,
    get_or_build_resolved_records,
    prepared_input_cache_active,
    prepared_record_membership,
    prepared_resource_identity,
    prepared_resource_value_identity,
    register_prepared_resource_value,
    resolve_feature_inputs,
)
from .record_planning import (
    ResolvedRecordCollection,
    ResolvedRecordProvenance,
    _load_source_records,
    _prepared_source_cache_spec,
    resolve_circular_batch_outputs,
    resolve_circular_options,
    resolve_implicit_record_output_prefix,
    resolve_linear_options,
    resolve_record_inputs,
)
from .render import preflight_output_paths, save_figure_to
from .requests import (
    CircularBatchRequest,
    CircularDiagramRequest,
    DiagramRequest,
    GenBankInputSource,
    GffFastaInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordCardinality,
    RecordCollectionOptions,
    RecordInput,
    RecordPresentation,
    RenderOutputRequest,
)

logger = logging.getLogger(__name__)


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


def _resolve_request_option_tables(
    options: CircularDiagramOptions | LinearDiagramOptions,
    *,
    mode: Literal["circular", "linear"],
    load_comparison_colors: bool,
) -> CircularDiagramOptions | LinearDiagramOptions:
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
        default_colors = colors.default_colors
        if default_colors is None and colors.default_colors_file is not None:
            default_colors = register_prepared_resource_value(
                load_default_colors(
                    user_defined_default_colors=colors.default_colors_file,
                    palette=colors.default_colors_palette,
                    load_comparison=load_comparison_colors,
                ),
                colors.default_colors_file,
            )
        if colors.color_table_file is not None or colors.default_colors_file is not None:
            colors = replace(
                colors,
                color_table=color_table,
                color_table_file=None,
                default_colors=default_colors,
                default_colors_file=None,
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
            AnnotationOptions(
                sets=tuple(
                    read_annotation_table(
                        annotations.table_file,
                        mode=mode,
                    )
                )
            ),
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
    """Memoized Circular companion FASTA records shared by a request batch."""

    paths: tuple[str | None, ...]
    _records: tuple[tuple[SeqRecord, ...], ...] | None = None

    def cache_specs(
        self,
    ) -> tuple[
        tuple[Hashable, frozenset[PreparedResourceIdentity]] | None,
        ...,
    ]:
        specs: list[
            tuple[Hashable, frozenset[PreparedResourceIdentity]] | None
        ] = []
        for path in self.paths:
            if path is None:
                specs.append(None)
                continue
            identity = prepared_resource_identity(path)
            specs.append(
                None
                if identity is None
                else (
                    ("parsed-source-v1", "comparison-fasta", identity),
                    frozenset({identity}),
                )
            )
        return tuple(specs)

    def load(self) -> tuple[tuple[SeqRecord, ...], ...]:
        if self._records is None:
            loaded: list[tuple[SeqRecord, ...]] = []
            for path, cache_spec in zip(
                self.paths,
                self.cache_specs(),
                strict=True,
            ):
                if not path:
                    loaded.append(())
                    continue

                def parse(path: str = path) -> tuple[SeqRecord, ...]:
                    return tuple(SeqIO.parse(path, "fasta"))

                if cache_spec is None:
                    loaded.append(parse())
                    continue
                key, identities = cache_spec
                loaded.append(
                    get_or_build_parsed_source(
                        key,
                        identities,
                        parse,
                        publish=lambda value: bool(value),
                    )
                )
            self._records = tuple(loaded)
        return self._records


@dataclass(frozen=True)
class PreparedDiagramInputs:
    """Resolved values reused by record loading, drawing, and metadata."""

    features: ResolvedFeatureInputs
    gff_candidate_features: tuple[str, ...]
    gff_keep_all_features: bool
    comparison_sequences: _ComparisonSequenceSources | None = None


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
            detached_manifest: Mapping[str, Any] | None = None
            if manifest is not None:
                if not isinstance(manifest, Mapping):
                    raise ValidationError("protein_identity_manifest must be an object.")
                detached_manifest = copy.deepcopy(dict(manifest))
            if raw_entries or derived_entries or detached_manifest is not None:
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
                validate_protein_identity_manifest(detached_manifest)
            if protein_entries and detached_manifest is None:
                raise ValidationError(
                    "Current protein LOSATP artifacts require protein_identity_manifest."
                )
            if detached_manifest is not None:
                for index, entry in enumerate(protein_entries):
                    _record_request_render_diagnostic_metric(
                        "proteinRawReferenceValidationCount"
                    )
                    if not validate_protein_raw_entry_references(
                        entry,
                        detached_manifest,
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
                    detached_manifest,
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
        object.__setattr__(self, "protein_identity_manifest", detached_manifest)


@dataclass(frozen=True)
class PreparedDiagramRequest:
    """A validated request with normalized records and its SVG drawing."""

    mode: Literal["circular", "linear"]
    request: DiagramRequest
    records: tuple[SeqRecord, ...]
    drawing: Drawing
    inputs: PreparedDiagramInputs | None = None
    linear_metadata: LinearDiagramMetadata | None = None
    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()
    losat_derived_cache_entries: tuple[Mapping[str, Any], ...] = ()
    protein_identity_manifest: Mapping[str, Any] | None = None


@dataclass(frozen=True)
class RequestRenderResult:
    """Files and normalized inputs produced by one request render."""

    mode: Literal["circular", "linear"]
    request: DiagramRequest
    records: tuple[SeqRecord, ...]
    drawing: Drawing
    output_paths: tuple[Path, ...]
    interactive_context: InteractiveSvgContext | None = None
    linear_metadata: LinearDiagramMetadata | None = None
    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()
    losat_derived_cache_entries: tuple[Mapping[str, Any], ...] = ()
    protein_identity_manifest: Mapping[str, Any] | None = None


@dataclass(frozen=True)
class PreparedCircularBatchRequest:
    """A validated Circular batch with one prepared diagram per record."""

    request: CircularBatchRequest
    records: tuple[SeqRecord, ...]
    items: tuple[PreparedDiagramRequest, ...]
    inputs: PreparedDiagramInputs | None = None

    @property
    def mode(self) -> Literal["circular"]:
        return "circular"


@dataclass(frozen=True)
class CircularBatchRenderResult:
    """Saved outputs for an explicit separate-diagram Circular batch."""

    request: CircularBatchRequest
    records: tuple[SeqRecord, ...]
    items: tuple[RequestRenderResult, ...]

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

    @property
    def mode(self) -> Literal["circular"]:
        return "circular"

    def preflight_outputs(self) -> None:
        """Validate every materialized output before diagram construction."""

        _preflight_render_output(self.request.output)

    def build(self) -> Drawing:
        shared_kwargs = (
            {"_resolved_feature_inputs": self.inputs.features}
            if self.inputs is not None
            else {}
        )
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
                **depth_kwargs,
            )
        return build_circular_multi_diagram(
            self.records,
            options=self.request.options,
            layout=self.layout,
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

    @property
    def mode(self) -> Literal["circular"]:
        return "circular"

    def preflight_outputs(self) -> None:
        """Validate every materialized batch output before diagram construction."""

        _preflight_circular_batch_outputs(self.request)

    def item_plans(self) -> tuple[CircularRequestPlan, ...]:
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
                        source=InMemoryRecordSource(record),
                        record_key=self.request.records[index].record_key,
                    ),
                ),
                options=item_options,
                output=output,
                grouping="single",
            )
            plans.append(
                CircularRequestPlan(
                    request=item_request,
                    records=(record,),
                    layout=None,
                    precomputed_depth_track_specs=(
                        tuple(normalized_depth[index])
                        if normalized_depth is not None
                        else None
                    ),
                    precomputed_depth_track_count=(
                        logical_depth_count if normalized_depth is not None else None
                    ),
                    inputs=self.inputs,
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
        kwargs: dict[str, Any] = {"options": self.request.options}
        if self.layout is not None:
            kwargs["layout"] = self.layout
        if losatp_cache is not None:
            kwargs["losatp_cache"] = losatp_cache
        if protein_extraction is not None:
            kwargs["protein_extraction"] = protein_extraction
        if self.inputs is not None:
            kwargs["_resolved_feature_inputs"] = self.inputs.features
        built = build_linear_diagram_result(self.records, **kwargs)
        if isinstance(built, LinearDiagramBuildResult):
            return built
        if not isinstance(built, Drawing):
            raise ValidationError(
                "The Linear diagram builder returned an unsupported result."
            )
        options = self.request.options
        return LinearDiagramBuildResult(
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
        or options.protein_blastp_mode != "none"
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
    if default_colors is None:
        default_colors = load_default_colors(
            user_defined_default_colors=(
                colors.default_colors_file
                if colors is not None and colors.default_colors_file
                else ""
            ),
            palette=(
                colors.default_colors_palette
                if colors is not None
                else "default"
            ),
            load_comparison=(
                isinstance(request, LinearDiagramRequest)
                and _linear_request_uses_comparisons(request)
            ),
        )
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
        _ComparisonSequenceSources(
            tuple(options.conservation_fasta_files or ())
        )
        if isinstance(options, CircularDiagramOptions)
        and options.conservation_fasta_files
        else None
    )
    return PreparedDiagramInputs(
        features=resolve_feature_inputs(
            color_table=color_table,
            default_colors=default_colors,
            feature_visibility_table=feature_visibility_table,
        ),
        gff_candidate_features=tuple(sorted(candidate_features)),
        gff_keep_all_features=keep_all_features,
        comparison_sequences=comparison_sequences,
    )


def _normalize_request_records(
    request: DiagramRequest,
    inputs: PreparedDiagramInputs,
) -> ResolvedRecordCollection:
    def resolve() -> ResolvedRecordCollection:
        return resolve_record_inputs(
            request.records,
            record_options=request.record_options,
            gff_candidate_features=inputs.gff_candidate_features,
            gff_keep_all_features=inputs.gff_keep_all_features,
            genbank_loader=load_gbks,
            gff_loader=load_gff_fasta,
        )

    if not prepared_input_cache_active():
        return resolve()
    source_specs = tuple(
        _prepared_source_cache_spec(
            record_input.source,
            gff_candidate_features=inputs.gff_candidate_features,
            gff_keep_all_features=inputs.gff_keep_all_features,
        )
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
        _load_source_records(
            record_input.source,
            gff_candidate_features=inputs.gff_candidate_features,
            gff_keep_all_features=inputs.gff_keep_all_features,
            genbank_loader=load_gbks,
            gff_loader=load_gff_fasta,
        )
    resource_identities = frozenset(
        identity
        for _key, identities in typed_source_specs
        for identity in identities
    )
    key = (
        "resolved-records-v1",
        tuple(source_key for source_key, _identities in typed_source_specs),
        tuple(_record_input_preparation_key(item) for item in request.records),
        _record_collection_preparation_key(request.record_options),
    )
    return get_or_build_resolved_records(key, resource_identities, resolve)


def _selector_preparation_key(value: Any) -> Hashable:
    if value is None:
        return None
    return (value.record_id, value.record_index)


def _region_preparation_key(value: Any) -> Hashable:
    if value is None:
        return None
    return (
        value.record_id,
        value.record_index,
        int(value.start),
        int(value.end),
        bool(value.reverse_complement),
    )


def _presentation_preparation_key(value: RecordPresentation) -> Hashable:
    return (
        value.label,
        value.subtitle,
        bool(value.reverse_complement),
        value.grid_row,
        value.grid_column,
    )


def _record_input_preparation_key(value: RecordInput) -> Hashable:
    return (
        _selector_preparation_key(value.selector),
        value.cardinality.value,
        _region_preparation_key(value.region),
        _presentation_preparation_key(value.presentation),
        value.record_key,
    )


def _record_collection_preparation_key(
    value: RecordCollectionOptions,
) -> Hashable:
    return (
        tuple(_region_preparation_key(region) for region in value.regions),
        tuple(value.labels),
        tuple(value.subtitles),
    )


def _coerce_resolved_collection(
    request: DiagramRequest,
    value: ResolvedRecordCollection | Sequence[SeqRecord],
) -> ResolvedRecordCollection:
    """Keep private test seams compatible while planners consume provenance."""

    if isinstance(value, ResolvedRecordCollection):
        return value
    records = tuple(value)
    if len(records) != len(request.records):
        raise ValidationError(
            "A record resolver without provenance must return one record per input."
        )
    provenance: list[ResolvedRecordProvenance] = []
    for index, (record, record_input) in enumerate(
        zip(records, request.records, strict=True)
    ):
        source = record_input.source
        if isinstance(source, GenBankInputSource):
            source_kind: Literal["genbank", "gff_fasta", "memory"] = "genbank"
            source_paths = (str(source.path),)
        elif isinstance(source, GffFastaInputSource):
            source_kind = "gff_fasta"
            source_paths = (str(source.gff_path), str(source.fasta_path))
        else:
            source_kind = "memory"
            source_paths = ()
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
                region=record_input.region,
                presentation=record_input.presentation,
            )
        )
    return ResolvedRecordCollection(records, tuple(provenance))


def normalize_request_records(request: DiagramRequest) -> tuple[SeqRecord, ...]:
    """Resolve typed record inputs according to their explicit cardinality."""

    inputs = _prepare_diagram_inputs(request)
    return _coerce_resolved_collection(
        request,
        _normalize_request_records(request, inputs),
    ).records


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
        annotations = getattr(record, "annotations", {})
        label = annotations.get("gbdraw_record_label")
        subtitle = annotations.get("gbdraw_record_subtitle")
        source_presentation = provenance.presentation
        materialized.append(
            RecordInput(
                source=InMemoryRecordSource(record),
                presentation=RecordPresentation(
                    label=str(label) if label is not None else None,
                    subtitle=str(subtitle) if subtitle is not None else None,
                    reverse_complement=False,
                    grid_row=source_presentation.grid_row,
                    grid_column=source_presentation.grid_column,
                ),
                record_key=provenance.record_key,
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
            and record_input.region is None
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
                "WARNING: The annotation indicates that record %s is linear. "
                "Are you sure you want to visualize it as circular?",
                record.id,
            )
        elif topology is not None and topology != "circular":
            logger.warning(
                "WARNING: Topology information not available for %s.",
                record.id,
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
    positioned.sort(
        key=lambda item: (
            int(item[1].grid_row or 0),
            int(item[1].grid_column)
            if item[1].grid_column is not None
            else item[0],
            item[0],
        )
    )
    positions = tuple(
        f"#{index + 1}@{presentation.grid_row}"
        for index, presentation in positioned
    )
    return replace(layout, multi_record_positions=positions)


def _linear_layout_with_record_placements(
    request: LinearDiagramRequest,
) -> LinearMultiRecordOptions | None:
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
    positioned.sort(
        key=lambda item: (
            int(item[1].grid_row or 0),
            int(item[1].grid_column)
            if item[1].grid_column is not None
            else item[0],
            item[0],
        )
    )
    positions = tuple(
        f"#{index + 1}@{presentation.grid_row}"
        for index, presentation in positioned
    )
    return replace(layout, multi_record_positions=positions)


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
            load_comparison_colors=False,
        )
        unresolved_request = (
            request
            if resolved_options is request.options
            else replace(request, options=resolved_options)
        )
        inputs = _prepare_diagram_inputs(unresolved_request)
    with _request_render_diagnostic_phase("recordLoad"):
        collection = _coerce_resolved_collection(
            unresolved_request,
            _normalize_request_records(unresolved_request, inputs),
        )
    with _request_render_diagnostic_phase("preparation"):
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
            if _is_materialized_exact_one_request(unresolved_request)
            and not unresolved_request.output.resolve_prefix_from_first_record
            and resolved_layout == unresolved_request.layout
            else projected_request
        )
        _warn_circular_topologies(records)
    return CircularRequestPlan(
        request=materialized_request,
        records=records,
        layout=resolved_layout,
        inputs=inputs,
        provenance=collection.provenance,
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
            load_comparison_colors=False,
        )
        unresolved_request = (
            request
            if resolved_options is request.options
            else replace(request, options=resolved_options)
        )
        inputs = _prepare_diagram_inputs(unresolved_request)
    with _request_render_diagnostic_phase("recordLoad"):
        collection = _coerce_resolved_collection(
            unresolved_request,
            _normalize_request_records(unresolved_request, inputs),
        )
    with _request_render_diagnostic_phase("preparation"):
        records = collection.records
        outputs = (
            unresolved_request.outputs
            if unresolved_request.outputs
            else resolve_circular_batch_outputs(
                unresolved_request.output_policy,
                records,
            )
        )
        materialized_request = (
            unresolved_request
            if _is_materialized_exact_one_request(unresolved_request)
            and unresolved_request.output_policy is None
            else replace(
                unresolved_request,
                records=_materialized_record_inputs(collection),
                outputs=outputs,
                output_policy=None,
                record_options=RecordCollectionOptions(),
            )
        )
        _warn_circular_topologies(records)
    return CircularBatchRequestPlan(
        request=materialized_request,
        records=records,
        inputs=inputs,
        provenance=collection.provenance,
    )


def plan_linear_request(
    request: LinearDiagramRequest,
) -> LinearRequestPlan:
    """Normalize a Linear request into one explicit builder plan."""

    if not isinstance(request, LinearDiagramRequest):
        raise ValidationError("request must be LinearDiagramRequest.")
    with _request_render_diagnostic_phase("preparation"):
        unresolved_request = replace(
            request,
            options=_resolve_request_option_tables(
                request.options,
                mode="linear",
                load_comparison_colors=_linear_request_uses_comparisons(request),
            ),
        )
        inputs = _prepare_diagram_inputs(unresolved_request)
    with _request_render_diagnostic_phase("recordLoad"):
        collection = _coerce_resolved_collection(
            unresolved_request,
            _normalize_request_records(unresolved_request, inputs),
        )
    with _request_render_diagnostic_phase("preparation"):
        projected_request = replace(
            unresolved_request,
            records=_materialized_record_inputs(collection),
            record_options=RecordCollectionOptions(),
        )
        resolved_layout = _linear_layout_with_record_placements(projected_request)
        resolved_options = resolve_linear_options(
            projected_request.options,
            records=collection.records,
            layout=resolved_layout,
        )
        materialized_request = (
            unresolved_request
            if _is_materialized_exact_one_request(unresolved_request)
            and resolved_layout == unresolved_request.layout
            and resolved_options is unresolved_request.options
            else replace(
                projected_request,
                options=resolved_options,
            )
        )
    return LinearRequestPlan(
        request=materialized_request,
        records=collection.records,
        layout=resolved_layout,
        inputs=inputs,
        provenance=collection.provenance,
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
    """Resolve an unresolved typed request to an exact-one projection."""

    return plan_request(request).request


def _source_protein_mode(
    request: LinearDiagramRequest,
    artifacts: CurrentRequestArtifacts,
) -> str:
    requested = str(request.options.protein_blastp_mode or "none")
    if requested != "none":
        return requested
    if artifacts.protein_source_mode in {"pairwise", "orthogroup", "collinear"}:
        return str(artifacts.protein_source_mode)
    for entry in artifacts.losat_derived_cache_entries:
        mode = str(entry.get("mode") or "")
        if mode in {"pairwise", "orthogroup", "collinear"}:
            return mode
    return "orthogroup" if request.options.orthogroups is not None else "pairwise"


def _extract_linear_request_proteins(
    request: LinearDiagramRequest,
    records: tuple[SeqRecord, ...],
    inputs: PreparedDiagramInputs,
) -> ProteinExtractionResult:
    record_instance_keys = tuple(
        record_input.record_key or f"record-{index + 1}"
        for index, record_input in enumerate(request.records)
    )
    extraction = extract_web_stable_cds_proteins(
        records,
        record_instance_keys=record_instance_keys,
        record_source_ids=tuple(str(record.id) for record in records),
        record_selectors=tuple(
            record_input.selector.raw if record_input.selector is not None else None
            for record_input in request.records
        ),
        regions=tuple(
            record_input.region.raw if record_input.region is not None else None
            for record_input in request.records
        ),
        feature_visibility_rules=inputs.features.feature_visibility_rules,
    )
    if extraction.identity_manifest is None:
        raise ValidationError("Protein extraction did not produce an identity manifest.")
    return extraction


def _prepare_linear_artifacts(
    request: LinearDiagramRequest,
    records: tuple[SeqRecord, ...],
    artifacts: CurrentRequestArtifacts,
    inputs: PreparedDiagramInputs,
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
        or request.options.protein_blastp_mode != "none"
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

    extraction = _extract_linear_request_proteins(request, records, inputs)
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
        threads_per_job=request.options.losatp_threads or "auto",
    )
    return _PreparedLinearArtifacts(
        cache=cache,
        extraction=extraction,
        nucleotide_entries=nucleotide_entries,
        passthrough_derived_entries=passthrough_derived,
        source_mode=_source_protein_mode(request, artifacts),
    )


def _comparison_frame_payload(
    frame: DataFrame,
    *,
    pair_index: int,
    query_index: int,
    subject_index: int,
) -> dict[str, Any]:
    missing = [column for column in COMPARISON_COLUMNS if column not in frame.columns]
    if missing:
        raise ValidationError(
            "Protein comparison is missing canonical column(s): "
            + ", ".join(missing)
        )
    tsv = frame.loc[:, list(COMPARISON_COLUMNS)].to_csv(
        sep="\t",
        header=False,
        index=False,
        lineterminator="\n",
    )
    rows = json.loads(frame.to_json(orient="records"))
    return {
        "pair_index": int(pair_index),
        "query_index": int(query_index),
        "subject_index": int(subject_index),
        "tsv": tsv,
        "rows": rows,
        "hit_count": len(frame),
    }


def _resolved_protein_pair_payloads(
    metadata: LinearDiagramMetadata | None,
    request: LinearDiagramRequest,
) -> list[dict[str, Any]]:
    direct = metadata.protein_comparisons if metadata is not None else None
    if direct is None:
        direct = request.options.protein_comparisons
    explicit = metadata.linear_comparisons if metadata is not None else None

    result: list[dict[str, Any]] = []
    for pair_index, frame in enumerate(direct or ()):
        if not isinstance(frame, DataFrame):
            continue
        result.append(
            _comparison_frame_payload(
                frame,
                pair_index=pair_index,
                query_index=pair_index,
                subject_index=pair_index + 1,
            )
        )
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
        result.append(
            _comparison_frame_payload(
                frame,
                pair_index=len(result),
                query_index=int(comparison.query_record_index),
                subject_index=int(comparison.subject_record_index),
            )
        )
    return result


def _collinearity_parameter_identity(
    params: LosslessCollinearityParameters | None,
) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    resolved = params or LosslessCollinearityParameters()
    snapshot = {
        "model": "lossless",
        "parameters": {
            "minAnchors": int(resolved.min_anchors),
            "maxUnitGap": int(resolved.max_unit_gap),
            "maxDiagonalDrift": int(resolved.max_diagonal_drift),
            "maxConflicts": int(resolved.max_conflicts),
            "mergeOrientation": str(resolved.merge_orientation),
        },
    }
    effective = {
        "minAnchors": int(resolved.min_anchors),
        "maxGeneGap": int(resolved.max_unit_gap),
        "maxDiagonalDrift": int(resolved.max_diagonal_drift),
        "maxConflictsInMergeGap": int(resolved.max_conflicts),
        "mergeOrientation": str(resolved.merge_orientation),
    }
    return snapshot, effective


def _build_current_derived_entries(
    metadata: LinearDiagramMetadata | None,
    request: LinearDiagramRequest,
    records: tuple[SeqRecord, ...],
    artifacts: _PreparedLinearArtifacts,
    raw_entries: Sequence[Mapping[str, Any]],
    *,
    requested_values: Mapping[str, Any] | None = None,
) -> tuple[Mapping[str, Any], ...]:
    if artifacts.extraction is None:
        return artifacts.passthrough_derived_entries
    mode = artifacts.source_mode
    if mode not in {"pairwise", "orthogroup", "collinear"}:
        return ()
    manifest = artifacts.extraction.identity_manifest
    if manifest is None:
        raise ValidationError("Derived protein payload requires an identity manifest.")
    pair_payloads = _resolved_protein_pair_payloads(metadata, request)
    orthogroups = (
        metadata.orthogroups
        if metadata is not None
        else request.options.orthogroups
    )
    orthogroup_payload = serialize_orthogroups_payload(
        orthogroups,
        records=records,
    )
    if not pair_payloads and not orthogroup_payload:
        return ()

    record_instances = manifest.record_instances
    instance_keys = tuple(
        record_input.record_key or f"record-{index + 1}"
        for index, record_input in enumerate(request.records)
    )
    protein_raw_entries = [
        entry for entry in raw_entries if is_protein_losat_cache_entry(entry)
    ]
    raw_cache_keys = sorted(
        {
            str(entry.get("key") or "")
            for entry in protein_raw_entries
            if str(entry.get("key") or "")
        }
    )
    identity = {
        "cacheSchema": 3,
        "idEncoding": "runtime-handle-v1",
        "mode": mode,
        # Orthogroup and collinearity inference can consume self, reverse, and
        # hidden non-adjacent comparisons that are not represented by the
        # displayed pair payloads below.  Bind the derived key to every raw
        # schema-4 entry available to that inference run.
        "rawCacheKeys": raw_cache_keys,
        "thresholds": {
            "bitscore": str(request.options.bitscore),
            "evalue": str(request.options.evalue),
            "identity": str(request.options.identity),
            "alignmentLength": str(request.options.alignment_length),
        },
        "records": [
            {
                "recordIndex": index,
                "recordInstanceKey": instance_key,
                "runtimeBindingHash": str(
                    record_instances.get(instance_key, {}).get(
                        "runtimeBindingHash"
                    )
                    or ""
                ),
                "displayBindingHash": str(
                    record_instances.get(instance_key, {}).get(
                        "displayBindingHash"
                    )
                    or ""
                ),
                "viewTransform": {
                    "length": len(records[index].seq),
                    "reverse": bool(
                        request.records[index].presentation.reverse_complement
                        or (
                            request.records[index].region is not None
                            and request.records[index].region.reverse_complement
                        )
                    ),
                },
            }
            for index, instance_key in enumerate(instance_keys)
        ],
        "pairs": [
            {
                "pairIndex": int(pair["pair_index"]),
                "queryIndex": int(pair["query_index"]),
                "subjectIndex": int(pair["subject_index"]),
                "rawCacheKeys": sorted(
                    str(entry.get("key") or "")
                    for entry in protein_raw_entries
                    if entry.get("queryRecordInstanceKey")
                    == instance_keys[int(pair["query_index"])]
                    and entry.get("subjectRecordInstanceKey")
                    == instance_keys[int(pair["subject_index"])]
                ),
            }
            for pair in pair_payloads
        ],
    }
    if mode == "pairwise":
        identity["pairwise"] = {
            "displayMaxHits": int(request.options.protein_blastp_max_hits),
        }
    if mode in {"orthogroup", "collinear"}:
        identity["orthogroup"] = {
            "membershipMode": str(request.options.orthogroup_membership_mode),
            "memberMaxHits": int(request.options.orthogroup_member_max_hits),
            "maxParalogLinksPerOrthogroup": int(
                request.options.collinear_max_paralog_links_per_orthogroup
            ),
        }
    if mode == "collinear":
        parameter_identity, effective_collinearity_params = (
            _collinearity_parameter_identity(
                request.options.collinearity_params
            )
        )
        identity["collinear"] = {
            **effective_collinearity_params,
            "unitMode": str(request.options.collinearity_unit_mode),
            "anchorMode": str(request.options.collinearity_anchor_mode),
            "searchScope": str(request.options.collinearity_search_scope),
            "parameterIdentity": parameter_identity,
        }
    key = hashlib.sha256(
        json.dumps(
            identity,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()
    effective_values: dict[str, Any] = {
        "mode": mode,
        "evalue": request.options.evalue,
        "bitscore": request.options.bitscore,
        "identity": request.options.identity,
        "alignmentLength": request.options.alignment_length,
        "proteinBlastpCandidateLimit": (
            request.options.protein_blastp_candidate_limit
        ),
    }
    if mode == "pairwise":
        effective_values["proteinBlastpMaxHits"] = (
            request.options.protein_blastp_max_hits
        )
    if mode in {"orthogroup", "collinear"}:
        effective_values.update(
            {
                "orthogroupMembershipMode": (
                    request.options.orthogroup_membership_mode
                ),
                "orthogroupMemberMaxHits": (
                    request.options.orthogroup_member_max_hits
                ),
                "collinearMaxParalogLinksPerOrthogroup": (
                    request.options.collinear_max_paralog_links_per_orthogroup
                ),
            }
        )
    if mode == "collinear":
        effective_values.update(
            {
                "collinearityParams": identity["collinear"]["parameterIdentity"],
                "collinearityUnitMode": request.options.collinearity_unit_mode,
                "collinearityAnchorMode": request.options.collinearity_anchor_mode,
                "collinearitySearchScope": request.options.collinearity_search_scope,
                "collinearityColorMode": request.options.collinearity_color_mode,
            }
        )
    requested_snapshot = copy.deepcopy(
        dict(requested_values)
        if requested_values is not None
        else effective_values
    )
    explicit_names = set(requested_snapshot)
    provenance = {
        "requested": requested_snapshot,
        "effective": copy.deepcopy(effective_values),
        "reasons": {
            name: "explicit" if name in explicit_names else "default"
            for name in effective_values
        },
        "actualInvocations": [
            {
                "rawStageIdentity": str(entry.get("key") or ""),
                "program": str(entry.get("program") or "blastp"),
                "outfmt": str(entry.get("outfmt") or "6"),
                "args": [str(arg) for arg in entry.get("args") or ()],
                **(
                    {"toolIdentity": str(entry["toolIdentity"])}
                    if entry.get("toolIdentity") is not None
                    else {}
                ),
            }
            for entry in sorted(
                protein_raw_entries,
                key=lambda item: str(item.get("key") or ""),
            )
        ],
        "upstreamStageIdentities": raw_cache_keys,
        "versions": {
            "gbdraw": __version__,
            "proteinRawCacheSchema": 4,
            "proteinDerivedArtifactSchema": 3,
            "proteinIdentityManifestSchema": 2,
        },
    }
    entry: Mapping[str, Any] = {
        "schema": 3,
        "kind": "derived-losatp-payload",
        "idEncoding": "runtime-handle-v1",
        "key": key,
        "mode": mode,
        "payload": {
            "identity": identity,
            "provenance": provenance,
            "pairs": pair_payloads,
            "orthogroups": orthogroup_payload,
        },
    }
    if orthogroups is not None:
        from gbdraw.session_request_codec import encode_canonical_typed_resource

        entry["payload"]["orthogroupResult"] = json.loads(
            encode_canonical_typed_resource(
                "orthogroupResult",
                orthogroups,
            ).decode("utf-8")
        )
    collinearity_result = (
        metadata.collinearity_result
        if metadata is not None
        else (
            request.options.collinearity_blocks
            if isinstance(request.options.collinearity_blocks, CollinearityResult)
            else None
        )
    )
    if collinearity_result is not None:
        from gbdraw.session_request_codec import encode_canonical_typed_resource

        entry["payload"]["collinearityResult"] = json.loads(
            encode_canonical_typed_resource(
                "result",
                collinearity_result,
            ).decode("utf-8")
        )
    for cached in artifacts.passthrough_derived_entries:
        if str(cached.get("key") or "") != key:
            continue
        cached_payload = cached.get("payload")
        cached_identity = (
            cached_payload.get("identity")
            if isinstance(cached_payload, Mapping)
            else None
        )
        if cached_identity != identity:
            raise ValidationError(
                "Current derived LOSATP artifact identity does not match its key."
            )
        admitted = copy.deepcopy(dict(cached))
        admitted_payload = admitted.get("payload")
        assert isinstance(admitted_payload, dict)
        admitted_payload["provenance"] = provenance
        if mode == "collinear":
            admitted_payload["pairs"] = pair_payloads
        return (admitted,)
    return (entry,)


def build_request_plan_diagram(
    plan: DiagramRequestPlan,
    *,
    artifacts: CurrentRequestArtifacts | None = None,
) -> PreparedDiagramRequest | PreparedCircularBatchRequest:
    """Build a previously resolved plan without loading or planning it again."""

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
                    drawing=item_plan.build(),
                    inputs=item_plan.inputs,
                )
                for item_plan in plan.item_plans()
            )
        return PreparedCircularBatchRequest(
            request=plan.request,
            records=plan.records,
            items=items,
            inputs=plan.inputs,
        )
    request = plan.request
    records = plan.records
    losat_cache_entries: tuple[Mapping[str, Any], ...] = ()
    losat_derived_cache_entries: tuple[Mapping[str, Any], ...] = ()
    protein_identity_manifest: Mapping[str, Any] | None = None
    linear_metadata: LinearDiagramMetadata | None = None
    if isinstance(plan, CircularRequestPlan):
        with _request_render_diagnostic_phase("drawing"):
            drawing = plan.build()
        losat_cache_entries = current_artifacts.losat_cache_entries
        losat_derived_cache_entries = current_artifacts.losat_derived_cache_entries
        protein_identity_manifest = current_artifacts.protein_identity_manifest
    else:
        if plan.inputs is None:  # pragma: no cover - planners always prepare inputs.
            raise ValidationError("Linear request plan has no prepared diagram inputs.")
        with _request_render_diagnostic_phase("comparisonPreparation"):
            linear_artifacts = _prepare_linear_artifacts(
                request,
                records,
                current_artifacts,
                plan.inputs,
            )
        with _request_render_diagnostic_phase("drawing"):
            linear_build = plan.build(
                losatp_cache=linear_artifacts.cache,
                protein_extraction=linear_artifacts.extraction,
            )
        drawing = linear_build.drawing
        linear_metadata = linear_build.metadata
        with _request_render_diagnostic_phase("comparisonPreparation"):
            protein_entries = (
                tuple(linear_artifacts.cache.session_entries())
                if linear_artifacts.cache is not None
                else ()
            )
            losat_cache_entries = (
                *protein_entries,
                *linear_artifacts.nucleotide_entries,
            )
            losat_derived_cache_entries = _build_current_derived_entries(
                linear_metadata,
                request,
                records,
                linear_artifacts,
                losat_cache_entries,
            )
            protein_identity_manifest = (
                linear_artifacts.extraction.identity_manifest.to_dict()
                if linear_artifacts.extraction is not None
                and linear_artifacts.extraction.identity_manifest is not None
                else current_artifacts.protein_identity_manifest
            )
    return PreparedDiagramRequest(
        mode=plan.mode,
        request=request,
        records=records,
        drawing=drawing,
        inputs=plan.inputs,
        linear_metadata=linear_metadata,
        losat_cache_entries=losat_cache_entries,
        losat_derived_cache_entries=losat_derived_cache_entries,
        protein_identity_manifest=protein_identity_manifest,
    )


def build_request_diagram(
    request: DiagramRequest,
    *,
    artifacts: CurrentRequestArtifacts | None = None,
) -> PreparedDiagramRequest | PreparedCircularBatchRequest:
    """Normalize inputs and build a drawing from current typed artifacts."""

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
        return build_interactive_svg_context(
            prepared.records,
            selected_features_set=options.selected_features_set,
            feature_visibility_rules=inputs.features.feature_visibility_rules,
            specific_color_rules=inputs.features.specific_color_rules,
            orthogroups=computed_orthogroups,
            linear_rendered_feature_ids=prepared.mode == "linear",
            annotations=options.annotations,
            mode=prepared.mode,
            comparison_sequence_records=comparison_sequence_records,
            collinearity_search_scope=collinearity_search_scope,
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


def _prepared_key_resource_identities(
    value: object,
) -> frozenset[PreparedResourceIdentity]:
    identities: set[PreparedResourceIdentity] = set()

    def collect(item: object) -> None:
        if isinstance(item, PreparedResourceIdentity):
            identities.add(item)
            return
        if isinstance(item, (tuple, list, frozenset)):
            for child in item:
                collect(child)

    collect(value)
    return frozenset(identities)


def _resource_backed_context_key(value: object) -> Hashable | None:
    if value is None:
        return ("none",)
    identity = prepared_resource_value_identity(value)
    if identity is None:
        return None
    return ("resource", identity)


def _annotation_style_context_key(value: object | None) -> Hashable:
    if value is None:
        return None
    hatch = getattr(value, "hatch", None)
    hatch_key = None
    if hatch is not None:
        hatch_key = (
            hatch.angle,
            hatch.spacing,
            hatch.color,
            hatch.width,
            hatch.cross,
        )
    return (
        value.stroke,
        value.stroke_width,
        tuple(value.stroke_dasharray),
        value.line_cap,
        value.fill,
        value.fill_opacity,
        hatch_key,
        value.label_color,
        value.label_font_size,
        value.label_orientation,
        value.label_position,
        value.label_offset,
    )


def _annotation_target_context_key(value: object) -> Hashable:
    record = getattr(value, "record", None)
    record_key = _selector_preparation_key(record)
    if hasattr(value, "start") and hasattr(value, "end"):
        return (
            "coordinates",
            record_key,
            int(value.start),
            int(value.end),
            value.coordinate_space,
            bool(value.wraps_origin),
            value.out_of_bounds,
        )
    return (
        "features",
        record_key,
        tuple((selector.key, selector.value) for selector in value.selectors),
        value.envelope,
        value.circular_path,
    )


def _annotation_context_key(value: AnnotationOptions | None) -> Hashable | None:
    if value is None:
        return ("none",)
    identity = prepared_resource_value_identity(value)
    if identity is not None:
        return ("resource", identity)
    if value.table is not None:
        table_identity = prepared_resource_value_identity(value.table)
        return (
            None
            if table_identity is None
            else ("resource-table", table_identity)
        )
    if value.table_file is not None:
        return None
    return (
        "inline-annotation-sets-v1",
        tuple(
            (
                annotation_set.id,
                annotation_set.legend_label,
                _annotation_style_context_key(annotation_set.default_style),
                tuple(
                    (
                        annotation.id,
                        _annotation_target_context_key(annotation.target),
                        annotation.label,
                        annotation.mark,
                        annotation.lane,
                        _annotation_style_context_key(annotation.style),
                        annotation.legend_label,
                        tuple(sorted(annotation.metadata.items())),
                    )
                    for annotation in annotation_set.annotations
                ),
            )
            for annotation_set in value.sets
        ),
    )


def _interactive_context_cache_spec(
    prepared: PreparedDiagramRequest,
    *,
    include_feature_catalog: bool,
) -> tuple[Hashable, frozenset[PreparedResourceIdentity]] | None:
    records = tuple(prepared.records)
    membership = prepared_record_membership(records)
    if membership is None:
        return None
    options = prepared.request.options
    inputs = prepared.inputs
    if inputs is None:
        return None
    visibility_key = _resource_backed_context_key(
        inputs.features.feature_visibility_table
    )
    color_key = _resource_backed_context_key(inputs.features.color_table)
    annotation_key = _annotation_context_key(options.annotations)
    if None in {visibility_key, color_key, annotation_key}:
        return None
    orthogroups = _prepared_orthogroups(prepared)
    orthogroup_key = _resource_backed_context_key(orthogroups)
    if (
        orthogroups is not None
        and orthogroup_key is None
        and prepared.linear_metadata is not None
        and prepared.linear_metadata.collinearity_result is not None
    ):
        collinearity_key = _resource_backed_context_key(
            prepared.linear_metadata.collinearity_result
        )
        if collinearity_key is not None:
            orthogroup_key = ("collinearity", collinearity_key)
    derived_keys = tuple(
        (
            str(entry.get("kind") or ""),
            str(entry.get("key") or ""),
            str(entry.get("mode") or ""),
        )
        for entry in prepared.losat_derived_cache_entries
        if str(entry.get("key") or "")
    )
    if orthogroups is not None and orthogroup_key is None:
        if not derived_keys:
            return None
        orthogroup_key = ("derived-artifacts", derived_keys)
    comparison_keys: list[Hashable] = []
    if inputs.comparison_sequences is not None:
        for path, spec in zip(
            inputs.comparison_sequences.paths,
            inputs.comparison_sequences.cache_specs(),
            strict=True,
        ):
            if path is None:
                comparison_keys.append(("none",))
            elif spec is None:
                return None
            else:
                comparison_keys.append(spec[0])
    record_keys = tuple(
        str(
            (getattr(record, "annotations", {}) or {}).get(
                "gbdraw_record_key",
                f"record-{index + 1}",
            )
        )
        for index, record in enumerate(records)
    )
    collinearity_scope = None
    if (
        isinstance(options, LinearDiagramOptions)
        and prepared.linear_metadata is not None
        and prepared.linear_metadata.collinearity_result is not None
        and prepared.linear_metadata.collinearity_result.orthogroups is not None
    ):
        collinearity_scope = str(options.collinearity_search_scope)
    popup_policy = (
        str(prepared.request.output.interactive_metadata_policy),
        bool(include_feature_catalog),
    )
    key = (
        "interactive-context-v1",
        prepared.mode,
        membership,
        tuple(options.selected_features_set or DEFAULT_SELECTED_FEATURES),
        visibility_key,
        color_key,
        annotation_key,
        orthogroup_key,
        derived_keys,
        tuple(comparison_keys),
        record_keys,
        collinearity_scope,
        popup_policy,
    )
    return key, _prepared_key_resource_identities(key)


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
    cache_spec = _interactive_context_cache_spec(
        prepared,
        include_feature_catalog=include_feature_catalog,
    )

    def build() -> InteractiveSvgContext:
        return build_prepared_interactive_context(
            prepared,
            comparison_sequence_records=comparison_sequence_records,
        )

    if cache_spec is None:
        return build()
    key, identities = cache_spec
    return get_or_build_interactive_context(key, identities, build)


def render_request(
    request: DiagramRequest,
    *,
    artifacts: CurrentRequestArtifacts | None = None,
    include_feature_catalog: bool = False,
) -> RequestRenderResult | CircularBatchRenderResult:
    """Build and save one typed request from current typed artifacts."""

    plan = plan_request(request)
    batch_outputs_preflighted = isinstance(plan, CircularBatchRequestPlan)
    plan.preflight_outputs()
    prepared = build_request_plan_diagram(plan, artifacts=artifacts)
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
            _preflight_circular_batch_outputs(prepared.request)
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
    if isinstance(prepared, PreparedCircularBatchRequest):
        return CircularBatchRenderResult(
            request=prepared.request,
            records=prepared.records,
            items=tuple(
                _render_prepared_request(
                    item,
                    comparison_sequence_records=comparison_sequence_records,
                    include_feature_catalog=include_feature_catalog,
                )
                for item in prepared.items
            ),
        )
    return _render_prepared_request(
        prepared,
        comparison_sequence_records=comparison_sequence_records,
        include_feature_catalog=include_feature_catalog,
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

    inputs = prepared.inputs
    if inputs is not None and inputs.comparison_sequences is not None:
        return require_interactive_svg_metadata(inputs.comparison_sequences.load)

    def load_unprepared() -> tuple[tuple[SeqRecord, ...], ...]:
        return tuple(
            tuple(SeqIO.parse(path, "fasta")) if path else ()
            for path in options.conservation_fasta_files or ()
        )

    return require_interactive_svg_metadata(load_unprepared)


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


def _preflight_circular_batch_outputs(
    request: CircularBatchRequest,
) -> None:
    """Reject batch collisions before writing any item."""

    path_groups = tuple(
        (output, _render_output_paths(output))
        for output in request.outputs
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
            "Could not resolve one or more Circular batch output paths."
        ) from exc
    if len(set(path_identities)) != len(path_identities):
        raise ValidationError(
            "Circular batch output requests resolve to duplicate file paths."
        )


def _render_prepared_request(
    prepared: PreparedDiagramRequest,
    *,
    comparison_sequence_records: Sequence[Sequence[SeqRecord]] = (),
    include_feature_catalog: bool = False,
) -> RequestRenderResult:
    """Write one already-planned diagram using its canonical output request."""

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
        )
    return RequestRenderResult(
        mode=prepared.mode,
        request=prepared.request,
        records=prepared.records,
        drawing=prepared.drawing,
        output_paths=tuple(Path(path) for path in paths),
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
    "normalize_request_records",
    "plan_request",
    "plan_circular_batch_request",
    "plan_circular_request",
    "plan_linear_request",
    "render_request",
    "render_prepared_request",
    "resolve_request",
]
