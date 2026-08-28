"""Internal codec for canonical, CLI-independent render request payloads.

Canonical request schemas are intentionally strict: unknown fields are rejected instead of being
silently discarded.  File content is represented by stable resource IDs; this
module does not own session envelopes, base64 encoding, or temporary files.
"""

from __future__ import annotations

from collections.abc import Mapping as MappingABC, Sequence as SequenceABC
from copy import deepcopy
from dataclasses import MISSING, asdict, dataclass, fields, is_dataclass
from io import StringIO
import json
import math
from pathlib import Path
import re
import types
from typing import Any, Literal, Mapping, Sequence, Union, get_args, get_origin, get_type_hints

from Bio import SeqIO  # type: ignore[reportMissingImports]
from pandas import DataFrame, read_csv  # type: ignore[reportMissingImports]

from gbdraw.analysis.collinearity import (  # type: ignore[reportMissingImports]
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityResult,
    LosslessCollinearityParameters,
)
from gbdraw.analysis.protein_colinearity import (  # type: ignore[reportMissingImports]
    OrthogroupMember,
    OrthogroupNameCandidate,
    OrthogroupResult,
    OrthologEdge,
    OrthologPath,
)
from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import RecordSelector
from gbdraw.io.regions import RegionSpec
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.linear_comparison import LinearComparison
from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot, ScalarSpec
from gbdraw.tracks.circular import (
    _InternalCircularTrackSlot,
    _internal_circular_track_slot,
    parse_circular_track_slot,
)
from gbdraw.annotations import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    FeatureSelector,
    FeatureSpan,
    HatchStyle,
    RegionAnnotation,
    RegionAnnotationStyle,
)

from .api.options import (
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    CircularOutputOptions,
    CircularRequestTrackOptions,
    ColorOptions,
    DepthTrackInput,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    LinearOutputOptions,
    LinearRequestTrackOptions,
    resolve_linear_diagram_options,
)
from .api.requests import (
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


CANONICAL_REQUEST_SCHEMA = 6
SUPPORTED_CANONICAL_REQUEST_SCHEMAS = frozenset({1, 2, 5, CANONICAL_REQUEST_SCHEMA})
UNKNOWN_FIELD_POLICY = "reject"


@dataclass(frozen=True)
class _LegacyStandardCollinearityPayload:
    """Reader-only representation of persisted ``kind=standard`` parameters."""

    min_anchors: int
    max_gene_gap: int
    block_merge_gap: int
    singleton_merge_gap: int
    max_diagonal_drift: int
    max_conflicts_in_merge_gap: int
    max_paralog_links_per_orthogroup: int
    gap_penalty: float
    nearby_duplicate_window: int
    score_mode: Literal["constant", "bitscore"]
    constant_anchor_score: float
    min_block_score: float | None
    block_evalue: float | None

    def validate(self) -> None:
        if self.min_anchors <= 0:
            raise ValidationError("collinear_min_anchors must be > 0")
        if self.max_gene_gap < 0:
            raise ValidationError("collinear_max_gene_gap must be >= 0")
        if self.block_merge_gap < 0:
            raise ValidationError("collinear_block_merge_gap must be >= 0")
        if self.singleton_merge_gap < 0:
            raise ValidationError("collinear_singleton_merge_gap must be >= 0")
        if self.max_diagonal_drift < 0:
            raise ValidationError("collinear_max_diagonal_drift must be >= 0")
        if self.max_conflicts_in_merge_gap < 0:
            raise ValidationError(
                "collinear_max_conflicts_in_merge_gap must be >= 0"
            )
        if self.max_paralog_links_per_orthogroup <= 0:
            raise ValidationError(
                "collinear_max_paralog_links_per_orthogroup must be > 0"
            )
        if self.gap_penalty < 0:
            raise ValidationError("collinear_gap_penalty must be >= 0")
        if self.nearby_duplicate_window < 0:
            raise ValidationError("collinear_nearby_duplicate_window must be >= 0")
        if str(self.score_mode) not in {"constant", "bitscore"}:
            raise ValidationError(
                "collinear_score_mode must be one of: constant, bitscore"
            )
        if self.constant_anchor_score <= 0:
            raise ValidationError("collinear_constant_anchor_score must be > 0")
        if self.block_evalue is not None:
            block_evalue = float(self.block_evalue)
            if not math.isfinite(block_evalue) or block_evalue < 0:
                raise ValidationError(
                    "collinear_block_evalue must be a finite value >= 0 or None"
                )

    def to_lossless(self) -> LosslessCollinearityParameters:
        return LosslessCollinearityParameters(
            min_anchors=self.min_anchors,
            max_unit_gap=self.max_gene_gap,
            max_diagonal_drift=self.max_diagonal_drift,
            max_conflicts=self.max_conflicts_in_merge_gap,
            merge_orientation="either",
        )


_RESOURCE_ID_RE = re.compile(r"^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$")
_TOP_LEVEL_FIELDS = frozenset(
    {"schema", "mode", "records", "diagramOptions", "layout", "comparisons", "output"}
)
_TOP_LEVEL_FIELDS_V5 = _TOP_LEVEL_FIELDS | {"grouping"}
_DEFAULT_CIRCULAR_OPTIONS = CircularDiagramOptions()
_DEFAULT_LINEAR_OPTIONS = LinearDiagramOptions()
_SHARED_OPTION_WRONG_MODE_DEFAULTS = {
    "circular": {
        "depth_track_heights": None,
        "pairwise_match_style": "ribbon",
    },
    "linear": {
        "conservation_blast_files": None,
        "conservation_fasta_files": None,
        "conservation_dataframes": None,
        "conservation_reference": "auto",
        "conservation_labels": None,
        "conservation_colors": None,
        "conservation_ring_width": None,
        "conservation_ring_gap": None,
        "keep_full_definition_with_plot_title": False,
        "species": None,
        "strain": None,
    },
}
_LEGACY_SPARSE_COMPARISON_DEFAULTS = {
    "evalue": 1e-5,
    "bitscore": 50.0,
    "identity": 70.0,
    "alignment_length": 0,
}
_LEGACY_SPARSE_FEATURE_TYPES = (
    "CDS",
    "rRNA",
    "tRNA",
    "tmRNA",
    "ncRNA",
    "misc_RNA",
    "repeat_region",
)
_LEGACY_SPARSE_CONFIG_OVERRIDES = {
    "circular": {
        "canvas.show_gc": True,
        "canvas.show_skew": True,
    },
    "linear": {
        "canvas.show_gc": True,
        "canvas.show_skew": True,
        "objects.axis.linear.stroke_color": "gray",
    },
}
_LEGACY_CONFIG_OVERRIDE_KEYS = {
    "depth_tick_interval": "depth_large_tick_interval",
}
_LEGACY_LINEAR_CONFIG_OVERRIDE_VALUES = {
    "label_placement": {
        "on_feature": "above_feature",
    },
    "linear_track_layout": {
        "spreadout": "above",
        "tuckin": "below",
    },
}
_LEGACY_FLAT_CONFIG_OVERRIDE_PATHS = {
    "block_stroke_width": (
        "objects.features.block_stroke_width.short",
        "objects.features.block_stroke_width.long",
    ),
    "block_stroke_color": ("objects.features.block_stroke_color",),
    "circular_axis_stroke_color": ("objects.axis.circular.stroke_color",),
    "circular_axis_stroke_width": (
        "objects.axis.circular.stroke_width.short",
        "objects.axis.circular.stroke_width.long",
    ),
    "linear_axis_stroke_color": ("objects.axis.linear.stroke_color",),
    "linear_axis_stroke_width": (
        "objects.axis.linear.stroke_width.short",
        "objects.axis.linear.stroke_width.long",
    ),
    "line_stroke_color": ("objects.features.line_stroke_color",),
    "line_stroke_width": (
        "objects.features.line_stroke_width.short",
        "objects.features.line_stroke_width.long",
    ),
    "gc_stroke_color": ("objects.gc_content.stroke_color",),
    "gc_content_mode": ("objects.gc_content.mode",),
    "gc_content_min_percent": ("objects.gc_content.min_percent",),
    "gc_content_max_percent": ("objects.gc_content.max_percent",),
    "gc_content_show_axis": ("objects.gc_content.show_axis",),
    "gc_content_show_ticks": ("objects.gc_content.show_ticks",),
    "gc_content_tick_interval": ("objects.gc_content.large_tick_interval",),
    "gc_content_large_tick_interval": ("objects.gc_content.large_tick_interval",),
    "gc_content_small_tick_interval": ("objects.gc_content.small_tick_interval",),
    "gc_content_tick_font_size": ("objects.gc_content.tick_font_size",),
    "gc_content_percent_background_color": (
        "objects.gc_content.percent_background_color",
    ),
    "gc_content_percent_background_opacity": (
        "objects.gc_content.percent_background_opacity",
    ),
    "gc_content_percent_border_color": (
        "objects.gc_content.percent_border_color",
    ),
    "gc_content_percent_border_width": (
        "objects.gc_content.percent_border_width",
    ),
    "show_gc": ("canvas.show_gc",),
    "show_skew": ("canvas.show_skew",),
    "show_depth": ("canvas.show_depth",),
    "depth_color": ("objects.depth.fill_color",),
    "depth_height": ("canvas.linear.depth_height",),
    "depth_min": ("objects.depth.min_depth",),
    "depth_max": ("objects.depth.max_depth",),
    "depth_normalize": ("objects.depth.normalize",),
    "depth_show_axis": ("objects.depth.show_axis",),
    "depth_show_ticks": ("objects.depth.show_ticks",),
    "depth_large_tick_interval": ("objects.depth.large_tick_interval",),
    "depth_small_tick_interval": ("objects.depth.small_tick_interval",),
    "depth_tick_font_size": ("objects.depth.tick_font_size",),
    "depth_share_axis": ("objects.depth.share_axis",),
    "align_center": ("canvas.linear.align_center",),
    "keep_definition_left_aligned": (
        "canvas.linear.keep_definition_left_aligned",
    ),
    "linear_track_layout": ("canvas.linear.track_layout",),
    "linear_track_axis_gap": ("canvas.linear.track_axis_gap",),
    "linear_ruler_on_axis": ("canvas.linear.ruler_on_axis",),
    "track_type": ("canvas.circular.track_type",),
    "strandedness": ("canvas.strandedness",),
    "resolve_overlaps": ("canvas.resolve_overlaps",),
    "label_radius_offset": ("labels.radius_factor",),
    "label_blacklist": ("labels.filtering.blacklist_keywords",),
    "outer_label_x_radius_offset": (
        "labels.unified_adjustment.outer_labels.x_radius_offset",
    ),
    "outer_label_y_radius_offset": (
        "labels.unified_adjustment.outer_labels.y_radius_offset",
    ),
    "inner_label_x_radius_offset": (
        "labels.unified_adjustment.inner_labels.x_radius_offset",
    ),
    "inner_label_y_radius_offset": (
        "labels.unified_adjustment.inner_labels.y_radius_offset",
    ),
    "comparison_height": ("canvas.linear.comparison_height",),
    "font_family": ("objects.text.font_family",),
    "default_cds_height": (
        "canvas.linear.default_cds_height.short",
        "canvas.linear.default_cds_height.long",
    ),
    "gc_height": ("canvas.linear.default_gc_height",),
    "scale_style": ("objects.scale.style",),
    "scale_stroke_color": ("objects.scale.stroke_color",),
    "scale_label_color": ("objects.scale.label_color",),
    "scale_stroke_width": ("objects.scale.stroke_width",),
    "scale_font_size": (
        "objects.scale.font_size.short",
        "objects.scale.font_size.long",
    ),
    "ruler_label_font_size": (
        "objects.scale.ruler_label_font_size.short",
        "objects.scale.ruler_label_font_size.long",
    ),
    "scale_interval": ("objects.scale.interval",),
    "tick_label_font_size": ("objects.ticks.tick_labels.font_size",),
    "pairwise_match_style": ("objects.blast_match.style",),
    "legend_box_size": (
        "objects.legends.color_rect_size.short",
        "objects.legends.color_rect_size.long",
    ),
    "legend_font_size": (
        "objects.legends.font_size.short",
        "objects.legends.font_size.long",
    ),
    "circular_label_spacing": ("labels.spacing.circular",),
    "circular_label_placement": ("labels.circular.placement",),
    "linear_label_spacing": ("labels.spacing.linear",),
    "label_rendering": ("labels.rendering",),
    "label_placement": ("labels.linear.placement",),
    "label_rotation": ("labels.linear.rotation",),
    "linear_definition_font_size": (
        "objects.definition.linear.font_size.short",
        "objects.definition.linear.font_size.long",
    ),
    "linear_definition_show_replicon": (
        "objects.definition.linear.show_replicon",
    ),
    "linear_definition_show_accession": (
        "objects.definition.linear.show_accession",
    ),
    "linear_definition_show_length": ("objects.definition.linear.show_length",),
    "plot_title_font_size": ("objects.definition.circular.plot_title_font_size",),
    "normalize_length": ("canvas.linear.normalize_length",),
}
_LEGACY_FILTERING_RAW_KEYS = {
    "label_whitelist": "whitelist_df",
    "qualifier_priority": "qualifier_priority_df",
    "label_table": "label_override_df",
}

_COMPARISON_SOURCE_FIELDS = frozenset(
    {
        "blast_files",
        "linear_comparisons",
        "protein_comparisons",
        "protein_comparison_pairs",
        "orthogroups",
        "collinearity_blocks",
    }
)
_PIPELINE_FIELDS = (
    "protein_blastp_mode",
    "collinearity_params",
    "collinearity_unit_mode",
    "collinearity_anchor_mode",
    "collinearity_search_scope",
    "collinearity_color_mode",
    "losatp_bin",
    "ncbi_blastp_bin",
    "losatp_threads",
    "protein_blastp_max_hits",
    "protein_blastp_candidate_limit",
    "orthogroup_membership_mode",
    "orthogroup_member_max_hits",
    "collinear_max_paralog_links_per_orthogroup",
    "align_orthogroup_feature",
)
_COMPARISON_FIELDS = _COMPARISON_SOURCE_FIELDS | frozenset(_PIPELINE_FIELDS)

_TABLE_FIELDS = frozenset(
    {
        "feature_visibility_table",
        "label_whitelist_table",
        "qualifier_priority_table",
        "label_override_table",
        "depth_table",
    }
)
_FILE_FIELDS = frozenset(
    {
        "feature_visibility_table_file",
        "label_whitelist_file",
        "qualifier_priority_file",
        "label_override_file",
        "depth_file",
    }
)
_TABLE_SEQUENCE_FIELDS = frozenset({"depth_tables", "conservation_dataframes"})
_FILE_SEQUENCE_FIELDS = frozenset({"depth_files", "conservation_blast_files"})
_OPTIONAL_FILE_SEQUENCE_FIELDS = frozenset({"conservation_fasta_files"})
_TABLE_MATRIX_FIELDS = frozenset({"depth_track_tables"})
_FILE_MATRIX_FIELDS = frozenset({"depth_track_files"})
_DEPTH_COMPATIBILITY_FIELDS = frozenset(
    {
        "depth_table",
        "depth_file",
        "depth_tables",
        "depth_files",
        "depth_track_tables",
        "depth_track_files",
        "depth_track_labels",
        "depth_track_colors",
        "depth_track_heights",
        "depth_track_large_tick_intervals",
        "depth_track_small_tick_intervals",
        "depth_track_tick_font_sizes",
    }
)
_ALL_DEPTH_INPUT_FIELDS = _DEPTH_COMPATIBILITY_FIELDS | {"depth_tracks"}

_TYPED_TREE_CLASSES = {
    cls.__name__: cls
    for cls in (
        CollinearityAnchor,
        CollinearityBlock,
        CollinearityResult,
        OrthogroupMember,
        OrthogroupNameCandidate,
        OrthogroupResult,
        OrthologEdge,
        OrthologPath,
    )
}
_SCHEMA1_TYPED_OPTIONAL_FIELDS = {
    "CollinearityAnchor": frozenset(
        {
            "queryViewFeatureSvgId",
            "subjectViewFeatureSvgId",
            "queryFeatureIndex",
            "subjectFeatureIndex",
        }
    ),
}


class CanonicalRequestCodecError(ValidationError):
    """Base error for canonical request conversion."""


class CanonicalRequestEncodingError(CanonicalRequestCodecError):
    """Raised when a typed request cannot be represented canonically."""


class CanonicalRequestDecodingError(CanonicalRequestCodecError):
    """Raised when a canonical payload cannot produce a typed request."""


@dataclass(frozen=True)
class CanonicalRequestResource:
    """One resource referenced by a canonical request payload."""

    resource_id: str
    kind: str
    name: str
    source_path: Path | None = None
    content: bytes | None = None

    def __post_init__(self) -> None:
        if not _RESOURCE_ID_RE.fullmatch(self.resource_id):
            raise CanonicalRequestEncodingError(
                f"Invalid canonical resource ID: {self.resource_id!r}."
            )
        if not isinstance(self.kind, str) or not self.kind:
            raise CanonicalRequestEncodingError("Canonical resource kind is required.")
        if not isinstance(self.name, str) or not self.name:
            raise CanonicalRequestEncodingError("Canonical resource name is required.")
        if (self.source_path is None) == (self.content is None):
            raise CanonicalRequestEncodingError(
                "A canonical resource must have exactly one path or byte content source."
            )


@dataclass(frozen=True)
class EncodedCanonicalRequest:
    """JSON-compatible request payload plus its out-of-band resources."""

    payload: dict[str, Any]
    resources: tuple[CanonicalRequestResource, ...]


class _ResourceBuilder:
    def __init__(self) -> None:
        self._resources: list[CanonicalRequestResource] = []
        self._ids: set[str] = set()

    def add_path(self, resource_id: str, *, kind: str, value: object) -> str:
        if not isinstance(value, (str, Path)) or not str(value).strip():
            raise CanonicalRequestEncodingError(
                f"Resource {resource_id!r} must identify a materialized file."
            )
        path = Path(str(value))
        if not path.is_file():
            raise CanonicalRequestEncodingError(
                f"Canonical request resource is not a file: {path}."
            )
        self._add(
            CanonicalRequestResource(
                resource_id=resource_id,
                kind=kind,
                name=path.name,
                source_path=path,
            )
        )
        return resource_id

    def add_bytes(
        self,
        resource_id: str,
        *,
        kind: str,
        name: str,
        content: bytes,
    ) -> str:
        self._add(
            CanonicalRequestResource(
                resource_id=resource_id,
                kind=kind,
                name=name,
                content=content,
            )
        )
        return resource_id

    def _add(self, resource: CanonicalRequestResource) -> None:
        if resource.resource_id in self._ids:
            raise CanonicalRequestEncodingError(
                f"Duplicate canonical resource ID: {resource.resource_id}."
            )
        self._ids.add(resource.resource_id)
        self._resources.append(resource)

    def result(self) -> tuple[CanonicalRequestResource, ...]:
        return tuple(self._resources)


def encode_canonical_request(request: DiagramRequest) -> EncodedCanonicalRequest:
    """Encode one typed request and normalize all conversion failures."""

    try:
        return _encode_canonical_request(request)
    except CanonicalRequestCodecError:
        raise
    except (TypeError, ValueError, ValidationError) as exc:
        raise CanonicalRequestEncodingError(
            f"Typed request could not be converted: {exc}"
        ) from exc


def _encode_canonical_request(request: DiagramRequest) -> EncodedCanonicalRequest:
    """Encode without embedding files or changing the session version."""

    if not isinstance(
        request,
        (CircularDiagramRequest, CircularBatchRequest, LinearDiagramRequest),
    ):
        raise CanonicalRequestEncodingError("Unsupported typed diagram request.")
    unresolved_reasons: list[str] = []
    if request.record_options != RecordCollectionOptions():
        unresolved_reasons.append("collection-level record transforms")
    if isinstance(request, CircularBatchRequest) and request.output_policy is not None:
        unresolved_reasons.append("Circular batch output policy")
    if isinstance(request, CircularBatchRequest) and any(
        output.resolve_prefix_from_first_record for output in request.outputs
    ):
        unresolved_reasons.append("record-derived output prefix")
    if isinstance(request, CircularDiagramRequest):
        if request.output.resolve_prefix_from_first_record:
            unresolved_reasons.append("record-derived output prefix")
        if request.options.conservation_table_file is not None:
            unresolved_reasons.append("Circular comparison table")
    if isinstance(request, LinearDiagramRequest):
        if request.output.resolve_prefix_from_first_record:
            unresolved_reasons.append("record-derived output prefix")
        if request.options.comparison_table_file is not None:
            unresolved_reasons.append("Linear comparison table")
    if unresolved_reasons:
        raise CanonicalRequestEncodingError(
            "Canonical schema 6 cannot encode unresolved request transforms; "
            "call gbdraw.api.resolve_request() before encoding (unresolved: "
            + ", ".join(unresolved_reasons)
            + ")."
        )
    _validate_dataclass_contract(request.options, path="diagramOptions", error="encode")
    if isinstance(request, CircularBatchRequest):
        for index, output in enumerate(request.outputs, start=1):
            _validate_dataclass_contract(
                output,
                path=f"output[{index - 1}]",
                error="encode",
            )
    else:
        _validate_dataclass_contract(request.output, path="output", error="encode")
    if getattr(request, "layout", None) is not None:
        _validate_dataclass_contract(request.layout, path="layout", error="encode")

    resources = _ResourceBuilder()
    mode: Literal["circular", "linear"] = (
        "circular"
        if isinstance(request, (CircularDiagramRequest, CircularBatchRequest))
        else "linear"
    )
    grouping = (
        request.grouping
        if isinstance(request, (CircularDiagramRequest, CircularBatchRequest))
        else "single"
    )
    payload = {
        "schema": CANONICAL_REQUEST_SCHEMA,
        "mode": mode,
        "grouping": grouping,
        "records": [
            _encode_record(record, index=index, resources=resources)
            for index, record in enumerate(request.records, start=1)
        ],
        "diagramOptions": _encode_diagram_options(
            request.options,
            record_count=len(request.records),
            resources=resources,
        ),
        "layout": _encode_layout(request),
        "comparisons": _encode_comparisons(request.options, mode=mode, resources=resources),
        "output": (
            [_encode_output(output) for output in request.outputs]
            if isinstance(request, CircularBatchRequest)
            else _encode_output(request.output)
        ),
    }
    return EncodedCanonicalRequest(payload=payload, resources=resources.result())


def decode_canonical_request(
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    output_directory: str | Path,
) -> DiagramRequest:
    """Decode every supported canonical schema and normalize validation failures."""

    try:
        return _decode_canonical_request(
            payload,
            resource_paths=resource_paths,
            output_directory=output_directory,
        )
    except CanonicalRequestCodecError:
        raise
    except (TypeError, ValueError, ValidationError) as exc:
        raise CanonicalRequestDecodingError(
            f"Canonical request could not be converted: {exc}"
        ) from exc


def _decode_canonical_request(
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    output_directory: str | Path,
) -> DiagramRequest:
    """Decode a supported schema using caller-owned materialized resources."""

    top = _object(
        payload,
        path="renderRequest",
        required={"schema"},
        exact=False,
    )
    schema = top["schema"]
    if not isinstance(schema, int) or isinstance(schema, bool):
        raise CanonicalRequestDecodingError("renderRequest.schema must be an integer.")
    if schema not in SUPPORTED_CANONICAL_REQUEST_SCHEMAS:
        raise CanonicalRequestDecodingError(
            f"Unsupported canonical request schema: {schema!r}."
        )
    _require_exact_fields(
        top,
        path="renderRequest",
        required=set(
            _TOP_LEVEL_FIELDS_V5
            if schema >= 5
            else _TOP_LEVEL_FIELDS
        ),
    )
    mode = top["mode"]
    if mode not in {"circular", "linear"}:
        raise CanonicalRequestDecodingError(
            f"Unsupported canonical request mode: {mode!r}."
        )
    if output_directory is None or not isinstance(output_directory, (str, Path)):
        raise CanonicalRequestDecodingError(
            "A replay output directory must be supplied by the caller."
        )
    if not str(output_directory).strip():
        raise CanonicalRequestDecodingError(
            "A replay output directory must be supplied by the caller."
        )

    raw_records = top["records"]
    if not isinstance(raw_records, list) or not raw_records:
        raise CanonicalRequestDecodingError(
            "renderRequest.records must be a non-empty array."
        )
    records = tuple(
        _decode_record(value, index=index, schema=schema, resource_paths=resource_paths)
        for index, value in enumerate(raw_records, start=1)
    )
    options_kwargs = _decode_diagram_options(
        top["diagramOptions"],
        mode=mode,
        schema=schema,
        resource_paths=resource_paths,
    )
    comparison_kwargs = _decode_comparisons(
        top["comparisons"], mode=mode, schema=schema, resource_paths=resource_paths
    )
    options_type = (
        CircularDiagramOptions if mode == "circular" else LinearDiagramOptions
    )
    options = options_type(**options_kwargs, **comparison_kwargs)
    _validate_dataclass_contract(options, path="diagramOptions", error="decode")
    grouping = top.get("grouping")
    if schema >= 5:
        allowed_groupings = (
            {"single", "grid", "batch"}
            if mode == "circular"
            else {"single"}
        )
        if grouping not in allowed_groupings:
            raise CanonicalRequestDecodingError(
                f"Unsupported {mode} request grouping: {grouping!r}."
            )

    if mode == "linear":
        output = _decode_output(top["output"], output_directory=output_directory)
        if schema == 1:
            layout_payload = _object(top["layout"], path="renderRequest.layout")
            if layout_payload:
                raise CanonicalRequestDecodingError(
                    "A schema 1 Linear canonical request must have an empty layout object."
                )
            linear_layout = None
        else:
            linear_layout = _decode_linear_layout(top["layout"], schema=schema)
        return LinearDiagramRequest(
            records=records,
            options=options,
            layout=linear_layout,
            output=output,
        )

    layout = _decode_circular_layout(top["layout"], schema=schema)
    if schema in {1, 2}:
        grouping = "grid" if layout is not None or len(records) > 1 else "single"
    if grouping == "batch":
        if layout is not None:
            raise CanonicalRequestDecodingError(
                "A Circular batch canonical request cannot define a grid layout."
            )
        return CircularBatchRequest(
            records=records,
            options=options,
            outputs=_decode_batch_outputs(
                top["output"],
                output_directory=output_directory,
                record_count=len(records),
            ),
        )
    output = _decode_output(top["output"], output_directory=output_directory)
    return CircularDiagramRequest(
        records=records,
        options=options,
        layout=layout,
        output=output,
        grouping=grouping,
    )


def _encode_record(
    record: RecordInput,
    *,
    index: int,
    resources: _ResourceBuilder,
) -> dict[str, Any]:
    source = record.source
    if isinstance(source, GenBankInputSource):
        resource_id = resources.add_path(
            f"record-{index}-genbank", kind="genbank", value=source.path
        )
        source_payload = {"kind": "genbank", "resourceId": resource_id}
    elif isinstance(source, GffFastaInputSource):
        gff_id = resources.add_path(
            f"record-{index}-gff3", kind="gff3", value=source.gff_path
        )
        fasta_id = resources.add_path(
            f"record-{index}-fasta", kind="fasta", value=source.fasta_path
        )
        source_payload = {
            "kind": "gffFasta",
            "gffResourceId": gff_id,
            "fastaResourceId": fasta_id,
        }
    elif isinstance(source, InMemoryRecordSource):
        stream = StringIO()
        try:
            serializable_record = deepcopy(source.record)
            serializable_record.annotations.setdefault("molecule_type", "DNA")
            SeqIO.write((serializable_record,), stream, "genbank")
        except Exception as exc:
            raise CanonicalRequestEncodingError(
                f"Could not serialize in-memory record {index} as GenBank."
            ) from exc
        resource_id = resources.add_bytes(
            f"record-{index}-genbank",
            kind="genbank",
            name=f"record-{index}.gbk",
            content=stream.getvalue().encode("utf-8"),
        )
        source_payload = {"kind": "genbank", "resourceId": resource_id}
    else:  # pragma: no cover - RecordInput validates its source union.
        raise CanonicalRequestEncodingError(f"Unsupported record source at index {index}.")

    presentation = record.presentation
    return {
        "recordKey": record.record_key or f"record-{index}",
        "cardinality": record.cardinality.value,
        "source": source_payload,
        "selector": _encode_selector(record.selector),
        "region": _encode_region(record.region),
        "presentation": {
            "label": presentation.label,
            "subtitle": presentation.subtitle,
            "reverseComplement": presentation.reverse_complement,
            "gridRow": presentation.grid_row,
            "gridColumn": presentation.grid_column,
        },
    }


def _decode_record(
    value: object,
    *,
    index: int,
    schema: int,
    resource_paths: Mapping[str, str | Path],
) -> RecordInput:
    path = f"renderRequest.records[{index - 1}]"
    required = {"source", "selector", "region", "presentation"}
    if schema >= 2:
        required.add("recordKey")
    if schema >= 6:
        required.add("cardinality")
    item = _object(value, path=path, required=required)
    source_payload = _object(
        item["source"], path=f"{path}.source", required={"kind"}, exact=False
    )
    kind = source_payload["kind"]
    if kind == "genbank":
        _require_exact_fields(
            source_payload, path=f"{path}.source", required={"kind", "resourceId"}
        )
        source = GenBankInputSource(
            _resolve_resource(
                source_payload["resourceId"],
                path=f"{path}.source.resourceId",
                resource_paths=resource_paths,
            )
        )
    elif kind == "gffFasta":
        _require_exact_fields(
            source_payload,
            path=f"{path}.source",
            required={"kind", "gffResourceId", "fastaResourceId"},
        )
        source = GffFastaInputSource(
            _resolve_resource(
                source_payload["gffResourceId"],
                path=f"{path}.source.gffResourceId",
                resource_paths=resource_paths,
            ),
            _resolve_resource(
                source_payload["fastaResourceId"],
                path=f"{path}.source.fastaResourceId",
                resource_paths=resource_paths,
            ),
        )
    else:
        raise CanonicalRequestDecodingError(
            f"Unsupported record source kind at {path}.source: {kind!r}."
        )

    presentation_payload = _object(
        item["presentation"],
        path=f"{path}.presentation",
        required={"label", "subtitle", "reverseComplement", "gridRow", "gridColumn"},
    )
    presentation = RecordPresentation(
        label=_optional_string(presentation_payload["label"], f"{path}.presentation.label"),
        subtitle=_optional_string(
            presentation_payload["subtitle"], f"{path}.presentation.subtitle"
        ),
        reverse_complement=_boolean(
            presentation_payload["reverseComplement"],
            f"{path}.presentation.reverseComplement",
        ),
        grid_row=_optional_integer(
            presentation_payload["gridRow"], f"{path}.presentation.gridRow"
        ),
        grid_column=_optional_integer(
            presentation_payload["gridColumn"], f"{path}.presentation.gridColumn"
        ),
    )
    return RecordInput(
        source=source,
        cardinality=(
            RecordCardinality(item["cardinality"])
            if schema >= 6
            else RecordCardinality.EXACTLY_ONE
        ),
        selector=_decode_selector(item["selector"], path=f"{path}.selector"),
        region=_decode_region(item["region"], path=f"{path}.region"),
        presentation=presentation,
        record_key=(
            _required_string(item["recordKey"], f"{path}.recordKey")
            if schema >= 2
            else f"record-{index}"
        ),
    )


def _encode_selector(selector: RecordSelector | None) -> dict[str, Any] | None:
    if selector is None:
        return None
    if selector.record_id is not None:
        return {"kind": "recordId", "value": selector.record_id}
    return {"kind": "recordIndex", "index": selector.record_index}


def _decode_selector(value: object, *, path: str) -> RecordSelector | None:
    if value is None:
        return None
    selector = _object(value, path=path, required={"kind"}, exact=False)
    kind = selector["kind"]
    if kind == "recordId":
        _require_exact_fields(selector, path=path, required={"kind", "value"})
        record_id = selector["value"]
        if not isinstance(record_id, str) or not record_id:
            raise CanonicalRequestDecodingError(f"{path}.value must be a non-empty string.")
        return RecordSelector(raw=record_id, record_id=record_id, record_index=None)
    if kind == "recordIndex":
        _require_exact_fields(selector, path=path, required={"kind", "index"})
        record_index = _integer(selector["index"], f"{path}.index")
        if record_index < 0:
            raise CanonicalRequestDecodingError(f"{path}.index must be non-negative.")
        return RecordSelector(
            raw=f"#{record_index + 1}", record_id=None, record_index=record_index
        )
    raise CanonicalRequestDecodingError(f"Unsupported selector kind at {path}: {kind!r}.")


def _encode_region(region: RegionSpec | None) -> dict[str, Any] | None:
    if region is None:
        return None
    selector: dict[str, Any] | None
    if region.record_id is not None:
        selector = {"kind": "recordId", "value": region.record_id}
    elif region.record_index is not None:
        selector = {"kind": "recordIndex", "index": region.record_index}
    else:
        selector = None
    return {
        "selector": selector,
        "start": region.start,
        "end": region.end,
        "reverseComplement": region.reverse_complement,
    }


def _decode_region(value: object, *, path: str) -> RegionSpec | None:
    if value is None:
        return None
    region = _object(
        value,
        path=path,
        required={"selector", "start", "end", "reverseComplement"},
    )
    selector = _decode_selector(region["selector"], path=f"{path}.selector")
    start = _integer(region["start"], f"{path}.start")
    end = _integer(region["end"], f"{path}.end")
    reverse = _boolean(region["reverseComplement"], f"{path}.reverseComplement")
    selector_text = selector.raw if selector is not None else ""
    coordinate_text = f"{start}-{end}" + (":rc" if reverse else "")
    raw = f"{selector_text}:{coordinate_text}" if selector_text else coordinate_text
    return RegionSpec(
        raw=raw,
        file_selector=None,
        record_id=selector.record_id if selector else None,
        record_index=selector.record_index if selector else None,
        start=start,
        end=end,
        reverse_complement=reverse,
    )


def _encode_layout(request: DiagramRequest) -> dict[str, Any]:
    if isinstance(request, CircularBatchRequest):
        return {}
    if isinstance(request, LinearDiagramRequest):
        if request.layout is None:
            return {}
        layout: dict[str, Any] = {"recordGapPx": request.layout.record_gap_px}
        if request.layout.multi_record_positions is not None:
            layout["multiRecordPositions"] = list(
                request.layout.multi_record_positions
            )
        return layout
    if not isinstance(request, CircularDiagramRequest) or request.layout is None:
        return {}
    layout = request.layout
    return {
        "multiRecordSizeMode": layout.multi_record_size_mode,
        "multiRecordMinRadiusRatio": layout.multi_record_min_radius_ratio,
        "multiRecordColumnGapRatio": layout.multi_record_column_gap_ratio,
        "multiRecordRowGapRatio": layout.multi_record_row_gap_ratio,
        "multiRecordPositions": (
            list(layout.multi_record_positions)
            if layout.multi_record_positions is not None
            else None
        ),
    }


def _encode_output(output: RenderOutputRequest) -> dict[str, Any]:
    return {
        "prefix": output.output_prefix,
        "formats": list(output.formats),
        "overwrite": False,
        "interactiveMetadataPolicy": output.interactive_metadata_policy,
    }


def _decode_circular_layout(
    value: object,
    *,
    schema: int,
) -> CircularMultiRecordOptions | None:
    layout = _object(value, path="renderRequest.layout")
    if not layout:
        return None
    required = {
        "multiRecordSizeMode",
        "multiRecordMinRadiusRatio",
        "multiRecordColumnGapRatio",
        "multiRecordRowGapRatio",
        "multiRecordPositions",
    }
    _require_exact_fields(layout, path="renderRequest.layout", required=required)
    positions = layout["multiRecordPositions"]
    if positions is not None and (
        not isinstance(positions, list)
        or not all(isinstance(item, str) for item in positions)
    ):
        raise CanonicalRequestDecodingError(
            "renderRequest.layout.multiRecordPositions must be a string array or null."
        )
    size_mode = layout["multiRecordSizeMode"]
    if schema in {1, 2} and size_mode == "sqrt":
        size_mode = "auto"
    result = CircularMultiRecordOptions(
        multi_record_size_mode=size_mode,
        multi_record_min_radius_ratio=layout["multiRecordMinRadiusRatio"],
        multi_record_column_gap_ratio=layout["multiRecordColumnGapRatio"],
        multi_record_row_gap_ratio=layout["multiRecordRowGapRatio"],
        multi_record_positions=tuple(positions) if positions is not None else None,
    )
    _validate_dataclass_contract(result, path="layout", error="decode")
    return result


def _decode_linear_layout(
    value: object,
    *,
    schema: int,
) -> LinearMultiRecordOptions | None:
    layout = _object(value, path="renderRequest.layout")
    if not layout:
        return None
    _require_exact_fields(
        layout,
        path="renderRequest.layout",
        required=(
            {"recordGapPx"}
            if schema >= 6
            else {"recordGapPx", "multiRecordPositions"}
        ),
        optional={"multiRecordPositions"} if schema >= 6 else frozenset(),
    )
    positions = layout.get("multiRecordPositions")
    if positions is not None and (
        not isinstance(positions, list)
        or not all(isinstance(item, str) for item in positions)
    ):
        raise CanonicalRequestDecodingError(
            "renderRequest.layout.multiRecordPositions must be a string array or null."
        )
    result = LinearMultiRecordOptions(
        record_gap_px=layout["recordGapPx"],
        multi_record_positions=tuple(positions) if positions is not None else None,
    )
    _validate_dataclass_contract(result, path="layout", error="decode")
    return result


def _encode_diagram_options(
    options: CircularDiagramOptions | LinearDiagramOptions,
    *,
    record_count: int,
    resources: _ResourceBuilder,
) -> dict[str, Any]:
    if isinstance(options, CircularDiagramOptions):
        default_options = _DEFAULT_CIRCULAR_OPTIONS
    elif isinstance(options, LinearDiagramOptions):
        default_options = _DEFAULT_LINEAR_OPTIONS
    else:
        raise CanonicalRequestEncodingError(
            "diagramOptions must use mode-specific options."
        )
    depth_tracks = _canonical_depth_tracks_for_encoding(
        options,
        record_count=record_count,
    )
    result: dict[str, Any] = {}
    for item in fields(options):
        name = item.name
        if name in _COMPARISON_FIELDS or name in _ALL_DEPTH_INPUT_FIELDS:
            continue
        value = getattr(options, name)
        default = getattr(default_options, name)
        if _same_default(value, default):
            continue
        result[_camel(name)] = _encode_option_value(name, value, resources=resources)
    if depth_tracks is not None:
        result["depthTracks"] = _encode_depth_tracks(
            depth_tracks,
            resources=resources,
        )
    return result


def _decode_diagram_options(
    value: object,
    *,
    mode: Literal["circular", "linear"],
    schema: int,
    resource_paths: Mapping[str, str | Path],
) -> dict[str, Any]:
    payload = _object(value, path="renderRequest.diagramOptions")
    if schema in {1, 2}:
        payload = _migrate_legacy_feature_visibility_fields(payload)
    payload = dict(payload)
    for name, default in _SHARED_OPTION_WRONG_MODE_DEFAULTS[mode].items():
        key = _camel(name)
        if key in payload and payload[key] == default:
            payload.pop(key)
    options_type = (
        CircularDiagramOptions if mode == "circular" else LinearDiagramOptions
    )
    known = {
        _camel(item.name): item.name
        for item in fields(options_type)
        if item.name not in _COMPARISON_FIELDS
    }
    unknown = set(payload) - set(known)
    if unknown:
        raise CanonicalRequestDecodingError(
            "Unknown field(s) at renderRequest.diagramOptions: "
            + ", ".join(sorted(unknown))
            + "."
        )
    decoded = {
        known[key]: _decode_option_value(
            known[key],
            raw,
            mode=mode,
            schema=schema,
            resource_paths=resource_paths,
        )
        for key, raw in payload.items()
    }
    if decoded.get("config_overrides") is not None:
        decoded["config_overrides"] = _decode_config_overrides(
            decoded["config_overrides"],
            mode=mode,
            schema=schema,
        )
    if schema in {1, 2}:
        feature_shapes = dict(decoded.get("feature_shapes") or {})
        feature_shapes.setdefault("repeat_region", "rectangle")
        decoded["feature_shapes"] = feature_shapes
    if schema in {1, 2}:
        _restore_legacy_sparse_defaults(decoded, mode=mode)
    return decoded


def _migrate_legacy_feature_visibility_fields(
    payload: Mapping[str, Any],
) -> dict[str, Any]:
    """Move removed feature-table aliases at the legacy schema boundary."""

    migrated = dict(payload)
    for old_name, current_name in (
        ("featureTable", "featureVisibilityTable"),
        ("featureTableFile", "featureVisibilityTableFile"),
    ):
        if old_name not in migrated:
            continue
        old_value = migrated.pop(old_name)
        current_value = migrated.get(current_name)
        if old_value is not None and current_value is not None:
            raise CanonicalRequestDecodingError(
                f"renderRequest.diagramOptions cannot contain both "
                f"{old_name} and {current_name}."
            )
        if current_name not in migrated or current_value is None:
            migrated[current_name] = old_value
    return migrated


def _decode_config_overrides(
    value: object,
    *,
    mode: Literal["circular", "linear"],
    schema: int,
) -> dict[str, Any]:
    path = "renderRequest.diagramOptions.configOverrides"
    overrides = dict(_object(value, path=path))
    if schema in {1, 2}:
        for old_name, current_name in _LEGACY_CONFIG_OVERRIDE_KEYS.items():
            if old_name not in overrides:
                continue
            old_value = overrides.pop(old_name)
            current_value = overrides.get(current_name)
            if old_value is not None and current_value is not None:
                raise CanonicalRequestDecodingError(
                    f"{path} cannot contain both {old_name} and {current_name}."
                )
            if current_name not in overrides or current_value is None:
                overrides[current_name] = old_value
        if mode == "linear":
            for name, replacements in _LEGACY_LINEAR_CONFIG_OVERRIDE_VALUES.items():
                raw = overrides.get(name)
                if isinstance(raw, str):
                    overrides[name] = replacements.get(
                        raw.strip().lower(),
                        raw,
                    )
    else:
        retired = _retired_config_override(overrides)
        if retired is not None:
            raise CanonicalRequestDecodingError(
                f"{path} contains retired value {retired!r}; schema {schema} "
                "requires canonical config overrides."
            )
    return _migrate_flat_config_overrides(overrides, mode=mode, path=path)


def _migrate_flat_config_overrides(
    overrides: Mapping[str, Any],
    *,
    mode: Literal["circular", "linear"],
    path: str,
) -> dict[str, Any]:
    """Project persisted flat aliases onto schema-derived dotted leaves."""

    legacy_label_paths = {
        "canvas.show_labels",
        "canvas.circular.show_labels",
        "canvas.linear.show_labels",
        "canvas.circular.allow_inner_labels",
    }
    migrated = {
        name: deepcopy(value)
        for name, value in overrides.items()
        if "." in name and name not in legacy_label_paths
    }
    migrated_sources = {name: name for name in migrated}
    consumed: set[str] = set()

    def assign(canonical_path: str, value: object, *, source: str) -> None:
        existing_source = migrated_sources.get(canonical_path)
        if existing_source is not None:
            if (
                existing_source == "gc_content_tick_interval"
                and source == "gc_content_large_tick_interval"
            ):
                migrated[canonical_path] = deepcopy(value)
                migrated_sources[canonical_path] = source
                return
            if migrated[canonical_path] == value:
                return
            raise CanonicalRequestDecodingError(
                f"{path} cannot contain both {existing_source} and {source}."
            )
        migrated[canonical_path] = deepcopy(value)
        migrated_sources[canonical_path] = source

    active_show_names = (
        ("show_labels", "canvas.show_labels", "canvas.circular.show_labels")
        if mode == "circular"
        else ("show_labels", "canvas.show_labels", "canvas.linear.show_labels")
    )
    wrong_show_name = (
        "canvas.linear.show_labels"
        if mode == "circular"
        else "canvas.circular.show_labels"
    )
    if overrides.get(wrong_show_name) is not None:
        raise CanonicalRequestDecodingError(
            f"{path}.{wrong_show_name} is for the other drawing mode."
        )
    consumed.add(wrong_show_name)
    legacy_scope: str | None = None
    legacy_scope_source: str | None = None
    for legacy_name in active_show_names:
        if legacy_name not in overrides:
            continue
        consumed.add(legacy_name)
        value = overrides[legacy_name]
        if value is None:
            continue
        circular_scope, linear_scope = _legacy_mode_label_scopes(
            value,
            mode=mode,
            path=f"{path}.{legacy_name}",
        )
        scope = circular_scope if mode == "circular" else linear_scope
        if legacy_scope is not None and legacy_scope != scope:
            raise CanonicalRequestDecodingError(
                f"{path} cannot contain conflicting {legacy_scope_source} "
                f"and {legacy_name} values."
            )
        legacy_scope = scope
        legacy_scope_source = legacy_name

    inner_label_names = (
        "allow_inner_labels",
        "canvas.circular.allow_inner_labels",
    )
    inner_labels: bool | None = None
    inner_source: str | None = None
    for legacy_name in inner_label_names:
        if legacy_name not in overrides:
            continue
        consumed.add(legacy_name)
        value = overrides[legacy_name]
        if value is None:
            continue
        if not isinstance(value, bool):
            raise CanonicalRequestDecodingError(
                f"{path}.{legacy_name} must be a boolean."
            )
        if inner_labels is not None and inner_labels != value:
            raise CanonicalRequestDecodingError(
                f"{path} cannot contain conflicting {inner_source} "
                f"and {legacy_name} values."
            )
        inner_labels = value
        inner_source = legacy_name

    if mode == "circular":
        if legacy_scope == "outer" and inner_labels:
            legacy_scope = "both"
        elif legacy_scope is None and inner_labels:
            legacy_scope = "both"
        if legacy_scope is not None:
            assign(
                "labels.circular.scope",
                legacy_scope,
                source=legacy_scope_source or inner_source or "legacy labels",
            )
    elif inner_labels:
        raise CanonicalRequestDecodingError(
            f"{path}.{inner_source} is Circular-only."
        )
    elif legacy_scope is not None:
        assign(
            "labels.linear.scope",
            legacy_scope,
            source=legacy_scope_source or "legacy labels",
        )

    for legacy_name, canonical_paths in _LEGACY_FLAT_CONFIG_OVERRIDE_PATHS.items():
        if legacy_name not in overrides:
            continue
        consumed.add(legacy_name)
        value = overrides[legacy_name]
        if value is None:
            continue
        if legacy_name == "label_blacklist" and isinstance(value, str):
            value = [item.strip() for item in value.split(",") if item.strip()]
        for canonical_path in canonical_paths:
            assign(canonical_path, value, source=legacy_name)

    if "label_font_size" in overrides:
        consumed.add("label_font_size")
        value = overrides["label_font_size"]
        if value is not None:
            prefix = (
                "labels.font_size"
                if mode == "circular"
                else "labels.font_size.linear"
            )
            for suffix in ("short", "long"):
                assign(
                    f"{prefix}.{suffix}",
                    value,
                    source="label_font_size",
                )

    if "circular_definition_font_size" in overrides:
        consumed.add("circular_definition_font_size")
        value = overrides["circular_definition_font_size"]
        if value is not None:
            assign(
                "objects.definition.circular.font_size",
                value,
                source="circular_definition_font_size",
            )
            try:
                interval = int(float(value) + 2)
            except (TypeError, ValueError, OverflowError) as exc:
                raise CanonicalRequestDecodingError(
                    f"{path}.circular_definition_font_size must be numeric."
                ) from exc
            assign(
                "objects.definition.circular.interval",
                interval,
                source="circular_definition_font_size",
            )

    if "linear_definition_line_styles" in overrides:
        consumed.add("linear_definition_line_styles")
        line_styles = overrides["linear_definition_line_styles"]
        if line_styles is not None:
            if not isinstance(line_styles, MappingABC):
                raise CanonicalRequestDecodingError(
                    f"{path}.linear_definition_line_styles must be an object."
                )
            for line_kind, style in line_styles.items():
                if not isinstance(style, MappingABC):
                    raise CanonicalRequestDecodingError(
                        f"{path}.linear_definition_line_styles.{line_kind} "
                        "must be an object."
                    )
                for property_name, property_value in style.items():
                    if property_value is None:
                        continue
                    assign(
                        "objects.definition.linear.line_styles."
                        f"{line_kind}.{property_name}",
                        property_value,
                        source=(
                            "linear_definition_line_styles."
                            f"{line_kind}.{property_name}"
                        ),
                    )

    raw_filtering = {
        raw_name: deepcopy(overrides[legacy_name])
        for legacy_name, raw_name in _LEGACY_FILTERING_RAW_KEYS.items()
        if legacy_name in overrides and overrides[legacy_name] is not None
    }
    consumed.update(
        legacy_name
        for legacy_name in _LEGACY_FILTERING_RAW_KEYS
        if legacy_name in overrides
    )
    if raw_filtering:
        assign(
            "labels.filtering.raw",
            raw_filtering,
            source=", ".join(
                legacy_name
                for legacy_name in _LEGACY_FILTERING_RAW_KEYS
                if legacy_name in overrides
                and overrides[legacy_name] is not None
            ),
        )

    for name, value in overrides.items():
        if name not in consumed and "." not in name:
            migrated[name] = deepcopy(value)
    return migrated


def _legacy_mode_label_scopes(
    value: object,
    *,
    mode: Literal["circular", "linear"],
    path: str,
) -> tuple[
    Literal["none", "outer"],
    Literal["none", "all", "first", "orthogroup_top"],
]:
    if isinstance(value, bool):
        return ("outer" if value else "none"), ("all" if value else "none")
    if not isinstance(value, str):
        raise CanonicalRequestDecodingError(
            f"{path} must be a boolean or label policy."
        )
    normalized = value.strip().lower()
    aliases = {
        "true": "all",
        "yes": "all",
        "on": "all",
        "false": "none",
        "no": "none",
        "off": "none",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in {"all", "first", "orthogroup_top", "none"}:
        raise CanonicalRequestDecodingError(
            f"{path} contains unsupported label policy {value!r}."
        )
    if mode == "circular" and normalized in {"first", "orthogroup_top"}:
        raise CanonicalRequestDecodingError(
            f"{path} contains Linear-only label policy {normalized!r}."
        )
    circular_scope = "outer" if normalized == "all" else "none"
    return circular_scope, normalized  # type: ignore[return-value]


def _retired_config_override(overrides: Mapping[str, Any]) -> str | None:
    for old_name in _LEGACY_CONFIG_OVERRIDE_KEYS:
        if old_name in overrides:
            return old_name
    for name, replacements in _LEGACY_LINEAR_CONFIG_OVERRIDE_VALUES.items():
        raw = overrides.get(name)
        if isinstance(raw, str):
            normalized = raw.strip().lower()
            if normalized in replacements:
                return f"{name}={normalized}"
    return None


def _restore_legacy_sparse_defaults(
    decoded: dict[str, Any],
    *,
    mode: Literal["circular", "linear"],
) -> None:
    """Keep supported sparse payloads on the defaults they were saved with.

    Canonical schemas 1 and 2 omitted values equal to the
    shared Python defaults. Fresh requests now use mode-specific profiles, so
    the compatibility reader materializes the historical values before the
    typed request applies the current profile.
    """

    for name, default in _LEGACY_SPARSE_COMPARISON_DEFAULTS.items():
        decoded.setdefault(name, default)
    decoded.setdefault("selected_features_set", _LEGACY_SPARSE_FEATURE_TYPES)

    if decoded.get("config") is not None:
        return
    explicit_overrides = decoded.get("config_overrides")
    if explicit_overrides is not None and not isinstance(
        explicit_overrides, MappingABC
    ):
        return
    overrides = dict(_LEGACY_SPARSE_CONFIG_OVERRIDES[mode])
    overrides.update(
        {
            key: value
            for key, value in (explicit_overrides or {}).items()
            if value is not None
        }
    )
    decoded["config_overrides"] = overrides


def _encode_depth_source(
    value: object,
    *,
    name: str,
    resources: _ResourceBuilder,
) -> dict[str, str]:
    if isinstance(value, DataFrame):
        return _table_ref(name, value, resources=resources)
    if isinstance(value, (str, Path)):
        return _file_ref(name, value, resources=resources)
    raise CanonicalRequestEncodingError(
        f"{name} must be a path or pandas DataFrame."
    )


def _expanded_depth_text_metadata(
    values: Sequence[object] | None,
    *,
    track_count: int,
    field_name: str,
) -> list[str | None]:
    if values is None:
        return [None] * track_count
    items = [str(value).strip() or None for value in values]
    if len(items) == 1:
        return items * track_count
    if len(items) != track_count:
        raise CanonicalRequestEncodingError(
            f"{field_name} count must be one or equal to the number of depth "
            f"tracks ({track_count})."
        )
    return items


def _expanded_depth_numeric_metadata(
    values: Sequence[float | str | None] | None,
    *,
    track_count: int,
    field_name: str,
) -> list[float | None]:
    if values is None:
        return [None] * track_count
    items: list[float | None] = []
    for value in values:
        if value is None or (
            isinstance(value, str)
            and value.strip().lower() in {"", "auto", "none", "null", "-"}
        ):
            items.append(None)
            continue
        try:
            numeric = float(value)
        except (TypeError, ValueError, OverflowError) as exc:
            raise CanonicalRequestEncodingError(
                f"{field_name} values must be finite numbers > 0 or None."
            ) from exc
        if not math.isfinite(numeric) or numeric <= 0:
            raise CanonicalRequestEncodingError(
                f"{field_name} values must be finite numbers > 0 or None."
            )
        items.append(numeric)
    if len(items) == 1:
        return items * track_count
    if len(items) != track_count:
        raise CanonicalRequestEncodingError(
            f"{field_name} count must be one or equal to the number of depth "
            f"tracks ({track_count})."
        )
    return items


def _canonical_depth_tracks_for_encoding(
    options: CircularDiagramOptions | LinearDiagramOptions,
    *,
    record_count: int,
) -> tuple[DepthTrackInput, ...] | None:
    """Project compatibility depth inputs onto the canonical logical-track form."""

    if options.depth_tracks is not None:
        return tuple(options.depth_tracks)

    singular_plural_present = any(
        getattr(options, name) is not None
        for name in ("depth_table", "depth_file", "depth_tables", "depth_files")
    )
    matrix_present = any(
        getattr(options, name, None) is not None
        for name in (
            "depth_track_tables",
            "depth_track_files",
            "depth_track_labels",
            "depth_track_colors",
            "depth_track_heights",
            "depth_track_large_tick_intervals",
            "depth_track_small_tick_intervals",
            "depth_track_tick_font_sizes",
        )
    )
    if singular_plural_present and matrix_present:
        raise CanonicalRequestEncodingError(
            "Singular/plural depth inputs cannot be combined with depth_track_* "
            "compatibility inputs."
        )
    if not singular_plural_present and not matrix_present:
        return None

    if singular_plural_present:
        if options.depth_table is not None and options.depth_file is not None:
            raise CanonicalRequestEncodingError(
                "Pass either depth_table or depth_file, not both."
            )
        singular = (
            options.depth_table
            if options.depth_table is not None
            else options.depth_file
        )
        plural = (
            options.depth_tables
            if options.depth_tables is not None
            else options.depth_files
        )
        if singular is not None and plural is not None:
            raise CanonicalRequestEncodingError(
                "Use depth_table/depth_file or depth_tables/depth_files, not both."
            )
        if options.depth_tables is not None and options.depth_files is not None:
            raise CanonicalRequestEncodingError(
                "Pass either depth_tables or depth_files, not both."
            )
        if singular is not None:
            return (DepthTrackInput(source=singular),)
        sources = tuple(plural or ())
        if not sources:
            raise CanonicalRequestEncodingError(
                "depth_tables/depth_files must include at least one source."
            )
        if len(sources) not in {1, record_count}:
            raise CanonicalRequestEncodingError(
                "depth_tables/depth_files must contain one shared source or one "
                f"source per displayed record ({record_count}); got {len(sources)}."
            )
        source: object = sources[0] if len(sources) == 1 else sources
        return (DepthTrackInput(source=source),)

    table_rows = options.depth_track_tables
    file_rows = options.depth_track_files
    if table_rows is not None and file_rows is not None:
        raise CanonicalRequestEncodingError(
            "Pass either depth_track_tables or depth_track_files, not both."
        )
    raw_rows = table_rows if table_rows is not None else file_rows
    if raw_rows is None:
        raise CanonicalRequestEncodingError(
            "depth_track metadata requires depth_track_tables or depth_track_files."
        )
    rows = [tuple(row) for row in raw_rows]
    if not rows:
        raise CanonicalRequestEncodingError(
            "depth_track_tables/depth_track_files must include at least one row."
        )
    if len(rows) not in {1, record_count}:
        raise CanonicalRequestEncodingError(
            "depth_track_tables/depth_track_files must contain one shared row or "
            f"one row per displayed record ({record_count}); got {len(rows)}."
        )
    track_count = max((len(row) for row in rows), default=0)
    if track_count <= 0:
        raise CanonicalRequestEncodingError(
            "depth_track_tables/depth_track_files must include at least one track."
        )

    labels = _expanded_depth_text_metadata(
        options.depth_track_labels,
        track_count=track_count,
        field_name="depth_track_labels",
    )
    colors = _expanded_depth_text_metadata(
        options.depth_track_colors,
        track_count=track_count,
        field_name="depth_track_colors",
    )
    heights = _expanded_depth_numeric_metadata(
        getattr(options, "depth_track_heights", None),
        track_count=track_count,
        field_name="depth_track_heights",
    )
    large_ticks = _expanded_depth_numeric_metadata(
        options.depth_track_large_tick_intervals,
        track_count=track_count,
        field_name="depth_track_large_tick_intervals",
    )
    small_ticks = _expanded_depth_numeric_metadata(
        options.depth_track_small_tick_intervals,
        track_count=track_count,
        field_name="depth_track_small_tick_intervals",
    )
    tick_fonts = _expanded_depth_numeric_metadata(
        options.depth_track_tick_font_sizes,
        track_count=track_count,
        field_name="depth_track_tick_font_sizes",
    )

    tracks: list[DepthTrackInput] = []
    for track_index in range(track_count):
        per_record = tuple(
            row[track_index] if track_index < len(row) else None
            for row in rows
        )
        source = per_record[0] if len(rows) == 1 else per_record
        tracks.append(
            DepthTrackInput(
                source=source,
                label=labels[track_index],
                color=colors[track_index],
                height=heights[track_index],
                large_tick_interval=large_ticks[track_index],
                small_tick_interval=small_ticks[track_index],
                tick_font_size=tick_fonts[track_index],
            )
        )
    return tuple(tracks)


def _encode_depth_tracks(
    value: object,
    *,
    resources: _ResourceBuilder,
) -> list[dict[str, Any]]:
    tracks = _sequence(value, name="depth_tracks")
    encoded: list[dict[str, Any]] = []
    for track_index, track in enumerate(tracks, start=1):
        if not isinstance(track, DepthTrackInput):
            raise CanonicalRequestEncodingError(
                "depth_tracks must contain DepthTrackInput values."
            )
        source_name = f"depth-tracks-{track_index}-source"
        source = track.source
        if isinstance(source, SequenceABC) and not isinstance(
            source,
            (str, bytes, bytearray),
        ) and not isinstance(source, DataFrame):
            encoded_source: object = [
                (
                    _encode_depth_source(
                        item,
                        name=f"{source_name}-record-{record_index}",
                        resources=resources,
                    )
                    if item is not None
                    else None
                )
                for record_index, item in enumerate(source, start=1)
            ]
        else:
            encoded_source = _encode_depth_source(
                source,
                name=source_name,
                resources=resources,
            )
        encoded.append(
            {
                "source": encoded_source,
                "label": track.label,
                "color": track.color,
                "height": track.height,
                "largeTickInterval": track.large_tick_interval,
                "smallTickInterval": track.small_tick_interval,
                "tickFontSize": track.tick_font_size,
            }
        )
    return encoded


def _decode_depth_source(
    value: object,
    *,
    name: str,
    resource_paths: Mapping[str, str | Path],
) -> str | DataFrame:
    ref = _object(
        value,
        path=f"resource reference {name}",
        required={"resourceId", "representation"},
    )
    representation = ref["representation"]
    if representation == "canonicalTsv":
        return _decode_table_ref(
            ref,
            name=name,
            resource_paths=resource_paths,
        )
    if representation == "file":
        return str(
            _decode_file_ref(
                ref,
                name=name,
                resource_paths=resource_paths,
            )
        )
    raise CanonicalRequestDecodingError(
        f"Unsupported resource representation for {name}: {representation!r}."
    )


def _decode_depth_tracks(
    value: object,
    *,
    resource_paths: Mapping[str, str | Path],
) -> tuple[DepthTrackInput, ...]:
    raw_tracks = _array(
        value,
        path="renderRequest.diagramOptions.depthTracks",
    )
    decoded: list[DepthTrackInput] = []
    required = {
        "source",
        "label",
        "color",
        "height",
        "largeTickInterval",
        "smallTickInterval",
        "tickFontSize",
    }
    for track_index, raw_track in enumerate(raw_tracks, start=1):
        path = f"renderRequest.diagramOptions.depthTracks[{track_index - 1}]"
        track = _object(raw_track, path=path, required=required)
        raw_source = track["source"]
        source_name = f"depth-tracks-{track_index}-source"
        if isinstance(raw_source, list):
            source: object = tuple(
                (
                    _decode_depth_source(
                        item,
                        name=f"{source_name}-record-{record_index}",
                        resource_paths=resource_paths,
                    )
                    if item is not None
                    else None
                )
                for record_index, item in enumerate(raw_source, start=1)
            )
        else:
            source = _decode_depth_source(
                raw_source,
                name=source_name,
                resource_paths=resource_paths,
            )
        decoded.append(
            DepthTrackInput(
                source=source,
                label=_optional_string(track["label"], f"{path}.label"),
                color=_optional_string(track["color"], f"{path}.color"),
                height=track["height"],
                large_tick_interval=track["largeTickInterval"],
                small_tick_interval=track["smallTickInterval"],
                tick_font_size=track["tickFontSize"],
            )
        )
    return tuple(decoded)


def _encode_option_value(
    name: str,
    value: object,
    *,
    resources: _ResourceBuilder,
) -> Any:
    if name == "config":
        if isinstance(value, GbdrawConfig):
            config = value
            value = asdict(config)
            value["labels"]["filtering"] = deepcopy(
                config.labels.filtering.as_dict()
            )
        return _json_value(value, path="diagramOptions.config")
    if name == "config_overrides" and isinstance(value, MappingABC):
        retired = _retired_config_override(value)
        if retired is not None:
            raise CanonicalRequestEncodingError(
                "diagramOptions.configOverrides contains retired value "
                f"{retired!r}; use canonical config overrides."
            )
    if name == "colors":
        return _encode_colors(value, resources=resources)
    if name == "tracks":
        return _encode_tracks(value)
    if name == "annotations":
        return _encode_annotations(value, resources=resources)
    if name == "output":
        return _encode_assembly_output(value)
    if name == "depth_tracks":
        return _encode_depth_tracks(value, resources=resources)
    if name in _TABLE_FIELDS:
        return _table_ref(name, value, resources=resources)
    if name in _FILE_FIELDS:
        return _file_ref(name, value, resources=resources)
    if name in _TABLE_SEQUENCE_FIELDS:
        return [
            _table_ref(f"{name}-{index}", item, resources=resources)
            for index, item in enumerate(_sequence(value, name=name), start=1)
        ]
    if name in _FILE_SEQUENCE_FIELDS:
        return [
            _file_ref(f"{name}-{index}", item, resources=resources)
            for index, item in enumerate(_sequence(value, name=name), start=1)
        ]
    if name in _OPTIONAL_FILE_SEQUENCE_FIELDS:
        return [
            (
                _file_ref(f"{name}-{index}", item, resources=resources)
                if item is not None
                else None
            )
            for index, item in enumerate(_sequence(value, name=name), start=1)
        ]
    if name in _TABLE_MATRIX_FIELDS:
        return _encode_resource_matrix(name, value, table=True, resources=resources)
    if name in _FILE_MATRIX_FIELDS:
        return _encode_resource_matrix(name, value, table=False, resources=resources)
    return _json_value(value, path=f"diagramOptions.{_camel(name)}")


def _decode_option_value(
    name: str,
    value: object,
    *,
    mode: Literal["circular", "linear"],
    schema: int,
    resource_paths: Mapping[str, str | Path],
) -> Any:
    if name == "config":
        return _migrate_legacy_full_config(
            _object(value, path="renderRequest.diagramOptions.config"),
            mode=mode,
            path="renderRequest.diagramOptions.config",
        )
    if name == "colors":
        return _decode_colors(value, resource_paths=resource_paths)
    if name == "tracks":
        return _decode_tracks(value, mode=mode, schema=schema)
    if name == "annotations":
        return _decode_annotations(value, resource_paths=resource_paths)
    if name == "output":
        return _decode_assembly_output(value, mode=mode, schema=schema)
    if name == "depth_tracks":
        return _decode_depth_tracks(value, resource_paths=resource_paths)
    if name in _TABLE_FIELDS:
        return _decode_table_ref(value, name=name, resource_paths=resource_paths)
    if name in _FILE_FIELDS:
        return str(_decode_file_ref(value, name=name, resource_paths=resource_paths))
    if name in _TABLE_SEQUENCE_FIELDS:
        raw = _array(value, path=f"renderRequest.diagramOptions.{_camel(name)}")
        return tuple(
            _decode_table_ref(item, name=f"{name}-{index}", resource_paths=resource_paths)
            for index, item in enumerate(raw, start=1)
        )
    if name in _FILE_SEQUENCE_FIELDS:
        raw = _array(value, path=f"renderRequest.diagramOptions.{_camel(name)}")
        return tuple(
            str(
                _decode_file_ref(
                    item, name=f"{name}-{index}", resource_paths=resource_paths
                )
            )
            for index, item in enumerate(raw, start=1)
        )
    if name in _OPTIONAL_FILE_SEQUENCE_FIELDS:
        raw = _array(value, path=f"renderRequest.diagramOptions.{_camel(name)}")
        return tuple(
            (
                str(
                    _decode_file_ref(
                        item,
                        name=f"{name}-{index}",
                        resource_paths=resource_paths,
                    )
                )
                if item is not None
                else None
            )
            for index, item in enumerate(raw, start=1)
        )
    if name in _TABLE_MATRIX_FIELDS:
        return _decode_resource_matrix(
            value, name=name, table=True, resource_paths=resource_paths
        )
    if name in _FILE_MATRIX_FIELDS:
        return _decode_resource_matrix(
            value, name=name, table=False, resource_paths=resource_paths
        )
    return value


def _migrate_legacy_full_config(
    value: Mapping[str, Any],
    *,
    mode: Literal["circular", "linear"],
    path: str,
) -> dict[str, Any]:
    """Move persisted config aliases onto the current typed schema."""

    migrated = deepcopy(dict(value))
    canvas = migrated.get("canvas")
    if not isinstance(canvas, dict):
        return migrated

    circular_canvas = canvas.get("circular")
    linear_canvas = canvas.get("linear")
    if not isinstance(circular_canvas, dict) or not isinstance(linear_canvas, dict):
        return migrated

    legacy_track_layout = linear_canvas.get("track_layout")
    if isinstance(legacy_track_layout, str):
        normalized_track_layout = legacy_track_layout.strip().lower()
        linear_canvas["track_layout"] = {
            "spreadout": "above",
            "tuckin": "below",
        }.get(normalized_track_layout, legacy_track_layout)

    labels = migrated.get("labels")
    if isinstance(labels, dict):
        linear_labels = labels.get("linear")
        if isinstance(linear_labels, dict):
            legacy_placement = linear_labels.get("placement")
            if (
                isinstance(legacy_placement, str)
                and legacy_placement.strip().lower() == "on_feature"
            ):
                linear_labels["placement"] = "above_feature"

    shared_show = canvas.pop("show_labels", None)
    circular_show = circular_canvas.pop("show_labels", None)
    linear_show = linear_canvas.pop("show_labels", None)
    allow_inner = circular_canvas.pop("allow_inner_labels", None)
    if all(
        item is None
        for item in (shared_show, circular_show, linear_show, allow_inner)
    ):
        return migrated

    circular_scope: str | None = None
    linear_scope: str | None = None
    if shared_show is not None:
        circular_scope, linear_scope = _legacy_mode_label_scopes(
            shared_show,
            mode=mode,
            path=f"{path}.canvas.show_labels",
        )
    if circular_show is not None:
        nested_circular_scope, _unused_linear_scope = _legacy_mode_label_scopes(
            circular_show,
            mode="circular",
            path=f"{path}.canvas.circular.show_labels",
        )
        if (
            circular_scope is not None
            and circular_scope != nested_circular_scope
        ):
            raise CanonicalRequestDecodingError(
                f"{path}.canvas contains conflicting Circular label values."
            )
        circular_scope = nested_circular_scope
    if linear_show is not None:
        _unused_circular_scope, nested_linear_scope = _legacy_mode_label_scopes(
            linear_show,
            mode="linear",
            path=f"{path}.canvas.linear.show_labels",
        )
        if linear_scope is not None and linear_scope != nested_linear_scope:
            raise CanonicalRequestDecodingError(
                f"{path}.canvas contains conflicting Linear label values."
            )
        linear_scope = nested_linear_scope
    if allow_inner is not None:
        if not isinstance(allow_inner, bool):
            raise CanonicalRequestDecodingError(
                f"{path}.canvas.circular.allow_inner_labels must be a boolean."
            )
        if allow_inner and circular_scope in {None, "outer"}:
            circular_scope = "both"

    if not isinstance(labels, dict):
        raise CanonicalRequestDecodingError(f"{path}.labels must be an object.")
    for section_name, scope in (
        ("circular", circular_scope),
        ("linear", linear_scope),
    ):
        if scope is None:
            continue
        section = labels.get(section_name)
        if section is None:
            section = {}
            labels[section_name] = section
        elif not isinstance(section, dict):
            raise CanonicalRequestDecodingError(
                f"{path}.labels.{section_name} must be an object."
            )
        existing = section.get("scope")
        if existing is not None and existing != scope:
            raise CanonicalRequestDecodingError(
                f"{path} contains conflicting legacy and current "
                f"{section_name} label scopes."
            )
        section.setdefault("scope", scope)
    return migrated


def _encode_colors(value: object, *, resources: _ResourceBuilder) -> dict[str, Any]:
    if not isinstance(value, ColorOptions):
        raise CanonicalRequestEncodingError("diagramOptions.colors must be ColorOptions.")
    return {
        "colorTable": (
            _table_ref("colors-color-table", value.color_table, resources=resources)
            if value.color_table is not None
            else None
        ),
        "colorTableFile": (
            _file_ref("colors-color-table-file", value.color_table_file, resources=resources)
            if value.color_table_file is not None
            else None
        ),
        "defaultColors": (
            _table_ref("colors-default-colors", value.default_colors, resources=resources)
            if value.default_colors is not None
            else None
        ),
        "defaultColorsPalette": value.default_colors_palette,
        "defaultColorsFile": (
            _file_ref(
                "colors-default-colors-file", value.default_colors_file, resources=resources
            )
            if value.default_colors_file is not None
            else None
        ),
    }


def _decode_colors(
    value: object,
    *,
    resource_paths: Mapping[str, str | Path],
) -> ColorOptions:
    path = "renderRequest.diagramOptions.colors"
    payload = _object(
        value,
        path=path,
        required={
            "colorTable",
            "colorTableFile",
            "defaultColors",
            "defaultColorsPalette",
            "defaultColorsFile",
        },
    )
    result = ColorOptions(
        color_table=(
            _decode_table_ref(
                payload["colorTable"],
                name="colors-color-table",
                resource_paths=resource_paths,
            )
            if payload["colorTable"] is not None
            else None
        ),
        color_table_file=(
            str(
                _decode_file_ref(
                    payload["colorTableFile"],
                    name="colors-color-table-file",
                    resource_paths=resource_paths,
                )
            )
            if payload["colorTableFile"] is not None
            else None
        ),
        default_colors=(
            _decode_table_ref(
                payload["defaultColors"],
                name="colors-default-colors",
                resource_paths=resource_paths,
            )
            if payload["defaultColors"] is not None
            else None
        ),
        default_colors_palette=payload["defaultColorsPalette"],
        default_colors_file=(
            str(
                _decode_file_ref(
                    payload["defaultColorsFile"],
                    name="colors-default-colors-file",
                    resource_paths=resource_paths,
                )
            )
            if payload["defaultColorsFile"] is not None
            else None
        ),
    )
    _validate_dataclass_contract(result, path="diagramOptions.colors", error="decode")
    return result


def _encode_tracks(value: object) -> dict[str, Any]:
    if isinstance(value, CircularRequestTrackOptions):
        circular_slots = value.circular_track_slots
        circular_axis_index = value.circular_track_axis_index
        linear_slots = None
        linear_axis_index = None
        center_reserved_radius = value.center_reserved_radius
    elif isinstance(value, LinearRequestTrackOptions):
        circular_slots = None
        circular_axis_index = None
        linear_slots = value.linear_track_slots
        linear_axis_index = value.linear_track_axis_index
        center_reserved_radius = None
    else:
        raise CanonicalRequestEncodingError(
            "diagramOptions.tracks must use mode-specific track options."
        )
    return {
        "circularTrackSlots": _encode_track_slots(circular_slots),
        "circularTrackAxisIndex": circular_axis_index,
        "linearTrackSlots": _encode_track_slots(linear_slots),
        "linearTrackAxisIndex": linear_axis_index,
        "centerReservedRadius": center_reserved_radius,
    }


def _decode_tracks(
    value: object,
    *,
    mode: Literal["circular", "linear"],
    schema: int,
) -> CircularRequestTrackOptions | LinearRequestTrackOptions:
    path = "renderRequest.diagramOptions.tracks"
    payload = _object(
        value,
        path=path,
        required={
            "circularTrackSlots",
            "circularTrackAxisIndex",
            "linearTrackSlots",
            "linearTrackAxisIndex",
            "centerReservedRadius",
        },
    )
    circular_slots = _decode_track_slots(
        payload["circularTrackSlots"],
        mode="circular",
        schema=schema,
        path=f"{path}.circularTrackSlots",
    )
    linear_slots = _decode_track_slots(
        payload["linearTrackSlots"],
        mode="linear",
        schema=schema,
        path=f"{path}.linearTrackSlots",
    )
    if mode == "circular":
        if linear_slots is not None or payload["linearTrackAxisIndex"] is not None:
            raise CanonicalRequestDecodingError(
                "A Circular request cannot contain Linear track values."
            )
        result: CircularRequestTrackOptions | LinearRequestTrackOptions = (
            CircularRequestTrackOptions(
                circular_track_slots=circular_slots,
                circular_track_axis_index=payload["circularTrackAxisIndex"],
                center_reserved_radius=payload["centerReservedRadius"],
            )
        )
    else:
        if (
            circular_slots is not None
            or payload["circularTrackAxisIndex"] is not None
            or payload["centerReservedRadius"] is not None
        ):
            raise CanonicalRequestDecodingError(
                "A Linear request cannot contain Circular track values."
            )
        result = LinearRequestTrackOptions(
            linear_track_slots=linear_slots,
            linear_track_axis_index=payload["linearTrackAxisIndex"],
        )
    _validate_dataclass_contract(result, path="diagramOptions.tracks", error="decode")
    return result


def _encode_annotation_style(style: RegionAnnotationStyle) -> dict[str, Any]:
    return {
        "stroke": style.stroke,
        "strokeWidth": style.stroke_width,
        "strokeDasharray": list(style.stroke_dasharray),
        "lineCap": style.line_cap,
        "fill": style.fill,
        "fillOpacity": style.fill_opacity,
        "hatch": (
            {
                "angle": style.hatch.angle,
                "spacing": style.hatch.spacing,
                "color": style.hatch.color,
                "width": style.hatch.width,
                "cross": style.hatch.cross,
            }
            if style.hatch is not None
            else None
        ),
        "labelColor": style.label_color,
        "labelFontSize": style.label_font_size,
        "labelOrientation": style.label_orientation,
        "labelPosition": style.label_position,
        "labelOffset": style.label_offset,
    }


def _decode_annotation_style(value: object, *, path: str) -> RegionAnnotationStyle:
    payload = _object(
        value,
        path=path,
        required={
            "stroke", "strokeWidth", "strokeDasharray", "lineCap", "fill",
            "fillOpacity", "hatch", "labelColor", "labelFontSize",
            "labelOrientation", "labelPosition", "labelOffset",
        },
    )
    hatch_payload = payload["hatch"]
    hatch = None
    if hatch_payload is not None:
        raw_hatch = _object(
            hatch_payload,
            path=f"{path}.hatch",
            required={"angle", "spacing", "color", "width", "cross"},
        )
        hatch = HatchStyle(**raw_hatch)
    dasharray = _array(payload["strokeDasharray"], path=f"{path}.strokeDasharray")
    return RegionAnnotationStyle(
        stroke=payload["stroke"],
        stroke_width=payload["strokeWidth"],
        stroke_dasharray=tuple(dasharray),
        line_cap=payload["lineCap"],
        fill=payload["fill"],
        fill_opacity=payload["fillOpacity"],
        hatch=hatch,
        label_color=payload["labelColor"],
        label_font_size=payload["labelFontSize"],
        label_orientation=payload["labelOrientation"],
        label_position=payload["labelPosition"],
        label_offset=payload["labelOffset"],
    )


def _encode_annotations(value: object, *, resources: _ResourceBuilder) -> dict[str, Any]:
    if not isinstance(value, AnnotationOptions):
        raise CanonicalRequestEncodingError("diagramOptions.annotations must be AnnotationOptions.")
    sets = []
    for annotation_set in value.sets:
        items = []
        for annotation in annotation_set.annotations:
            target = annotation.target
            if isinstance(target, CoordinateSpan):
                target_payload = {
                    "kind": "coordinateSpan",
                    "record": _encode_selector(target.record),
                    "start": target.start,
                    "end": target.end,
                    "coordinateSpace": target.coordinate_space,
                    "wrapsOrigin": target.wraps_origin,
                    "outOfBounds": target.out_of_bounds,
                }
            elif isinstance(target, FeatureSpan):
                target_payload = {
                    "kind": "featureSpan",
                    "record": _encode_selector(target.record),
                    "selectors": [
                        {"key": selector.key, "value": selector.value}
                        for selector in target.selectors
                    ],
                    "envelope": target.envelope,
                    "circularPath": target.circular_path,
                }
            else:  # pragma: no cover - RegionAnnotation validates the union.
                raise CanonicalRequestEncodingError("Unsupported annotation target.")
            items.append(
                {
                    "id": annotation.id,
                    "target": target_payload,
                    "label": annotation.label,
                    "mark": annotation.mark,
                    "lane": annotation.lane,
                    "style": _encode_annotation_style(annotation.style) if annotation.style else None,
                    "legendLabel": annotation.legend_label,
                    "metadata": dict(annotation.metadata),
                }
            )
        sets.append(
            {
                "id": annotation_set.id,
                "annotations": items,
                "defaultStyle": _encode_annotation_style(annotation_set.default_style),
                "legendLabel": annotation_set.legend_label,
            }
        )
    return {
        "sets": sets,
        "table": (
            _table_ref("annotations-table", value.table, resources=resources)
            if value.table is not None else None
        ),
        "tableFile": (
            _file_ref("annotations-table-file", value.table_file, resources=resources)
            if value.table_file is not None else None
        ),
    }


def _decode_annotations(
    value: object,
    *,
    resource_paths: Mapping[str, str | Path],
) -> AnnotationOptions:
    path = "renderRequest.diagramOptions.annotations"
    payload = _object(value, path=path, required={"sets", "table", "tableFile"})
    raw_sets = _array(payload["sets"], path=f"{path}.sets")
    sets: list[AnnotationSet] = []
    for set_index, raw_set in enumerate(raw_sets):
        set_path = f"{path}.sets[{set_index}]"
        set_payload = _object(
            raw_set,
            path=set_path,
            required={"id", "annotations", "defaultStyle", "legendLabel"},
        )
        raw_items = _array(set_payload["annotations"], path=f"{set_path}.annotations")
        annotations: list[RegionAnnotation] = []
        for item_index, raw_item in enumerate(raw_items):
            item_path = f"{set_path}.annotations[{item_index}]"
            item = _object(
                raw_item,
                path=item_path,
                required={"id", "target", "label", "mark", "lane", "style", "legendLabel", "metadata"},
            )
            target_payload = _object(item["target"], path=f"{item_path}.target", required={"kind"}, exact=False)
            kind = target_payload["kind"]
            if kind == "coordinateSpan":
                _require_exact_fields(
                    target_payload,
                    path=f"{item_path}.target",
                    required={"kind", "record", "start", "end", "coordinateSpace", "wrapsOrigin", "outOfBounds"},
                )
                target = CoordinateSpan(
                    record=_decode_selector(target_payload["record"], path=f"{item_path}.target.record"),
                    start=target_payload["start"],
                    end=target_payload["end"],
                    coordinate_space=target_payload["coordinateSpace"],
                    wraps_origin=target_payload["wrapsOrigin"],
                    out_of_bounds=target_payload["outOfBounds"],
                )
            elif kind == "featureSpan":
                _require_exact_fields(
                    target_payload,
                    path=f"{item_path}.target",
                    required={"kind", "record", "selectors", "envelope", "circularPath"},
                )
                raw_selectors = _array(target_payload["selectors"], path=f"{item_path}.target.selectors")
                selectors = tuple(
                    FeatureSelector(**_object(raw, path=f"{item_path}.target.selectors[{index}]", required={"key", "value"}))
                    for index, raw in enumerate(raw_selectors)
                )
                target = FeatureSpan(
                    record=_decode_selector(target_payload["record"], path=f"{item_path}.target.record"),
                    selectors=selectors,
                    envelope=target_payload["envelope"],
                    circular_path=target_payload["circularPath"],
                )
            else:
                raise CanonicalRequestDecodingError(f"Unsupported annotation target kind at {item_path}: {kind!r}.")
            metadata = _object(item["metadata"], path=f"{item_path}.metadata")
            annotations.append(
                RegionAnnotation(
                    id=item["id"], target=target, label=item["label"], mark=item["mark"],
                    lane=item["lane"],
                    style=(_decode_annotation_style(item["style"], path=f"{item_path}.style") if item["style"] is not None else None),
                    legend_label=item["legendLabel"], metadata=metadata,
                )
            )
        sets.append(
            AnnotationSet(
                id=set_payload["id"],
                annotations=tuple(annotations),
                default_style=_decode_annotation_style(set_payload["defaultStyle"], path=f"{set_path}.defaultStyle"),
                legend_label=set_payload["legendLabel"],
            )
        )
    table = (
        _decode_table_ref(payload["table"], name="annotations-table", resource_paths=resource_paths)
        if payload["table"] is not None else None
    )
    table_file = (
        str(_decode_file_ref(payload["tableFile"], name="annotations-table-file", resource_paths=resource_paths))
        if payload["tableFile"] is not None else None
    )
    return AnnotationOptions(sets=tuple(sets), table=table, table_file=table_file)


def _encode_track_slots(value: object) -> list[Any] | None:
    if value is None:
        return None
    result: list[Any] = []
    for slot in _sequence(value, name="track slots"):
        if isinstance(slot, str):
            result.append(slot)
        elif isinstance(slot, CircularTrackSlot):
            result.append(_encode_track_slot(slot, kind="circularTrackSlot"))
        elif isinstance(slot, LinearTrackSlot):
            result.append(_encode_track_slot(slot, kind="linearTrackSlot"))
        else:
            raise CanonicalRequestEncodingError(
                f"Unsupported track slot value: {type(slot).__name__}."
            )
    return result


def _encode_track_slot(
    slot: CircularTrackSlot | LinearTrackSlot,
    *,
    kind: str,
) -> dict[str, Any]:
    result = {"kind": kind}
    legacy_spacing = (
        slot.legacy_spacing
        if isinstance(slot, _InternalCircularTrackSlot)
        else None
    )
    slot_fields = (
        fields(CircularTrackSlot)
        if isinstance(slot, CircularTrackSlot)
        else fields(LinearTrackSlot)
    )
    for item in slot_fields:
        raw = getattr(slot, item.name)
        if isinstance(raw, ScalarSpec):
            raw = {"value": raw.value, "unit": raw.unit}
        elif item.name == "params" and isinstance(raw, MappingABC):
            raw = dict(raw)
            private_keys = [
                str(key)
                for key in raw
                if str(key).strip().startswith("_")
            ]
            if private_keys:
                raise CanonicalRequestEncodingError(
                    "trackSlot.params cannot contain private key(s): "
                    + ", ".join(sorted(private_keys))
                    + "."
                )
            if isinstance(raw.get("style_override"), RegionAnnotationStyle):
                raw["style_override"] = _encode_annotation_style(raw["style_override"])
        result[_camel(item.name)] = _json_value(raw, path=f"trackSlot.{_camel(item.name)}")
    if legacy_spacing is not None:
        if legacy_spacing.unit != "px":
            raise CanonicalRequestEncodingError(
                "A legacy factor-based Circular spacing value can be replayed "
                "but cannot be written to the current canonical schema. Set explicit "
                "inner_gap_px and outer_gap_px values before saving."
            )
        result["innerGapPx"] = (
            slot.inner_gap_px
            if slot.inner_gap_px is not None
            else legacy_spacing.value
        )
        result["outerGapPx"] = (
            slot.outer_gap_px
            if slot.outer_gap_px is not None
            else legacy_spacing.value
        )
    return result


def _decode_track_slots(
    value: object,
    *,
    mode: Literal["circular", "linear"],
    schema: int,
    path: str,
) -> tuple[str | CircularTrackSlot | LinearTrackSlot, ...] | None:
    if value is None:
        return None
    raw_slots = _array(value, path=path)
    expected_kind = "circularTrackSlot" if mode == "circular" else "linearTrackSlot"
    cls = CircularTrackSlot if mode == "circular" else LinearTrackSlot
    result: list[str | CircularTrackSlot | LinearTrackSlot] = []
    for index, raw in enumerate(raw_slots):
        if isinstance(raw, str):
            if mode == "circular" and schema in {1, 2}:
                migrated = _migrate_legacy_circular_track_slot_spec(raw)
                try:
                    result.append(
                        parse_circular_track_slot(
                            migrated,
                            _allow_legacy_transport=True,
                        )
                    )
                except (TypeError, ValueError) as exc:
                    raise CanonicalRequestDecodingError(
                        f"Invalid legacy Circular slot at {path}[{index}]: {exc}"
                    ) from exc
            else:
                result.append(raw)
            continue
        slot_path = f"{path}[{index}]"
        slot = _object(raw, path=slot_path, required={"kind"}, exact=False)
        legacy_spacing: ScalarSpec | None = None
        if mode == "circular" and schema in {1, 2} and "spacing" in slot:
            slot = dict(slot)
            raw_legacy_spacing = slot.pop("spacing")
            if raw_legacy_spacing is not None:
                decoded_spacing = _decode_scalar_if_needed(
                    raw_legacy_spacing,
                    path=f"{slot_path}.spacing",
                )
                legacy_spacing = (
                    decoded_spacing
                    if isinstance(decoded_spacing, ScalarSpec)
                    else ScalarSpec.parse(decoded_spacing)
                )
        if slot["kind"] != expected_kind:
            raise CanonicalRequestDecodingError(
                f"Unsupported track slot kind at {slot_path}: {slot['kind']!r}."
            )
        field_map = {_camel(item.name): item.name for item in fields(cls)}
        _require_exact_fields(
            slot, path=slot_path, required={"kind", *field_map.keys()}
        )
        scalar_fields = (
            {"radius", "width"}
            if mode == "circular"
            else {"height", "spacing"}
        )
        kwargs = {
            name: (
                _decode_scalar_if_needed(slot[key], path=f"{slot_path}.{key}")
                if name in scalar_fields
                else slot[key]
            )
            for key, name in field_map.items()
        }
        if mode == "circular" and schema in {1, 2}:
            params = dict(kwargs.get("params") or {})
            private_spacing = params.pop("__gbdraw_legacy_spacing", None)
            if private_spacing is not None:
                if legacy_spacing is not None:
                    raise CanonicalRequestDecodingError(
                        f"{slot_path} contains two legacy spacing values."
                    )
                legacy_spacing = (
                    private_spacing
                    if isinstance(private_spacing, ScalarSpec)
                    else ScalarSpec.parse(private_spacing)
                )
            kwargs["params"] = params
        if str(kwargs.get("renderer", "")).strip().lower() == "annotations":
            params = dict(kwargs.get("params") or {})
            if isinstance(params.get("style_override"), MappingABC):
                params["style_override"] = _decode_annotation_style(
                    params["style_override"], path=f"{slot_path}.params.style_override"
                )
            kwargs["params"] = params
        decoded_slot = cls(**kwargs)
        if (
            isinstance(decoded_slot, CircularTrackSlot)
            and legacy_spacing is not None
        ):
            decoded_slot = _internal_circular_track_slot(
                decoded_slot,
                legacy_spacing=legacy_spacing,
            )
        result.append(decoded_slot)
    return tuple(result)


def _migrate_legacy_circular_track_slot_spec(spec: str) -> str:
    """Translate retired Circular slot fields at the persisted-data boundary."""

    head, separator, options = str(spec).partition("@")
    if not separator:
        return str(spec)
    migrated: list[str] = []
    for raw_part in options.split(","):
        part = raw_part.strip()
        if not part or "=" not in part:
            migrated.append(part)
            continue
        raw_key, raw_value = part.split("=", 1)
        key = raw_key.strip().lower()
        if key in {"strict", "compress", "reserve"}:
            continue
        if key == "spacing":
            migrated.append(
                f"__gbdraw_legacy_spacing={raw_value.strip()}"
            )
        else:
            migrated.append(part)
    return head if not migrated else f"{head}@{','.join(migrated)}"


def _decode_scalar_if_needed(value: object, *, path: str) -> object:
    if not isinstance(value, MappingABC) or "unit" not in value:
        return value
    scalar = _object(value, path=path, required={"value", "unit"})
    result = ScalarSpec(value=scalar["value"], unit=scalar["unit"])
    _validate_dataclass_contract(result, path=path, error="decode")
    return result


def _encode_assembly_output(value: object) -> dict[str, Any]:
    if not isinstance(value, (CircularOutputOptions, LinearOutputOptions)):
        raise CanonicalRequestEncodingError(
            "diagramOptions.output must use mode-specific output options."
        )
    return {
        "legend": value.legend,
        "plotTitlePosition": value.plot_title_position,
    }


def _decode_assembly_output(
    value: object,
    *,
    mode: Literal["circular", "linear"],
    schema: int,
) -> CircularOutputOptions | LinearOutputOptions:
    required = {"legend", "plotTitlePosition"}
    if schema in {1, 2}:
        required.add("outputPrefix")
    payload = _object(
        value,
        path="renderRequest.diagramOptions.output",
        required=required,
    )
    output_type = (
        CircularOutputOptions if mode == "circular" else LinearOutputOptions
    )
    result = output_type(
        legend=payload["legend"],
        plot_title_position=payload["plotTitlePosition"],
    )
    _validate_dataclass_contract(result, path="diagramOptions.output", error="decode")
    return result


def _encode_resource_matrix(
    name: str,
    value: object,
    *,
    table: bool,
    resources: _ResourceBuilder,
) -> list[list[dict[str, str] | None]]:
    result: list[list[dict[str, str] | None]] = []
    for row_index, row in enumerate(_sequence(value, name=name), start=1):
        encoded_row: list[dict[str, str] | None] = []
        for column_index, item in enumerate(_sequence(row, name=name), start=1):
            if item is None:
                encoded_row.append(None)
                continue
            resource_name = f"{name}-{row_index}-{column_index}"
            encoded_row.append(
                _table_ref(resource_name, item, resources=resources)
                if table
                else _file_ref(resource_name, item, resources=resources)
            )
        result.append(encoded_row)
    return result


def _decode_resource_matrix(
    value: object,
    *,
    name: str,
    table: bool,
    resource_paths: Mapping[str, str | Path],
) -> tuple[tuple[DataFrame | str | None, ...], ...]:
    rows = _array(value, path=f"renderRequest.diagramOptions.{_camel(name)}")
    result: list[tuple[DataFrame | str | None, ...]] = []
    for row_index, raw_row in enumerate(rows, start=1):
        row = _array(
            raw_row,
            path=f"renderRequest.diagramOptions.{_camel(name)}[{row_index - 1}]",
        )
        decoded_row: list[DataFrame | str | None] = []
        for column_index, item in enumerate(row, start=1):
            if item is None:
                decoded_row.append(None)
                continue
            resource_name = f"{name}-{row_index}-{column_index}"
            if table:
                decoded_row.append(
                    _decode_table_ref(
                        item, name=resource_name, resource_paths=resource_paths
                    )
                )
            else:
                decoded_row.append(
                    str(
                        _decode_file_ref(
                            item, name=resource_name, resource_paths=resource_paths
                        )
                    )
                )
        result.append(tuple(decoded_row))
    return tuple(result)


def _encode_comparisons(
    options: CircularDiagramOptions | LinearDiagramOptions,
    *,
    mode: Literal["circular", "linear"],
    resources: _ResourceBuilder,
) -> list[dict[str, Any]]:
    if mode == "circular":
        return []
    if not isinstance(options, LinearDiagramOptions):
        raise CanonicalRequestEncodingError(
            "Linear comparisons require LinearDiagramOptions."
        )
    result: list[dict[str, Any]] = []
    for index, comparison in enumerate(options.linear_comparisons or (), start=1):
        ref = _table_ref(
            f"comparison-explicit-{index}", comparison.matches, resources=resources
        )
        result.append(
            {
                "kind": "precomputedProteinComparison",
                "resourceId": ref["resourceId"],
                "encoding": "canonicalTsv",
                "queryRecordIndex": comparison.query_record_index,
                "subjectRecordIndex": comparison.subject_record_index,
            }
        )
    for index, path in enumerate(options.blast_files or (), start=1):
        resource_id = resources.add_path(
            f"comparison-nucleotide-{index}", kind="nucleotide-blast", value=path
        )
        result.append(
            {
                "kind": "nucleotideBlast",
                "resourceId": resource_id,
                "queryRecordIndex": index - 1,
                "subjectRecordIndex": index,
            }
        )
    for index, table in enumerate(options.protein_comparisons or (), start=1):
        ref = _table_ref(f"comparison-protein-{index}", table, resources=resources)
        result.append(
            {
                "kind": "precomputedProteinComparison",
                "resourceId": ref["resourceId"],
                "encoding": "canonicalTsv",
                "queryRecordIndex": index - 1,
                "subjectRecordIndex": index,
            }
        )
    if options.orthogroups is not None:
        resource_id = _typed_json_resource(
            "comparison-orthogroups",
            kind="orthogroup-result",
            value_kind="orthogroupResult",
            value=options.orthogroups,
            resources=resources,
        )
        result.append(
            {
                "kind": "orthogroupResult",
                "resourceId": resource_id,
                "encoding": "canonicalJson",
            }
        )
    if options.collinearity_blocks is not None:
        value_kind = (
            "result"
            if isinstance(options.collinearity_blocks, CollinearityResult)
            else "blocks"
        )
        resource_id = _typed_json_resource(
            "comparison-collinearity",
            kind="collinearity-result",
            value_kind=value_kind,
            value=options.collinearity_blocks,
            resources=resources,
        )
        result.append(
            {
                "kind": "collinearityResult",
                "resourceId": resource_id,
                "encoding": "canonicalJson",
                "valueKind": value_kind,
            }
        )
    if any(
        not _same_default(
            getattr(options, name),
            getattr(_DEFAULT_LINEAR_OPTIONS, name),
        )
        for name in _PIPELINE_FIELDS
    ) or options.protein_comparison_pairs is not None:
        settings = {
            _camel(name): _encode_pipeline_value(name, getattr(options, name))
            for name in _PIPELINE_FIELDS
            if name != "protein_blastp_mode"
        }
        result.append(
            {
                "kind": "generatedProteinComparison",
                "mode": options.protein_blastp_mode,
                "pairs": [
                    {
                        "queryRecordIndex": int(pair[0]),
                        "subjectRecordIndex": int(pair[1]),
                    }
                    for pair in (options.protein_comparison_pairs or ())
                ],
                "settings": settings,
            }
        )
    return result


def _decode_comparisons(
    value: object,
    *,
    mode: Literal["circular", "linear"],
    schema: int,
    resource_paths: Mapping[str, str | Path],
) -> dict[str, Any]:
    comparisons = _array(value, path="renderRequest.comparisons")
    if mode == "circular" and comparisons:
        raise CanonicalRequestDecodingError(
            "A circular canonical request cannot contain linear comparisons."
        )
    blast_files: list[str] = []
    protein_tables: list[DataFrame] = []
    explicit_comparisons: list[LinearComparison] = []
    result: dict[str, Any] = {}
    singleton_kinds: set[str] = set()
    for index, raw in enumerate(comparisons):
        path = f"renderRequest.comparisons[{index}]"
        item = _object(raw, path=path, required={"kind"}, exact=False)
        kind = item["kind"]
        if kind == "nucleotideBlast":
            required = {"kind", "resourceId"}
            if schema >= 2:
                required |= {"queryRecordIndex", "subjectRecordIndex"}
            _require_exact_fields(item, path=path, required=required)
            resource_path = _resolve_resource(
                item["resourceId"],
                path=f"{path}.resourceId",
                resource_paths=resource_paths,
            )
            if schema >= 2:
                try:
                    table = read_csv(
                        resource_path,
                        sep="\t",
                        comment="#",
                        names=COMPARISON_COLUMNS,
                    )
                except Exception as exc:
                    raise CanonicalRequestDecodingError(
                        f"Could not decode BLAST resource for {path}."
                    ) from exc
                query_index = _non_negative_index(
                    item["queryRecordIndex"], f"{path}.queryRecordIndex"
                )
                subject_index = _non_negative_index(
                    item["subjectRecordIndex"], f"{path}.subjectRecordIndex"
                )
                if query_index == len(blast_files) and subject_index == query_index + 1:
                    blast_files.append(str(resource_path))
                else:
                    explicit_comparisons.append(
                        LinearComparison(query_index, subject_index, table)
                    )
            else:
                blast_files.append(str(resource_path))
        elif kind == "precomputedProteinComparison":
            required = {"kind", "resourceId", "encoding"}
            if schema >= 2:
                required |= {"queryRecordIndex", "subjectRecordIndex"}
            _require_exact_fields(item, path=path, required=required)
            if item["encoding"] != "canonicalTsv":
                raise CanonicalRequestDecodingError(
                    f"Unsupported protein comparison encoding at {path}."
                )
            table = _read_canonical_table(
                _resolve_resource(
                    item["resourceId"],
                    path=f"{path}.resourceId",
                    resource_paths=resource_paths,
                ),
                context=path,
            )
            if schema >= 2:
                query_index = _non_negative_index(
                    item["queryRecordIndex"], f"{path}.queryRecordIndex"
                )
                subject_index = _non_negative_index(
                    item["subjectRecordIndex"], f"{path}.subjectRecordIndex"
                )
                if query_index == len(protein_tables) and subject_index == query_index + 1:
                    protein_tables.append(table)
                else:
                    explicit_comparisons.append(
                        LinearComparison(query_index, subject_index, table)
                    )
            else:
                protein_tables.append(table)
        elif kind == "orthogroupResult":
            _check_singleton_kind(kind, singleton_kinds, path=path)
            _require_exact_fields(
                item, path=path, required={"kind", "resourceId", "encoding"}
            )
            if item["encoding"] != "canonicalJson":
                raise CanonicalRequestDecodingError(
                    f"Unsupported orthogroup encoding at {path}."
                )
            result["orthogroups"] = _read_typed_json_resource(
                item["resourceId"],
                value_kind="orthogroupResult",
                expected=OrthogroupResult,
                path=path,
                resource_paths=resource_paths,
            )
        elif kind == "collinearityResult":
            _check_singleton_kind(kind, singleton_kinds, path=path)
            _require_exact_fields(
                item,
                path=path,
                required={"kind", "resourceId", "encoding", "valueKind"},
            )
            if item["encoding"] != "canonicalJson" or item["valueKind"] not in {
                "result",
                "blocks",
            }:
                raise CanonicalRequestDecodingError(
                    f"Unsupported collinearity representation at {path}."
                )
            expected = CollinearityResult if item["valueKind"] == "result" else tuple[CollinearityBlock, ...]
            result["collinearity_blocks"] = _read_typed_json_resource(
                item["resourceId"],
                value_kind=item["valueKind"],
                expected=expected,
                path=path,
                resource_paths=resource_paths,
            )
        elif kind == "generatedProteinComparison":
            _check_singleton_kind(kind, singleton_kinds, path=path)
            result.update(_decode_pipeline(item, path=path, schema=schema))
        else:
            raise CanonicalRequestDecodingError(
                f"Unsupported comparison kind at {path}: {kind!r}."
            )
    if blast_files:
        result["blast_files"] = tuple(blast_files)
    if protein_tables:
        result["protein_comparisons"] = tuple(protein_tables)
    if explicit_comparisons:
        result["linear_comparisons"] = tuple(explicit_comparisons)
    return result


def _encode_pipeline_value(name: str, value: object) -> Any:
    if name == "collinearity_params":
        if value is None:
            return None
        if not isinstance(value, LosslessCollinearityParameters):
            raise CanonicalRequestEncodingError(
                "Unsupported collinearity parameter object."
            )
        value.validate()
        return {
            "kind": "lossless",
            "parameters": {
                _camel(item.name): _json_value(
                    getattr(value, item.name), path=f"collinearityParams.{item.name}"
                )
                for item in fields(value)
            },
        }
    return _json_value(value, path=f"comparisons.settings.{_camel(name)}")


def _decode_pipeline(
    item: Mapping[str, Any], *, path: str, schema: int
) -> dict[str, Any]:
    required = {"kind", "mode", "settings"}
    if schema >= 2:
        required.add("pairs")
    _require_exact_fields(item, path=path, required=required)
    settings = _object(item["settings"], path=f"{path}.settings")
    setting_fields = tuple(name for name in _PIPELINE_FIELDS if name != "protein_blastp_mode")
    field_map = {_camel(name): name for name in setting_fields}
    _require_exact_fields(
        settings,
        path=f"{path}.settings",
        required=set() if schema == CANONICAL_REQUEST_SCHEMA else set(field_map),
        optional=set(field_map) if schema == CANONICAL_REQUEST_SCHEMA else set(),
    )
    mode = item["mode"]
    result = {"protein_blastp_mode": mode}
    legacy_max_paralog_links: int | None = None
    if schema >= 2:
        pairs = _array(item["pairs"], path=f"{path}.pairs")
        decoded_pairs: list[tuple[int, int]] = []
        for index, raw_pair in enumerate(pairs):
            pair_path = f"{path}.pairs[{index}]"
            pair = _object(
                raw_pair,
                path=pair_path,
                required={"queryRecordIndex", "subjectRecordIndex"},
            )
            decoded_pairs.append(
                (
                    _non_negative_index(
                        pair["queryRecordIndex"], f"{pair_path}.queryRecordIndex"
                    ),
                    _non_negative_index(
                        pair["subjectRecordIndex"], f"{pair_path}.subjectRecordIndex"
                    ),
                )
            )
        # Early schema-2 writers also stored the derived row-adjacent search pairs
        # used by collinear mode here.  The public option is pairwise-only; current
        # collinear rendering derives those pairs from the saved layout instead.
        result["protein_comparison_pairs"] = (
            tuple(decoded_pairs) if decoded_pairs and mode == "pairwise" else None
        )
    for key, name in field_map.items():
        if key not in settings:
            continue
        raw = settings[key]
        if name == "collinearity_params":
            decoded, legacy_max_paralog_links = _decode_collinearity_params(
                raw,
                path=f"{path}.settings.{key}",
                allow_standard=schema in {1, 2},
            )
            result[name] = decoded
        else:
            result[name] = raw
    if (
        mode == "collinear"
        and legacy_max_paralog_links is not None
        and result.get("collinear_max_paralog_links_per_orthogroup", 2) == 2
    ):
        result["collinear_max_paralog_links_per_orthogroup"] = (
            legacy_max_paralog_links
        )
    return result


def decode_canonical_protein_comparison_intent(
    comparison: Mapping[str, Any],
    *,
    diagram_options: Mapping[str, Any] | None = None,
) -> LinearDiagramOptions:
    """Decode the canonical generated-protein fields used by a Worker helper."""

    canonical_comparison = _object(
        comparison,
        path="generatedProteinComparison",
    )
    mode = canonical_comparison.get("mode")
    common_fields = {"proteinBlastpCandidateLimit"}
    mode_fields = {
        "pairwise": {"proteinBlastpMaxHits"},
        "orthogroup": {
            "orthogroupMembershipMode",
            "orthogroupMemberMaxHits",
            "collinearMaxParalogLinksPerOrthogroup",
        },
        "collinear": {
            "orthogroupMembershipMode",
            "orthogroupMemberMaxHits",
            "collinearMaxParalogLinksPerOrthogroup",
            "collinearityParams",
            "collinearityUnitMode",
            "collinearityAnchorMode",
            "collinearitySearchScope",
            "collinearityColorMode",
        },
    }
    active_fields = common_fields | mode_fields.get(str(mode), set())
    raw_settings = _object(
        canonical_comparison.get("settings"),
        path="generatedProteinComparison.settings",
    )
    active_comparison = {
        **canonical_comparison,
        "settings": {
            key: value
            for key, value in raw_settings.items()
            if key in active_fields
        },
    }
    decoded = _decode_pipeline(
        active_comparison,
        path="generatedProteinComparison",
        schema=CANONICAL_REQUEST_SCHEMA,
    )
    canonical_options = (
        _object(diagram_options, path="diagramOptions")
        if diagram_options is not None
        else {}
    )
    option_fields = {
        "evalue": "evalue",
        "bitscore": "bitscore",
        "identity": "identity",
        "alignmentLength": "alignment_length",
    }
    _require_exact_fields(
        canonical_options,
        path="diagramOptions",
        required=set(),
        optional=set(option_fields),
    )
    decoded.update(
        {
            name: canonical_options[key]
            for key, name in option_fields.items()
            if key in canonical_options
        }
    )
    return resolve_linear_diagram_options(LinearDiagramOptions(**decoded))


def _decode_collinearity_params(
    value: object,
    *,
    path: str,
    allow_standard: bool,
) -> tuple[LosslessCollinearityParameters | None, int | None]:
    if value is None:
        return None, None
    required = {"kind", "parameters"} if allow_standard else {"parameters"}
    payload = _object(value, path=path, required=required, exact=allow_standard)
    if not allow_standard:
        _require_exact_fields(payload, path=path, required=required, optional={"kind"})
    kind = payload.get("kind", "lossless")
    supported_kinds = {"standard", "lossless"} if allow_standard else {"lossless"}
    if kind not in supported_kinds:
        raise CanonicalRequestDecodingError(
            f"Unsupported collinearity parameter kind at {path}: {kind!r}."
        )
    cls = (
        _LegacyStandardCollinearityPayload
        if kind == "standard"
        else LosslessCollinearityParameters
    )
    parameters = _object(payload["parameters"], path=f"{path}.parameters")
    field_map = {_camel(item.name): item.name for item in fields(cls)}
    allow_omission = kind == "lossless" and not allow_standard
    _require_exact_fields(
        parameters,
        path=f"{path}.parameters",
        required=set() if allow_omission else set(field_map),
        optional=set(field_map) if allow_omission else set(),
    )
    result = cls(
        **{
            name: parameters[key]
            for key, name in field_map.items()
            if key in parameters
        }
    )
    _validate_dataclass_contract(result, path=path, error="decode")
    result.validate()
    if isinstance(result, _LegacyStandardCollinearityPayload):
        lossless = result.to_lossless()
        lossless.validate()
        return lossless, result.max_paralog_links_per_orthogroup
    return result, None


def _decode_output(
    value: object,
    *,
    output_directory: str | Path,
) -> RenderOutputRequest:
    payload = _object(
        value,
        path="renderRequest.output",
        required={"prefix", "formats", "overwrite", "interactiveMetadataPolicy"},
    )
    formats = payload["formats"]
    if not isinstance(formats, list):
        raise CanonicalRequestDecodingError(
            "renderRequest.output.formats must be an array."
        )
    _boolean(payload["overwrite"], "renderRequest.output.overwrite")
    return RenderOutputRequest(
        output_prefix=payload["prefix"],
        output_directory=output_directory,
        formats=tuple(formats),
        overwrite=False,
        interactive_metadata_policy=payload["interactiveMetadataPolicy"],
    )


def _decode_batch_outputs(
    value: object,
    *,
    output_directory: str | Path,
    record_count: int,
) -> tuple[RenderOutputRequest, ...]:
    payloads = _array(value, path="renderRequest.output")
    if len(payloads) != record_count:
        raise CanonicalRequestDecodingError(
            "A Circular batch canonical request requires one output per record."
        )
    return tuple(
        _decode_output(
            payload,
            output_directory=output_directory,
        )
        for payload in payloads
    )


def _table_ref(
    name: str,
    value: object,
    *,
    resources: _ResourceBuilder,
) -> dict[str, str]:
    if not isinstance(value, DataFrame):
        raise CanonicalRequestEncodingError(f"{name} must be a pandas DataFrame.")
    if len(value.columns) == 0:
        raise CanonicalRequestEncodingError(
            f"DataFrame resource {name!r} must have at least one column."
        )
    try:
        text = value.to_csv(sep="\t", index=False, lineterminator="\n")
    except Exception as exc:
        raise CanonicalRequestEncodingError(
            f"Could not encode DataFrame resource {name!r} as canonical TSV."
        ) from exc
    resource_id = _resource_id(name)
    resources.add_bytes(
        resource_id,
        kind="canonical-tsv",
        name=f"{resource_id}.tsv",
        content=text.encode("utf-8"),
    )
    return {"resourceId": resource_id, "representation": "canonicalTsv"}


def _file_ref(
    name: str,
    value: object,
    *,
    resources: _ResourceBuilder,
) -> dict[str, str]:
    resource_id = _resource_id(name)
    resources.add_path(resource_id, kind=_resource_id(name), value=value)
    return {"resourceId": resource_id, "representation": "file"}


def _decode_table_ref(
    value: object,
    *,
    name: str,
    resource_paths: Mapping[str, str | Path],
) -> DataFrame:
    ref = _resource_ref(value, path=f"resource reference {name}", representation="canonicalTsv")
    return _read_canonical_table(
        _resolve_resource(
            ref["resourceId"],
            path=f"resource reference {name}.resourceId",
            resource_paths=resource_paths,
        ),
        context=name,
    )


def _decode_file_ref(
    value: object,
    *,
    name: str,
    resource_paths: Mapping[str, str | Path],
) -> Path:
    ref = _resource_ref(value, path=f"resource reference {name}", representation="file")
    return _resolve_resource(
        ref["resourceId"],
        path=f"resource reference {name}.resourceId",
        resource_paths=resource_paths,
    )


def _resource_ref(
    value: object,
    *,
    path: str,
    representation: str,
) -> Mapping[str, Any]:
    ref = _object(
        value, path=path, required={"resourceId", "representation"}
    )
    if ref["representation"] != representation:
        raise CanonicalRequestDecodingError(
            f"Unsupported resource representation at {path}: {ref['representation']!r}."
        )
    return ref


def _read_canonical_table(file_path: Path, *, context: str) -> DataFrame:
    try:
        table = read_csv(file_path, sep="\t")
    except Exception as exc:
        raise CanonicalRequestDecodingError(
            f"Could not decode canonical TSV resource for {context}."
        ) from exc
    from gbdraw.api.prepared import register_prepared_resource_value

    return register_prepared_resource_value(table, file_path)


def _typed_json_resource(
    resource_id: str,
    *,
    kind: str,
    value_kind: str,
    value: object,
    resources: _ResourceBuilder,
) -> str:
    content = encode_canonical_typed_resource(value_kind, value)
    resources.add_bytes(
        resource_id,
        kind=kind,
        name=f"{resource_id}.json",
        content=content,
    )
    return resource_id


def encode_canonical_typed_resource(value_kind: str, value: object) -> bytes:
    """Encode one value using the canonical request's typed JSON resource contract."""

    if not isinstance(value_kind, str) or not value_kind.strip():
        raise CanonicalRequestEncodingError(
            "Canonical typed resources require a non-empty value kind."
        )
    body = {
        "schema": 2,
        "kind": value_kind,
        "value": _encode_typed_tree(value),
    }
    try:
        return json.dumps(
            body,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise CanonicalRequestEncodingError(
            f"Could not encode canonical typed resource {value_kind!r}."
        ) from exc


def _read_typed_json_resource(
    resource_id: object,
    *,
    value_kind: str,
    expected: object,
    path: str,
    resource_paths: Mapping[str, str | Path],
) -> Any:
    resource_path = _resolve_resource(
        resource_id, path=f"{path}.resourceId", resource_paths=resource_paths
    )

    def decode() -> Any:
        try:
            raw = json.loads(resource_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise CanonicalRequestDecodingError(
                f"Could not decode canonical JSON resource for {path}."
            ) from exc
        payload = _object(
            raw, path=f"{path} resource", required={"schema", "kind", "value"}
        )
        resource_schema = payload["schema"]
        if resource_schema not in {1, 2} or payload["kind"] != value_kind:
            raise CanonicalRequestDecodingError(
                f"Canonical JSON resource metadata does not match {path}."
            )
        return _decode_typed_tree(
            payload["value"],
            expected,
            path=f"{path}.value",
            resource_schema=resource_schema,
        )

    from gbdraw.api.prepared import get_or_build_decoded_resource

    return get_or_build_decoded_resource(
        resource_path,
        ("canonical-typed-json", value_kind, repr(expected)),
        decode,
    )


def _encode_typed_tree(value: object) -> Any:
    if is_dataclass(value) and not isinstance(value, type):
        cls = type(value)
        if cls.__name__ not in _TYPED_TREE_CLASSES:
            raise CanonicalRequestEncodingError(
                f"Unsupported typed resource value: {cls.__name__}."
            )
        return {
            "type": cls.__name__,
            "fields": {
                _camel(item.name): _encode_typed_tree(getattr(value, item.name))
                for item in fields(value)
            },
        }
    if isinstance(value, MappingABC):
        if not all(isinstance(key, str) for key in value):
            raise CanonicalRequestEncodingError(
                "Canonical typed resource mappings require string keys."
            )
        return {key: _encode_typed_tree(item) for key, item in value.items()}
    if isinstance(value, SequenceABC) and not isinstance(value, (str, bytes, bytearray)):
        return [_encode_typed_tree(item) for item in value]
    return _json_value(value, path="typed resource")


def _decode_typed_tree(
    value: object,
    hint: object,
    *,
    path: str,
    resource_schema: int,
) -> Any:
    if hint is Any or hint is object:
        return value
    origin = get_origin(hint)
    args = get_args(hint)
    if origin in {types.UnionType, Union}:
        if value is None and type(None) in args:
            return None
        tagged_type = value.get("type") if isinstance(value, MappingABC) else None
        for branch in args:
            if branch is type(None):
                continue
            if isinstance(branch, type) and branch.__name__ == tagged_type:
                return _decode_typed_tree(
                    value, branch, path=path, resource_schema=resource_schema
                )
            branch_origin = get_origin(branch)
            if branch_origin in {tuple, list, Sequence, SequenceABC} and isinstance(value, list):
                return _decode_typed_tree(
                    value, branch, path=path, resource_schema=resource_schema
                )
            if _matches_type(value, branch):
                return _decode_typed_tree(
                    value, branch, path=path, resource_schema=resource_schema
                )
        raise CanonicalRequestDecodingError(
            f"Typed resource value at {path} does not match its union contract."
        )
    if isinstance(hint, type) and is_dataclass(hint):
        tagged = _object(value, path=path, required={"type", "fields"})
        if tagged["type"] != hint.__name__ or hint.__name__ not in _TYPED_TREE_CLASSES:
            raise CanonicalRequestDecodingError(
                f"Unsupported typed resource object at {path}: {tagged['type']!r}."
            )
        raw_fields = _object(tagged["fields"], path=f"{path}.fields")
        hints = get_type_hints(hint)
        field_defs = {_camel(item.name): item for item in fields(hint)}
        optional_fields = (
            _SCHEMA1_TYPED_OPTIONAL_FIELDS.get(hint.__name__, frozenset())
            if resource_schema == 1
            else frozenset()
        )
        _require_exact_fields(
            raw_fields,
            path=f"{path}.fields",
            required=set(field_defs).difference(optional_fields),
            optional=set(optional_fields),
        )
        kwargs: dict[str, object] = {}
        for key, item in field_defs.items():
            if key in raw_fields:
                kwargs[item.name] = _decode_typed_tree(
                    raw_fields[key],
                    hints[item.name],
                    path=f"{path}.fields.{key}",
                    resource_schema=resource_schema,
                )
            elif item.default is not MISSING:
                kwargs[item.name] = item.default
            elif item.default_factory is not MISSING:
                kwargs[item.name] = item.default_factory()
            else:  # pragma: no cover - guarded by the required-field check
                raise CanonicalRequestDecodingError(
                    f"Missing required field at {path}.fields.{key}."
                )
        return hint(**kwargs)
    if origin in {tuple, list, Sequence, SequenceABC}:
        raw = _array(value, path=path)
        item_hint = args[0] if args else Any
        decoded = [
            _decode_typed_tree(
                item,
                item_hint,
                path=f"{path}[{index}]",
                resource_schema=resource_schema,
            )
            for index, item in enumerate(raw)
        ]
        return tuple(decoded) if origin is tuple else decoded
    if origin in {dict, Mapping, MappingABC}:
        raw = _object(value, path=path)
        key_hint, value_hint = args if len(args) == 2 else (str, Any)
        if key_hint is not str:
            raise CanonicalRequestDecodingError(
                f"Unsupported typed resource mapping key contract at {path}."
            )
        return {
            key: _decode_typed_tree(
                item,
                value_hint,
                path=f"{path}.{key}",
                resource_schema=resource_schema,
            )
            for key, item in raw.items()
        }
    if origin is Literal:
        if value not in args:
            raise CanonicalRequestDecodingError(
                f"Typed resource literal at {path} has unsupported value {value!r}."
            )
        return value
    if not _matches_type(value, hint):
        raise CanonicalRequestDecodingError(
            f"Typed resource value at {path} has an invalid type."
        )
    return value


def _resolve_resource(
    resource_id: object,
    *,
    path: str,
    resource_paths: Mapping[str, str | Path],
) -> Path:
    if not isinstance(resource_id, str) or not _RESOURCE_ID_RE.fullmatch(resource_id):
        raise CanonicalRequestDecodingError(f"Invalid resource ID at {path}.")
    if resource_id not in resource_paths:
        raise CanonicalRequestDecodingError(
            f"Canonical resource {resource_id!r} referenced at {path} is missing."
        )
    raw_path = resource_paths[resource_id]
    if not isinstance(raw_path, (str, Path)) or not str(raw_path).strip():
        raise CanonicalRequestDecodingError(
            f"Materialized path for resource {resource_id!r} is invalid."
        )
    materialized = Path(str(raw_path))
    if not materialized.is_file():
        raise CanonicalRequestDecodingError(
            f"Materialized resource is not a file: {resource_id!r}."
        )
    return materialized


def _object(
    value: object,
    *,
    path: str,
    required: set[str] | frozenset[str] | None = None,
    exact: bool = True,
) -> Mapping[str, Any]:
    if not isinstance(value, MappingABC) or not all(
        isinstance(key, str) for key in value
    ):
        raise CanonicalRequestDecodingError(f"{path} must be an object.")
    if required is not None:
        missing = set(required) - set(value)
        if missing:
            raise CanonicalRequestDecodingError(
                f"Missing required field(s) at {path}: {', '.join(sorted(missing))}."
            )
        if exact:
            _require_exact_fields(value, path=path, required=set(required))
    return value


def _require_exact_fields(
    value: Mapping[str, Any],
    *,
    path: str,
    required: set[str],
    optional: set[str] | frozenset[str] = frozenset(),
) -> None:
    missing = required - set(value)
    if missing:
        raise CanonicalRequestDecodingError(
            f"Missing required field(s) at {path}: {', '.join(sorted(missing))}."
        )
    unknown = set(value) - required - set(optional)
    if unknown:
        raise CanonicalRequestDecodingError(
            f"Unknown field(s) at {path}: {', '.join(sorted(unknown))}."
        )


def _array(value: object, *, path: str) -> list[Any]:
    if not isinstance(value, list):
        raise CanonicalRequestDecodingError(f"{path} must be an array.")
    return value


def _sequence(value: object, *, name: str) -> Sequence[Any]:
    if not isinstance(value, SequenceABC) or isinstance(value, (str, bytes, bytearray)):
        raise CanonicalRequestEncodingError(f"{name} must be a sequence.")
    return value


def _json_value(value: object, *, path: str) -> Any:
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise CanonicalRequestEncodingError(f"{path} must be finite.")
        return value
    if isinstance(value, MappingABC):
        if not all(isinstance(key, str) for key in value):
            raise CanonicalRequestEncodingError(f"{path} requires string mapping keys.")
        return {
            key: _json_value(item, path=f"{path}.{key}")
            for key, item in value.items()
        }
    if isinstance(value, SequenceABC) and not isinstance(value, (str, bytes, bytearray)):
        return [
            _json_value(item, path=f"{path}[{index}]")
            for index, item in enumerate(value)
        ]
    raise CanonicalRequestEncodingError(
        f"{path} contains unsupported value type {type(value).__name__}."
    )


def _validate_dataclass_contract(value: object, *, path: str, error: str) -> None:
    if not is_dataclass(value) or isinstance(value, type):
        _raise_contract_error(error, f"{path} must be a dataclass value.")
    hints = get_type_hints(type(value))
    for item in fields(value):
        raw = getattr(value, item.name)
        hint = hints.get(item.name, Any)
        if not _matches_type(raw, hint):
            _raise_contract_error(
                error,
                f"{path}.{_camel(item.name)} does not match its typed contract.",
            )
        if is_dataclass(raw) and not isinstance(raw, type):
            _validate_dataclass_contract(
                raw, path=f"{path}.{_camel(item.name)}", error=error
            )


def _matches_type(value: object, hint: object) -> bool:
    if hint is Any or hint is object:
        return True
    origin = get_origin(hint)
    args = get_args(hint)
    if origin in {types.UnionType, Union}:
        return any(_matches_type(value, branch) for branch in args)
    if origin is Literal:
        return value in args and any(type(value) is type(item) for item in args if item == value)
    if origin in {tuple, list, Sequence, SequenceABC}:
        if not isinstance(value, SequenceABC) or isinstance(value, (str, bytes, bytearray)):
            return False
        if not args:
            return True
        item_hint = args[0]
        return all(_matches_type(item, item_hint) for item in value)
    if origin in {dict, Mapping, MappingABC}:
        if not isinstance(value, MappingABC):
            return False
        if len(args) != 2:
            return True
        return all(
            _matches_type(key, args[0]) and _matches_type(item, args[1])
            for key, item in value.items()
        )
    if hint is float:
        return (
            isinstance(value, (int, float))
            and not isinstance(value, bool)
            and math.isfinite(float(value))
        )
    if hint is int:
        return isinstance(value, int) and not isinstance(value, bool)
    if hint is bool:
        return isinstance(value, bool)
    if hint is str:
        return isinstance(value, str)
    if hint is type(None):
        return value is None
    return isinstance(hint, type) and isinstance(value, hint)


def _raise_contract_error(kind: str, message: str) -> None:
    if kind == "encode":
        raise CanonicalRequestEncodingError(message)
    raise CanonicalRequestDecodingError(message)


def _same_default(value: object, default: object) -> bool:
    if value is default:
        return True
    if isinstance(value, DataFrame) or isinstance(default, DataFrame):
        return False
    try:
        return bool(value == default)
    except (TypeError, ValueError):
        return False


def _resource_id(value: str) -> str:
    normalized = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    if not normalized or not normalized[0].isalpha():
        normalized = f"resource-{normalized}".rstrip("-")
    return normalized


def _camel(value: str) -> str:
    head, *tail = value.split("_")
    return head + "".join(part[:1].upper() + part[1:] for part in tail)


def _integer(value: object, path: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise CanonicalRequestDecodingError(f"{path} must be an integer.")
    return value


def _optional_integer(value: object, path: str) -> int | None:
    return None if value is None else _integer(value, path)


def _non_negative_index(value: object, path: str) -> int:
    index = _integer(value, path)
    if index < 0:
        raise CanonicalRequestDecodingError(f"{path} must be non-negative.")
    return index


def _boolean(value: object, path: str) -> bool:
    if not isinstance(value, bool):
        raise CanonicalRequestDecodingError(f"{path} must be a boolean.")
    return value


def _optional_string(value: object, path: str) -> str | None:
    if value is not None and not isinstance(value, str):
        raise CanonicalRequestDecodingError(f"{path} must be a string or null.")
    return value


def _required_string(value: object, path: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise CanonicalRequestDecodingError(f"{path} must be a non-empty string.")
    return value.strip()


def _check_singleton_kind(kind: str, seen: set[str], *, path: str) -> None:
    if kind in seen:
        raise CanonicalRequestDecodingError(
            f"Duplicate singleton comparison kind at {path}: {kind}."
        )
    seen.add(kind)


__all__ = [
    "CANONICAL_REQUEST_SCHEMA",
    "SUPPORTED_CANONICAL_REQUEST_SCHEMAS",
    "UNKNOWN_FIELD_POLICY",
    "CanonicalRequestCodecError",
    "CanonicalRequestDecodingError",
    "CanonicalRequestEncodingError",
    "CanonicalRequestResource",
    "EncodedCanonicalRequest",
    "decode_canonical_request",
    "encode_canonical_typed_resource",
    "encode_canonical_request",
]
