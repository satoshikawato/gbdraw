"""Internal codec for canonical, CLI-independent render request payloads.

Canonical request schemas are intentionally strict: unknown fields are rejected instead of being
silently discarded.  File content is represented by stable resource IDs; this
module does not own session envelopes, base64 encoding, or temporary files.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping as MappingABC, Sequence as SequenceABC
from copy import deepcopy
from dataclasses import MISSING, Field, asdict, dataclass, fields, is_dataclass, replace
from functools import lru_cache
from io import StringIO
import json
import logging
import math
from pathlib import Path
import re
import types
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    Mapping,
    NamedTuple,
    NoReturn,
    Sequence,
    Union,
    cast,
    get_args,
    get_origin,
    get_type_hints,
    overload,
)

from Bio import SeqIO
from pandas import DataFrame, read_csv

from gbdraw.analysis.collinearity import (
    COLLINEARITY_ANCHOR_MODE_SPELLINGS,
    COLLINEARITY_ANCHOR_MODES,
    COLLINEARITY_COLOR_MODE_SPELLINGS,
    COLLINEARITY_COLOR_MODES,
    COLLINEARITY_MERGE_ORIENTATION_SPELLINGS,
    COLLINEARITY_MERGE_ORIENTATIONS,
    COLLINEARITY_SEARCH_SCOPE_SPELLINGS,
    COLLINEARITY_SEARCH_SCOPES,
    RETIRED_COLLINEARITY_ANCHOR_MODES,
    CollinearityAnchor,
    CollinearityBlock,
    CollinearityResult,
    LosslessCollinearityParameters,
)
from gbdraw.analysis.collinearity_units import _COLLINEARITY_UNIT_MODES, COLLINEARITY_UNIT_MODE_SPELLINGS
from gbdraw.analysis.conservation import _default_label as _default_conservation_label
from gbdraw.analysis.ortholog_paths import OrthologPathCollection
from gbdraw.analysis.protein_colinearity import (
    ORTHOGROUP_MEMBERSHIP_MODE_SPELLINGS,
    ORTHOGROUP_MEMBERSHIP_MODES,
    OrthogroupMember,
    OrthogroupNameCandidate,
    OrthogroupResult,
    OrthogroupGraphResult,
    compact_ortholog_paths,
    OrthologEdge,
    OrthologPath,
    normalize_protein_blastp_mode,
)
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.models.root import read_stored_font_size, style_leaf_domain
from gbdraw.config.toml import load_config_toml
from gbdraw.comparisons.losat_jobs import record_source_paths
from gbdraw.core.numeric import ascii_float, parse_ascii_nonnegative_integer
from gbdraw.session_legacy_values import (
    LEGACY_FEATURE_SHAPES,
    LEGACY_FEATURE_TYPES,
    LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES,
)
from gbdraw.exceptions import ValidationError
from gbdraw.session_contract import SESSION_ID_RE
from gbdraw.features.overrides import FeatureOverride
from gbdraw.features.placement import (
    _LANE_SIDES,
    FeaturePlacementOverride,
    FeaturePlacementTarget,
)
from gbdraw.io.record_select import RecordSelector, parse_record_selector
from gbdraw.io.regions import RegionSpec
from gbdraw.io.comparisons import read_comparison_table
from gbdraw.io.table_text import EXACT_FLOAT_READ_OPTIONS
from gbdraw.core.record_metadata import _read_coord_map
from gbdraw.linear_comparison import (
    LinearComparison,
    reverse_endpoint_table_text,
    reverse_unbound_endpoint_rows,
)
from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot, ScalarSpec
from gbdraw.tracks.circular import (
    _LEGACY_SPACING_PARAM,
    _InternalCircularTrackSlot,
    _migrate_legacy_circular_track_slot_spec,
    _internal_circular_track_slot,
    _normalize_renderer as _circular_renderer_name,
    parse_circular_track_slot,
)
from gbdraw.tracks.linear import _normalize_renderer as _linear_renderer_name
from gbdraw.annotations import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    FeatureIdentitySpan,
    FeatureSelector,
    FeatureSpan,
    HatchStyle,
    RegionAnnotation,
    RegionAnnotationStyle,
    annotation_track_params_from_mapping,
)

from .api.options import (
    DEPTH_TRACK_POSITIVE_FIELDS,
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    CircularOutputOptions,
    CircularRequestTrackOptions,
    ColorOptions,
    DepthTrackInput,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
    LinearRecordTranslation,
    LinearOutputOptions,
    LinearRequestTrackOptions,
    LOSATP_MODE_WIRE,
    LosatRuntimeOptions,
    LosatSearchOptions,
    _invalid_input,
    losatp_analysis_mode,
)
from .comparisons.losat_runtime import AUTOMATIC_LOSAT_BIN
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
    RecordDisplayOptions,
    RecordInput,
    RecordPresentation,
    RenderOutputRequest,
    _LegacySimilarityAlignment,
    _with_legacy_similarity_alignment,
)
from .layout.similarity_alignment import (
    SIMILARITY_ALIGNMENT_PLAN_SCHEMA,
    AlignmentAnchorIdentity,
    AlignmentDecisionStatus,
    AlignmentRecordDecision,
    AlignmentResolutionRationale,
    SimilarityAlignmentPlan,
    SimilarityAlignmentReference,
)


if TYPE_CHECKING:
    from gbdraw.annotations.models import RegionTarget
    from gbdraw.session_resources import RequestResources, SessionResourceTable

    from .api.options import LosatpMode


CANONICAL_REQUEST_SCHEMA = 9
DISPLAY_PLACEMENT_SCHEMA = 7
FEATURE_OVERRIDE_SCHEMA = 9
SUPPORTED_CANONICAL_REQUEST_SCHEMAS = frozenset(
    {1, 2, 5, 6, 7, 8, CANONICAL_REQUEST_SCHEMA}
)
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
            **{
                lossless: getattr(self, standard)
                for standard, lossless in _LEGACY_STANDARD_COLLINEARITY_FIELDS.items()
            },
            merge_orientation=_LEGACY_STANDARD_MERGE_ORIENTATION,
        )


_TOP_LEVEL_FIELDS = frozenset(
    {"schema", "mode", "records", "diagramOptions", "layout", "comparisons", "output"}
)
_TOP_LEVEL_FIELDS_V5 = _TOP_LEVEL_FIELDS | {"grouping"}
_DEFAULT_CIRCULAR_OPTIONS = CircularDiagramOptions()
_DEFAULT_LINEAR_OPTIONS = LinearDiagramOptions()
logger = logging.getLogger(__name__)
_SHARED_OPTION_WRONG_MODE_DEFAULTS: dict[str, dict[str, object]] = {
    "circular": {
        "depth_track_heights": None,
        "pairwise_match_style": "ribbon",
    },
    "linear": {
        "conservation_blast_files": None,
        "conservation_sequence_files": None,
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
_LEGACY_SPARSE_CONFIG_OVERRIDES: dict[str, dict[str, object]] = {
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
# Schemas 1-2 diagramOptions names of the feature visibility table.
_LEGACY_FEATURE_VISIBILITY_RENAMES = (
    ("featureTable", "featureVisibilityTable"),
    ("featureTableFile", "featureVisibilityTableFile"),
)
# A ``kind=standard`` collinearity payload (schemas 1-2) read as lossless
# parameters: standard field -> lossless field, and the merge orientation.
_LEGACY_STANDARD_COLLINEARITY_FIELDS = {
    "min_anchors": "min_anchors",
    "max_gene_gap": "max_unit_gap",
    "max_diagonal_drift": "max_diagonal_drift",
    "max_conflicts_in_merge_gap": "max_conflicts",
}
# The standard payload's merge rule, which the lossless "strand" keeps (OV-2301).
_LEGACY_STANDARD_MERGE_ORIENTATION: Literal["strand"] = "strand"
# The paralog link setting a ``kind=standard`` payload's own value replaces
# in collinear mode, when the setting holds this default.
_LEGACY_PARALOG_LINKS_DEFAULT = 2
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
        "orthogroups",
        "collinearity_blocks",
    }
)
# ``generatedProteinComparison.settings``: persisted wire key -> typed owner
# and field (design D6). Wire keys never change with typed names. "options"
# is LinearDiagramOptions, "search" LosatSearchOptions, "runtime" its
# LosatRuntimeOptions. The order is the persisted key order.
_PROTEIN_SETTINGS_WIRE: tuple[tuple[str, str, str], ...] = (
    ("collinearityParams", "options", "collinearity_params"),
    ("collinearityUnitMode", "options", "collinearity_unit_mode"),
    ("collinearityAnchorMode", "options", "collinearity_anchor_mode"),
    ("collinearitySearchScope", "options", "collinearity_search_scope"),
    ("collinearityColorMode", "options", "collinearity_color_mode"),
    ("losatpBin", "runtime", "losat_executable"),
    ("ncbiBlastpBin", "runtime", "ncbi_blast_executable"),
    ("losatpThreads", "runtime", "threads"),
    ("proteinBlastpMaxHits", "search", "losatp_max_hits"),
    ("proteinBlastpCandidateLimit", "search", "losatp_max_target_seqs"),
    ("orthogroupMembershipMode", "options", "orthogroup_membership_mode"),
    ("orthogroupMemberMaxHits", "search", "losatp_member_max_hits"),
    ("collinearInferOrthogroups", "options", "collinear_infer_orthogroups"),
    (
        "collinearMaxParalogLinksPerOrthogroup",
        "options",
        "collinear_max_paralog_links_per_orthogroup",
    ),
)
# From this schema a comparison mode setting holds one of its names; older
# schemas hold any spelling the option's normalizer reads (OV-651).
_SETTING_NAMES_SCHEMA = 9
# Schemas 1-7 also stored the retired alignment target in the settings.
_LEGACY_ALIGNMENT_SETTING = "alignOrthogroupFeature"
_OPTIONAL_PROTEIN_SETTINGS = frozenset({"collinearInferOrthogroups"})
# What a request without an optional setting held. A request written before
# `collinearInferOrthogroups` existed always inferred Collinear groups; the
# reading stays when the option's default changes.
_ABSENT_PROTEIN_SETTINGS: dict[str, object] = {"collinearInferOrthogroups": True}
_WIRE_LOSATP_MODES = cast(
    "dict[str, LosatpMode | Literal['none']]",
    {wire: typed for typed, wire in LOSATP_MODE_WIRE.items()},
)
# Settings of a LOSATP search that has not been configured.
_UNSET_LOSAT_SEARCH = LosatSearchOptions(program="losatp", losatp_mode="none")
# Circular ring LOSAT intent (encoded only after the planner resolves it) and
# its resolved rows (encoded as conservationBlastFiles).
_RING_SEARCH_FIELDS = frozenset({"conservation_losat_gencodes", "conservation_search_results"})
_COMPARISON_FIELDS = _COMPARISON_SOURCE_FIELDS | frozenset(
    {"losat_search"}
    | {name for _key, owner, name in _PROTEIN_SETTINGS_WIRE if owner == "options"}
)

_TABLE_FIELDS = frozenset(
    {
        "feature_visibility_table",
        "label_whitelist_table",
        "qualifier_priority_table",
        "label_override_table",
        "depth_table",
    }
)
# The text styling tables. Every cell is text, as `read_literal_table` reads a
# table file: an empty cell is "" and "NA" or "nan" is that text.
_TEXT_TABLE_FIELDS = frozenset(
    {
        "feature_visibility_table",
        "label_whitelist_table",
        "qualifier_priority_table",
        "label_override_table",
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
_OPTIONAL_FILE_SEQUENCE_FIELDS = frozenset({"conservation_sequence_files"})
# Typed option -> persisted snake name, where the typed name changed and the
# wire name did not (design D6, D18). Keys and resource IDs use the wire name.
_OPTION_WIRE_NAMES = {"conservation_sequence_files": "conservation_fasta_files"}


def _option_wire_name(name: str) -> str:
    return _OPTION_WIRE_NAMES.get(name, name)


def _option_wire_key(name: str) -> str:
    return _camel(_option_wire_name(name))
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
# A comparison table path: encode_canonical_request refuses it (call
# resolve_request first) and no request of any schema holds it (OV-324).
# A request field of the wrong kind or name.
_REQUEST_FIELDS_INVALID = _invalid_input("schema", "FIELDS")
_TABLE_PATH_FIELDS = frozenset({"conservation_table_file", "comparison_table_file"})
# Options that hold a name. Their constructors read null (and an empty array)
# as the default; every writer wrote a name (OV-528, OV-596).
_NAME_OPTIONS = frozenset({"conservation_reference", "pairwise_match_style"})

_TYPED_TREE_CLASSES = {
    cls.__name__: cls
    for cls in (
        CollinearityAnchor,
        CollinearityBlock,
        CollinearityResult,
        OrthogroupMember,
        OrthogroupNameCandidate,
        OrthogroupResult,
        OrthogroupGraphResult,
        OrthologEdge,
        OrthologPath,
    )
}
# The resource schemas a typed JSON resource may name; a class may need one of them.
_TYPED_RESOURCE_SCHEMAS = (1, 2, 3)
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
        if not SESSION_ID_RE.fullmatch(self.resource_id):
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
    """JSON-compatible request payload plus the resources it added."""

    payload: dict[str, Any]
    resources: tuple[CanonicalRequestResource, ...]


def encode_canonical_request(
    request: DiagramRequest,
    *,
    table: SessionResourceTable | None = None,
) -> EncodedCanonicalRequest:
    """Encode one typed request and normalize all conversion failures.

    The request's resources are allocated in ``table``, a Session's resource
    table, or alone without one; ``resources`` lists the ones it added.
    """

    try:
        return _encode_canonical_request(request, table=table)
    except CanonicalRequestCodecError:
        raise
    except (TypeError, ValueError, ValidationError) as exc:
        raise CanonicalRequestEncodingError(
            f"Typed request could not be converted: {exc}"
        ) from exc


def _encode_canonical_request(
    request: DiagramRequest,
    *,
    table: SessionResourceTable | None,
) -> EncodedCanonicalRequest:
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
    if isinstance(request, (CircularDiagramRequest, CircularBatchRequest)):
        if request.options.conservation_table_file is not None:
            unresolved_reasons.append("Circular comparison table")
        ring_search = request.options.losat_search
        if ring_search is not None:
            unresolved_reasons.append(f"{ring_search.program} ring search")
    if isinstance(request, LinearDiagramRequest):
        if request.output.resolve_prefix_from_first_record:
            unresolved_reasons.append("record-derived output prefix")
        if request.options.comparison_table_file is not None:
            unresolved_reasons.append("Linear comparison table")
        if isinstance(request.similarity_alignment, SimilarityAlignmentReference):
            unresolved_reasons.append("similarity alignment reference")
        search = request.options.losat_search
        if search is not None and search.program != "losatp":
            unresolved_reasons.append(f"{search.program} search")
        if request._legacy_similarity_alignment is not None:
            raise CanonicalRequestEncodingError(
                "A legacy similarity alignment must be materialized before current encoding."
            )
    if unresolved_reasons:
        raise CanonicalRequestEncodingError(
            "The current canonical schema cannot encode unresolved request transforms; "
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
    layout = getattr(request, "layout", None)
    if layout is not None:
        _validate_dataclass_contract(layout, path="layout", error="encode")

    from gbdraw.session_resources import SessionResourceTable

    resources = (table if table is not None else SessionResourceTable()).request()
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
    # A Circular batch or grid draft has no per-record selection, crop, or
    # orientation.
    records, frames = _encode_records(
        request.records,
        resources=resources,
        per_record=mode == "linear" or grouping == "single",
    )
    payload = {
        "schema": CANONICAL_REQUEST_SCHEMA,
        "mode": mode,
        "grouping": grouping,
        "records": records,
        "diagramOptions": _encode_diagram_options(
            request.options,
            record_count=len(request.records),
            resources=resources,
        ),
        "layout": _encode_layout(request),
        "comparisons": _encode_comparisons(
            request.options,
            mode=mode,
            resources=resources,
            frames=frames,
        ),
        "output": (
            [_encode_output(output) for output in request.outputs]
            if isinstance(request, CircularBatchRequest)
            else _encode_output(request.output)
        ),
    }
    return EncodedCanonicalRequest(payload=payload, resources=resources.added())


def decode_canonical_request(
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    output_directory: str | Path,
) -> DiagramRequest:
    """Decode every supported canonical schema and normalize validation failures.

    A current-schema request must also match its description
    (:func:`describe_canonical_request`): it is checked before the typed
    constructors convert any value, and refused after the decoder's own
    checks, whose messages are more specific.
    """

    try:
        problems = _undescribed_parts(payload)
        request = _decode_canonical_request(
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
    if problems:
        shown = "; ".join(problems[:5]) + (f"; and {len(problems) - 5} more" if len(problems) > 5 else "")
        raise CanonicalRequestDecodingError(
            f"The schema {CANONICAL_REQUEST_SCHEMA} request does not match its description: {shown}.",
            diagnostic=_REQUEST_FIELDS_INVALID,
        )
    return request


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
        _decode_record(
            value, index=index, mode=mode, schema=schema, resource_paths=resource_paths
        )
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
    legacy_similarity_alignment = comparison_kwargs.pop(
        "_legacy_similarity_alignment",
        None,
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
        assert isinstance(options, LinearDiagramOptions)
        output = _decode_output(top["output"], output_directory=output_directory)
        if schema == 1:
            layout_payload = _object(top["layout"], path="renderRequest.layout")
            if layout_payload:
                raise CanonicalRequestDecodingError(
                    "A schema 1 Linear canonical request must have an empty layout object."
                )
            linear_layout = None
        else:
            linear_layout, similarity_alignment = _decode_linear_layout(
                top["layout"],
                schema=schema,
            )
        if schema == 1:
            similarity_alignment = None
        if schema < 6 and linear_layout is not None and linear_layout.multi_record_positions:
            # Schemas 2-5 wrote one ``#<card>@<row>`` per record card; as Web
            # Load promotes them, each card's row becomes its gridRow, which
            # the records an ALL card expands to inherit.
            rows = [
                parse_ascii_nonnegative_integer(token.rpartition("@")[2])
                for token in linear_layout.multi_record_positions
            ]
            records = tuple(
                replace(record, presentation=replace(record.presentation, grid_row=row))
                if row is not None and row > 0
                else record
                for record, row in zip(records, [*rows, *[None] * len(records)])
            )
            linear_layout = replace(linear_layout, multi_record_positions=None)
        request = LinearDiagramRequest(
            records=records,
            options=options,
            layout=linear_layout,
            similarity_alignment=similarity_alignment,
            output=output,
        )
        if legacy_similarity_alignment is not None:
            _with_legacy_similarity_alignment(
                request,
                legacy_similarity_alignment,
            )
        return request

    assert isinstance(options, CircularDiagramOptions)
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


def _genbank_bytes(records: Sequence[Any], *, index: int) -> bytes:
    stream = StringIO()
    try:
        serializable = [deepcopy(record) for record in records]
        for record in serializable:
            record.annotations.setdefault("molecule_type", "DNA")
        SeqIO.write(serializable, stream, "genbank")
    except Exception as exc:
        raise CanonicalRequestEncodingError(
            f"Could not serialize in-memory record {index} as GenBank."
        ) from exc
    return stream.getvalue().encode("utf-8")


def _shared_record_selector(
    record_id: str, *, record_index: int, unique: bool
) -> RecordSelector:
    """A record of a shared resource: its ID when unique and usable, else its index."""

    if unique and not record_id.startswith("#"):
        selector = parse_record_selector(record_id)
        if selector is not None and selector.record_id == record_id:
            return selector
    return RecordSelector(
        raw=f"#{record_index + 1}", record_id=None, record_index=record_index
    )


# The resources of a source file read by the planner: (resource kind, source field).
_SOURCE_FILE_RESOURCES: Mapping[str, tuple[tuple[str, str], ...]] = {
    "genbank": (("genbank", "resourceId"),),
    "gff_fasta": (("gff3", "gffResourceId"), ("fasta", "fastaResourceId")),
}


def _source_file_transform(record: Any) -> tuple[tuple[int, int] | None, bool] | None:
    """The crop and orientation that redraw ``record`` from its source file record.

    Each crop and reverse complement the planner applies moves the record's
    coordinate map (the source base of display base 1 and its step), and a
    crop marks the record; nothing else changes the bases it draws. Returns the
    cropped source span (``None`` when uncropped) and whether the drawn record
    is reverse-complemented, or ``None`` when the map is not a span of the
    source record.
    """

    annotations = record.annotations
    base, step = _read_coord_map(record)
    length = len(record)
    source_length = annotations.get("gbdraw_source_length", length)
    first, last = (base, base + length - 1) if step > 0 else (base - length + 1, base)
    if not isinstance(source_length, int) or not 1 <= first <= last <= source_length:
        return None
    if annotations.get("gbdraw_region_applied"):
        return (first, last), step < 0
    if (first, last) != (1, source_length):
        return None
    return None, step < 0


def _source_file_transforms(
    records: Sequence[Any], paths: tuple[str, ...], *, per_record: bool
) -> list[tuple[tuple[int, int] | None, bool]] | None:
    """Each record's crop and orientation when the records are drawn from their files.

    The planner annotates each record it reads with its source kind and its
    index among the file's records. Returns ``None`` unless each record of
    the files is drawn at most once and each is the file's record. Only when
    ``per_record`` allows it may some records of the files be left undrawn,
    or a record be cropped or reverse-complemented.
    """

    annotations = [record.annotations for record in records]
    kind = annotations[0].get("gbdraw_source_kind")
    count = annotations[0].get("gbdraw_source_record_count")
    indexes = sorted(item.get("gbdraw_source_record_index", -1) for item in annotations)
    if not (
        len(paths) == len(_SOURCE_FILE_RESOURCES.get(str(kind), ()))
        and isinstance(count, int)
        and all(Path(path).is_file() for path in paths)
        and (indexes == list(range(count)) or (
            per_record
            and len(set(indexes)) == len(indexes)
            and 0 <= indexes[0]
            and indexes[-1] < count
        ))
        and all(
            item.get("gbdraw_source_kind") == kind
            and item.get("gbdraw_source_record_count") == count
            for item in annotations
        )
    ):
        return None
    found = [_source_file_transform(record) for record in records]
    if None in found or not (per_record or all(item == (None, False) for item in found)):
        return None
    # None is not in found here.
    return cast("list[tuple[tuple[int, int] | None, bool]]", found)


def _source_region(
    span: tuple[int, int], reverse: bool, selector: RecordSelector | None
) -> RegionSpec:
    """A source-coordinate crop of the selected record, reverse-complemented when drawn so."""

    coordinates = f"{span[0]}-{span[1]}" + (":rc" if reverse else "")
    return RegionSpec(
        raw=f"{selector.raw}:{coordinates}" if selector is not None else coordinates,
        file_selector=None,
        record_id=selector.record_id if selector is not None else None,
        record_index=selector.record_index if selector is not None else None,
        start=span[0],
        end=span[1],
        reverse_complement=reverse,
    )


def _encode_records(
    records: Sequence[RecordInput],
    *,
    resources: RequestResources,
    per_record: bool,
) -> tuple[list[dict[str, Any]], tuple[tuple[int, bool], ...]]:
    """Encode records; the in-memory records of one source file share a resource.

    The planner projects each displayed record to an in-memory record. When
    they are records of their source file, each drawn once, that file's own
    bytes are the resource, and each record's crop (in source coordinates) and
    orientation are its region and reverse complement, as the Web writes them:
    the CLI binds the same file for the Web draft, so Session Load sees the
    request draw from its bound inputs. ``per_record`` says whether the Web
    draft selects, crops, and orients each record (every Linear record, and a
    single Circular record); otherwise only every record of the file, each
    unchanged, reads the file. Other records whose provenance names one source
    file are written as drawn to one GenBank resource named after that file.
    Records of a multi-record file are selected by record ID (by index when an
    ID is repeated or unusable), the shape the Web writes for the records it
    draws of one uploaded multi-record file, so a replay or Web load keeps one
    file as one File and one LOSAT source
    (:func:`gbdraw.comparisons.losat_jobs.losat_source_ids`) and draws only
    the selected records. Any other in-memory record is its own source and
    keeps its own resource.

    Also returns each record's persisted comparison frame ``(L, reversed)``
    (PD-OI-073): a record written as its reverse-complemented sequence has the
    displayed search frame, so comparison rows that touch it are written in
    that frame; a source file with a reverse complement keeps its rows.
    """

    groups: dict[tuple[str, ...], list[int]] = {}
    for position, record in enumerate(records):
        if isinstance(record.source, InMemoryRecordSource):
            paths = record_source_paths(record.source.record)
            if paths:
                groups.setdefault(paths, []).append(position)
    shared: dict[int, tuple[dict[str, Any], RecordSelector | None]] = {}
    transformed: dict[int, RecordInput] = {}
    for paths, members in groups.items():
        first = members[0] + 1
        # Only in-memory sources were grouped above.
        members_records = [
            cast(InMemoryRecordSource, records[position].source).record for position in members
        ]
        found = _source_file_transforms(members_records, paths, per_record=per_record)
        if found is not None:
            file_resources = _SOURCE_FILE_RESOURCES[
                members_records[0].annotations["gbdraw_source_kind"]
            ]
            source_payload: dict[str, Any] = {
                "kind": "genbank" if len(file_resources) == 1 else "gffFasta"
            }
            for path, (kind, field) in zip(paths, file_resources, strict=True):
                source_payload[field] = resources.add_bytes(
                    f"record-{first}-{kind}",
                    kind=kind,
                    name=Path(path).name,
                    content=Path(path).read_bytes(),
                )
            ids = [str(record.annotations.get("gbdraw_source_record_id", record.id))
                   for record in members_records]
            record_indexes = [record.annotations["gbdraw_source_record_index"]
                              for record in members_records]
            # The file's other records may share an ID with a drawn record.
            unique = [bool(record.annotations.get("gbdraw_source_record_id_unique"))
                      for record in members_records]
            single = members_records[0].annotations["gbdraw_source_record_count"] == 1
        elif len(members) > 1:
            source_payload = {"kind": "genbank", "resourceId": resources.add_bytes(
                f"record-{first}-genbank",
                kind="genbank",
                name=(
                    Path(paths[0]).name
                    if len(paths) == 1 and Path(paths[0]).name
                    else f"record-{first}.gbk"
                ),
                content=_genbank_bytes(members_records, index=first),
            )}
            ids = [str(record.id) for record in members_records]
            record_indexes = list(range(len(members)))
            unique = [ids.count(record_id) == 1 for record_id in ids]
            single = False
        else:
            continue
        for offset, (position, record_id, record_index, is_unique) in enumerate(
            zip(members, ids, record_indexes, unique, strict=True)
        ):
            selector = None if single else _shared_record_selector(
                record_id,
                record_index=record_index,
                unique=is_unique,
            )
            if found is not None:
                span, reverse = found[offset]
                record = records[position]
                transformed[position] = replace(
                    record,
                    region=None if span is None else _source_region(span, reverse, selector),
                    presentation=replace(
                        record.presentation, reverse_complement=reverse and span is None
                    ),
                )
                # A crop names its record in the region, as the Web writes it.
                selector = selector if span is None else None
            shared[position] = (source_payload, selector)
    payloads = [
        _encode_record(
            transformed.get(position, record),
            index=position + 1,
            resources=resources,
            shared_source=shared.get(position),
        )
        for position, record in enumerate(records)
    ]
    frames = tuple(
        (len(record.source.record), True)
        if position not in transformed
        and isinstance(record.source, InMemoryRecordSource)
        and _read_coord_map(record.source.record)[1] == -1
        else (0, False)
        for position, record in enumerate(records)
    )
    return payloads, frames


def _encode_record(
    record: RecordInput,
    *,
    index: int,
    resources: RequestResources,
    shared_source: tuple[dict[str, Any], RecordSelector | None] | None = None,
) -> dict[str, Any]:
    source = record.source
    selector = record.selector
    if shared_source is not None:
        source_payload, selector = shared_source
    elif isinstance(source, GenBankInputSource):
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
        resource_id = resources.add_bytes(
            f"record-{index}-genbank",
            kind="genbank",
            name=f"record-{index}.gbk",
            content=_genbank_bytes((source.record,), index=index),
        )
        source_payload = {"kind": "genbank", "resourceId": resource_id}
    else:  # pragma: no cover - RecordInput validates its source union.
        raise CanonicalRequestEncodingError(f"Unsupported record source at index {index}.")

    return {
        "recordKey": record.record_key or f"record-{index}",
        "cardinality": record.cardinality.value,
        "display": _encode_record_display(record.display),
        "source": source_payload,
        "selector": _encode_selector(selector),
        "region": _encode_region(record.region),
        "presentation": _encode_record_presentation(record.presentation),
    }


def _encode_record_display(display: RecordDisplayOptions) -> dict[str, Any]:
    return {"isCircular": display.is_circular, "startCoordinate": display.start_coordinate}


def _encode_record_presentation(presentation: RecordPresentation) -> dict[str, Any]:
    return {
        "label": presentation.label,
        "subtitle": presentation.subtitle,
        "reverseComplement": presentation.reverse_complement,
        "gridRow": presentation.grid_row,
        "gridColumn": presentation.grid_column,
    }


def legacy_record_cardinality(mode: object, record: Mapping[str, Any]) -> RecordCardinality:
    """The cardinality of a request schema 1-5 record, which saved none.

    As Web Load promotes it (``promoteRenderRequestToSchema8``), a Linear
    record without selector and region draws every record of its source; any
    other record draws one.
    """

    if mode == "linear" and record.get("selector") is None and record.get("region") is None:
        return RecordCardinality.ALL
    return RecordCardinality.EXACTLY_ONE


def _decode_record(
    value: object,
    *,
    index: int,
    mode: str,
    schema: int,
    resource_paths: Mapping[str, str | Path],
) -> RecordInput:
    path = f"renderRequest.records[{index - 1}]"
    required = {"source", "selector", "region", "presentation"}
    if schema >= 2:
        required.add("recordKey")
    if schema >= 6:
        required.add("cardinality")
    if schema >= DISPLAY_PLACEMENT_SCHEMA:
        required.add("display")
    item = _object(value, path=path, required=required)
    display = RecordDisplayOptions()
    if schema >= DISPLAY_PLACEMENT_SCHEMA:
        raw = _object(item["display"], path=f"{path}.display",
                      required={"isCircular", "startCoordinate"})
        display = RecordDisplayOptions(raw["isCircular"], raw["startCoordinate"])
    source_payload = _object(
        item["source"], path=f"{path}.source", required={"kind"}, exact=False
    )
    kind = source_payload["kind"]
    source: GenBankInputSource | GffFastaInputSource
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
            else legacy_record_cardinality(mode, item)
        ),
        selector=_decode_selector(item["selector"], path=f"{path}.selector"),
        region=_decode_region(item["region"], path=f"{path}.region"),
        presentation=presentation,
        display=display,
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


def _encode_alignment_anchor(anchor: AlignmentAnchorIdentity) -> dict[str, Any]:
    return {
        "recordKey": anchor.record_key,
        "biologicalFeatureId": anchor.biological_feature_id,
        "sourceFeatureIndex": anchor.source_feature_index,
        "stableFeatureSvgId": anchor.stable_feature_svg_id,
    }


def _encode_similarity_alignment_plan(
    plan: SimilarityAlignmentPlan | None,
) -> dict[str, Any] | None:
    if plan is None:
        return None
    return {
        "schema": plan.schema,
        "groupId": plan.group_id,
        "reference": _encode_alignment_anchor(plan.reference),
        "records": [
            {
                "recordKey": decision.record_key,
                "status": decision.status.value,
                "rationale": decision.rationale.value,
                "anchor": (
                    _encode_alignment_anchor(decision.anchor)
                    if decision.anchor is not None
                    else None
                ),
            }
            for decision in plan.records
        ],
    }


def _decode_alignment_anchor(
    value: object,
    *,
    path: str,
) -> AlignmentAnchorIdentity:
    payload = _object(
        value,
        path=path,
        required={
            "recordKey",
            "biologicalFeatureId",
            "sourceFeatureIndex",
            "stableFeatureSvgId",
        },
    )
    return AlignmentAnchorIdentity(
        record_key=payload["recordKey"],
        biological_feature_id=payload["biologicalFeatureId"],
        source_feature_index=payload["sourceFeatureIndex"],
        stable_feature_svg_id=payload["stableFeatureSvgId"],
    )


def _decode_similarity_alignment_plan(
    value: object,
    *,
    path: str,
) -> SimilarityAlignmentPlan | None:
    if value is None:
        return None
    payload = _object(
        value,
        path=path,
        required={"schema", "groupId", "reference", "records"},
    )
    raw_records = _array(payload["records"], path=f"{path}.records")
    decisions: list[AlignmentRecordDecision] = []
    for index, raw in enumerate(raw_records):
        decision_path = f"{path}.records[{index}]"
        decision = _object(
            raw,
            path=decision_path,
            required={
                "recordKey",
                "status",
                "rationale",
                "anchor",
            },
        )
        decisions.append(
            AlignmentRecordDecision(
                record_key=decision["recordKey"],
                status=AlignmentDecisionStatus(decision["status"]),
                rationale=AlignmentResolutionRationale(decision["rationale"]),
                anchor=(
                    _decode_alignment_anchor(
                        decision["anchor"],
                        path=f"{decision_path}.anchor",
                    )
                    if decision["anchor"] is not None
                    else None
                ),
            )
        )
    return SimilarityAlignmentPlan(
        schema=payload["schema"],
        group_id=payload["groupId"],
        reference=_decode_alignment_anchor(
            payload["reference"],
            path=f"{path}.reference",
        ),
        records=tuple(decisions),
    )


def _encode_layout(request: DiagramRequest) -> dict[str, Any]:
    if isinstance(request, CircularBatchRequest):
        return {}
    if isinstance(request, LinearDiagramRequest):
        if request.layout is None:
            return {}
        return {
            "recordGapPx": request.layout.record_gap_px,
            "multiRecordPositions": (
                list(request.layout.multi_record_positions)
                if request.layout.multi_record_positions is not None
                else None
            ),
            "recordTranslations": [
                {
                    "recordKey": translation.record_key,
                    "x": translation.x,
                    "y": translation.y,
                }
                for translation in request.layout.record_translations
            ],
            # encode_canonical_request rejects an unresolved reference first.
            "similarityAlignment": _encode_similarity_alignment_plan(
                cast("SimilarityAlignmentPlan | None", request.similarity_alignment)
            ),
        }
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
    if schema in {1, 2} and isinstance(size_mode, str):
        size_mode = LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES.get(size_mode, size_mode)
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
) -> tuple[LinearMultiRecordOptions | None, SimilarityAlignmentPlan | None]:
    layout = _object(value, path="renderRequest.layout")
    if not layout:
        return None, None
    if schema >= 8:
        _require_exact_fields(
            layout,
            path="renderRequest.layout",
            required={
                "recordGapPx",
                "multiRecordPositions",
                "recordTranslations",
                "similarityAlignment",
            },
        )
        raw_positions = layout["multiRecordPositions"]
        if raw_positions is not None:
            positions = _array(
                raw_positions,
                path="renderRequest.layout.multiRecordPositions",
            )
            if not all(isinstance(item, str) and item.strip() for item in positions):
                raise CanonicalRequestDecodingError(
                    "renderRequest.layout.multiRecordPositions must contain non-empty text."
                )
        else:
            positions = None
        raw_translations = _array(
            layout["recordTranslations"],
            path="renderRequest.layout.recordTranslations",
        )
        translations: list[LinearRecordTranslation] = []
        for index, raw in enumerate(raw_translations):
            path = f"renderRequest.layout.recordTranslations[{index}]"
            item = _object(
                raw,
                path=path,
                required={"recordKey", "x", "y"},
            )
            translations.append(
                LinearRecordTranslation(
                    record_key=item["recordKey"],
                    x=item["x"],
                    y=item["y"],
                )
            )
        result = LinearMultiRecordOptions(
            record_gap_px=layout["recordGapPx"],
            multi_record_positions=(
                tuple(positions) if positions is not None else None
            ),
            record_translations=tuple(translations),
        )
        _validate_dataclass_contract(result, path="layout", error="decode")
        return result, _decode_similarity_alignment_plan(
            layout["similarityAlignment"],
            path="renderRequest.layout.similarityAlignment",
        )
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
    return result, None


# The planner consumes it: a Session stores what the figure drew with.
_PLANNER_ONLY_FIELDS = frozenset({"wrap_defaults"})
_PLACEMENT_INPUT_FIELDS = frozenset({
    "feature_placements", "feature_placement_table", "feature_placement_table_file",
    "feature_overrides", "feature_override_table", "feature_override_table_file",
})


def _encode_diagram_options(
    options: CircularDiagramOptions | LinearDiagramOptions,
    *,
    record_count: int,
    resources: RequestResources,
) -> dict[str, Any]:
    default_options: CircularDiagramOptions | LinearDiagramOptions
    if isinstance(options, CircularDiagramOptions):
        default_options = _DEFAULT_CIRCULAR_OPTIONS
    elif isinstance(options, LinearDiagramOptions):
        default_options = _DEFAULT_LINEAR_OPTIONS
    else:
        raise CanonicalRequestEncodingError(
            "diagramOptions must use mode-specific options."
        )
    if (options.feature_placement_table is not None
            or options.feature_placement_table_file is not None):
        raise CanonicalRequestEncodingError(
            "Feature placement tables must be materialized before canonical encoding."
        )
    if (options.feature_override_table is not None
            or options.feature_override_table_file is not None):
        raise CanonicalRequestEncodingError(
            "Feature override tables must be materialized before canonical encoding."
        )
    depth_tracks = _canonical_depth_tracks_for_encoding(
        options,
        record_count=record_count,
    )
    result: dict[str, Any] = {"featurePlacements": [
        {
            "recordKey": row.record_key,
            "biologicalFeatureId": row.biological_feature_id,
            "placement": ({"kind": "main"} if row.target.kind == "main" else {
                "kind": "lane", "side": row.target.side, "level": row.target.level,
            }),
        }
        for row in options.feature_placements
    ], "featureOverrides": [row.to_mapping() for row in options.feature_overrides]}
    for item in fields(options):
        name = item.name
        if (
            name in _COMPARISON_FIELDS
            or name in _ALL_DEPTH_INPUT_FIELDS
            or name in _PLACEMENT_INPUT_FIELDS
            or name in _PLANNER_ONLY_FIELDS
        ):
            continue
        if name == "conservation_search_results":
            assert isinstance(options, CircularDiagramOptions)
            if options.conservation_search_results:
                # Planner-resolved ring rows are stored as the Web stores them.
                result["conservationBlastFiles"] = [
                    {
                        "resourceId": resources.add_bytes(
                            f"conservation-blast-files-{index}",
                            kind=f"conservation-blast-files-{index}",
                            name=ring.name,
                            content=ring.text.encode("utf-8"),
                        ),
                        "representation": "file",
                    }
                    for index, ring in enumerate(options.conservation_search_results, start=1)
                ]
            continue
        value = getattr(options, name)
        default = getattr(default_options, name)
        if _same_default(value, default):
            continue
        result[_option_wire_key(name)] = _encode_option_value(name, value, resources=resources)
    labels = _precomputed_ring_labels(options)
    if labels is not None:
        result["conservationLabels"] = _encode_option_value(
            "conservation_labels", labels, resources=resources
        )
    if depth_tracks is not None:
        result["depthTracks"] = _encode_depth_tracks(
            depth_tracks,
            resources=resources,
        )
    return result


def _precomputed_ring_labels(
    options: CircularDiagramOptions | LinearDiagramOptions,
) -> tuple[str, ...] | None:
    """The drawn label of every precomputed ring, stored so a replay never derives one.

    A Session without them predates D-03 and drew the full file name
    (``_decode_diagram_options``); a renamed resource would also change it.
    """

    files = tuple(getattr(options, "conservation_blast_files", None) or ())
    if not files:
        return None
    given = getattr(options, "conservation_labels", None)
    count = max(len(files), len(getattr(options, "conservation_dataframes", None) or ()))
    labels = []
    for index in range(count):
        label = str(given[index]) if given is not None and index < len(given) else ""
        path = files[index] if index < len(files) else None
        labels.append(label if label.strip() else _default_conservation_label(index, path))
    return tuple(labels)


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
    payload = with_automatic_scale_interval(payload)
    placements: tuple[FeaturePlacementOverride, ...] = ()
    if schema >= DISPLAY_PLACEMENT_SCHEMA:
        if not isinstance(payload.get("featurePlacements"), list):
            raise CanonicalRequestDecodingError("diagramOptions.featurePlacements must be an array.")
        placements = tuple(FeaturePlacementOverride.from_mapping(row)
                           for row in payload.pop("featurePlacements"))
        for row in placements:
            row.target.validate_mode(mode)
    overrides: tuple[FeatureOverride, ...] = ()
    if schema >= FEATURE_OVERRIDE_SCHEMA:
        if not isinstance(payload.get("featureOverrides"), list):
            raise CanonicalRequestDecodingError("diagramOptions.featureOverrides must be an array.")
        overrides = tuple(FeatureOverride.from_mapping(row) for row in payload.pop("featureOverrides"))
    for name, default in _SHARED_OPTION_WRONG_MODE_DEFAULTS[mode].items():
        key = _option_wire_key(name)
        if key in payload and payload[key] == default:
            payload.pop(key)
    options_type = (
        CircularDiagramOptions if mode == "circular" else LinearDiagramOptions
    )
    known = {
        _option_wire_key(item.name): item.name
        for item in fields(options_type)
        if item.name not in _COMPARISON_FIELDS
        and item.name not in _PLACEMENT_INPUT_FIELDS
        and item.name not in _RING_SEARCH_FIELDS
        and item.name not in _TABLE_PATH_FIELDS
        and item.name not in _PLANNER_ONLY_FIELDS
    }
    unknown = set(payload) - set(known)
    if unknown:
        raise CanonicalRequestDecodingError(
            "Unknown field(s) at renderRequest.diagramOptions: "
            + ", ".join(sorted(unknown))
            + ".",
            diagnostic=_REQUEST_FIELDS_INVALID,
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
    decoded["feature_placements"] = placements
    decoded["feature_overrides"] = overrides
    if "conservationBlastFiles" in payload and "conservationLabels" not in payload:
        # Writers before D-03 stored no label for a precomputed ring and drew its file name.
        files = decoded.get("conservation_blast_files") or ()
        count = max(len(files), len(decoded.get("conservation_dataframes") or ()))
        decoded["conservation_labels"] = tuple(
            Path(files[index]).name if index < len(files) else f"Conservation {index + 1}"
            for index in range(count)
        )
    if decoded.get("config_overrides") is not None:
        decoded["config_overrides"] = _decode_config_overrides(
            decoded["config_overrides"],
            mode=mode,
            schema=schema,
        )
    if schema in {1, 2}:
        feature_shapes = dict(decoded.get("feature_shapes") or {})
        for feature_type, shape in LEGACY_FEATURE_SHAPES.items():
            feature_shapes.setdefault(feature_type, shape)
        decoded["feature_shapes"] = feature_shapes
    if schema in {1, 2}:
        _restore_legacy_sparse_defaults(decoded, mode=mode)
    return decoded


def with_stored_font_sizes_as_auto(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Return a stored request with each font size <= 0 read as Auto, as Web Load does (Owner Q-13 a, OV-1282).

    The Session reader (``decode_session_request``) calls it; a fresh request,
    such as a Web Generate, keeps refusing such a size (PD-OI-081).

    A font size is a ``POSITIVE`` style leaf (``style_leaf_domain``), in
    ``config`` or ``configOverrides``. Auto is the value the request's base
    gives that path: an override leaf is dropped, and inside ``config`` (or an
    override that holds a table) the leaf takes the stored ``config`` value
    (for an override), else the packaged one, else ``None`` (a Linear
    definition line style then inherits its size). Changed containers are
    copies; the stored request stays as saved.
    """

    packaged: list[Mapping[str, object]] = []  # Read on the first replacement only.
    packaged_base = object()

    def lookup(tree: object, path: str) -> object:
        if tree is packaged_base:
            if not packaged:
                packaged.append(load_config_toml("gbdraw.data", "config.toml"))
            tree = packaged[0]
        for key in path.split("."):
            if not isinstance(tree, MappingABC) or key not in tree:
                return None
            tree = tree[key]
        return tree

    def read(value: object, path: str, where: str, bases: tuple[object, ...]) -> object:
        if isinstance(value, MappingABC):
            children = {
                key: read(child, f"{path}.{key}", f"{where}.{key}", bases) for key, child in value.items()
            }
            return children if any(children[key] is not value[key] for key in value) else value
        domain = style_leaf_domain(path)
        if (
            value is None  # Unset stays unset; only a stored number <= 0 is replaced.
            or domain is None
            or domain[0] != "POSITIVE"
            or read_stored_font_size(value, path=where) is not None
        ):
            return value
        return next((default for base in bases if (default := lookup(base, path)) is not None), None)

    result = dict(payload)
    stored = result.get("diagramOptions")
    if not isinstance(stored, MappingABC):
        return result
    options = dict(stored)
    result["diagramOptions"] = options
    config = options.get("config")
    if isinstance(config, MappingABC):
        options["config"] = {
            key: read(child, key, f"renderRequest.diagramOptions.config.{key}", (packaged_base,))
            for key, child in config.items()
        }
    overrides = options.get("configOverrides")
    if isinstance(overrides, MappingABC):
        bases = (options.get("config"), packaged_base)
        kept: dict[str, Any] = {}
        for key, value in overrides.items():
            where = f'renderRequest.diagramOptions.configOverrides["{key}"]'
            if isinstance(value, MappingABC):
                kept[key] = read(value, key, where, bases)
            elif read(value, key, where, ()) is value:
                kept[key] = value
            # Else dropped, as Web Load drops it: the base value applies.
        options["configOverrides"] = kept
    return result


def with_automatic_scale_interval(payload: Mapping[str, Any]) -> dict[str, Any]:
    """Return request ``diagramOptions`` with a scale interval <= 0 read as automatic.

    Every request writer drew it so; the CLI, the option objects, and the Web
    field refuse such a value as input (D-04), but a request keeps its meaning.
    The decoder and the Web Load check of the stored configuration use it.
    """

    def automatic(holder: object, key: str) -> object:
        if not isinstance(holder, Mapping):
            return holder
        value = holder.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or value > 0:
            return holder
        return {**holder, key: None}

    result = dict(payload)
    config = result.get("config")
    if isinstance(config, Mapping) and isinstance(objects := config.get("objects"), Mapping):
        scale = automatic(objects.get("scale"), "interval")
        if scale is not objects.get("scale"):
            result["config"] = {**config, "objects": {**objects, "scale": scale}}
    overrides = result.get("configOverrides")
    for key in ("objects.scale.interval", "scale_interval"):
        overrides = automatic(overrides, key)
    if overrides is not result.get("configOverrides"):
        result["configOverrides"] = overrides
    return result


def _migrate_legacy_feature_visibility_fields(
    payload: Mapping[str, Any],
) -> dict[str, Any]:
    """Move removed feature-table aliases at the legacy schema boundary."""

    migrated = dict(payload)
    for old_name, current_name in _LEGACY_FEATURE_VISIBILITY_RENAMES:
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
    return migrate_flat_config_overrides(overrides, mode=mode, path=path)


def migrate_flat_config_overrides(
    overrides: Mapping[str, Any],
    *,
    mode: Literal["circular", "linear"],
    path: str,
) -> dict[str, Any]:
    """Project persisted flat aliases onto schema-derived dotted leaves.

    The decoder reads the overrides of every stored request with it, and so
    does the Web Load check of a stored request (OV-799).
    """

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
        # The Web writer stored a null blacklist for a filter mode other than
        # Blacklist: no blacklist, not the packaged one.
        if value is None and legacy_name == "label_blacklist":
            value = []
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
            # The interval follows this font when the override is applied.
            assign(
                "objects.definition.circular.font_size",
                value,
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
    circular_scope: Literal["none", "outer"] = "outer" if normalized == "all" else "none"
    return circular_scope, cast(
        'Literal["none", "all", "first", "orthogroup_top"]', normalized
    )


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
    decoded.setdefault("selected_features_set", LEGACY_FEATURE_TYPES)

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
    resources: RequestResources,
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
            numeric = ascii_float(value)
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
    resources: RequestResources,
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


def _without_unwritten_feature_tolerance(name: str, value: object, error_type):
    if name not in {"config", "config_overrides"} or not isinstance(value, MappingABC):
        return value
    owner = value.get("canvas", {}) if name == "config" else value
    key = "feature_overlap_tolerance_bp" if name == "config" else "canvas.feature_overlap_tolerance_bp"
    if not isinstance(owner, MappingABC) or key not in owner:
        return value
    result = deepcopy(dict(value))
    owner = result["canvas"] if name == "config" else result
    tolerance = owner[key]
    if isinstance(tolerance, bool) or not isinstance(tolerance, int) or tolerance != 0:
        raise error_type(
            "Feature overlap tolerance cannot be represented by canonical schema 6; "
            "persistence requires the shared display/placement writer integration."
        )
    owner.pop(key)  # Zero is exactly the existing missing-value behavior.
    return result


def _encode_option_value(
    name: str,
    value: object,
    *,
    resources: RequestResources,
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
                _file_ref(f"{_option_wire_name(name)}-{index}", item, resources=resources)
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
    if schema < DISPLAY_PLACEMENT_SCHEMA:
        value = _without_unwritten_feature_tolerance(name, value, CanonicalRequestDecodingError)
    if name == "config":
        return migrate_legacy_full_config(
            _object(value, path="renderRequest.diagramOptions.config"),
            mode=mode,
            path="renderRequest.diagramOptions.config",
        )
    if name == "colors":
        return _decode_colors(value, resource_paths=resource_paths)
    if name == "tracks":
        return _decode_tracks(value, mode=mode, schema=schema)
    if name == "annotations":
        return _decode_annotations(value, schema=schema, resource_paths=resource_paths)
    if name == "output":
        return _decode_assembly_output(value, mode=mode, schema=schema)
    if name == "depth_tracks":
        return _decode_depth_tracks(value, resource_paths=resource_paths)
    if name in _TABLE_FIELDS:
        return _decode_table_ref(
            value,
            name=name,
            resource_paths=resource_paths,
            text_cells=name in _TEXT_TABLE_FIELDS,
        )
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
        raw = _array(value, path=f"renderRequest.diagramOptions.{_option_wire_key(name)}")
        return tuple(
            (
                str(
                    _decode_file_ref(
                        item,
                        name=f"{_option_wire_name(name)}-{index}",
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
    if name in _NAME_OPTIONS:
        return _string(value, f"renderRequest.diagramOptions.{_option_wire_key(name)}")
    return value


def migrate_legacy_full_config(
    value: Mapping[str, Any],
    *,
    mode: Literal["circular", "linear"],
    path: str,
) -> dict[str, Any]:
    """Move persisted config aliases onto the current typed schema.

    The decoder reads every stored full ``config`` with it, and so does the Web
    Load check of a stored request (OV-790).
    """

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


def _encode_colors(value: object, *, resources: RequestResources) -> dict[str, Any]:
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
                text_cells=True,
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
                text_cells=True,
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
        hatch = HatchStyle(**{**raw_hatch, "cross": _boolean(raw_hatch["cross"], f"{path}.hatch.cross")})
    dasharray = _array(payload["strokeDasharray"], path=f"{path}.strokeDasharray")
    return RegionAnnotationStyle(
        stroke=payload["stroke"],
        stroke_width=payload["strokeWidth"],
        stroke_dasharray=tuple(dasharray),
        line_cap=payload["lineCap"],
        fill=_optional_string(payload["fill"], f"{path}.fill"),
        fill_opacity=payload["fillOpacity"],
        hatch=hatch,
        label_color=payload["labelColor"],
        label_font_size=payload["labelFontSize"],
        label_orientation=payload["labelOrientation"],
        label_position=payload["labelPosition"],
        label_offset=payload["labelOffset"],
    )


def _encode_annotations(value: object, *, resources: RequestResources) -> dict[str, Any]:
    if not isinstance(value, AnnotationOptions):
        raise CanonicalRequestEncodingError("diagramOptions.annotations must be AnnotationOptions.")
    sets = []
    for annotation_set in value.sets:
        items = []
        for annotation in annotation_set.annotations:
            target = annotation.target
            target_payload: dict[str, Any]
            if isinstance(target, CoordinateSpan):
                target_payload = {
                    "kind": "coordinateSpan",
                    "record": _encode_selector(target.record),
                    "start": target.start,
                    "end": target.end,
                    "coordinateSpace": target.coordinate_space,
                    # CoordinateSpan reads wraps_origin by truth (OV-632).
                    "wrapsOrigin": bool(target.wraps_origin),
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
            elif isinstance(target, FeatureIdentitySpan):
                target_payload = {
                    "kind": "featureIdentity",
                    "recordKey": target.record_key,
                    "biologicalFeatureId": target.biological_feature_id,
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
    schema: int,
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
            target: RegionTarget
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
                    wraps_origin=_boolean(target_payload["wrapsOrigin"], f"{item_path}.target.wrapsOrigin"),
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
                    _decode_feature_selector(raw, path=f"{item_path}.target.selectors[{index}]")
                    for index, raw in enumerate(raw_selectors)
                )
                target = FeatureSpan(
                    record=_decode_selector(target_payload["record"], path=f"{item_path}.target.record"),
                    selectors=selectors,
                    envelope=target_payload["envelope"],
                    circular_path=target_payload["circularPath"],
                )
            elif kind == "featureIdentity" and schema >= FEATURE_OVERRIDE_SCHEMA:
                _require_exact_fields(
                    target_payload,
                    path=f"{item_path}.target",
                    required={"kind", "recordKey", "biologicalFeatureId", "envelope", "circularPath"},
                )
                target = FeatureIdentitySpan(
                    record_key=target_payload["recordKey"],
                    biological_feature_id=target_payload["biologicalFeatureId"],
                    envelope=target_payload["envelope"],
                    circular_path=target_payload["circularPath"],
                )
            else:
                raise CanonicalRequestDecodingError(f"Unsupported annotation target kind at {item_path}: {kind!r}.")
            metadata = _object(item["metadata"], path=f"{item_path}.metadata")
            # A Web writer wrote a fraction at lane (OV-630): the model reads its integer part.
            annotations.append(
                RegionAnnotation(
                    id=_string(item["id"], f"{item_path}.id"), target=target,
                    label=_string(item["label"], f"{item_path}.label"), mark=item["mark"],
                    lane=item["lane"],
                    style=(_decode_annotation_style(item["style"], path=f"{item_path}.style") if item["style"] is not None else None),
                    legend_label=_optional_string(item["legendLabel"], f"{item_path}.legendLabel"), metadata=metadata,
                )
            )
        sets.append(
            AnnotationSet(
                id=_string(set_payload["id"], f"{set_path}.id"),
                annotations=tuple(annotations),
                default_style=_decode_annotation_style(set_payload["defaultStyle"], path=f"{set_path}.defaultStyle"),
                legend_label=_optional_string(set_payload["legendLabel"], f"{set_path}.legendLabel"),
            )
        )
    table = (
        _decode_table_ref(
            payload["table"],
            name="annotations-table",
            resource_paths=resource_paths,
            text_cells=True,
        )
        if payload["table"] is not None else None
    )
    table_file = (
        str(_decode_file_ref(payload["tableFile"], name="annotations-table-file", resource_paths=resource_paths))
        if payload["tableFile"] is not None else None
    )
    return AnnotationOptions(sets=tuple(sets), table=table, table_file=table_file)


def _decode_feature_selector(value: object, *, path: str) -> FeatureSelector:
    selector = _object(value, path=path, required={"key", "value"})
    return FeatureSelector(
        key=_optional_string(selector["key"], f"{path}.key"),
        value=_string(selector["value"], f"{path}.value"),
    )


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
    result: dict[str, Any] = {"kind": kind}
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
    # The names the track planner reads a renderer and a side as (OV-653).
    renderer_name = (
        _circular_renderer_name if isinstance(slot, CircularTrackSlot) else _linear_renderer_name
    )
    names = {
        "renderer": renderer_name(str(slot.renderer)),
        "side": None if slot.side is None else str(slot.side).strip().lower(),
    }
    for item in slot_fields:
        raw = names[item.name] if item.name in names else getattr(slot, item.name)
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
            # A spec string keeps `marks=a|b` as text; the Session holds the
            # mark list the Web reads (OV-1298). Read before the style
            # override below becomes its camelCase Session form.
            if isinstance(raw.get("marks"), str):
                raw["marks"] = list(annotation_track_params_from_mapping(raw).marks or ())
            if isinstance(raw.get("style_override"), RegionAnnotationStyle):
                raw["style_override"] = _encode_annotation_style(raw["style_override"])
        result[_camel(item.name)] = _json_value(raw, path=f"trackSlot.{_camel(item.name)}")
    if legacy_spacing is not None:
        assert isinstance(slot, CircularTrackSlot)
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


@overload
def _decode_track_slots(
    value: object,
    *,
    mode: Literal["circular"],
    schema: int,
    path: str,
) -> tuple[str | CircularTrackSlot, ...] | None: ...


@overload
def _decode_track_slots(
    value: object,
    *,
    mode: Literal["linear"],
    schema: int,
    path: str,
) -> tuple[str | LinearTrackSlot, ...] | None: ...


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
                    else ScalarSpec.parse(str(decoded_spacing))
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
        # Wire values; _validate_dataclass_contract checks the constructed slot.
        kwargs: dict[str, Any] = {
            name: (
                _decode_scalar_if_needed(slot[key], path=f"{slot_path}.{key}")
                if name in scalar_fields
                else slot[key]
            )
            for key, name in field_map.items()
        }
        if mode == "circular" and schema in {1, 2}:
            params = dict(kwargs.get("params") or {})
            private_spacing = params.pop(_LEGACY_SPACING_PARAM, None)
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
        _validate_dataclass_contract(decoded_slot, path=slot_path, error="decode")
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
    resources: RequestResources,
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


def _record_frame(frames: tuple[tuple[int, bool], ...], index: int) -> tuple[int, bool]:
    return frames[index] if 0 <= index < len(frames) else (0, False)


def _encode_comparisons(
    options: CircularDiagramOptions | LinearDiagramOptions,
    *,
    mode: Literal["circular", "linear"],
    resources: RequestResources,
    frames: tuple[tuple[int, bool], ...],
) -> list[dict[str, Any]]:
    if mode == "circular":
        return []
    if not isinstance(options, LinearDiagramOptions):
        raise CanonicalRequestEncodingError(
            "Linear comparisons require LinearDiagramOptions."
        )
    result: list[dict[str, Any]] = []
    for index, comparison in enumerate(options.linear_comparisons or (), start=1):
        endpoint_frames = (
            _record_frame(frames, comparison.query_record_index),
            _record_frame(frames, comparison.subject_record_index),
        )
        if comparison.search_frame_text is not None:
            # A planner-resolved LOSAT edge or a decoded nucleotideBlast item
            # keeps its raw search-frame rows, as the Web writes them (design
            # 3.7); only a reversed persisted endpoint rewrites them.
            content = comparison.search_frame_text
            if endpoint_frames[0][1] or endpoint_frames[1][1]:
                content = reverse_endpoint_table_text(content, *endpoint_frames)
            resource_id = resources.add_bytes(
                _resource_id(f"comparison-losat-{index}"),
                kind="nucleotide-blast",
                name=f"comparison-losat-{index}.tsv",
                content=content.encode("utf-8"),
            )
            result.append(
                {
                    "kind": "nucleotideBlast",
                    "resourceId": resource_id,
                    "queryRecordIndex": comparison.query_record_index,
                    "subjectRecordIndex": comparison.subject_record_index,
                }
            )
            continue
        if comparison.drawn_feature_ids:
            raise CanonicalRequestEncodingError(
                f"Linear comparison {index} still holds the drawn-frame feature IDs of a "
                "Session up to 40; write it with gbdraw.session.save_session_document, "
                "which moves it to the search frame (OV-1183)."
            )
        comparison = reverse_unbound_endpoint_rows(comparison, *endpoint_frames)
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
        resource_id = _nucleotide_blast_resource(
            f"comparison-nucleotide-{index}",
            path,
            frames=(_record_frame(frames, index - 1), _record_frame(frames, index)),
            resources=resources,
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
    search = options.losat_search
    pairs = search.pairs if search is not None else None
    protein_mode = losatp_analysis_mode(search)
    settings = _encode_protein_settings(options)
    if (
        protein_mode != "none"
        or pairs is not None
        or settings != _encode_protein_settings(_DEFAULT_LINEAR_OPTIONS)
    ):
        result.append(
            {
                "kind": "generatedProteinComparison",
                "mode": protein_mode,
                "pairs": [
                    {
                        "queryRecordIndex": int(pair[0]),
                        "subjectRecordIndex": int(pair[1]),
                    }
                    for pair in (pairs or ())
                ],
                "settings": settings,
            }
        )
    return result


def _encode_protein_settings(options: LinearDiagramOptions) -> dict[str, Any]:
    search = options.losat_search or _UNSET_LOSAT_SEARCH
    owners = {"options": options, "search": search, "runtime": search.runtime}
    settings: dict[str, Any] = {}
    for key, owner, name in _PROTEIN_SETTINGS_WIRE:
        value = getattr(owners[owner], name)
        if key == "losatpBin" and value is None:
            value = AUTOMATIC_LOSAT_BIN
        settings[key] = _encode_pipeline_value(key, name, value)
    return settings


def _nucleotide_blast_resource(
    resource_id: str,
    path: object,
    *,
    frames: tuple[tuple[int, bool], tuple[int, bool]],
    resources: RequestResources,
) -> str:
    """Persist one ``-b`` table, rewritten for a reverse-complemented record.

    The table is in the search frame of the source records; a record persisted
    as its reverse-complemented sequence needs its rows in that sequence.
    """

    if not (frames[0][1] or frames[1][1]):
        return resources.add_path(resource_id, kind="nucleotide-blast", value=path)
    source = Path(str(path))
    try:
        text = source.read_text(encoding="utf-8-sig")
    except (OSError, UnicodeError) as exc:
        raise CanonicalRequestEncodingError(
            f"Could not read comparison table {source} for the Session."
        ) from exc
    content = reverse_endpoint_table_text(text, frames[0], frames[1])
    return resources.add_bytes(
        resource_id, kind="nucleotide-blast", name=source.name, content=content.encode("utf-8")
    )


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
                query_index = _non_negative_index(
                    item["queryRecordIndex"], f"{path}.queryRecordIndex"
                )
                subject_index = _non_negative_index(
                    item["subjectRecordIndex"], f"{path}.subjectRecordIndex"
                )
                from gbdraw.api.prepared import (
                    get_or_build_decoded_resource,
                    register_prepared_resource_value,
                )

                # The resource text stays on the edge, so a re-saved Session
                # writes the same nucleotideBlast resource. Both are read once
                # per Worker resource.
                def read_edge(
                    resource_path: Path = resource_path,
                    label: str = f"BLAST resource {item['resourceId']!r}",
                ) -> tuple[str, DataFrame]:
                    text = resource_path.read_bytes().decode("utf-8")
                    return text, read_comparison_table(StringIO(text), label=label)

                try:
                    text, table = get_or_build_decoded_resource(
                        resource_path, ("nucleotide-blast-edge-v1",), read_edge
                    )
                except (OSError, UnicodeError, ValidationError) as exc:
                    # The reader's diagnostic plus the record pair the Web names (CI-04).
                    cause = getattr(exc, "diagnostic", None) or (
                        None if isinstance(exc, ValidationError) else {"code": "INPUT_UNREADABLE", "field": "comparison"}
                    )
                    raise CanonicalRequestDecodingError(
                        f"Could not decode {path}: {exc}",
                        diagnostic=cause and {**cause, "queryRecordIndex": query_index, "subjectRecordIndex": subject_index},
                    ) from exc
                explicit_comparisons.append(
                    LinearComparison(
                        query_index,
                        subject_index,
                        register_prepared_resource_value(table, resource_path),
                        search_frame_text=text,
                    )
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
            protein_path = _resolve_resource(
                item["resourceId"],
                path=f"{path}.resourceId",
                resource_paths=resource_paths,
            )
            from gbdraw.api.prepared import get_or_build_decoded_resource

            # Read once per Worker resource, like the nucleotide edges above.
            table = get_or_build_decoded_resource(
                protein_path,
                ("canonical-protein-comparison-v1",),
                lambda: _read_canonical_table(protein_path, context=path),
            )
            if schema >= 2:
                query_index = _non_negative_index(
                    item["queryRecordIndex"], f"{path}.queryRecordIndex"
                )
                subject_index = _non_negative_index(
                    item["subjectRecordIndex"], f"{path}.subjectRecordIndex"
                )
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
                expected=OrthogroupResult | OrthogroupGraphResult,
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
            pipeline = _decode_pipeline(item, path=path, schema=schema)
            legacy_alignment = pipeline.pop("_legacy_similarity_alignment", None)
            if legacy_alignment is not None:
                result["_legacy_similarity_alignment"] = legacy_alignment
            result.update(pipeline)
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


def _encode_pipeline_value(key: str, name: str, value: object) -> Any:
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
    return _json_value(value, path=f"comparisons.settings.{key}")


def _decode_pipeline(
    item: Mapping[str, Any], *, path: str, schema: int
) -> dict[str, Any]:
    required = {"kind", "mode", "settings"}
    if schema >= 2:
        required.add("pairs")
    _require_exact_fields(item, path=path, required=required)
    settings = _object(item["settings"], path=f"{path}.settings")
    setting_keys = {key for key, _owner, _name in _PROTEIN_SETTINGS_WIRE}
    if schema < 8:
        setting_keys.add(_LEGACY_ALIGNMENT_SETTING)
    _require_exact_fields(
        settings, path=f"{path}.settings",
        required=setting_keys - _OPTIONAL_PROTEIN_SETTINGS,
        optional=set(_OPTIONAL_PROTEIN_SETTINGS),
    )
    mode = normalize_protein_blastp_mode(_string(item["mode"], f"{path}.mode"))
    result: dict[str, Any] = {}
    search_values: dict[str, Any] = {}
    runtime_values: dict[str, Any] = {}
    decoded_pairs: list[tuple[int, int]] = []
    legacy_max_paralog_links: int | None = None
    if schema >= 2:
        pairs = _array(item["pairs"], path=f"{path}.pairs")
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
    search_pairs = tuple(decoded_pairs) if decoded_pairs and mode == "pairwise" else None
    owned_settings = [(key, owner, name) for key, owner, name in _PROTEIN_SETTINGS_WIRE]
    if schema < 8:
        owned_settings.append((_LEGACY_ALIGNMENT_SETTING, "legacy", "align_orthogroup_feature"))
    for key, owner, name in owned_settings:
        raw = settings[key] if key in settings else _absent_protein_setting(key)
        if owner == "search":
            search_values[name] = raw
            continue
        if owner == "runtime":
            runtime_values[name] = raw
            continue
        if name == "align_orthogroup_feature":
            if raw is None:
                continue
            if not isinstance(raw, str) or not raw.strip() or "\0" in raw:
                raise CanonicalRequestDecodingError(
                    f"{path}.settings.{key} must be null or non-empty text without NUL."
                )
            result["_legacy_similarity_alignment"] = _LegacySimilarityAlignment(
                target=raw,
                source_schema=schema,
            )
            continue
        if name == "collinearity_anchor_mode":
            raw = _read_retired_anchor_mode(raw, path=f"{path}.settings.{key}")
        names = _setting_names(name)
        if names and schema >= _SETTING_NAMES_SCHEMA and raw not in names:
            raise CanonicalRequestDecodingError(
                f"{path}.settings.{key} must be one of: {', '.join(names)}.",
                diagnostic=_REQUEST_FIELDS_INVALID,
            )
        if name == "collinearity_params":
            decoded, legacy_max_paralog_links = _decode_collinearity_params(
                raw,
                path=f"{path}.settings.{key}",
                allow_standard=schema in {1, 2},
            )
            result[name] = decoded
        else:
            result[name] = raw
    search = LosatSearchOptions(
        program="losatp",
        losatp_mode=_WIRE_LOSATP_MODES[mode],
        pairs=search_pairs,
        runtime=LosatRuntimeOptions(**runtime_values),
        **search_values,
    )
    result["losat_search"] = None if search == _UNSET_LOSAT_SEARCH else search
    if (
        mode == "collinear"
        and legacy_max_paralog_links is not None
        and result["collinear_max_paralog_links_per_orthogroup"]
        == _LEGACY_PARALOG_LINKS_DEFAULT
    ):
        result["collinear_max_paralog_links_per_orthogroup"] = (
            legacy_max_paralog_links
        )
    return result


def _absent_protein_setting(key: str) -> object:
    """The value of an optional protein setting (``_OPTIONAL_PROTEIN_SETTINGS``) that a request does not hold."""

    return _ABSENT_PROTEIN_SETTINGS[key]


def _read_retired_anchor_mode(raw: object, *, path: str) -> object:
    """A retired anchor mode, read as its replacement with a warning."""

    if not isinstance(raw, str):
        return raw
    replacement = RETIRED_COLLINEARITY_ANCHOR_MODES.get(raw.strip().lower().replace("-", "_"))
    if replacement is None:
        return raw
    logger.warning(
        "%s: the Collinear anchor mode %r was removed; it is read as %r.",
        path,
        raw,
        replacement,
    )
    return replacement


def _setting_names(name: str) -> tuple[str, ...]:
    """The names a LinearDiagramOptions setting stores: the Literal part of its type."""

    return tuple(
        value
        for branch in get_args(_class_type_hints(LinearDiagramOptions)[name]) or ()
        if get_origin(branch) is Literal
        for value in get_args(branch)
    )


def _decode_collinearity_params(
    value: object,
    *,
    path: str,
    allow_standard: bool,
) -> tuple[LosslessCollinearityParameters | None, int | None]:
    if value is None:
        return None, None
    payload = _object(value, path=path, required={"kind", "parameters"})
    kind = payload["kind"]
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
    _require_exact_fields(
        parameters, path=f"{path}.parameters", required=set(field_map)
    )
    result = cls(**{name: parameters[key] for key, name in field_map.items()})
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


def encode_canonical_table(frame: object, *, name: str = "table") -> str:
    """Encode one DataFrame as the canonical request's TSV table text.

    This is the only encoder of canonical TSV tables: request resources and the
    Web protein comparison payload both use it, so equal frames give equal bytes.
    """

    if not isinstance(frame, DataFrame):
        raise CanonicalRequestEncodingError(f"{name} must be a pandas DataFrame.")
    if len(frame.columns) == 0:
        raise CanonicalRequestEncodingError(
            f"DataFrame resource {name!r} must have at least one column."
        )
    try:
        return frame.to_csv(sep="\t", index=False, lineterminator="\n")
    except Exception as exc:
        raise CanonicalRequestEncodingError(
            f"Could not encode DataFrame resource {name!r} as canonical TSV."
        ) from exc


def _table_ref(
    name: str,
    value: object,
    *,
    resources: RequestResources,
) -> dict[str, str]:
    text = encode_canonical_table(value, name=name)
    preferred_id = _resource_id(name)
    resource_id = resources.add_bytes(
        preferred_id,
        kind="canonical-tsv",
        name=f"{preferred_id}.tsv",
        content=text.encode("utf-8"),
    )
    return {"resourceId": resource_id, "representation": "canonicalTsv"}


def _file_ref(
    name: str,
    value: object,
    *,
    resources: RequestResources,
) -> dict[str, str]:
    preferred_id = _resource_id(name)
    resource_id = resources.add_path(preferred_id, kind=preferred_id, value=value)
    return {"resourceId": resource_id, "representation": "file"}


def _decode_table_ref(
    value: object,
    *,
    name: str,
    resource_paths: Mapping[str, str | Path],
    text_cells: bool = False,
) -> DataFrame:
    ref = _resource_ref(value, path=f"resource reference {name}", representation="canonicalTsv")
    return _read_canonical_table(
        _resolve_resource(
            ref["resourceId"],
            path=f"resource reference {name}.resourceId",
            resource_paths=resource_paths,
        ),
        context=name,
        text_cells=text_cells,
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


def _read_canonical_table(
    file_path: Path, *, context: str, text_cells: bool = False
) -> DataFrame:
    # A text table keeps every cell as written (an empty label text stays "");
    # the numeric tables (depth, comparisons) take pandas' inferred dtypes and
    # read floats exactly, so an encoded table survives a read unchanged.
    text_options: dict[str, Any] = (
        {"dtype": str, "keep_default_na": False}
        if text_cells
        else EXACT_FLOAT_READ_OPTIONS
    )
    try:
        table = read_csv(file_path, sep="\t", **text_options)
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
    resources: RequestResources,
) -> str:
    content = encode_canonical_typed_resource(value_kind, value)
    return resources.add_bytes(
        resource_id,
        kind=kind,
        name=f"{resource_id}.json",
        content=content,
    )


def encode_canonical_typed_resource(value_kind: str, value: object) -> bytes:
    """Encode one value using the canonical request's typed JSON resource contract."""

    if not isinstance(value_kind, str) or not value_kind.strip():
        raise CanonicalRequestEncodingError(
            "Canonical typed resources require a non-empty value kind."
        )
    body = {
        "schema": 3,
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


def decode_canonical_typed_resource(
    content: bytes,
    *,
    value_kind: str,
    expected: object,
    path: str = "typed resource",
) -> Any:
    """Decode one canonical typed JSON resource from its published bytes."""

    if not isinstance(content, bytes):
        raise CanonicalRequestDecodingError(
            f"Canonical JSON resource for {path} must be bytes."
        )
    try:
        raw = json.loads(content.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise CanonicalRequestDecodingError(
            f"Could not decode canonical JSON resource for {path}."
        ) from exc
    payload = _object(
        raw, path=f"{path} resource", required={"schema", "kind", "value"}
    )
    resource_schema = payload["schema"]
    # `True == 1` and `1.0 == 1`: only an integer names a schema (OV-1930).
    if type(resource_schema) is not int or resource_schema not in _TYPED_RESOURCE_SCHEMAS or payload["kind"] != value_kind:
        raise CanonicalRequestDecodingError(
            f"Canonical JSON resource metadata does not match {path}."
        )
    return _decode_typed_tree(
        payload["value"],
        expected,
        path=f"{path}.value",
        resource_schema=resource_schema,
    )


def describe_canonical_typed_resources(resources: Mapping[str, tuple[str, type]]) -> dict[str, Any]:
    """Describe the JSON that :func:`decode_canonical_typed_resource` decodes,
    for each ``key: (value_kind, expected)`` of ``resources``.

    Returns ``{"roots": {key: name}, "types": {name: shape}}`` with the shapes
    of :mod:`gbdraw.contract_shapes`: per key the resource object, one option
    per resource schema whose tree decoder can accept ``expected``, and each
    dataclass as the ``{"type", "fields"}`` object its decoder requires,
    named ``<class>Schema<n>`` after the schema it is first described for. The description is structural: a value inside
    it is still refused where a class checks its values (``__post_init__``)
    or where an ``OrthologPathCollection`` reads its own fields
    (``from_payload``); the shape of a path collection holds only the tagged
    object its decoder requires before that read.
    """

    from .contract_shapes import (
        BOOLEAN, INTEGER, NULL, NUMBER, STRING, UNKNOWN, Shape, array, literal, obj, record, ref, union,
    )

    scalars: dict[object, Shape] = {str: STRING, int: INTEGER, float: NUMBER, bool: BOOLEAN, type(None): NULL}
    named: dict[str, Shape] = {}

    def tree(hint: object, schema: int) -> Shape | None:
        """The shape of a value of ``hint`` under ``schema``; ``None`` where
        its decoder refuses every value (as :func:`_build_typed_tree_decoder`)."""

        if hint is OrthologPathCollection:
            if schema != 3:
                return None
            return obj({"type": literal("OrthologPathCollection"), "fields": record(UNKNOWN)})
        if hint is Any or hint is object:
            return UNKNOWN
        origin = get_origin(hint)
        args = get_args(hint)
        if origin in {types.UnionType, Union}:
            options = [shape for shape in (tree(branch, schema) for branch in args) if shape is not None]
            return union(*options) if options else None
        if (hint is OrthogroupGraphResult and schema != 3) or (hint is OrthogroupResult and schema == 3):
            return None
        if isinstance(hint, type) and is_dataclass(hint):
            return dataclass_shape(hint, schema)
        if origin in {tuple, list, Sequence, SequenceABC} or hint in {tuple, list}:
            return array(required_tree(args[0] if args else Any, schema, hint))
        if origin in {dict, Mapping, MappingABC} or hint is dict:
            key_hint, value_hint = args if len(args) == 2 else (str, Any)
            return record(required_tree(value_hint, schema, hint)) if key_hint is str else None
        if origin is Literal:
            return literal(*args)
        if hint in scalars:
            return scalars[hint]
        raise TypeError(f"No typed resource shape is defined for the type {hint!r}.")

    def required_tree(hint: object, schema: int, owner: object) -> Shape:
        shape = tree(hint, schema)
        if shape is None:
            raise TypeError(f"An item of {owner!r} has no value its decoder accepts.")
        return shape

    def dataclass_shape(hint: type, schema: int) -> Shape | None:
        name = f"{hint.__name__}Schema{schema}"
        if name in named:
            return ref(name)
        if hint.__name__ not in _TYPED_TREE_CLASSES:
            return None
        optional = _SCHEMA1_TYPED_OPTIONAL_FIELDS.get(hint.__name__, frozenset()) if schema == 1 else frozenset()
        hints = _class_type_hints(hint)
        named[name] = UNKNOWN  # a class that holds itself refers to its name
        fields_: dict[str, Shape] = {}
        for key, item in _class_field_defs(hint).items():
            shape = tree(hints[item.name], schema)
            if shape is not None:
                fields_[key] = shape
            elif key not in optional:
                del named[name]
                return None
        shape = obj({
            "type": literal(hint.__name__),
            "fields": obj(fields_, optional=sorted(optional.intersection(fields_))),
        })
        # An equal shape described for another schema keeps its name, so each shape is written once.
        same = next(
            (f"{hint.__name__}Schema{n}" for n in _TYPED_RESOURCE_SCHEMAS
             if n != schema and named.get(f"{hint.__name__}Schema{n}") == shape),
            None,
        )
        if same is not None:
            del named[name]
            return ref(same)
        named[name] = shape
        return ref(name)

    roots: dict[str, str] = {}
    for key, (value_kind, expected) in resources.items():
        options = [
            obj({"schema": literal(schema), "kind": literal(value_kind), "value": value})
            for schema in _TYPED_RESOURCE_SCHEMAS
            if (value := tree(expected, schema)) is not None
        ]
        if not options:
            raise TypeError(f"No typed resource schema decodes {expected!r}.")
        roots[key] = f"{expected.__name__}Resource"
        named[roots[key]] = union(*options)
    return {"roots": roots, "types": named}


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
            content = resource_path.read_bytes()
        except OSError as exc:
            raise CanonicalRequestDecodingError(
                f"Could not decode canonical JSON resource for {path}."
            ) from exc
        return decode_canonical_typed_resource(
            content,
            value_kind=value_kind,
            expected=expected,
            path=path,
        )

    from gbdraw.api.prepared import get_or_build_decoded_resource

    return get_or_build_decoded_resource(
        resource_path,
        ("canonical-typed-json", value_kind, repr(expected)),
        decode,
    )


def _encode_typed_tree(value: object) -> Any:
    value_type = type(value)
    # Exact built-in types first: none of them is a result, a path collection
    # or a dataclass, so they take the branches below without the checks.
    if value_type in _TYPED_TREE_JSON_SCALARS:
        return _json_value(value, path="typed resource")
    if value_type is list or value_type is tuple:
        return [_encode_typed_tree(item) for item in cast("Sequence[object]", value)]
    if isinstance(value, OrthogroupResult):
        value = compact_ortholog_paths(value)
    if isinstance(value, OrthologPathCollection):
        body = value.to_payload()
        if value.kind == "explicit":
            body["paths"] = _encode_typed_tree(body["paths"])
        return {"type": "OrthologPathCollection", "fields": body}
    if is_dataclass(value) and not isinstance(value, type):
        cls = type(value)
        if cls.__name__ not in _TYPED_TREE_CLASSES:
            raise CanonicalRequestEncodingError(
                f"Unsupported typed resource value: {cls.__name__}."
            )
        return {
            "type": cls.__name__,
            "fields": {
                key: _encode_typed_tree(getattr(value, item.name))
                for key, item in _class_field_defs(cls).items()
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


_TYPED_TREE_JSON_SCALARS = frozenset({str, int, float, bool, type(None)})


def _decode_typed_tree(
    value: object,
    hint: object,
    *,
    path: str,
    resource_schema: int,
) -> Any:
    return _typed_tree_decoder(hint, resource_schema).decode(value, path)


class _TreeDecoder(NamedTuple):
    """The decoder of one type hint. ``accepts(value)`` is true only where
    ``decode(value, path)`` returns ``value`` itself, so a container skips the
    call and the item's path for such a value."""

    decode: Callable[[object, str], object]
    accepts: Callable[[object], bool]


def _typed_tree_decoder(hint: object, resource_schema: int) -> _TreeDecoder:
    # ``A | B == B | A``, but a union tries its branches in order, so the key
    # also holds the repr, which keeps the order at every depth.
    return _build_typed_tree_decoder(hint, repr(hint), resource_schema)


@lru_cache(maxsize=1024)
def _build_typed_tree_decoder(hint: object, hint_repr: str, resource_schema: int) -> _TreeDecoder:
    """Build once per type hint and resource schema the decoder of a typed tree.

    A dataclass reads its field hints at its first decode, as when each node
    interpreted its hint, so a hint that cannot be read fails at the same value.
    """

    if hint is OrthologPathCollection:
        return _TreeDecoder(_path_collection_decoder(resource_schema), _reject)
    if hint is Any or hint is object:
        return _TreeDecoder(lambda value, path: value, _accept_any)
    origin = get_origin(hint)
    args = get_args(hint)
    if origin in {types.UnionType, Union}:
        return _union_decoder(args, resource_schema)
    if hint is OrthogroupGraphResult and resource_schema != 3:
        return _TreeDecoder(
            _raise_decoding("Graph orthogroups require typed resource schema 3."), _reject
        )
    if hint is OrthogroupResult and resource_schema == 3:
        return _TreeDecoder(
            _raise_decoding("Legacy orthogroups require typed resource schema 1 or 2."), _reject
        )
    if isinstance(hint, type) and is_dataclass(hint):
        return _TreeDecoder(_dataclass_decoder(hint, resource_schema), _reject)
    if origin in {tuple, list, Sequence, SequenceABC}:
        item = _typed_tree_decoder(args[0] if args else Any, resource_schema)
        return _TreeDecoder(_sequence_decoder(item, as_tuple=origin is tuple), _reject)
    if origin in {dict, Mapping, MappingABC}:
        key_hint, value_hint = args if len(args) == 2 else (str, Any)
        return _TreeDecoder(
            _mapping_decoder(
                _typed_tree_decoder(value_hint, resource_schema) if key_hint is str else None
            ),
            _reject,
        )
    if origin is Literal:

        def accepts_literal(value: object) -> bool:
            return value in args

        def decode_literal(value: object, path: str) -> object:
            if value not in args:
                raise CanonicalRequestDecodingError(
                    f"Typed resource literal at {path} has unsupported value {value!r}."
                )
            return value

        return _TreeDecoder(decode_literal, accepts_literal)
    matches = _type_checker(hint)

    def decode_value(value: object, path: str) -> object:
        if not matches(value):
            raise CanonicalRequestDecodingError(
                f"Typed resource value at {path} has an invalid type."
            )
        return value

    return _TreeDecoder(decode_value, matches)


def _raise_decoding(message: str) -> Callable[[object, str], NoReturn]:
    def decode(value: object, path: str) -> NoReturn:
        raise CanonicalRequestDecodingError(message)

    return decode


def _path_collection_decoder(resource_schema: int) -> Callable[[object, str], object]:
    def decode(value: object, path: str) -> OrthologPathCollection:
        tagged = _object(value, path=path, required={"type", "fields"})
        if resource_schema != 3 or tagged["type"] != "OrthologPathCollection":
            raise CanonicalRequestDecodingError(f"Invalid path collection at {path}.")
        body = dict(_object(tagged["fields"], path=f"{path}.fields"))
        if body.get("kind") == "explicit" and "paths" in body:
            body["paths"] = _decode_typed_tree(body["paths"], tuple[OrthologPath, ...],
                                             path=f"{path}.paths", resource_schema=resource_schema)
        try:
            return OrthologPathCollection.from_payload(body)
        except (ValidationError, ValueError, TypeError, KeyError) as exc:
            raise CanonicalRequestDecodingError(f"Invalid path collection at {path}: {exc}") from exc

    return decode


def _union_decoder(args: tuple[object, ...], resource_schema: int) -> _TreeDecoder:
    has_none = type(None) in args
    branches = tuple(
        (
            branch.__name__ if isinstance(branch, type) else None,
            get_origin(branch) in {tuple, list, Sequence, SequenceABC},
            _type_checker(branch),
            _typed_tree_decoder(branch, resource_schema),
        )
        for branch in args
        if branch is not type(None)
    )

    def decode_union(value: object, path: str) -> object:
        if value is None and has_none:
            return None
        tagged_type = value.get("type") if isinstance(value, MappingABC) else None
        for name, sequence, matches, branch in branches:
            if (
                (name is not None and name == tagged_type)
                or (sequence and isinstance(value, list))
                or matches(value)
            ):
                return branch.decode(value, path)
        raise CanonicalRequestDecodingError(
            f"Typed resource value at {path} does not match its union contract."
        )

    if any(branch.accepts is _reject for *_rest, branch in branches):
        return _TreeDecoder(decode_union, _reject)
    branch_checks = tuple(matches for _name, _sequence, matches, _branch in branches)

    def accepts_union(value: object) -> bool:
        # Every branch returns the value it matches; a mapping may name a branch.
        if value is None and has_none:
            return True
        return not isinstance(value, MappingABC) and any(
            matches(value) for matches in branch_checks
        )

    return _TreeDecoder(decode_union, accepts_union)


def _sequence_decoder(item: _TreeDecoder, *, as_tuple: bool) -> Callable[[object, str], object]:
    decode_item, accepts = item

    def decode(value: object, path: str) -> list[object] | tuple[object, ...]:
        raw = _array(value, path=path)
        decoded = [
            entry if accepts(entry) else decode_item(entry, f"{path}[{index}]")
            for index, entry in enumerate(raw)
        ]
        return tuple(decoded) if as_tuple else decoded

    return decode


def _mapping_decoder(item: _TreeDecoder | None) -> Callable[[object, str], object]:
    def decode(value: object, path: str) -> dict[str, object]:
        raw_mapping = _object(value, path=path)
        if item is None:
            raise CanonicalRequestDecodingError(
                f"Unsupported typed resource mapping key contract at {path}."
            )
        decode_item, accepts = item
        return {
            key: entry if accepts(entry) else decode_item(entry, f"{path}.{key}")
            for key, entry in raw_mapping.items()
        }

    return decode


class _DataclassDecodePlan(NamedTuple):
    required: set[str]
    optional: set[str]
    entries: tuple[tuple[str, Field[object], _TreeDecoder], ...]


def _dataclass_decoder(hint: type, resource_schema: int) -> Callable[[object, str], object]:
    plans: list[_DataclassDecodePlan] = []

    def plan() -> _DataclassDecodePlan:
        if not plans:
            hints = _class_type_hints(hint)
            field_defs = _class_field_defs(hint)
            optional_fields = (
                _SCHEMA1_TYPED_OPTIONAL_FIELDS.get(hint.__name__, frozenset())
                if resource_schema == 1
                else frozenset()
            )
            plans.append(
                _DataclassDecodePlan(
                    required=set(field_defs).difference(optional_fields),
                    optional=set(optional_fields),
                    entries=tuple(
                        (key, item, _typed_tree_decoder(hints[item.name], resource_schema))
                        for key, item in field_defs.items()
                    ),
                )
            )
        return plans[0]

    def decode(value: object, path: str) -> object:
        tagged = _object(value, path=path, required={"type", "fields"})
        if tagged["type"] != hint.__name__ or hint.__name__ not in _TYPED_TREE_CLASSES:
            raise CanonicalRequestDecodingError(
                f"Unsupported typed resource object at {path}: {tagged['type']!r}."
            )
        raw_fields = _object(tagged["fields"], path=f"{path}.fields")
        required, optional, field_plans = plan()
        _require_exact_fields(
            raw_fields,
            path=f"{path}.fields",
            required=required,
            optional=optional,
        )
        kwargs: dict[str, object] = {}
        for key, item, (decode_field, accepts) in field_plans:
            if key in raw_fields:
                entry = raw_fields[key]
                kwargs[item.name] = (
                    entry if accepts(entry) else decode_field(entry, f"{path}.fields.{key}")
                )
            elif item.default is not MISSING:
                kwargs[item.name] = item.default
            elif item.default_factory is not MISSING:
                kwargs[item.name] = item.default_factory()
            else:  # pragma: no cover - guarded by the required-field check
                raise CanonicalRequestDecodingError(
                    f"Missing required field at {path}.fields.{key}."
                )
        try:
            result = hint(**kwargs)
        except ValidationError as exc:
            raise CanonicalRequestDecodingError(f"Invalid typed resource at {path}: {exc}") from exc
        return compact_ortholog_paths(result) if isinstance(result, OrthogroupResult) else result

    return decode


def _resolve_resource(
    resource_id: object,
    *,
    path: str,
    resource_paths: Mapping[str, str | Path],
) -> Path:
    if not isinstance(resource_id, str) or not SESSION_ID_RE.fullmatch(resource_id):
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
    hints = _class_type_hints(type(value))
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
    return _type_checker(hint)(value)


def _accept_any(value: object) -> bool:
    return True


def _reject(value: object) -> bool:
    return False


def _is_finite_float(value: object) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
    )


def _is_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _is_bool(value: object) -> bool:
    return isinstance(value, bool)


def _is_str(value: object) -> bool:
    return isinstance(value, str)


def _is_none(value: object) -> bool:
    return value is None


def _is_sequence(value: object) -> bool:
    return isinstance(value, SequenceABC) and not isinstance(value, (str, bytes, bytearray))


def _is_mapping(value: object) -> bool:
    return isinstance(value, MappingABC)


_SCALAR_CHECKERS: dict[object, Callable[[object], bool]] = {
    float: _is_finite_float,
    int: _is_int,
    bool: _is_bool,
    str: _is_str,
    type(None): _is_none,
}


@lru_cache(maxsize=1024)
def _type_checker(hint: object) -> Callable[[object], bool]:
    """Build once per type hint the check that a value matches its typed contract."""
    if hint is Any or hint is object:
        return _accept_any
    origin = get_origin(hint)
    args = get_args(hint)
    if origin in {types.UnionType, Union}:
        branches = tuple(_type_checker(branch) for branch in args)
        return lambda value: any(branch(value) for branch in branches)
    if origin is Literal:
        return lambda value: value in args and any(
            type(value) is type(item) for item in args if item == value
        )
    if origin in {tuple, list, Sequence, SequenceABC}:
        if not args:
            return _is_sequence
        item_check = _type_checker(args[0])
        return lambda value: _is_sequence(value) and all(
            item_check(item) for item in cast(SequenceABC[object], value)
        )
    if origin in {dict, Mapping, MappingABC}:
        if len(args) != 2:
            return _is_mapping
        key_check = _type_checker(args[0])
        item_check = _type_checker(args[1])
        return lambda value: _is_mapping(value) and all(
            key_check(key) and item_check(item)
            for key, item in cast(MappingABC[object, object], value).items()
        )
    scalar = _SCALAR_CHECKERS.get(hint) if isinstance(hint, type) else None
    if scalar is not None:
        return scalar
    if isinstance(hint, type):
        cls = hint
        return lambda value: isinstance(value, cls)
    return _reject


def _raise_contract_error(kind: str, message: str) -> NoReturn:
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


@lru_cache(maxsize=256)
def _class_type_hints(cls: type) -> dict[str, Any]:
    return get_type_hints(cls)


@lru_cache(maxsize=256)
def _class_field_defs(cls: type) -> dict[str, Field[Any]]:
    return {_camel(item.name): item for item in fields(cls)}


@lru_cache(maxsize=1024)
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
        raise CanonicalRequestDecodingError(f"{path} must be a boolean.", diagnostic=_REQUEST_FIELDS_INVALID)
    return value


def _string(value: object, path: str) -> str:
    if not isinstance(value, str):
        raise CanonicalRequestDecodingError(f"{path} must be a string.", diagnostic=_REQUEST_FIELDS_INVALID)
    return value


def _optional_string(value: object, path: str) -> str | None:
    if value is not None and not isinstance(value, str):
        raise CanonicalRequestDecodingError(f"{path} must be a string or null.", diagnostic=_REQUEST_FIELDS_INVALID)
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


@lru_cache(maxsize=1)
def _current_request_description() -> dict[str, Any]:
    return describe_canonical_request()


def _undescribed_parts(payload: object) -> list[str]:
    """Where a current-schema request leaves its description (BA-D09); older
    schemas are not described."""

    from .contract_shapes import ref, violations

    schema = payload.get("schema") if isinstance(payload, MappingABC) else None
    if type(schema) is not int or schema != CANONICAL_REQUEST_SCHEMA:
        return []
    description = _current_request_description()
    return violations(
        _as_decoded_settings(payload), ref(description["root"]), description["types"], "renderRequest"
    )


def _as_decoded_settings(payload: object) -> object:
    """``payload`` with the retired anchor mode "one_to_one", which a released
    schema 9 writer stored, read as "rbh" (``_read_retired_anchor_mode``); the
    Web repair of the request reads it as the description's default, with a
    row."""

    if not isinstance(payload, MappingABC):
        return payload
    comparisons = payload.get("comparisons")
    if not isinstance(comparisons, list):
        return payload
    read: list[object] = []
    for item in comparisons:
        settings = item.get("settings") if isinstance(item, MappingABC) else None
        if not isinstance(settings, MappingABC):
            read.append(item)
            continue
        settings = dict(settings)
        mode = settings.get("collinearityAnchorMode")
        if isinstance(mode, str) and mode in RETIRED_COLLINEARITY_ANCHOR_MODES:
            settings["collinearityAnchorMode"] = RETIRED_COLLINEARITY_ANCHOR_MODES[mode]
        read.append({**item, "settings": settings})
    return {**payload, "comparisons": read}


def describe_legacy_request_tables() -> dict[str, Any]:
    """The tables of the request schemas 1-2 reading, with wire names.

    The Web lift of a stored schema 1 or 2 request into schema 5
    (``gbdraw/web/js/services/legacy-render-request.js``) reads them from the
    generated ``LEGACY_REQUEST_TABLES``; this module's tables stay their owner.
    """

    paralog_key = next(
        key
        for key, _owner, name in _PROTEIN_SETTINGS_WIRE
        if name == "collinear_max_paralog_links_per_orthogroup"
    )
    sparse_paths = {path for overrides in _LEGACY_SPARSE_CONFIG_OVERRIDES.values() for path in overrides}
    return {
        "comparisonDefaults": {
            _option_wire_key(name): value
            for name, value in _LEGACY_SPARSE_COMPARISON_DEFAULTS.items()
        },
        "featureTypes": list(LEGACY_FEATURE_TYPES),
        "featureShapes": dict(LEGACY_FEATURE_SHAPES),
        "multiRecordSizeModeRenames": dict(LEGACY_MULTI_RECORD_SIZE_MODE_RENAMES),
        # The slot parameter the decoder reads a slot spacing from (A8, N3).
        "privateSpacingParam": _LEGACY_SPACING_PARAM,
        "configOverrides": {
            mode: dict(overrides)
            for mode, overrides in _LEGACY_SPARSE_CONFIG_OVERRIDES.items()
        },
        "configOverrideRenames": dict(_LEGACY_CONFIG_OVERRIDE_KEYS),
        "linearConfigOverrideValues": {
            name: dict(values)
            for name, values in _LEGACY_LINEAR_CONFIG_OVERRIDE_VALUES.items()
        },
        # The canonical paths of the flat names the lift reads: the Linear
        # values above, a comparison height (B17, Web only), and the flat
        # names of the sparse config defaults (a default yields to them).
        "configOverridePaths": {
            name: list(paths)
            for name, paths in _LEGACY_FLAT_CONFIG_OVERRIDE_PATHS.items()
            if name in (*_LEGACY_LINEAR_CONFIG_OVERRIDE_VALUES, "comparison_height")
            or sparse_paths.intersection(paths)
        },
        # What replaces an invalid comparison height of a stored full config
        # (B17): the packaged default.
        "comparisonHeightDefault": load_config_toml("gbdraw.data", "config.toml")["canvas"]["linear"][
            "comparison_height"
        ],
        "diagramOptionRenames": dict(_LEGACY_FEATURE_VISIBILITY_RENAMES),
        "collinearityStandardParameters": [
            _camel(item.name) for item in fields(_LegacyStandardCollinearityPayload)
        ],
        "collinearityStandardFields": {
            _camel(standard): _camel(lossless)
            for standard, lossless in _LEGACY_STANDARD_COLLINEARITY_FIELDS.items()
        },
        "collinearityMergeOrientation": _LEGACY_STANDARD_MERGE_ORIENTATION,
        "collinearityStandardParalogLinks": _camel("max_paralog_links_per_orthogroup"),
        "paralogLinksSetting": paralog_key,
        "paralogLinksDefault": _LEGACY_PARALOG_LINKS_DEFAULT,
    }


def describe_canonical_request() -> dict[str, Any]:
    """Describe the JSON of a current-schema request as neutral data.

    Returns ``{"root": name, "types": {name: shape}}`` with the shapes of
    :mod:`gbdraw.contract_shapes`. Each object comes from the dataclass whose
    fields the encoder writes, and from the tables the encoder and the decoder
    route option fields with; an object with no dataclass names its encoder.
    Old schemas are not described.
    """

    from .contract_shapes import (
        INTEGER, NULL, NUMBER, STRING, UNKNOWN, Shape,
        array, bounded, literal, non_null, nullable, obj, record, ref, shape_of_hint, shape_violations, tuple_of,
        union,
    )
    from .mode_profiles import DINUCLEOTIDE_PATTERN
    from .annotations.models import BLANK_TEXT_PATTERN, COLOR_NAME_TEXT_PATTERN, HEX_COLOR_TEXT_PATTERN
    from .session_resources import SessionResourceTable

    types_: dict[str, Shape] = {}

    def named(name: str, shape: Shape) -> Shape:
        types_[f"Canonical{name}"] = shape
        return ref(f"Canonical{name}")

    def kinds(name: str, suffix: str, variants: Mapping[str, Mapping[str, Shape]]) -> Shape:
        """An object tagged by ``kind``; each variant is named ``<Kind><suffix>``."""

        return named(name, union(*(
            named(
                kind[0].upper() + kind[1:] + ("" if kind.endswith(suffix) else suffix),
                obj({"kind": literal(kind), **variant}),
            )
            for kind, variant in variants.items()
        )))

    # _table_ref and _file_ref.
    table_ref = named("TableRef", obj({"resourceId": STRING, "representation": literal("canonicalTsv")}))
    file_ref = named("FileRef", obj({"resourceId": STRING, "representation": literal("file")}))
    # _encode_selector.
    selector = kinds("RecordSelector", "Selector", {
        "recordId": {"value": STRING},
        "recordIndex": {"index": INTEGER},
    })

    # The shape of a class the encoder writes its own way, by class.
    wire: dict[object, Shape] = {
        DataFrame: table_ref,
        RecordSelector: selector,
        # asdict(GbdrawConfig) or a caller's mapping; the decoder migrates it.
        GbdrawConfig: record(UNKNOWN),
    }

    def of(hint: object) -> Shape:
        return shape_of_hint(hint, wire.get)

    def resource(hint: object) -> Shape:
        """A path or a table that the encoder writes as a resource reference."""

        return shape_of_hint(hint, lambda item: file_ref if item in {str, Path} else wire.get(item))

    def data(
        cls: type,
        name: str | None = None,
        *,
        files: frozenset[str] = frozenset(),
        replace: Mapping[str, Shape] | None = None,
        extra: Mapping[str, Shape] | None = None,
        kind: str | None = None,
        defaults: Mapping[str, object] | None = None,
    ) -> Shape:
        """The object the encoder writes for a dataclass: one camelCase key per field.

        ``files`` are the fields written as resource references, ``replace``
        the fields the encoder writes its own way, ``extra`` the keys it adds,
        ``kind`` the tag it adds, and ``defaults`` the values a repair of a
        stored request puts at a field (the encoder's form of the default).
        """

        hints = _class_type_hints(cls)
        members: dict[str, Shape] = {"kind": literal(kind)} if kind is not None else {}
        for item in fields(cls):
            members[_camel(item.name)] = (replace or {}).get(item.name) or (
                resource if item.name in files else of
            )(hints[item.name])
        wire[cls] = named(name or cls.__name__, obj({**members, **(extra or {})}, defaults=defaults))
        return wire[cls]

    # _encode_record: RecordInput; _SOURCE_FILE_RESOURCES names each source's resources.
    source = kinds("RecordSource", "RecordSource", {
        _camel(kind): {field: STRING for _resource_kind, field in resource_fields}
        for kind, resource_fields in _SOURCE_FILE_RESOURCES.items()
    })
    # The value domains the owners check beyond the kind: a field's bounds sit
    # beside the class it belongs to, each named with its owner's check.
    def bound(cls: type, name: str, **bounds: object) -> Shape:
        return bounded(of(_class_type_hints(cls)[name]), **bounds)

    # RecordDisplayOptions.__post_init__: a positive integer.
    data(RecordDisplayOptions, replace={"start_coordinate": bound(RecordDisplayOptions, "start_coordinate", minimum=1)})
    # _encode_record: every presentation value is the writer's default when unset.
    data(RecordPresentation, defaults={
        _camel(item.name): getattr(RecordPresentation(), item.name) for item in fields(RecordPresentation)
    })
    region_hints = _class_type_hints(RegionSpec)
    region = named("RecordRegion", obj({  # _encode_region
        "selector": nullable(selector),
        "start": of(region_hints["start"]),
        "end": of(region_hints["end"]),
        "reverseComplement": of(region_hints["reverse_complement"]),
    }))
    # A request has a record (validate_render_request); RecordInput.__post_init__
    # refuses a blank record key.
    records = bounded(array(data(RecordInput, "Record", replace={
        "source": source,
        "region": nullable(region),
        "record_key": bounded(STRING, pattern=r"\S"),  # written as given or as record-<n>
    }, defaults={
        "display": _encode_record_display(RecordDisplayOptions()),
        "presentation": _encode_record_presentation(RecordPresentation()),
        "cardinality": RecordInput.__dataclass_fields__["cardinality"].default.value,
    })), minItems=1)

    # _encode_output: output_directory is the caller's, overwrite is read and ignored.
    output_hints = _class_type_hints(RenderOutputRequest)
    output = named("RequestOutput", obj({
        "prefix": of(output_hints["output_prefix"]),
        "formats": bounded(array(STRING), minItems=1),  # RenderOutputRequest.__post_init__
        "overwrite": of(output_hints["overwrite"]),
        "interactiveMetadataPolicy": of(output_hints["interactive_metadata_policy"]),
    }))

    # _encode_layout: an empty object when the request has no layout.
    empty = named("EmptyLayout", obj({}))
    data(AlignmentAnchorIdentity)
    data(AlignmentRecordDecision)
    # SimilarityAlignmentPlan.__post_init__ refuses every other schema.
    data(SimilarityAlignmentPlan, replace={
        "schema": literal(SIMILARITY_ALIGNMENT_PLAN_SCHEMA),
        "records": bound(SimilarityAlignmentPlan, "records", minItems=1),
    })
    data(LinearRecordTranslation)
    layouts = {
        # CircularMultiRecordOptions.__post_init__: a radius ratio in (0, 1], gaps >= 0.
        "circular": union(empty, data(CircularMultiRecordOptions, "CircularLayout", replace={
            "multi_record_min_radius_ratio": bound(
                CircularMultiRecordOptions, "multi_record_min_radius_ratio", exclusiveMinimum=0, maximum=1
            ),
            "multi_record_column_gap_ratio": bound(CircularMultiRecordOptions, "multi_record_column_gap_ratio", minimum=0),
            "multi_record_row_gap_ratio": bound(CircularMultiRecordOptions, "multi_record_row_gap_ratio", minimum=0),
        }, defaults={
            _camel(item.name): getattr(CircularMultiRecordOptions(), item.name)
            for item in fields(CircularMultiRecordOptions) if item.name != "multi_record_positions"
        })),
        "linear": union(empty, data(LinearMultiRecordOptions, "LinearLayout", extra={
            "similarityAlignment": nullable(wire[SimilarityAlignmentPlan]),
        }, defaults={
            "recordGapPx": LinearMultiRecordOptions().record_gap_px,
            "recordTranslations": list(LinearMultiRecordOptions().record_translations),
            "similarityAlignment": None,
        })),
    }

    # _encode_diagram_options and _encode_option_value.
    data(ColorOptions, files=frozenset({"color_table_file", "default_colors_file"}))
    # annotations.models._color: a hex of 3 or 6 digits or a color name in any
    # case, with outer spaces; an unset fill is null or a blank text.
    color = named("AnnotationColor", union(
        bounded(STRING, pattern=HEX_COLOR_TEXT_PATTERN), bounded(STRING, pattern=COLOR_NAME_TEXT_PATTERN),
    ))
    hatch_default, style_default = HatchStyle(), RegionAnnotationStyle()
    data(HatchStyle, replace={"color": color}, defaults={"color": hatch_default.color})
    data(RegionAnnotationStyle, replace={  # RegionAnnotationStyle.__post_init__
        "fill_opacity": bound(RegionAnnotationStyle, "fill_opacity", minimum=0, maximum=1),
        "stroke": color, "label_color": color,
        "fill": union(color, bounded(STRING, pattern=BLANK_TEXT_PATTERN), NULL),
    }, defaults={"stroke": style_default.stroke, "labelColor": style_default.label_color, "fill": style_default.fill})
    data(FeatureSelector)
    # A Web writer wrote a fraction at lane, start and end (OV-630, BA-OD74);
    # the decoder reads its integer part.
    # CoordinateSpan.__post_init__: 1-based.
    data(CoordinateSpan, kind="coordinateSpan", replace={
        "start": bounded(NUMBER, minimum=1), "end": bounded(NUMBER, minimum=1),
    })
    # FeatureSpan.__post_init__: at least one selector.
    data(FeatureSpan, kind="featureSpan", replace={"selectors": bound(FeatureSpan, "selectors", minItems=1)})
    data(FeatureIdentitySpan, kind="featureIdentity")
    data(RegionAnnotation, replace={"lane": nullable(NUMBER)})
    data(AnnotationSet)
    data(AnnotationOptions, files=frozenset({"table_file"}))
    scalar = data(ScalarSpec)
    data(CircularTrackSlot, kind="circularTrackSlot")
    # A Linear slot's height and spacing are pixels (gbdraw.tracks.linear._validate_px_scalar).
    pixels = nullable(named("LinearScalarSpec", obj(
        {**types_[scalar["name"]]["fields"], "unit": literal("px")}
    )))
    data(LinearTrackSlot, kind="linearTrackSlot", replace={"height": pixels, "spacing": pixels})
    data(FeatureOverride)  # FeatureOverride.to_mapping
    target_hints = _class_type_hints(FeaturePlacementTarget)
    placement_hints = _class_type_hints(FeaturePlacementOverride)
    modes: dict[Literal["circular", "linear"], tuple[type, type, type]] = {
        "circular": (CircularDiagramOptions, CircularRequestTrackOptions, CircularOutputOptions),
        "linear": (LinearDiagramOptions, LinearRequestTrackOptions, LinearOutputOptions),
    }
    option_files = (
        _TABLE_FIELDS | _FILE_FIELDS | _TABLE_SEQUENCE_FIELDS | _FILE_SEQUENCE_FIELDS
        | _OPTIONAL_FILE_SEQUENCE_FIELDS
    )
    # The decoder also reads the depth inputs of older writers (the writer
    # stores every depth input as depthTracks).
    unwritten = (
        _COMPARISON_FIELDS | _PLACEMENT_INPUT_FIELDS | _RING_SEARCH_FIELDS | _DEPTH_COMPATIBILITY_FIELDS
        | _TABLE_PATH_FIELDS | _PLANNER_ONLY_FIELDS
    )

    def option(name: str, hint: object, mode: Literal["circular", "linear"]) -> Shape:
        """The hint's shape, without null where the option's decoder refuses null."""

        shape = (resource if name in option_files else of)(hint)
        try:
            _decode_option_value(
                name, None, mode=mode, schema=CANONICAL_REQUEST_SCHEMA, resource_paths={}
            )
        except CanonicalRequestDecodingError:
            return non_null(shape)
        return shape

    diagram_options: dict[str, Shape] = {}
    for mode, (options_class, track_class, output_class) in modes.items():
        title = mode.capitalize()
        # A depth track has a height only in a Linear diagram
        # (CircularDiagramOptions.__post_init__).
        # DepthTrackInput.__post_init__: at least one source, and positive numbers.
        depth_hints = _class_type_hints(DepthTrackInput)
        data(DepthTrackInput, f"{title}DepthTrack", files=frozenset({"source"}), replace={
            "source": bounded(resource(depth_hints["source"]), minItems=1),
            **{name: bounded(of(depth_hints[name]), exclusiveMinimum=0) for name in DEPTH_TRACK_POSITIVE_FIELDS},
            **({"height": NULL} if mode == "circular" else {}),
        })
        # FeaturePlacementOverride.from_mapping: a lane has a level and a side
        # of its mode (FeaturePlacementTarget.validate_mode), main neither.
        wire[FeaturePlacementOverride] = named(f"{title}FeaturePlacement", obj({
            "recordKey": of(placement_hints["record_key"]),
            "biologicalFeatureId": of(placement_hints["biological_feature_id"]),
            "placement": union(obj({"kind": literal("main")}), obj({
                "kind": literal("lane"),
                "side": literal(*_LANE_SIDES[mode]),
                # FeaturePlacementTarget: level 1 only.
                "level": bounded(non_null(of(target_hints["level"])), minimum=1, maximum=1),
            })),
        }))
        # _encode_tracks writes the track fields of both modes, the other mode's as null.
        # _validate_track_configuration: a Linear slot list holds a slot.
        track_hints = _class_type_hints(track_class)
        track_bounds: dict[str, dict[str, object]] = {"linear_track_slots": {"minItems": 1}}
        wire[track_class] = named(f"{title}TrackOptions", obj({
            _camel(item.name): (
                bounded(of(track_hints[item.name]), **track_bounds[item.name])
                if item.name in track_bounds else of(track_hints[item.name])
            ) if owner is track_class else NULL
            for _options, owner, _output in modes.values()
            for item in fields(owner)
        }))
        data(output_class)
        hints = _class_type_hints(options_class)
        members: dict[str, Shape] = {
            "featurePlacements": array(wire[FeaturePlacementOverride]),
            "featureOverrides": array(wire[FeatureOverride]),
        }
        required = sorted(members)
        # The diagram options' own checks: depth tracks hold a track, and a
        # dinucleotide is two bases (validate_dinucleotide).
        option_bounds: dict[str, dict[str, object]] = {
            "depth_tracks": {"minItems": 1}, "dinucleotide": {"pattern": DINUCLEOTIDE_PATTERN},
        }
        for item in fields(options_class):
            if item.name not in unwritten:
                members[_option_wire_key(item.name)] = bounded(
                    option(item.name, hints[item.name], mode), **option_bounds.get(item.name, {})
                )
        # A shared option of the other mode is accepted at its default only.
        for name, default in _SHARED_OPTION_WRONG_MODE_DEFAULTS[mode].items():
            members.setdefault(_option_wire_key(name), literal(default))
        # The decoder leaves an absent option at the dataclass default; the
        # default is written as the encoder writes it, and None is "unset".
        defaults = _DEFAULT_CIRCULAR_OPTIONS if mode == "circular" else _DEFAULT_LINEAR_OPTIONS
        diagram_options[mode] = named(f"{title}DiagramOptions", obj(
            members, optional=sorted(set(members) - set(required)), defaults={
                # _encode_diagram_options writes no row as the empty list.
                **{name: [] for name in required},
                **{
                    _option_wire_key(item.name): _encode_option_value(
                        item.name, getattr(defaults, item.name), resources=SessionResourceTable().request()
                    )
                    for item in fields(options_class)
                    if item.name not in unwritten and getattr(defaults, item.name) is not None
                },
            },
        ))

    # _encode_comparisons; _PROTEIN_SETTINGS_WIRE names each setting's typed owner.
    wire[LosslessCollinearityParameters] = named("CollinearityParameters", obj({
        "kind": literal("lossless"),
        "parameters": data(LosslessCollinearityParameters, replace={
            # Released schema 9 writers stored "either" (the Web default),
            # which every reader reads as "strand" (OV-2301); writers write
            # only the typed names.
            "merge_orientation": literal(
                *COLLINEARITY_MERGE_ORIENTATIONS,
                *COLLINEARITY_MERGE_ORIENTATION_SPELLINGS["aliases"],
                spellings=COLLINEARITY_MERGE_ORIENTATION_SPELLINGS,
            ),
        }),
    }))
    setting_hints = {
        "options": _class_type_hints(LinearDiagramOptions),
        "search": _class_type_hints(LosatSearchOptions),
        "runtime": _class_type_hints(LosatRuntimeOptions),
    }
    # LinearDiagramOptions.__post_init__ normalizes these settings, so the
    # encoder writes one of the names their normalizers check. The decoder
    # reads other spellings of a name as well (OV-651).
    # The description names the other spellings they read, which a reader of a
    # stored request reads as the name.
    setting_names: dict[str, Shape] = {
        "collinearity_unit_mode": literal(*_COLLINEARITY_UNIT_MODES, spellings=COLLINEARITY_UNIT_MODE_SPELLINGS),
        "collinearity_anchor_mode": literal(
            *COLLINEARITY_ANCHOR_MODES, spellings=COLLINEARITY_ANCHOR_MODE_SPELLINGS
        ),
        "collinearity_search_scope": literal(
            *COLLINEARITY_SEARCH_SCOPES, spellings=COLLINEARITY_SEARCH_SCOPE_SPELLINGS
        ),
        "collinearity_color_mode": literal(*COLLINEARITY_COLOR_MODES, spellings=COLLINEARITY_COLOR_MODE_SPELLINGS),
        "orthogroup_membership_mode": literal(
            *ORTHOGROUP_MEMBERSHIP_MODES, spellings=ORTHOGROUP_MEMBERSHIP_MODE_SPELLINGS
        ),
    }
    pair = {"queryRecordIndex": INTEGER, "subjectRecordIndex": INTEGER}
    comparison = kinds("Comparison", "Comparison", {
        "nucleotideBlast": {"resourceId": STRING, **pair},
        "precomputedProteinComparison": {
            "resourceId": STRING, "encoding": literal("canonicalTsv"), **pair,
        },
        "orthogroupResult": {"resourceId": STRING, "encoding": literal("canonicalJson")},
        "collinearityResult": {
            "resourceId": STRING,
            "encoding": literal("canonicalJson"),
            "valueKind": literal("result", "blocks"),
        },
        "generatedProteinComparison": {
            "mode": literal(*LOSATP_MODE_WIRE.values()),
            "pairs": array(obj(pair)),
            "settings": named("ProteinComparisonSettings", obj(
                {
                    key: setting_names[name] if name in setting_names else of(setting_hints[owner][name])
                    for key, owner, name in _PROTEIN_SETTINGS_WIRE
                },
                optional=sorted(_OPTIONAL_PROTEIN_SETTINGS),
                # A stored retired anchor mode reads as its replacement.
                defaults={
                    **{key: _absent_protein_setting(key) for key in _OPTIONAL_PROTEIN_SETTINGS},
                    "collinearityAnchorMode": RETIRED_COLLINEARITY_ANCHOR_MODES["one_to_one"],
                },
            )),
        },
    })

    # _encode_canonical_request and decode_canonical_request.
    # A CircularDiagramRequest resolves its grouping; a CircularBatchRequest is "batch".
    grouping = union(
        non_null(of(_class_type_hints(CircularDiagramRequest)["grouping"])), literal("batch")
    )
    requests = []
    for mode in modes:
        top = {
            "schema": literal(CANONICAL_REQUEST_SCHEMA),
            "mode": literal(mode),
            "grouping": grouping if mode == "circular" else literal("single"),
            "records": records,
            "diagramOptions": diagram_options[mode],
            "layout": layouts[mode],
            # A Circular request has no comparison item.
            "comparisons": array(comparison) if mode == "linear" else tuple_of(),
            # A Circular batch has one output per record.
            "output": bounded(union(output, array(output)), minItems=1) if mode == "circular" else output,
        }
        if set(top) != _TOP_LEVEL_FIELDS_V5:
            raise RuntimeError("The request description and _TOP_LEVEL_FIELDS_V5 differ.")
        # What a repair of a stored request puts at a part it cannot read:
        # what the encoder writes for the default options and output.
        options_class = modes[mode][0]
        requests.append(named(f"{mode.capitalize()}RenderRequest", obj(top, defaults={
            "diagramOptions": _encode_diagram_options(
                options_class(), record_count=1, resources=SessionResourceTable().request()
            ),
            "output": _encode_output(RenderOutputRequest()),
        })))
    named("RenderRequest", union(*requests))
    # A required field that holds null takes null where its stored value
    # cannot be read: the encoder writes null for an unset value (OV-777).
    # Null is no empty form where another value decides whether it may be
    # null: a decision's anchor and the anchor's identity
    # (AlignmentRecordDecision.__post_init__; the plan finds its reference by
    # the identity), a slot list an axis index names
    # (_validate_track_configuration), and the edits of a feature override
    # (FeatureOverride: at least one edit).
    null_decided_elsewhere = {
        "CanonicalAlignmentRecordDecision": {"anchor"},
        "CanonicalAlignmentAnchorIdentity": {"sourceFeatureIndex", "stableFeatureSvgId"},
        "CanonicalCircularTrackOptions": {"circularTrackSlots"},
        "CanonicalLinearTrackOptions": {"linearTrackSlots"},
        "CanonicalFeatureOverride": {"featureVisibility", "labelText", "labelVisibility"},
    }
    for name, shape in types_.items():
        if shape["kind"] == "object":
            nullable_fields = {
                key: None
                for key, field in shape["fields"].items()
                if key not in shape["optional"] and key not in null_decided_elsewhere.get(name, set())
                and field != UNKNOWN and not shape_violations(None, field, types_, key)
            }
            shape_defaults = {**nullable_fields, **shape.get("defaults", {})}
            if shape_defaults:
                shape["defaults"] = dict(sorted(shape_defaults.items()))
    # Only a default that fits its field is described.
    for shape in types_.values():
        if shape["kind"] == "object" and "defaults" in shape:
            shape["defaults"] = {
                key: value
                for key, value in shape["defaults"].items()
                if not shape_violations(value, shape["fields"][key], types_, key)
            }
            if not shape["defaults"]:
                del shape["defaults"]
    return {"root": "CanonicalRenderRequest", "types": types_}


__all__ = [
    "CANONICAL_REQUEST_SCHEMA",
    "SUPPORTED_CANONICAL_REQUEST_SCHEMAS",
    "UNKNOWN_FIELD_POLICY",
    "CanonicalRequestCodecError",
    "CanonicalRequestDecodingError",
    "CanonicalRequestEncodingError",
    "CanonicalRequestResource",
    "EncodedCanonicalRequest",
    "decode_canonical_typed_resource",
    "decode_canonical_request",
    "describe_canonical_request",
    "describe_canonical_typed_resources",
    "describe_legacy_request_tables",
    "encode_canonical_table",
    "encode_canonical_typed_resource",
    "encode_canonical_request",
    "migrate_flat_config_overrides",
    "migrate_legacy_full_config",
    "with_automatic_scale_interval",
    "with_stored_font_sizes_as_auto",
]
