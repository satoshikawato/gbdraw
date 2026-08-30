"""Small, typed, beginner-facing Python API for gbdraw."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from functools import wraps
from inspect import Parameter, signature
from os import PathLike
from pathlib import Path
from typing import Literal, Mapping, Sequence, TypeAlias

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from Bio import SeqIO  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]

from gbdraw.analysis.collinearity import (
    CollinearityAnchorMode,
    CollinearityBlock,
    CollinearityColorMode,
    CollinearityResult,
    CollinearitySearchScope,
    LosslessCollinearityParameters,
)
from gbdraw.analysis.collinearity_units import CollinearityUnitMode
from gbdraw.analysis.protein_colinearity import (
    OrthogroupMembershipMode,
    OrthogroupResult,
)
from gbdraw.annotations import AnnotationOptions
from gbdraw.api.io import load_gbks as _load_gbks, load_gff_fasta as _load_gff_fasta
from gbdraw.api.options import (
    CircularDiagramOptions as _CircularDiagramOptions,
    CircularMultiRecordOptions as _CircularLayout,
    CircularOutputOptions as _CircularOutputOptions,
    CircularRequestTrackOptions as _CircularRequestTrackOptions,
    ColorOptions as _ColorOptions,
    LinearDiagramOptions as _LinearDiagramOptions,
    LinearMultiRecordOptions as _LinearLayout,
    LinearOutputOptions as _LinearOutputOptions,
    LinearRequestTrackOptions as _LinearRequestTrackOptions,
    _validate_center_reserved_radius,
    _validate_track_configuration,
)
from gbdraw.render.output_paths import preflight_output_paths
from gbdraw.api.request_render import (
    PreparedDiagramRequest as _PreparedDiagramRequest,
    build_prepared_interactive_context as _build_prepared_interactive_context,
    build_request_diagram as _build_request_diagram,
)
from gbdraw.api.requests import (
    CircularDiagramRequest as _CircularDiagramRequest,
    InMemoryRecordSource as _InMemoryRecordSource,
    LinearDiagramRequest as _LinearDiagramRequest,
    RecordInput as _RecordInput,
)
from gbdraw.api.render import render_to_bytes
from gbdraw.exceptions import ExportError, ValidationError
from gbdraw.linear_comparison import LinearComparison
from gbdraw.config.models import GbdrawConfig
from gbdraw.mode_profiles import (
    CIRCULAR_MODE_PROFILE,
    ComparisonThresholds,
    DEFAULT_FEATURE_TYPES,
    LINEAR_MODE_PROFILE,
)
from gbdraw.render.formats import INTERACTIVE_SVG_FORMAT, normalize_format_token
from gbdraw.render.interactive_context import require_interactive_svg_metadata
from gbdraw.render.interactive_svg import InteractiveSvgContext
from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot


TableSource: TypeAlias = DataFrame | str | PathLike[str]
RecordCollection: TypeAlias = SeqRecord | Sequence[SeqRecord]


@dataclass(frozen=True)
class FeatureOptions:
    """Feature selection, colors, visibility, and shape overrides."""

    types: Sequence[str] | None = DEFAULT_FEATURE_TYPES
    color_table: TableSource | None = None
    default_colors: TableSource | None = None
    palette: str = "default"
    visibility: TableSource | None = None
    shapes: Mapping[str, str] | None = None


@dataclass(frozen=True)
class LabelOptions:
    """Optional label filtering and replacement tables."""

    whitelist: TableSource | None = None
    qualifier_priority: TableSource | None = None
    overrides: TableSource | None = None


@dataclass(frozen=True)
class TitleOptions:
    """Diagram title text and placement."""

    text: str | None = None
    position: Literal["none", "center", "top", "bottom"] | None = None
    font_size: float | None = None


@dataclass(frozen=True)
class Thresholds:
    """Optional comparison filters resolved by the owning diagram mode."""

    evalue: float | None = None
    bitscore: float | None = None
    identity: float | None = None
    alignment_length: int | None = None

    def __post_init__(self) -> None:
        validated = ComparisonThresholds(
            evalue=0.0 if self.evalue is None else self.evalue,
            bitscore=0.0 if self.bitscore is None else self.bitscore,
            identity=0.0 if self.identity is None else self.identity,
            alignment_length=(
                0 if self.alignment_length is None else self.alignment_length
            ),
        )
        for field_name in (
            "evalue",
            "bitscore",
            "identity",
            "alignment_length",
        ):
            if getattr(self, field_name) is not None:
                object.__setattr__(
                    self,
                    field_name,
                    getattr(validated, field_name),
                )


def _resolve_thresholds(
    thresholds: Thresholds,
    defaults: ComparisonThresholds,
) -> Thresholds:
    return Thresholds(
        evalue=defaults.evalue if thresholds.evalue is None else thresholds.evalue,
        bitscore=(
            defaults.bitscore
            if thresholds.bitscore is None
            else thresholds.bitscore
        ),
        identity=(
            defaults.identity
            if thresholds.identity is None
            else thresholds.identity
        ),
        alignment_length=(
            defaults.alignment_length
            if thresholds.alignment_length is None
            else thresholds.alignment_length
        ),
    )


@dataclass(frozen=True)
class DepthTrackOptions:
    """One depth track, with one source per displayed record when needed."""

    source: TableSource | Sequence[TableSource | None]
    label: str | None = None
    color: str | None = None
    height: float | str | None = None
    large_tick_interval: float | str | None = None
    small_tick_interval: float | str | None = None
    tick_font_size: float | str | None = None

    def __post_init__(self) -> None:
        source = self.source
        if isinstance(source, (DataFrame, str, PathLike)):
            _source(source, name="depth track source")
            return
        if isinstance(source, (str, bytes)) or not isinstance(source, Sequence):
            raise ValidationError(
                "DepthTrackOptions.source must be a table, path, or sequence "
                "of tables, paths, and None values."
            )
        for index, item in enumerate(source):
            _source(item, name=f"DepthTrackOptions.source[{index}]")


@dataclass(frozen=True)
class CircularTrackOptions:
    """Circular track order and axis placement."""

    slots: Sequence[str | CircularTrackSlot] | None = None
    axis_index: int | None = None
    center_reserved_radius: float | None = None

    def __post_init__(self) -> None:
        axis_index = _validate_track_configuration(
            self.slots,
            axis_index=self.axis_index,
            mode="circular",
            expected_type=CircularTrackSlot,
            slots_field_name="CircularTrackOptions.slots",
            axis_field_name="CircularTrackOptions.axis_index",
        )
        object.__setattr__(self, "axis_index", axis_index)
        object.__setattr__(
            self,
            "center_reserved_radius",
            _validate_center_reserved_radius(
                self.center_reserved_radius,
                field_name="CircularTrackOptions.center_reserved_radius",
            ),
        )


@dataclass(frozen=True)
class LinearTrackOptions:
    """Linear track order and axis placement."""

    slots: Sequence[str | LinearTrackSlot] | None = None
    axis_index: int | None = None

    def __post_init__(self) -> None:
        axis_index = _validate_track_configuration(
            self.slots,
            axis_index=self.axis_index,
            mode="linear",
            expected_type=LinearTrackSlot,
            slots_field_name="LinearTrackOptions.slots",
            axis_field_name="LinearTrackOptions.axis_index",
        )
        object.__setattr__(self, "axis_index", axis_index)


@dataclass(frozen=True)
class CircularLayout:
    """Grid layout used automatically when a circular diagram has multiple records."""

    size: Literal["linear", "auto", "equal"] = "auto"
    min_radius_ratio: float = 0.55
    column_gap_ratio: float = 0.10
    row_gap_ratio: float = 0.05
    positions: Sequence[str] | None = None

    def __post_init__(self) -> None:
        normalized = self._legacy()
        object.__setattr__(self, "size", normalized.multi_record_size_mode)
        object.__setattr__(
            self,
            "min_radius_ratio",
            normalized.multi_record_min_radius_ratio,
        )
        object.__setattr__(
            self,
            "column_gap_ratio",
            normalized.multi_record_column_gap_ratio,
        )
        object.__setattr__(
            self,
            "row_gap_ratio",
            normalized.multi_record_row_gap_ratio,
        )
        object.__setattr__(self, "positions", normalized.multi_record_positions)

    def _legacy(self) -> _CircularLayout:
        return _CircularLayout(
            multi_record_size_mode=self.size,
            multi_record_min_radius_ratio=self.min_radius_ratio,
            multi_record_column_gap_ratio=self.column_gap_ratio,
            multi_record_row_gap_ratio=self.row_gap_ratio,
            multi_record_positions=self.positions,
        )


@dataclass(frozen=True)
class LinearLayout:
    """Optional layout for linear records arranged in multiple rows."""

    record_gap: float = 24.0
    positions: Sequence[str] | None = None

    def __post_init__(self) -> None:
        normalized = self._legacy()
        object.__setattr__(self, "record_gap", normalized.record_gap_px)
        object.__setattr__(self, "positions", normalized.multi_record_positions)

    def _legacy(self) -> _LinearLayout:
        return _LinearLayout(
            record_gap_px=self.record_gap,
            multi_record_positions=self.positions,
        )


@dataclass(frozen=True)
class ComparisonRingTrackOptions:
    """One circular comparison ring drawn from sequence-similarity hits."""

    source: TableSource
    label: str | None = None
    color: str | None = None
    comparison_sequence_source: RecordCollection | str | PathLike[str] | None = None


@dataclass(frozen=True)
class ComparisonRingOptions:
    """Circular sequence-similarity comparison rings and their shared geometry."""

    tracks: Sequence[ComparisonRingTrackOptions] = ()
    reference: Literal["query", "subject", "auto"] = "auto"
    ring_width: float | None = None
    ring_gap: float | None = None


# Compatibility aliases for the original package-root names.
ConservationTrackOptions = ComparisonRingTrackOptions
ConservationOptions = ComparisonRingOptions


_LINEAR_DIAGRAM_DEFAULTS = _LinearDiagramOptions()


@dataclass(frozen=True)
class LinearComparisonOptions:
    """Precomputed or in-process comparison inputs for a linear diagram."""

    blast_files: Sequence[str] | None = None
    comparisons: Sequence[LinearComparison] | None = None
    protein_comparisons: Sequence[DataFrame] | None = None
    orthogroups: OrthogroupResult | None = None
    protein_mode: Literal["none", "pairwise", "orthogroup", "collinear"] = (
        _LINEAR_DIAGRAM_DEFAULTS.protein_blastp_mode
    )
    pairs: Sequence[tuple[int, int]] | None = None
    match_style: Literal["ribbon", "curve"] = "ribbon"
    collinearity_blocks: CollinearityResult | Sequence[CollinearityBlock] | None = None
    collinearity_params: LosslessCollinearityParameters | None = None
    collinearity_unit: CollinearityUnitMode | str = (
        _LINEAR_DIAGRAM_DEFAULTS.collinearity_unit_mode
    )
    collinearity_anchor: CollinearityAnchorMode | str = (
        _LINEAR_DIAGRAM_DEFAULTS.collinearity_anchor_mode
    )
    collinearity_scope: CollinearitySearchScope | str = (
        _LINEAR_DIAGRAM_DEFAULTS.collinearity_search_scope
    )
    collinearity_color: CollinearityColorMode | str = (
        _LINEAR_DIAGRAM_DEFAULTS.collinearity_color_mode
    )
    losat_executable: str = _LINEAR_DIAGRAM_DEFAULTS.losatp_bin
    blastp_executable: str | None = _LINEAR_DIAGRAM_DEFAULTS.ncbi_blastp_bin
    threads: int | None = _LINEAR_DIAGRAM_DEFAULTS.losatp_threads
    max_hits: int = _LINEAR_DIAGRAM_DEFAULTS.protein_blastp_max_hits
    candidate_limit: int | None = (
        _LINEAR_DIAGRAM_DEFAULTS.protein_blastp_candidate_limit
    )
    orthogroup_membership: OrthogroupMembershipMode | str = (
        _LINEAR_DIAGRAM_DEFAULTS.orthogroup_membership_mode
    )
    orthogroup_member_max_hits: int = (
        _LINEAR_DIAGRAM_DEFAULTS.orthogroup_member_max_hits
    )
    max_paralog_links: int = (
        _LINEAR_DIAGRAM_DEFAULTS.collinear_max_paralog_links_per_orthogroup
    )
    align_feature: str | None = None


@dataclass(frozen=True)
class _CommonOptions:
    features: FeatureOptions = field(default_factory=FeatureOptions)
    labels: LabelOptions = field(default_factory=LabelOptions)
    title: TitleOptions = field(default_factory=TitleOptions)
    legend: str = "right"
    annotations: AnnotationOptions | None = None
    config: GbdrawConfig | dict[str, object] | None = None
    config_overrides: Mapping[str, object] | None = None
    dinucleotide: str = "GC"
    window: int | None = None
    step: int | None = None
    depth_window: int | None = None
    depth_step: int | None = None
    depth_tracks: Sequence[DepthTrackOptions] = ()
    thresholds: Thresholds = field(default_factory=Thresholds)


def _validate_common_options(options: _CommonOptions) -> None:
    expected_types = (
        ("features", options.features, FeatureOptions),
        ("labels", options.labels, LabelOptions),
        ("title", options.title, TitleOptions),
        ("thresholds", options.thresholds, Thresholds),
    )
    for field_name, value, expected_type in expected_types:
        if not isinstance(value, expected_type):
            raise ValidationError(
                f"{field_name} must be {expected_type.__name__}."
            )
    if (
        isinstance(options.depth_tracks, (str, bytes))
        or not isinstance(options.depth_tracks, Sequence)
        or not all(
            isinstance(track, DepthTrackOptions) for track in options.depth_tracks
        )
    ):
        raise ValidationError(
            "depth_tracks must contain DepthTrackOptions values."
        )


@dataclass(frozen=True)
class CircularOptions(_CommonOptions):
    """Options for :func:`draw_circular`.

    ``comparison_rings`` is the canonical sequence-similarity ring field.
    ``conservation`` remains a runtime compatibility constructor and attribute alias.
    """

    tracks: CircularTrackOptions = field(default_factory=CircularTrackOptions)
    comparison_rings: ComparisonRingOptions | None = None
    species: str | None = None
    strain: str | None = None
    keep_full_definition_with_title: bool = False
    thresholds: Thresholds = field(default_factory=Thresholds)

    def __post_init__(self) -> None:
        _validate_common_options(self)
        object.__setattr__(
            self,
            "thresholds",
            _resolve_thresholds(
                self.thresholds,
                CIRCULAR_MODE_PROFILE.comparison,
            ),
        )
        if not isinstance(self.tracks, CircularTrackOptions):
            raise ValidationError("tracks must be CircularTrackOptions.")
        comparison_rings = self.comparison_rings
        if comparison_rings is None:
            comparison_rings = ComparisonRingOptions()
        if not isinstance(comparison_rings, ComparisonRingOptions):
            raise ValidationError(
                "comparison_rings must be ComparisonRingOptions "
                "(ConservationOptions is a compatibility alias)."
            )
        object.__setattr__(self, "comparison_rings", comparison_rings)
        if self.title.position not in {None, "none", "top", "bottom"}:
            raise ValidationError(
                "Circular title position must be one of: none, top, bottom."
            )
        for track in self.depth_tracks:
            if track.height is not None:
                raise ValidationError(
                    "Depth track height is supported only by linear diagrams."
                )

    @property
    def conservation(self) -> ComparisonRingOptions:
        """Compatibility attribute alias for :attr:`comparison_rings`."""

        comparison_rings = self.comparison_rings
        assert comparison_rings is not None
        return comparison_rings


_CIRCULAR_OPTIONS_GENERATED_INIT = CircularOptions.__init__
_CIRCULAR_OPTIONS_INIT_SIGNATURE = signature(_CIRCULAR_OPTIONS_GENERATED_INIT)


@wraps(_CIRCULAR_OPTIONS_GENERATED_INIT)
def _circular_options_init(
    self: CircularOptions,
    *args: object,
    conservation: ComparisonRingOptions | None = None,
    **kwargs: object,
) -> None:
    bound = _CIRCULAR_OPTIONS_INIT_SIGNATURE.bind_partial(self, *args, **kwargs)
    if conservation is not None:
        if bound.arguments.get("comparison_rings") is not None:
            raise ValidationError(
                "Pass comparison_rings or the conservation compatibility alias, "
                "not both."
            )
        bound.arguments["comparison_rings"] = conservation
    _CIRCULAR_OPTIONS_GENERATED_INIT(*bound.args, **bound.kwargs)


_circular_options_init.__signature__ = _CIRCULAR_OPTIONS_INIT_SIGNATURE.replace(  # type: ignore[attr-defined]
    parameters=(
        *_CIRCULAR_OPTIONS_INIT_SIGNATURE.parameters.values(),
        Parameter(
            "conservation",
            kind=Parameter.KEYWORD_ONLY,
            default=None,
            annotation=ComparisonRingOptions | None,
        ),
    ),
)
CircularOptions.__init__ = _circular_options_init  # type: ignore[method-assign]


@dataclass(frozen=True)
class LinearOptions(_CommonOptions):
    """Options accepted only by :func:`draw_linear`."""

    tracks: LinearTrackOptions = field(default_factory=LinearTrackOptions)
    comparisons: LinearComparisonOptions = field(default_factory=LinearComparisonOptions)
    thresholds: Thresholds = field(default_factory=Thresholds)

    def __post_init__(self) -> None:
        _validate_common_options(self)
        object.__setattr__(
            self,
            "thresholds",
            _resolve_thresholds(
                self.thresholds,
                LINEAR_MODE_PROFILE.comparison,
            ),
        )
        if not isinstance(self.tracks, LinearTrackOptions):
            raise ValidationError("tracks must be LinearTrackOptions.")
        if not isinstance(self.comparisons, LinearComparisonOptions):
            raise ValidationError("comparisons must be LinearComparisonOptions.")
        if self.title.position not in {None, "center", "top", "bottom"}:
            raise ValidationError(
                "Linear title position must be one of: center, top, bottom."
            )


class Diagram:
    """A rendered diagram that can be serialized or saved without SVG internals."""

    def __init__(
        self,
        drawing: Drawing,
        *,
        mode: Literal["circular", "linear"],
        records: Sequence[SeqRecord],
        interactive_context: (
            InteractiveSvgContext
            | Callable[[], InteractiveSvgContext]
            | None
        ),
    ) -> None:
        self._drawing = drawing
        self.mode = mode
        self.records = tuple(records)
        self._interactive_context = interactive_context

    def _resolve_interactive_context(self) -> InteractiveSvgContext | None:
        if callable(self._interactive_context):
            self._interactive_context = require_interactive_svg_metadata(
                self._interactive_context
            )
        return self._interactive_context

    def to_bytes(self, format: str = "svg") -> bytes:
        """Return the diagram in SVG or a supported binary format."""

        return render_to_bytes(
            self._drawing,
            format,
            interactive_context=(
                self._resolve_interactive_context()
                if normalize_format_token(format) == INTERACTIVE_SVG_FORMAT
                else None
            ),
        )

    def to_svg(self, *, interactive: bool = False) -> str:
        """Return SVG text, optionally with interactive feature metadata."""

        format = "interactive_svg" if interactive else "svg"
        return self.to_bytes(format).decode("utf-8")

    def save(
        self,
        path: str | PathLike[str],
        *,
        format: str | None = None,
        overwrite: bool = False,
    ) -> Path:
        """Write exactly one output file and return its path."""

        output_path = Path(path)
        resolved_format = format or _format_from_path(output_path)
        preflight_output_paths((output_path,), overwrite=True)
        if output_path.exists() and not overwrite:
            raise ValidationError(
                f"Output file already exists: {output_path}. "
                "Use overwrite=True to replace it."
            )
        try:
            output_path.parent.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise ExportError(
                f"Could not prepare output directory: {output_path.parent}"
            ) from exc
        payload = self.to_bytes(resolved_format)
        if overwrite and (output_path.exists() or output_path.is_symlink()):
            try:
                output_path.unlink()
            except OSError as exc:
                raise ExportError(
                    f"Could not replace existing output file: {output_path}"
                ) from exc
        try:
            with output_path.open("xb") as output_file:
                output_file.write(payload)
        except FileExistsError as exc:
            message = (
                f"Output target appeared during export: {output_path}. "
                "Refusing to replace it."
                if overwrite
                else (
                    f"Output file already exists: {output_path}. "
                    "Use overwrite=True to replace it."
                )
            )
            raise ValidationError(message) from exc
        except OSError as exc:
            raise ExportError(f"Could not write output file: {output_path}") from exc
        return output_path


def _format_from_path(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".interactive.svg"):
        return "interactive_svg"
    suffix = path.suffix.lower().lstrip(".")
    if suffix in {"svg", "png", "pdf", "eps", "ps"}:
        return suffix
    raise ValidationError(
        "Could not infer the output format from the path; pass format explicitly."
    )


def _records(value: RecordCollection) -> tuple[SeqRecord, ...]:
    if isinstance(value, SeqRecord):
        return (value,)
    try:
        records = tuple(value)
    except TypeError as exc:
        raise ValidationError("records must be a SeqRecord or a sequence of SeqRecord values.") from exc
    if not records:
        raise ValidationError("At least one record is required.")
    if not all(isinstance(record, SeqRecord) for record in records):
        raise ValidationError("Every record must be a Bio.SeqRecord.SeqRecord.")
    return records


def _paths(value: str | PathLike[str] | Sequence[str | PathLike[str]]) -> list[str]:
    if isinstance(value, (str, PathLike)):
        normalized = [str(value)]
    else:
        normalized = [str(path) for path in value]
    if not normalized or any(not path.strip() for path in normalized):
        raise ValidationError("At least one non-empty input path is required.")
    return normalized


def read_genbank(
    paths: str | PathLike[str] | Sequence[str | PathLike[str]],
) -> list[SeqRecord]:
    """Read all records from one or more GenBank files."""

    return _load_gbks(_paths(paths))


def read_gff(
    gff_paths: str | PathLike[str] | Sequence[str | PathLike[str]],
    fasta_paths: str | PathLike[str] | Sequence[str | PathLike[str]],
    *,
    features: Sequence[str] | None = None,
) -> list[SeqRecord]:
    """Read records from paired GFF3 and FASTA files."""

    normalized_gff = _paths(gff_paths)
    normalized_fasta = _paths(fasta_paths)
    if len(normalized_gff) != len(normalized_fasta):
        raise ValidationError("GFF3 and FASTA inputs must contain the same number of paths.")
    return _load_gff_fasta(
        normalized_gff,
        normalized_fasta,
        selected_features_set=features,
    )


def _source(value: TableSource | None, *, name: str) -> tuple[DataFrame | None, str | None]:
    if value is None:
        return None, None
    if isinstance(value, DataFrame):
        return value, None
    if isinstance(value, (str, PathLike)):
        path = str(value).strip()
        if not path:
            raise ValidationError(f"{name} path must not be empty.")
        return None, path
    raise ValidationError(f"{name} must be a pandas DataFrame or a path.")


def _depth_kwargs(
    tracks: Sequence[DepthTrackOptions],
    *,
    record_count: int,
    mode: Literal["circular", "linear"],
) -> dict[str, object]:
    if not tracks:
        return {}
    table_rows: list[list[DataFrame | None]] = [list() for _ in range(record_count)]
    file_rows: list[list[str | None]] = [list() for _ in range(record_count)]
    for track_index, track in enumerate(tracks, start=1):
        source = track.source
        if isinstance(source, (DataFrame, str, PathLike)):
            if record_count != 1:
                raise ValidationError(
                    f"Depth track {track_index} needs one source per record."
                )
            sources: Sequence[TableSource | None] = (source,)
        else:
            sources = tuple(source)
            if len(sources) != record_count:
                raise ValidationError(
                    f"Depth track {track_index} has {len(sources)} source(s); "
                    f"expected {record_count}."
                )
        for record_index, item in enumerate(sources):
            table, file_path = _source(item, name=f"depth track {track_index}")
            table_rows[record_index].append(table)
            file_rows[record_index].append(file_path)
        if mode == "circular" and track.height is not None:
            raise ValidationError("Depth track height is supported only by linear diagrams.")
    result: dict[str, object] = {
        "depth_track_labels": [track.label for track in tracks],
        "depth_track_colors": [track.color for track in tracks],
        "depth_track_large_tick_intervals": [track.large_tick_interval for track in tracks],
        "depth_track_small_tick_intervals": [track.small_tick_interval for track in tracks],
        "depth_track_tick_font_sizes": [track.tick_font_size for track in tracks],
    }
    if any(table is not None for row in table_rows for table in row):
        result["depth_track_tables"] = table_rows
    if any(file_path is not None for row in file_rows for file_path in row):
        result["depth_track_files"] = file_rows
    if mode == "linear":
        result["depth_track_heights"] = [track.height for track in tracks]
    return result


def _base_options(options: _CommonOptions, *, record_count: int, mode: Literal["circular", "linear"]):
    features = options.features
    labels = options.labels
    color_table, color_table_file = _source(features.color_table, name="feature color table")
    default_colors, default_colors_file = _source(
        features.default_colors,
        name="default color table",
    )
    visibility_table, visibility_file = _source(features.visibility, name="feature visibility")
    whitelist_table, whitelist_file = _source(labels.whitelist, name="label whitelist")
    priority_table, priority_file = _source(
        labels.qualifier_priority,
        name="qualifier priority",
    )
    override_table, override_file = _source(labels.overrides, name="label overrides")
    values: dict[str, object] = {
        "config": options.config,
        "config_overrides": options.config_overrides,
        "colors": _ColorOptions(
            color_table=color_table,
            color_table_file=color_table_file,
            default_colors=default_colors,
            default_colors_palette=features.palette,
            default_colors_file=default_colors_file,
        ),
        "annotations": options.annotations,
        "selected_features_set": features.types,
        "feature_visibility_table": visibility_table,
        "feature_visibility_table_file": visibility_file,
        "label_whitelist_table": whitelist_table,
        "label_whitelist_file": whitelist_file,
        "qualifier_priority_table": priority_table,
        "qualifier_priority_file": priority_file,
        "label_override_table": override_table,
        "label_override_file": override_file,
        "feature_shapes": features.shapes,
        "dinucleotide": options.dinucleotide,
        "window": options.window,
        "step": options.step,
        "depth_window": options.depth_window,
        "depth_step": options.depth_step,
        "plot_title": options.title.text,
        "plot_title_font_size": options.title.font_size,
        "evalue": options.thresholds.evalue,
        "bitscore": options.thresholds.bitscore,
        "identity": options.thresholds.identity,
        "alignment_length": options.thresholds.alignment_length,
    }
    values.update(_depth_kwargs(options.depth_tracks, record_count=record_count, mode=mode))
    return values


def _circular_options(
    options: CircularOptions,
    *,
    record_count: int,
) -> _CircularDiagramOptions:
    values = _base_options(options, record_count=record_count, mode="circular")
    values["tracks"] = _CircularRequestTrackOptions(
        circular_track_slots=options.tracks.slots,
        circular_track_axis_index=options.tracks.axis_index,
        center_reserved_radius=options.tracks.center_reserved_radius,
    )
    values["output"] = _CircularOutputOptions(
        legend=options.legend,
        plot_title_position=options.title.position,
    )
    comparison_rings = options.comparison_rings
    assert comparison_rings is not None
    if comparison_rings.tracks:
        kinds = {
            "table" if isinstance(track.source, DataFrame) else "file"
            for track in comparison_rings.tracks
        }
        if len(kinds) != 1:
            raise ValidationError(
                "Comparison rings must use either DataFrames or paths in one "
                "diagram, not both."
            )
        sources = [track.source for track in comparison_rings.tracks]
        if "table" in kinds:
            values["conservation_dataframes"] = sources
        else:
            values["conservation_blast_files"] = [str(source) for source in sources]
        labels = [track.label for track in comparison_rings.tracks]
        colors = [track.color for track in comparison_rings.tracks]
        if any(label is not None for label in labels):
            if any(label is None for label in labels):
                raise ValidationError(
                    "Set labels for every comparison ring or for none of them."
                )
            values["conservation_labels"] = labels
        if any(color is not None for color in colors):
            if any(color is None for color in colors):
                raise ValidationError(
                    "Set colors for every comparison ring or for none of them."
                )
            values["conservation_colors"] = colors
    values.update(
        conservation_reference=comparison_rings.reference,
        conservation_ring_width=comparison_rings.ring_width,
        conservation_ring_gap=comparison_rings.ring_gap,
        species=options.species,
        strain=options.strain,
        keep_full_definition_with_plot_title=options.keep_full_definition_with_title,
    )
    return _CircularDiagramOptions(**values)


def _linear_options(
    options: LinearOptions,
    *,
    record_count: int,
) -> _LinearDiagramOptions:
    values = _base_options(options, record_count=record_count, mode="linear")
    values["tracks"] = _LinearRequestTrackOptions(
        linear_track_slots=options.tracks.slots,
        linear_track_axis_index=options.tracks.axis_index,
    )
    values["output"] = _LinearOutputOptions(
        legend=options.legend,
        plot_title_position=options.title.position,
    )
    comparisons = options.comparisons
    values.update(
        blast_files=comparisons.blast_files,
        linear_comparisons=comparisons.comparisons,
        protein_comparisons=comparisons.protein_comparisons,
        orthogroups=comparisons.orthogroups,
        protein_blastp_mode=comparisons.protein_mode,
        protein_comparison_pairs=comparisons.pairs,
        pairwise_match_style=comparisons.match_style,
        collinearity_blocks=comparisons.collinearity_blocks,
        collinearity_params=comparisons.collinearity_params,
        collinearity_unit_mode=comparisons.collinearity_unit,
        collinearity_anchor_mode=comparisons.collinearity_anchor,
        collinearity_search_scope=comparisons.collinearity_scope,
        collinearity_color_mode=comparisons.collinearity_color,
        losatp_bin=comparisons.losat_executable,
        ncbi_blastp_bin=comparisons.blastp_executable,
        losatp_threads=comparisons.threads,
        protein_blastp_max_hits=comparisons.max_hits,
        protein_blastp_candidate_limit=comparisons.candidate_limit,
        orthogroup_membership_mode=comparisons.orthogroup_membership,
        orthogroup_member_max_hits=comparisons.orthogroup_member_max_hits,
        collinear_max_paralog_links_per_orthogroup=comparisons.max_paralog_links,
        align_orthogroup_feature=comparisons.align_feature,
    )
    return _LinearDiagramOptions(**values)


def _interactive_context(
    prepared: _PreparedDiagramRequest,
    *,
    options: _CommonOptions,
) -> InteractiveSvgContext:
    comparison_sequence_records: list[list[SeqRecord]] = []
    if isinstance(options, CircularOptions):
        comparison_rings = options.comparison_rings
        assert comparison_rings is not None
        for track in comparison_rings.tracks:
            source = track.comparison_sequence_source
            if source is None:
                comparison_sequence_records.append([])
            elif isinstance(source, SeqRecord):
                comparison_sequence_records.append([source])
            elif isinstance(source, Sequence) and not isinstance(source, (str, bytes, PathLike)):
                comparison_sequence_records.append(list(source))
            else:
                comparison_sequence_records.append(list(SeqIO.parse(str(source), "fasta")))
    context = _build_prepared_interactive_context(
        prepared,
        comparison_sequence_records=comparison_sequence_records,
    )
    return context


def draw_circular(
    records: RecordCollection,
    *,
    options: CircularOptions | None = None,
    layout: CircularLayout | None = None,
) -> Diagram:
    """Draw one circular record or a multi-record circular grid."""

    normalized = _records(records)
    options = options or CircularOptions()
    if not isinstance(options, CircularOptions):
        raise ValidationError("draw_circular options must be CircularOptions.")
    if layout is not None and not isinstance(layout, CircularLayout):
        raise ValidationError("draw_circular layout must be CircularLayout.")
    compiled = _circular_options(options, record_count=len(normalized))
    prepared = _build_request_diagram(
        _CircularDiagramRequest(
            records=tuple(
                _RecordInput(source=_InMemoryRecordSource(record))
                for record in normalized
            ),
            options=compiled,
            layout=layout._legacy() if layout is not None else None,
        )
    )
    return Diagram(
        prepared.drawing,
        mode="circular",
        records=normalized,
        interactive_context=lambda: _interactive_context(
            prepared,
            options=options,
        ),
    )


def draw_linear(
    records: RecordCollection,
    *,
    options: LinearOptions | None = None,
    layout: LinearLayout | None = None,
) -> Diagram:
    """Draw one or more records as a linear diagram."""

    normalized = _records(records)
    options = options or LinearOptions()
    if not isinstance(options, LinearOptions):
        raise ValidationError("draw_linear options must be LinearOptions.")
    if layout is not None and not isinstance(layout, LinearLayout):
        raise ValidationError("draw_linear layout must be LinearLayout.")
    compiled = _linear_options(options, record_count=len(normalized))
    prepared = _build_request_diagram(
        _LinearDiagramRequest(
            records=tuple(
                _RecordInput(source=_InMemoryRecordSource(record))
                for record in normalized
            ),
            options=compiled,
            layout=layout._legacy() if layout is not None else None,
        )
    )
    return Diagram(
        prepared.drawing,
        mode="linear",
        records=normalized,
        interactive_context=lambda: _interactive_context(
            prepared,
            options=options,
        ),
    )


__all__ = [
    "CircularLayout",
    "CircularOptions",
    "CircularTrackOptions",
    "ComparisonRingOptions",
    "ComparisonRingTrackOptions",
    "ConservationOptions",
    "ConservationTrackOptions",
    "DepthTrackOptions",
    "Diagram",
    "FeatureOptions",
    "LabelOptions",
    "LinearComparisonOptions",
    "LinearLayout",
    "LinearOptions",
    "LinearTrackOptions",
    "Thresholds",
    "TitleOptions",
    "draw_circular",
    "draw_linear",
    "read_genbank",
    "read_gff",
]
