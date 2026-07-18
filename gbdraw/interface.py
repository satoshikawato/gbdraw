"""Small, typed, beginner-facing Python API for gbdraw."""

from __future__ import annotations

from dataclasses import dataclass, field
from os import PathLike
from pathlib import Path
from typing import Literal, Mapping, Sequence, TypeAlias

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]

from gbdraw.analysis.collinearity import (
    CollinearityAnchorMode,
    CollinearityBlock,
    CollinearityColorMode,
    CollinearityParameters,
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
from gbdraw.api.diagram import (
    build_circular_diagram as _build_circular_diagram,
    build_circular_multi_diagram as _build_circular_multi_diagram,
    build_linear_diagram as _build_linear_diagram,
)
from gbdraw.api.io import load_gbks as _load_gbks, load_gff_fasta as _load_gff_fasta
from gbdraw.api.options import (
    CircularMultiRecordOptions as _CircularLayout,
    ColorOptions as _ColorOptions,
    DiagramOptions as _DiagramOptions,
    LinearMultiRecordOptions as _LinearLayout,
    OutputOptions as _OutputOptions,
    TrackOptions as _TrackOptions,
)
from gbdraw.api.render import render_to_bytes
from gbdraw.exceptions import ExportError, ValidationError
from gbdraw.features.visibility import read_feature_visibility_file
from gbdraw.io.colors import load_default_colors, read_color_table
from gbdraw.linear_comparison import LinearComparison
from gbdraw.config.models import GbdrawConfig
from gbdraw.render.interactive_context import build_interactive_svg_context
from gbdraw.render.interactive_svg import InteractiveSvgContext
from gbdraw.tracks import CircularTrackSlot, LinearTrackSlot


TableSource: TypeAlias = DataFrame | str | PathLike[str]
RecordCollection: TypeAlias = SeqRecord | Sequence[SeqRecord]


@dataclass(frozen=True)
class FeatureOptions:
    """Feature selection, colors, visibility, and shape overrides."""

    types: Sequence[str] | None = None
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
    """Filters shared by comparison and conservation inputs."""

    evalue: float = 1e-5
    bitscore: float = 50.0
    identity: float = 70.0
    alignment_length: int = 0


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


@dataclass(frozen=True)
class CircularTrackOptions:
    """Circular track order and axis placement."""

    slots: Sequence[str | CircularTrackSlot] | None = None
    axis_index: int | None = None
    center_reserved_radius: float | None = None


@dataclass(frozen=True)
class LinearTrackOptions:
    """Linear track order and axis placement."""

    slots: Sequence[str | LinearTrackSlot] | None = None
    axis_index: int | None = None


@dataclass(frozen=True)
class CircularLayout:
    """Grid layout used automatically when a circular diagram has multiple records."""

    size: Literal["linear", "auto", "equal", "sqrt"] = "auto"
    min_radius_ratio: float = 0.55
    column_gap_ratio: float = 0.10
    row_gap_ratio: float = 0.05
    positions: Sequence[str] | None = None

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

    def _legacy(self) -> _LinearLayout:
        return _LinearLayout(
            record_gap_px=self.record_gap,
            multi_record_positions=self.positions,
        )


@dataclass(frozen=True)
class ConservationTrackOptions:
    """One circular sequence-conservation ring."""

    source: TableSource
    label: str | None = None
    color: str | None = None


@dataclass(frozen=True)
class ConservationOptions:
    """Circular conservation rings and their shared geometry."""

    tracks: Sequence[ConservationTrackOptions] = ()
    reference: Literal["query", "subject", "auto"] = "auto"
    ring_width: float | None = None
    ring_gap: float | None = None


@dataclass(frozen=True)
class LinearComparisonOptions:
    """Precomputed or in-process comparison inputs for a linear diagram."""

    blast_files: Sequence[str] | None = None
    comparisons: Sequence[LinearComparison] | None = None
    protein_comparisons: Sequence[DataFrame] | None = None
    orthogroups: OrthogroupResult | None = None
    protein_mode: Literal["none", "pairwise", "orthogroup", "collinear"] = "none"
    pairs: Sequence[tuple[int, int]] | None = None
    match_style: Literal["ribbon", "curve"] = "ribbon"
    collinearity_blocks: CollinearityResult | Sequence[CollinearityBlock] | None = None
    collinearity_params: CollinearityParameters | LosslessCollinearityParameters | None = None
    collinearity_unit: CollinearityUnitMode | str = "auto"
    collinearity_anchor: CollinearityAnchorMode | str = "rbh"
    collinearity_scope: CollinearitySearchScope | str = "adjacent"
    collinearity_color: CollinearityColorMode | str = "orientation"
    losat_executable: str = "losat"
    blastp_executable: str | None = None
    threads: int | None = None
    max_hits: int = 5
    candidate_limit: int | None = None
    orthogroup_membership: OrthogroupMembershipMode | str = "anchor_core_v1"
    orthogroup_member_max_hits: int = 5
    max_paralog_links: int = 2
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


@dataclass(frozen=True)
class CircularOptions(_CommonOptions):
    """Options accepted only by :func:`draw_circular`."""

    tracks: CircularTrackOptions = field(default_factory=CircularTrackOptions)
    conservation: ConservationOptions = field(default_factory=ConservationOptions)
    species: str | None = None
    strain: str | None = None
    keep_full_definition_with_title: bool = False


@dataclass(frozen=True)
class LinearOptions(_CommonOptions):
    """Options accepted only by :func:`draw_linear`."""

    tracks: LinearTrackOptions = field(default_factory=LinearTrackOptions)
    comparisons: LinearComparisonOptions = field(default_factory=LinearComparisonOptions)


class Diagram:
    """A rendered diagram that can be serialized or saved without SVG internals."""

    def __init__(
        self,
        drawing: Drawing,
        *,
        mode: Literal["circular", "linear"],
        records: Sequence[SeqRecord],
        interactive_context: InteractiveSvgContext | None,
    ) -> None:
        self._drawing = drawing
        self.mode = mode
        self.records = tuple(records)
        self._interactive_context = interactive_context

    def to_bytes(self, format: str = "svg") -> bytes:
        """Return the diagram in SVG or a supported binary format."""

        return render_to_bytes(
            self._drawing,
            format,
            interactive_context=self._interactive_context,
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
        if output_path.exists() and not overwrite:
            raise ValidationError(
                f"Output file already exists: {output_path}. Use overwrite=True to replace it."
            )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_bytes(resolved_format)
        try:
            output_path.write_bytes(payload)
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

    return _load_gbks(_paths(paths), mode="linear")


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
        mode="linear",
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
        "output": _OutputOptions(
            legend=options.legend,
            plot_title_position=options.title.position,
        ),
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


def _circular_options(options: CircularOptions, *, record_count: int) -> _DiagramOptions:
    values = _base_options(options, record_count=record_count, mode="circular")
    values["tracks"] = _TrackOptions(
        circular_track_slots=options.tracks.slots,
        circular_track_axis_index=options.tracks.axis_index,
        center_reserved_radius=options.tracks.center_reserved_radius,
    )
    conservation = options.conservation
    if conservation.tracks:
        kinds = {
            "table" if isinstance(track.source, DataFrame) else "file"
            for track in conservation.tracks
        }
        if len(kinds) != 1:
            raise ValidationError(
                "Conservation tracks must use either DataFrames or paths in one diagram, not both."
            )
        sources = [track.source for track in conservation.tracks]
        if "table" in kinds:
            values["conservation_dataframes"] = sources
        else:
            values["conservation_blast_files"] = [str(source) for source in sources]
        labels = [track.label for track in conservation.tracks]
        colors = [track.color for track in conservation.tracks]
        if any(label is not None for label in labels):
            if any(label is None for label in labels):
                raise ValidationError("Set labels for every conservation track or for none of them.")
            values["conservation_labels"] = labels
        if any(color is not None for color in colors):
            if any(color is None for color in colors):
                raise ValidationError("Set colors for every conservation track or for none of them.")
            values["conservation_colors"] = colors
    values.update(
        conservation_reference=conservation.reference,
        conservation_ring_width=conservation.ring_width,
        conservation_ring_gap=conservation.ring_gap,
        species=options.species,
        strain=options.strain,
        keep_full_definition_with_plot_title=options.keep_full_definition_with_title,
    )
    return _DiagramOptions(**values)


def _linear_options(options: LinearOptions, *, record_count: int) -> _DiagramOptions:
    values = _base_options(options, record_count=record_count, mode="linear")
    values["tracks"] = _TrackOptions(
        linear_track_slots=options.tracks.slots,
        linear_track_axis_index=options.tracks.axis_index,
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
    return _DiagramOptions(**values)


def _interactive_context(
    records: Sequence[SeqRecord],
    *,
    options: _CommonOptions,
    legacy: _DiagramOptions,
    mode: Literal["circular", "linear"],
) -> InteractiveSvgContext:
    visibility_table, visibility_file = _source(
        options.features.visibility,
        name="feature visibility",
    )
    if visibility_table is None and visibility_file is not None:
        visibility_table = read_feature_visibility_file(visibility_file)
    color_table, color_file = _source(options.features.color_table, name="feature color table")
    if color_table is None and color_file is not None:
        color_table = read_color_table(color_file)
    default_colors, default_file = _source(
        options.features.default_colors,
        name="default color table",
    )
    if default_colors is None:
        default_colors = load_default_colors(default_file or "", options.features.palette)
    return build_interactive_svg_context(
        records,
        selected_features_set=options.features.types,
        feature_table=visibility_table,
        color_table=color_table,
        default_colors=default_colors,
        orthogroups=legacy.orthogroups,
        linear_rendered_feature_ids=mode == "linear",
        annotations=options.annotations,
        mode=mode,
    )


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
    if len(normalized) == 1 and layout is not None:
        raise ValidationError("CircularLayout requires at least two records.")
    legacy = _circular_options(options, record_count=len(normalized))
    if len(normalized) == 1:
        drawing = _build_circular_diagram(normalized[0], options=legacy)
    else:
        drawing = _build_circular_multi_diagram(
            normalized,
            options=legacy,
            layout=(layout or CircularLayout())._legacy(),
        )
    return Diagram(
        drawing,
        mode="circular",
        records=normalized,
        interactive_context=_interactive_context(
            normalized,
            options=options,
            legacy=legacy,
            mode="circular",
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
    legacy = _linear_options(options, record_count=len(normalized))
    drawing = _build_linear_diagram(
        normalized,
        options=legacy,
        layout=layout._legacy() if layout is not None else None,
    )
    return Diagram(
        drawing,
        mode="linear",
        records=normalized,
        interactive_context=_interactive_context(
            normalized,
            options=options,
            legacy=legacy,
            mode="linear",
        ),
    )


__all__ = [
    "CircularLayout",
    "CircularOptions",
    "CircularTrackOptions",
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
