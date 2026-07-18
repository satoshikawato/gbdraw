"""Normalize and render typed diagram requests without CLI/session orchestration."""

from __future__ import annotations

import copy
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal, Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite import Drawing  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError
from gbdraw.features.visibility import (
    read_feature_visibility_file,
    resolve_candidate_feature_types,
)
from gbdraw.io.colors import load_default_colors, read_color_table
from gbdraw.io.record_select import reverse_records, select_record
from gbdraw.render.interactive_context import build_interactive_svg_context
from gbdraw.render.interactive_svg import InteractiveSvgContext

from .diagram import (
    DEFAULT_SELECTED_FEATURES,
    build_circular_diagram,
    build_circular_multi_diagram,
    build_linear_diagram,
)
from .io import apply_region_specs, load_gbks, load_gff_fasta
from .options import CircularMultiRecordOptions, DiagramOptions, LinearMultiRecordOptions
from .render import save_figure_to
from .requests import (
    CircularDiagramRequest,
    DiagramRequest,
    GenBankInputSource,
    GffFastaInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
)


@dataclass(frozen=True)
class PreparedDiagramRequest:
    """A validated request with normalized records and its SVG drawing."""

    mode: Literal["circular", "linear"]
    request: DiagramRequest
    records: tuple[SeqRecord, ...]
    drawing: Drawing


@dataclass(frozen=True)
class RequestRenderResult:
    """Files and normalized inputs produced by one request render."""

    mode: Literal["circular", "linear"]
    request: DiagramRequest
    records: tuple[SeqRecord, ...]
    drawing: Drawing
    output_paths: tuple[Path, ...]
    warnings: tuple[str, ...] = ()


def _mode(request: DiagramRequest) -> Literal["circular", "linear"]:
    return "circular" if isinstance(request, CircularDiagramRequest) else "linear"


def _load_record_input(
    record_input: RecordInput,
    *,
    mode: Literal["circular", "linear"],
    gff_candidate_features: Sequence[str] | None,
    gff_keep_all_features: bool,
    color_table: DataFrame | None,
    feature_visibility_table: DataFrame | None,
) -> SeqRecord:
    selector = record_input.selector
    selector_values = [selector.raw] if selector is not None else None
    reverse = record_input.presentation.reverse_complement
    source = record_input.source

    if isinstance(source, GenBankInputSource):
        records = load_gbks(
            [str(source.path)],
            mode=mode,
            record_selectors=selector_values,
            reverse_flags=[reverse],
        )
    elif isinstance(source, GffFastaInputSource):
        records = load_gff_fasta(
            [str(source.gff_path)],
            [str(source.fasta_path)],
            mode=mode,
            selected_features_set=gff_candidate_features,
            keep_all_features=gff_keep_all_features,
            record_selectors=selector_values,
            reverse_flags=[reverse],
            color_table=color_table,
            feature_visibility_table=feature_visibility_table,
        )
    elif isinstance(source, InMemoryRecordSource):
        records = [copy.deepcopy(source.record)]
        try:
            records = select_record(records, selector)
            records = reverse_records(records, reverse)
        except ValueError as exc:
            raise ValidationError(str(exc)) from exc
    else:  # pragma: no cover - RecordInput validates this union.
        raise ValidationError("Unsupported record input source.")

    if record_input.region is not None:
        records = apply_region_specs(records, [record_input.region])
    if len(records) != 1:
        raise ValidationError(
            "Each RecordInput must resolve to exactly one record; add a selector or region."
        )

    record = records[0]
    presentation = record_input.presentation
    if getattr(record, "annotations", None) is None:
        record.annotations = {}
    if presentation.label:
        record.annotations["gbdraw_record_label"] = presentation.label
    if presentation.subtitle:
        record.annotations["gbdraw_record_subtitle"] = presentation.subtitle
    if record_input.record_key:
        record.annotations["gbdraw_record_key"] = record_input.record_key
    return record


def normalize_request_records(request: DiagramRequest) -> tuple[SeqRecord, ...]:
    """Load/copy every RecordInput and return exactly one record per input."""

    if not isinstance(request, (CircularDiagramRequest, LinearDiagramRequest)):
        raise ValidationError("Unsupported diagram request type.")
    mode = _mode(request)
    has_gff_source = any(
        isinstance(record_input.source, GffFastaInputSource)
        for record_input in request.records
    )
    color_table = _color_table(request.options) if has_gff_source else None
    feature_visibility_table = (
        _visibility_table(request.options) if has_gff_source else None
    )
    candidate_features, keep_all_features = (
        resolve_candidate_feature_types(
            request.options.selected_features_set or DEFAULT_SELECTED_FEATURES,
            color_table=color_table,
            feature_visibility_table=feature_visibility_table,
        )
        if has_gff_source
        else (set(request.options.selected_features_set or DEFAULT_SELECTED_FEATURES), False)
    )
    return tuple(
        _load_record_input(
            record_input,
            mode=mode,
            gff_candidate_features=tuple(sorted(candidate_features)),
            gff_keep_all_features=keep_all_features,
            color_table=color_table,
            feature_visibility_table=feature_visibility_table,
        )
        for record_input in request.records
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


def build_request_diagram(request: DiagramRequest) -> PreparedDiagramRequest:
    """Normalize inputs and build a drawing through the high-level API owners."""

    records = normalize_request_records(request)
    if isinstance(request, CircularDiagramRequest):
        layout = _layout_with_record_placements(request)
        if layout is None:
            drawing = build_circular_diagram(records[0], options=request.options)
        else:
            drawing = build_circular_multi_diagram(
                records,
                options=request.options,
                layout=layout,
            )
        mode: Literal["circular", "linear"] = "circular"
    elif isinstance(request, LinearDiagramRequest):
        linear_layout = _linear_layout_with_record_placements(request)
        if linear_layout is None:
            drawing = build_linear_diagram(records, options=request.options)
        else:
            drawing = build_linear_diagram(
                records,
                options=request.options,
                layout=linear_layout,
            )
        mode = "linear"
    else:  # pragma: no cover - normalize_request_records rejects this first.
        raise ValidationError("Unsupported diagram request type.")
    return PreparedDiagramRequest(
        mode=mode,
        request=request,
        records=records,
        drawing=drawing,
    )


def _visibility_table(options: DiagramOptions) -> DataFrame | None:
    table = (
        options.feature_visibility_table
        if options.feature_visibility_table is not None
        else options.feature_table
    )
    file_path = (
        options.feature_visibility_table_file
        if options.feature_visibility_table_file is not None
        else options.feature_table_file
    )
    if table is None and file_path is not None:
        return read_feature_visibility_file(file_path)
    return table


def _color_table(options: DiagramOptions) -> DataFrame | None:
    colors = options.colors
    if colors is None or colors.color_table is not None:
        return colors.color_table if colors is not None else None
    if colors.color_table_file is not None:
        return read_color_table(colors.color_table_file)
    return None


def _interactive_context(
    prepared: PreparedDiagramRequest,
) -> InteractiveSvgContext | None:
    output = prepared.request.output
    if "interactive_svg" not in output.formats or output.interactive_metadata_policy == "omit":
        return None

    options = prepared.request.options
    colors = options.colors
    color_table = _color_table(options)
    default_colors = colors.default_colors if colors is not None else None
    if default_colors is None:
        default_colors = load_default_colors(
            colors.default_colors_file if colors is not None and colors.default_colors_file else "",
            colors.default_colors_palette if colors is not None else "default",
        )
    return build_interactive_svg_context(
        prepared.records,
        selected_features_set=options.selected_features_set,
        feature_table=_visibility_table(options),
        color_table=color_table,
        default_colors=default_colors,
        orthogroups=options.orthogroups,
        linear_rendered_feature_ids=prepared.mode == "linear",
        annotations=options.annotations,
        mode=prepared.mode,
    )


def render_request(request: DiagramRequest) -> RequestRenderResult:
    """Build and save one typed request, returning only paths that were created."""

    prepared = build_request_diagram(request)
    output = request.output
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
        interactive_context=_interactive_context(prepared),
    )
    return RequestRenderResult(
        mode=prepared.mode,
        request=prepared.request,
        records=prepared.records,
        drawing=prepared.drawing,
        output_paths=tuple(Path(path) for path in paths),
    )


__all__ = [
    "PreparedDiagramRequest",
    "RequestRenderResult",
    "build_request_diagram",
    "normalize_request_records",
    "render_request",
]
