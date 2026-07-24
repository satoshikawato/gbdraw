"""CLI-independent builders for rich interactive SVG metadata."""

from __future__ import annotations

from typing import Any, Sequence

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.visibility import compile_feature_visibility_rules
from gbdraw.exceptions import ValidationError
from gbdraw.render.interactive_svg import InteractiveSvgContext
from gbdraw.annotations import AnnotationOptions, resolve_annotations
from gbdraw.web_support.feature_metadata import extract_features_from_records_payload
from gbdraw.web_support.orthogroup_metadata import (
    enrich_features_with_orthogroups,
    serialize_orthogroups_payload,
)


def build_interactive_svg_context(
    records: Sequence[SeqRecord],
    *,
    selected_features_set: Sequence[str] | None = None,
    feature_table: DataFrame | None = None,
    feature_visibility_rules: list[dict[str, Any]] | None = None,
    color_table: DataFrame | None = None,
    default_colors: DataFrame | None = None,
    orthogroups: object | None = None,
    linear_rendered_feature_ids: bool = False,
    annotations: AnnotationOptions | None = None,
    mode: str | None = None,
    comparison_sequence_records: Sequence[Sequence[SeqRecord]] | None = None,
) -> InteractiveSvgContext:
    """Build rich popup metadata from rendered records.

    ``linear_rendered_feature_ids`` identifies the diagram mode, but feature
    metadata remains keyed by stable biological IDs. The final SVG enricher maps
    those IDs to record-specific rendered IDs after inspecting the actual DOM.
    Passing precompiled visibility rules is useful to CLI callers; public callers
    normally pass ``feature_table``.
    """

    if feature_table is not None and feature_visibility_rules is not None:
        raise ValidationError(
            "Pass either feature_table or feature_visibility_rules, not both."
        )
    resolved_visibility_rules = (
        feature_visibility_rules
        if feature_visibility_rules is not None
        else compile_feature_visibility_rules(feature_table)
    )
    specific_color_rules = None
    if color_table is not None and default_colors is not None:
        specific_color_rules, _ = preprocess_color_tables(color_table, default_colors)

    record_list = list(records)
    payload = extract_features_from_records_payload(
        record_list,
        selected_features=selected_features_set,
        feature_visibility_rules=resolved_visibility_rules,
        specific_color_rules=specific_color_rules,
        linear_rendered_feature_ids=linear_rendered_feature_ids,
        include_biological_features=True,
    )
    features = payload.get("features", [])
    biological_features = payload.get("biological_features", [])
    orthogroup_payload: list[dict[str, object]] = []
    if orthogroups is not None:
        orthogroup_payload = serialize_orthogroups_payload(
            orthogroups,
            records=record_list,
        )
        if orthogroup_payload:
            features = enrich_features_with_orthogroups(features, orthogroup_payload)
            biological_features = enrich_features_with_orthogroups(
                biological_features,
                orthogroup_payload,
            )

    annotation_payload: list[dict[str, object]] = []
    if annotations is not None:
        resolved = resolve_annotations(
            annotations,
            record_list,
            mode=mode or ("linear" if linear_rendered_feature_ids else "circular"),
        )
        annotation_payload = [
            {
                "id": item.id,
                "set_id": item.set_id,
                "track_id": "",
                "record_index": item.record_index,
                "segments": [list(segment) for segment in item.segments],
                "midpoint_bp": item.midpoint_bp,
                "span_bp": item.span_bp,
                "label": item.label,
                "mark": item.mark,
                "lane": item.lane,
                "legend_label": item.legend_label,
                "metadata": dict(item.metadata),
            }
            for item in resolved.annotations
        ]

    resolved_mode = mode or ("linear" if linear_rendered_feature_ids else "circular")
    sequence_sources: list[dict[str, object]] = []
    for record_index, record in enumerate(record_list):
        aliases = {
            str(record.id or "").strip(),
            str(record.name or "").strip(),
        }
        accessions = (getattr(record, "annotations", {}) or {}).get("accessions", [])
        if isinstance(accessions, str):
            aliases.add(accessions.strip())
        else:
            aliases.update(str(value or "").strip() for value in accessions)
        sequence_sources.append(
            {
                "key": f"linear:record:{record_index}" if resolved_mode == "linear" else f"circular:record:{record_index}",
                "recordId": str(record.id),
                "aliases": sorted(value for value in aliases if value),
                "sequence": str(record.seq).upper(),
                "origin": "linear-record" if resolved_mode == "linear" else "circular-reference",
                "recordIndex": record_index,
            }
        )
    for source_index, source_records in enumerate(comparison_sequence_records or ()):
        for record in source_records:
            sequence_sources.append(
                {
                    "key": f"homology:comparison:{source_index}:{record.id}",
                    "recordId": str(record.id),
                    "aliases": [str(record.id), str(record.name)],
                    "sequence": str(record.seq).upper(),
                    "origin": "homology-comparison",
                    "sourceIndex": source_index,
                }
            )

    return InteractiveSvgContext(
        features=features,
        biological_features=biological_features,
        orthogroups=orthogroup_payload,
        annotations=annotation_payload,
        sequence_sources=sequence_sources,
    )


__all__ = ["build_interactive_svg_context"]
