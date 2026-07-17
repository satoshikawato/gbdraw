"""Resolve public annotation targets to local half-open record segments."""

from __future__ import annotations

from collections.abc import Sequence

from Bio.SeqFeature import SeqFeature  # type: ignore[reportMissingImports]
from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from gbdraw.core.record_metadata import _read_coord_map
from gbdraw.exceptions import ValidationError
from gbdraw.features.selector_values import (
    get_feature_hash,
    get_feature_location_str,
    get_feature_record_location_str,
    get_feature_qualifiers,
    get_feature_type,
    normalize_qualifier_values,
)
from gbdraw.io.record_select import RecordSelector

from .io import annotation_sets_from_dataframe, read_annotation_table
from .models import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    FeatureSelector,
    FeatureSpan,
    ResolvedAnnotationBundle,
    ResolvedRegionAnnotation,
    ResolutionWarning,
)


def _bind_record(records: Sequence[SeqRecord], selector: RecordSelector | None) -> int:
    if selector is None:
        if len(records) != 1:
            raise ValidationError(
                "Annotation target without a record selector is ambiguous for multiple records."
            )
        return 0
    if selector.record_index is not None:
        index = int(selector.record_index)
        if index < 0 or index >= len(records):
            raise ValidationError(
                f"Annotation record selector {selector.label()} is out of range for {len(records)} records."
            )
        return index
    matches = [index for index, record in enumerate(records) if record.id == selector.record_id]
    if not matches:
        raise ValidationError(f"Annotation record selector {selector.record_id!r} did not match any record.")
    if len(matches) > 1:
        raise ValidationError(
            f"Annotation record selector {selector.record_id!r} matched multiple records; use #index."
        )
    return matches[0]


def _merge_segments(segments: Sequence[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    ordered = sorted((int(start), int(end)) for start, end in segments if int(end) > int(start))
    if not ordered:
        return ()
    merged: list[list[int]] = [[ordered[0][0], ordered[0][1]]]
    for start, end in ordered[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return tuple((start, end) for start, end in merged)


def _clip_local_interval(
    start_1: int,
    end_1: int,
    record_length: int,
) -> tuple[tuple[int, int], bool]:
    raw_start, raw_end = start_1 - 1, end_1
    start = max(0, raw_start)
    end = min(record_length, raw_end)
    return ((start, end) if end > start else (0, 0)), (start != raw_start or end != raw_end)


def _source_interval_to_local(
    source_start: int,
    source_end: int,
    *,
    coord_base: int,
    coord_step: int,
    record_length: int,
) -> tuple[tuple[int, int], bool]:
    if record_length <= 0:
        return (0, 0), True
    last_source = coord_base + coord_step * (record_length - 1)
    visible_low, visible_high = sorted((coord_base, last_source))
    low = max(int(source_start), int(visible_low))
    high = min(int(source_end), int(visible_high))
    clipped = low != source_start or high != source_end
    if high < low:
        return (0, 0), True
    first_local = (low - coord_base) * coord_step
    last_local = (high - coord_base) * coord_step
    return (min(first_local, last_local), max(first_local, last_local) + 1), clipped


def _coordinate_segments(
    target: CoordinateSpan,
    record: SeqRecord,
    *,
    mode: str,
) -> tuple[tuple[tuple[int, int], ...], bool]:
    length = len(record.seq)
    if target.wraps_origin and mode != "circular":
        raise ValidationError("Origin-spanning annotation targets are not supported in linear diagrams.")
    if target.coordinate_space == "local":
        ranges = (
            ((target.start, length), (1, target.end))
            if target.wraps_origin
            else ((target.start, target.end),)
        )
        resolved = [_clip_local_interval(start, end, length) for start, end in ranges]
    else:
        coord_base, coord_step = _read_coord_map(record)
        last_source = coord_base + coord_step * max(0, length - 1)
        source_high = max(coord_base, last_source, target.start, target.end)
        ranges = (
            ((target.start, source_high), (1, target.end))
            if target.wraps_origin
            else ((target.start, target.end),)
        )
        resolved = [
            _source_interval_to_local(
                start,
                end,
                coord_base=coord_base,
                coord_step=coord_step,
                record_length=length,
            )
            for start, end in ranges
        ]
    segments = tuple(segment for segment, _ in resolved if segment[1] > segment[0])
    return _merge_segments(segments), any(clipped for _, clipped in resolved)


def _seqfeature_segments(feature: object) -> tuple[tuple[int, int], ...]:
    location = getattr(feature, "location", None)
    if location is None:
        return ()
    parts = getattr(location, "parts", None) or (location,)
    segments: list[tuple[int, int]] = []
    for part in parts:
        try:
            start, end = int(part.start), int(part.end)
        except (TypeError, ValueError, AttributeError):
            continue
        if end > start:
            segments.append((start, end))
    return tuple(segments)


def _feature_matches(feature: object, selector: FeatureSelector, record_id: str) -> bool:
    key = selector.key.lower() if selector.key else None
    expected = selector.value
    if key == "hash":
        return get_feature_hash(feature, record_id) == expected
    if key in {"location", "position"}:
        return get_feature_location_str(feature) == expected
    if key == "record_location":
        return get_feature_record_location_str(feature, record_id) == expected
    if key in {"type", "feature_type"}:
        return get_feature_type(feature) == expected
    qualifiers = get_feature_qualifiers(feature)
    if key:
        return any(
            expected == value
            for qualifier_key, raw_values in qualifiers.items()
            if str(qualifier_key).lower() == key
            for value in normalize_qualifier_values(raw_values)
        )
    if get_feature_hash(feature, record_id) == expected:
        return True
    return any(
        expected == value
        for raw_values in qualifiers.values()
        for value in normalize_qualifier_values(raw_values)
    )


def _feature_segments(target: FeatureSpan, record: SeqRecord, *, mode: str) -> tuple[tuple[int, int], ...]:
    matched: list[object] = []
    unmatched: list[str] = []
    for selector in target.selectors:
        selector_matches = [
            feature
            for feature in getattr(record, "features", ())
            if isinstance(feature, SeqFeature) and _feature_matches(feature, selector, record.id)
        ]
        if not selector_matches:
            unmatched.append(
                f"{selector.key}={selector.value}" if selector.key else selector.value
            )
        matched.extend(selector_matches)
    if unmatched:
        raise ValidationError(
            f"Feature annotation selector(s) did not match record {record.id!r}: {', '.join(unmatched)}."
        )
    segments = _merge_segments(
        [segment for feature in matched for segment in _seqfeature_segments(feature)]
    )
    if not segments:
        return ()
    if target.envelope == "segments":
        return segments

    length = len(record.seq)
    normal = ((segments[0][0], segments[-1][1]),)
    if mode != "circular" or len(segments) == 1:
        return normal
    normal_span = normal[0][1] - normal[0][0]
    wrap_span = (length - segments[-1][0]) + segments[0][1]
    use_wrap = target.circular_path == "reverse" or (
        target.circular_path == "shortest" and wrap_span < normal_span
    )
    if target.circular_path == "forward":
        use_wrap = False
    return ((segments[-1][0], length), (0, segments[0][1])) if use_wrap else normal


def _midpoint(segments: Sequence[tuple[int, int]], record_length: int) -> float:
    if len(segments) == 2 and segments[0][0] == 0 and segments[-1][1] == record_length:
        ordered = (segments[-1], segments[0])
        span = sum(end - start for start, end in ordered)
        return (ordered[0][0] + span / 2.0) % max(1, record_length)
    total = sum(end - start for start, end in segments)
    cursor = total / 2.0
    for start, end in segments:
        width = end - start
        if cursor <= width:
            return start + cursor
        cursor -= width
    return float(segments[-1][1])


def resolve_annotation_set(
    annotation_set: AnnotationSet,
    records: Sequence[SeqRecord],
    *,
    mode: str,
) -> ResolvedAnnotationBundle:
    """Resolve one set against already materialized records."""

    if mode not in {"circular", "linear"}:
        raise ValidationError("Annotation resolver mode must be 'circular' or 'linear'.")
    resolved: list[ResolvedRegionAnnotation] = []
    warnings: list[ResolutionWarning] = []
    for annotation in annotation_set.annotations:
        record_index = _bind_record(records, annotation.target.record)
        record = records[record_index]
        clipped = False
        if isinstance(annotation.target, CoordinateSpan):
            segments, clipped = _coordinate_segments(annotation.target, record, mode=mode)
            policy = annotation.target.out_of_bounds
            if clipped and policy == "error":
                raise ValidationError(
                    f"Annotation {annotation_set.id}/{annotation.id} extends outside record {record.id!r}."
                )
            if clipped and policy == "skip":
                warnings.append(
                    ResolutionWarning(
                        code="out_of_bounds_skipped",
                        set_id=annotation_set.id,
                        annotation_id=annotation.id,
                        message=f"Skipped annotation outside record {record.id!r}.",
                    )
                )
                continue
            if clipped:
                warnings.append(
                    ResolutionWarning(
                        code="out_of_bounds_clipped",
                        set_id=annotation_set.id,
                        annotation_id=annotation.id,
                        message=f"Clipped annotation to record {record.id!r}.",
                    )
                )
        else:
            segments = _feature_segments(annotation.target, record, mode=mode)
        if not segments:
            warnings.append(
                ResolutionWarning(
                    code="empty_span",
                    set_id=annotation_set.id,
                    annotation_id=annotation.id,
                    message=f"Annotation resolved to an empty span on record {record.id!r}.",
                )
            )
            continue
        style = annotation.style or annotation_set.default_style
        if style.line_cap == "arrow" and isinstance(annotation.target, CoordinateSpan):
            raise ValidationError(
                f"Annotation {annotation_set.id}/{annotation.id} uses line_cap=arrow without a directional feature target."
            )
        span_bp = sum(end - start for start, end in segments)
        resolved.append(
            ResolvedRegionAnnotation(
                id=annotation.id,
                set_id=annotation_set.id,
                record_index=record_index,
                segments=tuple(segments),
                midpoint_bp=_midpoint(segments, len(record.seq)),
                span_bp=span_bp,
                label=annotation.label,
                mark=annotation.mark,
                lane=annotation.lane,
                style=style,
                style_is_annotation_override=annotation.style is not None,
                legend_label=annotation.legend_label or annotation_set.legend_label,
                metadata=annotation.metadata,
            )
        )
    return ResolvedAnnotationBundle(tuple(resolved), tuple(warnings), (annotation_set.id,))


def _materialize_sets(options: AnnotationOptions, *, mode: str) -> tuple[AnnotationSet, ...]:
    if options.sets:
        return options.sets
    if options.table is not None:
        return annotation_sets_from_dataframe(options.table, mode=mode)  # type: ignore[arg-type]
    if options.table_file:
        return read_annotation_table(options.table_file, mode=mode)  # type: ignore[arg-type]
    return ()


def resolve_annotations(
    annotations: AnnotationOptions | Sequence[AnnotationSet] | None,
    records: Sequence[SeqRecord],
    *,
    mode: str,
) -> ResolvedAnnotationBundle:
    """Resolve all configured sets into one deterministic bundle."""

    if annotations is None:
        return ResolvedAnnotationBundle(())
    sets = (
        _materialize_sets(annotations, mode=mode)
        if isinstance(annotations, AnnotationOptions)
        else tuple(annotations)
    )
    resolved = [resolve_annotation_set(item, records, mode=mode) for item in sets]
    return ResolvedAnnotationBundle(
        tuple(annotation for bundle in resolved for annotation in bundle.annotations),
        tuple(warning for bundle in resolved for warning in bundle.warnings),
        tuple(item.id for item in sets),
    )


__all__ = ["resolve_annotation_set", "resolve_annotations"]
