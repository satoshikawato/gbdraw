"""Resolve public annotation targets to local half-open record segments."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal, cast

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.core.record_metadata import _read_coord_map, _source_to_local_index, record_identity_key
from gbdraw.exceptions import ValidationError
from gbdraw.features.selector_values import (
    get_feature_hash,
    get_feature_location_str,
    get_feature_record_location_str,
    get_feature_qualifiers,
    get_feature_type,
    matches_feature_selector,
)
from gbdraw.io.record_select import RecordSelector, drawn_record_units, select_record_units
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.layout.record_transform import bind_record_transforms, spans_may_cross_origin

from .io import annotation_sets_from_dataframe, read_annotation_table
from .models import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    FeatureIdentitySpan,
    FeatureSelector,
    FeatureSpan,
    ResolvedAnnotationBundle,
    ResolvedRegionAnnotation,
    ResolutionWarning,
    region_run,
)

if TYPE_CHECKING:
    from gbdraw.features.source import IdentityBinding, IdentityBindings

# Planner bindings of FeatureIdentitySpan targets, keyed by (record key, feature ID), one per unit that draws it.
FeatureBindings = Mapping[tuple[str, str], "IdentityBindings"]


def _bind_record_key(records: Sequence[SeqRecord], record_key: str) -> tuple[int, ...]:
    """The drawn units of the record ``record_key`` names: all windows of a record with several (D-01)."""

    units = tuple(index for index, record in enumerate(records) if record_identity_key(record) == record_key)
    if units:
        return units
    raise ValidationError(
        f"Unknown feature identity record key {record_key!r}.",
        diagnostic={"code": "FEATURE_IDENTITY"},
    )


def _bind_record(records: Sequence[SeqRecord], selector: RecordSelector | None) -> tuple[int, ...]:
    """The drawn units of the record ``selector`` names; ``#index`` counts records (``select_record_units``)."""

    if selector is None and len(drawn_record_units(records)) != 1:
        raise ValidationError(
            "Annotation target without a record selector is ambiguous for multiple records."
        )
    try:
        return select_record_units(records, selector)
    except ValueError as exc:
        raise ValidationError(f"Annotation record selector: {exc}") from exc


def merge_annotation_segments(
    segments: Sequence[tuple[int, int]],
) -> tuple[tuple[int, int], ...]:
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
    first_local = _source_to_local_index(low, coord_base, coord_step)
    last_local = _source_to_local_index(high, coord_base, coord_step)
    return (min(first_local, last_local), max(first_local, last_local) + 1), clipped


def _coordinate_segments(
    target: CoordinateSpan,
    record: SeqRecord,
    *,
    is_circular: bool,
) -> tuple[tuple[tuple[int, int], ...], bool]:
    length = len(record.seq)
    if target.wraps_origin and not is_circular:
        raise ValidationError("Origin-spanning annotation targets require a complete circular record.")
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
    return merge_annotation_segments(segments), any(clipped for _, clipped in resolved)


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
    return matches_feature_selector(
        key=selector.key, value=selector.value,
        feature_type=get_feature_type(feature),
        feature_hash=get_feature_hash(feature, record_id),
        location=get_feature_location_str(feature),
        record_location=get_feature_record_location_str(feature, record_id),
        qualifiers=get_feature_qualifiers(feature),
    )


@dataclass(frozen=True)
class _FeatureSegments:
    segments: tuple[tuple[int, int], ...]
    missing_count: int = 0


def _feature_segments(target: FeatureSpan, record: SeqRecord, *, is_circular: bool) -> _FeatureSegments:
    matched: list[object] = []
    unmatched: list[FeatureSelector] = []
    for selector in dict.fromkeys(target.selectors):
        selector_matches = [
            feature
            for feature in getattr(record, "features", ())
            if isinstance(feature, SeqFeature) and _feature_matches(feature, selector, record.id)
        ]
        if not selector_matches:
            unmatched.append(selector)
        matched.extend(selector_matches)
    if unmatched:
        return _FeatureSegments((), len(unmatched))
    return _envelope_segments(matched, target, record, is_circular=is_circular)


def _identity_segments(
    target: FeatureIdentitySpan,
    record: SeqRecord,
    binding: IdentityBinding | None,
    *,
    is_circular: bool,
) -> _FeatureSegments:
    """The bound feature's span; a feature that is not drawn is one unmatched selector."""
    if binding is None or binding.status != "present":
        return _FeatureSegments((), 1)
    return _envelope_segments((binding.feature,), target, record, is_circular=is_circular)


def _envelope_segments(
    matched: Sequence[object],
    target: FeatureSpan | FeatureIdentitySpan,
    record: SeqRecord,
    *,
    is_circular: bool,
) -> _FeatureSegments:
    segments = merge_annotation_segments(
        [segment for feature in matched for segment in _seqfeature_segments(feature)]
    )
    if not segments:
        return _FeatureSegments(())
    if target.envelope == "segments":
        return _FeatureSegments(segments)

    length = len(record.seq)
    normal = ((segments[0][0], segments[-1][1]),)
    if not is_circular or len(segments) == 1:
        return _FeatureSegments(normal)
    normal_span = normal[0][1] - normal[0][0]
    wrap_span = (length - segments[-1][0]) + segments[0][1]
    use_wrap = target.circular_path == "reverse" or (
        target.circular_path == "shortest" and wrap_span < normal_span
    )
    if target.circular_path == "forward":
        use_wrap = False
    return _FeatureSegments(((segments[-1][0], length), (0, segments[0][1])) if use_wrap else normal)


def segments_cross_origin(segments: Sequence[tuple[int, int]], record_length: int, *, is_circular: bool) -> bool:
    """Whether a span's segments go on across the origin: on a record whose spans may cross it, one segment starts
    at 0 and another ends at the record's end, so the two are one stretch of the circle."""

    return (
        is_circular
        and len(segments) > 1
        and min(start for start, _end in segments) == 0
        and max(end for _start, end in segments) == record_length
    )


def annotation_midpoint(
    segments: Sequence[tuple[int, int]],
    record_length: int,
    *,
    crosses_origin: bool,
) -> float:
    """Halfway along a span's segments in the order the region runs (``region_run``), on its record."""

    walk = region_run(segments, crosses_origin=crosses_origin)
    cursor = sum(end - start for start, end in walk) / 2.0
    for start, end in walk:
        if cursor <= end - start:
            return (start + cursor) % max(1, record_length)
        cursor -= end - start
    return float(walk[-1][1]) % max(1, record_length)


@dataclass(frozen=True)
class _UnitPiece:
    """An annotation target resolved on one drawn unit of its record."""

    record_index: int
    segments: tuple[tuple[int, int], ...]
    clipped: bool
    missing_count: int
    is_circular: bool

    @property
    def span_bp(self) -> int:
        return sum(end - start for start, end in self.segments)


def _unit_piece(
    target: CoordinateSpan | FeatureSpan | FeatureIdentitySpan,
    record: SeqRecord,
    record_index: int,
    binding: IdentityBinding | None,
    transform: RecordDisplayTransform,
) -> _UnitPiece:
    is_circular = spans_may_cross_origin(record, transform)
    if isinstance(target, CoordinateSpan):
        segments, clipped = _coordinate_segments(target, record, is_circular=is_circular)
        return _UnitPiece(record_index, tuple(segments), clipped, 0, is_circular)
    feature_segments = (
        _identity_segments(target, record, binding, is_circular=is_circular)
        if isinstance(target, FeatureIdentitySpan)
        else _feature_segments(target, record, is_circular=is_circular)
    )
    return _UnitPiece(
        record_index, tuple(feature_segments.segments), False, feature_segments.missing_count, is_circular
    )


def resolve_annotation_set(
    annotation_set: AnnotationSet,
    records: Sequence[SeqRecord],
    *,
    mode: str,
    record_transforms: Sequence[RecordDisplayTransform] | None = None,
    feature_bindings: FeatureBindings | None = None,
) -> ResolvedAnnotationBundle:
    """Resolve one set against already materialized records."""

    if mode not in {"circular", "linear"}:
        raise ValidationError("Annotation resolver mode must be 'circular' or 'linear'.")
    resolved: list[ResolvedRegionAnnotation] = []
    warnings: list[ResolutionWarning] = []
    transforms = bind_record_transforms(records, record_transforms)
    for annotation in annotation_set.annotations:
        target = annotation.target
        # The drawn units of the target's record: one, or each window of a record with several (D-01).
        units: tuple[tuple[int, IdentityBinding | None], ...]
        if isinstance(target, FeatureIdentitySpan):
            bound = (feature_bindings or {}).get((target.record_key, target.biological_feature_id))
            units = (
                tuple((binding.record_index, binding) for binding in bound)
                if bound is not None
                else tuple((index, None) for index in _bind_record_key(records, target.record_key))
            )
        else:
            units = tuple((index, None) for index in _bind_record(records, target.record))
        pieces = [
            _unit_piece(target, records[record_index], record_index, binding, transforms[record_index])
            for record_index, binding in units
        ]
        drawn = [piece for piece in pieces if piece.missing_count == 0 and piece.segments]
        reference = drawn[0] if drawn else pieces[0]
        record, record_index = records[reference.record_index], reference.record_index
        if isinstance(target, CoordinateSpan):
            # A span is clipped when no unit draws all of it.
            clipped = all(piece.clipped for piece in pieces)
            policy = target.out_of_bounds
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
                        record_id=record.id,
                        record_index=record_index,
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
                        record_id=record.id,
                        record_index=record_index,
                        message=f"Clipped annotation to record {record.id!r}.",
                    )
                )
        elif all(piece.missing_count for piece in pieces):
            missing_count = pieces[0].missing_count
            warnings.append(
                ResolutionWarning(
                    code="feature_selector_unmatched",
                    set_id=annotation_set.id,
                    annotation_id=annotation.id,
                    record_id=record.id,
                    record_index=record_index,
                    missing_count=missing_count,
                    message=f"Skipped annotation: {missing_count} feature selector(s) unmatched.",
                )
            )
            continue
        if not drawn:
            warnings.append(
                ResolutionWarning(
                    code="empty_span",
                    set_id=annotation_set.id,
                    annotation_id=annotation.id,
                    record_id=record.id,
                    record_index=record_index,
                    message=f"Annotation resolved to an empty span on record {record.id!r}.",
                )
            )
            continue
        style = annotation.style or annotation_set.default_style
        if style.line_cap == "arrow" and isinstance(annotation.target, CoordinateSpan):
            raise ValidationError(
                f"Annotation {annotation_set.id}/{annotation.id} uses line_cap=arrow without a directional feature target."
            )
        # One annotation drawn in pieces: the longest piece carries its label.
        labelled = max(drawn, key=lambda piece: piece.span_bp)
        for piece in drawn:
            piece_record = records[piece.record_index]
            segments = piece.segments
            crosses_origin = segments_cross_origin(segments, len(piece_record.seq), is_circular=piece.is_circular)
            resolved.append(
                ResolvedRegionAnnotation(
                    id=annotation.id,
                    set_id=annotation_set.id,
                    record_index=piece.record_index,
                    segments=tuple(segments),
                    midpoint_bp=annotation_midpoint(segments, len(piece_record.seq), crosses_origin=crosses_origin),
                    span_bp=piece.span_bp,
                    label=annotation.label if piece is labelled else "",
                    mark=annotation.mark,
                    lane=annotation.lane,
                    style=style,
                    style_is_annotation_override=annotation.style is not None,
                    legend_label=annotation.legend_label or annotation_set.legend_label,
                    metadata=annotation.metadata,
                    crosses_origin=crosses_origin,
                    source_step=-1 if _read_coord_map(piece_record)[1] < 0 else 1,
                    strand=annotation.strand,
                )
            )
    return ResolvedAnnotationBundle(
        tuple(resolved), tuple(warnings), (annotation_set.id,),
        (annotation_set.id,) if annotation_set.keep_highlight_overlay else (),
    )


def _materialize_sets(options: AnnotationOptions, *, mode: str) -> tuple[AnnotationSet, ...]:
    if options.sets:
        return options.sets
    # `mode` is checked by the resolver; the table readers accept the same two values.
    _mode = cast('Literal["circular", "linear"]', mode)
    if options.table is not None:
        return annotation_sets_from_dataframe(options.table, mode=_mode)
    if options.table_file:
        return read_annotation_table(options.table_file, mode=_mode)
    return ()


def read_annotation_options(
    options: AnnotationOptions | None, *, mode: Literal["circular", "linear"]
) -> AnnotationOptions | None:
    """The options with their table file read into sets, so no later resolution reads it."""

    if options is None or not options.table_file:
        return options
    return AnnotationOptions(sets=_materialize_sets(options, mode=mode))


def resolve_annotations(
    annotations: AnnotationOptions | Sequence[AnnotationSet] | None,
    records: Sequence[SeqRecord],
    *,
    mode: str,
    record_transforms: Sequence[RecordDisplayTransform] | None = None,
    feature_bindings: FeatureBindings | None = None,
) -> ResolvedAnnotationBundle:
    """Resolve all configured sets into one deterministic bundle.

    ``feature_bindings`` are the planner's bindings of FeatureIdentitySpan
    targets; without one, such a target is skipped as unmatched.
    """

    if annotations is None:
        return ResolvedAnnotationBundle(())
    sets = (
        _materialize_sets(annotations, mode=mode)
        if isinstance(annotations, AnnotationOptions)
        else tuple(annotations)
    )
    resolved = [
        resolve_annotation_set(
            item, records, mode=mode, record_transforms=record_transforms,
            feature_bindings=feature_bindings,
        )
        for item in sets
    ]
    return ResolvedAnnotationBundle(
        tuple(annotation for bundle in resolved for annotation in bundle.annotations),
        tuple(warning for bundle in resolved for warning in bundle.warnings),
        tuple(item.id for item in sets),
        tuple(item.id for item in sets if item.keep_highlight_overlay),
    )


__all__ = [
    "annotation_midpoint",
    "merge_annotation_segments",
    "read_annotation_options",
    "resolve_annotation_set",
    "resolve_annotations",
]
