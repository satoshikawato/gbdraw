#!/usr/bin/env python
# coding: utf-8

"""Circular feature paths: introns, arrows and blocks between explicit radii.

Every arc goes through :func:`gbdraw.svg.circular_arcs.circle_arc_to`, which
splits it so that no drawn arc is near 180 or 360 degrees. The factor-based
functions resolve the three radii of a feature track and draw through the
explicit-radii functions.
"""

import math
from typing import TYPE_CHECKING, Dict, Union

from .arrows import (
    LinearBlockGlyph,
    calculate_arrow_shaft_bounds,
    block_path_with_open_caps,
    cap_arrow_head_length,
    has_arrow_shaft,
    is_legacy_arrow_short,
    linear_block_glyph,
    set_arrow_shoulder,
)
from ..layout.circular import calculate_feature_position_factors_circular
from .circular_arcs import circle_arc_to, circle_point, position_angle, span_angle
from .path_data import path_point

if TYPE_CHECKING:
    from ..features.objects import FeatureObject

CircularCoordDict = Dict[str, Union[str, int]]


def _point(radius: float, position: float, total_length: int) -> str:
    return path_point(*circle_point(radius, position_angle(position, total_length)))


def _arc_to(radius: float, from_position: float, span_bp: float, total_length: int) -> str:
    return circle_arc_to(radius, position_angle(from_position, total_length), span_angle(span_bp, total_length))


def _coord_len_bp(start: int, end: int, total_length: int, *, strand: str = "positive") -> int:
    if strand == "negative":
        return (int(start) - int(end)) % int(total_length)
    return (int(end) - int(start)) % int(total_length)


def _strand_ends(coord_dict: Dict[str, Union[str, int]], total_length: int) -> tuple[str, int, int, float]:
    """The strand, the tail and head positions of a directional block, and the direction of travel."""
    coord_strand = str(coord_dict["coord_strand"])
    coord_start = int(coord_dict["coord_start"]) % total_length
    coord_end = int(coord_dict["coord_end"]) % total_length
    if coord_strand == "positive":
        return coord_strand, coord_start, coord_end, 1.0
    return coord_strand, coord_end, coord_start, -1.0


def _arrow_span(coord_dict: Dict[str, Union[str, int]], total_length: int) -> tuple[str, int, int, float, float]:
    """``_strand_ends`` and the block's length in bp, from tail to tip (a display fragment keeps its own span)."""
    coord_strand, arrow_start, arrow_end, direction = _strand_ends(coord_dict, total_length)
    coord_len_bp: float = _coord_len_bp(arrow_start, arrow_end, total_length, strand=coord_strand)
    if coord_dict.get("display_fragment"):
        coord_len_bp = int(coord_dict["coord_end"]) - int(coord_dict["coord_start"])
    return coord_strand, arrow_start, arrow_end, direction, float(coord_len_bp)


def _open_tail(coord_dict: Dict[str, Union[str, int]], coord_strand: str) -> bool:
    return bool(coord_dict.get("open_start" if coord_strand == "positive" else "open_end"))


def _factor_radii(
    radius: float,
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    track_ratio: float,
    cds_ratio: float,
    offset: float,
    track_type: str,
    strandedness: bool,
    track_id: int,
) -> tuple[float, float, float]:
    factors = calculate_feature_position_factors_circular(
        total_length, str(coord_dict["coord_strand"]), track_ratio, cds_ratio, offset, track_type, strandedness, track_id
    )
    return radius * factors[0], radius * factors[1], radius * factors[2]


def generate_circular_intron_path_with_radii(
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    center_radius_px: float,
) -> list[str]:
    raw_start = int(coord_dict["coord_start"])
    raw_end = int(coord_dict["coord_end"])
    arc_start = raw_start % total_length
    coord_len_bp = (raw_end % total_length - arc_start) % total_length
    # BioPython encodes an effective 1 bp intron as start = end + 1; it is not a near-full circle.
    if raw_start == raw_end + 1:
        coord_len_bp = 1
    # Introns follow the axis direction (clockwise in this coordinate system).
    feature_path = (
        f"M {_point(center_radius_px, arc_start, total_length)} "
        f"{_arc_to(center_radius_px, arc_start, coord_len_bp, total_length)}"
    )
    return ["line", feature_path]


def generate_circular_arrow_path_with_radii(
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    cds_arrow_length: float,
    inner_radius_px: float,
    center_radius_px: float,
    outer_radius_px: float,
) -> list[str]:
    """Generate the legacy five-vertex circular arrow with explicit radii."""
    coord_strand, arrow_start, arrow_end, direction, coord_len_bp = _arrow_span(coord_dict, total_length)

    point = _point(center_radius_px, arrow_end, total_length)
    start_inner = _point(inner_radius_px, arrow_start, total_length)
    start_outer = _point(outer_radius_px, arrow_start, total_length)
    open_tail = _open_tail(coord_dict, coord_strand)
    if is_legacy_arrow_short(coord_len_bp, cds_arrow_length):
        return block_path_with_open_caps(f"M {start_inner} L {point} L {start_outer} z", open_tail=open_tail)

    shoulder = set_arrow_shoulder(coord_strand, arrow_end, cds_arrow_length) % total_length
    arc_len_bp = direction * (coord_len_bp - cds_arrow_length)
    feature_path = (
        f"M {start_inner} {_arc_to(inner_radius_px, arrow_start, arc_len_bp, total_length)} L {point} "
        f"L {_point(outer_radius_px, shoulder, total_length)} "
        f"{_arc_to(outer_radius_px, shoulder, -arc_len_bp, total_length)} z"
    )
    return block_path_with_open_caps(feature_path, open_tail=open_tail)


def generate_circular_narrow_arrow_path_with_radii(
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    head_length_bp: float,
    inner_radius_px: float,
    center_radius_px: float,
    outer_radius_px: float,
    shaft_width_ratio: float,
) -> list[str]:
    """Generate a narrow-shaft arrow using explicit feature radii."""
    coord_strand, arrow_start, arrow_end, direction, coord_len_bp = _arrow_span(coord_dict, total_length)
    resolved_head_length_bp = cap_arrow_head_length(coord_len_bp, head_length_bp)

    point = _point(center_radius_px, arrow_end, total_length)
    open_tail = _open_tail(coord_dict, coord_strand)
    if not has_arrow_shaft(coord_len_bp, resolved_head_length_bp):
        feature_path = (
            f"M {_point(inner_radius_px, arrow_start, total_length)} L {point} "
            f"L {_point(outer_radius_px, arrow_start, total_length)} z"
        )
        return block_path_with_open_caps(feature_path, open_tail=open_tail)

    shoulder = set_arrow_shoulder(coord_strand, arrow_end, resolved_head_length_bp) % total_length
    shaft_inner_radius_px, shaft_outer_radius_px = calculate_arrow_shaft_bounds(
        inner_radius_px,
        center_radius_px,
        outer_radius_px,
        shaft_width_ratio,
    )
    shaft_len_bp = direction * (coord_len_bp - resolved_head_length_bp)
    feature_path = (
        f"M {_point(shaft_inner_radius_px, arrow_start, total_length)} "
        f"{_arc_to(shaft_inner_radius_px, arrow_start, shaft_len_bp, total_length)} "
        f"L {_point(inner_radius_px, shoulder, total_length)} L {point} "
        f"L {_point(outer_radius_px, shoulder, total_length)} L {_point(shaft_outer_radius_px, shoulder, total_length)} "
        f"{_arc_to(shaft_outer_radius_px, shoulder, -shaft_len_bp, total_length)} z"
    )
    return block_path_with_open_caps(feature_path, open_tail=open_tail)


def generate_circular_rectangle_path_with_radii(
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    inner_radius_px: float,
    center_radius_px: float,
    outer_radius_px: float,
) -> list[str]:
    del center_radius_px
    coord_strand = str(coord_dict["coord_strand"])
    coord_start = int(coord_dict["coord_start"])
    coord_end = int(coord_dict["coord_end"])
    coord_len_bp = coord_end - coord_start if coord_end >= coord_start else (coord_end - coord_start) + total_length
    if coord_strand == "positive":
        rect_start, rect_end, span_bp = coord_start, coord_end, float(coord_len_bp)
    else:
        rect_start, rect_end, span_bp = coord_end, coord_start, -float(coord_len_bp)

    start_inner = _point(inner_radius_px, rect_start, total_length)
    end_outer = f"L {_point(outer_radius_px, rect_end, total_length)}"
    feature_path = (
        f"M {start_inner} {_arc_to(inner_radius_px, rect_start, span_bp, total_length)} {end_outer} "
        f"{_arc_to(outer_radius_px, rect_end, -span_bp, total_length)} z"
    )
    return block_path_with_open_caps(
        feature_path,
        open_tail=_open_tail(coord_dict, coord_strand),
        tail_cap=f"L {start_inner}",
        end_cap=end_outer if coord_dict.get("open_end" if coord_strand == "positive" else "open_start") else None,
    )


def generate_circular_narrow_arrow_shaft_path_with_radii(
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    inner_radius_px: float,
    center_radius_px: float,
    outer_radius_px: float,
    shaft_width_ratio: float,
) -> list[str]:
    """Generate a reduced-width body block for a multipart arrow."""
    shaft_inner_radius_px, shaft_outer_radius_px = calculate_arrow_shaft_bounds(
        inner_radius_px,
        center_radius_px,
        outer_radius_px,
        shaft_width_ratio,
    )
    return generate_circular_rectangle_path_with_radii(
        coord_dict,
        total_length,
        shaft_inner_radius_px,
        center_radius_px,
        shaft_outer_radius_px,
    )


def circular_feature_parts(feature_object: "FeatureObject", total_length: int) -> list[tuple[CircularCoordDict, bool]]:
    """The parts a Circular feature is drawn with, each with its ``is_last``: the record's display fragments when
    it is projected; one block across the origin for a two-block feature that wraps it, so it is drawn without a
    split; else its location parts. The feature path and the leader anchor (``circular_feature_edge_radii_at``)
    read these same parts."""
    if feature_object.display_parts is None:
        merged = _merged_origin_spanning_block(feature_object, total_length)
        if merged is not None:
            return [merged]
    coords = (feature_object.location if feature_object.display_parts is None
              else [part.location for part in feature_object.display_parts])
    parts: list[tuple[CircularCoordDict, bool]] = []
    for coord in coords:
        coord_dict: CircularCoordDict = {
            "coord_type": str(coord.kind),
            "coord_strand": coord.strand,
            "coord_start": coord.start,
            "coord_end": coord.end,
            "open_start": coord.open_start,
            "open_end": coord.open_end,
        }
        if feature_object.display_parts is not None:
            coord_dict["display_fragment"] = True
        parts.append((coord_dict, bool(coord.is_last)))
    return parts


def _merged_origin_spanning_block(
    feature_object: "FeatureObject", total_length: int
) -> tuple[CircularCoordDict, bool] | None:
    """A two-block feature across the origin (e.g. a D-loop as ``join(end..genome,1..start)``) as one block."""
    block_coords = [coord for coord in feature_object.location if coord.kind == "block"]
    if len(block_coords) != 2:
        return None
    starts = [int(coord.start) for coord in block_coords]
    ends = [int(coord.end) for coord in block_coords]
    if min(starts) > 1 or max(ends) < total_length:
        return None
    left_blocks = [coord for coord in block_coords if int(coord.start) <= 1]
    right_blocks = [coord for coord in block_coords if int(coord.end) >= total_length]
    if len(left_blocks) != 1 or len(right_blocks) != 1 or left_blocks[0] is right_blocks[0]:
        return None
    merged_start = max(starts)
    merged_end = min(ends)
    if merged_start <= merged_end:
        return None
    # The block that runs to the record end starts the merged block; the one from 1 ends it.
    merged: CircularCoordDict = {
        "coord_type": "block",
        "coord_strand": block_coords[0].strand,
        "coord_start": merged_start,
        "coord_end": merged_end,
        "open_start": right_blocks[0].open_start,
        "open_end": left_blocks[0].open_end,
    }
    return merged, any(coord.is_last for coord in block_coords)


def _ray_radius_on_chord(angle_deg: float, start: tuple[float, float], end: tuple[float, float]) -> float:
    """Where the ray from the centre at ``angle_deg`` meets the straight edge from ``start`` to ``end``."""
    ux, uy = circle_point(1.0, angle_deg)
    dx, dy = end[0] - start[0], end[1] - start[1]
    across = ux * dy - uy * dx
    if abs(across) <= 1e-12:
        return math.hypot(*start)
    return (start[0] * dy - start[1] * dx) / across


def circular_block_edge_radii_at(
    coord_dict: CircularCoordDict,
    glyph: LinearBlockGlyph,
    position_bp: float,
    total_length: int,
    head_length_bp: float,
    inner_radius_px: float,
    center_radius_px: float,
    outer_radius_px: float,
    shaft_width_ratio: float,
) -> tuple[float, float]:
    """The inner and outer radius, at ``position_bp`` (held to the block), of the outline the block's path draws as
    ``glyph`` (``generate_circular_*_path_with_radii``): the lane's edges on a rectangle, the shaft's on a narrowed
    shaft, and in an arrowhead the straight edges from the shoulder's corners to the tip on the centre radius."""
    if glyph == "rectangle":
        return inner_radius_px, outer_radius_px
    full_width = float(shaft_width_ratio) == 1.0
    shaft_inner_px, shaft_outer_px = (
        (inner_radius_px, outer_radius_px)
        if full_width
        else calculate_arrow_shaft_bounds(inner_radius_px, center_radius_px, outer_radius_px, shaft_width_ratio)
    )
    if glyph == "shaft":
        return shaft_inner_px, shaft_outer_px
    coord_strand, tail, tip, direction, coord_len_bp = _arrow_span(coord_dict, total_length)
    if full_width:
        head_bp, triangle = float(head_length_bp), is_legacy_arrow_short(coord_len_bp, head_length_bp)
    else:
        head_bp = cap_arrow_head_length(coord_len_bp, head_length_bp)
        triangle = not has_arrow_shaft(coord_len_bp, head_bp)
    shoulder_from_tail_bp = 0.0 if triangle else coord_len_bp - head_bp
    from_tail_bp = (direction * (float(position_bp) - tail)) % total_length
    if from_tail_bp > coord_len_bp:
        # Past the block: hold to the nearer end.
        from_tail_bp = coord_len_bp if from_tail_bp - coord_len_bp < total_length - from_tail_bp else 0.0
    if from_tail_bp < shoulder_from_tail_bp:
        return shaft_inner_px, shaft_outer_px
    shoulder_angle = position_angle(tail + direction * shoulder_from_tail_bp, total_length)
    tip_point = circle_point(center_radius_px, position_angle(tip, total_length))
    angle = position_angle(tail + direction * from_tail_bp, total_length)
    return (
        _ray_radius_on_chord(angle, circle_point(inner_radius_px, shoulder_angle), tip_point),
        _ray_radius_on_chord(angle, circle_point(outer_radius_px, shoulder_angle), tip_point),
    )


def circular_feature_edge_radii_at(
    feature_object: "FeatureObject",
    position_bp: float,
    total_length: int,
    head_length_bp: float,
    radii_px: tuple[float, float, float],
    shaft_width_ratio: float,
) -> tuple[float, float]:
    """The inner and outer radius of the outline ``feature_object`` is drawn with at ``position_bp``, in its lane of
    ``radii_px`` (inner, centre, outer): read on the drawn block that holds the position, or the nearest one. A
    leader that starts there starts on the feature's ink, also beside an arrowhead that narrows toward its tip."""
    best: tuple[float, CircularCoordDict, bool] | None = None
    for coord_dict, is_last in circular_feature_parts(feature_object, total_length):
        if coord_dict["coord_type"] != "block":
            continue
        start, end = int(coord_dict["coord_start"]), int(coord_dict["coord_end"])
        span_bp = (end - start) if coord_dict.get("display_fragment") else (end - start) % total_length
        offset_bp = (float(position_bp) - start) % total_length
        miss_bp = 0.0 if offset_bp <= span_bp else min(offset_bp - span_bp, total_length - offset_bp)
        if best is None or miss_bp < best[0]:
            best = (miss_bp, coord_dict, is_last)
    if best is None:
        return radii_px[0], radii_px[2]
    _, coord_dict, is_last = best
    glyph = linear_block_glyph(feature_object.glyph_kind, str(coord_dict["coord_strand"]), is_last, shaft_width_ratio)
    return circular_block_edge_radii_at(
        coord_dict, glyph, position_bp, total_length, head_length_bp, *radii_px, shaft_width_ratio
    )


def generate_circular_intron_path(
    radius: float,
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    track_ratio: float,
    cds_ratio: float,
    offset: float,
    track_type: str,
    strandedness: bool,
    track_id: int = 0,
) -> list[str]:
    """An intron on the centre radius of its feature track (factor-based layout); ``["line", path]``."""
    _, center_radius_px, _ = _factor_radii(
        radius, coord_dict, total_length, track_ratio, cds_ratio, offset, track_type, strandedness, track_id
    )
    return generate_circular_intron_path_with_radii(coord_dict, total_length, center_radius_px)


def generate_circular_arrow_path(
    radius: float,
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    cds_arrow_length: float,
    track_ratio: float,
    cds_ratio: float,
    offset: float,
    track_type: str,
    strandedness: bool,
    track_id: int = 0,
) -> list[str]:
    """The legacy five-vertex arrow of a feature track (factor-based layout); ``["block", path, ...]``."""
    return generate_circular_arrow_path_with_radii(
        coord_dict,
        total_length,
        cds_arrow_length,
        *_factor_radii(radius, coord_dict, total_length, track_ratio, cds_ratio, offset, track_type, strandedness, track_id),
    )


def generate_circular_narrow_arrow_path(
    radius: float,
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    head_length_bp: float,
    track_ratio: float,
    cds_ratio: float,
    offset: float,
    track_type: str,
    strandedness: bool,
    shaft_width_ratio: float,
    track_id: int = 0,
) -> list[str]:
    """Generate a narrow-shaft arrow using legacy factor-based layout."""
    return generate_circular_narrow_arrow_path_with_radii(
        coord_dict,
        total_length,
        head_length_bp,
        *_factor_radii(radius, coord_dict, total_length, track_ratio, cds_ratio, offset, track_type, strandedness, track_id),
        shaft_width_ratio,
    )


def generate_circular_narrow_arrow_shaft_path(
    radius: float,
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    track_ratio: float,
    cds_ratio: float,
    offset: float,
    track_type: str,
    strandedness: bool,
    shaft_width_ratio: float,
    track_id: int = 0,
) -> list[str]:
    """Generate a multipart narrow arrow body using legacy factor-based layout."""
    return generate_circular_narrow_arrow_shaft_path_with_radii(
        coord_dict,
        total_length,
        *_factor_radii(radius, coord_dict, total_length, track_ratio, cds_ratio, offset, track_type, strandedness, track_id),
        shaft_width_ratio,
    )


def generate_circular_rectangle_path(
    radius: float,
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    track_ratio: float,
    cds_ratio,
    offset,
    track_type: str,
    strandedness: bool,
    track_id: int = 0,
) -> list[str]:
    """A rectangular block of a feature track (factor-based layout); ``["block", path, ...]``."""
    return generate_circular_rectangle_path_with_radii(
        coord_dict,
        total_length,
        *_factor_radii(radius, coord_dict, total_length, track_ratio, cds_ratio, offset, track_type, strandedness, track_id),
    )


def generate_circular_arrowhead_path_with_radii(
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    cds_arrow_length: float,
    inner_radius_px: float,
    center_radius_px: float,
    outer_radius_px: float,
) -> list[str]:
    """Compatibility wrapper for the legacy five-vertex arrow helper."""
    return generate_circular_arrow_path_with_radii(
        coord_dict,
        total_length,
        cds_arrow_length,
        inner_radius_px,
        center_radius_px,
        outer_radius_px,
    )


def generate_circular_arrowhead_path(
    radius: float,
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    cds_arrow_length: float,
    track_ratio: float,
    cds_ratio: float,
    offset: float,
    track_type: str,
    strandedness: bool,
    track_id: int = 0,
) -> list[str]:
    """Compatibility wrapper for the legacy five-vertex arrow helper."""
    return generate_circular_arrow_path(
        radius,
        coord_dict,
        total_length,
        cds_arrow_length,
        track_ratio,
        cds_ratio,
        offset,
        track_type,
        strandedness,
        track_id,
    )


__all__ = [
    "CircularCoordDict",
    "circular_block_edge_radii_at",
    "circular_feature_edge_radii_at",
    "circular_feature_parts",
    "generate_circular_arrow_path",
    "generate_circular_arrow_path_with_radii",
    "generate_circular_arrowhead_path",
    "generate_circular_arrowhead_path_with_radii",
    "generate_circular_intron_path",
    "generate_circular_intron_path_with_radii",
    "generate_circular_narrow_arrow_path",
    "generate_circular_narrow_arrow_path_with_radii",
    "generate_circular_narrow_arrow_shaft_path",
    "generate_circular_narrow_arrow_shaft_path_with_radii",
    "generate_circular_rectangle_path",
    "generate_circular_rectangle_path_with_radii",
]
