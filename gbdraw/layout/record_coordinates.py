"""Pure source/local/display coordinates, independent of diagram geometry.

Bases are 1-based integers; intervals and samples use 0-based source boundaries.
``length`` always describes the complete source. The existing affine map is
source_base + source_step * local_base_index. Unset local projection returns its
input directly; source projection supplies provenance when fragments are needed.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from math import isfinite
from typing import Literal, Sequence, TypeVar

from gbdraw.exceptions import ValidationError


_Coordinate = TypeVar("_Coordinate", int, float)


def _integer(value: int, name: str, lower: int, upper: int) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or not lower <= value <= upper
    ):
        raise ValidationError(f"{name} must be an integer in {lower}..{upper}.")


def _step(value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value not in (-1, 1):
        raise ValidationError("Orientation step must be +1 or -1.")


@dataclass(frozen=True)
class SourceInterval:
    """Ascending half-open boundaries, with a direction of biological traversal.

    Parts must be supplied in biological order; part_index is the original part
    identity (or the input ordinal when omitted). Terminal flags denote real
    biological ends, not every exon boundary. Empty spans cover no bases.
    """

    start: int
    end: int
    strand: Literal[-1, 1] = 1
    part_index: int | None = None
    biological_start: bool = True
    biological_end: bool = True

    def __post_init__(self) -> None:
        if isinstance(self.end, bool) or not isinstance(self.end, int) or self.end < 0:
            raise ValidationError("Interval end must be a non-negative integer.")
        _integer(self.start, "Interval start", 0, self.end)
        _step(self.strand)
        if self.part_index is not None:
            if (
                isinstance(self.part_index, bool)
                or not isinstance(self.part_index, int)
                or self.part_index < 0
            ):
                raise ValidationError(
                    "Part index must be a non-negative integer or None."
                )
        if not isinstance(self.biological_start, bool) or not isinstance(
            self.biological_end, bool
        ):
            raise ValidationError("Biological terminal flags must be booleans.")


@dataclass(frozen=True)
class DisplayFragment:
    """One continuous span; all coordinate pairs are ascending half-open bounds.

    orientation is biological strand * display step. biological_start/end refer
    to the biological traversal; artificial_start/end refer to the geometric
    display_start/end. An artificial seam can never create a biological end.
    """

    source_start: int
    source_end: int
    local_start: int
    local_end: int
    display_start: int
    display_end: int
    orientation: Literal[-1, 1]
    biological_start: bool
    biological_end: bool
    artificial_start: bool
    artificial_end: bool
    part_index: int


@dataclass(frozen=True)
class SeriesPoint:
    """A boundary-domain sample position and value (never a 1-based base)."""

    position: float
    value: float

    def __post_init__(self) -> None:
        for value in (self.position, self.value):
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not isfinite(value)
            ):
                raise ValidationError(
                    "Series positions and values must be finite numbers."
                )


@dataclass(frozen=True)
class DisplaySeries:
    """A continuous segment in display order, opened at the display seam.

    source_indices identify original samples; None denotes interpolation.
    seam_sampled is None on the unset affine path, otherwise it distinguishes
    a duplicated existing sample from an interpolated seam endpoint.
    """

    points: tuple[SeriesPoint, ...]
    source_indices: tuple[int | None, ...]
    seam_sampled: bool | None


@dataclass(frozen=True)
class RecordDisplayTransform:
    length: int
    source_base: int
    source_step: Literal[-1, 1]
    start_coordinate: int | None
    is_circular: bool
    _anchor_boundary: int = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if (
            isinstance(self.length, bool)
            or not isinstance(self.length, int)
            or self.length <= 0
        ):
            raise ValidationError("Complete source length must be a positive integer.")
        _step(self.source_step)
        _integer(self.source_base, "Source base", 1, self.length)
        if not isinstance(self.is_circular, bool):
            raise ValidationError("Effective circular topology must be a boolean.")
        anchor = self.source_base
        if self.start_coordinate is not None:
            _integer(self.start_coordinate, "Display start", 1, self.length)
            if not self.is_circular:
                raise ValidationError(
                    "An explicit display start requires a circular record."
                )
            anchor = self.start_coordinate
        object.__setattr__(self, "_anchor_boundary", anchor - (self.source_step == 1))

    def source_base_to_display_index(self, base: int) -> int:
        _integer(base, "Source base", 1, self.length)
        if self.start_coordinate is None:
            return self.source_step * (base - self.source_base)
        return self.source_step * (base - self.start_coordinate) % self.length

    def source_boundary_to_display_offset(self, boundary: int) -> int:
        _integer(boundary, "Source boundary", 0, self.length)
        return self._boundary_offset(boundary)

    def _boundary_offset(self, boundary: _Coordinate) -> _Coordinate:
        offset = self.source_step * (boundary - self._anchor_boundary)
        return offset if self.start_coordinate is None else offset % self.length

    def display_index_to_source_base(self, index: int) -> int:
        _integer(index, "Display base index", 0, self.length - 1)
        if self.start_coordinate is None:
            base = self.source_base + self.source_step * index
            _integer(base, "Source base", 1, self.length)
            return base
        return (self.start_coordinate - 1 + self.source_step * index) % self.length + 1

    def _local_boundary_to_source(self, boundary: _Coordinate) -> _Coordinate:
        return self.source_base - (self.source_step == 1) + self.source_step * boundary

    def local_boundary_to_display_offset(self, boundary: int) -> int:
        _integer(boundary, "Local boundary", 0, self.length)
        if self.start_coordinate is None:
            return boundary
        return self.source_boundary_to_display_offset(
            self._local_boundary_to_source(boundary)
        )

    def _validate_interval(self, interval: SourceInterval) -> None:
        if not isinstance(interval, SourceInterval):
            raise ValidationError("Expected a SourceInterval.")
        _integer(interval.end, "Source interval end", 0, self.length)

    def project_interval(self, interval: SourceInterval) -> tuple[DisplayFragment, ...]:
        """Project source coverage, preserving biological traversal order."""
        self._validate_interval(interval)
        lo, hi = interval.start, interval.end
        if lo == hi:
            return ()
        cuts = [lo, hi]
        if self.start_coordinate is not None and lo < self._anchor_boundary < hi:
            cuts.insert(1, self._anchor_boundary)
        spans = list(zip(cuts, cuts[1:]))
        if interval.strand == -1:
            spans.reverse()
        local_anchor = self.source_base - (self.source_step == 1)
        fragments = []
        for index, (start, end) in enumerate(spans):
            leading = start if self.source_step == 1 else end
            display_start = self.source_boundary_to_display_offset(leading)
            local_start, local_end = sorted(
                (
                    self.source_step * (start - local_anchor),
                    self.source_step * (end - local_anchor),
                )
            )
            artificial_low, artificial_high = start != lo, end != hi
            fragments.append(
                DisplayFragment(
                    start,
                    end,
                    local_start,
                    local_end,
                    display_start,
                    display_start + end - start,
                    interval.strand * self.source_step,
                    interval.biological_start and index == 0,
                    interval.biological_end and index == len(spans) - 1,
                    artificial_low if self.source_step == 1 else artificial_high,
                    artificial_high if self.source_step == 1 else artificial_low,
                    interval.part_index if interval.part_index is not None else 0,
                )
            )
        return tuple(fragments)

    def project_parts(
        self, parts: Sequence[SourceInterval]
    ) -> tuple[DisplayFragment, ...]:
        """Preserve input part identities/order; only outer parts have terminals."""
        for part in parts:
            self._validate_interval(part)
        occupied = [i for i, part in enumerate(parts) if part.start != part.end]
        if not occupied:
            return ()
        return tuple(
            fragment
            for i in occupied
            for fragment in self.project_interval(
                replace(
                    parts[i],
                    part_index=parts[i].part_index
                    if parts[i].part_index is not None
                    else i,
                    biological_start=parts[i].biological_start and i == occupied[0],
                    biological_end=parts[i].biological_end and i == occupied[-1],
                )
            )
        )

    def project_local_parts(
        self,
        parts: Sequence[SourceInterval],
    ) -> Sequence[SourceInterval] | tuple[DisplayFragment, ...]:
        """Pass through existing validated local geometry when start is unset."""
        if self.start_coordinate is None:
            return parts
        source_parts = []
        for part in parts:
            start, end = sorted(
                (
                    self._local_boundary_to_source(part.start),
                    self._local_boundary_to_source(part.end),
                )
            )
            source_parts.append(
                replace(
                    part, start=start, end=end, strand=part.strand * self.source_step
                )
            )
        return self.project_parts(source_parts)

    def project_series(
        self, points: Sequence[SeriesPoint]
    ) -> tuple[DisplaySeries, ...]:
        """Project one periodic track sampled at increasing boundaries in [0,L).

        Explicit projection opens the periodic polyline at the display cut. Both
        edge values match: use the cut sample if present, otherwise interpolate
        linearly between its neighbors, including the source-origin gap. No
        path joins display L back to 0. Unset projection neither closes the track
        nor adds samples. This primitive accepts one continuous periodic track;
        it does not infer missing-data gaps.
        """
        previous = -1
        for point in points:
            if (
                not isinstance(point, SeriesPoint)
                or not 0 <= point.position < self.length
            ):
                raise ValidationError(
                    "Source series samples must be in [0, source length)."
                )
            if point.position <= previous:
                raise ValidationError(
                    "Source series samples must be strictly increasing."
                )
            previous = point.position
        if not points:
            return ()
        if (
            self.start_coordinate is None
            and self.source_base == 1
            and self.source_step == 1
        ):
            return (DisplaySeries(tuple(points), tuple(range(len(points))), None),)
        mapped = sorted(
            (self._boundary_offset(point.position), point.value, i)
            for i, point in enumerate(points)
        )
        seam_sampled = None
        if self.start_coordinate is not None:
            seam_sampled = mapped[0][0] == 0
            if seam_sampled:
                seam_value, seam_index = mapped[0][1:]
            else:
                left, right = mapped[-1], mapped[0]
                left_position = left[0] - self.length
                fraction = -left_position / (right[0] - left_position)
                seam_value = left[1] + fraction * (right[1] - left[1])
                seam_index = None
                mapped.insert(0, (0, seam_value, seam_index))
            mapped.append((self.length, seam_value, seam_index))
        return (
            DisplaySeries(
                tuple(SeriesPoint(position, value) for position, value, _ in mapped),
                tuple(index for _, _, index in mapped),
                seam_sampled,
            ),
        )

    def project_local_series(
        self,
        points: Sequence[SeriesPoint],
    ) -> Sequence[SeriesPoint] | tuple[DisplaySeries, ...]:
        """Keep unset local samples unchanged; otherwise adapt boundaries once."""
        if self.start_coordinate is None:
            return points
        source_points = []
        for point in points:
            if not 0 <= point.position < self.length:
                raise ValidationError(
                    "Local series samples must be in [0, source length)."
                )
            position = self._local_boundary_to_source(point.position) % self.length
            source_points.append(SeriesPoint(position, point.value))
        order = sorted(
            range(len(source_points)), key=lambda i: source_points[i].position
        )
        segments = self.project_series(tuple(source_points[i] for i in order))
        return tuple(
            replace(
                segment,
                source_indices=tuple(
                    order[i] if i is not None else None for i in segment.source_indices
                ),
            )
            for segment in segments
        )

    def alignment_cut_breakpoints(self, span: SourceInterval) -> tuple[float, ...]:
        """Interior cuts in directed endpoint-linear t, never exact gapped bases."""
        self._validate_interval(span)
        if (
            self.start_coordinate is None
            or not span.start < self._anchor_boundary < span.end
        ):
            return ()
        distance = (
            self._anchor_boundary - span.start
            if span.strand == 1
            else span.end - self._anchor_boundary
        )
        return (distance / (span.end - span.start),)


def alignment_cut_breakpoints(
    *spans: tuple[RecordDisplayTransform, SourceInterval],
) -> tuple[float, ...]:
    """Union query/subject interior cuts on the same endpoint-linear parameter."""
    return tuple(
        sorted(
            {
                cut
                for transform, span in spans
                for cut in transform.alignment_cut_breakpoints(span)
            }
        )
    )
