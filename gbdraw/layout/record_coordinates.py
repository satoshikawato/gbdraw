"""Pure source/local/display coordinates, independent of diagram geometry.

Bases are 1-based integers; intervals and samples use 0-based source boundaries.
``length`` always describes the complete source. The existing affine map is
source_base + source_step * local_base_index. Unset local projection returns its
input directly; source projection supplies provenance when fragments are needed.

A wrapped transform draws one record over rows: row r holds the display offsets
[r * wrap_length, (r + 1) * wrap_length). The display offset already holds the
rotation, the strand step, and the source base, so a row is not a source interval.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from math import ceil, isfinite
from typing import Literal, Mapping, Sequence, TypeVar, cast

from gbdraw.exceptions import ValidationError


_Coordinate = TypeVar("_Coordinate", int, float)


MAX_WRAP_LENGTH = 2**53 - 1

_DISPLAY_START_DIAGNOSTIC = {
    "code": "INPUT_INVALID",
    "field": "start",
    "reason": "DISPLAY_START_BOUNDS",
}


def _integer(
    value: int,
    name: str,
    lower: int,
    upper: int,
    diagnostic: Mapping[str, object] | None = None,
) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or not lower <= value <= upper
    ):
        raise ValidationError(
            f"{name} must be an integer in {lower}..{upper}.", diagnostic=diagnostic
        )


def _step(value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value not in (-1, 1):
        raise ValidationError("Orientation step must be +1 or -1.")


def strand_ordered(first: bool, second: bool, strand: int | None) -> tuple[bool, bool]:
    """Swap a pair of end flags between biological (5', 3') and ascending (low, high) order.

    On the minus strand the 5' end is the high one; any other strand reads low to high.
    The map is its own inverse, so it converts in both directions.
    """
    return (second, first) if strand == -1 else (first, second)


@dataclass(frozen=True)
class SourceInterval:
    """Ascending half-open boundaries, with a direction of biological traversal.

    Parts must be supplied in biological order; part_index is the original part
    identity (or the input ordinal when omitted). Terminal flags denote real
    biological ends, not every exon boundary. Partial flags mark a biological end
    that the source location writes as partial (``<``/``>``): the feature
    continues past it, so it is not a terminal either. Empty spans cover no bases.
    """

    start: int
    end: int
    strand: Literal[-1, 1] = 1
    part_index: int | None = None
    biological_start: bool = True
    biological_end: bool = True
    partial_start: bool = False
    partial_end: bool = False

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
        if not isinstance(self.partial_start, bool) or not isinstance(
            self.partial_end, bool
        ):
            raise ValidationError(
                "Partial end flags must be booleans.",
                diagnostic={"code": "INPUT_INVALID"},
            )


@dataclass(frozen=True)
class DisplayFragment:
    """One continuous span; all coordinate pairs are ascending half-open bounds.

    orientation is biological strand * display step. biological_start/end refer
    to the biological traversal; artificial_start/end refer to the geometric
    display_start/end. An artificial seam can never create a biological end.
    open_start/open_end also refer to display_start/end: the drawing continues
    past that edge, at an artificial seam or at a partial end of the source
    location, so a glyph draws no arrowhead, cap, or edge stroke there.
    On a wrapped record display_start/end are local to ``row`` and a row boundary
    is one more artificial, open seam; only the true outer ends stay biological.
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
    open_start: bool
    open_end: bool
    row: int = 0


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
    a duplicated existing sample from an interpolated seam endpoint. On a wrapped
    record positions are local to ``row``; a row boundary opens the segment the same
    way, and seam_sampled still describes the display cut.
    """

    points: tuple[SeriesPoint, ...]
    source_indices: tuple[int | None, ...]
    seam_sampled: bool | None
    row: int = 0


def _line_value(left: tuple[float, float], right: tuple[float, float], position: float) -> float:
    """The value at ``position`` on the line between two (position, value) samples."""
    fraction = (position - left[0]) / (right[0] - left[0])
    return left[1] + fraction * (right[1] - left[1])


@dataclass(frozen=True)
class RecordDisplayTransform:
    length: int
    source_base: int
    source_step: Literal[-1, 1]
    start_coordinate: int | None
    is_circular: bool
    # Set at the draw boundary (``with_wrap``), never by the planner.
    wrap_length: int | None = None
    # The bases on the display axis when a crop shows fewer than ``length``.
    display_length: int | None = None
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
            _integer(
                self.start_coordinate,
                "Display start",
                1,
                self.length,
                _DISPLAY_START_DIAGNOSTIC,
            )
            if not self.is_circular:
                raise ValidationError(
                    "An explicit display start requires a circular record."
                )
            anchor = self.start_coordinate
        object.__setattr__(self, "_anchor_boundary", anchor - (self.source_step == 1))
        if self.wrap_length is not None:
            _integer(self.wrap_length, "Wrap length", 1, MAX_WRAP_LENGTH)
        if self.display_length is not None:
            _integer(self.display_length, "Display length", 1, self.length)

    def with_wrap(self, wrap_length: int, display_length: int) -> RecordDisplayTransform:
        """This record drawn over rows of ``wrap_length`` bp along its ``display_length`` bases."""
        return replace(self, wrap_length=wrap_length, display_length=display_length)

    @property
    def projects(self) -> bool:
        """Whether geometry goes through display fragments, not the unset local path."""
        return self.start_coordinate is not None or self.wrap_length is not None

    @property
    def drawn_length(self) -> int:
        """The bases on the display axis, over every row."""
        return self.length if self.display_length is None else self.display_length

    @property
    def axis_length(self) -> int:
        """The width of one drawn row in bp: the whole display axis when not wrapped."""
        drawn = self.drawn_length
        return drawn if self.wrap_length is None else min(drawn, self.wrap_length)

    @property
    def row_count(self) -> int:
        return ceil(self.drawn_length / self.axis_length)

    def row_length(self, row: int) -> int:
        """The bases on one row: the axis length, or the remainder on the last row."""
        return min(self.axis_length, self.drawn_length - row * self.axis_length)

    def display_index_row(self, index: int) -> tuple[int, int]:
        """The row of a whole-axis display base index, and its index in that row."""
        return (0, index) if self.wrap_length is None else divmod(index, self.wrap_length)

    def source_base_to_display_index(self, base: int) -> int:
        _integer(base, "Source base", 1, self.length)
        if self.start_coordinate is None:
            return self.source_step * (base - self.source_base)
        return self.source_step * (base - self.start_coordinate) % self.length

    def source_boundary_to_display_offset(self, boundary: int) -> int:
        _integer(boundary, "Source boundary", 0, self.length)
        return self._boundary_offset(boundary)

    def source_position_to_display_offset(self, position: float) -> float:
        """Project a finite source-boundary position into display coordinates."""
        if (
            isinstance(position, bool)
            or not isinstance(position, (int, float))
            or not isfinite(float(position))
            or not 0.0 <= float(position) <= float(self.length)
        ):
            raise ValidationError(
                "Source position must be finite and inside [0, source length]."
            )
        return float(self._boundary_offset(float(position)))

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
        # Biological ends lie at the low source boundary on +1 and the high one on -1.
        partial_lo, partial_hi = strand_ordered(
            interval.partial_start, interval.partial_end, interval.strand
        )
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
            open_low = artificial_low or partial_lo
            open_high = artificial_high or partial_hi
            fragments.append(
                DisplayFragment(
                    start,
                    end,
                    local_start,
                    local_end,
                    display_start,
                    display_start + end - start,
                    # Product of two +/-1 steps is +/-1.
                    cast(Literal[-1, 1], interval.strand * self.source_step),
                    interval.biological_start and index == 0,
                    interval.biological_end and index == len(spans) - 1,
                    artificial_low if self.source_step == 1 else artificial_high,
                    artificial_high if self.source_step == 1 else artificial_low,
                    interval.part_index if interval.part_index is not None else 0,
                    open_low if self.source_step == 1 else open_high,
                    open_high if self.source_step == 1 else open_low,
                )
            )
        return self._split_rows(fragments)

    def _split_rows(self, fragments: Sequence[DisplayFragment]) -> tuple[DisplayFragment, ...]:
        """Cut whole-axis fragments at the row boundaries, keeping their traversal order."""
        width = self.wrap_length
        if width is None:
            return tuple(fragments)
        result: list[DisplayFragment] = []
        for fragment in fragments:
            pieces = [
                self._row_piece(fragment, row, width)
                for row in range(fragment.display_start // width, (fragment.display_end - 1) // width + 1)
            ]
            # Biological traversal runs down the rows on a fragment reading against the display.
            result.extend(pieces if fragment.orientation == 1 else reversed(pieces))
        return tuple(result)

    def _row_piece(self, fragment: DisplayFragment, row: int, width: int) -> DisplayFragment:
        """The part of a whole-axis fragment on one row: only its own outer edges keep their flags."""
        origin = row * width
        start = max(fragment.display_start, origin)
        end = min(fragment.display_end, origin + width)
        at_low, at_high = start == fragment.display_start, end == fragment.display_end
        offset, size = start - fragment.display_start, end - start
        source_start = (
            fragment.source_start + offset if self.source_step == 1 else fragment.source_end - offset - size
        )
        at_biological_start, at_biological_end = (
            (at_low, at_high) if fragment.orientation == 1 else (at_high, at_low)
        )
        return DisplayFragment(
            source_start,
            source_start + size,
            fragment.local_start + offset,
            fragment.local_start + offset + size,
            start - origin,
            end - origin,
            fragment.orientation,
            fragment.biological_start and at_biological_start,
            fragment.biological_end and at_biological_end,
            fragment.artificial_start if at_low else True,
            fragment.artificial_end if at_high else True,
            fragment.part_index,
            fragment.open_start if at_low else True,
            fragment.open_end if at_high else True,
            row,
        )

    def project_parts(
        self, parts: Sequence[SourceInterval]
    ) -> tuple[DisplayFragment, ...]:
        """Preserve input part identities/order; only outer parts have terminals or partial ends."""
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
                    partial_start=parts[i].partial_start and i == occupied[0],
                    partial_end=parts[i].partial_end and i == occupied[-1],
                )
            )
        )

    def project_local_parts(
        self,
        parts: Sequence[SourceInterval],
    ) -> Sequence[SourceInterval] | tuple[DisplayFragment, ...]:
        """Pass through existing validated local geometry when nothing projects it."""
        if not self.projects:
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
                    part,
                    start=start,
                    end=end,
                    # Product of two +/-1 steps is +/-1.
                    strand=cast(Literal[-1, 1], part.strand * self.source_step),
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
        previous: float = -1
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
        if not self.projects and self.source_base == 1 and self.source_step == 1:
            return (DisplaySeries(tuple(points), tuple(range(len(points))), None),)
        mapped: list[tuple[float, float, int | None]] = sorted(
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
                seam_value = _line_value((left[0] - self.length, left[1]), right[:2], 0)
                seam_index = None
                mapped.insert(0, (0, seam_value, seam_index))
            mapped.append((self.length, seam_value, seam_index))
        return self._split_series_rows(
            DisplaySeries(
                tuple(SeriesPoint(position, value) for position, value, _ in mapped),
                tuple(index for _, _, index in mapped),
                seam_sampled,
            )
        )

    def _split_series_rows(self, series: DisplaySeries) -> tuple[DisplaySeries, ...]:
        """One segment per row that has samples, opened at each row boundary on one shared edge.

        The edge is the sample on the boundary when there is one, otherwise the line between
        its neighbours, as at the display cut. A last sample on a boundary closes its row.
        """
        width = self.wrap_length
        if width is None or not series.points:
            return (series,)
        rows: list[DisplaySeries] = []
        row = int(series.points[0].position // width)
        points: list[SeriesPoint] = []
        indices: list[int | None] = []

        def close(value: float, index: int | None) -> None:
            nonlocal row, points, indices
            points.append(SeriesPoint(width, value))
            indices.append(index)
            rows.append(DisplaySeries(tuple(points), tuple(indices), series.seam_sampled, row))
            row += 1
            points, indices = [SeriesPoint(0, value)], [index]

        last = series.points[-1].position
        previous: tuple[float, float] | None = None
        for point, index in zip(series.points, series.source_indices, strict=True):
            sample = (point.position, point.value)
            while previous is not None and point.position > (row + 1) * width:
                boundary = (row + 1) * width
                close(_line_value(previous, sample, boundary), None)
            if point.position == (row + 1) * width and point.position < last:
                close(point.value, index)
            else:
                points.append(SeriesPoint(point.position - row * width, point.value))
                indices.append(index)
            previous = sample
        rows.append(DisplaySeries(tuple(points), tuple(indices), series.seam_sampled, row))
        return tuple(rows)

    def project_local_series(
        self,
        points: Sequence[SeriesPoint],
    ) -> Sequence[SeriesPoint] | tuple[DisplaySeries, ...]:
        """Keep unset local samples unchanged; otherwise adapt boundaries once."""
        if not self.projects:
            return points
        if self.start_coordinate is None:
            # A local position already is the display offset of an affine axis.
            order = sorted(range(len(points)), key=lambda i: points[i].position)
            segments = self._split_series_rows(
                DisplaySeries(tuple(points[i] for i in order), tuple(range(len(order))), None)
            )
        else:
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
