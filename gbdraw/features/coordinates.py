#!/usr/bin/env python
# coding: utf-8

import logging
from typing import List, NamedTuple, Protocol, Sequence, Tuple, cast

from Bio.SeqFeature import AfterPosition, BeforePosition, SimpleLocation

from .objects import FeatureDisplayPart, FeatureLocation, FeatureLocationPart, Strand
from ..layout.record_coordinates import DisplayFragment, RecordDisplayTransform, SourceInterval, strand_ordered


logger = logging.getLogger(__name__)


class PartialEnds(NamedTuple):
    """Whether a location part's biological start (5') and end (3') are partial."""

    start: bool
    end: bool


NOT_PARTIAL = PartialEnds(False, False)


def location_partial_ends(exons: Sequence[SimpleLocation]) -> tuple[PartialEnds, ...]:
    """The partial ends of each part, in biological order, for feature drawing.

    The source writes ``<`` on a low boundary and ``>`` on a high one, whether the
    file or a gbdraw crop wrote it. Only the outer ends of a location can be
    partial: the first part's start and the last part's end. A part without a
    strand is read low to high. CDS translation reads the 5' marker on its own
    (``gbdraw.core.sequence._cds_five_prime_is_complete``, which also reads the
    GFF3 ``start_range``/``end_range`` attributes).
    """
    last = len(exons) - 1
    ends = []
    for index, exon in enumerate(exons):
        low = isinstance(exon.start, BeforePosition)
        high = isinstance(exon.end, AfterPosition)
        start, end = strand_ordered(low, high, exon.strand)
        ends.append(PartialEnds(start and index == 0, end and index == last))
    return tuple(ends)


class _Span(Protocol):
    """Half-open, 0-based bounds of a feature part."""

    @property
    def start(self) -> int: ...

    @property
    def end(self) -> int: ...


def location_parts_text(parts: Sequence[tuple[int, int, int | None]]) -> str:
    """A location's parts (0-based, half-open, with Biopython strands) as a GenBank file writes them: 1-based, in
    ``join`` when there are several, in ``complement`` when the first part is on the minus strand. Messages that name
    a feature write its location this way."""

    text = ",".join(f"{start + 1}..{end}" for start, end, _strand in parts)
    text = text if len(parts) == 1 else f"join({text})"
    return f"complement({text})" if parts and parts[0][2] == -1 else text


def location_fits_record(parts: Sequence[_Span], length: int) -> bool:
    """Whether a location fits in a record of ``length`` bases: each part starts in the record and covers at most
    all of it. A part that starts in the record may end past its end, the GFF3 way to write a span over a circular
    record's origin; a part that starts outside the record, or is longer than it, cannot be drawn on it (OV-3678)."""

    return all(0 <= int(part.start) < length and int(part.end) - int(part.start) <= length for part in parts)


class OriginCrossing(NamedTuple):
    """A location that crosses its record's origin: the record's length and the gap between its parts that the
    location does not cover (``uncovered_gap``)."""

    length: int
    gap_start: int
    gap_end: int


def uncovered_gap(parts: Sequence[_Span], length: int) -> tuple[int, int] | None:
    """The gap a set of parts leaves out when the shortest arc that covers them crosses the record's origin: the
    largest empty gap between the parts, when it is larger than the gap over the origin; else None.

    The arc follows the positions of the parts, not their order, so a location whose parts come in either order
    (a GFF3 reader sorts its rows; a trans-spliced feature's parts come from separate transcripts) is read alike.
    """

    spans = sorted((int(p.start), int(p.end)) for p in parts if int(p.start) < int(p.end))
    if len(spans) < 2:
        return None
    best: tuple[int, int] | None = None
    reach = spans[0][1]
    for start, end in spans[1:]:
        if start > reach and (best is None or start - reach > best[1] - best[0]):
            best = (reach, start)
        reach = max(reach, end)
    over_origin = length - reach + spans[0][0]
    return best if best is not None and best[1] - best[0] > over_origin else None


def origin_crossing(parts: Sequence[_Span], origin_length: int | None) -> OriginCrossing | None:
    """Where a location crosses its record's origin, on a record whose spans may cross it (``origin_length``, its
    length; ``spans_may_cross_origin``), else None. Connector lines (``connector_gaps``) read it."""

    gap = uncovered_gap(parts, origin_length) if origin_length is not None else None
    return OriginCrossing(origin_length, *gap) if gap is not None and origin_length is not None else None


def connector_gaps(
    previous: _Span, part: _Span, previous_strand: object, strand: object, crossing: OriginCrossing | None
) -> list[tuple[int, int]]:
    """The half-open gaps a connector line spans from one part to the next (Biopython strands).

    The gap along the first part's strand, as the location reads; parts that touch or overlap have none. Parts
    of one strand in the other order say nothing by their order: they connect over the origin, in two pieces,
    only when they lie on either side of the gap a location over the origin leaves out (``crossing``), and
    otherwise not at all. The location path and the display path both draw these gaps; a display start or
    wrapping adds and drops no connection.
    """

    if previous_strand == 1 and previous.end < part.start:
        return [(previous.end, part.start)]
    if previous_strand == -1 and part.end < previous.start:
        return [(part.end, previous.start)]
    if crossing is None or previous_strand != strand or strand not in (1, -1):
        return []
    low, high = (previous, part) if previous.start <= part.start else (part, previous)
    if not (low.end <= crossing.gap_start and crossing.gap_end <= high.start):
        return []
    pieces = [(high.end, crossing.length), (0, low.start)]
    return [(start, end) for start, end in (pieces if strand == 1 else pieces[::-1]) if start < end]


def project_feature_parts(
    exons: List[SimpleLocation],
    transform: RecordDisplayTransform,
    *,
    origin_length: int | None,
    is_trans_spliced: bool = False,
) -> tuple[FeatureDisplayPart, ...]:
    """Adapt true local boundaries once; retain unknown-strand rectangle glyphs.

    ``origin_length`` is the record's length when its spans may cross its origin, else None (``connector_gaps``).

    The +1 traversal for an undirected block orders geometry only. It does not
    assert biological strand or terminals, and never changes the source feature.
    """
    if not transform.projects:
        raise ValueError("Unset feature geometry must use the existing local path")
    crossing = None if is_trans_spliced else origin_crossing(exons, origin_length)
    parts = [SourceInterval(
        int(exon.start), int(exon.end),
        exon.strand if exon.strand in (-1, 1) else 1,
        index,
        biological_start=exon.strand in (-1, 1) and not partial.start,
        biological_end=exon.strand in (-1, 1) and not partial.end,
        partial_start=partial.start,
        partial_end=partial.end,
    ) for index, (exon, partial) in enumerate(zip(exons, location_partial_ends(exons), strict=True))]
    # A projecting transform (checked above) makes project_local_parts return DisplayFragments.
    by_part: dict[int, list[DisplayFragment]] = {}
    for fragment in cast("tuple[DisplayFragment, ...]", transform.project_local_parts(parts)):
        by_part.setdefault(fragment.part_index, []).append(fragment)
    result: list[FeatureDisplayPart] = []
    previous_index = None
    for index, part in enumerate(parts):
        if part.start == part.end:
            continue
        if previous_index is not None and not is_trans_spliced:
            previous = parts[previous_index]
            gaps = connector_gaps(previous, part, exons[previous_index].strand, exons[index].strand, crossing)
            for start, end in gaps:
                connector = SourceInterval(start, end, part.strand, index, False, False)
                result.extend(FeatureDisplayPart("line", get_strand(exons[index].strand), f)
                              for f in cast(
                                  "tuple[DisplayFragment, ...]",
                                  transform.project_local_parts((connector,)),
                              ))
        result.extend(FeatureDisplayPart("block", get_strand(exons[index].strand), fragment)
                      for fragment in by_part.get(index, ()))
        previous_index = index
    return tuple(result)


def get_exon_coordinate(
    exon_line: SimpleLocation,
    previous_exon_count: int,
    last_or_not: bool,
    partial: PartialEnds = NOT_PARTIAL,
) -> Tuple[int, FeatureLocationPart]:
    """
    Calculates the coordinate of an exon.

    A partial 3' end is not a known biological end, so the block is not the
    arrowhead block; a partial end opens the block's edge on that side.

    Returns:
        (exon_count, exon_coordinate)
    """
    exon_strand: Strand = get_strand(exon_line.strand)
    exon_count: int = previous_exon_count + 1
    exon_id: str = str(exon_count).zfill(3)
    # Half-open, 0-based bounds, as the display path projects them: a part at the first base starts at 0.
    exon_start = int(exon_line.start)
    exon_end = int(exon_line.end)
    open_start, open_end = strand_ordered(partial.start, partial.end, exon_line.strand)
    exon_coordinate = FeatureLocationPart(
        "block", exon_id, exon_strand, exon_start, exon_end, last_or_not and not partial.end,
        open_start, open_end,
    )
    return exon_count, exon_coordinate


def get_intron_coordinate(
    previous_exon: FeatureLocationPart,
    current_exon: FeatureLocationPart,
    previous_intron_count: int,
    strands: tuple[object, object],
    crossing: OriginCrossing | None,
) -> tuple[int, list[FeatureLocationPart]]:
    """
    The connector line from one part to the next, on the first part's strand (``strands``: the two parts'
    Biopython strands).

    The line spans the half-open gaps ``connector_gaps`` gives. Parts that touch or overlap (e.g. an RNA editing
    frameshift such as join(1850..2340,2339..3026)) have none.
    """
    intron_count: int = previous_intron_count + 1
    intron_id: str = str(intron_count).zfill(3)
    return intron_count, [
        FeatureLocationPart("line", intron_id, previous_exon.strand, start, end, False)
        for start, end in connector_gaps(previous_exon, current_exon, *strands, crossing)
    ]


def get_exon_and_intron_coordinates(
    exons: List[SimpleLocation], is_trans_spliced: bool = False, *, origin_length: int | None
) -> FeatureLocation:
    """
    Processes a list of exons to determine the coordinates of both exons and introns.
    Handles trans-spliced features by drawing only exons without introns. ``origin_length`` is the record's
    length when its spans may cross its origin, else None (``connector_gaps``).
    """
    coordinates_list: FeatureLocation = []
    exon_count: int = 0
    intron_count: int = 0
    num_exons: int = len(exons)

    if num_exons == 0:
        return []
    partial_ends = location_partial_ends(exons)

    if num_exons == 1:
        exon_count, exon_coord = get_exon_coordinate(exons[0], exon_count, True, partial_ends[0])
        coordinates_list.append(exon_coord)

    elif num_exons > 1 and not is_trans_spliced:
        crossing = origin_crossing(exons, origin_length)
        for i, exon in enumerate(exons):
            is_last_exon = i == num_exons - 1
            exon_count, exon_coord = get_exon_coordinate(exon, i, is_last_exon, partial_ends[i])

            if i > 0:
                previous_exon_coord = coordinates_list[-1]
                intron_count, intron_coords = get_intron_coordinate(
                    previous_exon_coord, exon_coord, intron_count, (exons[i - 1].strand, exon.strand), crossing
                )
                coordinates_list.extend(intron_coords)

            coordinates_list.append(exon_coord)

    elif num_exons > 1 and is_trans_spliced:
        for i, exon in enumerate(exons):
            is_last = i == num_exons - 1
            exon_count, exon_coord = get_exon_coordinate(exon, i, is_last, partial_ends[i])
            coordinates_list.append(exon_coord)

    return coordinates_list


def get_coordinate(coord: FeatureLocationPart) -> Tuple[str, Strand, int, int]:
    """
    Extracts and returns basic feature coordinates.
    """
    feat_type: str = coord.kind
    feat_strand: Strand = coord.strand
    feat_start: int = coord.start
    feat_end: int = coord.end
    return feat_type, feat_strand, feat_start, feat_end


def get_strand(strand_value: int) -> Strand:
    """
    Converts strand value to a string representation.
    """
    if strand_value == 1:
        strand: Strand = "positive"
    elif strand_value == -1:
        strand = "negative"
    else:
        strand = "undefined"
    return strand


__all__ = [
    "PartialEnds",
    "OriginCrossing",
    "connector_gaps",
    "origin_crossing",
    "uncovered_gap",
    "get_coordinate",
    "get_exon_and_intron_coordinates",
    "get_exon_coordinate",
    "get_intron_coordinate",
    "get_strand",
    "location_partial_ends",
]


