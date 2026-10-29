#!/usr/bin/env python
# coding: utf-8

"""Comparison intervals that cross a circular record's origin.

An arc is a 0-based half-open ``(start, end)`` pair. ``start > end`` means the
arc crosses the record origin: it covers ``[start, L)`` and then ``[0, end)``
for a record of length ``L``. Comparison rows use directed 1-based endpoints,
whose order gives the strand (``start > end`` on the minus strand); a span that
crosses the origin reverses that order, so its endpoints disagree with its
strand (:func:`directed_arc`). A row cannot hold such a span, so :func:`split_link_at_origin` cuts a
link into rows that each stay inside the record.
"""

from __future__ import annotations

from math import inf
from typing import NamedTuple, Sequence

Arc = tuple[int, int]


class LinkSide(NamedTuple):
    """One side of a comparison link: directed 1-based endpoints on a record."""

    start: int
    end: int
    crosses_origin: bool
    record_length: int


def location_origin_arc(
    parts: Sequence[tuple[int, int, int | None]], record_length: int,
) -> Arc | None:
    """The arc of a location whose parts, in biological order, cross the origin.

    ``parts`` are 0-based half-open ``(start, end, strand)`` triples on a
    record of ``record_length`` bases. A location crosses the origin when its
    parts all share one strand and step back exactly once, from a part that
    ends at the record's end to a part that starts at 0 (``join(L-x..L,1..y)``
    and its complement). Any other location, such as a minus-strand location
    listed in ascending order or a trans-spliced one, returns None and keeps
    its linear extent.
    """

    if len(parts) < 2:
        return None
    strands = {strand for _start, _end, strand in parts}
    if len(strands) != 1:
        return None
    (strand,) = strands
    if strand not in (-1, 1):
        return None
    # Genomic (5' to 3' on the plus strand) order: a step back there is a wrap.
    ordered = list(parts) if strand == 1 else list(reversed(parts))
    steps = [
        (previous, current)
        for previous, current in zip(ordered, ordered[1:])
        if int(current[0]) < int(previous[0])
    ]
    if len(steps) != 1:
        return None
    ((last, first),) = steps
    if int(last[1]) != int(record_length) or int(first[0]) != 0:
        return None
    start, end = int(ordered[0][0]), int(ordered[-1][1])
    return (start, end) if start > end else None


def location_center(parts: Sequence[tuple[int, int, int | None]], record_length: int) -> float:
    """The source boundary at the middle of a location.

    A location across the origin centres on the middle of its arc
    (:func:`location_origin_arc`), wrapped into ``[0, record_length]``; any
    other location centres on the middle of its linear extent.
    """

    arc = location_origin_arc(parts, record_length)
    if arc is None:
        return (min(part[0] for part in parts) + max(part[1] for part in parts)) / 2.0
    center = (arc[0] + arc[1] + int(record_length)) / 2.0
    return center - record_length if center > record_length else center


def cover_arcs(arcs: Sequence[Arc]) -> Arc | None:
    """The smallest arc that covers ``arcs``; None when they cover the whole record.

    Without an arc across the origin the cover is the linear extent, as before.
    With one, the cover crosses the origin too and leaves out the largest gap
    between the arcs.
    """

    if not any(start > end for start, end in arcs):
        return min(start for start, _end in arcs), max(end for _start, end in arcs)
    pieces = sorted(
        piece
        for start, end in arcs
        for piece in (((start, inf), (0, end)) if start > end else ((start, end),))
    )
    gap: tuple[float, float] | None = None
    reach = pieces[0][1]
    for start, end in pieces[1:]:
        if start > reach and (gap is None or start - reach > gap[1] - gap[0]):
            gap = (reach, start)
        reach = max(reach, end)
    if gap is None:
        return None
    return int(gap[1]), int(gap[0])


def directed_arc(
    start: int, end: int, strand: int | None, record_length: int | None,
) -> tuple[Arc, bool]:
    """The arc of directed 1-based endpoints and whether it crosses the origin.

    Endpoints that disagree with their strand cross the origin when the record
    length is known and that arc is the shorter of the two the endpoints bound;
    otherwise they keep their linear extent, as every reader did before.
    """

    linear = abs(int(end) - int(start)) + 1
    crosses = (
        ((strand == 1 and start > end) or (strand == -1 and start < end))
        and record_length is not None
        and int(record_length) - linear + 2 < linear
    )
    forward = (start > end) if crosses else (start <= end)
    if forward:
        return (int(start) - 1, int(end)), crosses
    return (int(end) - 1, int(start)), crosses


def arc_link_endpoints(arc: Arc, forward: bool) -> tuple[int, int]:
    """Directed 1-based endpoints of an arc, traversed forward or backward."""

    start, end = arc
    return (start + 1, end) if forward else (end, start + 1)


def _unwrapped(side: LinkSide) -> tuple[float, float]:
    """Boundaries where the traversal begins and ends, past ``L`` when it crosses."""

    length = side.record_length if side.crosses_origin else 0
    forward = (side.start > side.end) if side.crosses_origin else (side.start <= side.end)
    if forward:
        return float(side.start - 1), float(side.end + length)
    return float(side.start + length), float(side.end - 1)


def split_link_at_origin(query: LinkSide, subject: LinkSide) -> tuple[tuple[int, int, int, int], ...]:
    """Directed ``(qstart, qend, sstart, send)`` rows of a link, cut at each origin.

    Both sides advance along one endpoint-linear parameter, the geometry of a
    drawn ribbon; each side that crosses its origin cuts the link there, and the
    other side is cut at the matching interpolated base.
    """

    sides = (query, subject)
    bounds = [_unwrapped(side) for side in sides]
    cuts = {0.0, 1.0}
    for side, (begin, finish) in zip(sides, bounds):
        if side.crosses_origin:
            cuts.add((side.record_length - begin) / (finish - begin))
    ordered = sorted(cuts)
    rows: list[tuple[int, int, int, int]] = []
    for left, right in zip(ordered, ordered[1:]):
        row: list[int] = []
        for side, (begin, finish) in zip(sides, bounds):
            first = begin + (finish - begin) * left
            last = begin + (finish - begin) * right
            low, high = sorted((first, last))
            if side.crosses_origin and (low + high) / 2 > side.record_length:
                low, high = low - side.record_length, high - side.record_length
                first, last = first - side.record_length, last - side.record_length
            low_base = int(round(low))
            if side.crosses_origin:
                low_base = min(low_base, side.record_length - 1)
            high_base = max(int(round(high)), low_base + 1)
            row.extend((low_base + 1, high_base) if first <= last else (high_base, low_base + 1))
        rows.append((row[0], row[1], row[2], row[3]))
    return tuple(rows)


def rows_inside_records(
    row: dict[str, object], query: LinkSide, subject: LinkSide,
) -> list[dict[str, object]]:
    """``row`` with the link's endpoints, as one row per piece inside the records.

    A link that crosses no origin keeps one row; every piece of one that does
    carries the same metadata (block, proteins, scores) as the link.
    """

    if not (query.crosses_origin or subject.crosses_origin):
        return [row]
    return [
        {**row, "qstart": qstart, "qend": qend, "sstart": sstart, "send": send}
        for qstart, qend, sstart, send in split_link_at_origin(query, subject)
    ]


__all__ = [
    "Arc",
    "LinkSide",
    "arc_link_endpoints",
    "cover_arcs",
    "directed_arc",
    "location_center",
    "location_origin_arc",
    "rows_inside_records",
    "split_link_at_origin",
]
