"""A whole sequence read a span at a time: the one loop of every reader that must not copy a record sequence.

A count, a digest or a comparison of a whole record sequence reads it here, so it
holds one span of its bases at a time. A slice of a record sequence read from its
input file reads only the bases of the slice (:mod:`gbdraw.io.record_sequences`).
"""

from __future__ import annotations

from collections.abc import Iterator

from Bio.Seq import MutableSeq, Seq

SPAN_BASES = 1 << 20  # the bases a whole-sequence reader takes at a time, unless it names its own size


def sequence_spans(seq: Seq | MutableSeq | bytes, size: int | None = None) -> Iterator[bytes]:
    """``seq``'s bases ``size`` (``SPAN_BASES`` when not given) at a time, in order."""

    step = size or SPAN_BASES
    for start in range(0, len(seq), step):
        yield bytes(seq[start:start + step])


__all__ = ["SPAN_BASES", "sequence_spans"]
