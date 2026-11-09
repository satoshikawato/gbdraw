"""The sliding windows of the GC content, GC skew and Depth tracks: where each window lies and which bases it reads.

The one owner of the window rule. The value of the sample at local position ``p`` is that of the ``window`` bases
centred on it, ``[p - window // 2, p - window // 2 + window)`` (OV-2870, OV-2919). Past an end of the drawn record
the window reads what lies there on the molecule (OV-2890):

- a circular record's windows wrap round its origin;
- a crop (a region, or one window of a record) reads its flanking bases from the record it was cut from, which it
  keeps as its :class:`WindowSource`; past that record's ends it wraps or stops as the record does;
- a linear record's windows stop at its ends, so a window there reads fewer bases.

The GC tracks count bases (``gbdraw.analysis.skew``) and the Depth track sums depths (``gbdraw.analysis.depth``)
over the same windows, read through the same :class:`SourcePiece` map.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord
from numpy.typing import NDArray

from gbdraw.core.record_metadata import _read_coord_map

WINDOW_SOURCE_ANNOTATION = "gbdraw_window_source"


@dataclass(frozen=True)
class WindowSource:
    """The whole sequence a record's windows read: its own, or for a crop the record it was cut from.

    Index ``i`` of ``sequence`` is the source coordinate ``base + step * i`` (the coordinate map of the record
    the sequence belongs to, ``gbdraw.core.record_metadata``).
    """

    sequence: Seq | MutableSeq
    base: int
    step: Literal[-1, 1]


def keep_window_source(crop: SeqRecord, record: SeqRecord) -> None:
    """Let ``crop``, cut from ``record``, read its flanking bases from ``record``.

    A crop of a crop already holds the first record's source (its annotations are copied), and keeps it.
    """

    if not isinstance(crop.annotations.get(WINDOW_SOURCE_ANNOTATION), WindowSource):
        crop.annotations[WINDOW_SOURCE_ANNOTATION] = window_source(record)


def window_source(record: SeqRecord) -> WindowSource:
    """The sequence ``record``'s windows read: the source a crop keeps, else the record's own sequence."""

    stored = (getattr(record, "annotations", None) or {}).get(WINDOW_SOURCE_ANNOTATION)
    if isinstance(stored, WindowSource):
        return stored
    base, step = _read_coord_map(record)
    return WindowSource(record.seq, base, 1 if step > 0 else -1)


@dataclass(frozen=True)
class SourcePiece:
    """Local positions ``local .. local + length - 1`` of a record read source indexes ``index``, ``index +
    direction``, ...: one run of the source, without a turn round its origin."""

    local: int
    index: int
    length: int
    direction: Literal[-1, 1]

    @property
    def low(self) -> int:
        """The first source index of the piece in source order."""
        return self.index if self.direction > 0 else self.index - self.length + 1


@dataclass(frozen=True)
class SlidingWindows:
    """The windows of one record's samples on its local axis, which runs on past its ends into its source.

    ``samples`` are 0, ``step``, ... below the record length: where each value is drawn. The window of
    ``samples[k]`` reads local positions ``starts[k] .. ends[k] - 1``; ``pieces`` map the local positions
    ``starts.min() .. ends.max() - 1`` to the source, in local order.
    """

    samples: NDArray[np.int64]
    starts: NDArray[np.int64]
    ends: NDArray[np.int64]
    source: WindowSource
    pieces: tuple[SourcePiece, ...]

    @property
    def bases(self) -> NDArray[np.int64]:
        """How many bases each window reads: ``window``, or fewer at the ends of a linear record."""
        return self.ends - self.starts


def sliding_windows(record: SeqRecord, window: int, step: int, *, circular: bool) -> SlidingWindows:
    """The window of each sample of ``record`` (``window`` > 0, ``step`` > 0, a record of one base or more).

    ``circular`` says whether the molecule is circular: the record drawn as circular (its bound
    ``RecordDisplayTransform.is_circular``), for a crop the record it was cut from.
    """

    length = len(record.seq)
    source = window_source(record)
    size = len(source.sequence)
    base, record_step = _read_coord_map(record)
    # Local position x reads source index first + direction * x.
    first = (base - source.base) * source.step
    direction: Literal[-1, 1] = 1 if record_step * source.step > 0 else -1
    samples = np.arange(0, length, int(step), dtype=np.int64)
    starts = samples - int(window) // 2
    ends = starts + int(window)
    if not circular:
        # A linear molecule has bases at source indexes 0 .. size - 1 only.
        low, high = (-first, size - first) if direction > 0 else (first - size + 1, first + 1)
        starts = np.clip(starts, low, high)
        ends = np.clip(ends, low, high)
    pieces: list[SourcePiece] = []
    local, stop = int(starts.min()), int(ends.max())
    while local < stop:
        index = (first + direction * local) % size
        run = min(stop - local, size - index if direction > 0 else index + 1)
        pieces.append(SourcePiece(local, index, run, direction))
        local += run
    return SlidingWindows(samples, starts, ends, source, tuple(pieces))


__all__ = [
    "WINDOW_SOURCE_ANNOTATION",
    "SlidingWindows",
    "SourcePiece",
    "WindowSource",
    "keep_window_source",
    "sliding_windows",
    "window_source",
]
