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

from dataclasses import dataclass, replace
from typing import Literal

import numpy as np
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord
from numpy.typing import NDArray

from gbdraw.core.record_metadata import (
    _iter_source_features,
    _read_coord_map,
    RecordContent,
    _source_feature_index,
)

WINDOW_SOURCE_ANNOTATION = "gbdraw_window_source"


@dataclass(frozen=True, eq=False)
class WindowSource:
    """The whole record a record's windows read: itself, or for a crop the record it was cut from.

    Index ``i`` of ``sequence`` is the source coordinate ``base + step * i`` (the coordinate map of ``record``,
    ``gbdraw.core.record_metadata``). A crop holds its source by reference: a copy of the crop shares it, and a
    Session writes a crop of an in-memory record as that record and the crop's window
    (``gbdraw.session_request_codec``), so a replay reads the same flanking bases.
    """

    record: SeqRecord
    base: int
    step: Literal[-1, 1]
    # The crop's content when it was cut, while it is a crop of ``record`` as cut: ``None`` when it was cut from a
    # crop changed since that crop was cut.
    crop_content: RecordContent | None = None
    # The source feature indexes of the features the crop was cut with (``keep_window_source``): the one owner of which
    # features the crop cut. A feature that carries any other index (one copied from another crop of ``record``, say)
    # is a feature added to the crop, live and in a Session (``record_metadata._holder_of``, OV-3256). Empty for a
    # record that is not a crop.
    cut_indexes: frozenset[int] = frozenset()

    @property
    def sequence(self) -> Seq | MutableSeq:
        return self.record.seq

    def holds_for(self, crop: SeqRecord) -> bool:
        """Whether ``crop`` is still ``record``'s crop as it was cut, so cutting ``record`` again redraws it."""
        return self.crop_content is not None and self.crop_content.same_as(RecordContent(crop))

    def __deepcopy__(self, _memo: object) -> WindowSource:
        # Nothing writes a source record after it is cut: copies of a crop share it, as they share its sequence.
        return self


def keep_window_source(crop: SeqRecord, record: SeqRecord) -> None:
    """Let ``crop``, cut from ``record`` (and reverse-complemented when the cut says so), read its flanking bases
    from ``record``, and mark it as ``record``'s crop as cut.

    A crop of a crop reads the first record's bases; it is that record's crop as cut only while the crop it was
    cut from is, and it cut only the features that crop cut.
    """

    stored = record.annotations.get(WINDOW_SOURCE_ANNOTATION)
    as_cut = not isinstance(stored, WindowSource) or stored.holds_for(record)
    carried = (_source_feature_index(feature) for feature in _iter_source_features(crop.features))
    cut = frozenset(index for index in carried if index is not None)
    crop.annotations[WINDOW_SOURCE_ANNOTATION] = replace(
        window_source(record),
        crop_content=RecordContent(crop) if as_cut else None,
        cut_indexes=cut & stored.cut_indexes if isinstance(stored, WindowSource) else cut,
    )


def reversed_window_source(reversed_record: SeqRecord, record: SeqRecord) -> None:
    """Keep ``reversed_record``, the reverse complement of the crop ``record``, a crop as cut while ``record`` is."""

    stored = record.annotations.get(WINDOW_SOURCE_ANNOTATION)
    if isinstance(stored, WindowSource):
        reversed_record.annotations[WINDOW_SOURCE_ANNOTATION] = replace(
            stored, crop_content=RecordContent(reversed_record) if stored.holds_for(record) else None
        )


def cut_content(crop: SeqRecord) -> RecordContent | None:
    """The content ``crop`` had when it was cut (or reverse-complemented as a crop), else ``None``.

    Only a caller that knows the crop changed since in nothing but gbdraw's annotations, which the content leaves
    out, may take it as the crop's content now: it saves taking the crop's text again.
    """

    stored = (getattr(crop, "annotations", None) or {}).get(WINDOW_SOURCE_ANNOTATION)
    return stored.crop_content if isinstance(stored, WindowSource) else None


def window_source(record: SeqRecord) -> WindowSource:
    """The sequence ``record``'s windows read: the source a crop keeps, else the record's own sequence."""

    stored = (getattr(record, "annotations", None) or {}).get(WINDOW_SOURCE_ANNOTATION)
    if isinstance(stored, WindowSource):
        return stored
    base, step = _read_coord_map(record)
    return WindowSource(record, base, 1 if step > 0 else -1)


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
    ``min(starts.min(), 0) .. max(ends.max(), length) - 1`` (the windows and the drawn record) to the source, in
    local order.
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
    # The pieces cover the windows and the whole drawn record: a step wider than half the window leaves bases of
    # the record in no window, and a reader of the drawn bases (Depth's count) still maps them.
    local, stop = min(int(starts.min()), 0), max(int(ends.max()), length)
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
    "cut_content",
    "keep_window_source",
    "sliding_windows",
    "window_source",
]
