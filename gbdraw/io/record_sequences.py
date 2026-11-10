"""Record sequences read from their input file: the one owner of record sequence access.

Every sequence input (a GenBank, DDBJ or EMBL flat file, a FASTA file, the FASTA
part of a GFF3 file) is read in one pass over its bytes. The pass hands the lines
that are not sequence to Biopython or BCBio-GFF, which build the record headers and
features as before, and records, for each record, where its sequence lies in the
file (``_SequenceTable``) instead of building it. A record's ``seq`` is then a
``Seq`` over a :class:`RecordSequence`, which reads bases on demand from the file
through one block cache (``_BlockCache``). Its bytes equal the ``Seq``
Biopython builds from the same lines: a flat file's sequence in upper case, without
the coordinate column and the spaces; a FASTA sequence as written, without
whitespace.

A gzip input is expanded as the pass reads it (:class:`gbdraw.io.input_text.GzipExpansion`,
with its cap), and nothing of it is written: the pass keeps a point of the expansion
every ``_POINT_BYTES`` expanded bytes, and a read at an offset expands from the
point before it, or goes on from where the last read of that file stopped. The pass
counts bytes itself and never asks the file for its position, and every read below
a record sequence goes through ``_SequenceFile.read_at``, the one reader, so the
Web can put its own range reader there.

The sequence lines are checked in runs of lines of one length. A run laid out as
writers lay out sequence lines (coordinates, then letters and spaces) is read
with one ``bytes.translate``; any other line is read with the rule Biopython
applies to each line, so both give Biopython's bytes. A cache budget decides how
many decoded bases stay in memory, never what a read returns: a short read (a
feature's bases) keeps the chunks it decoded, and a read of ``_BLOCK_BYTES`` or
more (a pass over a span) is decoded in bulk and keeps nothing.

The same pass counts the bases (``gbdraw.core.base_counts``): the one producer of
a record's base counts. The letter check of a FASTA entry, the GC percent of a
record or a crop and the GC and skew windows read those counts
(:func:`known_base_counts_at`, :func:`first_base_beyond`), reading only the bases
at the edges of a span again.
"""

from __future__ import annotations

import bisect
import contextlib
import copy
import functools
import itertools
import os
import re
import threading
import warnings
import weakref
from array import array
from collections import OrderedDict, deque
from collections.abc import Callable, Generator, Iterable, Iterator
from dataclasses import dataclass, field
from os import PathLike, fspath
from typing import IO, Literal, NamedTuple, TypeVar

import numpy as np
from numpy.typing import NDArray

from Bio import BiopythonParserWarning
from Bio.Seq import MutableSeq, Seq, SequenceDataAbstractBaseClass, _dna_complement_table, _rna_complement_table
from Bio.SeqRecord import SeqRecord

from ..core.base_counts import (
    BASES,
    BEYOND_NUCLEOTIDE,
    BEYOND_PROTEIN,
    KNOWN,
    NUCLEOTIDE,
    OTHER,
    PROTEIN,
    ROW_BASES,
    SLICE_BYTES,
    BASES_MATRIX,
    CountedRows,
    RowCounts,
    byte_counts,
    count_matrix,
    count_rows,
    known_prefix_counts,
)
from ..core.sequence_spans import sequence_spans
from ..exceptions import InputFileError
from . import input_text
from .input_text import (
    GZIP_ERRORS,
    GZIP_START,
    FlatFileFormat,
    GzipExpansion,
    GzipPoint,
    in_sequence_encoding,
    is_gzip_file,
    open_input_bytes,
)

SequenceFormat = Literal["genbank", "embl", "fasta"]

_BLOCK_BYTES = 1 << 20  # sequential reads of the parse pass; a read this long or longer is decoded in bulk
_CHUNK_BYTES = 64 << 10  # most raw bytes of one indexed chunk (whole lines)
_CACHE_BYTES = 256 << 20  # decoded bases the block cache keeps, all input files together
_OPEN_FILES = 16  # input files kept open for reads at offsets, all records together
_POINT_BYTES = 4 << 20  # expanded bytes of a gzip input between two points a read may start from
_AHEAD_BLOCKS = 4  # blocks a gzip input is expanded ahead of its parse pass, on a thread of their own
_BYTE_ORDER_MARK = "﻿"
_DIGITS = b"0123456789"
_LETTERS = bytes(range(ord("A"), ord("Z") + 1)) + bytes(range(ord("a"), ord("z") + 1)) + b"*-."
_UPPER = bytes.maketrans(bytes(range(ord("a"), ord("z") + 1)), bytes(range(ord("A"), ord("Z") + 1)))
_FLAT_DELETE = b" \r\n" + _DIGITS
_FASTA_DELETE = b" \t\r\n"

# Chunk kinds: how the raw bytes of a chunk become bases.
_FLAT_RUN = 0  # coordinates, letters and spaces: delete digits and whitespace, upper case
_FASTA_RUN = 1  # ASCII FASTA lines: delete whitespace
_GENBANK_LINES = 2  # Biopython's GenBank rule, line by line
_EMBL_LINES = 3  # Biopython's EMBL rule, line by line
_FASTA_LINES = 4  # FASTA lines that are not ASCII: Biopython's text rule
# A run is counted from its raw bytes, a chunk of it per histogram call (any split of a run decodes alone); lines
# read line by line are counted from their bases.
_RUN_MATRICES = {_FLAT_RUN: count_matrix(_FLAT_DELETE), _FASTA_RUN: count_matrix(_FASTA_DELETE)}
assert _CHUNK_BYTES == SLICE_BYTES and _CHUNK_BYTES <= 1 << 16  # a row's offset in its chunk is a ``uint16``
_FLAT_KEPT = np.ones(256, dtype=bool)
_FLAT_KEPT[list(_FLAT_DELETE)] = False
_FLAT_NOT_LETTER = _FLAT_KEPT.copy()
_FLAT_NOT_LETTER[list(_LETTERS)] = False
_EDGE_BATCH = 1 << 20  # most bases of edge rows counted together


class _Rows(NamedTuple):
    """Count rows of a sequence's bases (``gbdraw.core.base_counts``): ``counts[j]`` (A, C, G and T, U as T,
    either case) counts row ``j``, bases ``j * ROW_BASES`` up to the next row or ``end``.
    They hold no reader of the bases, so a holder that keeps its rows is not kept alive by them."""

    counts: NDArray[np.uint16]
    end: int


def _alphabets(totals: NDArray[np.int64]) -> NDArray[np.uint8]:
    """What the bases counted in each row of ``totals`` (7 columns) hold (``gbdraw.core.base_counts``)."""

    return np.where(totals[:, BEYOND_PROTEIN] > 0, OTHER,
                    np.where(totals[:, BEYOND_NUCLEOTIDE] > 0, PROTEIN, NUCLEOTIDE)).astype(np.uint8)

# Biopython's record start and sequence headers (Bio.GenBank.Scanner).
_GENBANK_RECORD_START = "LOCUS       "
_GENBANK_SEQUENCE_HEADERS = ("CONTIG", "ORIGIN", "BASE COUNT", "WGS", "TSA", "TLS")
_EMBL_RECORD_START = "ID   "
_EMBL_SEQUENCE_HEADERS = ("SQ", "CO")

_Read = TypeVar("_Read")


class _BlockCache:
    """The decoded bases read so far, least recently used dropped first past ``_CACHE_BYTES``.

    One cache serves every input file, so the budget bounds the memory of all of
    them. A block is one chunk.
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._blocks: OrderedDict[tuple[int, int, int], bytes] = OrderedDict()
        self._bytes = 0

    def get(self, key: tuple[int, int, int], build: Callable[[], bytes]) -> bytes:
        _release_pending()
        with self._lock:
            hit = self._blocks.get(key)
            if hit is not None:
                self._blocks.move_to_end(key)
                return hit
        value = build()
        dropped: list[bytes] = []  # freed after the lock is let go
        with self._lock:
            if key not in self._blocks:
                self._blocks[key] = value
                self._bytes += len(value)
                while self._bytes > _CACHE_BYTES and len(self._blocks) > 1:
                    dropped.append(self._blocks.popitem(last=False)[1])
                    self._bytes -= len(dropped[-1])
        return value

    def drop_file(self, file_key: int) -> None:
        """Forget the blocks of a file no record sequence reads any more."""

        dropped: list[bytes] = []
        with self._lock:
            for key in [key for key in self._blocks if key[0] == file_key]:
                dropped.append(self._blocks.pop(key))
                self._bytes -= len(dropped[-1])


_CACHE = _BlockCache()
_FILE_KEYS = iter(range(1, 1 << 62))

# What identifies the bytes a file holds: size, modification time (ns) and inode.
_Stamp = tuple[int, int, int]


def _stamp(handle: IO[bytes]) -> _Stamp:
    return _stamp_of(os.fstat(handle.fileno()))


def _stamp_of(status: os.stat_result) -> _Stamp:
    return status.st_size, status.st_mtime_ns, status.st_ino


class _Reader:
    """One input file open for reads at offsets: the bytes a record sequence reads, from the file as indexed.

    A plain file is read where its bytes lie. A gzip file is expanded from the point
    of its parse pass before the offset (``_SequenceFile.points``), or on from where
    the last read stopped when that is nearer: the last block expanded is kept, so
    reads in order (window edges, a sequence written out) expand each byte once.
    It refers to its input weakly (a read comes through a record, which holds it), so
    the input's finalizer runs when its last record goes, while the reader still
    sits among the open files (OV-3355).
    """

    def __init__(self, source: _SequenceFile) -> None:
        self._source = weakref.ref(source)
        self.handle = open_input_bytes(source.path)
        self._expansion: GzipExpansion | None = None
        self._at = 0  # the expanded offset of ``_block``
        self._block = b""

    def read_at(self, offset: int, size: int) -> bytes:
        source = self._source()
        assert source is not None  # the record reading holds its input
        if not source.gzip:
            self.handle.seek(offset)
            return self.handle.read(size)
        points = source.points
        point = points[bisect.bisect_right(points, offset, key=_expanded_at) - 1]
        if self._expansion is None or not (self._at <= offset and point.expanded <= self._at + len(self._block)):
            self.handle.seek(point.compressed)
            self._expansion = GzipExpansion(self.handle.read, source.name, start=point, verify=False)
            self._at, self._block = point.expanded, b""
        pieces: list[bytes] = []
        position, end = offset, offset + size
        while position < end:
            block_end = self._at + len(self._block)
            if position < block_end:
                pieces.append(self._block[position - self._at: end - self._at])
                position += len(pieces[-1])
            elif block := self._expansion.read(_BLOCK_BYTES):
                self._at, self._block = block_end, block
            else:
                break  # the file ends before the bytes indexed: the caller refuses it
        return pieces[0] if len(pieces) == 1 else b"".join(pieces)

    def close(self) -> None:
        self.handle.close()


def _expanded_at(point: GzipPoint) -> int:
    return point.expanded


class _Ahead:
    """What the expansion thread and the parse pass share: the blocks made and not yet read, and which side waits.

    Each side wakes the other only when that one waits (a flag set under the condition before it looks), so a
    block passes without a handoff while neither waits.
    """

    def __init__(self) -> None:
        self.made: deque[bytes] = deque()
        self.turn = threading.Condition()
        self.reader_waits = False
        self.maker_waits = False
        self.done = False
        self.stop = False
        self.error: BaseException | None = None

    def wake(self) -> None:
        with self.turn:
            self.turn.notify()


def _ahead(blocks: Iterator[bytes]) -> Iterator[bytes]:
    """``blocks``, made on a thread of their own while it is fewer than ``_AHEAD_BLOCKS`` ahead of the reader.

    zlib lets go of the interpreter lock while it expands, so a gzip input expands
    beside the parse pass's own work. Where no thread can start (Pyodide), the
    blocks are made in turn. A reader that stops early, or a refusal on the way,
    stops the thread before the file is closed.
    """

    shared = _Ahead()

    def make() -> None:
        try:
            for block in blocks:
                shared.made.append(block)
                if shared.reader_waits:
                    shared.wake()
                if len(shared.made) >= _AHEAD_BLOCKS:
                    with shared.turn:
                        shared.maker_waits = True
                        while len(shared.made) >= _AHEAD_BLOCKS and not shared.stop:
                            shared.turn.wait()
                        shared.maker_waits = False
                if shared.stop:
                    break
        except BaseException as error:  # handed to the reader, which raises it after the blocks before it
            shared.error = error
        finally:
            with shared.turn:
                shared.done = True
                shared.turn.notify_all()

    thread = threading.Thread(target=make, name="gbdraw-gzip-ahead", daemon=True)
    try:
        thread.start()
    except RuntimeError:  # no threads here
        yield from blocks
        return
    try:
        while True:
            if shared.made:
                block = shared.made.popleft()
                if shared.maker_waits:
                    shared.wake()
                yield block
                continue
            with shared.turn:
                shared.reader_waits = True
                while not shared.made and not shared.done:
                    shared.turn.wait()
                shared.reader_waits = False
                if not shared.made:
                    break
        if shared.error is not None:
            raise shared.error
    finally:
        with shared.turn:
            shared.stop = True
            shared.turn.notify_all()
        thread.join()


class _OpenFiles:
    """The input files open for reads at offsets: at most ``_OPEN_FILES``, least recently used closed first.

    A file is opened by its path when a read needs it, so the number of inputs is
    not bounded by the number of files a process may hold open. Every read checks
    that the file still holds the bytes its records were indexed from. The lock
    guards only the set of readers: a reader is taken out of it while it reads,
    and opening, reading and closing happen outside the lock.
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._handles: OrderedDict[int, _Reader] = OrderedDict()

    def read(self, source: _SequenceFile, offset: int, size: int) -> bytes:
        _release_pending()
        with self._lock:
            handle = self._handles.pop(source.key, None)
        try:
            if handle is None:
                handle = _Reader(source)
        except OSError as error:
            raise source.changed() from error
        try:
            if _stamp(handle.handle) != source.stamp:
                raise source.changed()
            data = handle.read_at(offset, size)
        except InputFileError:
            handle.close()
            raise
        except (OSError, *GZIP_ERRORS) as error:  # a gzip file that no longer expands as indexed changed too
            handle.close()
            raise source.changed() from error
        closing: list[_Reader] = []
        with self._lock:
            other = self._handles.pop(source.key, None)  # another thread's handle of the same file
            if other is not None:
                closing.append(other)
            self._handles[source.key] = handle
            while len(self._handles) > _OPEN_FILES:
                closing.append(self._handles.popitem(last=False)[1])
        for stale in closing:
            stale.close()
        if len(data) != size:
            raise source.changed()
        return data

    def close(self, file_key: int) -> None:
        with self._lock:
            handle = self._handles.pop(file_key, None)
        if handle is not None:
            handle.close()


_OPEN = _OpenFiles()

# Inputs whose last record sequence has gone. A finalizer runs at any allocation, also inside
# the readers' locked sections, so it only queues the input here (deque.append takes no lock);
# the readers release the queue before they take their locks.
_RELEASED: deque[int] = deque()


def _release(file_key: int) -> None:
    _RELEASED.append(file_key)


def _release_pending() -> None:
    """Close the files of the queued inputs and forget their blocks. Called outside every lock of this module."""

    while True:
        try:
            file_key = _RELEASED.popleft()
        except IndexError:
            return
        _OPEN.close(file_key)
        _CACHE.drop_file(file_key)


class _SequenceFile:
    """One input file's bytes: read in order by the parse pass, and at offsets by its record sequences.

    The bytes of a gzip input are its expansion: the parse pass expands it and keeps
    ``points`` of the expansion, from which a read at an offset expands again
    (:class:`_Reader`). A link is followed once: the records read the file it names,
    so they keep reading after the link goes (the Web Worker links each input into
    the folder of one request and removes the folder when the request ends).
    ``read_at`` is the one read below a record sequence. The first parse pass records
    the file's stamp; a later pass (the encoding retry) of a file that changed since,
    and a read of one, are refused with an ``InputFileError``, never answered with
    other bases.
    """

    def __init__(self, path: str | PathLike[str]) -> None:
        _release_pending()
        self.name = fspath(path)
        self.path = os.path.realpath(self.name)
        self.encoding = "utf-8"
        self.key = next(_FILE_KEYS)
        self.stamp: _Stamp | None = None
        self.gzip = is_gzip_file(self.path)
        # Where a read of a gzip input may start expanding, in expanded order: the start, then a point at least
        # every ``_POINT_BYTES`` the parse pass expanded.
        self.points: list[GzipPoint] = [GZIP_START]
        weakref.finalize(self, _release, self.key)
        self._tables = 0

    def table_key(self) -> int:
        """A key no other sequence table of this file has (its blocks' cache keys hold it)."""

        self._tables += 1
        return self._tables

    def blocks(self) -> Generator[bytes, None, None]:
        """The file's bytes from its start, in order, from a handle of the parse pass's own.

        A gzip input is expanded as the pass reads it, under the cap of
        :mod:`gbdraw.io.input_text`, which refuses only a byte past it that the pass
        asks for; a pass that stops early expands no more.
        """

        with open_input_bytes(self.path) as handle:
            stamp = _stamp(handle)
            if self.stamp is None:
                self.stamp = stamp
            elif stamp != self.stamp:
                raise self.changed()
            if self.gzip:
                yield from self._expanded(handle)
            else:
                while block := handle.read(_BLOCK_BYTES):
                    yield block
            if _stamp(handle) != self.stamp:
                raise self.changed()

    def _expanded(self, handle: IO[bytes]) -> Iterator[bytes]:
        """A gzip input's expansion, a point of it kept every ``_POINT_BYTES`` (the points of a later pass replace
        those of an earlier one)."""

        expansion = GzipExpansion(handle.read, self.name, limit=input_text.MAX_EXPANDED_INPUT_BYTES)
        points = [GZIP_START]
        self.points = points

        def expanded() -> Iterator[bytes]:
            while block := expansion.read(_BLOCK_BYTES):
                if expansion.expanded - points[-1].expanded >= _POINT_BYTES and (point := expansion.point()) is not None:
                    points.append(point)
                yield block

        yield from _ahead(expanded())

    def read_at(self, offset: int, size: int) -> bytes:
        """``size`` bytes from ``offset``, all within the file as it was indexed."""

        return _OPEN.read(self, offset, size)

    def unchanged(self) -> bool:
        """Whether the file still holds the bytes the parse pass indexed (one stat, no read)."""

        try:
            return _stamp_of(os.stat(self.path)) == self.stamp
        except OSError:
            return False

    def changed(self) -> InputFileError:
        return InputFileError(
            f"{self.name} changed after gbdraw read it, so its sequence can no longer be read from it. "
            "Read the file again.",
            diagnostic={"code": "INPUT_UNREADABLE"},
        )


@dataclass(eq=False)
class _SequenceTable:
    """Where one record's sequence lies in its file: chunks of whole lines, in order.

    ``starts[i]`` is the number of bases before chunk ``i``; ``offsets[i]`` and
    ``sizes[i]`` its raw bytes in the file; ``kinds[i]`` how they become bases;
    ``alphabets[i]`` what its bases hold (``gbdraw.core.base_counts``).

    The parse pass counts the bases as it indexes them (the one producer of a
    record's base counts), in rows of ``ROW_BASES`` bases from the first, and keeps
    where in its chunk each row's first base lies. A count of any span is then
    whole rows summed and the bases of its edge rows read (:func:`known_base_counts_at`).
    """

    source: _SequenceFile
    key: int
    offsets: array[int] = field(default_factory=lambda: array("q"))
    sizes: array[int] = field(default_factory=lambda: array("q"))
    starts: array[int] = field(default_factory=lambda: array("q"))
    kinds: bytearray = field(default_factory=bytearray)
    alphabets: bytearray = field(default_factory=bytearray)
    length: int = 0
    _counts: RowCounts = field(default_factory=RowCounts)
    _row_offsets: list[NDArray[np.uint16]] = field(default_factory=list)
    _rows: _Rows | None = None  # built once, when the parse pass has added every line (``finish``)

    def add(self, offset: int, size: int, kind: int, counted: CountedRows) -> None:
        """Add ``size`` raw bytes from file offset ``offset``, of one kind, with what was counted in them
        (:func:`gbdraw.core.base_counts.count_rows`: in the bytes for a run, in its bases for lines read line by line).

        A run is cut into chunks of ``_CHUNK_BYTES`` raw bytes; lines read line by line join the chunk before them
        up to the usual size.
        """

        self._counts.add(counted.rows)
        run = kind in _RUN_MATRICES
        # Where each row's first base lies in its chunk; lines read line by line are decoded whole.
        self._row_offsets.append(
            (counted.starts % _CHUNK_BYTES).astype(np.uint16) if run else np.zeros(len(counted.starts), dtype=np.uint16))
        totals = counted.histograms @ (_RUN_MATRICES[kind] if run else BASES_MATRIX)
        last = len(self.kinds) - 1
        if (not run and last >= 0 and self.kinds[last] == kind
                and self.offsets[last] + self.sizes[last] == offset and self.sizes[last] + size <= _CHUNK_BYTES):
            self.sizes[last] += size  # lines read line by line stay in chunks of the usual size
            self.alphabets[last] = max(self.alphabets[last], int(_alphabets(totals.sum(axis=0, keepdims=True))[0]))
            self.length += counted.bases
            return
        if not run:
            totals = totals.sum(axis=0, keepdims=True)
        # File offsets stay Python ints: numpy's default integer is 32 bits in Pyodide, and a file passes 2 GiB.
        chunk_offsets = range(offset, offset + size, _CHUNK_BYTES) if run else range(offset, offset + 1)
        self.offsets.extend(chunk_offsets)
        self.sizes.extend(min(_CHUNK_BYTES, offset + size - start) if run else size for start in chunk_offsets)
        self.starts.extend((self.length + np.cumsum(totals[:, BASES]) - totals[:, BASES]).tolist())
        self.kinds.extend(bytes((kind,)) * len(chunk_offsets))
        self.alphabets.extend(_alphabets(totals).tobytes())
        self.length += counted.bases

    def finish(self) -> None:
        """Join the count rows, once the parse pass has added every line, before any reader can ask for them."""

        self._row_offsets = [np.concatenate(self._row_offsets) if self._row_offsets else np.zeros(0, np.uint16)]
        self._rows = _Rows(self._counts.joined(), self.length)

    def rows(self) -> _Rows:
        """The count rows of the whole sequence."""

        assert self._rows is not None, "a record's rows are read after its parse pass"
        return self._rows

    def row_bases(self, row: int, size: int) -> bytes:
        """The first ``size`` bases of count row ``row`` (fewer than a row), read from the file and decoded alone:
        from the row's first base up to the next row's at most."""

        (row_offsets,) = self._row_offsets
        position = row * ROW_BASES
        stop, following = min(position + size, self.length), position + ROW_BASES
        index = bisect.bisect_right(self.starts, position) - 1
        skip = int(row_offsets[row])
        pieces: list[bytes] = []
        while position < stop:
            kind, begin, offset = self.kinds[index], self.starts[index], self.offsets[index]
            end = offset + self.sizes[index]
            if kind in _RUN_MATRICES:
                if following < (self.starts[index + 1] if index + 1 < len(self.starts) else self.length):
                    end = offset + int(row_offsets[row + 1])  # the next row begins in this chunk
                raw = self.source.read_at(offset + skip, end - offset - skip)
                bases = _decode_chunk(kind, raw, self.source.encoding, warn=False)[: stop - position]
            else:
                bases = _decode_chunk(kind, self.source.read_at(offset, end - offset), self.source.encoding,
                                      warn=False)[position - begin: stop - begin]
            pieces.append(bases)
            position += len(bases)
            index, skip = index + 1, 0
        return pieces[0] if len(pieces) == 1 else b"".join(pieces)

    def _block(self, index: int) -> bytes:
        """Chunk ``index``'s bases, through the block cache."""

        return _CACHE.get((self.source.key, self.key, index), lambda: self._decode(index, index + 1))

    def _decode(self, first: int, stop: int) -> bytes:
        """The bases of chunks ``first``..``stop - 1``, which are of one kind and adjacent in the file."""

        start = self.offsets[first]
        raw = self.source.read_at(start, self.offsets[stop - 1] + self.sizes[stop - 1] - start)
        return _decode_chunk(self.kinds[first], raw, self.source.encoding, warn=False)

    def read(self, start: int, stop: int) -> bytes:
        """Bases ``start``..``stop - 1`` (``0 <= start <= stop <= length``).

        A read of ``_BLOCK_BYTES`` bases or more is decoded in bulk, around the cache:
        it is a pass over a span (a digest, a crop written out), which keeps nothing.
        A shorter read (a feature's bases) reads its chunks through the cache.
        """

        if start >= stop:
            return b""
        if stop - start >= _BLOCK_BYTES:
            return self._bulk(start, stop)
        index = bisect.bisect_right(self.starts, start) - 1
        pieces: list[bytes] = []
        position = start
        while position < stop:
            bases = self._block(index)
            begin = self.starts[index]
            pieces.append(bases[position - begin: stop - begin])
            position = begin + len(bases)
            index += 1
        return pieces[0] if len(pieces) == 1 else b"".join(pieces)

    def _bulk(self, start: int, stop: int) -> bytes:
        """Bases ``start``..``stop - 1``, adjacent chunks of one kind decoded together, around the cache."""

        index = bisect.bisect_right(self.starts, start) - 1
        pieces: list[bytes] = []
        position = start
        while position < stop:
            last = index + 1
            while (last < len(self.starts) and self.starts[last] < stop and self.kinds[last] == self.kinds[index]
                   and self.offsets[last] == self.offsets[last - 1] + self.sizes[last - 1]
                   and self.offsets[last] + self.sizes[last] - self.offsets[index] <= _BLOCK_BYTES):
                last += 1
            bases = self._decode(index, last)
            begin = self.starts[index]
            pieces.append(bases[position - begin: stop - begin] if position > begin or stop - begin < len(bases) else bases)
            position = begin + len(bases)
            index = last
        return pieces[0] if len(pieces) == 1 else b"".join(pieces)


class _HeldBases:
    """Bases held in memory, read as a record's sequence table is read: their count rows are counted once, the
    first time a count is asked for (``gbdraw.core.base_counts``, on the same rows as a parse pass)."""

    __slots__ = ("bases", "length", "_rows", "_lock")

    def __init__(self, bases: bytes, rows: _Rows | None = None) -> None:
        self.bases = bases
        self.length = len(bases)
        self._rows = rows  # given: the rows a parse pass counted in these bases
        self._lock = threading.Lock()

    def read(self, start: int, stop: int) -> bytes:
        return self.bases[start:stop]

    def rows(self) -> _Rows:
        with self._lock:
            if self._rows is None:
                self._rows = _sliced_rows(self.read, self.length)
            return self._rows

    def row_bases(self, row: int, size: int) -> bytes:
        return self.bases[row * ROW_BASES: min(self.length, row * ROW_BASES + size)]


# What a record sequence reads: the sequence table of a record in its input file, or bases held in memory.
_Bases = _SequenceTable | _HeldBases


class RecordSequence(SequenceDataAbstractBaseClass):
    """A record's sequence, or a view of one, read on demand (``Seq(RecordSequence(...))``).

    It reads ``pieces`` of its bases (a record's sequence table in its input file,
    or bases held in memory) one after the other: a crop is one piece, two turns
    of a circular record two. On the reverse strand it reads them backwards; every
    byte read goes through ``translate`` (the complement table of a reverse
    complement), when it has one. ``len`` is known without reading; a slice reads
    only its bases and returns them (a plain ``Seq`` slice); ``bytes()`` of the
    whole reads them all. A copy is the same object: the bases never change.
    Views are made by :func:`sequence_view` and :func:`two_turns`, never here.
    """

    __slots__ = ("_bases", "_pieces", "_ends", "_reverse", "_translate")

    def __init__(
        self, bases: _Bases, pieces: tuple[tuple[int, int], ...] | None = None, *,
        reverse: bool = False, translate: bytes | None = None,
    ) -> None:
        self._bases = bases
        self._pieces = ((0, self._bases.length),) if pieces is None else pieces
        self._ends = tuple(itertools.accumulate(stop - start for start, stop in self._pieces))
        self._reverse = reverse
        self._translate = translate
        super().__init__()

    def __len__(self) -> int:
        return self._ends[-1] if self._ends else 0

    def _forward_pieces(self, start: int, stop: int) -> tuple[tuple[int, int], ...]:
        """The pieces of the bases that local positions ``start``..``stop - 1`` read, in forward order."""

        if self._reverse:
            length = len(self)
            start, stop = length - stop, length - start
        pieces: list[tuple[int, int]] = []
        index = bisect.bisect_right(self._ends, start)
        position = start
        while position < stop:
            low = self._ends[index - 1] if index else 0
            first, last = self._pieces[index]
            end = min(stop, self._ends[index])
            pieces.append((first + position - low, first + end - low))
            position = end
            index += 1
        return tuple(pieces)

    def _read(self, start: int, stop: int) -> bytes:
        """Local positions ``start``..``stop - 1`` (``0 <= start <= stop <= len``)."""

        if start >= stop:
            return b""
        bases = self._bases
        reads = [bases.read(first, last) for first, last in self._forward_pieces(start, stop)]
        data = reads[0] if len(reads) == 1 else b"".join(reads)
        del reads  # each step below holds the bases read and their next form only
        if self._translate is not None:
            data = data.translate(self._translate)
        return data[::-1] if self._reverse else data

    def __getitem__(self, key: int | slice) -> int | bytes:
        length = len(self)
        if isinstance(key, slice):
            start, stop, step = key.indices(length)
            if step == 1:
                return self._read(start, stop)
            indexes = range(start, stop, step)
            if not indexes:
                return b""
            low = min(indexes[0], indexes[-1])
            data = self._read(low, max(indexes[0], indexes[-1]) + 1)
            end = indexes.stop - low
            return data[indexes.start - low: end if end >= 0 else None: step]
        index = key + length if key < 0 else key
        if not 0 <= index < length:
            raise IndexError("sequence index out of range")
        return self._read(index, index + 1)[0]

    def __copy__(self) -> RecordSequence:
        return self

    def __deepcopy__(self, memo: dict[int, object]) -> RecordSequence:
        return self

    def __reduce__(self) -> tuple[type[bytes], tuple[bytes]]:
        # A pickle holds the bases: the file it reads from may not exist where it is loaded.
        return bytes, (bytes(self),)

    def __bytes__(self) -> bytes:
        return self._read(0, len(self))


def _viewed(seq: Seq | MutableSeq) -> RecordSequence | None:
    """``seq``'s bases as a record sequence a view can read, or ``None`` when they cannot be read so.

    A record sequence is itself; bases held as ``bytes`` are read in place, through
    one holder of those bytes, so their count rows are counted once however often
    the sequence (or a copy) is viewed or counted. Bases that may change (a ``MutableSeq``) and a
    sequence that is not defined everywhere are left to Biopython's own slices and
    reverse complements.
    """

    data = getattr(seq, "_data", None)
    if isinstance(data, RecordSequence):
        return data
    if isinstance(seq, Seq) and type(data) is bytes:
        return RecordSequence(_held(seq, data))
    return None


# The holder of the bytes of each ``Seq`` a view or a count has read, by the identity of the bytes (a copy of a
# ``Seq`` shares them): kept while the ``Seq`` that first read them lives, which keeps the bytes and so their
# identity. The weakref callback takes no lock: it may run inside any allocation.
_HELD: dict[int, tuple[_HeldRef, _HeldBases]] = {}
_HELD_LOCK = threading.Lock()  # one holder per bytes object, so their rows are counted once


class _HeldRef(weakref.ref[Seq]):
    """The ``Seq`` that first read some bytes, with the key of their holder: one object per held record."""

    __slots__ = ("key",)

    key: int

    def __new__(cls, seq: Seq, key: int) -> _HeldRef:
        ref = super().__new__(cls, seq, _forget)
        ref.key = key
        return ref


def _forget(ref: weakref.ref[Seq]) -> None:  # takes no lock: it runs as a finalizer
    assert isinstance(ref, _HeldRef)
    if _HELD.get(ref.key, (None,))[0] is ref:
        del _HELD[ref.key]


def _held(seq: Seq, data: bytes, rows: _Rows | None = None) -> _HeldBases:
    key = id(data)
    with _HELD_LOCK:
        entry = _HELD.get(key)
        if entry is not None and entry[1].bases is data:
            return entry[1]
        held = _HeldBases(data, rows)
        _HELD[key] = (_HeldRef(seq, key), held)
        return held


def _then(first: bytes | None, second: bytes | None) -> bytes | None:
    """The translate table that maps as ``first`` and then ``second`` (``None``: each byte as itself)."""

    if first is None or second is None:
        return first if second is None else second
    composed = first.translate(second)
    return None if composed == _IDENTITY else composed


_IDENTITY = bytes(range(256))


def sequence_view(
    seq: Seq | MutableSeq, start: int = 0, stop: int | None = None, *, strand: Literal[1, -1] = 1, rna: bool = False,
) -> Seq:
    """Bases ``start``..``stop - 1`` of ``seq``, on the reverse strand when ``strand`` is -1, read from ``seq``.

    The one way gbdraw cuts, reverse-complements and reads round a sequence
    without copying it: a crop, a reverse-complemented record, a window piece.
    The result reads ``seq``'s bases on demand; a view of a view reads the first
    sequence's bases. Its bytes equal ``seq[start:stop]``, and on the reverse
    strand Biopython's ``reverse_complement()`` of them (``reverse_complement_rna()``
    when ``rna``), byte for byte; ``start`` and ``stop`` are read as that slice reads them.
    """

    length = len(seq)
    start, stop, _ = slice(start, stop).indices(length)  # read as ``seq[start:stop]`` reads them
    stop = max(start, stop)
    data = _viewed(seq)
    if data is None:
        piece = seq[start:stop]
        if strand < 0:
            piece = piece.reverse_complement_rna() if rna else piece.reverse_complement()
        return piece if isinstance(piece, Seq) else Seq(piece)
    if strand > 0 and start == 0 and stop == length and isinstance(seq, Seq):
        return seq
    complement = (_rna_complement_table if rna else _dna_complement_table) if strand < 0 else None
    return Seq(RecordSequence(
        data._bases, data._forward_pieces(start, stop), reverse=data._reverse != (strand < 0),
        translate=_then(data._translate, complement),
    ))


def two_turns(seq: Seq | MutableSeq) -> Seq:
    """``seq`` read twice round (a circular record across its origin), read from ``seq``: its bytes equal ``seq + seq``."""

    data = _viewed(seq)
    if data is None:
        turns = seq + seq
        return turns if isinstance(turns, Seq) else Seq(turns)
    return Seq(RecordSequence(data._bases, data._pieces * 2, reverse=data._reverse, translate=data._translate))


def same_bases(left: Seq | MutableSeq, right: Seq | MutableSeq) -> bool:
    """Whether two sequences hold the same bases, read a span at a time."""

    return len(left) == len(right) and all(
        a == b for a, b in zip(sequence_spans(left), sequence_spans(right), strict=True))


def _sums_before(counts: NDArray[np.uint16], row: NDArray[np.int64]) -> NDArray[np.int64]:
    """Column ``k``: the counts of A, C, G and T in rows ``0..row[k] - 1``, summed between the rows asked for."""

    marks = np.unique(np.concatenate(([0], row))).astype(np.intp)  # reduceat takes intp: int32 in Pyodide (OV-3356)
    sums = np.zeros((len(marks), 4), dtype=np.int64)
    if len(marks) > 1:
        np.cumsum(np.add.reduceat(counts[: marks[-1]], marks[:-1], axis=0, dtype=np.int64), axis=0, out=sums[1:])
    return sums[np.searchsorted(marks, row)].T


def _counts_at(rows: _Rows, at: NDArray[np.int64], read: Callable[[int, int], bytes]) -> NDArray[np.int64]:
    """Column ``k``: the counts of A, C, G and T in the rows' bases before ``at[k]`` (within them): the whole rows
    before it summed, and the bases of its own row before it counted, every edge row read once (``read(row,
    size)``: the first ``size`` bases of a row)."""

    grid = at // ROW_BASES
    ended = at >= rows.end
    row = np.where(ended, len(rows.counts), grid)
    within = np.where(ended, 0, at - grid * ROW_BASES)
    result = _sums_before(rows.counts, row)
    need = np.flatnonzero(within > 0)
    if not len(need):
        return result
    need = need[np.argsort(row[need], kind="stable")]
    edge_rows, first, inverse = np.unique(row[need], return_index=True, return_inverse=True)
    lengths = np.maximum.reduceat(within[need], first)
    group = 0
    while group < len(edge_rows):
        end, size = group + 1, int(lengths[group])
        while end < len(edge_rows) and size + lengths[end] <= _EDGE_BATCH:
            size += int(lengths[end])
            end += 1
        bases = b"".join(read(int(edge), int(length))
                         for edge, length in zip(edge_rows[group:end], lengths[group:end]))
        begins = np.concatenate(([0], np.cumsum(lengths[group:end])[:-1]))
        low, high = first[group], first[end] if end < len(first) else len(need)
        taken = need[low:high]
        begin = begins[inverse[low:high] - group]
        counts = known_prefix_counts(bases, np.concatenate((begin + within[taken], begin)))
        result[:, taken] += counts[:, : len(taken)] - counts[:, len(taken):]
        group = end
    return result


def _sliced_rows(read: Callable[[int, int], bytes], stop: int) -> _Rows:
    """The count rows of the bases ``0..stop - 1`` that ``read`` gives (bases held in memory), counted now on the
    rows a record's parse pass counts."""

    counts = RowCounts()
    for low in range(0, stop, SLICE_BYTES):
        counts.add(count_rows(read(low, min(stop, low + SLICE_BYTES)), low).rows)
    return _Rows(counts.joined(), stop)


@functools.lru_cache(maxsize=8)
def _counted_order(translate: bytes | None) -> tuple[int, ...] | None:
    """The forward rows (A, C, G, T) whose counts are the counts of each base read through ``translate``, or
    ``None`` when ``translate`` does not map each counted base, and the unknown bases, as one."""

    if translate is None:
        return (0, 1, 2, 3)
    mapped: dict[int, int] = {}
    for byte in range(256):
        before, after = _count_class(byte), _count_class(translate[byte])
        if mapped.setdefault(before, after) != after:
            return None
    if sorted(mapped.get(known, -1) for known in range(4)) != [0, 1, 2, 3] or mapped.get(4, 4) != 4:
        return None
    order = [0] * 4
    for forward in range(4):
        order[mapped[forward]] = forward
    return tuple(order)


def _count_class(byte: int) -> int:
    """The row (0-3: A, C, G, T) a byte is counted in, or 4 for an unknown base."""

    counted = BASES_MATRIX[byte, KNOWN]
    return int(np.argmax(counted)) if counted.any() else 4


def known_base_counts_at(seq: Seq | MutableSeq, positions: NDArray[np.int64]) -> NDArray[np.int64]:
    """Row ``b`` (A, C, G, T), column ``k``: the count of base ``b`` in ``seq[:positions[k]]`` (clipped to it).

    The one reader of base counts (``gbdraw.core.base_counts``: either case, U as T). A record
    read from a file reads the rows its parse pass counted: whole rows summed, each edge row read once,
    so a count reads no whole sequence again. Bases held in memory (a ``Seq`` of bytes) are counted once, on the
    same rows, the first time. A view reads the rows of the sequence it views, piece by piece, on the reverse
    strand counted from the other end with the complement's bases. Bases that may change (a ``MutableSeq``) are
    counted at each call, a slice at a time, never translated or copied whole.
    """

    positions = np.clip(np.asarray(positions, dtype=np.int64), 0, len(seq))
    data = _viewed(seq)
    order = None if data is None else _counted_order(data._translate)
    if data is None or order is None:
        def read(low: int, high: int) -> bytes:
            return bytes(seq[low:high])

        return _counts_at(_sliced_rows(read, len(seq)), positions,
                          lambda row, size: read(row * ROW_BASES, min(len(seq), row * ROW_BASES + size)))
    length = len(data)
    # On the reverse strand the first p bases are the last p read forward: all of them less the first length - p.
    local = np.append(length - positions, length) if data._reverse else positions
    counts = np.zeros((4, len(local)), dtype=np.int64)
    bases = data._bases
    rows = bases.rows()
    low = 0
    for (first, last), end in zip(data._pieces, data._ends):
        inside = first + np.clip(local - low, 0, last - first)
        piece = _counts_at(rows, np.concatenate((inside, [first])), bases.row_bases)
        counts += piece[:, :-1] - piece[:, -1:]
        low = end
    if data._reverse:
        counts = counts[:, -1:] - counts[:, :-1]
    return counts[list(order)]


def known_base_counts(seq: Seq | MutableSeq) -> NDArray[np.int64]:
    """The counts of A, C, G and T in the whole of ``seq`` (:func:`known_base_counts_at`)."""

    return known_base_counts_at(seq, np.array([len(seq)], dtype=np.int64))[:, 0]


def first_base_beyond(seq: Seq | MutableSeq, alphabet: int) -> int | None:
    """A position at or before the first base of ``seq`` outside ``alphabet`` (``NUCLEOTIDE`` or ``PROTEIN``,
    ``gbdraw.core.base_counts``), or ``None`` when every base is within it.

    A record read from a file answers from what its parse pass counted, reading no base; other bases are counted
    a slice at a time.
    """

    table = _whole_table(seq)
    if table is not None:
        beyond = np.flatnonzero(np.frombuffer(bytes(table.alphabets), dtype=np.uint8) > alphabet)
        return int(table.starts[beyond[0]]) if len(beyond) else None
    column = BEYOND_NUCLEOTIDE if alphabet == NUCLEOTIDE else BEYOND_PROTEIN
    for start, span in zip(itertools.count(0, SLICE_BYTES), sequence_spans(seq, SLICE_BYTES)):
        if byte_counts(span)[column]:
            return start
    return None


def _uniform_run(buf: bytes, start: int, end: int, most: int) -> tuple[int, int, bytes]:
    """The longest run of lines of one length from ``start`` (a line start), within ``end`` and ``most`` bytes.

    Returns the length of each line with its ending, the number of lines, and the
    ending. ``buf[start:end]`` ends at a line end.
    """

    newline = buf.find(b"\n", start, end)
    carriage = buf.find(b"\r", start, newline if newline >= 0 else end)
    if carriage >= 0:
        stop = carriage + 2 if buf[carriage + 1: carriage + 2] == b"\n" else carriage + 1
    else:
        stop = newline + 1 if newline >= 0 else end
    width = stop - start
    ending = buf[start:stop].removeprefix(buf[start:stop].rstrip(b"\r\n"))
    if not ending or width > _CHUNK_BYTES:
        return width, 1, ending
    limit = min(end, start + most)
    column = buf[start + width - 1: limit: width]
    count = len(column) - len(column.lstrip(ending[-1:]))
    if count * width > limit - start:
        count = (limit - start) // width
    # Each line ends where the first does; a line end inside a line is found by the layout check
    # (``_flat_run_bases``), and a run that fails it is read line by line, whatever its lines are.
    region = buf[start: start + count * width]
    if count > 1 and any(region[width - len(ending) + j::width] != ending[j:j + 1] * count
                         for j in range(len(ending) - 1)):
        count = 1
    return width, max(count, 1), ending


def _flat_run_bases(
    raw: bytes, start: int, width: int, count: int, ending: bytes, histogram: NDArray[np.int64],
) -> int | None:
    """The number of bases of a run of GenBank sequence lines from ``start`` of ``raw`` laid out as writers do, or
    ``None``.

    The layout: up to nine (ten past position 999,999,999) digits or spaces, a
    space, then letters and spaces. Deleting digits and whitespace from such lines
    gives the bases Biopython takes from each line (everything after the
    coordinate column, spaces removed). ``histogram`` counts each byte value of the run.
    """

    content = width - len(ending)
    if content < 10:
        return None
    stop = start + width * count
    prefix = 10 if raw[start + 9:stop:width] == b" " * count else 11
    if prefix == 11 and (content < 11 or raw[start + 10:stop:width] != b" " * count):
        return None
    prefix_spaces = 0
    for column in range(prefix - 1):
        values = raw[start + column:stop:width]
        if values.translate(None, b" " + _DIGITS):
            return None
        prefix_spaces += values.count(b" ")
    if histogram[_FLAT_NOT_LETTER].any():  # a byte left after deleting digits and whitespace is not a letter
        return None
    kept = int(histogram[_FLAT_KEPT].sum())
    expected = count * (content - prefix) - (int(histogram[ord(" ")]) - prefix_spaces - count)
    if kept != expected:
        return None
    if prefix == 11:  # Biopython reads such a line one column to the left, and says so
        warnings.warn("Invalid indentation for sequence line", BiopythonParserWarning)
    return kept


def _genbank_lines(raw: bytes, encoding: str, *, warn: bool) -> bytes:
    """Biopython's GenBank sequence rule (``GenBankScanner.parse_footer``), line by line."""

    pieces: list[str] = []
    for raw_line in raw.splitlines():
        line = raw_line.decode(encoding).rstrip()
        if not line:
            if warn:
                warnings.warn("Blank line in sequence data", BiopythonParserWarning)
            continue
        if len(line) > 9 and line[9:10] != " ":
            if warn:
                warnings.warn("Invalid indentation for sequence line", BiopythonParserWarning)
            line = line[1:]
            if len(line) > 9 and line[9:10] != " ":
                raise ValueError(f"Sequence line mal-formed, '{line}'")
        pieces.append(line[10:].replace(" ", "").upper())
    return "".join(pieces).encode("ascii")


def _embl_lines(raw: bytes, encoding: str, *, warn: bool) -> bytes:
    """Biopython's EMBL sequence rule (``EmblScanner.parse_footer``), line by line."""

    pieces: list[str] = []
    for raw_line in raw.splitlines():
        line = raw_line.decode(encoding).strip()
        if not line:
            raise ValueError("Blank line in sequence data")
        split = line.rsplit(None, 1)
        if len(split) == 2 and split[1].isdigit():
            pieces.append(split[0])
        elif not line.isdigit():
            if warn:
                warnings.warn("EMBL sequence line missing coordinates", BiopythonParserWarning)
            pieces.append(line)
    return "".join(pieces).replace(" ", "").upper().encode("ascii")


def _decode_chunk(kind: int, raw: bytes, encoding: str, *, warn: bool) -> bytes:
    if kind == _FLAT_RUN:
        return raw.translate(_UPPER, _FLAT_DELETE)
    if kind == _FASTA_RUN:
        return raw.translate(None, _FASTA_DELETE)
    if kind == _GENBANK_LINES:
        return _genbank_lines(raw, encoding, warn=warn)
    if kind == _EMBL_LINES:
        return _embl_lines(raw, encoding, warn=warn)
    # Biopython joins the decoded lines and encodes them as UTF-8.
    return raw.decode(encoding).encode("utf-8").translate(None, _FASTA_DELETE)


def _index_lines(table: _SequenceTable, offset: int, raw: bytes, fmt: SequenceFormat) -> None:
    """Add the sequence lines ``raw`` (whole lines from file offset ``offset``) to ``table``, counted as they are read."""

    if fmt == "fasta" and raw.isascii():  # one run: chunks from its start, counted in one go
        table.add(offset, len(raw), _FASTA_RUN, count_rows(raw, table.length, _FASTA_DELETE))
        return
    position = 0
    encoding = table.source.encoding
    view = memoryview(raw)
    lines_until = 0  # a run whose layout failed is taken again in runs of one chunk at most, up to here
    while position < len(raw):
        counted = None
        if fmt == "embl":
            stop = _line_aligned(raw, position, position + _CHUNK_BYTES)
            kind = _EMBL_LINES
        elif fmt == "fasta":
            stop = min(len(raw), position + _CHUNK_BYTES)
            kind = _FASTA_RUN
            if not raw[position:stop].isascii():  # a character of several bytes is decoded whole, with its line
                stop = _line_aligned(raw, position, stop)
                kind = _FASTA_LINES
            else:
                counted = count_rows(view[position:stop], table.length, _FASTA_DELETE)
        else:
            most = _CHUNK_BYTES if position < lines_until else len(raw)
            width, count, ending = _uniform_run(raw, position, len(raw), most)
            stop = position + width * count
            counted = count_rows(view[position:stop], table.length, _FLAT_DELETE)
            kind = _FLAT_RUN
            if _flat_run_bases(raw, position, width, count, ending, counted.histograms.sum(axis=0)) is None:
                if position >= lines_until and stop - position > _CHUNK_BYTES:
                    lines_until = stop
                    continue
                kind, counted = _GENBANK_LINES, None
        if counted is None:
            counted = count_rows(_decode_chunk(kind, raw[position:stop], encoding, warn=True), table.length)
        table.add(offset + position, stop - position, kind, counted)
        position = stop


def _line_aligned(raw: bytes, start: int, limit: int) -> int:
    """The end of the last line that ends within ``limit`` (at least one line from ``start``)."""

    if limit >= len(raw):
        return len(raw)
    end = max(raw.rfind(b"\n", start, limit), raw.rfind(b"\r", start, limit))
    if end < 0:
        end = raw.find(b"\n", limit)
        return len(raw) if end < 0 else end + 1
    return end + 2 if raw[end:end + 2] == b"\r\n" else end + 1


# A line end followed by a line that does not start with a space.
_NEXT_UNINDENTED = re.compile(rb"(?:\r\n|\r(?!\n)|\n)(?! )")


class _Lines:
    """The lines of a file's bytes (``\\n``, ``\\r\\n`` or ``\\r`` ends one, as text mode reads them).

    ``offset`` is counted here, from the bytes read, never asked of the file.
    """

    def __init__(self, blocks: Iterator[bytes]) -> None:
        self._blocks = blocks
        self._buf = b""
        self._pos = 0
        self._base = 0
        self._eof = False
        self.ended = False
        self._carriage = -1  # the next "\r" at or after the position, or len(buf) when the buffer has none
        self._line_start = 0

    @property
    def offset(self) -> int:
        return self._base + self._pos

    def _fill(self) -> bool:
        if self._eof:
            return False
        more = next(self._blocks, b"")
        if not more:
            self._eof = True
            return False
        self._base += self._pos
        self._buf = self._buf[self._pos:] + more
        self._pos = 0
        self._carriage = -1
        return True

    def readline(self) -> bytes:
        """The next line, its ending read as ``\\n`` (as text mode reads it); ``b""`` at the end."""

        while True:
            buf, pos = self._buf, self._pos
            self._line_start = pos
            newline = buf.find(b"\n", pos)
            if self._carriage < pos:
                carriage = buf.find(b"\r", pos)
                self._carriage = len(buf) if carriage < 0 else carriage
            carriage = self._carriage
            if 0 <= newline < carriage:
                self._pos = newline + 1
                return buf[pos:self._pos]
            if carriage < len(buf) and (carriage + 1 < len(buf) or self._eof):
                self._pos = carriage + 2 if buf[carriage + 1: carriage + 2] == b"\n" else carriage + 1
                return buf[pos:carriage] + b"\n"
            if not self._fill():
                self._pos = len(buf)
                return buf[pos:]

    def indented_block(self) -> bytes:
        """The whole lines from here that start with a space, as one block (``b""`` when the next line does not)."""

        while True:
            buf, pos = self._buf, self._pos
            if buf[pos:pos + 1] != b" ":
                if pos < len(buf) or not self._fill():
                    return b""
                continue
            found = _NEXT_UNINDENTED.search(buf, pos)
            if found is not None:
                end = found.end()
            else:
                end = max(buf.rfind(b"\n"), buf.rfind(b"\r", 0, len(buf) - 1)) + 1
                if end <= pos:
                    if self._fill():
                        continue
                    end = len(buf)
            self._pos = end
            self._carriage = -1
            return buf[pos:end]

    def unread(self) -> None:
        """Step back over the line ``readline`` returned last."""

        self._pos = self._line_start
        self._carriage = -1

    def region(self, end: _RegionEnd) -> Iterator[tuple[int, bytes]]:
        """The whole lines from here to the first end line (left unread), in pieces, with their offsets.

        ``ended`` then tells whether an end line was found before the end of the file.
        """

        while True:
            buf, pos = self._buf, self._pos
            limit = len(buf) if self._eof else max(buf.rfind(b"\n"), buf.rfind(b"\r")) + 1
            if limit > pos:
                found = end.find(buf, pos, limit)
                stop = limit if found < 0 else found
                self._pos = stop
                self._carriage = -1
                if stop > pos:
                    yield self._base + pos, buf[pos:stop]
                if found >= 0:
                    self.ended = True
                    return
            if not self._fill():
                if self._pos < len(self._buf):
                    rest = self._buf[self._pos:]
                    self._pos = len(self._buf)
                    yield self._base + len(self._buf) - len(rest), rest
                self.ended = False
                return


def _line_end(buf: bytes, start: int) -> int:
    """The end of the line that starts at ``start`` (``\\n``, ``\\r\\n`` or ``\\r``), or of ``buf``."""

    newline = buf.find(b"\n", start)
    carriage = buf.find(b"\r", start, newline if newline >= 0 else len(buf))
    if carriage >= 0:
        return carriage + 2 if buf[carriage + 1: carriage + 2] == b"\n" else carriage + 1
    return newline + 1 if newline >= 0 else len(buf)


@dataclass(frozen=True)
class _RegionEnd:
    """The line that ends a sequence: one holding ``mark`` (at its start when ``anchored``) that ``accepts`` takes."""

    marks: tuple[bytes, ...]
    accepts: Callable[[bytes], bool]
    anchored: bool = True

    def find(self, buf: bytes, pos: int, limit: int) -> int:
        """The start of the first end line in ``buf[pos:limit]`` (whole lines from a line start), else -1."""

        best = -1
        for mark in self.marks:
            index = buf.find(mark, pos, limit)
            while index >= 0 and (best < 0 or index < best):
                if self.anchored:
                    start = index if index == pos or buf[index - 1] in b"\r\n" else -1
                else:
                    start = max(pos, buf.rfind(b"\n", pos, index) + 1, buf.rfind(b"\r", pos, index) + 1)
                if start >= 0 and self.accepts(buf[start:_line_end(buf, start)]):
                    best = start
                    break
                index = buf.find(mark, index + 1, limit)
        return best


@dataclass
class _Slot:
    """One record the pass met, and the table of its sequence when it has one."""

    table: _SequenceTable | None = None


class _Scan:
    """One pass over a sequence file: the lines for the parser, and each record's sequence table."""

    def __init__(self, source: _SequenceFile, fmt: SequenceFormat, blocks: Iterator[bytes]) -> None:
        self.source = source
        self.fmt = fmt
        self.slots: list[_Slot] = []
        self._lines = _Lines(blocks)
        self._first = True

    def _decode(self, line: bytes) -> str:
        text = line.decode(self.source.encoding)
        if self._first:
            self._first = False
            text = text.removeprefix(_BYTE_ORDER_MARK)
        return text

    def _table(self) -> _SequenceTable:
        table = _SequenceTable(self.source, key=self.source.table_key())
        self.slots[-1].table = table
        return table

    def flat_lines(self) -> Iterator[str]:
        """A flat file's lines without its sequence lines (Biopython reads an empty sequence): a line at a time, or
        a block of indented lines (header continuations, a feature table) as one piece."""

        genbank = self.fmt == "genbank"
        width = 12 if genbank else 5
        record_start = _GENBANK_RECORD_START if genbank else _EMBL_RECORD_START
        headers = _GENBANK_SEQUENCE_HEADERS if genbank else _EMBL_SEQUENCE_HEADERS
        footer = False
        while True:
            if not footer and not self._first:  # header continuations and the feature table, as blocks
                block = self._lines.indented_block()
                if block:
                    yield _text_mode_endings(block.decode(self.source.encoding))
                    continue
            line = self._lines.readline()
            if not line:
                return
            text = self._decode(line)
            if text[:width] == record_start:
                self.slots.append(_Slot())
                footer = False
            elif not self.slots:  # Biopython skips the lines before the first record
                pass
            elif text[:width].rstrip() in headers:
                footer = True
            elif footer and not (genbank and (text[:width] == " " * width or text[:3] == "WGS")):
                footer = False
                self._lines.unread()
                if genbank:
                    if not self._index_region(_RegionEnd((b"//", b"CONTIG"), self._genbank_end)):
                        # Biopython ends a record cut off in its sequence as if "//" followed.
                        warnings.warn("Premature end of file in sequence data", BiopythonParserWarning)
                        yield "//\n"
                    elif self._lines.readline().startswith(b"CONTIG"):
                        # Biopython stops at a CONTIG line after sequence lines and refuses the record there.
                        raise ValueError("A CONTIG line follows the sequence lines of a record")
                    else:
                        self._lines.unread()
                    continue
                _check_embl_first_line(text)
                if not self._index_region(_RegionEnd((b"//",), _embl_end, anchored=False)):
                    raise ValueError("Premature end of file in sequence data")
                continue
            yield text

    def _genbank_end(self, line: bytes) -> bool:
        return line.startswith(b"CONTIG") or line.decode(self.source.encoding).rstrip() == "//"

    def _index_region(self, end: _RegionEnd) -> bool:
        """Index the sequence lines from here to ``end`` into a new table of the current record."""

        table = self._table()
        for offset, raw in self._lines.region(end):
            _index_lines(table, offset, raw, self.fmt)
        return self._lines.ended

    def fasta_lines(self) -> Iterator[str]:
        """A FASTA file's header lines from the first one (Biopython reads each entry empty)."""

        while line := self._lines.readline():
            text = self._decode(line)
            if text.startswith(">"):
                self.slots.append(_Slot())
                yield text
                self._index_region(_FASTA_ENTRY)

    def gff3_lines(self) -> Iterator[str]:
        """A GFF3 file's lines, its ``##FASTA`` part as header lines only (BCBio reads each entry empty)."""

        while line := self._lines.readline():
            text = self._decode(line)
            yield text
            if text.startswith("##FASTA"):
                break
        while line := self._lines.readline():
            text = self._decode(line)
            if text.startswith(">"):
                yield text
            for _offset, raw in self._lines.region(_FASTA_ENTRY):
                raw.decode(self.source.encoding)  # text that is not in the encoding stops the read, as before


def _text_mode_endings(text: str) -> str:
    """Whole lines with their endings read as ``\\n``, as text mode reads them."""

    return text.replace("\r\n", "\n").replace("\r", "\n") if "\r" in text else text


def _starts_entry(line: bytes) -> bool:
    return line.startswith(b">")


_FASTA_ENTRY = _RegionEnd((b">",), _starts_entry)


def _embl_end(line: bytes) -> bool:
    return line.strip() == b"//"


def _check_embl_first_line(text: str) -> None:
    line = text.rstrip()
    if not (line[:5] == " " * 5 or line.strip() == "//"):
        raise ValueError(f"Unexpected content after SQ or CO line: {line!r}")


def _with_encoding(source: _SequenceFile, read: Callable[[Iterator[bytes]], _Read]) -> _Read:
    """``read`` of one pass over the file's blocks, in the encoding of every sequence file (:func:`in_sequence_encoding`).

    The pass is closed when ``read`` returns or raises, so a gzip expansion ahead of it stops then, not when the
    blocks are collected (an error's traceback can hold them).
    """

    def read_in(encoding: str) -> _Read:
        source.encoding = encoding
        with contextlib.closing(source.blocks()) as blocks:
            return read(blocks)

    return in_sequence_encoding(read_in)


def _attach(records: list[SeqRecord], slots: list[_Slot]) -> list[SeqRecord]:
    """Give each record the sequence its slot found, in file order."""

    if len(records) != len(slots):
        raise ValueError(f"{len(slots)} record starts were found, but the parser read {len(records)} records")
    for record, slot in zip(records, slots):
        table = slot.table
        if table is None or not table.length:
            continue  # Biopython's own sequence: undefined of the stated length, or empty
        expected = len(record.seq)
        if expected and expected != table.length:
            warnings.warn(f"Expected sequence length {expected}, found {table.length} ({record.id}).",
                          BiopythonParserWarning)
        table.finish()
        record.seq = Seq(RecordSequence(table))
    return records


def read_flat_records(
    path: str | PathLike[str], parse: Callable[[Iterator[str], FlatFileFormat], list[SeqRecord]],
    fmt_of: Callable[[str], FlatFileFormat],
) -> list[SeqRecord]:
    """The records ``parse`` builds from a flat file's lines, each with its sequence read from the file.

    ``fmt_of`` names the format from the first non-blank line, read from the
    blocks the parse pass reads (the file is read once).
    """

    source = _SequenceFile(path)

    def read(pass_blocks: Iterator[bytes]) -> list[SeqRecord]:
        first, blocks = _first_text_line(source, pass_blocks)
        fmt = fmt_of(first)
        scan = _Scan(source, fmt, blocks)
        return _attach(parse(scan.flat_lines(), fmt), scan.slots)

    return _with_encoding(source, read)


def _first_text_line(source: _SequenceFile, blocks: Iterator[bytes]) -> tuple[str, Iterator[bytes]]:
    """The first line that is not blank, and the file's blocks from its start again.

    The blocks read to find the line (usually one) are kept until the pass reads
    them again, so the file is read once.
    """

    head: list[bytes] = []
    scan = _Scan(source, "genbank", _recorded(blocks, head))
    lines = (scan._decode(line) for line in iter(scan._lines.readline, b""))
    return next((text for text in lines if text.strip()), ""), itertools.chain(head, blocks)


def _recorded(blocks: Iterator[bytes], into: list[bytes]) -> Iterator[bytes]:
    for block in blocks:
        into.append(block)
        yield block


def read_fasta_records(
    path: str | PathLike[str], parse: Callable[[Iterator[str]], list[SeqRecord]],
) -> list[SeqRecord] | None:
    """The entries ``parse`` builds from a FASTA file's header lines, each with its sequence, or ``None`` without one."""

    source = _SequenceFile(path)

    def read(blocks: Iterator[bytes]) -> list[SeqRecord] | None:
        scan = _Scan(source, "fasta", blocks)
        lines = scan.fasta_lines()
        first = next(lines, None)
        if first is None:
            return None
        records = parse(_prepend(first, lines))
        return _attach(records, scan.slots)

    return _with_encoding(source, read)


def _prepend(first: str, rest: Iterator[str]) -> Iterator[str]:
    yield first
    yield from rest


def _whole_table(seq: object) -> _SequenceTable | None:
    """The sequence table of a record in an input file when ``seq`` reads the whole of it, forward, else ``None``."""

    data = getattr(seq, "_data", None)
    if (isinstance(data, RecordSequence) and isinstance(data._bases, _SequenceTable) and not data._reverse
            and data._translate is None and data._pieces == ((0, data._bases.length),)):
        return data._bases
    return None


def _file_table(seq: object) -> _SequenceTable | None:
    """The sequence table in an input file that ``seq`` reads, or ``None`` when it holds its bases in memory."""

    data = getattr(seq, "_data", None)
    if isinstance(data, RecordSequence) and isinstance(data._bases, _SequenceTable):
        return data._bases
    return None


def detached_records(records: Iterable[SeqRecord]) -> list[SeqRecord]:
    """The records, each holding its sequence in memory instead of reading it from its input file.

    Records handed to a caller (the public readers) are detached: the caller may
    rewrite or delete the file, which a record sequence would then refuse to read.
    A record that reads its file is copied (a shallow copy, with its own
    sequence), so records a cache owns stay as they are. A whole record keeps the
    count rows its parse pass counted, so its bases are not counted again.
    """

    detached: list[SeqRecord] = []
    for record in records:
        if _file_table(record.seq) is not None:
            table = _whole_table(record.seq)
            bases = bytes(record.seq)
            record = copy.copy(record)
            record.seq = Seq(bases)
            if table is not None:
                _held(record.seq, bases, table.rows())
        detached.append(record)
    return detached


def sequence_files_unchanged(records: Iterable[SeqRecord]) -> bool:
    """Whether every file the records' sequences read still holds the bytes it held when they were read.

    One stat per distinct file and no read: a holder of records kept across
    requests (the prepared input cache) asks before it reuses them, and builds
    them again from the file when a file was removed or written again.
    Records that hold their sequence in memory read no file.
    """

    checked: set[int] = set()
    for record in records:
        table = _file_table(record.seq)
        if table is None:
            continue
        source = table.source
        if source.key in checked:
            continue
        if not source.unchanged():
            return False
        checked.add(source.key)
    return True


def read_gff3_text(path: str | PathLike[str], read: Callable[[Iterator[str]], _Read]) -> _Read:
    """What ``read`` makes of a GFF3 file's lines; its ``##FASTA`` part keeps only its header lines.

    The FASTA partner of the pair (the same file, or another) gives the sequences.
    """

    source = _SequenceFile(path)
    return _with_encoding(source, lambda blocks: read(_Scan(source, "fasta", blocks).gff3_lines()))


__all__ = [
    "RecordSequence",
    "SequenceFormat",
    "detached_records",
    "first_base_beyond",
    "known_base_counts",
    "known_base_counts_at",
    "read_fasta_records",
    "read_flat_records",
    "read_gff3_text",
    "same_bases",
    "sequence_files_unchanged",
    "sequence_view",
    "two_turns",
]
