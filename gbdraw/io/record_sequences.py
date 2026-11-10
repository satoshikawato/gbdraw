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

A gzip input is expanded into a temporary file as the pass reads it, because gzip
cannot be read at an offset; a pass that stops early expands no more, and the
expansion keeps the cap of :mod:`gbdraw.io.input_text`. The pass
counts bytes itself and never asks the file for its position, and every read below
a record sequence goes through ``_SequenceFile.read_at``, the one reader, so the
Web can put its own range reader there.

The sequence lines are checked in runs of lines of one length. A run laid out as
writers lay out sequence lines (coordinates, then letters and spaces) is read
with one ``bytes.translate``; any other line is read with the rule Biopython
applies to each line, so both give Biopython's bytes. A cache budget decides how
many decoded bases stay in memory, never what a read returns.

The same pass counts the bases (``gbdraw.core.base_counts``): the one producer of
a record's base counts. The letter check of a FASTA entry, the GC percent of a
record or a crop and the GC and skew windows read those counts
(:func:`known_base_counts_at`, :func:`first_base_beyond`), reading only the bases
at the edges of a span again.
"""

from __future__ import annotations

import bisect
import copy
import functools
import itertools
import os
import re
import tempfile
import threading
import warnings
import weakref
from array import array
from collections import OrderedDict, deque
from collections.abc import Callable, Iterable, Iterator
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
    ROW_BYTES,
    SLICE_BYTES,
    byte_histograms,
    count_matrix,
    known_prefix_counts,
    row_counts,
)
from ..core.sequence_spans import sequence_spans
from ..exceptions import InputFileError
from . import input_text
from .input_text import FlatFileFormat, in_sequence_encoding, is_gzip_file, open_expanded_input, open_input_bytes

SequenceFormat = Literal["genbank", "embl", "fasta"]

_BLOCK_BYTES = 8 << 20  # sequential reads of the parse pass
_CHUNK_BYTES = 64 << 10  # most raw bytes of one indexed chunk (whole lines)
_CACHE_BYTES = 256 << 20  # decoded bases the block cache keeps, all input files together
_OPEN_FILES = 16  # input files kept open for reads at offsets, all records together
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
# A chunk of a run kind is counted from its raw bytes, in rows of ``ROW_BYTES`` (any split of a run decodes alone);
# a chunk read line by line is counted from its bases, as one row.
_RUN_MATRICES = {_FLAT_RUN: count_matrix(_FLAT_DELETE), _FASTA_RUN: count_matrix(_FASTA_DELETE)}
assert _CHUNK_BYTES % SLICE_BYTES == 0 and SLICE_BYTES % ROW_BYTES == 0
_FLAT_KEPT = np.ones(256, dtype=bool)
_FLAT_KEPT[list(_FLAT_DELETE)] = False
_FLAT_NOT_LETTER = _FLAT_KEPT.copy()
_FLAT_NOT_LETTER[list(_LETTERS)] = False
_EDGE_BATCH = 1 << 20  # most bases of edge rows counted together


class _Rows(NamedTuple):
    """Count rows of a sequence's bases: ``starts`` (each row's first base, then the end), ``counts`` (rows x 4: A,
    C, G and T, U as T, either case), and ``read``, which gives the bases of one row."""

    starts: NDArray[np.int64]
    counts: NDArray[np.int32] | NDArray[np.int64]
    read: Callable[[int], bytes]


def _alphabet(rows: NDArray[np.int64]) -> int:
    """What the bases counted in ``rows`` hold (``gbdraw.core.base_counts``)."""

    beyond = rows[:, (BEYOND_NUCLEOTIDE, BEYOND_PROTEIN)].any(axis=0)
    return OTHER if beyond[1] else PROTEIN if beyond[0] else NUCLEOTIDE

# Biopython's record start and sequence headers (Bio.GenBank.Scanner).
_GENBANK_RECORD_START = "LOCUS       "
_GENBANK_SEQUENCE_HEADERS = ("CONTIG", "ORIGIN", "BASE COUNT", "WGS", "TSA", "TLS")
_EMBL_RECORD_START = "ID   "
_EMBL_SEQUENCE_HEADERS = ("SQ", "CO")

_Read = TypeVar("_Read")


class _BlockCache:
    """The decoded bases read so far, least recently used dropped first past ``_CACHE_BYTES``.

    One cache serves every input file, so the budget bounds the memory of all of
    them. A block is a whole record sequence that fits the budget, else one chunk.
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


class _OpenFiles:
    """The input files open for reads at offsets: at most ``_OPEN_FILES``, least recently used closed first.

    A file is opened by its path when a read needs it, so the number of inputs is
    not bounded by the number of files a process may hold open. Every read checks
    that the file still holds the bytes its records were indexed from. The lock
    guards only the set of handles: a handle is taken out of it while it is read,
    and opening, reading and closing happen outside the lock.
    """

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._handles: OrderedDict[int, IO[bytes]] = OrderedDict()

    def read(self, source: _SequenceFile, offset: int, size: int) -> bytes:
        _release_pending()
        with self._lock:
            handle = self._handles.pop(source.key, None)
        try:
            if handle is None:
                handle = open_input_bytes(source.path)
        except OSError as error:
            raise source.changed() from error
        try:
            if _stamp(handle) != source.stamp:
                raise source.changed()
            handle.seek(offset)
            data = handle.read(size)
        except InputFileError:
            handle.close()
            raise
        except OSError as error:
            handle.close()
            raise source.changed() from error
        closing: list[IO[bytes]] = []
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
_EXPANSIONS: list[tempfile.TemporaryDirectory[str]] = []
_EXPANSIONS_LOCK = threading.Lock()


def _expansion_directory() -> str:
    """The one directory of this process's expanded gzip inputs, removed when the process ends."""

    with _EXPANSIONS_LOCK:
        if not _EXPANSIONS:
            _EXPANSIONS.append(tempfile.TemporaryDirectory(prefix="gbdraw-inputs-"))
        return _EXPANSIONS[0].name


# Inputs whose last record sequence has gone: (file key, gzip copy). A finalizer runs at any
# allocation, also inside the readers' locked sections, so it only queues the input here
# (deque.append takes no lock); the readers release the queue before they take their locks.
_RELEASED: deque[tuple[int, str | None]] = deque()


def _release(file_key: int, copy: str | None) -> None:
    _RELEASED.append((file_key, copy))


def _release_pending() -> None:
    """Close the files of the queued inputs, forget their blocks and remove their gzip copies.

    Called outside every lock of this module.
    """

    while True:
        try:
            file_key, copy = _RELEASED.popleft()
        except IndexError:
            return
        _OPEN.close(file_key)
        _CACHE.drop_file(file_key)
        if copy is not None:
            try:
                os.unlink(copy)
            except OSError:
                pass  # removed with its directory already


class _SequenceFile:
    """One input file's bytes: read in order by the parse pass, and at offsets by its record sequences.

    A gzip input is expanded into a file of its own (gzip cannot be read at an
    offset) as the parse pass reads it, removed when the last record sequence of
    the input goes. A link
    is followed once: the records read the file it names, so they keep reading
    after the link goes (the Web Worker links each input into the folder of one
    request and removes the folder when the request ends).
    ``read_at`` is the one read below a record sequence. The parse pass records
    the file's stamp, and a read of a file that changed since is refused with an
    ``InputFileError``, never answered with other bases.
    """

    def __init__(self, path: str | PathLike[str]) -> None:
        _release_pending()
        self.name = fspath(path)
        self.encoding = "utf-8"
        self.key = next(_FILE_KEYS)
        self.stamp: _Stamp = (-1, -1, -1)
        copy = None
        if is_gzip_file(self.name):
            handle, copy = tempfile.mkstemp(dir=_expansion_directory(), prefix="input-")
            self.stamp = _stamp_of(os.fstat(handle))
            os.close(handle)
        self.path = copy or os.path.realpath(self.name)
        # Bytes of a gzip input written to its copy so far, and whether that is all of them (None: a plain file).
        self.expanded: int | None = None if copy is None else 0
        self.expanded_all = copy is None
        self.original: _Stamp | None = None  # a gzip input's own stamp when its first pass expanded it
        weakref.finalize(self, _release, self.key, copy)
        self._tables = 0

    def table_key(self) -> int:
        """A key no other sequence table of this file has (its blocks' cache keys hold it)."""

        self._tables += 1
        return self._tables

    def blocks(self) -> Iterator[bytes]:
        """The file's bytes from its start, in order, from a handle of the parse pass's own.

        A gzip input's copy is read as far as it is written, then expanded on: each
        block is in the copy before the pass sees it.
        """

        with open_input_bytes(self.path) as handle:
            if self.expanded is None:
                self.stamp = _stamp(handle)
            while block := handle.read(_BLOCK_BYTES):
                yield block
            if _stamp(handle) != self.stamp:
                raise self.changed()
        if not self.expanded_all:
            yield from self._expand()

    def _expand(self) -> Iterator[bytes]:
        """A gzip input's bytes past its copy, each written to the copy before it is handed on.

        A later pass (the encoding retry) expands the input again and passes over what the copy holds, so the
        input must be the file the first pass expanded. A read stops at the cap (its last block is cut there), so
        the cap refuses only a byte past it that a reader asks for. A pass refused while it expands leaves no byte
        of its copy behind.
        """

        assert self.expanded is not None
        try:
            original = _stamp_of(os.stat(self.name))
            if self.original is not None and original != self.original:
                raise self.changed()
            self.original = original
            with open_expanded_input(self.name) as expanded, open(self.path, "ab") as copy:
                skip = self.expanded  # in the copy already: expanded again (the cap counts them) and passed over
                while skip:
                    if not (passed := len(expanded.read(min(skip, _BLOCK_BYTES)))):
                        raise self.changed()
                    skip -= passed
                while block := expanded.read(max(1, min(_BLOCK_BYTES, input_text.MAX_EXPANDED_INPUT_BYTES - self.expanded))):
                    copy.write(block)
                    copy.flush()
                    self.expanded += len(block)
                    self.stamp = _stamp(copy)
                    yield block
        except GeneratorExit:
            raise
        except BaseException:
            with open(self.path, "wb"):
                pass  # the copy of a refused pass: emptied, so a later pass expands the input again
            self.expanded = 0
            self.stamp = _stamp_of(os.stat(self.path))
            raise
        self.expanded_all = True

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
    record's base counts): chunk ``i`` holds the count rows from ``first_rows[i]``,
    each row the counts of A, C, G and T in at most ``ROW_BYTES`` raw bytes of a
    run, or in a whole chunk read line by line. A count of any span is then whole
    rows summed and its two edge rows counted (:func:`known_base_counts_at`).
    """

    source: _SequenceFile
    key: int
    offsets: array[int] = field(default_factory=lambda: array("q"))
    sizes: array[int] = field(default_factory=lambda: array("q"))
    starts: array[int] = field(default_factory=lambda: array("q"))
    kinds: bytearray = field(default_factory=bytearray)
    alphabets: bytearray = field(default_factory=bytearray)
    first_rows: array[int] = field(default_factory=lambda: array("q"))
    length: int = 0
    _row_starts: list[NDArray[np.int64]] = field(default_factory=list)
    _row_counts: list[NDArray[np.int32]] = field(default_factory=list)
    _rows: _Rows | None = None
    _row_total: int = 0

    def add(self, offset: int, raw: bytes, kind: int, rows: NDArray[np.int64]) -> None:
        """Add the raw bytes ``raw`` from file offset ``offset``, of one kind, with the counts of their rows.

        A run is cut into chunks of ``_CHUNK_BYTES`` raw bytes, rows of ``ROW_BYTES``; lines read line by line
        join the chunk before them up to the usual size, their counts its one row.
        """

        last = len(self.kinds) - 1
        if (kind not in _RUN_MATRICES and last >= 0 and self.kinds[last] == kind
                and self.offsets[last] + self.sizes[last] == offset and self.sizes[last] + len(raw) <= _CHUNK_BYTES):
            self.sizes[last] += len(raw)  # lines read line by line stay in chunks of the usual size
            self.alphabets[last] = max(self.alphabets[last], _alphabet(rows))
            self._row_counts[-1][-1] += rows[0, KNOWN].astype(np.int32)
            self.length += int(rows[0, BASES])
            return
        per_chunk = _CHUNK_BYTES // ROW_BYTES if kind in _RUN_MATRICES else len(rows)
        bases = rows[:, BASES]
        row_starts = self.length + np.concatenate(([0], np.cumsum(bases)))
        for first in range(0, len(rows), per_chunk):
            start = first * ROW_BYTES
            self.offsets.append(offset + start)
            self.sizes.append(len(raw) - start if first + per_chunk >= len(rows) else per_chunk * ROW_BYTES)
            self.starts.append(int(row_starts[first]))
            self.kinds.append(kind)
            self.alphabets.append(_alphabet(rows[first:first + per_chunk]))
            self.first_rows.append(self._row_total + first)
        self._row_starts.append(row_starts[:-1])
        self._row_counts.append(rows[:, KNOWN].astype(np.int32))
        self._row_total += len(rows)
        self.length = int(row_starts[-1])

    def rows(self) -> _Rows:
        """The count rows of the whole sequence (joined once, when the parse pass has added them all)."""

        if self._rows is None:
            starts = np.concatenate([*self._row_starts, [self.length]]).astype(np.int64)
            counts = np.concatenate(self._row_counts) if self._row_counts else np.zeros((0, 4), dtype=np.int32)
            self._rows = _Rows(starts, counts, self._row_bases)
            self._row_starts, self._row_counts = [], []
        return self._rows

    def _row_bases(self, row: int) -> bytes:
        """The bases of count row ``row``: read from the file and decoded alone."""

        index = bisect.bisect_right(self.first_rows, row) - 1
        kind = self.kinds[index]
        offset, size = self.offsets[index], self.sizes[index]
        if kind in _RUN_MATRICES:
            skip = (row - self.first_rows[index]) * ROW_BYTES
            offset, size = offset + skip, min(ROW_BYTES, size - skip)
        return _decode_chunk(kind, self.source.read_at(offset, size), self.source.encoding, warn=False)

    def _block(self, index: int) -> bytes:
        """Chunk ``index``'s bases (``-1``: the whole sequence), through the block cache."""

        if index < 0:
            return _CACHE.get((self.source.key, self.key, -1), lambda: self._bulk(0, self.length))
        return _CACHE.get((self.source.key, self.key, index), lambda: self._decode(index, index + 1))

    def _decode(self, first: int, stop: int) -> bytes:
        """The bases of chunks ``first``..``stop - 1``, which are of one kind and adjacent in the file."""

        start = self.offsets[first]
        raw = self.source.read_at(start, self.offsets[stop - 1] + self.sizes[stop - 1] - start)
        return _decode_chunk(self.kinds[first], raw, self.source.encoding, warn=False)

    def read(self, start: int, stop: int) -> bytes:
        """Bases ``start``..``stop - 1`` (``0 <= start <= stop <= length``).

        A sequence that fits the cache budget is one block: read whole once, then
        every read is a slice of it (the whole read the block itself). A longer one
        is read a chunk at a time through the cache, or past half the budget in
        bulk, around the cache.
        """

        if start >= stop:
            return b""
        if self.length <= _CACHE_BYTES:
            whole = self._block(-1)
            return whole if stop - start == self.length else whole[start:stop]
        if stop - start > _CACHE_BYTES // 2:
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


# What a record sequence reads: the sequence table of a record in its input file, or bases held in memory.
_Bases = _SequenceTable | bytes


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
        self._pieces = ((0, len(bases) if isinstance(bases, bytes) else bases.length),) if pieces is None else pieces
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
        reads = [bases[first:last] if isinstance(bases, bytes) else bases.read(first, last)
                 for first, last in self._forward_pieces(start, stop)]
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

    A record sequence is itself; bases held as ``bytes`` are read in place.
    Bases that may change (a ``MutableSeq``) and a sequence that is not defined
    everywhere are left to Biopython's own slices and reverse complements.
    """

    data = getattr(seq, "_data", None)
    if isinstance(data, RecordSequence):
        return data
    if isinstance(seq, Seq) and type(data) is bytes:
        return RecordSequence(data)
    return None


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


def _counts_at(rows: _Rows, at: NDArray[np.int64]) -> NDArray[np.int64]:
    """Column ``k``: the counts of A, C, G and T in the rows' bases before ``at[k]`` (within them): the whole rows
    before it summed, and the bases of its own row before it counted, every edge row read once."""

    cumulative = np.zeros((len(rows.counts) + 1, 4), dtype=np.int64)
    np.cumsum(rows.counts, axis=0, out=cumulative[1:])
    row = np.searchsorted(rows.starts, at, side="right") - 1
    row = np.minimum(row, len(rows.counts))
    result = cumulative[row].T.copy()
    within = at - rows.starts[row]
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
        bases = b"".join(rows.read(int(edge))[: int(length)] for edge, length in zip(edge_rows[group:end], lengths[group:end]))
        begins = np.concatenate(([0], np.cumsum(lengths[group:end])[:-1]))
        low, high = first[group], first[end] if end < len(first) else len(need)
        taken = need[low:high]
        begin = begins[inverse[low:high] - group]
        counts = known_prefix_counts(bases, np.concatenate((begin + within[taken], begin)))
        result[:, taken] += counts[:, : len(taken)] - counts[:, len(taken):]
        group = end
    return result


def _sliced_rows(read: Callable[[int, int], bytes], start: int, stop: int) -> _Rows:
    """The count rows of bases ``start``..``stop - 1`` that ``read`` gives (bases held in memory), counted now."""

    counts = [row_counts(read(low, min(stop, low + SLICE_BYTES)))[:, KNOWN] for low in range(start, stop, SLICE_BYTES)]
    starts = np.append(np.arange(start, stop, ROW_BYTES, dtype=np.int64), stop)
    return _Rows(
        starts, np.concatenate(counts) if counts else np.zeros((0, 4), dtype=np.int64),
        lambda row: read(int(starts[row]), int(starts[row + 1])),
    )


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

    counted = row_counts(bytes((byte,)))[0, KNOWN]
    return int(np.argmax(counted)) if counted.any() else 4


def known_base_counts_at(seq: Seq | MutableSeq | bytes, positions: NDArray[np.int64]) -> NDArray[np.int64]:
    """Row ``b`` (A, C, G, T), column ``k``: the count of base ``b`` in ``seq[:positions[k]]`` (clipped to it).

    The one reader of base counts (``gbdraw.core.base_counts``: either case, U as T). A record
    read from a file reads the rows its parse pass counted: whole rows summed, each edge row read once,
    so a count reads no whole sequence again. A view reads the rows of the sequence it views, piece by piece, on
    the reverse strand counted from the other end with the complement's bases. Bases held in memory are counted
    here, a slice at a time, never translated or copied whole.
    """

    positions = np.clip(np.asarray(positions, dtype=np.int64), 0, len(seq))
    data = _viewed(seq) if not isinstance(seq, bytes) else RecordSequence(seq)
    order = None if data is None else _counted_order(data._translate)
    if data is None or order is None:
        return _counts_at(_sliced_rows(lambda low, high: bytes(seq[low:high]), 0, len(seq)), positions)
    length = len(data)
    # On the reverse strand the first p bases are the last p read forward: all of them less the first length - p.
    local = np.append(length - positions, length) if data._reverse else positions
    counts = np.zeros((4, len(local)), dtype=np.int64)
    bases = data._bases
    table_rows = bases.rows() if isinstance(bases, _SequenceTable) else None
    sliced: dict[tuple[int, int], _Rows] = {}
    low = 0
    for (first, last), end in zip(data._pieces, data._ends):
        inside = first + np.clip(local - low, 0, last - first)
        if table_rows is not None:
            piece = _counts_at(table_rows, np.concatenate((inside, [first])))
            counts += piece[:, :-1] - piece[:, -1:]
        else:
            assert isinstance(bases, bytes)
            rows = sliced.get((first, last))
            if rows is None:
                rows = sliced[(first, last)] = _sliced_rows(lambda a, b: bases[a:b], first, last)
            counts += _counts_at(rows, inside)
        low = end
    if data._reverse:
        counts = counts[:, -1:] - counts[:, :-1]
    return counts[list(order)]


def known_base_counts(seq: Seq | MutableSeq | bytes) -> NDArray[np.int64]:
    """The counts of A, C, G and T in the whole of ``seq`` (:func:`known_base_counts_at`)."""

    return known_base_counts_at(seq, np.array([len(seq)], dtype=np.int64))[:, 0]


def first_base_beyond(seq: Seq | MutableSeq, alphabet: int) -> int | None:
    """A position at or before the first base of ``seq`` outside ``alphabet`` (``NUCLEOTIDE`` or ``PROTEIN``,
    ``gbdraw.core.base_counts``), or ``None`` when every base is within it.

    A record read from a file answers from what its parse pass counted, reading no base; other bases are counted
    a slice at a time.
    """

    data = _viewed(seq)
    if (data is not None and isinstance(data._bases, _SequenceTable) and not data._reverse
            and data._translate is None and data._pieces == ((0, data._bases.length),)):
        table = data._bases
        beyond = np.flatnonzero(np.frombuffer(bytes(table.alphabets), dtype=np.uint8) > alphabet)
        return int(table.starts[beyond[0]]) if len(beyond) else None
    column = BEYOND_NUCLEOTIDE if alphabet == NUCLEOTIDE else BEYOND_PROTEIN
    for start, span in zip(itertools.count(0, SLICE_BYTES), sequence_spans(seq, SLICE_BYTES)):
        beyond = np.flatnonzero(row_counts(span)[:, column])
        if len(beyond):
            return start + int(beyond[0]) * ROW_BYTES
    return None


def _uniform_run(buf: bytes, start: int, end: int) -> tuple[int, int, bytes]:
    """The longest run of lines of one length from ``start`` (a line start), within ``end``.

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
    limit = min(end, start + _CHUNK_BYTES)
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


def _flat_run_bases(raw: bytes, width: int, count: int, ending: bytes, histogram: NDArray[np.int64]) -> int | None:
    """The number of bases of a run of GenBank sequence lines laid out as writers do, or ``None``.

    The layout: up to nine (ten past position 999,999,999) digits or spaces, a
    space, then letters and spaces. Deleting digits and whitespace from such lines
    gives the bases Biopython takes from each line (everything after the
    coordinate column, spaces removed). ``histogram`` counts each byte value of ``raw``.
    """

    content = width - len(ending)
    if content < 10:
        return None
    prefix = 10 if raw[9::width] == b" " * count else 11
    if prefix == 11 and (content < 11 or raw[10::width] != b" " * count):
        return None
    prefix_spaces = 0
    for column in range(prefix - 1):
        values = raw[column::width]
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

    if fmt == "fasta" and raw.isascii():  # one run: chunks and rows from its start, counted in one go
        table.add(offset, raw, _FASTA_RUN, row_counts(raw, _RUN_MATRICES[_FASTA_RUN]))
        return
    position = 0
    encoding = table.source.encoding
    while position < len(raw):
        rows = None
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
                rows = row_counts(raw[position:stop], _RUN_MATRICES[kind])
        else:
            width, count, ending = _uniform_run(raw, position, len(raw))
            stop = position + width * count
            histograms = byte_histograms(raw[position:stop])
            kind = _GENBANK_LINES
            if _flat_run_bases(raw[position:stop], width, count, ending, histograms.sum(axis=0)) is not None:
                kind, rows = _FLAT_RUN, histograms @ _RUN_MATRICES[_FLAT_RUN]
        piece = raw[position:stop]
        if rows is None:
            rows = row_counts(_decode_chunk(kind, piece, encoding, warn=True)).sum(axis=0, keepdims=True)
        table.add(offset + position, piece, kind, rows)
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

    def __init__(self, source: _SequenceFile, fmt: SequenceFormat, blocks: Iterator[bytes] | None = None) -> None:
        self.source = source
        self.fmt = fmt
        self.slots: list[_Slot] = []
        self._lines = _Lines(source.blocks() if blocks is None else blocks)
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
        """A flat file's lines without its sequence lines (Biopython reads an empty sequence)."""

        genbank = self.fmt == "genbank"
        width = 12 if genbank else 5
        record_start = _GENBANK_RECORD_START if genbank else _EMBL_RECORD_START
        headers = _GENBANK_SEQUENCE_HEADERS if genbank else _EMBL_SEQUENCE_HEADERS
        footer = False
        while True:
            if not footer and not self._first:  # header continuations and the feature table, as blocks
                block = self._lines.indented_block()
                if block:
                    yield from _block_lines(block.decode(self.source.encoding))
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


def _block_lines(text: str) -> Iterator[str]:
    """The lines of a block of whole lines, each ending in ``\\n`` as text mode reads them."""

    if "\r" in text:
        text = text.replace("\r\n", "\n").replace("\r", "\n")
    lines = text.split("\n")
    lines.pop()
    for line in lines:
        yield line + "\n"


def _starts_entry(line: bytes) -> bool:
    return line.startswith(b">")


_FASTA_ENTRY = _RegionEnd((b">",), _starts_entry)


def _embl_end(line: bytes) -> bool:
    return line.strip() == b"//"


def _check_embl_first_line(text: str) -> None:
    line = text.rstrip()
    if not (line[:5] == " " * 5 or line.strip() == "//"):
        raise ValueError(f"Unexpected content after SQ or CO line: {line!r}")


def _with_encoding(source: _SequenceFile, read: Callable[[], _Read]) -> _Read:
    """``read()`` with the file read in the encoding of every sequence file (:func:`in_sequence_encoding`)."""

    def read_in(encoding: str) -> _Read:
        source.encoding = encoding
        return read()

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

    def read() -> list[SeqRecord]:
        first, blocks = _first_text_line(source, source.blocks())
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

    def read() -> list[SeqRecord] | None:
        scan = _Scan(source, "fasta")
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
    sequence), so records a cache owns stay as they are.
    """

    detached: list[SeqRecord] = []
    for record in records:
        if _file_table(record.seq) is not None:
            record = copy.copy(record)
            record.seq = Seq(bytes(record.seq))
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
    return _with_encoding(source, lambda: read(_Scan(source, "fasta").gff3_lines()))


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
