"""Gzip output whose compressed bytes resynchronize after a local change ("rsyncable").

The compressor resets its dictionary (a deflate full flush) at cut points chosen
from the uncompressed bytes alone, as ``gzip --rsyncable`` does: a cut follows the
first byte, at least ``MIN_CHUNK`` bytes after the previous cut, where the sum of
the last ``WINDOW`` bytes is a multiple of ``WINDOW``. Bytes between two cuts that a
change leaves alone compress to the same bytes, so git stores a regenerated file
as a small delta of its previous version. The stream stays plain gzip (RFC 1952),
with no file name and mtime 0: equal input gives equal bytes, however the writes
split it, and every gzip reader expands it.
"""

from __future__ import annotations

import io
import struct
import zlib
from typing import TYPE_CHECKING, BinaryIO

if TYPE_CHECKING:
    from _typeshed import ReadableBuffer

WINDOW = 4096
MIN_CHUNK = 64 * 1024
# A run with no cut is compressed as it arrives once it is this long, so memory
# stays bounded; only the last WINDOW bytes are kept for the rolling sum.
_MAX_PENDING = 1024 * 1024
# ID1 ID2 CM=deflate FLG=0 MTIME=0 XFL=0 OS=255 (unknown), as gzip.compress(mtime=0) writes at level 6.
_HEADER = b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff"


class RsyncableGzipWriter(io.BufferedIOBase):
    """A binary stream that writes rsyncable gzip to ``fileobj``; closing it ends the
    gzip stream and leaves ``fileobj`` open."""

    def __init__(self, fileobj: BinaryIO, *, compresslevel: int = 6) -> None:
        super().__init__()
        self._fileobj = fileobj
        self._deflate = zlib.compressobj(compresslevel, zlib.DEFLATED, -zlib.MAX_WBITS)
        self._crc = 0
        self._size = 0
        # The bytes since the last cut that are not compressed yet.
        self._pending = bytearray()
        # The last candidate cut tested in ``_pending`` and the sum of the WINDOW
        # bytes before it; None until ``_pending`` holds MIN_CHUNK bytes.
        self._tested: int | None = None
        self._sum = 0
        fileobj.write(_HEADER)

    @property
    def name(self) -> str:
        """No file name, as the gzip header records none (``io.TextIOWrapper`` reads it)."""

        return ""

    def writable(self) -> bool:
        return True

    def write(self, data: ReadableBuffer) -> int:
        if self.closed:
            raise ValueError("write to a closed RsyncableGzipWriter")
        view = memoryview(data).cast("B")
        self._crc = zlib.crc32(view, self._crc)
        self._size += len(view)
        # In pieces, so one large write is not copied whole into ``_pending``.
        for start in range(0, len(view), _MAX_PENDING):
            self._pending += view[start : start + _MAX_PENDING]
            self._cut_pending()
        return len(view)

    def _cut_pending(self) -> None:
        pending = self._pending
        while True:
            if self._tested is None:
                if len(pending) < MIN_CHUNK:
                    return
                position = MIN_CHUNK - 1
                window_sum = sum(pending[position - WINDOW : position])
            else:
                position, window_sum = self._tested, self._sum
            end = len(pending)
            cut = None
            with memoryview(pending) as view:
                for added, dropped in zip(view[position:end], view[position - WINDOW : end - WINDOW]):
                    window_sum += added - dropped
                    position += 1
                    if not window_sum % WINDOW:
                        cut = position
                        break
            if cut is None:
                if end > _MAX_PENDING:
                    done = end - WINDOW
                    self._fileobj.write(self._deflate.compress(pending[:done]))
                    del pending[:done]
                    position -= done
                self._tested, self._sum = position, window_sum
                return
            self._fileobj.write(self._deflate.compress(pending[:cut]))
            self._fileobj.write(self._deflate.flush(zlib.Z_FULL_FLUSH))
            del pending[:cut]
            self._tested = None

    def close(self) -> None:
        if self.closed:
            return
        try:
            self._fileobj.write(self._deflate.compress(self._pending))
            self._fileobj.write(self._deflate.flush(zlib.Z_FINISH))
            self._fileobj.write(struct.pack("<II", self._crc, self._size & 0xFFFFFFFF))
            self._pending.clear()
        finally:
            super().close()
