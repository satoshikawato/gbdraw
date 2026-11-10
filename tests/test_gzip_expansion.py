"""The one gzip expansion (:class:`gbdraw.io.input_text.GzipExpansion`) reads as Python's ``gzip`` module reads.

Every gzip input is expanded there: the text readers through ``open_input_text`` and
``read_input_bytes``, the record sequence reader in its parse pass and at offsets.
So it gives the module's bytes and refuses what the module refuses, with the same
exception and message, and it starts again at any point it gives.
"""

from __future__ import annotations

import gzip
import io
import random
import struct
import zlib

import pytest

from gbdraw.io.input_text import GZIP_START, GzipExpansion, read_input_bytes


def _member(data: bytes, *, flags: int = 0, extra: bytes = b"", name: bytes = b"", comment: bytes = b"",
            crc: int | None = None, size: int | None = None, method: int = 8) -> bytes:
    header = struct.pack("<BBBBIBB", 0x1F, 0x8B, method, flags, 0, 0, 255)
    if flags & 4:
        header += struct.pack("<H", len(extra)) + extra
    if flags & 8:
        header += name + b"\0"
    if flags & 16:
        header += comment + b"\0"
    if flags & 2:
        header += b"\x12\x34"  # a header CRC: read past, never checked
    compressor = zlib.compressobj(6, zlib.DEFLATED, -zlib.MAX_WBITS)
    body = compressor.compress(data) + compressor.flush()
    trailer = struct.pack("<II", zlib.crc32(data) if crc is None else crc, (len(data) if size is None else size) & 0xFFFFFFFF)
    return header + body + trailer


_RANDOM = random.Random(7)
_TEXT = b"".join(b">r%d\n" % i + bytes(_RANDOM.choices(b"ACGTN", k=3000)) + b"\n" for i in range(600))
_RUNS = b"N" * 2_000_000 + _TEXT  # a run that one compressed chunk expands past any read size
_SHORT = _TEXT[:20000]
_ONE = gzip.compress(_SHORT)

CASES = {
    "one member": gzip.compress(_TEXT),
    "runs": gzip.compress(_RUNS),
    "empty member": gzip.compress(b""),
    "members": _member(_TEXT[:5000]) + _member(b"") + _member(_TEXT[5000:]),
    "bgzip-like extra field": b"".join(
        _member(part, flags=4, extra=b"BC\x02\x00\x00\x00") for part in (_TEXT[:70000], _TEXT[70000:], b"")),
    "name, comment and header crc": _member(_SHORT, flags=8 | 16 | 2, name=b"a.fa", comment=b"note"),
    "zero padding": _ONE + b"\0" * 100,
    "padding then a member": _ONE + b"\0" * 3 + _member(b"ACGT"),
    "trailing garbage": _ONE + b"garbage",
    "one trailing byte": _ONE + b"x",
    "cut in the stream": _ONE[: len(_ONE) // 2],
    "cut in the trailer": _ONE[:-3],
    "cut in a header": _ONE + _ONE[:6],
    "cut in a name": _ONE + _member(b"A", flags=8, name=b"long name")[:13],
    "wrong crc": _member(_SHORT, crc=1),
    "wrong length": _member(_SHORT, size=3),
    "unknown method": _member(_SHORT, method=7),
    "damaged stream": b"\x1f\x8b\x08\x00" + bytes(6) + bytes(range(128, 256)),
}


def _outcome(read):
    try:
        return read()
    except (EOFError, zlib.error, gzip.BadGzipFile) as error:
        return type(error), str(error)


def _module(data: bytes):
    return _outcome(lambda: gzip.GzipFile(fileobj=io.BytesIO(data)).read())


def _expanded(data: bytes, size: int, start=GZIP_START):
    handle = io.BytesIO(data)
    handle.seek(start.compressed)
    expansion = GzipExpansion(handle.read, "input.gz", start=start)

    def read() -> bytes:
        pieces = []
        while piece := expansion.read(size):
            assert len(piece) <= size
            pieces.append(piece)
        return b"".join(pieces)

    return _outcome(read)


@pytest.mark.parametrize("case", sorted(CASES))
@pytest.mark.parametrize("size", [1, 4093, 1 << 20])
def test_the_expansion_gives_and_refuses_what_the_gzip_module_does(case: str, size: int) -> None:
    data = CASES[case]
    if size == 1 and len(data) > 50_000:
        pytest.skip("byte at a time over a long input: the short inputs read so")
    assert _expanded(data, size) == _module(data)


@pytest.mark.parametrize("case", ["one member", "runs", "members", "bgzip-like extra field"])
def test_reading_starts_again_at_every_point_the_expansion_gives(case: str) -> None:
    data = CASES[case]
    expected = gzip.decompress(data)
    handle = io.BytesIO(data)
    expansion = GzipExpansion(handle.read, "input.gz")
    points = [GZIP_START]
    while expansion.read(_RANDOM.choice((97, 65536, 1 << 20, 1 << 20))):
        if (point := expansion.point()) is not None:
            points.append(point)
    assert len(points) > 3
    for point in points:
        assert _expanded(data, 1 << 20, point) == expected[point.expanded:]
    restarted = _expanded(data, 1 << 20, points[2])
    assert _expanded(data, 1 << 20, points[2]) == restarted  # a point stays: each start reads from a copy of it


def test_the_text_readers_expand_through_it(tmp_path) -> None:
    for case in ("members", "trailing garbage", "wrong crc"):
        path = tmp_path / "input.gz"
        path.write_bytes(CASES[case])
        assert _outcome(lambda: read_input_bytes(path)) == _module(CASES[case])
