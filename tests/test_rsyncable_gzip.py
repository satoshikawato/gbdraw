"""The rsyncable gzip writer every Session ``.gz`` goes through (gbdraw/io/rsyncable_gzip.py)."""

from __future__ import annotations

import gzip
import io
import random
import struct
import zlib

import pytest

from gbdraw.io import rsyncable_gzip
from gbdraw.io.rsyncable_gzip import MIN_CHUNK, WINDOW, RsyncableGzipWriter

# No file name, mtime 0, OS=255 on every Python version. gzip.compress(mtime=0) is not a
# reference for it: Python 3.11 and 3.12 write OS=3 there, 3.10 and 3.13+ write 255.
HEADER = b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\xff"


def _compress(pieces: list[bytes]) -> bytes:
    output = io.BytesIO()
    with RsyncableGzipWriter(output) as writer:
        for piece in pieces:
            writer.write(piece)
    return output.getvalue()


def _reference(data: bytes) -> bytes:
    """The cut rule written out position by position, as the measured prototype applies it."""

    deflate = zlib.compressobj(6, zlib.DEFLATED, -zlib.MAX_WBITS)
    parts = [HEADER]
    previous = 0
    window_sum = sum(data[:WINDOW])
    for end in range(WINDOW, len(data) + 1):
        if window_sum % WINDOW == 0 and end - previous >= MIN_CHUNK:
            parts += [deflate.compress(data[previous:end]), deflate.flush(zlib.Z_FULL_FLUSH)]
            previous = end
        if end < len(data):
            window_sum += data[end] - data[end - WINDOW]
    parts += [deflate.compress(data[previous:]), deflate.flush()]
    parts.append(struct.pack("<II", zlib.crc32(data), len(data) & 0xFFFFFFFF))
    return b"".join(parts)


def _session_like_text(size: int, seed: int) -> bytes:
    """JSON-like text with the byte mix of a Session: keys, numbers and base64."""

    rng = random.Random(seed)
    alphabet = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
    parts: list[bytes] = []
    total = 0
    while total < size:
        number = rng.randrange(10**6)
        blob = bytes(rng.choices(alphabet, k=rng.randrange(20, 400)))
        part = b'{"start":%d,"end":%d,"data":"%s"},' % (number, number + 900, blob)
        parts.append(part)
        total += len(part)
    return b"".join(parts)[:size]


@pytest.mark.parametrize(
    "data",
    [b"", b"{}", b"x" * (MIN_CHUNK - 1), _session_like_text(700_000, 1)],
    ids=["empty", "tiny", "below-min-chunk", "several-chunks"],
)
def test_output_is_gzip_that_expands_to_the_input(data: bytes) -> None:
    compressed = _compress([data])

    assert gzip.decompress(compressed) == data
    assert compressed[:10] == HEADER


def test_cuts_follow_the_rule() -> None:
    zeros = bytes(3 * MIN_CHUNK + 5)  # every window sums to 0: a cut every MIN_CHUNK bytes
    text = _session_like_text(700_000, 5)

    assert _compress([zeros]) == _reference(zeros)
    assert _compress([zeros]).count(b"\x00\x00\xff\xff") >= 3  # three full flushes
    assert _compress([text]) == _reference(text)


def test_equal_input_gives_equal_bytes_however_the_writes_split_it() -> None:
    data = _session_like_text(700_000, 2)
    rng = random.Random(3)
    pieces = []
    offset = 0
    while offset < len(data):
        size = rng.choice([1, 7, 4095, 8192, MIN_CHUNK + 1, 300_000])
        pieces.append(data[offset : offset + size])
        offset += size

    assert _compress(pieces) == _compress([data])


def test_bytes_after_a_local_edit_compress_to_the_same_bytes() -> None:
    data = _session_like_text(2_000_000, 4)
    edited = data[:100_000] + b'"inserted":true,' + data[100_000:]
    before = _compress([data])
    after = _compress([edited])

    shared_tail = 0
    while shared_tail < min(len(before), len(after)) - 8 and (
        before[-9 - shared_tail] == after[-9 - shared_tail]
    ):
        shared_tail += 1
    # Everything after the first cut past the edit is identical; the gzip
    # trailer (CRC-32 and size, 8 bytes) is not compared. Plain gzip would share
    # almost nothing after the edit.
    assert shared_tail > 0.8 * len(before)


def test_a_run_with_no_cut_is_compressed_as_it_arrives(monkeypatch: pytest.MonkeyPatch) -> None:
    # Any 4096-byte window of "AB" sums to 2048 * 131, never a multiple of 4096;
    # a cut follows in the zeros after it.
    monkeypatch.setattr(rsyncable_gzip, "_MAX_PENDING", 200_000)
    data = b"AB" * 300_000 + bytes(2 * MIN_CHUNK) + b"AB" * 50_000
    output = io.BytesIO()
    writer = RsyncableGzipWriter(output)
    for offset in range(0, len(data), 8192):
        writer.write(data[offset : offset + 8192])
        assert len(writer._pending) <= 200_000 + 8192
    writer.close()

    assert output.getvalue() == _reference(data)
