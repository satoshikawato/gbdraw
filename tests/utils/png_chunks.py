"""Read the chunks of a PNG, for tests of its size and resolution."""

from __future__ import annotations

import struct
import zlib

PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def png_chunks(payload: bytes) -> list[tuple[bytes, bytes]]:
    """The ``(type, data)`` of every chunk, each CRC checked."""

    assert payload[:8] == PNG_SIGNATURE
    chunks = []
    offset = 8
    while offset < len(payload):
        (length,) = struct.unpack(">I", payload[offset : offset + 4])
        kind = payload[offset + 4 : offset + 8]
        data = payload[offset + 8 : offset + 8 + length]
        (crc,) = struct.unpack(">I", payload[offset + 8 + length : offset + 12 + length])
        assert crc == zlib.crc32(kind + data), kind
        chunks.append((kind, data))
        offset += 12 + length
    return chunks


def png_size(payload: bytes) -> tuple[int, int]:
    """The pixel width and height from the IHDR chunk."""

    return struct.unpack(">II", png_chunks(payload)[0][1][:8])


def png_phys(payload: bytes) -> tuple[int, int, int] | None:
    """``(pixels per metre x, y, unit)`` of the one pHYs chunk, before the first IDAT; ``None`` without one."""

    kinds = [kind for kind, _data in png_chunks(payload)]
    assert kinds.count(b"pHYs") <= 1
    if b"pHYs" not in kinds:
        return None
    assert kinds.index(b"pHYs") < kinds.index(b"IDAT")
    data = dict(png_chunks(payload))[b"pHYs"]
    return struct.unpack(">IIB", data)


def tiny_png() -> bytes:
    """A valid 1 x 1 PNG without a pHYs chunk, for fakes of CairoSVG's ``svg2png``."""

    def chunk(kind: bytes, data: bytes) -> bytes:
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data))

    return (
        PNG_SIGNATURE
        + chunk(b"IHDR", struct.pack(">IIBBBBB", 1, 1, 8, 2, 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(b"\x00\x00\x00\x00"))
        + chunk(b"IEND", b"")
    )
