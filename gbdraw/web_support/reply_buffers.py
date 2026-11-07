"""The byte buffers of a Generate reply, written in chunks.

A Generate reply carries UTF-8 buffers: the metadata JSON, the feature catalog
JSON, and each Result's sent children texts. Python writes each one into a
buffer the Worker owns (``new_reply_buffer`` makes one), in chunks of about
1 MiB, so neither the whole text nor its whole bytes exist in the Python heap,
whose WebAssembly memory never shrinks. Joined, the JSON chunks are the bytes
``json.dumps(value, ensure_ascii=False, separators=(",", ":")).encode("utf-8")``
returns.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Iterable, Iterator
from typing import Protocol, TypeVar

# The characters of JSON text gathered before they are encoded and written.
REPLY_CHUNK_CHARACTERS = 1 << 20
# Dicts and lists above this depth are walked; deeper values, and lists at this
# depth in slices of _SLICE_ITEMS items, go to the C encoder whole. For the
# feature catalog the slices are runs of features: {"items": [{"features": [...]}]}.
_WALKED_DEPTH = 3
_SLICE_ITEMS = 256
_ENCODER = json.JSONEncoder(ensure_ascii=False, separators=(",", ":"))

_Finished = TypeVar("_Finished", covariant=True)


class ReplyBuffer(Protocol[_Finished]):
    """One reply buffer: the Worker copies each chunk out of the Python heap."""

    def write(self, chunk: bytes) -> None: ...

    def finish(self) -> _Finished: ...


_Buffer = TypeVar("_Buffer")


def write_reply_buffer(
    new_reply_buffer: Callable[[], ReplyBuffer[_Buffer]],
    chunks: Iterable[bytes],
) -> _Buffer:
    """Write ``chunks`` in order into a new reply buffer and return what it holds."""

    buffer = new_reply_buffer()
    for chunk in chunks:
        buffer.write(chunk)
    return buffer.finish()


def _json_texts(value: object, depth: int) -> Iterator[str]:
    encode = _ENCODER.encode
    if type(value) is dict and value and depth < _WALKED_DEPTH and all(type(key) is str for key in value):
        opener = "{"
        for key, member in value.items():
            yield opener + encode(key) + ":"
            opener = ","
            yield from _json_texts(member, depth + 1)
        yield "}"
    elif type(value) is list and value and depth < _WALKED_DEPTH:
        opener = "["
        for member in value:
            yield opener
            opener = ","
            yield from _json_texts(member, depth + 1)
        yield "]"
    elif type(value) is list and value:
        # A slice encodes as its items joined by "," inside brackets.
        opener = "["
        for start in range(0, len(value), _SLICE_ITEMS):
            yield opener
            yield encode(value[start:start + _SLICE_ITEMS])[1:-1]
            opener = ","
        yield "]"
    else:
        yield encode(value)


def json_utf8_chunks(value: object) -> Iterator[bytes]:
    """The compact UTF-8 JSON of ``value`` in chunks of about 1 MiB."""

    pending: list[str] = []
    characters = 0
    for text in _json_texts(value, 0):
        pending.append(text)
        characters += len(text)
        if characters >= REPLY_CHUNK_CHARACTERS:
            yield "".join(pending).encode("utf-8")
            pending.clear()
            characters = 0
    if pending:
        yield "".join(pending).encode("utf-8")
