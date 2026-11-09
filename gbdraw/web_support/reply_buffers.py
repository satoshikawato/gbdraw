"""The byte buffers of a Generate reply or a sequence read's reply, written in chunks.

A Generate reply carries UTF-8 buffers: the metadata JSON, the feature catalog
JSON, and each Result's sent children texts; a sequence read (the feature,
match and export reads) replies with its JSON in one buffer, since Pyodide
makes a JavaScript string of a long str with many times its size of Worker
heap. Python writes each one into a
buffer the Worker owns (``new_reply_buffer`` makes one) in chunks of about
1 MiB, so the Python heap, whose WebAssembly memory never shrinks, never holds a
whole text's bytes (or, for JSON, the whole text). Each chunk crosses through a
file: Python writes it to the buffer's file and the Worker takes it out of the
file system. A file write reads the chunk's heap address as unsigned, while
Pyodide's buffer views read it as signed and fail above 2 GiB (OV-2800). The
Worker joins each buffer's chunks after Python returns. Joined, the JSON chunks
are the bytes
``json.dumps(value, ensure_ascii=False, separators=(",", ":")).encode("utf-8")``
returns.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Iterable, Iterator
from typing import Protocol, TypeVar

# The characters gathered before they are encoded and written; also the slice
# length of a longer text.
REPLY_CHUNK_CHARACTERS = 1 << 20
# Dicts and lists above this depth are walked; deeper values, and lists at this
# depth in slices of _SLICE_ITEMS items, go to the C encoder whole. For the
# feature catalog the slices are runs of features: {"items": [{"features": [...]}]}.
_WALKED_DEPTH = 3
_SLICE_ITEMS = 256
_ENCODER = json.JSONEncoder(ensure_ascii=False, separators=(",", ":"))


class ReplyBuffer(Protocol):
    """One reply buffer the Worker owns: Python writes each chunk to the file at
    ``path``, then ``take`` moves it out of the file system into the buffer."""

    @property
    def path(self) -> str: ...

    def take(self) -> None: ...


_Buffer = TypeVar("_Buffer", bound=ReplyBuffer)


def write_reply_buffer(new_reply_buffer: Callable[[], _Buffer], chunks: Iterable[bytes]) -> _Buffer:
    """Write ``chunks`` in order into a new reply buffer and return the buffer.

    The Worker joins its chunks after Python returns (``finishReplyBuffers``).
    """

    buffer = new_reply_buffer()
    path = buffer.path
    for chunk in chunks:
        with open(path, "wb") as file:
            file.write(chunk)
        buffer.take()
    return buffer


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


def _utf8_chunks(texts: Iterable[str]) -> Iterator[bytes]:
    # Short texts are gathered up to REPLY_CHUNK_CHARACTERS; a longer one is cut
    # into slices of that many code points (a Python str slice never splits a
    # character). A chunk is at most twice REPLY_CHUNK_CHARACTERS characters.
    pending: list[str] = []
    characters = 0
    for text in texts:
        if len(text) > REPLY_CHUNK_CHARACTERS:
            if pending:
                yield "".join(pending).encode("utf-8")
                pending.clear()
                characters = 0
            for start in range(0, len(text), REPLY_CHUNK_CHARACTERS):
                yield text[start:start + REPLY_CHUNK_CHARACTERS].encode("utf-8")
            continue
        pending.append(text)
        characters += len(text)
        if characters >= REPLY_CHUNK_CHARACTERS:
            yield "".join(pending).encode("utf-8")
            pending.clear()
            characters = 0
    if pending:
        yield "".join(pending).encode("utf-8")


def json_utf8_chunks(value: object) -> Iterator[bytes]:
    """The compact UTF-8 JSON of ``value`` in chunks of about 1 MiB."""

    return _utf8_chunks(_json_texts(value, 0))


def text_utf8_chunks(text: str) -> Iterator[bytes]:
    """The UTF-8 of ``text`` in chunks of about 1 MiB (none for an empty text)."""

    return _utf8_chunks((text,) if text else ())
