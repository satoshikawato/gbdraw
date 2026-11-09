"""Generate reply buffers: chunked UTF-8 JSON with the bytes json.dumps writes."""

from __future__ import annotations

import json
import math
import os
import tempfile
import tracemalloc
from pathlib import Path

import pytest

from gbdraw.web_support.reply_buffers import (
    REPLY_CHUNK_CHARACTERS,
    json_utf8_chunks,
    text_utf8_chunks,
    write_reply_buffer,
)


class _Buffer:
    """A reply buffer as the Worker's: it takes each chunk out of its file."""

    def __init__(self) -> None:
        self.chunks: list[bytes] = []
        descriptor, self.path = tempfile.mkstemp(prefix="gbdraw-reply-")
        os.close(descriptor)

    def take(self) -> None:
        file = Path(self.path)
        self.chunks.append(file.read_bytes())
        file.unlink()


def _dumped(value: object) -> bytes:
    return json.dumps(value, ensure_ascii=False, separators=(",", ":")).encode("utf-8")


def _catalog(feature_count: int) -> dict[str, object]:
    feature = {
        "svg_id": "f0001",
        "qualifiers": {"product": ["µ-protein \"x\"   日本"], "note": []},
        "nucleotide_sequence": "ACGT" * 200,
        "location": [{"start": 1, "end": 800, "strand": -1}],
        "score": 0.5,
        "flags": (True, None),
    }
    return {
        "schema": 5,
        "items": [
            {
                "resultIndex": index,
                "resultName": f"r{index}.svg",
                "features": [dict(feature, svg_id=f"f{i}") for i in range(feature_count)],
                "biologicalFeatures": [],
                "anchorProfile": {"1": "a", "b": {}},
            }
            for index in range(2)
        ],
    }


@pytest.mark.parametrize(
    "value",
    [
        {},
        [],
        "µ",
        3,
        math.nan,
        {"a": {1: "int key", "b": [1, 2.5, -math.inf]}},
        [[[[["deep"]]]], ("tuple", "items"), {}],
        {"list": list(range(1000)), "nested": [{"x": [[i] * 3 for i in range(600)]}]},
        _catalog(10),
    ],
)
def test_chunks_join_to_the_json_dumps_bytes(value: object) -> None:
    assert b"".join(json_utf8_chunks(value)) == _dumped(value)


def test_a_large_catalog_is_written_in_chunks_of_about_one_mebibyte() -> None:
    catalog = _catalog(4000)
    buffer = _Buffer()

    assert write_reply_buffer(lambda: buffer, json_utf8_chunks(catalog)) is buffer

    assert b"".join(buffer.chunks) == _dumped(catalog)
    assert not os.path.exists(buffer.path)
    assert len(buffer.chunks) > 4
    # A chunk closes once it passes the size; one slice of features is the overshoot.
    assert max(len(chunk) for chunk in buffer.chunks) < 2 * REPLY_CHUNK_CHARACTERS


def test_writing_a_catalog_never_holds_its_whole_text_in_the_python_heap() -> None:
    catalog = _catalog(16000)
    size = len(_dumped(catalog))

    class _Discard:
        path = os.devnull

        def take(self) -> None:
            pass

    tracemalloc.start()
    try:
        write_reply_buffer(_Discard, json_utf8_chunks(catalog))
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()

    # json.dumps(...).encode() holds the text and its bytes at once: about twice
    # the size. The chunks hold a few chunks' worth, whatever the size.
    assert size > 24 * REPLY_CHUNK_CHARACTERS
    assert peak < 12 * REPLY_CHUNK_CHARACTERS


def test_a_long_text_is_cut_by_code_point_and_never_splits_a_character() -> None:
    # Characters of one to four UTF-8 bytes; a slice boundary falls inside the run.
    text = ("a" + "µ" + "日" + "\U0001F9EC") * (REPLY_CHUNK_CHARACTERS // 2 + 3)
    chunks = list(text_utf8_chunks(text))

    assert len(chunks) == 3
    assert b"".join(chunks) == text.encode("utf-8")
    assert all(chunk.decode("utf-8") for chunk in chunks)
    assert max(len(chunk.decode("utf-8")) for chunk in chunks) == REPLY_CHUNK_CHARACTERS
    assert list(text_utf8_chunks("")) == []


def test_a_json_string_longer_than_a_chunk_is_cut_too() -> None:
    value = {"items": [{"features": [{"nucleotide_sequence": "ACGT" * REPLY_CHUNK_CHARACTERS}]}]}
    chunks = list(json_utf8_chunks(value))

    assert b"".join(chunks) == _dumped(value)
    assert max(len(chunk) for chunk in chunks) <= 2 * REPLY_CHUNK_CHARACTERS
