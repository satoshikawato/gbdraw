"""The SVG text of a Web render Result, joined from its reply pieces (tests only).

The Web reply carries a drawing's text as pieces and never joins it (OV-3027); the page assembles them. A test that
reads the whole text joins a Result whose every content piece carries its children text (no held revision).
"""

from __future__ import annotations

from collections.abc import Mapping


def web_result_svg(result: Mapping[str, object]) -> str:
    """The root element of ``result`` (a ``WebRenderResult``): head, pieces in order, tail."""

    document = result["document"]
    assert isinstance(document, Mapping)
    texts: list[str] = [str(document["head"])]
    for piece in document["pieces"]:  # type: ignore[attr-defined]
        if piece["kind"] == "content":
            assert "children" in piece, "a held revision carries no text to join"
            texts.extend((piece["start"], piece["children"], "</g>"))
        else:
            texts.append(piece["text"])
    texts.append(str(document["tail"]))
    return "".join(texts)
