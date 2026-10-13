"""SVG document text: the one serialization of a drawing per request."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TextIO

from svgwrite import Drawing


@dataclass(frozen=True)
class SvgDocument:
    """The text of one serialized drawing.

    ``prolog`` is what svgwrite writes before the root element of a file: the
    XML declaration and any stylesheet instructions. ``body`` is the root
    element alone; in-memory SVG, interactive SVG and the CairoSVG input use it
    without the prolog, as ``Drawing.tostring()`` gave them.
    """

    prolog: str
    body: str

    def write(self, handle: TextIO) -> None:
        """Write the file form, prolog then body, to an open text file."""

        handle.write(self.prolog)
        handle.write(self.body)


class _TextParts:
    """A text sink that keeps each write svgwrite makes as one part."""

    def __init__(self) -> None:
        self.parts: list[str] = []

    def write(self, text: str) -> int:
        self.parts.append(text)
        return len(text)


def write_svg_document(drawing: Drawing) -> SvgDocument:
    """Serialize ``drawing`` once with svgwrite's own file writer.

    ``Drawing.write`` writes the prolog and then the root element as its last
    write, so the file text and the body come from the same serialization.
    """

    sink = _TextParts()
    drawing.write(sink)
    *prolog, body = sink.parts
    return SvgDocument(prolog="".join(prolog), body=body)


__all__ = ["SvgDocument", "write_svg_document"]
