"""SVG document text: the one serialization of a drawing per request."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator, TextIO, TypeAlias

from svgwrite import Drawing
from svgwrite.base import BaseElement

from gbdraw.render.group_content import (
    GroupContent,
    element_text,
    placed_content,
    require_element_markup,
    start_tag,
)


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
        write_text(handle, self.body)


# The characters handed to a text file at a time: its encoder then holds a
# chunk's bytes, not a second copy of the whole document.
_WRITE_CHUNK = 1 << 20


def write_text(handle: TextIO, text: str) -> None:
    """Write ``text`` to an open text file a chunk at a time; the file gets the same bytes."""

    for start in range(0, len(text), _WRITE_CHUNK):
        handle.write(text[start : start + _WRITE_CHUNK])


_DECLARATION = '<?xml version="1.0" encoding="utf-8" ?>\n'
_STYLESHEET = (
    '<?xml-stylesheet href="%s" type="text/css" title="%s" alternate="%s" media="%s"?>\n'
)


@dataclass(frozen=True)
class ContentPiece:
    """A placed content group with children: its start tag and its content's children.

    Its end tag is always ``</g>``. ``revision`` and ``byte_length`` name the
    children text (its sha256 and UTF-8 length), which the Web reply sends only
    when the page does not hold that revision.
    """

    start: str
    children: str
    content: GroupContent = field(repr=False, compare=False)

    @property
    def revision(self) -> str:
        return self.content.revision

    @property
    def byte_length(self) -> int:
        return self.content.children_byte_length


@dataclass(frozen=True)
class RawPiece:
    """Any other top-level element, by svgwrite's serialization."""

    text: str


Piece: TypeAlias = ContentPiece | RawPiece


@dataclass(frozen=True)
class DocumentParts:
    """A drawing's text as ordered pieces: ``prolog + head + pieces + tail``.

    ``head`` is the root start tag and ``tail`` its end tag; ``pieces`` are the
    top-level elements in document order. Every text but the prolog holds only
    whole elements with escaped text and no ``<!`` or ``<?`` (precondition P).
    """

    prolog: str
    head: str
    pieces: tuple[Piece, ...]
    tail: str

    def body(self) -> str:
        """The root element: head, pieces and tail joined in order."""

        return "".join((self.head, *_piece_texts(self.pieces), self.tail))


def _piece_texts(pieces: tuple[Piece, ...]) -> Iterator[str]:
    for piece in pieces:
        if isinstance(piece, ContentPiece):
            yield piece.start
            yield piece.children
            yield "</g>"
        else:
            yield piece.text


def _piece(element: BaseElement) -> Piece:
    content = placed_content(element)
    if content is not None:
        return ContentPiece(
            start=require_element_markup(start_tag(element)),
            children=content.children_text,
            content=content,
        )
    # defs, a group without content, and an empty group (``<g ... />``).
    return RawPiece(require_element_markup(element_text(element)))


def document_parts(drawing: Drawing) -> DocumentParts:
    """Serialize ``drawing`` once, as ordered pieces (the one order of serialization).

    The prolog is svgwrite's: the XML declaration and the stylesheet
    instructions. A placed content group with children is a ``ContentPiece``
    over the stored text of its children; every other top-level element is a
    ``RawPiece`` by svgwrite's own serialization. Raises ``ExportError`` when a
    piece holds a comment, CDATA or a processing instruction.
    """

    prolog = _DECLARATION + "".join(
        _STYLESHEET % stylesheet for stylesheet in drawing._stylesheets
    )
    if not drawing.elements:
        return DocumentParts(
            prolog=prolog, head=require_element_markup(drawing.tostring()), pieces=(), tail=""
        )
    return DocumentParts(
        prolog=prolog,
        head=require_element_markup(start_tag(drawing)),
        pieces=tuple(_piece(element) for element in drawing.elements),
        tail=f"</{drawing.elementname}>",
    )


def write_svg_document(drawing: Drawing) -> SvgDocument:
    """Serialize ``drawing`` once, as ``Drawing.write`` would write it (``document_parts``)."""

    parts = document_parts(drawing)
    return SvgDocument(prolog=parts.prolog, body=parts.body())


__all__ = [
    "ContentPiece",
    "DocumentParts",
    "Piece",
    "RawPiece",
    "SvgDocument",
    "document_parts",
    "write_svg_document",
    "write_text",
]
